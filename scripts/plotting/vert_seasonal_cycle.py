
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

from datetime import datetime
from dateutil.relativedelta import relativedelta
import pandas as pd

import xarray as xr
from pathlib import Path
import glob
from itertools import chain



#Set seasonal ranges:
seasons = {"DJF": [12, 1, 2],
            "JJA": [6, 7, 8],
            "MAM": [3, 4, 5],
            "SON": [9, 10, 11]}

#Set monthly codes:
month_dict = {1:'JAN',
    		  2:'FEB',
    		  3:'MAR',
    		  4:'APR',
    		  5:'MAY',
    		  6:'JUN',
    		  7:'JUL',
    		  8:'AUG',
    		  9:'SEP',
    		  10:'OCT',
    		  11:'NOV',
    		  12:'DEC'}

delta_symbol = r'$\Delta$'

temp_levs = np.arange(140, 300, 10)
temp_diff_levs = np.arange(-40, 41, 4)

#contour_levels = np.arange(-30, 31, 3)
#orig_contour_levels = np.arange(-120, 121, 10)
wind_levs = np.arange(-120, 121, 10)
wind_diff_levs = np.arange(-30, 31, 3)
cont_ranges = {"U":{"levs":wind_levs,"diff_levs":wind_diff_levs,"units":"m/s"},
               "T":{"levs":temp_levs,"diff_levs":temp_diff_levs,"units":"K"}}


calc_var_list = ['U','T','V','UU','VU','VT','OMEGAU','OMEGA','O3','Q','UTEND_VDIFF']
#calc_var_list = ['U','T','V','O3','Q','UTEND_VDIFF']
calc_var_list = ['U','T']

merra2_vars = ['U','T','V','lat','lev']


obs_cam_vars={"saber":{"U":"u", "T":"temp"},
              "merra":{"U":"U", "T":"T"}}

#
# --- Main Function Shares Name with Module: regional_map_multicase ---
#
def vert_seasonal_cycle(adfobj):
    """
    Chemistry Map main function
        * Initially start with  Zonal maps
            - This probably can be expanded to LatLon if given single pressure levels?

        
    """

    #CAM simulation variables (this is always assumed to be a list):
    case_names = adfobj.get_cam_info("cam_case_name", required=True)
    #Extract cam history files location:
    cam_hist_locs = adfobj.get_cam_info('cam_hist_loc')

    #Special ADF variable which contains the output paths for
    #all generated plots and tables:
    #plot_locations = adfobj.cam_diag_plot_loc

    #Grab case years
    syear_cases = adfobj.climo_yrs["syears"]
    eyear_cases = adfobj.climo_yrs["eyears"]

    if not adfobj.get_basic_info("compare_obs"):
        obs = False
        data_name = adfobj.get_baseline_info("cam_case_name", required=True) # does not get used, is just here as a placemarker
        data_list = [data_name] # gets used as just the name to search for climo files HAS TO BE LIST
    

        #Grab baseline years (which may be empty strings if using Obs):
        syear_baseline = adfobj.climo_yrs["syear_baseline"]
        syear_cases = syear_cases + [syear_baseline]
        eyear_baseline = adfobj.climo_yrs["eyear_baseline"]
        eyear_cases = eyear_cases + [eyear_baseline]

        #Grab all case nickname(s)
        test_nicknames = adfobj.case_nicknames["test_nicknames"]
        base_nickname = adfobj.case_nicknames["base_nickname"]

        case_names = case_names + data_list

        #Get climo years for verification or assignment if missing
        baseline_hist_locs = adfobj.get_baseline_info('cam_hist_loc')
        cam_hist_locs = cam_hist_locs + [baseline_hist_locs]
    #End if

    # Notify user that script has started:
    print("\n  Generating zonal vertical seasonal cycle plots plots ...")


    #var_list = adfobj.diag_var_list
    var_list = calc_var_list + ['lat','lev','time']

    #case_names = case_deets["case_names"]["cases"] + case_deets["case_names"]["baseline"]
    #runs = [adfobj.get_cam_info("cam_case_name")[0],adfobj.get_cam_baseline_info("cam_case_name")[0]+"_fake",

    #runs = [adfobj.get_cam_info("cam_case_name")[0],adfobj.get_cam_baseline_info("cam_case_name"),
       #adfobj.get_cam_info("cam_case_name")[0]+"_fake2",
       #adfobj.get_cam_info("cam_case_name")[0]+"_fake3",
       #adfobj.get_cam_info("cam_case_name")[0]+"_fake4",
       #adfobj.get_cam_info("cam_case_name")[0]+"_fake5",
       #adfobj.get_cam_info("cam_case_name")[0]+"_fake6",
       #adfobj.get_cam_info("cam_case_name")[0]+"_fake7",
     #  ]
    


    


    #Get MERRA2 data and seasonal and monthly averages
    
    cases_coords = {}
    cases_seasonal = {}
    cases_monthly = {}
    for idx,case_name in enumerate(case_names):
        #run = run[0]
        print(case_name)
        #zmdir = adfobj.basic_info_dict['diag_loc'] + f"zm/{run}/" #+ "/atm/hist/"
        #file = rundir + "waccm_135.nc"
        #if idx == 0:
        hist_loc = cam_hist_locs[idx]
        #case_name = case[idx]
        syr = syear_cases[idx]
        eyr = eyear_cases[idx]

        #make_zm_files(hist_loc,case_name,calc_var_list,syr,eyr,return_ds=True):
        #ncfile = make_zm_files(adfobj,hist_loc,case_name,calc_var_list,syr,eyr,return_ds=True)
        file = "/glade/work/richling/ADF/ADF_dev/notebooks/chem-diags/new_tests/waccm_135_acom_ne16pg3_ne16pg3_mg17_1536_long2.nc"
        ncfile = xr.open_dataset(file, decode_times=True, use_cftime=True)

        #zmdir = run
        #file = f"waccm_135_{zmdir}.nc"
        #print(file)
        #ncfile = nc.Dataset(file, 'r')
        #ncfile = xr.open_dataset(file, decode_times=True, use_cftime=True)

        #Check if plot output directory exists, and if not, then create it:
        #if not plot_loc.is_dir():
        #    print("    {} not found, making new directory".format(plot_loc))
        #    plot_loc.mkdir(parents=True)

        case_coords = {}
        case_seasonal = {}
        case_monthly = {}
        for index, var in enumerate(var_list):

            if var not in case_seasonal:
                case_seasonal[var] = {}
            if var not in case_monthly:
                case_monthly[var] = {}
            case_coords[var] =  ncfile[var]

            #TODO: clean this up,
            if var in calc_var_list:
                for season in seasons:
                    if season not in case_seasonal[var]:
                        case_seasonal[var][season] = time_mean(ncfile, case_coords[var],
                                                               time_avg="season",
                                                               interval=season,
                                                               is_climo=None)

                #Set months number to reflect actual month num, ie 1:Jan, 2:Feb, etc
                for month in np.arange(1,13,1):
                    case_monthly[var][month_dict[month]] = time_mean(ncfile, case_coords[var],
                                                                    time_avg="month",
                                                                    interval=month,
                                                                    is_climo=None)

        cases_coords[case_name] = case_coords
        cases_monthly[case_name] = case_monthly
        cases_seasonal[case_name] = case_seasonal
    
    #Make nested dictionary of all case data
    case_ds_dict = {"coords":cases_coords,
                    "monthly":cases_monthly,
                    "seasonal":cases_seasonal}


    #Get Obs and seasonal and monthly averages
    saber_monthly, saber_seasonal = saber_data(filename = "/glade/work/richling/ADF/ADF_dev/notebooks/chem-diags/SABER_monthly_2002-2014.nc")
    merra2_monthly, merra2_seasonal = merra_data(filename = "/glade/work/richling/ADF/ADF_dev/notebooks/chem-diags/MERRA2_met.nc")

    obs_seas_dict = {"saber":saber_seasonal, "merra":merra2_seasonal}
    obs_month_dict = {"saber":saber_monthly, "merra":merra2_monthly}
    obs_ds_dict = {"monthly":obs_month_dict,
                   "seasonal":obs_seas_dict}
    
    
    for cam_var in calc_var_list:
        for month in [6,12]:
            comparison_plots(adfobj, cam_var, case_names, case_ds_dict, obs_ds_dict, "month", month)
        for season in ["DJF", "JJA"]:
            comparison_plots(adfobj, cam_var, case_names, case_ds_dict, obs_ds_dict, "season", season)
            #comparison_plots(adfobj, cam_var, case_names, case_ds_dict, obs_ds_dict, time_avg, interval):






"""

case_deets = {"years":{"syears":syear_cases,"eyears":eyear_cases,
                                  "syear_baseline":syear_baseline,"eyear_baseline":eyear_baseline},
                 "nicknames":{"cases":test_nicknames,
                             "baseline":base_nickname},
                 "case_names":{"cases":case_names,
                               "baseline":data_name},
                 "ptype":plot_type
                }
"""

def comparison_plots(adfobj, cam_var, case_names, case_ds_dict, obs_ds_dict, time_avg, interval):
    """

    """

    #Get plotting details for variable
    levs = cont_ranges[cam_var]["levs"]
    diff_levs = cont_ranges[cam_var]["diff_levs"]
    units = cont_ranges[cam_var]["units"]

    #Grab obs variable corresponding to CAM variable
    saber_var = obs_cam_vars['saber'][cam_var]
    merra_var = obs_cam_vars['merra'][cam_var]

    font_size = 8

    #Get number of test cases (number of columns)
    casenum = len(case_names)

    #Number of obs to compare
    #Currently, just compared to MERRA2 and SABER
    obsnum = 2
    nrows = obsnum+1

    #Set up plot
    fig = plt.figure(figsize=(casenum*4,nrows*5))

    for idx,case_name in enumerate(case_names):
        """
        cases_coords[run] = case_coords
        cases_monthly[run] = case_monthly
        cases_seasonal[run] = case_seasonal
    
    #Make nested dictionary of all case data
    case_ds_dict = {"coords":cases_coords,
                    "monthly":cases_monthly,
                    "seasonal":cases_seasonal}
        """
        data_coords = case_ds_dict["coords"][case_name]
        #data_coords = case_runs[case_name]
        data_lev = data_coords['lev']
        data_lat = data_coords['lat']

        #Set lat/lev grid for plotting
        [lat_grid, lev_grid] = np.meshgrid(data_lev,data_lat)

        if time_avg == "season":
            #case_ds_dict["seasonal"]
            data_array = case_ds_dict["seasonal"][case_name][cam_var][interval]
            #data_array = case_runs_seasonal[case_name][cam_var][interval]

            #Make Obs interpolated field from case
            #merra_ds = merra2_seasonal[merra_var][interval]
            merra_ds = obs_ds_dict["seasonal"]["merra"][merra_var][interval]
            merra_rfield = merra_ds.interp(lat=data_lat, lev=data_lev, method='linear')

            #saber_ds = saber_seasonal[saber_var][interval]
            saber_ds = obs_ds_dict["seasonal"]["saber"][saber_var][interval]
            saber_rfield = saber_ds.interp(lat=data_lat, lev=data_lev, method='linear')

        if time_avg == "month":
            case_ds_dict["monthly"]
            str_month = month_dict[interval]
            data_array = case_ds_dict["monthly"][case_name][cam_var][str_month]
            #data_array = case_runs_monthly[case_name][cam_var][month_dict[interval]]

            #Make Obs interpolated fields from case
            #merra_ds = merra2_monthly[merra_var][str_month]
            merra_ds = obs_ds_dict["monthly"]["merra"][merra_var][str_month]
            merra_rfield = merra_ds.interp(lat=data_lat, lev=data_lev, method='linear')

            #saber_ds = saber_monthly[saber_var][str_month]
            saber_ds = obs_ds_dict["monthly"]["saber"][saber_var][str_month]
            saber_rfield = saber_ds.interp(lat=data_lat, lev=data_lev, method='linear')

        #Case plots (contours and contour fill)
        #######################################

        #Set up set of axes for first row
        ax = fig.add_subplot(nrows, casenum, idx+1)

        #Plot case contour fill
        cf=plt.contourf(lev_grid, lat_grid,
                        data_array.transpose(transpose_coords=True),
                        levels=levs, cmap='RdBu_r')

        #Plot case contours (for highlighting)
        plt.contour(lev_grid, lat_grid,
                        data_array.transpose(transpose_coords=True),
                    colors="black",linewidths=0.5,levels=levs,zorder=100)

        #Format axes
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        ax.set_xticklabels([])
        if idx > 0:
            plt.yticks([])
        else:
            plt.ylabel('hPa')

        #Set individual plot title
        plt.title(case_name, fontsize=font_size)

        #Make colorbar on last plot only
        if idx == casenum-1:
            axins = inset_axes(ax,
                        width="5%",
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                       )
            fig.colorbar(cf, cax=axins, orientation="vertical", label=units)

        #Difference with MERRA2 and MERRA2 contours
        ###########################################

        #Set up new set of axes for second row
        ax = fig.add_subplot(nrows, casenum, casenum+idx+1)

        #Plot interpolated contour
        contour = plt.contour(lev_grid, lat_grid, merra_rfield.transpose(transpose_coords=True),
                    colors='black', levels=levs,
                    negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        if idx == 0:
            #Add a legend for the contour lines for first plot only
            legend_elements = [Line2D([0], [0],
                               color=contour.collections[0].get_edgecolor(),
                               label='MERRA2 interp')]

            ax.legend(handles=legend_elements, loc='upper right' )
        #End if

        #Plot difference contour fill
        cf=plt.contourf(lev_grid, lat_grid,
                        (data_array-merra_rfield).transpose(transpose_coords=True),
                        levels=diff_levs, cmap='RdBu_r')
        #Format axes
        plt.yscale("log")
        ax.set_ylim(1000,0.1)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        ax.set_xticklabels([])
        if idx > 0:
            plt.yticks([])
        else:
            plt.ylabel('hPa')

        #Set individual plot title
        local_title = f'{case_name}\n {delta_symbol} from MERRA2'
        plt.title(local_title, fontsize=font_size)

        #Make colorbar on last plot only
        if idx == casenum-1:
            axins = inset_axes(ax,
                        width="5%",
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                       )
            fig.colorbar(cf, cax=axins, orientation="vertical", label=units)

        #Difference with SABER and SABER contours
        #########################################

        #Set up new set of axes for third row
        ax = fig.add_subplot(nrows, casenum, (casenum*2)+idx+1)

        #Plot interpolated contour
        contour = plt.contour(lev_grid, lat_grid, saber_rfield.transpose(transpose_coords=True),
                    colors='black', levels=levs,
                    negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        if idx == 0:
            #Add a legend for the contour lines for first plot only
            legend_elements = [Line2D([0], [0],
                               color=contour.collections[0].get_edgecolor(),
                               label='SABER interp')]

            ax.legend(handles=legend_elements, loc='upper right')
        #End if

        #Plot difference contour fill
        cf=plt.contourf(lev_grid, lat_grid,
                        (data_array-saber_rfield).transpose(transpose_coords=True),
                        levels=diff_levs, cmap='RdBu_r')

        #Format axes
        plt.yscale("log")
        ax.set_ylim(1000,0.1)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if idx > 0:
            plt.yticks([])
        else:
            plt.ylabel('hPa')

        #Set individual plot title
        local_title = f'{case_name}\n {delta_symbol} from SABER'
        plt.title(local_title, fontsize=font_size)

        #Make colorbar on last plot only
        if idx == casenum-1:
            axins = inset_axes(ax,
                        width="5%",
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                       )
            fig.colorbar(cf, cax=axins, orientation="vertical", label=units)

    #Set up main plot title
    if time_avg == "month":
        str_interval = month_dict[interval].lower().capitalize()
    else:
        str_interval = interval
    fig.suptitle(f"Zonal Mean {cam_var} - {str_interval}",fontsize=16,y=0.93)

    #Add plot to website (if enabled):
    #Special ADF variable which contains the output paths for
    #all generated plots and tables:
    #plot_locations = adfobj.get_basic_info('cam_diag_plot_loc', required=True)
    #plot_loc = Path(plot_locations) / case_names[0]
    plot_locations = adfobj.plot_location
    #print("plot_locations",plot_locations)
    plot_loc = Path(plot_locations[0])
    plot_type = "png"
    plot_name = plot_loc / f"{cam_var}_{str_interval}_Zonal_Mean_scycle.{plot_type}"

    #adfobj.add_website_data(plot_name, cam_var, case_name, season=str_interval, plot_type="Zonal", category="SeasonalCycle")
    #Write the figure to provided workspace/file:
    fig.savefig(plot_name, bbox_inches='tight', dpi=300)
    #plt.savefig(plot_name,dpi=300)
    



def polar_car_temp(hemi, case_names, case_runs, cases_monthly, merra2_monthly):
    """
    """

    if hemi == "s":
        slat = -90
        nlat = -60
    if hemi == "n":
        slat = 60
        nlat = 90

    nplots = len(case_names)
    if nplots > 4:
        ncols = 4
    else:
        ncols = nplots
    #End if
    ncols = 4
    nrows = int(np.ceil(nplots/ncols))

    fig = plt.figure(figsize=(2*7,nrows*5))

    #for run in range(len(runs)):
    for idx,case_name in enumerate(case_names):
        ds = case_runs[case_name]
        ds_month = cases_monthly[case_name]

        rfield_seas = np.zeros((12,len(ds['lev']),len(ds['lat'])))
        rfield_seas = xr.DataArray(rfield_seas, dims=['month','lev', 'lat'],
                                            coords={'month': np.arange(1,13,1),
                                                    'lev': ds['lev'],
                                                    'lat': ds['lat']})

        case_seas = np.zeros((12,len(ds['lev']),len(ds['lat'])))
        case_seas = xr.DataArray(case_seas, dims=['month','lev', 'lat'],
                                 coords={'month': np.arange(1,13,1),
                                         'lev': ds['lev'],
                                         'lat': ds['lat']})
        #Make array of monthly temp data
        for m in range(0,12):
            rfield_seas[m] = merra2_monthly['T'][month_dict[m+1]].interp(lat=ds['lat'], lev=ds['lev'],
                                                                method='linear')
            case_seas[m] = ds_month['T'][month_dict[m+1]]

        #Average over set of latitudes
        merra2_pcap = coslat_average(rfield_seas,slat,nlat)
        case_pcap = coslat_average(case_seas,slat,nlat)

        #
        [time_grid, lev_grid] = np.meshgrid(ds['lev'],np.arange(0,12))

        #Set up plot
        ax = fig.add_subplot(nrows, ncols, idx+1)

        cf=plt.contourf(lev_grid, time_grid, (case_pcap-merra2_pcap),#.transpose(transpose_coords=True),
                        levels=temp_diff_levs,cmap='RdBu_r'
                       )

        c=plt.contour(lev_grid, time_grid, merra2_pcap, colors='black',
                           levels=temp_levs,
                           negative_linestyles='dashed',
                           linewidths=.5, alpha=0.5)
        #Format the axes
        plt.yscale("log")
        ax.set_ylim(300,1)
        ax.set_yticks([300,100,30,10])
        ax.set_xticks(np.arange(0,12,2),rotation=40)
        ax.set_xticklabels(('Jan','Mar','May','Jul','Sep','Nov'),rotation=40)
        if idx > 0:
            plt.yticks([])
        else:
            plt.ylabel('hPa')

        #Set title
        local_title=case_names[idx]+"\n minus MERRA2"
        plt.title(local_title, fontsize=font_size)

        # Calculate the required wspace based on the length of titles
        #title_lengths = [len(ax.get_title()) for ax in axs]
        #max_title_length = max(title_lengths)
        #required_wspace = .003 * max_title_length  # Adjust the multiplier as needed
        #required_wspace = 0.
        # Adjust the wspace dynamically
        #plt.subplots_adjust(wspace=required_wspace)

        #Check for start of new row
        if idx % 4 == 0:
            row = idx // 4 + 1

        #Check to see where the colorbar will go
        #The idea is to if the plots fill up each row, put the colorbar on last plot of row
        #If the row isn't filled up, put the color bar on last possible plot of row
        if ((4*(row-1) < idx < 4*(row+1)) and (idx == nplots-1)) or ((idx+1) % 4 == 0):
                #if idx == nplots-1:
                axins = inset_axes(ax,
                                width="5%",
                                height="80%",
                                loc='center right',
                                borderpad=-2.5
                               )
                fig.colorbar(cf, cax=axins, orientation="vertical", label='K')

    plt.savefig(f'temp_{hemi}pcap_merra2.png',dpi=300)

########

# Helper functions
##################
def make_zm_files(adfobj,hist_loc,case_name,calc_var_list,syr,eyr,return_ds=True):
    """
    Make zonal mean files from history monthly files

    args:
    -----
       * hist_loc: Path object
          - place to find history files
       * case_name: str
          - name fo current case
       * calc_var_list: list
          - list of variables to compute and save zonal means
       * syr, eyr
          - start and end desired climo years
       * return_ds: boolean
          - return the dataset to xarray DataSet object   

    output: netcdf file
    ------
       - case specific file name with case data, saved to where???????
    """

    save_path = adfobj.get_basic_info('diag_loc', required=True)
    if not Path(f"{save_path}/waccm_135_{case_name}.nc").exists():
        h0_lists = []

        for yr in np.arange(int(syr),int(eyr)+1):
            h0_lists.append(sorted(glob.glob(f'{hist_loc}*cam.h0.{yr}-*')))

        h0_list = list(chain(*h0_lists))

        waccm_135 = xr.open_mfdataset(h0_list, use_cftime=True, data_vars=calc_var_list)
        waccm_135 = waccm_135[calc_var_list].mean(dim='lon')
        
        waccm_135.to_netcdf(f"{save_path}/waccm_135_{case_name}.nc")
    else:
        waccm_135 = xr.open_mfdataset(f"{save_path}/waccm_135_{case_name}.nc")
    if return_ds:
        return waccm_135
########

def saber_data(filename = "../SABER_monthly_2002-2014.nc"):
    """

    """
    saber = {}
    saber_seasonal = {}
    saber_monthly = {}
    saber_vars = ['u','temp','lat','lev']

    saber_ncfile = xr.open_dataset(filename, decode_times=True, use_cftime=True)
    saber_ncfile = saber_ncfile.rename({"latitude":"lat"})
    saber_ncfile = saber_ncfile.rename({"pressure":"lev"})

    #WARNING: there is no actual time information in the `time` coordinate!
    # - !! The assigned times are strictly from the file name !!
    start_date = datetime(2002, 1, 1)

    #Grab number of months
    num_months = len(saber_ncfile.time)

    # List to store datetime objects
    datetime_list = []

    # Generate datetime objects incrementally by month
    for i in range(num_months):
        new_date = start_date + relativedelta(months=i)

        # Set the day to the first day of the month
        datetime_list.append(new_date.replace(day=1))

    for index, var in enumerate(saber_vars):
        if var not in saber_seasonal:
            saber_seasonal[var] = {}
        if var not in saber_monthly:
            saber_monthly[var] = {}
        saber[var] = saber_ncfile[var]
        if index < len(saber_vars)-2:
            saber_ncfile[var] = saber_ncfile[var].assign_coords({"time": datetime_list})
            saber[var] = saber_ncfile[var]
            for season in seasons:
                saber_seasonal[var][season] = time_mean(saber_ncfile, saber_ncfile[var],
                                                        time_avg="season", interval=season,
                                                        is_climo=None, obs=True)
            for month in np.arange(1,13,1):
                saber_monthly[var][month_dict[month]] = time_mean(saber_ncfile, saber_ncfile[var],
                                                                  time_avg="month", interval=month,
                                                                  is_climo=None, obs=True)
    return saber_monthly, saber_seasonal
########

def merra_data(filename = "/glade/work/richling/ADF/ADF_dev/notebooks/chem-diags/MERRA2_met.nc"):
    """
    """

    merra2 = {}
    merra2_seasonal = {}
    merra2_monthly = {}
    merra2_vars = ['U','T','V','lat','lev']

    merra_ncfile = xr.open_dataset(filename, decode_times=True, use_cftime=True)
    merra_ncfile = merra_ncfile.sel(time=merra_ncfile.time.values[0])
    merra_ncfile = merra_ncfile.rename({"time":"first-time"})
    merra_ncfile = merra_ncfile.rename({"record":"time"})

    for index, var in enumerate(merra2_vars):

        merra2[var] = merra_ncfile[var]
        if index < len(merra2_vars)-2:

            start_date = datetime(1980, 1, 1)

            # Number of months to generate
            num_months = len(merra_ncfile[var].time)

            # List to store datetime objects
            datetime_list = []

            # Generate datetime objects incrementally by month
            for i in range(num_months):
                new_date = start_date + relativedelta(months=i)
                datetime_list.append(new_date.replace(day=1))  # Set the day to the first day of the month

            merra_ncfile[var] = merra_ncfile[var].assign_coords({"time": datetime_list})
            if var not in merra2_seasonal:
                merra2_seasonal[var] = {}
            if var not in merra2_monthly:
                merra2_monthly[var] = {}
            merra2[var] = merra_ncfile[var]

            for season in seasons:
                merra2_seasonal[var][season] = time_mean(merra_ncfile, merra2[var], time_avg="season", interval=season, is_climo=None, obs=True)
            for month in np.arange(1,13,1):
                merra2_monthly[var][month_dict[month]] = time_mean(merra_ncfile, merra2[var], time_avg="month", interval=month, is_climo=None, obs=True)

    return merra2_monthly, merra2_seasonal
########

def coslat_average(darray, lat1, lat2):
    """
    Calculate the weighted average for an [:,lat] array over the region
    lat1 to lat2
    """

    # flip latitudes if they are decreasing
    if (darray.lat[0] > darray.lat[darray.lat.size -1]):
        print("flipping latitudes")
        darray = darray.sortby('lat')

    region = darray.sel(lat=slice(lat1, lat2))
    weights=np.cos(np.deg2rad(region.lat))
    regionw = region.weighted(weights)
    regionm = regionw.mean("lat")

    return regionm
########

def time_mean(ncfile, data, time_avg, interval, is_climo=None, obs=False):
    """Calculates the time-weighted seasonal average (or average over all time).

    Parameters
    ----------
    data : xarray.DataArray or xarray.Dataset
        data to be averaged
    season : str, optional
        the season to extract from `data`
        If season is `ANN` or None, average all available time.
    is_climo : bool, optional
        If True, expects data to have time or month dimenion of size 12.
        If False, then 'time' must be a coordinate,
        and the `time.dt.days_in_month` attribute must be available.

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        the average of `data` in season `season`

    Notes
    -----
    If the data is a climatology, the code will make an attempt to understand the time or month
    dimension, but will assume that it is ordered from January to December.
    If the data is a climatology and is just a numpy array with one dimension that is size 12,
    it will assume that dimension is time running from January to December.
    """
    if time_avg == "season":
        #seasons = ["ANN", "DJF", "JJA", "MAM", "SON"]
        if interval is not None:
            assert interval in seasons, f"Unrecognized season string provided: '{season}'"
        elif interval is None:
            interval = "ANN"

    #data = data.drop("time")
    #data = data.rename({"record":"times"})
    #data.assign_coords(timez=datetime_list)

    # Can't remove the 'time' dimension, so I will build
    # a `times' dimension to replace the `record` dimension and 
    # update those values to date time objects for month/year

    try:
        month_length = data.time.dt.days_in_month
    except (AttributeError, TypeError):
        print("Nah, nope workingn nope")
        # do our best to determine the temporal dimension and assign weights
        if not is_climo:
            raise ValueError("Non-climo file provided, but without a decoded time dimension.")
        else:
            # CLIMO file: try to determine which dimension is month

            has_time = False
            #print("not has_time")z
            if isinstance(data, xr.DataArray):
                has_time = 'time' in data.dims
                if not has_time:
                    if "month" in data.dims:
                        data = data.rename({"month":"time"})
                        has_time = True
            if not has_time:
                print("not has_time")
                # this might happen if a pure numpy array gets passed in
                # --> assumes ordered January to December.
                assert ((12 in data.shape) and (data.shape.count(12) == 1)), f"Sorry, {data.shape.count(12)} dimensions have size 12, making determination of which dimension is month ambiguous. Please provide a `time` or `month` dimension."
                time_dim_num = data.shape.index(12)
                fakedims = [f"dim{n}" for n in range(len(data.shape))]
                fakedims[time_dim_num] = "time"
                data = xr.DataArray(data, dims=fakedims, attrs=data.attrs)

            timefix = pd.date_range(start='1/1/1999', end='12/1/1999', freq='MS') # generic time coordinate from a non-leap-year
            data = data.assign_coords({"time":timefix})
        month_length = data.time.dt.days_in_month
    #End try/except
    
    if not obs:
        syr = ncfile.time.dt.year.values[0]
        eyr = ncfile.time.dt.year.values[-2]
        #if season == "DJF":
        #    print(ncfile.time.values)
        #    print(np.unique(ncfile.time.dt.year.values))
    
        #    print(syr,eyr)
        data.attrs["year_range"] = f"{syr}-{eyr}"
        timefix = pd.date_range(start=f'1/1/{syr}', end=f'12/1/{eyr}', freq='MS')
        data['time']=timefix
    data.attrs[time_avg] = interval
    if time_avg == "season":
        #data.attrs[time_avg] = interval
        data = data.sel(time=data.time.dt.month.isin(seasons[interval])) # directly take the months we want based on season kwarg
    if time_avg == "month":
        #data.attrs[time_avg] = season
        data = data.sel(time=data.time.dt.month.isin(interval)) # directly take the months we want based on season kwarg

        #try:
        #    del data.attrs['season']
        #except KeyError:
        #   pass

    return data.weighted(data.time.dt.daysinmonth).mean(dim='time', keep_attrs=True)
########










