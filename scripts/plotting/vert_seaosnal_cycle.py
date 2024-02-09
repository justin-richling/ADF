
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy
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



calc_var_list = ['U','T','V','UU','VU','VT','OMEGAU','OMEGA','O3','Q','UTEND_VDIFF']
#calc_var_list = ['U','T','V','O3','Q','UTEND_VDIFF']

merra2_vars = ['U','T','V','lat','lev']

#
# --- Main Function Shares Name with Module: regional_map_multicase ---
#
def make_scycle_maps(adfobj, diag, data_dict, case_deets):
    """
    Chemistry Map main function
        * Initially start with Aerosol Zonal maps
            - This probably can be expanded to LatLon if given single pressure levels?

        
    """

    # Notify user that script has started:
    print("\n  Generating zonal vertical seasonal cycle plots plots ...")

    var_list = adfobj.diag_var_list

    case_names = case_deets["case_names"]["cases"] + case_deets["case_names"]["baseline"]




    cases_coords = {}
    cases_seasonal = {}
    cases_monthly = {}


    #Get MERRA2 data and seasonal and monthly averages
    merra2 = {}
    merra2_seasonal = {}
    merra2_monthly = {}
    
    merra_ncfile = xr.open_dataset("/glade/work/richling/ADF/ADF_dev/notebooks/chem-diags/MERRA2_met.nc",
                                decode_times=True, use_cftime=True)
    #Change MERRA2 coords to match CESM
    #Drop the `time` coordinate and rename to first-time
    # - this is a single index coordinate, really want time, lat, lev only
    merra_ncfile = merra_ncfile.sel(time=merra_ncfile.time.values[0])
    merra_ncfile = merra_ncfile.rename({"time":"first-time"})
    #Rename `record` coordinate to `time` to match CESM
    merra_ncfile = merra_ncfile.rename({"record":"time"})

    for index, var in enumerate(merra2_vars):

        merra2[var] = merra_ncfile[var]
        if index < len(merra2_vars)-2:


            start_date = datetime(1980, 1, 1)

            # Number of months to generate
            num_months = len(merra_ncfile[var].time)#record

            # List to store datetime objects
            datetime_list = []

            # Generate datetime objects incrementally by month
            for i in range(num_months):
                new_date = start_date + relativedelta(months=i)
                datetime_list.append(new_date.replace(day=1))  # Set the day to the first day of the month

            #Make datetime objects for the `time` coordinate
            merra_ncfile[var] = merra_ncfile[var].assign_coords({"time": datetime_list})

            #Initialize nested data dictionaries
            if var not in merra2_seasonal:
                merra2_seasonal[var] = {}
            if var not in merra2_monthly:
                merra2_monthly[var] = {}
            merra2[var] = merra_ncfile[var]

            #Grab all seasonal mean data
            for season in seasons:
                merra2_seasonal[var][season] = time_mean(merra_ncfile, merra2[var],
                                                         time_avg="season",
                                                         interval=season,
                                                         is_climo=None, obs=True)
            #Grab all monthly mean data
            for month in np.arange(1,13,1):
                merra2_monthly[var][month_dict[month]] = time_mean(merra_ncfile, merra2[var],
                                                                   time_avg="month",
                                                                   interval=month,
                                                                   is_climo=None,
                                                                   obs=True)

    #Load all data, process TEM vars, seasonal averaging
    for idx,run in enumerate(runs):
        #run = run[0]
        print(run)
        #zmdir = adfobj.basic_info_dict['diag_loc'] + f"zm/{run}/" #+ "/atm/hist/"
        #file = rundir + "waccm_135.nc"

        hist_loc = adfobj.get_cam_info("cam_hist_loc")[0]
        case_name = adfobj.get_cam_info("cam_case_name")[0]
        syr = adfobj.climo_yrs['syears'][0]
        eyr = adfobj.climo_yrs['eyears'][0]

        ncfile = make_zm_files(hist_loc,case_name,syr,eyr,return_ds=True)


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

        cases_coords[run] = case_coords
        cases_monthly[run] = case_monthly
        cases_seasonal[run] = case_seasonal





    zonal_wind(data_dict, case_names, obs="MERRA2")






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

def seasonal_comparison():
    



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
def make_zm_files(hist_loc,case_name,calc_var_list,syr,eyr,return_ds=True):
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
    h0_lists = []

    for yr in np.arange(int(syr),int(eyr)+1):
        h0_lists.append(sorted(glob.glob(f'{hist_loc}*cam.h0.{yr}-*')))

    h0_list = list(chain(*h0_lists))

    waccm_135 = xr.open_mfdataset(h0_list, use_cftime=True, data_vars=calc_var_list)
    waccm_135 = waccm_135[calc_var_list].mean(dim='lon')
    waccm_135.to_netcdf(f"waccm_135_{case_name}.nc")
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
        if interval is not None:
            assert interval in seasons, f"Unrecognized season string provided: '{interval}'"
        elif interval is None:
            interval = "ANN"

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
    
    #Check if it's CESM and (for now) correct the time
    if not obs:
        syr = ncfile.time.dt.year.values[0]
        eyr = ncfile.time.dt.year.values[-2]

        data.attrs["year_range"] = f"{syr}-{eyr}"
        timefix = pd.date_range(start=f'1/1/{syr}', end=f'12/1/{eyr}', freq='MS')
        data['time']=timefix

    return data.weighted(data.time.dt.daysinmonth).mean(dim='time', keep_attrs=True)
########










def zonal_wind(data_dict, case_names, obs="MERRA2"):
    # CAM vars needed
    # - U

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 
    season="MAM"


    #Do 001 - 001 diff, etc, in addition to direct MERRA2 comparison

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 
    #for run in range(len(runs)):
    for case_name in case_names:
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],
                                        data_runs[runs[run]]['lat'])
        
        contour_levels = np.arange(-120, 121, 10)
        #ax = fig.add_subplot(2, 4, run+1)
        ax = fig.add_subplot(2, len(runs), run+1)
        
        cf=plt.contourf(lev_grid, lat_grid,
                        data_runs_seasonal[runs[run]]['U'][season].transpose(transpose_coords=True),
                        levels=contour_levels)
        
        #cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_seasonal[runs[run]]['U'][season,:,:]), levels=contour_levels)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        plt.title(run_names[run])
        ax.set_xticklabels([])
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='m/s')
    run=0
    for run in range(len(runs)):
        #reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],
        #                                                       merra2_seasonal['U'],
        #                                                       #merra2_seasonal['U'][season],
        #                                                       data_runs[runs[run]]['lat'],
        #                                                       data_runs[runs[run]]['lev'])

        reference_field = merra2_seasonal['U'].interp(lat=data_runs[runs[run]]['lat'],
                                                    lev=data_runs[runs[run]]['lev'], method='linear')

        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(-30, 31, 3)
        ax = fig.add_subplot(2, len(runs), len(runs)+run+1)
        #ax = fig.add_subplot(2, 4, 4+run+1)
        
        cf=plt.contourf(lev_grid, lat_grid,
                        (data_runs_seasonal[runs[run]]['U'][season]-reference_field).transpose(transpose_coords=True),
                        levels=contour_levels)
        
        plt.set_cmap('RdBu_r')
        orig_contour_levels = np.arange(-120, 121, 10)
        plt.contour(lev_grid, lat_grid, reference_field.transpose(transpose_coords=True),
                    colors='black', levels=orig_contour_levels,
                    negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        
        plt.yscale("log")
        ax.set_ylim(1000,0.1)
        plt.yticks([1000,100,10,1,0.1])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        local_title = run_names[run] + " - MERRA2"
        plt.title(local_title) 
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", ticks=np.arange(-30,31,10), label='m/s')
    fig.suptitle(f"Xarray Created Zonal Wind\nXarray package read in file\nADF Seasonal Means - {season}")
    plt.savefig('output/zonal_wind_mam_ADF_seasonal_xarray_zm_xarray_read.png',dpi=300)



"""
def zonal_temp():
    # CAM vars needed
    # - T

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 
    season=2;

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 
    for run in range(len(runs)):
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(160, 301, 10)
        ax = fig.add_subplot(2, 4, run+1)
        cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_seasonal[runs[run]]['T'][season,:,:]), levels=contour_levels)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        plt.title(run_names[run])
        ax.set_xticklabels([])
        if run == 3:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')
    run=0
    #reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_seasonal['T'][season,:,:],data_runs[runs[0]]['lat'],data_runs[runs[0]]['lev'])
    for run in range(len(runs)):
        #reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_seasonal['U'][season,:,:],data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
        reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_seasonal['T'][season,:,:],data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
        
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(-20, 21, 2)
        ax = fig.add_subplot(2, 4, 4+run+1)
        cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_seasonal[runs[run]]['T'][season,:,:]-reference_field), levels=contour_levels)
        plt.set_cmap('RdBu_r')
        orig_contour_levels = np.arange(150, 300, 10)
        plt.contour(lev_grid, lat_grid, np.transpose(reference_field), colors='black', levels=orig_contour_levels, negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        local_title = run_names[run] + " - MERRA2"
        plt.title(local_title) 
        if run == 3:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", ticks=np.arange(-20,21,4), label='K')
    plt.savefig('output/temp_jja_mine.png',dpi=300)



def temp_vs_saber():
    # CAM vars needed
    # - T

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 
    month=11

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 
    for run in range(len(runs)):
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(140, 301, 10)
        ax = fig.add_subplot(2, len(runs), run+1)
        cf=plt.contourf(lev_grid, lat_grid, data_runs_monthly[runs[run]]['T'].transpose(transpose_coords=True), levels=contour_levels)
        #cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_monthly[runs[run]]['T'][month,:,:]), levels=contour_levels)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        plt.title(run_names[run])
        ax.set_xticklabels([])
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')
    run=0

    for run in range(len(runs)):
        #reference_field=interpolating_functions.interpolate_2d(saber['latitude'],saber['pressure'],saber_monthly['temp'][month,:,:],data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
        reference_field = saber_monthly['temp'].interp(latitude=data_runs[runs[run]]['lat'], pressure=data_runs[runs[run]]['lev'], method='linear')
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(-40, 41, 4)
        ax = fig.add_subplot(2, len(runs), len(runs)+run+1)
        cf=plt.contourf(lev_grid, lat_grid, (data_runs_monthly[runs[run]]['T']-reference_field).transpose(transpose_coords=True), levels=contour_levels)
        #cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_monthly[runs[run]]['T'][month,:,:]-reference_field), levels=contour_levels)
        plt.set_cmap('RdBu_r')
        orig_contour_levels = np.arange(140, 300, 10)
        ctours=plt.contour(lev_grid, lat_grid, reference_field.transpose(transpose_coords=True), colors='black', levels=orig_contour_levels, negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        local_title = run_names[run] + " - SABER"
        plt.title(local_title) 
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')
    plt.savefig('output/temp_dec_saber_ADF_monthly_xarray_zm_xarray_read.png',dpi=300)

def qbo():
    # CAM vars needed
    # - U


    plt.rc('font', size=8) 
    run=0
    #Plot zonal wind data
    #fig, axes = plt.subplot_mosaic('AAAABC;DDDDEF;GGGGHI;JJJJKL;MMMMNO')
    #main_key=['A','D','G','J','M']
    #side1_key=['B','E','H','K','N']
    #side2_key=['C','F','I','L','O']

    #fig, axes = plt.subplot_mosaic('AAAABC;DDDDEF;GGGGHI;JJJJKL',figsize=(15,10))
    #main_key=['A','D','G','J']
    #side1_key=['B','E','H','K']
    #side2_key=['C','F','I','L']


    fig, axes = plt.subplot_mosaic('AAAABC;DDDDEF;GGGGHI',figsize=(15,10))
    main_key=['A','D','G']
    side1_key=['B','E','H']
    side2_key=['C','F','I']



    runs_local={}
    #runs_local[0]=fv
    runs_local[0]=data_runs[runs[0]]
    runs_local[1]=data_runs[runs[1]]
    #runs_local[2]=data_runs[runs[2]]

    y = 1.00

    contour_levels = np.arange(-35, 35, 2.5)

    # MERRA2 plot
    plot_num = len(runs)
    merra_plot = plot_num
    nt=len(data_runs[runs[run]]['date'])
    if nt > 120:
        nt = 120
    #else:
    #    nt = nt
    print("MERRA2",nt)

    #plotdata=averaging_functions.coslat_average(merra2['U'],merra2['lat'],-10,10)
    plotdata = cosweightlat(merra2['U'],-10,10)

    plotdata_clip = np.clip(np.abs(plotdata), None, 35)
    plotdata=np.sign(plotdata)*plotdata_clip
    [time_grid, lev_grid] = np.meshgrid(merra2['lev'],np.arange(1,nt+1,1))
    start_ind=252
    end_ind=start_ind+nt
    axes[main_key[merra_plot]].contourf(lev_grid, time_grid, plotdata[start_ind:end_ind,0,:], levels=contour_levels, cmap='RdBu_r')
    axes[main_key[merra_plot]].set_ylim(100,3)
    axes[main_key[merra_plot]].set_yscale("log")
    axes[main_key[merra_plot]].set_ylabel('hPa')
    axes[main_key[merra_plot]].set_title("MERRA2",y=y)
    #axes[main_key[plot_num]].set_xticks(np.arange(1,nt+1,12),rotation=40)
    #axes[main_key[plot_num]].set_xticks(np.arange(1,121,12),rotation=40)
    #axes[main_key[plot_num]].set_xticks(np.arange(int(merra_start_yr),int(merra_end_yr)+1,5),rotation=40)
    #axes[main_key[plot_num]].set_xticklabels(['2000','2001','2002','2003','2004','2005','2006','2007','2008','2009'])
    #axes[main_key[merra_plot]].set_xticklabels(np.arange(int(merra_start_yr),int(merra_end_yr)+1,5))

    axes[main_key[plot_num]].set_xticks(np.arange(1,nt+1,12),rotation=40)
    #axes[main_key[run]].set_xticklabels(np.arange(int(yrs[0]),int(yrs[0])+10,1))
    #axes[main_key[plot_num]].set_xticklabels(['2000','2001','2002','2003','2004','2005','2006','2007','2008','2009',"2010",
    #                                          '2011','2012','2013'])

    start_year = int(str(plotdata[252].record.values)[0:4])
    axes[main_key[plot_num]].set_xticklabels(np.arange(start_year,start_year+(nt/12),1).astype(int))

    amp_m = qbo_functions.qbo_amplitude(plotdata[:,0,:])
    axes[side1_key[merra_plot]].plot(amp_m,merra2['lev'],color='k',linewidth=1)
    axes[side1_key[merra_plot]].set_ylim(100,3)
    axes[side1_key[merra_plot]].set_yscale("log")
    axes[side1_key[merra_plot]].set_xlim(0,20)
    axes[side1_key[merra_plot]].set_xticks(np.arange(0,21,5))
    axes[side1_key[merra_plot]].set_xlabel('m/s')
    axes[side1_key[merra_plot]].set_yticks([])

    period_m = qbo_functions.qbo_frequency(plotdata[:,0,:])
    axes[side2_key[merra_plot]].plot(period_m,merra2['lev'],color='k',linewidth=1)
    axes[side2_key[merra_plot]].set_ylim(100,3)
    axes[side2_key[merra_plot]].set_yscale("log")
    axes[side2_key[merra_plot]].set_xlim(0,40)
    axes[side2_key[merra_plot]].set_xticks(np.arange(0,41,10))
    axes[side2_key[merra_plot]].set_yticks([])
    axes[side2_key[merra_plot]].set_xlabel('months')

    for run in range(len(runs)):
        yrs = wowsa[runs[run]]
        nt=len(runs_local[run]['date'])
        if nt > 120:
            nt_sub = 120
        else:
            nt_sub = nt
        print(run_names[run],nt,"\n")
        
        if run ==0:
            [time_grid, lev_grid] = np.meshgrid(runs_local[run]['lev'],np.arange(1,nt+1+12,1))
        else:
            [time_grid, lev_grid] = np.meshgrid(runs_local[run]['lev'],np.arange(1,nt+1,1))
        contour_levels = np.arange(-35, 35, 2.5)
        #ax = fig.add_subplot(4, 4, run,colspan=3)
        #plotdata=averaging_functions.coslat_average(runs_local[run]['U'],runs_local[run]['lat'],-10,10)
        plotdata = cosweightlat(runs_local[run]['U'],-10,10)
        plotdata_clip = np.clip(np.abs(plotdata), None, 35)
        plotdata=np.sign(plotdata)*plotdata_clip
        if run == 0:
            start_idx = 11
            axes[main_key[run]].contourf(lev_grid[start_idx:nt_sub+12,:], time_grid[start_idx:nt_sub+12,:], plotdata[start_idx:nt_sub+12,:], levels=contour_levels, cmap='RdBu_r')
        else:
            start_idx = 0
            axes[main_key[run]].contourf(lev_grid[start_idx:nt_sub,:], time_grid[start_idx:nt_sub,:], plotdata[start_idx:nt_sub,:], levels=contour_levels, cmap='RdBu_r')
        axes[main_key[run]].set_ylim(100,3)
        axes[main_key[run]].set_yscale("log")
        axes[main_key[run]].set_ylabel('hPa')
        axes[main_key[run]].set_title(run_names[run],y=y)
        if run == 0:
            axes[main_key[run]].set_xticks(np.arange(12,nt_sub+12,12),rotation=40)
        else:
            axes[main_key[run]].set_xticks(np.arange(1,nt_sub,12),rotation=40)
        #axes[main_key[run]].set_xticks([])
        #if nt > 120:
        if run == 0:
            axes[main_key[run]].set_xticklabels(np.arange(int(yrs[0]+1),int(yrs[0])+10+1,1))
        else:
            axes[main_key[run]].set_xticklabels(np.arange(int(yrs[0]),int(yrs[0])+10,1))
        
        amp = qbo_functions.qbo_amplitude(plotdata)
        axes[side1_key[run]].plot(amp,runs_local[run]['lev'],linewidth=1)
        axes[side1_key[run]].plot(amp_m,merra2['lev'],color='k',linewidth=1)
        axes[side1_key[run]].set_ylim(100,3)
        axes[side1_key[run]].set_yscale("log")
        axes[side1_key[run]].set_xlim(0,20)
        axes[side1_key[run]].set_xticks(np.arange(0,20,5))
        #axes[side1_key[run]].set_xticks([])
        axes[side1_key[run]].set_yticks([])
        if run==0:
            axes[side1_key[run]].set_title('Amplitude',y=y)
            
        
        period = qbo_functions.qbo_frequency(plotdata)
        axes[side2_key[run]].plot(period,runs_local[run]['lev'],linewidth=1)
        axes[side2_key[run]].plot(period_m,merra2['lev'],color='k',linewidth=1)
        axes[side2_key[run]].set_ylim(100,3)
        axes[side2_key[run]].set_yscale("log")
        axes[side2_key[run]].set_xlim(0,40)
        axes[side2_key[run]].set_xticks(np.arange(0,40,10))
        #axes[side2_key[run]].set_xticks([])
        axes[side2_key[run]].set_yticks([])
        if run==0:
            axes[side2_key[run]].set_title('Period',y=y)

    plt.savefig('output/qbo_xarray_zm_xarray_read.png',dpi=300)


def polar_cap_temp():
    # CAM vars needed
    # - T

    month=11

    #Pcap seasonal cycle
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 

    #reference_field_seas=np.zeros((12,len(data_runs[runs[0]]['lev']),len(data_runs[runs[0]]['lat'])))
    #for month in range(0,12):
    #    reference_field_seas[month,:,::]=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_monthly['T'][month,:,:],data_runs[runs[0]]['lat'],data_runs[runs[0]]['lev'])
    #merra2_pcap=averaging_functions.coslat_average(reference_field_seas,data_runs[runs[0]]['lat'],-90,-60)
        
    for run in range(len(runs)):
        reference_field_seas=np.zeros((12,len(data_runs[runs[run]]['lev']),len(data_runs[runs[run]]['lat'])))
        for month in range(0,12):
            reference_field_seas[month,:,::]=interpolating_functions.interpolate_2d(merra2['lat'],
                                                                                    merra2['lev'],
                                                                                    merra2_monthly['T'][month,:,:],
                                                                                    data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
            
        merra2_pcap=averaging_functions.coslat_average(reference_field_seas,data_runs[runs[run]]['lat'],-90,-60)

        [time_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],np.arange(0,12))
        contour_levels = np.arange(-10,11,1)
        #contour_levels = np.arange(-35,16,1)
        
        ax = fig.add_subplot(2, 4, run+1)
        
        run_pcap=averaging_functions.coslat_average(data_runs_monthly[runs[run]]['T'],data_runs[runs[run]]['lat'],-90,-60)

        #, levels=contour_levels
        cf=plt.contourf(lev_grid, time_grid, run_pcap-merra2_pcap,
                        levels=contour_levels
                    )
        plt.yscale("log")
        ax.set_ylim(300,1)
        ax.set_yticks([300,100,30,10])
        #ax.set_xticks(np.arange(0,12,2),rotation=40)
        ax.set_xticklabels(('Jan','Mar','May','Jul','Sep','Nov'),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        local_title=run_names[run]+" - MERRA2"
        plt.title(local_title)
        ctours=plt.contour(lev_grid, time_grid, merra2_pcap, colors='black', levels=orig_contour_levels, negative_linestyles='dashed', linewidths=.5, alpha=0.5)

        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-2.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')

    plt.savefig('output/temp_spcap_mine.png',dpi=300)
"""