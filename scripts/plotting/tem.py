#Import standard modules:
from pathlib import Path
import numpy as np
import xarray as xr
import warnings  # use to warn user about missing files.
from datetime import date
import matplotlib.pyplot as plt

from scipy import integrate
from numpy import ma

#Format warning messages:
def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = my_formatwarning

def tem(adf):
    """
    This script/function is designed to generate global
    2-D lat/lon maps of model fields with continental overlays.
    Description of needed inputs:
    case_name         -> Name of CAM case provided by "cam_case_name".
    model_rgrid_loc   -> Location of re-gridded CAM climo files provided by "cam_regrid_loc".
    data_name         -> Name of data set CAM case is being compared against,
                         which is always either "obs" or the baseline CAM case name,
                         depending on whether "compare_obs" is true or false.
    data_loc          -> Location of comparison data, which is either the location listed
                         in each variable's ""obs_file", or the same as "model_rgrid_loc",
                         depending on whether "compare_obs" is true or false.
    var_list          -> List of CAM output variables provided by "diag_var_list"
    data_list         -> List of data sets CAM will be compared against, which
                         is simply the baseline case name in situations when
                         "compare_obs" is false.
    plot_location     -> Location where plot files will be written to, which is
                         specified by "cam_diag_plot_loc".
    climo_yrs         -> Dictionary containing the start and end years of the test
                        and baseline model data (if applicable).
    variable_defaults -> optional,
                         Dict that has keys that are variable names and values that are plotting preferences/defaults.
    """

    #Import necessary modules:
    #------------------------
    import pandas as pd

    #CAM diagnostic plotting functions:
    #import plotting_functions as pf
    #-------------------------

    # Steps:
    # - load regridded climo files for model and obs
    # - calculate all-time and seasonal fields (from individual months)
    # - Take difference, calculate statistics
    # - make plot

    #
    # Use ADF api to get all necessary information
    #

    #Special ADF variable which contains the output paths for
    #all generated plots and tables for each case:
    plot_location = Path(adf.plot_location[0])
    #Check if plot output directory exists, and if not, then create it:
    if not plot_location.is_dir():
        print("    {} not found, making new directory".format(plot_location))
        plot_location.mkdir(parents=True)

    #CAM simulation variables (this is always assumed to be a list):
    case_names = adf.get_cam_info("cam_case_name", required=True)

    # CAUTION:
    # "data" here refers to either obs or a baseline simulation,
    # Until those are both treated the same (via intake-esm or similar)
    # we will do a simple check and switch options as needed:
    if adf.get_basic_info("compare_obs"):

        #Extract variable-obs dictionary:
        var_obs_dict = adf.var_obs_dict

        #If dictionary is empty, then  there are no observations to regrid to,
        #so quit here:
        if not var_obs_dict:
            print("No observations found to plot against, so TEM maps won't be generated.")
            return
    else:
        base_name = adf.get_baseline_info("cam_case_name", required=True) # does not get used, is just here as a placemarker
    #End if

    #Extract test case years
    syear_cases = adf.climo_yrs["syears"]
    eyear_cases = adf.climo_yrs["eyears"]

    #Extract baseline years (which may be empty strings if using Obs):
    syear_baseline = adf.climo_yrs["syear_baseline"]
    eyear_baseline = adf.climo_yrs["eyear_baseline"]

    #Grab all case nickname(s)
    test_nicknames = adf.case_nicknames["test_nicknames"]
    base_nickname = adf.case_nicknames["base_nickname"]
    case_nicknames = test_nicknames.append(base_nickname)

    #case_names = case_names + [base_name]
    if base_name:
        case_names.append(base_name)
        syear_cases.append(syear_baseline)
        eyear_cases.append(eyear_baseline)

    '''if syear_baseline:
        syear_cases.append(syear_baseline)
    if eyear_baseline:
        eyear_cases.append(eyear_baseline)'''

    res = adf.variable_defaults # will be dict of variable-specific plot preferences
    # or an empty dictionary if use_defaults was not specified in YAML.

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adf.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')
    print(f"\t NOTE: Plot type is set to {plot_type}")

    # check if existing plots need to be redone
    redo_plot = adf.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")
    #-----------------------------------------


    #Location to saved TEM netCDF files
    output_loc = adf.get_basic_info("tem_loc")
    #If path not specified, skip TEM calculation?
    if output_loc is None:
        return
    else:
        #Notify user that script has started:
        print("\n  Generating TEM plots...")
    
    #Set seasonal ranges:
    seasons = {"ANN": np.arange(1,13,1),
               "DJF": [12, 1, 2],
               "JJA": [6, 7, 8],
               "MAM": [3, 4, 5],
               "SON": [9, 10, 11]
               }

    var_list_og = ['uzm','vzm','epfy','epfz','vtem','wtem',
     'psitem','utendepfd','utendvtem','utendwtem']

    if "qbo" in adf._AdfDiag__plotting_scripts:
        var_list = ['uzm','epfy','epfz','vtem','wtem',
                    'psitem','utendepfd','utendvtem','utendwtem']
    else:
        var_list = ['uzm','epfy','epfz','vtem','wtem','psitem','utendepfd']

    print(var_list)
    
    #Setup TEM plots
    nrows = len(var_list)
    ncols = 2
    fig_width = 15
    fig_height = 15+(3*nrows) #try and dynamically create size of fig based off number of cases (therefore rows)

    #Loop over season dictionary:
    for s in seasons:
        
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,fig_height),
                                facecolor='w', edgecolor='k', sharex=True)

        #Loop over model cases:
        for idx,case_name in enumerate(case_names):

        #Loop over season dictionary:
        #for s in seasons:

            #Extract start and end year values:
            start_year = syear_cases[idx]
            end_year   = eyear_cases[idx]

            #Ope the TEM file
            output_loc_idx = Path(output_loc) / case_name
            tem = output_loc_idx / f'{case_name}.TEMdiag_{start_year}-{end_year}.nc'

            ds = xr.open_dataset(tem)

            #Location to save plots
            plot_name = plot_location / f"{case_name}_{s}_TEM_Mean.png"

            """#Plot TEM
            nrows = 5
            ncols = 2
            fig_width = 15
            fig_height = 15+(3*nrows) #try and dynamically create size of fig based off number of cases (therefore rows)
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,fig_height),
                                    facecolor='w', edgecolor='k', sharex=True)"""

            #Setup and plot the sub-plots
            if len(case_names) > 1:
                tem_plot(ds, idx, axs, s, var_list, res)

            if len(case_names) == 1:
                tem_plot_single(ds, axs, s, var_list, res)


            #Set figure title
            yrs = f"{syear_cases[idx]} - {eyear_cases[idx]}"
            plt.suptitle(f'TEM Diagnostics: {case_nicknames[idx]} - {s}\nyrs: {yrs}', 
                            fontsize=16, y=.91)

            #Write the figure to provided workspace/file:
            fig.savefig(plot_name, bbox_inches='tight', dpi=300)

            #Add plot to website (if enabled):
            adf.add_website_data(plot_name, "TEM", case_name, season=s)

def tem_plot_single(ds, axs, s, var_list, res):

    for var in var_list:
        mdata = ds[var].squeeze()
        vres = res[var]

        #Create array to avoid weighting missing values:
        md_ones = xr.where(mdata.isnull(), 0.0, 1.0)

        #Calculate monthly weights based on number of days:
        month_length = mdata.time.dt.days_in_month
        weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())

        #Create array to avoid weighting missing values:
        md_ones = xr.where(mdata.isnull(), 0.0, 1.0)

        #Calculate monthly weights based on number of days:
        month_length = mdata.time.dt.days_in_month
        weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())


        #Calculate monthly-weighted seasonal averages:
        mseasons = {}
        if s == 'ANN':

            #Calculate annual weights (i.e. don't group by season):
            weights_ann = month_length / month_length.sum()

            mseasons[s] = (mdata * weights_ann).sum(dim='time')
            mseasons[s] = mseasons[s] / (md_ones*weights_ann).sum(dim='time')


        else:
            #this is inefficient because we do same calc over and over
            mseasons[s] = (mdata * weights).groupby("time.season").sum(dim="time").sel(season=s)
            wgt_denom = (md_ones*weights).groupby("time.season").sum(dim="time").sel(season=s)
            mseasons[s] = mseasons[s] / wgt_denom

        """if vres["vmin"]:
            vmin = vres["vmin"]
        else:
            vmin = None
        if vres["vmax"]:
            vmax = vres["vmax"]
        else:
            vmax = None
        if vres["cmap"]:
            cmap = vres["cmap"]
        else:
            cmap = None

        mseasons[s].plot(ax=axs[vres["axs"][0],vres["axs"][1]], y='lev', yscale='log',
                         vmin=vmin, vmax=vmax, ylim=vres["ylim"],
                         cmap=cmap, cbar_kwargs={'label': vres["units"]})
        
        axs[vres["axs"][0],vres["axs"][1]].set_title(vres["long_name"])"""

        # Row 1
        #------------------------------------------------------------------------------------------
        if var == "uzm":
            mseasons[s].plot(ax=axs[0,0], y='lev', yscale='log',ylim=[1e3,1],
                                    cbar_kwargs={'label': ds[var].units})
            axs[0,0].set_title(ds[var].long_name)
        
        if var == "vzm":
            mseasons[s].plot(ax=axs[0,1], y='lev', yscale='log',ylim=[1e3,1],
                                    cbar_kwargs={'label': ds[var].units})
            axs[0,1].set_title(ds[var].long_name)

        # Row 2
        #------------------------------------------------------------------------------------------
        if var == "epfy":
            mseasons[s].plot(ax=axs[1,0], y='lev', yscale='log',vmax=1e6,ylim=[1e2,1],
                                    cbar_kwargs={'label': ds[var].units})
            axs[1,0].set_title(ds[var].long_name)
        
        if var == "epfz":
            mseasons[s].plot(ax=axs[1,1], y='lev', yscale='log',vmax=1e5,ylim=[1e2,1],
                                    cbar_kwargs={'label': ds[var].units})
            axs[1,1].set_title(ds[var].long_name)

        # Row 3
        #------------------------------------------------------------------------------------------
        if var == "vtem":
            mseasons[s].plot.contourf(ax=axs[2,0], levels = 21, y='lev', yscale='log',
                                                vmax=3,vmin=-3,ylim=[1e2,1], cmap='RdBu_r',
                                                cbar_kwargs={'label': ds[var].units})
            mseasons[s].plot.contour(ax=axs[2,0], levels = 11, y='lev', yscale='log',
                                                vmax=3,vmin=-3,ylim=[1e2,1],
                                                colors='black', linestyles=None)
            axs[2,0].set_title(ds[var].long_name)

        if var == "wtem":
            mseasons[s].plot.contourf(ax=axs[2,1], levels = 21, y='lev', yscale='log',
                                                vmax=0.005, vmin=-0.005, ylim=[1e2,1], cmap='RdBu_r',
                                                cbar_kwargs={'label': ds[var].units})
            mseasons[s].plot.contour(ax=axs[2,1], levels = 7, y='lev', yscale='log',
                                            vmax=0.03, vmin=-0.03, ylim=[1e2,1],
                                            colors='black', linestyles=None)
            axs[2,1].set_title(ds[var].long_name)

        # Row 4
        #------------------------------------------------------------------------------------------
        if var == "psitem":
            mseasons[s].plot.contourf(ax=axs[3,0], levels = 21, y='lev', yscale='log',
                                                vmax=5e9, ylim=[1e2,2],
                                                cbar_kwargs={'label': ds[var].units})
            axs[3,0].set_title(ds[var].long_name)

        if var == "utendepfd":
            mseasons[s].plot(ax=axs[3,1], y='lev', yscale='log',
                                            vmax=0.0001, vmin=-0.0001, ylim=[1e2,2],
                                            cbar_kwargs={'label': ds[var].units})
            axs[3,1].set_title(ds[var].long_name)

        # Row 5
        #------------------------------------------------------------------------------------------
        if var == "utendvtem":
            mseasons[s].plot(ax=axs[4,0], y='lev', yscale='log',vmax=0.001, ylim=[1e3,1],
                                            cbar_kwargs={'label': ds[var].units})
            axs[4,0].set_title(ds[var].long_name)

        if var == "utendwtem":
            mseasons[s].plot(ax=axs[4,1], y='lev', yscale='log',vmax=0.0001, ylim=[1e3,1],
                                            cbar_kwargs={'label': ds[var].units})
            axs[4,1].set_title(ds[var].long_name)
            

    """# Row 1
    #------------------------------------------------------------------------------------------
    ds.uzm.isel(time=-1).plot(ax=axs[0,0], y='lev', yscale='log',ylim=[1e3,1],
                            cbar_kwargs={'label': ds.uzm.units})
    axs[0,0].set_title(ds.uzm.long_name)
    ds.vzm.isel(time=-1).plot(ax=axs[0,1], y='lev', yscale='log',ylim=[1e3,1],
                            cbar_kwargs={'label': ds.vzm.units})
    axs[0,1].set_title(ds.vzm.long_name)

    # Row 2
    #------------------------------------------------------------------------------------------
    ds.epfy.isel(time=-1).plot(ax=axs[1,0], y='lev', yscale='log',vmax=1e6,ylim=[1e2,1],
                            cbar_kwargs={'label': ds.epfy.units})
    axs[1,0].set_title(ds.epfy.long_name)
    ds.epfz.isel(time=-1).plot(ax=axs[1,1], y='lev', yscale='log',vmax=1e5,ylim=[1e2,1],
                            cbar_kwargs={'label': ds.epfz.units})
    axs[1,1].set_title(ds.epfz.long_name)

    # Row 3
    #------------------------------------------------------------------------------------------
    ds.vtem.isel(time=-1).plot.contourf(ax=axs[2,0], levels = 21, y='lev', yscale='log',
                                        vmax=3,vmin=-3,ylim=[1e2,1], cmap='RdBu_r',
                                        cbar_kwargs={'label': ds.vtem.units})
    ds.vtem.isel(time=-1).plot.contour(ax=axs[2,0], levels = 11, y='lev', yscale='log',
                                        vmax=3,vmin=-3,ylim=[1e2,1],
                                        colors='black', linestyles=None)
    axs[2,0].set_title(ds.vtem.long_name)

    ds.wtem.isel(time=-1).plot.contourf(ax=axs[2,1], levels = 21, y='lev', yscale='log',
                                        vmax=0.005, vmin=-0.005, ylim=[1e2,1], cmap='RdBu_r',
                                        cbar_kwargs={'label': ds.wtem.units})
    ds.wtem.isel(time=-1).plot.contour(ax=axs[2,1], levels = 7, y='lev', yscale='log',
                                    vmax=0.03, vmin=-0.03, ylim=[1e2,1],
                                    colors='black', linestyles=None)
    axs[2,1].set_title(ds.wtem.long_name)

    # Row 4
    #------------------------------------------------------------------------------------------
    ds.psitem.isel(time=-1).plot.contourf(ax=axs[3,0], levels = 21, y='lev', yscale='log',
                                        vmax=5e9, ylim=[1e2,2],
                                        cbar_kwargs={'label': ds.psitem.units})
    axs[3,0].set_title(ds.psitem.long_name)

    ds.utendepfd.isel(time=-1).plot(ax=axs[3,1], y='lev', yscale='log',
                                    vmax=0.0001, vmin=-0.0001, ylim=[1e2,2],
                                    cbar_kwargs={'label': ds.utendepfd.units})
    axs[3,1].set_title(ds.utendepfd.long_name)

    # Row 5
    #------------------------------------------------------------------------------------------
    ds.utendvtem.isel(time=-1).plot(ax=axs[4,0], y='lev', yscale='log',vmax=0.001, ylim=[1e3,1],
                                    cbar_kwargs={'label': ds.utendvtem.units})
    axs[4,0].set_title(ds.utendvtem.long_name)
    ds.utendwtem.isel(time=-1).plot(ax=axs[4,1], y='lev', yscale='log',vmax=0.0001, ylim=[1e3,1],
                                    cbar_kwargs={'label': ds.utendwtem.units})
    axs[4,1].set_title(ds.utendwtem.long_name)"""

    #Adjust subplots
    hspace = 0.3
    plt.subplots_adjust(wspace=0.3, hspace=hspace)

    return axs


def tem_plot(ds, idx, axs, s, var_list, res):

    for var in var_list:
        mdata = ds[var].squeeze()
        vres = res[var]

        #Create array to avoid weighting missing values:
        md_ones = xr.where(mdata.isnull(), 0.0, 1.0)

        #Calculate monthly weights based on number of days:
        month_length = mdata.time.dt.days_in_month
        weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())

        #Create array to avoid weighting missing values:
        md_ones = xr.where(mdata.isnull(), 0.0, 1.0)

        #Calculate monthly weights based on number of days:
        month_length = mdata.time.dt.days_in_month
        weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())

        #Calculate monthly-weighted seasonal averages:
        mseasons = {}
        if s == 'ANN':

            #Calculate annual weights (i.e. don't group by season):
            weights_ann = month_length / month_length.sum()

            mseasons[s] = (mdata * weights_ann).sum(dim='time')
            mseasons[s] = mseasons[s] / (md_ones*weights_ann).sum(dim='time')


        else:
            #this is inefficient because we do same calc over and over
            mseasons[s] = (mdata * weights).groupby("time.season").sum(dim="time").sel(season=s)
            wgt_denom = (md_ones*weights).groupby("time.season").sum(dim="time").sel(season=s)
            mseasons[s] = mseasons[s] / wgt_denom

        #Run through vars and plot each against the baseline on same row
        #Each column will be a case, ie (test, base), or (test, test, base) , ...
        #                         column: 0  ,   1         0,     1,    2     ...

        # Var 1
        #------------------------------------------------------------------------------------------
        if var == "uzm":
            mseasons[s].plot(ax=axs[0,idx], y='lev', yscale='log',ylim=[1e3,1],
                                    cbar_kwargs={'label': ds[var].units})
            axs[0,idx].set_title(ds[var].long_name)

        # Var 2
        #------------------------------------------------------------------------------------------
        if var == "epfy":
            mseasons[s].plot(ax=axs[1,idx], y='lev', yscale='log',vmax=1e6,ylim=[1e2,1],
                                    cbar_kwargs={'label': ds[var].units})
            axs[1,idx].set_title(ds[var].long_name)
        
        # Var 3
        #------------------------------------------------------------------------------------------
        if var == "epfz":
            mseasons[s].plot(ax=axs[2,idx], y='lev', yscale='log',vmax=1e5,ylim=[1e2,1],
                                    cbar_kwargs={'label': ds[var].units})
            axs[1,idx].set_title(ds[var].long_name)

        # Var 4
        #------------------------------------------------------------------------------------------
        if var == "vtem":
            mseasons[s].plot.contourf(ax=axs[3,idx], levels = 21, y='lev', yscale='log',
                                                vmax=3,vmin=-3,ylim=[1e2,1], cmap='RdBu_r',
                                                cbar_kwargs={'label': ds[var].units})
            mseasons[s].plot.contour(ax=axs[3,idx], levels = 11, y='lev', yscale='log',
                                                vmax=3,vmin=-3,ylim=[1e2,1],
                                                colors='black', linestyles=None)
            axs[3,idx].set_title(ds[var].long_name)

        # Var 5
        #------------------------------------------------------------------------------------------
        if var == "wtem":
            mseasons[s].plot.contourf(ax=axs[4,idx], levels = 21, y='lev', yscale='log',
                                                vmax=0.005, vmin=-0.005, ylim=[1e2,1], cmap='RdBu_r',
                                                cbar_kwargs={'label': ds[var].units})
            mseasons[s].plot.contour(ax=axs[4,idx], levels = 7, y='lev', yscale='log',
                                            vmax=0.03, vmin=-0.03, ylim=[1e2,1],
                                            colors='black', linestyles=None)
            axs[4,idx].set_title(ds[var].long_name)

        # Var 6
        #------------------------------------------------------------------------------------------
        if var == "psitem":
            mseasons[s].plot.contourf(ax=axs[5,idx], levels = 21, y='lev', yscale='log',
                                                vmax=5e9, ylim=[1e2,2],
                                                cbar_kwargs={'label': ds[var].units})
            axs[5,idx].set_title(ds[var].long_name)

        # Var 7
        #------------------------------------------------------------------------------------------
        if var == "utendepfd":
            mseasons[s].plot(ax=axs[6,idx], y='lev', yscale='log',
                                            vmax=0.0001, vmin=-0.0001, ylim=[1e2,2],
                                            cbar_kwargs={'label': ds[var].units})
            axs[6,idx].set_title(ds[var].long_name)

        # Var 8
        #------------------------------------------------------------------------------------------
        if var == "utendvtem":
            mseasons[s].plot(ax=axs[7,idx], y='lev', yscale='log',vmax=0.001, ylim=[1e3,1],
                                            cbar_kwargs={'label': ds[var].units})
            axs[7,idx].set_title(ds[var].long_name)

        # Var 9
        #------------------------------------------------------------------------------------------

        if var == "utendwtem":
            mseasons[s].plot(ax=axs[8,idx], y='lev', yscale='log',vmax=0.0001, ylim=[1e3,1],
                                            cbar_kwargs={'label': ds[var].units})
            axs[8,idx].set_title(ds[var].long_name)
##############
#END OF SCRIPT