from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
import plotting_functions as pf
import matplotlib.pyplot as plt

import warnings  # use to warn user about missing files.

def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'

warnings.formatwarning = my_formatwarning

def time_series(adfobj):

    """
    This script plots zonal averages.
    Compare CAM climatologies against
    other climatological data (observations or baseline runs).
    Description of needed inputs from ADF:
    case_name        -> Name of CAM case provided by "cam_case_name".
    model_rgrid_loc  -> Location of re-gridded CAM climo files provided by "cam_regrid_loc".
    data_name        -> Name of data set CAM case is being compared against,
                        which is always either "obs" or the baseline CAM case name,
                        depending on whether "compare_obs" is true or false.
    data_loc         -> Location of comparison data, which is either "obs_climo_loc"
                        or "cam_baseline_climo_loc", depending on whether
                        "compare_obs" is true or false.
    var_list         -> List of CAM output variables provided by "diag_var_list"
    data_list        -> List of data sets CAM will be compared against, which
                        is simply the baseline case name in situations when
                        "compare_obs" is false.
    plot_location    -> Location where plot files will be written to, which is
                        specified by "cam_diag_plot_loc".
    Notes:
        The script produces plots of 2-D and 3-D variables,
        but needs to determine which type along the way.
        For 3-D variables, the default behavior is to interpolate
        climo files to pressure levels, which requires the hybrid-sigma
        coefficients and surface pressure. That ASSUMES that the climo
        files are using native hybrid-sigma levels rather than being
        transformed to pressure levels.
    """

    #Notify user that script has started:
    print("\n  Generating time series plots...")

    # Extract needed quantities from ADF object:
    # -----------------------------------------
    
    case_names = adfobj.get_cam_info('cam_case_name', required=True)
    data_name = adfobj.get_baseline_info('cam_case_name', required=True)

    case_ts_loc = adfobj.get_cam_info("cam_ts_loc", required=True)
    data_ts_loc = adfobj.get_baseline_info("cam_ts_loc", required=True)

    plot_locations = adfobj.plot_location
    
    res = adfobj.variable_defaults # dict of variable-specific plot preferences
    # or an empty dictionary if use_defaults was not specified in YAML.

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')
    print(f"\t NOTE: Plot type is set to {plot_type}")

    #Check if existing plots need to be redone
    redo_plot = adfobj.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")

    #Set seasonal ranges:
    seasons = {"ANN": np.arange(1,13,1),}
    '''           "DJF": [12, 1, 2],
               "JJA": [6, 7, 8],
               "MAM": [3, 4, 5],
               "SON": [9, 10, 11]}'''

    #Time Series variables:
    #Might be best as a list of variables in the config file?
    #var_list = ['TS', 'SST', 'FSNT', 'FLNT', "RESTOM"]
    var_list = ['TS', 'FSNT', 'FLNT',"RESTOM"]

    #Loop over model cases:
    for case_idx, case_name in enumerate(case_names):
        #Loop over seasons
        for s in seasons:
            #Loop over variables:
            for var in var_list:
                #Set output plot location:
                plot_loc = Path(plot_locations[case_idx])
                plot_name = plot_loc / f"{var}_{s}_TimeSeries_Mean.{plot_type}"

                #Check if plot output directory exists, and if not, then create it:
                if not plot_loc.is_dir():
                    print(f"    {plot_loc} not found, making new directory")
                    plot_loc.mkdir(parents=True)
                print(f"\t - Plotting Time Series, {var} {s}")

                # Check redo_plot. If set to True: remove old plot, if it already exists:
                if (not redo_plot) and plot_name.is_file():
                    continue
                elif (redo_plot) and plot_name.is_file():
                    plot_name.unlink()
                    
                if var == "RESTOM":
                    ave_case,ave_base,unit,yrs_case,yrs_base = get_restom_data(case_ts_loc,data_ts_loc)
                else:
                    ave_case,ave_base,unit,yrs_case,yrs_base = get_data(case_ts_loc,data_ts_loc,var)

                vals_case = [ave_case.sel(time=i).mean() for i in yrs_case]
                vals_base = [ave_base.sel(time=i).mean() for i in yrs_base]
                    
                fig, ax = plt.subplots(figsize=(15,6))
                plt.subplots_adjust(wspace=0.52,hspace=0.2)
                        
                ax = ts_plot(ax,case_names,data_name,vals_case,vals_base,unit,yrs_case,yrs_base)    

                #Set legend for case plots
                fig.legend(bbox_to_anchor=(-0.1, .88, 1., .102), loc="right",
                    ncol=1, borderaxespad=0.0)

                #Set Main title for subplots:
                fig.suptitle(f"Time Series: {var} - {s}", fontsize=15,horizontalalignment = 'left', x=fig.subplotpars.left,y=0.92)
                fig.savefig(plot_name, bbox_inches='tight', dpi=300)
                #print(f"\t Time Series: completed {s}. \n\t File: {plot_name}")
            
            #End for (variables loop)
        #End for (seasons loop)
    #End for (case names loop)

    #Notify user that script has ended:
    print("  ...Time Series have been generated successfully.")

#
# Helpers
#
def _load_dataset(fils):
    if len(fils) == 0:
        warnings.warn(f"Input file list is empty.")
        return None
    elif len(fils) > 1:
        return xr.open_mfdataset(fils, combine='by_coords')
    else:
        sfil = str(fils[0])
        return xr.open_dataset(sfil)

def get_data(case_ts_loc,data_ts_loc,var):
    ave_case,unit,yrs_case = _data_calcs(case_ts_loc[0],var)
    ave_base,unit,yrs_base = _data_calcs(data_ts_loc,var)

    return ave_case,ave_base,unit,yrs_case,yrs_base

def get_restom_data(case_ts_loc,data_ts_loc):
    ave_case_FSNT,unit,yrs_case = _data_calcs(case_ts_loc[0],'FSNT')
    ave_base_FSNT,unit,yrs_base = _data_calcs(data_ts_loc,'FSNT')

    ave_case_FLNT,_,_ = _data_calcs(case_ts_loc[0],"FLNT")
    ave_base_FLNT,_,_ = _data_calcs(data_ts_loc,"FLNT")

    restom_case = ave_case_FSNT - ave_case_FLNT
    restom_base = ave_base_FSNT - ave_base_FLNT

    return restom_case,restom_base,unit,yrs_case,yrs_base

def ts_plot(ax, case_names, data_name, vals_case, vals_base, unit, yrs_case, yrs_base):
    yrs_lists = []
    for _,val in enumerate(case_names):
        ax.plot(yrs_case, vals_case, c="b", label=val)
        yrs_lists.append(yrs_case)
    ax.plot(yrs_base, vals_base, "--", c="orange",label=data_name)
    ax.set_ylabel(unit,fontsize=15,labelpad=20)
    ax.set_xlabel("Years",fontsize=15,labelpad=20)

    yrs_lists.append(yrs_base)
    #list1 = [1,2,4,5,7,6,5]
    #ist2 = [1,4,5,6]
    #list3 = [2,5,4,6,7,8,5,3]
    #list4 = [5,7,3,7,5,2,3,8]

    #lists = [list1, list2, list3, list4]

    list_max = filter(lambda i: len(i) == max([len(l) for l in yrs_lists]), yrs_lists)
    #print(list(list_max))
    ax.set_xticks(list_max, list_max, rotation=45, ha='right', rotation_mode='anchor')
    
    
    
    ax.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
    
    return ax

def _data_calcs(ts_loc,var):
    fils = sorted(list(Path(ts_loc).glob(f"*{var}*.nc")))
    ts = xr.open_mfdataset(fils)[var].compute()
    w = np.cos(np.radians(ts.lat))  # area weighting
    ave  = ts.weighted(w).mean(dim=("lat","lon")) # global averaging
    unit = ts.units
    
    dtimes = [val.item().strftime()[0:4] for _,val in enumerate(ts["time"])]
    yrs = np.unique(dtimes)

    return ave,unit,yrs