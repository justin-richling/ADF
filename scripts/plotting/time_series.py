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
    
    """
    # Case names:
    # NOTE: "baseline" == "reference" == "observations" will be called `base`
    #       test case(s) == case(s) to be diagnosed  will be called `case` (assumes a list)
    case_names = adfobj.get_cam_info('cam_case_name', required=True)  # Loop over these
    case_climo_loc = adfobj.get_cam_info('cam_climo_loc', required=True)

    # ADF variable which contains the output path for plots and tables:
    plot_location = adfobj.plot_location
    if not plot_location:
        plot_location = adfobj.get_basic_info("cam_diag_plot_loc")
    if isinstance(plot_location, list):
        for pl in plot_location:
            plpth = Path(pl)
            #Check if plot output directory exists, and if not, then create it:
            if not plpth.is_dir():
                print(f"\t    {pl} not found, making new directory")
                plpth.mkdir(parents=True)
        if len(plot_location) == 1:
            plot_loc = Path(plot_location[0])
        else:
            print(f"Ambiguous plotting location since all cases go on same plot. Will put them in first location: {plot_location[0]}")
            plot_loc = Path(plot_location[0])
    else:
        plot_loc = Path(plot_location)

    # CAUTION:
    # "data" here refers to either obs or a baseline simulation,
    # Until those are both treated the same (via intake-esm or similar)
    # we will do a simple check and switch options as needed:
    if adfobj.get_basic_info("compare_obs"):
        data_name = "obs"  # does not get used, is just here as a placemarker
        data_list = adfobj.read_config_var("obs_type_list")  # Double caution!
        data_loc = adfobj.get_basic_info("obs_climo_loc", required=True)
    else:
        data_name = adfobj.get_baseline_info('cam_case_name', required=True)
        data_list = data_name # should not be needed (?)
        data_loc = adfobj.get_baseline_info("cam_climo_loc", required=True)
    """
    
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

    #Check if the variables needed for the Taylor diags are present,
    #If not then skip this script:
    taylor_var_set = {'TS', 'SST', 'FSNT', 'FLNT'}
    """if not taylor_var_set.issubset(adfobj.diag_var_list) or \
       (not ('PRECT' in adfobj.diag_var_list) and (not ('PRECL' in adfobj.diag_var_list) or not ('PRECC' in adfobj.diag_var_list))):
        print("\tThe Taylor Diagrams require the variables: ")
        print("\tU, PSL, SWCF, LWCF, PRECT (or PRECL and PRECC), LANDFRAC, TREFHT, TAUX, RELHUM,T")
        print("\tSome variables are missing so Taylor diagrams will be skipped.")
        return
    #End if"""

    #Set seasonal ranges:
    seasons = {"ANN": np.arange(1,13,1),}
    '''           "DJF": [12, 1, 2],
               "JJA": [6, 7, 8],
               "MAM": [3, 4, 5],
               "SON": [9, 10, 11]}'''

    # TAYLOR PLOT VARIABLES:
    var_list = ['TS', 'SST', 'FSNT', 'FLNT','RESTOM']

    #Loop over variables:
    for var in var_list:

        #Loop over model cases:
        for case_idx, case_name in enumerate(case_names):

            # LOOP OVER SEASON
            #
            for s in seasons:
                #Set output plot location:
                plot_loc = Path(plot_locations[case_idx])
                plot_name = plot_loc / f"cam_{s}_TimeSeries_Mean.{plot_type}"
                print(f"\t - Plotting Time Series, {s}")

                # Check redo_plot. If set to True: remove old plot, if it already exists:
                if (not redo_plot) and plot_name.is_file():
                    continue
                elif (redo_plot) and plot_name.is_file():
                    plot_name.unlink()

                """ # hold the data in a DataFrame for each case
                # variable | correlation | stddev ratio | bias
                df_template = pd.DataFrame(index=var_list, columns=['corr', 'ratio', 'bias'])
                result_by_case = {cname: df_template.copy() for cname in case_names}
                #"""
                """# LOOP OVER VARIABLES
                #
                for v in var_list:
                    base_x = _retrieve(adfobj, v, data_name, data_loc) # get the baseline field
                    for casenumber, case in enumerate(case_names):     # LOOP THROUGH CASES
                        case_x = _retrieve(adfobj, v, case, case_climo_loc[casenumber])
                        # ASSUMING `time` is 1-12, get the current season:
                        case_x = case_x.sel(time=seasons[s]).mean(dim='time')
                        result_by_case[case].loc[v] = taylor_stats_single(case_x, base_x)"""
                    
                
                ts_aave_case,ts_aave_base,yrs = get_data(case_ts_loc,data_ts_loc,var)

                vals_case = [ts_aave_case.sel(time=i).mean() for i in yrs]
                vals_base = [ts_aave_base.sel(time=i).mean() for i in yrs]
                    
                #
                # -- PLOTTING (one per season) --
                #
                fig, ax = plt.subplots(figsize=(15,6))
                plt.subplots_adjust(wspace=0.52,hspace=0.2)
                        
                ax = ts_plot(ax,case_names,data_name,vals_case,vals_base)    

                # add text with variable names:
                #txtstrs = [f"{i+1} - {v}" for i, v in enumerate(var_list)]
                fig.legend(loc='upper left', bbox_to_anchor=(0.33, 0.02))
                #fig.text(0.9, 0.9, "\n".join(txtstrs), va='top')
                fig.savefig(plot_name, bbox_inches='tight')
                print(f"\t Time Series: completed {s}. \n\t File: {plot_name}")

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

def get_data(case_ts_loc,data_ts_loc,var):#var='TS'
    fils_case = sorted(list(Path(case_ts_loc[0]).glob(f"*v{var}*.nc")))
    ts_case = xr.open_mfdataset(fils_case)[f'{var}'].compute()
    w_case = np.cos(np.radians(ts_case.lat))  # area weighting
    ts_aave_case = ts_case.weighted(w_case).mean(dim=("lat","lon"))

    fils_base = sorted(list(Path(data_ts_loc).glob(f"*{var}*.nc")))
    ts_base = xr.open_mfdataset(fils_base)[var].compute()
    w_base = np.cos(np.radians(ts_base.lat))  # area weighting
    ts_aave_base = ts_base.weighted(w_base).mean(dim=("lat","lon"))

    woo = [val.item().strftime()[0:4] for _,val in enumerate(ts_case["time"])]
    yrs = np.unique(woo)
    return ts_aave_case,ts_aave_base,yrs

def ts_plot(ax, case_name, data_name, vals_case, vals_base, yrs):
    for i, c in enumerate(case_name):
        labels = [data_name, c]
    
    ax.plot(yrs, vals_case, c="b", label=case_name)
    ax.scatter(yrs, vals_base, marker="*",c="orange", s=100,label=data_name)

    ax.set_xticks(yrs)

    ax.set_xlabel("Time")
    return ax