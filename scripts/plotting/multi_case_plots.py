"""
Module: regional_map_multicase

Provides a plot with regional maps of specified variables for all cases (up to 10) in a row.

Since this is a specialty plot, it looks for several custom options to be provided in the YAML configuration file. For example, one might include this block in the YAML:

region_multicase:
    region_spec: [slat, nlat, wlon, elon]
    region_time_option: <calendar | zeroanchor>  If calendar, will look for specified years. If zeroanchor will use a nyears starting from year_offset from the beginning of timeseries
    region_start_year:
    region_end_year:
    region_nyear:
    region_year_offset:
    region_month: <NULL means look for season>
    region_season: <NULL means use annual mean>
    region_variables: <list of variables to try to use; allows for a subset of the total diag variables>

"""
#
# --- imports and configuration ---
#
from pathlib import Path
import warnings  # use to warn user about missing files.

import numpy as np
import xarray as xr
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from plotting_functions import pres_from_hybrid, prep_contour_plot


def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + "\n"


warnings.formatwarning = my_formatwarning
#
# --- Main Function Shares Name with Module: regional_map_multicase ---
#
def regional_map_multicase(adfobj):
    """
    regional_map_multicase
    input -> ADF object

    Sketch of workflow:
    - check for regional options (see module docstring),
    - get case names/locations (time series are used),
    - determine plot location and type
    - detect per-variable plot options
    - Loop through variables: make one multi-panel plot per variable
    """

    # Notify user that script has started:
    print("\n  Generating regional contour plots ...")

    # We need to know:
    # - Variable to plot
    # - Region to plot
    # - Years to include in average
    # - Months/Season to include in average
    # - Cases ... reference on left, test(s) to the right --> IMPOSE UPPER LIMIT OF 10 PANELS

    #
    # Check if regional options were specified... can not proceed without them
    #
    multi_case_plots = adfobj.read_config_var("multi_case_plots")


    # case information
    case_names = adfobj.get_cam_info("cam_case_name", required=True)
    if len(case_names) > 10:
        print("ERROR: multi_case_plots is limited to <= 10 cases.")
        return

    #Grab test case nickname(s)
    test_nicknames = adfobj.get_cam_info('case_nickname')
    if test_nicknames == None:
        test_nicknames = case_names

    # CAUTION:
    # "data" here refers to either obs or a baseline simulation,
    # Until those are both treated the same (via intake-esm or similar)
    # we will do a simple check and switch options as needed:
    if adfobj.get_basic_info("compare_obs"):

        #Extract variable-obs dictionary:
        var_obs_dict = adfobj.var_obs_dict
        base_nickname = "Obs"

        #If dictionary is empty, then  there are no observations to regrid to,
        #so quit here:
        if not var_obs_dict:
            print("No observations found to plot against, so no lat/lon maps will be generated.")
            return

    else:
        data_name = adfobj.get_baseline_info("cam_case_name", required=True) # does not get used, is just here as a placemarker
        data_list = [data_name] # gets used as just the name to search for climo files HAS TO BE LIST
        #data_loc  = model_rgrid_loc #Just use the re-gridded model data path

        syear_baseline = adfobj.climo_yrs["syear_baseline"]
        eyear_baseline = adfobj.climo_yrs["eyear_baseline"]

        #Grab baseline case nickname
        base_nickname = adfobj.get_baseline_info('case_nickname')
        if base_nickname == None:
            base_nickname = data_name

    nicknames = test_nicknames.append(base_nickname)
    case_loc = adfobj.get_cam_info("cam_ts_loc", required=True)


    #
    # Determine input "reference case"
    #
    if adfobj.get_basic_info("compare_obs"):
        print("NotImplementedError: the multi case plots do not have an observation option yet.")
        return
        # TODO: add observation (need time series), use `var_obs_dict`.
    else:
        #Extract required baseline run info:
        data_name = adfobj.get_baseline_info("cam_case_name", required=True)
        data_loc = adfobj.get_baseline_info("cam_ts_loc", required=True)
    #
    # Set plot options
    #
    res = adfobj.variable_defaults  # dict of variable-specific plot preferences
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get("plot_type", "png")
    redo_plot = adfobj.get_basic_info("redo_plot")
    plot_loc = _get_plot_location(adfobj)

    ptype_order_dict = {'global_latlon_map': ["LatLon"],
                            'zonal_mean': ["Zonal"],
                            'global_latlon_vect_map': ["LatLon_Vector"],
                            'polar_map': ["NHPolar","SHPolar"],
                            'cam_taylor_diagram': ["TaylorDiag"],
                            'time_series':['TimeSeries'],
                            'top_10':['Top10']}

    seasons = ["ANN","DJF","MAM","JJA","SON"]








#############

def _get_plot_location(adfobj):
    """
    Determine plot location based on ADF object.
    Create location if it is not a directory.
    Return a Path object to the directory.
    """
    plot_location = adfobj.plot_location  # ADF output path for plots and tables
    if not plot_location:
        plot_location = adfobj.get_basic_info("cam_diag_plot_loc")
    if isinstance(plot_location, list):
        for pl in plot_location:
            plpth = Path(pl)
            # Check if plot output directory exists, and if not, then create it:
            if not plpth.is_dir():
                print(f"\t    {pl} not found, making new directory")
                plpth.mkdir(parents=True)
        if len(plot_location) == 1:
            plot_loc = Path(plot_location[0])
        else:
            print(
                f"Ambiguous plotting location since all cases go on same plot. Will put them in first location: {plot_location[0]}"
            )
            plot_loc = Path(plot_location[0])
    else:
        plot_loc = Path(plot_location)
    return plot_loc

####

def make_plots(nicknames):


    #hspace values for subplots based off number of cases (plots) with figsize=(15,15)
    #These will need to be changed if the figsize changes!
    hspace_dict = {
    **dict.fromkeys([2, 3], 0), 
    **dict.fromkeys([4,5,6], -0.72),
    **dict.fromkeys([7,8,9], -0.6), 
    **dict.fromkeys([10,11,12], -0.35),
    **dict.fromkeys([13,14,15], 0.4),
    }

    titles = []
    #ncols = int(np.sqrt(nplots)) + 1
    ncols = 3
    nplots = len(nicknames)
    nrows = int(np.ceil(nplots/ncols))
    if nrows < 2:
        nrows = 2

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,figsize=(15,15), facecolor='w', edgecolor='k',
                            sharex=True,
                            sharey=True,
                            subplot_kw={"projection": ccrs.PlateCarree()})

    count = 0
    for l in range(0,nrows):
        for c in range(0,ncols):
                
                if count < nplots:            
                    axs[l,c].coastlines()
                    gl = axs[l,c].gridlines(draw_labels=True, linewidth=1)
                    gl.top_labels = False
                    gl.right_labels = False
                    if (axs[l,c].get_subplotspec().is_first_col()) and (axs[l,c].get_subplotspec().is_first_row()):
                        titles.append(axs[l,c].set_title(nicknames[-1],loc='left'))
                    else:
                        titles.append(axs[l,c].set_title(nicknames[count],loc='left'))
                        
                    if axs[l,c].get_subplotspec().is_first_col():
                        gl.left_labels = True
                    else:
                        gl.left_labels = False
                else:
                    axs[l,c].set_visible(False)
                count = count + 1
                
    plt.subplots_adjust(wspace=0.1, hspace=hspace_dict[nplots])
    plt.savefig(f"multi_case_plots_{nplots}_cases.png",bbox_inches="tight")


