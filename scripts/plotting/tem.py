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

    syear_cases = adf.climo_yrs["syears"]
    eyear_cases = adf.climo_yrs["eyears"]

    #Grab all case nickname(s)
    test_nicknames = adf.case_nicknames["test_nicknames"]
    base_nickname = adf.case_nicknames["base_nickname"]

    # CAUTION:
    # "data" here refers to either obs or a baseline simulation,
    # Until those are both treated the same (via intake-esm or similar)
    # we will do a simple check and switch options as needed:
    if adf.get_basic_info("compare_obs"):

        #Extract variable-obs dictionary:
        var_obs_dict = adf.var_obs_dict
        base_nickname = "Obs"

        #If dictionary is empty, then  there are no observations to regrid to,
        #so quit here:
        if not var_obs_dict:
            print("No observations found to plot against, so no lat/lon maps will be generated.")
            return

    else:
        data_name = adf.get_baseline_info("cam_case_name", required=True) # does not get used, is just here as a placemarker

    #End if

    #case_names = case_names + [data_name]

    #Extract baseline years (which may be empty strings if using Obs):
    syear_baseline = adf.climo_yrs["syear_baseline"]
    eyear_baseline = adf.climo_yrs["eyear_baseline"]

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
    
    #Loop over model cases:
    for idx,case_name in enumerate(case_names):
        #Ope the TEM file
        output_loc_idx = Path(output_loc) / case_name
        tem = output_loc_idx / f'{case_name}.TEMdiag.nc'

        ds = xr.open_dataset(tem)

        #Location to save plots
        plot_name = str(plot_location)+"/"+case_name+"_tem.png"

        #Plot TEM
        nrows = 5
        ncols = 2
        fig_width = 15
        fig_height = 15+(3*nrows) #try and dynamically create size of fig based off number of cases (therefore rows)
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,fig_height),
                                facecolor='w', edgecolor='k', sharex=True)
        
        # Row 1
        axs[0,0].set_title('First Plot')
        ds.uzm.isel(time=-1).plot(ax=axs[0,0], y='lev', yscale='log',ylim=[1e3,1])
        axs[0,0].set_title('uzm')
        ds.vzm.isel(time=-1).plot(ax=axs[0,1], y='lev', yscale='log',ylim=[1e3,1])
        axs[0,1].set_title('vzm')

        # Row 2
        ds.epfy.isel(time=-1).plot(ax=axs[1,0], y='lev', yscale='log',vmax=1e6,ylim=[1e2,1])
        axs[1,0].set_title('epfy')
        ds.epfz.isel(time=-1).plot(ax=axs[1,1], y='lev', yscale='log',vmax=1e5,ylim=[1e2,1])
        axs[1,1].set_title('epfz')

        # Row 3
        ds.vtem.isel(time=-1).plot.contourf(ax=axs[2,0], levels = 21, y='lev', yscale='log',
                                            vmax=3,vmin=-3,ylim=[1e2,1], cmap='RdBu_r')
        ds.vtem.isel(time=-1).plot.contour(ax=axs[2,0], levels = 11, y='lev', yscale='log',
                                            vmax=3,vmin=-3,ylim=[1e2,1],
                                            colors='black', linestyles=None)
        axs[2,0].set_title('vtem')

        ds.wtem.isel(time=-1).plot.contourf(ax=axs[2,1], levels = 21, y='lev', yscale='log',
                                            vmax=0.005,vmin=-0.005,ylim=[1e2,1], cmap='RdBu_r')
        ds.wtem.isel(time=-1).plot.contour(ax=axs[2,1], levels = 7, y='lev', yscale='log',
                                           vmax=0.03,vmin=-0.03,ylim=[1e2,1], 
                                           colors='black', linestyles=None)
        axs[2,1].set_title('wtem')

        # Row 4
        ds.psitem.isel(time=-1).plot.contourf(ax=axs[3,0], levels = 21, y='lev', yscale='log',
                                              vmax=5e9,ylim=[1e2,2])
        axs[3,0].set_title('psitem')

        ds.utendepfd.isel(time=-1).plot(ax=axs[3,1], y='lev', yscale='log',
                                        vmax=0.0001,vmin=-0.0001,ylim=[1e2,2])
        axs[3,1].set_title('utendepfd')

        # Row 5
        ds.utendvtem.isel(time=-1).plot(ax=axs[4,0], y='lev', yscale='log',vmax=0.001, ylim=[1e3,1])
        axs[4,0].set_title('utendvtem')
        ds.utendwtem.isel(time=-1).plot(ax=axs[4,1], y='lev', yscale='log',vmax=0.0001, ylim=[1e3,1])
        axs[4,1].set_title('utendwtem')
        
        #Adjust subplots
        hspace = 0.3
        plt.subplots_adjust(wspace=0.3, hspace=hspace)

        #Set figure title
        plt.suptitle(f'TEM Diagnostics: {test_nicknames[idx]} - ANN\nyrs: {syear_cases[idx]} - {eyear_cases[idx]}', fontsize=16, y=.91)

        #Write the figure to provided workspace/file:
        fig.savefig(plot_name, bbox_inches='tight', dpi=300)

        #Add plot to website (if enabled):
        adf.add_website_data(plot_name, "TEM", case_name, season="ANN")


##############
#END OF SCRIPT