#Import standard modules:
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import timedelta

import warnings  # use to warn user about missing files.

#Format warning messages:
def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = my_formatwarning

def time_series(adfobj):
    """
    This script plots time series.
    Compare CAM climatologies against other
    climatological data (observations or baseline runs).

    Description of needed inputs from ADF:
    case_name        -> Name of CAM case provided by "cam_case_name".

    ts_loc           -> Location of CAM time series files provided by "cam_ts_loc".

    data_name        -> Name of data set CAM case is being compared against,
                        which is always either "obs" or the baseline CAM case name,
                        depending on whether "compare_obs" is true or false.

    ts_var_list      -> List of CAM output variables provided by "timeseries_var_list".

    data_list        -> List of data sets CAM will be compared against, which
                        is simply the baseline case name in situations when
                        "compare_obs" is false.

    plot_location    -> Location where plot files will be written to, which is
                        specified by "cam_diag_plot_loc".
    Notes:
        * This script runs annual/seasonal and global weighting.
        * It will be pretty flexible for the variables plotted and layout of figure.
        * This currently only works for single case comparison
            - multi-case comparison is in the works. 02/2023 - JR
    """

    #CAM diagnostic plotting functions:
    import plotting_functions as pf
    #-------------------------

    #Notify user that script has started:
    print("\n  Generating time series plots...")

    #Extract needed quantities from ADF object:
    #-----------------------------------------
    
    #DUMMY wont need eventually:
    main_site_assets_path = "./"

    #List of desired (if available) CAM variables from config file
    var_list = adfobj.diag_var_list
    
    #Check if ocean or land fraction exist
    #in the variable list:
    for var in ["OCNFRAC", "LANDFRAC"]:
        if var in var_list:
            #If so, then move them to the front of variable list so
            #that they can be used to mask or vertically interpolate
            #other model variables if need be:
            var_idx = var_list.index(var)
            var_list.pop(var_idx)
            var_list.insert(0,var)
        else:
            #Since the masking is important, add these just in case the user 
            #forgot to specifiy in the config file
            #NOTE: this may still break if these aren't in the CAM output history files!!
            var_list = [var] + var_list
        #End if
    #End if

    #pressure levels:
    pres_levs = adfobj.get_basic_info("plot_press_levels")

    case_names = adfobj.get_cam_info('cam_case_name', required=True)
    var_defaults = adfobj.variable_defaults

    #Check for multi-case diagnostics
    if len(case_names) > 1:
        #Grab all multi-case diagnostic directories
        #main_site_path = adfobj.main_site_paths["main_site_path"]
        main_site_assets_path = adfobj.main_site_paths["main_site_assets_path"]
        #main_site_img_path = adfobj.main_site_paths["main_site_img_path"]
        case = None
        multi_case = True
    else:
        case = case_names[0]
        multi_case = False

    case_ts_loc = adfobj.get_cam_info("cam_ts_loc", required=True)

    #Grab all case nickname(s)
    test_nicknames = adfobj.case_nicknames["test_nicknames"]
    if not test_nicknames:
        test_nicknames = case_names
    base_nickname = adfobj.case_nicknames["base_nickname"]

    #CAUTION:
    #"data" here refers to either obs or a baseline simulation,
    #Until those are both treated the same (via intake-esm or similar)
    #we will do a simple check and switch options as needed:
    if adfobj.compare_obs:
        print("NOTE: the ADF currently can't plot observational time series, so only test case will plot!")

        base_nickname = "Obs"
        data_name = "Obs"
        all_case_names = case_names
        all_nicknames = test_nicknames
        case_ts_locs = case_ts_loc
        
        #Currently we don't have time series files for obs, so just skip
        # the plotting of obs data.
        #NOTE: the ts plots will still be created, just the test case will plot
        
        #Code commented out below: for obs if we get time series files
        '''
        #Extract variable-obs dictionary:
        var_obs_dict = adfobj.var_obs_dict
        
        #If dictionary is empty, then there are no observations to compare against,
        #so quit here:
        if not var_obs_dict:
            print("No observations found to plot against. So just the test case will be plotted. I pitty the fool who doesn't have obs values! RRRRRah!!\n")
            
            #Just keep the test case info
            all_case_names = case_names
            all_nicknames = test_nicknames
            case_ts_locs = case_ts_loc
            all_nicknames = test_nicknames
        else:
            #Bundle all case names
            all_case_names = case_names + [base_nickname]
            #Gather all nicknames
            all_nicknames = test_nicknames + [base_nickname]
            case_ts_locs = case_ts_loc + [data_ts_loc]
        #End if var_obs_dict
        
        '''
    else:
        data_name = adfobj.get_baseline_info("cam_case_name")
        data_ts_loc = adfobj.get_baseline_info("cam_ts_loc")
        #Bundle all case names
        all_case_names = case_names + [data_name]
        case_ts_locs = case_ts_loc + [data_ts_loc]
        all_nicknames = test_nicknames + [base_nickname]
    #End if

    #Get number of cases
    case_num = len(all_case_names)

    #Special ADF variable which contains the output paths for plots:
    plot_location = adfobj.plot_location
    plot_loc = Path(plot_location[0])   

    res = adfobj.variable_defaults #dict of variable-specific plot preferences
    #or an empty dictionary if use_defaults was not specified in config (YAML) file.

    #Set plot file type:
    #-- this should be set in basic_info_dict, but is not required
    #-- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info", required=True)
    plot_type = basic_info_dict.get('plot_type', 'png')
    print(f"\t NOTE: Plot type is set to {plot_type}")

    #Check if existing plots need to be redone
    redo_plot = adfobj.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to '{redo_plot}'")

    #Set seasons:
    seasons = ["ANN","DJF","MAM","JJA","SON"]
    
    """
    syear_cases = adfobj.climo_yrs["syears"]
    eyear_cases = adfobj.climo_yrs["eyears"]
    #Extract baseline years (which may be empty strings if using Obs):
    syear_baseline = adfobj.climo_yrs["syear_baseline"]
    eyear_baseline = adfobj.climo_yrs["eyear_baseline"]
    
    syears = syear_cases + [syear_baseline]
    eyears = eyear_cases + [eyear_baseline]
    """

    #Set up the plots
    #################

    #Add more colors as needed for number of test cases
    #** Baseline is already added as green dashed line in plotting function **
    #matplotlib colors here: https://matplotlib.org/stable/gallery/color/named_colors.html
    colors = ["k", "aqua", "r", "b", "magenta",
              "orange", "slategrey", "rosybrown"]
    
    #Make a list for vars to skip plotting if desired
    #Check the variable defaults yaml file to add
    # example
    #VAR:
    #  timeseries: 
    #    skip_plot: True
    skip_list = []

    #Create/reset new variable that potentially stores the re-gridded
    #ocean fraction xarray data-array:
    ocn_frc_da = None
    ocn_frc_da = {case_names[0]:None,data_name:None}
    
    #Dictionary for vars with vertical levels specified by user
    var_lev_dict = {}
    
    #Loop over variables:
    #--------------------
    for var in var_list:
        
        #Initialize nested dictionary for each variable
        var_lev_dict[var] = {}
        
        #TODO: Add regional subset so this can change when implemented - JR
        title_var = "Global"
        
        #Extract defaults for variable:
        var_default_dict = var_defaults.get(var, {})
        
        #Check res for any variable specific options that need to be used BEFORE going to the plot:
        if var in res:
            vres = res[var]
            #If found then notify user, assuming debug log is enabled:
            adfobj.debug_log(f"time_series: Found variable defaults for {var}")

            #Extract category (if available):
            web_category = vres.get("category", None)
        else:
            vres = {}
            web_category = None
        #End if

        #Add variables that user doesn't want to plot
        if vres.get('timeseries'):
            if vres["timeseries"].get('skip_plot'):
                skip_list.append(var)
        #End if

        if var not in skip_list:
            print(f"\t - time series plots for {var}")

        #Set plotting parameters based off whether the user wants rolling average
        #Currently only RESTOM is defaulted to 5-yr rolling avg
        #Check the variable defaults yaml file
        #Example to add
        #VAR:
        #  timeseries: 
        #    rolling:
        #      years: 5
        rolling,roll = check_rolling(vres)
        
        #Loop over seasons:
        #------------------
        for season in seasons:
            
            #Initialize nested dictionary for each season
            var_lev_dict[var][season] = {}
            
            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot(111)

            #Set up list to gather whether the var exists for each case (or if var has vertical levs for now)
            #This will close the fig if neither case has a variable to plot
            #TODO: there might be a better way of doing this but becasue of the nested for-loops it's a work around for now - JR
            bad = []
            
            #Initialize dictionary to keep track of climo years for all cases
            yrs = {}

            #Loop over test cases:
            #----------------------
            for case_idx, case_name in enumerate(all_case_names):
                
                #Initialize nested dictionary for each case
                var_lev_dict[var][season][case_name] = {}                
                
                #Locate the time series files
                input_location = Path(case_ts_locs[case_idx])
                ts_filenames = f'{case_name}.*.{var}.*nc'
                ts_files = sorted(input_location.glob(ts_filenames))
                
                # If no files exist, try to move to next variable. 
                # --> Means we can not proceed with this variable, and it'll be problematic later.
                if not ts_files:
                    if season == seasons[0]:
                        errmsg = f"Time series files for variable '{var}' not found.  Script will continue to next variable."
                        warnings.warn(errmsg)
                    bad.append(True)
                    continue
                #End if

                #TEMPORARY:  For now, make sure only one file exists:
                if len(ts_files) != 1:
                    errmsg =  "Currently the time series plotting script can only handle one time series file per variable."
                    errmsg += f" Multiple files were found for the variable '{var}'"
                    raise AdfError(errmsg)
                #End if

                # Load the data
                data = _load_data(ts_files[0], var)

                #Extract units string, if available:
                if hasattr(data, 'units'):
                    unit_str = data.units
                else:
                    unit_str = '--'
                #End if

                #Check if variable has a vertical coordinate:
                if 'lev' in data.coords or 'ilev' in data.coords:                                            
                    #Skip this variable and move to the next variable in var_list
                    # during 2-d plotting. 
                    bad.append(True)
                    continue
                #End if

                #Check if variable should be masked:
                if 'mask' in var_default_dict:
                    if var_default_dict['mask'].lower() == 'ocean':
                        #Check if the ocean fraction has already been regridded
                        #and saved:
                        if ocn_frc_da[case_name] is not None:
                            ofrac = ocn_frc_da[case_name]
                            # set the bounds of regridded ocnfrac to 0 to 1
                            ofrac = xr.where(ofrac>1,1,ofrac)
                            ofrac = xr.where(ofrac<0,0,ofrac)

                            # apply ocean fraction mask to variable
                            #print(data.time,ofrac.time)
                            data = pf.mask_land_or_ocean(data, ofrac, use_nan=True)
                        else:
                            print(f"OCNFRAC not found, unable to apply mask to '{var}'")
                        #End if
                    else:
                        #Currently only an ocean mask is supported, so print warning here:
                        wmsg = "Currently the only variable mask option is 'ocean',"
                        wmsg += f"not '{var_default_dict['mask'].lower()}'"
                        print(wmsg)
                    #End if
                #End if

                #If the variable is ocean fraction, then save the dataset for use later:
                if var == 'OCNFRAC':
                    ocn_frc_da[case_name] = data
                #End if

                # we should check if we need to do area averaging:
                if len(data.dims) > 1:
                    # flags that we have spatial dimensions
                    # Note: that could be 'lev' which should trigger different behavior
                    # Note: we should be able to handle (lat, lon) or (ncol,) cases, at least
                    data_sp_avg = pf.spatial_average(data)  # changes data "in place"
                #End if

                #Nicknames for plot legend
                name = all_nicknames[case_idx]

                #Set the baseline plot line as green dashed
                #TODO: change for color deficiency - JR
                if case_name == data_name:
                    color_dict = {"color":'g',"marker":"--.",
                                  "label":f"{name} (baseline)"}
                else:
                    color_dict = {"color":colors[case_idx],"marker":"-.",
                                  "label":f"{name}"}
                #End if

                if season == "ANN":
                    ds = pf.annual_mean(data_sp_avg, whole_years=True, time_name='time')
                else:
                    ds = pf.seasonal_mean(data_sp_avg, season=season, is_climo=False)
                    ds = ds.groupby('time.year').mean(dim='time')
                #End if

                if rolling:
                    if len(ds.year) > roll:
                        ds = ds.rolling(year=roll,center=True).mean().dropna("year")
                    else:
                        msg = f"'{var}' doesn't have enough years for the chosen rolling average."
                        msg += " No rolling average applied."
                        adfobj.debug_log(msg)
                        rolling = False

                #Gather years to plot
                #yrs = ds.year
                yrs[case_name] = ds.year

                #Add case to plot (ax)
                ax.plot(yrs[case_name], ds, color_dict["marker"], c=color_dict["color"],label=color_dict["label"])

                #For the minor ticks, use no labels; default NullFormatter.
                ax.tick_params(which='major', length=7)
                ax.tick_params(which='minor', length=5)
                ds.close()
            #End for (case names)

            #Set up plots
            plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"

            if rolling:
                #Add rolling to file name for extra info?
                ax.set_title(f"{roll}-yr rolling average",loc="right")
            if multi_case:
                #Add multi_plot to file name for ADF multi-case web generation
                plot_name = plot_name.replace(f".{plot_type}",f"_multi_plot.{plot_type}")

            #Check if any cases were flagged
            #If there is at least one case that has the var to plot, add it
            #Check if against obs (case_num = 1)
            if (case_num == 1 and len(bad) < case_num) or (case_num > 1 and len(bad) < case_num):
                if var not in skip_list:
                    #Set Main title for subplots:
                    ax.set_title(f"Time Series {title_var}: {var} - {season}",loc="left")
                    
                    #Grab first and last years from each case to determine range
                    #of x-axis -> want to encompass all possible years
                    first_yrs = []
                    last_yrs = []
                    for case,years in yrs.items():
                        first_yrs.append(min(years.values))
                        last_yrs.append(max(years.values))
                    
                    #Set range based off earliest and latest years of all cases involved
                    yrs_cleaned = np.arange(min(first_yrs),max(last_yrs)+1,1)
                    
                    #Format axes
                    ax = _format_yaxis(ax, case_num, unit_str, **vres)
                    ax = _format_xaxis(ax, yrs_cleaned)

                    #Set up legend
                    fig = _make_fig_legend(case_num, fig)
                    plt.savefig(plot_name, facecolor='w')

                    #Add plot to website (if enabled):
                    adfobj.add_website_data(plot_name, f"{var}", all_case_names[0], category=web_category,
                                            season=season, plot_type="TimeSeries")
            #End if (plotting for good vars - vs obs)

            #Close the figure
            plt.close()

        #End for (seasons)
    #End for (variables)
    
    #Notify user the plots are finished
    print("  ...time series plots have been generated successfully.")



#Helper functions
#----------------

def save_to_nc(tosave, outname, attrs=None, proc=None):
    """Saves xarray variable to new netCDF file"""

    xo = tosave  # used to have more stuff here.
    # deal with getting non-nan fill values.
    if isinstance(xo, xr.Dataset):
        enc_dv = {xname: {'_FillValue': None} for xname in xo.data_vars}
    else:
        enc_dv = {}
    #End if
    enc_c = {xname: {'_FillValue': None} for xname in xo.coords}
    enc = {**enc_c, **enc_dv}
    if attrs is not None:
        xo.attrs = attrs
    if proc is not None:
        xo.attrs['Processing_info'] = f"Start from file {origname}. " + proc
    xo.to_netcdf(outname, format='NETCDF4', encoding=enc)

#####

def _set_ymargin(ax, top, bottom):
    """
    Allow for custom padding of plot lines and axes borders
    -----
    """
    ax.set_ymargin(0)
    ax.autoscale_view()
    lim = ax.get_ylim()
    delta = np.diff(lim)
    top = lim[1] + delta*top
    bottom = lim[0] - delta*bottom
    ax.set_ylim(bottom,top)
    return ax

########

def _format_yaxis(ax, case_num, unit, **kwargs):
    """
    Gather variable data and format y-axis
    -----
        - Set the y-axis plot limits to guarantee data range from all cases (including baseline)
        - Pad the top of plot to allow for flexible-sized legend in top left corner
            -> For multi-case, this will pad the plot according to number of cases
    """

    #Set up plot details, if applicable from the adf_variable_defaults.yaml file
    if 'ts' in kwargs:
        if "units" in kwargs['ts']:
            print("Looks like desired units are different than from raw file...\n")
            unit = kwargs['ts']["units"]

    ax.set_ylabel(unit,fontsize=20,labelpad=12)

    #Attempt flexible pad based on number of cases for both single and
    #multi-case scenarios, too
    pad = 0.075*case_num
    ax = _set_ymargin(ax, top=pad, bottom=0.1)
    
    #ax.yaxis.set_major_locator(MultipleLocator(1))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    return ax

########

def _format_xaxis(ax, yrs):
    """
    Gather climo year data and format x-axis
    -----
        - Set the x-axis plot limits to guarantee data range from all cases (including baseline)
        - Set minor and major locators based on number of years
        - Round the range to the nearest 5-year interval for cleaner appearance
    """

    #Grab all unique years and find min/max years
    uniq_yrs = sorted(yrs)
    first_year = int(uniq_yrs[0])
    last_year = int(uniq_yrs[-1])
    
    # Pad the last year by one -> just to add space on the plot?
    #last_year = max_year# + 1
    #last_year = max_year - max_year % 5
    #if (max_year > 5) and (last_year < max_year):
    #    last_year += 5
    
    # Pad the first year by one -> just to add space on the plot?
    #first_year = min_year# - 1
    #first_year = min_year - min_year % 5
    #if min_year < 5:
    #    first_year = 0

    #print(first_year, last_year)
    ax.set_xlim(first_year, last_year)
    ax.set_xlabel("Years",fontsize=15,labelpad=20)

    #x-axis ticks and numbers
    if len(uniq_yrs) > 120:
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
    if 50 <= len(uniq_yrs) <= 120:
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
    if 10 <= len(uniq_yrs) < 50:
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
    if 0 < len(uniq_yrs) < 10:
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
    
    return ax

########

def _make_fig_legend(case_num, fig):
    """
    Defualt matplotlib legend
    -----
        - This will just plot the colored lines and case names as given by the adf obj
          Function to generate legend and labels for all plots
    """
    
    #Gather labels based on case names and plotted line format (color, style, etc)
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        
    #Make height based on number of cases
    #NOTE: I think this works for multi-case, just need to test - JR
    h = 0.05*(case_num-1)
    
    fig.legend(lines[:case_num+1], labels[:case_num+1],loc="upper left",
                bbox_to_anchor=(0.12, 0.885-h, 0.042, h) #bbox_to_anchor(x0, y0, width, height)
                ) 

    return fig

########

def _load_data(dataloc, varname):
    ds = xr.open_dataset(dataloc)
    time_index_shifted  = ds.time.get_index('time') - timedelta(days=1)
    ds['time'] = time_index_shifted
    return ds[varname]

########

def check_rolling(vres):
    """
    Search variable defaults file for any rolling mean desired for variable

    Currently, this is only set up for rolling average of years.
    """

    rolling = False
    roll = None

    if vres.get('timeseries'):
        if 'rolling' in vres["timeseries"]:
            # roll_interval -> years, months, etc
            #NOTE: start with years only for time series for now
            if 'years' in vres["timeseries"]['rolling']:
                rolling = True
                roll = vres['timeseries']["rolling"]["years"]
                print(f"\t   rolling interval: {roll} years\n")
            #End if
        #End if
    #End if
    
    return rolling,roll

########

def seasonal_mean(data, season=None, is_climo=None):
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
    seasons = {
            "DJF": [12, 1, 2],
            "JJA": [6, 7, 8],
            "MAM": [3, 4, 5],
            "SON": [9, 10, 11]}
    if season is not None:
        assert season in ["ANN", "DJF", "JJA", "MAM", "SON"], f"Unrecognized season string provided: '{season}'"
    elif season is None:
        season = "ANN"

    try:
        month_length = data.time.dt.days_in_month
    except (AttributeError, TypeError):
        # do our best to determine the temporal dimension and assign weights
        if not is_climo:
            raise ValueError("Non-climo file provided, but without a decoded time dimension.")
        else:
            # CLIMO file: try to determine which dimension is month
            has_time = False
            if isinstance(data, xr.DataArray):
                has_time = 'time' in data.dims
                if not has_time:
                    if "month" in data.dims:
                        data = data.rename({"month":"time"})
                        has_time = True
            if not has_time:
                # this might happen if a pure numpy array gets passed in
                # --> assumes ordered January to December.
                assert ((12 in data.shape) and (data.shape.count(12) == 1)), f"Sorry, {data.shape.count(12)} dimensions have size 12, making determination of which dimension is month ambiguous. Please provide a `time` or `month` dimension."
                time_dim_num = data.shape.index(12)
                fakedims = [f"dim{n}" for n in range(len(data.shape))]
                fakedims[time_dim_num] = "time"
                data = xr.DataArray(data, dims=fakedims)
            timefix = pd.date_range(start='1/1/1999', end='12/1/1999', freq='MS') # generic time coordinate from a non-leap-year
            data = data.assign_coords({"time":timefix})
        month_length = data.time.dt.days_in_month
    #End try/except

    data = data.sel(time=data.time.dt.month.isin(seasons[season])) # directly take the months we want based on season kwarg
    
    if not is_climo: #ie time series
        #For time series this will keep the time coordinates for calculation later
        return data
    else:
        return data.weighted(data.time.dt.daysinmonth).mean(dim='time')
