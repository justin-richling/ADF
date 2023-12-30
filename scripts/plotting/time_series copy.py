#Import standard modules:
from pathlib import Path
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker

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

     #Import necessary modules:
    #------------------------
    #CAM diagnostic plotting functions:
    import plotting_functions as pf

    #ADF warning format
    from adf_base import AdfError
    #-------------------------

    #Notify user that script has started:
    print("\n  Generating time series plots...")

    #Extract needed quantities from ADF object:
    #-----------------------------------------
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
    base_nickname = adfobj.case_nicknames["base_nickname"]

    #CAUTION:
    #"data" here refers to either obs or a baseline simulation,
    #Until those are both treated the same (via intake-esm or similar)
    #we will do a simple check and switch options as needed:
    if adfobj.get_basic_info("compare_obs"):

        #Extract variable-obs dictionary:
        var_obs_dict = adfobj.var_obs_dict
        base_nickname = "Obs"

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
            all_case_names = case_names + [data_name]
            #Gather all nicknames
            all_nicknames = test_nicknames + [base_nickname]
            case_ts_locs = case_ts_loc + [data_ts_loc]
        
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

    #ADF variable which contains the output path for plots and tables:
    plot_location = adfobj.plot_location
    if not plot_location:
        plot_location = adfobj.get_basic_info("cam_diag_plot_loc")
    if isinstance(plot_location, list):
        for pl in plot_location:
            plpth = Path(pl)
            #Check if plot output directory exists, and if not, then create it:
            if not plpth.is_dir():
                print(f"\t    '{pl}' not found, making new directory")
                plpth.mkdir(parents=True)
        if len(plot_location) == 2:
            plot_loc = Path(plot_location[0])
        else:
            print(f"Ambiguous plotting location since all cases go on same plot. Will put them in first location: {main_site_assets_path}")
            plot_loc = main_site_assets_path
    else:
        plot_loc = Path(plot_location)    

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
    #seasons = ["ANN"]

    #Set up the plots
    #################

    #Add more colors as needed for number of test cases
    #** Baseline is already added as green dashed line in plotting function **
    #matplotlib colors here: https://matplotlib.org/stable/gallery/color/named_colors.html
    colors = ["k", "aqua", "r", "b", "magenta",
              "orange", "slategrey", "rosybrown"]
    
    #Create/reset new variable that potentially stores the re-gridded
    #ocean fraction xarray data-array:
    ocn_frc_da = None
    
    #Loop over variables:
    for var in var_list[2:5]:
        
        #Extract defaults for variable:
        var_default_dict = var_defaults.get(var, {}) #pull this put and put above!
        
        #Check res for any variable specific options that need to be used BEFORE going to the plot:
        if var in res:
            vres = res[var]
            #If found then notify user, assuming debug log is enabled:
            adfobj.debug_log(f"time_series: Found variable defaults for {var}")
        else:
            vres = {}
        #End if

        title_var = "Global"
        print(f"\t - time series for {var}")

        #Set plotting parameters based off whether the user wants
        #5-yr rolling average
        #Currently RESTOM is defaulted to 5-yr rolling avg
        rolling = False
        if 'timeseries' in vres:
            if "rolling" in vres['timeseries']:
                rolling = True
                rolling_yrs = vres['timeseries']["rolling"]["years"]
                print("rolling_months",rolling_yrs,"\n")

        #Loop over test cases:
        #----------------------
        for season in seasons:
            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot(111)
            bad = []
            for case_idx, case_name in enumerate(all_case_names):
                input_location = Path(case_ts_locs[case_idx])
                ts_filenames = f'{case_name}.*.{var}.*nc'
                ts_files = sorted(input_location.glob(ts_filenames))
                
                
                # If no files exist, try to move to next variable. --> Means we can not proceed with this variable, and it'll be problematic later.
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

                #Check if variable has a vertical coordinate:
                if 'lev' in data.coords or 'ilev' in data.coords:
                    #TODO: calc for the desired level in config file if var has vertical levels - JR
                    #
                    #
                    #
                    #
                    
                    # this warning will appear for each case that the variable is missing...
                    if season == seasons[0]:
                        print(f"\t   Variable '{var}' has a vertical dimension, "+\
                              "which is currently not supported for the time series plot. Skipping...")
                    #Skip this variable and move to the next variable in var_list:
                    plt.close()
                    bad.append(True)
                    continue
                #End if

                #Extract defaults for variable:
                var_default_dict = var_defaults.get(var, {}) #pull this put and put above!

                #Check if variable should be masked:
                if 'mask' in var_default_dict:
                    if var_default_dict['mask'].lower() == 'ocean':
                        #Check if the ocean fraction has already been regridded
                        #and saved:
                        if ocn_frc_da is not None:
                            ofrac = ocn_frc_da
                            # set the bounds of regridded ocnfrac to 0 to 1
                            ofrac = xr.where(ofrac>1,1,ofrac)
                            ofrac = xr.where(ofrac<0,0,ofrac)

                            # apply ocean fraction mask to variable
                            data = mask_land_or_ocean(data, ofrac, use_nan=True)
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
                    ocn_frc_da = data
                #End if

                # we should check if we need to do area averaging:
                if len(data.dims) > 1:
                    # flags that we have spatial dimensions
                    # Note: that could be 'lev' which should trigger different behavior
                    # Note: we should be able to handle (lat, lon) or (ncol,) cases, at least
                    data_sp_avg = spatial_average(data)  # changes data "in place"

                if case_name == data_name:
                    color_dict = {"color":'g',"marker":"--"}
                else:
                    color_dict = {"color":colors[case_idx],"marker":"-"}
                #End if

                if season == "ANN":
                    ds = annual_mean(data_sp_avg, whole_years=True, time_name='time')
                else:
                    ds = seasonal_mean(data_sp_avg, season=season, is_climo=False)
                    ds = ds.groupby('time.year').mean(dim='time')
                    
                if rolling:
                    print("ROLLING IS WORKING:",rolling_months,"months")
                    #avg_case = avg_case.rolling(time=rolling_months,center=True).mean()
                    ds = ds.rolling(time=rolling_months,center=True).mean()

                name = all_nicknames[case_idx]
                if case_idx == (case_num-1):
                    name = f"{name} (baseline)"
                
                #Gather years to plot
                yrs = ds.year
                
                #Add case to plot (ax)
                ax.plot(yrs, ds, color_dict["marker"], c=color_dict["color"],label=name)

                #For the minor ticks, use no labels; default NullFormatter.
                ax.tick_params(which='major', length=7)
                ax.tick_params(which='minor', length=5)
                ds.close()
            #End for (case names)
            

            #print(var,bad)
            if len(bad) < 2:
                #Set Main title for subplots:
                ax.set_title(f"Time Series {title_var}: {var} - {season}",loc="left")

                if rolling:
                    ax.set_title(f"5-yr rolling average",loc="right")
                unit = unit_str

                ax = _format_yaxis(ax, case_num, unit, **vres)
                ax = _format_xaxis(ax, yrs)

                #Set up legend
                fig = _make_fig_legend(case_num, fig)

                #Save plot
                #plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"

                #Save plot
                #plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"
                #if multi_case:
                if 0==1:
                    plot_name = plot_loc / f"{var}_{season}_TimeSeries_multi_plot.{plot_type}"
                else:
                    #plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"
                    plot_name = f"./{var}_{season}_TimeSeries_Mean.{plot_type}"

                plt.savefig(plot_name, facecolor='w')
            #    plt.close()
            """
            else:
                # if no case has a plotted variable then close the figure so 
                #plt.close()
                if season == seasons[0]: 
                    print(f"Time series plot for {var} not created, sorry :(")
            """
            plt.close()
        #End for (seasons)
    #End for (variables)


# Helper functions
##################

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
    #uniq_yrs = sorted({x for v in yrs.values() for x in v})
    #uniq_yrs = sorted(x for v in yrs.values() for x in v)
    uniq_yrs = yrs #sorted(x for v in yrs for x in v)
    max_year = int(max(uniq_yrs))
    min_year = int(min(uniq_yrs))
    
    # Pad the last year by one -> just to add space on the plot?
    last_year = max_year# + 1
    #last_year = max_year - max_year % 5
    #if (max_year > 5) and (last_year < max_year):
    #    last_year += 5
    
    # Pad the first year by one -> just to add space on the plot?
    first_year = min_year# - 1
    #first_year = min_year - min_year % 5
    #if min_year < 5:
    #    first_year = 0

    #print(first_year, last_year)
    ax.set_xlim(first_year, last_year)
    ax.set_xlabel("Years",fontsize=15,labelpad=20)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    
    """
    #x-axis ticks and numbers
    if max_year > 120:
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
    if 10 <= max_year <= 120:
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
    if 0 < max_year < 10:
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
    """
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

    fig.legend(lines[:case_num+1], labels[:case_num+1],loc="upper left",
                bbox_to_anchor=(0.12, 1.025,.042,.05*(case_num-1))
                #bbox_to_anchor=(0.12, 0.835,.042,.05)
                ) #bbox_to_anchor(x0, y0, width, height)

    return fig

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
    
    if not is_climo:
        return data
    else:
        return data.weighted(data.time.dt.daysinmonth).mean(dim='time')

########

def annual_mean(data, whole_years=False, time_name='time'):
    """Calculate annual averages from monthly time series data.

    Parameters
    ----------
    data : xr.DataArray or xr.Dataset
        monthly data values with temporal dimension
    whole_years : bool, optional
        whether to restrict endpoints of the average to
        start at first January and end at last December
    time_name : str, optional
        name of the time dimension, defaults to `time`

    Returns
    -------
    result : xr.DataArray or xr.Dataset
        `data` reduced to annual averages

    Notes
    -----
    This function assumes monthly data, and weights the average by the
    number of days in each month.

    `result` includes an attribute that reports the date range used for the average.
    """
    assert time_name in data.coords, f"Did not find the expected time coordinate '{time_name}' in the data"
    if whole_years:
        first_january = np.argwhere((data.time.dt.month == 1).values)[0].item()
        last_december = np.argwhere((data.time.dt.month == 12).values)[-1].item()
        data_to_avg = data.isel(time=slice(first_january,last_december+1)) # PLUS 1 BECAUSE SLICE DOES NOT INCLUDE END POINT
    else:
        data_to_avg = data
    date_range_string = f"{data_to_avg['time'][0]} -- {data_to_avg['time'][-1]}"

    # this provides the normalized monthly weights in each year
    # -- do it for each year to allow for non-standard calendars (360-day)
    # -- and also to provision for data with leap years
    days_gb = data_to_avg.time.dt.daysinmonth.groupby('time.year').map(lambda x: x / x.sum())
    # weighted average with normalized weights: <x> = SUM x_i * w_i  (implied division by SUM w_i)
    result =  (data_to_avg * days_gb).groupby('time.year').sum(dim='time')
    result.attrs['averaging_period'] = date_range_string
    return result

########

def _load_data(dataloc, varname):
    import xarray as xr
    ds = xr.open_dataset(dataloc)
    time_index_shifted  = ds.time.get_index('time') - timedelta(days=1)
    ds['time'] = time_index_shifted
    return ds[varname]

########

def mask_land_or_ocean(arr, msk, use_nan=False):
    """Apply a land or ocean mask to provided variable.

    Parameters
    ----------
    arr : xarray.DataArray
        the xarray variable to apply the mask to.
    msk : xarray.DataArray
        the xarray variable that contains the land or ocean mask,
        assumed to be the same shape as "arr".
    use_nan : bool, optional
        argument for whether to set the missing values
        to np.nan values instead of the defaul "-999." values.

    Returns
    -------
    arr : xarray.DataArray
        Same as input `arr` but masked as specified.
    """

    if use_nan:
        missing_value = np.nan
    else:
        missing_value = -999.
    #End if

    arr = xr.where(msk>=0.9,arr,missing_value)
    arr.attrs["missing_value"] = missing_value
    return(arr)

########

def spatial_average(indata, weights=None, spatial_dims=None):
    """Compute spatial average.

    Parameters
    ----------
    indata : xr.DataArray
        input data
    weights : np.ndarray or xr.DataArray, optional
        the weights to apply, see Notes for default behavior
    spatial_dims : list, optional
        list of dimensions to average, see Notes for default behavior

    Returns
    -------
    xr.DataArray
        weighted average of `indata`

    Notes
    -----
    When `weights` is not provided, tries to find sensible values.
    If there is a 'lat' dimension, use `cos(lat)`.
    If there is a 'ncol' dimension, looks for `area` in `indata`.
    Otherwise, set to equal weights.

    Makes an attempt to identify the spatial variables when `spatial_dims` is None.
    Will average over `ncol` if present, and then will check for `lat` and `lon`.
    When none of those three are found, raise an AdfError.
    """
    import warnings

    if weights is None:
        #Calculate spatial weights:
        if 'lat' in indata.coords:
            weights = np.cos(np.deg2rad(indata.lat))
            weights.name = "weights"
        elif 'ncol' in indata.dims:
            if 'area' in indata:
                warnings.warn("area variable being used to generated normalized weights.")
                weights = indata['area'] / indata['area'].sum()
            else:
                warnings.warn("We need a way to get area variable. Using equal weights.")
                weights = xr.DataArray(1.)
            weights.name = "weights"
        else:
            weights = xr.DataArray(1.)
            weights.name = "weights"
            warnings.warn("Un-recognized spatial dimensions: using equal weights for all grid points.")
        #End if
    #End if

    #Apply weights to input data:
    weighted = indata.weighted(weights)

    # we want to average over all non-time dimensions
    if spatial_dims is None:
        if 'ncol' in indata.dims:
            spatial_dims = ['ncol']
        else:
            spatial_dims = [dimname for dimname in indata.dims if (('lat' in dimname.lower()) or ('lon' in dimname.lower()))]

    if not spatial_dims:
        #Scripts using this function likely expect the horizontal dimensions
        #to be removed via the application of the mean. So in order to avoid
        #possibly unexpected behavior due to arrays being incorrectly dimensioned
        #(which could be difficult to debug) the ADF should die here:
        emsg = "spatial_average: No spatial dimensions were identified,"
        emsg += " so can not perform average."
        raise AdfError(emsg)

    return weighted.mean(dim=spatial_dims)

########

def global_average(fld, wgt, verbose=False):
    """A simple, pure numpy global average.

    Parameters
    ----------
    fld : np.ndarray
        an input ndarray
    wgt : np.ndarray
        a 1-dimensional array of weights, should be same size as one dimension of `fld`
    verbose : bool, optional
        prints information when `True`

    Returns
    -------
    weighted average of `fld`
    """

    s = fld.shape
    for i in range(len(s)):
        if np.size(fld, i) == len(wgt):
            a = i
            break
    fld2 = np.ma.masked_invalid(fld)
    if verbose:
        print("(global_average)-- fraction of mask that is True: {}".format(np.count_nonzero(fld2.mask) / np.size(fld2)))
        print("(global_average)-- apply ma.average along axis = {} // validate: {}".format(a, fld2.shape))
    avg1, sofw = np.ma.average(fld2, axis=a, weights=wgt, returned=True) # sofw is sum of weights

    return np.ma.average(avg1)