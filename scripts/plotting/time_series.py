from pathlib import Path
from typing import OrderedDict
import numpy as np
import xarray as xr
import pandas as pd
import plotting_functions as pf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
import warnings  #use to warn user about missing files.

def my_formatwarning(msg, *args, **kwargs):
    #ignore everything except the message
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

    #Notify user that script has started:
    print("\n  Generating time series plots...")

    #Extract needed quantities from ADF object:
    #-----------------------------------------

    #Timeseries yaml file options
    ts_opts = adfobj.read_config_var("time_series_var_list")
    if not ts_opts:
        print(
            "Time series options were not specified, so time_series can not run."+\
            " See documentation or config_cam_baseline_example.yaml for options to add to configuration file."
        )
        return

    ann_var_list = ts_opts.get("ann_var_list", adfobj.diag_var_list)
    season_var_list = ts_opts.get("season_var_list", adfobj.diag_var_list)

    case_names = adfobj.get_cam_info('cam_case_name', required=True)

    #Check for multi-case diagnostics
    if len(case_names) > 1:
        #Grab all multi-case diagnostic directory
        main_site_assets_path = adfobj.main_site_paths["main_site_assets_path"]
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
            print("No observations found to plot against, so no time series plots will be generated.")
            return

    else:
        data_name = adfobj.get_baseline_info("cam_case_name")
        data_ts_loc = adfobj.get_baseline_info("cam_ts_loc")
    #End if

    #Bundle all case names
    all_case_names = case_names + [data_name]
    case_num = len(all_case_names)

    #Gather all nicknames
    all_nicknames = test_nicknames + [base_nickname]

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

    #Set seasonal ranges
    seasons = ["DJF","MAM","JJA","SON"]

    #Set up the plots
    #################

    #Add more colors as needed for number of test cases
    #** Baseline is already added as green dashed line in plotting function **
    #matplotlib colors here: https://matplotlib.org/stable/gallery/color/named_colors.html
    colors = ["k", "aqua", "r", "b", "magenta",
              "orange", "slategrey", "rosybrown"]

    case_ts_locs = case_ts_loc + [data_ts_loc]

    #Grab all the seasonally weighted data up front
    print("\n  Grabbing seasonally weighted data...")
    vals, yrs, units, season_var_list = _get_seasonal_data(season_var_list, all_case_names, case_ts_locs)
    
    #Annual global weighted
    #######################
    #Treat ANN differently than seasonally weighted (below)
    season = "ANN"
    print(f"\n  Generating time series for {season}...")
    #Loop over variables:
    for var in ann_var_list:
        #Check res for any variable specific options that need to be used BEFORE going to the plot:
        if var in res:
            vres = res[var]
            #If found then notify user, assuming debug log is enabled:
            adfobj.debug_log(f"time_series: Found variable defaults for {var}")
        else:
            vres = {}
        #End if

        fig = plt.figure(figsize=(12,8))
        ax = fig.add_subplot(111)

        title_var = "Global"
        print(f"\t - time series for {var}")

        #Set plotting parameters based off whether the user wants
        #5-yr rolling average
        #Currently RESTOM is defaulted to 5-yr rolling avg
        rolling = False
        if 'ts' in vres:
            if "rolling" in vres['ts']:
                rolling = True
                rolling_months = vres['ts']["rolling"]["months"]
                print("rolling_months",rolling_months,"\n")

        #Loop over test cases:
        #----------------------
        #Create lists to hold all min/max values for var data (for each case)
        mins = []
        maxs = []

        for case_idx, case_name in enumerate(all_case_names):
            if var == "RESTOM":
                fils = sorted(Path(case_ts_locs[case_idx]).glob(f"*FSNT*.nc"))
                ts_ds = _load_dataset(fils)
                avg_case_FSNT,_,yrs_case,unit = _data_calcs('FSNT',ts_ds=ts_ds,subset=None)
                fils = sorted(Path(case_ts_locs[case_idx]).glob(f"*FLNT*.nc"))
                ts_ds = _load_dataset(fils)
                avg_case_FLNT,_,_,_ = _data_calcs("FLNT",ts_ds=ts_ds,subset=None)
                if len(yrs_case) < 5:
                    print(f"Not a lot of climo years for {case_name}, only doing 1-yr avg for RESTOM...")
                    FSNT_case = avg_case_FSNT
                    FLNT_case = avg_case_FLNT
                    if case_name == data_name:
                        color_dict = {"color":'g',"marker":"--*"}
                    else:
                        color_dict = {"color":colors[case_idx],"marker":"-*"}
                else:
                    FSNT_case = avg_case_FSNT.rolling(time=rolling_months,center=True).mean()
                    FLNT_case = avg_case_FLNT.rolling(time=rolling_months,center=True).mean()
                    if case_name == data_name:
                        color_dict = {"color":'g',"marker":"--"}
                    else:
                        color_dict = {"color":colors[case_idx],"marker":"-"}
                #End if
                avg_case = FSNT_case - FLNT_case
            else:
                if case_name == data_name:
                    color_dict = {"color":'g',"marker":"--"}
                else:
                    color_dict = {"color":colors[case_idx],"marker":"-"}
                #End if

                fils = sorted(Path(case_ts_locs[case_idx]).glob(f"*{var}.*.nc"))
                ts_ds = _load_dataset(fils)

                #Check if variable has a vertical coordinate:
                if 'lev' in ts_ds.coords or 'ilev' in ts_ds.coords:
                    print(f"\t   Variable '{var}' has a vertical dimension, "+\
                        "which is currently not supported for the time series plot. Skipping...")
                    #Skip this variable and move to the next variable in var_list:
                    continue
                avg_case,_,yrs_case,unit = _data_calcs(var,ts_ds=ts_ds,subset=None)
                if rolling:
                    avg_case = avg_case.rolling(time=rolling_months,center=True).mean()

            #End if (RESTOM)

            #Get yearly averages for all available years
            vals_case = [avg_case.sel(time=i).mean() for i in yrs_case]

            #Grab min and max vals from each test case
            mins.append(np.nanmin(vals_case))
            maxs.append(np.nanmax(vals_case))

            #Get int of years for plotting on x-axis
            yrs_case_int = yrs_case.astype(int)

            name = all_nicknames[case_idx]
            if case_idx == (case_num-1):
                name = f"{name} (baseline)"

            #Add case to plot (ax)
            ax.plot(yrs_case_int, vals_case, color_dict["marker"], c=color_dict["color"],label=name)

            #For the minor ticks, use no labels; default NullFormatter.
            ax.tick_params(which='major', length=7)
            ax.tick_params(which='minor', length=5)
        #End for (case names)

        #Set Main title for subplots:
        ax.set_title(f"Time Series {title_var}: {var} - {season}",loc="left")
        
        if rolling:
            ax.set_title(f"5-yr rolling average",loc="right")

        #Minor tweak to not plot variables that have vertical levels.
        #TODO: Clean up this check - JR
        #mins and maxs are blank lists when trying to ignore variables...
        if (mins) or (maxs):
            ax = _format_yaxis(ax, case_num, unit, **vres)
            ax = _format_xaxis(ax, yrs)

            #Set up legend
            lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
            lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

            fig.legend(lines[:case_num+1], labels[:case_num+1],loc="upper left",
                        bbox_to_anchor=(0.12, 0.835,.042,.05)) # bbox_to_anchor(x0, y0, width, height)

            #Save plot
            if multi_case:
                plot_name = plot_loc / f"{var}_{season}_TimeSeries_multi_plot.{plot_type}"
            else:
                plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"

            plt.savefig(plot_name, facecolor='w')
            plt.close()
            #Add plot to website (if enabled):
            adfobj.add_website_data(plot_name, var, case,
                                    season=season,
                                    plot_type="TimeSeries",
                                    multi_case=multi_case)
        #End if (min/max)

    #End for (vars)
    #End ANN plots

    #Seasonally weighted plots
    # - DJF, MAM, JJA, SON
    ##########################
    print("\n  Generating seasonally weighted time series...")
    for var in season_var_list:
        print(f"\t - time series for {var}")

        vres = res[var]
        for season in seasons:
            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot(111)
                
            #make list to capture max/min for each case
            y_mins = []
            y_maxs = []
            for case_idx, case_name in enumerate(all_case_names):
                    
                #Check for baseline, and set linestyle to dashed
                if case_name == data_name:
                    label=f"{base_nickname} (baseline)"
                    marker = "--"
                    ax.plot(yrs[case_name].astype(int),
                            vals[var][case_name][season],
                            marker, c='g', label=label)
                #Set linestyle for test cases to solid line
                else:
                    label=f"{test_nicknames[case_idx]}"
                    marker = "-"
                    ax.plot(yrs[case_name].astype(int),
                            vals[var][case_name][season],
                            marker, c=colors[case_idx],label=label)
                #End if

                #For the minor ticks, use no labels; default NullFormatter.
                ax.tick_params(which='major', length=7)
                ax.tick_params(which='minor', length=5)

                #Attempt to set custom y-ranges
                #Grab mins/maxes
                y_mins.append(np.nanmin(vals[var][case_name][season]))
                y_maxs.append(np.nanmax(vals[var][case_name][season]))
            #End for (cases)

            #Set Main title for subplots:
            ax.set_title(f"Time Series {title_var}: {var} - {season}",loc="left")

            #Format axes
            ax = _format_xaxis(ax, yrs)
            ax = _format_yaxis(ax, case_num, units[var], **vres)

            #Set up legend
            lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
            lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

            fig.legend(lines[:case_num+1], labels[:case_num+1],loc="upper left",
                        bbox_to_anchor=(0.12, 0.835,.042,.05)) # bbox_to_anchor(x0, y0, width, height)

            #Save plot
            if multi_case:
                plot_name = plot_loc / f"{var}_{season}_TimeSeries_multi_plot.{plot_type}"
            else:
                plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"

            plt.savefig(plot_name, facecolor='w')
            plt.close()

            #Add plot to website (if enabled):
            adfobj.add_website_data(plot_name, var, case, season=season,
                                    plot_type="TimeSeries",
                                    multi_case=multi_case)
            #End for (cases)
        #End for (season)
    #End for (var)
    #End seasonally weighted

    #Notify user that script has ended:
    print("  ...time series have been generated successfully.")

#Helper functions:
##################

def _load_dataset(fils):
    if len(fils) == 0:
        print("Input file list is empty.")
        return None
    elif len(fils) > 1:
        return xr.open_mfdataset(fils, combine='by_coords')
    else:
        sfil = str(fils[0])
        return xr.open_dataset(sfil)
    #End if
#End def

########

def _data_calcs(var, ts_ds, subset=None):
    """
    Function to calculate global weighted data
    ----
        - weight global data by latitude
        - option for lat/lon subset
        - returns
            * month lengths
            * climo years
            * units
    """

    time = ts_ds['time']
    time = xr.DataArray(ts_ds['time_bnds'].load().mean(dim='nbnd').values, dims=time.dims, attrs=time.attrs)
    ts_ds['time'] = time
    ts_ds.assign_coords(time=time)
    ts_ds = xr.decode_cf(ts_ds)

    if subset != None:
        ts_ds = ts_ds.sel(lat=slice(subset["s"],subset["n"]), lon=slice(subset["w"],subset["e"])) 

    data = ts_ds[var].squeeze()
    month_length = data.time.dt.days_in_month
    unit = data.units

    #global weighting
    w = np.cos(np.radians(data.lat))
    avg = data.weighted(w).mean(dim=("lat","lon"))

    yrs = np.unique([str(val.item().timetuple().tm_year).zfill(4) for val in ts_ds["time"]])

    return avg,month_length,yrs,unit

########

def seasonal_data(data, month_length):
    """
    Function to grab seasonal weighted data, grouped by season
    -----
        - DJF, MAM, JJA, and SON
    """

    weighted_mean = ((data * month_length).resample(time='QS-DEC').sum() /
                          month_length.resample(time='QS-DEC').sum())

    mdata_seasonal_mean = weighted_mean.groupby('time.season')

    return mdata_seasonal_mean

########

def _get_seasonal_data(season_var_list, all_case_names, case_ts_locs):
    """
    Gather seasonally weighted data
    -----
    - DFJ, MAM, JJA, SON 
        ** No ANN **
    """

    #Ignore variables that have vertical levels for now
    season_var_list_mod = []

    vals = OrderedDict()
    yrs = {}
    units = {}
    for var in season_var_list:
        #print("var seaonal_data:",var,"\n")
        if var not in vals:
            vals[var] = OrderedDict()
        for case_idx, case_name in enumerate(all_case_names):
            fils = sorted(Path(case_ts_locs[case_idx]).glob(f"*{var}.*.nc"))
            ts_ds = _load_dataset(fils)

            if ts_ds:
                if ('lev' not in ts_ds.coords) or ('ilev' not in ts_ds.coords):
                    season_var_list_mod.append(var)
                    data,month_length,_,unit =_data_calcs(var,ts_ds=ts_ds,subset=None)
                    units[var] = unit
                    mdata_seasonal_mean = seasonal_data(data, month_length)

                    if case_name not in vals[var]:
                        vals[var][case_name] = OrderedDict()

                    for season, arr in mdata_seasonal_mean:
                        #Weight DJF differently:
                        if season == "DJF":
                            #Get years
                            yrs_s = []
                            #Grab first year avialable for December
                            for val in arr.time[1:]:
                                yr = val.item().timetuple().tm_year
                                #Need to put year in 4-digit expression
                                #so grab leading zeros (ie yr 31 would be 0031)
                                yrs_s.append(str(yr).zfill(4))
                            yrs_s = np.unique(yrs_s)
                            yrs[case_name] = yrs_s

                            yrs_mean = [arr.sel(time=i).mean().values for i in yrs[case_name]]
                            vals[var][case_name][season] = yrs_mean

                        else:
                            #Get years
                            yrs_s = []
                            for val in arr.time:
                                yr = val.item().timetuple().tm_year
                                yrs_s.append(str(yr).zfill(4))
                            yrs_s = np.unique(yrs_s)
                            yrs[case_name] = yrs_s

                            yrs_mean = [arr.sel(time=i).mean().values for i in yrs[case_name]]
                            vals[var][case_name][season] = yrs_mean
                else:
                    print(f"{var} has vertical levels, skipping for now...")    
            else:
                print(f"{var} not found, skipping...")
    return vals, yrs, units, season_var_list_mod

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
            unit = kwargs['ts']["units"]

    ax.set_ylabel(unit,fontsize=20,labelpad=12)

    #Attempt flexible pad based on number of cases for both single and
    #multi-case scenarios, too
    pad = 0.075*case_num
    ax.set_ymargin(0)
    ax.autoscale_view()
    lim = ax.get_ylim()
    delta = np.diff(lim)
    top = lim[1] + delta*pad
    bottom = lim[0] - delta*0.1
    ax.set_ylim(bottom,top)

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
    uniq_yrs = sorted(x for v in yrs.values() for x in v)
    max_year = int(max(uniq_yrs))
    min_year = int(min(uniq_yrs))

    last_year = max_year - max_year % 5
    if (max_year > 5) and (last_year < max_year):
        last_year += 5

    first_year = min_year - min_year % 5
    if min_year < 5:
        first_year = 0

    ax.set_xlim(first_year, last_year)
    ax.set_xlabel("Years",fontsize=15,labelpad=20)

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

    return ax

########

##############
#END OF SCRIPT