from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
import plotting_functions as pf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D

import time

import warnings  # use to warn user about missing files.

def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'

warnings.formatwarning = my_formatwarning

#Set seasonal ranges:
seasons = {"ANN": np.arange(1,13,1),
            "DJF": [12, 1, 2],
            "JJA": [6, 7, 8],
            "MAM": [3, 4, 5],
            "SON": [9, 10, 11]}

def time_series(adfobj):
    """
    This script plots time series.
    Compare all CAM runs against
    other climatological data (observations or baseline runs).
    Description of needed inputs from ADF:

    adfobj
    ------
        * ADF object

    Notes:
        This script can have a limited option for annual/seasonal weighting.
        It will be pretty flexible for the variables plotted and layout of figure.
    """

    #Notify user that script has started:
    print("\n  Generating time series plots...")

    # Extract needed quantities from ADF object:
    # -----------------------------------------
    
    case_names = adfobj.get_cam_info('cam_case_name', required=True)
    data_name = adfobj.get_baseline_info('cam_case_name', required=True)

    if len(case_names) > 1:
        multi_path = Path(adfobj.get_basic_info('cam_diag_plot_loc', required=True))
        main_site_path = multi_path / "main_website"
        main_site_path.mkdir(exist_ok=True)
        main_site_assets_path = main_site_path / "assets"
        main_site_assets_path.mkdir(exist_ok=True)

    case_ts_loc = adfobj.get_cam_info("cam_ts_loc", required=True)
    data_ts_loc = adfobj.get_baseline_info("cam_ts_loc", required=True)

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
    #season = "ANN"
    plot_type = "png"
    
    
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

    #Get list of plotted variables
    #ts_var_list = adfobj.diag_var_list
    ts_var_list = adfobj.timeseries_var_list

    #Subset - optional
    subsetz = {}

    """subset_vars = ["ICEFRAC"]
    for i in subset_vars:
        #w = res[i]["LabSea"]["w"]
        #e = res[i]["LabSea"]["e"]
        #s = res[i]["LabSea"]["s"]
        #n = res[i]["LabSea"]["n"]
        subset_dict = {"s":res[i]["LabSea"]["s"], "n":res[i]["LabSea"]["n"], 
                       "e":res[i]["LabSea"]["e"], "w":res[i]["LabSea"]["w"]}

        subsetz[i] = [subset_dict,"Subset"]"""

    custom_leg = True
    labels = ["022c: ice/ocn spunup + large RESTOM asegdosjdvp",
              "case",
              "othe4r case",
              "adkncakdsc",
              "022: Baseline"]

    # Add more colors as needed for number of test cases
    # ** Baseline is already added as green dashed line in plotting function **
    # matplotlib colors here: https://matplotlib.org/stable/gallery/color/named_colors.html
    colors = ["k", "aqua", "r", "b", "magenta", "orange", "slategrey", "rosybrown"]
    
    case_ts_locs = [i for i in case_ts_loc]
    case_ts_locs.append(data_ts_loc)
    
    print("Gathering data, this may take a while...")
    ticz = time.perf_counter()
    vals,yrs,units = _get_seasonal_data(ts_var_list, case_ts_locs, subsetz)
    tocz = time.perf_counter()
    print("... Did this take long? I bet it did...\n")
    print(f"\t...Seaosnal weighted calcs take {(tocz-ticz)/60:0.4f} minutes\n")

    """#Remove FSNT and FLNT from being plotted
    if "FSNT" in ts_var_list:
        ts_var_list.remove("FSNT")
    if "FLNT" in ts_var_list:
        ts_var_list.remove("FLNT")"""
    #Add derived variable names into plotting list
    #ts_var_list += ["RESTOM"]
    
    case_base_names = case_names + [data_name]

    for season in seasons:

        # Loop over variables:
        for var in ts_var_list:
            # Check res for any variable specific options that need to be used BEFORE going to the plot:
            if var in res:
                vres = res[var]
                #If found then notify user, assuming debug log is enabled:
                adfobj.debug_log(f"time_series: Found variable defaults for {var}")

                #Extract category (if available):
                web_category = vres.get("category", None)

            else:
                vres = {}
            #End if

            print(f"\t - Plotting Time Series, {season}")
                
            print("Plotting variable:",var)
                
            title_var = "Global"

            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot(111)

            ax = ts_plot(ax, var, season, units, title_var)

            # Create lists to hold all sets of years (for each case) and
            # sets of var data (for each case)
            vals_cases = []
            yrs_cases = []

            if len(case_names) > 1:
                plot_name = plot_loc / f"{var}_{season}_TimeSeries_multi_plot.{plot_type}"
            else:
                plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"

            for case_idx, case_name in enumerate(case_base_names):
                
                if case_idx == len(case_base_names)-1:
                    if custom_leg == True:
                        label=f"{labels[case_idx]} (baseline)"
                    else:
                        label=f"{case_name} (baseline)"
                    marker = "--"

                    ax.plot(yrs[case_name].astype(int), vals[var][case_name][season], marker, c='g',
                                    label=label)
                else:
                    if custom_leg == True:
                        label=f"{labels[case_idx]}"
                    else:
                        label=f"{case_name}"
                    marker = "-"
                    ax.plot(yrs[case_name].astype(int), vals[var][case_name][season], marker, c=colors[case_idx],
                                    label=label)

                vals_cases.append(vals[var][case_name][season])
                yrs_cases.append(yrs[case_name])

            ax = plot_var_details(ax, var, vals_cases, units[var], title_var, **vres)

            ax = _format_xaxis(ax, yrs_cases)

            # Set up legend
            # If custom_legend = True, change the code in make_fig_legend() function for custom legend
            #fig = make_fig_legend(case_names_len, fig, custom_legend=False)
            fig.legend(loc="center left",fontsize=12,
                                    bbox_to_anchor=(0.122, 0.82,.05,.05))
            plt.savefig(plot_name, facecolor='w')
            #Add plot to website (if enabled):
            adfobj.add_website_data(plot_name, var, case_name=None, category=web_category, season=season, plot_type="TimeSeries",multi_case=True)
            #Close plots:
            plt.close()
        # End for (variables loop)

        if multi_path:
            for season in seasons:

                # Loop over variables:
                for var in ts_var_list:
                    # Check res for any variable specific options that need to be used BEFORE going to the plot:
                    if var in res:
                        vres = res[var]
                        #If found then notify user, assuming debug log is enabled:
                        adfobj.debug_log(f"time_series: Found variable defaults for {var}")
                        
                        #Extract category (if available):
                        web_category = vres.get("category", None)

                    else:
                        vres = {}
                    #End if

                    print(f"\t - Plotting Time Series, {season}")
                        
                    print("Plotting variable:",var)
                        
                    title_var = "Global"

                    """fig = plt.figure(figsize=(12,8))
                    ax = fig.add_subplot(111)

                    ax = ts_plot(ax, var, season, units, title_var)

                    # Create lists to hold all sets of years (for each case) and
                    # sets of var data (for each case)
                    vals_cases = []
                    yrs_cases = []"""

                    for case_idx, case_name in enumerate(case_names):
                        plot_name = plot_location[case_idx] / f"{var}_{season}_TimeSeries_Mean.{plot_type}"

                        fig = plt.figure(figsize=(12,8))
                        ax = fig.add_subplot(111)

                        ax = ts_plot(ax, var, season, units, title_var)

                        # Create lists to hold all sets of years (for each case) and
                        # sets of var data (for each case)
                        vals_cases = []
                        yrs_cases = []
                        
                        #if case_idx == len(case_base_names)-1:
                        if custom_leg == True:
                            label=f"{labels[case_idx]} (baseline)"
                        else:
                            label=f"{case_name} (baseline)"
                        marker = "--"

                        ax.plot(yrs[case_name].astype(int), vals[var][data_name][season], marker, c='g',
                                            label=label)
                        #else:
                        if custom_leg == True:
                            label=f"{labels[case_idx]}"
                        else:
                            label=f"{case_name}"
                        marker = "-"
                        ax.plot(yrs[case_name].astype(int), vals[var][case_name][season], marker, c=colors[case_idx],
                                            label=label)

                        vals_cases.append(vals[var][case_name][season])
                        yrs_cases.append(yrs[case_name])

                        ax = plot_var_details(ax, var, vals_cases, units[var], title_var, **vres)

                        ax = _format_xaxis(ax, yrs_cases)

                        # Set up legend
                        # If custom_legend = True, change the code in make_fig_legend() function for custom legend
                        #fig = make_fig_legend(case_names_len, fig, custom_legend=False)
                        fig.legend(loc="center left",fontsize=12,
                                            bbox_to_anchor=(0.122, 0.82,.05,.05))
                        plt.savefig(plot_name, facecolor='w')
                        #Add plot to website (if enabled):
                        adfobj.add_website_data(plot_name, var, case_name=case_name, category=web_category, season=season, plot_type="TimeSeries")
                        #Close plots:
                        plt.close()
                    #End for (case loop)
                #End for (variables loop)
            #End for (season loop)
        #End if (multi_path)

        #Derived quantities:
        #-------------------
        
        # - RESTOM
        #if ("FSNT" in ts_var_list) and ("FLNT" in ts_var_list):
        if all(value in ts_var_list for value in ["FSNT","FLNT"]):
            print("DID IT MAKE IT TO THIS POINT, AHHHHHHHHHHHH\n")
            #if any((match := item) in ts_var_list for item in ["FSNT","FLNT"]):
   
            if season == "ANN":
                var = "RESTOM"
                if var in res:
                    vres = res[var]
                    #If found then notify user, assuming debug log is enabled:
                    adfobj.debug_log(f"time_series: Found variable defaults for {var}")

                else:
                    vres = {}
                #End if

                if len(case_names) > 1:
                    plot_name = plot_loc / f"{var}_{season}_TimeSeries_multi_plot.{plot_type}"
                else:
                    plot_name = plot_loc / f"{var}_{season}_TimeSeries_Mean.{plot_type}"
                print(f"\t - Plotting Time Series, {season}")

                print("Plotting variable:",var)

                fig = plt.figure(figsize=(12,8))
                ax = fig.add_subplot(111)

                ax = ts_plot(ax, var, season, units, title_var)
                                
                # Create lists to hold all sets of years (for each case) and
                # sets of var data (for each case)
                vals_cases = []
                yrs_cases = []

                case_base_names = [i for i in case_names]
                case_base_names.append(data_name)

                for case_idx, case_name in enumerate(case_base_names):
                    if case_idx == len(case_base_names)-1:
                        if custom_leg == True:
                            label=f"{labels[case_idx]} (baseline)"
                        else:
                            label=f"{case_name} (baseline)"
                                    
                        if len(yrs[case_name]) < 5:
                            marker = "--*"
                        else:
                            marker = "--"
                        ax.plot(yrs[case_name].astype(int), vals[var][case_name][season], marker, c='g',
                                            label=label)
                    else:
                        if custom_leg == True:
                            label=f"{labels[case_idx]}"
                        else:
                            label=f"{case_name}"
                        if len(yrs[case_name]) < 5:
                            marker = "-*"
                        else:
                            marker = "-"

                        ax.plot(yrs[case_name].astype(int), vals[var][case_name][season], marker, c=colors[case_idx],
                                            label=label)

                    vals_cases.append(vals[var][case_name][season])
                    yrs_cases.append(yrs[case_name])

                ax = plot_var_details(ax, var, vals_cases, units[var], title_var, **vres)

                ax = _format_xaxis(ax, yrs_cases)

                # Set up legend
                # If custom_legend = True, change the code in make_fig_legend() function for custom legend
                #fig = make_fig_legend(case_names_len, fig, custom_legend=False)
                fig.legend(loc="center left",fontsize=12,
                                            bbox_to_anchor=(0.122, 0.82,.05,.05))
                plt.savefig(plot_name, facecolor='w')
                #Add plot to website (if enabled):
                adfobj.add_website_data(plot_name, var, case_name=None, category=web_category, season=season, plot_type="TimeSeries",multi_case=True)
                #Close plots:
                plt.close()
    # End for (seasons loop)

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

def _data_calcs(ts_loc,var,subset=None):
    """
    args
    ----
    - ts_loc: Path
        path to time series file
            
    - var: str
        name of variable
            
    - subset (optional): dict 
        lat/lon extents (south, north, east, west)
    """
    #print("ts_loc",ts_loc,"\n")
    fils = sorted(list(Path(ts_loc).glob(f"*{var}*.nc")))

    ts_ds = _load_dataset(fils)

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

    # global weighting
    w = np.cos(np.radians(data.lat))
    avg = data.weighted(w).mean(dim=("lat","lon"))

    yrs = np.unique([str(val.item().timetuple().tm_year).zfill(4) for _,val in enumerate(ts_ds["time"])])

    return avg,month_length,yrs,unit

########

def seasonal_data(data, month_length):
    """
    Function to grab seasonal weighted data, grouped by season
    
    -> DJF, MAM, JJA, and SON
    """

    weighted_mean = ((data * month_length).resample(time='QS-DEC').sum() /
                          month_length.resample(time='QS-DEC').sum()) #.resample(time='QS-DEC')

    mdata_season_mean_all_years = weighted_mean.resample(time='QS-DEC').mean()
    mdata_seasonal_mean = mdata_season_mean_all_years.groupby('time.season')
    
    return mdata_seasonal_mean

########

def _get_seasonal_data(ts_var_list, case_ts_locs, subsetz):
    vals = {}
    yrs = {}
    units = {}
    restom_dict = {}
    ts_var_list2 = ts_var_list#+["RESTOM"]
    restom_dict["RESTOM"] = {}
    vals["RESTOM"] = {}
    #ts_var_list2.remove("RESTOM")
    for var in ts_var_list2:
        print(f"Getting {var} data for all cases...\n")
        vals[var] = {}
        if (var == "FSNT") or (var == "FLNT"):
            restom_dict[var] = {}

        for case_loc in case_ts_locs:
            case = Path(case_loc).parts[-2]

            if var in subsetz: #subsetz["ICEFRAC"] = [subset_dict,"Subset"]
                print(f"{var} ({case}) has subset call...\n")
                data,month_length,_,unit = _data_calcs(case_loc,var,subsetz[var][0])
            else:
                data,month_length,_,unit = _data_calcs(case_loc,var)
                if var == "FSNT":
                    restom_dict[var][case] = data
                    #Copy units for RESTOM as well:
                    units["RESTOM"] = unit
                if var == "FLNT":
                    restom_dict[var][case] = data

            mdata_seasonal_mean = seasonal_data(data, month_length)

            vals[var][case] = {}
            vals["RESTOM"][case] = {}

            # Do ANN separately
            #Grab the years from the first variable (per case) only, no need to repeat over all variables
            if var == ts_var_list[0]:
                yrs_s = np.unique([str(val.item().timetuple().tm_year).zfill(4) for _,val in enumerate(data.time[1:])])
                yrs[case] = yrs_s
            vals[var][case]["ANN"] = [data.sel(time=i).mean().values for i in yrs[case]]

            #Loop over seasons (not including ANN)
            for season, vardata in mdata_seasonal_mean:
                if season == "DJF":
                    #Grab the years from the first case only, no need to repeat over all variables
                    if var == ts_var_list[0]:
                        yrs_s = np.unique([str(val.item().timetuple().tm_year).zfill(4) for _,val in enumerate(vardata.time[1:])])
                        yrs[case] = yrs_s
                    vals[var][case][season] = [vardata.sel(time=i).mean().values for i in yrs[case]]

                else:
                    #Grab the years from the first case only, no need to repeat over all variables
                    if var == ts_var_list[0]:
                        yrs_s = np.unique([str(val.item().timetuple().tm_year).zfill(4) for _,val in enumerate(vardata.time)])
                        yrs[case] = yrs_s
                    vals[var][case][season] = [vardata.sel(time=i).mean().values for i in yrs[case]]
        units[var] = unit

    #Section for derived quantities??
    #--------------------------------

    # - RESTOM
    if "FSNT" and "FLNT" in ts_var_list2:
        for case_loc in case_ts_locs:
            case = Path(case_loc).parts[-2]
            if len(yrs[case]) < 5:
                restom_data = restom_dict["FSNT"][case] - restom_dict["FLNT"][case]
                vals["RESTOM"][case]["ANN"] = [restom_data.sel(time=i).mean().values for i in yrs[case]]
            else:
                FSNT_avg = restom_dict["FSNT"][case].rolling(time=60,center=True).mean()
                FLNT_avg = restom_dict["FLNT"][case].rolling(time=60,center=True).mean()
                restom_data = FSNT_avg - FLNT_avg
                vals["RESTOM"][case]["ANN"] = [restom_data.sel(time=i).mean().values for i in yrs[case]]

    return vals, yrs, units

########

def ts_plot(ax, var, season, units, title_var):
    #Set Main title for subplots:
    ax.set_title(f"Time Series {title_var}: {var} - {season}",loc="left")
    # set axes titles
    ax.set_ylabel(units[var],fontsize=20,labelpad=12)
    ax.set_xlabel("Years",fontsize=15,labelpad=20)

    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=5)

    return ax

########

def plot_var_details(ax, var, vals_cases, unit, title_var, **kwargs):

    mins = []
    maxs = []
    for i,_ in enumerate(vals_cases):
        mins.append(np.nanmin(vals_cases[i]))
        maxs.append(np.nanmax(vals_cases[i]))

    #Set up plot details, if applicable from the adf_variable_defaults.yaml file
    if 'ts' in kwargs:
        print("Checking if desired units are different than raw from file...\n")
        if "units" in kwargs['ts']:
            unit = kwargs['ts']["units"]

        if 'major_locator' in kwargs['ts']:
            major_locator = kwargs['ts']['major_locator']

        if 'minor_locator' in kwargs['ts']:
            minor_locator = kwargs['ts']['minor_locator']
    else:
        minor_locator = None
        major_locator = None

    if minor_locator:
        ax.yaxis.set_minor_locator(MultipleLocator(minor_locator))
    if major_locator:
        ax.yaxis.set_major_locator(MultipleLocator(major_locator))

    #print(f"for {var}: major_locator: {major_locator} and minor_locator: {minor_locator}\n")
    ax.set_ylabel(unit,fontsize=20,labelpad=12)

    if major_locator:
        ax.set_ylim(np.floor(min(mins))-major_locator,max(maxs)+major_locator)
    #ax.set_ylim(np.floor(min(mins)),np.ceil(max(maxs)))

    if var == "RESTOM":
        # Set label to show if RESTOM is 1 or 5-yr avg
        line_1yr = Line2D([], [], label='1-yr avg', color='k', linewidth=1,marker='*',)       
        line_5yr = Line2D([], [], label='5-yr avg', color='k', linewidth=1,)
        ax.legend(handles=[line_1yr,line_5yr], bbox_to_anchor=(0.965, 0.96,.042,.05),) # 0.083,0.885 fontsize=14) # bbox_to_anchor=(0.99, 0.9))
        if major_locator:
            ax.set_ylim(np.floor(min(mins))-major_locator,max(maxs)+(major_locator*1.5))

    # Add extra space on the y-axis, except for ICEFRAC
    if var == "ICEFRAC":
        if major_locator:
            ax.set_ylim(np.floor(min(mins)),np.ceil(max(maxs)))

    return ax

########

def _format_xaxis(ax, yrs):
    """
    Set the x-axis plot limits to guarantee data range from all cases (including baseline)
    """

   # Grab the earliest and latest climo years of all cases excluding baseline
    min_year = int(min([min(i) for i in yrs]))
    max_year = int(max([max(i) for i in yrs]))

    # Grab the earliest and latest climo years of baseline
    #max_year = max([yrs_cases_max,max(yrs_base_int)])
    #min_year = min(yrs_cases_min,min(yrs_base_int))

    # Set the x-axis plot limits
    # to guarantee data range from all cases (including baseline)
    # Also set x-range to round to nearest five integer
    #last_val = max([yrs_cases_max,max(yrs_base_int)])
    #last_val = 
    last_year = int(max_year) - int(max_year) % 5
    if (max_year > 5) and (last_year < max_year):
        last_year += 5
   
    #first_val = min(yrs_cases_min,min(yrs_base_int))
    first_year = int(min_year) - int(min_year) % 5
    if min_year < 5:
        min_year = 0

    #print(first_year, last_year)
    ax.set_xlim(first_year, last_year)

    # x-axis ticks and numbers
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


