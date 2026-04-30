"""Module to make polar stereographic maps."""
from pathlib import Path
from collections import OrderedDict
import numpy as np

# ADF library
import plotting_functions as pf
import adf_utils as utils

def get_hemisphere(hemi_type):
    """Helper function to convert plot type to hemisphere code.
    
    Parameters
    ----------
    hemi_type : str
        if `NHPolar` set NH, otherwise SH
        
    Returns
    -------
    str
        NH or SH
    """
    return "NH" if hemi_type == "NHPolar" else "SH"

def process_seasonal_data(mdata, odata, season):
    """Helper function to calculate seasonal means and differences.
    Parameters
    ----------
    mdata : xarray.DataArray
        test case data
    odata : xarray.DataArray
        reference case data
    season : str
        season (JJA, DJF, MAM, SON)

    Returns
    -------
    mseason : xarray.DataArray
    oseason : xarray.DataArray
    dseason : xarray.DataArray
    pseason : xarray.DataArray
        Seasonal means for test, reference, difference, and percent difference    
    """
    mseason = utils.seasonal_mean(mdata, season=season, is_climo=True)
    oseason = utils.seasonal_mean(odata, season=season, is_climo=True)
    
    # Calculate differences
    dseason = mseason - oseason
    dseason.attrs['units'] = mseason.attrs['units']
    
    # Calculate percent change
    pseason = (mseason - oseason) / np.abs(oseason) * 100.0
    pseason.attrs['units'] = '%'
    pseason = pseason.where(np.isfinite(pseason), np.nan)
    pseason = pseason.fillna(0.0)
    
    return mseason, oseason, dseason, pseason

def polar_map(adfobj):
    """Generate polar maps of model fields with continental overlays."""
    #Notify user that script has started:
    msg = "\n  Generating polar maps..."
    print(f"{msg}\n  {'-' * (len(msg)-3)}")

    var_list = adfobj.diag_var_list

    #Special ADF variable which contains the output paths for
    #all generated plots and tables for each case:
    plot_locations = adfobj.plot_location

    #CAM simulation variables (this is always assumed to be a list):
    case_names = adfobj.get_cam_info("cam_case_name", required=True)

    """if "polar_map" in adfobj.multi_case_plots:
        if isinstance(adfobj.multi_case_plots,dict):
            multi_case_latlon = adfobj.multi_case_plots.get("polar_map",[])
        else:
            multi_case_latlon = True

    multi_plots = False
    if len(case_names) > 1:
        #Check if multi-plots are desired from yaml file
        if multi_case_latlon:
            multi_plots = True
            multi_dict = OrderedDict()
            #for multi_var in adfobj.get_multi_case_info("polar_map"):
            #    if multi_var not in multi_dict:
            #        multi_dict[multi_var] = OrderedDict()
    #End if (check for multiple cases)"""

    """if multi_plots:
        #if not adfobj.get_multi_case_info("global_latlon_map"):
        #        multi_dict[var] = OrderedDict()
        if adfobj.get_multi_case_info("polar_map"):
            for multi_var in adfobj.get_multi_case_info("polar_map"):
                if multi_var not in multi_dict:
                    multi_dict[multi_var] = OrderedDict()"""

    multi_plots = False
    if len(case_names) > 1:
        if adfobj.get_multi_case_info and "polar_map" in adfobj.get_multi_case_info:
            print("did it come here?")
            multi_plots = True
            multi_dict = OrderedDict()

    #Grab case years
    syear_cases = adfobj.climo_yrs["syears"]
    eyear_cases = adfobj.climo_yrs["eyears"]

    # if doing comparison to obs, but no observations are found, quit
    if adfobj.get_basic_info("compare_obs"):
        var_obs_dict = adfobj.var_obs_dict
        if not var_obs_dict:
            print("\t No observations found to plot against, so no polar maps will be generated.")
            return

    #Grab baseline years (which may be empty strings if using Obs):
    syear_baseline = adfobj.climo_yrs["syear_baseline"]
    eyear_baseline = adfobj.climo_yrs["eyear_baseline"]

    #Grab all case nickname(s)
    test_nicknames = adfobj.case_nicknames["test_nicknames"]
    base_nickname = adfobj.case_nicknames["base_nickname"]

    res = adfobj.variable_defaults # will be dict of variable-specific plot preferences
    # or an empty dictionary if use_defaults was not specified in YAML.

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')
    print(f"\t NOTE: Plot type is set to {plot_type}")

    # check if existing plots need to be redone
    redo_plot = adfobj.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")
    #-----------------------------------------

    #Determine if user wants to plot 3-D variables on
    #pressure levels:
    pres_levs = adfobj.get_basic_info("plot_press_levels")

    #Set seasonal ranges:
    seasons = {"ANN": np.arange(1,13,1),
               "DJF": [12, 1, 2],
               "JJA": [6, 7, 8],
               "MAM": [3, 4, 5],
               "SON": [9, 10, 11]
               }
    results = None
    has_levs = {}
    # probably want to do this one variable at a time:
    for var in var_list:
        print(f"\t - polar maps for {var}")

        if multi_plots:
            #for multi_var in adfobj.get_multi_case_info("polar_map"):
            if var not in multi_dict:
                multi_dict[var] = OrderedDict()

        if var not in adfobj.data.ref_var_nam:
            dmsg = f"\t    WARNING: No reference data found for variable `{var}`, polar lat/lon mean plotting skipped."
            adfobj.debug_log(dmsg)
            print(dmsg)
            continue

        if not adfobj.compare_obs:
            base_name = adfobj.data.ref_labels[var]
        else:
            base_name = adfobj.data.ref_case_label


        # Get variable-specific settings
        vres = res.get(var, {})
        web_category = vres.get("category", None)

        if multi_plots:
            vres["multi_plots"] = True

        # Get all plot info and check existence
        plot_info = []
        all_plots_exist = True
        
        for case_idx, case_name in enumerate(case_names):

            #Grab data for desired multi-plots (from yaml file)
            if multi_plots:
                multi_dict[var][case_name] = OrderedDict()
                """if var in adfobj.get_multi_case_info("polar_map"):
                    multi_dict[var][case_name] = OrderedDict()
                else:
                    continue"""

            plot_loc = Path(plot_locations[case_idx])

            tmp_ds = adfobj.data.load_regrid_dataset(case_name, var)
            if tmp_ds is None:
                continue

            has_lev = "lev" in tmp_ds.dims

            if has_lev:
                has_levs[var] = adfobj.get_basic_info("plot_press_levels")
            else:
                has_levs[var] = None

            for s in seasons:
                #if var in adfobj.get_multi_case_info("polar_map"):
                #    multi_dict[var][case_name][s] = {}
                multi_dict[var][case_name][s] = {}
                for hemi_type in ["NHPolar", "SHPolar"]:
                    #if var in adfobj.get_multi_case_info("polar_map"):
                    #    multi_dict[var][case_name][s][hemi_type] = {}
                    multi_dict[var][case_name][s][hemi_type] = {}
                    if pres_levs and has_lev: # 3-D variable & pressure levels specified
                        for pres in pres_levs:
                            plot_name = plot_loc / f"{var}_{pres}hpa_{s}_{hemi_type}_Mean.{plot_type}"
                            info = {
                                'path': plot_name,
                                'var': f"{var}_{pres}hpa",
                                'case': case_name,
                                'case_idx': case_idx,
                                'season': s,
                                'type': hemi_type,
                                'pressure': pres,
                                'exists': plot_name.is_file()
                            }

                            plot_info.append(info)
                            if (redo_plot is False) and info['exists']:
                                adfobj.add_website_data(info['path'], info['var'],
                                                    info['case'], category=web_category,
                                                    season=s, plot_type=hemi_type,script=__file__)
                            else:
                                all_plots_exist = False
                    elif (not has_lev): # 2-D variable
                        plot_name = plot_loc / f"{var}_{s}_{hemi_type}_Mean.{plot_type}"
                        info = {
                            'path': plot_name,
                            'var': var,
                            'case': case_name,
                            'case_idx': case_idx,
                            'season': s,
                            'type': hemi_type,
                            'exists': plot_name.is_file()
                        }
                        plot_info.append(info)
                        if (redo_plot is False) and info['exists']:
                            adfobj.add_website_data(info['path'], info['var'],
                                                  info['case'], category=web_category,
                                                  season=s, plot_type=hemi_type,script=__file__)
                        else:
                            all_plots_exist = False

        if all_plots_exist:
            print(f"\t    Skipping {var} - all plots already exist")
            continue

        odata = adfobj.data.load_reference_regrid_da(base_name, var)
        if odata is None:
            print(f"\t    WARNING: No reference data found for {var}")
            continue

        # Process each case
        for plot in plot_info:
            if plot['exists'] and not redo_plot:
                print("plot['exists'] and not redo_plot - DTFYTGUYHIUJOK")
                continue

            case_name = plot['case']
            case_idx = plot['case_idx']
            plot_loc = Path(plot_locations[case_idx])

            # Ensure plot directory exists
            plot_loc.mkdir(parents=True, exist_ok=True)

            # Load and validate model data (units transformation included in load_regrid_da)
            mdata = adfobj.data.load_regrid_da(case_name, var)
            if mdata is None:
                continue

            # Process data based on dimensionality
            if "lev" in mdata.dims:
                has_lev = True
            else:
                has_lev = False

            if has_lev and pres_levs and plot.get('pressure'):
                if not all(dim in mdata.dims for dim in ['lat', 'lev']):
                    continue
                mdata = mdata.sel(lev=plot['pressure'])
                odata_level = odata.sel(lev=plot['pressure'])
            else:
                if not utils.lat_lon_validate_dims(mdata):
                    continue

            # Calculate seasonal means and differences
            use_odata = odata_level if has_lev else odata
            mseason, oseason, dseason, pseason = process_seasonal_data(
                mdata, 
                use_odata,
                plot['season']
            )

            # Create plot
            if plot['path'].exists():
                plot['path'].unlink()

            result = pf.make_polar_plot(
                plot['path'], test_nicknames[case_idx], base_nickname,
                [syear_cases[case_idx], eyear_cases[case_idx]],
                [syear_baseline, eyear_baseline],
                mseason, oseason, dseason, pseason,
                hemisphere=get_hemisphere(plot['type']),
                obs=adfobj.compare_obs, **vres
            )


            test_str = f'{case_name} - test'
            base_str = f'{case_name} - base'
            diff_str = f'{case_name} - diff'
            press = plot.get('pressure', None)
            if press:
                if press not in multi_dict[var][case_name][plot['season']][plot['type']]:
                    multi_dict[var][case_name][plot['season']][plot['type']][plot['pressure']] = {}
                multi_dict[var][case_name][plot['season']][plot['type']][plot['pressure']]["m_data"] = result[test_str]
                multi_dict[var][case_name][plot['season']][plot['type']][plot['pressure']]["o_data"] = result[base_str]
                multi_dict[var][case_name][plot['season']][plot['type']][plot['pressure']]["diff_data"] = result[diff_str]
            else:
                multi_dict[var][case_name][plot['season']][plot['type']]["m_data"] = result[test_str]
                multi_dict[var][case_name][plot['season']][plot['type']]["o_data"] = result[base_str]
                multi_dict[var][case_name][plot['season']][plot['type']]["diff_data"] = result[diff_str]


            """test_str = f'{case_name} - test'
            base_str = f'{case_name} - base'
            diff_str = f'{case_name} - diff'

            if plot['pressure']:
                multi_dict[var][case_name][plot['season']][plot['type']][plot['pressure']]["m_data"] = result[test_str]
            else:
                multi_dict[var][case_name][plot['season']][plot['type']]["m_data"] = result[test_str]

            if plot['pressure']:
                multi_dict[var][case_name][plot['season']][plot['type']][plot['pressure']]["o_data"] = result[base_str]
            else:
                multi_dict[var][case_name][plot['season']][plot['type']]["o_data"] = result[base_str]

            if plot['pressure']:
                multi_dict[var][case_name][plot['season']][plot['type']][plot['pressure']]["diff_data"] = result[diff_str]
            else:
                multi_dict[var][case_name][plot['season']][plot['type']]["diff_data"] = result[diff_str]"""

            # Add to website
            adfobj.add_website_data(
                plot['path'], plot['var'], case_name,
                category=web_category, season=plot['season'],
                plot_type=plot['type'],script=__file__
            )

    #This will be a list of variables for multi-case plotting based off LatLon plot type
    if multi_plots:
        #Notify user that script has started:
        print("\n     Generating polar lat/lon multi-case plots...")

        hemis = ["NHPolar", "SHPolar"]
        main_site_assets_path = adfobj.main_site_paths["main_site_assets_path"]
        for var in multi_dict.keys():
            print("polar plotting sscript VAR:",var)
            
            vres = res.get(var, {})
            #print("vres",vres)
            web_category = vres.get("category", None)

            if has_levs[var]:
                print(f"DOES {var} have dims????")
                vres["levs"] = adfobj.get_basic_info("plot_press_levels")

            for hemi in hemis:
                #pf.multi_polar_plots(main_site_assets_path, var, hemi, case_names,
                #                    [test_nicknames,base_nickname], 
                #                    [syear_cases,eyear_cases], [syear_baseline,eyear_baseline], multi_dict[var],
                #                    web_category, adfobj, multi_case=True, **vres)
                pf.multi_map_plots(
                                main_site_assets_path,
                                var,
                                hemi,
                                case_names,
                                [test_nicknames,base_nickname],
                                [syear_cases,eyear_cases],
                                [syear_baseline,eyear_baseline],
                                multi_dict[var],
                                web_category,
                                adfobj,
                                **vres
                            )

        print("     ...polar lat/lon multi-case plots have been generated successfully.")
    
    """if not all(x is None for x in results):
        #Notify user that script has started:
        print("\n     Generating polar lat/lon multi-case plots...")

        hemis = ["NHPolar", "SHPolar"]
        main_site_assets_path = adfobj.main_site_paths["main_site_assets_path"]
        for var in multi_dict.keys():
            print("VAR:",var)
            print("multi_dict",multi_dict[var].keys(),"\n\n")
            
            vres = res.get(var, {})
            print("vres",vres)
            web_category = vres.get("category", None)
            for hemi in hemis:
                
                pf.multi_polar_plots(main_site_assets_path, var, hemi, case_names,
                                    [test_nicknames,base_nickname], 
                                    [syear_cases,eyear_cases], [syear_baseline,eyear_baseline], multi_dict[var],
                                    web_category, adfobj, multi_case=True, **vres)

        print("     ...polar lat/lon multi-case plots have been generated successfully.")"""

    print("  ...polar maps have been generated successfully.")

##############
#END OF `polar_map` function

##############
# END OF FILE
