"""
Generate global maps of 2-D fields

Functions
---------
global_latlon_map(adfobj)
    use ADF object to make maps
my_formatwarning(msg, *args, **kwargs)
    format warning messages
    (private method)
plot_file_op
    Check on status of output plot file.
"""
# Import standard modules:
from pathlib import Path
import numpy as np
from collections import OrderedDict

# Import local modules:
import plotting_functions as pf
from aod_latlon import aod_latlon 
import adf_utils as utils
import plotting_utils as plot_utils

# Warnings
import warnings  # use to warn user about missing files.
warnings.formatwarning = utils.my_formatwarning
results = []
#########

def global_latlon_map(adfobj):
    """
    This script/function is designed to generate global
    2-D lat/lon maps of model fields with continental overlays.

    Parameters
    ----------
    adfobj : AdfDiag
        The diagnostics object that contains all the configuration information

    Returns
    -------
    Does not return a value; produces plots and saves files.

    Notes
    -----

    It uses the AdfDiag object's methods to get necessary information.
    Makes use of AdfDiag's data sub-class.
    Explicitly accesses:
    adfobj.diag_var_list
        List of variables
    adfobj.plot_location
        output plot path
    adfobj.climo_yrs
        start and end climo years of the case(s), `syears` & `eyears`
        start and end climo years of the reference, `syear_baseline` & `eyear_baseline`
    adfobj.variable_defaults 
        dict of variable-specific plot preferences
    adfobj.read_config_var
        dict of basic info, `diag_basic_info`
        Then use to check `plot_type`
    adfobj.debug_log
        Issues debug message
    adfobj.add_website_data
        Communicates information to the website generator
    adfobj.compare_obs
        Logical to determine if comparing to observations

        
    The `plotting_functions` module is needed for:
    plot_utils.get_central_longitude()
        determine central longitude for global plots
    utils.lat_lon_validate_dims()
        makes sure latitude and longitude are valid
    utils.seasonal_mean()
        calculate seasonal mean
    plot_utils.plot_map_and_save()
        send information to make the plot and save the file
    utils.zm_validate_dims()
        Checks on pressure level dimension
    """

    msg = "\n  Generating lat/lon maps..."
    print(f"{msg}\n  {'-' * (len(msg)-3)}")

    # Get configuration
    config = get_plot_config(adfobj)
    
    multi_plots = False
    multi_dict = {}

    #Grab case years
    syear_cases = adfobj.climo_yrs["syears"]
    eyear_cases = adfobj.climo_yrs["eyears"]

    case_names = adfobj.data.case_names
    """#print("adfobj.get_multi_case_info",adfobj.get_multi_case_info)
    #read_config_var('multi_case_plots')
    #multi_case_latlon = adfobj.multi_case_plots.get("global_latlon_map",[])
    multi_case_latlon = False
    if "global_latlon_map" in adfobj.multi_case_plots:
        if isinstance(adfobj.multi_case_plots,dict):
            multi_case_latlon = adfobj.multi_case_plots.get("global_latlon_map",[])
        else:
            multi_case_latlon = True
    #if multi_case:
    if len(case_names) > 1:
        #Check if multi-plots are desired from yaml file
        if multi_case_latlon: #adfobj.get_multi_case_info("global_latlon_map"):
            multi_plots = True
            multi_dict = OrderedDict()
            if multi_case_latlon != True:
                for multi_var in multi_case_latlon: #adfobj.get_multi_case_info("global_latlon_map"):
                    if multi_var not in multi_dict:
                        multi_dict[multi_var] = OrderedDict()
            #else:
            #    if multi_var not in multi_dict:
            #           multi_dict[multi_var] = OrderedDict()"""
    multi_plots = False
    if len(case_names) > 1:
        multi_plots = True
        multi_dict = OrderedDict()

    # Process regular variables
    config["multi_dict"] = multi_dict
    #print("INITIAL: multi_dict",multi_dict)

    if multi_plots:
        
        print("\n     Generating lat/lon multi-case plots...")
        for var in adfobj.diag_var_list:
            if var not in multi_dict:
                multi_dict[var] = OrderedDict()

            multi_dict, has_dims = process_variable(adfobj, var, **config)
            #This will be a list of variables for multi-case plotting based off LatLon plot type
            if multi_plots and multi_dict:
                #Notify user that script has started:
                

                test_nicknames = adfobj.data.test_nicknames
                base_nickname = adfobj.data.ref_nickname
                #test_climo_yrs = [adfobj.climo_yrs["syears"], adfobj.climo_yrs["eyears"]]
                test_climo_yrs = []
                for i in range(len(adfobj.climo_yrs["syears"])):
                    test_climo_yrs.append([adfobj.climo_yrs["syears"][i], adfobj.climo_yrs["eyears"][i]])
                base_climo_yrs = [adfobj.climo_yrs["syear_baseline"], adfobj.climo_yrs["eyear_baseline"]]
                res = adfobj.variable_defaults

                main_site_assets_path = adfobj.main_site_paths["main_site_assets_path"]
                #for var in multi_dict.keys():
                if 1==1:
                    print("VAR:",var)
                    vres = res.get(var, {})
                    if has_dims['has_lev']:
                        print(f"DOES {var} have dims????")
                        vres["levs"] = adfobj.get_basic_info("plot_press_levels")
                    web_category = vres.get("category", None)
                    
                    #pf.multi_latlon_plots(main_site_assets_path, var, "LatLon", case_names,
                    #                    [test_nicknames,base_nickname],
                    #                    test_climo_yrs, base_climo_yrs, multi_dict[var],
                    #                    web_category, adfobj, multi_case=True, **vres)
                    pf.multi_map_plots(
                                    main_site_assets_path,
                                    var,
                                    "LatLon",
                                    case_names,
                                    [test_nicknames,base_nickname],
                                    #test_climo_yrs,
                                    [syear_cases,eyear_cases],
                                    base_climo_yrs,
                                    multi_dict[var],
                                    web_category,
                                    adfobj,
                                    **vres
                                )

        print("     ...lat/lon multi-case plots have been generated successfully.")
    
    #print("FINAL: multi_dict",multi_dict)
    # Handle AOD special case
    if "AODVISdn" in adfobj.diag_var_list:
        print("\tRunning AOD panel diagnostics against MERRA and MODIS...")
        aod_latlon(adfobj)
        
    print("  ...lat/lon maps have been generated successfully.")

    #results.append(result)

    '''#This will be a list of variables for multi-case plotting based off LatLon plot type
    if multi_plots and multi_dict:
        #Notify user that script has started:
        print("\n     Generating lat/lon multi-case plots...")

        test_nicknames = adfobj.data.test_nicknames
        base_nickname = adfobj.data.ref_nickname
        #test_climo_yrs = [adfobj.climo_yrs["syears"], adfobj.climo_yrs["eyears"]]
        test_climo_yrs = []
        for i in range(len(adfobj.climo_yrs["syears"])):
            test_climo_yrs.append([adfobj.climo_yrs["syears"][i], adfobj.climo_yrs["eyears"][i]])
        base_climo_yrs = [adfobj.climo_yrs["syear_baseline"], adfobj.climo_yrs["eyear_baseline"]]
        res = adfobj.variable_defaults

        main_site_assets_path = adfobj.main_site_paths["main_site_assets_path"]
        for var in multi_dict.keys():
            print("VAR:",var)
            vres = res.get(var, {})
            web_category = vres.get("category", None)
            
            pf.multi_latlon_plots(main_site_assets_path, var, "LatLon", case_names,
                                [test_nicknames,base_nickname],
                                test_climo_yrs, base_climo_yrs, multi_dict[var],
                                web_category, adfobj, multi_case=True, **vres)

        print("     ...lat/lon multi-case plots have been generated successfully.")'''


def process_variable(adfobj, var, seasons, pres_levs, plot_type, redo_plot, multi_dict):
    vres = adfobj.variable_defaults.get(var, {})
    web_category = vres.get("category", None)

    # For global maps, also set the central longitude:
    # can be specified in adfobj basic info as 'central_longitude' or supplied as a number,
    # otherwise defaults to 180
    vres['central_longitude'] = plot_utils.get_central_longitude(adfobj)

    # Load reference data
    odata = load_reference_data(adfobj, var)
    if odata is None:
        print(f"[global_latlon_map][process_variable] finds no reference data.")
        return

    #Loop over model cases:
    for case_idx, case_name in enumerate(adfobj.data.case_names):
        multi_dict, has_dims = process_case(adfobj, case_name, case_idx, var, odata, 
                    seasons, pres_levs, plot_type, redo_plot,
                    vres, web_category, multi_dict)
    #print("\t process_variable: multi_dict",multi_dict)
    return multi_dict, has_dims

def load_reference_data(adfobj, var):
    """Load and validate reference data."""
    if not adfobj.compare_obs:
        base_name = adfobj.data.ref_case_label
    else:
        if var not in adfobj.data.ref_var_nam:
            dmsg = f"\t    WARNING: No obs data found for variable `{var}`, global lat/lon mean plotting skipped."
            adfobj.debug_log(dmsg)
            print(dmsg)
            return None
        base_name = adfobj.data.ref_labels[var]

    odata = adfobj.data.load_reference_regrid_da(base_name, var)
    if odata is None:
        print(f"\t    WARNING: No reference data found for {var}")
        return None

    o_has_dims = utils.validate_dims(odata, ["lat", "lon", "lev"])
    if (not o_has_dims['has_lat']) or (not o_has_dims['has_lon']):
        print(f"\t    WARNING: Reference data missing lat/lon dimensions")
        return None
        
    return odata


def process_case(adfobj, case_name, case_idx, var, odata, seasons, 
                pres_levs, plot_type, redo_plot, vres, web_category, multi_dict):
    """Process individual case data and generate plots."""
    plot_loc = Path(adfobj.plot_location[case_idx])
    plot_loc.mkdir(parents=True, exist_ok=True)

    mdata = adfobj.data.load_regrid_da(case_name, var)
    if mdata is None:
        return

    has_dims = utils.validate_dims(mdata, ["lat", "lon", "lev"])
    if (not has_dims['has_lat']) or (not has_dims['has_lon']):
        print(f"\t    WARNING: Model data missing lat/lon dimensions")
        return

    # Check pressure levels if 3D data
    if has_dims['has_lev'] and not pres_levs:
        print(f"\t    WARNING: 3D variable found but no pressure levels specified")
        return

    #multi_case = adfobj.multi_case_plots.get("global_latlon_map",[])
    if multi_dict:
        if var in multi_dict: #adfobj.get_multi_case_info("global_latlon_map"):
            multi_dict[var][case_name] = OrderedDict()
            for s in seasons.keys():
                if s not in multi_dict[var][case_name]:
                    multi_dict[var][case_name][s] = OrderedDict()
            #print("\n999999999999999999999999\n",multi_dict[var][case_name].keys(),"\n999999999999999999999999\n")
            multi_dict = process_plots(adfobj, mdata, odata, case_name, case_idx,
                        var, seasons, pres_levs, plot_loc, plot_type,
                        redo_plot, vres, web_category, has_dims, multi_dict)
            #print("\n999999990000000000000999999999\n",multi_dict[var][case_name].keys(),"\n999999990000000000000999999999\n")
            return multi_dict, has_dims


def get_plot_config(adfobj):
    """Get plotting configuration from ADF object."""
    return {
        'seasons': {
            "ANN": np.arange(1,13,1),
            "DJF": [12, 1, 2],
            "JJA": [6, 7, 8],
            "MAM": [3, 4, 5],
            "SON": [9, 10, 11]
        },
        'plot_type': adfobj.read_config_var("diag_basic_info").get('plot_type', 'png'),
        'redo_plot': adfobj.get_basic_info('redo_plot'),
        'pres_levs': adfobj.get_basic_info("plot_press_levels")
    }


def process_seasonal_data(mdata, odata, season, weight_season=True):
    """Helper function to calculate seasonal means and differences."""
    if weight_season:
        mseason = utils.seasonal_mean(mdata, season=season, is_climo=True)
        oseason = utils.seasonal_mean(odata, season=season, is_climo=True)
    else:
        mseason = mdata.sel(time=season).mean(dim='time')
        oseason = odata.sel(time=season).mean(dim='time')
    
    # Calculate differences
    dseason = mseason - oseason
    
    # Calculate percent change
    pseason = (mseason - oseason) / np.abs(oseason) * 100.0
    pseason = pseason.where(np.isfinite(pseason), np.nan)
    
    return mseason, oseason, dseason, pseason


def plot_file_op(adfobj, plot_name, var, case_name, season, web_category, redo_plot, plot_type):
    """Check if output plot needs to be made or remade.
    
    Parameters
    ----------
    adfobj : AdfDiag
        The diagnostics object that contains all the configuration information

    plot_name : Path
        path of the output plot

    var : str
        name of variable

    case_name : str
        case name
    
    season : str
        season being plotted

    web_category : str
        the category for this variable

    redo_plot : bool
        whether to overwrite existing plot with this file name

    plot_type : str
        the file type for the output plot

    Returns
    -------
    bool
        Returns True if existing file is removed or no existing file, i.e. make the plot.
        Returns False if file exists and redo_plot is False

    Notes
    -----
    The long list of parameters is because add_website_data is called
    when the file exists and will not be overwritten.
    
    """
    # Check redo_plot. If set to True: remove old plot, if it already exists:
    if plot_name.is_file():
        if redo_plot:
            plot_name.unlink()
            return True
        else:
            #Add already-existing plot to website (if enabled):
            adfobj.add_website_data(plot_name, var, case_name, category=web_category,
                                    season=season, plot_type=plot_type,script=__file__)
            return False  # False tells caller that file exists and not to overwrite
    else:
        return True


def process_plots(adfobj, mdata, odata, case_name, case_idx, var, seasons, 
                 pres_levs, plot_loc, plot_type, redo_plot, vres, web_category, has_dims, multi_dict):
    """Process and generate plots for different seasons and pressure levels.
    
    Parameters
    ----------
    adfobj : AdfDiag
        The diagnostics object containing configuration
    mdata : xarray.DataArray  
        Model data
    odata : xarray.DataArray
        Reference/observation data
    case_name : str
        Name of current case
    case_idx : int
        Index of current case
    var : str
        Variable name
    seasons : dict
        Dictionary of season definitions
    pres_levs : list
        Pressure levels to plot
    plot_loc : Path
        Output plot directory
    plot_type : str
        Plot file type (e.g. 'png')
    redo_plot : bool
        Whether to regenerate existing plots
    vres : dict
        Variable-specific plot settings
    web_category : str
        Category for website organization
    has_dims : dict
        Dictionary indicating which dimensions exist in data
        
    Returns
    -------
    None
    """
    def get_key_paths(d, parent=""):
        paths = []
        for key, value in d.items():
            full_key = f"{parent}.{key}" if parent else key
            paths.append(full_key)
            if isinstance(value, OrderedDict):
                paths.extend(get_key_paths(value, full_key))
        return paths

    #print("multi_dict at start of process_plots",get_key_paths(multi_dict),"\n-----------------------------------\n")
    # Get case nickname and years
    case_nickname = adfobj.data.test_nicknames[case_idx]
    syear_cases = adfobj.climo_yrs["syears"]
    eyear_cases = adfobj.climo_yrs["eyears"]
    syear_baseline = adfobj.climo_yrs["syear_baseline"]
    eyear_baseline = adfobj.climo_yrs["eyear_baseline"]
    
    # Check if files exist and build doplot dict
    doplot = check_existing_plots(adfobj, var, plot_loc, plot_type, 
                                case_name, seasons, pres_levs, 
                                has_dims, web_category, redo_plot)
    
    if not any(value for value in doplot.values()):
        print(f"\t    INFO: All plots exist for {var}. Redo is {redo_plot}. Existing plots added to website data.")
        return

    # Initialize seasonal data dictionaries
    mseasons = {}
    oseasons = {}
    dseasons = {} 
    pseasons = {}

    if not has_dims['has_lev']:
        # Process 2D data
        multi_dict = process_2d_plots(adfobj, mdata, odata, case_name, case_nickname,
                        var, seasons, plot_loc, plot_type, doplot,
                        mseasons, oseasons, dseasons, pseasons,
                        syear_cases[case_idx], eyear_cases[case_idx],
                        syear_baseline, eyear_baseline,
                        web_category, vres, multi_dict)
    else:
        # Process 3D data with pressure levels
        multi_dict = process_3d_plots(adfobj, mdata, odata, case_name, case_nickname, 
                        var, seasons, pres_levs, plot_loc, plot_type, doplot,
                        mseasons, oseasons, dseasons, pseasons,
                        syear_cases[case_idx], eyear_cases[case_idx],
                        syear_baseline, eyear_baseline,
                        web_category, vres, multi_dict)
    #print("\t process_plots: multi_dict",multi_dict)
    #print("multi_dict at END of process_plots",get_key_paths(multi_dict),"\n-----------------------------------\n")
    #print("\n999999---------999999999\n",multi_dict[var][case_name]["ANN"].keys(),"\n999999---------999999999\n")
    return multi_dict

def check_existing_plots(adfobj, var, plot_loc, plot_type, case_name, 
                        seasons, pres_levs, has_dims, web_category, redo_plot):
    """Check which plots need to be generated."""
    doplot = {}
    
    if not has_dims['has_lev']:
        for s in seasons:
            plot_name = plot_loc / f"{var}_{s}_LatLon_Mean.{plot_type}"
            doplot[plot_name] = plot_file_op(adfobj, plot_name, var, 
                                           case_name, s, web_category, 
                                           redo_plot, "LatLon")
    else:
        for pres in pres_levs:
            for s in seasons:
                plot_name = plot_loc / f"{var}_{pres}hpa_{s}_LatLon_Mean.{plot_type}"
                doplot[plot_name] = plot_file_op(adfobj, plot_name, 
                                               f"{var}_{pres}hpa",
                                               case_name, s, web_category, 
                                               redo_plot, "LatLon")
    return doplot


def process_2d_plots(adfobj, mdata, odata, case_name, case_nickname,
                    var, seasons, plot_loc, plot_type, doplot,
                    mseasons, oseasons, dseasons, pseasons,
                    syear_case, eyear_case, syear_baseline, eyear_baseline,
                    web_category, vres, multi_dict=None):
    """Process and generate 2D plots."""
    for s in seasons.keys():
        plot_name = plot_loc / f"{var}_{s}_LatLon_Mean.{plot_type}"
        if doplot[plot_name] is None:
            continue
            
        # Calculate seasonal means and differences
        mseasons[s], oseasons[s], dseasons[s], pseasons[s] = \
            process_seasonal_data(mdata, odata, s)

        if multi_dict:
            #if 1==1:
            if "LatLon" not in multi_dict[var][case_name][s]:
                multi_dict[var][case_name][s]["LatLon"] = {}

            multi_dict[var][case_name][s]["LatLon"]["m_data"] = {}
            multi_dict[var][case_name][s]["LatLon"]["o_data"] = {}
            multi_dict[var][case_name][s]["LatLon"]["diff_data"] = {}

        # Generate plot
        result = pf.plot_map_and_save(plot_name, case_nickname, adfobj.data.ref_nickname,
                            [syear_case, eyear_case],
                            [syear_baseline, eyear_baseline],
                            mseasons[s], oseasons[s], dseasons[s], pseasons[s],
                            obs=adfobj.compare_obs, multi_plot_dict=multi_dict, season=s, **vres)

        # Add to website
        adfobj.add_website_data(plot_name, var, case_name, 
                               category=web_category,
                               season=s, plot_type="LatLon",script=__file__)

        check_str = f'{case_name} - test'
        multi_dict[var][case_name][s]["LatLon"]["m_data"] = result[check_str]

        check_str = f'{case_name} - base'
        multi_dict[var][case_name][s]["LatLon"]["o_data"] = result[check_str]

        check_str = f'{case_name} - diff'
        multi_dict[var][case_name][s]["LatLon"]["diff_data"] = result[check_str]

    return multi_dict

def process_3d_plots(adfobj, mdata, odata, case_name, case_nickname,
                    var, seasons, pres_levs, plot_loc, plot_type, doplot,
                    mseasons, oseasons, dseasons, pseasons, 
                    syear_case, eyear_case, syear_baseline, eyear_baseline,
                    web_category, vres, multi_dict=None):
    """Process and generate 3D plots with pressure levels."""
    for pres in pres_levs:
        # Validate pressure level exists
        if (not (pres in mdata['lev'])) or (not (pres in odata['lev'])):
            print(f"\t    WARNING: plot_press_levels value '{pres}' not present " 
                  f"in {var} [test: {(pres in mdata['lev'])}, "
                  f"ref: {pres in odata['lev']}], so skipping.")
            continue

        for s in seasons.keys():
            plot_name = plot_loc / f"{var}_{pres}hpa_{s}_LatLon_Mean.{plot_type}"
            if doplot[plot_name] is None:
                continue

            # Calculate seasonal means and differences
            mseasons[s], oseasons[s], dseasons[s], pseasons[s] = \
                process_seasonal_data(mdata, odata, s)
            
            if multi_dict:
                if "LatLon" not in multi_dict[var][case_name][s]:
                    multi_dict[var][case_name][s]["LatLon"] = {}
                #print("process 3d plots pres:",pres,"\n\n")
                if pres not in multi_dict[var][case_name][s]["LatLon"]:
                    multi_dict[var][case_name][s]["LatLon"][pres] = {}

                multi_dict[var][case_name][s]["LatLon"][pres]["m_data"] = {}
                multi_dict[var][case_name][s]["LatLon"][pres]["o_data"] = {}
                multi_dict[var][case_name][s]["LatLon"][pres]["diff_data"] = {}
            #print("\nZZZZZZZZZZZ---------nZZZZZZZZZZZ\n",multi_dict[var][case_name][s].keys(),"\nZZZZZZZZZZZ---------nZZZZZZZZZZZ\n")
            # Generate plot
            vres["pres"] = pres
            vres["s"] = s
            #vres["case_name"] = case_name
            result = pf.plot_map_and_save(plot_name, case_nickname, adfobj.data.ref_nickname,
                                [syear_case, eyear_case],
                                [syear_baseline, eyear_baseline],
                                mseasons[s].sel(lev=pres), 
                                oseasons[s].sel(lev=pres),
                                dseasons[s].sel(lev=pres),
                                pseasons[s].sel(lev=pres),
                                obs=adfobj.compare_obs, multi_plot_dict=multi_dict, season=s,**vres)

            # Add to website
            adfobj.add_website_data(plot_name, f"{var}_{pres}hpa",
                                   case_name, category=web_category,
                                   season=s, plot_type="LatLon",script=__file__)

            check_str = f'{case_name} - test'
            multi_dict[var][case_name][s]["LatLon"][pres]["m_data"] = result[check_str]

            check_str = f'{case_name} - base'
            multi_dict[var][case_name][s]["LatLon"][pres]["o_data"] = result[check_str]

            check_str = f'{case_name} - diff'
            multi_dict[var][case_name][s]["LatLon"][pres]["diff_data"] = result[check_str]
    #print("multi_dict_var[case_names[r]][season][ptype].keys()",multi_dict[var][case_name][s]["LatLon"].keys())
    return multi_dict