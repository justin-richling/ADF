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
#Import standard modules:
import os
import subprocess
from pathlib import Path
import numpy as np
import xarray as xr
import warnings  # use to warn user about missing files.

# Import plotting modules:
import matplotlib as mpl
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import plotting_functions as pf

# Warnings
import warnings  # use to warn user about missing files.
#     - Format warning messages:
def my_formatwarning(msg, *args, **kwargs):
    """Issue `msg` as warning."""
    return str(msg) + '\n'
warnings.formatwarning = my_formatwarning

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
    pf.get_central_longitude()
        determine central longitude for global plots
    pf.lat_lon_validate_dims()
        makes sure latitude and longitude are valid
    pf.seasonal_mean()
        calculate seasonal mean
    pf.plot_map_and_save()
        send information to make the plot and save the file
    pf.zm_validate_dims()
        Checks on pressure level dimension
    """

    #Notify user that script has started:
    print("\n  Generating lat/lon maps...")

    #
    # Use ADF api to get all necessary information
    #
    var_list = adfobj.diag_var_list
    #Special ADF variable which contains the output paths for
    #all generated plots and tables for each case:
    plot_locations = adfobj.plot_location

    #Grab case years
    syear_cases = adfobj.climo_yrs["syears"]
    eyear_cases = adfobj.climo_yrs["eyears"]

    #Grab baseline years (which may be empty strings if using Obs):
    syear_baseline = adfobj.climo_yrs["syear_baseline"]
    eyear_baseline = adfobj.climo_yrs["eyear_baseline"]

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

    weight_season = True  #always do seasonal weighting

    #Set seasonal ranges:
    seasons = {"ANN": np.arange(1,13,1),
               "DJF": [12, 1, 2],
               "JJA": [6, 7, 8],
               "MAM": [3, 4, 5],
               "SON": [9, 10, 11]
               }

    # probably want to do this one variable at a time:
    for var in var_list:
        if var not in adfobj.data.ref_var_nam:
            dmsg = f"No reference data found for variable `{var}`, global lat/lon mean plotting skipped."
            adfobj.debug_log(dmsg)
            print(dmsg)
            continue        

        #Notify user of variable being plotted:
        print("\t - lat/lon maps for {}".format(var))

        # Check res for any variable specific options that need to be used BEFORE going to the plot:
        if var in res:
            vres = res[var]
            #If found then notify user, assuming debug log is enabled:
            adfobj.debug_log(f"global_latlon_map: Found variable defaults for {var}")

            #Extract category (if available):
            web_category = vres.get("category", None)

        else:
            vres = {}
            web_category = None
        #End if

        # For global maps, also set the central longitude:
        # can be specified in adfobj basic info as 'central_longitude' or supplied as a number,
        # otherwise defaults to 180
        vres['central_longitude'] = pf.get_central_longitude(adfobj)

        # load reference data (observational or baseline)
        if not adfobj.compare_obs:
            base_name = adfobj.data.ref_case_label
        else:
            base_name = adfobj.data.ref_labels[var]

        # Gather reference variable data
        odata = adfobj.data.load_reference_regrid_da(base_name, var)

        if odata is None:
            dmsg = f"No regridded test file for {base_name} for variable `{var}`, global lat/lon mean plotting skipped."
            adfobj.debug_log(dmsg)
            continue

        o_has_dims = pf.validate_dims(odata, ["lat", "lon", "lev"]) # T iff dims are (lat,lon) -- can't plot unless we have both
        if (not o_has_dims['has_lat']) or (not o_has_dims['has_lon']):
            print(f"\t = skipping global map for {var} as REFERENCE does not have both lat and lon")
            continue

        #Loop over model cases:
        for case_idx, case_name in enumerate(adfobj.data.case_names):

            #Set case nickname:
            case_nickname = adfobj.data.test_nicknames[case_idx]

            #Set output plot location:
            plot_loc = Path(plot_locations[case_idx])

            #Check if plot output directory exists, and if not, then create it:
            if not plot_loc.is_dir():
                print("    {} not found, making new directory".format(plot_loc))
                plot_loc.mkdir(parents=True)

            #Load re-gridded model files:
            mdata = adfobj.data.load_regrid_da(case_name, var)

            #Skip this variable/case if the regridded climo file doesn't exist:
            if mdata is None:
                dmsg = f"No regridded test file for {case_name} for variable `{var}`, global lat/lon mean plotting skipped."
                adfobj.debug_log(dmsg)
                continue

            #Determine dimensions of variable:
            has_dims = pf.validate_dims(mdata, ["lat", "lon", "lev"])
            if (not has_dims['has_lat']) or (not has_dims['has_lon']):
                print(f"\t = skipping global map for {var} for case {case_name} as it does not have both lat and lon")
                continue
            else: # i.e., has lat&lon
                if (has_dims['has_lev']) and (not pres_levs):
                    print(f"\t - skipping global map for {var} as it has more than lev dimension, but no pressure levels were provided")
                    continue

            # Check output file. If file does not exist, proceed.
            # If file exists:
            #   if redo_plot is true: delete it now and make plot
            #   if redo_plot is false: add to website and move on
            doplot = {}

            if not has_dims['has_lev']:
                for s in seasons:
                    plot_name = plot_loc / f"{var}_{s}_LatLon_Mean.{plot_type}"
                    doplot[plot_name] = plot_file_op(adfobj, plot_name, var, case_name, s, web_category, redo_plot, "LatLon")
            else:
                for pres in pres_levs:
                    for s in seasons:
                        plot_name = plot_loc / f"{var}_{pres}hpa_{s}_LatLon_Mean.{plot_type}"
                        doplot[plot_name] = plot_file_op(adfobj, plot_name, f"{var}_{pres}hpa", case_name, s, web_category, redo_plot, "LatLon")
            if all(value is None for value in doplot.values()):
                print(f"All plots exist for {var}. Redo is {redo_plot}. Existing plots added to website data. Continue.")
                continue

            #Create new dictionaries:
            mseasons = {}
            oseasons = {}
            dseasons = {} # hold the differences

            if not has_dims['has_lev']:  # strictly 2-d data          

                #Loop over season dictionary:
                for s in seasons:
                    plot_name = plot_loc / f"{var}_{s}_LatLon_Mean.{plot_type}"
                    if doplot[plot_name] is None:
                        continue

                    if weight_season:
                        mseasons[s] = pf.seasonal_mean(mdata, season=s, is_climo=True)
                        oseasons[s] = pf.seasonal_mean(odata, season=s, is_climo=True)
                    else:
                        #Just average months as-is:
                        mseasons[s] = mdata.sel(time=seasons[s]).mean(dim='time')
                        oseasons[s] = odata.sel(time=seasons[s]).mean(dim='time')
                    #End if

                    # difference: each entry should be (lat, lon)
                    dseasons[s] = mseasons[s] - oseasons[s]

                    pf.plot_map_and_save(plot_name, case_nickname, adfobj.data.ref_nickname,
                                            [syear_cases[case_idx],eyear_cases[case_idx]],
                                            [syear_baseline,eyear_baseline],
                                            mseasons[s], oseasons[s], dseasons[s],
                                            obs=adfobj.compare_obs, **vres)

                    #Add plot to website (if enabled):
                    adfobj.add_website_data(plot_name, var, case_name, category=web_category,
                                            season=s, plot_type="LatLon")

            else: # => pres_levs has values, & we already checked that lev is in mdata (has_lev)

                for pres in pres_levs:

                    #Check that the user-requested pressure level
                    #exists in the model data, which should already
                    #have been interpolated to the standard reference
                    #pressure levels:
                    if (not (pres in mdata['lev'])) or (not (pres in odata['lev'])):
                        print(f"plot_press_levels value '{pres}' not present in {var} [test: {(pres in mdata['lev'])}, ref: {pres in odata['lev']}], so skipping.")
                        continue

                    #Loop over seasons:
                    for s in seasons:
                        plot_name = plot_loc / f"{var}_{pres}hpa_{s}_LatLon_Mean.{plot_type}"
                        if doplot[plot_name] is None:
                            continue

                        if weight_season:
                            mseasons[s] = pf.seasonal_mean(mdata, season=s, is_climo=True)
                            oseasons[s] = pf.seasonal_mean(odata, season=s, is_climo=True)
                        else:
                            #Just average months as-is:
                            mseasons[s] = mdata.sel(time=seasons[s]).mean(dim='time')
                            oseasons[s] = odata.sel(time=seasons[s]).mean(dim='time')
                        #End if

                        # difference: each entry should be (lat, lon)
                        dseasons[s] = mseasons[s] - oseasons[s]

                        pf.plot_map_and_save(plot_name, case_nickname, adfobj.data.ref_nickname,
                                                [syear_cases[case_idx],eyear_cases[case_idx]],
                                                [syear_baseline,eyear_baseline],
                                                mseasons[s].sel(lev=pres), oseasons[s].sel(lev=pres), dseasons[s].sel(lev=pres),
                                                obs=adfobj.compare_obs, **vres)

                        #Add plot to website (if enabled):
                        adfobj.add_website_data(plot_name, f"{var}_{pres}hpa", case_name, category=web_category,
                                                season=s, plot_type="LatLon")
                    #End for (seasons)
                #End for (pressure levels)
            #End if (plotting pressure levels)
        #End for (case loop)
    #End for (variable loop)

    # Check for AOD, and run the 4-panel diagnostics against MERRA and MODIS
    if "AODVISdn" in var_list:
        print("\tRunning AOD panel diagnostics...")
        aod_latlon(adfobj)

    #Notify user that script has ended:
    print("  ...lat/lon maps have been generated successfully.")


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
    int, None
        Returns 1 if existing file is removed or no existing file.
        Returns None if file exists and redo_plot is False

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
                                    season=season, plot_type=plot_type)
            return False  # False tells caller that file exists and not to overwrite
    else:
        return True
########


def aod_latlon(adfobj):
    var = "AODVISdn"
    season_abbr = ['Dec-Jan-Feb', 'Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov']
    # Define a list of season labels
    seasons = ['DJF', 'MAM', 'JJA', 'SON']

    test_case_names = adfobj.get_cam_info('cam_case_name', required=True)
    # load reference data (observational or baseline)
    if not adfobj.compare_obs:
        base_name = adfobj.data.ref_case_label
        case_names = test_case_names + [base_name]
    else:
        case_names = test_case_names

    #Grab all case nickname(s)
    test_nicknames = adfobj.case_nicknames["test_nicknames"]
    base_nickname = adfobj.case_nicknames["base_nickname"]
    case_nicknames = test_nicknames + [base_nickname]

    res = adfobj.variable_defaults # will be dict of variable-specific plot preferences
    # or an empty dictionary if use_defaults was not specified in YAML.
    res_aod_diags = res["aod_diags"]
    plot_params = res_aod_diags["plot_params"]
    plot_params_relerr = res_aod_diags["plot_params_relerr"]

    # Observational Datasets
    #-----------------------
    obs_dir = adfobj.get_basic_info("obs_data_loc")
    file_merra2 = os.path.join(obs_dir, 'MERRA2_192x288_AOD_2001-2020_climo.nc')
    file_mod08_m3 = os.path.join(obs_dir, 'MOD08_M3_192x288_AOD_2001-2020_climo.nc')

    ds_merra2 = xr.open_dataset(file_merra2)
    ds_merra2 = ds_merra2['TOTEXTTAU']
    ds_merra2['lon'] = ds_merra2['lon'].round(5)
    ds_merra2['lat'] = ds_merra2['lat'].round(5)

    ds_mod08_m3 = xr.open_dataset(file_mod08_m3)
    ds_mod08_m3 = ds_mod08_m3['AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean']
    ds_mod08_m3['lon'] = ds_mod08_m3['lon'].round(5)
    ds_mod08_m3['lat'] = ds_mod08_m3['lat'].round(5)

    ds_merra2_season = monthly_to_seasonal(ds_merra2)
    ds_merra2_season['lon'] = ds_merra2_season['lon'].round(5)
    ds_merra2_season['lat'] = ds_merra2_season['lat'].round(5)

    ds_mod08_m3_season = monthly_to_seasonal(ds_mod08_m3)
    ds_mod08_m3_season['lon'] = ds_mod08_m3_season['lon'].round(5)
    ds_mod08_m3_season['lat'] = ds_mod08_m3_season['lat'].round(5)

    ds_obs = [ds_mod08_m3_season, ds_merra2_season]
    obs_titles = ["TERRA MODIS", "MERRA2"]

    # Model Case Datasets
    #-----------------------
    ds_cases = []

    for case in test_case_names:
        """
        TODO: Need to check grid of test data in case they are on a different
                 grid than these particular observational data sets!
        
        """
        #Load re-gridded model files:
        ds_case = adfobj.data.load_climo_da(case, var)

        #Skip this variable/case if the regridded climo file doesn't exist:
        if ds_case is None:
            dmsg = f"No regridded test file for {case} for variable `{var}`, global lat/lon plots skipped."
            adfobj.debug_log(dmsg)
            continue
        else:
            ds_case['lon'] = ds_case['lon'].round(5)
            ds_case['lat'] = ds_case['lat'].round(5)

            # Calculate seasonal means
            ds_case_season = monthly_to_seasonal(ds_case)
            ds_case_season['lon'] = ds_case_season['lon'].round(5)
            ds_case_season['lat'] = ds_case_season['lat'].round(5)
            ds_cases.append(ds_case_season)

    # load reference data (observational or baseline)
    if not adfobj.compare_obs:
        base_name = adfobj.data.ref_case_label
    
        # Gather reference variable data
        ds_base = adfobj.data.load_reference_climo_da(base_name, var)
        if ds_base is None:
            dmsg = f"No regridded test file for {base_name} for variable `{var}`, global lat/lon plots skipped."
            adfobj.debug_log(dmsg)
        else:
            ds_base['lon'] = ds_base['lon'].round(5)
            ds_base['lat'] = ds_base['lat'].round(5)

            # Calculate seasonal means
            ds_base_season = monthly_to_seasonal(ds_base)
            ds_base_season['lon'] = ds_base_season['lon'].round(5)
            ds_base_season['lat'] = ds_base_season['lat'].round(5)
            ds_cases.append(ds_base_season)

    # Number of relevant cases        
    case_num = len(ds_cases)

    
    # 4-Panel global lat/lon plots
    #-----------------------------
    for i_obs,ds_ob in enumerate(ds_obs):
        for i_s,season in enumerate(seasons):
            plotnames = []
            fields = []
            params = []
            types = []
            case_namez = []

            obs_name = obs_titles[i_obs]
            chem_season = season_abbr[i_s]

            for i_case,ds_case in enumerate(ds_cases):
                case_nickname = case_nicknames[i_case]

                case_field = ds_case.sel(season=season) - ds_ob.sel(season=season)
                plotnames.append(f'{case_nickname} - {obs_name}\nAOD 550 nm - ' + chem_season)
                fields.append(case_field)
                params.append(plot_params)
                types.append("Diff")
                case_namez.append(case_names[i_case])

                field_relerr = 100 * case_field / ds_ob.sel(season=season)
                field_relerr = np.clip(field_relerr, -100, 100)
                plotnames.append(f'Percent Diff {case_nickname} - {obs_name}\nAOD 550 nm - ' + chem_season)
                fields.append(field_relerr)
                params.append(plot_params_relerr)
                types.append("Percent Diff")
                case_namez.append(case_names[i_case])
            # End for

            # Create 4-panel plot for season
            aod_panel_latlon(adfobj, plotnames, params, fields, season, obs_name, case_namez, case_num, types, symmetric=True)
        # End for
    # End for


#######################################
# Helper functions for AOD 4-panel plts
#######################################

def monthly_to_seasonal(ds,obs=False):
    ds_season = xr.Dataset(
        coords={'lat': ds.coords['lat'], 'lon': ds.coords['lon'],
                'season': np.arange(4)})
    da_season = xr.DataArray(
         coords=ds_season.coords, dims=['lat', 'lon', 'season'])
    
    # Create a list of DataArrays
    dataarrays = []
    # Define a list of season labels
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    
    if obs:
        for varname in ds:
            if '_n' not in varname:
                ds_season = xr.zeros_like(da_season)
                for s in seasons:
                    dataarrays.append(pf.seasonal_mean(ds, season=s, is_climo=True))
    else:
        for s in seasons:
            dataarrays.append(pf.seasonal_mean(ds, season=s, is_climo=True))

    # Use xr.concat to combine along a new 'season' dimension
    ds_season = xr.concat(dataarrays, dim='season')

    # Assign the 'season' labels to the new 'season' dimension
    ds_season['season'] = seasons
    ds_season = ds_season.transpose('lat', 'lon', 'season')

    return ds_season
#######


def aod_panel_latlon(adfobj, plotnames, plot_params, fields, season, obs_name, case_name, case_num, types, symmetric=False):

    #Set plot details:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    file_type = basic_info_dict.get('plot_type', 'png')
    plot_dir = adfobj.plot_location[0]

    # create figure:
    fig = plt.figure(figsize=(7*case_num,10))
    proj = ccrs.PlateCarree()

    # LAYOUT WITH GRIDSPEC
    plot_len = int(3*case_num)
    gs = mpl.gridspec.GridSpec(2*case_num, plot_len, wspace=0.5, hspace=0.0)
    gs.tight_layout(fig)

    axs = []
    for i in range(case_num):
        start = i * 3
        end = (i + 1) * 3
        axs.append(plt.subplot(gs[0:case_num, start:end], projection=proj))
        axs.append(plt.subplot(gs[case_num:, start:end], projection=proj))

    # formatting for tick labels
    lon_formatter = LongitudeFormatter(number_format='0.0f',
                                        degree_symbol='',
                                        dateline_direction_label=False)
    lat_formatter = LatitudeFormatter(number_format='0.0f',
                                        degree_symbol='')

    # 
    for i,field in enumerate(fields):
        # Set up sub plots for main panel plot
        ind_fig, ind_ax = plt.subplots(1, 1, figsize=((7*case_num)/2,10/2),subplot_kw={'projection': proj})

        lon_values = field.lon.values
        lat_values = field.lat.values

        # Define plot levels
        levels = np.linspace(
            plot_params[i]['range_min'], plot_params[i]['range_max'],
            plot_params[i]['nlevel'], endpoint=True)
        if 'augment_levels' in plot_params[i]:
            levels = sorted(np.append(
                levels, np.array(plot_params[i]['augment_levels'])))

        if field.ndim > 2:
            field_values = field.values[0,:,:]
            print(f"Required 2d lat/lon coordinates, got {field.ndim}d")
            emg = "AOD panel plot:\n"
            emg += f"\t Too many dimensions for {case_name}. Needs 2 (lat/lon) but got {field.ndim}"
            adfobj.debug_log(emg)
            return

        # Get data
        field_values = field.values[:,:]
        field_values, lon_values  = add_cyclic_point(field_values, coord=lon_values)
        lon_mesh, lat_mesh = np.meshgrid(lon_values, lat_values)
        field_mean = np.nanmean(field_values)

        extend_option = 'both' if symmetric else 'max'
        cmap_option = plt.cm.bwr if symmetric else plt.cm.turbo

        img = axs[i].contourf(lon_mesh, lat_mesh, field_values,
            levels, cmap=cmap_option, extend=extend_option,
                              transform_first=True,
            transform=ccrs.PlateCarree())
        ind_img = ind_ax.contourf(lon_mesh, lat_mesh, field_values,
            levels, cmap=cmap_option, extend=extend_option,
                              transform_first=True,
            transform=ccrs.PlateCarree())

        # ax.gridlines()
        axs[i].set_facecolor('gray')
        ind_ax.set_facecolor('gray')
        axs[i].coastlines()
        ind_ax.coastlines()

        # Averages plot titles
        axs[i].set_title(plotnames[i] + ('  Mean %.2g' % field_mean),fontsize=10)
        ind_ax.set_title(plotnames[i] + ('  Mean %.2g' % field_mean),fontsize=10)

        cbar = plt.colorbar(img, orientation='horizontal', pad=0.05)
        ind_cbar = plt.colorbar(ind_img, orientation='horizontal', pad=0.05)

        if 'ticks' in plot_params[i]:
            cbar.set_ticks(plot_params[i]['ticks'])
            ind_cbar.set_ticks(plot_params[i]['ticks'])
        if 'tick_labels' in plot_params[i]:
            cbar.ax.set_xticklabels(plot_params[i]['tick_labels'])
            ind_cbar.ax.set_xticklabels(plot_params[i]['tick_labels'])
        cbar.ax.tick_params(labelsize=6)

        # Save the individual figure
        pbase = f'AOD_{case_name[i]}_vs_{obs_name.replace(" ","_")}_{types[i].replace(" ","_")}'
        ind_plotfile = f'{pbase}_{season}_Chemistry_Mean.{file_type}'
        ind_png_file = Path(plot_dir) / ind_plotfile
        ind_fig.savefig(f'{ind_png_file}', bbox_inches='tight', dpi=300)
        plt.close(ind_fig)
    # End for

    # Save the panel figure
    plotfile = f'AOD_diff_panel_plot_{obs_name.replace(" ","_")}_{season}_LatLon_Mean'
    png_file = f'{plotfile}.{file_type}'
    png_file = Path(plot_dir) / png_file
    fig.savefig(png_file, bbox_inches='tight', dpi=300)
    adfobj.add_website_data(png_file, f'AOD_diff_panel_plot_{obs_name.replace(" ","_")}', None, season=season, multi_case=True, plot_type="LatLon", category="4-Panel AOD Diags")

    #pdf_file = f'{plotfile}.pdf'
    #ps_file = f'{plotfile}.ps'
    
    #fig.savefig(pdf_file, bbox_inches='tight')

    # Create and save PDF file
    #cmd = ["pdf2ps", pdf_file, ps_file]
    #subprocess.run(cmd, shell=False)
    plt.close(fig)


##############
#END OF SCRIPT