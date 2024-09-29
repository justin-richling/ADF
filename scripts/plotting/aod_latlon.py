import os
import logging

import numpy as np
from scipy import stats

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import xarray as xr
import pandas as pd

import warnings  # use to warn user about missing files.

import plotting_functions as pf
from pathlib import Path

#Format warning messages:
def my_formatwarning(msg, *args, **kwargs):
    """Issue `msg` as warning."""
    return str(msg) + '\n'


def aod_latlon(adfobj):
    var = "AODVISdn"
    season_abbr = ['Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov', 'Dec-Jan-Feb']
    # Define a list of season labels
    seasons = ['MAM', 'JJA', 'SON','DJF']

    test_case_names = adfobj.get_cam_info('cam_case_name', required=True)
    case_names = test_case_names + [adfobj.get_baseline_info('cam_case_name')]

    base_name = adfobj.get_baseline_info('cam_case_name')

    #Grab all case nickname(s)
    test_nicknames = adfobj.case_nicknames["test_nicknames"]
    base_nickname = adfobj.case_nicknames["base_nickname"]
    case_nicknames = test_nicknames + [base_nickname]

    #Grab case years
    syears_case = adfobj.climo_yrs["syears"]
    eyears_case = adfobj.climo_yrs["eyears"]

    #Grab baseline years (which may be empty strings if using Obs):
    #syear_baseline = adfobj.climo_yrs["syear_baseline"]
    #eyear_baseline = adfobj.climo_yrs["eyear_baseline"]

    syears = syears_case + [adfobj.climo_yrs["syear_baseline"]]
    eyears = eyears_case + [adfobj.climo_yrs["eyear_baseline"]]

    plot_dir = adfobj.plot_location

    mam_dir = "/glade/derecho/scratch/richling/adf-output/tests/aod/climos/"


    # check if existing plots need to be redone
    redo_plot = adfobj.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")

    


    """# Check redo_plot. If set to True: remove old plots, if they already exist:
    if (not redo_plot) and plot_file.is_file():
        #Add already-existing plot to website (if enabled):
        adfobj.debug_log(f"'{plot_file}' and '{plot_file}' exist and clobber is false.")
        adfobj.add_website_data(plot_file, "QBO", None, season="TimeSeries", multi_case=True, non_season=True)
    """
    """plot_params = dict()
    plot_params['outdir'] = plot_dir
    plot_params['range_min'] = -0.4
    plot_params['range_max'] = 0.4
    plot_params['nlevel'] = 17

    plot_params_relerr = dict()
    plot_params_relerr['outdir'] = plot_dir
    plot_params_relerr['range_min'] = -100
    plot_params_relerr['range_max'] = 100
    plot_params_relerr['nlevel'] = 21
    """
    res = adfobj.variable_defaults # will be dict of variable-specific plot preferences
    # or an empty dictionary if use_defaults was not specified in YAML.
    res_aod_diags = res[var]["aod_diags"]
    plot_params = res_aod_diags["plot_params"]
    plot_params_relerr = res_aod_diags["plot_params_relerr"]

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')
    print(f"\t NOTE: Plot type is set to {plot_type}")

    # check if existing plots need to be redone
    redo_plot = adfobj.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")


    # Model Case Datasets
    #-----------------------
    ds_cases = []
    #mam_dir = "/glade/derecho/scratch/richling/adf-output/tests/aod/climos/"
    #mam4_1_dir = f"{mam_dir}/{case_names[0]}/yrs_1997-2000/"
    #mam4_2_dir = f"{mam_dir}/{case_names[1]}/yrs_1997-2000/"


    #o_has_dims = pf.validate_dims(ds_base, ["lat", "lon", "lev"]) # T iff dims are (lat,lon) -- can't plot unless we have both
    #if (not o_has_dims['has_lat']) or (not o_has_dims['has_lon']):
    #    print(f"\t = skipping global map for {var} as REFERENCE does not have both lat and lon")
    #   #continue

    #casedir = f"{mam_dir}/{case}/yrs_{syears[idx]}-{eyears[idx]}/"
    #case_path = os.path.join(casedir, f'{case}_{var}_climo.nc')
    #ds_case = xr.open_dataset(case_path)
    """ds_base['lon'] = ds_base['lon'].round(5)
    ds_base['lat'] = ds_base['lat'].round(5)

    ds_base_season = monthly_to_seasonal(ds_base)
    ds_base_season['lon'] = ds_base_season['lon'].round(5)
    ds_base_season['lat'] = ds_base_season['lat'].round(5)
    ds_base_season = ds_base_season[var]
    ds_cases.append(ds_base_season)"""
    #file_mam4_1 = os.path.join(mam4_1_dir, f'{case_names[0]}_AODVISdn_climo.nc')
    #file_mam4_2 = os.path.join(mam4_2_dir, f'{case_names[1]}_AODVISdn_climo.nc')

    for idx,case in enumerate(case_names):
        """
        TODO: Need to check grid of test data in case they are on a different
                 grid than these particular observational data sets!
        
        """
        #mam_dir = 
        #Load re-gridded model files:

        ds_case = adfobj.data.load_regrid_da(case, var)

        #Skip this variable/case if the regridded climo file doesn't exist:
        if ds_case is None:
            dmsg = f"No regridded test file for {case} for variable `{var}`, global lat/lon plots skipped."
            adfobj.debug_log(dmsg)
            continue
        #casedir = f"{mam_dir}/{case}/yrs_{syears[idx]}-{eyears[idx]}/"
        #case_path = os.path.join(casedir, f'{case}_{var}_climo.nc')
        #ds_case = xr.open_dataset(case_path)
        ds_case['lon'] = ds_case['lon'].round(5)
        ds_case['lat'] = ds_case['lat'].round(5)

        ds_case_season = monthly_to_seasonal(ds_case)
        ds_case_season['lon'] = ds_case_season['lon'].round(5)
        ds_case_season['lat'] = ds_case_season['lat'].round(5)
        #ds_case_season = ds_case_season[var]
        if idx == 0:
            print("\n",ds_case_season.shape,"\n")
        ds_cases.append(ds_case_season)
    
    # Gather reference variable data
    ds_base = adfobj.data.load_reference_regrid_da(base_name, var)
    if ds_base is None:
        dmsg = f"No regridded test file for {base_name} for variable `{var}`, global lat/lon plots skipped."
        adfobj.debug_log(dmsg)
        #continue

    ds_base['lon'] = ds_base['lon'].round(5)
    ds_base['lat'] = ds_base['lat'].round(5)

    ds_base_season = monthly_to_seasonal(ds_base)
    ds_base_season['lon'] = ds_base_season['lon'].round(5)
    ds_base_season['lat'] = ds_base_season['lat'].round(5)
    #ds_base_season = ds_base_season[var]
    ds_cases.append(ds_base_season)
    
    case_num = len(ds_cases)

    # Observational Datasets
    #-----------------------
    climo_dir = "/glade/work/fillmore/Data/AOD_climos/"
    file_merra2 = os.path.join(climo_dir, 'MERRA2_192x288_AOD_2001-2020_climo.nc')
    file_mod08_m3 = os.path.join(climo_dir, 'MOD08_M3_192x288_AOD_2001-2020_climo.nc')

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
    #ds_merra2_season.to_netcdf('MERRA2_192x288_AOD_2001-2020_seasonal_climo.nc')

    ds_mod08_m3_season = monthly_to_seasonal(ds_mod08_m3)
    ds_mod08_m3_season['lon'] = ds_mod08_m3_season['lon'].round(5)
    ds_mod08_m3_season['lat'] = ds_mod08_m3_season['lat'].round(5)
    #ds_mod08_m3_season.to_netcdf('MOD08_M3_192x288_AOD_2001-2020_seasonal_climo.nc')


    #ds_mod08_m3_season = ds_mod08_m3_season['AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean']
    #ds_merra2_season = ds_merra2_season['TOTEXTTAU']

    ds_obs = [ds_mod08_m3_season, ds_merra2_season]
    obs_titles = ["TERRA MODIS", "MERRA2"]


    # Individual global lat/lon plots
    #--------------------------------
    #for i_season in range(1):
    #    # Loop over all test cases - could be multi-case scenario
    #    for i_case,ds_case in enumerate(ds_cases):
    #
    #        case_name = str(case_names[i_case])
    #        varname = ds_case.name
    #        #print("case_name",case_name)
    #        plotfile = case_name + '_' + var + '_' + season_abbr[i_season]
    #        field = ds_case[:,:,i_season]
    #        plot_lon_lat(adfobj, plotfile, plot_dir, case_name, f'{case_name}\nAOD 550 nm  1997-2000' + ' ' + season_abbr[i_season], plot_params, field, i_season)
    """
        print("OBSIES")
        # Loop over supplied obs datasets
        for i_obs,ds_ob in enumerate(ds_obs):
            obs_name = obs_titles[i_obs]
            plotfile = obs_name + '_' + varname + '_' + season_abbr[i_season]
            field = ds_ob[:,:,i_season]
            plot_lon_lat(adfobj, plotfile, plot_dir, obs_name, f'{obs_name}\nAOD 550 nm  2001-2020' + ' ' + season_abbr[i_season], plot_params, field, i_season)
    """
    # 4-Panel global lat/lon plots
    #-----------------------------
    for i_obs,ds_ob in enumerate(ds_obs):
        for i_s,season in enumerate(seasons):
            plotnames = []
            fields = []
            params = []
            types = []
            case_namez = []
            season_abbr = ['Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov', 'Dec-Jan-Feb']
            print("ds_ob",ds_ob,"\n")

            obs_name = obs_titles[i_obs]
            chem_season = season_abbr[i_s]
            print(season)

            for i_case,ds_case in enumerate(ds_cases):
                case_nickname = case_nicknames[i_case]
                print(f"{case_nickname} minus {obs_name}")
                #case_field = ds_case[:,:,season]- ds_ob[:,:,season]
                case_field = ds_case.sel(season=season) - ds_ob.sel(season=season)
                plotnames.append(f'{case_nickname} - {obs_name}\nAOD 550 nm - ' + chem_season)
                fields.append(case_field)
                params.append(plot_params)
                types.append("Diff")
                case_namez.append(case_names[i_case])

                print(f"{case_nickname} minus {obs_name} % Diff")
                field_relerr = 100 * case_field / ds_ob.sel(season=season)
                #field_relerr = 100 * case_field / ds_ob[:,:,season]
                field_relerr = np.clip(field_relerr, -100, 100)
                plotnames.append(f'Percent Diff {case_nickname} - {obs_name}\nAOD 550 nm - ' + chem_season)
                fields.append(field_relerr)
                params.append(plot_params_relerr)
                types.append("Percent Diff")
                case_namez.append(case_names[i_case])

            yeah_boi(adfobj, plotnames, params, fields, season, obs_name, case_namez, case_num, types, symmetric=True)
            #yeah_boi(adfobj, plotnames, plot_params, fields, season, obs_name, case_name, case_num, symmetric=False)
























def monthly_to_seasonal(ds,obs=False):
    da = xr.DataArray(
        coords={'lat': ds.coords['lat'], 'lon': ds.coords['lon']},
        dims=['lat', 'lon'])
    ds_season = xr.Dataset(
        coords={'lat': ds.coords['lat'], 'lon': ds.coords['lon'],
                'season': np.arange(4)})
    da_season = xr.DataArray(
         coords=ds_season.coords, dims=['lat', 'lon', 'season'])
    dataarrays = []
    # Define a list of season labels
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    if obs:
        for varname in ds:
            if '_n' not in varname:
                #print(varname)
                # MAM, JJA, SON, DJF
                ds_season = xr.zeros_like(da_season)
                for i,s in enumerate(seasons):
                    #ds_season = pf.seasonal_mean(ds, season=s, is_climo=True)
                    dataarrays.append(pf.seasonal_mean(ds, season=s, is_climo=True))
    else:
        for i,s in enumerate(seasons):
            #ds_season = pf.seasonal_mean(ds, season=s, is_climo=True)
            dataarrays.append(pf.seasonal_mean(ds, season=s, is_climo=True))

    # Create a list of DataArrays
    #dataarrays = [da_DJF, da_MAM, da_JJA, da_SON]


    # Use xr.concat to combine along a new 'season' dimension
    ds_season = xr.concat(dataarrays, dim='season')

    # Assign the 'season' labels to the new 'season' dimension
    ds_season['season'] = seasons
    ds_season = ds_season.transpose('lat', 'lon', 'season')

    # The new DataArray now has dimensions ('season', 'lat', 'lon')
    #print(ds_season)
    print(ds_season.dims)
    return ds_season



def plot_lon_lat(adfobj, plotfile, plot_dir, case_name, plotname, plot_params, field, season, symmetric=False):
    import numpy as np
    logging.info(plotfile)


    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')

    plot_dir = adfobj.plot_location[0]

    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    ax = plt.axes(projection=ccrs.PlateCarree())

    lon_values = field.lon.values
    lat_values = field.lat.values

    levels = np.linspace(
        plot_params['range_min'], plot_params['range_max'],
        plot_params['nlevel'], endpoint=True)
    if 'augment_levels' in plot_params:
        levels = sorted(np.append(
            levels, np.array(plot_params['augment_levels'])))

    # ****** QUESTION: Would this ever be the case???? ******
    if field.ndim == 3:
        # field_values = np.clip(field.values[0,:,:], levels[0], levels[-1])
        field_values = field.values[0,:,:]
    else:
        # field_values = np.clip(field.values[:,:], levels[0], levels[-1])
        field_values = field.values[:,:]
        # field_values = field

    #import numpy as np

    #diffs = np.diff(lon_values)
    #print("diffies",diffs)

    field_values, lon_values \
        = add_cyclic_point(field_values, coord=lon_values)

    #print(lat_values.shape)
    #print(lon_values.shape)
    #print(field_values.shape)

    lon_mesh, lat_mesh = np.meshgrid(lon_values, lat_values)

    # print(np.nanmin(field_values), np.nanmax(field_values))
    field_mean = np.nanmean(field_values)

    extend_option = 'both' if symmetric else 'max'
    cmap_option = plt.cm.bwr if symmetric else plt.cm.turbo
    #print("I'm guessing this is where it's going to brizzeak??")

    cp = ax.contourf(lon_mesh, lat_mesh, field_values,
        levels, cmap=cmap_option, extend=extend_option,
        transform=ccrs.PlateCarree())

    # ax.gridlines()
    ax.set_facecolor('gray')
    ax.coastlines()
    # ax.add_feature(cfeature.BORDERS)
    # ax.add_feature(states_provinces)

    plt.title(plotname + ('  Mean %.2g' % field_mean))

    cbar = plt.colorbar(cp, orientation='horizontal', pad=0.05)

    if 'ticks' in plot_params:
        cbar.set_ticks(plot_params['ticks'])
    if 'tick_labels' in plot_params:
        cbar.ax.set_xticklabels(plot_params['tick_labels'])
    cbar.ax.tick_params(labelsize=6)

    #png_file = os.path.join(plot_dir, plotfile) + plot_type
    #pdf_file = os.path.join(plot_dir, plotfile) + '.pdf'
    #ps_file = os.path.join(plot_dir, plotfile) + '.ps'
    #plt.savefig(png_file, bbox_inches='tight', dpi=300)
    #plt.savefig(pdf_file, bbox_inches='tight')

    #plotfile = obs_name + '_' + varname + '_' + season_abbr[i_season]
    plotfile = f"AOD_{case_name}_{season}_Chemistry_Mean.{plot_type}"
    png_file = Path(plot_dir) / plotfile
    print("png_file")
    #png_file = Path(plot_dir) / f'QBO_Amplitude_Special_Mean.{plot_type}'
    #adfobj.add_website_data(png_file, f"AOD_{case_name}", None, season=season, multi_case=True, plot_type="Chemistry")
    #print()
    # Write final figure to file
    plt.savefig(png_file, bbox_inches='tight', dpi=300)
    # plot_loc_amp = Path(plot_locations[0]) / f'QBO_Amplitude_Special_Mean.{plot_type}'
    # adfobj.add_website_data(plot_loc_amp, "QBO", None, season="Amplitude", multi_case=True, non_season=True)
    #command = 'pdf2ps ' + pdf_file + ' ' + ps_file
    #os.system(command)
    plt.clf()




#yeah_boi(adfobj, plotnames, plot_params, fields, season, obs_name, case_name, case_num, types, symmetric=False)
def yeah_boi(adfobj, plotnames, plot_params, fields, season, obs_name, case_name, case_num, types, symmetric=False):

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')

    plot_dir = adfobj.plot_location[0]
    #plotfile = f'aod_output2/cases_vs_{obs_name.replace(" ","_")}_{season}'
    plotfile = f'AOD_diff_{obs_name.replace(" ","_")}_{season}_Chemistry_Mean.{plot_type}'
    png_file = Path(plot_dir) / plotfile
    #png_file = plotfile + '.png'
    #pdf_file = plotfile + '.pdf'
    #ps_file = plotfile + '.ps'

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    file_type = basic_info_dict.get('plot_type', 'png')

    # create figure:
    #fig = plt.figure(figsize=(14,10))
    fig = plt.figure(figsize=(7*case_num,10))

    proj = ccrs.PlateCarree()

    # LAYOUT WITH GRIDSPEC
    plot_len = int(3*case_num)
    print(plot_len)
    gs = mpl.gridspec.GridSpec(4, plot_len, wspace=0.5, hspace=0.0)
    gs.tight_layout(fig)

    axs = []
    for i in range(case_num):
        start = i * 3
        end = (i + 1) * 3
        print(f"{start}:{end}")
        axs.append(plt.subplot(gs[0:2, start:end], projection=proj))
        axs.append(plt.subplot(gs[2:, start:end], projection=proj))

    # formatting for tick labels
    lon_formatter = LongitudeFormatter(number_format='0.0f',
                                        degree_symbol='',
                                        dateline_direction_label=False)
    lat_formatter = LatitudeFormatter(number_format='0.0f',
                                        degree_symbol='')

    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')


    for i,field in enumerate(fields):

        ind_fig, ind_ax = plt.subplots(1, 1, figsize=((7*case_num)/2,10/2),subplot_kw={'projection': proj})

        lon_values = field.lon.values
        lat_values = field.lat.values

        levels = np.linspace(
            plot_params[i]['range_min'], plot_params[i]['range_max'],
            plot_params[i]['nlevel'], endpoint=True)
        if 'augment_levels' in plot_params[i]:
            levels = sorted(np.append(
                levels, np.array(plot_params[i]['augment_levels'])))

        if field.ndim == 3:
            field_values = field.values[0,:,:]
        else:
            field_values = field.values[:,:]

        #field_values, lon_values  = add_cyclic_point(field_values, coord=lon_values)
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
        # ax.add_feature(cfeature.BORDERS)
        # ax.add_feature(states_provinces)

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
        ind_cbar

        """#imgs.append(img)# = [img1,img2,img3,img4]
        axs[i].tick_params('both', length=5, width=1.5, which='major')
        ind_ax
        axs[i].tick_params('both', length=5, width=1.5, which='minor')
        ind_ax
        axs[i].xaxis.set_major_formatter(lon_formatter)
        ind_ax
        axs[i].yaxis.set_major_formatter(lat_formatter)
        ind_ax"""

        # Save the individual figure
        #ind_plotfile = f'aod_output2/{case_name[i]}_vs_{obs_name.replace(" ","_")}_{season}_{types[i]}'
        ind_plotfile = f'AOD_{case_name[i]}_vs_{obs_name.replace(" ","_")}_{season}_Chemistry_{types[i]}.{file_type}'
        #adfobj.add_website_data(ind_plotfile, "AOD", None, season=season, multi_case=False, plot_type="Chemistry")
        #adfobj.add_website_data(png_file, f'AOD_diff_{obs_name.replace(" ","_")}', None, season=season, multi_case=True, plot_type="Chemistry")
        print(ind_plotfile,"\n")
        ind_fig.savefig(f'{ind_plotfile}.{plot_type}', bbox_inches='tight', dpi=300)
        plt.close(ind_fig)

        # plot_loc_amp = Path(plot_locations[0]) / f'QBO_Amplitude_Special_Mean.{plot_type}'
        # adfobj.add_website_data(plot_loc_amp, "QBO", None, season="Amplitude", multi_case=True, non_season=True)

    # Save the 4-panel figure
    fig.savefig(png_file, bbox_inches='tight', dpi=300)
    adfobj.add_website_data(png_file, f'AOD_diff_{obs_name.replace(" ","_")}', None, season=season, multi_case=True, plot_type="Chemistry")
    if season == "MAM":
        png_file2 = png_file.replace(season,"ANN")
        adfobj.add_website_data(png_file2, f'AOD_diff_{obs_name.replace(" ","_")}', None, season="ANN", multi_case=True, plot_type="Chemistry")
    #fig.savefig(pdf_file, bbox_inches='tight')
    #command = 'pdf2ps ' + pdf_file + ' ' + ps_file
    #os.system(command)
    plt.close(fig)