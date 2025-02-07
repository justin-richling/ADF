#Import standard modules:
from pathlib import Path
import numpy as np
import xarray as xr
import warnings  # use to warn user about missing files.
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib as mpl
import matplotlib.cm as cm
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import metpy.calc.thermo as thermo
from metpy.units import units
#import metpy.constants as mconst

import plotting_functions as pf

#Format warning messages:
def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = my_formatwarning

def tem(adf):
    """
    Plot the contents of the TEM dignostic ouput of 2-D latitude vs vertical pressure maps.
    
    Steps:
     - loop through TEM variables
     - calculate all-time fields (from individual months)
     - take difference, calculate statistics
     - make plots

    Notes:
     - If any of the TEM cases are missing, the ADF skips this plotting script and moves on.

    """

    # Notify user that script has started:
    print("\n  Generating TEM plots ...")

    #Special ADF variable which contains the output paths for
    #all generated plots and tables for each case:
    plot_location = Path(adf.plot_location[0])

    #Check if plot output directory exists, and if not, then create it:
    if not plot_location.is_dir():
        print(f"    {plot_location} not found, making new directory")
        plot_location.mkdir(parents=True)

    #CAM simulation variables (this is always assumed to be a list):
    case_names = adf.get_cam_info("cam_case_name", required=True)

    res = adf.variable_defaults # will be dict of variable-specific plot preferences

    #Check if comparing against observations
    if adf.compare_obs:
        obs = True
        base_name = "Obs"
    else:
        obs = False
        base_name = adf.get_baseline_info("cam_case_name", required=True)
    #End if

    #Extract test case years
    syear_cases = adf.climo_yrs["syears"]
    eyear_cases = adf.climo_yrs["eyears"]

    #Extract baseline years (which may be empty strings if using Obs):
    syear_baseline = adf.climo_yrs["syear_baseline"]
    eyear_baseline = adf.climo_yrs["eyear_baseline"]

    #Grab all case nickname(s)
    test_nicknames = adf.case_nicknames["test_nicknames"]
    base_nickname = adf.case_nicknames["base_nickname"]
    case_nicknames = test_nicknames + [base_nickname]
 
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
    
    #Initialize list of input TEM file locations
    tem_locs = []

    #Extract TEM file save locations
    tem_case_locs = adf.get_cam_info("cam_tem_loc",required=True)
    tem_base_loc = adf.get_baseline_info("cam_tem_loc")
    output_loc       = adf.get_basic_info("cam_regrid_loc", required=True)
    tem_case_locs = [f"{output_loc}/tem"]
    tem_base_loc = f"{output_loc}/tem"

    #If path not specified, skip TEM calculation?
    if tem_case_locs is None:
        print("\t 'cam_tem_loc' not found for test case(s) in config file, so no TEM plots will be generated.")
        return
    else:
        for tem_case_loc in tem_case_locs:
            tem_case_loc = Path(tem_case_loc)
            #Check if TEM directory exists, and if not, then create it:
            if not tem_case_loc.is_dir():
                print(f"    {tem_case_loc} not found, making new directory")
                tem_case_loc.mkdir(parents=True)
            #End if
            tem_locs.append(tem_case_loc)
        #End for

    #Set seasonal ranges:
    seasons = {"ANN": np.arange(1,13,1),
               "DJF": [12, 1, 2],
               "JJA": [6, 7, 8],
               "MAM": [3, 4, 5],
               "SON": [9, 10, 11]
               }

    #Suggestion from Rolando, if QBO is being produced, add utendvtem and utendwtem?
    if "qbo" in adf.plotting_scripts:
        var_list = ['uzm','thzm','epfy','epfz','vtem','wtem',
                    'psitem','utendepfd','utendvtem','utendwtem']
        #var_list = ['uzm','epfy','epfz','vtem','wtem',
        #            'psitem','utendepfd','utendvtem','utendwtem']
    #Otherwise keep it simple
    else:
        var_list = ['uzm','thzm','epfy','epfz','vtem','wtem','psitem','utendepfd']
        #var_list = ['uzm','epfy','epfz','vtem','wtem','psitem','utendepfd']

    #Baseline TEM location
    #input_loc_idx = Path(tem_base_loc)

    #Check if comparing against obs
    if adf.compare_obs:
        obs = True
        #Set TEM file for observations
        base_file_name = 'Obs.TEMdiag.nc'
        input_loc_idx = Path(tem_locs[0])
    else:
        #Set TEM file for baseline
        input_loc_idx = Path(tem_base_loc)
        base_file_name = f'{base_name}.TEMdiag_{syear_baseline}-{eyear_baseline}.nc'
    
    #Set full path for baseline/obs file
    tem_base = input_loc_idx / base_file_name

    #Check to see if baseline/obs TEM file exists    
    if tem_base.is_file():
        ds_base = xr.open_dataset(tem_base)
    else:
        print(f"\t'{base_file_name}' does not exist. TEM plots will be skipped.")
        return

    input_ts_locs = adf.get_cam_info("cam_ts_loc", required=True)

    #Loop over variables:
    for var in var_list:
        """
        if adf.compare_obs:
            #Check if obs exist for the variable:
            if var in var_obs_dict:
                #Note: In the future these may all be lists, but for
                #now just convert the target_list.
                #Extract target file:
                dclimo_loc = var_obs_dict[var]["obs_file"]
                #Extract target list (eventually will be a list, for now need to convert):
                data_list = [var_obs_dict[var]["obs_name"]]
                #Extract target variable name:
                data_var = var_obs_dict[var]["obs_var"]
            else:
                dmsg = f"No obs found for variable `{var}`, zonal mean plotting skipped."
                adfobj.debug_log(dmsg)
                continue
            #End if
        else:
            #Set "data_var" for consistent use below:
            data_var = var
        #End if
        """
        if (adf.compare_obs) and (var == "thzm"):
            print("Obs case is missing potential temperature, so this variable will be skipped.")
            continue

        #Notify user of variable being plotted:
        print(f"\t - TEM plots for {var}")

        #Loop over model cases:
        for idx,case_name in enumerate(case_names):

            """# Check redo_plot. If set to True: remove old plot, if it already exists:
            if (not redo_plot) and plot_name.is_file():
                #Add already-existing plot to website (if enabled):
                adf.debug_log(f"'{plot_name}' exists and clobber is false.")
                adf.add_website_data(plot_name, "TEM", case_name, season=s, plot_type="WACCM",ext="Mean",category="Seasonal Cycle")

                #Continue to next iteration:
                continue
            elif (redo_plot) and plot_name.is_file():
                plot_name.unlink()"""

            tem_loc = tem_case_locs[idx]

            #Extract start and end year values:
            start_year = syear_cases[idx]
            end_year   = eyear_cases[idx]

            #Open the TEM file
            output_loc_idx = Path(tem_loc)
            case_file_name = f'{case_name}.TEMdiag_{start_year}-{end_year}.nc'
            tem_case = output_loc_idx / case_file_name

            #Grab the data for the TEM netCDF files
            if tem_case.is_file():
                ds = xr.open_dataset(tem_case)
            else:
                print(f"\t'{tem_case}' does not exist. TEM plots will be skipped.")
                return

            #Loop over season dictionary:
            for s in seasons:

                #Location to save plots
                plot_name = plot_location / f"{var}_{s}_WACCM_SeasonalCycle_Mean.png"

                # Check redo_plot. If set to True: remove old plot, if it already exists:
                if (not redo_plot) and plot_name.is_file():
                    #Add already-existing plot to website (if enabled):
                    adf.debug_log(f"'{plot_name}' exists and clobber is false.")
                    adf.add_website_data(plot_name, var, None, season=s, plot_type="WACCM",ext="SeasonalCycle_Mean",category="TEM",multi_case=True)

                #plot_name = plot_loc / f"CPT_ANN_WACCM_SeasonalCycle_Mean.{plot_type}"
                elif ((redo_plot) and plot_name.is_file()) or (not plot_name.is_file()):
                    if plot_name.is_file():
                        plot_name.unlink()

                #Grab variable defaults for this variable
                vres = res[var]

                #Gather data for both cases
                mdata = ds[var].squeeze()
                odata = ds_base[var].squeeze()

                # APPLY UNITS TRANSFORMATION IF SPECIFIED:
                # NOTE: looks like our climo files don't have all their metadata
                mdata = mdata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
                # update units
                mdata.attrs['units'] = vres.get("new_unit", mdata.attrs.get('units', 'none'))

                # Do the same for the baseline case if need be:
                if not obs:
                    odata = odata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
                    # update units
                    odata.attrs['units'] = vres.get("new_unit", odata.attrs.get('units', 'none'))
                # Or for observations
                else:
                    odata = odata * vres.get("obs_scale_factor",1) + vres.get("obs_add_offset", 0)
                    # Note: we are going to assume that the specification ensures the conversion makes the units the same. Doesn't make sense to add a different unit.

                #Create array to avoid weighting missing values:
                md_ones = xr.where(mdata.isnull(), 0.0, 1.0)
                od_ones = xr.where(odata.isnull(), 0.0, 1.0)

                month_length = mdata.time.dt.days_in_month
                weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())

                #Calculate monthly-weighted seasonal averages:
                if s == 'ANN':

                    #Calculate annual weights (i.e. don't group by season):
                    weights_ann = month_length / month_length.sum()

                    mseasons = (mdata * weights_ann).sum(dim='time')
                    mseasons = mseasons / (md_ones*weights_ann).sum(dim='time')

                    #Calculate monthly weights based on number of days:
                    if obs:
                        month_length_obs = odata.time.dt.days_in_month
                        weights_ann_obs = month_length_obs / month_length_obs.sum()
                        oseasons = (odata * weights_ann_obs).sum(dim='time')
                        oseasons = oseasons / (od_ones*weights_ann_obs).sum(dim='time')
                    else:
                        month_length_base = odata.time.dt.days_in_month
                        weights_ann_base = month_length_base / month_length_base.sum()
                        oseasons = (odata * weights_ann_base).sum(dim='time')
                        oseasons = oseasons / (od_ones*weights_ann_base).sum(dim='time')

                else:
                    #this is inefficient because we do same calc over and over
                    mseasons = (mdata * weights).groupby("time.season").sum(dim="time").sel(season=s)
                    wgt_denom = (md_ones*weights).groupby("time.season").sum(dim="time").sel(season=s)
                    mseasons = mseasons / wgt_denom

                    if obs:
                        month_length_obs = odata.time.dt.days_in_month
                        weights_obs = (month_length_obs.groupby("time.season") / month_length_obs.groupby("time.season").sum())
                        oseasons = (odata * weights_obs).groupby("time.season").sum(dim="time").sel(season=s)
                        wgt_denom = (od_ones*weights_obs).groupby("time.season").sum(dim="time").sel(season=s)
                        oseasons = oseasons / wgt_denom
                    else:
                        month_length_base = odata.time.dt.days_in_month
                        weights_base = (month_length_base.groupby("time.season") / month_length_base.groupby("time.season").sum())
                        oseasons = (odata * weights_base).groupby("time.season").sum(dim="time").sel(season=s)
                        wgt_denom_base = (od_ones*weights_base).groupby("time.season").sum(dim="time").sel(season=s)
                        oseasons = oseasons / wgt_denom_base

                if var == "thzm":
                    print(f"\t       INFO: deriving zonal mean temperature from potential temperature")


                    """
                    from metpy.calc import temperature_from_potential_temperature
                    # potential temperature
                    theta = np.array([ 286.12859679, 288.22362587]) * units.kelvin
                    p = 850 * units.mbar
                    T = temperature_from_potential_temperature(p, theta)
                    """


                    #path = "/glade/derecho/scratch/richling/adf-output/ADF-data/timeseries/"
                    #path += "f.cam6_3_132.FMTHIST_ne30.sponge.001/1996-2005/"
                    #ds_pmid = xr.open_dataset(path+"f.cam6_3_132.FMTHIST_ne30.sponge.001.cam.h0.PMID.199601-200512.nc")

                    #path = "/glade/derecho/scratch/richling/adf-output/ADF-data/timeseries/"
                    #path += "f.cam6_3_132.FMTHIST_ne30.sponge.001/1996-2005/"
                    path = input_ts_locs[idx]
                    ds_pmid = xr.open_dataset(f"{path}{case_name}.cam.h0.PMID.{start_year}01-{end_year}12.nc")


                    ds_pmid_interp = ds_pmid.interp(lat=mseasons.zalat,method="nearest")
                    pmid = ds_pmid_interp["PMID"]
                    pmid.attrs['units'] = 'Pa'
                    #print(pmid)

                    #Create array to avoid weighting missing values:
                    pmid_ones = xr.where(pmid.isnull(), 0.0, 1.0)

                    #month_length = pmid.time.dt.days_in_month
                    #weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())
                    if s == 'ANN':

                        #Calculate annual weights (i.e. don't group by season):
                        weights_ann = month_length / month_length.sum()

                        pmid = (pmid * weights_ann).sum(dim='time')
                        pmid = pmid / (pmid_ones*weights_ann).sum(dim='time')
                    else:
                        #this is inefficient because we do same calc over and over
                        pmid = (pmid * weights).groupby("time.season").sum(dim="time").sel(season=s)
                        wgt_denom = (pmid_ones*weights).groupby("time.season").sum(dim="time").sel(season=s)
                        pmid = pmid / wgt_denom


                    mseasons.attrs['units'] = "K"
                    oseasons.attrs['units'] = "K"
                    pmid = pmid.mean(dim="lon")
                    #mseasons = thermo.temperature_from_potential_temperature(pmid* units.mbar,mseasons* units.kelvin)
                    #print("AHHH",np.max(mseasons.values))
                    #oseasons = thermo.temperature_from_potential_temperature(pmid* units.mbar,oseasons* units.kelvin)

                    mseasons = thermo.temperature_from_potential_temperature(pmid* units.Pa,mseasons* units.kelvin)
                    #mseasons_metpy = thermo.temperature_from_potential_temperature(pmid* units.Pa,mseasons* units.kelvin)
                    #print("AHHH",np.max(mseasons.values))
                    #oseasons_metpy = thermo.temperature_from_potential_temperature(pmid* units.Pa,oseasons* units.kelvin)
                    oseasons = thermo.temperature_from_potential_temperature(pmid* units.Pa,oseasons* units.kelvin)

                    # exner_function(pressure, reference_pressure=mpconsts.P0)
                    # potential_temperature * exner_function(pressure)
                    #exner_function = (pmid / mconst.P0*100)**mconst.kappa
                    #(pmid / mconst.P0)**mconst.kappa

                    #mseasons = (mseasons * exner_function)/13.894954
                    #oseasons = (oseasons * exner_function)/13.894954


                    #print("mseasons",np.max(mseasons_metpy)==np.max(mseasons))
                    #print("oseasons",np.max(oseasons_metpy)==np.max(oseasons))



                if var == "utendepfd":
                    mseasons = mseasons*1000
                    oseasons = oseasons*1000

                if s ==list(seasons.keys())[0] and var==var_list[0]:
                    print("\n\nmseasons",mseasons,"\n\n")

                    print("mseasons.shape",mseasons.shape)
                    print("oseasons.shape",oseasons.shape,"\n\n")

                test_lons = mseasons.lev
                test_lats = mseasons.zalat

                obs_lons = oseasons.lev
                obs_lats = oseasons.zalat

                """if obs_lons.shape == test_lons.shape:
                    try:
                        xr.testing.assert_equal(test_lons, obs_lons)
                        same_lats = True
                        if s ==list(seasons.keys())[0] and var==var_list[0]:
                            print("the lons ARE the same")
                    except AssertionError as e:
                        same_lons = False
                        if s ==list(seasons.keys())[0] and var==var_list[0]:
                            print("the lons aren't the same")
                    try:
                        xr.testing.assert_equal(test_lats, obs_lats)
                        same_lats = True
                        if s ==list(seasons.keys())[0] and var==var_list[0]:
                            print("the lats ARE the same")
                    except AssertionError as e:
                        same_lats = False
                        if s ==list(seasons.keys())[0] and var==var_list[0]:
                            print("the lats aren't the same")
                else:
                    same_lats = False
                    same_lons = False
                    if s ==list(seasons.keys())[0] and var==var_list[0]:
                        print("The ensemble array lat/lon shape does not match the " \
                            "obs mask array.\nRegridding to ensemble lats and lons")"""



                """#if (not same_lats) and (not same_lons):
                if 0 == 1:

                    # Define standard pressure levels (vertical target levels)
                    standard_lev = np.array([1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50,
                                            30, 20, 10, 7, 5, 3, 2, 1])

                    # Determine which dataset has more vertical levels
                    if len(mseasons.lev) > len(oseasons.lev):
                        if s == list(seasons.keys())[0] and var == var_list[0]:
                            print("Source data is oseasons")
                        source_data = oseasons
                        target_data = mseasons
                    else:
                        if s == list(seasons.keys())[0] and var == var_list[0]:
                            print("Source data is mseasons")
                        source_data = mseasons
                        target_data = oseasons

                    source_data = oseasons
                    target_data = mseasons

                    # Extract source and target coordinates
                    source_lat = source_data.zalat.values
                    source_lev = source_data.lev.values
                    source_values = source_data.values

                    target_lat = target_data.zalat.values
                    target_lev = target_data.lev.values
                    target_values = target_data.values

                    ### Step 1: Interpolate Source Data to Standard Vertical Levels ###
                    source_interpolator = RegularGridInterpolator((source_lev, source_lat), source_values, bounds_error=False, fill_value=np.nan)

                    # Generate meshgrid points for interpolation
                    source_points = np.array(np.meshgrid(standard_lev, source_lat, indexing='ij')).reshape(2, -1).T
                    source_regridded_values_vert = source_interpolator(source_points).reshape(len(standard_lev), len(source_lat))

                    ### Step 2: Interpolate Target Data to Standard Vertical Levels ###
                    target_interpolator = RegularGridInterpolator((target_lev, target_lat), target_values, bounds_error=False, fill_value=np.nan)

                    # Generate meshgrid points for interpolation
                    target_points = np.array(np.meshgrid(standard_lev, target_lat, indexing='ij')).reshape(2, -1).T
                    target_regridded_values_vert = target_interpolator(target_points).reshape(len(standard_lev), len(target_lat))

                    ### Step 3: Interpolate Source Data Horizontally ###
                    source_regridded_values_horiz = []
                    for i, lev in enumerate(standard_lev):
                        interpolator = RegularGridInterpolator((source_lat,), source_regridded_values_vert[i, :], bounds_error=False, fill_value=np.nan)
                        regridded_values = interpolator(source_lat)
                        source_regridded_values_horiz.append(regridded_values)
                    source_regridded_values_horiz = np.array(source_regridded_values_horiz)

                    ### Step 4: Interpolate Target Data Horizontally ###
                    target_regridded_values_horiz = []
                    for i, lev in enumerate(standard_lev):
                        interpolator = RegularGridInterpolator((target_lat,), target_regridded_values_vert[i, :], bounds_error=False, fill_value=np.nan)
                        regridded_values = interpolator(source_lat)
                        target_regridded_values_horiz.append(regridded_values)
                    target_regridded_values_horiz = np.array(target_regridded_values_horiz)

                    ### Step 5: Convert Regridded Data Back into xarray.DataArray ###
                    source_regridded_data = xr.DataArray(
                        data=source_regridded_values_horiz,
                        dims=["lev", "zalat"],
                        coords={"lev": standard_lev, "zalat": source_lat},
                        name="source_regridded_data"
                    )

                    target_regridded_data = xr.DataArray(
                        data=target_regridded_values_horiz,
                        dims=["lev", "zalat"],
                        coords={"lev": standard_lev, "zalat": source_lat},
                        name="target_regridded_data"
                    )

                    ### Step 6: Assign Back to Variables ###
                    oseasons = source_regridded_data
                    mseasons = target_regridded_data
                    lat = mseasons["zalat"]
                    lev = mseasons["lev"]

                    ### Debug Output ###
                    if s == list(seasons.keys())[0] and var == var_list[0]:
                        print("Source Regridded Data:")
                        print(source_regridded_data, "\n\n")

                        print("Target Regridded Data:")
                        print(target_regridded_data, "\n\n")
                else:
                    lat = mseasons['zalat']
                    lev = mseasons['lev']"""
                    
                lat = mseasons['zalat']
                lev = mseasons['lev']
                #difference: each entry should be (lat, lon)
                dseasons = mseasons-oseasons
                
                #Gather contour plot options
                cp_info = pf.prep_contour_plot(mseasons, oseasons, dseasons, **vres)
                clevs = np.unique(np.array(cp_info['levels1']))

                norm = cp_info['norm1']
                cmap = cp_info['cmap1']
                clevs_diff = np.unique(np.array(cp_info['levelsdiff']))

                # mesh for plots:
                lats, levs = np.meshgrid(lat, lev)

                # Find the next value below highest vertical level
                prev_major_tick = 10 ** (np.floor(np.log10(np.min(levs))))
                prev_major_tick

                # Set padding for colorbar form axis
                cmap_pad = 0.005

                # create figure object
                fig = plt.figure(figsize=(14,10))
                # LAYOUT WITH GRIDSPEC
                # 4 rows, 8 columns, but each map will take up 4 columns and 2 rows
                gs = mpl.gridspec.GridSpec(4, 8, wspace=0.75,hspace=0.5)
                ax1 = plt.subplot(gs[0:2, :4], **cp_info['subplots_opt'])
                ax2 = plt.subplot(gs[0:2, 4:], **cp_info['subplots_opt'])
                ax3 = plt.subplot(gs[2:, 2:6], **cp_info['subplots_opt'])
                ax = [ax1,ax2,ax3]

                #Contour fill
                img0 = ax[0].contourf(lats, levs,mseasons, levels=clevs, norm=norm, cmap=cmap)
                img1 = ax[1].contourf(lats, levs,oseasons, levels=clevs, norm=norm, cmap=cmap)
                    
                #Add contours for highlighting
                c0 = ax[0].contour(lats,levs,mseasons,levels=clevs[::2], norm=norm,
                                    colors="k", linewidths=0.5)

                #Check if contour labels need to be adjusted
                #ie if the values are large and/or in scientific notation, just label the 
                #contours with the leading numbers.
                #EXAMPLE: plot values are 200000; plot the contours as 2.0 and let the colorbar
                #         indicate that it is e5.
                fmt = {}
                if 'contour_adjust' in vres:
                    test_strs = c0.levels/float(vres['contour_adjust'])
                    for l, str0 in zip(c0.levels, test_strs):
                        fmt[l] = str0

                    # Add contour labels
                    plt.clabel(c0, inline=True, fontsize=8, levels=c0.levels, fmt=fmt)
                else:
                    # Add contour labels
                    plt.clabel(c0, inline=True, fontsize=8, levels=c0.levels)

                #Add contours for highlighting
                c1 = ax[1].contour(lats,levs,oseasons,levels=clevs[::2], norm=norm,
                                    colors="k", linewidths=0.5)

                #Check if contour labels need to be adjusted
                #ie if the values are large and/or in scientific notation, just label the 
                #contours with the leading numbers.
                #EXAMPLE: plot values are 200000; plot the contours as 2.0 and let the colorbar
                #         indicate that it is e5.
                fmt = {}
                if 'contour_adjust' in vres:
                    base_strs = c1.levels/float(vres['contour_adjust'])
                    for l, str0 in zip(c1.levels, base_strs):
                        fmt[l] = str0

                    # Add contour labels
                    plt.clabel(c1, inline=True, fontsize=8, levels=c1.levels, fmt=fmt)
                else:
                    # Add contour labels
                    plt.clabel(c1, inline=True, fontsize=8, levels=c1.levels)


                #Check if difference plot has contour levels, if not print notification
                if len(dseasons.lev) == 0:
                    #Set empty message for comparison of cases with different vertical levels
                    #TODO: Work towards getting the vertical and horizontal interpolations!! - JR
                    empty_message = "These have different vertical levels\nCan't compare cases currently"
                    props = {'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.9}
                    prop_x = 0.18
                    prop_y = 0.42
                    ax[2].text(prop_x, prop_y, empty_message,
                                    transform=ax[2].transAxes, bbox=props)
                else:
                    img2 = ax[2].contourf(lats, levs, dseasons,
                                            #cmap="BrBG",
                                            cmap=cp_info['cmapdiff'],
                                            levels=clevs_diff,
                                            norm=cp_info['normdiff'])
                    ax[2].contour(lats, levs, dseasons, colors="k", linewidths=0.5,
                                    levels=clevs_diff[::2], norm=cp_info['normdiff'])
                    cp_info['diff_colorbar_opt']["label"] = cp_info['colorbar_opt']["label"]
                    plt.colorbar(img2, ax=ax[2], location='right', pad=cmap_pad,**cp_info['diff_colorbar_opt'])

                #Format y-axis
                #for a in ax[:]:
                for i,a in enumerate(ax[:]):
                    a.set_yscale("log")
                    a.set_xlabel("Latitude")
                    # Only plot y-axis label for test case
                    if i == 0:
                        a.set_ylabel('Pressure [hPa]', va='center', rotation='vertical')
                    if 'ylim' in vres:
                        y_lims = [float(lim) for lim in vres['ylim']]

                        #print("y_lims",y_lims,"\n")
                        np.min(levs)
                        y_lims[-1]=prev_major_tick #np.min(levs)
                        #print("y_lims",y_lims,"\n")
                        a.set_ylim(y_lims)
                    else:
                        a.set_ylim(a.get_ylim()[::-1])

                # Format color bars
                #print("cp_info['colorbar_opt']",cp_info['colorbar_opt'],"\n")
                plt.colorbar(img1, ax=ax[1], location='right', pad=cmap_pad,**cp_info['colorbar_opt'])
                # Remove the colorbar label for baseline
                cp_info['colorbar_opt'].pop("label", None)
                plt.colorbar(img0, ax=ax[0], location='right', pad=cmap_pad,**cp_info['colorbar_opt'])

                #Set titles of subplots
                #Set figure title
                #plt.suptitle(f'TEM Diagnostics: {s}', fontsize=20, y=.98)

                #Variable plot title name
                longname = vres["long_name"]
                #plt.text(0.5, 0.915, f"{longname}\n", fontsize=12, ha='center',
                #            transform=fig.transFigure)

                plt.suptitle(f'{longname}: {s}', fontsize=20, y=.97)

                test_yrs = f"{start_year}-{end_year}"
                #ax[0].set_title(f"{test_nicknames[idx]}\n{test_yrs}",fontsize=10)

                
                plot_title = "$\mathbf{Test}:$"+f"{test_nicknames[idx]}\nyears: {test_yrs}"
                ax[0].set_title(plot_title, loc='left', fontsize=10)
                #ax[idx].set_title(plot_title, loc='left', fontsize=10)

                if obs:
                    obs_title = Path(vres["obs_name"]).stem
                    ax[1].set_title(f"{obs_title}\n",fontsize=10)

                else:
                    base_yrs = f"{syear_baseline}-{eyear_baseline}"
                    plot_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {base_yrs}"
                    ax[1].set_title(plot_title, loc='left', fontsize=10)
                
                #Set main title for difference plots column
                ax[2].set_title("$\mathbf{Test} - \mathbf{Baseline}$",fontsize=10)

                #Write the figure to provided workspace/file:
                fig.savefig(plot_name, bbox_inches='tight', dpi=300)

                #Add plot to website (if enabled):
                adf.add_website_data(plot_name, var, case_name, season=s, plot_type="WACCM",
                                     ext="SeasonalCycle_Mean",category="TEM")

                plt.close()
    print("  ...TEM plots have been generated successfully.")

# Helper functions
##################

import xesmf as xe



def interp_tem(arr_anom1, arr_anom2):
    """
    Check if the Obs array needs to be interpolated
    to the ensemble file

    Most likely the input array will need to be interpolated!
    """

    #arr_anom1 = arr_anom1.rename(lev="lon", zalat="lat")
    #arr_anom2 = arr_anom2.rename(lev="lon", zalat="lat")
    import numpy as np

    # Ensure the source data is contiguous
    #if not arr_anom1.values.flags['C_CONTIGUOUS']:
    #    arr_anom1 = arr_anom1.copy(data=np.ascontiguousarray(arr_anom1.values))

    test_lons = arr_anom1.lon
    test_lats = arr_anom1.lat

    obs_lons = arr_anom2.lon
    obs_lats = arr_anom2.lat

    # Just set these to true for now
    same_lats = True
    same_lons = True

    arr_prime = None

    if obs_lons.shape == test_lons.shape:
        try:
            xr.testing.assert_equal(test_lons, obs_lons)
            print("the lons ARE the same")
        except AssertionError as e:
            same_lons = False
            print("the lons aren't the same")
        try:
            xr.testing.assert_equal(test_lats, obs_lats)
            print("the lats ARE the same")
        except AssertionError as e:
            same_lats = False
            print("the lats aren't the same")
    else:
        same_lats = False
        same_lons = False
        print("The ensemble array lat/lon shape does not match the " \
             "obs mask array.\nRegridding to ensemble lats and lons")
        return arr_anom1

    if (not same_lons) and (not same_lats):

        ds_out = xr.Dataset(
            {   
                "lon": (["lon"], obs_lons.values, {"units": "degrees_east"}),
                "lat": (["lat"], obs_lats.values, {"units": "degrees_north"}),
            }
        )

        # Regrid to the ensemble grid to make altered obs grid
        regridder = xe.Regridder(arr_anom1, ds_out, "bilinear", periodic=True)
        arr_prime = regridder(arr_anom1, keep_attrs=True)

    # Return the new interpolated obs array
    return arr_prime






'''def interp_tem(arr_anom1, arr_anom2):
    """
    Calculate seasonal averages and regrid the seasonal data to match the target grid,
    only if the levels (lev) or latitudes (zalat) are different between the two datasets.
    """
    # Step 1: Compute seasonal averages for arr_anom1 (ensemble) and arr_anom2 (observations)
    
    # Resample to seasonal averages: We assume monthly data, resampling to DJF, MAM, JJA, SON
    #seasonal_anom1 = arr_anom1.resample(time='QS-DEC').mean()  # 'QS-DEC' starts the season in December (for DJF)
    #seasonal_anom2 = arr_anom2.resample(time='QS-DEC').mean()

    # Step 2: Check if levs and zalats are the same between the two datasets
    same_levs = xr.DataArray.equals(arr_anom1.lev, arr_anom2.lev)
    same_zalats = xr.DataArray.equals(arr_anom1.zalat, arr_anom2.zalat)

    # Step 3: If both levs and zalats are the same, no regridding is needed
    if same_levs and same_zalats:
        print("The levels (lev) and zalats are the same. No regridding required.")
        return arr_anom1  # Return the seasonal data as is

    # Step 4: If levs or zalats are different, create a new output grid (target grid) from arr_anom2
    print("The levels (lev) or zalats are different. Proceeding with regridding.")
    ds_out = xr.Dataset(
        {
            "zalat": (["zalat"], arr_anom2.zalat.values, {"units": "degrees_north"}),
            "lev": (["lev"], arr_anom2.lev.values, {"units": "hPa"}),
        }
    )
    arr_anom1 = arr_anom1.rename(lev="lat", zalat="lon")
    arr_anom2 = arr_anom2.rename(lev="lat", zalat="lon")
    print("arr_anom1.shape",arr_anom1.shape)
    print("arr_anom2.shape",arr_anom2.shape)

    # Step 5: Apply the regridding
    regridder = xe.Regridder(arr_anom1, arr_anom2, method="bilinear")

    # Step 6: Regrid the seasonal data from arr_anom1 to the target grid
    anom1_prime = regridder(arr_anom1, keep_attrs=True)
    
    # Step 7: Return the regridded seasonal anomalies
    return anom1_prime'''







'''def interp_tem(arr_anom1, arr_anom2):
    """
    Check if the Obs array needs to be interpolated
    to the ensemble file

    Most likely the input array will need to be interpolated!
    """
    import xesmf as xe
    test_levs = arr_anom1.lev
    test_zalats = arr_anom1.zalat

    obs_levs = arr_anom2.lev
    obs_zalats = arr_anom2.zalat

    # Just set these to true for now
    same_zalats = True
    same_levs = True

    arr_prime = None

    if obs_levs.shape == test_levs.shape:
        try:
            xr.testing.assert_equal(test_levs, obs_levs)
            print("the lons ARE the same")
        except AssertionError as e:
            same_levs = False
            print("the lons aren't the same")
        try:
            xr.testing.assert_equal(test_zalats, obs_zalats)
            print("the zalats ARE the same")
        except AssertionError as e:
            same_zalats = False
            print("the zalats aren't the same")
    else:
        same_zalats = False
        same_levs = False
        print("The ensemble array lev/zalats shape does not match the " \
             "obs mask array.\nRegridding to ensemble lats and lons")

    if (not same_levs) and (not same_zalats):

        ds_out = xr.Dataset(
            {
                "zalat": (["lat"], obs_zalats.values, {"units": "degrees_north"}),
                "lev": (["lon"], obs_levs.values, {"units": "hPa"}),
            }
        )

        # Regrid to the ensemble grid to make altered obs grid
        regridder = xe.Regridder(arr_anom1, ds_out, "bilinear")
        arr_prime = regridder(arr_anom1, keep_attrs=True)

    # Return the new interpolated obs array
    return arr_prime

#######'''