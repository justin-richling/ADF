#Import standard modules:
import xarray as xr

def regrid_and_vert_interp(adf):

    """
    This funtion regrids the test cases to the same horizontal
    grid as the observations or baseline climatology.  It then
    vertically interpolates the test case (and baseline case
    if need be) to match a default set of pressure levels, which
    are (in hPa):

    1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50,
    30, 20, 10, 7, 5, 3, 2, 1

    Currently any 3-D observations file needs to have equivalent pressure
    levels in order to work properly, although in the future it is hoped
    to enable the vertical interpolation of observations as well.

    Description of needed inputs from ADF:

    case_name        -> Name of CAM case provided by "cam_case_name"
    input_climo_loc  -> Location of CAM climo files provided by "cam_climo_loc"
    output_loc       -> Location to write re-gridded CAM files, specified by "cam_regrid_loc"
    var_list         -> List of CAM output variables provided by "diag_var_list"
    var_defaults     -> Dict that has keys that are variable names and values that are plotting preferences/defaults.
    target_list      -> List of target data sets CAM could be regridded to
    taget_loc        -> Location of target files that CAM will be regridded to
    overwrite_regrid -> Logical to determine if already existing re-gridded
                        files will be overwritten. Specified by "cam_overwrite_regrid"
    """

    #Import necessary modules:
    import numpy as np
    import plotting_functions as pf

    from pathlib import Path
    import os

    # regridding
    # Try just using the xarray method
    # import xesmf as xe  # This package is for regridding, and is just one potential solution.

    # Steps:
    # - load climo files for model and obs
    # - calculate all-time and seasonal fields (from individual months)
    # - regrid one to the other (probably should be a choice)

    #Notify user that script has started:
    print("\n  Regridding CAM climatologies...")

    #Extract needed quantities from ADF object:
    #-----------------------------------------
    overwrite_regrid = adf.get_basic_info("cam_overwrite_regrid", required=True)
    output_loc       = adf.get_basic_info("cam_regrid_loc", required=True)

    #CAM simulation variables (these quantities are always lists):
    case_names = adf.get_cam_info("cam_case_name", required=True)
    input_climo_locs = adf.get_cam_info("cam_climo_loc", required=True)

    #Grab case years
    syear_cases = adf.climo_yrs["syears"]
    eyear_cases = adf.climo_yrs["eyears"]

    #Extract model baseline variables:
    target_loc = adf.get_baseline_info("cam_climo_loc", required=True)
    target_list = [adf.get_baseline_info("cam_case_name", required=True)]

    obs_dir = adf.get_basic_info("obs_data_loc")
    file_merra2 = os.path.join(obs_dir, 'MERRA2_192x288_AOD_2001-2020_climo.nc')
    file_mod08_m3 = os.path.join(obs_dir, 'MOD08_M3_192x288_AOD_2001-2020_climo.nc')

    target_list = []

    #Grab baseline years (which may be empty strings if using Obs):
    syear_baseline = adf.climo_yrs["syear_baseline"]
    eyear_baseline = adf.climo_yrs["eyear_baseline"]

    #Set attributes dictionary for climo years to save in the file attributes
    base_climo_yrs_attr = f"{target_list[0]}: {syear_baseline}-{eyear_baseline}"

    #-----------------------------------------

    #Set output/target data path variables:
    #------------------------------------
    rgclimo_loc = Path(output_loc)
    if not adf.compare_obs:
        tclimo_loc  = Path(target_loc)
    #------------------------------------

    #Check if re-gridded directory exists, and if not, then create it:
    if not rgclimo_loc.is_dir():
        print(f"    {rgclimo_loc} not found, making new directory")
        rgclimo_loc.mkdir(parents=True)
    #End if

    #Loop over CAM cases:
    for case_idx, case_name in enumerate(case_names):

        #Notify user of model case being processed:
        print(f"\t Regridding case '{case_name}' :")

        #Set case climo data path:
        mclimo_loc  = Path(input_climo_locs[case_idx])


        #Get climo years for case
        syear = syear_cases[case_idx]
        eyear = eyear_cases[case_idx]

        # probably want to do this one variable at a time:
        for var in ["AODVISdn"]:

            #Notify user of variable being regridded:
            print(f"\t - regridding {var} (known targets: {target_list})")

            ds_case = adf.data.load_climo_da(case_name, var)

            #loop over regridding targets:
            for target in target_list:

                #Write to debug log if enabled:
                adf.debug_log(f"regrid_example: regrid target = {target}")

                #Determine regridded variable file name:
                regridded_file_loc = rgclimo_loc / f'{target}_{case_name}_{var}_regridded.nc'


                #Check if re-gridded file already exists and over-writing is allowed:
                if regridded_file_loc.is_file() and overwrite_regrid:
                    #If so, then delete current file:
                    regridded_file_loc.unlink()
                #End if

                #Check again if re-gridded file already exists:
                if not regridded_file_loc.is_file():

                    #Create list of regridding target files (we should explore intake as an alternative to having this kind of repeated code)
                    # NOTE: This breaks if you have files from different cases in same directory!
                    if adf.compare_obs:
                        #For now, only grab one file (but convert to list for use below):
                        tclim_fils = [tclimo_loc]
                    else:
                       tclim_fils = sorted(tclimo_loc.glob(f"{target}*_{var}_climo.nc"))
                    #End if

                    #Write to debug log if enabled:
                    adf.debug_log(f"regrid_example: tclim_fils (n={len(tclim_fils)}): {tclim_fils}")

                    if len(tclim_fils) > 1:
                        #Combine all target files together into a single data set:
                        tclim_ds = xr.open_mfdataset(tclim_fils, combine='by_coords')
                    elif len(tclim_fils) == 0:
                        print(f"\t - regridding {var} failed, no file. Continuing to next variable.")
                        continue
                    else:
                        #Open single file as new xarray dataset:
                        tclim_ds = xr.open_dataset(tclim_fils[0])
                    #End if

                    #Generate CAM climatology (climo) file list:
                    mclim_fils = sorted(mclimo_loc.glob(f"{case_name}_{var}_*.nc"))

                    if len(mclim_fils) > 1:
                        #Combine all cam files together into a single data set:
                        mclim_ds = xr.open_mfdataset(mclim_fils, combine='by_coords')
                    elif len(mclim_fils) == 0:
                        wmsg = f"\t - Unable to find climo file for '{var}'."
                        wmsg += " Continuing to next variable."
                        print(wmsg)
                        continue
                    else:
                        #Open single file as new xarray dataset:
                        mclim_ds = xr.open_dataset(mclim_fils[0])
                    #End if

                    #Create keyword arguments dictionary for regridding function:
                    regrid_kwargs = {}


                    #Perform regridding and interpolation of variable:
                    rgdata_interp = _regrid_and_interpolate_levs(mclim_ds, var,
                                                                 regrid_dataset=tclim_ds,
                                                                 **regrid_kwargs)

                    #Finally, write re-gridded data to output file:
                    #Convert the list of Path objects to a list of strings
                    climatology_files_str = [str(path) for path in mclim_fils]
                    climatology_files_str = ', '.join(climatology_files_str)
                    test_attrs_dict = {
                            "adf_user": adf.user,
                            "climo_yrs": f"{case_name}: {syear}-{eyear}",
                            "climatology_files": climatology_files_str,
                        }
                    rgdata_interp = rgdata_interp.assign_attrs(test_attrs_dict)
                    save_to_nc(rgdata_interp, regridded_file_loc)
                    rgdata_interp.close()  # bpm: we are completely done with this data

                    #Now vertically interpolate baseline (target) climatology,
                    #if applicable:

                    #Set interpolated baseline file name:
                    interp_bl_file = rgclimo_loc / f'{target}_{var}_baseline.nc'

                    if not adf.compare_obs and not interp_bl_file.is_file():

                        #Look for a baseline climo file for surface pressure (PS):
                        bl_ps_fil = tclimo_loc / f'{target}_PS_climo.nc'

                        #Also look for a baseline climo file for mid-level pressure (PMID):
                        bl_pmid_fil = tclimo_loc / f'{target}_PMID_climo.nc'

                        #Create new keyword arguments dictionary for regridding function:
                        regrid_kwargs = {}

                        #Check if PS and PMID files exist:
                        if bl_ps_fil.is_file():
                            regrid_kwargs.update({'ps_file': bl_ps_fil})
                        #End if
                        if bl_pmid_fil.is_file():
                            regrid_kwargs.update({'pmid_file': bl_pmid_fil})
                        #End if

                        #Generate vertically-interpolated baseline dataset:
                        tgdata_interp = _regrid_and_interpolate_levs(tclim_ds, var,
                                                                     **regrid_kwargs)

                        if tgdata_interp is None:
                            #Something went wrong during interpolation, so just cycle through
                            #for now:
                            continue
                        #End if

                        # Convert the list to a string (join with commas or another separator)
                        climatology_files_str = [str(path) for path in tclim_fils]
                        climatology_files_str = ', '.join(climatology_files_str)
                        # Create a dictionary of attributes
                        base_attrs_dict = {
                            "adf_user": adf.user,
                            "climo_yrs": f"{case_name}: {syear}-{eyear}; {base_climo_yrs_attr}",
                            "climatology_files": climatology_files_str,
                        }
                        tgdata_interp = tgdata_interp.assign_attrs(base_attrs_dict)

                        #Write interpolated baseline climatology to file:
                        save_to_nc(tgdata_interp, interp_bl_file)
                    #End if
                else:
                    print("\t Regridded file already exists, so skipping...")
                #End if (file check)
            #End do (target list)
        #End do (variable list)
    #End do (case list)

    #Notify user that script has ended:
    print("  ...CAM climatologies have been regridded successfully.")

#################
#Helper functions
#################

def _regrid_and_interpolate_levs(model_dataset, var_name, regrid_dataset=None, regrid_ofrac=False, **kwargs):

    """
    Function that takes a variable from a model xarray
    dataset, regrids it to another dataset's lat/lon
    coordinates (if applicable), and then interpolates
    it vertically to a set of pre-defined pressure levels.
    ----------
    model_dataset -> The xarray dataset which contains the model variable data
    var_name      -> The name of the variable to be regridded/interpolated.

    Optional inputs:

    ps_file        -> A NetCDF file containing already re-gridded surface pressure
    regrid_dataset -> The xarray dataset that contains the lat/lon grid that
                      "var_name" will be regridded to.  If not present then
                      only the vertical interpolation will be done.

    kwargs         -> Keyword arguments that contain paths to surface pressure
                      and mid-level pressure files, which are necessary for
                      certain types of vertical interpolation.

    This function returns a new xarray dataset that contains the regridded
    and/or vertically-interpolated model variable.
    """

    #Import ADF-specific functions:
    import plotting_functions as pf


    #Extract variable info from model data (and remove any degenerate dimensions):
    mdata = model_dataset[var_name].squeeze()
   

    #Regrid variable to target dataset (if available):
    if regrid_dataset:

        #Extract grid info from target data:
        if 'time' in regrid_dataset.coords:
            tgrid = regrid_dataset.isel(time=0).squeeze()
        #End if

        #End if
    else:
        #Just rename variables:
        rgdata = mdata
    #End if


   
    #Just rename variable:
    rgdata_interp = rgdata


    #Convert to xarray dataset:
    rgdata_interp = rgdata_interp.to_dataset()

    #Return dataset:
    return rgdata_interp

#####

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

def regrid_data(fromthis, tothis, method=1):
    """Regrid data using various different methods"""

    if method == 1:
        # kludgy: spatial regridding only, seems like can't automatically deal with time
        if 'time' in fromthis.coords:
            result = [fromthis.isel(time=t).interp_like(tothis) for t,time in enumerate(fromthis['time'])]
            result = xr.concat(result, 'time')
            return result
        else:
            return fromthis.interp_like(tothis)
    elif method == 2:
        newlat = tothis['lat']
        newlon = tothis['lon']
        coords = dict(fromthis.coords)
        coords['lat'] = newlat
        coords['lon'] = newlon
        return fromthis.interp(coords)
    elif method == 3:
        newlat = tothis['lat']
        newlon = tothis['lon']
        ds_out = xr.Dataset({'lat': newlat, 'lon': newlon})
        regridder = xe.Regridder(fromthis, ds_out, 'bilinear')
        return regridder(fromthis)
    elif method==4:
        # geocat
        newlat = tothis['lat']
        newlon = tothis['lon']
        result = geocat.comp.linint2(fromthis, newlon, newlat, False)
        result.name = fromthis.name
        return result
    #End if

#####