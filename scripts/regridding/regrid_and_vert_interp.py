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

    # regridding
    # Try just using the xarray method
    # import xesmf as xe  # This package is for regridding, and is just one potential solution.

    # Steps:
    # - load climo files for model and obs
    # - calculate all-time and seasonal fields (from individual months)
    # - regrid one to the other (probably should be a choice)

    #Notify user that script has started:
    msg = "\n  Regridding CAM climatologies..."
    print(f"{msg}\n  {'-' * (len(msg)-3)}")

    #Extract needed quantities from ADF object:
    #-----------------------------------------
    overwrite_regrid_locs = adf.get_cam_info("cam_overwrite_climo_regrid", required=True)
    test_output_loc       = adf.get_cam_info("cam_climo_regrid_loc", required=True)
    var_list         = adf.diag_var_list
    var_defaults     = adf.variable_defaults

    comp = adf.model_component
    if comp == "atm":
        spec_vars = ["PMID", "OCNFRAC", "LANDFRAC"]
    if comp == "lnd":
        spec_vars = ["LANDFRAC"]

    #CAM simulation variables (these quantities are always lists):
    case_names = adf.get_cam_info("cam_case_name", required=True)
    input_climo_locs = adf.get_cam_info("cam_climo_loc", required=True)
    unstruct_cases = adf.unstructs['unstruct_tests']

    case_latlon_files   = adf.latlon_files["test_latlon_file"]
    #print("case_latlon_file",case_latlon_file,"\n")
    #Check if any a weights file exists if using native grid, OPTIONAL
    case_wgts_files   = adf.latlon_wgt_files["test_wgts_file"]
    case_methods = adf.latlon_regrid_method["test_regrid_method"]

    #Grab case years
    syear_cases = adf.climo_yrs["syears"]
    eyear_cases = adf.climo_yrs["eyears"]

    #Check if mid-level pressure, ocean fraction or land fraction exist
    #in the variable list:
    for var in spec_vars:
        if var in var_list:
            #If so, then move them to the front of variable list so
            #that they can be used to mask or vertically interpolate
            #other model variables if need be:
            var_idx = var_list.index(var)
            var_list.pop(var_idx)
            var_list.insert(0,var)
        #End if
    #End for

    #Create new variables that potentially stores the re-gridded
    #ocean/land fraction dataset:
    frc_ds = None
    tgt_frc_ds = None

    #Check if surface pressure exists in variable list:
    if "PS" in var_list:
        #If so, then move it to front of variable list so that
        #it can be used to vertically interpolate model variables
        #if need be.  This should be done after PMID so that the order
        #is PS, PMID, other variables:
        ps_idx = var_list.index("PS")
        var_list.pop(ps_idx)
        var_list.insert(0,"PS")
    #End if

    #Regrid target variables (either obs or a baseline run):
    if adf.compare_obs:

        #Set obs name to match baseline (non-obs)
        target_list = ["Obs"]

        #Extract variable-obs dictionary:
        var_obs_dict = adf.var_obs_dict

        #If dictionary is empty, then  there are no observations to regrid to,
        #so quit here:
        if not var_obs_dict:
            print("\t No observations found to regrid to, so no re-gridding will be done.")
            return
        #End if

    else:

        #Extract model baseline variables:
        target_loc = adf.get_baseline_info("cam_climo_loc", required=True)
        target_list = [adf.get_baseline_info("cam_case_name", required=True)]
        trgclimo_loc = Path(adf.get_baseline_info("cam_climo_regrid_loc", required=True))
        unstruct_base = adf.unstructs['unstruct_base']
        #Check if re-gridded directory exists, and if not, then create it:
        if not trgclimo_loc.is_dir():
            print(f"    {trgclimo_loc} not found, making new directory")
            trgclimo_loc.mkdir(parents=True)
        #End if
    #End if

    #Grab baseline years (which may be empty strings if using Obs):
    syear_baseline = adf.climo_yrs["syear_baseline"]
    eyear_baseline = adf.climo_yrs["eyear_baseline"]

    #Set attributes dictionary for climo years to save in the file attributes
    base_climo_yrs_attr = f"{target_list[0]}: {syear_baseline}-{eyear_baseline}"

    #-----------------------------------------

    #Set output/target data path variables:
    #------------------------------------
    if not adf.compare_obs:
        tclimo_loc  = Path(target_loc)
    #------------------------------------


    #Loop over CAM cases:
    for case_idx, case_name in enumerate(case_names):

        #Notify user of model case being processed:
        print(f"\t Regridding case '{case_name}' :")

        rgclimo_loc = Path(test_output_loc[case_idx])
        #Check if re-gridded directory exists, and if not, then create it:
        if not rgclimo_loc.is_dir():
            print(f"    {rgclimo_loc} not found, making new directory")
            rgclimo_loc.mkdir(parents=True)
        #End if

        overwrite_mregrid = overwrite_regrid_locs[case_idx]

        #Set case climo data path:
        mclimo_loc  = Path(input_climo_locs[case_idx])

        #Create empty dictionaries which store the locations of regridded surface
        #pressure and mid-level pressure fields if needed:
        ps_loc_dict = {}
        pmid_loc_dict = {}

        #Get climo years for case
        syear = syear_cases[case_idx]
        eyear = eyear_cases[case_idx]

        

        # probably want to do this one variable at a time:
        for var in var_list:

            if adf.compare_obs:
                #Check if obs exist for the variable:
                if var in var_obs_dict:
                    #Note: In the future these may all be lists, but for
                    #now just convert the target_list.
                    #Extract target file:
                    tclimo_loc = var_obs_dict[var]["obs_file"]
                    #Extract target list (eventually will be a list, for now need to convert):
                    target_list = [var_obs_dict[var]["obs_name"]]
                else:
                    dmsg = f"No obs found for variable `{var}`, regridding skipped."
                    adf.debug_log(dmsg)
                    continue
                #End if
            #End if

            #Notify user of variable being regridded:
            print(f"\t - regridding {var} (known targets: {target_list})")

            #loop over regridding targets:
            for target in target_list:

                #Write to debug log if enabled:
                adf.debug_log(f"regrid_example: regrid target = {target}")

                #Determine regridded variable file name:
                print("target",target,"\n")
                regridded_file_loc = rgclimo_loc / f'{target}_{case_name}_{var}_regridded.nc'

                if comp == "atm":
                    #If surface or mid-level pressure, then save for potential use by other variables:
                    if var == "PS":
                        ps_loc_dict[target] = regridded_file_loc
                    elif var == "PMID":
                        pmid_loc_dict[target] = regridded_file_loc
                    #End if

                #Check if re-gridded file already exists and over-writing is allowed:
                if regridded_file_loc.is_file() and overwrite_mregrid:
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
                        print(f"\t    WARNING: regridding {var} failed, no climo file for case '{target}'. Continuing to next variable.")
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
                        #wmsg = f"\t    WARNING: Unable to find climo file for '{var}'."
                        #wmsg += " Continuing to next variable."
                        wmsg= f"\t    WARNING: regridding {var} failed, no climo file for case '{case_name}'. Continuing to next variable."
                        print(wmsg)
                        continue
                    else:
                        #Open single file as new xarray dataset:
                        mclim_ds = xr.open_dataset(mclim_fils[0])
                    #End if

                    #Create keyword arguments dictionary for regridding function:
                    regrid_kwargs = {}
                    native_regrid_kwargs = {}

                    if comp == "atm":
                        #Check if target in relevant pressure variable dictionaries:
                        if target in ps_loc_dict:
                            regrid_kwargs.update({'ps_file': ps_loc_dict[target]})
                        #End if
                        if target in pmid_loc_dict:
                            regrid_kwargs.update({'pmid_file': pmid_loc_dict[target]})
                        #End if

                        ##Perform regridding and interpolation of variable:
                        #rgdata_interp = _regrid_and_interpolate_levs(mclim_ds, var,
                        #                                         regrid_dataset=tclim_ds,
                        #                                         **regrid_kwargs)
                    #if comp == "lnd":
                    #    #Perform regridding of variable:
                    #    rgdata_interp = _regrid(mclim_ds, var,
                    #                        regrid_dataset=tclim_ds,
                    #                        **regrid_kwargs)
                    #if unstruct_cases[case_idx]:
                    if ('lat' not in mclim_ds.dims) and ('lat' not in mclim_ds.dims):
                        if ('ncol' in mclim_ds.dims) or ('lndgrid' in mclim_ds.dims):
                            #mclim_ds
                            print(f"Looks like test case '{case_name}' is unstructured, eh?")
                            #Check if any a FV file exists if using native grid
                            #case_latlon_file   = adf.latlon_files("test_latlon_file")
                            case_latlon_file = case_latlon_files[case_idx]
                            print("case_latlon_file",case_latlon_file,"\n")
                            #Check if any a weights file exists if using native grid, OPTIONAL
                            case_wgts_file   = case_wgts_files[case_idx]
                            case_method = case_methods[case_idx]
                            if case_wgts_file:
                                native_regrid_kwargs["wgt_file"] = case_wgts_file
                            else:
                                print("This looks like an unstructured case, but missing weights file")
                                case_wgts_file = None
                            if case_latlon_file:
                                native_regrid_kwargs["latlon_file"] = case_latlon_file
                            else:
                                print("This looks like an unstructured case, but missing lat/lon file")
                                #adf error thingy
                            #rgdata_interp = _regrid(mclim_ds, var,
                            #                    comp=comp,
                            #                    method=case_method,
                            #                    **native_regrid_kwargs)

                            rgdata_interp = _regrid(mclim_ds, var, comp, case_method, 
                                    wgt_file=case_wgts_file, 
                                    latlon_file=case_latlon_file)
                            #case_latlon_file
                            #fv_ds = xr.open_dataset(case_latlon_file)
                            #rgdata_interp = regrid_unstructured_to_latlon(mclim_ds, fv_ds.lat, fv_ds.lon, fv_ds)
                        else:
                            print("Trying everything and nothing is working. I guess this really is a problem!")
                    else:
                        rgdata_interp = mclim_ds
                    #else:
                    rgdata_interp = _regrid_and_interpolate_levs(rgdata_interp, var,
                                                                 regrid_dataset=tclim_ds,
                                                                 **regrid_kwargs)
                    #Extract defaults for variable:
                    var_default_dict = var_defaults.get(var, {})
                    if comp == "atm":
                        mask_var = "OCNFRAC"
                        mask_name = "ocean"
                    if comp == "lnd":
                        mask_var = "LANDFRAC"
                        mask_name = "land"
                    if 'mask' in var_default_dict:
                        if var_default_dict['mask'].lower() == mask_name:
                            #Check if the ocean fraction has already been regridded
                            #and saved:
                            if frc_ds:
                                frac = frc_ds[mask_var]
                                # set the bounds of regridded ocnfrac to 0 to 1
                                frac = xr.where(frac>1,1,frac)
                                frac = xr.where(frac<0,0,frac)

                                # apply ocean fraction mask to variable
                                rgdata_interp[mask_var] = frac
                                var_tmp = rgdata_interp[var]
                                var_tmp = pf.mask_land_or_ocean(var_tmp,frac)
                                rgdata_interp[var] = var_tmp
                            else:
                                print(f"\t    WARNING: {mask_var} not found, unable to apply mask to '{var}'")
                            #End if
                        else:
                            #Currently only ocean or land masks are supported, so print warning here:
                            wmsg = f"\t    WARNING: Currently the only variable mask option is '{mask_name}',"
                            wmsg += f"not '{var_default_dict['mask'].lower()}'"
                            print(wmsg)
                        #End if
                    #End if

                    #If the variable is the mask fraction, then save the dataset for use later:
                    if var == mask_var:
                        frc_ds = rgdata_interp
                    #End if

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
                    interp_bl_file = trgclimo_loc / f'{target}_{var}_baseline.nc'

                    if not adf.compare_obs and not interp_bl_file.is_file():
                        if comp == "atm":
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
                            #tgdata_interp = _regrid_and_interpolate_levs(tclim_ds, var,
                            #                                         **regrid_kwargs)
                        #if comp == "lnd":
                        #    #Perform regridding of variable:
                        #    tgdata_interp = _regrid(mclim_ds, var,
                        #                        regrid_dataset=tclim_ds,
                        #                        **regrid_kwargs)
            
                        #if unstruct_base:
                        if ('lat' not in tclim_ds.dims) and ('lat' not in tclim_ds.dims):
                            if ('ncol' in tclim_ds.dims) or ('lndgrid' in tclim_ds.dims):
                                print(f"Looks like baseline case '{target}' is unstructured, eh?")
                                native_regrid_kwargs = {}
                                #Check if any a FV file exists if using native grid
                                baseline_latlon_file   = adf.latlon_files["base_latlon_file"]

                                #Check if any a weights file exists if using native grid, OPTIONAL
                                baseline_wgts_file   = adf.latlon_wgt_files["base_wgts_file"]
                                base_method = adf.latlon_regrid_method["base_regrid_method"]
                                if baseline_wgts_file:
                                    native_regrid_kwargs["wgt_file"] = baseline_wgts_file
                                else:
                                    print("This looks like an unstructured case, but missing weights file")
                                if baseline_latlon_file:
                                    native_regrid_kwargs["latlon_file"] = baseline_latlon_file
                                else:
                                    print("This looks like an unstructured case, but missing lat/lon file")
                                    #adf error thingy
                                tgdata_interp = _regrid(tclim_ds, var,
                                                comp=comp,
                                                method=base_method,
                                                **native_regrid_kwargs)
                            else:
                                print("Trying everything and nothing is working. I guess this really is a problem!")
                        else:
                            tgdata_interp = tclim_ds

                        tgdata_interp = _regrid_and_interpolate_levs(tgdata_interp, var,
                                                                    regrid_dataset=tclim_ds,
                                                                    **regrid_kwargs)
                        if tgdata_interp is None:
                            #Something went wrong during interpolation, so just cycle through
                            #for now:
                            continue
                        #End if

                        if comp == "atm":
                            mask_var = "OCNFRAC"
                            mask_name = "ocean"
                        if comp == "lnd":
                            mask_var = "LANDFRAC"
                            mask_name = "land"
                        #If the variable is ocean fraction, then save the dataset for use later:
                        if var == mask_var:
                            frc_ds = tgdata_interp
                        #End if
                        if 'mask' in var_default_dict:
                            if var_default_dict['mask'].lower() == mask_name:
                                #Check if the ocean fraction has already been regridded
                                #and saved:
                                if frc_ds:
                                    frac = frc_ds[mask_var]
                                    # set the bounds of regridded ocnfrac to 0 to 1
                                    frac = xr.where(frac>1,1,frac)
                                    frac = xr.where(frac<0,0,frac)

                                    # apply ocean fraction mask to variable
                                    rgdata_interp[mask_var] = frac
                                    var_tmp = rgdata_interp[var]
                                    var_tmp = pf.mask_land_or_ocean(var_tmp,frac)
                                    rgdata_interp[var] = var_tmp
                                else:
                                    print(f"\t    WARNING: {mask_var} not found, unable to apply mask to '{var}'")
                                #End if
                            else:
                                #Currently only ocean or land masks are supported, so print warning here:
                                wmsg = f"\t    WARNING: Currently the only variable mask option is '{mask_name}',"
                                wmsg += f"not '{var_default_dict['mask'].lower()}'"
                                print(wmsg)
                            #End if
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
                    print("\t    INFO: Regridded file already exists, so skipping...")
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
    from regrid_se_to_fv import make_se_regridder, regrid_se_data_conservative

    #Extract keyword arguments:
    if 'ps_file' in kwargs:
        ps_file = kwargs['ps_file']
    else:
        ps_file = None
    #End if
    if 'pmid_file' in kwargs:
        pmid_file = kwargs['pmid_file']
    else:
        pmid_file = None
    #End if

    #Extract variable info from model data (and remove any degenerate dimensions):
    #mdata = model_dataset[var_name].squeeze()
    mdata = model_dataset[var_name].squeeze()
    mdat_ofrac = None
    #if regrid_ofrac:
    #    if 'OCNFRAC' in model_dataset:
    #        mdat_ofrac = model_dataset['OCNFRAC'].squeeze()

    #Check if variable has a vertical component:
    if 'lev' in mdata.dims or 'ilev' in mdata.dims:
        has_lev = True

        #If lev exists, then determine what kind of vertical coordinate
        #is being used:
        if 'lev' in mdata.dims:
            lev_attrs = model_dataset['lev'].attrs
        elif 'ilev' in mdata.dims:
            lev_attrs = model_dataset['ilev'].attrs

        #First check if there is a "vert_coord" attribute:
        if 'vert_coord' in lev_attrs:
            vert_coord_type = lev_attrs['vert_coord']
        else:
            #Next check that the "long_name" attribute exists:
            if 'long_name' in lev_attrs:
                #Extract long name:
                lev_long_name = lev_attrs['long_name']

                #Check for "keywords" in the long name:
                if 'hybrid level' in lev_long_name:
                    #Set model to hybrid vertical levels:
                    vert_coord_type = "hybrid"
                elif 'zeta level' in lev_long_name:
                    #Set model to height (z) vertical levels:
                    vert_coord_type = "height"
                else:
                    #Print a warning, and skip variable re-gridding/interpolation:
                    wmsg = "WARNING! Unable to determine the vertical coordinate"
                    wmsg +=f" type from the 'lev' long name, which is:\n'{lev_long_name}'"
                    print(wmsg)
                    return None
                #End if

            else:
                #Print a warning, and assume hybrid levels (for now):
                wmsg = "WARNING!  No long name found for the 'lev' dimension,"
                wmsg += f" so no re-gridding/interpolation will be done."
                print(wmsg)
                return None
            #End if
        #End if

    else:
        has_lev = False
    #End if

    #Check if variable has a vertical levels dimension:
    if has_lev:

        if vert_coord_type == "hybrid":
            # Need hyam, hybm, and P0 for vertical interpolation of hybrid levels:
            if 'lev' in mdata.dims:
                if ('hyam' not in model_dataset) or ('hybm' not in model_dataset):
                    print(f"!! PROBLEM -- NO hyam or hybm for 3-D variable {var_name}, so it will not be re-gridded.")
                    return None #Return None to skip to next variable.
                #End if
                mhya = model_dataset['hyam']
                mhyb = model_dataset['hybm']
            elif 'ilev' in mdata.dims:
                if ('hyai' not in model_dataset) or ('hybi' not in model_dataset):
                    print(f"!! PROBLEM -- NO hyai or hybi for 3-D variable {var_name}, so it will not be re-gridded.")
                    return None #Return None to skip to next variable.
                #End if
                mhya = model_dataset['hyai']
                mhyb = model_dataset['hybi']
            if 'time' in mhya.dims:
                mhya = mhya.isel(time=0).squeeze()
            if 'time' in mhyb.dims:
                mhyb = mhyb.isel(time=0).squeeze()
            if 'P0' in model_dataset:
                P0_tmp = model_dataset['P0']
                if isinstance(P0_tmp, xr.DataArray):
                    #All of these value should be the same,
                    #so just grab the first one:
                    P0 = P0_tmp[0]
                else:
                    #Just rename variable:
                    P0 = P0_tmp
                #End if
            else:
                P0 = 100000.0  # Pa
            #End if

        elif vert_coord_type == "height":
            #Initialize already-regridded PMID logical:
            regridded_pmid = False

            #Need mid-level pressure for vertical interpolation of height levels:
            if 'PMID' in model_dataset:
                mpmid = model_dataset['PMID']
            else:
                #Check if target has an associated surface pressure field:
                if pmid_file:
                    mpmid_ds = xr.open_dataset(pmid_file)
                    mpmid = mpmid_ds['PMID']
                    #This mid-level pressure field has already been regridded:
                    regridded_pmid = True
                else:
                    print(f"!! PROBLEM -- NO PMID for 3-D variable {var_name}, so it will not be re-gridded.")
                    return None
                #End if
            #End if
        #End if (vert_coord_type)

        #It is probably good to try and acquire PS for all vertical coordinate types, so try here:
        regridded_ps = False
        if 'PS' in model_dataset:
            mps = model_dataset['PS']
        else:
            #Check if target has an associated surface pressure field:
            if ps_file:
                mps_ds = xr.open_dataset(ps_file)
                mps = mps_ds['PS']
                #This surface pressure field has already been regridded:
                regridded_ps = True
            else:
                print(f"!! PROBLEM -- NO PS for 3-D variable {var_name}, so it will not be re-gridded.")
                return None
            #End if
        #End if
    #End if (has_lev)

    #Regrid variable to target dataset (if available):
    if regrid_dataset:

        #Extract grid info from target data:
        if 'time' in regrid_dataset.coords:
            if 'lev' in regrid_dataset.coords:
                tgrid = regrid_dataset.isel(time=0, lev=0).squeeze()
            else:
                tgrid = regrid_dataset.isel(time=0).squeeze()
            #End if
        #End if

        #Regrid model data to match target grid:
        rgdata = regrid_data(mdata, tgrid, method=1)
        if mdat_ofrac:
            rgofrac = regrid_data(mdat_ofrac, tgrid, method=1)
        #Regrid surface pressure if need be:
        if has_lev:
            if not regridded_ps:
                rg_ps = regrid_data(mps, tgrid, method=1)
            else:
                rg_ps = mps
            #End if

            #Also regrid mid-level pressure if need be:
            if vert_coord_type == "height":
                if not regridded_pmid:
                    rg_pmid = regrid_data(mpmid, tgrid, method=1)
                else:
                    rg_pmid = mpmid
                #End if
            #End if
        #End if
    else:
        #Just rename variables:
        rgdata = mdata
        if has_lev:
            rg_ps = mps
            if vert_coord_type == "height":
                rg_pmid = mpmid
            #End if
        #End if
    #End if

    #Vertical interpolation:

    #Interpolate variable to default pressure levels:
    if has_lev:

        if vert_coord_type == "hybrid":
            #Interpolate from hybrid sigma-pressure to the standard pressure levels:
            rgdata_interp = pf.lev_to_plev(rgdata, rg_ps, mhya, mhyb, P0=P0, \
                                           convert_to_mb=True)
        elif vert_coord_type == "height":
            #Interpolate variable using mid-level pressure (PMID):
            rgdata_interp = pf.pmid_to_plev(rgdata, rg_pmid, convert_to_mb=True)
        else:
            #The vertical coordinate type is un-recognized, so print warning and
            #skip vertical interpolation:
            wmsg = f"WARNING! Un-recognized vertical coordinate type: '{vert_coord_type}',"
            wmsg += f" for variable '{var_name}'.  Skipping vertical interpolation."
            print(wmsg)
            #Don't process variable:
            return None
        #End if
    else:
        #Just rename variable:
        rgdata_interp = rgdata
    #End if

    #Convert to xarray dataset:
    rgdata_interp = rgdata_interp.to_dataset()
    if mdat_ofrac:
        rgdata_interp['OCNFRAC'] = rgofrac

    #Add surface pressure to variable if a hybrid (just in case):
    if has_lev:
        rgdata_interp['PS'] = rg_ps

        #Update "vert_coord" attribute for variable "lev":
        rgdata_interp['lev'].attrs.update({"vert_coord": "pressure"})
    #End if

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







'''def _regrid(model_dataset, var_name, comp, method, **kwargs):
    """
    Function that takes a variable from a model xarray
    dataset, regrids it to another dataset's lat/lon
    coordinates (if applicable)
    ----------
    model_dataset -> The xarray dataset which contains the model variable data
    var_name      -> The name of the variable to be regridded/interpolated.

    Optional inputs:

    ps_file        -> NOT APPLICABLE: A NetCDF file containing already re-gridded surface pressure
    regrid_dataset -> The xarray dataset that contains the lat/lon grid that
                      "var_name" will be regridded to.  If not present then
                      only the vertical interpolation will be done.

    kwargs         -> Keyword arguments that contain paths to THE REST IS NOT APPLICABLE: surface pressure
                      and mid-level pressure files, which are necessary for
                      certain types of vertical interpolation.

    This function returns a new xarray dataset that contains the regridded
    model variable.
    """

    #Import ADF-specific functions:
    import numpy as np
    import plotting_functions as pf

    #comp = adf.model_component
    if comp == "atm":
        comp_grid = "ncol"
    if comp == "lnd":
        comp_grid = "lndgrid"

    #Extract variable info from model data (and remove any degenerate dimensions):
    if var_name:
        mdata = model_dataset[var_name].squeeze()
    else:
        #if isinstance(xo, xr.Dataset):
        mdata = model_dataset
    mdat_ofrac = None
    #if regrid_lfrac:
    #    if 'LANDFRAC' in model_dataset:
    #        mdat_lfrac = model_dataset['LANDFRAC'].squeeze()

    #Regrid variable to target dataset (if available):
    #if regrid_dataset:
    if 1==1:

        # Hardwiring for now
        #con_weight_file = "/glade/work/wwieder/map_ne30pg3_to_fv0.9x1.25_scripgrids_conserve_nomask_c250108.nc"
        if "wgt_file" in kwargs:
            weight_file = kwargs["wgt_file"]
        #latlon_file = '/glade/derecho/scratch/wwieder/ctsm5.3.018_SP_f09_t232_mask/run/ctsm5.3.018_SP_f09_t232_mask.clm2.h0.0001-01.nc'
        if "latlon_file" in kwargs:
            latlon_file = kwargs["latlon_file"]
        else:
            print("Well, it looks like you're missing a target grid file for regridding!")
            #adferror thing

        print("\nlatlon_file",latlon_file,"\n")
        fv_ds = xr.open_dataset(latlon_file)

        model_dataset[var_name] = model_dataset[var_name].fillna(0)
        if comp == "lnd":
            model_dataset['landfrac']= model_dataset['landfrac'].fillna(0)
            if var_name:
                model_dataset[var_name] = model_dataset[var_name] * model_dataset.landfrac  # weight flux by land frac
            s_data = model_dataset.landmask.isel(time=0)
            d_data = fv_ds.landmask
        else:
            if var_name:
                s_data = model_dataset[var_name].isel(time=0)
                d_data = fv_ds[var_name]
            else:
                s_data = model_dataset.isel(time=0)
                d_data = fv_ds

        #Regrid model data to match target grid:
        # These two functions come with import regrid_se_to_fv
        regridder = make_se_regridder(weight_file=weight_file,
                                      s_data = s_data, #model_dataset.landmask.isel(time=0),
                                      d_data = d_data, #fv_ds.landmask,
                                      Method = method,  # Bug in xesmf needs this without "n"
                                      )
        if method == 'coservative':
            rgdata = regrid_se_data_conservative(regridder, model_dataset, comp_grid)

        if comp == "lnd":
            rgdata[var_name] = (rgdata[var_name] / rgdata.landfrac)
            rgdata['landmask'] = fv_ds.landmask
            rgdata['landfrac'] = rgdata.landfrac.isel(time=0)

        rgdata['lat'] = fv_ds.lat
        #rgdata['landmask'] = fv_t232.landmask
        #rgdata['landfrac'] = rgdata.landfrac.isel(time=0)

        # calculate area
        area_km2 = np.zeros(shape=(len(rgdata['lat']), len(rgdata['lon'])))
        earth_radius_km = 6.37122e3  # in meters

        yres_degN = np.abs(np.diff(rgdata['lat'].data))  # distances between gridcell centers...
        xres_degE = np.abs(np.diff(rgdata['lon']))  # ...end up with one less element, so...
        yres_degN = np.append(yres_degN, yres_degN[-1])  # shift left (edges <-- centers); assume...
        xres_degE = np.append(xres_degE, xres_degE[-1])  # ...last 2 distances bet. edges are equal

        dy_km = yres_degN * earth_radius_km * np.pi / 180  # distance in m
        phi_rad = rgdata['lat'].data * np.pi / 180  # degrees to radians

        # grid cell area
        for j in range(len(rgdata['lat'])):
            for i in range(len(rgdata['lon'])):
              dx_km = xres_degE[i] * np.cos(phi_rad[j]) * earth_radius_km * np.pi / 180  # distance in m
              area_km2[j,i] = dy_km[j] * dx_km

        rgdata['area'] = xr.DataArray(area_km2,
                                      coords={'lat': rgdata.lat, 'lon': rgdata.lon},
                                      dims=["lat", "lon"])
        rgdata['area'].attrs['units'] = 'km2'
        rgdata['area'].attrs['long_name'] = 'Grid cell area'
    else:
        #Just rename variables:
        rgdata = mdata
    #End if

    #Return dataset:
    return rgdata


import xarray as xr
import xesmf
import numpy as np

def make_se_regridder(weight_file, s_data, d_data,
                      Method='coservative'
                      ):
    # Intialize dict for xesmf.Regridder
    regridder_kwargs = {}

    if weight_file:
        weights = xr.open_dataset(weight_file)
        regridder_kwargs['weights'] = weights
    else:
        print("No weights file given, so I'm gonna need to make one. Please have a seat and the next associate will be with you shortly. Please don't tap the glass!")
        regridder_kwargs['method'] = 'coservative'
    
    in_shape = weights.src_grid_dims.load().data

    # Since xESMF expects 2D vars, we'll insert a dummy dimension of size-1
    if len(in_shape) == 1:
        in_shape = [1, in_shape.item()]

    # output variable shape
    out_shape = weights.dst_grid_dims.load().data.tolist()[::-1]

    dummy_in = xr.Dataset(
        {
            "lat": ("lat", np.empty((in_shape[0],))),
            "lon": ("lon", np.empty((in_shape[1],))),
        }
    )
    dummy_out = xr.Dataset(
        {
            "lat": ("lat", weights.yc_b.data.reshape(out_shape)[:, 0]),
            "lon": ("lon", weights.xc_b.data.reshape(out_shape)[0, :]),
        }
    )
    # Hard code masks for now, not sure this does anything?
    s_mask = xr.DataArray(s_data.data.reshape(in_shape[0],in_shape[1]), dims=("lat", "lon"))
    dummy_in['mask']= s_mask
    
    d_mask = xr.DataArray(d_data.values, dims=("lat", "lon"))  
    dummy_out['mask']= d_mask                

    # do source and destination grids need masks here?
    # See xesmf docs https://xesmf.readthedocs.io/en/stable/notebooks/Masking.html#Regridding-with-a-mask
    regridder = xesmf.Regridder(
        dummy_in,
        dummy_out,
        weights=weight_file,
        # results seem insensitive to this method choice
        # choices are coservative_normed, coservative, and bilinear
        method=Method,
        reuse_weights=True,
        periodic=True,
        **regridder_kwargs
    )
    return regridder






'''
def regrid_se_data_bilinear(regridder, data_to_regrid, comp_grid):
    updated = data_to_regrid.copy().transpose(..., comp_grid).expand_dims("dummy", axis=-2)
    regridded = regridder(updated.rename({"dummy": "lat", comp_grid: "lon"}),
                         skipna=True, na_thres=1,
                         )
    return regridded

def regrid_se_data_conservative(regridder, data_to_regrid, comp_grid):
    updated = data_to_regrid.copy().transpose(..., comp_grid).expand_dims("dummy", axis=-2)
    regridded = regridder(updated.rename({"dummy": "lat", comp_grid: "lon"}) )
    return regridded






import xarray as xr
import numpy as np
import xesmf

def _regrid(model_dataset, var_name, comp, method, **kwargs):
    """
    Regrid model data (with optional vertical levels) to a lat/lon grid.
    """
    if comp == "atm":
        comp_grid = "ncol"
    if comp == "lnd":
        comp_grid = "lndgrid"

    # Extract and squeeze variable
    mdata = model_dataset[var_name].squeeze()

    if "wgt_file" in kwargs:
        weight_file = kwargs["wgt_file"]
    else:
        raise ValueError("Missing regridder weight file (wgt_file).")

    if "latlon_file" in kwargs:
        latlon_file = kwargs["latlon_file"]
    else:
        raise ValueError("Missing target lat/lon grid file (latlon_file).")

    # Load target grid (lat/lon) from the provided dataset
    fv_ds = xr.open_dataset(latlon_file)

    model_dataset[var_name] = model_dataset[var_name].fillna(0)

    # Identify source and destination data for regridding
    if comp == "lnd":
        model_dataset['landfrac']= model_dataset['landfrac'].fillna(0)
        model_dataset[var_name] = model_dataset[var_name] * model_dataset.landfrac  # weight flux by land frac
        s_data = model_dataset.landmask.isel(time=0)
        d_data = fv_ds.landmask
    else:
        s_data = mdata.isel(time=0)
        d_data = fv_ds[var_name] if var_name in fv_ds else fv_ds

    # Create regridder
    regridder = make_se_regridder(weight_file, s_data, d_data, method)

    # Handle 2D vs 3D data (with or without 'lev')
    if "lev" in mdata.dims:
        # Iterate through each level, regrid separately
        regridded_data = []
        for lev in mdata.lev.values:
            lev_slice = mdata.sel(lev=lev)
            rgdata = regrid_se_data_conservative(regridder, lev_slice, comp_grid)
            rgdata = rgdata.expand_dims("lev")
            rgdata["lev"] = lev
            regridded_data.append(rgdata)

        # Combine all levels
        rgdata = xr.concat(regridded_data, dim="lev")
    else:
        # 2D regridding (no vertical levels)
        rgdata = regrid_se_data_conservative(regridder, mdata, comp_grid)

    # Ensure output matches the target grid dimensions
    rgdata["lat"] = fv_ds.lat
    rgdata["lon"] = fv_ds.lon

    # Compute grid cell area (optional but useful for post-processing)
    rgdata["area"] = _calculate_area(rgdata)

    return rgdata

def make_se_regridder(weight_file, s_data, d_data, Method='conservative'):
    """
    Create xESMF regridder for spectral element grids.
    """
    regridder_kwargs = {}

    # Load weights if available
    weights = xr.open_dataset(weight_file)
    regridder_kwargs['weights'] = weights

    in_shape = weights.src_grid_dims.load().data

    # Ensure 2D compatibility for xESMF (reshape if needed)
    if len(in_shape) == 1:
        in_shape = [1, in_shape.item()]

    out_shape = weights.dst_grid_dims.load().data.tolist()[::-1]

    dummy_in = xr.Dataset({
        "lat": ("lat", np.empty((in_shape[0],))),
        "lon": ("lon", np.empty((in_shape[1],))),
    })
    dummy_out = xr.Dataset({
        "lat": ("lat", weights.yc_b.data.reshape(out_shape)[:, 0]),
        "lon": ("lon", weights.xc_b.data.reshape(out_shape)[0, :]),
    })

    # Handle source and destination masks
    s_mask = xr.DataArray(s_data.data.reshape(in_shape[0], in_shape[1]), dims=("lat", "lon"))
    dummy_in['mask'] = s_mask

    d_mask = xr.DataArray(d_data.values, dims=("lat", "lon"))
    dummy_out['mask'] = d_mask

    # Create xESMF regridder
    regridder = xesmf.Regridder(
        dummy_in,
        dummy_out,
        weights=weight_file,
        method=Method,
        reuse_weights=True,
        periodic=True,
        #**regridder_kwargs
    )

    return regridder

'''def regrid_se_data_conservative(regridder, model_dataset, comp_grid):
    """
    Apply conservative regridding to the dataset.
    """
    return regridder(model_dataset)
'''

def _calculate_area(rgdata):
    """
    Compute grid cell area for regridded dataset.
    """
    area_km2 = np.zeros((len(rgdata['lat']), len(rgdata['lon'])))
    earth_radius_km = 6.37122e3

    yres_degN = np.abs(np.diff(rgdata['lat'].data))
    xres_degE = np.abs(np.diff(rgdata['lon']))

    yres_degN = np.append(yres_degN, yres_degN[-1])
    xres_degE = np.append(xres_degE, xres_degE[-1])

    dy_km = yres_degN * earth_radius_km * np.pi / 180
    phi_rad = rgdata['lat'].data * np.pi / 180

    for j in range(len(rgdata['lat'])):
        for i in range(len(rgdata['lon'])):
            dx_km = xres_degE[i] * np.cos(phi_rad[j]) * earth_radius_km * np.pi / 180
            area_km2[j, i] = dy_km[j] * dx_km

    return xr.DataArray(area_km2, coords={'lat': rgdata.lat, 'lon': rgdata.lon}, dims=["lat", "lon"], attrs={'units': 'km2', 'long_name': 'Grid cell area'})





