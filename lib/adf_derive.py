import glob
import os
from pathlib import Path
import xarray as xr


def check_derive(self, res, var, case_name, diag_var_list, constit_dict, hist_file_ds, hist0):
    """
    For incoming variable, look for list of constituents if available
     - as a list in variable defaults file

     If the variable does not have the argument `derivable_from` or `derivable_from_cam_chem`,
     then it will be assumed not to be a derivable variable, just missing from history file

     If the variable does have the argument `derivable_from` or `derivable_from_cam_chem`,
     first check cam-chem, then regular cam.

    Arguments
    ---------
        self: self object?
            - 
        res: dict
            - 
        var: str
            - 
        case_name: str
            - 
        diag_var_list: list
            - 
        constit_dict: dict
            - 
        hist_file_ds: xarray.DataSet
            - 
        hist0: str
            - 
    
    Returns
    -------
        constit_list: list
           - list of declared consituents from the variable defaults yaml file
           - empty list:
             * if missing `derived_from` argument(s)
             * if `derived_from` argument(s) exist but not declared
        
        diag_var_list: list
           - updated list (if applicable) of ADF variables for time series creation
    """

    # Aerosol Calcs
    #--------------

    # Always make sure PMID is made if aerosols are desired in config file
    # Since there's no requirement for `aerosol_zonal_list`, allow it to be absent:
    azl = res.get("aerosol_zonal_list", [])
    if azl:
        if "PMID" not in diag_var_list:
            if any(item in azl for item in diag_var_list):
                diag_var_list += ["PMID"]
        if "T" not in diag_var_list:
            if any(item in azl for item in diag_var_list):
                diag_var_list += ["T"]
    #End aerosol calcs

    # Set error messages for printing/debugging
    # Derived variable, but missing constituent list
    constit_errmsg = f"create time series for {case_name}:"
    constit_errmsg += f"\n Can't create time series for {var}. \n\tThis variable"
    constit_errmsg += " is flagged for derivation, but is missing list of constiuents."
    constit_errmsg += "\n\tPlease add list of constituents to 'derivable_from' "
    constit_errmsg += f"for {var} in variable defaults yaml file."

    # No time series creation
    exit_msg = f"WARNING: {var} is not in the file {hist0} and can't be derived."
    exit_msg += "\n\t  ** No time series will be generated. **\n"

    # Initialiaze list for constituents
    # NOTE: This is if the variable is NOT derivable but needs
    #       an empty list as a check later
    constit_list = []

    try_cam_constits = True
    try:
        vres = res[var]
    except KeyError:
        #if verbose: # make this a wrapper!
        #    print(exit_msg)
        print(exit_msg)
        self.debug_log(exit_msg)
        return diag_var_list, constit_dict

    # Check first if variable is potentially part of a CAM-CHEM run
    if "derivable_from_cam_chem" in vres:
        constit_list = vres["derivable_from_cam_chem"]

        if constit_list:
            if all(item in hist_file_ds.data_vars for item in constit_list):
                # Set check to look for regular CAM constituents in variable defaults
                try_cam_constits = False
                msg = f"derive time series for {case_name}:"
                msg += "\n\tLooks like this a CAM-CHEM run, "
                msg += f"checking constituents for '{var}'"
                self.debug_log(msg)
        else:
            self.debug_log(constit_errmsg)
        # End if
    # End if
    
    # If not CAM-CHEM, check regular CAM runs
    if try_cam_constits:
        if "derivable_from" in vres:
            constit_list = vres["derivable_from"]
        else:
            # Missing variable or missing derivable_from argument
            der_from_msg = f"derive time series for {case_name}:"
            der_from_msg += f"\n Can't create time series for {var}.\n\tEither "
            der_from_msg += "the variable is missing from CAM output or it is a "
            der_from_msg += "derived quantity and is missing the 'derivable_from' "
            der_from_msg += "config argument.\n\tPlease add variable to CAM run "
            der_from_msg += "or set appropriate argument in variable "
            der_from_msg += "defaults yaml file."
            self.debug_log(der_from_msg)
        # End if
    # End if

    # Log if this variable can be derived but is missing list of constituents
    if isinstance(constit_list, list) and not constit_list:
        self.debug_log(constit_errmsg)

    # Check if any constituents were found
    if constit_list:
        # Add variable and constituent list to dictionary
        constit_dict[var] = constit_list

        # Aadd constituents to ADF diag variable list for time series generation
        for constit in constit_list:
            if constit not in diag_var_list:
                diag_var_list.append(constit)
    else:
        #if verbose: # make this a wrapper!
        #    print(exit_msg)
        print(exit_msg)
        self.debug_log(exit_msg)
    # End if

    return diag_var_list, constit_dict

########

def derive_variable(self, case_name, res=None, vars_to_derive=None, ts_dir=None,
                         constit_list=None, overwrite=None):
    """
    Derive variables acccording to steps given here.  Since derivations will depend on the
    variable, each variable to derive will need its own set of steps below.

    Caution: this method assumes that there will be one time series file per variable

    If the file for the derived variable exists, the kwarg `overwrite` determines
    whether to overwrite the file (true) or exit with a warning message.

    """

    # Loop through derived variables
    for var in vars_to_derive:
        print(f"\t - deriving time series for {var}")

        #Check whether there are parts to derive from and if there is an associated equation
        vres = res.get(var, {})
        #Check if appropriate config variable
        #if "derive" in vres:
        #if "from" in vres["derive"]:
                    
        #First check if it is just a simple derivation from constituents
        constit_list = vres["derivable_from"]
        flag = "derivable_from"

        # Grab all required time series files for derived variable
        constit_files = []
        constit_files_dict = {}

        for constit in constit_list:
            if len(sorted(glob.glob(os.path.join(ts_dir, f"*.{constit}.*.nc")))) > 1:
                mutli_ts = True
                constit_files_dict[constit] = sorted(glob.glob(os.path.join(ts_dir, f"*.{constit}.*.nc")))
            else:
                mutli_ts = False
                if glob.glob(os.path.join(ts_dir, f"*.{constit}.*.nc")):
                    print("single time series HERE??")
                    constit_files.append(glob.glob(os.path.join(ts_dir, f"*.{constit}.*"))[0])

        if mutli_ts:
            for i in range(len(constit_files_dict[constit_list[0]])):

                ahh = []
                for cons in constit_files_dict.keys():
                    ahh.append(constit_files_dict[cons][i])

                #Check if all the constituent files were found
                if len(ahh) != len(constit_list):
                    ermsg = f"Not all constituent files present; {var} cannot be calculated."
                    ermsg += f" Please remove {var} from diag_var_list or find the relevant CAM files."
                    print(ermsg)
                    continue
                #Open a new dataset with all the constituent files/variables
                ds = xr.open_mfdataset(ahh, compat='override')
                # create new file name for derived variable

                #print("fdghjk",constit_files_dict[constit_list[0]][i])

                derived_file = constit_files_dict[constit_list[0]][i].replace(list(constit_files_dict.keys())[0], var)
                #print("derived_file",derived_file,"\n")
                #Check if clobber is true for file
                if Path(derived_file).is_file():
                    if overwrite:
                        Path(derived_file).unlink()
                    else:
                        print(
                            f"[{__name__}] Warning: '{var}' file was found and overwrite is False. Will use existing file."
                        )
                        continue

                ###



        else:
            #Check if all the constituent files were found
            print(f"{var}: matchies??",len(constit_files) == len(constit_list))
            if len(constit_files) != len(constit_list):
                ermsg = f"Not all constituent files present; {var} cannot be calculated."
                ermsg += f" Please remove {var} from diag_var_list or find the relevant CAM files."
                print(ermsg)
                continue
            #Open a new dataset with all the constituent files/variables
            ds = xr.open_mfdataset(constit_files, compat='override')
        
            # create new file name for derived variable
            derived_file = constit_files[0].replace(constit_list[0], var)


                
            #print("derived_file",derived_file,"\n")
            #Check if clobber is true for file
            if Path(derived_file).is_file():
                if overwrite:
                    Path(derived_file).unlink()
                else:
                    print(
                        f"[{__name__}] Warning: '{var}' file was found and overwrite is False. Will use existing file."
                    )
                    continue

            """if flag == "derive_interp":
                derive_interp(ds, vres, constit_list, var, derived_file)
                    
            if flag == "derive_mask":
                derive_masked(ds, vres, constit_list, var, ts_dir, i, derived_file)
      
            if flag == "derivable_from":
                derive_from_constits(ds, constit_list, var, derived_file)"""

        if flag == "derivable_from":
            derive_from_constits(ds, constit_list, var, derived_file)


        
########

#Helper Function(s)
def derive_from_constits(ds, ts_dir, res, constit_list, var, derived_file):
    #NOTE: this will need to be changed when derived equations are more complex! - JR
    if var == "RESTOM":
        print("RESTOM\n")
        der_val = ds["FSNT"]-ds["FLNT"]
    else:
        #Loop through all constituents and sum
        der_val = 0
        for v in constit_list:
            der_val += ds[v]

    #Set derived variable name and add to dataset
    der_val.name = var
    ds[var] = der_val

    #If variable is an aerosol, do extra processesing
    if var in res["aerosol_zonal_list"]:
        ds = calc_aerosol(var, ts_dir, res, ds)

    #Drop all constituents from final saved dataset
    #These are not necessary because they have their own time series files
    ds_final = ds.drop_vars(constit_list)
    ds_final.to_netcdf(derived_file, unlimited_dims='time', mode='w')

def derive_interp(ds, vres, constit_list, var, derived_file):
    for dim in ["time","lat","lon","lev","ilev"]:
        if dim in vres["derive"].keys():
            #ts_exist = glob.glob(os.path.join(ts_case_dir, f"*.{der_from}.*"))
            #der_from_ds = xr.open_dataset(ts_exist[0])
            der_from_ds = ds
            #Grab variable to derive from
            #constit_list
            der_from_var = der_from_ds[constit_list[0]]
            #der_from_var = der_from_ds[der_from]

            # Interpolate the data to the nearest requested value: vres["derive"][dim]
            der_var = der_from_var.interp({dim: vres["derive"][dim]}, method='nearest')
                                                                
            #Set derived variable in dataset and remove the original variable
            der_from_ds[var] = der_var
            ds_final = der_from_ds.drop_vars(constit_list)
            ds_final.to_netcdf(derived_file, unlimited_dims='time', mode='w')

def derive_masked(ds, vres, constit_list, var, ts_dir, index, derived_file):
    der_from_ds = ds
    der_from_var = der_from_ds[constit_list[0]]
    #Derive variables that come from other means
    #EXAMPLE: derive SST's from TS if not in CAM output
    #if 'SST' in diag_var_list and not glob.glob(os.path.join(ts_dir, f"*SST*")):
    #if var in diag_var_list and not glob.glob(os.path.join(ts_dir, f"*{var}*")):

    if 'mask' in vres:
        if vres['mask'].lower() == 'ocean':
            #Check if the ocean fraction has already been regridded
            #and saved:
            #if ts_ds:
            if ds:
                #ofrac_ds = xr.open_dataset(glob.glob(os.path.join(ts_dir, f"*OCNFRAC*"))[0])
                ocnfrac_file = sorted(glob.glob(os.path.join(ts_dir, "*OCNFRAC*")))
                ofrac_ds = xr.open_mfdataset(ocnfrac_file[index], compat='override')
                if ofrac_ds:
                    ofrac = ofrac_ds['OCNFRAC']
                    # set the bounds of regridded ocnfrac to 0 to 1
                    ofrac = xr.where(ofrac>1,1,ofrac)
                    ofrac = xr.where(ofrac<0,0,ofrac)
                    # mask the land in TS for global means
                    #ts_ds['OCNFRAC'] = ofrac
                    ds['OCNFRAC'] = ofrac
                    #der_from_var = der_from_ds[constit_list[0]]
                    #ts_tmp = ts_ds[constit_list[0]]
                    ts_tmp = ds[constit_list[0]]
                    #Import ADF-specific modules:
                    import plotting_functions as pf
                    ts_tmp = pf.mask_land_or_ocean(ts_tmp,ofrac)
                    #ts_ds['SST'] = ts_tmp
                    #ts_ds[var] = ts_tmp
                    ds[var] = ts_tmp
                    #Set derived variable in dataset and remove the original variable
                    der_from_ds[var] = der_from_var
                    ds_final = der_from_ds.drop_vars(constit_list)
                    ds_final.to_netcdf(derived_file, unlimited_dims='time', mode='w')
                else:
                    wmsg = "OCNFRAC not found in CAM output,"
                    wmsg += f" unable to apply mask to '{var}'"
                    print(wmsg)
            else:
                wmsg = f"{der_from_var} not found in CAM output,"
                wmsg += f" unable to apply mask to '{var}'"
                print(wmsg)
            #End if
        #End if
    #End if



def calc_aerosol(var, ts_dir, res, ds):
    #These will be multiplied by rho (density of dry air)
    ds_pmid_done = False
    ds_t_done = False

    #Only calculate once for all aerosol vars
    if not ds_pmid_done:
        ds_pmid = _load_dataset(glob.glob(os.path.join(ts_dir, "*.PMID.*"))[0])
        ds_pmid_done = True
        if not ds_pmid:
            errmsg = "Missing necessary files for dry air density (rho) "
            errmsg += "calculation.\nPlease make sure 'PMID' is in the CAM "
            errmsg += "run for aerosol calculations"
            print(errmsg)
            #continue
    if not ds_t_done:
        ds_t = _load_dataset(glob.glob(os.path.join(ts_dir, "*.T.*"))[0])
        ds_t_done = True
        if not ds_t:
            errmsg = "Missing necessary files for dry air density (rho) "
            errmsg += "calculation.\nPlease make sure 'T' is in the CAM "
            errmsg += "run for aerosol calculations"
            print(errmsg)
            #continue

    #Multiply aerosol by dry air density (rho): (P/Rd*T)
    ds[var] = ds[var]*(ds_pmid["PMID"]/(res["Rgas"]*ds_t["T"]))

    #Sulfate conversion factor
    if var == "SO4":
        ds[var] = ds[var]*(96./115.)
    return ds

def _load_dataset(fils):
    """
    This method exists to get an xarray Dataset from input file information that
    can be passed into the plotting methods.

    Parameters
    ----------
    fils : list
        strings or paths to input file(s)

    Returns
    -------
    xr.Dataset

    Notes
    -----
    When just one entry is provided, use `open_dataset`, otherwise `open_mfdatset`
    """
    import warnings  # use to warn user about missing files.

    #Format warning messages:
    def my_formatwarning(msg, *args, **kwargs):
        """Issue `msg` as warning."""
        return str(msg) + '\n'
    warnings.formatwarning = my_formatwarning

    if len(fils) == 0:
        warnings.warn("Input file list is empty.")
        return None
    if len(fils) > 1:
        return xr.open_mfdataset(fils, combine='by_coords')
    else:
        return xr.open_dataset(fils[0])
    #End if
#End def

########
