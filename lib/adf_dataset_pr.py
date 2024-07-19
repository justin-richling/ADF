from pathlib import Path
import xarray as xr

import warnings # use to warn user about missing files

def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = my_formatwarning

# "reference data"
# It is often just a "baseline case", 
# but could be some totally external data (reanalysis or observation or other model)
# When it is another simulation, it gets treated like another "case"
# When it is external data expect:
# - "climo" files (12 monthly climos in the file)
# - one variable per "climo"
# - source can differ for each variable, requires label
# - resolution can differ for each variable, requires regridded file(s)
# - the variable name and units in the file may differ from CAM; use defaults.yaml to set conversion
# - there could be multiple instances of a variable from different sources (e.g. different observations)

# NOTE: the last item (multiple instances of a variable) is not allowed in AdfObs.var_obs_dict
#       Since ADF is not able to handle this case, for now it is excluded the AdfData class.

# NOTE: To make the "baseline case" vs "external data" cases as similar as possible,
#       below construct the "baseline case" version to be similar to "external data".
#       - provide a dictionary of (variable: file-path)
#         + For external data, that dictionay is from AdfObs.var_obs_dict,
#           which provides a dict of all the available variables. 
#         + For reference simulation, look for files that match the diag_var_list

# NOTE: There is currently a "base_nickname" allowed from AdfInfo. 
#       Set AdfData.ref_nickname to that.
#       Could be altered from "Obs" to be the data source label.

class AdfData:
    """A class instantiated with an AdfDiag object. 
       Methods provide means to load data. 
       This class does not interact with plotting, 
       just provides access to data locations and loading data.

       A future need is to add some kind of frequency/sampling
       parameters to allow for non-h0 files. 

    """
    def __init__(self, adfobj):
        self.adf = adfobj  # provides quick access to the AdfDiag object
        # paths 
        self.model_rgrid_loc = adfobj.get_basic_info("cam_regrid_loc", required=True)

        # variables (and info for unit transform)
        # use self.adf.diag_var_list and self.adf.self.adf.variable_defaults

        # case names and nicknames
        self.case_names = adfobj.get_cam_info("cam_case_name", required=True)
        self.test_nicknames = adfobj.case_nicknames["test_nicknames"]
        self.base_nickname = adfobj.case_nicknames["base_nickname"]
        self.ref_nickname = self.base_nickname

        # define reference data
        self.set_reference(init=True) # specify "ref_labels" -> called "data_list" in zonal_mean (name of data source)

    # Reference case setup (baseline/obs)
    #------------------------------------
    def set_reference(self, init=False):
        """Set attributes for reference (aka baseline) data location, names, and variables."""
        if self.adf.compare_obs:
            self.ref_var_loc = {v: self.adf.var_obs_dict[v]['obs_file'] for v in self.adf.var_obs_dict}
            self.ref_labels = {v: self.adf.var_obs_dict[v]['obs_name'] for v in self.adf.var_obs_dict}
            self.ref_var_nam = {v: self.adf.var_obs_dict[v]['obs_var'] for v in self.adf.var_obs_dict}
            if not self.adf.var_obs_dict:
                warnings.warn("\t WARNING: reference is observations, but no observations found to plot against.")
        else:
            self.ref_var_loc = {}
            self.ref_var_nam = {}
            self.ref_labels = {}
            # when using a reference simulation, allow a "special" attribute with the case name:
            self.ref_case_label = self.adf.get_baseline_info("cam_case_name", required=True)
            for v in self.adf.diag_var_list:
                self.ref_var_nam[v] = v
                self.ref_labels[v] = self.adf.get_baseline_info("cam_case_name", required=True)
                f = self.get_reference_climo_file(v)
                if f is None:
                    if not init:
                        warnings.warn(f"\t WARNING: ADFData found no reference climo file for {v}")
                    continue
                else:
                    self.ref_labels[v] = self.adf.get_baseline_info("cam_case_name", required=True)


    # Time series files
    #------------------
    def get_timeseries_file(self, case, field):
        ts_locs = self.adf.get_cam_info("cam_ts_loc", required=True) # list of paths (could be multiple cases)
        caseindex = (self.case_names).index(case)
        ts_loc = Path(ts_locs[caseindex])
        ts_filenames = f'{case}.*.{field}.*nc'
        ts_files = sorted(ts_loc.glob(ts_filenames))
        return ts_files


    def get_ref_timeseries_file(self, field):
        if self.adf.compare_obs:
            return None
        else:
            ts_loc = Path(self.adf.get_baseline_info("cam_ts_loc", required=True))
            ts_filenames = f'{self.ref_case_label}.*.{field}.*nc'
            ts_files = sorted(ts_loc.glob(ts_filenames))
            return ts_files


    def load_timeseries_dataset(self, fils):
        if (len(fils) == 0):
            warnings.warn("Input file list is empty.")
            return None
        elif (len(fils) > 1):
            ds = xr.open_mfdataset(fils, decode_times=False)
        else:
            sfil = str(fils[0])
            if not Path(sfil).is_file():
                warnings.warn(f"Expecting to find file: {sfil}")
                return None
            ds = xr.open_dataset(sfil, decode_times=False)
        if ds is None:
            warnings.warn(f"invalid data on load_dataset")
        # assign time to midpoint of interval (even if it is already)
        if 'time_bnds' in ds:
            t = ds['time_bnds'].mean(dim='nbnd')
            t.attrs = ds['time'].attrs
            ds = ds.assign_coords({'time':t})
        elif 'time_bounds' in ds:
            t = ds['time_bounds'].mean(dim='nbnd')
            t.attrs = ds['time'].attrs
            ds = ds.assign_coords({'time':t})
        else:
            warnings.warn("Timeseries file does not have time bounds info.")
        return xr.decode_cf(ds)

    #----------------

    
    # Climatology files
    #------------------
    def get_climo_file(self, case, variablename):
        """Retrieve the climo file path(s) for variablename for a specific case."""
        a = self.adf.get_cam_info("cam_climo_loc", required=True) # list of paths (could be multiple cases)
        caseindex = (self.case_names).index(case) # the entry for specified case
        model_cl_loc = Path(a[caseindex])
        return sorted(model_cl_loc.glob(f"{case}_{variablename}_climo.nc"))


    def get_reference_climo_file(self, var):
        """Return a list of files to be used as reference (aka baseline) for variable var."""
        if self.adf.compare_obs:
            fils = self.ref_var_loc.get(var, None)
            return [fils] if fils is not None else None
        ref_loc = self.adf.get_baseline_info("cam_climo_loc")
        # NOTE: originally had this looking for *_baseline.nc
        fils = sorted(Path(ref_loc).glob(f"{self.ref_case_label}_{var}_climo.nc"))
        if fils:
            return fils
        return None


    def load_climo_da(self, case, variablename):
        """Return DataArray from climo file"""
        fils = self.get_climo_file(case, variablename)
        return self.load_da(case, fils, variablename)


    def load_climo_file(self, case, variablename):
        """Return Dataset for climo of variablename"""
        fils = self.get_climo_file(case, variablename)
        if not fils:
            warnings.warn(f"ERROR: Did not find climo file for variable: {variablename}. Will try to skip.")
            return None
        return self.load_dataset(fils)


    def load_obs_climo_dataset(self, variablename):
        """Return Dataset for observation climo of variablename"""
        fils = self.get_reference_climo_file(variablename)
        if not fils:
            warnings.warn(f"ERROR: Did not find any reference climo files for variable: {variablename}. Will try to skip.")
            return None
        return self.load_dataset(fils)

    #----------------


    # Regridded files
    #----------------
    def load_reference_regrid_dataset(self, case, field):
        fils = self.get_ref_regrid_file(case, field)
        if not fils:
            warnings.warn(f"ERROR: Did not find regrid file(s) for case: {case}, variable: {field}")
            return None
        return self.load_dataset(fils)


    def load_reference_regrid_da(self, case, field):
        fils = self.get_ref_regrid_file(case, field)
        if not fils:
            warnings.warn(f"ERROR: Did not find regrid file(s) for case: {case}, variable: {field}")
            return None
        return self.load_da(case, fils, field)


    def get_ref_regrid_file(self, case, field):
        if self.adf.compare_obs:
            obs_loc = self.ref_var_loc.get(field, None)
            fils = [str(obs_loc)]
        else:
            model_rg_loc = Path(self.adf.get_basic_info("cam_regrid_loc", required=True))
            fils = sorted(model_rg_loc.glob(f"{case}_{field}_*.nc"))
        return fils
    

    def get_regrid_file(self, case, field):
        model_rg_loc = Path(self.adf.get_basic_info("cam_regrid_loc", required=True))
        # rlbl = "reference label" = the name of the reference data that defines target grid
        rlbl = self.ref_labels[field]
        return sorted(model_rg_loc.glob(f"{rlbl}_{case}_{field}_*.nc"))

    
    def load_regrid_dataset(self, case, field):
        fils = self.get_regrid_file(case, field)
        if not fils:
            warnings.warn(f"ERROR: Did not find regrid file(s) for case: {case}, variable: {field}")
            return None
        return self.load_dataset(fils)

    
    def load_regrid_da(self, case, field):
        fils = self.get_regrid_file(case, field)
        if not fils:
            warnings.warn(f"ERROR: Did not find regrid file(s) for case: {case}, variable: {field}")
            return None
        return self.load_da(case, fils, field)
    
    #----------------


    # DataSet and DataArray load
    #---------------------------

    # Load DataSet
    def load_dataset(self, fils):
        if (len(fils) == 0):
            warnings.warn("Input file list is empty.")
            return None
        elif (len(fils) > 1):
            ds = xr.open_mfdataset(fils, combine='by_coords')
        else:
            sfil = str(fils[0])
            if not Path(sfil).is_file():
                warnings.warn(f"Expecting to find file: {sfil}")
                return None
            ds = xr.open_dataset(sfil)
        if ds is None:
            warnings.warn(f"invalid data on load_dataset")
        return ds


    # Load DataArray
    def load_da(self, case, fils, variablename):
        #Check if case is baseline and if it is, check if comparing against obs
        if (case == self.ref_labels[variablename]) and (self.adf.compare_obs):
            da = self.load_obs_climo_dataset(variablename)[self.ref_var_nam[variablename]]
        #Else, its either a test case, or baseline and NOT comparing against obs
        else:
            ds = self.load_dataset(fils)
            da = (ds[variablename]).squeeze()

        if variablename in self.adf.variable_defaults:
            vres = self.adf.variable_defaults[variablename]
            if (case == self.ref_labels[variablename]) and (self.adf.compare_obs):
                scale_factor = vres.get("obs_scale_factor",1)
                add_offset = vres.get("obs_add_offset", 0)
            else:
                scale_factor = vres.get("scale_factor",1)
                add_offset = vres.get("add_offset", 0)
            da = da * scale_factor + add_offset
            da.attrs['units'] = vres.get("new_unit", da.attrs.get('units', 'none'))
        return da

 
    #----------------

#End Script
##########
    
    
