



'''def unstructure_regrid(model_dataset, var_name, comp="atm", **kwargs):

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
        """#Extract keyword arguments:
        if 'ps_file' in kwargs:
            ps_file = kwargs['ps_file']
        else:
            ps_file = None
        #End if"""
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

        """#Extract grid info from target data:
        if 'time' in regrid_dataset.coords:
            if 'lev' in regrid_dataset.coords:
                tgrid = regrid_dataset.isel(time=0, lev=0).squeeze()
            else:
                tgrid = regrid_dataset.isel(time=0).squeeze()
            #End if
        #End if"""

        # Hardwiring for now
        #con_weight_file = "/glade/work/wwieder/map_ne30pg3_to_fv0.9x1.25_scripgrids_conserve_nomask_c250108.nc"
        if "wgt_file" in kwargs:
            weight_file = kwargs["wgt_file"]


        #fv_file = '/glade/derecho/scratch/wwieder/ctsm5.3.018_SP_f09_t232_mask/run/ctsm5.3.018_SP_f09_t232_mask.clm2.h0.0001-01.nc'
        if "latlon_file" in kwargs:
            fv_file = kwargs["latlon_file"]
        else:
            print("Well, it looks like you're missing a target grid file for regridding!")
            #adferror thing
        fv_ds = xr.open_dataset(fv_file)

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
                                      Method = 'coservative',  # Bug in xesmf needs this without "n"
                                      )

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




# Regrids unstructured SE grid to regular lat-lon
# Shamelessly borrowed from @maritsandstad with NorESM who deserves credit for this work
# https://github.com/NorESMhub/xesmf_clm_fates_diagnostic/blob/main/src/xesmf_clm_fates_diagnostic/plotting_methods.py

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
        #weights=weight_file,
        # results seem insensitive to this method choice
        # choices are coservative_normed, coservative, and bilinear
        #method=Method,
        reuse_weights=True,
        periodic=True,
        **regridder_kwargs
    )
    return regridder

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
'''
import xarray as xr
import xesmf
import numpy as np

def unstructure_regrid(model_dataset, var_name, comp, weight_file, latlon_file, method):
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
    #print("\n\nmodel_dataset and var:",var_name,model_dataset,"\n\n")
    mdata = model_dataset[var_name].squeeze()

    # Load target grid (lat/lon) from the provided dataset
    fv_ds = xr.open_dataset(latlon_file)

    mdata = mdata.fillna(0)
    if comp == "lnd":
        model_dataset['landfrac']= model_dataset['landfrac'].fillna(0)
        mdata = mdata * model_dataset.landfrac  # weight flux by land frac
        #print("\n\nmodel_dataset.landmask:",model_dataset.landmask,"\n\n")
        if 'time' in model_dataset.landmask.variable.dims:
            s_data = model_dataset.landmask.isel(time=0)
        else:
            s_data = model_dataset.landmask
        d_data = fv_ds.landmask
    else:
        print("comp",comp)
        s_data = mdata.isel(time=0)
        # Get the first data variable with 'time' and 'ncol' dimensions
        matching_vars = [var for var in fv_ds.data_vars if set(fv_ds[var].dims) == {'time', 'lndgrid'}]

        if matching_vars:
            first_matching_var = matching_vars[0]
            print(f"First matching variable: {first_matching_var}")
        d_data = fv_ds[first_matching_var]

    print("\nWeights file? ",weight_file,"\n")

    #Regrid model data to match target grid:
    regridder = make_se_regridder(weight_file=weight_file,
                                    s_data = s_data,
                                    d_data = d_data,
                                    Method = method,
                                    )
    if method == 'coservative':
        rgdata = regrid_se_data_conservative(regridder, model_dataset, comp_grid)
    if method == 'bilinear':
        rgdata = regrid_se_data_bilinear(regridder, model_dataset, comp_grid)

    if comp == "lnd":
        rgdata[var_name] = (rgdata[var_name] / rgdata.landfrac)
        rgdata['landmask'] = fv_ds.landmask
        if 'time' in rgdata.landfrac.variable.dims:
            rgdata['landfrac'] = rgdata.landfrac.isel(time=0)
        else:
            rgdata['landfrac'] = rgdata.landfrac

    rgdata['lat'] = fv_ds.lat
    rgdata['lon'] = fv_ds.lon

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

    #Return dataset:
    return rgdata



def make_se_regridder(weight_file, s_data, d_data,
                      Method='coservative'
                      ):
    """
    Create xESMF regridder for spectral element grids.
    """

    if weight_file:
        weights = xr.open_dataset(weight_file)
    
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
    # Handle source and destination masks
    s_mask = xr.DataArray(s_data.data.reshape(in_shape[0],in_shape[1]), dims=("lat", "lon"))
    dummy_in['mask']= s_mask
    
    d_mask = xr.DataArray(d_data.values, dims=("lat", "lon"))  
    dummy_out['mask']= d_mask                

    # QUESTION: Do source and destination grids need masks here?
    # See xesmf docs https://xesmf.readthedocs.io/en/stable/notebooks/Masking.html#Regridding-with-a-mask
    regridder = xesmf.Regridder(
        dummy_in,
        dummy_out,
        weights=weight_file,
        # NOTE: results seem insensitive to this method choice
        # INFO: choices are coservative_normed, coservative, and bilinear
        method=Method,
        reuse_weights=True,
        periodic=True,
    )
    return regridder


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