# Regrids unstructured SE grid to regular lat-lon
# Shamelessly borrowed from @maritsandstad with NorESM who deserves credit for this work
# https://github.com/NorESMhub/xesmf_clm_fates_diagnostic/blob/main/src/xesmf_clm_fates_diagnostic/plotting_methods.py

import xarray as xr
import xesmf
import numpy as np

def make_se_regridder(weight_file, s_data, d_data,
                      Method='coservative'
                      ):
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
    )
    return regridder

def regrid_se_data_bilinear(regridder, data_to_regrid):
    updated = data_to_regrid.copy().transpose(..., "lndgrid").expand_dims("dummy", axis=-2)
    regridded = regridder(updated.rename({"dummy": "lat", "lndgrid": "lon"}),
                         skipna=True, na_thres=1,
                         )
    return regridded

def regrid_se_data_conservative(regridder, data_to_regrid):
    updated = data_to_regrid.copy().transpose(..., "lndgrid").expand_dims("dummy", axis=-2)
    regridded = regridder(updated.rename({"dummy": "lat", "lndgrid": "lon"}) )
    return regridded