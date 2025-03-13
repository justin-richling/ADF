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
        #**regridder_kwargs
    )
    return regridder