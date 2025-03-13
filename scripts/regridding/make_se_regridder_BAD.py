def make_se_regridder_BAD(weight_file, s_data, d_data, Method='conservative'):
    """
    Create xESMF regridder for spectral element grids.
    """
    regridder_kwargs = {}

    # Load weights if available
    if weight_file:
        weights = xr.open_dataset(weight_file)
        regridder_kwargs['weights'] = weights
    else:
        print("No weights file given, so I'm gonna need to make one. Please have a seat and the next associate will be with you shortly. Please don't tap the glass!")
        regridder_kwargs['method'] = 'coservative'

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