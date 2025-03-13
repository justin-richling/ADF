def _regrid_BAD(model_dataset, var_name, comp, method, **kwargs):
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
    if comp == "atm":
        comp_grid = "ncol"
    if comp == "lnd":
        comp_grid = "lndgrid"

    #Extract variable info from model data (and remove any degenerate dimensions):
    mdata = model_dataset[var_name].squeeze()

    if "wgt_file" in kwargs:
        weight_file = kwargs["wgt_file"]
    if "latlon_file" in kwargs:
        latlon_file = kwargs["latlon_file"]
    else:
        print("Well, it looks like you're missing a target grid file for regridding!")
        #adferror thing

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
        d_data = fv_ds[var_name]
        #d_data = fv_ds[var_name] if var_name in fv_ds else fv_ds

    # Create regridder
    regridder = make_se_regridder_BAD(weight_file=weight_file,
                                      s_data = s_data, #model_dataset.landmask.isel(time=0),
                                      d_data = d_data, #fv_ds.landmask,
                                      Method = method,  # Bug in xesmf needs this without "n"
                                      )
    if method == 'coservative':
        rgdata = regrid_se_data_conservative(regridder, model_dataset, comp_grid)


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

    if comp == "lnd":
        rgdata[var_name] = (rgdata[var_name] / rgdata.landfrac)
        rgdata['landmask'] = fv_ds.landmask
        rgdata['landfrac'] = rgdata.landfrac.isel(time=0)

    # Ensure output matches the target grid dimensions
    rgdata["lat"] = fv_ds.lat
    rgdata["lon"] = fv_ds.lon

    # Compute grid cell area (optional but useful for post-processing)
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