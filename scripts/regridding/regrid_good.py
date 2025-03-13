def _regrid(model_dataset, var_name, comp, method, **kwargs):
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