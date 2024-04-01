"""                                                                    .
Generic computation and plotting helper functions

Functions
---------
use_this_norm()
    switches matplotlib color normalization method
get_difference_colors(values)
    Provide a color norm and colormap assuming `values` is a difference field.
mask_land_or_ocean(arr, msk, use_nan=False)
    Apply a land or ocean mask to provided variable.
get_central_longitude(*args)
    Determine central longitude for maps.
global_average(fld, wgt, verbose=False)
    pure numpy global average.
spatial_average(indata, weights=None, spatial_dims=None)
    Compute spatial average
wgt_rmse(fld1, fld2, wgt):
    Calculate the area-weighted RMSE.
annual_mean(data, whole_years=False, time_name='time'):
    Calculate annual averages from time series data.
seasonal_mean(data, season=None, is_climo=None):
    Calculates the time-weighted seasonal average (or average over all time).
domain_stats(data, domain):
    Provides statistics in specified region.
make_polar_plot(wks, case_nickname, base_nickname,
                    case_climo_yrs, baseline_climo_yrs,
                    d1:xr.DataArray, d2:xr.DataArray, difference:Optional[xr.DataArray]=None,
                    domain:Optional[list]=None, hemisphere:Optional[str]=None, **kwargs):
    Make a stereographic polar plot for the given data and hemisphere.
plot_map_vect_and_save(wks, case_nickname, base_nickname,
                           case_climo_yrs, baseline_climo_yrs,
                           plev, umdlfld_nowrap, vmdlfld_nowrap,
                           uobsfld_nowrap, vobsfld_nowrap,
                           udiffld_nowrap, vdiffld_nowrap, **kwargs):
    Plots a vector field on a map.
plot_map_and_save(wks, case_nickname, base_nickname,
                      case_climo_yrs, baseline_climo_yrs,
                      mdlfld, obsfld, diffld, **kwargs):
    Map plots of `mdlfld`, `obsfld`, and their difference, `diffld`.
pres_from_hybrid(psfc, hya, hyb, p0=100000.):
    Converts a hybrid level to a pressure
vert_remap(x_mdl, p_mdl, plev)
    Interpolates to specified pressure levels.
lev_to_plev(data, ps, hyam, hybm, P0=100000., new_levels=None, convert_to_mb=False)
    Interpolate model hybrid levels to specified pressure levels.
pmid_to_plev(data, pmid, new_levels=None, convert_to_mb=False)
    Interpolate `data` from hybrid-sigma levels to isobaric levels using provided mid-level pressures.
zonal_mean_xr(fld)
    Average over all dimensions except `lev` and `lat`.
validate_dims(fld, list_of_dims)
    Checks if specified dimensions are in a DataArray
lat_lon_validate_dims(fld)
    Check if input field has lat and lon.
zm_validate_dims(fld)
    Check for dimensions for zonal average.
zonal_plot(lat, data, ax=None, color=None, **kwargs)
    Make a line plot or pressure-latitude plot of `data`.
meridional_plot(lon, data, ax=None, color=None, **kwargs)
    Make a line plot or pressure-longitude plot of `data`.
prep_contour_plot
    Preparation for making contour plots.
plot_zonal_mean_and_save
    zonal mean plot
plot_meridional_mean_and_save
    meridioanl mean plot
square_contour_difference
    Produce filled contours of fld1, fld2, and their difference with square axes.

Notes
-----
This module includes several "private" methods intended for internal use only.

_plot_line(axobject, xdata, ydata, color, **kwargs)
    Create a generic line plot
_meridional_plot_line

_zonal_plot_line

_zonal_plot_preslat

_meridional_plot_preslon

"""

#import statements:
from signal import SIG_DFL
from typing import Optional
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import cartopy.crs as ccrs
#nice formatting for tick labels
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
import geocat.comp as gcomp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D

from adf_diag import AdfDiag
from adf_base import AdfError

#Set non-X-window backend for matplotlib:
mpl.use('Agg')

#Now import pyplot:
import matplotlib.pyplot as plt

empty_message = "No Valid\nData Points"
props = {'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.9}


#Set seasonal ranges:
seasons = {"ANN": np.arange(1,13,1),
            "DJF": [12, 1, 2],
            "JJA": [6, 7, 8],
            "MAM": [3, 4, 5],
            "SON": [9, 10, 11]
            }


#################
#HELPER FUNCTIONS
#################

def use_this_norm():
    """Just use the right normalization; avoids a deprecation warning."""

    mplversion = [int(x) for x in mpl.__version__.split('.')]
    if mplversion[0] < 3:
        return mpl.colors.Normalize, mplversion[0]
    else:
        if mplversion[1] < 2:
            return mpl.colors.DivergingNorm, mplversion[0]
        else:
            return mpl.colors.TwoSlopeNorm, mplversion[0]


def get_difference_colors(values):
    """Provide a color norm and colormap assuming this is a difference field.

    Parameters
    ----------
    values : array-like
        can be either the data field or a set of specified contour levels.

    Returns
    -------
    dnorm
        Matplotlib color nomalization
    cmap
        Matplotlib colormap

    Notes
    -----
    Uses 'OrRd' colormap for positive definite, 'BuPu_r' for negative definite,
    and 'RdBu_r' centered on zero if there are values of both signs.
    """
    normfunc, mplv = use_this_norm()
    dmin = np.min(values)
    dmax = np.max(values)
    # color normalization for difference
    if ((dmin < 0) and (0 < dmax)):
        dnorm = normfunc(vmin=np.min(values), vmax=np.max(values), vcenter=0.0)
        cmap = mpl.cm.RdBu_r
    else:
        dnorm = mpl.colors.Normalize(vmin=np.min(values), vmax=np.max(values))
        if dmin >= 0:
            cmap = mpl.cm.OrRd
        elif dmax <= 0:
            cmap = mpl.cm.BuPu_r
        else:
            dnorm = mpl.colors.TwoSlopeNorm(vmin=dmin, vcenter=0, vmax=dmax)
    return dnorm, cmap


def mask_land_or_ocean(arr, msk, use_nan=False):
    """Apply a land or ocean mask to provided variable.

    Parameters
    ----------
    arr : xarray.DataArray
        the xarray variable to apply the mask to.
    msk : xarray.DataArray
        the xarray variable that contains the land or ocean mask,
        assumed to be the same shape as "arr".
    use_nan : bool, optional
        argument for whether to set the missing values
        to np.nan values instead of the defaul "-999." values.

    Returns
    -------
    arr : xarray.DataArray
        Same as input `arr` but masked as specified.
    """

    if use_nan:
        missing_value = np.nan
    else:
        missing_value = -999.
    #End if

    arr = xr.where(msk>=0.9,arr,missing_value)
    arr.attrs["missing_value"] = missing_value
    return(arr)


def get_central_longitude(*args):
    """Determine central longitude for maps.

    Allows an arbitrary number of arguments.
    If any of the arguments is an instance of `AdfDiag`, then check
    whether it has a `central_longitude` in `diag_basic_info`.
    _This case takes precedence._
    _Else_, if any of the arguments are scalars in [-180, 360],
    assumes the FIRST ONE is the central longitude.
    There are no other possible conditions, so if none of those are met,
    returns the default value of 180.

    Parameters
    ----------
    *args : tuple
        Any number of objects to check for `central_longitude`.
        After that, looks for the first number between -180 and 360 in the args.

    Notes
    -----
    This allows a script to, for example, allow a config file to specify, but also have a preference:
    `get_central_longitude(AdfObj, 30.0)`
    """
    chk_for_adf = [isinstance(arg, AdfDiag) for arg in args]
    # preference is to get value from AdfDiag:
    if any(chk_for_adf):
        for arg in args:
            if isinstance(arg, AdfDiag):
                result = arg.get_basic_info('central_longitude', required=False)
                if (isinstance(result, int) or isinstance(result, float)) and \
                   (result >= -180) and (result <= 360):
                    return result
                else:
                    #If result exists, then write info to debug log:
                    if result:
                        msg = f"central_lngitude of type '{type(result).__name__}'"
                        msg += f" and value '{result}', which is not a valid longitude"
                        msg += " for the ADF."
                        arg.debug_log(msg)
                    #End if

                    #There is only one ADF object per ADF run, so if its
                    #not present or configured correctly then no
                    #reason to keep looking:
                    break
                #End if
            #End if
        #End for
    #End if

    # 2nd pass through arguments, look for numbers:
    for arg in args:
        if (isinstance(arg, float) or isinstance(arg, int)) and ((arg >= -180) and (arg <= 360)):
            return arg
        #End if
    else:
        # this is the `else` on the for loop --> if non of the arguments meet the criteria, do this.
        print("No valid central longitude specified. Defaults to 180.")
        return 180
    #End if

#######

def coslat_average(darray, lat1, lat2):
    """
    Calculate the weighted average for an [:,lat] array over the region
    lat1 to lat2
    """

    # flip latitudes if they are decreasing
    if (darray.lat[0] > darray.lat[darray.lat.size -1]):
        print("flipping latitudes")
        darray = darray.sortby('lat')

    region = darray.sel(lat=slice(lat1, lat2))
    weights=np.cos(np.deg2rad(region.lat))
    regionw = region.weighted(weights)
    regionm = regionw.mean("lat")

    return regionm
#######

def global_average(fld, wgt, verbose=False):
    """A simple, pure numpy global average.

    Parameters
    ----------
    fld : np.ndarray
        an input ndarray
    wgt : np.ndarray
        a 1-dimensional array of weights, should be same size as one dimension of `fld`
    verbose : bool, optional
        prints information when `True`

    Returns
    -------
    weighted average of `fld`
    """

    s = fld.shape
    for i in range(len(s)):
        if np.size(fld, i) == len(wgt):
            a = i
            break
    fld2 = np.ma.masked_invalid(fld)
    if verbose:
        print("(global_average)-- fraction of mask that is True: {}".format(np.count_nonzero(fld2.mask) / np.size(fld2)))
        print("(global_average)-- apply ma.average along axis = {} // validate: {}".format(a, fld2.shape))
    avg1, sofw = np.ma.average(fld2, axis=a, weights=wgt, returned=True) # sofw is sum of weights

    return np.ma.average(avg1)


def spatial_average(indata, weights=None, spatial_dims=None):
    """Compute spatial average.

    Parameters
    ----------
    indata : xr.DataArray
        input data
    weights : np.ndarray or xr.DataArray, optional
        the weights to apply, see Notes for default behavior
    spatial_dims : list, optional
        list of dimensions to average, see Notes for default behavior

    Returns
    -------
    xr.DataArray
        weighted average of `indata`

    Notes
    -----
    When `weights` is not provided, tries to find sensible values.
    If there is a 'lat' dimension, use `cos(lat)`.
    If there is a 'ncol' dimension, looks for `area` in `indata`.
    Otherwise, set to equal weights.

    Makes an attempt to identify the spatial variables when `spatial_dims` is None.
    Will average over `ncol` if present, and then will check for `lat` and `lon`.
    When none of those three are found, raise an AdfError.
    """
    import warnings

    if weights is None:
        #Calculate spatial weights:
        if 'lat' in indata.coords:
            weights = np.cos(np.deg2rad(indata.lat))
            weights.name = "weights"
        elif 'ncol' in indata.dims:
            if 'area' in indata:
                warnings.warn("area variable being used to generated normalized weights.")
                weights = indata['area'] / indata['area'].sum()
            else:
                warnings.warn("We need a way to get area variable. Using equal weights.")
                weights = xr.DataArray(1.)
            weights.name = "weights"
        else:
            weights = xr.DataArray(1.)
            weights.name = "weights"
            warnings.warn("Un-recognized spatial dimensions: using equal weights for all grid points.")
        #End if
    #End if

    #Apply weights to input data:
    weighted = indata.weighted(weights)

    # we want to average over all non-time dimensions
    if spatial_dims is None:
        if 'ncol' in indata.dims:
            spatial_dims = ['ncol']
        else:
            spatial_dims = [dimname for dimname in indata.dims if (('lat' in dimname.lower()) or ('lon' in dimname.lower()))]

    if not spatial_dims:
        #Scripts using this function likely expect the horizontal dimensions
        #to be removed via the application of the mean. So in order to avoid
        #possibly unexpected behavior due to arrays being incorrectly dimensioned
        #(which could be difficult to debug) the ADF should die here:
        emsg = "spatial_average: No spatial dimensions were identified,"
        emsg += " so can not perform average."
        raise AdfError(emsg)

    return weighted.mean(dim=spatial_dims)


def wgt_rmse(fld1, fld2, wgt):
    """Calculate the area-weighted RMSE.

    Parameters
    ----------
    fld1, fld2 : array-like
        2-dimensional spatial fields with the same shape.
        They can be xarray DataArray or numpy arrays.
    wgt : array-like
        the weight vector, expected to be 1-dimensional,
        matching length of one dimension of the data.

    Returns
    -------
    float
        root mean squared error

    Notes:
    ```rmse = sqrt( mean( (fld1 - fld2)**2 ) )```
    """
    assert len(fld1.shape) == 2,     "Input fields must have exactly two dimensions."
    assert fld1.shape == fld2.shape, "Input fields must have the same array shape."
    # in case these fields are in dask arrays, compute them now.
    if hasattr(fld1, "compute"):
        fld1 = fld1.compute()
    if hasattr(fld2, "compute"):
        fld2 = fld2.compute()
    if isinstance(fld1, xr.DataArray) and isinstance(fld2, xr.DataArray):
        return (np.sqrt(((fld1 - fld2)**2).weighted(wgt).mean())).values.item()
    else:
        check = [len(wgt) == s for s in fld1.shape]
        if ~np.any(check):
            raise IOError(f"Sorry, weight array has shape {wgt.shape} which is not compatible with data of shape {fld1.shape}")
        check = [len(wgt) != s for s in fld1.shape]
        dimsize = fld1.shape[np.argwhere(check).item()]  # want to get the dimension length for the dim that does not match the size of wgt
        warray = np.tile(wgt, (dimsize, 1)).transpose()   # May need more logic to ensure shape is correct.
        warray = warray / np.sum(warray) # normalize
        wmse = np.sum(warray * (fld1 - fld2)**2)
        return np.sqrt( wmse ).item()


#######
# Time-weighted averaging

def annual_mean(data, whole_years=False, time_name='time'):
    """Calculate annual averages from monthly time series data.

    Parameters
    ----------
    data : xr.DataArray or xr.Dataset
        monthly data values with temporal dimension
    whole_years : bool, optional
        whether to restrict endpoints of the average to
        start at first January and end at last December
    time_name : str, optional
        name of the time dimension, defaults to `time`

    Returns
    -------
    result : xr.DataArray or xr.Dataset
        `data` reduced to annual averages

    Notes
    -----
    This function assumes monthly data, and weights the average by the
    number of days in each month.

    `result` includes an attribute that reports the date range used for the average.
    """
    assert time_name in data.coords, f"Did not find the expected time coordinate '{time_name}' in the data"
    if whole_years:
        first_january = np.argwhere((data.time.dt.month == 1).values)[0].item()
        last_december = np.argwhere((data.time.dt.month == 12).values)[-1].item()
        data_to_avg = data.isel(time=slice(first_january,last_december+1)) # PLUS 1 BECAUSE SLICE DOES NOT INCLUDE END POINT
    else:
        data_to_avg = data
    date_range_string = f"{data_to_avg['time'][0]} -- {data_to_avg['time'][-1]}"

    # this provides the normalized monthly weights in each year
    # -- do it for each year to allow for non-standard calendars (360-day)
    # -- and also to provision for data with leap years
    days_gb = data_to_avg.time.dt.daysinmonth.groupby('time.year').map(lambda x: x / x.sum())
    # weighted average with normalized weights: <x> = SUM x_i * w_i  (implied division by SUM w_i)
    result =  (data_to_avg * days_gb).groupby('time.year').sum(dim='time')
    result.attrs['averaging_period'] = date_range_string
    return result


def seasonal_mean(data, season=None, is_climo=None):
    """Calculates the time-weighted seasonal average (or average over all time).

    Parameters
    ----------
    data : xarray.DataArray or xarray.Dataset
        data to be averaged
    season : str, optional
        the season to extract from `data`
        If season is `ANN` or None, average all available time.
    is_climo : bool, optional
        If True, expects data to have time or month dimenion of size 12.
        If False, then 'time' must be a coordinate,
        and the `time.dt.days_in_month` attribute must be available.

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        the average of `data` in season `season`

    Notes
    -----
    If the data is a climatology, the code will make an attempt to understand the time or month
    dimension, but will assume that it is ordered from January to December.
    If the data is a climatology and is just a numpy array with one dimension that is size 12,
    it will assume that dimension is time running from January to December.
    """
    if season is not None:
        assert season in ["ANN", "DJF", "JJA", "MAM", "SON"], f"Unrecognized season string provided: '{season}'"
    elif season is None:
        season = "ANN"

    try:
        month_length = data.time.dt.days_in_month
    except (AttributeError, TypeError):
        # do our best to determine the temporal dimension and assign weights
        if not is_climo:
            raise ValueError("Non-climo file provided, but without a decoded time dimension.")
        else:
            # CLIMO file: try to determine which dimension is month
            has_time = False
            if isinstance(data, xr.DataArray):
                has_time = 'time' in data.dims
                if not has_time:
                    if "month" in data.dims:
                        data = data.rename({"month":"time"})
                        has_time = True
            if not has_time:
                # this might happen if a pure numpy array gets passed in
                # --> assumes ordered January to December.
                assert ((12 in data.shape) and (data.shape.count(12) == 1)), f"Sorry, {data.shape.count(12)} dimensions have size 12, making determination of which dimension is month ambiguous. Please provide a `time` or `month` dimension."
                time_dim_num = data.shape.index(12)
                fakedims = [f"dim{n}" for n in range(len(data.shape))]
                fakedims[time_dim_num] = "time"
                data = xr.DataArray(data, dims=fakedims, attrs=data.attrs)
            timefix = pd.date_range(start='1/1/1999', end='12/1/1999', freq='MS') # generic time coordinate from a non-leap-year
            data = data.assign_coords({"time":timefix})
        month_length = data.time.dt.days_in_month
    #End try/except

    data = data.sel(time=data.time.dt.month.isin(seasons[season])) # directly take the months we want based on season kwarg
    return data.weighted(data.time.dt.daysinmonth).mean(dim='time', keep_attrs=True)



#######

#Polar Plot functions

def domain_stats(data, domain):
    """Provides statistics in specified region.

    Parameters
    ----------
    data : xarray.DataArray
        data values
    domain : list or tuple or numpy.ndarray
        the domain specification as:
        [west_longitude, east_longitude, south_latitude, north_latitude]

    Returns
    -------
    x_region_mean : float
        the regional area-weighted average
    x_region_max : float
        the maximum value in the region
    x_region_min : float
        the minimum value in the region

    Notes
    -----
    Currently assumes 'lat' is a dimension and uses `cos(lat)` as weight.
    Should use `spatial_average`

    See Also
    --------
    spatial_average

    """
    x_region = data.sel(lat=slice(domain[2],domain[3]), lon=slice(domain[0],domain[1]))
    x_region_mean = x_region.weighted(np.cos(np.deg2rad(x_region['lat']))).mean().item()
    x_region_min = x_region.min().item()
    x_region_max = x_region.max().item()
    return x_region_mean, x_region_max, x_region_min

def make_polar_plot(wks, case_nickname, base_nickname,
                    case_climo_yrs, baseline_climo_yrs,
                    d1:xr.DataArray, d2:xr.DataArray, difference:Optional[xr.DataArray]=None,
                    domain:Optional[list]=None, hemisphere:Optional[str]=None, obs=False, **kwargs):

    """Make a stereographic polar plot for the given data and hemisphere.

    Parameters
    ----------
    wks : str or Path
        output file path
    case_nickname : str
        short case name for `d1`
    base_nickname : str
        short case name for `d2`
    case_climo_yrs : list
        years for case `d1`, used for annotation
    baseline_climo_yrs : list
        years for case `d2`, used for annotation
    d1, d2 : xr.DataArray
        input data, must contain dimensions `lat` and `lon`
    difference : xr.DataArray, optional
        data to use as the difference, otherwise `d2 - d1`
    domain : list, optional
        the domain to plot, specified as
        [west_longitude, east_longitude, south_latitude, north_latitude]
        Defaults to pole to 45-degrees, all longitudes
    hemisphere : {'NH', 'SH'}, optional
        Hemsiphere to plot
    kwargs : dict, optional
        variable-dependent options for plots, See Notes.

    Notes
    -----
    - Uses contourf. No contour lines (yet).
    - kwargs is checked for:
        + `colormap`
        + `contour_levels`
        + `contour_levels_range`
        + `diff_contour_levels`
        + `diff_contour_range`
        + `diff_colormap`
        + `units`
    """
    if difference is None:
        dif = d2 - d1
    else:
        dif = difference

    if hemisphere.upper() == "NH":
        proj = ccrs.NorthPolarStereo()
    elif hemisphere.upper() == "SH":
        proj = ccrs.SouthPolarStereo()
    else:
        raise AdfError(f'[make_polar_plot] hemisphere not specified, must be NH or SH; hemisphere set as {hemisphere}')

    if domain is None:
        if hemisphere.upper() == "NH":
            domain = [-180, 180, 45, 90]
        else:
            domain = [-180, 180, -90, -45]

    # statistics for annotation (these are scalars):
    d1_region_mean, d1_region_max, d1_region_min = domain_stats(d1, domain)
    d2_region_mean, d2_region_max, d2_region_min = domain_stats(d2, domain)
    dif_region_mean, dif_region_max, dif_region_min = domain_stats(dif, domain)

    #downsize to the specified region; makes plotting/rendering/saving much faster
    d1 = d1.sel(lat=slice(domain[2],domain[3]))
    d2 = d2.sel(lat=slice(domain[2],domain[3]))
    dif = dif.sel(lat=slice(domain[2],domain[3]))

    # add cyclic point to the data for better-looking plot
    d1_cyclic, lon_cyclic = add_cyclic_point(d1, coord=d1.lon)
    d2_cyclic, _ = add_cyclic_point(d2, coord=d2.lon)  # since we can take difference, assume same longitude coord.
    dif_cyclic, _ = add_cyclic_point(dif, coord=dif.lon)

    # -- deal with optional plotting arguments that might provide variable-dependent choices

    # determine levels & color normalization:
    minval    = np.min([np.min(d1), np.min(d2)])
    maxval    = np.max([np.max(d1), np.max(d2)])
    absmaxdif = np.max(np.abs(dif))

    if 'colormap' in kwargs:
        cmap1 = kwargs['colormap']
    else:
        cmap1 = 'coolwarm'

    if 'contour_levels' in kwargs:
        levels1 = kwargs['contour_levels']
        norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
    elif 'contour_levels_range' in kwargs:
        assert len(kwargs['contour_levels_range']) == 3, "contour_levels_range must have exactly three entries: min, max, step"
        levels1 = np.arange(*kwargs['contour_levels_range'])
        norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
    else:
        levels1 = np.linspace(minval, maxval, 12)
        norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)

    if ('colormap' not in kwargs) and ('contour_levels' not in kwargs):
        norm1, cmap1 = get_difference_colors(levels1)  # maybe these are better defaults if nothing else is known.

    if "diff_contour_levels" in kwargs:
        levelsdiff = kwargs["diff_contour_levels"]  # a list of explicit contour levels
    elif "diff_contour_range" in kwargs:
            assert len(kwargs['diff_contour_range']) == 3, "diff_contour_range must have exactly three entries: min, max, step"
            levelsdiff = np.arange(*kwargs['diff_contour_range'])
    else:
        # set levels for difference plot (with a symmetric color bar):
        levelsdiff = np.linspace(-1*absmaxdif, absmaxdif, 12)
    #End if

    #NOTE: Sometimes the contour levels chosen in the defaults file
    #can result in the "contourf" software stack generating a
    #'TypologyException', which should manifest itself as a
    #"PredicateError", but due to bugs in the stack itself
    #will also sometimes raise an AttributeError.

    #To prevent this from happening, the polar max and min values
    #are calculated, and if the default contour values are significantly
    #larger then the min-max values, then the min-max values are used instead:
    #-------------------------------
    if max(levels1) > 10*maxval:
        levels1 = np.linspace(minval, maxval, 12)
        norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    elif minval < 0 and min(levels1) < 10*minval:
        levels1 = np.linspace(minval, maxval, 12)
        norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    #End if

    if max(abs(levelsdiff)) > 10*absmaxdif:
        levelsdiff = np.linspace(-1*absmaxdif, absmaxdif, 12)
    #End if
    #-------------------------------

    # Difference options -- Check in kwargs for colormap and levels
    if "diff_colormap" in kwargs:
        cmapdiff = kwargs["diff_colormap"]
        dnorm, _ = get_difference_colors(levelsdiff)  # color map output ignored
    else:
        dnorm, cmapdiff = get_difference_colors(levelsdiff)
    #End if

    # -- end options

    lons, lats = np.meshgrid(lon_cyclic, d1.lat)

    fig = plt.figure(figsize=(10,10))
    gs = mpl.gridspec.GridSpec(2, 4, wspace=0.9)

    ax1 = plt.subplot(gs[0, :2], projection=proj)
    ax2 = plt.subplot(gs[0, 2:], projection=proj)
    ax3 = plt.subplot(gs[1, 1:3], projection=proj)

    levs = np.unique(np.array(levels1))
    levs_diff = np.unique(np.array(levelsdiff))

    print(lons.shape)
    print(lats.shape)
    print(d1_cyclic.shape)

    if len(levs) < 2:
        img1 = ax1.contourf(lons, lats, d1_cyclic, transform=ccrs.PlateCarree(), colors="w", norm=norm1)
        ax1.text(0.4, 0.4, empty_message, transform=ax1.transAxes, bbox=props)

        img2 = ax2.contourf(lons, lats, d2_cyclic, transform=ccrs.PlateCarree(), colors="w", norm=norm1)
        ax2.text(0.4, 0.4, empty_message, transform=ax2.transAxes, bbox=props)
    else:
        img1 = ax1.contourf(lons, lats, d1_cyclic, transform=ccrs.PlateCarree(), cmap=cmap1, norm=norm1, levels=levels1)
        img2 = ax2.contourf(lons, lats, d2_cyclic, transform=ccrs.PlateCarree(), cmap=cmap1, norm=norm1, levels=levels1)

    if len(levs_diff) < 2:
        img3 = ax3.contourf(lons, lats, dif_cyclic, transform=ccrs.PlateCarree(), colors="w", norm=dnorm)
        ax3.text(0.4, 0.4, empty_message, transform=ax3.transAxes, bbox=props)
    else:
        img3 = ax3.contourf(lons, lats, dif_cyclic, transform=ccrs.PlateCarree(), cmap=cmapdiff, norm=dnorm, levels=levelsdiff)

    #Set Main title for subplots:
    st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=18)
    st.set_y(0.95)

    #Set plot titles
    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"
    ax1.set_title(case_title, loc='left', fontsize=6) #fontsize=tiFontSize

    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
        ax2.set_title(base_title, loc='left', fontsize=6) #fontsize=tiFontSize
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"
        ax2.set_title(base_title, loc='left', fontsize=6)

    ax1.text(-0.2, -0.10, f"Mean: {d1_region_mean:5.2f}\nMax: {d1_region_max:5.2f}\nMin: {d1_region_min:5.2f}", transform=ax1.transAxes)

    ax2.text(-0.2, -0.10, f"Mean: {d2_region_mean:5.2f}\nMax: {d2_region_max:5.2f}\nMin: {d2_region_min:5.2f}", transform=ax2.transAxes)

    ax3.text(-0.2, -0.10, f"Mean: {dif_region_mean:5.2f}\nMax: {dif_region_max:5.2f}\nMin: {dif_region_min:5.2f}", transform=ax3.transAxes)
    ax3.set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=8)

    if "units" in kwargs:
        ax2.set_ylabel(kwargs["units"])
        ax3.set_ylabel(kwargs["units"])
    else:
        ax2.set_ylabel(f"{dif.units}")
        ax3.set_ylabel(f"{dif.units}")


    [a.set_extent(domain, ccrs.PlateCarree()) for a in [ax1, ax2, ax3]]
    [a.coastlines() for a in [ax1, ax2, ax3]]

    # __Follow the cartopy gallery example to make circular__:
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpl.path.Path(verts * radius + center)
    [a.set_boundary(circle, transform=a.transAxes) for a in [ax1, ax2, ax3]]

    # __COLORBARS__
    cb_mean_ax = inset_axes(ax2,
                    width="5%",  # width = 5% of parent_bbox width
                    height="90%",  # height : 90%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0.05, 1, 1),
                    bbox_transform=ax2.transAxes,
                    borderpad=0,
                    )
    fig.colorbar(img1, cax=cb_mean_ax)

    cb_diff_ax = inset_axes(ax3,
                    width="5%",  # width = 5% of parent_bbox width
                    height="90%",  # height : 90%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0.05, 1, 1),
                    bbox_transform=ax3.transAxes,
                    borderpad=0,
                    )
    fig.colorbar(img3, cax=cb_diff_ax)

    # Save files
    fig.savefig(wks, bbox_inches='tight', dpi=300)

    # Close figures to avoid memory issues:
    plt.close(fig)

#######

def plot_map_vect_and_save(wks, case_nickname, base_nickname,
                           case_climo_yrs, baseline_climo_yrs,
                           plev, umdlfld_nowrap, vmdlfld_nowrap,
                           uobsfld_nowrap, vobsfld_nowrap,
                           udiffld_nowrap, vdiffld_nowrap, obs=False, **kwargs):

    """This plots a vector plot.

    Vector fields constructed from x- and y-components (u, v).

    Parameters
    ----------
    wks : str or Path
        output file path
    case_nickname : str
        short name for case
    base_nickname : str
        short name for base case
    case_climo_yrs : list
        list of years in case climatology, used for annotation
    baseline_climo_yrs : list
        list of years in base case climatology, used for annotation
    plev : str or float or None
        if not None, label denoting the pressure level
    umdlfld_nowrap, vmdlfld_nowrap : xarray.DataArray
        input data for case, the x- and y- components of the vectors
    uobsfld_nowrap, vobsfld_nowrap : xarray.DataArray
        input data for base case, the x- and y- components of the vectors
    udiffld_nowrap, vdiffld_nowrap : xarray.DataArray
        input difference data, the x- and y- components of the vectors
    kwargs : dict, optional
        variable-specific options, See Notes

    Notes
    -----
    kwargs expected to be a variable-specific section,
    possibly provided by an ADF Variable Defaults YAML file.
    Currently it is inspected for:
    - `central_longitude`
    - `var_name`
    - `case_name`
    - `baseline`
    - `tiString`
    - `tiFontSize`
    - `units`

    _Note_ The title string constructed by kwargs appears to not be used.
    """

    # specify the central longitude for the plot:
    cent_long = kwargs.get('central_longitude', 180)

    # generate projection:
    proj = ccrs.PlateCarree(central_longitude=cent_long)
    lat = umdlfld_nowrap['lat']
    wgt = np.cos(np.radians(lat))

    # add cyclic longitude:
    umdlfld, lon = add_cyclic_point(umdlfld_nowrap, coord=umdlfld_nowrap['lon'])
    vmdlfld, _   = add_cyclic_point(vmdlfld_nowrap, coord=vmdlfld_nowrap['lon'])
    uobsfld, _   = add_cyclic_point(uobsfld_nowrap, coord=uobsfld_nowrap['lon'])
    vobsfld, _   = add_cyclic_point(vobsfld_nowrap, coord=vobsfld_nowrap['lon'])
    udiffld, _   = add_cyclic_point(udiffld_nowrap, coord=udiffld_nowrap['lon'])
    vdiffld, _   = add_cyclic_point(vdiffld_nowrap, coord=vdiffld_nowrap['lon'])

    # create mesh for plots:
    lons, lats = np.meshgrid(lon, lat)

    # create figure:
    fig = plt.figure(figsize=(14,10))

    # LAYOUT WITH GRIDSPEC
    gs = mpl.gridspec.GridSpec(3, 6, wspace=0.5, hspace=0.0)
    gs.tight_layout(fig)
    ax1 = plt.subplot(gs[0:2, :3], projection=proj)
    ax2 = plt.subplot(gs[0:2, 3:], projection=proj)
    ax3 = plt.subplot(gs[2, 1:5], projection=proj)
    ax = [ax1,ax2,ax3]

    # formatting for tick labels
    lon_formatter = LongitudeFormatter(number_format='0.0f',
                                        degree_symbol='',
                                        dateline_direction_label=False)
    lat_formatter = LatitudeFormatter(number_format='0.0f',
                                        degree_symbol='')

    # too many vectors to see well, so prune by striding through data:
    skip=(slice(None,None,5),slice(None,None,8))

    title_string = "Missing title!"
    title_string_base = title_string
    if "var_name" in kwargs:
        var_name = kwargs["var_name"]
    else:
        var_name = "missing VAR name"
    #End if

    if "case_name" in kwargs:
        case_name = kwargs["case_name"]
        if plev:
            title_string = f"{case_name} {var_name} [{plev} hPa]"
        else:
            title_string = f"{case_name} {var_name}"
        #End if
    #End if
    if "baseline" in kwargs:
        data_name = kwargs["baseline"]
        if plev:
            title_string_base = f"{data_name} {var_name} [{plev} hPa]"
        else:
            title_string_base = f"{data_name} {var_name}"
        #End if
    #End if

    # Calculate vector magnitudes.
    # Please note that the difference field needs
    # to be calculated from the model and obs fields
    # in order to get the correct sign:
    mdl_mag_ma  = np.sqrt(umdlfld**2 + vmdlfld**2)
    obs_mag_ma  = np.sqrt(uobsfld**2 + vobsfld**2)

    #Convert vector magnitudes to xarray DataArrays:
    mdl_mag  = xr.DataArray(mdl_mag_ma)
    obs_mag  = xr.DataArray(obs_mag_ma)
    diff_mag = mdl_mag - obs_mag

    # Get difference limits, in order to plot the correct range:
    min_diff_val = np.min(diff_mag)
    max_diff_val = np.max(diff_mag)

    # Color normalization for difference
    if (min_diff_val < 0) and (0 < max_diff_val):
        normdiff = mpl.colors.TwoSlopeNorm(vmin=min_diff_val, vmax=max_diff_val, vcenter=0.0)
    else:
        normdiff = mpl.colors.Normalize(vmin=min_diff_val, vmax=max_diff_val)
    #End if

    # Generate vector plot:
    #  - contourf to show magnitude w/ colorbar
    #  - vectors (colored or not) to show flow --> subjective (?) choice for how to thin out vectors to be legible
    img1 = ax1.contourf(lons, lats, mdl_mag, cmap='Greys', transform=ccrs.PlateCarree())
    ax1.quiver(lons[skip], lats[skip], umdlfld[skip], vmdlfld[skip], mdl_mag.values[skip], transform=ccrs.PlateCarree(), cmap='Reds')

    img2 = ax2.contourf(lons, lats, obs_mag, cmap='Greys', transform=ccrs.PlateCarree())
    ax2.quiver(lons[skip], lats[skip], uobsfld[skip], vobsfld[skip], obs_mag.values[skip], transform=ccrs.PlateCarree(), cmap='Reds')

    # We should think about how to do plot customization and defaults.
    # Here I'll just pop off a few custom ones, and then pass the rest into mpl.
    if 'tiString' in kwargs:
        tiString = kwargs.pop("tiString")
    else:
        tiString = ''
    if 'tiFontSize' in kwargs:
        tiFontSize = kwargs.pop('tiFontSize')
    else:
        tiFontSize = 8
    #End if

    #Set Main title for subplots:
    st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=18)
    st.set_y(0.85)

    #Set plot titles
    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"
    ax[0].set_title(case_title, loc='left', fontsize=tiFontSize)

    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)

    #Set stats: area_avg
    ax[0].set_title(f"Mean: {mdl_mag.weighted(wgt).mean().item():5.2f}\nMax: {mdl_mag.max():5.2f}\nMin: {mdl_mag.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    ax[1].set_title(f"Mean: {obs_mag.weighted(wgt).mean().item():5.2f}\nMax: {obs_mag.max():5.2f}\nMin: {mdl_mag.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    ax[-1].set_title(f"Mean: {diff_mag.weighted(wgt).mean().item():5.2f}\nMax: {diff_mag.max():5.2f}\nMin: {mdl_mag.min():5.2f}", loc='right',
                       fontsize=tiFontSize)

    # set rmse title:
    ax[-1].set_title(f"RMSE: ", fontsize=tiFontSize)
    ax[-1].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=tiFontSize)

    if "units" in kwargs:
        ax[1].set_ylabel(f"[{kwargs['units']}]")
        ax[-1].set_ylabel(f"[{kwargs['units']}]")
    #End if

    # Add cosmetic plot features:
    for a in ax:
        a.spines['geo'].set_linewidth(1.5) #cartopy's recommended method
        a.coastlines()
        a.set_xticks(np.linspace(-180, 120, 6), crs=ccrs.PlateCarree())
        a.set_yticks(np.linspace(-90, 90, 7), crs=ccrs.PlateCarree())
        a.tick_params('both', length=5, width=1.5, which='major')
        a.tick_params('both', length=5, width=1.5, which='minor')
        a.xaxis.set_major_formatter(lon_formatter)
        a.yaxis.set_major_formatter(lat_formatter)
    #End for

    # Add colorbar to vector plot:
    cb_c2_ax = inset_axes(ax2,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 100%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0, 1, 1),
                   bbox_transform=ax2.transAxes,
                   borderpad=0,
                   )
    fig.colorbar(img2, cax=cb_c2_ax)

    # Plot vector differences:
    img3 = ax3.contourf(lons, lats, diff_mag, transform=ccrs.PlateCarree(), norm=normdiff, cmap='PuOr', alpha=0.5)
    ax3.quiver(lons[skip], lats[skip], udiffld[skip], vdiffld[skip], transform=ccrs.PlateCarree())

    # Add color bar to difference plot:
    cb_d_ax = inset_axes(ax3,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 100%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0, 1, 1),
                   bbox_transform=ax3.transAxes,
                   borderpad=0
                   )
    fig.colorbar(img3, cax=cb_d_ax)

    # Write final figure to file
    fig.savefig(wks, bbox_inches='tight', dpi=300)

    #Close plots:
    plt.close()


#######

def plot_map_and_save(wks, case_nickname, base_nickname,
                      case_climo_yrs, baseline_climo_yrs,
                      mdlfld, obsfld, diffld, obs=False, **kwargs):
    """This plots mdlfld, obsfld, diffld in a 3-row panel plot of maps.

    Parameters
    ----------
    wks : str or Path
        output file path
    case_nickname : str
        short name for case
    base_nickname : str
        short name for base case
    case_climo_yrs : list
        list of years in case climatology, used for annotation
    baseline_climo_yrs : list
        list of years in base case climatology, used for annotation
    mdlfld : xarray.DataArray
        input data for case
    obsfld : xarray.DataArray
        input data for base case
    diffld : xarray.DataArray
        input difference data
    kwargs : dict, optional
        variable-specific options, See Notes

    Notes
    -----
    kwargs expected to be a variable-specific section,
    possibly provided by an ADF Variable Defaults YAML file.
    Currently it is inspected for:
    - colormap -> str, name of matplotlib colormap
    - contour_levels -> list of explict values or a tuple: (min, max, step)
    - diff_colormap
    - diff_contour_levels
    - tiString -> str, Title String
    - tiFontSize -> int, Title Font Size
    - mpl -> dict, This should be any matplotlib kwargs that should be passed along. Keep reading:
        + Organize these by the mpl function. In this function (`plot_map_and_save`)
          we will check for an entry called `subplots`, `contourf`, and `colorbar`. So the YAML might looks something like:
          ```
           mpl:
             subplots:
               figsize: (3, 9)
             contourf:
               levels: 15
               cmap: Blues
             colorbar:
               shrink: 0.4
          ```
        + This is experimental, and if you find yourself doing much with this, you probably should write a new plotting script that does not rely on this module.
    When these are not provided, colormap is set to 'coolwarm' and limits/levels are set by data range.
    """

    #nice formatting for tick labels
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    # preprocess
    # - assume all three fields have same lat/lon
    lat = obsfld['lat']
    wgt = np.cos(np.radians(lat))
    mwrap, lon = add_cyclic_point(mdlfld, coord=mdlfld['lon'])
    owrap, _ = add_cyclic_point(obsfld, coord=obsfld['lon'])
    dwrap, _ = add_cyclic_point(diffld, coord=diffld['lon'])
    wrap_fields = (mwrap, owrap, dwrap)
    # mesh for plots:
    lons, lats = np.meshgrid(lon, lat)
    # Note: using wrapped data makes spurious lines across plot (maybe coordinate dependent)
    lon2, lat2 = np.meshgrid(mdlfld['lon'], mdlfld['lat'])

    # get statistics (from non-wrapped)
    fields = (mdlfld, obsfld, diffld)
    area_avg = [global_average(x, wgt) for x in fields]

    d_rmse = wgt_rmse(mdlfld, obsfld, wgt)  # correct weighted RMSE for (lat,lon) fields.

    # We should think about how to do plot customization and defaults.
    # Here I'll just pop off a few custom ones, and then pass the rest into mpl.
    if 'tiString' in kwargs:
        tiString = kwargs.pop("tiString")
    else:
        tiString = ''
    #End if

    if 'tiFontSize' in kwargs:
        tiFontSize = kwargs.pop('tiFontSize')
    else:
        tiFontSize = 8
    #End if

    # generate dictionary of contour plot settings:
    cp_info = prep_contour_plot(mdlfld, obsfld, diffld, **kwargs)

    # specify the central longitude for the plot
    central_longitude = kwargs.get('central_longitude', 180)

    # create figure object
    fig = plt.figure(figsize=(14,10))

    # LAYOUT WITH GRIDSPEC
    gs = mpl.gridspec.GridSpec(3, 6, wspace=0.5,hspace=0.0) # 2 rows, 4 columns, but each map will take up 2 columns
    #gs.tight_layout(fig)
    proj = ccrs.PlateCarree(central_longitude=central_longitude)
    ax1 = plt.subplot(gs[0:2, :3], projection=proj, **cp_info['subplots_opt'])
    ax2 = plt.subplot(gs[0:2, 3:], projection=proj, **cp_info['subplots_opt'])
    ax3 = plt.subplot(gs[2, 1:5], projection=proj,  **cp_info['subplots_opt'])
    ax = [ax1,ax2,ax3]

    img = [] # contour plots
    cs = []  # contour lines
    cb = []  # color bars

    # formatting for tick labels
    lon_formatter = LongitudeFormatter(number_format='0.0f',
                                        degree_symbol='',
                                        dateline_direction_label=False)
    lat_formatter = LatitudeFormatter(number_format='0.0f',
                                        degree_symbol='')

    for i, a in enumerate(wrap_fields):

        if i == len(wrap_fields)-1:
            levels = cp_info['levelsdiff']
            cmap = cp_info['cmapdiff']
            norm = cp_info['normdiff']
        else:
            levels = cp_info['levels1']
            cmap = cp_info['cmap1']
            norm = cp_info['norm1']

        levs = np.unique(np.array(levels))
        if len(levs) < 2:
            img.append(ax[i].contourf(lons,lats,a,colors="w",transform=ccrs.PlateCarree()))
            ax[i].text(0.4, 0.4, empty_message, transform=ax[i].transAxes, bbox=props)
        else:
            img.append(ax[i].contourf(lons, lats, a, levels=levels, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), **cp_info['contourf_opt']))
        #End if
        ax[i].set_title("AVG: {0:.3f}".format(area_avg[i]), loc='right', fontsize=11)

        # add contour lines <- Unused for now -JN
        # TODO: add an option to turn this on -BM
        #cs.append(ax[i].contour(lon2, lat2, fields[i], transform=ccrs.PlateCarree(), colors='k', linewidths=1))
        #ax[i].clabel(cs[i], cs[i].levels, inline=True, fontsize=tiFontSize-2, fmt='%1.1f')
        #ax[i].text( 10, -140, "CONTOUR FROM {} to {} by {}".format(min(cs[i].levels), max(cs[i].levels), cs[i].levels[1]-cs[i].levels[0]),
        #bbox=dict(facecolor='none', edgecolor='black'), fontsize=tiFontSize-2)

    st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=18)
    st.set_y(0.85)

    #Set plot titles
    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"
    ax[0].set_title(case_title, loc='left', fontsize=tiFontSize)

    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)

    #Set stats: area_avg
    ax[0].set_title(f"Mean: {mdlfld.weighted(wgt).mean().item():5.2f}\nMax: {mdlfld.max():5.2f}\nMin: {mdlfld.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    ax[1].set_title(f"Mean: {obsfld.weighted(wgt).mean().item():5.2f}\nMax: {obsfld.max():5.2f}\nMin: {obsfld.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    ax[-1].set_title(f"Mean: {diffld.weighted(wgt).mean().item():5.2f}\nMax: {diffld.max():5.2f}\nMin: {diffld.min():5.2f}", loc='right',
                       fontsize=tiFontSize)

    # set rmse title:
    ax[-1].set_title(f"RMSE: {d_rmse:.3f}", fontsize=tiFontSize)
    ax[-1].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=tiFontSize)

    if "units" in kwargs:
        ax[1].set_ylabel(f"[{kwargs['units']}]")
        ax[-1].set_ylabel(f"[{kwargs['units']}]")
    else:
        ax[1].set_ylabel(f"{mdlfld.units}")
        ax[-1].set_ylabel(f"{mdlfld.units}")

    for a in ax:
        a.spines['geo'].set_linewidth(1.5) #cartopy's recommended method
        a.coastlines()
        a.set_xticks(np.linspace(-180, 120, 6), crs=ccrs.PlateCarree())
        a.set_yticks(np.linspace(-90, 90, 7), crs=ccrs.PlateCarree())
        a.tick_params('both', length=5, width=1.5, which='major')
        a.tick_params('both', length=5, width=1.5, which='minor')
        a.xaxis.set_major_formatter(lon_formatter)
        a.yaxis.set_major_formatter(lat_formatter)

    # __COLORBARS__
    cb_mean_ax = inset_axes(ax2,
                    width="5%",  # width = 5% of parent_bbox width
                    height="100%",  # height : 100%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0, 1, 1),
                    bbox_transform=ax2.transAxes,
                    borderpad=0,
                    )
    fig.colorbar(img[1], cax=cb_mean_ax, **cp_info['colorbar_opt'])

    cb_diff_ax = inset_axes(ax3,
                    width="5%",  # width = 5% of parent_bbox width
                    height="100%",  # height : 100%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0, 1, 1),
                    bbox_transform=ax3.transAxes,
                    borderpad=0,
                    )
    fig.colorbar(img[2], cax=cb_diff_ax, **cp_info['colorbar_opt'])

    # Write final figure to file
    fig.savefig(wks, bbox_inches='tight', dpi=300)

    #Close plots:
    plt.close()

#
#  -- vertical interpolation code --
#

def pres_from_hybrid(psfc, hya, hyb, p0=100000.):
    """Calculates pressure field

    pressure derived with the formula:
    ```p = a(k)*p0 + b(k)*ps```

    Parameters
    ----------
    psfc
        surface pressure
    hya, hyb
        hybrid-sigma A and B coefficients
    p0 : optional
        reference pressure, defaults to 100000 Pa

    Returns
    -------
    pressure, size is same as `psfc` with `len(hya)` levels
    """
    return hya*p0 + hyb*psfc

#####

def vert_remap(x_mdl, p_mdl, plev):
    """Apply simple 1-d interpolation to a field

    Parameters
    ----------
    x_mdl : xarray.DataArray or numpy.ndarray
        input data
    p_mdl : xarray.DataArray or numpy.ndarray
        pressure field, same shape as `x_mdl`
    plev : xarray.DataArray or numpy.ndarray
        the new pressures

    Returns
    -------
    output
        `x_mdl` interpolated to `plev`

    Notes
    -----
    Interpolation done in log pressure
    """

    #Determine array shape of output array:
    out_shape = (plev.shape[0], x_mdl.shape[1])

    #Initialize interpolated output numpy array:
    output = np.full(out_shape, np.nan)

    #Perform 1-D interpolation in log-space:
    for i in range(out_shape[1]):
        output[:,i] = np.interp(np.log(plev), np.log(p_mdl[:,i]), x_mdl[:,i])
    #End for

    #Return interpolated output:
    return output

#####

def lev_to_plev(data, ps, hyam, hybm, P0=100000., new_levels=None,
                convert_to_mb=False):
    """Interpolate model hybrid levels to specified pressure levels.

    Parameters
    ----------
    data :
    ps :
        surface pressure
    hyam, hybm :
        hybrid-sigma A and B coefficients
    P0 : float, optional
        reference pressure, defaults to 100000 Pa
    new_levels : numpy.ndarray, optional
        1-D array containing pressure levels in Pascals (Pa).
        If not specified, then the levels will be set
        to the GeoCAT defaults, which are (in hPa):
        `1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50,
        30, 20, 10, 7, 5, 3, 2, 1`
    convert_to_mb : bool, optional
        If True, then vertical (lev) dimension will have
        values of mb/hPa, otherwise the units are Pa.

    Returns
    -------
    data_interp_rename
        data interpolated to new pressure levels

    Notes
    -----
    The function `interp_hybrid_to_pressure` used here is dask-enabled,
    and so can potentially be sped-up via the use of a DASK cluster.
    """

    #Temporary print statement to notify users to ignore warning messages.
    #This should be replaced by a debug-log stdout filter at some point:
    print("Please ignore the interpolation warnings that follow!")

    #Apply GeoCAT hybrid->pressure interpolation:
    if new_levels is not None:
        data_interp = gcomp.interpolation.interp_hybrid_to_pressure(data, ps,
                                                                    hyam,
                                                                    hybm,
                                                                    p0=P0,
                                                                    new_levels=new_levels
                                                                   )
    else:
        data_interp = gcomp.interpolation.interp_hybrid_to_pressure(data, ps,
                                                                    hyam,
                                                                    hybm,
                                                                    p0=P0
                                                                   )

    # data_interp may contain a dask array, which can cause
    # trouble downstream with numpy functions, so call compute() here.
    if hasattr(data_interp, "compute"):
        data_interp = data_interp.compute()

    #Rename vertical dimension back to "lev" in order to work with
    #the ADF plotting functions:
    data_interp_rename = data_interp.rename({"plev": "lev"})

    #Convert vertical dimension to mb/hPa, if requested:
    if convert_to_mb:
        data_interp_rename["lev"] = data_interp_rename["lev"] / 100.0

    return data_interp_rename

#####

def pmid_to_plev(data, pmid, new_levels=None, convert_to_mb=False):
    """Interpolate data from hybrid-sigma levels to isobaric levels.

    Parameters
    ----------
    data : xarray.DataArray
        field with a 'lev' coordinate
    pmid : xarray.DataArray
        the pressure field (Pa), same shape as `data`
    new_levels : optional
        the output pressure levels (Pa), defaults to standard levels
    convert_to_mb : bool, optional
        flag to convert output to mb (i.e., hPa), defaults to False

    Returns
    -------
    output : xarray.DataArray
        `data` interpolated onto `new_levels`
    """

    # determine pressure levels to interpolate to:
    if new_levels is None:
        pnew = 100.0 * np.array([1000, 925, 850, 700, 500, 400,
                                 300, 250, 200, 150, 100, 70, 50,
                                 30, 20, 10, 7, 5, 3, 2, 1])  # mandatory levels, converted to Pa
    else:
        pnew = new_levels
    #End if

    # save name of DataArray:
    data_name = data.name

    # reshape data and pressure assuming "lev" is the name of the coordinate
    zdims = [i for i in data.dims if i != 'lev']
    dstack = data.stack(z=zdims)
    pstack = pmid.stack(z=zdims)
    output = vert_remap(dstack.values, pstack.values, pnew)
    output = xr.DataArray(output, name=data_name, dims=("lev", "z"),
                          coords={"lev":pnew, "z":pstack['z']})
    output = output.unstack()

    # convert vertical dimension to mb/hPa, if requested:
    if convert_to_mb:
        output["lev"] = output["lev"] / 100.0
    #End if

    #Return interpolated output:
    return output

#
#  -- zonal & meridional mean code --
#

def zonal_mean_xr(fld):
    """Average over all dimensions except `lev` and `lat`."""
    if isinstance(fld, xr.DataArray):
        d = fld.dims
        davgovr = [dim for dim in d if dim not in ('lev','lat')]
    else:
        raise IOError("zonal_mean_xr requires Xarray DataArray input.")
    return fld.mean(dim=davgovr)


def validate_dims(fld, list_of_dims):
    """Check if specified dimensions are in a DataArray.

    Parameters
    ----------
    fld : xarray.DataArray
        field to check for named dimensions
    list_of_dims : list
        list of strings that specifiy the dimensions to check for

    Returns
    -------
    dict
        dict with keys that are "has_{x}" where x is the name from
        `list_of_dims` and values that are boolean

    """
    if not isinstance(list_of_dims, list):
        list_of_dims = list(list_of_dims)
    return { "_".join(["has",f"{v}"]):(v in fld.dims) for v in list_of_dims}


def lat_lon_validate_dims(fld):
    """Check if input field has lat and lon.

    Parameters
    ----------
    fld : xarray.DataArray
        data with named dimensions

    Returns
    -------
    bool
        True if lat and lon are both dimensions, False otherwise.

    See Also
    --------
    validate_dims
    """
    # note: we can only handle variables that reduce to (lat,lon)
    if len(fld.dims) > 3:
        return False
    validate = validate_dims(fld, ['lat','lon'])
    if not all(validate.values()):
        return  False
    else:
        return True


def zm_validate_dims(fld):
    """Check for dimensions for zonal average.

    Looks for dimensions called 'lev' and 'lat'.


    Parameters
    ----------
    fld : xarray.DataArray
        field to check for lat and/or lev dimensions
    Returns
    -------
    tuple
        (has_lat, has_lev) each are bool
    None
        If 'lat' is not in dimensions, returns None.
    """
    # note: we can only handle variables that reduce to (lev, lat) or (lat,)
    if len(fld.dims) > 4:
        print(f"Sorry, too many dimensions: {fld.dims}")
        return None
    validate = validate_dims(fld, ['lev','lat'])
    has_lev, has_lat = validate['has_lev'], validate['has_lat']
    if not has_lat:
        return None
    else:
        return has_lat, has_lev

def _plot_line(axobject, xdata, ydata, color, **kwargs):
    """Create a generic line plot and check for some ways to annotate."""

    if color != None:
        axobject.plot(xdata, ydata, c=color, **kwargs)
    else:
        axobject.plot(xdata, ydata, **kwargs)

    #Set Y-axis label:
    if hasattr(ydata, "units"):
        axobject.set_ylabel("[{units}]".format(units=getattr(ydata,"units")))
    elif "units" in kwargs:
        axobject.set_ylabel("[{units}]".format(kwargs["units"]))
    #End if

    return axobject

def _meridional_plot_line(ax, lon, data, color, **kwargs):
    """Create line plot with longitude as the X-axis."""

    ax = _plot_line(ax, lon, data, color, **kwargs)
    ax.set_xlim([lon.min(), lon.max()])
    #
    # annotate
    #
    ax.set_xlabel("LONGITUDE")
    if hasattr(data, "units"):
        ax.set_ylabel("{units}".format(units=getattr(data,"units")))
    elif "units" in kwargs:
        ax.set_ylabel("{units}".format(kwargs["units"]))
    return ax

def _zonal_plot_line(ax, lat, data, color, **kwargs):
    """Create line plot with latitude as the X-axis."""
    ax = _plot_line(ax, lat, data, color, **kwargs)
    ax.set_xlim([max([lat.min(), -90.]), min([lat.max(), 90.])])
    #
    # annotate
    #
    ax.set_xlabel("LATITUDE")
    if hasattr(data, "units"):
        ax.set_ylabel("{units}".format(units=getattr(data,"units")))
    elif "units" in kwargs:
        ax.set_ylabel("{units}".format(kwargs["units"]))
    return ax

def _zonal_plot_preslat(ax, lat, lev, data, **kwargs):
    """Create plot with latitude as the X-axis, and pressure as the Y-axis."""
    mlev, mlat = np.meshgrid(lev, lat)
    if 'cmap' in kwargs:
        cmap = kwargs.pop('cmap')
    else:
        cmap = 'Spectral_r'

    img = ax.contourf(mlat, mlev, data.transpose('lat', 'lev'), cmap=cmap, **kwargs)

    minor_locator = mpl.ticker.FixedLocator(lev)
    ax.yaxis.set_minor_locator(minor_locator)
    ax.tick_params(which='minor', length=4, color='r')
    ax.set_ylim([np.max(lev), np.min(lev)])
    return img, ax


def _meridional_plot_preslon(ax, lon, lev, data, **kwargs):
    """Create plot with longitude as the X-axis, and pressure as the Y-axis."""

    mlev, mlon = np.meshgrid(lev, lon)
    if 'cmap' in kwargs:
        cmap = kwargs.pop('cmap')
    else:
        cmap = 'Spectral_r'

    img = ax.contourf(mlon, mlev, data.transpose('lon', 'lev'), cmap=cmap, **kwargs)

    minor_locator = mpl.ticker.FixedLocator(lev)
    ax.yaxis.set_minor_locator(minor_locator)
    ax.tick_params(which='minor', length=4, color='r')
    ax.set_ylim([np.max(lev), np.min(lev)])
    return img, ax

def zonal_plot(lat, data, ax=None, color=None, **kwargs):
    """Make zonal plot

    Determine which kind of zonal plot is needed based
    on the input variable's dimensions.

    Parameters
    ----------
    lat
        latitude
    data
        input data
    ax : Axes, optional
        axes object to use
    color : str or mpl color specification
        color for the curve
    kwargs : dict, optional
        plotting options

    Notes
    -----
    Checks if there is a `lev` dimension to determine if
    it is a lat-pres plot or a line plot.
    """
    if ax is None:
        ax = plt.gca()
    if 'lev' in data.dims:
        img, ax = _zonal_plot_preslat(ax, lat, data['lev'], data, **kwargs)
        return img, ax
    else:
        ax = _zonal_plot_line(ax, lat, data, color, **kwargs)
        return ax

def meridional_plot(lon, data, ax=None, color=None, **kwargs):
    """Make meridional plot

    Determine which kind of meridional plot is needed based
    on the input variable's dimensions.


    Parameters
    ----------
    lon
        longitude
    data
        input data
    ax : Axes, optional
        axes object to use
    color : str or mpl color specification
        color for the curve
    kwargs : dict, optional
        plotting options

    Notes
    -----
    Checks if there is a `lev` dimension to determine if
    it is a lon-pres plot or a line plot.
    """
    if ax is None:
        ax = plt.gca()
    if 'lev' in data.dims:
        img, ax = _meridional_plot_preslon(ax, lon, data['lev'], data, **kwargs)
        return img, ax
    else:
        ax = _meridional_plot_line(ax, lon,  data, color, **kwargs)
        return ax

def prep_contour_plot(adata, bdata, diffdata, **kwargs):
    """Preparation for making contour plots.

    Prepares for making contour plots of adata, bdata, and diffdata, which is
    presumably the difference between adata and bdata.
    - set colormap from kwargs or defaults to coolwarm
    - set contour levels from kwargs or 12 evenly spaced levels to span the data
    - normalize colors based on specified contour levels or data range
    - set option for linear or log pressure when applicable
    - similar settings for difference, defaults to symmetric about zero
    - separates Matplotlib kwargs into their own dicts

    Parameters
    ----------
    adata, bdata, diffdata
        the data to be plotted
    kwargs : dict, optional
        plotting options

    Returns
    -------
    dict
        a dict with the following:
        - 'subplots_opt': mpl kwargs for subplots
        - 'contourf_opt': mpl kwargs for contourf
        - 'colorbar_opt': mpl kwargs for colorbar
        - 'diff_colorbar_opt' : mpl kwargs for difference colorbar
        - 'normdiff': color normalization for difference panel
        - 'cmapdiff': colormap for difference panel
        - 'levelsdiff': contour levels for difference panel
        - 'cmap1': color map for a and b panels
        - 'norm1': color normalization for a and b panels
        - 'levels1' : contour levels for a and b panels
        - 'plot_log_p' : true/false whether to plot log(pressure) axis
    """
    # determine levels & color normalization:
    # Replace inf values with NaN
    #adata = adata.where(~np.isinf(adata))
    #bdata = bdata.where(~np.isinf(bdata))
    minval = np.min([np.nanmin(adata), np.nanmin(bdata)])
    maxval = np.max([np.nanmax(adata), np.nanmax(bdata)])

    # determine norm to use (deprecate this once minimum MPL version is high enough)
    normfunc, mplv = use_this_norm()

    if 'colormap' in kwargs:
        cmap1 = kwargs['colormap']
    else:
        cmap1 = 'coolwarm'
    #End if

    if 'contour_levels' in kwargs:
        #Make these floats in case the contour levels are in scientific notation
        levels1 = [float(x) for x in kwargs['contour_levels']]
        norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
    elif 'contour_levels_range' in kwargs:
        assert len(kwargs['contour_levels_range']) == 3, \
        "contour_levels_range must have exactly three entries: min, max, step"
        #Make these floats in case the contour levels are in scientific notation
        lev_range = [float(x) for x in kwargs['contour_levels_range']]
        levels1 = np.arange(*lev_range)
        norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
    else:
        levels1 = np.linspace(minval, maxval, 12)
        norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    #End if

    if ('non_linear_levels' in kwargs) and (kwargs['non_linear_levels']):
        cmap_obj = cm.get_cmap(cmap1)
        norm1 = mpl.colors.BoundaryNorm(levels1, cmap_obj.N)

    #Check if the minval and maxval are actually different.  If not,
    #then set "levels1" to be an empty list, which will cause the
    #plotting scripts to add a label instead of trying to plot a variable
    #with no contours:
    if minval == maxval:
        levels1 = []
    #End if

    if ('colormap' not in kwargs) and ('contour_levels' not in kwargs):
        if ((minval < 0) and (0 < maxval)) and mplv > 2:
            norm1 = normfunc(vmin=minval, vmax=maxval, vcenter=0.0)
        else:
            norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)
        #End if
    #End if

# Difference options -- Check in kwargs for colormap and levels
    if "diff_colormap" in kwargs:
        cmapdiff = kwargs["diff_colormap"]
    else:
        cmapdiff = 'coolwarm'
    #End if

    if "diff_contour_levels" in kwargs:
        #Make these floats in case the contour levels are in scientific notation
        levelsdiff = [float(x) for x in kwargs['diff_contour_range']]
    elif "diff_contour_range" in kwargs:
        assert len(kwargs['diff_contour_range']) == 3, \
        "diff_contour_range must have exactly three entries: min, max, step"
        #Make these floats in case the contour levels are in scientific notation
        lev_range = [float(x) for x in kwargs['diff_contour_range']]
        levelsdiff = np.arange(*lev_range)
    else:
        # set a symmetric color bar for diff:
        absmaxdif = np.max(np.abs(diffdata))
        # set levels for difference plot:
        levelsdiff = np.linspace(-1*absmaxdif, absmaxdif, 12)
    #End if

    if "plot_log_pressure" in kwargs:
        plot_log_p = kwargs["plot_log_pressure"]
    else:
        plot_log_p = False

    # color normalization for difference
    if ((np.min(levelsdiff) < 0) and (0 < np.max(levelsdiff))) and mplv > 2:
        normdiff = normfunc(vmin=np.min(levelsdiff), vmax=np.max(levelsdiff), vcenter=0.0)
    else:
        normdiff = mpl.colors.Normalize(vmin=np.min(levelsdiff), vmax=np.max(levelsdiff))
    #End if

    subplots_opt = {}
    contourf_opt = {}
    colorbar_opt = {}
    diff_colorbar_opt = {}

    # extract any MPL kwargs that should be passed on:
    if 'mpl' in kwargs:
        subplots_opt.update(kwargs['mpl'].get('subplots',{}))
        contourf_opt.update(kwargs['mpl'].get('contourf',{}))
        colorbar_opt.update(kwargs['mpl'].get('colorbar',{}))
        diff_colorbar_opt.update(kwargs['mpl'].get('diff_colorbar',{}))
    #End if
    return {'subplots_opt': subplots_opt,
            'contourf_opt': contourf_opt,
            'colorbar_opt': colorbar_opt,
            'diff_colorbar_opt': diff_colorbar_opt,
            'normdiff': normdiff,
            'cmapdiff': cmapdiff,
            'levelsdiff': levelsdiff,
            'cmap1': cmap1,
            'norm1': norm1,
            'levels1': levels1,
            'plot_log_p': plot_log_p
            }


def plot_zonal_mean_and_save(wks, case_nickname, base_nickname,
                             case_climo_yrs, baseline_climo_yrs,
                             adata, bdata, has_lev, log_p, obs=False, **kwargs):

    """This is the default zonal mean plot

    Parameters
    ----------
    adata : data to plot ([lev], lat, [lon]).
            The vertical coordinate (lev) must be pressure levels.
    bdata : baseline or observations to plot adata against.

        - For 2-d variables (reduced to (lat,)):
          + 2 panels: (top) zonal mean, (bottom) difference
        - For 3-D variables (reduced to (lev,lat)):
          + 3 panels: (top) zonal mean adata, (middle) zonal mean bdata, (bottom) difference
          + pcolormesh/contour plot
    kwargs -> optional dictionary of plotting options
             ** Expecting this to be variable-specific section, possibly provided by ADF Variable Defaults YAML file.**
    - colormap -> str, name of matplotlib colormap
    - contour_levels -> list of explict values or a tuple: (min, max, step)
    - diff_colormap
    - diff_contour_levels
    - tiString -> str, Title String
    - tiFontSize -> int, Title Font Size
    - mpl -> dict, This should be any matplotlib kwargs that should be passed along. Keep reading:
        + Organize these by the mpl function. In this function (`plot_map_and_save`)
          we will check for an entry called `subplots`, `contourf`, and `colorbar`. So the YAML might looks something like:
          ```
           mpl:
             subplots:
               figsize: (3, 9)
             contourf:
               levels: 15
               cmap: Blues
             colorbar:
               shrink: 0.4
          ```
    """

    # style the plot:
    # We should think about how to do plot customization and defaults.
    # Here I'll just pop off a few custom ones, and then pass the rest into mpl.
    if 'tiFontSize' in kwargs:
        tiFontSize = kwargs.pop('tiFontSize')
    else:
        tiFontSize = 8
    #End if

    #Set plot titles
    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"

    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"
    if has_lev:

        # calculate zonal average:
        azm = zonal_mean_xr(adata)
        bzm = zonal_mean_xr(bdata)

        # calculate difference:
        diff = azm - bzm

        # generate dictionary of contour plot settings:
        cp_info = prep_contour_plot(azm, bzm, diff, **kwargs)

        # Generate zonal plot:
        fig, ax = plt.subplots(figsize=(10,10),nrows=3, constrained_layout=True, sharex=True, sharey=True,**cp_info['subplots_opt'])
        levs = np.unique(np.array(cp_info['levels1']))

        levs_diff = np.unique(np.array(cp_info['levelsdiff']))


        if len(levs) < 2:
            img0, ax[0] = zonal_plot(adata['lat'], azm, ax=ax[0])
            ax[0].text(0.4, 0.4, empty_message, transform=ax[0].transAxes, bbox=props)
            img1, ax[1] = zonal_plot(bdata['lat'], bzm, ax=ax[1])
            ax[1].text(0.4, 0.4, empty_message, transform=ax[1].transAxes, bbox=props)
        else:
            img0, ax[0] = zonal_plot(adata['lat'], azm, ax=ax[0], norm=cp_info['norm1'],cmap=cp_info['cmap1'],levels=cp_info['levels1'],**cp_info['contourf_opt'])
            img1, ax[1] = zonal_plot(bdata['lat'], bzm, ax=ax[1], norm=cp_info['norm1'],cmap=cp_info['cmap1'],levels=cp_info['levels1'],**cp_info['contourf_opt'])
            fig.colorbar(img0, ax=ax[0], location='right',**cp_info['colorbar_opt'])
            fig.colorbar(img1, ax=ax[1], location='right',**cp_info['colorbar_opt'])
        #End if

        if len(levs_diff) < 2:
            img2, ax[2] = zonal_plot(adata['lat'], diff, ax=ax[2])
            ax[2].text(0.4, 0.4, empty_message, transform=ax[2].transAxes, bbox=props)
        else:
            img2, ax[2] = zonal_plot(adata['lat'], diff, ax=ax[2], norm=cp_info['normdiff'],cmap=cp_info['cmapdiff'],levels=cp_info['levelsdiff'],**cp_info['contourf_opt'])
            fig.colorbar(img2, ax=ax[2], location='right',**cp_info['colorbar_opt'])

        ax[0].set_title(case_title, loc='left', fontsize=tiFontSize)
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)
        ax[2].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=tiFontSize)


        # style the plot:
        #Set Main title for subplots:
        st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=15)
        st.set_y(0.85)
        ax[-1].set_xlabel("LATITUDE")

        if log_p:
            [a.set_yscale("log") for a in ax]

        fig.text(-0.03, 0.5, 'PRESSURE [hPa]', va='center', rotation='vertical')
    else:
        line = Line2D([0], [0], label="$\mathbf{Test}:$"+f"{case_nickname} - years: {case_climo_yrs[0]}-{case_climo_yrs[-1]}",
                        color="#1f77b4") # #1f77b4 -> matplotlib standard blue

        line2 = Line2D([0], [0], label=base_title,
                        color="#ff7f0e") # #ff7f0e -> matplotlib standard orange

        azm = zonal_mean_xr(adata)
        bzm = zonal_mean_xr(bdata)
        diff = azm - bzm
        fig, ax = plt.subplots(nrows=2)
        ax = [ax[0],ax[1]]

        #Set Main title for subplots:
        st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=15)
        st.set_y(1.02)

        zonal_plot(adata['lat'], azm, ax=ax[0],color="#1f77b4") # #1f77b4 -> matplotlib standard blue
        zonal_plot(bdata['lat'], bzm, ax=ax[0],color="#ff7f0e") # #ff7f0e -> matplotlib standard orange

        fig.legend(handles=[line,line2],bbox_to_anchor=(-0.15, 0.87, 1.05, .102),loc="right",
                   borderaxespad=0.0,fontsize=6,frameon=False)

        zonal_plot(adata['lat'], diff, ax=ax[1], color="k")
        ax[1].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=10)

        for a in ax:
            try:
                a.label_outer()
            except:
                pass
            #End except
        #End for
    #End if

    #Write the figure to provided workspace/file:
    fig.savefig(wks, bbox_inches='tight', dpi=300)

    #Close plots:
    plt.close()



def plot_meridional_mean_and_save(wks, case_nickname, base_nickname,
                             case_climo_yrs, baseline_climo_yrs,
                             adata, bdata, has_lev, latbounds=None, obs=False, **kwargs):

    """Default meridional mean plot

    Parameters
    ----------
    wks :
        the figure object to plot in
    case_nickname : str
        short name of `adata` case, use for annotation
    base_nickname : str
        short name of `bdata` case, use for annotation
    case_climo_yrs : list
        years in the `adata` case, use for annotation
    baseline_climo_yrs : list:
        years in the `bdata` case, use for annotation
    adata : xarray.DataArray
        data to plot ([lev], [lat], lon).
        The vertical coordinate (lev) must be pressure levels.
    bdata : xarray.DataArray
        baseline or observations to plot adata against.
        It must have the same dimensions and vertical levels as adata.
    has_lev : bool
        whether lev dimension is present
    latbounds : numbers.Number or slice, optional
        indicates latitude bounds to average over
        if it is a number, assume symmetric about equator,
        otherwise expects `slice(south, north)`
        defaults to `slice(-5,5)`
    kwargs : dict, optional
        optional dictionary of plotting options, See Notes

    Notes
    -----

    - For 2-d variables (reduced to (lon,)):
        + 2 panels: (top) meridional mean, (bottom) difference
    - For 3-D variables (reduced to (lev,lon)):
        + 3 panels: (top) meridonal mean adata, (middle) meridional mean bdata, (bottom) difference
        + pcolormesh/contour plot

    - kwargs -> optional dictionary of plotting options
        ** Expecting this to be variable-specific section, possibly
        provided by ADF Variable Defaults YAML file.**
        - colormap             -> str, name of matplotlib colormap
        - contour_levels       -> list of explicit values or a tuple: (min, max, step)
        - diff_colormap        -> str, name of matplotlib colormap used for different plot
        - diff_contour_levels  -> list of explicit values or a tuple: (min, max, step)
        - tiString             -> str, Title String
        - tiFontSize           -> int, Title Font Size
        - mpl -> dict, This should be any matplotlib kwargs that should be passed along. Keep reading:
            + Organize these by the mpl function. In this function (`plot_meridional_mean_and_save`)
            we will check for an entry called `subplots`, `contourf`, and `colorbar`.
            So the YAML might looks something like:
            ```
            mpl:
                subplots:
                figsize: (3, 9)
                contourf:
                levels: 15
                cmap: Blues
                colorbar:
                shrink: 0.4
            ```
        """
    # apply averaging:
    import numbers  # built-in; just checking on the latbounds input
    if latbounds is None:
        latbounds = slice(-5, 5)
    elif isinstance(latbounds, numbers.Number):
        latbounds = slice(-1*np.absolute(latbounds), np.absolute(latbounds))
    elif not isinstance(latbounds, slice):  #If not a slice object, then quit this routine.
        print(f"ERROR: plot_meridonal_mean_and_save - received an invalid value for latbounds ({latbounds}). Must be a number or a slice.")
        return None
    #End if

    # style the plot:
    # We should think about how to do plot customization and defaults.
    # Here I'll just pop off a few custom ones, and then pass the rest into mpl.
    if 'tiFontSize' in kwargs:
        tiFontSize = kwargs.pop('tiFontSize')
    else:
        tiFontSize = 8
    #End if

    # possible that the data has time, but usually it won't
    if len(adata.dims) > 4:
        print(f"ERROR: plot_meridonal_mean_and_save - too many dimensions: {adata.dims}")
        return None

    if 'time' in adata.dims:
        adata = adata.mean(dim='time', keep_attrs=True)
    if 'lat' in adata.dims:
        latweight = np.cos(np.radians(adata.lat))
        adata = adata.weighted(latweight).mean(dim='lat', keep_attrs=True)
    if 'time' in bdata.dims:
        adata = bdata.mean(dim='time', keep_attrs=True)
    if 'lat' in bdata.dims:
        latweight = np.cos(np.radians(bdata.lat))
        bdata = bdata.weighted(latweight).mean(dim='lat', keep_attrs=True)
    # If there are other dimensions, they are still going to be there:
    if len(adata.dims) > 2:
        print(f"ERROR: plot_meridonal_mean_and_save - AFTER averaging, there are too many dimensions: {adata.dims}")
        return None

    diff = adata - bdata

    # plot-controlling parameters:
    xdim = 'lon' # the name used for the x-axis dimension
    pltfunc = meridional_plot  # the plotting function ... maybe we can generalize to get zonal/meridional into one function (?)

    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"

    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"

    if has_lev:
        # generate dictionary of contour plot settings:
        cp_info = prep_contour_plot(adata, bdata, diff, **kwargs)

        # generate plot objects:
        fig, ax = plt.subplots(figsize=(10,10),nrows=3, constrained_layout=True, sharex=True, sharey=True,**cp_info['subplots_opt'])
        levs = np.unique(np.array(cp_info['levels1']))
        levs_diff = np.unique(np.array(cp_info['levelsdiff']))


        if len(levs) < 2:
            img0, ax[0] = pltfunc(adata[xdim], adata, ax=ax[0])
            ax[0].text(0.4, 0.4, empty_message, transform=ax[0].transAxes, bbox=props)
            img1, ax[1] = pltfunc(bdata[xdim], bdata, ax=ax[1])
            ax[1].text(0.4, 0.4, empty_message, transform=ax[1].transAxes, bbox=props)
        else:
            img0, ax[0] = pltfunc(adata[xdim], adata, ax=ax[0], norm=cp_info['norm1'],cmap=cp_info['cmap1'],levels=cp_info['levels1'],**cp_info['contourf_opt'])
            img1, ax[1] = pltfunc(bdata[xdim], bdata, ax=ax[1], norm=cp_info['norm1'],cmap=cp_info['cmap1'],levels=cp_info['levels1'],**cp_info['contourf_opt'])
            cb0 = fig.colorbar(img0, ax=ax[0], location='right',**cp_info['colorbar_opt'])
            cb1 = fig.colorbar(img1, ax=ax[1], location='right',**cp_info['colorbar_opt'])
        #End if

        if len(levs_diff) < 2:
            img2, ax[2] = pltfunc(adata[xdim], diff, ax=ax[2])
            ax[2].text(0.4, 0.4, empty_message, transform=ax[2].transAxes, bbox=props)
        else:
            img2, ax[2] = pltfunc(adata[xdim], diff, ax=ax[2], norm=cp_info['normdiff'],cmap=cp_info['cmapdiff'],levels=cp_info['levelsdiff'],**cp_info['contourf_opt'])
            cb2 = fig.colorbar(img2, ax=ax[2], location='right',**cp_info['colorbar_opt'])

        #Set plot titles
        ax[0].set_title(case_title, loc='left', fontsize=tiFontSize)
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)
        ax[2].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=tiFontSize)

        # style the plot:
        #Set Main title for subplots:
        st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=15)
        st.set_y(0.85)
        ax[-1].set_xlabel("LONGITUDE")
        if cp_info['plot_log_p']:
            [a.set_yscale("log") for a in ax]
        fig.text(-0.03, 0.5, 'PRESSURE [hPa]', va='center', rotation='vertical')

    else:
        line = Line2D([0], [0], label="$\mathbf{Test}:$"+f"{case_nickname} - years: {case_climo_yrs[0]}-{case_climo_yrs[-1]}",
                        color="#1f77b4") # #1f77b4 -> matplotlib standard blue

        line2 = Line2D([0], [0], label=base_title,
                        color="#ff7f0e") # #ff7f0e -> matplotlib standard orange



        fig, ax = plt.subplots(nrows=2)
        ax = [ax[0],ax[1]]

        pltfunc(adata[xdim], adata, ax=ax[0],color="#1f77b4") # #1f77b4 -> matplotlib standard blue
        pltfunc(bdata[xdim], bdata, ax=ax[0],color="#ff7f0e") # #ff7f0e -> matplotlib standard orange
        pltfunc(adata[xdim], diff, ax=ax[1], color="k")

        ax[1].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=10)

        #Set Main title for subplots:
        st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=15)
        st.set_y(1.02)

        fig.legend(handles=[line,line2],bbox_to_anchor=(-0.15, 0.87, 1.05, .102),loc="right",
                borderaxespad=0.0,fontsize=6,frameon=False)

        for a in ax:
            try:
                a.label_outer()
            except:
                pass
            #End except
        #End for
    #End if

    #Write the figure to provided workspace/file:
    fig.savefig(wks, bbox_inches='tight', dpi=300)

    #Close plots:
    plt.close()

#
#  -- zonal mean annual cycle --
#

def square_contour_difference(fld1, fld2, **kwargs):
    """Produce filled contours of fld1, fld2, and their difference with square axes.

    Intended use is latitude-by-month to show the annual cycle.
    Example use case: use climo files to get data, take zonal averages,
    rename "time" to "month" if desired,
    and then provide resulting DataArrays to this function.

    Parameters
    ----------
        fld1, fld2 : xarray.DataArray
            2-dimensional DataArrays with same shape
        **kwargs : dict, optional
            optional keyword arguments
            this function _only checks_ `kwargs` for `case1name`, `case2name`

    Returns
    -------
    fig
        figure object

    Notes
    -----
    Assumes `fld1.shape == fld2.shape` and `len(fld1.shape) == 2`

    Will try to label the cases by looking for
    `case1name` and `case2name` in `kwargs`,
    and then `fld1['case']` and `fld2['case']` (i.e., attributes)
    If no case name is found proceeds with empty strings.
    **IF THERE IS A BETTER CONVENTION WE SHOULD USE IT.**

    Each panel also puts the Min/Max values into the title string.

    Axis labels are upper-cased names of the coordinates of `fld1`.
    Ticks are automatic with the exception that if the
    first coordinate is "month" and is length 12, use `np.arange(1,13)`.
    """

    assert len(fld1.shape) == 2,     "Input fields must have exactly two dimensions."  # 2-Dimension => plot contourf
    assert fld1.shape == fld2.shape, "Input fields must have the same array shape."    # Same shape => allows difference


    if "case1name" in kwargs:
        case1name = kwargs.pop("case1name")
    elif hasattr(fld1, "case"):
        case1name = getattr(fld1, "case")
    else:
        case1name = ""

    if "case2name" in kwargs:
        case2name = kwargs.pop("case2name")
    elif hasattr(fld2, "case"):
        case2name = getattr(fld2, "case")
    else:
        case2name = ""

    # Geometry of the figure is hard-coded
    fig = plt.figure(figsize=(10,10))

    rows = 5
    columns = 5
    grid = mpl.gridspec.GridSpec(rows, columns, wspace=1, hspace=1,
                            width_ratios=[1,1,1,1,0.2],
                            height_ratios=[1,1,1,1,0.2])
    # plt.subplots_adjust(wspace= 0.01, hspace= 0.01)
    ax1 = plt.subplot(grid[0:2, 0:2])
    ax2 = plt.subplot(grid[0:2, 2:4])
    ax3 = plt.subplot(grid[2:4, 1:3])
    # color bars / means share top bar.
    cbax_top = plt.subplot(grid[0:2, -1])
    cbax_bot = plt.subplot(grid[-1, 1:3])

    # determine color normalization for means:
    mx = np.max([fld1.max(), fld2.max()])
    mn = np.min([fld1.min(), fld2.min()])
    mnorm = mpl.colors.Normalize(mn, mx)

    coord1, coord2 = fld1.coords  # ASSUMES xarray WITH coords AND 2-dimensions
    print(f"{coord1}, {coord2}")
    xx, yy = np.meshgrid(fld1[coord2], fld1[coord1])
    print(f"shape of meshgrid: {xx.shape}")

    img1 = ax1.contourf(xx, yy, fld1.transpose())
    if (coord1 == 'month') and (fld1.shape[0] ==12):
        ax1.set_xticks(np.arange(1,13))
    ax1.set_ylabel(coord2.upper())
    ax1.set_xlabel(coord1.upper())
    ax1.set_title(f"{case1name}\nMIN:{fld1.min().values:.2f}  MAX:{fld1.max().values:.2f}")

    ax2.contourf(xx, yy, fld2.transpose())
    if (coord1 == 'month') and (fld1.shape[0] ==12):
        ax2.set_xticks(np.arange(1,13))
    ax2.set_xlabel(coord1.upper())
    ax2.set_title(f"{case2name}\nMIN:{fld2.min().values:.2f}  MAX:{fld2.max().values:.2f}")


    diff = fld1 - fld2
    ## USE A DIVERGING COLORMAP CENTERED AT ZERO
    ## Special case is when min > 0 or max < 0
    dmin = diff.min()
    dmax = diff.max()
    if dmin > 0:
        dnorm = mpl.colors.Normalize(dmin, dmax)
        cmap = mpl.cm.OrRd
    elif dmax < 0:
        dnorm = mpl.colors.Normalize(dmin, dmax)
        cmap = mpl.cm.BuPu_r
    else:
        dnorm = mpl.colors.TwoSlopeNorm(vmin=dmin, vcenter=0, vmax=dmax)
        cmap = mpl.cm.RdBu_r

    img3 = ax3.contourf(xx, yy, diff.transpose(), cmap=cmap, norm=dnorm)
    if (coord1 == 'month') and (fld1.shape[0] ==12):
        ax3.set_xticks(np.arange(1,13))
    ax3.set_ylabel(coord2.upper())
    ax3.set_xlabel(coord1.upper())
    ax3.set_title(f"DIFFERENCE (= a - b)\nMIN:{diff.min().values:.2f}  MAX:{diff.max().values:.2f}")


    # Try to construct the title:
    if hasattr(fld1, "long_name"):
        tstr = getattr(fld1, "long_name")
    elif hasattr(fld2, "long_name"):
        tstr = getattr(fld2, "long_name")
    elif hasattr(fld1, "short_name"):
        tstr = getattr(fld1, "short_name")
    elif hasattr(fld2, "short_name"):
        tstr = getattr(fld2, "short_name")
    elif hasattr(fld1, "name"):
        tstr = getattr(fld1, "name")
    elif hasattr(fld2, "name"):
        tstr = getattr(fld2, "name")
    else:
        tstr = ""
    if hasattr(fld1, "units"):
        tstr = tstr + f" [{getattr(fld1, 'units')}]"
    elif hasattr(fld2, "units"):
        tstr = tstr + f" [{getattr(fld2, 'units')}]"
    else:
        tstr = tstr + "[-]"

    fig.suptitle(tstr, fontsize=18)

    cb1 = fig.colorbar(img1, cax=cbax_top)
    cb2 = fig.colorbar(img3, cax=cbax_bot, orientation='horizontal')
    return fig

########


# WACCM Plots
#############

#Set monthly codes:
month_dict = {1:'JAN',
    		  2:'FEB',
    		  3:'MAR',
    		  4:'APR',
    		  5:'MAY',
    		  6:'JUN',
    		  7:'JUL',
    		  8:'AUG',
    		  9:'SEP',
    		  10:'OCT',
    		  11:'NOV',
    		  12:'DEC'}

delta_symbol = r'$\Delta$'

#temp_levs = np.arange(140, 300, 10)
#temp_diff_levs = np.arange(-40, 41, 4)

#wind_levs = np.arange(-120, 121, 10)
#wind_diff_levs = np.arange(-30, 31, 3)
#cont_ranges = {"U":{"levs":wind_levs,"diff_levs":wind_diff_levs,"units":"m/s"},
#               "T":{"levs":temp_levs,"diff_levs":temp_diff_levs,"units":"K"}}

#obs_cam_vars={"saber":{"U":"u", "T":"temp"},
#              "merra":{"U":"U", "T":"T"}}




def comparison_plots(plot_name, cam_var, case_names, case_ds_dict, obs_ds_dict, time_avg, interval, comp_plots_dict, obs_cam_vars):
    """

    """

    #Get plotting details for variable
    levs = np.arange(*comp_plots_dict[cam_var]["levs"])
    diff_levs = np.arange(*comp_plots_dict[cam_var]["diff_levs"])
    units = comp_plots_dict[cam_var]["units"]

    #Grab obs variable corresponding to CAM variable
    saber_var = obs_cam_vars['saber'][cam_var]
    merra_var = obs_cam_vars['merra'][cam_var]

    font_size = 8

    #Get number of test cases (number of columns)
    casenum = len(case_names)

    #Number of obs to compare
    #Currently, just compared to MERRA2 and SABER
    obsnum = 2
    nrows = obsnum+1

    #Set up plot
    fig = plt.figure(figsize=(casenum*4,nrows*5))

    for idx,case_name in enumerate(case_names):

        data_coords = case_ds_dict["coords"][case_name]

        data_lev = data_coords['lev']
        data_lat = data_coords['lat']

        #Set lat/lev grid for plotting
        [lat_grid, lev_grid] = np.meshgrid(data_lev,data_lat)

        if time_avg == "season":

            data_array = case_ds_dict["seasonal"][case_name][cam_var][interval]

            #Make Obs interpolated field from case
            merra_ds = obs_ds_dict["seasonal"]["merra"][merra_var][interval]
            merra_rfield = merra_ds.interp(lat=data_lat, lev=data_lev, method='linear')

            saber_ds = obs_ds_dict["seasonal"]["saber"][saber_var][interval]
            saber_rfield = saber_ds.interp(lat=data_lat, lev=data_lev, method='linear')

        if time_avg == "month":
            case_ds_dict["monthly"]
            str_month = interval #month_dict[interval]
            data_array = case_ds_dict["monthly"][case_name][cam_var][str_month]

            #Make Obs interpolated fields from case
            merra_ds = obs_ds_dict["monthly"]["merra"][merra_var][str_month]
            merra_rfield = merra_ds.interp(lat=data_lat, lev=data_lev, method='linear')

            saber_ds = obs_ds_dict["monthly"]["saber"][saber_var][str_month]
            saber_rfield = saber_ds.interp(lat=data_lat, lev=data_lev, method='linear')

        #Case plots (contours and contour fill)
        #######################################

        #Set up set of axes for first row
        ax = fig.add_subplot(nrows, casenum, idx+1)

        #Plot case contour fill
        cf=plt.contourf(lev_grid, lat_grid,
                        data_array.transpose(transpose_coords=True),
                        levels=levs, cmap='RdYlBu_r')

        #Plot case contours (for highlighting)
        contour = plt.contour(lev_grid, lat_grid,
                        data_array.transpose(transpose_coords=True),
                    colors="black",linewidths=0.5,levels=levs,zorder=100)
        fmt = {lev: '{:.0f}'.format(lev) for lev in contour.levels}
        ax.clabel(contour, contour.levels[::2], inline=True, fmt=fmt, fontsize=8)

        #Format axes
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        ax.tick_params(axis='x', labelsize=8)
        if idx > 0:
            plt.yticks([])
        else:
            plt.ylabel('hPa',fontsize=10)
        ax.tick_params(axis='y', labelsize=8)

        #Set individual plot title
        plt.title(case_name, fontsize=font_size)

        #Make colorbar on last plot only
        if idx == casenum-1:
            axins = inset_axes(ax,
                        width="5%",
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                       )
            cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label=units)
            cbar.ax.tick_params(axis='y', labelsize=8)
            cbar.set_label(units, fontsize=10, labelpad=1)

        #Difference with MERRA2 and MERRA2 contours
        ###########################################

        #Set up new set of axes for second row
        ax = fig.add_subplot(nrows, casenum, casenum+idx+1)

        #Plot interpolated contour
        contour = plt.contour(lev_grid, lat_grid, merra_rfield.transpose(transpose_coords=True),
                    colors='black', levels=levs,
                    negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        #fmt = {lev: '{:.0f}'.format(lev) for lev in contour.levels}
        #ax.clabel(contour, contour.levels[::4], inline=True, fmt=fmt, fontsize=8)
        #if idx == 0:
        #Add a legend for the contour lines for first plot only
        legend_elements = [Line2D([0], [0],
                               color=contour.collections[0].get_edgecolor(),
                               label=f'MERRA2 interp {cam_var}')]

        ax.legend(handles=legend_elements, loc='upper right', fontsize=5, bbox_to_anchor=(1., 1.))
        #End if

        #Plot difference contour fill
        cf=plt.contourf(lev_grid, lat_grid,
                        (data_array-merra_rfield).transpose(transpose_coords=True),
                        levels=diff_levs, cmap='RdYlBu_r')
        #fmt = {lev: '{:.0f}'.format(lev) for lev in cf.levels}
        #ax.clabel(cf, cf.levels[::3], inline=True, fmt=fmt, fontsize=8)

        #Plot case contours (for highlighting)
        contour = plt.contour(lev_grid, lat_grid,
                        (data_array-merra_rfield).transpose(transpose_coords=True),
                    colors="black",linewidths=0.5,levels=diff_levs[::2],zorder=100)
        fmt = {lev: '{:.0f}'.format(lev) for lev in contour.levels}
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=8)

        #Format axes
        plt.yscale("log")
        ax.set_ylim(1000,0.1)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        ax.tick_params(axis='x', labelsize=8)
        if idx > 0:
            plt.yticks([])
        else:
            plt.ylabel('hPa',fontsize=10)

        ax.tick_params(axis='y', labelsize=8)

        #Set individual plot title
        local_title = f'{case_name}\n {delta_symbol} from MERRA2'
        plt.title(local_title, fontsize=font_size)

        #Make colorbar on last plot only
        if idx == casenum-1:
            axins = inset_axes(ax,
                        width="5%",
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                       )
            cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label=units)
            cbar.ax.tick_params(axis='y', labelsize=8)
            cbar.set_label(units, fontsize=10, labelpad=1)

        #Difference with SABER and SABER contours (Temp only)
        #####################################################
        if cam_var == "T":
            #Set up new set of axes for third row
            ax = fig.add_subplot(nrows, casenum, (casenum*2)+idx+1)

            #Plot interpolated contour
            contour = plt.contour(lev_grid, lat_grid, saber_rfield.transpose(transpose_coords=True),
                        colors='black', levels=levs,
                        negative_linestyles='dashed', linewidths=.5, alpha=0.5)
            #if idx == 0:
            #Add a legend for the contour lines for first plot only
            legend_elements = [Line2D([0], [0],
                                color=contour.collections[0].get_edgecolor(),
                                label='SABER interp T')]

            ax.legend(handles=legend_elements, loc='upper right', fontsize=5, bbox_to_anchor=(1., 1.))
            #End if

            #Plot difference contour fill
            cf=plt.contourf(lev_grid, lat_grid,
                            (data_array-saber_rfield).transpose(transpose_coords=True),
                            levels=diff_levs, cmap='RdYlBu_r')

            #Plot case contours (for highlighting)
            contour = plt.contour(lev_grid, lat_grid,
                            (data_array-saber_rfield).transpose(transpose_coords=True),
                        colors="black",linewidths=0.5,levels=diff_levs,zorder=100)
            fmt = {lev: '{:.0f}'.format(lev) for lev in contour.levels}
            ax.clabel(contour, contour.levels[::2], inline=True, fmt=fmt, fontsize=8)

            #Format axes
            plt.yscale("log")
            ax.set_ylim(1000,0.1)
            plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
            plt.xticks(np.arange(-90,91,45),rotation=40)
            ax.tick_params(axis='x', labelsize=8)
            if idx > 0:
                plt.yticks([])
            else:
                plt.ylabel('hPa',fontsize=10)
            ax.tick_params(axis='y', labelsize=8)

            #Set individual plot title
            local_title = f'{case_name}\n {delta_symbol} from SABER'
            plt.title(local_title, fontsize=font_size)

            #Make colorbar on last plot only
            if idx == casenum-1:
                axins = inset_axes(ax,
                            width="5%",
                            height="80%",
                            loc='center right',
                            borderpad=-1.5
                        )
                cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label=units)
                cbar.ax.tick_params(axis='y', labelsize=8)
                # Set the font size for the colorbar label
                cbar.set_label(units, fontsize=10, labelpad=1)

    #Set up main plot title
    if time_avg == "month":
        str_interval = interval.lower().capitalize()
    else:
        str_interval = interval
    fig.suptitle(f"Zonal Mean {cam_var} - {str_interval}",fontsize=12,y=0.91)
    
    fig.savefig(plot_name, bbox_inches='tight', dpi=300)

    plt.close()



def polar_cap_temp(plot_name, hemi, case_names, cases_coords, cases_monthly, merra2_monthly, pcap_dict):
    """
    """

    levs = np.arange(*pcap_dict["T"]["levs"])

    #Get number of test cases (number of columns)
    casenum = len(case_names)

    font_size = 8
    if hemi == "s":
        slat = -90
        nlat = -60
        title_ext = f"{np.abs(nlat)}-{np.abs(slat)}\u00b0S"

    if hemi == "n":
        slat = 60
        nlat = 90
        title_ext = f"{slat}-{nlat}\u00b0N"

    """
    nplots = len(case_names)
    if nplots > 4:
        ncols = 4
    else:
        ncols = nplots
    #End if
    """
    ncols = 4
    nrows = 2

    #fig = plt.figure(figsize=(ncols*7,nrows*5))
    fig = plt.figure(figsize=(casenum*4,nrows*5))

    for idx,case_name in enumerate(case_names):
        ds = cases_coords[case_name]
        ds_month = cases_monthly[case_name]

        rfield_seas = np.zeros((12,len(ds['lev']),len(ds['lat'])))
        rfield_seas = xr.DataArray(rfield_seas, dims=['month','lev', 'lat'],
                                            coords={'month': np.arange(1,13,1),
                                                    'lev': ds['lev'],
                                                    'lat': ds['lat']})

        case_seas = np.zeros((12,len(ds['lev']),len(ds['lat'])))
        case_seas = xr.DataArray(case_seas, dims=['month','lev', 'lat'],
                                 coords={'month': np.arange(1,13,1),
                                         'lev': ds['lev'],
                                         'lat': ds['lat']})
        #Make array of monthly temp data
        for m in range(0,12):
            rfield_seas[m] = merra2_monthly['T'][month_dict[m+1]].interp(lat=ds['lat'], lev=ds['lev'],
                                                                method='linear')
            case_seas[m] = ds_month['T'][month_dict[m+1]]

        #Average over set of latitudes
        merra2_pcap = coslat_average(rfield_seas,slat,nlat)
        case_pcap = coslat_average(case_seas,slat,nlat)

        #
        [time_grid, lev_grid] = np.meshgrid(ds['lev'],np.arange(0,12))

        #Set up first row - Temps
        ax = fig.add_subplot(nrows, casenum, idx+1)
        cf=plt.contourf(lev_grid, time_grid, case_pcap,
                        levels=levs,cmap='RdYlBu_r'
                       ) #np.arange(-10,11,1)
        c0=plt.contour(lev_grid, time_grid, case_pcap, colors='grey',
                           levels=levs[::2],
                           negative_linestyles='dashed',
                           linewidths=.5, alpha=1)
        fmt = {lev: '{:.0f}'.format(lev) for lev in c0.levels}
        ax.clabel(c0, c0.levels, inline=True, fmt=fmt, fontsize=8)

        #Format the axes
        plt.yscale("log")
        ax.set_ylim(300,1)
        ax.set_yticks([300,100,30,10])
        ax.set_xticks(np.arange(0,12,2),rotation=40)
        ax.set_xticklabels(('Jan','Mar','May','Jul','Sep','Nov'),rotation=40,fontsize=8)
        if idx > 0:
            plt.yticks([])
        else:
            ax.set_yticklabels(["","$10^{2}$","","$10^{1}$"],fontsize=10)
            plt.ylabel('hPa',fontsize=10)
        
        #Set title
        local_title=f"{case_names[idx]}"
        plt.title(local_title, fontsize=font_size)

        #Make colorbar on last plot only
        if idx == casenum-1:
            axins = inset_axes(ax,
                        width="5%",
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                       )
            cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label="K")
            cbar.ax.tick_params(axis='y', labelsize=8)
            cbar.set_label("K", fontsize=10, labelpad=1)


        #Set up second row - Temp anomlies and Merra2 contours
        ax = fig.add_subplot(nrows, casenum, casenum+idx+1)
        clevs = np.arange(-10,11,1)
        cf=plt.contourf(lev_grid, time_grid, (case_pcap-merra2_pcap),
                        levels=np.arange(-10,11,1),cmap='RdYlBu_r'
                       ) #np.arange(-10,11,1)
        c0=plt.contour(lev_grid, time_grid, (case_pcap-merra2_pcap), colors='grey',
                           levels=clevs[::3],
                           negative_linestyles='dashed',
                           linewidths=.5, alpha=1)
        fmt = {lev: '{:.0f}'.format(lev) for lev in c0.levels}
        ax.clabel(c0, c0.levels, inline=True, fmt=fmt, fontsize=8)

        c=plt.contour(lev_grid, time_grid, merra2_pcap, colors='black',
                           levels=levs,
                           negative_linestyles='dashed',
                           linewidths=.5, alpha=0.5)

        #fmt = {lev: '{:.0f}'.format(lev) for lev in c.levels}
        #ax.clabel(c, c.levels, inline=True, fmt=fmt, fontsize=8)
        #if idx == 0:
        #Add a legend for the contour lines for first plot only
        legend_elements = [Line2D([0], [0],
                               color=c.collections[0].get_edgecolor(),
                               label='MERRA2 interp T')]

        ax.legend(handles=legend_elements, loc='upper right', fontsize=5, bbox_to_anchor=(1., 1.))
        #Format the axes
        plt.yscale("log")
        ax.set_ylim(300,1)
        ax.set_yticks([300,100,30,10])
        ax.set_xticks(np.arange(0,12,2),rotation=40)
        ax.set_xticklabels(('Jan','Mar','May','Jul','Sep','Nov'),rotation=40,fontsize=8)
        if idx > 0:
            plt.yticks([])
        else:
            ax.set_yticklabels(["","$10^{2}$","","$10^{1}$"],fontsize=10)
            plt.ylabel('hPa',fontsize=10)

        #Set title
        local_title=f"{case_names[idx]}\n {delta_symbol} from MERRA2"
        plt.title(local_title, fontsize=font_size)

        # Calculate the required wspace based on the length of titles
        #title_lengths = [len(ax.get_title()) for ax in axs]
        #max_title_length = max(title_lengths)
        #required_wspace = .003 * max_title_length  # Adjust the multiplier as needed
        #required_wspace = 0.
        # Adjust the wspace dynamically
        #plt.subplots_adjust(wspace=required_wspace)
        #Make colorbar on last plot only
        if idx == casenum-1:
            axins = inset_axes(ax,
                                width="5%",
                                height="80%",
                                loc='center right',
                                borderpad=-1.5
                               )
            cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label='K', ticks=np.arange(-9,10,3))
            cbar.ax.tick_params(axis='y', labelsize=8)
            # Set the font size for the colorbar label
            cbar.set_label("K", fontsize=10, labelpad=1)

        """
        #Check for start of new row
        if idx % 4 == 0:
            row = idx // 4 + 1

        #Check to see where the colorbar will go
        #The idea is to if the plots fill up each row, put the colorbar on last plot of row
        #If the row isn't filled up, put the color bar on last possible plot of row
        if ((4*(row-1) < idx < 4*(row+1)) and (idx == casenum-1)) or ((idx+1) % 4 == 0):
                axins = inset_axes(ax,
                                width="5%",
                                height="80%",
                                loc='center right',
                                borderpad=-1.5
                               )
                cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label='K', ticks=np.arange(-9,10,3))
                cbar.ax.tick_params(axis='y', labelsize=8)
                # Set the font size for the colorbar label
                cbar.set_label("K", fontsize=10, labelpad=1)
        """
    
    
    if hemi == "s":
        ptype = "SHPolar"
    if hemi == "n":
        ptype = "NHPolar"

    fig.suptitle(f"{hemi.upper()}H Polar Cap Temp Anomolies - {title_ext}",fontsize=12,y=0.93) #,horizontalalignment="center"
 
    fig.savefig(plot_name, bbox_inches='tight', dpi=300)
    
    #Close plots:
    plt.close()
########


#def cold_point_temp(var, var_dict, plot_name, case_names, case_runs, cases_monthly):
def month_vs_lat_plot(var, var_dict, plot_name, case_names, case_nicknames, climo_yrs, case_runs, cases_monthly, vert_lev):
    """
    """

    ahh = []
    for i in list(month_dict.keys())[::3]:
        ahh.append(month_dict[i].lower().capitalize())

    #Grab values for the month vs lat plot in variable defaults yaml file
    slat = var_dict[var]["slat"]
    nlat = var_dict[var]["nlat"]
    cmap = var_dict[var]["cmap"]
    diff_cmap = var_dict[var]["diff_cmap"]
    levs = np.arange(*var_dict[var]["levels"])
    diff_levs = np.arange(*var_dict[var]["diff_levels"])
    units = var_dict[var]["units"]
    title = var_dict[var]["title"]
    y_labels = var_dict[var]["y_labels"]
    tick_inter = var_dict[var]["tick_inter"]


    

    nplots = len(case_names)
    if nplots > 4:
        ncols = 4
    else:
        ncols = nplots
    #End if
    ncols = 2
    nrows = int(np.ceil(nplots/ncols))

    #fig = plt.figure(figsize=(2*8,nrows*5))
    #Gather contour plot options
    #cp_info = prep_contour_plot(mseasons, oseasons, dseasons, **vres)

    # create figure:
    fig = plt.figure(figsize=(16,10))

    # LAYOUT WITH GRIDSPEC
    gs = mpl.gridspec.GridSpec(4, 8, wspace=0.95,hspace=0.5)
    ax1 = plt.subplot(gs[0:2, :4])#, **cp_info['subplots_opt'])
    ax2 = plt.subplot(gs[0:2, 4:])#, **cp_info['subplots_opt'])
    ax3 = plt.subplot(gs[2:, 2:6])#, **cp_info['subplots_opt'])
    ax = [ax1,ax2,ax3]

    #
    pcap_vals = {}

    #for run in range(len(runs)):
    for idx,case_name in enumerate(case_names):
        ds = case_runs[case_name]
        ds_month = cases_monthly[case_name]

        #if idx == len(case_names)-1:
        #    levs = diff_levs
        #    cmap = "BrBG"

        #Make 24 months so we can have Jan-Dec repeated twice
        case_seas = np.zeros((25,len(ds['lev']),len(ds['lat'])))
        case_seas = xr.DataArray(case_seas, dims=['month','lev', 'lat'],
                                 coords={'month': np.arange(1,26,1),
                                         'lev': ds['lev'],
                                         'lat': ds['lat']})
        #Make array of monthly temp data
        for m in range(0,25):
            month = m
            if m > 11:
                month = m-12
            if month == 12:
                month = 0

            case_seas[m] = ds_month[var][month_dict[month+1]]

        #Average over set of latitudes
        case_pcap = coslat_average(case_seas,slat,nlat)
        case_pcap = case_seas.sel(lev=vert_lev,method="nearest").sel(lat=slice(slat, nlat))
        if var == "Q":
            case_pcap = case_pcap*1e6

        pcap_vals[case_name] = case_pcap
        #cases_monthly["diff"] = cases_monthly[case_names[0]] - cases_monthly[case_names[1]]
        #case_runs["diff"]

        #
        [time_grid, lat_grid] = np.meshgrid(ds['lat'].sel(lat=slice(slat, nlat)),
                                            np.arange(0,25))
        #Set up plot
        #ax = fig.add_subplot(nrows, ncols, idx+1)

        """cf=plt.contourf(lat_grid, time_grid, (case_pcap),
                        levels=levs,
                        cmap=cmap,#zorder=100
                      )
        c=plt.contour(lat_grid, time_grid, (case_pcap),
                        levels=levs,
                        colors='k',linewidths=0.5,alpha=0.5
                      )"""

        cf=ax[idx].contourf(lat_grid, time_grid, (case_pcap),
                        levels=levs,
                        cmap=cmap,#zorder=100
                      )
        c=ax[idx].contour(lat_grid, time_grid, (case_pcap),
                        levels=levs,
                        colors='k',linewidths=0.5,alpha=0.5
                      )
        
        # add contour labels
        #lb = plt.clabel(c, fontsize=6, inline=True, fmt='%r')

        # Format contour labels
        if var == "T":
            fmt = {lev: '{:.0f}'.format(lev) for lev in c.levels}
        else:
            fmt = {lev: '{:.1f}'.format(lev) for lev in c.levels}
        ax[idx].clabel(c, c.levels[::2], inline=True, fmt=fmt, fontsize=8)

        #Add a horizontal line at 0 degrees latitude
        #plt.axhline(0, color='grey', linestyle='-',zorder=200,alpha=0.7)
        ax[idx].axhline(0, color='grey', linestyle='-',zorder=200,alpha=0.7)

        #Format the x-axis
        ax[idx].set_xticks(np.arange(0,25,3),rotation=40)
        ax[idx].set_xticklabels(ahh+ahh+["Jan"],rotation=40)

        #Set title
        if idx == 0:
            #eyear_cases = eyear_cases + [eyear_baseline]
            #climo_yrs = [syear_cases, eyear_cases]

            plot_title = "$\mathbf{Test}:$"+f"{case_nicknames[0]}\nyears: {climo_yrs[0][0]}-{climo_yrs[1][0]}"
        if idx == 1:
            plot_title = "$\mathbf{Baseline}:$"+f"{case_nicknames[1]}\nyears: {climo_yrs[0][1]}-{climo_yrs[1][1]}"
        ax[idx].set_title(plot_title, loc='left', fontsize=8)
        #local_title=case_names[idx]
        #plt.title(local_title, fontsize=8)

        #Check for start of new row
        if idx % 2 == 0:
            row = idx // 2 + 1

        #Add latitude label to first column of each row
        #if idx==2*(row-1):
        #   plt.ylabel('Latitude',fontsize=10)
        
        #plt.ylabel('Latitude',fontsize=10)
        ax[idx].set_ylabel('Latitude',fontsize=10)

        #Format the y-axis
        ax[idx].set_yticks(np.arange(slat,nlat+1,tick_inter))
        ax[idx].set_yticklabels(y_labels,fontsize=10)

        """#Check to see where the colorbar will go
        if ((idx==2*(row-1)) and (idx == nplots-1)) or ((idx+1) % 2 == 0):
                axins = inset_axes(ax,
                                width="3%",
                                height="80%",
                                loc='center right',
                                borderpad=-1.5
                               )
                cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label=units,
                                    #ticks=levs
                                   )
                cbar.add_lines(c)
                cbar.ax.tick_params(axis='y', labelsize=8)
                # Set the font size for the colorbar label
                cbar.set_label(units, fontsize=10, labelpad=1)"""

        axins = inset_axes(ax[idx],
                                width="3%",
                                height="80%",
                                loc='center right',
                                borderpad=-1.5
                               )
        cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label=units,
                                    #ticks=levs
                                   )
        cbar.add_lines(c)
        cbar.ax.tick_params(axis='y', labelsize=8)
        # Set the font size for the colorbar label
        cbar.set_label(units, fontsize=10, labelpad=1)
    #End cases


    #Difference Plots
    #----------------

    idx = 2

    #ds = case_runs[case_name]
    #ds_month = cases_monthly[case_name]

    #cases_monthly["diff"] = cases_monthly[case_names[0]] - cases_monthly[case_names[1]]
    #case_runs["diff"] = case_runs[case_names[0]] - case_runs[case_names[1]]

    diff_pcap = pcap_vals[case_names[0]] - pcap_vals[case_names[1]]
    #ds_month = case_runs[case_names[0]] - case_runs[case_names[1]]

    #if idx == len(case_names)-1:
    levs = diff_levs
    cmap = diff_cmap
    """
    #Make 24 months so we can have Jan-Dec repeated twice
    case_seas = np.zeros((25,len(ds['lev']),len(ds['lat'])))
    case_seas = xr.DataArray(case_seas, dims=['month','lev', 'lat'],
                                 coords={'month': np.arange(1,26,1),
                                         'lev': ds['lev'],
                                         'lat': ds['lat']})
    #Make array of monthly temp data
    for m in range(0,25):
        month = m
        if m > 11:
            month = m-12
        if month == 12:
            month = 0

        case_seas[m] = ds_month[var][month_dict[month+1]]

    #Average over set of latitudes
    case_pcap = coslat_average(case_seas,slat,nlat)
    case_pcap = case_seas.sel(lev=vert_lev,method="nearest").sel(lat=slice(slat, nlat))
    if var == "Q":
        case_pcap = case_pcap*1e6"""

    #
    [time_grid, lat_grid] = np.meshgrid(ds['lat'].sel(lat=slice(slat, nlat)),
                                            np.arange(0,25))
    #Set up plot
    #ax = fig.add_subplot(nrows, ncols, idx+1)

       
    cf=ax[idx].contourf(lat_grid, time_grid, (diff_pcap),
                        levels=levs,
                        cmap=cmap,#zorder=100
                      )
    c=ax[idx].contour(lat_grid, time_grid, (diff_pcap),
                        levels=levs,
                        colors='k',linewidths=0.5,alpha=0.5
                      )
        
    # add contour labels
    #lb = plt.clabel(c, fontsize=6, inline=True, fmt='%r')

    # Format contour labels
    if var == "T":
        fmt = {lev: '{:.0f}'.format(lev) for lev in c.levels}
    else:
        fmt = {lev: '{:.1f}'.format(lev) for lev in c.levels}
    ax[idx].clabel(c, c.levels[::2], inline=True, fmt=fmt, fontsize=8)

    #Add a horizontal line at 0 degrees latitude
    #plt.axhline(0, color='grey', linestyle='-',zorder=200,alpha=0.7)
    ax[idx].axhline(0, color='grey', linestyle='-',zorder=200,alpha=0.7)

    #Format the x-axis
    ax[idx].set_xticks(np.arange(0,25,3),rotation=40)
    ax[idx].set_xticklabels(ahh+ahh+["Jan"],rotation=40)

    #Set title
    local_title="$\mathbf{Test} - \mathbf{Baseline}$" #"Test - Baseline"#case_names[idx]
    #plt.title(local_title, fontsize=8)
    ax[idx].set_title(local_title, fontsize=8)

    #Check for start of new row
    #if idx % 2 == 0:
    #    row = idx // 2 + 1

    #Add latitude label to first column of each row
    #if idx==2*(row-1):
    #   plt.ylabel('Latitude',fontsize=10)
        
    #plt.ylabel('Latitude',fontsize=10)
    ax[idx].set_ylabel('Latitude',fontsize=10)

    #Format the y-axis
    ax[idx].set_yticks(np.arange(slat,nlat+1,tick_inter))
    ax[idx].set_yticklabels(y_labels,fontsize=10)

       

    axins = inset_axes(ax[idx],
                                width="3%",
                                height="80%",
                                loc='center right',
                                borderpad=-1.5
                               )
    cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label=units,
                                    #ticks=levs
                                   )
    cbar.add_lines(c)
    cbar.ax.tick_params(axis='y', labelsize=8)
    # Set the font size for the colorbar label
    cbar.set_label(units, fontsize=10, labelpad=1)





    fig.suptitle(f"{title} - {vert_lev}hPa",fontsize=16,y=0.99,horizontalalignment="center")

    fig.savefig(plot_name, bbox_inches='tight', dpi=300)


########

#WACCM QBO
import numpy as np

def qbo_amplitude(data):
    """
    Calculate the QBO amplitude
    """
    
    from scipy.signal import convolve
    
    boxcar = np.ones((6, 1)) / 6
    filtered_data = convolve(data, boxcar, mode='valid')
    amplitude=np.std(filtered_data, axis=0)
    
    return amplitude


def qbo_frequency(data):
    """
    Calculate the QBO frequency
    """
    
    [dt,dx]=data.shape
    dt2=int(dt/2)
    f=1*np.arange(0,dt2+1,1)/dt
    f=f[1:]
    f = np.tile(f, (dx, 1)).swapaxes(0,1)
    
    fft_data = np.fft.fft(data, axis=0)
    fft_data = fft_data[1:dt2+1,:]
    
    power_spectrum = np.abs(fft_data)**2
    
    period=np.sum(power_spectrum*(1/f),axis=0)/np.sum(power_spectrum,axis=0)
    
    return period
  

def waccm_qbo(plot_name, case_names, nicknames, case_runs, merra2, syear_cases, eyear_cases):
    """
    
    """

    def format_side_axes(axes, side_axis, x, merra_data, data=None, case_lev=None, merra=False):
        """
        Format the period and amplitiude side axes
        """
        axes[side_axis].plot(merra_data,merra2['lev'],color='k')# s=1
        axes[side_axis].set_ylim(y_lims[0],y_lims[1])
        axes[side_axis].set_yscale("log")

        if merra==False:
            axes[side_axis].plot(data,case_lev)#linewidths=1
        if x == "period":
            axes[side_axis].set_xlim(0,40)
            axes[side_axis].set_xticks(np.arange(0,41,10))
            axes[side_axis].set_xticklabels(np.arange(0,41,10),fontsize=8)
            axes[side_axis].set_xlabel('months',fontsize=10)
        if x == "amplitude":
            axes[side_axis].set_xlim(0,20)
            axes[side_axis].set_xticks(np.arange(0,21,5))
            axes[side_axis].set_xticklabels(np.arange(0,21,5),fontsize=8)
            axes[side_axis].set_xlabel('m/s',fontsize=10)
        axes[side_axis].set_yticks([])
        return axes

    def advance_string(input_string):
        advanced_chars = [chr(ord(char) + 3) if char.isalpha() else char for char in input_string]
        advanced_string = ''.join(advanced_chars)
        return advanced_string


    #Build subplot mosiac based off number of CAM cases
    input_string0 = 'AAAABC'
    ahh = []
    ahh.append(input_string0)
    for i in range(len(case_names)):
        if i ==0:
            input_string = advance_string(input_string0)
        else:
            input_string = advance_string(input_string)
            input_string = f"{input_string}"
        ahh.append(input_string)

    main_key = []
    side1_key = []
    side2_key = []
    mos_str = input_string0
    for idx,i in enumerate(ahh):
        if idx != 0:
            mos_str += f";{i}"
        main_key.append(i[0])
        side1_key.append(i[-2])
        side2_key.append(i[-1])

    fig, axes = plt.subplot_mosaic(mos_str,figsize=(12,5*len(case_names)))

    y = 1.00
    y_lims = [100,0.1]

    contour_levels = np.arange(-35, 36, 2.5)

    #Plot MERRA2 last; this will be based on number of CAM cases
    merra_idx = len(case_names)

    nt = 108
    plotdata = coslat_average(merra2['U'],-10,10)

    plotdata_clip = np.clip(np.abs(plotdata), None, 35)
    plotdata=np.sign(plotdata)*plotdata_clip
    [time_grid, lev_grid] = np.meshgrid(merra2['lev'],np.arange(1,nt+1,1))
    start_ind=252-12
    end_ind=start_ind+nt

    data = plotdata[start_ind:end_ind,:]
    cf = axes[main_key[merra_idx]].contourf(lev_grid, time_grid, data,
                                        levels=contour_levels, cmap='RdBu_r')

    c = axes[main_key[merra_idx]].contour(lev_grid, time_grid, data, alpha=0.5,linewidths=0.3,
                                        levels=contour_levels[::5], colors='k',linestyles=['dashed' if val < 0 else 'solid' for val in np.unique(data)])
    # add contour labels
    lb = plt.clabel(c, fontsize=6, inline=True, fmt='%r')
    """
    axins = inset_axes(axes[main_key[merra_plot]], width="3%", height="80%", loc='center right', borderpad=-0.5)
    cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label="m/s",
                                    ticks=contour_levels)
    cbar.add_lines(c)
    """
    axins = inset_axes(axes[main_key[merra_idx]], width="100%", height="5%", loc='lower center', borderpad = -3.5)
    cbar = fig.colorbar(cf, cax=axins, orientation="horizontal", label="m/s",
                                        ticks=contour_levels[::2])
    cbar.ax.tick_params(axis='x', labelsize=8)
    # Set the font size for the colorbar label
    cbar.set_label("m/s", fontsize=10, labelpad=1)

    axes[main_key[merra_idx]].set_ylim(y_lims[0],y_lims[1])
    axes[main_key[merra_idx]].set_yscale("log")
    axes[main_key[merra_idx]].set_ylabel('hPa',fontsize=10)
    axes[main_key[merra_idx]].tick_params(axis='y', labelsize=8)
    axes[main_key[merra_idx]].set_title("MERRA2",y=y,fontsize=10)
    axes[main_key[merra_idx]].set_xticks(np.arange(1,nt+1,12),rotation=40)

    start_year = int(str(plotdata[start_ind].time.values)[0:4])
    axes[main_key[merra_idx]].set_xticklabels(np.arange(start_year,start_year+(nt/12),1).astype(int),fontsize=8)

    #MERRA QBO Amplitude side axis
    amp_m = qbo_amplitude(plotdata)
    axes = format_side_axes(axes,side1_key[merra_idx],"amplitude",amp_m,merra=True)

    #MERRA QBO Period side axis
    period_m = qbo_frequency(plotdata)
    axes = format_side_axes(axes,side2_key[merra_idx],"period",period_m,merra=True)

    #Loop over CAM case data
    for idx,case_name in enumerate(case_names):
        case_data = case_runs[case_name]
        nickname = nicknames[idx]
        yrs = syear_cases[idx]
        
        #Get number of time steps
        nt = len(case_data['time'])
        #If the number is greater than 10 years, clip it to 10 years?
        if nt > 120:
            nt_sub = 120
            nt_sub = 108
        else:
            nt_sub = nt

        [time_grid, lev_grid] = np.meshgrid(case_data['lev'],np.arange(0,nt_sub+1,1))

        contour_levels = np.arange(-35, 35, 2.5)

        plotdata = coslat_average(case_data['U'],-10,10)
        plotdata_clip = np.clip(np.abs(plotdata), None, 35)
        plotdata=np.sign(plotdata)*plotdata_clip

        #TODO: this will need to be adjusted??
        #Curently this is finding (start_idx)th month and then going out 9 years
        #QUESTION: what if the data doesn't have 9 years? - we will need to clip this...
        start_idx = 0 #119-24
        #print(plotdata[start_idx:start_idx+(12*9),:].shape)
        end_idx = start_idx+(12*9)+1
        cf = axes[main_key[idx]].contourf(lev_grid[start_idx:end_idx,:], time_grid[start_idx:end_idx,:], plotdata[start_idx:end_idx,:],
                                    levels=contour_levels, cmap='RdBu_r')

        c = axes[main_key[idx]].contour(lev_grid[start_idx:end_idx,:], time_grid[start_idx:end_idx,:], plotdata[start_idx:end_idx,:],
                                    levels=contour_levels[::5], colors='k',alpha=0.5)
        # add contour labels
        lb = plt.clabel(c, fontsize=6, inline=True, fmt='%r')
        
        """
        axins = inset_axes(axes[main_key[idx]], width="3%", height="80%", loc='center right', borderpad=-0.5)
        cbar = fig.colorbar(cf, cax=axins, orientation="vertical", label="m/s",
                                        ticks=contour_levels)
        cbar.add_lines(c)
        cbar.ax.tick_params(axis='y', labelsize=8)
        # Set the font size for the colorbar label
        cbar.set_label("m/s", fontsize=10, labelpad=1)
        """
        axes[main_key[idx]].set_ylim(y_lims[0],y_lims[1])
        axes[main_key[idx]].set_yscale("log")
        axes[main_key[idx]].set_ylabel('hPa',fontsize=10)
        axes[main_key[idx]].tick_params(axis='y', labelsize=8)
        axes[main_key[idx]].set_title(nickname,y=y,fontsize=10)
        #print((nt_sub/12)+1)
        
        """
        if idx == 0:
            axes[main_key[idx]].set_xticks(np.arange(0,(nt_sub)+1,12),rotation=40)
        else:
            axes[main_key[idx]].set_xticks(np.arange(0,(nt_sub)+1,12),rotation=40)
        
        if idx == 0:
            axes[main_key[idx]].set_xticklabels(np.arange(int(yrs[0]+int(nt_sub/12)),int(yrs[0]+int(nt_sub/12))+int(nt_sub/12)+1,1))
        else:
            axes[main_key[idx]].set_xticklabels(np.arange(int(yrs[0]+int(nt_sub/12)),int(yrs[0]+int(nt_sub/12))+int(nt_sub/12)+1,1))
        """
        #print("nt_sub",nt_sub)
        axes[main_key[idx]].set_xticks(np.arange(0,(nt_sub)+1,12),rotation=40)
        #axes[main_key[idx]].set_xticklabels(np.arange(int(yrs+int(nt_sub/12)),int(yrs+int(nt_sub/12))+int(nt_sub/12)+1,1))
        yr0 = int(yrs+int(start_idx/12))
        axes[main_key[idx]].set_xticklabels(np.arange(yr0, yr0+int(nt_sub/12)+1, 1), fontsize=8)
        #axes[main_key[merra_plot]].tick_params(axis='y', labelsize=10)

        #Case QBO Amplitude side axis
        amp = qbo_amplitude(plotdata)
        axes = format_side_axes(axes, side1_key[idx], "amplitude",amp_m,amp,case_data['lev'])

        #Case QBO Period side axis
        period = qbo_frequency(plotdata)
        axes = format_side_axes(axes, side2_key[idx], "period",period_m,period,case_data['lev'])

        #Label first row of side axes only
        if idx==0:
            axes[side1_key[idx]].set_title('Amplitude',y=y,fontsize=12)
            axes[side2_key[idx]].set_title('Period',y=y,fontsize=12)

    # Adjust the vertical spacing (hspace)
    plt.subplots_adjust(hspace=0.35)

    fig.suptitle(f"QBO Diagnostics",fontsize=16,y=0.93,horizontalalignment="center")

    fig.savefig(plot_name, bbox_inches='tight', dpi=300)




def tem_plot(wks, case_nickname, base_nickname,
             case_climo_yrs, baseline_climo_yrs,
             mdlfld, obsfld, diffld, obs=False, **kwargs):

    # generate dictionary of contour plot settings:
    cp_info = prep_contour_plot(mdlfld, obsfld, diffld, **kwargs)

    # create figure object
    fig = plt.figure(figsize=(14,10))
    img = []

    # LAYOUT WITH GRIDSPEC
    gs = mpl.gridspec.GridSpec(3, 6, wspace=0.5,hspace=0.0) # 2 rows, 4 columns, but each map will take up 2 columns
    #gs.tight_layout(fig)
    ax1 = plt.subplot(gs[0:2, :3], **cp_info['subplots_opt'])
    ax2 = plt.subplot(gs[0:2, 3:], **cp_info['subplots_opt'])
    ax3 = plt.subplot(gs[2, 1:5], **cp_info['subplots_opt'])
    ax = [ax1,ax2,ax3]

    fields = (mdlfld, obsfld, diffld)
    for i, a in enumerate(fields):

        if i == len(fields)-1:
            levels = cp_info['levelsdiff']
            cmap = cp_info['cmapdiff']
            norm = cp_info['normdiff']
        else:
            levels = cp_info['levels1']
            cmap = cp_info['cmap1']
            norm = cp_info['norm1']

        levs = np.unique(np.array(levels))
        if len(levs) < 2:
            img.append(ax[i].contourf(lons,lats,a,colors="w",transform=ccrs.PlateCarree()))
            ax[i].text(0.4, 0.4, empty_message, transform=ax[i].transAxes, bbox=props)
        else:
            img.append(ax[i].contourf(lons, lats, a, levels=levels, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), **cp_info['contourf_opt']))
        #End if
        #ax[i].set_title("AVG: {0:.3f}".format(area_avg[i]), loc='right', fontsize=11)

        # add contour lines <- Unused for now -JN
        # TODO: add an option to turn this on -BM
        #cs.append(ax[i].contour(lon2, lat2, fields[i], transform=ccrs.PlateCarree(), colors='k', linewidths=1))
        #ax[i].clabel(cs[i], cs[i].levels, inline=True, fontsize=tiFontSize-2, fmt='%1.1f')
        #ax[i].text( 10, -140, "CONTOUR FROM {} to {} by {}".format(min(cs[i].levels), max(cs[i].levels), cs[i].levels[1]-cs[i].levels[0]),
        #bbox=dict(facecolor='none', edgecolor='black'), fontsize=tiFontSize-2)

#####################
#END HELPER FUNCTIONS
