"""                                                                    .
Generic computation helper functions

Functions
---------
load_dataset()
    generalized load dataset method used for plotting/analysis functions
mask_land_or_ocean(arr, msk, use_nan=False)
    Apply a land or ocean mask to provided variable.
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

Notes
-----

"""

import netCDF4
import pathlib as path
import pandas as pd
import numpy as np
import os
from datetime import datetime, timedelta
from fnmatch import fnmatch
from typing import Iterable
import xarray as xr

def find_flight_fnames(dir_path: str) -> list[str]:
    """
    find_flight_fnames just searches a directory for all *.nc files and returns a list of them.

    :param dir_path: a path to the directory containing flight netcdf files

    :return: Returns a list of flight netcdf files.
    """
    flight_paths=[]
    flight_fnames = sorted([fname for fname in os.listdir(dir_path) if fnmatch(fname, "*.nc")])
    for i in range(len(flight_fnames)):
        flight_paths.append(dir_path + '/' + flight_fnames[i])
    
    return flight_paths

def find_nc_fnames(dir_path: str) -> list[str]:
    """
    find_flight_fnames just searches a directory for all *.nc files and returns a list of them.
    
    :param dir_path: a path to the directory containing flight netcdf files
    
    :return: Returns a list of flight netcdf files.
    """
    nc_paths=[]
    nc_fnames = sorted([fname for fname in os.listdir(dir_path) if fnmatch(fname, "*.nc")])
    for i in range(len(nc_fnames)):
        nc_paths.append(dir_path + '/' + nc_fnames[i])
        
        nudg_path = [file for file in nc_paths if ".hs." in file]
        free_path = [file for file in nc_paths if ".h0." in file]
        # save dictionary with the paths for 
        paths = {'Free': free_path,'Nudg': nudg_path}
        
    return paths

def open_nc(flight_paths: str) -> netCDF4._netCDF4.Dataset:
    """
    open_flight_nc simply checks to see if the file at the provided path string exists and opens it.

    :param file_path: A path string to a flight data file, e.g. "./test/test_flight.nc"

    :return: Returns xr.open_dataset object.
    """
    fp_path = path.Path(flight_paths)
    if not fp_path.is_file():
        raise FileNotFoundError('testing excptions')

    return xr.open_dataset(flight_paths)

def read_flight_nc_1hz(nc: xr.open_dataset, read_vars) -> pd.DataFrame:
    """
    read_flight_nc reads a set of variables into memory.

    NOTE: a low-rate, 1 Hz, flight data file is assumed

    :param nc: netCDF4._netCDF4.Dataset object opened by open_flight_nc.
    :param read_vars: An list of strings of variable names to be read into memory.

    :return: Returns a pandas data frame.
    """
    long_names = [nc[var].long_name if 'long_name' in nc[var].attrs else None for var in read_vars]
    data = [] # an empty list to accumulate Dataframes of each variable to be read in
    for var in read_vars:
        try:
            if var == "Time":
                # df = xr.open_dataset(nc)
                time = np.array(nc.Time)
                data.append(pd.DataFrame({var: time}))
                # dt_list = sfm_to_datetime(time, tunits)
                # data.append(pd.DataFrame({'datetime': time}))
            else:
                output = nc[var][:]
                data.append(pd.DataFrame({var: output}))
        except Exception as e:
            print(f"Issue reading {var}: {e}")
            pass
    
    dataframe = pd.concat(data, axis=1, ignore_index=False)
    dataframe.attrs['long_names'] = long_names
    # concatenate the list of dataframes into a single dataframe and return it
    return dataframe

def read_flight_nc_25hz(nc: xr.open_dataset, read_vars) -> pd.DataFrame:
    """
    read_flight_nc reads a set of variables into memory.
    
    NOTE: a high-rate, usually 25 Hz, flight data file is assumed.
    
    :param nc: netCDF4._netCDF4.Dataset object opened by open_flight_nc.
    :param read_vars: An optional list of strings of variable names to be read into memory. A default
                      list, vars_to_read, is specified above. Passing in a similar list will read in those variables
                      instead.
    
    :return: Returns a pandas data frame.
    """
    data = []
    sub_seconds = np.arange(0, 25, 1)/25.
    hz = 25
    for var in read_vars:
        try:
            if var == "Time":
                time = nc[var].values  # Get NumPy array from Xarray
                # Convert sub_seconds into timedelta in nanoseconds
                sub_seconds_ns = (sub_seconds * 1e9).astype('timedelta64[ns]')
                # Expand time into 2D, add sub-second offsets
                time_25hz = time[:, None] + sub_seconds_ns
                output = time_25hz.ravel()  # Flatten to 1D
                data.append(pd.DataFrame({var: output}))
            else:
                ndims = len(np.shape(nc[var][:]))
                if ndims == 2:
                    # 2-D, 25 Hz variables can just be raveled into 1-D time series
                    output = np.ravel(nc[var].values)
                    data.append(pd.DataFrame({var: output}))
                elif ndims == 1:
                    values = nc[var].values  # Extract as NumPy array
                    if values.shape[0] != len(time):  # Interpolation case (e.g., GGALT-style)
                        print(f"Skipping {var} due to shape mismatch: {values.shape[0]} != {len(time)}")
                        continue
                    # Interpolate to 25 Hz (fudged interpolation)
                    output_2d = np.full((len(values), hz), np.nan)
                    for i in range(len(values) - 1):
                        output_2d[i, :] = values[i] + sub_seconds * (values[i+1] - values[i])
                    output = output_2d[:-1].ravel()  # remove the last NaN row
                    data.append(pd.DataFrame({var: output}))
        except Exception as e:
            print(f"Issue reading {var}: {e}")
            pass
    # concatenate the list of dataframes into a single dataframe and return it
    dataframe = pd.concat(data, axis=1, ignore_index=False)
    return dataframe

def read_flight_nc(nc: xr.open_dataset, vars2read: list[str]) -> pd.DataFrame:
    """
    read_flight_nc simply figures out if the flight netcdf object is 1 hz or 25 hz and calls the appropriate reader.

    :param nc: A netcdf object for a flight netcdf file.
    :param read_vars: A list of variable names to be read in the netcdf object. Optional. Default is "vars_to_read" specified
                      above.

    :return: Returns Pandas DataFrame
    """
    dim_names = list(nc.dims)
    if 'sps25' in dim_names:
        df = read_flight_nc_25hz(nc, vars2read)
    else:
        df = read_flight_nc_1hz(nc, vars2read)
    return df

# Function to read in all the relevant variables from the NSF aircraft datasets
def read_vars(nc):

    var_list = nc.data_vars
    time = 'Time'
    # Spatial variables
    lat, lon, alt = 'GGLAT', 'GGLON', 'GGALT'
    
    # state variables
    temp = 'ATX'
    dwpt = 'DPXC'
    u = 'UIC' if 'UIC' in var_list else 'UIX'
    v = 'VIC' if 'VIC' in var_list else 'VIX'
    w = 'WIC' if 'WIC' in var_list else 'WIX'
    p = 'PSXC'
    ew = 'EWX'
    rh = 'RHUM'
    vars_to_read = [time, lat, lon, alt, temp, dwpt, u,  w, p, ew, rh]
    # Thermodynamic data
    if any('THETA' in var for var in var_list): 
        theta_vars = [var for var in var_list if 'THETA' in var and ('_GP' not in var)]
        vars_to_read.extend(theta_vars)
    # Cloud microphysical
    if any('CONC' in var for var in var_list): # cloud concentrations
        conc_vars = [var for var in var_list if 'CONC' in var and 'D' in var and ('R_' not in var and 'CN' not in var and \
                    'CV' not in var and '0_' not in var and 'UD' not in var)]
        # print(conc_vars)
        vars_to_read.extend(conc_vars)
    if any('PLW' in var for var in var_list): # Liquid/Ice water contents
        # v = [var for var in var_list if '2' not in var]
        wc_vars = [var for var in var_list if 'PLW' in var and ('2V' not in var)]
        vars_to_read.extend(wc_vars)
    # Aerosol data
    if any('UHSAS' in var for var in var_list) or any('CONCN' in var for var in var_list):
        aer_var = [var for var in var_list if ('UHSAS' in var or 'CONCU' in var or 'CONCN' in var) and ('AU' not in var and 'UD'  not in var and 'CUH' not in var and 'CFDC' not in var)]
        # uhsas_cells = var_list['CUHSAS_LWII'].CellSizes
        vars_to_read.extend(aer_var)
    # print("Loaded variables:")
    # print(vars_to_read)
    return vars_to_read

def read_sizedist_vars(nc):
    # Ensure we’re iterating over plain strings (variable names)
    names = list(nc.data_vars.keys())

    out = []
        # Include Time if present (coord or data var)
    if 'Time' in nc:
        out.append('Time')
        
    def add_prefix(prefix, exclude_substr=None):
        for n in names:
            if n.startswith(prefix) and (exclude_substr is None or exclude_substr not in n):
                out.append(n)

    # Cloud probe size distributions
    add_prefix('CCDP')
    add_prefix('C2DCA')
    add_prefix('C2DSA')          # <-- startswith enforces "at the start"

    # Aerosol size distributions
    add_prefix('CUHSAS', exclude_substr='CVI')   # exclude any CUHSAS* containing CVI
    add_prefix('CS200') # PCASP

    # Deduplicate while preserving original order
    out = list(dict.fromkeys(out))
    
    return out

def _prep_probe(nc, varname):
    da = nc[varname]
    # collapse any sps* to 1 Hz
    sps_dims = [d for d in da.dims if d.lower().startswith('sps')]
    if sps_dims:
        da = da.mean(dim=sps_dims, keep_attrs=True)
    # find bin dim and order (Time, Bin)
    bin_dim = next(d for d in da.dims if d.lower().startswith(('vector','bin','cell')))
    time_name = 'Time' if 'Time' in da.dims else 'time'
    da = da.transpose(time_name, bin_dim)

    # restrict to used bins
    first_bin = int(da.attrs.get('FirstBin', 0))
    last_bin  = int(da.attrs.get('LastBin', da.sizes[bin_dim]-1))
    da = da.isel({bin_dim: slice(first_bin, last_bin+1)})
    nbins = da.sizes[bin_dim]

    # upper edges for used bins
    cells_all = np.asarray(da.attrs.get('CellSizes', []), dtype=float)
    if cells_all.size == 0:
        raise ValueError(f"{varname} missing CellSizes attr")
    cells_used = cells_all[first_bin:last_bin+1]  # length == nbins

    return da, bin_dim, time_name, cells_used, nbins

def _sum_range_by_upper_edge(da, bin_dim, cells_used, lower_um=None, upper_um=None):
    """Sum across bins chosen by upper-edge thresholds."""
    nbins = da.sizes[bin_dim]

    if lower_um is None and upper_um is None:
        raise ValueError("Provide at least lower_um or upper_um")

    # choose start
    if lower_um is None:
        i0 = 0
    else:
        i0 = int(np.searchsorted(cells_used, lower_um, side='left'))

    # choose end (inclusive)
    if upper_um is None:
        i1 = nbins - 1
    else:
        i1 = int(np.searchsorted(cells_used, upper_um, side='right')) - 1

    i0 = np.clip(i0, 0, nbins - 1)
    i1 = np.clip(i1, 0, nbins - 1)

    if i1 < i0:
        # empty selection → NaNs (shape preserves time axis)
        return da.isel({bin_dim: slice(0, 0)}).sum(dim=bin_dim) * np.nan

    return da.isel({bin_dim: slice(i0, i1 + 1)}).sum(dim=bin_dim, skipna=True)

def calc_concs_from_sd(sizedist_vars, nc):
    cols = []

    # --- C2DC branch (always compute Ndriz + Nprecip if present) ---
    var_2dc = next((v for v in sizedist_vars if v.startswith('C2DC')), None)
    if var_2dc is not None:
        da, bin_dim, time_name, cells_used, _ = _prep_probe(nc, var_2dc)
        # drizzle: 100–500 µm
        ndriz_2dc = _sum_range_by_upper_edge(da, bin_dim, cells_used, 100.0, 500.0)
        # precip: ≥1000 µm
        nprecip_2dc = _sum_range_by_upper_edge(da, bin_dim, cells_used, 1000.0, None)
        t = pd.to_datetime(da[time_name].values)
        df_2dc = pd.DataFrame({"Ndriz_2DC": ndriz_2dc.values, "Nprecip_2DC": nprecip_2dc.values},index=t)
        df_2dc.index.name = "time"
        cols.append(df_2dc)

    # --- C2DS branch (optional; only Nprecip requested) ---
    var_2ds = next((v for v in sizedist_vars if v.startswith('C2DS') and v.endswith('2H')), None)
    if var_2ds is not None:
        da, bin_dim, time_name, cells_used, _ = _prep_probe(nc, var_2ds)
        nprecip_2ds = _sum_range_by_upper_edge(da, bin_dim, cells_used, 1000.0, None)
        ndriz_2ds = _sum_range_by_upper_edge(da, bin_dim, cells_used, 100.0, 500.0)
        t = pd.to_datetime(da[time_name].values)
        df_2ds = pd.DataFrame({"Ndriz_2DS": ndriz_2ds, "Nprecip_2DS": nprecip_2ds.values}, index=t)
        df_2ds.index.name = "time"
        cols.append(df_2ds)

    # # --- C2DS branch (optional; only Nprecip requested) ---
    var_uhsas = next((v for v in sizedist_vars if v.startswith('CUH')), None)
    if var_uhsas is not None:
        da, bin_dim, time_name, cells_used, _ = _prep_probe(nc, var_uhsas)
        n_accum = _sum_range_by_upper_edge(da, bin_dim, cells_used, 100.0, None)
        n_ait = _sum_range_by_upper_edge(da, bin_dim, cells_used, 70.0, 100.0)
        t = pd.to_datetime(da[time_name].values)
        df_uhs = pd.DataFrame({"Naitk_UH": n_ait, "Naccum_UH": n_accum.values}, index=t)
        df_uhs.index.name = "time"
        cols.append(df_uhs)

    if not cols:
        return pd.DataFrame()

    # time-align and return
    return pd.concat(cols, axis=1).sort_index()

def load_flight_data(dir_path: str, idx: int = 0, add_sizedist: bool = True,
                     asof: bool = False, tol: str = "1s") -> pd.DataFrame:
    """
    High-level loader: base 1 Hz vars + (optional) drizzle/precip from sizedists.
    """    """
    High-Level Function for finding and reading in flight data.

    This function searches a directory for NetCDF (*.nc) flight data files, selects one based on the provided 
    index, opens it, identifies relevant variables based on the dataset contents, and reads those variables 
    into a Pandas DataFrame.

    :param dir_path: Path to the directory containing NetCDF flight data files.
    :param idx: Index of the file to load from the sorted list of *.nc files in the directory.

    :return: A Pandas DataFrame containing the extracted flight data variables.
    """
    flight_dat_paths = find_flight_fnames(dir_path)
    nc = open_nc(flight_dat_paths[idx])
    vars2read = read_vars(nc)
    df = read_flight_nc(nc,vars2read)
    
    # 2) derived from size distributions (returns time-indexed DF)
    sd_vars = read_sizedist_vars(nc)          # your prefix-based picker
    conc_df = calc_concs_from_sd(sd_vars, nc) # columns like Ndriz_2DC, Nprecip_2DC, Nprecip_2DS, ...

    if conc_df is None or conc_df.empty:
        return df

    # 3) join on time
    df2 = df.copy()
    df2["Time"] = pd.to_datetime(df2["Time"]).dt.tz_localize(None).round("S")

    conc = conc_df.copy()
    conc.index = pd.to_datetime(conc.index).tz_localize(None).round("S")

    if asof:
        # nearest match within tolerance (useful if clocks are off by <1s)
        out = pd.merge_asof(
            df2.sort_values("Time"),
            conc.reset_index().rename(columns={"index": "Time"}).sort_values("Time"),
            on="Time",
            direction="nearest",
            tolerance=pd.Timedelta(tol),
        )
    else:
        # exact join on second-resolution timestamps
        out = df2.set_index("Time").join(conc, how="left").reset_index()

    return out


def find_sondes(dir_path: str) -> list[str]:
    """
    find_flight_fnames just searches a directory for all *.nc files and returns a list of them.

    :param dir_path: a path to the directory containing flight netcdf files

    :return: Returns a list of flight netcdf files.
    """
    flight_paths=[]
    flight_fnames = sorted([fname for fname in os.listdir(dir_path) if fnmatch(fname, "*.cls")])
    for i in range(len(flight_fnames)):
        flight_paths.append(dir_path + '/' + flight_fnames[i])
    
    return flight_paths

def read_sonde2df(file_path):
    """
    Reads a `.cls` radiosonde file and extracts multiple datasets with their nominal release times.

    :param file_path: Path to the `.cls` file containing radiosonde data.

    :return: A tuple containing:
        - A list of Pandas DataFrames, each representing an individual radiosonde dataset.
        - A list of corresponding nominal release times as Pandas Timestamps.
    
    :raises FileNotFoundError: If the file does not exist.
    :raises ValueError: If the file is incorrectly formatted or missing essential data.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File '{file_path}' not found.")

    # Read the entire file into memory
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    # Identify all "Nominal Release Time" occurrences
    start_indices = [i for i, line in enumerate(lines) if "Nominal Release Time" in line]
    
    if not start_indices:
        raise ValueError(f"No 'Nominal Release Time' entries found in file: {file_path}")

    datasets = []
    drop_times = []
    # Process each radiosonde dataset
    for idx, start in enumerate(start_indices):
        # Find the start of the tabular data
        data_start = None
        for i in range(start, len(lines) - 2):
            if lines[i].strip().startswith("Time") and "Press" in lines[i]:  # Detect header row
                data_start = i + 3  # Data starts 2 lines after the column names
                break

        if data_start is None:
            print(f"Warning: No data start found for entry at line {start}")
            continue  # Skip this dataset

        # Extract and convert nominal release time
        try:
            date_time_str = lines[start].split("):")[1].strip()
            drop_time = pd.to_datetime(date_time_str, format='%Y, %m, %d, %H:%M:%S')
        except (IndexError, ValueError) as e:
            print(f"Warning: Failed to parse drop time at line {start}: {e}")
            continue  # Skip this dataset

        drop_times.append(drop_time)

        # Extract column names
        columns = lines[data_start - 3].strip().split()

        # Determine dataset end (next "Nominal Release Time" or end of file)
        end = start_indices[idx + 1] if idx + 1 < len(start_indices) else len(lines)

        # Extract and clean data
        data_lines = [line.strip().split() for line in lines[data_start:end]]
        data = [row for row in data_lines if len(row) == len(columns)]

        # Convert to DataFrame
        df = pd.DataFrame(data, columns=columns)

        if df.empty:
            print(f"Warning: Empty dataset at line {start}")
            continue  # Skip empty datasets

        # Remove rows containing 9999.0 in any column
        df = df[(df != "9999.0").all(axis=1)]
        # Convert numeric columns where possible
        df = df.apply(pd.to_numeric, errors='coerce')
        # Store drop time in DataFrame metadata
        df.attrs["drop_time"] = drop_time
        # Append dataset
        datasets.append(df)

    return datasets

def load_nc_cldrgme(file_paths):

    combined_blocks = []   
    for path in file_paths:
        rf = path.split("/")[-1].split(".")[0]  # e.g., "RF01"
        ds = xr.open_dataset(path)
        
        # Get unique combinations
        labels = ds["block_label"].values
        indices = ds["block_index"].values
        
        # Convert to DataFrame for convenient filtering
        df = ds.to_dataframe().reset_index().drop(columns="index")  # remove redundant index
        
        # Get all unique (label, index) pairs
        unique_blocks = df[["block_label", "block_index"]].drop_duplicates()
        
        # Loop through each unique block
        for _, row in unique_blocks.iterrows():
            label = row["block_label"]
            idx = row["block_index"]
        
            # Filter DataFrame
            df_block = df[(df["block_label"] == label) & (df["block_index"] == idx)].copy()
        
            # (Optional) add flight ID if you have it
            df_block["block_label"] = label
            df_block["block_index"] = idx
        
            combined_blocks.append(df_block)
        all_blocks = pd.concat(combined_blocks, ignore_index=True)
    return all_blocks




import pathlib as path
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation  


def grid_flight(cesm: xr.open_dataset, cesm_dat: xr.open_dataset, df: pd.DataFrame) -> dict:    
    
    # Step 1: Identify Variables Automatically
    lat_var = next((var for var in df.columns if 'GGLAT' in var), None)
    lon_var = next((var for var in df.columns if 'GGLON' in var), None)
    alt_var = next((var for var in df.columns if 'GGALT' in var or 'PSXC' in var), None)

    # Compute pressure altitude (palt) from CESM hybrid coordinates
    p0 = cesm.P0  # Reference pressure
    ps = cesm.PS  # Surface pressure [=] Pa
    hyai = cesm.hyai  # Hybrid A coefficient at layer interface
    hybi = cesm.hybi  # Hybrid B coefficient at layer interface
    
    P_dummy = p0 * hyai
    midP = (P_dummy + (hybi * ps))  # [=] Pa
    palt = midP * 0.01  # Convert to hPa
    
    df_vars = [col for col in df.columns if col.lower() != 'time']
    
    if not lat_var or not lon_var or not alt_var:
        raise ValueError("Missing essential latitude, longitude, or altitude variables.")
    
    # Step 2: Create a 3D grid based on CESM & flight data
    # Compute midpoints
    lat_mp = (cesm.lat[:-1] + cesm.lat[1:]) / 2
    lon_mp = (cesm.lon[:-1] + cesm.lon[1:]) / 2
    # Find lat & lon bounds based on aircraft min/max values
    lat_bounds = [int(np.abs(lat_mp - np.min(df[lat_var])).argmin()) - 1, 
                  int(np.abs(lat_mp - np.max(df[lat_var])).argmin()) + 1]
    lon_bounds = [int(np.abs(lon_mp - np.min(df[lon_var])).argmin()) - 1, 
                  int(np.abs(lon_mp - np.max(df[lon_var])).argmin()) + 1]

    # Ensure bounds are within valid range
    lat_bounds = [max(0, lat_bounds[0]), min(len(lat_mp) - 1, lat_bounds[1])]
    lon_bounds = [max(0, lon_bounds[0]), min(len(lon_mp) - 1, lon_bounds[1])]
    
    # Select the subset of `palt` corresponding to the lat/lon bounds
    palt_subset = palt.isel(lat=slice(lat_bounds[0], lat_bounds[1] + 1),
                            lon=slice(lon_bounds[0], lon_bounds[1] + 1))

    # Extract altitude values based on aircraft altitude range
    min_alt, max_alt = np.min(df[alt_var]), np.max(df[alt_var])

    # Find altitude bounds dynamically for each lat-lon grid cell
    alt_indices = []
    for lat_idx in range(lat_bounds[0], lat_bounds[1] + 1):
        for lon_idx in range(lon_bounds[0], lon_bounds[1] + 1):
            local_palt = palt.isel(lat=lat_idx, lon=lon_idx).values  # Get 1D alt profile for this grid cell

            min_alt_idx = np.abs(local_palt - min_alt).argmin()
            max_alt_idx = np.abs(local_palt - max_alt).argmin()

            alt_indices.append((min_alt_idx, max_alt_idx))

    # Determine overall altitude bounds from all lat-lon cells
    min_alt_bound = min(i[0] for i in alt_indices)
    max_alt_bound = max(i[1] for i in alt_indices)
    
    # Ensure altitude bounds are valid
    alt_bounds = [max(0, min_alt_bound - 1), min(palt.shape[0] - 1, max_alt_bound + 1)]
    # Define bounds dictionary
    bounds = {'lat': lat_bounds, 'lon': lon_bounds, 'palt': alt_bounds}

    # Create grid with the correct shape
    grid_shape = (alt_bounds[1] - alt_bounds[0] + 1, 
                  lat_bounds[1] - lat_bounds[0] + 1, 
                  lon_bounds[1] - lon_bounds[0] + 1)
        
        # return np.zeros(grid_shape), bounds
    grid = np.zeros(grid_shape)

    # Step 3: Match aircraft times with CESM times
    da = xr.DataArray(cesm_dat.time, dims="time")
    cesm_times = pd.to_datetime([pd.Timestamp(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
                        for dt in da.values])    
    aircraft_times = pd.to_datetime(df['Time'])
    
    times = np.array(cesm_times[np.isin(cesm_times, aircraft_times)])
    
    # Step 4: Initialize grid arrays
    mean_lat, mean_lon, mean_alt = np.zeros_like(grid, dtype=float), np.zeros_like(grid, dtype=float), np.zeros_like(grid, dtype=float)
    
    mean_values = {var: np.zeros_like(grid, dtype=float) for var in df_vars}
    mean_lat, mean_lon, mean_alt = np.zeros_like(grid, dtype=float), np.zeros_like(grid, dtype=float), np.zeros_like(grid, dtype=float)
    
    # Generate latitude, longitude grid values
    lats = np.array(lat_mp[bounds['lat'][0]:bounds['lat'][1] + 1])
    lons = np.array(lon_mp[bounds['lon'][0]:bounds['lon'][1] + 1])
    
    # Select the region of interest in `palt`
    palt_subset = palt.sel(lat=slice(lats.min(), lats.max()), lon=slice(lons.min(), lons.max()))
    # Adjust time range
    new_time = times[-1] + np.timedelta64(30, 'm')
    minus_time = times[0] - np.timedelta64(30, 'm')
    times = np.append(minus_time, np.append(times, new_time))
    
    mid_time = []
    for t in range(0, len(times) - 1):
        time_start, time_end = times[t], times[t + 1]
    
        # Select aircraft data within the time interval
        air_time_indices = (df['Time'] > time_start) & (df['Time'] <= time_end)
        sliced_df_time = df[air_time_indices]
    
        if not sliced_df_time.empty:
            # Digitize lat & lon into grid bins
            lat_bins = np.digitize(sliced_df_time['GGLAT'], lats) - 1
            lon_bins = np.digitize(sliced_df_time['GGLON'], lons) - 1
    
            # Ensure valid lat/lon indices
            valid_mask = (lat_bins >= 0) & (lat_bins < len(lats) - 1) & \
                         (lon_bins >= 0) & (lon_bins < len(lons) - 1)
    
            if valid_mask.any():
                # Filter valid rows
                sliced_df = sliced_df_time.loc[valid_mask].copy()
                sliced_df['lat_bin'] = lat_bins[valid_mask]
                sliced_df['lon_bin'] = lon_bins[valid_mask]
                # Compute altitude bins dynamically based on (lat_bin, lon_bin)
                alt_bins = np.full(len(sliced_df), -1, dtype=int)  # Initialize invalid bins

                for i, (lat_idx, lon_idx, psxc) in enumerate(zip(sliced_df['lat_bin'], sliced_df['lon_bin'], sliced_df['PSXC'])):
                    # Extract altitude levels for this (lat_idx, lon_idx) from `palt`
                    local_alt_profile = palt_subset.isel(lat=lat_idx, lon=lon_idx).values  # Keep natural order
                    # Ensure the pressure profile is in the correct shape
                    local_alt_profile = local_alt_profile.squeeze()
                
                    # Use np.searchsorted instead of np.digitize to find the correct bin
                    alt_bins[i] = np.searchsorted(local_alt_profile, psxc, side='right') - 1

                # Assign altitude bins
                sliced_df['alt_bin'] = alt_bins
    
                # Ensure valid altitude indices
                valid_alt_mask = (alt_bins >= 0) & (alt_bins < palt.shape[0] - 1)
                sliced_df = sliced_df[valid_alt_mask]
    
                # Group by grid cells and compute means                
                grouped = sliced_df.groupby(['alt_bin', 'lat_bin', 'lon_bin'])
    
                for var in df_vars:
                    grouped_mean = grouped[var].mean()
                    if not grouped_mean.empty:  # ✅ Skip empty bins
                        mean_values[var][tuple(zip(*grouped_mean.index.to_numpy()))] = grouped_mean.values
    
                        mean_lat[tuple(zip(*grouped[lat_var].mean().index.to_numpy()))] = grouped[lat_var].mean().values
                        mean_lon[tuple(zip(*grouped[lon_var].mean().index.to_numpy()))] = grouped[lon_var].mean().values
                        mean_alt[tuple(zip(*grouped[alt_var].mean().index.to_numpy()))] = grouped[alt_var].mean().values
    
                # Compute mid-time for each grid cell
                grouped_times = grouped['Time'].agg(lambda x: x.min() + (x.max() - x.min()) / 2)
                mid_time.extend(grouped_times.values)
    
    # Convert mid_time to NumPy array
    mid_time = np.array(mid_time, dtype='datetime64[ns]')
    
    # Filter valid data points
    mean_t = mean_values['ATX']
    valid_indices = np.argwhere(mean_t != 0)  # (N, 3) array of (alt, lat, lon)
    
    alt_indices, lat_indices, lon_indices = valid_indices[:, 0], valid_indices[:, 1], valid_indices[:, 2]
    
    selected_time = mid_time[:len(valid_indices)]
    
    grid_dict = {
        'Time': selected_time,
        'Latitude': mean_lat[mean_t != 0],
        'Longitude': mean_lon[mean_t != 0],
        'Altitude': mean_alt[mean_t != 0],
    }
    grid_dict.update({var: mean_values[var][mean_t != 0] for var in df_vars})
    
    # Sort the data by time
    sorted_indices = np.argsort(selected_time)
    grid_dict = {key: np.array(value)[sorted_indices] for key, value in grid_dict.items()}
    grid_dict.update({'long_names': df.attrs['long_names']})
    
    # Confirm data integrity
    if all(len(v) > 0 for v in grid_dict.values()):
        print("✅ Grid dictionary successfully populated with data!")
    else:
        print("⚠️ Warning: Some entries in grid_dict are empty!")
    
    return grid_dict, grid, bounds


def plot_3d_track(grid_data,df):

    # Create a figure
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    sc = ax.scatter(grid_data['GGLON'], grid_data['GGLAT'], grid_data['PSXC'], c=grid_data['ATX'], cmap='viridis', marker='^',label='grid-mean values',s=100)
    # Invert the Z-axis
    ax.scatter(df['GGLON'], df['GGLAT'],df['PSXC'], c=df['ATX'], label='3D Flight track',s=12)
    ax.invert_zaxis()

    ax.set_xlabel('latitude (deg)')
    ax.set_ylabel('longitude (deg)')
    ax.set_zlabel('pressure alt (hPa)') 

    ax.legend()
    # # Color bar to show the mapping of color to the fourth dimension
    plt.colorbar(sc, label='Mean Temperature (°C)')

    # Animation function to rotate the view
    def rotate(angle):
        ax.view_init(elev=30, azim=angle)

    # Create animation
    ani = FuncAnimation(fig, rotate, frames=np.arange(-180, 360, 20), interval=100)

    # Show the animation in Jupyter Notebook
    from IPython.display import HTML
    HTML(ani.to_jshtml())

    plt.show()

def write_nc(grid_data, filename="test_grid_data.nc"):
    """
    Automatically creates and saves a NetCDF file from the given grid data dictionary.
    
    :param grid_data: Dictionary containing time series data with "Time" and corresponding variables.
    :param filename: Name of the NetCDF file to be saved.
    """

    # Extract headers dynamically (excluding "Time")
    headers = [key for key in grid_data.keys() if key.lower() != "time"]

    # Create the xarray dataset dynamically
    ds = xr.Dataset(
        {var: (["time"], grid_data[var]) for var in headers},  # Assign all variables dynamically
        coords={"time": grid_data["Time"]},  # Set "Time" as the coordinate
    )
    # # Save to a NetCDF file
    ds.to_netcdf("test_grid_data.nc")
    print("NetCDF file 'grid_data.nc' saved successfully!")


import pandas as pd
import numpy as np

import glob
import xarray as xr
import datetime
from scipy.spatial import cKDTree

def assign_flight_type(df):
    """
    Assigns flight type ('level' or 'profile') to each row of the input DataFrame based on stable altitude blocks 
    and gaps between these blocks. The function uses rolling standard deviation of altitude to identify level legs
    and combines consecutive blocks of stable altitude with a specified time gap threshold. Additionally, it labels 
    flight segments as "level" for level legs and "profile" for the aircraft vertical profile.

    Parameters:
    -----------
    df : pandas.DataFrame
        The input DataFrame with at least the following columns:
        - 'Time' (timestamp)
        - 'GGALT' (altitude in meters)

    Returns:
    --------
    pandas.DataFrame
        The input DataFrame with a new column 'flight_type', where each row is assigned a flight type:
        - 'level' for stable altitude periods
        - 'profile' for gaps between stable altitude blocks

    Example:
    --------
    df = pd.read_csv('flight_data.csv')  # Assuming the CSV contains relevant columns
    df_with_flight_types = assign_flight_type(df)
    """

    #-----------------------------------------
    #----- Find profiles and level legs ------
    #-----------------------------------------
    
    # Define a time gap threshold to combine blocks (e.g., 120 seconds)
    time_gap_threshold = pd.Timedelta(seconds=120)
    
    # Compute rolling standard deviation of altitude to smooth noise
    df['rolling_std'] = df['GGALT'].rolling(window=10, center=True).std()
    
    # Identify where altitude remains stable within the threshold
    df['stable'] = df['rolling_std'] < 3  # You can adjust the threshold (meters)
    
    # Assign unique block IDs when stability changes
    df['block_id'] = (df['stable'] != df['stable'].shift()).cumsum()
    
    # Group by block_id and filter for long-duration stable blocks
    block_info = df[df['stable']].groupby('block_id').agg(
        start_time=('Time', 'first'),
        end_time=('Time', 'last'),
        lower_bound=('GGALT', 'min'),  # Minimum altitude (lower bound)
        upper_bound=('GGALT', 'max'),  # Maximum altitude (upper bound)
        duration=('Time', lambda x: x.max() - x.min())
    )
    
    # Filter out short-duration blocks
    valid_blocks = block_info[block_info['duration'] > pd.Timedelta(seconds=150)] ## EDIT?
    
    # Sort the blocks by start time
    valid_blocks = valid_blocks.sort_values(by='start_time')
    
    # Define a time gap threshold to combine blocks (e.g., 120 seconds)
    time_gap_threshold = pd.Timedelta(seconds=120)
    
    # Combine consecutive blocks that are less than the threshold apart
    combined_blocks = []
    previous_block = valid_blocks.iloc[0]
    
    for idx, current_block in valid_blocks.iloc[1:].iterrows():
        # Check if the gap between the end time of the previous block and start time of the current block is below the threshold
        if current_block['start_time'] - previous_block['end_time'] <= time_gap_threshold:
            # Extend the previous block's end time to the current block's end time
            previous_block['end_time'] = current_block['end_time']
        else:
            # If the gap is too large, append the previous block and update to the current block
            combined_blocks.append(previous_block)
            previous_block = current_block
    
    # Add the last block after the loop
    combined_blocks.append(previous_block)
    
    # Convert combined blocks back to DataFrame
    combined_blocks_df = pd.DataFrame(combined_blocks)
    
    # --- Identify and Label "Profiles" between "Level" (Stable) sections ---
    
    # Create a new column 'flight_type' to categorize the blocks as "level" or "profile"
    combined_blocks_df['flight_type'] = 'level'  # By default, label as 'level'
    
    # Now identify the gaps between "level" blocks and label as "profile"
    profile_blocks = []
    for i in range(len(combined_blocks_df) - 1):
        end_time_current = combined_blocks_df.iloc[i]['end_time']
        start_time_next = combined_blocks_df.iloc[i + 1]['start_time']
        
        # If there's a gap between two 'level' blocks, label the gap as 'profile'
        if start_time_next - end_time_current > time_gap_threshold:
            # Assign 'profile' to the gap between two level blocks and calculate duration
            profile_duration = start_time_next - end_time_current  # Duration of the profile block
            
            profile_blocks.append({
                'start_time': end_time_current,
                'end_time': start_time_next,
                'flight_type': 'profile',
                'duration': profile_duration
            })
    
    # Convert 'profile_blocks' to DataFrame
    profile_blocks_df = pd.DataFrame(profile_blocks)
    
    # Append profile blocks to the original combined blocks DataFrame
    combined_blocks_with_profiles = pd.concat([combined_blocks_df, profile_blocks_df], ignore_index=True)
    
    # Sort again by time
    combined_blocks_with_profiles = combined_blocks_with_profiles.sort_values(by='start_time')
    
    # Check for "Profile" after the Last Level Block
    last_end_time = combined_blocks_with_profiles.iloc[-1]['end_time']
    last_time_in_data = df['Time'].max()
    
    if last_time_in_data - last_end_time > time_gap_threshold:
        # If the gap is greater than the threshold, consider it a "profile" block
        profile_block = pd.DataFrame([{
            'start_time': last_end_time,
            'end_time': last_time_in_data,
            'flight_type': 'profile',
        }])
    
        # Concatenate the new profile block to the existing DataFrame
        combined_blocks_with_profiles = pd.concat([combined_blocks_with_profiles, profile_block], ignore_index=True)
    
    # Check for "Profile" before the First Level Block, used for takeoff/landing
    first_start_time = combined_blocks_with_profiles.iloc[0]['start_time']
    first_time_in_data = df['Time'].min()
    
    if first_start_time - first_time_in_data > time_gap_threshold:
        profile_block_before_first = pd.DataFrame([{
            'start_time': first_time_in_data,
            'end_time': first_start_time,
            'flight_type': 'profile',
        }])
    
        combined_blocks_with_profiles = pd.concat([profile_block_before_first, combined_blocks_with_profiles], ignore_index=True)
    
    # List of columns to remove
    columns_to_remove = ['rolling_std','stable','block_id']
    # Drop the specified columns from df2
    df = df.drop(columns=columns_to_remove)
    
    # Add new column "flight_type" as either "level" or "profile"
    for _, row in combined_blocks_with_profiles.iterrows():
        flight_type = row['flight_type']
        # Find rows in df2 where the time is between start_time and end_time
        mask = (df['Time'] >= row['start_time']) & (df['Time'] <= row['end_time'])
        df.loc[mask, 'flight_type'] = flight_type
    # Assign the flight_type to the first few rows that fall before the first start_time in df1
    df.loc[df['Time'] < first_start_time, 'flight_type'] = df.iloc[0]['flight_type']

    #------------------------------
    #----- Find cloud layers ------
    #------------------------------
    # Ensure 'Time' is in datetime format
    df = df.copy()  # Avoid modifying original DataFrame
    df['Time'] = pd.to_datetime(df['Time'])
    
    # Find best match for column names dynamically
    plwc_col = next((col for col in df.columns if 'PLWCD' in col), None) or \
           next((col for col in df.columns if 'PLWC' in col), None)
    concd_col = next((col for col in df.columns if 'CONCD' in col), None)
    # Add check if there are any cloudy periods
    if not plwc_col or not concd_col:
        print("Required columns not found. Skipping cloud detection.")
        final_cloud_blocks = pd.DataFrame(columns=['start_time', 'end_time', 'lower_bound', 'upper_bound', 'duration', 'Location'])
        df['cloud_status'] = 'Out-of-cloud'
        df['Location'] = 'Free'
    else:
        df['blocked'] = (df[plwc_col] > 0.001) & (df[concd_col] > 10)
        if not df['blocked'].any():
            print("No valid cloud blocks found. Skipping cloud layer logic.")
            final_cloud_blocks = pd.DataFrame(columns=['start_time', 'end_time', 'lower_bound', 'upper_bound', 'duration', 'Location'])
            df['cloud_status'] = 'Out-of-cloud'
            df['Location'] = 'Free'
        else:
            df['block_id'] = (df['blocked'] != df['blocked'].shift()).cumsum()
            block_info = df[df['blocked']].groupby('block_id').agg(
                start_time=('Time', 'first'),
                end_time=('Time', 'last'),
                lower_bound=('GGALT', 'min'),
                upper_bound=('GGALT', 'max'),
            )

            # Calculate duration directly by subtracting start_time from end_time
            block_info['duration'] = block_info['end_time'] - block_info['start_time']
            
            # Filter out short-duration blocks
            min_vertical = 30  # Adjust as needed (100 meters in your case)
            valid_blocks = block_info[(block_info['upper_bound'] - block_info['lower_bound']) > min_vertical].reset_index(drop=True)
            
            # Define the altitude difference and time gap thresholds
            altitude_gap_threshold = 200  # Increased altitude gap threshold
            # time_gap_threshold = pd.Timedelta(minutes=20)  # Time gap threshold for merging
            
            # Sort the valid blocks by their start time to process them in sequence
            valid_blocks = valid_blocks.sort_values(by='lower_bound')
            
            # Initialize a list to store combined blocks
            combined_blocks = []
            previous_block = valid_blocks.iloc[0].to_dict()
            
            # Iterate through the blocks and merge those that are within the thresholds
            for idx, current_block in valid_blocks.iloc[1:].iterrows():
                # Calculate the altitude gap between the current block's lower bound and the previous block's upper bound
                altitude_gap = abs(current_block['lower_bound'] - previous_block['upper_bound'])
                
                # Calculate the time gap between the current block's start time and the previous block's end time
                time_gap = current_block['start_time'] - previous_block['end_time']
                # Check if the altitude gap is within the threshold or if the time gap is within the allowed range for smaller altitudes
                if altitude_gap <= altitude_gap_threshold :
                    # If both criteria are met, merge the blocks
                    previous_block['end_time'] = max(previous_block['end_time'], current_block['end_time'])  # Get the latest end time
                    previous_block['start_time'] = min(previous_block['start_time'], current_block['start_time'])  # Get the earliest start time
                    previous_block['upper_bound'] = max(previous_block['upper_bound'], current_block['upper_bound'])  # Update upper bound
                    previous_block['lower_bound'] = min(previous_block['lower_bound'], current_block['lower_bound'])  # Update lower bound
                    
                    # Recalculate the duration for the merged block
                    previous_block['duration'] = previous_block['end_time'] - previous_block['start_time']
                else:
                    # If the blocks are far apart, save the previous block and move to the next one
                    combined_blocks.append(previous_block)
                    previous_block = current_block.to_dict()
            
            # Add the last block after the loop
            combined_blocks.append(previous_block)
            
            # Convert the merged blocks back into a DataFrame
            combined_blocks_df = pd.DataFrame(combined_blocks)
            
            # Second check for merging adjacent blocks in combined_blocks_df
            final_combined_blocks = []
            previous_block = combined_blocks_df.iloc[0].to_dict()
            
            # Apply additional check for merging based on both time and altitude gap
            for idx, current_block in combined_blocks_df.iloc[1:].iterrows():
                # Calculate the altitude gap and time gap
                altitude_gap = abs(current_block['lower_bound'] - previous_block['upper_bound'])
                time_gap = current_block['start_time'] - previous_block['end_time']
                
                # Check for overlap in the altitude ranges
                overlap_check = (current_block['lower_bound'] >= previous_block['lower_bound']) and (current_block['lower_bound'] <= previous_block['upper_bound'])
            
                # Check if both the altitude gap, time gap, or overlap condition is met
                if altitude_gap <= altitude_gap_threshold or overlap_check:
                    # Merge the blocks
                    previous_block['end_time'] = max(previous_block['end_time'], current_block['end_time'])
                    previous_block['start_time'] = min(previous_block['start_time'], current_block['start_time'])
                    previous_block['upper_bound'] = max(previous_block['upper_bound'], current_block['upper_bound'])
                    previous_block['lower_bound'] = min(previous_block['lower_bound'], current_block['lower_bound'])
                    
                    # Recalculate the duration for the merged block
                    previous_block['duration'] = previous_block['end_time'] - previous_block['start_time']
                else:
                    # Save the previous block and move to the next one
                    final_combined_blocks.append(previous_block)
                    previous_block = current_block.to_dict()
            
            # Add the last block after the loop
            final_combined_blocks.append(previous_block)
            
            # Convert the final combined blocks back into a DataFrame
            final_cloud_blocks = pd.DataFrame(final_combined_blocks)

            # Add 'cloud_status' based on whether altitude and time fall within any blocked region (in the cloud or out of cloud)
            df['cloud_status'] = 'Out-of-cloud'  # Default label
            # Loop through each block and label altitudes as "In-cloud" if they fall within the block's range
            for _, block in final_cloud_blocks.iterrows():
                # Create a mask that checks both altitude and time conditions
                mask = (
                    (df['GGALT'] >= block['lower_bound']) & (df['GGALT'] <= block['upper_bound']) &
                    (df['Time'] >= block['start_time']) & (df['Time'] <= block['end_time'])
                )
                
                # Apply the 'In-cloud' label where the mask is True
                df.loc[mask, 'cloud_status'] = 'In-cloud'
            
            # List of columns to remove
            columns_to_remove = ['blocked','block_id']
            # Drop the specified columns from df2
            df = df.drop(columns=columns_to_remove)
        
            df['Location'] = 'Free'
            
            # Find the minimum in-cloud altitude
            min_ic_alt = np.min(final_cloud_blocks['lower_bound'])-5
            mask = df.GGALT < min_ic_alt
            # Define 
            df.loc[mask, 'Location'] = 'BL'
        
            # Update the Location column based on the GGALT and cloud status
            df.loc[df['GGALT'] < min_ic_alt, 'Location'] = 'BL'
        
            # --- Add Location to final_cloud_blocks DataFrame ---
            # Add 'Location' based on the minimum in-cloud altitude
            final_cloud_blocks['Location'] = final_cloud_blocks['lower_bound'].apply(
                lambda x: 'BL' if x < min_ic_alt else 'Free'
            )
                # Add 'Location' based on the minimum in-cloud altitude
            combined_blocks_with_profiles['Location'] = combined_blocks_with_profiles['lower_bound'].apply(
                lambda x: 'BL' if x < min_ic_alt else 'Free'
            )

    # Sort the dataframe by Time for continuous time grouping
    df = df.sort_values(by='Time')
    # Remove rows where 'flight_type' is NaN
    df = df.dropna(subset=['flight_type'])
    # Create a new column 'block_id' to group continuous time periods based on flight_type, cloud_status, and Location
    df['block_id'] = (df['flight_type'] != df['flight_type'].shift()) | \
                      (df['cloud_status'] != df['cloud_status'].shift()) | \
                      (df['Location'] != df['Location'].shift())
    df['block_id'] = df['block_id'].cumsum()

    Final_ds = {'DataFrame': df,
                'flight_blocks': combined_blocks_with_profiles,
                'Cloud_blocks': final_cloud_blocks
               }
    
    return Final_ds

def block_flight(df):
    """
    Segments a flight dataset into different flight block categories based on cloud status, location, and flight type.

    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame containing flight data with at least the following columns:
        - 'block_id' (int): Identifies different flight segments.
        - 'Location' (str): Can be 'BL' (Boundary Layer) or 'Free' airspace.
        - 'flight_type' (str): Can be 'level' or 'profile'.
        - 'cloud_status' (str): Either 'In-cloud' or 'Out-of-cloud'.
        - 'GGALT' (float): Altitude dbata, used for filtering profile segments.
        - 'Time' (datetime): Used to filter level flight segments.

    Returns:
    --------
    Flight_blocks : dict
        A dictionary containing categorized flight data:
        - 'Level BL': List of DataFrames for level flight in the boundary layer.
        - 'In-Cloud Profiles': List of DataFrames for in-cloud profile flights with altitude variation > 30m.
        - 'In-Cloud Level FT': List of DataFrames for level flights in free airspace within clouds.
        - 'Out-of-cloud Level FT': List of DataFrames for level flights in free airspace, lasting at least 3 minutes.

    Notes:
    ------
    - The function removes the first and last 'Level BL' periods to exclude takeoff/landing effects.
    - Only level flights lasting more than 180 seconds are included in 'Out-of-cloud Level FT'.
    - Profile flights are only included if their altitude change is greater than 30 meters.
    """
    # ---------- Level BL periods (KEEP ALL) ----------
    out_of_cloud_bl = df[(df['Location'] == 'BL') & (df['flight_type'] == 'level')]
    bl_ids = sorted(out_of_cloud_bl['block_id'].unique())
    bl_blocks_ds = [df[df['block_id'] == i] for i in bl_ids]
    
    # Find In-Cloud profile periods
    in_cloud_prof = df[(df['cloud_status'] == 'In-cloud') & (df['flight_type'] == 'profile')]
    ic_prof_ids = sorted(in_cloud_prof['block_id'].unique())
    ic_pro_blocks_ds = [df[df['block_id'] == i] for i in ic_prof_ids if df[df['block_id'] == i]['GGALT'].max() - df[df['block_id'] == i]['GGALT'].min() > 30]
    
    # Find level FT periods out-of-cloud
    level_ft = df[(df['cloud_status'] == 'Out-of-cloud') & (df['flight_type'] == 'level') & (df['Location'] == 'Free')]
    out_of_cloud_ft_ids = sorted(level_ft['block_id'].unique())
    level_ft_out_blocks_ds = [df[df['block_id'] == i] for i in out_of_cloud_ft_ids if df[df['block_id'] == i]['Time'].iloc[-1] - df[df['block_id'] == i]['Time'].iloc[0] > pd.Timedelta(seconds=180)]

    # Find level FT periods in-cloud
    level_ft_ic = df[(df['cloud_status'] == 'In-cloud') & (df['flight_type'] == 'level') & (df['Location'] == 'Free')]
    in_cloud_ft_ids = sorted(level_ft_ic['block_id'].unique())
    level_ft_ic_blocks_ds = [df[df['block_id'] == i] for i in in_cloud_ft_ids]
    
    # Save blocks of flight as dictionary for output
    Flight_blocks = {
        'Level BL': bl_blocks_ds,
        'In-Cloud Profiles': ic_pro_blocks_ds,
        'In-Cloud Level FT': level_ft_ic_blocks_ds,
        'Out-of-cloud Level FT': level_ft_out_blocks_ds
    }

    return Flight_blocks

def assign_cloud_type_HCR(flight_blocks, dir, idx: int = 0):
    """
    Assigns cloud echo type classifications from HCR (HIAPER-Cloud Radar) data 
    to flight data blocks within the global Flight_blocks variable.

    Parameters:
    -----------
    dir : str
        Base directory path where the RF (Research Flight) subfolders are located.
    idx : int, optional
        Index of the flight number (used to construct the folder name as 'RF{idx}'), 
        by default 0.

    Returns:
    --------
    dict
        Updated Flight_blocks dictionary with an added 'Echo_Type' column in each block, 
        indicating the radar-derived cloud classification at each time step.

    Notes:
    ------
    - Requires global Flight_blocks to be defined externally.
    - Each time in the flight data is matched with HCR timestamps to assign echo types.
    - Echo type values are pulled from the 'HCR_ECHO_TYPE_1D' variable in netCDF files.
    """

    dir_fold = dir + f"RF{idx+1:02d}" + '/'
    flight_paths = find_flight_fnames(dir_fold)

    hcr_time = []
    echo_type_1D = []
    
    for file in flight_paths:
        nc = open_nc(file)
        hcr_time.extend(np.array(nc.time))
        echo_type_1D.append(np.array(nc.HCR_ECHO_TYPE_1D))
    
    echo_type_1D = np.concatenate(echo_type_1D)
    
    for val in flight_blocks:
        block_type = flight_blocks[val]
        for i in range(len(block_type)):
            # block = block_type[i]
            block = block_type[i].copy()  # Make a copy to avoid SettingWithCopyWarning
            # Extract start and end time from the block
            start_time = block['Time'].iloc[0]
            end_time = block['Time'].iloc[-1]
            
            # Convert start_time and end_time to the same format as hcr_time if necessary
            # Assuming hcr_time is a list of datetime objects, convert if needed
            if isinstance(start_time, str):  # If start_time is in string format, convert it
                start_time = pd.to_datetime(start_time)
            if isinstance(end_time, str):  # If end_time is in string format, convert it
                end_time = pd.to_datetime(end_time)
    
            hcr_time_array = np.array(hcr_time)
            echo_column = np.full(len(block), np.nan)
    
            for j in range(len(block)):
                time_point = pd.to_datetime(block['Time'].iloc[j])
                match_indices = np.where(hcr_time_array == time_point)[0]
                if len(match_indices) > 0:
                    echo_column[j] = echo_type_1D[match_indices[0]]
    
            block.loc[:, 'Echo_Type'] = echo_column  # <- Use loc for safe assignment
            block_type[i] = block
    
        flight_blocks[val] = block_type
        
    return flight_blocks


# High-Level function
def VAP_process_flight_data(df,i):
    """
    High-Level Function for Processing Flight Data in Value Added Products.

    This function serves as the main entry point for processing flight data. It first calls 
    `assign_flight_type` to assign flight types (e.g., 'level' or 'profile') to different segments of the flight 
    based on altitude stability and time gaps. After flight types are assigned, it proceeds to categorize the data 
    into different flight blocks (e.g., level flight in boundary layer, in-cloud profile flight, etc.) by calling 
    the `block_flight` function.

    Parameters:
    -----------
    df : pandas.DataFrame
        A DataFrame containing flight data with at least the following columns:
        - 'Time' (datetime): Time of each flight record.
        - 'GGALT' (float): Altitude of the aircraft.
        - 'PLWCD_' (float): Cloud Droplet Probe LWC.
        - 'CONCD_' (float): Cloud Droplet Probe Number Concentration.

    Returns:
    --------
    dict
        A dictionary containing:
        - 'DataFrame': A modified DataFrame with assigned flight types, cloud status, and location.
        - 'flight_blocks': A dictionary of flight blocks categorized by flight type and cloud status.
        - 'cloud_blocks': A dataframe of flight blocks including blocks of aircraft data inside cloud layers.

    Notes:
    ------
    - The `assign_flight_type` function is responsible for determining whether the flight segments are 'level' or 'profile'.
    - The `block_flight` function segments the flight data based on the assigned flight types and cloud status into specific blocks (e.g., 'Level BL', 'In-Cloud Profiles', etc.).
    - The function ensures proper labeling of different flight segments for further analysis, including cloud status and location (e.g., boundary layer or free airspace).
    """
    # Function to assign flight type "Level" and "Profile" when in/out of cloud
    dict_flight_type = assign_flight_type(df)

    # Extract dataframe that has been modified from the assign_flight_type function
    df_mod = dict_flight_type['DataFrame']
    # Plot time series of aircraft defined flight blocks
    # plot_block_ts(dict_flight_type,i)

    # Run block flight function to return list of Dataframes of "blocked" flight data
    flight_blocks = block_flight(df_mod)
    # Function to assign cloud type from the HCR data
    # flight_block_comp = assign_cloud_type_HCR(flight_blocks,dir,i)
    
    # Plot time series of HCR defined cloud types  
    # plot_hcr_cloud_type(df_mod,flight_block_comp,i)
    return flight_blocks

def select_ERA5_4flight(df,campaign):
    # Define function to filter ERA5 files based on time
    def get_matching_files(pattern, start_dt, end_dt):
        file_list = glob.glob(pattern)
        selected = []
        for file in file_list:
            time_strs = file.split('.')[-2].split('_')
            file_start = datetime.datetime.strptime(time_strs[0], "%Y%m%d%H")
            file_end = datetime.datetime.strptime(time_strs[1], "%Y%m%d%H")
            if file_start <= end_dt and file_end >= start_dt:
                selected.append(file)
        return selected
    
    filepath_sfc = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"
    filepath_pl = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/"
    # Extract the times of the research flight
    month, year = df.Time[0].month, df.Time[0].year
    day_start,day_end = df.Time[0].day, df.Time.iloc[-1].day
    start_hour, end_hour = df.Time[0].hour, df.Time.iloc[-1].hour
    
    # Select the latitude/longitude box to reduce size of era5 data
    lat_max, lat_min = np.floor(df.GGLAT.min()), np.ceil(df.GGLAT.max())
    if campaign == 'SOCRATES':
        lon_adj = 180
    elif campaign  == 'CSET':
        lon_adj = 360
    lon_min, lon_max = np.floor(df.GGLON.min())+lon_adj, np.ceil(df.GGLON.max())+lon_adj
    
    # Flight start and end times
    start_dt = datetime.datetime(year, month, day_start, start_hour) 
    end_dt = datetime.datetime(year, month, day_end, end_hour)+datetime.timedelta(hours=1)
    
    # Make the yearmonth string for file selection
    dir_date = f"{year}{month:02d}"
    
    # Load SST data
    filepath="/glade/campaign/collections/rda/data/ds633.0/e5.oper.an.sfc/"
    # # Construct the search pattern (SST)
    ds_sst = xr.open_mfdataset(get_matching_files(f"{filepath_sfc}{dir_date}/*_sstk.*.nc", start_dt, end_dt), combine='by_coords')
    ds_sst = ds_sst.sel(
    latitude=slice(lat_min, lat_max),  # Select latitudes within range
    longitude=slice(lon_min, lon_max),  # Select longitudes within range
    time=slice(start_dt, end_dt),
    )
    
    # Load 2m tempeature
    ds_t2m = xr.open_mfdataset(get_matching_files(f"{filepath_sfc}{dir_date}/*_2t.*.nc", start_dt, end_dt), combine='by_coords')
    ds_t2m = ds_t2m.sel(
    latitude=slice(lat_min, lat_max),  # Select latitudes within range
    longitude=slice(lon_min, lon_max),  # Select longitudes within range
    time=slice(start_dt, end_dt),
    )
    
    # Load 10m wind speed (u and v components)
    ds_u10 = xr.open_mfdataset(get_matching_files(f"{filepath_sfc}{dir_date}/*_10u.*.nc", start_dt, end_dt), combine='by_coords')[['VAR_10U']]
    ds_u10 = ds_u10.sel(
    latitude=slice(lat_min, lat_max),  # Select latitudes within range
    longitude=slice(lon_min, lon_max),  # Select longitudes within range
    time=slice(start_dt, end_dt),
    )
    ds_v10 = xr.open_mfdataset(get_matching_files(f"{filepath_sfc}{dir_date}/*_10v.*.nc", start_dt, end_dt), combine='by_coords')[['VAR_10V']]
    ds_v10 = ds_v10.sel(
    latitude=slice(lat_min, lat_max),  # Select latitudes within range
    longitude=slice(lon_min, lon_max),  # Select longitudes within range
    time=slice(start_dt, end_dt),
    )
    
    ws = np.sqrt(ds_u10.VAR_10U**2 + ds_v10.VAR_10V**2)
    wind_dir = (270 - np.degrees(np.arctan2(ds_v10.VAR_10V, ds_u10.VAR_10U))) % 360
    
    # Load w data for all pressure levels
    # # Construct the search pattern (SST)
    ds_w = xr.open_mfdataset(get_matching_files(f"{filepath_pl}{dir_date}/*_w.*.nc", start_dt, end_dt), combine='nested', concat_dim='time')
    w_700 = ds_w['W'].sel(level=700).squeeze()
    
    # Select w at 500 hPa
    w_700 = w_700.sortby('time').sel(
        latitude=slice(lat_min, lat_max),
        longitude=slice(lon_min, lon_max),
        time=slice(start_dt, end_dt)
    )
    
    ds_rh700 = xr.open_mfdataset(get_matching_files(f"{filepath_pl}{dir_date}/*_r.*.nc", start_dt, end_dt), combine='nested', concat_dim='time')[['R']].sel(level=700).drop_vars('level', errors='ignore')
    rh = ds_rh700['R'].rename("RH") if 'R' in ds_rh700 else None
    rh = rh.sortby('time').sel(
        latitude=slice(lat_min, lat_max),
        longitude=slice(lon_min, lon_max),
        time=slice(start_dt, end_dt)
    )
    
    # Load pressure level wind speed (u and v components)
    ds_u700 = xr.open_mfdataset(get_matching_files(f"{filepath_pl}{dir_date}/*_u.*.nc", start_dt, end_dt), combine='nested', concat_dim='time')[['U']].sel(level=700).drop_vars('level', errors='ignore')
    ds_u700 = ds_u700.sortby('time').sel(
    latitude=slice(lat_min, lat_max),  # Select latitudes within range
    longitude=slice(lon_min, lon_max),  # Select longitudes within range
    time=slice(start_dt, end_dt)
    )
    ds_v700 = xr.open_mfdataset(get_matching_files(f"{filepath_pl}{dir_date}/*_v.*.nc", start_dt, end_dt), combine='nested', concat_dim='time')[['V']].sel(level=700).drop_vars('level', errors='ignore')
    ds_v700 = ds_v700.sortby('time').sel(
        latitude=slice(lat_min, lat_max),
        longitude=slice(lon_min, lon_max),
        time=slice(start_dt, end_dt)
    )
    
    # Calcualte wind shear (SFC - 700mb)
    ws700 = np.sqrt(ds_u700.U**2 + ds_v700.V**2)  
    wind_shear = ws700-ws
    
    # Load upper temperature data with select pressure levels
    ds_t = xr.open_mfdataset(get_matching_files(f"{filepath_pl}{dir_date}/*_t.*.nc", start_dt, end_dt), combine='nested', concat_dim='time')[['T']].sel(level=800).drop_vars('level', errors='ignore')
    ds_t = ds_t.sortby('time').sel(
        latitude=slice(lat_min, lat_max),
        longitude=slice(lon_min, lon_max),
        time=slice(start_dt, end_dt)
    )
    # Load upper temperature data with select pressure levels
    ds_t700 = xr.open_mfdataset(get_matching_files(f"{filepath_pl}{dir_date}/*_t.*.nc", start_dt, end_dt), combine='nested', concat_dim='time')[['T']].sel(level=700).drop_vars('level', errors='ignore')
    ds_t700 = ds_t700.sortby('time').sel(
        latitude=slice(lat_min, lat_max),
        longitude=slice(lon_min, lon_max),
        time=slice(start_dt, end_dt)
    )
    
    # Calculate M-value
    Rd = 287
    Cp = 1005     
    theta_sfc = ds_sst.SSTK*(1)**(Rd/Cp)
    theta_800 = ds_t*(1013.25/800)**(Rd/Cp)
    
    M = theta_sfc.T - theta_800.T
    M = M.transpose("time", "latitude", "longitude")
    
    dt = ds_t2m.VAR_2T - ds_sst.SSTK
    # Constants
    Re = 6.371e6  # Earth radius in meters
    deg2rad = np.pi / 180
    phi = np.deg2rad(ds_sst.SSTK['latitude'])
    # meters per 1° at this latitude
    m_per_deg_lon = Re * np.cos(phi) * deg2rad
    m_per_deg_lat = Re * deg2rad
    # gradients in K/m  (NOTE the division by meters-per-degree)
    dT_dx = ds_sst.SSTK.differentiate("longitude") / m_per_deg_lon   # K/m
    dT_dy = ds_sst.SSTK.differentiate('latitude') / m_per_deg_lat   # K/m
    # advection: K/s -> K/day
    Tadv = -(ds_u10['VAR_10U'] * dT_dx + ds_v10['VAR_10V'] * dT_dy) * 86400.0
    # Convert to K/day
    Tadv = Tadv.rename("Tadv")

    # Calculate EIS following Wood and Bretherton (2006, J. Climate)

    cp = 1004.     # specific heat at constant pressure for dry air (J / kg / K)
    Rd = 287.         # gas constant for dry air (J / kg / K)
    kappa = Rd / cp
    Lhvap = 2.5e6    # Latent heat of vaporization (J / kg)
    g = 9.81 # m/s^2
    cp = 1004 # J/K/kg
    Lv = 2.5e6 # J/kg
    
    Rv = 461 # J/K/kg;
    Ra = 287 # J/K/kg
    
    def get_qsat(T,p):
        Tcel = T-273.15
        es=6.11*10**(7.5*Tcel/(Tcel+273.15))
        return 0.622*es/p
    
    # Calculate lower tropospheric stability (LTS)
    theta_700 = ds_t700.T*(1013.25/700)**kappa
    LTS = theta_700 - theta_sfc
    
    T850 = (ds_t2m.VAR_2T+ds_t700.T)/2
    
    Gammam = (g/cp*(1.0 - (1.0 + Lhvap*get_qsat(T850,850) / Rd / T850) /
                 (1.0 + Lhvap**2 * get_qsat(T850,850)/ cp/Rv/T850**2)))
    
    # Assume exponential decrease of pressure with scale height given by surface temperature
    z700 = (Rd * ds_t2m.VAR_2T / g) * np.log(1000 / 700)
    # Assume 80% relative humidity to compute LCL, appropriate for marine boundary layer
    Tadj = Tadj = ds_t2m.VAR_2T-55.  # in Kelvin
    LCL = cp/g*(Tadj - (1/Tadj - np.log(0.8)/2840.)**(-1))
    
    EIS = LTS - Gammam*(z700 - LCL)
    
    # Merge the dataset variables used later
    ds = {
    'deltaT': dt,
    'Tadv': Tadv,
    'M': M,
    'w_700': w_700,
    'SST': ds_sst.SSTK,
    'WS': ws,
    'Wind_shear': wind_shear,
    'RH700': rh,
    'EIS': EIS
     }

    return ds

def wrap180(lon):
    # Map any longitude to [-180, 180)
    return (lon + 180.0) % 360.0 - 180.0

def nearest_time_indices(era5_times_ns, flight_times_ns):
    # era5_times_ns: 1D int64 nanoseconds, sorted
    # flight_times_ns: 1D int64 nanoseconds
    idx_right = np.searchsorted(era5_times_ns, flight_times_ns, side="left")
    idx_left  = np.clip(idx_right - 1, 0, len(era5_times_ns) - 1)
    idx_right = np.clip(idx_right,       0, len(era5_times_ns) - 1)
    # choose whichever neighbor is closer
    choose_right = np.abs(era5_times_ns[idx_right] - flight_times_ns) < np.abs(era5_times_ns[idx_left] - flight_times_ns)
    return np.where(choose_right, idx_right, idx_left)

def collocate_ERA5_dat(ds, blocks):
    """
    Vectorized collocation of ERA5 fields onto flight blocks.
    Expects ds variables with dims ('time','latitude','longitude').
    Returns updated `blocks` in-place with columns:
    ERA5_SST, M, w700, deltaT, Wind_sp, Wind_shear, Tadv, RH700, EIS.
    """
    # --- Pull coords as numpy (keep ordering exactly as in ds)
    ds = xr.Dataset(ds)  # now ds has proper coords
    lat_vals = ds['latitude'].values
    lon_vals_ds = ds['longitude'].values
    lon_vals_wrapped = wrap180(lon_vals_ds)  # build KDTree in [-180,180)

    # --- KDTree on (lat, lon_wrapped)
    lon_grid, lat_grid = np.meshgrid(lon_vals_wrapped, lat_vals)
    tree = cKDTree(np.c_[lat_grid.ravel(), lon_grid.ravel()])

    # --- Era5 time (sorted)
    t_era = ds['time'].values.astype('datetime64[ns]')
    t_era_ns = t_era.view('int64')

    # --- Preload arrays once (-> NumPy) for super-fast indexing
    def arr3(name):
        return ds[name].transpose('time','latitude','longitude').compute().values  # (T,Y,X) numpy
    arr = {
        'ERA5_SST': arr3('SST'),
        'M':        arr3('M'),
        'w700':     arr3('w_700'),
        'deltaT':   arr3('deltaT'),
        'Wind_sp':  arr3('WS'),
        'Wind_shear': arr3('Wind_shear'),
        'Tadv':     arr3('Tadv'),
        'RH700':    arr3('RH700'),
        'EIS':      arr3('EIS'),
    }

    ny, nx = len(lat_vals), len(lon_vals_ds)

    # --- Iterate blocks, but do all lookups in one vectorized shot per block
    for val in blocks:
        block_list = blocks[val]
        for i, block in enumerate(block_list):
            block = block.copy()
            block = block.dropna(subset=['GGLAT','GGLON'])
            if len(block) == 0:
                block_list[i] = block
                continue

            # Flight coords/time
            flt_lat = block['GGLAT'].values.astype(float)
            flt_lon = wrap180(block['GGLON'].values.astype(float))  # match KDTree frame
            flt_t   = block['Time'].astype('datetime64[ns]').values
            flt_t_ns = flt_t.view('int64')

            # Nearest time indices (vectorized, no giant argmin)
            ti = nearest_time_indices(t_era_ns, flt_t_ns)  # (N,)

            # Nearest gridpoint via KDTree (vectorized)
            _, flat_idx = tree.query(np.c_[flt_lat, flt_lon])     # (N,)
            yi, xi = np.unravel_index(flat_idx, (ny, nx))         # (N,), (N,)

            # Gather all variables in one pass
            for out_name, A in arr.items():
                block[out_name] = A[ti, yi, xi]

            # Put back
            block_list[i] = block
        blocks[val] = block_list

    return blocks

def cloud_regime_old(fblks):
    # Iterate through blocks
    for val in fblks:
        block_type = fblks[val]
        for i in range(len(block_type)):
            block = block_type[i].copy()
            block['cloud_regime'] = pd.Series('Undetermined', index=block.index, dtype='object')
            condition_cum = ((block['M'] > -7) & (block['Wind_shear'] < 6)) | ((block['M'] > -7) & (block['Wind_sp'] > 10))
            condition_strcu = ((block['M'] <= -10) & (block['Wind_sp'] < 10)) | ((block['M'] <= -10) | (block['Wind_shear'] > 6))
            # Assign values based on conditions
            block.loc[condition_cum, 'cloud_regime'] = 'Open-Cell Cu'
            block.loc[condition_strcu, 'cloud_regime'] = 'Stratiform'
    
            block_type[i] = block
    
        fblks[val] = block_type

    return fblks

def cloud_regime(fblks, campaign):
    """
    Assign cloud_regime per block:
      - For CSET:
          Stratocumulus if any of:
              1) M < -10 & RH700 > 30
              2) M < -10 & ERA5_SST < 295
              3) M < -11 & Tadv < 0
          Open-cell Cumulus if any of:
              1) M >= -10 & RH700 <= 30
              2) M >= -10 & ERA5_SST >= 296
              3) M >= -10 & Tadv >= 0.0005   # assume same units as your Tadv column (ideally K/day)
        - For SOCRATES: (keeps your existing rules)
    """
    import numpy as np
    import pandas as pd

    for val in fblks:
        block_type = fblks[val]
        for i in range(len(block_type)):
            block = block_type[i].copy()

            # default class
            block['cloud_regime'] = pd.Series('Unknown', index=block.index, dtype='object')

            if campaign == 'SOCRATES':
                condition_cum = ((block['M'] > -7) & (block['Wind_shear'] < 6)) | \
                                ((block['M'] > -7) & (block['Wind_sp'] > 10))
                condition_strcu = ((block['M'] <= -10) & (block['Wind_sp'] < 10)) | \
                                  ((block['M'] <= -10) | (block['Wind_shear'] > 6))

                block.loc[condition_cum,  'cloud_regime'] = 'Open-Cell'
                block.loc[condition_strcu,'cloud_regime'] = 'Stratocumulus'  # or 'Closed-Cell' if you prefer

            elif campaign == 'CSET':
                # Safely get needed variables (treat missing as NaN so conditions become False)
                M    = block['M']
                RH   = block.get('RH700',     pd.Series(np.nan, index=block.index))
                SST  = block.get('ERA5_SST',  block.get('sst', pd.Series(np.nan, index=block.index)))
                Tadv = block.get('Tadv',      pd.Series(np.nan, index=block.index))

                # --- Stratocumulus (any of the three) ---
                cond_strat = (
                    ((M < -10) & (RH > 30)) |
                    ((M < -10) & (SST < 295)) |
                    ((M < -10) & (Tadv < 0))
                )

                # --- Open-cell Cumulus (any of the three) ---
                cond_opencu = (
                    ((M >= -10) & (RH <= 30)) |
                    ((M >= -10) & (SST >= 296)) |
                    ((M >= -10) & (Tadv >= 0.0005))  # assumes Tadv units match your data (recommended: K/day)
                )

                # Apply labels
                block.loc[cond_strat,  'cloud_regime'] = 'Stratocumulus'
                block.loc[cond_opencu, 'cloud_regime'] = 'Open-Cell'

            # write back
            block_type[i] = block
        fblks[val] = block_type

    return fblks
       
def write_RF_nc(fblks_cr, rf, campaign, save_path="./"):
    if not os.path.exists(save_path):
        print(f"Creating directory: {save_path}")
        os.makedirs(save_path)
    name = f"{save_path}/{campaign}_{rf}.nc"
    if os.path.exists(name):
        print(f"{name} exists, skipping write.")
        return
    combined = []
    if isinstance(fblks_cr, dict):
        for label, df_list in fblks_cr.items():
            for i, df in enumerate(df_list):
                df = df.copy()
                df["flight"] = rf
                df["block_label"] = label
                df["block_index"] = i
                combined.append(df)

        df_all = pd.concat(combined, ignore_index=True)
        df_all = df_all.set_index(["block_label", "block_index", "Time"])
        ds = df_all.reset_index().to_xarray()

        ds.to_netcdf(name)
        print(f"Wrote {name}")

def plot_block_ts(dict,idx):

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.dates as mdates
    import numpy as np
    import pandas as pd

    # Assuming the four DataFrames are already created
    # profile_in_cloud, level_in_cloud, level_out_cloud_bl, level_out_cloud_fr
    # Dict = assign_flight_type(df)
    df = dict['DataFrame']
    # print(df)
    blocks = dict['flight_blocks']
    incloud = dict['Cloud_blocks']
    # Creating the four DataFrames based on flight_type, cloud status, and Location
    profile_in_cloud = df[(df['flight_type'] == 'profile') & (df['cloud_status'] == 'In-cloud')]
    level_in_cloud = df[(df['flight_type'] == 'level') & (df['cloud_status'] == 'In-cloud') & (df['Location'] != 'BL')]
    level_out_cloud_bl = df[(df['flight_type'] == 'level') & (df['cloud_status'] != 'In-cloud') & (df['Location'] == 'BL')]
    level_out_cloud_fr = df[(df['flight_type'] == 'level') & (df['cloud_status'] != 'In-cloud') & (df['Location'] == 'Free')]
    # Create figure and GridSpec layout
    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.7, 1])  # First subplot is 3x the width of the second

    # First subplot (larger)
    ax1 = fig.add_subplot(gs[0])  # Assigning first subplot

    # Plot flight altitude for each DataFrame in different colors
    ax1.plot(df['Time'], df['GGALT'], color='k', linewidth=3)

    # Plot In-cloud, Level In-cloud, and other data
    ax1.scatter(profile_in_cloud['Time'], profile_in_cloud['GGALT'], color='red', label='Profile In-cloud', marker='s', s=1, zorder=2)
    ax1.scatter(level_in_cloud['Time'], level_in_cloud['GGALT'], color='blue', label='Level In-cloud', marker='s', s=1, zorder=2)
    ax1.scatter(level_out_cloud_bl['Time'], level_out_cloud_bl['GGALT'], color='green', label='Level Out-of-cloud (BL)', marker='s', s=1, zorder=2)
    ax1.scatter(level_out_cloud_fr['Time'], level_out_cloud_fr['GGALT'], color='purple', label='Level Out-of-cloud (Free)', marker='s', s=1, zorder=2)

    # Set x-axis limits
    start_limit = df['Time'].min()
    end_limit = df['Time'].max()
    ax1.set_xlim([start_limit, end_limit])

    # Format x-axis to show only hours, minutes, and seconds
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))  # Adjust format as needed
    fig.autofmt_xdate()

    # Shade blocks according to their type (level or profile)
    for _, row in blocks.iterrows():
        start_time, end_time, flight_type = row['start_time'], row['end_time'], row['flight_type']
        if flight_type == 'level':
            ax1.axvspan(start_time, end_time, color='blue', alpha=0.3, label="Level leg" if 'Level leg' not in ax1.get_legend_handles_labels()[1] else None)
        elif flight_type == 'profile':
            ax1.axvspan(start_time, end_time, color='goldenrod', alpha=0.5, label="Profiling" if 'Profiling' not in ax1.get_legend_handles_labels()[1] else None)

    # Labels and title
    ax1.set_xlabel('Time UTC (MM-dd HH:mm)')
    ax1.set_ylabel('Altitude (m)')
    ax1.set_title('Altitude Time Series separating into "Level legs" and "Profiles"')
    ax1.set_ylim(-200, np.max(df['GGALT'] + 600))
    ax1.legend(loc='upper right', ncol=6, markerscale=6,fontsize=8)
    ax1.grid(True)

    # Second subplot (smaller)
    ax2 = fig.add_subplot(gs[1])  # Assigning second subplot

    # Find best match for column names dynamically
    plwc_col = next((col for col in df.columns if 'PLWCD' in col), None) or \
           next((col for col in df.columns if 'PLWC' in col), None)
    concd_col = next((col for col in df.columns if 'CONCD' in col), None)

    # Scatter plot for concentration and altitude
    ax2.scatter(df[concd_col], df.GGALT, color='b', label='CDP', alpha=0.5, marker='^', s=4)

    # Log scale for x-axis
    ax2.set_xscale('log')
    ax2.set_xlabel('Conc (#/cm3)')
    ax2.set_ylabel('Altitude (meters)')
    ax2.grid(True)
    ax2.set_xlim(.01, 1000)

    # Add second x-axis on top
    ax2_top = ax2.twiny()
    ax2_top.plot(df[plwc_col], df.GGALT, color='orange', alpha=.7, linestyle='--', label='CDP LWC')
    ax2_top.set_xlabel('g/m3')
    ax2_top.set_xscale('log')
    ax2_top.set_xlim(0.0001, 10)
    ax2.set_ylim(-200, np.max(df['GGALT'] + 600))

    # Shade blocked altitude regions in ax2
    for i, row in incloud.iterrows():
        ax2.fill_betweenx(
            y=[row['lower_bound'], row['upper_bound']],  # Altitude range for shading
            x1=0.01,  # Left bound (min x-value)
            x2=1000,  # Right bound (max x-value)
            color='red', alpha=0.3, label="Cloud layer" if 'Cloud layer' not in ax2.get_legend_handles_labels()[1] else None
        )

    # Merge legends from both axes
    lines, labels = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_top.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper right')  # Merge legends from both axes

    # Reorder the legends if needed
    all_lines = lines + lines2
    all_labels = labels + labels2
    if len(all_labels) >= 3:
        new_order = [0, 2, 1]  # Modify based on the desired order
        all_lines = [all_lines[i] for i in new_order]
        all_labels = [all_labels[i] for i in new_order]

    # Apply reordered legend
    ax2.legend(all_lines, all_labels, loc='upper right', markerscale=3)

    # Adjust layout to prevent overlap
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.12)  # Increase spacing between subplots

    # Save the figure
    rf_id = f"RF_{idx+1:02d}"
    filename = f"CSET_Altitude_Flight_Cloud_type{rf_id}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')  # Save as PNG with high resolution

   
def plot_hcr_cloud_type(df,Flight_blocks,idx):
    
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import ListedColormap, BoundaryNorm
    import matplotlib.dates as mdates
    
    # Initialize figure and gridspec for plotting
    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(2, 1, height_ratios=[10, 1])
    
    ax1 = fig.add_subplot(gs[0, 0])  # Main altitude plot
    ax2 = fig.add_subplot(gs[1, 0])  # Echo Type plot
    
    tick_values = [14, 16, 18, 25, 30, 32, 34, 36, 38]
    
    # Update to use the new colormap interface
    spectral_cmap = plt.colormaps['Set1']  # Access colormap directly
    tick_to_color = {tick: spectral_cmap(i / len(tick_values)) for i, tick in enumerate(tick_values)}  # Map each tick to a color
       
    # Plot flight altitude for each DataFrame in different colors
    ax1.plot(df['Time'], df['GGALT'], color='k', linewidth=3)
    
    start_limit = df['Time'].min()
    end_limit = df['Time'].max()
    ax1.set_xlim([start_limit, end_limit])
    
    for val in Flight_blocks:
        block_type = Flight_blocks[val]
        for i in range(len(block_type)):
            block = block_type[i]
            start_time = block['Time'].iloc[0]
            end_time = block['Time'].iloc[-1]
            mean_echo_type = np.nanmean(block['Echo_Type'])
        
            # Find the closest tick value and corresponding color
            closest_tick = tick_values[np.argmin(np.ceil(np.abs(np.array(tick_values) - mean_echo_type)))]
            color = tick_to_color[closest_tick]
        
            # Shade region on ax1
            ax1.axvspan(start_time, end_time, color=color, alpha=0.8)
    
            # Plot scatter on ax2
            ax2.scatter(
                block.Time,
                np.zeros(len(block.Time)),
                c=block.Echo_Type,
                cmap='Set1',
                marker='s',
                s=4,
                vmin=min(tick_values),
                vmax=max(tick_values)
            )
    
    # start_time = pd.to_datetime("2018-01-16 01:30:00")
    # end_time = pd.to_datetime("2018-01-16 03:00:00")
    # ax1.set_xlim(start_time, end_time)
    ax1.set_ylabel('Altitude (m)')
    # Clean up ax2 to make it look like a color strip
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_yticks([])
    ax2.set_xlabel('Time UTC (MM-dd HH:mm)')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    fig.autofmt_xdate()
    # Set x-axis limits
    
    ax2.set_xlim([start_limit, end_limit])
    
    # Define the tick values (bin labels)
    tick_values = [14, 16, 18, 25, 30, 32, 34, 36, 38]
    tick_labels = [
        "stratiform low",
        "stratiform mid",
        "stratiform high",
        "mixed",
        "convective",
        "conv. elevated",
        "conv. shallow",
        "conv. mid",
        "conv. deep"
    ]
    # We need to define edges for each bin; to get N blocks, we need N+1 boundaries
    bounds = list(range(len(tick_values) + 1))  # e.g., 0, 1, 2, ..., 9
    
    # Create a colormap with N colors
    colors = plt.cm.Set1(np.linspace(0, 1, len(tick_values)))
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)
    
    # Add vertical colorbar on the right
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=[ax1, ax2],
        orientation='vertical',
        ticks=np.arange(len(tick_values)) + 0.5,  # Tick in the center of each block
        pad=0.012
    )
    # Set category labels instead of numbers
    cbar.ax.set_yticklabels(tick_labels)
    # Set colorbar labels and title
    cbar.ax.set_yticklabels(tick_labels)
    cbar.set_label("HCR Echo Type")
    # Turn ticks inside for all three axes
    for ax in [ax1, ax2]:
        ax.tick_params(direction='in', which='both', top=True, right=True)
    # # Show the plot
    # plt.show()

    # Save the figure
    rf_id = f"RF_{idx+1:02d}"
    filename = f"SOCRATES_Altitude_Flight_HCR_Cloud_Echo_{rf_id}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')  # Save as PNG with high resolution