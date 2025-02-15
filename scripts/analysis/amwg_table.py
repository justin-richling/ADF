import numpy as np
import xarray as xr
import sys
from pathlib import Path
import warnings  # use to warn user about missing files.
import shutil


#Import "special" modules:
try:
    import scipy.stats as stats # for easy linear regression and testing
except ImportError:
    print("Scipy module does not exist in python path, but is needed for amwg_table.")
    print("Please install module, e.g. 'pip install scipy'.")
    sys.exit(1)
#End except

try:
    import pandas as pd
except ImportError:
    print("Pandas module does not exist in python path, but is needed for amwg_table.")
    print("Please install module, e.g. 'pip install pandas'.")
    sys.exit(1)
#End except

#Import ADF-specific modules:
import plotting_functions as pf

def amwg_table(adf):

    """
    Main function goes through series of steps:
    - load the variable data
    - Determine whether there are spatial dims; if yes, do global average (TODO: regional option)
    - Apply annual average (TODO: add seasonal here)
    - calculates the statistics
      + mean
      + sample size
      + standard deviation
      + standard error of the mean
      + 5/95% confidence interval of the mean
      + linear trend
      + p-value of linear trend
    - puts statistics into a CSV file
    - generates simple HTML that can display the data

    Description of needed inputs from ADF:

    case_names      -> Name(s) of CAM case provided by "cam_case_name"
    input_ts_locs   -> Location(s) of CAM time series files provided by "cam_ts_loc"
    output_loc      -> Location to write AMWG table files to, provided by "cam_diag_plot_loc"
    var_list        -> List of CAM output variables provided by "diag_var_list"
    var_defaults    -> Dict that has keys that are variable names and values that are plotting preferences/defaults.

    and if doing a CAM baseline comparison:

    baseline_name     -> Name of CAM baseline case provided by "cam_case_name"
    input_ts_baseline -> Location of CAM baseline time series files provied by "cam_ts_loc"

    """

    #Import necessary modules:
    from adf_base import AdfError
    #from adf_dataset import SuppressWarningsPrint
    

    #Additional information:
    #----------------------

    # GOAL: replace the "Tables" set in AMWG
    #       Set Description
    #   1 Tables of ANN, DJF, JJA, global and regional means and RMSE.
    #
    # STRATEGY:
    # I think the right solution is to generate one CSV (or other?) file that
    # contains all of the data.
    # So we need:
    # - a function that would produces the data, and
    # - then call a function that adds the data to a file
    # - another function(module?) that uses the file to produce a "web page"

    # IMPLEMENTATION:
    # - assume that we will have time series of global averages already ... that should be done ahead of time
    # - given a variable or file for a variable (equivalent), we will calculate the all-time, DJF, JJA, MAM, SON
    #   + mean
    #   + standard error of the mean
    #     -- 95% confidence interval for the mean, estimated by:
    #     ---- CI95 = mean + (SE * 1.96)
    #     ---- CI05 = mean - (SE * 1.96)
    #   + standard deviation
    # AMWG also includes the RMSE b/c it is comparing two things, but I will put that off for now.

    # DETAIL: we use python's type hinting as much as possible

    # in future, provide option to do multiple domains
    # They use 4 pre-defined domains:
    domains = {"global": (0, 360, -90, 90),
               "tropics": (0, 360, -20, 20),
               "southern": (0, 360, -90, -20),
               "northern": (0, 360, 20, 90)}

    # and then in time it is DJF JJA ANN

    # within each domain and season
    # the result is just a table of
    # VARIABLE-NAME, RUN VALUE, OBS VALUE, RUN-OBS, RMSE
    #----------------------

    #Notify user that script has started:
    msg = "\n  Calculating AMWG variable tables..."
    print(f"{msg}\n  {'-' * (len(msg)-3)}")

    #with adf.data.SuppressWarningsPrint(suppress=True):  # Suppress warnings inside this block
    with adf.data.SuppressWarningsPrint(suppress=adf.verbose):  # Suppress warnings inside this block

        #Extract needed quantities from ADF object:
        #-----------------------------------------
        var_list     = adf.diag_var_list
        var_defaults = adf.variable_defaults

        #Check if ocean or land fraction exist
        #in the variable list:
        for var in ["OCNFRAC", "LANDFRAC"]:
            if var in var_list:
                #If so, then move them to the front of variable list so
                #that they can be used to mask or vertically interpolate
                #other model variables if need be:
                var_idx = var_list.index(var)
                var_list.pop(var_idx)
                var_list.insert(0,var)
            #End if
        #End if

        #Special ADF variable which contains the output paths for
        #all generated plots and tables for each case:
        output_locs = adf.plot_location

        #CAM simulation variables (these quantities are always lists):
        #case_names    = adf.get_cam_info("cam_case_name", required=True)
        #CAM simulation variables (these quantities are always lists):
        test_case_names = adf.get_cam_info("cam_case_name", required=True)
        #input_ts_locs = adf.get_cam_info("cam_ts_loc", required=True)
        input_locs = adf.ts_locs["test"]

        #adf.get_baseline_info("cam_climo_loc")
        #input_climo_locs = adf.get_cam_info("cam_climo_loc")
        input_climo_locs = adf.climo_locs["test"]

        #Grab case years
        syear_cases = adf.climo_yrs["syears"]
        eyear_cases = adf.climo_yrs["eyears"]

        test_nicknames = adf.case_nicknames["test_nicknames"]
        base_nickname = adf.case_nicknames["base_nickname"]
        nicknames = test_nicknames + [base_nickname]
        





        """
        #Check if a baseline simulation is also being used:
        if not adf.get_basic_info("compare_obs"):
            #Extract CAM baseline variaables:
            baseline_name     = adf.get_baseline_info("cam_case_name", required=True)
            input_ts_baseline = adf.get_baseline_info("cam_ts_loc", required=True)

            case_names.append(baseline_name)
            input_ts_locs.append(input_ts_baseline)

            #Save the baseline to the test case's plots directory:
            if len(test_case_names) == 1:
                output_locs.append(output_locs[0])
        else:
            print("AMWG table doesn't currently work with obs, so obs table won't be created.")
        #End if
        """






        #Check if user wants to skip time series file creation
        '''if not input_locs:
            #print("User indicates no time series files will be used")
            #print()
            emsg = "\n  User indicates no time series files will be used."
            emsg += " Looking if table already exisits:"
            print(emsg)

            #if ah:
            for case_idx, case_name in enumerate(case_names):
                #Convert output location string to a Path object:
                output_location = Path(output_locs[case_idx])
                #Create output file name:
                output_csv_file = output_location / f"amwg_table_{case_name}.csv"
                if Path(output_csv_file).is_file():
                    print(f"\t - AMWG table for '{case_name}' exists, adding to website.")
                    table_df = pd.read_csv(output_csv_file)
                    # last step is to add table dataframe to website (if enabled):
                    adf.add_website_data(table_df, case_name, case_name, plot_type="Tables")
                else:
                    print(f"\t - AMWG table for '{case_name}' does not exist.")
                    print('\t  check here:',output_csv_file,"\n")
            #input_locs = []
            pass#return
        else:'''
        #if 1==1:
        #    input_locs = [None]*len(case_names)
        #End if
        print("\nTest input_locs",input_locs,"\n")
        csv_locs = {}
        #Check if a baseline simulation is also being used:
        if not adf.get_basic_info("compare_obs"):
            #Extract CAM baseline variaables:
            baseline_name     = adf.get_baseline_info("cam_case_name", required=True)
            #input_loc = adf.get_baseline_info("cam_ts_loc", required=True)
            input_loc = adf.ts_locs["baseline"]
            print("\nBaseline input_locs",input_loc,"\n")
            #input_climo_loc = adf.get_baseline_info("cam_climo_loc")
            input_climo_loc = adf.climo_locs["baseline"]
            input_climo_locs.append(input_climo_loc)

            #Grab baseline years (which may be empty strings if using Obs):
            syear_baseline = adf.climo_yrs["syear_baseline"]
            eyear_baseline = adf.climo_yrs["eyear_baseline"]

            syear_cases.append(syear_baseline)
            eyear_cases.append(eyear_baseline)

            #Convert output location string to a Path object:
            output_location = Path(output_locs[0])
            if not input_loc:
                #print("User indicates no time series files will be used")
                #print()
                emsg = "\n  User indicates no time series files will be used."
                emsg += " Looking if table already exisits:"
                print(emsg)

                #if ah:
                #for case_idx, case_name in enumerate(case_names):
                #Create output file name:
                output_csv_file = output_location / f"amwg_table_{baseline_name}.csv"
                csv_locs[baseline_name] = output_csv_file
                if Path(output_csv_file).is_file():
                    print(f"\t - AMWG table for '{baseline_name}' exists, adding to website.")
                    table_df = pd.read_csv(output_csv_file)
                    # last step is to add table dataframe to website (if enabled):
                    adf.add_website_data(table_df, baseline_name, baseline_name, plot_type="Tables")
                else:
                    print(f"\t - AMWG table for '{baseline_name}' does not exist.")
                    print('\t  check here:',output_csv_file,"\n")
                input_locs.append(None)
                pass#return
            else:
                #input_loc = adf.get_baseline_info("cam_climo_loc")
                input_locs.append(input_loc)

            #case_names.append(baseline_name)
            #if input_loc:
            #case_names.append(baseline_name)
            case_names = test_case_names + [baseline_name]
                #input_locs.append(input_loc)

            #Save the baseline to the first case's plots directory:
            output_locs.append(output_location)
        else:
            print("AMWG table doesn't currently work with obs, so obs table won't be created.")
        #End if

        #-----------------------------------------
        print("input_locs",input_locs,"\n")
        #Loop over CAM cases:
        #Initialize list of case name csv files for case comparison check later
        csv_list = []
        for case_idx, case_name in enumerate(case_names):
            syear = syear_cases[case_idx]
            eyear = eyear_cases[case_idx]

            #Convert output location string to a Path object:
            output_location = Path(output_locs[case_idx])

            """#Generate input file path:
            input_location = input_locs[case_idx]
            print("\n\tTS input_location",input_location)

            if not input_location:
                print(f"\t ** User supplied climo files for {case_name}, will make only global mean (no other stats) for each variable. Thanks and have a nice day.")
                is_climo = True
            else:
                is_climo = False

            #print("\n\tis_climo:",is_climo,"\n")

            #Generate input file path:
            if not is_climo:
                input_location = Path(input_locs[case_idx])
            if is_climo:
                input_location = Path(input_climo_locs[case_idx])
            print("\tinput_location",input_location)

            #Check that time series input directory actually exists:
            if not input_location.is_dir():
                errmsg = f"Directory '{input_location}' not found.  Script is exiting."
                raise AdfError(errmsg)
            #Write to debug log if enabled:
            adf.debug_log(f"DEBUG: location of files is {str(input_location)}")

            #Notify user as attempting table creation:
            print(f"\n  Calculating AMWG variable table for '{case_name}'...")"""

            #Create output file name:
            output_csv_file = output_location / f"amwg_table_{case_name}.csv"
            csv_locs[case_name] = output_csv_file

            #Given that this is a final, user-facing analysis, go ahead and re-do it every time:
            if Path(output_csv_file).is_file():
                Path.unlink(output_csv_file)
            #End if

            #Create/reset new variable that potentially stores the re-gridded
            #ocean fraction xarray data-array:
            ocn_frc_da = None

            #Notify user that script has started:
            print(f"\n  Calculating AMWG variable table for '{case_name}'...")
        
            #Loop over CAM output variables:
            for var in var_list:
                is_climo = False # default to time series
                #Generate input file path:
                if input_locs[case_idx]:
                    input_location = Path(input_locs[case_idx])
                    #print("\n\tTS input_location",input_location)

                    filenames = f'{case_name}.*.{var}.*nc'
                    files = sorted(input_location.glob(filenames))
                else:
                    files = None

                # If no files exist, try to move to next variable. --> Means we can not proceed with this variable, and it'll be problematic later.
                if not files:
                    # Try for climo files:
                    msg = f"\t    INFO: Time series files for variable '{var}' not found.  Checking on climo files."
                    print(msg)
                    filenames = f'{case_name}_{var}_climo.nc'
                    try_input_location = Path(input_climo_locs[case_idx])
                    try_files = sorted(try_input_location.glob(filenames))
                    if not try_files:
                        #set_warning_filter(enable=True)  # Suppress warnings
                        errmsg = f"\t    WARNING: Climo files for variable '{var}' not found.  Script will continue to next variable."
                        print(errmsg)
                        continue
                    else:
                        print(f"\t    INFO: User supplied climo files for {var} in {case_name}, will make only global mean (no other stats) for each variable. Thanks and have a nice day.")
                        files = try_files
                        input_location = try_input_location
                        is_climo = True
                #End if

                """if not input_location:
                    print(f"\t ** User supplied climo files for {var} in {case_name}, will make only global mean (no other stats) for each variable. Thanks and have a nice day.")
                    is_climo = True
                else:
                    is_climo = False"""

                #print("\n\tis_climo:",is_climo,"\n")

                """#Generate input file path:
                if not is_climo:
                    input_location = Path(input_locs[case_idx])
                if is_climo:
                    input_location = Path(input_climo_locs[case_idx])
                print("\tinput_location",input_location)"""

                #Check that time series input directory actually exists:
                if not input_location.is_dir():
                    errmsg = f"amwg_table: Time series directory '{input_location}' not found.  Script is exiting."
                    raise AdfError(errmsg)
                #Write to debug log if enabled:
                adf.debug_log(f"DEBUG: location of files is {str(input_location)}")

                #Notify users of variable being added to table:
                print(f"\t - Variable '{var}' being added to table")

                #Create list of time series files present for variable:
                #ts_filenames = f'{case_name}.*.{var}.*nc'
                #ts_files = sorted(input_location.glob(ts_filenames))


                """if is_climo:
                    #Create list of climo files present for variable:
                    filenames = f'{case_name}_{var}_climo.nc'
                else:
                    #Create list of time series files present for variable:
                    filenames = f'{case_name}.*.{var}.*nc'"""
                """files = sorted(input_location.glob(filenames))

                # If no files exist, try to move to next variable. --> Means we can not proceed with this variable, and it'll be problematic later.
                if not files:
                    errmsg = f"\t    ** Time series files for variable '{var}' not found.  Script will continue to next variable."
                    warnings.warn(errmsg)
                    continue
                #End if"""

                """#TEMPORARY:  For now, make sure only one file exists:
                if len(files) != 1:
                    errmsg =  "Currently the AMWG table script can only handle one time series file per variable."
                    errmsg += f" Multiple files were found for the variable '{var}', so it will be skipped."
                    print(errmsg)
                    continue
                #End if"""

                #Load model variable data from file:
                ds = pf.load_dataset(files)

                if not is_climo:
                    #Average time dimension over time bounds, if bounds exist:
                    if 'time_bnds' in ds:
                        time = ds['time']
                        # NOTE: force `load` here b/c if dask & time is cftime, throws a NotImplementedError:
                        time = xr.DataArray(ds['time_bnds'].load().mean(dim='nbnd').values, dims=time.dims, attrs=time.attrs)
                        ds['time'] = time
                        ds.assign_coords(time=time)
                        ds = xr.decode_cf(ds)

                #print("afdasdfs",ds.time,"\n")
                #data = ds[var]
                if len(files) > 1:
                    # Slice for years 0500 to 0521
                    # Slice using only the 4-digit year
                    time_slice = slice(str(syear).zfill(4), str(eyear).zfill(4))
                    ds = ds.sel(time=time_slice)
                    #print("afdasdfs",ds.time,"\n")
                    data = ds[var].compute()
                    #print(data.time)
                    #data = data.sel(time=slice())
                else:
                    data = ds[var]

                #Extract units string, if available:
                if hasattr(data, 'units'):
                    unit_str = data.units
                else:
                    unit_str = '--'

                #Check if variable has a vertical coordinate:
                if 'lev' in data.coords or 'ilev' in data.coords:
                    #set_warning_filter(enable=True)  # Suppress warnings
                    print(f"\t    WARNING: Variable '{var}' has a vertical dimension, "+\
                        "which is currently not supported for the AMWG Table. Skipping...")
                    #Skip this variable and move to the next variable in var_list:
                    continue
                #End if

                #Extract defaults for variable:
                var_default_dict = var_defaults.get(var, {})

                #Check if variable should be masked:
                if 'mask' in var_default_dict:
                    if var_default_dict['mask'].lower() == 'ocean':
                        #Check if the ocean fraction has already been regridded
                        #and saved:
                        if ocn_frc_da is not None:
                            ofrac = ocn_frc_da
                            # set the bounds of regridded ocnfrac to 0 to 1
                            ofrac = xr.where(ofrac>1,1,ofrac)
                            ofrac = xr.where(ofrac<0,0,ofrac)

                            # apply ocean fraction mask to variable
                            data = pf.mask_land_or_ocean(data, ofrac, use_nan=True)
                            #data = var_tmp
                        else:
                            #set_warning_filter(enable=True)  # Suppress warnings
                            print(f"\t    WARNING: OCNFRAC not found, unable to apply mask to '{var}'")
                        #End if
                    else:
                        #Currently only an ocean mask is supported, so print warning here:
                        wmsg = "\t    INFO: Currently the only variable mask option is 'ocean',"
                        wmsg += f"not '{var_default_dict['mask'].lower()}'"
                        print(wmsg)
                    #End if
                #End if

                #If the variable is ocean fraction, then save the dataset for use later:
                if var == 'OCNFRAC':
                    ocn_frc_da = data
                #End if

                # we should check if we need to do area averaging:
                if len(data.dims) > 1:
                    # flags that we have spatial dimensions
                    # Note: that could be 'lev' which should trigger different behavior
                    # Note: we should be able to handle (lat, lon) or (ncol,) cases, at least
                    data = pf.spatial_average(data)  # changes data "in place"

                """# In order to get correct statistics, average to annual or seasonal
                data = pf.annual_mean(data, whole_years=True, time_name='time')

                # create a dataframe:
                cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                        'standard error', '95% CI', 'trend', 'trend p-value']

                # These get written to our output file:
                stats_list = _get_row_vals(data)
                row_values = [var, unit_str] + stats_list"""
                #if var == "RESTOM":
                #    print("data before annual mean",data,"\n")
                #    print(len(data),"\n\n")
                

                if is_climo:
                    data = pf.seasonal_mean(data, season="ANN", is_climo=True)
                    #Conditional Formatting depending on type of float
                    if np.abs(data) < 1:
                        formatter = ".3g"
                    else:
                        formatter = ".3f"
                    mean_final = f'{data:{formatter}}'

                    # create a dataframe:
                    cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                                'standard error', '95% CI', 'trend', 'trend p-value']
                    row_values = [var, unit_str] + [mean_final] + ["-","-","-","-","-","-"]
                else:
                    # In order to get correct statistics, average to annual or seasonal
                    data = pf.annual_mean(data, whole_years=True, time_name='time')
                    #if var == "RESTOM":
                    #    print("data AFTER annual mean",data,"\n")
                    #    print(len(data),"\n\n")
                    # create a dataframe:
                    cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                                'standard error', '95% CI', 'trend', 'trend p-value']
                    stats_list = _get_row_vals(data)
                    row_values = [var, unit_str] + stats_list
                #End if

                # Format entries:
                dfentries = {c:[row_values[i]] for i,c in enumerate(cols)}

                # Add entries to Pandas structure:
                df = pd.DataFrame(dfentries)

                # Check if the output CSV file exists,
                # if so, then append to it:
                if output_csv_file.is_file():
                    df.to_csv(output_csv_file, mode='a', header=False, index=False)
                else:
                    df.to_csv(output_csv_file, header=cols, index=False)

            #End of var_list loop
            #--------------------

            # Move RESTOM to top of table (if applicable)
            #--------------------------------------------
            try:
                table_df = pd.read_csv(output_csv_file)
                if 'RESTOM' in table_df['variable'].values:
                    table_df = pd.concat([table_df[table_df['variable'] == 'RESTOM'], table_df]).reset_index(drop = True)
                    table_df = table_df.drop_duplicates()
                    table_df.to_csv(output_csv_file, header=cols, index=False)

                # last step is to add table dataframe to website (if enabled):
                adf.add_website_data(table_df, case_name, case_name, plot_type="Tables")
            except FileNotFoundError:
                print(f"\n\tAMWG table for '{case_name}' not created.\n")
            #End try/except

            #Keep track of case csv files for comparison table check later
            csv_list.extend(sorted(output_location.glob(f"amwg_table_{case_name}.csv")))

        #End of model case loop
        #----------------------

        """#Start case comparison tables
        #----------------------------
        #Check if observations are being compared to, if so skip table comparison...
        if not adf.get_basic_info("compare_obs"):
            #Check if all tables were created to compare against, if not, skip table comparison...
            if len(csv_list) != len(case_names):
                print("\tNot enough cases to compare, skipping comparison table...")
            else:
                #Create comparison table for both cases
                print("\n  Making comparison table...")
                _df_comp_table(adf, csv_locs, case_names)
                print("  ... Comparison table has been generated successfully")
            #End if
        else:
            print(" No comparison table will be generated due to running against obs.")
        #End if

        #Notify user that script has ended:
        print("  ...AMWG variable table(s) have been generated successfully.")"""



        #Start case comparison tables
        #----------------------------
        # Copy the file to all individual directories
        if len(test_case_names) > 1:
            base_csv = sorted(Path(output_locs[-1]).glob(f"amwg_table_{baseline_name}.csv"))
            for i,case in enumerate(test_case_names):
                shutil.copy(base_csv[0], output_locs[i])
            base_csv[0].unlink()

        #Check if observations are being compared to, if so skip table comparison...
        if not adf.get_basic_info("compare_obs"):
            #Check if all tables were created to compare against, if not, skip table comparison...
            if len(csv_list) != len(case_names):
                print("\tNot enough cases to compare, skipping comparison table...")
            else:
                if len(test_case_names) == 1:
                    #Create comparison table for both cases
                    print("\n  Making comparison table...")
                    _df_comp_table(adf, output_location, Path(output_locs[0]), case_names)
                    print("  ... Comparison table has been generated successfully")

                if len(test_case_names) > 1:
                    print("\n  Making comparison table for multiple cases...")
                    _df_multi_comp_table(adf, csv_locs, case_names, nicknames)
                    print("\n  Making comparison table for each case...")
                    for idx,case in enumerate(case_names[0:-1]):
                        _df_comp_table(adf, Path(output_locs[idx]), Path(output_locs[0]), [case,baseline_name])
                    print("  ... Multi-case comparison table has been generated successfully")
            #End if
        else:
            print(" No comparison table will be generated due to running against obs.")
        #End if


##################
# Helper functions
##################

def _get_row_vals(data):
    # Now that data is (time,), we can do our simple stats:

    data_mean = data.data.mean()
    #Conditional Formatting depending on type of float
    if np.abs(data_mean) < 1:
        formatter = ".3g"
    else:
        formatter = ".3f"

    data_sample = len(data)
    data_std = data.std()
    data_sem = data_std / data_sample
    data_ci = data_sem * 1.96  # https://en.wikipedia.org/wiki/Standard_error
    data_trend = stats.linregress(data.year, data.values)

    stdev = f'{data_std.data.item() : {formatter}}'
    sem = f'{data_sem.data.item() : {formatter}}'
    ci = f'{data_ci.data.item() : {formatter}}'
    slope_int = f'{data_trend.intercept : {formatter}} + {data_trend.slope : {formatter}} t'
    pval = f'{data_trend.pvalue : {formatter}}'

    return [f'{data_mean:{formatter}}', data_sample, stdev, sem, ci, slope_int, pval]

#####



def _df_comp_table(adf, output_location, case_names):
    import pandas as pd

    output_csv_file_comp = output_location / "amwg_table_comp.csv"

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #This will be for single-case for now (case_names[0]),
    #will need to change to loop as multi-case is introduced
    case = output_location/f"amwg_table_{case_names[0]}.csv"
    baseline = output_location/f"amwg_table_{case_names[-1]}.csv"

    #Read in test case and baseline dataframes:
    df_case = pd.read_csv(case)
    df_base = pd.read_csv(baseline)

    #Create a merged dataframe that contains only the variables
    #contained within both the test case and the baseline:
    df_merge = pd.merge(df_case, df_base, how='inner', on=['variable'])

    #Create the "comparison" dataframe:
    df_comp = pd.DataFrame(dtype=object)
    df_comp[['variable','unit','case']] = df_merge[['variable','unit_x','mean_x']]
    df_comp['baseline'] = df_merge[['mean_y']]

    diffs = df_comp['case'].values-df_comp['baseline'].values
    df_comp['diff'] = [f'{i:.3g}' if np.abs(i) < 1 else f'{i:.3f}' for i in diffs]

    #Write the comparison dataframe to a new CSV file:
    cols_comp = ['variable', 'unit', 'test', 'control', 'diff']
    df_comp.to_csv(output_csv_file_comp, header=cols_comp, index=False)

    #Add comparison table dataframe to website (if enabled):
    adf.add_website_data(df_comp, "Case Comparison", case_names[0], plot_type="Tables")






def _df_multi_comp_table(adf, csv_locs, case_names, test_nicknames):
    """
    Function to build comparison AMWG table for all cases
    ------
        - Read in each previously made table from file
        and compile full comparison.
    """

    #Create path to main website in mutli-case directory
    main_site_path = Path(adf.get_basic_info('cam_diag_plot_loc', required=True))
    output_csv_file_comp = main_site_path / "amwg_table_comp_all.csv"

    #Create the "comparison" dataframe:
    df_comp = pd.DataFrame(dtype=object)

    #Create new colummns
    cols_comp = ['variable', 'unit']

    #Read baseline case
    baseline = str(csv_locs[0])+f"/amwg_table_{case_names[-1]}.csv"
    df_base = pd.read_csv(baseline)

    #Read all test cases and add to table
    for i,val in enumerate(csv_locs[:-1]):
        case = str(val)+f"/amwg_table_{case_names[i]}.csv"
        df_case = pd.read_csv(case)

        #If no custom nicknames, shorten column name to case number
        if test_nicknames[i] == case_names[i]:
            df_comp[['variable','unit',f"case {i+1}"]] = df_case[['variable','unit','mean']]
            cols_comp.append(f"case {i+1}")
        #Else, name columns after nicknames
        else:
            df_comp[['variable','unit',f"{test_nicknames[i]}"]] = df_case[['variable','unit','mean']]
            cols_comp.append(test_nicknames[i])

    #Add baseline cases to end of the table
    if test_nicknames[-1] == case_names[-1]:
        df_comp["baseline"] = df_base[['mean']]
        cols_comp.append("baseline")
    else:
        df_comp[f"{test_nicknames[-1]} ( baseline )"] = df_base[['mean']]
        cols_comp.append(f"{test_nicknames[-1]} ( baseline )")

    #Format the floats:
    for col in df_comp.columns:
        #Ignore columns that don't contain floats
        if (col != 'variable') and (col != "unit"):
            if "baseline" not in col:

                #Iterate over rows and check magnitude of value
                for idx,row in enumerate(df_comp[col]):
                    #Check if value is less than one, keep 3 non-zero decimal values
                    #Else, keep 3 main digits, including decimal values
                    if np.abs(df_comp[col][idx]) < 1:
                        formatter = ".3g"
                    else:
                        formatter = ".3f"
                    #Replace value in dataframe
                    df_comp.at[idx,col]= f'{df_comp[col][idx]:{formatter}}   ({(df_comp[col][idx]-df_base["mean"][idx]):{formatter}})'

    #Finally, write data to csv
    df_comp.to_csv(output_csv_file_comp, header=cols_comp, index=False)

    #Add comparison table dataframe to website (if enabled):
    adf.add_website_data(df_comp, "all_case_comparison", case_names[0], plot_type="Tables")

'''import builtins

def set_warning_filter(enable=True):
    """Enable or disable filtering of print statements containing 'WARNING'."""
    original_print = builtins.print

    def filtered_print(*args, **kwargs):
        message = " ".join(map(str, args))
        if enable and "WARNING" in message:
            return  # Skip printing warnings
        original_print(*args, **kwargs)

    builtins.print = filtered_print if enable else original_print'''

##############
#END OF SCRIPT