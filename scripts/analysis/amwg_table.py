import numpy as np
import xarray as xr
import sys
from pathlib import Path
import warnings  # use to warn user about missing files.

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
    print("\n  Calculating AMWG variable table...")


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
    case_names    = adf.get_cam_info("cam_case_name", required=True)
    
    #Check if time series location was provided
    # * time series files need to be made from history files (save to this location)
    # * premade time series files were supplied (use this location)
    input_ts_locs = adf.get_cam_info("cam_ts_loc")
    ts_locs = {}
    if not input_ts_locs:
        for case in case_names:
            ts_locs[case] = None
    else:
        for i,case in enumerate(case_names):
            if input_ts_locs[i]:
                ts_locs[case] = input_ts_locs[i]
            else:
                ts_locs[case] = None

    """#Check if time series files need to be calculated
    calc_cam_ts   = adf.get_cam_info("calc_cam_ts")
    calc_ts = {}
    if not calc_cam_ts:
        for case in case_names:
            calc_ts[case] = False
    else:
        for i,case in enumerate(case_names):
            if calc_cam_ts[i]:
                calc_ts[case] = calc_cam_ts[i]
            else:
                calc_ts[case] = False"""

    #Check if climo location was provided
    # * climo files need to be made from time series (save to this location)
    # * premade climo files were supplied (use this location)
    input_climo_locs = adf.get_cam_info("cam_climo_loc")
    climo_locs = {}
    if not input_climo_locs:
        for case in case_names:
            climo_locs[case] = None
    else:
        for i,case in enumerate(case_names):
            if input_climo_locs[i]:
                climo_locs[case] = input_climo_locs[i]
            else:
                climo_locs[case] = None

    #Grab case years
    syear_cases = adf.climo_yrs["syears"]
    eyear_cases = adf.climo_yrs["eyears"]

    #Grab baseline years (which may be empty strings if using Obs):
    syear_baseline = adf.climo_yrs["syear_baseline"]
    eyear_baseline = adf.climo_yrs["eyear_baseline"]

    syear_cases.append(syear_baseline)
    eyear_cases.append(eyear_baseline)

    #Check if a baseline simulation is also being used:
    if not adf.get_basic_info("compare_obs"):
        #Extract CAM baseline variaables:
        baseline_name     = adf.get_baseline_info("cam_case_name", required=True)
        case_names.append(baseline_name)

        #Check if time series location was provided
        # * time series files need to be made from history files (save to this location)
        # * premade time series files were supplied (use this location)
        input_ts_baseline = adf.get_baseline_info("cam_ts_loc")
        if input_ts_baseline:
            ts_locs[baseline_name] = input_ts_baseline
        else:
            ts_locs[baseline_name] = None
        #End if

        """#Check if time series files need to be calculated
        calc_baseline_ts   = adf.get_baseline_info("calc_cam_ts")
        if calc_baseline_ts:
            calc_ts[baseline_name] = calc_baseline_ts
        else:
            calc_ts[baseline_name] = False
        #End if"""

        #Check if climo location was provided
        # * climo files need to be made from time series (save to this location)
        # * premade climo files were supplied (use this location)
        input_base_climo_loc = adf.get_baseline_info("cam_climo_loc")
        if input_base_climo_loc:
            climo_locs[baseline_name] = input_base_climo_loc
        else:
            climo_locs[baseline_name] = None
        #End if

        #Save the baseline to the first case's plots directory:
        output_locs.append(output_locs[0])
    else:
        print("AMWG table doesn't currently work with obs, so obs table won't be created.")
    #End if

    print("ts_locs",ts_locs,"\n")
    print("climo_locs",climo_locs)
    #print("calc_ts",calc_ts,"\n")

    #-----------------------------------------

    #Check if user wants to skip time series file creation
    #calc_cam_ts   = adf.get_cam_info("calc_cam_ts")
    #if not isinstance(calc_cam_ts, list):
    #    # If so, then check if any of the entries are "True":
    #    calc_cam_ts = list(calc_cam_ts)
    # End if

    #if calc_cam_ts is None:
    #    calc_cam_ts = [False]*len(case_names)

    #Loop over CAM cases:
    #Initialize list of case name csv files for case comparison check later
    csv_list = []
    for case_idx, case_name in enumerate(case_names):
        print(f"Making AMWG table for case'{case_name}'")
        """if calc_ts[case_name]:
            is_climo = False
        else:
            print("User supplied Climo files, will make only global mean for each variable. Thanks and have a nice day.")
            is_climo = True"""

        """if not calc_ts[case_name]:
            print("User supplied Climo files, will make only global mean for each variable. Thanks and have a nice day.")
            is_climo = True
        else:
            is_climo = False"""

        
        if not ts_locs[case_name]:
            print("User supplied Climo files, will make only global mean for each variable. Thanks and have a nice day.")
            is_climo = True
        else:
            is_climo = False

        print("\nis_climo:",is_climo,"\n")

        #Convert output location string to a Path object:
        output_location = Path(output_locs[case_idx])

        #Generate input file path:
        if not is_climo:
            input_location = Path(ts_locs[case_name])
        if is_climo:
            input_location = Path(climo_locs[case_name])
        #print("input_location",input_location,"\n")

        #Check that time series/climo input directory actually exists:
        if not input_location.is_dir():
            errmsg = f"Directory '{input_location}' not found.  Script is exiting."
            raise AdfError(errmsg)
        #Write to debug log if enabled:
        adf.debug_log(f"DEBUG: location of files is {str(input_location)}")
        #Check if analysis directory exists, and if not, then create it:
        if not output_location.is_dir():
            print(f"\t    {output_locs[case_idx]} not found, making new directory")
            output_location.mkdir(parents=True)

        #Create output file name:
        output_csv_file = output_location / f"amwg_table_{case_name}.csv"

        #Given that this is a final, user-facing analysis, go ahead and re-do it every time:
        if Path(output_csv_file).is_file():
            Path.unlink(output_csv_file)
        #End if
        
        
        #Make and save the table to CSV file. Keep track of the file too for comparison table
        #csv_list = make_table(adf, var_list, case_name, input_location, var_defaults, output_csv_file, output_location, csv_list, premade_climo=is_climo)

        #Create/reset new variable that potentially stores the re-gridded
        #ocean fraction xarray data-array:
        ocn_frc_da = None

        #Loop over CAM output variables:
        for var in var_list:

            #Notify users of variable being added to table:
            print(f"\t - Variable '{var}' being added to table")

            if is_climo:
                #Create list of climo files present for variable:
                filenames = f'{case_name}_{var}_climo.nc'
            else:
                #Create list of time series files present for variable:
                filenames = f'{case_name}.*.{var}.*nc'
            files = sorted(input_location.glob(filenames))
            #print(f"TABLES for {case_name}")
            #print("input_location",input_location)
            #print("filenames",filenames,"\n")

            # If no files exist, try to move to next variable. --> Means we can not proceed with this variable, and it'll be problematic later.
            if not files:
                errmsg = f"Time series files for variable '{var}' not found.  Script will continue to next variable."
                warnings.warn(errmsg)
                continue
            #End if

            #TEMPORARY:  For now, make sure only one file exists:
            if len(files) != 1:
                errmsg =  "Currently the AMWG table script can only handle one time series file per variable."
                errmsg += f" Multiple files were found for the variable '{var}', so it will be skipped."
                print(errmsg)
                continue
            #End if

            #Load model variable data from file:
            ds = pf.load_dataset(files)
            data = ds[var]

            if not is_climo:
                data = fixcesmtime(data,syear_cases[case_idx],eyear_cases[case_idx])

            #Extract units string, if available:
            if hasattr(data, 'units'):
                unit_str = data.units
            else:
                unit_str = '--'

            #Check if variable has a vertical coordinate:
            if 'lev' in data.coords or 'ilev' in data.coords:
                print(f"\t   Variable '{var}' has a vertical dimension, "+\
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
                        print(f"OCNFRAC not found, unable to apply mask to '{var}'")
                    #End if
                else:
                    #Currently only an ocean mask is supported, so print warning here:
                    wmsg = "Currently the only variable mask option is 'ocean',"
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
            #End if
            


            if is_climo:
                data = pf.seasonal_mean(data, season="ANN", is_climo=True)
                #Conditional Formatting depending on type of float
                if np.abs(data) < 1:
                    formatter = ".3g"
                else:
                    formatter = ".3f"
                mean_final = f'{data:{formatter}}'

                # create a dataframe:
                cols = ['variable', 'unit', 'mean']
                row_values = [var, unit_str] + [mean_final]
            else:
                # In order to get correct statistics, average to annual or seasonal
                data = pf.annual_mean(data, whole_years=True, time_name='time')
                # create a dataframe:
                cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                            'standard error', '95% CI', 'trend', 'trend p-value']
                stats_list = _get_row_vals(data)
                row_values = [var, unit_str] + stats_list
            #End if

            """# These get written to our output file:
            stats_list = _get_row_vals(data)
            row_values = [var, unit_str] + stats_list



            # In order to get correct statistics, average to annual or seasonal
            data = pf.annual_mean(data, whole_years=True, time_name='time')

            # create a dataframe:
            cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                    'standard error', '95% CI', 'trend', 'trend p-value']

            # These get written to our output file:
            stats_list = _get_row_vals(data)
            row_values = [var, unit_str] + stats_list"""

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

    #Start case comparison tables
    #----------------------------
    #Check if observations are being compared to, if so skip table comparison...
    if not adf.get_basic_info("compare_obs"):
        #Check if all tables were created to compare against, if not, skip table comparison...
        if len(csv_list) != len(case_names):
            print("\tNot enough cases to compare, skipping comparison table...")
        else:
            #Create comparison table for both cases
            print("\n  Making comparison table...")
            _df_comp_table(adf, output_location, case_names)
            print("  ... Comparison table has been generated successfully")
        #End if
    else:
        print(" No comparison table will be generated due to running against obs.")
    #End if

    #Notify user that script has ended:
    print("  ...AMWG variable table(s) have been generated successfully.")


##################
# Helper functions
##################

def make_table(adf, var_list, case_name, input_location, var_defaults, output_csv_file, output_location, csv_list, premade_climo=False):
    #Create/reset new variable that potentially stores the re-gridded
    #ocean fraction xarray data-array:
    ocn_frc_da = None

    #Loop over CAM output variables:
    for var in var_list:

        #Notify users of variable being added to table:
        print(f"\t - Variable '{var}' being added to table")

        if premade_climo:
            #Create list of climo files present for variable:
            filenames = f'{case_name}_{var}_climo.nc'
        else:
            #Create list of time series files present for variable:
            filenames = f'{case_name}.*.{var}.*nc'
        files = sorted(input_location.glob(filenames))
        #print(f"TABLES for {case_name}")
        #print("input_location",input_location)
        #print("filenames",filenames,"\n")

        # If no files exist, try to move to next variable. --> Means we can not proceed with this variable, and it'll be problematic later.
        if not files:
            errmsg = f"Time series files for variable '{var}' not found.  Script will continue to next variable."
            warnings.warn(errmsg)
            continue
        #End if

        #TEMPORARY:  For now, make sure only one file exists:
        if len(files) != 1:
            errmsg =  "Currently the AMWG table script can only handle one time series file per variable."
            errmsg += f" Multiple files were found for the variable '{var}', so it will be skipped."
            print(errmsg)
            continue
        #End if

        #Load model variable data from file:
        ds = pf.load_dataset(files)
        data = ds[var]

        data = fixcesmtime(data,start_years[idx],end_years[idx])

        #Extract units string, if available:
        if hasattr(data, 'units'):
            unit_str = data.units
        else:
            unit_str = '--'

        #Check if variable has a vertical coordinate:
        if 'lev' in data.coords or 'ilev' in data.coords:
            print(f"\t   Variable '{var}' has a vertical dimension, "+\
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
                    print(f"OCNFRAC not found, unable to apply mask to '{var}'")
                #End if
            else:
                #Currently only an ocean mask is supported, so print warning here:
                wmsg = "Currently the only variable mask option is 'ocean',"
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

        if premade_climo:
            data = pf.seasonal_mean(data, season="ANN", is_climo=True)
            #Conditional Formatting depending on type of float
            if np.abs(data) < 1:
                formatter = ".3g"
            else:
                formatter = ".3f"
            mean_final = f'{data:{formatter}}'

            # create a dataframe:
            cols = ['variable', 'unit', 'mean']
            row_values = [var, unit_str] + [mean_final]

        else:
            # In order to get correct statistics, average to annual or seasonal
            data = pf.annual_mean(data, whole_years=True, time_name='time')
            # create a dataframe:
            cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                        'standard error', '95% CI', 'trend', 'trend p-value']

            # These get written to our output file:
            stats_list = _get_row_vals(data)
            row_values = [var, unit_str] + stats_list

        """# create a dataframe:
        cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                    'standard error', '95% CI', 'trend', 'trend p-value']

        # These get written to our output file:
        stats_list = _get_row_vals(data)
        row_values = [var, unit_str] + stats_list"""

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
    return csv_list


def fixcesmtime(dat,syear,eyear):
    """
    Fix the CESM timestamp with a simple set of dates
    """
    timefix = pd.date_range(start=f'1/1/{syear}', end=f'12/1/{eyear}', freq='MS') # generic time coordinate from a non-leap-year
    dat = dat.assign_coords({"time":timefix})

    return dat


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

##############
#END OF SCRIPT