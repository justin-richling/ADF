import numpy as np
import xarray as xr
import warnings
import sys
from pathlib import Path
from collections import OrderedDict
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

    #Grab all case nickname(s)
    test_nicknames = adf.case_nicknames["test_nicknames"]
    base_nickname = adf.case_nicknames["base_nickname"]
    all_nicknames = test_nicknames + [base_nickname]
    
    input_ts_locs = adf.get_cam_info("cam_ts_loc", required=True)

    #Check if a baseline simulation is also being used:
    if not adf.get_basic_info("compare_obs"):
        #Extract CAM baseline variaables:
        baseline_name     = adf.get_baseline_info("cam_case_name", required=True)
        input_ts_baseline = adf.get_baseline_info("cam_ts_loc", required=True)

        if "CMIP" in baseline_name:
            print("CMIP files detected, skipping AMWG table (for now)...")

        else:
            #Append to case list:
            case_names.append(baseline_name)
            input_ts_locs.append(input_ts_baseline)

        #Save the baseline to the first case's plots directory:
        output_locs.append(output_locs[0])

    #Declare any derived quantities here:
    derived_vars = {}
    derived_var_list = adf.derived_var_list
    for derived_var in derived_var_list:
        derived_vars[derived_var] = var_defaults[derived_var]["constituents"]

    derived_consts_list = [item for sublist in derived_vars.values() for item in sublist]

    #Make list of all constituents of derived variables
    constituents = []
    for const_set in derived_vars.values():
        for consts in const_set:
            constituents.append(consts)

    #Create (empty) dictionary to use for the derived calculations:
    derived_dict = {}

    #Hold output paths for csv files
    csv_locs = []

    #Loop over CAM cases:
    for case_idx, case_name in enumerate(case_names):

        #Convert output location string to a Path object:
        output_location = Path(output_locs[case_idx])

        #Generate input file path:
        input_location = Path(input_ts_locs[case_idx])

        #Add output paths for csv files
        csv_locs.append(output_locs[case_idx])

        #Check that time series input directory actually exists:
        if not input_location.is_dir():
            errmsg = f"Time series directory '{input_location}' not found.  Script is exiting."
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

        #Save case name as a new key in the derived quantities dictonary:
        derived_dict[case_name] = {}

        #Create/reset new variable that potentially stores the re-gridded
        #ocean fraction xarray data-array:
        ocn_frc_da = None

        #Loop over CAM output variables:
        for var in var_list:

            #Notify users of variable being added to table:
            print(f"\t - Variable '{var}' being added to table")

            #Create list of time series files present for variable:
            ts_filenames = f'{case_name}.*.{var}.*nc'
            ts_files = sorted(input_location.glob(ts_filenames))                

            # If no files exist, try to move to next variable. --> Means we can not proceed with this variable, and it'll be problematic later.
            if not ts_files:
                errmsg = f"Time series files for variable '{var}' not found.  Script will continue to next variable."
                warnings.warn(errmsg)
                continue
            #End if

            #TEMPORARY:  For now, make sure only one file exists:
            if len(ts_files) != 1:
                errmsg =  "Currently the AMWG table script can only handle one time series file per variable."
                errmsg += f" Multiple files were found for the variable '{var}'"
                raise AdfError(errmsg)
            #End if

            #Load model data from file:
            data = _load_data(ts_files[0], var)

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
                data = _spatial_average(data)  # changes data "in place"

            #Add necessary data for derived calcs below
            if var in derived_consts_list:
                derived_dict[case_name][var] = [data, unit_str]

            # In order to get correct statistics, average to annual or seasonal
            data = data.groupby('time.year').mean(dim='time') # this should be fast b/c time series should be in memory
                                                              # NOTE: data will now have a 'year' dimension instead of 'time'

            # create a dataframe:
            cols = ['variable', 'unit', 'mean', 'sample size', 'standard dev.',
                    'standard error', '95% CI', 'trend', 'trend p-value']

            # These get written to our output file:
            stats_list = _get_row_vals(data)
            row_values = [var, unit_str] + stats_list

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
        
        #Space for derived quantities
        #----------------------------

        #Variable Difference derived quantaties (ie RESTOM, etc.)
        #-------
        _derive_diff_var(case_name, derived_dict, derived_vars, output_csv_file, cols)

        # last step is to add table dataframe to website (if enabled):
        table_df = pd.read_csv(output_csv_file)

        #Reorder RESTOM to top of tables
        idx = table_df.index[table_df['variable'] == 'RESTOM'].tolist()[0]
        table_df = pd.concat([table_df[table_df['variable'] == 'RESTOM'], table_df]).reset_index(drop = True)
        table_df = table_df.drop([idx+1]).reset_index(drop=True)
        table_df = table_df.drop_duplicates()

        adf.add_website_data(table_df, case_name, case_name, plot_type="Tables")
        #End derived quantities

        
    #End of model case loop
    #----------------------
    test_case_names = adf.get_cam_info("cam_case_name", required=True)

    #Notify user that script has ended:
    print("  ...AMWG variable table has been generated successfully.")

    #Check if observations are being comapred to, if so skip table comparison...
    if not adf.get_basic_info("compare_obs"):
        if "CMIP" in baseline_name:
            print("CMIP case detected, skipping comparison table...")
        else:

            if len(test_case_names) == 1:
                #Create comparison table for both cases
                print("\n  Making comparison table...")
                _df_comp_table(adf, output_location, Path(output_locs[0]), case_names)
                print("  ... Comparison table has been generated successfully")

            if len(test_case_names) > 1:
                print("\n  Making comparison table for multiple cases...")
                _df_multi_comp_table(adf, csv_locs, case_names, all_nicknames)
                print("\n  Making comparison table for each case...")
                for idx,case in enumerate(case_names[0:-1]):
                    _df_comp_table(adf, Path(output_locs[idx]), Path(output_locs[0]), [case,baseline_name])
                print("  ... Multi-case comparison table has been generated successfully")
        #End if
    else:
        print(" Comparison table currently doesn't work with obs, so skipping...")
    #End if


##################
# Helper functions
##################

def _load_data(dataloc, varname):
    ds = xr.open_dataset(dataloc)
    return ds[varname]

#####

def _spatial_average(indata):

    #Make sure there is no veritcal level dimension:
    assert 'lev' not in indata.coords
    assert 'ilev' not in indata.coords

    #Calculate spatial weights:
    if 'lat' in indata.coords:
        weights = np.cos(np.deg2rad(indata.lat))
        weights.name = "weights"
    elif 'ncol' in indata.coords:
        warnings.warn("We need a way to get area variable. Using equal weights.")
        weights = xr.DataArray(1.)
        weights.name = "weights"
    else:
        weights = xr.DataArray(1.)
        weights.name = "weights"
    #End if

    #Apply weights to input data:
    weighted = indata.weighted(weights)

    # we want to average over all non-time dimensions
    avgdims = [dim for dim in indata.dims if dim != 'time']
    return weighted.mean(dim=avgdims)

#####

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

def _df_comp_table(adf, output_location, base_output_location, case_names):
    """
    Function to build case vs baseline AMWG table
    -----
        - Read in table data and create side by side comaprison table
        - Write output to csv file and add to website
    """
    
    output_csv_file_comp = output_location / "amwg_table_comp.csv"

    case = output_location/f"amwg_table_{case_names[0]}.csv"
    baseline = base_output_location/f"amwg_table_{case_names[-1]}.csv"

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
    cols_comp = ['variable', 'unit', 'test', 'baseline', 'diff']

    #Reorder RESTOM to top of tables
    idx = df_comp.index[df_comp['variable'] == 'RESTOM'].tolist()[0]
    df_comp = pd.concat([df_comp[df_comp['variable'] == 'RESTOM'], df_comp]).reset_index(drop = True)
    df_comp = df_comp.drop([idx+1]).reset_index(drop=True)
    df_comp = df_comp.drop_duplicates()

    df_comp.to_csv(output_csv_file_comp, header=cols_comp, index=False)

    #Add comparison table dataframe to website (if enabled):
    adf.add_website_data(df_comp, "case_comparison", case_names[0], plot_type="Tables")

#####

def _df_multi_comp_table(adf, csv_locs, case_names, all_nicknames):
    """
    Function to build all case comparison AMWG table
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
    baseline = str(csv_locs[-1])+f"/amwg_table_{case_names[-1]}.csv"
    df_base = pd.read_csv(baseline)

    #Read all test cases and add to table
    for i,val in enumerate(csv_locs[:-1]): 
        case = str(val)+f"/amwg_table_{case_names[i]}.csv"
        df_case = pd.read_csv(case)
        
        #If no custom nicknames, shorten column name to case number
        if all_nicknames[i] == case_names[i]:
            df_comp[['variable','unit',f"case {i+1}"]] = df_case[['variable','unit','mean']]
            cols_comp.append(f"case {i+1}")
        #Else, name columns after nicknames
        else:
            df_comp[['variable','unit',f"{all_nicknames[i]}"]] = df_case[['variable','unit','mean']]
            cols_comp.append(all_nicknames[i])

    #Add baseline cases to end of the table
    if all_nicknames[-1] == case_names[-1]:
        df_comp["baseline"] = df_base[['mean']]
        cols_comp.append("baseline")
    else:
        df_comp[f"{all_nicknames[-1]} ( baseline )"] = df_base[['mean']]
        cols_comp.append(f"{all_nicknames[-1]} ( baseline )")

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

    #Reorder RESTOM to top of tables
    idx = df_comp.index[df_comp['variable'] == 'RESTOM'].tolist()[0]
    df_comp = pd.concat([df_comp[df_comp['variable'] == 'RESTOM'], df_comp]).reset_index(drop = True)
    df_comp = df_comp.drop([idx+1]).reset_index(drop=True)
    df_comp = df_comp.drop_duplicates()

    #Finally, write data to csv
    df_comp.to_csv(output_csv_file_comp, header=cols_comp, index=False)

    #Add comparison table dataframe to website (if enabled):
    adf.add_website_data(df_comp, "all_case_comparison", case_names[0], plot_type="Tables")

#####

#Derived quantity function space
################################


#_derive_diff_var(case_name, derived_dict, output_csv_file, cols)
#def _derive_diff_var(case_name, var, var_consts, derived_dict, output_csv_file, cols):
def _derive_diff_var(case_name, derived_dict, derived_vars, output_csv_file, cols):
    """
    derived_vars -> dictioanry that houses derived variable as key and constituents as values
    derived_dict -> dictionary that houses consituents (all) as keys and data as value
                    * each key has data
                    derived_dict[case_name][var] = [data, unit_str]
    """
    
    for der_var,consts in derived_vars.items():
        #var = "RESTOM" #RESTOM = FSNT-FLNT
        print(f"\t - Variable '{der_var}' being added to table")


        #print("YAHHOOO",derived_dict[case_name][consts[0]][0])
        data = derived_dict[case_name][consts[0]][0]
        for consts_var in consts[1:]:
            data -= derived_dict[case_name][consts_var][0]

        # In order to get correct statistics, average to annual or seasonal
        data = data.groupby('time.year').mean(dim='time') # this should be fast b/c time series should be in memory
                                                                    # NOTE: data will now have a 'year' dimension instead of 'time'
        # These get written to our output file:
        stats_list = _get_row_vals(data)
        #Extract units string, if available:
        if hasattr(data, 'units'):
            unit_str = data.units
        else:
            unit_str = '--'
        #End if
        
        row_values = [der_var, unit_str] + stats_list

        # Format entries:
        #NOTE: col (column) values were declared above
        dfentries = {c:[row_values[i]] for i,c in enumerate(cols)}

        # Add entries to Pandas structure:
        df = pd.DataFrame(dfentries)

        # Check if the output CSV file exists,
        # if so, then append to it:
        if output_csv_file.is_file():
            df.to_csv(output_csv_file, mode='a', header=False, index=False)
        else:
            df.to_csv(output_csv_file, header=cols, index=False)
        #End if
                
        table_df = pd.read_csv(output_csv_file)

        """#Reorder RESTOM to top of tables
        idx = table_df.index[table_df['variable'] == 'RESTOM'].tolist()[0]
        table_df = pd.concat([table_df[table_df['variable'] == 'RESTOM'], table_df]).reset_index(drop = True)
        table_df = table_df.drop([idx+1]).reset_index(drop=True)
        table_df = table_df.drop_duplicates()"""

        #Re-save the csv file
        table_df.to_csv(output_csv_file, header=cols, index=False)


# RESTOM
def _derive_restom(case_name, derived_dict, output_csv_file, cols):
    
    var = "RESTOM" #RESTOM = FSNT-FLNT
    print(f"\t - Variable '{var}' being added to table")
    data = derived_dict[case_name]["FSNT"][0] - derived_dict[case_name]["FLNT"][0]
     # In order to get correct statistics, average to annual or seasonal
    data = data.groupby('time.year').mean(dim='time') # this should be fast b/c time series should be in memory
                                                                # NOTE: data will now have a 'year' dimension instead of 'time'
    # These get written to our output file:
    stats_list = _get_row_vals(data)
    #Extract units string, if available:
    if hasattr(data, 'units'):
        unit_str = data.units
    else:
        unit_str = '--'
    #End if
    
    row_values = [var, unit_str] + stats_list

    # Format entries:
    #NOTE: col (column) values were declared above
    dfentries = {c:[row_values[i]] for i,c in enumerate(cols)}

    # Add entries to Pandas structure:
    df = pd.DataFrame(dfentries)

    # Check if the output CSV file exists,
    # if so, then append to it:
    if output_csv_file.is_file():
        df.to_csv(output_csv_file, mode='a', header=False, index=False)
    else:
        df.to_csv(output_csv_file, header=cols, index=False)
    #End if
            
    table_df = pd.read_csv(output_csv_file)

    #Reorder RESTOM to top of tables
    idx = table_df.index[table_df['variable'] == 'RESTOM'].tolist()[0]
    table_df = pd.concat([table_df[table_df['variable'] == 'RESTOM'], table_df]).reset_index(drop = True)
    table_df = table_df.drop([idx+1]).reset_index(drop=True)
    table_df = table_df.drop_duplicates()

    #Re-save the csv file
    table_df.to_csv(output_csv_file, header=cols, index=False)

######

##############
#END OF SCRIPT