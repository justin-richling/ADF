#Import "special" modules:
try:
    import pandas as pd
except ImportError:
    print("Pandas module does not exist in python path, but is needed for amwg_table.")
    print("Please install module, e.g. 'pip install pandas'.")
    sys.exit(1)
try:
    import scipy.stats as stats # for easy linear regression and testing
except ImportError:
    print("Scipy module does not exist in python path, but is needed for amwg_table.")
    print("Please install module, e.g. 'pip install scipy'.")

from datetime import datetime, timedelta
import time
import timeit
import numpy as np
#import pandas as pd
import xarray as xr
from pathlib import Path
import sys
import warnings  # use to warn user about missing files.

# list of the variables to be caculated. 
CHEMS =["CH4",
        "CH3CCL3",
        "CO",
        "O3", #Currently missing some components 
        "ISOP",
        "C10H16",
        "CH3OH",
        "CH3COCH3"]

AEROSOLS = ['SOA', 'SALT', 'DUST', 'POM', 'BC', 'SULF']

def amwg_chem_table(adf):

    """
    Main function goes through series of steps:
    - load the variable data
    - Determine whether there are spatial dims; if yes, do global average (TODO: regional option)
    - Apply annual average (TODO: add seasonal here)

    Description of needed inputs from ADF:

    case_names      -> Name(s) of CAM case provided by "cam_case_name"
    input_ts_locs   -> Location(s) of CAM time series files provided by "cam_ts_loc"
    output_loc      -> Location to write AMWG table files to, provided by "cam_diag_plot_loc"
    var_list        -> List of CAM output variables provided by "diag_var_list"

    and if doing a CAM baseline comparison:

    baseline_name     -> Name of CAM baseline case provided by "cam_case_name"
    input_ts_baseline -> Location of CAM baseline time series files provied by "cam_ts_loc"

    """

    #Import necessary modules:
    from adf_base import AdfError

    #Additional information:
    #----------------------

    # DETAIL: we use python's type hinting as much as possible

    # in future, provide option to do multiple domains
    # They use 4 pre-defined domains:
    domains = {"global": (0, 360, -90, 90),
               "tropics": (0, 360, -20, 20),
               "southern": (0, 360, -90, -20),
               "northern": (0, 360, 20, 90)}

    #Extract needed quantities from ADF object:
    #-----------------------------------------

    #Special ADF variable which contains the output paths for
    #all generated plots and tables for each case:
    output_locs = adf.plot_location

    #CAM simulation variables (these quantities are always lists):
    case_names    = adf.get_cam_info("cam_case_name", required=True)
    input_ts_locs = adf.get_cam_info("cam_ts_loc", required=True)

    start_year = adf.climo_yrs["syears"]
    end_year = adf.climo_yrs["eyears"]

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
    
    #CHANGE THIS TO READ FROM YAML!!!!
    h_case = "h0"

    #Convert output location string to a Path object:
    output_location = Path(output_locs[0])

    #Generate input file path:
    input_location = Path(input_ts_locs[0])

    #Check that time series input directory actually exists:
    if not input_location.is_dir():
        errmsg = f"Time series directory '{input_location}' not found.  Script is exiting."
        raise AdfError(errmsg)
    #Write to debug log if enabled:
    adf.debug_log(f"DEBUG: location of files is {str(input_location)}")
    #Check if analysis directory exists, and if not, then create it:
    if not output_location.is_dir():
        print(f"\t    {output_locs[0]} not found, making new directory")
        output_location.mkdir(parents=True)

    #Grab history file locations from config yaml file
    cam_hist_locs = adf.get_cam_info("cam_hist_loc", required=True)
    cam_hist_locs = cam_hist_locs + [adf.get_baseline_info("cam_hist_loc", required=True)]

    #Create path object for the CAM history file(s) location:
    data_dirs = []
    for case_idx,case in enumerate(case_names):
        data_dirs.append(cam_hist_locs[case_idx])

    #End gathering case, path, and data info
    #-----------------------------------------


    # Look for specific h-case    
    scenarios = [f'{ix}.cam.{h_case}' for ix in case_names]
    #scenarios = [f'{ix}.{h_case}' for ix in case_names]

    # In CAM-Chem (or MUSICA-v0), user can save the outputs for only a box region.
    # ext1_SE: string specifying if the files are for only a region, which changes to variable names.
    # ex: if you saved files for only a box region ($LL_lat$,$LL_lon$,$UR_lat$,$UR_lon$),
    #      the 'lat' variable will be saved as: 'lat_$LL_lon$e_to_$UR_lon$e_$LL_lat$n_to_$UR_lat$n'
    #      for instance: 'lat_65e_to_91e_20n_to_32n'
    ext1_SE=''
    #ext1_SE='_65e_to_91e_20n_to_32n'

    # Tropospheric Values
    # -------------------
    # if True, calculate only Tropospheric values
    # if False, all layers
    # tropopause is defiend as o3>150ppb. If needed, change accordingly.
    Tropospheric=True

    # Regional subset
    # ---------------
    regional=False
    #dir_shapefile="/Users/roozitalab/INDIA/Shapefile/Bangladesh/bgd_adm_bbs_20201113_shp/bgd_adm_bbs_20201113_SHP/"
    dir_shapefile="/Users/roozitalab/INDIA/Shapefile/Countries/world"
    
    # Lat/Lon extent
    limit=(20,20,40,120)

    # Periods of Interest
    # -------------------
    # choose the period of interest. Plots will be averaged within this period
    start_dates = [f"{start_year[0]}-1-1", f"{start_year[0]}-1-1"]
    end_dates = [f"{end_year[0]}-1-1", f"{end_year[0]}-1-1"]

    start_periods = []
    end_periods = []
    durations = []
    for i,val in enumerate(start_dates):
    # convert date strings to datetime format
        start_period = datetime.strptime(start_dates[i], "%Y-%m-%d")
        end_period = datetime.strptime(end_dates[i], "%Y-%m-%d")
        
        start_periods.append(start_period)
        end_periods.append(end_period)
        
        durations.append((end_period-start_period).days*86400)

    print(timeit.timeit(lambda: Get_files(data_dirs,scenarios,start_dates,end_dates,area=True), number=1),"\n")
    tic = time.perf_counter()
    #Get the files for each case and set of start and end years
    Files,Lats,Lons,areas= Get_files(data_dirs,scenarios,start_dates,end_dates,area=True)
    toc = time.perf_counter()
    print(f"Get_files took {toc - tic:0.4f} seconds")

    # Files (list) can have serveral values
    # It looks like tmp_file is only one case at a time???
    # -> just to get variable names, dont need above statement

    # find the name of all the variables in the file.
    # this will help the code to work for the variables that are not in the files (assingn 0s)
    tmp_file=xr.open_dataset(data_dirs[0]+Files[scenarios[0]][0])
    ListVars=tmp_file.variables
    tmp_file.close()
    
    # Chemistry tables
    #-----------------
    #Notify user that script has started:
    print("\n  Calculating AMWG chemistry variable table...")

    '''
    NOTE - can probably make this more efficient/correct
    '''
    #Get dict for critical values and dict for case/variables mean ANN values
    #Dic_crit,var_dict = make_var_dict(CHEMS)
    print(timeit.timeit(lambda: create_dic_SE(CHEMS, ListVars, ext1_SE), number=1),"\n")
    tic = time.perf_counter()
    dic_SE = create_dic_SE(CHEMS, ListVars, ext1_SE)
    toc = time.perf_counter()
    print(f"create_dic_SE took {toc - tic:0.4f} seconds")
        
    # extract all the data
    var_dict={}

    # this is for finding tropospheric values
    Dic_crit={}

    #Create output file name:
    output_csv_file = output_location / f"amwg_chem_table_{case_names[0]}.csv"

    """if output_csv_file.is_file():
        print(f"'{output_csv_file}' already exists, so skipping partner!\n")
        table_df = pd.read_csv(output_csv_file)
        adf.add_website_data(table_df, "Chemistry", case_names[0], plot_type="Tables")"""
    
    if 1==0:
        print("this should never run!")
        
    else:
        for i,scn in enumerate(scenarios):
            
            area=areas[scn]
            current_lat=Lats[scn]
            current_lon=Lons[scn]

            if regional:
                inside=Inside_SE(current_lat,current_lon,limit)
            else:
                if len(np.shape(area)) == 1:
                    inside=np.full((len(current_lon)),True)
                else:
                    inside=np.full((len(current_lat),len(current_lon)),True)

            current_dir=data_dirs[i]
            current_files=Files[scn] 

            var_dict[scn]={}
            Dic_var_comp={}
                
            for _,var in enumerate(CHEMS):

                # Components are: burden, chemical loss, chemical prod, dry deposition, 
                #                 surface emissions, elevated emissions, chemical tendency
                # I always add Lightning NOx production for gaseous species.
                components=[var+'_BURDEN',var+'_CHML',var+'_CHMP',
                            var+'_DDF',var+'_WDF', var+'_SF', var+'_CLXF',
                            var+'_TEND',var+'_LNO']         

                Dic_comp={}
                for comp in components:
                    current_data=SEbudget(dic_SE,current_dir,current_files,comp,level=50)
                        
                    Dic_comp[comp]=current_data
                Dic_var_comp[var]=Dic_comp
            var_dict[scn]= Dic_var_comp    
            
            print("here comes the fun cooker for SEbudget:")
            print(timeit.timeit(lambda: SEbudget(dic_SE,current_dir,current_files,'O3',level=50), number=1),"\n")

            #Critical threshholds????
            tic = time.perf_counter()
            current_crit=SEbudget(dic_SE,current_dir,current_files,'O3',level=50)
            toc = time.perf_counter()
            print(f"SEbudget took {toc - tic:0.4f} seconds")
            #print(timeit.timeit(lambda: SEbudget(dic_SE,current_dir,current_files,'O3',level=50), number=1))
            Dic_crit[scn]=current_crit

            print(f'\nCurrent Scenario: {scn}')
            print(len(f'Current Scenario: {scn}')*'-','\n')
            
            if Tropospheric:
                trop=np.where(current_crit>150,np.nan,current_crit)
                strat=np.where(current_crit>150,current_crit,np.nan)
            else:
                trop=current_crit

        comp_ext_full = {'_BURDEN':'_BURDEN',
                    '_CHML':'_CHEM_LOSS','_CHMP':'_CHEM_PROD','_NET':'_NET',
                    '_SF':'_EMIS',
                    '_DDF':'_DRYDEP','_WDF':'_WETDEP',
                    '_LIFETIME':'_LIFETIME',
                    '_STE':'_STE',
                    '_LNO':'_LNO'}

        #O3 has been causing issues, let's try and separate O3 from other variables...
        #It seems to be missing components:
        # 'BURDEN', 'CHML', 'SF','LIFETIME', and 'LNO'
        not_O3_ext = {k: v for k, v in comp_ext_full.items() if k in ['_BURDEN', '_CHML', '_SF','_LIFETIME','_LNO']}
        O3_ext = comp_ext_full.copy()

        #Create the table
        #----------------
        
        #Use this for multi-case --> down the road a bit, yeah?
        cols = ['variable']+[f"Test {i+1}" for i,_ in enumerate(case_names[0:-1])]+["Baseline"]
        
        for current_var in CHEMS:

            #Run O3 calcs
            #------------
            if current_var == "O3":
                for key,ext in O3_ext.items():
                    row_values = []
                
                    for i,scn in enumerate(scenarios):
                        tic = time.perf_counter()
                        my_val = calc_chem_data(scn,current_var,var_dict,trop,
                                                area,durations[i],inside)[key]
                        toc = time.perf_counter()
                        print(f"calc_chem_data for O3 took {toc - tic:0.4f} seconds")

                        if ext == "_BURDEN":
                            new_ext = ext+" (Tg)"
                        elif ext == "_LNO":
                            new_ext = ext+" (TgN/yr)"
                        elif ext == "_LIFETIME":
                            new_ext = ext+" (days)"
                            my_val = my_val*365
                        else:
                            new_ext = ext+" (Tg/yr)"

                        row_values.append(np.round(my_val,3))

                    row_values = [current_var+new_ext]+row_values
                    dfentries = {c:[row_values[idx]] for idx,c in enumerate(cols)}

                    # Add entries to Pandas structure:
                    df = pd.DataFrame(dfentries,columns=cols)
                    if output_csv_file.is_file():
                        df.to_csv(output_csv_file, mode='a', header=False, index=False)
                    else:
                        df.to_csv(output_csv_file, header=False, index=False)
            
            #Run most other variables
            #------------------------
            elif current_var not in ['C10H16', 'CH3OH', 'CH3COCH3', 'ISOP', "O3"]:
                for key,ext in not_O3_ext.items():
                    row_values = []
                    for i,scn in enumerate(scenarios):
                        tic = time.perf_counter()
                        my_val = calc_chem_data(scn,current_var,var_dict,trop,
                                                area,durations[i],inside)[key]
                        toc = time.perf_counter()
                        print(f"calc_chem_data for {current_var} took {toc - tic:0.4f} seconds")
                    
                        if ext == "_BURDEN":
                            new_ext = ext+" (Tg)"
                        elif ext == "_LNO":
                            new_ext = ext+" (TgN/yr)"
                        elif ext == "_LIFETIME":
                            if my_val < 1:
                                my_val = my_val*365
                                new_ext = ext+" (days)"
                            else:
                                new_ext = ext+" (yr)"
                        else:
                            new_ext = ext+" (Tg/yr)"

                        row_values.append(np.round(my_val,3))

                    row_values = [current_var+new_ext]+row_values

                    dfentries = {c:[row_values[idx]] for idx,c in enumerate(cols)}
                    # Add entries to Pandas structure:
                    df = pd.DataFrame(dfentries,columns=cols)
                    if output_csv_file.is_file():
                        df.to_csv(output_csv_file, mode='a', header=False, index=False)
                    else:
                        df.to_csv(output_csv_file, header=False, index=False)
                        
            #Run ISOP, Monoterpene, Methanol, and Acetone emmission calcs
            #------------------------------------------------------------
            elif current_var in ['C10H16', 'CH3OH', 'CH3COCH3', 'ISOP']:
                print("current_var",current_var,"\n")
                row_values = []
                new_ext = "_EMIS (Tg/yr)"

                for i,scn in enumerate(scenarios):
                    tic = time.perf_counter()
                    my_val = calc_chem_data(scn,current_var,var_dict,trop,
                                                area,durations[i],inside)['_SF']
                    toc = time.perf_counter()
                    print(f"calc_chem_data for {current_var} took {toc - tic:0.4f} seconds")
                    row_values.append(np.round(my_val,3))
                row_values = [current_var+new_ext]+row_values

                dfentries = {c:[row_values[idx]] for idx,c in enumerate(cols)}
                # Add entries to Pandas structure:
                df = pd.DataFrame(dfentries,columns=cols)
                if output_csv_file.is_file():
                    df.to_csv(output_csv_file, mode='a', header=False, index=False)
                else:
                    df.to_csv(output_csv_file, header=False, index=False)

        #Do some extracurricular work to clean up the tables
        # - mostly to try and match the old AMWG chem tables
        table_df = pd.read_csv(output_csv_file,names=cols)

        # Change some compounds to match old AMWG chem table names
        #table_df = table_df.replace('C10H16','Monoterpene', regex=True)
        #table_df = table_df.replace('CH3OH','Methanol', regex=True)
        #table_df = table_df.replace('CH3COCH3','Acetone', regex=True)

        #Looks like LNO calculated in each variable are the same value, so just use one LNO??
        # There's probably a better way to do this 
        drop_vals = ["CH3CCL3_LNO","CO_LNO","O3_LNO"]
        for val in drop_vals:
            table_df = table_df[table_df["variable"].str.contains(val) == False]
            table_df.reset_index(drop=True, inplace = True)

        # Grab one LNO value (from CH4) change to LNO_PROD and add units
        # Also, for no good reason, move to the bottom of table for completeness' sake
        table_df = table_df.replace('CH4_LNO','LNO_PROD', regex=True)
        idx = table_df.index[table_df['variable'] == 'LNO_PROD (TgN/yr)'].tolist()[0]
        table_df = table_df.append(table_df.iloc[idx], ignore_index=True)
        table_df = table_df.drop([idx]).reset_index(drop=True)

        table_df = table_df.drop_duplicates()

        table_df.to_csv(output_csv_file, index=False)
        adf.add_website_data(table_df, "Chemistry", case_names[0], plot_type="Tables")

        #Notify user that script has ended:
        print("  ...AMWG chemistry variable table has been generated successfully.")
    #End if chem table csv exists
    #End chemistry tables
    #--------------------




    """# Aerosol tables
    #-----------------
    #Notify user that script has started:
    print("\n  Calculating AMWG aerosol variable table...")

    aerosols_ext = {'_BURDEN':'_BURDEN','_CHMP':'_CHEM_PROD','_SF':'_EMIS',
                     '_DDF':'_DRYDEP','_WDF':'_WETDEP','_LIFETIME':'_LIFETIME'}

    #Create output file name:
    output_csv_file = output_location / f"amwg_aerosol_table_{case_names[0]}.csv"

    if output_csv_file.is_file():
        print(f"'{output_csv_file}' already exists, so skipping partner!\n")
        table_df = pd.read_csv(output_csv_file)
        adf.add_website_data(table_df, "Aerosols", case_names[0], plot_type="Tables")
        
    else:
        dic_SE = create_dic_SE(AEROSOLS,ListVars,ext1_SE)

        # extract all the data
        var_dict={}

        # this is for finding tropospheric values
        Dic_crit={}

        for i,scn in enumerate(scenarios):
            
            area=areas[scn]
            current_lat=Lats[scn]
            current_lon=Lons[scn]

            if regional:
                inside=Inside_SE(current_lat,current_lon,limit)
            else:
                if len(np.shape(area)) == 1:
                    inside=np.full((len(current_lon)),True)
                else:
                    inside=np.full((len(current_lat),len(current_lon)),True)

            current_dir=data_dirs[i]
            current_files=Files[scn] 

            var_dict[scn]={}
            Dic_var_comp={}

            for _,current_var in enumerate(AEROSOLS):
                # Components are: burden, chemical loss, chemical prod, dry deposition,
                #                 surface emissions, elevated emissions, wet deposition, gas-aerosol exchange

                if current_var=='SULF':
                    # For SULF we also have AQS and NUCLEATION
                    components=[current_var+'_BURDEN',current_var+'_CHML',current_var+'_CHMP',
                                current_var+'_DDF',current_var+'_WDF', current_var+'_SF', current_var+'_CLXF',
                                current_var+'_GAEX',current_var+'_DDFC',current_var+'_WDFC',current_var+'_AQS',
                                current_var+'_NUCL']
                else:
                    components=[current_var+'_BURDEN',current_var+'_CHML',current_var+'_CHMP',
                                current_var+'_DDF',current_var+'_WDF', current_var+'_SF', current_var+'_CLXF',
                                current_var+'_GAEX',current_var+'_DDFC',current_var+'_WDFC']

                Dic_comp={}
                for comp in components:
                    current_data=SEbudget(dic_SE,current_dir,current_files,comp,level=50)
                        
                    Dic_comp[comp]=current_data
                Dic_var_comp[current_var]=Dic_comp
            var_dict[scn]= Dic_var_comp

            #Critical threshholds????
            current_crit=SEbudget(dic_SE,current_dir,current_files,'O3',level=50)
            Dic_crit[scn]=current_crit

            if Tropospheric:
                trop=np.where(current_crit>150,np.nan,current_crit)
                strat=np.where(current_crit>150,current_crit,np.nan)
            else:
                trop=current_crit

        for current_var in AEROSOLS:
            for key,ext in aerosols_ext.items():
                row_values = []

                for i,scn in enumerate(scenarios):
                    my_val = calc_aerosol_data(scn,current_var,var_dict,trop,
                                                    area,durations[i],inside)[key]
                
                    if ext == "_BURDEN":
                        if current_var == "SULF":
                            new_ext = ext+" (TgS)"
                        else:
                            new_ext = ext+" (TgC)"
                    elif ext == "_LIFETIME":
                        if i == 0:
                            if 0 < my_val < 1:
                                my_val = my_val*365
                                new_ext = ext+" (days)"
                                
                            elif my_val > 1:
                                new_ext = ext+" (yr)"
                            elif int(my_val) == 0:
                                new_ext = ext+" (days)"
                        if i > 0:
                            if 0 < my_val < 1:
                                my_val = my_val*365
                                new_ext = ext+" (days)"
                    else:
                        if current_var == "SULF":
                            new_ext = ext+" (TgS/yr)"
                        else:
                            new_ext = ext+" (TgC/yr)"

                    row_values.append(np.round(my_val,3))
                row_values = [current_var+new_ext]+row_values
                
                dfentries = {c:[row_values[idx]] for idx,c in enumerate(cols)}
                # Add entries to Pandas structure:
                df = pd.DataFrame(dfentries,columns=cols)
                if output_csv_file.is_file():
                    df.to_csv(output_csv_file, mode='a', header=False, index=False)
                else:
                    df.to_csv(output_csv_file, header=False, index=False)
                #End scenarios
            #End keys()        
            
            # Add aqueous calc for SO4 only
            if current_var == "SULF":

                for key,ext in {'_AQS':'_AQ_PROD',}.items():
                    row_values = []

                    for i,scn in enumerate(scenarios):
                        my_val = calc_aerosol_data(scn,current_var,var_dict,trop,
                                                    area,durations[i],inside)[key]
                        row_values.append(np.round(my_val,3))
                    new_ext = ext+" (TgS/yr)"
                    row_values = [current_var+new_ext]+row_values
                
                    dfentries = {c:[row_values[idx]] for idx,c in enumerate(cols)}
                    # Add entries to Pandas structure:
                    df = pd.DataFrame(dfentries,columns=cols)
                    if output_csv_file.is_file():
                        df.to_csv(output_csv_file, mode='a', header=False, index=False)
                    else:
                        df.to_csv(output_csv_file, header=False, index=False)   
            # End if - SULF

        table_df = pd.read_csv(output_csv_file,names=cols)

        table_df = table_df.replace('SULF','SO4', regex=True)

        table_df.to_csv(output_csv_file, index=False)
        adf.add_website_data(table_df, "Aerosols", case_names[0], plot_type="Tables")

        #Notify user that script has ended:
        print("  ...AMWG aerosol variable table has been generated successfully.")
    #End if aerosol table exists"""

##################
# Helper functions
##################

def _load_data(dataloc, varname):
    ds = xr.open_dataset(dataloc)
    return ds[varname]

#####

def list_files(directory,scenario,start_date,end_date):

    """
        This function extracts the files in the directory that are within the chosen dates.
        
        directory: string showing he directory to look for files. always end with '/'
        
        scenario: string used to find which kind of files to find. 
                  Use below example to determine the scenario:
                      filename: 'f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_NudgeOutside_new.cam.h0.2017-11.nc'
                      scenario: 'f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_NudgeOutside_new.cam.h0'
        
        start_period: datetime foramt showing the starting date. ex: datetime.datetime(2017, 11, 1, 0, 0)
        
        end_period: datetime format showing the ending date . ex: datetime.datetime(2017, 12, 1, 0, 0)
    """
 
    # Get all the files within the directory
    # Can be a spot to liit to only the desired variables. Right now it is taking all available variables
    # from files. 
    #           *** Flag for possible upgrade/update ***
    #
    #for i in sorted(Path(directory).glob(f'*')):
    #    print(i)
    print("start_date[0:4]",start_date[0:4],"\n")
    start_filenames = sorted(Path(directory).glob(f'*.{start_date[0:4]}-*'))
    #print("start_filenames",start_filenames)
    all_start_filenames = [i.stem+".nc" for i in start_filenames]
    #print("all_start_filenames",all_start_filenames)

    end_filenames = sorted(Path(directory).glob(f'*.{end_date[0:4]}-*'))
    #print("end_filenames",end_filenames)
    all_end_filenames = [i.stem+".nc" for i in end_filenames]
    #print("all_end_filenames",all_end_filenames)
    
    all_filenames = sorted(all_start_filenames+all_end_filenames)

    if len(all_filenames)==0 : sys.exit(" Directory has no outputs ")
    if len(all_filenames)==0:
        return
    #all_filenames.sort()
    
    # this is used to discern what files to extract
    scenario_len=len(scenario)
    all_fileNames=[]

    for i in range(len(all_filenames)):
        if all_filenames[i][0:scenario_len]==scenario: # check if the file is relevant
            tmp_file=xr.open_dataset(directory+all_filenames[i])    
            # the times on filenames may not represent the exact time but time_bnds always does
            dim_time=tmp_file.dims['time']
            time_bounds=tmp_file['time_bnds'].data

            time_bounds0=time_bounds[0,0]
            time_bounds1=time_bounds[0,1]
           
            if dim_time==1:
                if time_bounds0==time_bounds1:
                    continue # initial file
            
            #I think the ADF will take care of the below time bounds check...
            # We need to only extract the files that are within the chosen dates.
            # first timestep is used for this purpose
            
            # For CAM data
            """
            filetime0=np.datetime64(time_bounds[0,0]) # beginning time of first timestep
            filetime1=np.datetime64(time_bounds[0,1]) # ending time of first timestep
            
            start_period = datetime.strptime(start_date, "%Y-%m-%d")
            end_period = datetime.strptime(end_date, "%Y-%m-%d")
            
            if '.h0' in scenario: # this is hard coded. User should change it (e.g. to ".h1") accordingly to reflect monthly files.
                if  (start_period<=filetime0<end_period) :
                    print ('list_files_SE Warning: "h0" is hard-coded to contain monthly files. If not, change it in the function.') 
                    all_fileNames.append(all_filenames[i])
                
            else:
                if (start_period<=filetime0<end_period) or (start_period<=filetime1<end_period):
                    all_fileNames.append(all_filenames[i])
            """

            if '.h0' in scenario:
                all_fileNames.append(all_filenames[i])
        
    return all_fileNames

#####

def Get_files(data_dirs, scenarios, start_periods, end_periods, **kwargs):
        
    """
        This function retrieves the files, latitude, and longitude information
        in all the directories within the chosen dates.
        
        data_dirs: list showing the directories to look for files. always end with '/'
        
        scenarios: list the scenarios used to find which kind of files to find. 
                  Use below example to determine each scenario:
                      filename: 'f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_NudgeOutside_new.cam.h0.2017-11.nc'
                      scenario: 'f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_NudgeOutside_new.cam.h0'
        
        start_date: string showing the starting date. ex: "$Y-$M-$D"
        
        end_period: string showing the ending date. ex: "$Y-$M-$D"
        
        ext1_SE: string specifying if the files are for only a region, which changes to variable names.
                ex: if you saved files for a only a box region ($LL_lat$,$LL_lon$,$UR_lat$,$UR_lon$),
                    the 'lat' variable will be saved as: 'lat_$LL_lon$e_to_$UR_lon$e_$LL_lat$n_to_$UR_lat$n'
                    for instance: 'lat_65e_to_91e_20n_to_32n'
    """
    ext1_SE=kwargs.pop('ext1_SE','')
    Area=kwargs.pop('area',False)

    files={}
    Lats={}
    Lons={}

    areas={}
    Earth_rad=6.371e6 # Earth Radius in m 

    
    for i,scn in enumerate(scenarios):

        current_dir=data_dirs[i]

        # find the needed the files
        print("start_periods[i],end_periods[i]",start_periods[i],end_periods[i],"\n")
        current_files=list_files(current_dir,scn,start_periods[i],end_periods[i])
        # get the Lat and Lons for each scenario
        tmp_file=xr.open_dataset(current_dir+current_files[i])
        lon=tmp_file['lon'+ext1_SE].data
        lon[lon > 180.] -= 360 # shift longitude from 0-360˚ to -180-180˚
        lat=tmp_file['lat'+ext1_SE].data

        if Area==True:
            try:
                tmp_area=tmp_file['area'+ext1_SE].data
                Earth_area= 4 * np.pi * Earth_rad**(2)
                areas[scn]=tmp_area*Earth_area/np.nansum(tmp_area)

            except KeyError:
                dlon= np.abs(lon[1]-lon[0])
                dlat= np.abs(lat[1]-lat[0])

                lon2d,lat2d=np.meshgrid(lon,lat)
                #area=np.zeros_like(lat2d)

                dy=Earth_rad*dlat*np.pi/180
                dx=Earth_rad*np.cos(lat2d*np.pi/180)*dlon*np.pi/180

                area=dx*dy
                areas[scn]=area

        files[scn]=current_files
        Lats[scn]=lat
        Lons[scn]=lon
        areas[scn]=area

    print("Got the files, lats, lons, and areas (hopefully)...")
    return files, Lats, Lons, areas

#####

def SEbudget(dic_SE,data_dir,files,var,**kwargs):     

    """
        This function is used for getting the data for the budget calculation.

        
        dic_SE: dictionary specyfing what variables to get. For example, 
                for precipitation you can define SE as:
                    dic_SE['PRECT']={'PRECC'+ext1_SE:8.64e7,'PRECL'+ext1_SE:8.64e7}
                    - It means to sum the file variables "PRECC" and "PRECL" 
                      for my arbitrary desired variable named "PRECT"
                      
                    - It also has the option to apply conversion factors. 
                      For instance, PRECL and PRECC are in m/s. 8.64e7 is used to convernt m/s to mm/day


        data_dir: string of the directory that contains the files. always end with '/'
        
        files: list of the files to be read 
        
        var: string showing the variable to be extracted.
 
    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # This currently is gathering and storing data with numpy arrays
    # Probably makes sense to transistion this to xarray data arrays for seasonal weighting, right?
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # gas constanct
    Rgas=287.04 #[J/K/Kg]=8.314/0.028965
        
    all_data=[]
    for file in files:
        ds=xr.open_dataset(data_dir+file)
        data=[]
    
        for i in dic_SE[var].keys():
            #Check to see if the product is in the actual dataset, if not, move on and set to 0
            if i in ds:
                data.append(np.array(ds[i].isel(time=0))*dic_SE[var][i])
            else:
                print(f"Looks like {var} is missing {i} for {file}, so skipping...\n")
        data=np.sum(data,axis=0)
            
        if ('CHML' in var) or ('CHMP' in var) : 

            if 'T' in ds:
                Temp=np.array(ds['T'].isel(time=0))
            else:
                Temp=0
            if 'PMID' in ds:
                Pres=np.array(ds['PMID'].isel(time=0))
            else:
                Pres=0
            rho= Pres/(Rgas*Temp)
            
            if 'PDELDRY' in ds:
                delP=np.array(ds['PDELDRY'].isel(time=0))
            else:
                delP=0
            data=data*delP/rho
        elif ('BURDEN' in var):

            if 'PDELDRY' in ds:
                delP=np.array(ds['PDELDRY'].isel(time=0))
            else:
                delP=0 
            data=data*delP
        else:
            data=data
        #End if - vars
        all_data.append(data)
    # End for - files
    
    #Flush out the nans and take mean
    all_data=np.nanmean(all_data,axis=0)
    
    return all_data 

#####

def create_dic_SE(variables, ListVars, ext1_SE):
    #Dictionary for Molecular weights. Keys must be consistent with variable name
    MW={'O3':48,
        'CH4':16,
        'CO':28,
        'ISOP':68,
        'C10H16':136,
        'SOA':144.132,
        'SALT':12.011,
        'SULF':115.11,
        'POM':12.011,
        'BC':12.011 ,
        'DUST':12.011,
        'CH3CCL3':133.4042,
        'C10H16':136.2340,
        'CH3OH':32.0419,
        'CH3COCH3':58.0791,
        'AIR':28.97}

    # Avogadro's Number
    avo = 6.022e23
    # gravity - replace by numpy or scipy for best guess of gravity
    gr = 9.80616
    
    #dic_SE = set_dic_SE(ext1_SE)

    dic_SE={}

    dic_SE['O3']={'O3'+ext1_SE:1e9} # covert to ppb for Tropopause calculation
    dic_SE['CH4']={'CH4'+ext1_SE:1}
    dic_SE['CO']={'CO'+ext1_SE:1}

    dic_SE['ISOP']={'ISOP'+ext1_SE:1}
    dic_SE['C10H16']={'MTERP'+ext1_SE:1}
    dic_SE['CH3OH']={'CH3OH'+ext1_SE:1}
    dic_SE['CH3COCH3']={'CH3COCH3'+ext1_SE:1}
    dic_SE['CH3CCL3']={'CH3CCL3'+ext1_SE:1}


    dic_SE['SOA']={'soa1_a1'+ext1_SE:1,
                  'soa2_a1'+ext1_SE:1,
                  'soa3_a1'+ext1_SE:1,
                  'soa4_a1'+ext1_SE:1,
                  'soa5_a1'+ext1_SE:1,
                  'soa1_a2'+ext1_SE:1,
                  'soa2_a2'+ext1_SE:1,
                  'soa3_a2'+ext1_SE:1,
                  'soa4_a2'+ext1_SE:1,
                  'soa5_a2'+ext1_SE:1}

    dic_SE['DUST']={'dst_a1'+ext1_SE:1,
                  'dst_a2'+ext1_SE:1,
                  'dst_a3'+ext1_SE:1}

    dic_SE['SALT']={'ncl_a1'+ext1_SE:1,
                  'ncl_a2'+ext1_SE:1,
                  'ncl_a3'+ext1_SE:1}

    dic_SE['POM']={'pom_a1'+ext1_SE:1,
                  'pom_a4'+ext1_SE:1}

    dic_SE['BC']={'bc_a1'+ext1_SE:1,
                  'bc_a4'+ext1_SE:1}

    dic_SE['SULF']={'so4_a1'+ext1_SE:1,
                  'so4_a2'+ext1_SE:1,
                  'so4_a3'+ext1_SE:1}
    
    
    # Here, we deal with conversion factors for different components
    # Some conversion factors need density or Layer's pressure, that will be
    # accounted for when reading the files. 
    # We convert everying to kg/m2/s or kg/m2 or kg/s, so that final Tg/yr or Tg results are consistent
    for var in variables:

        dic_SE[var+'_BURDEN']={}
        dic_SE[var+'_CHML']={}    
        dic_SE[var+'_CHMP']={}    

        dic_SE[var+'_SF']={}    
        dic_SE[var+'_CLXF']={}    

        dic_SE[var+'_DDF']={}     
        dic_SE[var+'_WDF']={}

        if var in AEROSOLS:
            dic_SE[var+'_GAEX']={}       
            dic_SE[var+'_DDFC']={}
            dic_SE[var+'_WDFC']={}
        else:
            dic_SE[var+'_TEND']={}               
            dic_SE[var+'_LNO']={}               

        # We have nucleation and aqueous chemistry for sulfate.
        if var=='SULF':
            dic_SE[var+'_NUCL']={}
            dic_SE[var+'_AQS']={}

        var_keys=dic_SE[var].keys()

        for key in var_keys:

            # for CHML and CHMP:
            # original unit : [molec/cm3/s]
            # following Tilmes code to convert to [kg/m2/s]
            # conversion: Mw*rho*delP*1e3/Avo/gr
            # rho and delP will be applied when reading the files in SEbudget function.   

            if key=='O3'+ext1_SE: 
                # for O3, we should not include fast cycling reactions
                # As a result, we use below diagnostics in the model
                dic_SE[var+'_CHML'][key+'_Loss'+ext1_SE]=MW[var]*1e3/avo/gr
                dic_SE[var+'_CHMP'][key+'_Prod'+ext1_SE]=MW[var]*1e3/avo/gr 

            else:
                if key+'_CHML' in ListVars:
                    dic_SE[var+'_CHML'][key+'_CHML'+ext1_SE]=MW[var]*1e3/avo/gr
                else:
                    dic_SE[var+'_CHML']['O3'+ext1_SE]=0.

                if key+'_CHMP' in ListVars:            
                    dic_SE[var+'_CHMP'][key+'_CHMP'+ext1_SE]=MW[var]*1e3/avo/gr        
                else:
                    dic_SE[var+'_CHMP']['O3'+ext1_SE]=0.        

            # for SF:
            # original unit: [kg/m2/s]        
            if 'SF'+key in ListVars:            
                dic_SE[var+'_SF']['SF'+key+ext1_SE]=1        
            else:
                dic_SE[var+'_SF']['SFCO'+ext1_SE]=0.   

            # for CLXF:
            # original unit: [molec/cm2/s]
            # conversion: Mw*10/Avo
            if key+'_CLXF' in ListVars:            
                dic_SE[var+'_CLXF'][key+'_CLXF'+ext1_SE]=MW[var]*10/avo  # convert [molec/cm2/s] to [kg/m2/s]        
            else:
                dic_SE[var+'_CLXF']['O3'+ext1_SE]=0. 

            if var in AEROSOLS:
                # for each species:
                # original unit : [kg/kg]  in dry air
                # convert to [kg/m2]
                # conversion: delP/gr     
                # delP will be applied when reading the files in SEbudget function. 
                if key in ListVars:
                    dic_SE[var+'_BURDEN'][key+ext1_SE]=1/gr
                else:
                    dic_SE[var+'_BURDEN']['O3'+ext1_SE]=0              

                # for DDF:
                # original unit: [kg/m2/s]
                if key+'DDF' in ListVars:            
                    dic_SE[var+'_DDF'][key+ext1_SE+'DDF']=1        
                else:
                    dic_SE[var+'_DDF']['DF_O3'+ext1_SE]=0.  
                    
                # for SFWET:
                # original unit: [kg/m2/s]
                if key+'SFWET' in ListVars:            
                    dic_SE[var+'_WDF'][key+ext1_SE+'SFWET']=1        
                else:
                    dic_SE[var+'_WDF']['DF_O3'+ext1_SE]=0.                  
                    #dic_SE[var+'_WDF']['O3'+ext1_SE]=0.                  

                # for sfgaex1:
                # original unit: [kg/m2/s]
                if key+'_sfgaex1' in ListVars:            
                    dic_SE[var+'_GAEX'][key+ext1_SE+'_sfgaex1']=1        
                else:
                    dic_SE[var+'_GAEX']['DF_O3'+ext1_SE]=0.                

                # for DDF in cloud water:
                # original unit: [kg/m2/s]
                cloud_key=key[:-2]+'c'+key[-1]
                if cloud_key+ext1_SE+'DDF' in ListVars:            
                    dic_SE[var+'_DDFC'][cloud_key+ext1_SE+'DDF']=1        
                else:
                    dic_SE[var+'_DDFC']['DF_O3'+ext1_SE]=0.  

                # for SFWET in cloud water:
                # original unit: [kg/m2/s]
                if cloud_key+ext1_SE+'SFWET' in ListVars:            
                    dic_SE[var+'_WDFC'][cloud_key+ext1_SE+'SFWET']=1        
                else:
                    dic_SE[var+'_WDFC']['DF_O3'+ext1_SE]=0.                  

                if var=='SULF':
                    # for Nucleation :
                    # original unit: [kg/m2/s]
                    if key+ext1_SE+'_sfnnuc1' in ListVars:            
                        dic_SE[var+'_NUCL'][key+ext1_SE+'_sfnnuc1']=1        
                    else:
                        dic_SE[var+'_NUCL']['DF_O3'+ext1_SE]=0.  
                        #dic_SE[var+'_NUCL']['O3'+ext1_SE]=0.  

                    # for Aqueous phase :
                    # original unit: [kg/m2/s]
                    if cloud_key+ext1_SE+'AQSO4' in ListVars:            
                        dic_SE[var+'_AQS'][cloud_key+ext1_SE+'AQSO4']=1        
                    else:
                        dic_SE[var+'_AQS']['DF_O3'+ext1_SE]=0.

                    if cloud_key+ext1_SE+'AQH2SO4' in ListVars:            
                        dic_SE[var+'_AQS'][cloud_key+ext1_SE+'AQH2SO4']=1        
                    else:
                        dic_SE[var+'_AQS']['DF_O3'+ext1_SE]=0.

            else:
                # for each species:
                # original unit : [mole/mole]  in dry air
                # convert to [kg/m2]
                # conversion: Mw*delP/Mwair/gr     Mwair=28.97 gr/mole
                # delP will be applied when reading the files in SEbudget function. 
                if key in ListVars:
                    dic_SE[var+'_BURDEN'][key+ext1_SE]=MW[var]/28.97/gr
                else:
                    dic_SE[var+'_BURDEN']['O3'+ext1_SE]=0            

                # for DF:
                # original unit: [kg/m2/s]
                if 'DF_'+key in ListVars:            
                    dic_SE[var+'_DDF']['DF_'+key+ext1_SE]=1        
                else:
                    dic_SE[var+'_DDF']['DF_O3'+ext1_SE]=0.              

                # for WD:
                # original unit: [kg/m2/s]
                if 'WD_'+key in ListVars:            
                    dic_SE[var+'_WDF']['WD_'+key+ext1_SE]=1        
                else:
                    dic_SE[var+'_WDF']['DF_O3'+ext1_SE]=0.                  

                # for Chem tendency:
                # original unit: [kg/s]
                # conversion: not needed
                if 'D'+key+'CHM' in ListVars:            
                    dic_SE[var+'_TEND']['D'+key+'CHM'+ext1_SE]=1  # convert [kg/s] to [kg/s]        
                else:
                    dic_SE[var+'_TEND']['DO3CHM'+ext1_SE]=0    

                # for Lightning NO production: (always in gas)
                # original unit: [Tg N/Yr]
                # conversion: not needed
                if 'LNO_COL_PROD' in ListVars:            
                    dic_SE[var+'_LNO']['LNO_COL_PROD'+ext1_SE]=1  # convert [Tg N/yr] to [Tg N /yr]        
                else:
                    dic_SE[var+'_LNO']['DF_O3'+ext1_SE]=0
    return dic_SE

#####

def calc_chem_data(scn, var, var_dict, trop, area, duration, inside):
    """
    Calcs for chem diagnostics (no aerosols)
    
    Meant to be a brute force idea of just adding code blocks for each new derived quantity
    Seems likely there is a better/proper way of doing this, will look into
    
        *** Flag for upgrade ^^^^^^ ***
    
     - Add to list each budget item
     - returns list of final budget values
    """
    
    chem_dict = {}
    
    # Burden      
    spc_burd=var_dict[scn][var][var+'_BURDEN']       
    spc_burd=np.where(np.isnan(trop),np.nan,spc_burd)
    tmp_burden=np.nansum(spc_burd*area,axis=0)
    burden=np.ma.masked_where(inside==False,tmp_burden)  #convert Kg/m2 to Tg
    BURDEN = np.ma.sum(burden*1e-9)
    chem_dict['_BURDEN'] = np.round(BURDEN,5)

    # Chemical Loss
    spc_chml=var_dict[scn][var][var+'_CHML'] 
    spc_chml=np.where(np.isnan(trop),np.nan,spc_chml)       
    tmp_chml=np.nansum(spc_chml*area,axis=0)
    chml=np.ma.masked_where(inside==False,tmp_chml)  #convert Kg/m2/s to Tg/yr
    CHML = np.ma.sum(chml*duration*1e-9)
    chem_dict['_CHML'] = np.round(CHML,5)

    # Chemical Production
    spc_chmp=var_dict[scn][var][var+'_CHMP'] 
    spc_chmp=np.where(np.isnan(trop),np.nan,spc_chmp)
    tmp_chmp=np.nansum(spc_chmp*area,axis=0)
    chmp=np.ma.masked_where(inside==False,tmp_chmp)  #convert Kg/m2/s to Tg/yr
    CHMP = np.ma.sum(chmp*duration*1e-9)
    chem_dict['_CHMP'] = np.round(CHMP,5)
        
    # Surface Emissions
    spc_sf=var_dict[scn][var][var+'_SF'] 
    tmp_sf=spc_sf
    sf=np.ma.masked_where(inside==False,tmp_sf*area)  #convert Kg/m2/s to Tg/yr
    SF = np.ma.sum(sf*duration*1e-9)
    chem_dict['_SF'] = np.round(SF,5)

    # Elevated Emissions
    if var == "CO":
        print(f"Smoethign is borken with {var}")
        CLXF = np.nan
    else:
        spc_clxf=var_dict[scn][var][var+'_CLXF'] 
        tmp_clxf=np.nansum(spc_clxf*area,axis=0)
        clxf=np.ma.masked_where(inside==False,tmp_clxf)  #convert Kg/m2/s to Tg/yr
        CLXF = np.ma.sum(clxf*duration*1e-9)
    chem_dict['_CLXF'] = np.round(CLXF,5)

    # Dry Deposition Flux 
    spc_ddf=var_dict[scn][var][var+'_DDF'] 
    tmp_ddf=spc_ddf
    ddf=np.ma.masked_where(inside==False,tmp_ddf*area)  #convert Kg/m2/s to Tg/yr
    DDF = np.ma.sum(ddf*duration*1e-9)
    chem_dict['_DDF'] = np.round(DDF,5)
            
    # Wet Deposition Flux   
    spc_wdf=var_dict[scn][var][var+'_WDF'] 
    tmp_wdf=spc_wdf
    wdf=np.ma.masked_where(inside==False,tmp_wdf*area)  #convert Kg/m2/s to Tg/yr
    WDF = np.ma.sum(wdf*duration*1e-9)
    chem_dict['_WDF'] = np.round(WDF,5)
             
    # Chemical Tendency
    spc_tnd=var_dict[scn][var][var+'_TEND'] 
    spc_tnd=np.where(np.isnan(trop),np.nan,spc_tnd)
    tmp_tnd=np.nansum(spc_tnd,axis=0)
    tnd=np.ma.masked_where(inside==False,tmp_tnd)  #convert Kg/s to Tg/yr
    TND = np.ma.sum(tnd*duration*1e-9)
    chem_dict['_TEND'] = np.round(TND,5)
    
    # Stratospheric-Tropospheric Exchange
    STE=DDF-TND
    chem_dict['_STE'] = np.round(STE,5)

    # LifeTime = Burden/(loss+deposition)
    LT=BURDEN/(CHML+DDF-WDF)*duration/86400/365 # days
    chem_dict['_LIFETIME'] = np.round(LT,5)

    # Lightning NOX production
    spc_lno=var_dict[scn][var][var+'_LNO']
    tmp_lno=np.ma.masked_where(inside==False,spc_lno)  
    LNO = np.ma.sum(tmp_lno)
    chem_dict['_LNO'] = np.round(LNO,5)
    
    NET = CHMP-CHML
    chem_dict['_NET'] = np.round(NET,5)
    
    return chem_dict

#####

def calc_aerosol_data(scn, var, var_dict, trop, area, duration, inside):
    
    aerosol_dict = {}
    
    # Burden      
    spc_burd=var_dict[scn][var][var+'_BURDEN']       
    spc_burd=np.where(np.isnan(trop),np.nan,spc_burd)
    tmp_burden=np.nansum(spc_burd*area,axis=0)
    burden=np.ma.masked_where(inside==False,tmp_burden)  #convert Kg/m2 to Tg
    BURDEN = np.ma.sum(burden*1e-9)
    aerosol_dict['_BURDEN'] = np.round(BURDEN,5)

    # Chemical Loss
    spc_chml=var_dict[scn][var][var+'_CHML'] 
    spc_chml=np.where(np.isnan(trop),np.nan,spc_chml)       
    tmp_chml=np.nansum(spc_chml*area,axis=0)
    chml=np.ma.masked_where(inside==False,tmp_chml)  #convert Kg/m2/s to Tg/yr
    CHML = np.ma.sum(chml*duration*1e-9)
    aerosol_dict['_CHML'] = np.round(CHML,5)
    
    # Chemical Production
    spc_chmp=var_dict[scn][var][var+'_CHMP'] 
    spc_chmp=np.where(np.isnan(trop),np.nan,spc_chmp)
    tmp_chmp=np.nansum(spc_chmp*area,axis=0)
    chmp=np.ma.masked_where(inside==False,tmp_chmp)  #convert Kg/m2/s to Tg/yr
    CHMP = np.ma.sum(chmp*duration*1e-9)
    aerosol_dict['_CHMP'] = np.round(CHMP,5)
        
    # Surface Emissions
    spc_sf=var_dict[scn][var][var+'_SF'] 
    tmp_sf=spc_sf
    sf=np.ma.masked_where(inside==False,tmp_sf*area)  #convert Kg/m2/s to Tg/yr
    SF = np.ma.sum(sf*duration*1e-9)
    aerosol_dict['_SF'] = np.round(SF,5)
 
    # Elevated Emissions
    spc_clxf=var_dict[scn][var][var+'_CLXF'] 
    tmp_clxf=np.nansum(spc_clxf*area,axis=0)
    clxf=np.ma.masked_where(inside==False,tmp_clxf)  #convert Kg/m2/s to Tg/yr
    CLXF = np.ma.sum(clxf*duration*1e-9)
    aerosol_dict['_CLXF'] = np.round(CLXF,5)
    
    # Dry Deposition Flux      
    spc_ddfa=var_dict[scn][var][var+'_DDF'] 
    spc_ddfc=var_dict[scn][var][var+'_DDFC']
    spc_ddf=spc_ddfa +spc_ddfc
    tmp_ddf=spc_ddf
    ddf=np.ma.masked_where(inside==False,tmp_ddf*area)  #convert Kg/m2/s to Tg/yr
    DDF = np.ma.sum(ddf*duration*1e-9)
    aerosol_dict['_DDF'] = np.round(DDF,5)
            
    # Wet deposition
    spc_wdfa=var_dict[scn][var][var+'_WDF'] 
    spc_wdfc=var_dict[scn][var][var+'_WDFC']
    spc_wdf=spc_wdfa +spc_wdfc            
    tmp_wdf=spc_wdf
    wdf=np.ma.masked_where(inside==False,tmp_wdf*area)  #convert Kg/m2/s to Tg/yr
    WDF = np.ma.sum(wdf*duration*1e-9)
    aerosol_dict['_WDF'] = np.round(WDF,5)
            
    # gas-aerosol Exchange
    spc_gaex=var_dict[scn][var][var+'_GAEX'] 
    tmp_gaex=spc_gaex
    gaex=np.ma.masked_where(inside==False,tmp_gaex*area)  #convert Kg/m2/s to Tg/yr
    GAEX = np.ma.sum(gaex*duration*1e-9)
    aerosol_dict['_GAEX'] = np.round(GAEX,5)      
            
    # LifeTime = Burden/(loss+deposition)
    LT=BURDEN/(CHML+DDF-WDF)* duration/86400 # days   
    aerosol_dict['_LIFETIME'] = np.round(LT,5)
            
    if var=='SULF':     
        # Aqueous Chemistry
        spc_aqs=var_dict[scn][var][var+'_AQS'] 
        tmp_aqs=spc_aqs
        aqs=np.ma.masked_where(inside==False,tmp_aqs*area)  #convert Kg/m2/s to Tg/yr
        AQS = np.ma.sum(aqs*duration*1e-9)
        aerosol_dict['_AQS'] = np.round(AQS,5)    
        
        # Nucleation
        spc_nucl=var_dict[scn][var][var+'_NUCL'] 
        tmp_nucl=spc_nucl
        nucl=np.ma.masked_where(inside==False,tmp_nucl*area)  #convert Kg/m2/s to Tg/yr
        NUCL = np.ma.sum(nucl*duration*1e-9)
        aerosol_dict['_NUCL'] = np.round(NUCL,5)

    return aerosol_dict

#####


##############
#END OF SCRIPT
