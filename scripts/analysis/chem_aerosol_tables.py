import numpy as np
import xarray as xr
import sys
from pathlib import Path
import warnings  # use to warn user about missing files.

from datetime import datetime, timedelta
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt # Core library for plotting
import numpy as np

from cftime import DatetimeNoLeap
from netCDF4 import Dataset

from glob import glob
import pickle
import json

import itertools

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

# Import necessary ADF modules:
from adf_diag import AdfDiag
from adf_base import AdfError


def chem_aerosol_tables(adfobj):


    # INPUTS
    #list of the gaseous variables to be caculated.
    GAS_VARIABLES = ["CH4",'CH3COCH3','CO']

    # list of the aerosol variables to be caculated.
    #AEROSOL_VARIABLES=['AOD','DAOD','SOA', 'SALT', 'DUST', 'POM', 'BC', 'SO4']
    AEROSOL_VARIABLES = ['AOD','SO4']

    #list of all the variables to be caculated.
    VARIABLES = GAS_VARIABLES + AEROSOL_VARIABLES

    # For the case that outputs are saved for a specific region.
    # i.e., when using fincllonlat in user_nl_cam
    ext1_SE = ''

    # Tropospheric Values
    # -------------------
    # if True, calculate only Tropospheric values
    # if False, all layers
    # tropopause is defiend as o3>150ppb. If needed, change accordingly.
    Tropospheric = True


    ### NOT WORKING FOR NOW
    # To calculate the budgets only for a region
    # Lat/Lon extent
    limit = (20,20,40,120)
    regional = False


    #Dictionary for Molecular weights. Keys must be consistent with variable name
    # For aerosols, the MW is used only for chemical loss, chemical production, and elevated emission calculations
    # For SO4, we report everything in terms of Sulfur, so we use Sulfur MW here
    MW={'O3':48,
        'CH4':16,
        'CO':28,
        'ISOP':68,
        'MTERP':136,
        'SOA':144.132,
        'SALT':58.4412,
        'SO4':32.066,
        'POM':12.011,
        'BC':12.011 ,
        'DUST':168.0456,
        'CH3CCL3':133.4042,
        'CH3OH':32,
        'CH3COCH3':58}

    # Avogadro's Number
    AVO=6.022e23
    # gravity
    gr=9.80616
    # Mw air
    Mwair=28.97


    # The variables in the list below must be aerosols - do not add AOD and DAOD
    # no need to change this list, unless for a specific need!
    AEROSOLS = ['SOA', 'SALT', 'DUST', 'POM', 'BC', 'SO4']




    #CAM simulation variables (these quantities are always lists):
    #case_names = adfobj.get_cam_info('cam_case_name', required=True)
    case_names = adfobj.get_cam_info('cam_case_name', required=True) + [adfobj.get_baseline_info("cam_case_name")]
    case_names_len = len(case_names)

    #Grab all case nickname(s)
    test_nicknames_list = adfobj.case_nicknames["test_nicknames"]
    base_nickname_list = adfobj.case_nicknames["base_nickname"]
    nicknames_list = test_nicknames_list + [base_nickname_list]
    nicknames = {}


    res = adfobj.variable_defaults # dict of variable-specific plot preferences

    start_yrs = adfobj.climo_yrs["syears"] + [adfobj.climo_yrs["syear_baseline"]]
    end_yrs = adfobj.climo_yrs["eyears"] + [adfobj.climo_yrs["eyear_baseline"]]

    #Histort file num
    #h_case = "h0a"
    #h_case = "h0"

    #Grab history strings:
    cam_hist_strs = adfobj.hist_string["test_hist_str"]
    hist_strs = cam_hist_strs + [adfobj.hist_string["base_hist_str"]]

    # Filter the list to include only strings that are possible h0 strings
    # - Search for either h0 or h0a
    substrings = {"cam.h0","cam.h0a"}
    case_hist_strs = []
    for cam_case_str in hist_strs:
        # Check each possible h0 string
        for string in cam_case_str:
            if string in substrings:
                case_hist_strs.append(string)
                break

    #Grab history file locations from config yaml file
    cam_hist_locs = adfobj.get_cam_info("cam_hist_loc", required=True)
    #cam_hist_locs = [adfobj.get_baseline_info("cam_hist_loc", required=True)]
    hist_locs = cam_hist_locs + [adfobj.get_baseline_info("cam_hist_loc")]

    #Create path object for the CAM history file(s) location:
    data_dirs = []
    for case_idx,case in enumerate(case_names):
        nicknames[case] = nicknames_list[case_idx]
        #Check that history file input directory actually exists:
        if not Path(hist_locs[case_idx]).is_dir():
            errmsg = f"History files directory '{hist_locs[case_idx]}' not found.  Script is exiting."
            raise AdfError(errmsg)

        #Write to debug log if enabled:
        adfobj.debug_log(f"DEBUG: location of files is {str(hist_locs[case_idx])}")

        data_dirs.append(hist_locs[case_idx])

    #End gathering case, path, and data info
    #-----------------------------------------

    print("\nhist_strs",hist_strs,"\n")
    # Look for specific h-case
    #scenarios = [f'{ix}.cam.{h_case}' for ix in case_names]
    scenarios = [f'{case}.{hist_strs[ix][0]}' for ix,case in enumerate(case_names)]

    print("scenarios info:",len(scenarios),scenarios)


    # Periods of Interest
    # -------------------
    # choose the period of interest. Plots will be averaged within this period
    #start_dates = [f"{start_yrs[0]}-1-1", f"{start_yrs[1]}-1-1"]
    #end_dates = [f"{end_yrs[0]}-1-1", f"{end_yrs[1]}-1-1"]

    start_dates = []
    for syr in start_yrs:
        start_dates.append(f"{syr}-1-1")
    end_dates = []
    for eyr in end_yrs:
        end_dates.append(f"{eyr}-1-1")

    #start_periods = []
    #end_periods = []
    durations = {}
    num_yrs = {}
    """for i,val in enumerate(start_dates):
    # convert date strings to datetime format
        start_period = datetime.strptime(start_dates[i], "%Y-%m-%d")
        end_period = datetime.strptime(end_dates[i], "%Y-%m-%d")

        #start_periods.append(start_period)
        #end_periods.append(end_period)

        #durations.append((end_period-start_period).days*86400+365*86400)
        durations[case_names[i]] = (end_period-start_period).days*86400+365*86400

        # Get number of years for calculations
        #num_yrs.append(int(end_dates[i])-int(start_dates[i])+1)
        num_yrs[case_names[i]] = int(end_yrs[i])-int(start_yrs[i])+1
        print(f"number of years: {int(end_yrs[i])-int(start_yrs[i])+1}")
    """





    # Main function
    #--------------

    # Set dictionary of components for each case
    Dic_scn_var_comp = {}
    areas = {}
    trops = {}
    insides = {}
    for i,case in enumerate(case_names):
        print(f'Current Scenario: {case}',"\n",len(f'Current Scenario: {case}')*'-','\n')
        #print(len(f'Current Scenario: {case}')*'-','\n')

        #print(f'Current Scenario: {case}\n{len(f'Current Scenario: {case}')*'-','\n')

        start_period = datetime.strptime(start_dates[i], "%Y-%m-%d")
        end_period = datetime.strptime(end_dates[i], "%Y-%m-%d")

        #start_periods.append(start_period)
        #end_periods.append(end_period)

        #durations.append((end_period-start_period).days*86400+365*86400)
        durations[case_names[i]] = (end_period-start_period).days*86400+365*86400

        # Get number of years for calculations
        #num_yrs.append(int(end_dates[i])-int(start_dates[i])+1)
        num_yrs[case_names[i]] = int(end_yrs[i])-int(start_yrs[i])+1
        print(f"number of years: {int(end_yrs[i])-int(start_yrs[i])+1}")

        Files,Lats,Lons,areas[case],ext1_SE = Get_files(data_dirs[i],start_yrs[i],end_yrs[i],case_hist_strs[i],area=True)

        # find the name of all the variables in the file.
        # this will help the code to work for the variables that are not in the files (assingn 0s)
        tmp_file = Dataset(Path(data_dirs[i]) / Files[0])
        ListVars = tmp_file.variables.keys()
        tmp_file.close()

        # Set up and fill dictionaries for components for current cases
        dic_SE = set_dic_SE(ListVars,ext1_SE)
        #dic_SE = fill_dic_SE(dic_SE,VARIABLES,ListVars)
        dic_SE = fill_dic_SE(dic_SE, VARIABLES, ListVars, ext1_SE, AEROSOLS, MW, AVO, gr, Mwair)

        try:
            with open(f'{case}_Dic_scn_var_comp.pickle', 'rb') as handle:
                Dic_scn_var_comp[case] = pickle.load(handle)

            with open(f'{case}_Dic_crit.pickle', 'rb') as handle2:
                Dic_crit = pickle.load(handle2)
        except:
            print("JSON file not found, need to create the files.")

            # Make dictionary of all data for each case
            print(f"\t Calculating values for {case}")
            #Dic_crit, Dic_scn_var_comp[case] = make_Dic_scn_var_comp(VARIABLES, data_dirs[i], dic_SE)
            Dic_crit, Dic_scn_var_comp[case] = make_Dic_scn_var_comp(adfobj, VARIABLES, data_dirs[i], dic_SE, Files, ext1_SE, AEROSOLS)

            with open(f'{case}_Dic_scn_var_comp.pickle', 'wb') as handle:
                pickle.dump(Dic_scn_var_comp[case], handle, protocol=pickle.HIGHEST_PROTOCOL)

            with open(f'{case}_Dic_crit.pickle', 'wb') as handle:
                pickle.dump(Dic_crit, handle, protocol=pickle.HIGHEST_PROTOCOL)

        if regional:
            #inside=Inside_SE_region(current_lat,current_lon,dir_shapefile)
            inside = Inside_SE(Lats,Lons,limit)
        else:
            if len(np.shape(areas)) == 1:
                inside = np.full((len(Lons)),True)
            else:
                inside = np.full((len(Lats),len(Lons)),True)

        current_crit = Dic_crit[0]

        if Tropospheric:
            trop = np.where(current_crit>150,np.nan,current_crit)
            #strat=np.where(current_crit>150,current_crit,np.nan)
        else:
            trop=current_crit
        trops[case] = trop
        insides[case] = inside


    if len(AEROSOL_VARIABLES) > 0:
        print("\tMaking table for aerosols")
        #aerosol_table = make_table(adfobj, AEROSOL_VARIABLES,'aerosols',Dic_scn_var_comp,areas,trops,case_names, durations, insides, AEROSOLS)
        aerosol_table = make_table(adfobj, AEROSOL_VARIABLES, 'aerosols', Dic_scn_var_comp, areas, trops, case_names, nicknames, durations, insides, num_yrs, AEROSOLS)
        #make_table(vars, chem_type, Dic_scn_var_comp, areas, trops, case_names, durations, insides, AEROSOLS)

    if len(GAS_VARIABLES) > 0:
        print("\tMaking table for gases")
        #gas_table = make_table(adfobj, GAS_VARIABLES,'gases',Dic_scn_var_comp,areas,trops,case_names, durations, insides, AEROSOLS)
        gas_table = make_table(adfobj, GAS_VARIABLES, 'gases', Dic_scn_var_comp, areas, trops, case_names, nicknames, durations, insides, num_yrs, AEROSOLS)
        #make_table(vars, chem_type, Dic_scn_var_comp, areas, trops, case_names, durations, insides, AEROSOLS)

    #return


##################
# Helper functions
##################



'''
SE_functions.py
this code is designed for compiling the functions used for processing SE files


MODIFICATION HISTORY:
    Behrooz Roozitalab, 02, NOV, 2022: VERSION 1.00
    - Initial version

    Justin Richling, 27 Nov, 2023
    - updated to fit to ADF and check with old AMWG chem/aerosol tables
    - fixed:
        * lifetime inconsitencies
        * added difference bewtween cases column to tables

    Behrooz Roozitalab, 8 Aug, 2024
    - fixed:
        * lifetime inconsitencies
        * Removed redundant calculations to improve the speed
        * Verified the results against the NCL script.
'''

def list_files(directory,start_year,end_year,h_case):

    """
    This function extracts the files in the directory that are within the chosen dates.

    directory: string showing he directory to look for files. always end with '/'

    scenario: string used to find which kind of files to find.
                Use below example to determine the scenario:
                    filename: 'f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_NudgeOutside_new.cam.h0.2017-11.nc'
                    scenario: 'f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_NudgeOutside_new.cam.h0'
                    casename: 'f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_NudgeOutside_new'

    start_period: datetime format showing the starting date. ex: datetime.datetime(2017, 11, 1, 0, 0)
    start_year

    end_period: datetime format showing the ending date . ex: datetime.datetime(2017, 12, 1, 0, 0)
    end_year
    """

    #History file year range
    yrs = np.arange(int(start_year), int(end_year)+1)

    all_filenames = []
    for i in yrs:
        all_filenames.append(sorted(Path(directory).glob(f'*.{h_case}.{i}-*')))

    # Flattening the list of lists
    filenames = list(itertools.chain.from_iterable(sorted(all_filenames)))
    if len(filenames)==0 : sys.exit(" Directory has no outputs ")

    return filenames
#####


def Get_files(data_dir,start_year,end_year,h_case,**kwargs):

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
    ext1_SE = kwargs.pop('ext1_SE','')
    area = kwargs.pop('area',False)

    Earth_rad=6.371e6 # Earth Radius in meters

    current_dir = data_dir

    current_files = list_files(data_dir,start_year,end_year,h_case)

    # get the Lat and Lons for each case
    tmp_file = xr.open_dataset(Path(data_dir) / current_files[0])
    lon = tmp_file['lon'+ext1_SE].data
    lon[lon > 180.] -= 360 # shift longitude from 0-360˚ to -180-180˚
    lat = tmp_file['lat'+ext1_SE].data

    if area == True:
        try:
            tmp_area = tmp_file['area'+ext1_SE].data
            Earth_area = 4 * np.pi * Earth_rad**(2)

            areas = tmp_area*Earth_area/np.nansum(tmp_area)

        except KeyError:
            try:
                tmp_area = tmp_file['AREA'+ext1_SE].isel(time=0).data
                Earth_area = 4 * np.pi * Earth_rad**(2)
                areas = tmp_area*Earth_area/np.nansum(tmp_area)

            except:
                dlon = np.abs(lon[1]-lon[0])
                dlat = np.abs(lat[1]-lat[0])

                lon2d,lat2d = np.meshgrid(lon,lat)
                #area=np.zeros_like(lat2d)

                dy = Earth_rad*dlat*np.pi/180
                dx = Earth_rad*np.cos(lat2d*np.pi/180)*dlon*np.pi/180

                tmp_area = dx*dy
                areas = tmp_area


    files = current_files
    Lats = lat
    Lons = lon

    return files,Lats,Lons,areas,ext1_SE
#####

def set_dic_SE(ListVars, ext1_SE):

    # Initialize dictionary
    #----------------------
    dic_SE={}

    # Chemistry
    #----------
    dic_SE['O3']={'O3'+ext1_SE:1e9} # covert to ppb for Tropopause calculation
    dic_SE['CH4']={'CH4'+ext1_SE:1}
    dic_SE['CO']={'CO'+ext1_SE:1}

    dic_SE['ISOP']={'ISOP'+ext1_SE:1}
    dic_SE['MTERP']={'MTERP'+ext1_SE:1}
    dic_SE['CH3OH']={'CH3OH'+ext1_SE:1}
    dic_SE['CH3COCH3']={'CH3COCH3'+ext1_SE:1}
    dic_SE['CH3CCL3']={'CH3CCL3'+ext1_SE:1}


    # Aerosols
    #---------

    dic_SE['DAOD']={'AODDUSTdn'+ext1_SE:1}
    dic_SE['AOD']={'AODVISdn'+ext1_SE:1}

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


    dic_SE['SO4']={'so4_a1'+ext1_SE:1,
                   'so4_a2'+ext1_SE:1,
                   'so4_a3'+ext1_SE:1,
                   'so4_a5'+ext1_SE:1}

    # FOR SOA, first check if the integrated bins are included
    if (('soa_a1'+ext1_SE in ListVars ) & ('soa_a1'+ext1_SE in ListVars )):
        dic_SE['SOA'] = {'soa_a1'+ext1_SE:1,
                       'soa_a2'+ext1_SE:1}
    else:
        dic_SE['SOA'] = {'soa1_a1'+ext1_SE:1,
                   'soa2_a1'+ext1_SE:1,
                   'soa3_a1'+ext1_SE:1,
                   'soa4_a1'+ext1_SE:1,
                   'soa5_a1'+ext1_SE:1,
                   'soa1_a2'+ext1_SE:1,
                   'soa2_a2'+ext1_SE:1,
                   'soa3_a2'+ext1_SE:1,
                   'soa4_a2'+ext1_SE:1,
                   'soa5_a2'+ext1_SE:1}

    return dic_SE
#####

def fill_dic_SE(dic_SE, variables, ListVars, ext1_SE, AEROSOLS, MW, AVO, gr, Mwair):
    # Here, we deal with conversion factors for different components
    # Some conversion factors need density or Layer's pressure, that will be
    # accounted for when reading the files.
    # We convert everying to kg/m2/s or kg/m2 or kg/s, so that final Tg/yr or Tg results are consistent
    for var in variables:

        print("Main variable:",var,"\n")

        if 'AOD' in var:
            dic_SE[var+'_AOD']={}
        else:
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
            if var=='SO4':
                dic_SE[var+'_NUCL']={}
                dic_SE[var+'_AQS']={}


        var_keys=dic_SE[var].keys()

        for key in var_keys:
            print("key",key,"\n")
            # for CHML and CHMP:
            # original unit : [molec/cm3/s]
            # following Tilmes code to convert to [kg/m2/s]
            # conversion: Mw*rho*delP*1e3/Avo/gr
            # rho and delP will be applied when reading the files in SEbudget function.


            if 0==1:
                print("huh, broke mathematics...")
            else:

                # for AOD and DAOD:
                if 'AOD' in var:
                    if key in ListVars:

                        dic_SE[var+'_AOD'][key+ext1_SE]=1 
                    else:
                        dic_SE[var+'_AOD']['PS'+ext1_SE]=0. 

                    continue # AOD doesn't need any other budget calculations

                # for CHML and CHMP:
                # original unit : [molec/cm3/s]
                # following Tilmes code to convert to [kg/m2/s]
                # conversion: Mw*rho*delP*1e3/Avo/gr
                # rho and delP will be applied when reading the files in SEbudget function.
                if key=='O3'+ext1_SE:
                    # for O3, we should not include fast cycling reactions
                    # As a result, we use below diagnostics in the model if available. If not, we use CHML and CHMP
                    if ((key+'_Loss' in ListVars) & (key+'_Prod' in ListVars)) :
                        dic_SE[var+'_CHML'][key+'_Loss'+ext1_SE]=MW[var]*1e3/AVO/gr
                        dic_SE[var+'_CHMP'][key+'_Prod'+ext1_SE]=MW[var]*1e3/AVO/gr
                    else:
                        if key+'_CHML' in ListVars:
                            dic_SE[var+'_CHML'][key+'_CHML'+ext1_SE]=MW[var]*1e3/AVO/gr
                        else:
                            dic_SE[var+'_CHML']['U'+ext1_SE]=0

                        if key+'_CHMP' in ListVars:
                            dic_SE[var+'_CHMP'][key+'_CHMP'+ext1_SE]=MW[var]*1e3/AVO/gr
                        else:
                            dic_SE[var+'_CHMP']['U'+ext1_SE]=0

                else:

                    if key+'_CHML' in ListVars:
                        dic_SE[var+'_CHML'][key+'_CHML'+ext1_SE]=MW[var]*1e3/AVO/gr
                    else:
                        dic_SE[var+'_CHML']['U'+ext1_SE]=0

                    if key+'_CHMP' in ListVars:
                        dic_SE[var+'_CHMP'][key+'_CHMP'+ext1_SE]=MW[var]*1e3/AVO/gr
                    else:
                        dic_SE[var+'_CHMP']['U'+ext1_SE]=0


                # for SF:
                # original unit: [kg/m2/s]
                if 'SF'+key in ListVars:
                    if var=='SO4':

                        dic_SE[var+'_SF']['SF'+key+ext1_SE]=32.066/115.11
                    else:
                        dic_SE[var+'_SF']['SF'+key+ext1_SE]=1
                elif key+'SF' in ListVars:
                    dic_SE[var+'_SF'][key+ext1_SE+'SF']=1
                else:
                    dic_SE[var+'_SF']['PS'+ext1_SE]=0.


                # for CLXF:
                # original unit: [molec/cm2/s]
                # conversion: Mw*10/Avo
                if key+'_CLXF' in ListVars:
                    dic_SE[var+'_CLXF'][key+'_CLXF'+ext1_SE]=MW[var]*10/AVO  # convert [molec/cm2/s] to [kg/m2/s]
                else:
                    dic_SE[var+'_CLXF']['PS'+ext1_SE]=0.


                if var in AEROSOLS:
                    # for each species:
                    # original unit : [kg/kg]  in dry air
                    # convert to [kg/m2]
                    # conversion: delP/gr
                    # delP will be applied when reading the files in SEbudget function.
                    if key in ListVars:
                        if var=='SO4': # For SO4, we report all the budget calculation for Sulfur.
                            dic_SE[var+'_BURDEN'][key+ext1_SE]=(32.066/115.11)/gr
                        else:
                            dic_SE[var+'_BURDEN'][key+ext1_SE]=1/gr

                    else:
                        dic_SE[var+'_BURDEN']['U'+ext1_SE]=0


                    # for DDF:
                    # original unit: [kg/m2/s]
                    if key+'DDF' in ListVars:
                        if var=='SO4':
                            dic_SE[var+'_DDF'][key+ext1_SE+'DDF']=32.066/115.11
                        else:
                            dic_SE[var+'_DDF'][key+ext1_SE+'DDF']=1

                    else:
                        dic_SE[var+'_DDF']['PS'+ext1_SE]=0.


                    # for SFWET:
                    # original unit: [kg/m2/s]
                    if key+'SFWET' in ListVars:
                        if var=='SO4':
                            dic_SE[var+'_WDF'][key+ext1_SE+'SFWET']=32.066/115.11
                        else:
                            dic_SE[var+'_WDF'][key+ext1_SE+'SFWET']=1

                    else:
                        dic_SE[var+'_WDF']['PS'+ext1_SE]=0.


                    # for sfgaex1:
                    # original unit: [kg/m2/s]
                    if key+'_sfgaex1' in ListVars:
                        if var=='SO4':
                            dic_SE[var+'_GAEX'][key+ext1_SE+'_sfgaex1']=32.066/115.11
                        else:
                            dic_SE[var+'_GAEX'][key+ext1_SE+'_sfgaex1']=1

                    else:
                        dic_SE[var+'_GAEX']['PS'+ext1_SE]=0.


                    # for DDF in cloud water:
                    # original unit: [kg/m2/s]
                    cloud_key=key[:-2]+'c'+key[-1]
                    if cloud_key+ext1_SE+'DDF' in ListVars:
                        if var=='SO4':
                            dic_SE[var+'_DDFC'][cloud_key+ext1_SE+'DDF']=32.066/115.11
                        else:
                            dic_SE[var+'_DDFC'][cloud_key+ext1_SE+'DDF']=1
                    else:
                        dic_SE[var+'_DDFC']['PS'+ext1_SE]=0.


                    # for SFWET in cloud water:
                    # original unit: [kg/m2/s]
                    if cloud_key+ext1_SE+'SFWET' in ListVars:
                        if var=='SO4':
                            dic_SE[var+'_WDFC'][cloud_key+ext1_SE+'SFWET']=32.066/115.11
                        else:
                            dic_SE[var+'_WDFC'][cloud_key+ext1_SE+'SFWET']=1
                    else:
                        dic_SE[var+'_WDFC']['PS'+ext1_SE]=0.


                    if var=='SO4':
                        # for Nucleation :
                        # original unit: [kg/m2/s]
                        if key+ext1_SE+'_sfnnuc1' in ListVars:
                            dic_SE[var+'_NUCL'][key+ext1_SE+'_sfnnuc1']=32.066/115.11
                        else:
                            dic_SE[var+'_NUCL']['PS'+ext1_SE]=0.


                        # for Aqueous phase :
                        # original unit: [kg/m2/s]

                        if (('AQSO4_H2O2'+ext1_SE in ListVars) & ('AQSO4_O3'+ext1_SE in ListVars)) :
                                dic_SE[var+'_AQS']['AQSO4_H2O2'+ext1_SE]=32.066/115.11
                                dic_SE[var+'_AQS']['AQSO4_O3'+ext1_SE]=32.066/115.11
                        else:
                            # original unit: [kg/m2/s]
                            if cloud_key+'AQSO4'+ext1_SE in ListVars:
                                dic_SE[var+'_AQS'][cloud_key+'AQSO4'+ext1_SE]=32.066/115.11
                            else:
                                dic_SE[var+'_AQS']['PS'+ext1_SE]=0.

                            if cloud_key+'AQH2SO4'+ext1_SE in ListVars:
                                dic_SE[var+'_AQS'][cloud_key+'AQH2SO4'+ext1_SE]=32.066/115.11
                            else:
                                dic_SE[var+'_AQS']['PS'+ext1_SE]=0.

                else:
                    # for each species:
                    # original unit : [mole/mole]  in dry air
                    # convert to [kg/m2]
                    # conversion: Mw*delP/Mwair/gr     Mwair=28.97 gr/mole
                    # delP will be applied when reading the files in SEbudget function.
                    if key in ListVars:
                        dic_SE[var+'_BURDEN'][key+ext1_SE]=MW[var]/Mwair/gr
                    else:
                        dic_SE[var+'_BURDEN']['U'+ext1_SE]=0

                    # for DF:
                    # original unit: [kg/m2/s]
                    if 'DF_'+key in ListVars:
                        dic_SE[var+'_DDF']['DF_'+key+ext1_SE]=1
                    else:
                        dic_SE[var+'_DDF']['PS'+ext1_SE]=0.

                    # for WD:
                    # original unit: [kg/m2/s]
                    if 'WD_'+key in ListVars:
                        dic_SE[var+'_WDF']['WD_'+key+ext1_SE]=1
                    else:
                        dic_SE[var+'_WDF']['PS'+ext1_SE]=0.

                    # for Chem tendency:
                    # original unit: [kg/s]
                    # conversion: not needed
                    if 'D'+key+'CHM' in ListVars:
                        dic_SE[var+'_TEND']['D'+key+'CHM'+ext1_SE]=1  # convert [kg/s] to [kg/s]
                    else:
                        dic_SE[var+'_TEND']['U'+ext1_SE]=0

                    # for Lightning NO production: (always in gas)
                    # original unit: [Tg N/Yr]
                    # conversion: not needed
                    if 'LNO_COL_PROD' in ListVars:
                        dic_SE[var+'_LNO']['LNO_COL_PROD'+ext1_SE]=1  # convert [Tg N/yr] to [Tg N /yr]
                    else:
                        dic_SE[var+'_LNO']['PS'+ext1_SE]=0
    return dic_SE
#####


def make_Dic_scn_var_comp(adfobj, variables, data_dir, dic_SE, Files, ext1_SE, AEROSOLS):
    """
    the LNO is lightning NOx, which should be reported explicitly rather as CO_LNO, O3_LNO, ...
    """

    # this is for finding tropospheric values
    Dic_crit={}

    missing_vars_tot = []
    needed_vars_tot = []

    current_dir=data_dir
    #current_scn=scenario
    current_files=Files#[current_scn]

    #Dic_scn_var_comp[current_scn]={}
    Dic_var_comp={}

    for ivar,current_var in enumerate(variables):
        if 'AOD' in current_var:
            components=[current_var+'_AOD']
        else:
            if current_var in AEROSOLS: # AEROSOLS
                print("\n current var:",current_var)

                # Components are: burden, chemical loss, chemical prod, dry deposition,
                #                 surface emissions, elevated emissions, wet deposition, gas-aerosol exchange
                components=[current_var+'_BURDEN',current_var+'_CHML',current_var+'_CHMP',
                            current_var+'_DDF',current_var+'_WDF', current_var+'_SF', current_var+'_CLXF',
                            current_var+'_DDFC',current_var+'_WDFC']

                if current_var=='SO4':
                    # For SULF we also have AQS, nucleation, and strat-trop gas exchange
                    components.append(current_var+'_AQS')
                    components.append(current_var+'_NUCL')
                    components.append(current_var+'_GAEX')
                    components.remove(current_var+'_CHMP')

                    #components.append(current_var+'_CLXF') # BRT -  CLXF is added above.
                if current_var == "SOA":
                    components.append(current_var+'_GAEX')
            #End if - AEROSOLS

            else: # CHEMS
                # Components are: burden, chemical loss, chemical prod, dry/wet deposition,
                #                 surface emissions, elevated emissions, chemical tendency
                # I always add Lightning NOx production when calculating O3 budget.

                components=[current_var+'_BURDEN',current_var+'_CHML',current_var+'_CHMP',
                            current_var+'_DDF',current_var+'_WDF', current_var+'_SF', current_var+'_CLXF',
                            current_var+'_TEND']

                if current_var =="O3":
                        components.append(current_var+'_LNO')


        print("components",components,"\n")
        Dic_comp={}
        for comp in components:
            print(f"Component: {comp} for main variable, {current_var}\n")
            current_data,missing_vars,needed_vars = SEbudget(dic_SE,current_dir,current_files,comp,ext1_SE)

            for var_m in missing_vars:
                if var_m not in missing_vars_tot:
                    missing_vars_tot.append(var_m)

            for var_n in needed_vars:
                if var_n not in needed_vars_tot:
                    needed_vars_tot.append(var_n)

    #TODO: check this section to see if it can't be better run
            Dic_comp[comp] = current_data
        Dic_var_comp[current_var] = Dic_comp
    #Dic_scn_var_comp[current_scn]= Dic_var_comp
    Dic_scn_var_comp = Dic_var_comp

    # Critical threshholds?
    # Just run this once
    current_crit=SEbudget(dic_SE,current_dir,current_files,'O3',ext1_SE)
    #current_crit=SEbudget_dask(dic_SE[i],current_dir,current_files,'O3')
    #Dic_crit[current_scn]=current_crit
    Dic_crit=current_crit

    print("missing_vars_tot: ",missing_vars_tot,"\n")
    msg = f"chem/aerosol tables:"
    msg += f"\n\t - missing variables from budget? {missing_vars_tot}"
    adfobj.debug_log(msg)

    print("needed_vars_tot: ",needed_vars_tot,"\n")
    msg = f"chem/aerosol tables:"
    msg += f"\n\t - needed variables for budget? {needed_vars_tot}"
    adfobj.debug_log(msg)

    return Dic_crit,Dic_scn_var_comp
#####


def SEbudget(dic_SE,data_dir,files,var,ext1_SE,**kwargs):

    """
    Function used for getting the data for the budget calculation. This is the
    chunk of code that takes the longest by far.

    Example:
    ~70/75 mins per case for 9 years
    ** This is for both chemistry and aeorosl calculations


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
     -> this will be the individual componnent, ie O3_CHMP, SOA_WDF, etc.
       => this could be a place to limit unneccessary calculations??

    """
    print(var," - ",dic_SE[var].keys())
    # gas constanct
    Rgas=287.04 #[J/K/Kg]=8.314/0.028965

    missing_vars = []
    needed_vars = []

    all_data=[]
    for file in range(len(files)):

        ds=xr.open_dataset(Path(data_dir) / files[file])


        if file==0:
            mock_2d=np.zeros_like(np.array(ds['PS'+ext1_SE].isel(time=0)))
            mock_3d=np.zeros_like(np.array(ds['U'+ext1_SE].isel(time=0)))


        data=[]
        for i in dic_SE[var].keys():
            if i not in needed_vars:
                needed_vars.append(i)
            if file == 0:
                print("Variable:",i)

            if ((i!='PS'+ext1_SE) and (i!='U'+ext1_SE) ) :

                data.append(np.array(ds[i].isel(time=0))*dic_SE[var][i])

            else:
                if i=='PS'+ext1_SE:
                    data.append(mock_2d)
                else:
                    data.append(mock_3d)

                if file == 0:
                #     print(f"An error occurred: {e}")
                    print(f"No variable was found for: {var}")

                if var not in missing_vars:
                    missing_vars.append(var)

        data=np.sum(data,axis=0)

        try:
            delP=np.array(ds['PDELDRY'+ext1_SE].isel(time=0))
        except:

            hyai=np.array(ds['hyai'])
            hybi=np.array(ds['hybi'])

            try:
                PS=np.array(ds['PSDRY'+ext1_SE].isel(time=0))
            except:
                PS=np.array(ds['PS'+ext1_SE].isel(time=0))

            P0=1e5
            Plevel=np.zeros_like(np.array(ds['U'+ext1_SE]))


            for i in range(len(Plevel)):
                Plevel[i]=hyai[i]*P0+hybi[i]*PS

            delP=Plevel[1:]-Plevel[:-1]


        if ('CHML' in var) or ('CHMP' in var) :
            Temp=np.array(ds['T'+ext1_SE].isel(time=0))
            Pres=np.array(ds['PMID'+ext1_SE].isel(time=0))
            rho= Pres/(Rgas*Temp)

            data=data*delP/rho

        elif ('BURDEN' in var):
            data=data*delP

        else:
            data=data

        all_data.append(data)

    all_data=np.nanmean(all_data,axis=0)

    return all_data,missing_vars,needed_vars
#####


def SEbudget_dask(dic_SE,data_dir,files,var,ext1_SE,**kwargs):

    """
    Function used for getting the data for the budget calculation. This is the
    chunk of code that takes the longest by far.

    Example:
    ~70/75 mins per case for 9 years
    ** This is for both chemistry and aeorosl calculations


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
     -> this will be the individual componnent, ie O3_CHMP, SOA_WDF, etc.
       => this could be a place to limit unneccessary calculations??
    """

    def preprocess(ds):
        return ds[variables]

    print("Ahhhhh: ",var," - ",dic_SE[var].keys())
    # gas constanct
    Rgas=287.04 #[J/K/Kg]=8.314/0.028965


    all_files=[]
    for file in files:
        all_files.append(data_dir+file)


    variables=list(dic_SE[var].keys())
    variables+=['U','hyai','hybi','PDELDRY','PS','PSDRY','T','PMID']
    ds = xr.open_mfdataset(
        all_files,
        engine="netcdf4",
        parallel=True,
        concat_dim="time",
        combine="nested",
        preprocess=preprocess,
     chunks={'time': 1})


    mock_2d=np.zeros_like(np.array(ds['PS'].isel(time=0)))
    mock_3d=np.zeros_like(np.array(ds['U'].isel(time=0)))

    missing_vars = []
    needed_vars = []


    data=[]
    for i in dic_SE[var].keys():
        if i not in needed_vars:
            needed_vars.append(i)
        if file == 0:
            print("Variable:",i)

        if ((i!='PS'+ext1_SE) and (i!='U'+ext1_SE) ) :

            data.append(np.array(ds[i].isel(time=0))*dic_SE[var][i])

        else:
            if i=='PS'+ext1_SE:
                data.append(mock_2d)
            else:
                data.append(mock_3d)

            if file == 0:
            #     print(f"An error occurred: {e}")
                print(f"No variable was found for: {var}")

            if var not in missing_vars:
                missing_vars.append(var)

    if len(data)==0:
        # to make sure there is an array with correct shape
        data=np.array(ds[i])*dic_SE[var][i]
    else:
        data=np.sum(data,axis=0)

    try:
        delP=np.array(ds['PDELDRY'+ext1_SE])
    except:

        hyai=np.array(ds['hyai'])
        hybi=np.array(ds['hybi'])

        try:
            PS=np.array(ds['PSDRY'+ext1_SE])
        except:
            PS=np.array(ds['PS'+ext1_SE])

        P0=1e5

        Plevel=np.zeros_like(np.array(ds['U'+ext1_SE]))
        for i in range(len(Plevel[0])):
            Plevel[:,i,:,:]=hyai[:,i]*P0+hybi[:,i]*PS

        delP=Plevel[:,1:,:,:]-Plevel[:,:-1,:,:]

    if ('CHML' in var) or ('CHMP' in var) :
        Temp=np.array(ds['T'+ext1_SE])
        Pres=np.array(ds['PMID'+ext1_SE])
        rho= Pres/(Rgas*Temp)
        data=data*delP/rho
    elif (('BURDEN' in var) or ('MEAN' in var)) :
        data=data*delP
    else:
        data=data

    all_data=np.nanmean(data,axis=0)

    return all_data,missing_vars,needed_vars
#####


def calc_budget_data(current_var, Dic_scn_var_comp, area, trop, inside, num_yrs, duration, AEROSOLS):

    #duration = durations[case]
    chem_dict = {}
    #chem_dict[case] = {}

    if current_var == 'SO4':
        specifier = ' S'
    else:
        specifier = ''


    if 'AOD' in current_var:
        # Burden
        spc_burd = Dic_scn_var_comp[current_var][current_var+'_AOD']
        burden = np.ma.masked_where(inside==False,spc_burd)  #convert Kg/m2 to Tg
        BURDEN = np.ma.sum(burden*area)/np.ma.sum(area)
        chem_dict[f"{current_var}_mean"] = np.round(BURDEN,5)
    else:
        # Surface Emissions
        #print("Surface Emissions")
        spc_sf = Dic_scn_var_comp[current_var][current_var+'_SF']
        tmp_sf = spc_sf
        sf = np.ma.masked_where(inside==False,tmp_sf*area)  #convert Kg/m2/s to Tg/yr
        SF = np.ma.sum(sf*duration*1e-9)/num_yrs
        chem_dict[f"{current_var}_EMIS (Tg{specifier}/yr)"] = np.round(SF,5)

        # Elevated Emissions
        spc_clxf = Dic_scn_var_comp[current_var][current_var+'_CLXF']
        tmp_clxf = spc_clxf
        clxf = np.ma.masked_where(inside==False,tmp_clxf*area)  #convert Kg/m2/s to Tg/yr
        CLXF = np.ma.sum(clxf*duration*1e-9)/num_yrs
        chem_dict[f"{current_var}_EMIS_elevated (Tg{specifier}/yr)"] = np.round(CLXF,5)

        # Burden
        spc_burd = Dic_scn_var_comp[current_var][current_var+'_BURDEN']
        spc_burd = np.where(np.isnan(trop),np.nan,spc_burd)
        tmp_burden = np.nansum(spc_burd*area,axis=0)
        burden = np.ma.masked_where(inside==False,tmp_burden)  #convert Kg/m2 to Tg
        BURDEN = np.ma.sum(burden*1e-9)
        chem_dict[f"{current_var}_BURDEN (Tg{specifier})"] = np.round(BURDEN,5)

        # Chemical Loss
        spc_chml = Dic_scn_var_comp[current_var][current_var+'_CHML']
        spc_chml = np.where(np.isnan(trop),np.nan,spc_chml)
        tmp_chml = np.nansum(spc_chml*area,axis=0)
        chml = np.ma.masked_where(inside==False,tmp_chml)  #convert Kg/m2/s to Tg/yr
        CHML = np.ma.sum(chml*duration*1e-9)/num_yrs
        chem_dict[f"{current_var}_CHEM_LOSS (Tg{specifier}/yr)"] = np.round(CHML,5)

        # Chemical Production
        if current_var == 'SO4': # chemical production is basically the elevated emissions.
                               # We have removed it for SO4 budget. and put 0 here, so, we don't report it
            chem_dict[f"{current_var}_CHEM_PROD (Tg{specifier}/yr)"] = 0
        else:
            spc_chmp = Dic_scn_var_comp[current_var][current_var+'_CHMP']
            spc_chmp = np.where(np.isnan(trop),np.nan,spc_chmp)
            tmp_chmp = np.nansum(spc_chmp*area,axis=0)
            chmp = np.ma.masked_where(inside==False,tmp_chmp)  #convert Kg/m2/s to Tg/yr
            CHMP = np.ma.sum(chmp*duration*1e-9)/num_yrs
            chem_dict[f"{current_var}_CHEM_PROD (Tg{specifier}/yr)"] = np.round(CHMP,5)


        if current_var in AEROSOLS:

           # Dry Deposition Flux
            spc_ddfa = Dic_scn_var_comp[current_var][current_var+'_DDF']
            spc_ddfc = Dic_scn_var_comp[current_var][current_var+'_DDFC']
            spc_ddf = spc_ddfa +spc_ddfc
            tmp_ddf = spc_ddf
            ddf = np.ma.masked_where(inside==False,tmp_ddf*area)  #convert Kg/m2/s to Tg/yr
            DDF = np.ma.sum(ddf*duration*1e-9)/num_yrs
            chem_dict[f"{current_var}_DRYDEP (Tg{specifier}/yr)"] = np.round(DDF,5)

            # Wet deposition
            spc_wdfa = Dic_scn_var_comp[current_var][current_var+'_WDF']
            spc_wdfc = Dic_scn_var_comp[current_var][current_var+'_WDFC']
            spc_wdf = spc_wdfa +spc_wdfc
            tmp_wdf = spc_wdf
            wdf = np.ma.masked_where(inside==False,tmp_wdf*area)  #convert Kg/m2/s to Tg/yr
            WDF = np.ma.sum(wdf*duration*1e-9)/num_yrs
            chem_dict[f"{current_var}_WETDEP (Tg{specifier}/yr)"] = np.round(WDF,5)

            if current_var in ["SOA",'SO4']:
                # gas-aerosol Exchange
                spc_gaex = Dic_scn_var_comp[current_var][current_var+'_GAEX']
                tmp_gaex = spc_gaex
                gaex = np.ma.masked_where(inside==False,tmp_gaex*area)  #convert Kg/m2/s to Tg/yr
                GAEX = np.ma.sum(gaex*duration*1e-9)/num_yrs
                chem_dict[f"{current_var}_GAEX (Tg{specifier}/yr)"] = np.round(GAEX,5)

            # LifeTime = Burden/(loss+deposition)
            LT = BURDEN/(CHML+DDF-WDF)* duration/86400/num_yrs # days
            chem_dict[f"{current_var}_LIFETIME (days)"] = np.round(LT,5)

            if current_var == 'SO4':
                # Aqueous Chemistry
                spc_aqs = Dic_scn_var_comp[current_var][current_var+'_AQS']
                tmp_aqs = spc_aqs
                aqs = np.ma.masked_where(inside==False,tmp_aqs*area)  #convert Kg/m2/s to Tg/yr
                AQS = np.ma.sum(aqs*duration*1e-9)/num_yrs
                chem_dict[f"{current_var}_AQUEOUS (Tg{specifier}/yr)"] = np.round(AQS,5)

                # Nucleation
                spc_nucl = Dic_scn_var_comp[current_var][current_var+'_NUCL']
                tmp_nucl = spc_nucl
                nucl = np.ma.masked_where(inside==False,tmp_nucl*area)  #convert Kg/m2/s to Tg/yr
                NUCL = np.ma.sum(nucl*duration*1e-9)/num_yrs
                chem_dict[f"{current_var}_NUCLEATION (Tg{specifier}/yr)"] = np.round(NUCL,5)

        else:
            # Dry Deposition Flux
            #print("Dry Deposition Flux")
            spc_ddf = Dic_scn_var_comp[current_var][current_var+'_DDF']
            tmp_ddf = spc_ddf
            ddf = np.ma.masked_where(inside==False,tmp_ddf*area)  #convert Kg/m2/s to Tg/yr
            DDF = np.ma.sum(ddf*duration*1e-9)/num_yrs
            chem_dict[f"{current_var}_DRYDEP (Tg/yr)"] = np.round(DDF,5)

            # Wet Deposition Flux
            #print("Wet Deposition Flux")
            spc_wdf = Dic_scn_var_comp[current_var][current_var+'_WDF']
            tmp_wdf = spc_wdf
            wdf = np.ma.masked_where(inside==False,tmp_wdf*area)  #convert Kg/m2/s to Tg/yr
            WDF = -1*np.ma.sum(wdf*duration*1e-9)/num_yrs
            chem_dict[f"{current_var}_WETDEP (Tg/yr)"] = np.round(WDF,5)

            # Total Deposition
            TDEP = DDF - WDF
            chem_dict[f"{current_var}_TDEP (Tg/yr)"] = np.round(TDEP,5)

            # LifeTime = Burden/(loss+deposition)
            #print("LifeTime")
            if current_var == "CH4":
                LT = BURDEN/(CHML+DDF-WDF) # years
                chem_dict[f"{current_var}_LIFETIME (years)"] = np.round(LT,5)
            else:
                if (CHML+DDF-WDF) > 0:
                    if CHML != 0:
                        LT = BURDEN/(CHML+DDF-WDF)*duration/86400/num_yrs # days
                        chem_dict[f"{current_var}_LIFETIME (days)"] = np.round(LT,5)
                    else:
                        # do not report lifetime if chemical loss (for gases) is not included in the model outputs
                        # and put 0 here, so, we don't report it
                        chem_dict[f"{current_var}_LIFETIME (days)"] = 0

            #NET = CHMP-CHML
            # Chemical Tendency
            #print("Chemical Tendency")
            spc_tnd = Dic_scn_var_comp[current_var][current_var+'_TEND']
            spc_tnd = np.where(np.isnan(trop),np.nan,spc_tnd)
            tmp_tnd = np.nansum(spc_tnd,axis=0)
            tnd = np.ma.masked_where(inside==False,tmp_tnd)  #convert Kg/s to Tg/yr
            TND = np.ma.sum(tnd*duration*1e-9)/num_yrs
            chem_dict[f"{current_var}_TEND (Tg/yr)"] = np.round(TND,5)


            if current_var == "O3":

                # Stratospheric-Tropospheric Exchange
                #print("Stratospheric-Tropospheric Exchange")
                STE = DDF-TND
                chem_dict[f"{current_var}_STE (Tg/yr)"] = np.round(STE,5)

                # Lightning NOX production
                #print("Lightning NOX production")
                spc_lno = Dic_scn_var_comp[current_var][current_var+'_LNO']
                tmp_lno = np.ma.masked_where(inside==False,spc_lno)
                LNO = np.ma.sum(tmp_lno)
                chem_dict[f"{current_var}_LNO (Tg N/yr)"] = np.round(LNO,5)

    return chem_dict
#####


def make_table(adfobj, vars, chem_type, Dic_scn_var_comp, areas, trops, case_names, nicknames, durations, insides, num_yrs, AEROSOLS):
    # Initialize an empty dictionary to store DataFrames
    dfs = {}

    for case in case_names:
        nickname = nicknames[case]
        # Collect row data in a list of dictionaries
        durations[case]
        rows = []
        for current_var in vars:
            chem_dict = calc_budget_data(current_var, Dic_scn_var_comp[case], areas[case], trops[case], insides[case], num_yrs[case], durations[case], AEROSOLS)
            #chem_dict = calc_budget_data(case, current_var, Dic_scn_var_comp, area, trop, inside, num_yrs, durations, AEROSOLS)

            for key, val in chem_dict.items():
                if val != 0:  # Skip variables with a value of 0
                    #rows.append({'variable': key, case: np.round(val, 3)})
                    rows.append({'variable': key, nickname: np.round(val, 3)})
                else:
                    msg = f"chem/aerosol tables:"
                    msg += f"\n\t - Variable '{key}' has value of 0, will not add to table"
                    adfobj.debug_log(msg)

        # Create the DataFrame for the current case
        table_df = pd.DataFrame(rows)

        if chem_type == 'gases':
            # Replace compound names directly in the DataFrame
            replacements = {
                'MTERP': 'Monoterpene',
                'CH3OH': 'Methanol',
                'CH3COCH3': 'Acetone',
                'O3_LNO': 'LNOx_PROD'
            }
            table_df['variable'] = table_df['variable'].replace(replacements)

        # Store the DataFrame in the dictionary
        dfs[nickname] = table_df

    # Merge the DataFrames on the 'variable' column
    merged_df = pd.merge(dfs[nicknames[case_names[0]]], dfs[nicknames[case_names[1]]], on='variable')
    print("merged_df.keys()",merged_df.keys(),"\n")

    # Calculate the differences between case columns
    merged_df['difference'] = merged_df[nicknames[case_names[0]]] - merged_df[nicknames[case_names[1]]]

    # Optional: Save the result to a new CSV file
    merged_df.to_csv(f'ADF_amwg_{chem_type}_table.csv', index=False)
    #Add comparison table dataframe to website (if enabled):
    adfobj.add_website_data(merged_df, case, case, plot_type=f"{chem_type} Tables")
    #adf.add_website_data(table_df, case_name, case_name, plot_type="Tables")

    return merged_df
#####