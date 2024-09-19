import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
#os.chdir('/glade/work/nadavis/waccm_tuning')
#os.getcwd()
import averaging_functions
import tem_functions
import interpolating_functions
import qbo_functions

import xarray as xr
from pathlib import Path
import glob
from itertools import chain

import plotting_functions as pf




def vert_seasonal_cycle(adfobj):

    #Special ADF variable which contains the output paths for
    #all generated plots and tables for each case:
    plot_locations = adfobj.plot_location

    diag_loc = adfobj.basic_info_dict['diag_loc']

    #CAM simulation variables (this is always assumed to be a list):
    case_names = adfobj.get_cam_info("cam_case_name", required=True)

    #Grab case years
    syear_cases = adfobj.climo_yrs["syears"]
    eyear_cases = adfobj.climo_yrs["eyears"]

    #Grab all case nickname(s)
    test_nicknames = adfobj.case_nicknames["test_nicknames"]

    #Grab h4 history files locations
    cam_hist_locs = adfobj.get_cam_info("cam_hist_loc", required=True)

    #Check if comparing to observations
    if adfobj.get_basic_info("compare_obs"):
        var_obs_dict = adfobj.var_obs_dict
        #If dictionary is empty, then there are no observations, so quit here:
        if not var_obs_dict:
            print("No observations found to plot against, so no TEM will be generated.")
            return

        base_name = "Obs"
    else:
        base_name = adfobj.get_baseline_info("cam_case_name", required=True)
        cam_hist_locs.append(adfobj.get_baseline_info("cam_hist_loc", required=True))

        #Extract baseline years (which may be empty strings if using Obs):
        syear_baseline = adfobj.climo_yrs["syear_baseline"]
        eyear_baseline = adfobj.climo_yrs["eyear_baseline"]

        base_nickname = adfobj.case_nicknames["base_nickname"]

        test_nicknames.append(base_nickname)
        case_names.append(base_name)
        syear_cases.append(syear_baseline)
        eyear_cases.append(eyear_baseline)
        #overwrite_tem_cases.append(tem_opts.get("overwrite_tem_base", False))
    #End if

    res = adfobj.variable_defaults # will be dict of variable-specific plot preferences
    # or an empty dictionary if use_defaults was not specified in YAML.

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')
    print(f"\t NOTE: Plot type is set to {plot_type}")

    # check if existing plots need to be redone
    redo_plot = adfobj.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")
    #-----------------------------------------


    load_dir = cam_hist_locs
    #load_dir = [adfobj.get_cam_info("cam_hist_loc"),adfobj.get_baseline_info("cam_hist_loc")]

    runs = case_names #+ [data_name]
    #runs = [adfobj.get_cam_info("cam_case_name")[0],adfobj.get_baseline_info("cam_case_name")]
    run_names = test_nicknames #+ [base_nickname]
    #run_names = [adfobj.get_cam_info("case_nickname")[0],adfobj.get_baseline_info("case_nickname")]
    ['taubgnd5.energy','nudged_dst11']
    #var_list = ['U','T','V','lat','lev','time','date']
    #var_list = ['U','T','V','UU','VU','VT','OMEGAU','OMEGA','O3','Q','UTEND_VDIFF','lat','lev','time','date']
    #calc_var_list = ['U','T','V','O3','Q','UTEND_VDIFF']
    calc_var_list = ['U','T','V','Q']
    var_list = calc_var_list + ['lat','lev','time', 'date']
    #tem_vars = ['vtem','omegatem','delF','Fy','Fz']
    #syr = adfobj.climo_yrs['syears'][0]
    #eyr = adfobj.climo_yrs['eyears'][0]

    #syr = adfobj.climo_yrs['syear_baseline']
    #eyr = adfobj.climo_yrs['eyear_baseline']

    wowsa = {}
    for idx,run in enumerate(runs[:-1]):
        wowsa[run] = (adfobj.climo_yrs['syears'][idx],adfobj.climo_yrs['eyears'][idx])
    wowsa[runs[-1]] = (adfobj.climo_yrs['syear_baseline'],adfobj.climo_yrs['eyear_baseline'])

    mrr = xr.open_dataset("MERRA2_met.nc")
    merra_dates_str = mrr.TemporalRange
    merra_start_yr = merra_dates_str[0:4]
    merra_end_yr = merra_dates_str[-10:-6]
    merra_start_yr,merra_end_yr


    data_runs = {}
    data_runs_seasonal = {}
    data_runs_monthly = {}

    #var_list = ['U','T','V','lat','lev','time']
    #var_list = var_list+["date"]

    #Load all data, process TEM vars, seasonal averaging
    for idx,case_name in enumerate(case_names):



        #hist_loc = adfobj.get_cam_info("cam_hist_loc")[0]
        hist_loc = cam_hist_locs[case_name]
        #case_name = adfobj.get_cam_info("cam_case_name")[0]
        #syr = adfobj.climo_yrs['syears'][0]
        syr = syear_cases
        #eyr = adfobj.climo_yrs['eyears'][0]
        eyr = syear_cases
        h0_lists = []
        for yr in np.arange(int(syr),int(eyr)+1):
            h0_lists.append(sorted(glob.glob(f'{hist_loc}*cam.h0.{yr}-*')))

        h0_list = list(chain(*h0_lists))




        #run = run[0]
        print(run)
        #rundir = adfobj.basic_info_dict['diag_loc'] + f"zm/{run}/" #+ "/atm/hist/"
        #file = rundir + "waccm_135.nc"
        ncfile = xr.open_mfdataset(h0_list)
        #ncfile = nc.Dataset(file, 'r')
        #ncfile = ahh
        data_run_local = {}
        data_run_local_seasonal = {}
        data_run_local_monthly = {}
        for index, var in enumerate(var_list):
            #data_run_local[var] =  np.squeeze(np.array(ncfile.variables[var][:]))
            data_run_local[var] =  np.squeeze(ncfile[var].values)
            #if index < len(var_list)-4:
            if var in calc_var_list:
                data_run_local_seasonal[var]=averaging_functions.calc_seasonal_average(data_run_local[var])
                data_run_local_monthly[var]=averaging_functions.calc_monthly_average(data_run_local[var])

        #data_run_local = tem_functions.process_tem(data_run_local,data_run_local['U'],data_run_local['V'],data_run_local['T'],
                                                #data_run_local['OMEGA'],data_run_local['OMEGAU'],data_run_local['VU'],
                                                #data_run_local['VT'],data_run_local['lat'],data_run_local['lev'])
        #for index, var in enumerate(tem_vars):
        #    data_run_local_seasonal[var]=averaging_functions.calc_seasonal_average(data_run_local[tem_vars[index]])
        data_runs[run] = data_run_local
        data_runs_monthly[run] = data_run_local_monthly
        data_runs_seasonal[run] = data_run_local_seasonal
        
    merra2 = {}
    merra2_seasonal = {}
    merra2_monthly = {}
    merra2_vars = ['U','T','V','lat','lev']
    ncfile = nc.Dataset("MERRA2_met.nc", 'r')
    for index, var in enumerate(merra2_vars):
        merra2[var] =  np.squeeze(np.array(ncfile.variables[merra2_vars[index]][:]))
        if index < len(merra2_vars)-2:
                merra2_seasonal[var] = averaging_functions.calc_seasonal_average(merra2[var])
                merra2_monthly[var] = averaging_functions.calc_monthly_average(merra2[var])

    saber = {}
    saber_seasonal = {}
    saber_monthly = {}
    saber_vars = ['u','temp','latitude','pressure']
    ncfile = nc.Dataset("SABER_monthly_2002-2014.nc", 'r')
    for index, var in enumerate(saber_vars):
        saber[var] =  np.squeeze(np.array(ncfile.variables[saber_vars[index]][:]))
        if index < len(saber_vars)-2:
                saber_seasonal[var] = averaging_functions.calc_seasonal_average(saber[var])
                saber_monthly[var] = averaging_functions.calc_monthly_average(saber[var])

    swoosh = {}
    swoosh_seasonal = {}
    swoosh_monthly = {}
    ncfile = nc.Dataset("swoosh-v02.7-198401-202301-latpress-2.5deg-L31.nc", 'r')
    swoosh_vars = ['combinedanomfillo3q','lat','level']
    for index, var in enumerate(swoosh_vars):
        swoosh[var] =  np.squeeze(np.array(ncfile.variables[swoosh_vars[index]][:]))
        if index < len(swoosh_vars)-2:
                swoosh_seasonal[var] = averaging_functions.calc_seasonal_average(swoosh[var])
                swoosh_monthly[var] = averaging_functions.calc_monthly_average(swoosh[var])   

        




    
    #Seasonal Averages
    mseasons[s] = pf.seasonal_mean(mdata, season=s, is_climo=True)
    oseasons[s] = pf.seasonal_mean(odata, season=s, is_climo=True)





def make_zm_files(hist_loc,diag_loc,case_name,syr,eyr,result_string):
    h0_lists = []
    for yr in np.arange(int(syr),int(eyr)+1):
        h0_lists.append(sorted(glob.glob(f'{hist_loc}*cam.h0.{yr}-*')))

    h0_list = list(chain(*h0_lists))

    zm_loc = Path(diag_loc + "zm/" + case_name)
    #Check if plot output directory exists, and if not, then create it:
    if not zm_loc.is_dir():
        print(f"    {zm_loc} not found, making new directory")
        zm_loc.mkdir(parents=True)

    waccm_135 = xr.open_mfdataset(h0_list)
    waccm_135.to_netcdf(f'{zm_loc}/waccm_135_2.nc')


    

    """for idx,fil in enumerate(h0_list):
        fil_name = Path(fil).parts[-1]
        #if (idx == 0) or (idx == len(h0_list)):
        #    print(fil,"\n")
        out_fil = f"{zm_loc}/{fil_name.replace('.nc','.zm.nc')}"
        #if (idx == 0) or (idx == len(h0_list)):
        #    print(out_fil)  
        #!ncwa -a lon '{fil}' '{out_fil}'
        !ncwa -a lon -v '{result_string}' '{fil}' '{out_fil}'
        #for var in var_list:
        #   !ncwa -a lon -v '{var}' '{fil}' '{out_fil}'
    """

    return zm_loc