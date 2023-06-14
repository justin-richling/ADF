from datetime import datetime, timedelta
import os
import xarray as xr
from scipy.interpolate import griddata
import matplotlib.pyplot as plt # Core library for plotting
import sys
import numpy as np
from cftime import DatetimeNoLeap
from netCDF4 import Dataset
import pandas as pd
from pathlib import Path



def amwg_chem_table(adf):
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
            return

        else:
            #Append to case list:
            case_names.append(baseline_name)
            input_ts_locs.append(input_ts_baseline)

        #Save the baseline to the first case's plots directory:
        output_locs.append(output_locs[0])
    
  

    #Convert output location string to a Path object:
    output_location = Path(output_locs[0])


    # List of directories (each item must end with "/")

    # data_dirs=['/glade/campaign/acom/acom-weather/behroozr/f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_new/atm/hist/',
    #             '/glade/campaign/acom/acom-weather/behroozr/f.e22.FCnudged.ne0np4.India07.ne30x1_ne30x1_mt12/atm/hist/']

    #data_dirs=['/glade/campaign/acom/acom-weather/behroozr/f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_new/atm/hist/']

    #data_dirs=['/glade/scratch/tilmes/archive/f.cesm3_cam058_mom_e.FCHIST.ne30_L58.26c_non-orogw_off.001/atm/hist/']

    data_dirs=['/glade/scratch/tilmes/archive/FCnudged_MAM4_f09.carma_trop_strat.aqchem.2001_2020.atom/atm/hist/']



    #data_dirs=['/glade/campaign/acom/acom-weather/behroozr/f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_new/atm/hist/']
    #data_dirs=['/glade/scratch/jzhan166/archive/f.e21.FWscHIST.ne30_L58_BL10_cam6_3_019_plus_CESM2.2.003_zm2_chemistry.006.hf/atm/hist/']
    #data_dirs=['/glade/campaign/acom/acom-climate/tilmes/CO_CONUS/f.e22.FCcotagsNudged.ne0CONUSne30x8.cesm220.2012-01/atm/hist/']



    # List of scenarios. Must be compatible with the files name. 
    # for example, the following scenario name is for a sample filename of 'f.e22.FCcotagsNudged.ne0CONUSne30x8.cesm220.2012-01.cam.h1.{YYYY}-{MM}.nc 
    # In other words, it must be the filename exluding ".{YYYY}-{MM}.nc"

    # scenarios=['f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_new.cam.h0',
    #            'f.e22.FCnudged.ne0np4.India07.ne30x1_ne30x1_mt12.cam.h0']

    #scenarios=['f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_new.cam.h0']

    #scenarios=['f.cesm3_cam058_mom_e.FCHIST.ne30_L58.26c_non-orogw_off.001.cam.h0']
    scenarios=['FCnudged_MAM4_f09.carma_trop_strat.aqchem.2001_2020.atom.cam.h0']

    #scenarios=['f.e22.FCnudged.ne0np4.India07.ne30x8_ne30x8_mt12_new.cam.h0']
    #scenarios=['f.e21.FWscHIST.ne30_L58_BL10_cam6_3_019_plus_CESM2.2.003_zm2_chemistry.006.hf.cam.h1']

    #scenarios=['f.e22.FCcotagsNudged.ne0CONUSne30x8.cesm220.2012-01.cam.h1']





    # List of labels for printing and plotting uses

    # labels=['ne30x8',
    #         'ne30x1']

    #labels=['ne30x8']

    labels=['camChem']

    #labels=['ne30x1']

    #labels=['ne30x8']



    # List of scrip files.
    scrip_files=['/glade/work/behroozr/VRM_files/ne0np4.India07.ne30x8/grids/India07.ne30x8_np4_SCRIP.nc',
                '/glade/work/behroozr/VRM_files/ne0np4.India07.ne30x1/grids/India07.ne30x1_np4_SCRIP.nc']

    #scrip_files=['/glade/work/behroozr/VRM_files/ne0np4.India07.ne30x8/grids/India07.ne30x8_np4_SCRIP.nc']
    #scrip_files=['/glade/p/acom/MUSICA/grids/ne30np4/ne30np4_091226_pentagons.nc']
    #scrip_files=['/glade/p/acom/MUSICA/grids/ne0CONUSne30x8/ne0CONUS_ne30x8_np4_SCRIP.nc']





    # In CAM-Chem (or MUSICA-v0), user can save the outputs for only a box region.
    # ext1_SE: string specifying if the files are for only a region, which changes to variable names.
    # ex: if you saved files for only a box region ($LL_lat$,$LL_lon$,$UR_lat$,$UR_lon$),
    #     the 'lat' variable will be saved as: 'lat_$LL_lon$e_to_$UR_lon$e_$LL_lat$n_to_$UR_lat$n'
    #     for instance: 'lat_65e_to_91e_20n_to_32n'
    ext1_SE=''
    #ext1_SE='_65e_to_91e_20n_to_32n'








    # list of the variables to be caculated. 
    #variables=["CH4","CH3CCL3","CO","O3","ISOP","C10H16","CH3OH","CH3COCH3"]
    #variables=["O3","CH4","ISOP"]
    #variables=["SOA",'O3','CH4']
    #variables=["SOA",'CH4','SALT','DUST','POM','BC','SULF','CH3CCL3','CO','ISOP']
    variables=["SALT",'CO']
    #variables=["O3",'CH4','ISOP']


    # if True, calculate only Tropospheric values
    # if False, all layers
    # tropopause is defiend as o3>150ppb. If neede, change accordingly.
    Tropospheric=True  

    # if True, calculate only Tropospheric values
    # if False, all layers
    # tropopause is defiend as o3>150ppb. If neede, change accordingly.
    regional=False
    #dir_shapefile="/Users/roozitalab/INDIA/Shapefile/Bangladesh/bgd_adm_bbs_20201113_shp/bgd_adm_bbs_20201113_SHP/"
    dir_shapefile="/Users/roozitalab/INDIA/Shapefile/Countries/world"

    limit=(20,20,40,120)


    # choose the period of interest. Plots will be averaged within this period
    start_date = "2016-1-1"
    end_date = "2018-1-1"




    # Dictionary of model variables to be used. 
    # User can use this to include a combination of different variables in the calculation.
    # For example, for precipitation you can define PRECT as:
    #      dic_SE['PRECT']={'PRECC'+ext1_SE:8.64e7,'PRECL'+ext1_SE:8.64e7}
    #      - It means to sum the file variables "PRECC" and "PRECL" 
    #        for this arbitrary desired variable named "PRECT"
                        
    #      - It also has the option to apply conversion factors. 
    #        For instance, PRECL and PRECC are in m/s. 8.64e7 is used to convert m/s to mm/day
    dic_SE={}



    dic_SE['O3']={'O3'+ext1_SE:1e9} # covert to ppb for Tropopause calculation
    dic_SE['CH3CCL3']={'CH3CCL3'+ext1_SE:1e9} # covert to ppb for Tropopause calculation

    dic_SE['CH4']={'CH4'+ext1_SE:1}
    dic_SE['CO']={'CO'+ext1_SE:1}

    dic_SE['ISOP']={'ISOP'+ext1_SE:1}
    # dic_SE['C10H16']={'MTERP'+ext1_SE:1}
    # dic_SE['CH3OH']={'CH3OH'+ext1_SE:1}
    # dic_SE['CH3COCH3']={'CH3COCH3'+ext1_SE:1}
    # dic_SE['CH3OH']={'CH3OH'+ext1_SE:1}


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






    #Dictionary for Molecular weights. Keys must be consistent with variable name
    MW={'O3':48,
        'CH4':16,
        'CO':28,
        'ISOP':68,
        'C10H16':136,
        'CH3CCL3':133.4,    
        'SOA':144.132,
        'SALT':32.066,
        'SULF':115.11,
        'POM':12.011,
        'BC':12.011 ,
        'DUST':12.011}

    # MW={'O3':48,
    #     'CH4':16,
    #     'CO':28,
    #     'ISOP':68,
    #     'C10H16':136,
    #     'CH3CCL3':133.4042,    
    #     'SOA':144.132,
    #     'SALT':12.011,
    #     'SULF':115.11,
    #     'POM':12.011,
    #     'BC':12.011 ,
    #     'DUST':12.011}



    # Avogadro's Number
    AVO=6.022e23
    # gravity 
    gr=9.80616
    # Mw air
    #Mwair=28.97


    # Collect all the filenames and coordinates
    Files,Lats,Lons,areas= Get_files(data_dirs,scenarios,start_date,end_date,area=True)





    # find the name of all the variables in the file.
    # this will help the code to work for the variables that are not in the files (assingn 0s)
    tmp_file=Dataset(data_dirs[0]+Files[scenarios[0]][0])
    ListVars=tmp_file.variables.keys()
    tmp_file.close()




    AEROSOLS=['SOA','SALT','DUST','POM','BC','SULF']

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
                dic_SE[var+'_CHML'][key+'_Loss'+ext1_SE]=MW[var]*1e3/AVO/gr
                dic_SE[var+'_CHMP'][key+'_Prod'+ext1_SE]=MW[var]*1e3/AVO/gr 
                
            # if key=='khar':
            #     print('key')
                
            else:

                if key+'_CHML' in ListVars:            
                    dic_SE[var+'_CHML'][key+'_CHML'+ext1_SE]=MW[var]*1e3/AVO/gr        
                else:
                    dic_SE[var+'_CHML']['O3'+ext1_SE]=0.

                if key+'_CHMP' in ListVars:            
                    dic_SE[var+'_CHMP'][key+'_CHMP'+ext1_SE]=MW[var]*1e3/AVO/gr        
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
                dic_SE[var+'_CLXF'][key+'_CLXF'+ext1_SE]=MW[var]*10/AVO  # convert [molec/cm2/s] to [kg/m2/s]        
            else:
                dic_SE[var+'_CLXF']['CO_CLXF'+ext1_SE]=0. 


                        
                    
                    
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

    # extract all the data
    Dic_scn_var_comp={}

    # this is for finding tropospheric values
    Dic_crit={}

    for i in range(len(scenarios)):
        
        current_dir=data_dirs[i]
        current_scn=scenarios[i]  
        current_files=Files[current_scn] 


        Dic_scn_var_comp[current_scn]={}
        
        
        

        Dic_var_comp={}
        for v in range(len(variables)):
            current_var=variables[v]
    

            if current_var in AEROSOLS:

                # Components are: burden, chemical loss, chemical prod, dry deposition, 
                #                 surface emissions, elevated emissions, wet deposition, gas-aerosol exchange
                components=[current_var+'_BURDEN',current_var+'_CHML',current_var+'_CHMP',
                            current_var+'_DDF',current_var+'_WDF', current_var+'_SF', current_var+'_CLXF',
                            current_var+'_GAEX',current_var+'_DDFC',current_var+'_WDFC']   
                
                if current_var=='SULF':
                    # For SULF we also have AQS and NUCLEATION
                    components=[current_var+'_BURDEN',current_var+'_CHML',current_var+'_CHMP',
                                current_var+'_DDF',current_var+'_WDF', current_var+'_SF', current_var+'_CLXF',
                                current_var+'_GAEX',current_var+'_DDFC',current_var+'_WDFC',current_var+'_AQS',
                            current_var+'_NUCL']                
                

            else:

                # Components are: burden, chemical loss, chemical prod, dry deposition, 
                #                 surface emissions, elevated emissions, chemical tendency
                # I always add Lightning NOx production for gaseous species.
                components=[current_var+'_BURDEN',current_var+'_CHML',current_var+'_CHMP',
                            current_var+'_DDF',current_var+'_WDF', current_var+'_SF', current_var+'_CLXF',
                            current_var+'_TEND',current_var+'_LNO']         
            
        

    
            Dic_comp={}
            for comp in components:
                print(comp)
                
                
        
                    
                current_data=SEbudget(dic_SE,current_dir,current_files,comp,level=50)
                print(np.shape(current_data))
        
                    
                Dic_comp[comp]=current_data
    
            
            Dic_var_comp[current_var]=Dic_comp
            
        
        Dic_scn_var_comp[current_scn]= Dic_var_comp    
        
        
        
        current_crit=SEbudget(dic_SE,current_dir,current_files,'O3',level=50)  
        Dic_crit[current_scn]=current_crit


    # convert date strings to datetime format
    start_period = datetime.strptime(start_date, "%Y-%m-%d")
    end_period = datetime.strptime(end_date, "%Y-%m-%d")


    # in seconds
    duration=(end_period-start_period).days*86400


    # Tropospheric values: 
    for i in range(len(scenarios)):

        
        
        current_scn=scenarios[i]      
        area=areas[current_scn]

        current_lat=Lats[current_scn]
        current_lon=Lons[current_scn]

        if regional:
            #inside=Inside_SE_region(current_lat,current_lon,dir_shapefile)
            inside=Inside_SE(current_lat,current_lon,limit)
        else:
            if len(np.shape(area)) ==1:
                inside=np.full((len(current_lon)),True)
            else:
                inside=np.full((len(current_lat),len(current_lon)),True)
                
        
        print('Current Scenario: '+current_scn)
        print('*********\n')

        current_crit=Dic_crit[current_scn]
    
        if Tropospheric:
        
            trop=np.where(current_crit>150,np.nan,current_crit)
            strat=np.where(current_crit>150,current_crit,np.nan)
        else:
            trop=current_crit
        


        for v in range(len(variables)):
            current_var=variables[v]
            print(current_var)


            # Burden        
            spc_burd=Dic_scn_var_comp[current_scn][current_var][current_var+'_BURDEN']       
            spc_burd=np.where(np.isnan(trop),np.nan,spc_burd)
            tmp_burden=np.nansum(spc_burd*area,axis=0)
            burden=np.ma.masked_where(inside==False,tmp_burden)  #convert Kg/m2 to Tg
            BURDEN = np.ma.sum(burden*1e-9)

            # Chemical Loss
            spc_chml=Dic_scn_var_comp[current_scn][current_var][current_var+'_CHML'] 
            spc_chml=np.where(np.isnan(trop),np.nan,spc_chml)       
            tmp_chml=np.nansum(spc_chml*area,axis=0)
            chml=np.ma.masked_where(inside==False,tmp_chml)  #convert Kg/m2/s to Tg/yr
            CHML = np.ma.sum(chml*duration*1e-9)

            # Chemical Production
            spc_chmp=Dic_scn_var_comp[current_scn][current_var][current_var+'_CHMP'] 
            spc_chmp=np.where(np.isnan(trop),np.nan,spc_chmp)
            tmp_chmp=np.nansum(spc_chmp*area,axis=0)
            chmp=np.ma.masked_where(inside==False,tmp_chmp)  #convert Kg/m2/s to Tg/yr
            CHMP = np.ma.sum(chmp*duration*1e-9)
            
            # Surface Emissions
            spc_sf=Dic_scn_var_comp[current_scn][current_var][current_var+'_SF'] 
            tmp_sf=spc_sf
            sf=np.ma.masked_where(inside==False,tmp_sf*area)  #convert Kg/m2/s to Tg/yr
            SF = np.ma.sum(sf*duration*1e-9)

            
            # Elevated Emissions
            spc_clxf=Dic_scn_var_comp[current_scn][current_var][current_var+'_CLXF'] 
            print(np.shape(spc_clxf))

            tmp_clxf=spc_clxf
            clxf=np.ma.masked_where(inside==False,tmp_clxf)  #convert Kg/m2/s to Tg/yr
            CLXF = np.ma.sum(clxf*duration*1e-9)
            
            if current_var in AEROSOLS:

                # Dry Deposition Flux      
                spc_ddfa=Dic_scn_var_comp[current_scn][current_var][current_var+'_DDF'] 
                spc_ddfc=Dic_scn_var_comp[current_scn][current_var][current_var+'_DDFC']
                spc_ddf=spc_ddfa +spc_ddfc
                tmp_ddf=spc_ddf
                ddf=np.ma.masked_where(inside==False,tmp_ddf*area)  #convert Kg/m2/s to Tg/yr
                DDF = np.ma.sum(ddf*duration*1e-9)
                
                # Wet deposition
                spc_wdfa=Dic_scn_var_comp[current_scn][current_var][current_var+'_WDF'] 
                spc_wdfc=Dic_scn_var_comp[current_scn][current_var][current_var+'_WDFC']
                spc_wdf=spc_wdfa +spc_wdfc            
                tmp_wdf=spc_wdf
                wdf=np.ma.masked_where(inside==False,tmp_wdf*area)  #convert Kg/m2/s to Tg/yr
                WDF = np.ma.sum(wdf*duration*1e-9)
                
                
                # gas-aerosol Exchange
                spc_gaex=Dic_scn_var_comp[current_scn][current_var][current_var+'_GAEX'] 
                tmp_gaex=spc_gaex
                gaex=np.ma.masked_where(inside==False,tmp_gaex*area)  #convert Kg/m2/s to Tg/yr
                GAEX = np.ma.sum(gaex*duration*1e-9)
                
                
                # LifeTime = Burden/(loss+deposition)
                LT=BURDEN/(CHML+DDF-WDF)* duration/86400 # days 
                
                
                if current_var=='SULF':
                    
                    # Aqueous Chemistry
                    spc_aqs=Dic_scn_var_comp[current_scn][current_var][current_var+'_AQS'] 
                    tmp_aqs=spc_aqs
                    aqs=np.ma.masked_where(inside==False,tmp_aqs*area)  #convert Kg/m2/s to Tg/yr
                    AQS = np.ma.sum(aqs*duration*1e-9)
                
                
                    # Nucleation
                    spc_nucl=Dic_scn_var_comp[current_scn][current_var][current_var+'_NUCL'] 
                    tmp_nucl=spc_nucl
                    nucl=np.ma.masked_where(inside==False,tmp_nucl*area)  #convert Kg/m2/s to Tg/yr
                    NUCL = np.ma.sum(nucl*duration*1e-9)            
            
                
            else:

                # Dry Deposition Flux      
                spc_ddf=Dic_scn_var_comp[current_scn][current_var][current_var+'_DDF'] 
                tmp_ddf=spc_ddf
                ddf=np.ma.masked_where(inside==False,tmp_ddf*area)  #convert Kg/m2/s to Tg/yr
                DDF = np.ma.sum(ddf*duration*1e-9)

                
                # Wet Deposition Flux      
                spc_wdf=Dic_scn_var_comp[current_scn][current_var][current_var+'_WDF'] 
                tmp_wdf=spc_wdf
                wdf=np.ma.masked_where(inside==False,tmp_wdf*area)  #convert Kg/m2/s to Tg/yr
                WDF = np.ma.sum(wdf*duration*1e-9)
                
                
                # Chemical Tendency
                spc_tnd=Dic_scn_var_comp[current_scn][current_var][current_var+'_TEND'] 
                spc_tnd=np.where(np.isnan(trop),np.nan,spc_tnd)
                tmp_tnd=np.nansum(spc_tnd,axis=0)
                tnd=np.ma.masked_where(inside==False,tmp_tnd)  #convert Kg/s to Tg/yr
                TND = np.ma.sum(tnd*duration*1e-9)
                
                # Stratospheric-Tropospheric Exchange
                STE=DDF-TND

                # LifeTime = Burden/(loss+deposition)
                if current_var=='CO':
                    LT=BURDEN/(CHML)*duration/86400/365 # days
                    
                else:
                    LT=BURDEN/(CHML+DDF-WDF)*duration/86400/365 # days

                # Lightning NOX production
                spc_lno=Dic_scn_var_comp[current_scn][current_var][current_var+'_LNO']
                tmp_lno=np.ma.masked_where(inside==False,spc_lno)  
                LNO = np.ma.sum(tmp_lno)              
            row_values = []
            if current_var in AEROSOLS:
                
                
                print('Current Variable: '+current_var)
                print('Global Burden (Tg): '+str(np.round(BURDEN,3)))
                row_values.append(np.round(BURDEN,3))
                print('Global Chemical Loss (Tg/yr): '+str(np.round(CHML,2)))
                row_values.append(np.round(CHML,3))
                print('Global Chemical Prod (Tg/yr): '+str(np.round(CHMP,2)))
                row_values.append(np.round(CHMP,3))
                print('Global Chemical NET (Tg/yr): '+str(np.round(CHMP-CHML,2)))
                row_values.append(np.round(CHMP-CHML,3))        
                print('Global Dry Deposition (Tg/yr): '+str(np.round(DDF,2)))
                row_values.append(np.round(DDF,3))
                print('Global Wet Deposition (Tg/yr): '+str(np.round(WDF,2)))
                row_values.append(np.round(WDF,3))
                print('Global Surface Emis (Tg/yr): '+str(np.round(SF,2)))
                row_values.append(np.round(SF,3))
                print('Global Elevated Emis (Tg/yr): '+str(np.round(CLXF,2)))
                row_values.append(np.round(CLXF,3))
                print('Global Gas-Aerosol Exch (Tg/yr): '+str(np.round(GAEX,2)))
                row_values.append(np.round(GAEX,3))
                print('LifeTime (day): '+str(np.round(LT,0)))
                row_values.append(np.round(LT,0))

                if current_var=='SULF':
                    print('Global AQUEOUS Chem (Tg/yr): '+str(np.round(AQS,2)))
                    row_values.append(np.round(AQS,0))
                    print('Global Nucleation (Tg/yr): '+str(np.round(NUCL,2)))
                    row_values.append(np.round(NUCL,0))
                    
                
                print('****   *****')
                #Create output file name:
                output_csv_file = output_location / f"amwg_aerosol_table_{case_names[0]}.csv"
                cols = ['variable']+[f"Test {i+1}" for i,_ in enumerate(case_names[0:-1])]
                #cols = ['variable']+[f"Test {i+1}" for i,_ in enumerate(case_names[0:-1])]+["Baseline"]
                #row_values.append(np.round(my_val,3))
                dfentries = {c:[row_values[idx]] for idx,c in enumerate(cols)}
                df = pd.DataFrame(dfentries,columns=cols)
                df.to_csv(output_csv_file, header=False, index=False)
                
            else:
                print('Current Variable: '+current_var)
                print('Global Burden (Tg): '+str(np.round(BURDEN,2)))
                print('Global Chemical Loss (Tg/yr): '+str(np.round(CHML,2)))
                print('Global Chemical Prod (Tg/yr): '+str(np.round(CHMP,2))) 
                print('Global Chemical NET (Tg/yr): '+str(np.round(CHMP-CHML,2)))            
                print('Global Dry Deposition (Tg/yr): '+str(np.round(DDF,2)))
                print('Global Wet Deposition (Tg/yr): '+str(np.round(WDF,2)))
                print('Global Surface Emis (Tg/yr): '+str(np.round(SF,2)))
                print('Global Elevated Emis (Tg/yr): '+str(np.round(CLXF,2)))
                print('Global TND (Tg/yr): '+str(np.round(TND,2)))
                print('Global STE (Tg/yr): '+str(np.round(STE,2)))
                print('LifeTime (day): '+str(np.round(LT,2)))
                print('Global Lightning NO (Tg N/yr): '+str(np.round(LNO,2)))

                print('****   *****')        
        
        
        





# Helper fucntions
'''
SE_functions.py
this code is designed for compiling the functions used for processing SE files


MODIFICATION HISTORY:
    Behrooz Roozitalab, 02, NOV, 2022: VERSION 1.00
    - Initial version
   
'''

# import matplotlib as mpl

# from matplotlib.collections import PolyCollection
# from scipy.interpolate import griddata
# import cartopy.crs as ccrs
# import cartopy
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# import matplotlib.colors as colors
# from mpl_toolkits.basemap import Basemap
# from shapely.geometry import Polygon,Point

# from cartopy.io.shapereader import Reader
# from cartopy.feature import ShapelyFeature

def list_files(directory,scenario,start_period,end_period):

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
    all_filenames =list (file for file in os.listdir(directory) 
         if os.path.isfile(os.path.join(directory, file)))
    if len(all_filenames)==0 : sys.exit(" Directory has no outputs ")
    all_filenames.sort()


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
            
            # We need to only extract the files that are within the chosen dates.
            # first timestep is used for this purpose
            print(time_bounds[0,0])
            #filetime0=datetime.utcfromtimestamp(time_bounds[0,0].tolist()/1e9) # beginning time of first timestep
            #filetime1=datetime.utcfromtimestamp(time_bounds[0,1].tolist()/1e9) # ending time of first timestep


            # For Jun Zhang data
            filetime0=np.datetime64(time_bounds[0,0]) # beginning time of first timestep
            filetime1=np.datetime64(time_bounds[0,1]) # ending time of first timestep

            """if '.h0' in scenario: # this is hard coded. User should change it (e.g. to ".h1") accordingly to reflect monthly files.
                if  (start_period<=filetime0<end_period) :
                    print ('list_files_SE Warning: "h0" is hard-coded to contain monthly files. If not, change it in the function.') 
                    all_fileNames.append(all_filenames[i])
                
            else:
                
                if  (start_period<=filetime0<end_period) or (start_period<=filetime1<end_period):
     
                    all_fileNames.append(all_filenames[i])"""
            if '.h0' in scenario:
                all_fileNames.append(all_filenames[i])
                    

    return all_fileNames



def Get_files(data_dirs,scenarios,start_date,end_date,**kwargs):
        
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
    area=kwargs.pop('area',False)
    
    
    # convert date strings to datetime format
    start_period = datetime.strptime(start_date, "%Y-%m-%d")
    end_period = datetime.strptime(end_date, "%Y-%m-%d")



    files={}
    Lats={}
    Lons={}

    areas={}
    Earth_rad=6.371e6 # Earth Radius in m 



    
    for i in range(len(scenarios)):

        current_dir=data_dirs[i]
        current_scn=scenarios[i]


        # find the needed the files
        current_files=list_files(current_dir,current_scn,start_period,end_period)

        # get the Lat and Lons for each scenario
        tmp_file=xr.open_dataset(current_dir+current_files[0])
        lon=tmp_file['lon'+ext1_SE].data
        lon[lon > 180.] -= 360 # shift longitude from 0-360˚ to -180-180˚
        lat=tmp_file['lat'+ext1_SE].data

        if area==True:
            try:
                tmp_area=tmp_file['area'+ext1_SE].data
                Earth_area= 4 * np.pi * Earth_rad**(2)

                areas[current_scn]=tmp_area*Earth_area/np.nansum(tmp_area)


            except KeyError:
                dlon= np.abs(lon[1]-lon[0])
                dlat= np.abs(lat[1]-lat[0])

                lon2d,lat2d=np.meshgrid(lon,lat)
                #area=np.zeros_like(lat2d)

                dy=Earth_rad*dlat*np.pi/180
                dx=Earth_rad*np.cos(lat2d*np.pi/180)*dlon*np.pi/180

                area=dx*dy
                areas[current_scn]=area


        files[current_scn]=current_files
        Lats[current_scn]=lat
        Lons[current_scn]=lon



    return files,Lats,Lons,areas



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



    # gas constanct
    Rgas=287.04 #[J/K/Kg]=8.314/0.028965
        
    all_data=[]
    for file in range(len(files)):
        ds=xr.open_dataset(data_dir+files[file])


        data=[]
        
    
        for i in dic_SE[var].keys():
            data.append(np.array(ds[i].isel(time=0))*dic_SE[var][i])
  
        data=np.sum(data,axis=0)
  
    
            
        if ('CHML' in var) or ('CHMP' in var) : 

            
            Temp=np.array(ds['T'].isel(time=0))
            Pres=np.array(ds['PMID'].isel(time=0))
            rho= Pres/(Rgas*Temp)
            
            
            delP=np.array(ds['PDELDRY'].isel(time=0))
#             hyai=np.array(ds['hyai'])
#             hybi=np.array(ds['hybi'])
#             PSD=np.array(ds['PSDRY'][0])
#             P0=1e+5

#             PSI=np.zeros((len(hyai),len(PSD)))

#             for i in range(len(hyai)):           
#                 PSI[i]=hyai[i]*P0 + hybi[i]*PSD            

            
#             delP=PSI[1:]-PSI[:-1]
            
                        
            data=data*delP/rho
        elif ('BURDEN' in var):
        
            delP=np.array(ds['PDELDRY'].isel(time=0))
            
                        
            data=data*delP

        else:
            data=data
            

        
        all_data.append(data)



    all_data=np.nanmean(all_data,axis=0)
        

            
    return all_data



def Inside_SE_region(lat_center,lon_center,dir_shapefile):

    bsmap=Basemap(llcrnrlon=-180.,llcrnrlat=-90.,urcrnrlon=180.\
                        ,urcrnrlat=90.,projection='cyl')
    
    #dir_shapefile="/Users/roozitalab/INDIA/Shapefile/india_administrative_outline_boundary/india_administrative_outline_boundary"
 
    #bsmap.readshapefile(dir_shapefile+"bgd_admbnda_adm1_bbs_20201113",'states',drawbounds=True,linewidth=0.5)
    bsmap.readshapefile(dir_shapefile,'states',drawbounds=True,linewidth=0.5)

    inside=np.full((len(lon_center)),False)

    states=[]
    info=bsmap.states_info
    for i in range(len(info)):
        #if info[i]['ST_NM']+'_'+info[i][shape_level] not in states:       # ST_NM for states
            #if info[i]['ST_NM']=='NCT of Delhi':
            states.append(info[i]['CNTRY_NAME']+'_'+str(info[i]['RINGNUM']))   
    
    states_unique=[]
    for i in range(len(info)):
        if info[i]['CNTRY_NAME'] not in states_unique:       # ST_NM for states
            #if info[i]['ST_NM']=='NCT of Delhi':
            states_unique.append(info[i]['CNTRY_NAME'])

    for item in range(len(states)):
    

        
             if 'India' in states[item]:  
     
                 polygon=State_idx(states[item],bsmap,'CNTRY_NAME')
                 if len(polygon)>0:
                     #print(len(polygon))
                     
              
                     for i in range(len(lat_center)):
                         if (40.<lon_center[i]<105.) & (0.<lat_center[i]<55.) :
                         
                         
                             point=Point(lon_center[i],lat_center[i])
                             if point.within(Polygon(polygon[0])):
                                 inside[i]=True
    return inside


 


def Inside_SE(lat_center,lon_center,limit):
        
        inside=np.full((len(lon_center)),False)
        
        a=np.where(((limit[1]<lon_center) & (lon_center<limit[3]) & (limit[0]<lat_center) & (lat_center<limit[2])))
        inside[a]=True
        
        return inside
    
def State_idx(state,map,shape_level):

    #i=0
    poly=[]
    for info, shp in zip(map.states_info, map.states):
        proid = info['CNTRY_NAME']+'_'+ str(info['RINGNUM'])     # ST_NM for states and DITSRICT for counties
        if proid == state:
            #print(i)
            #for i in range(len(shp)):
                
            poly.append(shp)
            break
            #i+=1
    return poly    