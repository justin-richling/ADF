
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

from datetime import datetime
from dateutil.relativedelta import relativedelta
import pandas as pd


#Set seasonal ranges:
seasons = {"DJF": [12, 1, 2],
            "JJA": [6, 7, 8],
            "MAM": [3, 4, 5],
            "SON": [9, 10, 11]}

import xarray as xr
from pathlib import Path
import glob
from itertools import chain


#
# --- Main Function Shares Name with Module: regional_map_multicase ---
#
def make_scycle_maps(adfobj, diag, data_dict, case_deets):
    """
    Chemistry Map main function
        * Initially start with Aerosol Zonal maps
            - This probably can be expanded to LatLon if given single pressure levels?

        
    """

    # Notify user that script has started:
    print("\n  Generating zonal aerosol plots ...")

    var_list = adfobj.diag_var_list

    #Aerosol Calculations
    if diag == "aerosol":

        for var,constits in aerosol_dict.items():
            if all(elem in var_list  for elem in constits):
                print(f"\t - zonal mean aerosol plots for {var}")
                        
                #If found then notify user, assuming debug log is enabled:
                adfobj.debug_log(f"zonal_mean: Found variable defaults for {var}")
                aerosol_plot(adfobj, var, data_dict, case_deets)

            else:
                print(f"No constituents for {var}, moving on ...")




merra2 = {}
merra2_seasonal = {}
merra2_monthly = {}
merra2_vars = ['U','T','V','lat','lev']
#ncfile = nc.Dataset("MERRA2_met.nc", 'r')
merra_ncfile = xr.open_dataset("MERRA2_met.nc", decode_times=True, use_cftime=True)

            

for index, var in enumerate(merra2_vars):
    
    merra2[var] = merra_ncfile[var]
    if index < len(merra2_vars)-2:


        #start_date = datetime(merra_ncfile.time) #datetime(1980, 1, 1)
        time0 = merra_ncfile.time.dt.strftime("%Y-%m-%d")
        start_time = time0.values[0]
        start_date = datetime.strptime(str(start_time), '%Y-%m-%d')

        # Number of months to generate
        # Number of months to generate
        num_months = len(merra_ncfile[var].record)#record
                    
        # List to store datetime objects
        datetime_list = []
                    
        # Generate datetime objects incrementally by month
        for i in range(num_months):
            new_date = start_date + relativedelta(months=i)
            datetime_list.append(new_date.replace(day=1))  # Set the day to the first day of the month
    
        merra_ncfile[var] = merra_ncfile[var].assign_coords({"record": datetime_list})
        merra2[var] = merra_ncfile[var]

        
        #for season in seasons:
            #merra2_seasonal[var] = averaging_functions.calc_seasonal_average(merra2[var])
        merra2_seasonal[var] = obs_seasonal_mean(merra2[var], season="DJF", is_climo=None)
        #merra2_monthly[var] = averaging_functions.calc_monthly_average(merra2[var])
        #merra2_monthly[var] = averaging_functions.calc_monthly_average(merra2[var])
        #merra2_monthly[var] = monthly_mean(merra2[var], month=11, is_climo=None)












def zonal_wind(data_dict, obs="MERRA2"):
    # CAM vars needed
    # - U

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 
    season="MAM"


    #Do 001 - 001 diff, etc, in addition to direct MERRA2 comparison

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 
    for run in range(len(runs)):
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],
                                        data_runs[runs[run]]['lat'])
        
        contour_levels = np.arange(-120, 121, 10)
        #ax = fig.add_subplot(2, 4, run+1)
        ax = fig.add_subplot(2, len(runs), run+1)
        
        cf=plt.contourf(lev_grid, lat_grid,
                        data_runs_seasonal[runs[run]]['U'][season].transpose(transpose_coords=True),
                        levels=contour_levels)
        
        #cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_seasonal[runs[run]]['U'][season,:,:]), levels=contour_levels)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        plt.title(run_names[run])
        ax.set_xticklabels([])
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='m/s')
    run=0
    for run in range(len(runs)):
        #reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],
        #                                                       merra2_seasonal['U'],
        #                                                       #merra2_seasonal['U'][season],
        #                                                       data_runs[runs[run]]['lat'],
        #                                                       data_runs[runs[run]]['lev'])

        reference_field = merra2_seasonal['U'].interp(lat=data_runs[runs[run]]['lat'],
                                                    lev=data_runs[runs[run]]['lev'], method='linear')

        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(-30, 31, 3)
        ax = fig.add_subplot(2, len(runs), len(runs)+run+1)
        #ax = fig.add_subplot(2, 4, 4+run+1)
        
        cf=plt.contourf(lev_grid, lat_grid,
                        (data_runs_seasonal[runs[run]]['U'][season]-reference_field).transpose(transpose_coords=True),
                        levels=contour_levels)
        
        plt.set_cmap('RdBu_r')
        orig_contour_levels = np.arange(-120, 121, 10)
        plt.contour(lev_grid, lat_grid, reference_field.transpose(transpose_coords=True),
                    colors='black', levels=orig_contour_levels,
                    negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        
        plt.yscale("log")
        ax.set_ylim(1000,0.1)
        plt.yticks([1000,100,10,1,0.1])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        local_title = run_names[run] + " - MERRA2"
        plt.title(local_title) 
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", ticks=np.arange(-30,31,10), label='m/s')
    fig.suptitle(f"Xarray Created Zonal Wind\nXarray package read in file\nADF Seasonal Means - {season}")
    plt.savefig('output/zonal_wind_mam_ADF_seasonal_xarray_zm_xarray_read.png',dpi=300)




def zonal_temp():
    # CAM vars needed
    # - T

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 
    season=2;

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 
    for run in range(len(runs)):
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(160, 301, 10)
        ax = fig.add_subplot(2, 4, run+1)
        cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_seasonal[runs[run]]['T'][season,:,:]), levels=contour_levels)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        plt.title(run_names[run])
        ax.set_xticklabels([])
        if run == 3:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')
    run=0
    #reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_seasonal['T'][season,:,:],data_runs[runs[0]]['lat'],data_runs[runs[0]]['lev'])
    for run in range(len(runs)):
        #reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_seasonal['U'][season,:,:],data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
        reference_field=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_seasonal['T'][season,:,:],data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
        
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(-20, 21, 2)
        ax = fig.add_subplot(2, 4, 4+run+1)
        cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_seasonal[runs[run]]['T'][season,:,:]-reference_field), levels=contour_levels)
        plt.set_cmap('RdBu_r')
        orig_contour_levels = np.arange(150, 300, 10)
        plt.contour(lev_grid, lat_grid, np.transpose(reference_field), colors='black', levels=orig_contour_levels, negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        local_title = run_names[run] + " - MERRA2"
        plt.title(local_title) 
        if run == 3:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", ticks=np.arange(-20,21,4), label='K')
    plt.savefig('output/temp_jja_mine.png',dpi=300)



def temp_vs_saber():
    # CAM vars needed
    # - T

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 
    month=11

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 
    for run in range(len(runs)):
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(140, 301, 10)
        ax = fig.add_subplot(2, len(runs), run+1)
        cf=plt.contourf(lev_grid, lat_grid, data_runs_monthly[runs[run]]['T'].transpose(transpose_coords=True), levels=contour_levels)
        #cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_monthly[runs[run]]['T'][month,:,:]), levels=contour_levels)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        plt.title(run_names[run])
        ax.set_xticklabels([])
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')
    run=0

    for run in range(len(runs)):
        #reference_field=interpolating_functions.interpolate_2d(saber['latitude'],saber['pressure'],saber_monthly['temp'][month,:,:],data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
        reference_field = saber_monthly['temp'].interp(latitude=data_runs[runs[run]]['lat'], pressure=data_runs[runs[run]]['lev'], method='linear')
        [lat_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],data_runs[runs[run]]['lat'])
        contour_levels = np.arange(-40, 41, 4)
        ax = fig.add_subplot(2, len(runs), len(runs)+run+1)
        cf=plt.contourf(lev_grid, lat_grid, (data_runs_monthly[runs[run]]['T']-reference_field).transpose(transpose_coords=True), levels=contour_levels)
        #cf=plt.contourf(lev_grid, lat_grid, np.transpose(data_runs_monthly[runs[run]]['T'][month,:,:]-reference_field), levels=contour_levels)
        plt.set_cmap('RdBu_r')
        orig_contour_levels = np.arange(140, 300, 10)
        ctours=plt.contour(lev_grid, lat_grid, reference_field.transpose(transpose_coords=True), colors='black', levels=orig_contour_levels, negative_linestyles='dashed', linewidths=.5, alpha=0.5)
        plt.yscale("log")
        ax.set_ylim(1000,0.0001)
        plt.yticks([1000,100,10,1,0.1,.01,.001,.0001])
        plt.xticks(np.arange(-90,91,45),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        local_title = run_names[run] + " - SABER"
        plt.title(local_title) 
        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-1.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')
    plt.savefig('output/temp_dec_saber_ADF_monthly_xarray_zm_xarray_read.png',dpi=300)

def qbo():
    # CAM vars needed
    # - U


    plt.rc('font', size=8) 
    run=0
    #Plot zonal wind data
    #fig, axes = plt.subplot_mosaic('AAAABC;DDDDEF;GGGGHI;JJJJKL;MMMMNO')
    #main_key=['A','D','G','J','M']
    #side1_key=['B','E','H','K','N']
    #side2_key=['C','F','I','L','O']

    #fig, axes = plt.subplot_mosaic('AAAABC;DDDDEF;GGGGHI;JJJJKL',figsize=(15,10))
    #main_key=['A','D','G','J']
    #side1_key=['B','E','H','K']
    #side2_key=['C','F','I','L']


    fig, axes = plt.subplot_mosaic('AAAABC;DDDDEF;GGGGHI',figsize=(15,10))
    main_key=['A','D','G']
    side1_key=['B','E','H']
    side2_key=['C','F','I']



    runs_local={}
    #runs_local[0]=fv
    runs_local[0]=data_runs[runs[0]]
    runs_local[1]=data_runs[runs[1]]
    #runs_local[2]=data_runs[runs[2]]

    y = 1.00

    contour_levels = np.arange(-35, 35, 2.5)

    # MERRA2 plot
    plot_num = len(runs)
    merra_plot = plot_num
    nt=len(data_runs[runs[run]]['date'])
    if nt > 120:
        nt = 120
    #else:
    #    nt = nt
    print("MERRA2",nt)

    #plotdata=averaging_functions.coslat_average(merra2['U'],merra2['lat'],-10,10)
    plotdata = cosweightlat(merra2['U'],-10,10)

    plotdata_clip = np.clip(np.abs(plotdata), None, 35)
    plotdata=np.sign(plotdata)*plotdata_clip
    [time_grid, lev_grid] = np.meshgrid(merra2['lev'],np.arange(1,nt+1,1))
    start_ind=252
    end_ind=start_ind+nt
    axes[main_key[merra_plot]].contourf(lev_grid, time_grid, plotdata[start_ind:end_ind,0,:], levels=contour_levels, cmap='RdBu_r')
    axes[main_key[merra_plot]].set_ylim(100,3)
    axes[main_key[merra_plot]].set_yscale("log")
    axes[main_key[merra_plot]].set_ylabel('hPa')
    axes[main_key[merra_plot]].set_title("MERRA2",y=y)
    #axes[main_key[plot_num]].set_xticks(np.arange(1,nt+1,12),rotation=40)
    #axes[main_key[plot_num]].set_xticks(np.arange(1,121,12),rotation=40)
    #axes[main_key[plot_num]].set_xticks(np.arange(int(merra_start_yr),int(merra_end_yr)+1,5),rotation=40)
    #axes[main_key[plot_num]].set_xticklabels(['2000','2001','2002','2003','2004','2005','2006','2007','2008','2009'])
    #axes[main_key[merra_plot]].set_xticklabels(np.arange(int(merra_start_yr),int(merra_end_yr)+1,5))

    axes[main_key[plot_num]].set_xticks(np.arange(1,nt+1,12),rotation=40)
    #axes[main_key[run]].set_xticklabels(np.arange(int(yrs[0]),int(yrs[0])+10,1))
    #axes[main_key[plot_num]].set_xticklabels(['2000','2001','2002','2003','2004','2005','2006','2007','2008','2009',"2010",
    #                                          '2011','2012','2013'])

    start_year = int(str(plotdata[252].record.values)[0:4])
    axes[main_key[plot_num]].set_xticklabels(np.arange(start_year,start_year+(nt/12),1).astype(int))

    amp_m = qbo_functions.qbo_amplitude(plotdata[:,0,:])
    axes[side1_key[merra_plot]].plot(amp_m,merra2['lev'],color='k',linewidth=1)
    axes[side1_key[merra_plot]].set_ylim(100,3)
    axes[side1_key[merra_plot]].set_yscale("log")
    axes[side1_key[merra_plot]].set_xlim(0,20)
    axes[side1_key[merra_plot]].set_xticks(np.arange(0,21,5))
    axes[side1_key[merra_plot]].set_xlabel('m/s')
    axes[side1_key[merra_plot]].set_yticks([])

    period_m = qbo_functions.qbo_frequency(plotdata[:,0,:])
    axes[side2_key[merra_plot]].plot(period_m,merra2['lev'],color='k',linewidth=1)
    axes[side2_key[merra_plot]].set_ylim(100,3)
    axes[side2_key[merra_plot]].set_yscale("log")
    axes[side2_key[merra_plot]].set_xlim(0,40)
    axes[side2_key[merra_plot]].set_xticks(np.arange(0,41,10))
    axes[side2_key[merra_plot]].set_yticks([])
    axes[side2_key[merra_plot]].set_xlabel('months')

    for run in range(len(runs)):
        yrs = wowsa[runs[run]]
        nt=len(runs_local[run]['date'])
        if nt > 120:
            nt_sub = 120
        else:
            nt_sub = nt
        print(run_names[run],nt,"\n")
        
        if run ==0:
            [time_grid, lev_grid] = np.meshgrid(runs_local[run]['lev'],np.arange(1,nt+1+12,1))
        else:
            [time_grid, lev_grid] = np.meshgrid(runs_local[run]['lev'],np.arange(1,nt+1,1))
        contour_levels = np.arange(-35, 35, 2.5)
        #ax = fig.add_subplot(4, 4, run,colspan=3)
        #plotdata=averaging_functions.coslat_average(runs_local[run]['U'],runs_local[run]['lat'],-10,10)
        plotdata = cosweightlat(runs_local[run]['U'],-10,10)
        plotdata_clip = np.clip(np.abs(plotdata), None, 35)
        plotdata=np.sign(plotdata)*plotdata_clip
        if run == 0:
            start_idx = 11
            axes[main_key[run]].contourf(lev_grid[start_idx:nt_sub+12,:], time_grid[start_idx:nt_sub+12,:], plotdata[start_idx:nt_sub+12,:], levels=contour_levels, cmap='RdBu_r')
        else:
            start_idx = 0
            axes[main_key[run]].contourf(lev_grid[start_idx:nt_sub,:], time_grid[start_idx:nt_sub,:], plotdata[start_idx:nt_sub,:], levels=contour_levels, cmap='RdBu_r')
        axes[main_key[run]].set_ylim(100,3)
        axes[main_key[run]].set_yscale("log")
        axes[main_key[run]].set_ylabel('hPa')
        axes[main_key[run]].set_title(run_names[run],y=y)
        if run == 0:
            axes[main_key[run]].set_xticks(np.arange(12,nt_sub+12,12),rotation=40)
        else:
            axes[main_key[run]].set_xticks(np.arange(1,nt_sub,12),rotation=40)
        #axes[main_key[run]].set_xticks([])
        #if nt > 120:
        if run == 0:
            axes[main_key[run]].set_xticklabels(np.arange(int(yrs[0]+1),int(yrs[0])+10+1,1))
        else:
            axes[main_key[run]].set_xticklabels(np.arange(int(yrs[0]),int(yrs[0])+10,1))
        
        amp = qbo_functions.qbo_amplitude(plotdata)
        axes[side1_key[run]].plot(amp,runs_local[run]['lev'],linewidth=1)
        axes[side1_key[run]].plot(amp_m,merra2['lev'],color='k',linewidth=1)
        axes[side1_key[run]].set_ylim(100,3)
        axes[side1_key[run]].set_yscale("log")
        axes[side1_key[run]].set_xlim(0,20)
        axes[side1_key[run]].set_xticks(np.arange(0,20,5))
        #axes[side1_key[run]].set_xticks([])
        axes[side1_key[run]].set_yticks([])
        if run==0:
            axes[side1_key[run]].set_title('Amplitude',y=y)
            
        
        period = qbo_functions.qbo_frequency(plotdata)
        axes[side2_key[run]].plot(period,runs_local[run]['lev'],linewidth=1)
        axes[side2_key[run]].plot(period_m,merra2['lev'],color='k',linewidth=1)
        axes[side2_key[run]].set_ylim(100,3)
        axes[side2_key[run]].set_yscale("log")
        axes[side2_key[run]].set_xlim(0,40)
        axes[side2_key[run]].set_xticks(np.arange(0,40,10))
        #axes[side2_key[run]].set_xticks([])
        axes[side2_key[run]].set_yticks([])
        if run==0:
            axes[side2_key[run]].set_title('Period',y=y)

    plt.savefig('output/qbo_xarray_zm_xarray_read.png',dpi=300)


def polar_cap_temp():
    # CAM vars needed
    # - T

    month=11

    #Pcap seasonal cycle
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica'] + plt.rcParams['font.sans-serif']
    plt.rc('font', size=8) 

    #Plot zonal wind data
    fig = plt.figure(figsize=(len(runs)*4,10)) 

    #reference_field_seas=np.zeros((12,len(data_runs[runs[0]]['lev']),len(data_runs[runs[0]]['lat'])))
    #for month in range(0,12):
    #    reference_field_seas[month,:,::]=interpolating_functions.interpolate_2d(merra2['lat'],merra2['lev'],merra2_monthly['T'][month,:,:],data_runs[runs[0]]['lat'],data_runs[runs[0]]['lev'])
    #merra2_pcap=averaging_functions.coslat_average(reference_field_seas,data_runs[runs[0]]['lat'],-90,-60)
        
    for run in range(len(runs)):
        reference_field_seas=np.zeros((12,len(data_runs[runs[run]]['lev']),len(data_runs[runs[run]]['lat'])))
        for month in range(0,12):
            reference_field_seas[month,:,::]=interpolating_functions.interpolate_2d(merra2['lat'],
                                                                                    merra2['lev'],
                                                                                    merra2_monthly['T'][month,:,:],
                                                                                    data_runs[runs[run]]['lat'],data_runs[runs[run]]['lev'])
            
        merra2_pcap=averaging_functions.coslat_average(reference_field_seas,data_runs[runs[run]]['lat'],-90,-60)

        [time_grid, lev_grid] = np.meshgrid(data_runs[runs[run]]['lev'],np.arange(0,12))
        contour_levels = np.arange(-10,11,1)
        #contour_levels = np.arange(-35,16,1)
        
        ax = fig.add_subplot(2, 4, run+1)
        
        run_pcap=averaging_functions.coslat_average(data_runs_monthly[runs[run]]['T'],data_runs[runs[run]]['lat'],-90,-60)

        #, levels=contour_levels
        cf=plt.contourf(lev_grid, time_grid, run_pcap-merra2_pcap,
                        levels=contour_levels
                    )
        plt.yscale("log")
        ax.set_ylim(300,1)
        ax.set_yticks([300,100,30,10])
        #ax.set_xticks(np.arange(0,12,2),rotation=40)
        ax.set_xticklabels(('Jan','Mar','May','Jul','Sep','Nov'),rotation=40)
        if run > 0:
            plt.yticks([]) 
        else:
            plt.ylabel('hPa')
        plt.set_cmap('RdBu_r')
        local_title=run_names[run]+" - MERRA2"
        plt.title(local_title)
        ctours=plt.contour(lev_grid, time_grid, merra2_pcap, colors='black', levels=orig_contour_levels, negative_linestyles='dashed', linewidths=.5, alpha=0.5)

        if run == len(runs)-1:
            axins = inset_axes(ax,
                        width="5%",  
                        height="80%",
                        loc='center right',
                        borderpad=-2.5
                    )
            fig.colorbar(cf, cax=axins, orientation="vertical", label='K')

    plt.savefig('output/temp_spcap_mine.png',dpi=300)