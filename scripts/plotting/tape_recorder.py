# Import necessary packages for the new script
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import xarray as xr
import pandas as pd

from dateutil.relativedelta import relativedelta
import glob
from pathlib import Path

def tape_recorder(adfobj):
    """
    Calculate the weighted latitude average for the simulations and 
    plot the values of Q against two sets of obseravations, MLS and ERA5, for the tropics
    between 10S and 10N.

    MLS h2o data is for 09/2004-11/2021
    ERA5 Q data is for 01/1980-12/2020

    Optional built-in colormaps:
      - blue2red
      - precip
      - precip_nowhite -> default cmap
      - red2blue

    NOTE: If the baseline case is observations, it will be ignored
        since a defualt set of obs are already being compared against in the tape recorder.
    """

    #CAM diagnostic plotting functions:
    import plotting_functions as pf

    #Notify user that script has started:
    print("\n  Generating tape recorder plots...")

    #Special ADF variable which contains the output paths for plots:
    plot_location = adfobj.plot_location
    plot_loc = Path(plot_location[0])

    #Grab test case name(s)
    case_names = adfobj.get_cam_info('cam_case_name', required=True)

    #Grab test case time series locs(s)
    case_ts_locs = adfobj.get_cam_info("cam_ts_loc", required=True)

    #Grab test case climo years
    start_years = adfobj.climo_yrs["syears"]
    end_years = adfobj.climo_yrs["eyears"]

    #Grab test case nickname(s)
    test_nicknames = adfobj.get_cam_info('case_nickname')
    if test_nicknames == None:
        test_nicknames = case_names

    # CAUTION:
    # "data" here refers to either obs or a baseline simulation,
    # Until those are both treated the same (via intake-esm or similar)
    # we will do a simple check and switch options as needed:
    if not adfobj.get_basic_info("compare_obs"):

        #Append all baseline objects to test case lists
        data_name = adfobj.get_baseline_info("cam_case_name", required=True)
        case_names = case_names + [data_name]
        
        data_ts_loc = adfobj.get_baseline_info("cam_ts_loc", required=True)
        case_ts_locs = case_ts_locs+[data_ts_loc]

        base_nickname = adfobj.get_baseline_info('case_nickname')
        if base_nickname == None:
            base_nickname = data_name
        test_nicknames = test_nicknames+[base_nickname]

        data_start_year = adfobj.climo_yrs["syear_baseline"]
        data_end_year = adfobj.climo_yrs["eyear_baseline"]
        start_years = start_years+[data_start_year]
        end_years = end_years+[data_end_year]
    #End if

    res = adfobj.variable_defaults # will be dict of variable-specific plot preferences
    # or an empty dictionary if use_defaults was not specified in the config YAML file.

    #Grab location of ADF default obs files
    adf_obs_loc = Path(adfobj.get_basic_info("obs_data_loc"))

    # Default colormap
    cmap='precip_nowhite'

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

    #This may have to change if other variables are desired in this plot type?
    plot_name = plot_loc / f"Q_TapeRecorder_ANN_WACCM_SeasonalCycle_Mean.{plot_type}"
    print(f"\t - Plotting annual tape recorder for Q")

    # Check redo_plot. If set to True: remove old plot, if it already exists:
    if (not redo_plot) and plot_name.is_file():
        #Add already-existing plot to website (if enabled):
        adfobj.debug_log(f"'{plot_name}' exists and clobber is false.")
        adfobj.add_website_data(plot_name, "Q_TapeRecorder", None, season="ANN", plot_type="WACCM", ext="SeasonalCycle_Mean",multi_case=True,category="Seasonal Cycle")
        return

    elif ((redo_plot) and plot_name.is_file()) or (not plot_name.is_file()):
        if plot_name.is_file():
            plot_name.unlink()
    
    #Make dictionary for case names and associated timeseries file locations
    runs_LT2={}
    for i,val in enumerate(test_nicknames):
        runs_LT2[val] = case_ts_locs[i]

    # MLS data
    #mls_filename = res['tape_recorder']['mls']['obs_file']
    mls_filename = "mls_h2o_latNpressNtime_3d_monthly_v5.nc"
    mls_file = adf_obs_loc / mls_filename
    #saber_filename = seas_cyc_res['saber_file']
    #saber_file = adf_obs_loc / saber_filename
    #merra_filename = seas_cyc_res['merra2_file']
    #merra_file = adf_obs_loc / merra_filename



    #mls_file = pf.check_obs_file(adfobj, Path(mls_file))

    mls = pf.load_dataset(str(mls_file))
    if mls:
        #mls = xr.open_dataset("/glade/campaign/cgd/cas/islas/CAM7validation/MLS/mls_h2o_latNpressNtime_3d_monthly_v5.nc")
        mls = mls.rename(x='lat', y='lev', t='time')
        time = pd.date_range("2004-09","2021-11",freq='M')
        mls['time'] = time
        mls = pf.coslat_average(mls.H2O,-10,10)
        mls = mls.groupby('time.month').mean('time')
        # Convert mixing ratio values from ppmv to kg/kg
        mls = mls*18.015280/(1e6*28.964)
    else:
        no_mls = True
        print("Incorrect MLS file/path provided, so MLS won't be plotted. ")
        print("Please check your location in the 'tape_recorder' section of the variable defaults yaml file.\n")


    # ERA5 data
    #era5_filename = res['tape_recorder']['era5']['obs_file']
    era5_filename = "ERA5_Q_10Sto10N_1980to2020.nc"
    era5_file = adf_obs_loc / era5_filename
    #era5_file = pf.check_obs_file(adfobj, Path(era5_file))


    era5 = pf.load_dataset([str(era5_file)])
    if era5:
        #era5 = xr.open_dataset("/glade/campaign/cgd/cas/islas/CAM7validation/ERA5/ERA5_Q_10Sto10N_1980to2020.nc")
        era5 = era5.groupby('time.month').mean('time')
    else:
        no_era5 = True
        print("Incorrect ERA5 file/path provided, so ERA5 won't be plotted. ")
        print("Please check your location in the 'tape_recorder' section of the variable defaults yaml file.\n")


    alldat=[]
    runname_LT=[]
    for idx,key in enumerate(runs_LT2):
        fils= sorted(Path(runs_LT2[key]).glob(f'*{adfobj.hist_str}*.Q.*.nc'))
        print(fils,"\n")
        print(len(fils))
        dat = pf.load_dataset(str(fils))

        #Check if data files exist, skip current case if not
        if not dat:
            errmsg = f"No files for '{key}'\n"
            errmsg += "Please make sure Q is in the CAM output"
            print(errmsg)
            continue

        dat = fixcesmtime(dat,start_years[idx],end_years[idx])
        datzm = dat.mean('lon')
        dat_tropics = pf.coslat_average(datzm.Q, -10, 10)
        dat_mon = dat_tropics.groupby('time.month').mean('time').load()
        alldat.append(dat_mon)
        runname_LT.append(key)

    runname_LT=xr.DataArray(runname_LT, dims='run', coords=[np.arange(0,len(runname_LT),1)], name='run')
    
    #Check if any CAM cases were made, if none, kill this script and have ADF continue on
    if len(alldat) == 0:
        print("No CAM cases for Q, so tape recorder plots will not be made. Moving on.")
        return

    #Combine all case data arrays
    alldat_concat_LT = xr.concat(alldat, dim=runname_LT)

    #Total number of plots
    case_num = alldat_concat_LT.run.size+2

    #Calculate the number of rows needed
    rows = (case_num + 5 - 1) // 5

    #Setup plots
    fig = plt.figure(figsize=(25,rows*16))
    x1, x2, y1, y2 = get5by5coords_zmplots()

    plot_step = 0.5e-7
    plot_min = 1.5e-6
    plot_max = 3e-6

    if not no_mls:
        ax = plot_pre_mon(fig, mls, plot_step,plot_min,plot_max,'MLS',
                        x1[0],x2[0],y1[0],y2[0],cmap=cmap, paxis='lev',
                        taxis='month',climo_yrs="2004-2021")

    if not no_era5:
        ax = plot_pre_mon(fig, era5.Q, plot_step,plot_min,plot_max,
                        'ERA5',x1[1],x2[1],y1[1],y2[1], cmap=cmap, paxis='pre',
                        taxis='month',climo_yrs="1980-2020")

    #Start count at 2 to account for MLS and ERA5 plots above
    count=2
    for irun in np.arange(0,alldat_concat_LT.run.size,1):
        title = f"{alldat_concat_LT.run.isel(run=irun).values}"
        ax = plot_pre_mon(fig, alldat_concat_LT.isel(run=irun),
                          plot_step, plot_min, plot_max, title,
                          x1[count],x2[count],y1[count],y2[count],cmap=cmap, paxis='lev',
                          taxis='month',climo_yrs=f"{start_years[irun]}-{end_years[irun]}")
        count=count+1
    
    #Shift colorbar if there are less than 5 subplots
    # There will always be at least 2 (MLS and ERA5)
    if len(case_ts_locs) == 1:
        x1_loc = (x1[1]-x1[0])/2
        x2_loc = ((x2[2]-x2[1])/2)+x2[1]
    elif len(case_ts_locs) == 2:
        x1_loc = (x1[1]-x1[0])/2
        x2_loc = ((x2[3]-x2[2])/2)+x2[2]
    else:
        x1_loc = x1[1]
        x2_loc = x2[3]

    y1_loc = y1[count]-0.03
    y2_loc = y1[count]-0.02

    ax = plotcolorbar(fig, plot_step, plot_min, plot_max, 'Q', #'Q (vmr)'
                      x1_loc, x2_loc, y1_loc, y2_loc,
                      cmap=cmap)

    #Save image
    fig.savefig(plot_name, bbox_inches='tight', facecolor='white')

    #Add plot to website (if enabled):
    adfobj.add_website_data(plot_name, "Q_TapeRecorder", None, season="ANN", plot_type="WACCM",
                            ext="SeasonalCycle_Mean", multi_case=True, category="Seasonal Cycle")

    #Notify user that script has ended:
    print("  ...Tape recorder plots have been generated successfully.")

    #End tape recorder plotting script:
    return


# Helper Functions
###################

def blue2red_cmap(n, nowhite = False):
    """
    combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals
    nowhite = choice of white separating the diverging colors or not
    """

    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    if (nowhite):
        nwhite=0

    colors1 = plt.cm.Blues_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.YlOrRd(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))

    colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap 

#########

def red2blue_cmap(n, nowhite = False):
    """ 
    combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals
    nowhite = choice of white separating the diverging colors or not
    """

    if (int(n/2) == n/2):
        #even number of contours
        nwhite=1
        nneg = n/2
        npos = n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    if (nowhite):
        nwhite=0

    colors1 = plt.cm.YlOrRd_r(np.linspace(0.1,1,int(npos)))
    colors2 = plt.cm.Blues(np.linspace(0.1,1,int(nneg)))
    colorsw = np.ones((nwhite,4))

    if (nowhite):
        colors = np.vstack( (colors1, colors2))
    else:
        colors = np.vstack((colors1, colorsw, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
  
    return mymap

#########

def precip_cmap(n, nowhite=False):
    """
    combine two existing color maps to create a diverging color map with white in the middle.
    browns for negative, blues for positive
    n = the number of contour intervals
    nowhite = choice of white separating the diverging colors or not
    """
    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2

    if nowhite:
        colors1 = plt.cm.YlOrBr_r(np.linspace(0,0.8, int(nneg)))
        colors2 = plt.cm.GnBu(np.linspace(0.2,1, int(npos)))
        colors = np.vstack((colors1, colors2))
    else:
        colors1 = plt.cm.YlOrBr_r(np.linspace(0,1, int(nneg)))
        colors2 = plt.cm.GnBu(np.linspace(0,1, int(npos)))
        colorsw = np.ones((nwhite,4))
        colors = np.vstack((colors1, colorsw, colors2))

    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    return mymap

#########

def fixcesmtime(dat,syear,eyear):
    """
    Fix the CESM timestamp with a simple set of dates
    """
    timefix = pd.date_range(start=f'1/1/{syear}', end=f'12/1/{eyear}', freq='MS') # generic time coordinate from a non-leap-year
    dat = dat.assign_coords({"time":timefix})

    return dat

#########

def get5by5coords_zmplots():
    """
    positioning for 5x5 plots
    """
    x1 = [0.02,0.225,0.43,0.635,0.84,
          0.02,0.225,0.43,0.635,0.84,
          0.02,0.225,0.43,0.635,0.84,
          0.02,0.225,0.43,0.635,0.84,
          0.02,0.225,0.43,0.635,0.84] 
    x2 = [0.18,0.385,0.59,0.795,1,
          0.18,0.385,0.59,0.795,1,
          0.18,0.385,0.59,0.795,1,
          0.18,0.385,0.59,0.795,1,
          0.18,0.385,0.59,0.795,1]
    y1 = [0.8,0.8,0.8,0.8,0.8,
          0.6,0.6,0.6,0.6,0.6,
          0.4,0.4,0.4,0.4,0.4,
          0.2,0.2,0.2,0.2,0.2,
          0.,0.,0.,0.,0.]
    y2 = [0.95,0.95,0.95,0.95,0.95,
          0.75,0.75,0.75,0.75,0.75,
          0.55,0.55,0.55,0.55,0.55,
          0.35,0.35,0.35,0.35,0.35,
          0.15,0.15,0.15,0.15,0.15]
    
    return x1, x2, y1, y2

#########

def plotcolorbar(fig, ci, cmin, cmax, titlestr, x1, x2, y1, y2, 
   cmap='blue2red', orient='horizontal', posneg='both', ticks=None, fsize=14, nowhite=False,
   contourlines=False, contourlinescale=1):
    """
    plot a color bar
       Input:
           fig = the figure identified
           ci = the contour interval for the color map
           cmin = the minimum extent of the contour range
           cmax = the maximum extent of the contour range
           titlestr = the label for the color bar
           x1 = the location of the left edge of the color bar
           x2 = the location of the right edge of the color bar
           y1 = the location of the bottom edge of the color bar
           y2 = the location of the top edge of the color bar
           cmap = the color map to be used (only set up for blue2red at the moment)
           orient = the orientation (horizontal or vertical)
           posneg = if "both", both positive and negative sides are plotted
                    if "pos", only the positive side is plotted
                    if "neg", only the negative side is plotted
           ticks = user specified ticklabels
           fsize = user specified font size
           contourlines = used to overplot contour lines
           contourlinescale = scale factor for contour lines to be overplotted
           nowhite = choice of white separating the diverging colors or not
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = ci * np.arange(cmin/ci, (cmax+ci)/ci, 1) 

    if (cmap == "blue2red"):
        mymap = blue2red_cmap(nlevs, nowhite)

    if (cmap == "precip"):
        mymap = precip_cmap(nlevs, nowhite)

    if (cmap == "precip_nowhite"):
        mymap = precip_cmap(nlevs, nowhite=True)

    if (cmap == 'red2blue'):
        mymap = red2blue_cmap(nlevs, nowhite)

    clevplot=clevs
    if (posneg == "pos"):
        clevplot = clevs[clevs >= 0]
    if (posneg == "neg"):
        clevplot = clevs[clevs <= 0]

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    norm = mcolors.Normalize(vmin=cmin, vmax=cmax)
    
    if (ticks):
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
           orientation=orient, norm=norm, values=clevplot, ticks=ticks)
    else:
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap, 
           orientation=orient, norm=norm, values=clevplot)

    clb.ax.tick_params(labelsize=fsize)
    clb.set_label(titlestr, fontsize=fsize+2)

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci/2.]
        if (orient=='horizontal'):
            ax.vlines(clevlines[clevlines > 0],-5,5, colors='black', linestyle='solid')
            ax.vlines(clevlines[clevlines < 0],-5,5, colors='black', linestyle='dashed')
        if (orient=='vertical'):
            ax.hlines(clevlines[clevlines > 0],-10,15, colors='black', linestyle='solid')
            ax.hlines(clevlines[clevlines < 0],-10,15, colors='black', linestyle='dashed')

    return ax

#########

def cosweightlat(darray, lat1, lat2):
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

#########

def plot_pre_mon(fig, data, ci, cmin, cmax, expname, x1=None, x2=None, y1=None, y2=None, 
                 oplot=False, ax=None, cmap='precip', taxis='time', paxis='lev', climo_yrs=None):
    """
    Plot seasonal cycle, pressure versus time.
    """

    # move the time axis to the first
    if (data.dims[1] != taxis):
        data = data.transpose(..., taxis)

    #Make 24 months so we can have Jan-Dec repeated twice
    case_seas = np.zeros((25,len(data[paxis])))
    case_seas = xr.DataArray(case_seas, dims=[taxis,paxis],
                                     coords={taxis: np.arange(1,26,1),
                                             paxis: data[paxis]})
    #Make array of monthly temp data
    for m in range(0,25):
        month = m
    
        if m > 11:
            month = m-12
        if month == 12:
            month = 0    
        case_seas[m] = data.sel(month=month+1)

    nlevs = (cmax - cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = precip_cmap(nlevs)

    if (cmap == "precip_nowhite"):
        mymap = precip_cmap(nlevs, nowhite=True)

    # if overplotting, check for axis input
    if (oplot and (not ax)):
        print("This isn't going to work.  If overplotting, specify axis")
        return

    plt.rcParams['font.size'] = '14'

    case_seas = case_seas.transpose(..., taxis)

    monticks_temp = np.arange(0,25,1)
    monticks = monticks_temp

    if not oplot:
        if (x1):
            ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
        else:
            ax = fig.add_axes()
    ax.xaxis.set_label_position('top')
    if climo_yrs:
        ax.set_xlabel(f"{climo_yrs}", loc='center',
                           fontsize=8)
    
    #ax.contourf(monticks_temp, -np.log10(case_seas[paxis]), case_seas*(29/18), levels=clevs*(29/18), cmap=mymap, extend='max')
    #c= ax.contour(monticks_temp, -np.log10(case_seas[paxis]), case_seas*(29/18), levels=clevs[::3]*(29/18), colors="k", extend='max',linewidths=0.25)

    ax.contourf(monticks_temp, -np.log10(case_seas[paxis]), case_seas, levels=clevs, cmap=mymap, extend='max')
    c= ax.contour(monticks_temp, -np.log10(case_seas[paxis]), case_seas, levels=clevs, colors="k", extend='max',linewidths=0.25) #clevs[::3]
    fmt = {lev: '{:.1f}'.format(lev) for lev in c.levels}
    ax.clabel(c, c.levels, inline=True, fmt=fmt, fontsize=8)
    ax.set_ylim(-np.log10(100),-np.log10(3))
    ax.set_yticks([-np.log10(100),-np.log10(30),-np.log10(10),-np.log10(3)])
    ax.set_yticklabels(['100','30','10','3'])
    ax.set_xticks(monticks[0:25:3])
    ax.set_xticklabels(['Jan','Apr','Jul','Oct']*2+["Jan"], fontsize=10)
    ax.set_title(expname, fontsize=16)

    return ax

#########

###############
