"""
Generate plots that compare ENSO characteristics across various versions of CESM development
"""

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import warnings # use to warn user about missing files
from pathlib import Path

def enso_comparison_plots(adfobj):
    """
    This script/function is designed to generate ENSO-related plots across
    various CESM simulations

    Parameters
    ----------
    adfobj : AdfDiag
        The diagnostics object that contains all the configuration information

    Returns
    -------
    Does not return a value; produces plots and saves files.
    """
    
    # Notify user that script has started:
    msg = "\n  Generating ENSO plots to compare against all runs..."
    print(f"{msg}\n  {'-' * (len(msg)-3)}")

    plot_locations = adfobj.plot_location
    plot_type = adfobj.get_basic_info('plot_type')
    if not plot_type:
        plot_type = 'png'

    # check if existing plots need to be redone
    redo_plot = adfobj.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")


    #Grab saved files
    obs_ds   = xr.open_dataset('/glade/derecho/scratch/mdfowler/ENSOmetrics_Obs.nc')
    cesm1_ds = xr.open_dataset("/glade/derecho/scratch/mdfowler/ENSOmetrics_CESM1.nc")
    cesm2_ds = xr.open_dataset("/glade/derecho/scratch/mdfowler/ENSOmetrics_CESM2.nc")
    dev_ds   = xr.open_dataset("/glade/derecho/scratch/mdfowler/ENSOmetrics_CESM3dev.nc")

    #+++++++++++++
    # Make Nino Variance Plot 
    #+++++++++++++

    #Set path for variance figures:
    plot_loc_ts  = Path(plot_locations[0]) / f'Special_Nino34_Variance_Special_Mean.{plot_type}'

    # Check redo_plot. If set to True: remove old plots, if they already exist:
    if (not redo_plot) and plot_loc_ts.is_file():
        #Add already-existing plot to website (if enabled):
        adfobj.debug_log(f"'{plot_loc_ts}' exists and clobber is false.")
        adfobj.add_website_data(plot_loc_ts, "Special", None, season="Nino34_Variance", multi_case=True, non_season=True)

        #Continue to next iteration:
        return
    elif (redo_plot):
        if plot_loc_ts.is_file():
            plot_loc_ts.unlink()
    #End if



    #                          J   F   M    A  M   J   J   A    S   O    N  D
    daysPerMonth = np.asarray([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

    fig,axs=plt.subplots(1,1,figsize=(15,5))

    for iCase in range(len(dev_ds.case.values)):
        case_max = np.nanmax(dev_ds.nino34_variance.values[iCase,:])
        case_min = np.nanmin(dev_ds.nino34_variance.values[iCase,:])
        # Weighted mean
        weights = ( daysPerMonth / daysPerMonth.sum() )
        weighted_mean = (dev_ds.nino34_variance.values[iCase,:] * weights).sum() / weights.sum()
        
        axs.plot(iCase+np.ones(2), [case_min, case_max], 'k-')
        axs.plot(iCase+1, weighted_mean,'o', color='k')
        axs.plot(iCase+1, case_min,'^',color='k')
        axs.plot(iCase+1, case_max,'v',color='k')


    cesm1_xCenter = len(dev_ds.case.values)+1.5
    cesm2_xCenter = len(dev_ds.case.values)+3
    obs_xCenter   = len(dev_ds.case.values)+4

    offset = np.linspace(cesm1_xCenter-0.5, cesm1_xCenter+0.5, len(cesm1_ds.event.values))
    for iEvent in range(len(cesm1_ds.event.values)):
        case_max = np.nanmax(cesm1_ds.nino34_variance.values[iEvent,:])
        case_min = np.nanmin(cesm1_ds.nino34_variance.values[iEvent,:])
        # Weighted mean
        weights = ( daysPerMonth / daysPerMonth.sum() )
        weighted_mean = (cesm1_ds.nino34_variance.values[iEvent,:] * weights).sum() / weights.sum()
        
        axs.plot(offset[iEvent]*np.ones(2), [case_min, case_max],
                '-', color='mediumpurple', alpha=0.2)
        axs.plot(offset[iEvent], weighted_mean,'o', color='mediumpurple',alpha=0.4)
        axs.plot(offset[iEvent], case_min,'^',color='mediumpurple',alpha=0.4)
        axs.plot(offset[iEvent],case_max,'v',color='mediumpurple',alpha=0.4)


    offset = np.linspace(cesm2_xCenter-0.5, cesm2_xCenter+0.5, len(cesm2_ds.event.values))
    for iEvent in range(len(cesm2_ds.event.values)):
        case_max = np.nanmax(cesm2_ds.nino34_variance.values[iEvent,:])
        case_min = np.nanmin(cesm2_ds.nino34_variance.values[iEvent,:])
        # Weighted mean
        weights = ( daysPerMonth / daysPerMonth.sum() )
        weighted_mean = (cesm2_ds.nino34_variance.values[iEvent,:] * weights).sum() / weights.sum()
        
        axs.plot(offset[iEvent]*np.ones(2), [case_min, case_max],
                '-', color='orange', alpha=0.2)
        axs.plot(offset[iEvent], weighted_mean,'o', color='orange',alpha=0.4)
        axs.plot(offset[iEvent], case_min,'^',color='orange',alpha=0.4)
        axs.plot(offset[iEvent],case_max,'v',color='orange',alpha=0.4)


    # Get obs 
    obs_max = np.nanmax(obs_ds.nino34_variance.values)
    obs_min = np.nanmin(obs_ds.nino34_variance.values)
    # Weighted mean
    weights = ( daysPerMonth / daysPerMonth.sum() )
    weighted_mean = (obs_ds.nino34_variance.values * weights).sum() / weights.sum()

    axs.plot(obs_xCenter*np.ones(2), [obs_min, obs_max], '-', color='firebrick')
    axs.plot(obs_xCenter, weighted_mean,'o', color='firebrick')
    axs.plot(obs_xCenter, obs_min,'^',color='firebrick')
    axs.plot(obs_xCenter, obs_max,'v',color='firebrick')

    ## General plot settings 
    ticks = np.append(1+np.arange(len(dev_ds.case.values)), cesm1_xCenter)
    ticks = np.append(ticks, cesm2_xCenter)
    ticks = np.append(ticks, obs_xCenter)

    tickLabels = np.append(dev_ds.case.values, 'CESM1')
    tickLabels = np.append(tickLabels, 'CESM2')
    tickLabels = np.append(tickLabels, 'HadiSST')

    axs.set_xticks(ticks)
    axs.set_xticklabels(tickLabels)
    plt.setp( axs.xaxis.get_majorticklabels(), rotation=45 )

    axs.axhline(obs_min, color='firebrick',alpha=0.3)
    axs.axhline(obs_max, color='firebrick',alpha=0.3)

    axs.set_title('Monthly nino3.4 variance')

    #Save figure to file:
    fig.savefig(plot_loc_ts, bbox_inches='tight', facecolor='white')

    #Add plot to website (if enabled):
    adfobj.add_website_data(plot_loc_ts, "Special", None, season="Nino34_Variance", multi_case=True, non_season=True)

    return