#Import standard modules:
from pathlib import Path
import numpy as np
import xarray as xr
import warnings  # use to warn user about missing files.
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib as mpl
import matplotlib.cm as cm
import pandas as pd

import plotting_functions as pf

#Format warning messages:
def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = my_formatwarning

def tem(adf):
    """
    Plot the contents of the TEM dignostic ouput of 2-D latitude vs vertical pressure maps.
    
    Steps:
     - loop through TEM variables
     - calculate all-time fields (from individual months)
     - Take difference, calculate statistics
     - make plot

    """

    #Special ADF variable which contains the output paths for
    #all generated plots and tables for each case:
    plot_location = Path(adf.plot_location[0])

    #Check if plot output directory exists, and if not, then create it:
    if not plot_location.is_dir():
        print(f"    {plot_location} not found, making new directory")
        plot_location.mkdir(parents=True)

    #CAM simulation variables (this is always assumed to be a list):
    case_names = adf.get_cam_info("cam_case_name", required=True)

    res = adf.variable_defaults # will be dict of variable-specific plot preferences

    #Check if comparing against observations
    if adf.compare_obs:
        obs = True
        base_name = "Obs"
    else:
        obs = False
        base_name = adf.get_baseline_info("cam_case_name", required=True)
    #End if

    #Extract test case years
    syear_cases = adf.climo_yrs["syears"]
    eyear_cases = adf.climo_yrs["eyears"]

    #Extract baseline years (which may be empty strings if using Obs):
    syear_baseline = adf.climo_yrs["syear_baseline"]
    eyear_baseline = adf.climo_yrs["eyear_baseline"]

    #Grab all case nickname(s)
    test_nicknames = adf.case_nicknames["test_nicknames"]
    base_nickname = adf.case_nicknames["base_nickname"]
    case_nicknames = test_nicknames + [base_nickname]
 
    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adf.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')
    print(f"\t NOTE: Plot type is set to {plot_type}")

    # check if existing plots need to be redone
    redo_plot = adf.get_basic_info('redo_plot')
    print(f"\t NOTE: redo_plot is set to {redo_plot}")
    #-----------------------------------------

    #Grab TEM diagnostics options
    tem_opts = adf.read_config_var("tem_info")

    if not tem_opts:
        print("\n  No TEM options provided, skipping TEM plots." \
        "\nSee documentation or config_cam_baseline_example.yaml for options to add to configuration file.")
        return

    #Location of saved TEM netCDF files
    tem_loc = tem_opts.get("tem_loc")

    #If path not specified, skip TEM calculation
    if tem_loc is None:
        print("'tem_loc' not found in config file, so TEM plots will be skipped.")
        return
    else:
        #Notify user that script has started:
        print("\n  Generating TEM plots...")
    
    #Set seasonal ranges:
    seasons = {"ANN": np.arange(1,13,1),
               "DJF": [12, 1, 2],
               "JJA": [6, 7, 8],
               "MAM": [3, 4, 5],
               "SON": [9, 10, 11]
               }

    #Suggestion from Rolando, if QBO is being produced, add utendvtem and utendwtem?
    if "qbo" in adf.plotting_scripts:
        var_list = ['uzm','thzm','epfy','epfz','vtem','wtem',
                    'psitem','utendepfd','utendvtem','utendwtem']
    #Otherwise keep it simple
    else:
        var_list = ['uzm','thzm','epfy','epfz','vtem','wtem','psitem','utendepfd']

    #Baseline TEM location
    input_loc_idx = Path(tem_loc) / base_name

    #Check if comparing against obs
    if obs:
        #Set TEM file for observations
        base_file_name = f'{base_name}.TEMdiag.nc'
    else:
        #Set TEM file for baseline
        base_file_name = f'{base_name}.TEMdiag_{syear_baseline}-{eyear_baseline}.nc'
    
    #Set full path for baseline/obs file
    tem_base = input_loc_idx / base_file_name

    #Check to see if baseline/obs TEM file exists    
    if tem_base.is_file():
        ds_base = xr.open_dataset(tem_base)
    else:
        print(f"\t'{base_file_name}' does not exist. TEM plots will be skipped.")
        return

    """
    #Setup TEM plots
    nrows = len(var_list)
    ncols = len(case_nicknames)+1
    fig_width = 20

    #try and dynamically create size of fig based off number of cases
    fig_height = 15+(ncols*nrows)

    #Loop over season dictionary:
    for s in seasons:
        #Location to save plots
        plot_name = plot_location / f"TEM_{s}_WACCM_SeasonalCycle_Mean.png"

        # Check redo_plot. If set to True: remove old plot, if it already exists:
        if (not redo_plot) and plot_name.is_file():
            #Add already-existing plot to website (if enabled):
            adf.debug_log(f"'{plot_name}' exists and clobber is false.")
            adf.add_website_data(plot_name, "TEM", None, season=s, plot_type="WACCM",ext="SeasonalCycle_Mean",category="Seasonal Cycle",multi_case=True)

        #plot_name = plot_loc / f"CPT_ANN_WACCM_SeasonalCycle_Mean.{plot_type}"
        elif ((redo_plot) and plot_name.is_file()) or (not plot_name.is_file()):
            if plot_name.is_file():
                plot_name.unlink()
        
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,fig_height),
                                    facecolor='w', edgecolor='k')

            #Loop over model cases:
            for idx,case_name in enumerate(case_names):

                # Check redo_plot. If set to True: remove old plot, if it already exists:
                #if (not redo_plot) and plot_name.is_file():
                    #Add already-existing plot to website (if enabled):
                #    adf.debug_log(f"'{plot_name}' exists and clobber is false.")
                #    adf.add_website_data(plot_name, "TEM", case_name, season=s, plot_type="WACCM",ext="Mean",category="Seasonal Cycle")

                    #Continue to next iteration:
                    #continue
                #elif (redo_plot) and plot_name.is_file():
                #    plot_name.unlink()

                #Extract start and end year values:
                start_year = syear_cases[idx]
                end_year   = eyear_cases[idx]

                #Open the TEM file
                output_loc_idx = Path(tem_loc) / case_name
                case_file_name = f'{case_name}.TEMdiag_{start_year}-{end_year}.nc'
                tem = output_loc_idx / case_file_name

                #Grab the data for the TEM netCDF files
                if tem.is_file():
                    ds = xr.open_dataset(tem)
                else:
                    print(f"\t'{case_file_name}' does not exist. TEM plots will be skipped.")
                    return

                climo_yrs = {"test":[syear_cases[idx], eyear_cases[idx]],
                            "base":[syear_baseline, eyear_baseline]}

                #Setup and plot the sub-plots
                tem_plot(ds, ds_base, case_nicknames, axs, s, var_list, res, obs, climo_yrs)

            #Set figure title
            plt.suptitle(f'TEM Diagnostics: {s}', fontsize=20, y=.928)

            #Write the figure to provided workspace/file:
            fig.savefig(plot_name, bbox_inches='tight', dpi=300)

            #Add plot to website (if enabled):
            adf.add_website_data(plot_name, "TEM", case_name, season=s, plot_type="WACCM",ext="SeasonalCycle_Mean",category="Seasonal Cycle")
    """
    print("  ...TEM plots have been generated successfully.")


    #Loop over variables:
    for var in var_list:
        """
        if adf.compare_obs:
            #Check if obs exist for the variable:
            if var in var_obs_dict:
                #Note: In the future these may all be lists, but for
                #now just convert the target_list.
                #Extract target file:
                dclimo_loc = var_obs_dict[var]["obs_file"]
                #Extract target list (eventually will be a list, for now need to convert):
                data_list = [var_obs_dict[var]["obs_name"]]
                #Extract target variable name:
                data_var = var_obs_dict[var]["obs_var"]
            else:
                dmsg = f"No obs found for variable `{var}`, zonal mean plotting skipped."
                adfobj.debug_log(dmsg)
                continue
            #End if
        else:
            #Set "data_var" for consistent use below:
            data_var = var
        #End if
        """

        #Loop over model cases:
        for idx,case_name in enumerate(case_names):

            # Check redo_plot. If set to True: remove old plot, if it already exists:
            #if (not redo_plot) and plot_name.is_file():
                #Add already-existing plot to website (if enabled):
            #    adf.debug_log(f"'{plot_name}' exists and clobber is false.")
            #    adf.add_website_data(plot_name, "TEM", case_name, season=s, plot_type="WACCM",ext="Mean",category="Seasonal Cycle")

                #Continue to next iteration:
                #continue
            #elif (redo_plot) and plot_name.is_file():
            #    plot_name.unlink()

            #Extract start and end year values:
            start_year = syear_cases[idx]
            end_year   = eyear_cases[idx]

            #Open the TEM file
            output_loc_idx = Path(tem_loc) / case_name
            case_file_name = f'{case_name}.TEMdiag_{start_year}-{end_year}.nc'
            tem = output_loc_idx / case_file_name

            #Grab the data for the TEM netCDF files
            if tem.is_file():
                ds = xr.open_dataset(tem)
            else:
                print(f"\t'{case_file_name}' does not exist. TEM plots will be skipped.")
                return

            #Loop over season dictionary:
            for s in seasons:

                #Location to save plots
                plot_name = plot_location / f"{var}_{s}_WACCM_SeasonalCycle_Mean.png"

                # Check redo_plot. If set to True: remove old plot, if it already exists:
                if (not redo_plot) and plot_name.is_file():
                    #Add already-existing plot to website (if enabled):
                    adf.debug_log(f"'{plot_name}' exists and clobber is false.")
                    adf.add_website_data(plot_name, var, None, season=s, plot_type="WACCM",ext="SeasonalCycle_Mean",category="TEM",multi_case=True)

                #plot_name = plot_loc / f"CPT_ANN_WACCM_SeasonalCycle_Mean.{plot_type}"
                elif ((redo_plot) and plot_name.is_file()) or (not plot_name.is_file()):
                    if plot_name.is_file():
                        plot_name.unlink()

                print(var,s)
                #Grab variable defaults for this variable
                vres = res[var]

                #Gather data for both cases
                mdata = ds[var].squeeze()
                odata = ds_base[var].squeeze()

                # APPLY UNITS TRANSFORMATION IF SPECIFIED:
                # NOTE: looks like our climo files don't have all their metadata
                mdata = mdata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
                # update units
                mdata.attrs['units'] = vres.get("new_unit", mdata.attrs.get('units', 'none'))

                # Do the same for the baseline case if need be:
                if not obs:
                    odata = odata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
                    # update units
                    odata.attrs['units'] = vres.get("new_unit", odata.attrs.get('units', 'none'))
                # Or for observations
                else:
                    odata = odata * vres.get("obs_scale_factor",1) + vres.get("obs_add_offset", 0)
                    # Note: we are going to assume that the specification ensures the conversion makes the units the same. Doesn't make sense to add a different unit.

                #Create array to avoid weighting missing values:
                md_ones = xr.where(mdata.isnull(), 0.0, 1.0)
                od_ones = xr.where(odata.isnull(), 0.0, 1.0)

                month_length = mdata.time.dt.days_in_month
                weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())

                #Calculate monthly-weighted seasonal averages:
                if s == 'ANN':

                    #Calculate annual weights (i.e. don't group by season):
                    weights_ann = month_length / month_length.sum()

                    mseasons = (mdata * weights_ann).sum(dim='time')
                    mseasons = mseasons / (md_ones*weights_ann).sum(dim='time')

                    #Calculate monthly weights based on number of days:
                    if obs:
                        month_length_obs = odata.time.dt.days_in_month
                        weights_ann_obs = month_length_obs / month_length_obs.sum()
                        oseasons = (odata * weights_ann_obs).sum(dim='time')
                        oseasons = oseasons / (od_ones*weights_ann_obs).sum(dim='time')
                    else:
                        month_length_base = odata.time.dt.days_in_month
                        weights_ann_base = month_length_base / month_length_base.sum()
                        oseasons = (odata * weights_ann_base).sum(dim='time')
                        oseasons = oseasons / (od_ones*weights_ann_base).sum(dim='time')

                else:
                    #this is inefficient because we do same calc over and over
                    mseasons = (mdata * weights).groupby("time.season").sum(dim="time").sel(season=s)
                    wgt_denom = (md_ones*weights).groupby("time.season").sum(dim="time").sel(season=s)
                    mseasons = mseasons / wgt_denom

                    if obs:
                        month_length_obs = odata.time.dt.days_in_month
                        weights_obs = (month_length_obs.groupby("time.season") / month_length_obs.groupby("time.season").sum())
                        oseasons = (odata * weights_obs).groupby("time.season").sum(dim="time").sel(season=s)
                        wgt_denom = (od_ones*weights_obs).groupby("time.season").sum(dim="time").sel(season=s)
                        oseasons = oseasons / wgt_denom
                    else:
                        month_length_base = odata.time.dt.days_in_month
                        weights_base = (month_length_base.groupby("time.season") / month_length_base.groupby("time.season").sum())
                        oseasons = (odata * weights_base).groupby("time.season").sum(dim="time").sel(season=s)
                        wgt_denom_base = (od_ones*weights_base).groupby("time.season").sum(dim="time").sel(season=s)
                        oseasons = oseasons / wgt_denom_base

                #difference: each entry should be (lat, lon)
                dseasons = mseasons-oseasons
                #print(dseasons.min(),dseasons.max())

                #Run through variables and plot each against the baseline on each row
                #Each column will be a case, ie (test, base, difference)

                if 'units' in vres:
                    units = vres['units']
                else:
                    units = ''

            
                
                #
                cp_info = pf.prep_contour_plot(mseasons, oseasons, dseasons, **vres)
                clevs = np.unique(np.array(cp_info['levels1']))
                print(clevs)
                norm = cp_info['norm1']
                cmap = cp_info['cmap1']
                print(cmap,"\n")


                levs_diff = np.unique(np.array(cp_info['levelsdiff']))



                

                # mesh for plots:
                lat = mseasons['zalat']
                lev = mseasons['lev']
                lats, levs = np.meshgrid(lat, lev)

                # Generate zonal plot:
                #fig, ax = plt.subplots(figsize=(10,10),nrows=3, constrained_layout=True,
                #                    sharex=True, sharey=True,**cp_info['subplots_opt'])
                # create figure object
                fig = plt.figure(figsize=(14,10))
                # LAYOUT WITH GRIDSPEC
                gs = mpl.gridspec.GridSpec(4, 8, wspace=0.75,hspace=0.5) # 2 rows, 4 columns, but each map will take up 2 columns
                #gs.tight_layout(fig)
                """ax1 = plt.subplot(gs[0:2, :2], **cp_info['subplots_opt'])
                ax2 = plt.subplot(gs[0:2, 2:], **cp_info['subplots_opt'])
                ax3 = plt.subplot(gs[2:, 1:3], **cp_info['subplots_opt'])"""
                ax1 = plt.subplot(gs[0:2, :4], **cp_info['subplots_opt'])
                ax2 = plt.subplot(gs[0:2, 4:], **cp_info['subplots_opt'])
                ax3 = plt.subplot(gs[2:, 2:6], **cp_info['subplots_opt'])
                ax = [ax1,ax2,ax3]

                #Get axis number for variable
                #axs_id = var_axs[var]
                #print(axs_id)

                #x_filtered = x[~np.isnan(y)]

                #Contour fill
                img0 = ax[0].contourf(lats, levs,mseasons, levels=clevs, norm=norm, cmap=cmap)
                img1 = ax[1].contourf(lats, levs,oseasons, levels=clevs, norm=norm, cmap=cmap)

                # Define a custom formatting function
                from matplotlib.ticker import FuncFormatter
                def format_contour_label(x, pos):
                    if abs(x) > 1000:
                        #x = "{:e}".format(x)
                        return '{:.2f}'.format(x)
                        #return x.split('e')[0]
                    if x == 0:#else:
                        return x
                    
                #Add contours for highlighting
                c0 = ax[0].contour(lats,levs,mseasons,levels=clevs[::2], norm=norm, colors="k",linewidths=0.5)
                # Add contour labels every third contour line
                plt.clabel(c0, inline=True, fontsize=8, levels=c0.levels[::2],
                            #fmt=FuncFormatter(format_contour_label)
                            )

                c1 = ax[1].contour(lats,levs,oseasons,levels=clevs[::2], norm=norm, colors="k",linewidths=0.5)
                # Add contour labels every third contour line
                plt.clabel(c1, inline=True, fontsize=8, levels=c1.levels[::2],
                            #fmt=FuncFormatter(format_contour_label)
                            )


                #Check if difference plot has contour levels, if not print notification
                if len(dseasons.lev) == 0:
                    #Set empty message for comparison of cases with different vertical levels
                    #TODO: Work towards getting the vertical and horizontal interpolations!! - JR
                    empty_message = "These have different vertical levels\nCan't compare cases currently"
                    props = {'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.9}
                    prop_x = 0.18
                    prop_y = 0.42
                    ax[2].text(prop_x, prop_y, empty_message,
                                    transform=ax[2].transAxes, bbox=props)
                else:
                    img2 = ax[2].contourf(lats,levs,dseasons, cmap="BrBG",levels=levs_diff,norm=cp_info['normdiff'])#levels=levs_diff
                    ax[2].contour(lats,levs,dseasons, colors="k",linewidths=0.5,levels=levs_diff[::2],norm=cp_info['normdiff'])#levels=diff_levs[::2]
                    plt.colorbar(img2, ax=ax[2], location='right',**cp_info['diff_colorbar_opt'])#**cp_info['diff_colorbar_opt']

                #Format y-axis
                for a in ax[:]:
                    a.set_yscale("log")
                    a.set_ylim(ax[2].get_ylim()[::-1])

                plt.colorbar(img0, ax=ax[0], location='right',**cp_info['colorbar_opt'])
                plt.colorbar(img1, ax=ax[1], location='right',**cp_info['colorbar_opt'])

                #Set titles of subplots
                #Set figure title
                plt.suptitle(f'TEM Diagnostics: {s}', fontsize=20, y=.98)

                #Variable plot title name
                longname = vres["long_name"]
                #ax.text(0.5, 0.95, 'Title 1', transform=ax.transAxes, fontsize=14,
                #        verticalalignment='top', horizontalalignment='center')
                plt.text(0.5, 0.915, f"{longname}\n", fontsize=12, ha='center', transform=fig.transFigure)
                #ax[1].set_title(longname+"\n",fontsize=14,loc="center")

                test_yrs = f"{start_year}-{end_year}"
                #ax[0].set_title(f"\n\n"+"$\mathbf{Test}$"+f"  yrs: {test_yrs}\n"+f"{test_nicknames[idx]}\n",fontsize=10)
                ax[0].set_title(f"{test_nicknames[idx]}\n{test_yrs}",fontsize=10)

                if obs:
                    obs_title = Path(vres["obs_name"]).stem
                    #ax[1].set_title(f"\n\n"+"$\mathbf{Baseline}$\n"+f"{obs_title}"+longname+"\n",fontsize=14)
                    ax[1].set_title(f"{obs_title}\n",fontsize=10)

                else:
                    base_yrs = f"{syear_baseline}-{eyear_baseline}"
                    #ax[1].set_title(f"\n\n"+"$\mathbf{Baseline}$"+f"  yrs: {base_yrs}\n"+f"{base_nickname}\n",fontsize=10)
                    ax[1].set_title(f"{base_nickname}\n{base_yrs}",fontsize=10)
                
                #Set main title for difference plots column
                ax[2].set_title("$\mathbf{Test} - \mathbf{Baseline}$",fontsize=10)

                #Write the figure to provided workspace/file:
                fig.savefig(plot_name, bbox_inches='tight', dpi=300)

                #Add plot to website (if enabled):
                adf.add_website_data(plot_name, var, case_name, season=s, plot_type="WACCM",ext="SeasonalCycle_Mean",category="TEM")
                #Set variable name on center plot (except first plot, see above)
                #for i in range(1,len(var_list)):
                    #vres = res[var_list[i]]

                #Create new plot:
                #pf.plot_zonal_mean_and_save(plot_name, case_nickname, base_nickname,
                #                                            [syear_cases[case_idx],eyear_cases[case_idx]],
                #                                            [syear_baseline,eyear_baseline],
                #                                            mseasons[s], oseasons[s], has_lev=True, log_p=False, obs=obs, **vres)
                plt.close()
        

# Helper functions
##################

def tem_plot(ds, ds_base, case_names, axs, s, var_list, res, obs, climo_yrs):
    """
    TEM subplots
    
    """

    var_axs = {"uzm":0, "thzm":1, "epfy":2, "epfz":3, "vtem":4, "wtem":5,
                   "psitem":6, "utendepfd":7, "utendvtem":8, "utendwtem":9
                  }

    #Set empty message for comparison of cases with different vertical levels
    #TODO: Work towards getting the vertical and horizontal interpolations!! - JR
    empty_message = "These have different vertical levels\nCan't compare cases currently"
    props = {'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.9}
    prop_x = 0.18
    prop_y = 0.42

    for var in var_list:
        print(var,s)
        #Grab variable defaults for this variable
        vres = res[var]

        #Gather data for both cases
        mdata = ds[var].squeeze()
        odata = ds_base[var].squeeze()

        # APPLY UNITS TRANSFORMATION IF SPECIFIED:
        # NOTE: looks like our climo files don't have all their metadata
        mdata = mdata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
        # update units
        mdata.attrs['units'] = vres.get("new_unit", mdata.attrs.get('units', 'none'))

        # Do the same for the baseline case if need be:
        if not obs:
            odata = odata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
            # update units
            odata.attrs['units'] = vres.get("new_unit", odata.attrs.get('units', 'none'))
        # Or for observations
        else:
            odata = odata * vres.get("obs_scale_factor",1) + vres.get("obs_add_offset", 0)
            # Note: we are going to assume that the specification ensures the conversion makes the units the same. Doesn't make sense to add a different unit.

        #Create array to avoid weighting missing values:
        md_ones = xr.where(mdata.isnull(), 0.0, 1.0)
        od_ones = xr.where(odata.isnull(), 0.0, 1.0)

        month_length = mdata.time.dt.days_in_month
        weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())

        #Calculate monthly-weighted seasonal averages:
        if s == 'ANN':

            #Calculate annual weights (i.e. don't group by season):
            weights_ann = month_length / month_length.sum()

            mseasons = (mdata * weights_ann).sum(dim='time')
            mseasons = mseasons / (md_ones*weights_ann).sum(dim='time')

            #Calculate monthly weights based on number of days:
            if obs:
                month_length_obs = odata.time.dt.days_in_month
                weights_ann_obs = month_length_obs / month_length_obs.sum()
                oseasons = (odata * weights_ann_obs).sum(dim='time')
                oseasons = oseasons / (od_ones*weights_ann_obs).sum(dim='time')
            else:
                month_length_base = odata.time.dt.days_in_month
                weights_ann_base = month_length_base / month_length_base.sum()
                oseasons = (odata * weights_ann_base).sum(dim='time')
                oseasons = oseasons / (od_ones*weights_ann_base).sum(dim='time')

        else:
            #this is inefficient because we do same calc over and over
            mseasons = (mdata * weights).groupby("time.season").sum(dim="time").sel(season=s)
            wgt_denom = (md_ones*weights).groupby("time.season").sum(dim="time").sel(season=s)
            mseasons = mseasons / wgt_denom

            if obs:
                month_length_obs = odata.time.dt.days_in_month
                weights_obs = (month_length_obs.groupby("time.season") / month_length_obs.groupby("time.season").sum())
                oseasons = (odata * weights_obs).groupby("time.season").sum(dim="time").sel(season=s)
                wgt_denom = (od_ones*weights_obs).groupby("time.season").sum(dim="time").sel(season=s)
                oseasons = oseasons / wgt_denom
            else:
                month_length_base = odata.time.dt.days_in_month
                weights_base = (month_length_base.groupby("time.season") / month_length_base.groupby("time.season").sum())
                oseasons = (odata * weights_base).groupby("time.season").sum(dim="time").sel(season=s)
                wgt_denom_base = (od_ones*weights_base).groupby("time.season").sum(dim="time").sel(season=s)
                oseasons = oseasons / wgt_denom_base

        #difference: each entry should be (lat, lon)
        dseasons = mseasons-oseasons
        #print(dseasons.min(),dseasons.max())

        #Run through variables and plot each against the baseline on each row
        #Each column will be a case, ie (test, base, difference)
        """
        if i == len(wrap_fields)-1:
            levels = cp_info['levelsdiff']
            cmap = cp_info['cmapdiff']
            norm = cp_info['normdiff']
        else:
            levels = cp_info['levels1']
            cmap = cp_info['cmap1']
            norm = cp_info['norm1']
        cmap=cmap, norm=norm,levels=levels"""

        if 'units' in vres:
            units = vres['units']
        else:
            units = ''

        """
        minval = np.min([np.min(mseasons), np.min(oseasons)])
        maxval = np.max([np.max(mseasons), np.max(oseasons)])

        #Gather contour level data (if applicable)
        #if 'contour_levels_range' in vres:
        #    levs = vres['contour_levels_range']
        #else:
        #    levs = 20
        if 'diff_contour_range' in vres:
            diff_levs = vres['diff_contour_range']
            diff_levs = [float(x) for x in diff_levs]
            print(diff_levs,"\n")
            diff_levs = np.arange(*diff_levs)
        else:
            diff_levs = 20
        

        # extract any MPL kwargs that should be passed on:
        if 'mpl' in vres:
            #subplots_opt.update(kwargs['mpl'].get('subplots',{}))
            #contourf_opt.update(kwargs['mpl'].get('contourf',{}))
            #colorbar_opt.update(kwargs['mpl'].get('colorbar',{}))
            if vres['mpl'].get('colorbar',{}):
                if vres['mpl']['colorbar'].get('ticks',{}):
                    cbar_ticks = vres['mpl']['colorbar']['ticks']
        
        if 'colormap' in vres:
            cmap1 = vres['colormap']
        else:
            cmap1 = "RdYlBu_r"
        #End if

        if 'contour_levels' in vres:
            levels1 = vres['contour_levels']
            if ('non_linear' in vres) and (vres['non_linear']):
                cmap_obj = cm.get_cmap(cmap1)
                norm1 = mpl.colors.BoundaryNorm(levels1, cmap_obj.N)
            else:
                norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
        elif 'contour_levels_range' in vres:
            assert len(vres['contour_levels_range']) == 3, \
            "contour_levels_range must have exactly three entries: min, max, step"

            lev_range = [float(x) for x in vres['contour_levels_range']]

            levels1 = np.arange(*lev_range)
            
            if ('non_linear' in vres) and (vres['non_linear']):
                cmap_obj = cm.get_cmap(cmap1)
                norm1 = mpl.colors.BoundaryNorm(levels1, cmap_obj.N)
            else:
                norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
        else:
            levels1 = np.linspace(minval, maxval, 12)
            if ('non_linear' in vres) and (vres['non_linear']):
                cmap_obj = cm.get_cmap(cmap1)
                norm1 = mpl.colors.BoundaryNorm(levels1, cmap_obj.N)
            else:
                norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)
        #End if
        """
        
        #
        cp_info = pf.prep_contour_plot(mseasons, oseasons, dseasons, **vres)
        clevs = np.unique(np.array(cp_info['levels1']))
        print(clevs)
        norm = cp_info['norm1']
        cmap = cp_info['cmap1']
        print(cmap,"\n")


        levs_diff = np.unique(np.array(cp_info['levelsdiff']))



        

        # mesh for plots:
        lat = mseasons['zalat']
        lev = mseasons['lev']
        lats, levs = np.meshgrid(lat, lev)


        #Get axis number for variable
        axs_id = var_axs[var]
        #print(axs_id)

        #x_filtered = x[~np.isnan(y)]

        #Contour fill
        img0 = axs[axs_id,0].contourf(lats, levs,mseasons, levels=clevs, norm=norm, cmap=cmap)
        img1 = axs[axs_id,1].contourf(lats, levs,oseasons, levels=clevs, norm=norm, cmap=cmap)
            
        #Add contours for highlighting
        axs[axs_id,0].contour(lats,levs,mseasons,levels=clevs[::2], norm=norm, colors="k")
        axs[axs_id,1].contour(lats,levs,oseasons,levels=clevs[::2], norm=norm, colors="k")

        #Check if difference plot has contour levels, if not print notification
        if len(dseasons.lev) == 0:
            axs[axs_id,2].text(prop_x, prop_y, empty_message,
                               transform=axs[axs_id,2].transAxes, bbox=props)
        else:
            img2 = axs[axs_id,2].contourf(lats,levs,dseasons, cmap="BrBG",levels=levs_diff,norm=cp_info['normdiff'])#levels=levs_diff
            axs[axs_id,2].contour(lats,levs,dseasons, colors="k",levels=levs_diff[::2],norm=cp_info['normdiff'])#levels=diff_levs[::2]
            plt.colorbar(img2, ax=axs[axs_id,2], location='right',**cp_info['diff_colorbar_opt'])#**cp_info['diff_colorbar_opt']

        #Format y-axis
        for a in axs[axs_id,:]:
            a.set_yscale("log")
            a.set_ylim(axs[axs_id,2].get_ylim()[::-1])

        plt.colorbar(img0, ax=axs[axs_id,0], location='right',**cp_info['colorbar_opt'])
        plt.colorbar(img1, ax=axs[axs_id,1], location='right',**cp_info['colorbar_opt'])
        """















        # Generate zonal plot:
        fig, ax = plt.subplots(figsize=(10,10),nrows=3, constrained_layout=True, sharex=True, sharey=True,**cp_info['subplots_opt'])
        levs = np.unique(np.array(cp_info['levels1']))

        levs_diff = np.unique(np.array(cp_info['levelsdiff']))


        if len(levs) < 2:
            img0, = axs[axs_id,0].contourf(lats, levs, mseasons, ax=ax[0])
            axs[axs_id,1].text(0.4, 0.4, empty_message, transform=ax[0].transAxes, bbox=props)
            img1, = axs[axs_id,1].contourf(lats, levs, oseasons, ax=ax[1])
            axs[axs_id,1].text(0.4, 0.4, empty_message, transform=ax[1].transAxes, bbox=props)
        else:
            img0, = axs[axs_id,0].contourf(lats, levs, mseasons, ax=ax[0], levels=clevs, norm=norm, cmap=cmap,**cp_info['contourf_opt'])
            img1, = ax[1].contourf(lats, levs, oseasons, ax=ax[1], levels=clevs, norm=norm, cmap=cmap,**cp_info['contourf_opt'])
            #Add contours for highlighting
            axs[axs_id,0].contour(lats,levs,mseasons,levels=clevs[::2], norm=norm, colors="k")
            axsaxs[axs_id,1].contour(lats,levs,oseasons,levels=clevs[::2], norm=norm, colors="k")
            fig.colorbar(img0, ax=axs[axs_id,1], location='right',**cp_info['colorbar_opt'])
            fig.colorbar(img1, ax=axs[axs_id,1], location='right',**cp_info['colorbar_opt'])
        #End if

        if len(levs_diff) < 2:
            img2, = axs[axs_id,2].contourf(lats, levs, dseasons, ax=ax[2])
            axs[axs_id,2].text(0.4, 0.4, empty_message, transform=ax[2].transAxes, bbox=props)
        else:
            img2, axs[axs_id,2].contourf(lats, levs, dseasons, norm=cp_info['normdiff'],
                                 cmap=cp_info['cmapdiff'],levels=cp_info['levelsdiff'],
                                 **cp_info['contourf_opt'])
            ax[2].contourf(lats, levs, dseasons, norm=cp_info['normdiff'],
                           colors="k",levels=cp_info['levelsdiff'],**cp_info['contourf_opt'])
            #img2, ax[2].contourf(lats, levs, dseasons, ax=ax[2], )
            #ax[2].contourf(lats, levs, dseasons, ax=ax[2],)
            
            
            #**cp_info['diff_colorbar_opt']
            fig.colorbar(img2, ax=axs[axs_id,2], location='right',**cp_info['diff_colorbar_opt'])
        tiFontSize = 10
        ax[0].set_title("TEST", loc='left', fontsize=tiFontSize)
        ax[1].set_title("BASE", loc='left', fontsize=tiFontSize)
        ax[2].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=tiFontSize)


        # style the plot:
        #Set Main title for subplots:
        #st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=15)
        #st.set_y(0.85)
        ax[-1].set_xlabel("LATITUDE")

        #if log_p:
        [a.set_yscale("log") for a in ax]

        fig.text(-0.03, 0.5, 'PRESSURE [hPa]', va='center', rotation='vertical')
        """




















        """
        #Extract variable of interest
        odata = oclim_ds[data_var].squeeze()  # squeeze in case of degenerate dimensions
        mdata = mclim_ds[var].squeeze()

        # APPLY UNITS TRANSFORMATION IF SPECIFIED:
        # NOTE: looks like our climo files don't have all their metadata
        mdata = mdata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
        # update units
        mdata.attrs['units'] = vres.get("new_unit", mdata.attrs.get('units', 'none'))

        # Do the same for the baseline case if need be:
        if not adfobj.compare_obs:
            odata = odata * vres.get("scale_factor",1) + vres.get("add_offset", 0)
            # update units
            odata.attrs['units'] = vres.get("new_unit", odata.attrs.get('units', 'none'))
        # Or for observations
        else:
            odata = odata * vres.get("obs_scale_factor",1) + vres.get("obs_add_offset", 0)
            # Note: we are going to assume that the specification ensures the conversion makes the units the same. Doesn't make sense to add a different unit.

         # determine whether it's 2D or 3D
        # 3D triggers search for surface pressure
        has_lat, has_lev = pf.zm_validate_dims(mdata)  # assumes will work for both mdata & odata

        #Notify user of level dimension:
        if has_lev:
            print(f"\t   {var} has lev dimension.")

        #
        # Seasonal Averages
        #

        #Create new dictionaries:
        mseasons = {}
        oseasons = {}

        #Loop over season dictionary:
        for s in seasons:
                    
            # time to make plot; here we'd probably loop over whatever plots we want for this variable
            # I'll just call this one "Zonal_Mean"  ... would this work as a pattern [operation]_[AxesDescription] ?
            # NOTE: Up to this point, nothing really differs from global_latlon_map,
            #       so we could have made one script instead of two.
            #       Merging would make overall timing better because looping twice will double I/O steps.
            #
            plot_name = plot_loc / f"{var}_{s}_Zonal_Mean.{plot_type}"

            if plot_name not in zonal_skip:

                #Seasonal Averages
                mseasons[s] = pf.seasonal_mean(mdata, season=s, is_climo=True)
                oseasons[s] = pf.seasonal_mean(odata, season=s, is_climo=True)

                # difference: each entry should be (lat, lon) or (plev, lat, lon)
                # dseasons[s] = mseasons[s] - oseasons[s]
                # difference will be calculated in plot_zonal_mean_and_save;
                # because we can let any pressure-level interpolation happen there
                # This could be re-visited for efficiency or improved code structure.

                #Create new plot with log-p:
                if has_lev:
                    plot_name_log = plot_loc / f"{var}_{s}_Zonal_logp_Mean.{plot_type}"

                    if plot_name_log not in logp_zonal_skip:
                        pf.plot_zonal_mean_and_save(plot_name_log, case_nickname, base_nickname,
                                                            [syear_cases[case_idx],eyear_cases[case_idx]],
                                                            [syear_baseline,eyear_baseline],
                                                            mseasons[s], oseasons[s], has_lev, log_p=True, obs=obs, **vres)

                        #Add plot to website (if enabled):
                        adfobj.add_website_data(plot_name_log, f"{var}_logp", case_name, season=s, plot_type="Zonal", category="Log-P")

        """








        """
        # Zonal mean zonal wind
        #------------------------------------------------------------------------------------------
        if var == "uzm":    
            mseasons.plot(ax=axs[0,0], y='lev', yscale='log',ylim=[1e3,1],levels=levels1,
                          norm=norm1,
                                    cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[0,0], levels = levels1[::2], y='lev', yscale='log',
                                                ylim=[1e3,1],norm=norm1,
                                                colors='black', linestyles=None)

            oseasons.plot(ax=axs[0,1], y='lev', yscale='log',ylim=[1e3,1],levels=levels1,
                                    norm=norm1,cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[0,1], levels = levels1[::2], y='lev', yscale='log',
                                                ylim=[1e3,1],norm=norm1,
                                                colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[0,2].text(prop_x, prop_y, empty_message, transform=axs[0,2].transAxes, bbox=props)
            else:
                dseasons.plot.contourf(ax=axs[0,2], y='lev', yscale='log', ylim=[1e3,1],cmap="BrBG",
                                      levels=diff_levs, cbar_kwargs={'label': units})

        # Zonal mean temperature
        #------------------------------------------------------------------------------------------
        if var == "thzm":
            #mseasons = np.log(mseasons)
            #oseasons = np.log(oseasons)

            # Plot the logarithmic data
            #log_temperature.plot()
            mseasons.plot(ax=axs[1,0], y='lev', yscale='log',ylim=[1e3,1],levels=np.arange(260,500,10),
                                    cbar_kwargs={'label': units})

            oseasons.plot(ax=axs[1,1], y='lev', yscale='log',ylim=[1e3,1],levels=np.arange(260,500,10),
                                    cbar_kwargs={'label': units})

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[1,2].text(prop_x, prop_y, empty_message, transform=axs[1,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[1,2], y='lev', yscale='log', ylim=[1e3,1],cmap="BrBG",levels=11,
                                    cbar_kwargs={'label': units})

        # EP Flux - meridional component
        #------------------------------------------------------------------------------------------
        if var == "epfy":
            mseasons.plot(ax=axs[2,0], y='lev', yscale='log',ylim=[1e2,1],levels=levels1,
                                    norm=norm1,cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[2,0], y='lev', yscale='log',norm=norm1,
                                                ylim=[1e2,1],levels=levels1[::2],
                                                colors='black', linestyles=None)

            oseasons.plot(ax=axs[2,1], y='lev', yscale='log',ylim=[1e2,1],levels=levels1,
                                    norm=norm1,cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[2,1], y='lev', yscale='log',norm=norm1,
                                                ylim=[1e2,1],levels=levels1[::2],
                                                colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[2,2].text(prop_x, prop_y, empty_message, transform=axs[2,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[2,2], y='lev', yscale='log',
                            ylim=[1e2,1],cmap="BrBG",levels=diff_levs,
                                    cbar_kwargs={'label': units})
        
        # EP Flux - vertical component vmax=1e5
        #------------------------------------------------------------------------------------------
        if var == "epfz":
            mseasons.plot(ax=axs[3,0], y='lev', yscale='log',ylim=[1e2,1],levels=levels1,vmax=3e4,
                                    norm=norm1,cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[3,0],  y='lev', yscale='log',norm=norm1,
                                                levels = levels1[::2],ylim=[1e2,1],
                                                colors='black', linestyles=None)

            oseasons.plot(ax=axs[3,1], y='lev', yscale='log',ylim=[1e2,1],levels=levels1,vmax=3e4,
                                    norm=norm1,cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[3,1],  y='lev', yscale='log',norm=norm1,
                                                levels = levels1[::2],ylim=[1e2,1],
                                                colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[3,2].text(prop_x, prop_y, empty_message, transform=axs[3,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[3,2], y='lev', yscale='log',
                            ylim=[1e2,1],cmap="BrBG",levels=diff_levs,
                                    cbar_kwargs={'label': units})

        # TEM meridional wind 
        #------------------------------------------------------------------------------------------
        if var == "vtem":
            mseasons.plot.contourf(ax=axs[4,0], levels = levels1,y='lev', yscale='log',
                                                ylim=[1e2,1],norm=norm1,
                                                cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[4,0],  y='lev', yscale='log',levels = levels1[::2],
                                                ylim=[1e2,1],norm=norm1,
                                                colors='black', linestyles=None)

            oseasons.plot.contourf(ax=axs[4,1], levels = levels1, y='lev', yscale='log',
                                                ylim=[1e2,1],norm=norm1,
                                                cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[4,1],  y='lev', yscale='log',levels = levels1[::2],
                                               ylim=[1e2,1],norm=norm1,
                                                colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[4,2].text(prop_x, prop_y, empty_message, transform=axs[4,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[4,2], y='lev', yscale='log',
                            ylim=[1e2,1],cmap="BrBG",levels=diff_levs,
                                    cbar_kwargs={'label': units})

        # TEM vertical wind
        #------------------------------------------------------------------------------------------
        if var == "wtem":
            mseasons = mseasons*100
            oseasons = oseasons*100

            #Contour fill
            img0 = axs[5,0].contourf(lats,levs,mseasons,levels=levels1, norm=norm1,cmap=cmap1)
            img1 = axs[5,1].contourf(lats,levs,oseasons,levels=levels1, norm=norm1,cmap=cmap1)
            
            #Add contours for highlighting
            axs[5,0].contour(lats,levs,mseasons,levels=levels1[::2], norm=norm1,colors="k")
            axs[5,1].contour(lats,levs,oseasons,levels = levels1[::2], norm=norm1,colors="k")

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[5,2].text(prop_x, prop_y, empty_message, transform=axs[3,2].transAxes, bbox=props)
            else:
                img2 = axs[5,2].contourf(lats,levs,dseasons, cmap="BrBG",levels=diff_levs)
                axs[5,2].contour(lats,levs,dseasons, colors="k",)#levels=diff_levs[::2]

            for a in axs[5,:]:
                a.set_yscale("log")
                a.set_ylim(axs[5,2].get_ylim()[::-1])

            plt.colorbar(img0, ax=axs[5,0], location='right',ticks=cbar_ticks)
            plt.colorbar(img1, ax=axs[5,1], location='right',ticks=cbar_ticks)
            plt.colorbar(img2, ax=axs[5,2], location='right',)

        # TEM mass stream function
        #------------------------------------------------------------------------------------------
        if var == "psitem":
            mseasons.plot.contourf(ax=axs[6,0], levels = levels1, y='lev', yscale='log',
                                                 ylim=[1e2,2],norm=norm1,
                                                cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[6,0], y='lev', yscale='log',norm=norm1,
                                            ylim=[1e2,2],levels=levels1[::2],
                                            colors='black', linestyles=None)

            oseasons.plot.contourf(ax=axs[6,1], levels = levels1, y='lev', yscale='log',
                                                 ylim=[1e2,2],norm=norm1,
                                                cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[6,1], y='lev', yscale='log',norm=norm1,
                                             ylim=[1e2,2],levels=levels1[::2],
                                            colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[6,2].text(prop_x, prop_y, empty_message, transform=axs[6,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[6,2], y='lev', yscale='log',
                                    ylim=[1e2,2],cmap="BrBG",levels=diff_levs,
                                    cbar_kwargs={'label': units})

        # EP flux divergence
        #------------------------------------------------------------------------------------------
        if var == "utendepfd":
            mseasons.plot.contourf(ax=axs[7,0], y='lev', yscale='log',levels=levels1,
                                             ylim=[1e2,2],norm=norm1,
                                            cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[7,0], y='lev', yscale='log',norm=norm1,
                                            ylim=[1e2,2],levels=levels1[::2],
                                            colors='black', linestyles=None)

            oseasons.plot.contourf(ax=axs[7,1], y='lev', yscale='log',levels=levels1,
                                             ylim=[1e2,2],norm=norm1,
                                            cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[7,1], y='lev', yscale='log',norm=norm1,
                                            ylim=[1e2,2],levels=levels1[::2],
                                            colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[7,2].text(prop_x, prop_y, empty_message, transform=axs[7,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[7,2], y='lev', yscale='log',
                                    ylim=[1e2,2],cmap="BrBG",levels=diff_levs,
                                    cbar_kwargs={'label': units})

        # EP flux divergence - meridional component
        #------------------------------------------------------------------------------------------
        if var == "utendvtem":
            mseasons.plot(ax=axs[8,0], y='lev', yscale='log',vmax=0.001, ylim=[1e3,1],
                                            cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[8,0],  y='lev', yscale='log',
                                            ylim=[1e3,2],levels=11,
                                            colors='black', linestyles=None)

            oseasons.plot(ax=axs[8,1], y='lev', yscale='log',vmax=0.001, ylim=[1e3,1],
                                            cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[8,1],  y='lev', yscale='log',
                                            ylim=[1e3,2],levels=11,
                                            colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[8,2].text(prop_x, prop_y, empty_message, transform=axs[8,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[8,2], y='lev', yscale='log', ylim=[1e3,1],cmap="BrBG",levels=11,
                                    cbar_kwargs={'label': units})

        # EP flux divergence - vertical component
        #------------------------------------------------------------------------------------------
        if var == "utendwtem":
            mseasons.plot(ax=axs[9,0], y='lev', yscale='log',vmax=0.0001, ylim=[1e3,1],
                                            cbar_kwargs={'label': units},cmap=cmap1)
            mseasons.plot.contour(ax=axs[9,0],  y='lev', yscale='log',
                                            ylim=[1e3,1],levels=11,
                                            colors='black', linestyles=None)

            oseasons.plot(ax=axs[9,1], y='lev', yscale='log',vmax=0.0001, ylim=[1e3,1],
                                            cbar_kwargs={'label': units},cmap=cmap1)
            oseasons.plot.contour(ax=axs[9,1],  y='lev', yscale='log',
                                            ylim=[1e3,1],levels=11,
                                            colors='black', linestyles=None)

            #Check if difference plot has contour levels, if not print notification
            if len(dseasons.lev) == 0:
                axs[9,2].text(prop_x, prop_y, empty_message, transform=axs[9,2].transAxes, bbox=props)
            else:
                dseasons.plot(ax=axs[9,2], y='lev', yscale='log', ylim=[1e3,1],cmap="BrBG",levels=11,
                                    cbar_kwargs={'label': units})
        """
    # Set the ticks and ticklabels for all x-axes
    #NOTE: This has to come after all subplots have been done,
    #I am assuming this is because of the way xarray plots info automatically for labels and titles
    #This is to change the default xarray labels for each instance of the xarray plot method
    plt.setp(axs, xticks=np.arange(-80,81,20), xlabel='latitude', title="")

    #Set titles of subplots
    #Set case names in first subplot only (zonal mean zonal wind)    
    longname = res["uzm"]["long_name"]

    test_yrs = f"{climo_yrs['test'][0]}-{climo_yrs['test'][1]}"
    axs[0,0].set_title(f"\n\n"+"$\mathbf{Test}$"+f"  yrs: {test_yrs}\n"+f"{case_names[0]}\n\n\n",fontsize=14)

    if obs:
        obs_title = Path(vres["obs_name"]).stem
        axs[0,1].set_title(f"\n\n"+"$\mathbf{Baseline}$\n"+f"{obs_title}\n\n"+longname+"\n",fontsize=14)

    else:
        base_yrs = f"{climo_yrs['base'][0]}-{climo_yrs['base'][1]}"
        axs[0,1].set_title(f"\n\n"+"$\mathbf{Baseline}$"+f"  yrs: {base_yrs}\n"+f"{case_names[1]}\n\n"+longname+"\n",fontsize=14)
    
    #Set main title for difference plots column
    axs[0,2].set_title("$\mathbf{Test} - \mathbf{Baseline}$"+"\n\n\n",fontsize=14)
    
    #Set variable name on center plot (except first plot, see above)
    for i in range(1,len(var_list)):
        vres = res[var_list[i]]

        #Variable plot title name
        longname = vres["long_name"]
        axs[i,1].set_title(longname+"\n",fontsize=14)
    
    #Adjust subplots
    #May need to adjust hspace and wspace depending on if multi-case diagnostics ever happen for TEM diags
    hspace = 0.4
    plt.subplots_adjust(wspace=0.5, hspace=hspace)

    return axs

##############
#END OF SCRIPT