def plot_map_and_save(wks, case_nickname, base_nickname,
                      case_climo_yrs, baseline_climo_yrs,
                      mdlfld, obsfld, diffld, pctld, unstructured=False,
                      obs=False, **kwargs):
    """This plots mdlfld, obsfld, diffld in a 3-row panel plot of maps.

    Parameters
    ----------
    wks : str or Path
        output file path
    case_nickname : str
        short name for case
    base_nickname : str
        short name for base case
    case_climo_yrs : list
        list of years in case climatology, used for annotation
    baseline_climo_yrs : list
        list of years in base case climatology, used for annotation
    mdlfld : xarray.DataArray
        input data for case
    obsfld : xarray.DataArray
        input data for base case
    diffld : xarray.DataArray
        input difference data
    pctld : xarray.DataArray
        input percent difference data
    kwargs : dict, optional
        variable-specific options, See Notes

    Notes
    -----
    kwargs expected to be a variable-specific section,
    possibly provided by an ADF Variable Defaults YAML file.
    Currently it is inspected for:
    - colormap -> str, name of matplotlib colormap
    - contour_levels -> list of explict values or a tuple: (min, max, step)
    - diff_colormap
    - diff_contour_levels
    - tiString -> str, Title String
    - tiFontSize -> int, Title Font Size
    - mpl -> dict, This should be any matplotlib kwargs that should be passed along. Keep reading:
        + Organize these by the mpl function. In this function (`plot_map_and_save`)
          we will check for an entry called `subplots`, `contourf`, and `colorbar`. So the YAML might looks something like:
          ```
           mpl:
             subplots:
               figsize: (3, 9)
             contourf:
               levels: 15
               cmap: Blues
             colorbar:
               shrink: 0.4
          ```
        + This is experimental, and if you find yourself doing much with this, you probably should write a new plotting script that does not rely on this module.
    When these are not provided, colormap is set to 'coolwarm' and limits/levels are set by data range.
    """

    #nice formatting for tick labels
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    # preprocess
    if not unstructured:
        # - assume all three fields have same lat/lon
        lat = obsfld['lat']
        wgt = np.cos(np.radians(lat))
        mwrap, lon = add_cyclic_point(mdlfld, coord=mdlfld['lon'])
        owrap, _ = add_cyclic_point(obsfld, coord=obsfld['lon'])
        dwrap, _ = add_cyclic_point(diffld, coord=diffld['lon'])
        pwrap, _ = add_cyclic_point(pctld, coord=pctld['lon'])
        wrap_fields = (mwrap, owrap, pwrap, dwrap)
        # mesh for plots:
        lons, lats = np.meshgrid(lon, lat)
        # Note: using wrapped data makes spurious lines across plot (maybe coordinate dependent)
        lon2, lat2 = np.meshgrid(mdlfld['lon'], mdlfld['lat'])

        # get statistics (from non-wrapped)
        fields = (mdlfld, obsfld, diffld, pctld)
        area_avg = [spatial_average(x, weights=wgt, spatial_dims=None) for x in fields]

        d_rmse = wgt_rmse(mdlfld, obsfld, wgt)  # correct weighted RMSE for (lat,lon) fields.
        # specify the central longitude for the plot
        central_longitude = kwargs.get('central_longitude', 180)
    else:
        wgt = kwargs["wgt"]
        wrap_fields = (mdlfld, obsfld, diffld, pctld)
        area_avg = [global_average(x, wgt) for x in wrap_fields]

        # TODO Check this is correct, weighted rmse uses xarray weighted function
        #d_rmse = wgt_rmse(a, b, wgt)  
        d_rmse = (np.sqrt(((diffld**2)*wgt).sum())).values.item()

    # We should think about how to do plot customization and defaults.
    # Here I'll just pop off a few custom ones, and then pass the rest into mpl.
    if 'tiString' in kwargs:
        tiString = kwargs.pop("tiString")
    else:
        tiString = ''
    #End if

    if 'tiFontSize' in kwargs:
        tiFontSize = kwargs.pop('tiFontSize')
    else:
        tiFontSize = 8
    #End if

    central_longitude = kwargs.get('central_longitude', 180)

    # generate dictionary of contour plot settings:
    cp_info = prep_contour_plot(mdlfld, obsfld, diffld, pctld, **kwargs)

    # create figure object
    fig = plt.figure(figsize=(14,10))

    # LAYOUT WITH GRIDSPEC
    gs = mpl.gridspec.GridSpec(3, 6, wspace=2.0,hspace=0.0) # 2 rows, 4 columns, but each map will take up 2 columns
    proj = ccrs.PlateCarree(central_longitude=central_longitude)
    ax1 = plt.subplot(gs[0:2, :3], projection=proj, **cp_info['subplots_opt'])
    ax2 = plt.subplot(gs[0:2, 3:], projection=proj, **cp_info['subplots_opt'])
    ax3 = plt.subplot(gs[2, :3], projection=proj, **cp_info['subplots_opt'])
    ax4 = plt.subplot(gs[2, 3:], projection=proj, **cp_info['subplots_opt'])
    ax = [ax1,ax2,ax3,ax4]

    img = [] # contour plots
    cs = []  # contour lines
    cb = []  # color bars

    # formatting for tick labels
    lon_formatter = LongitudeFormatter(number_format='0.0f',
                                        degree_symbol='',
                                        dateline_direction_label=False)
    lat_formatter = LatitudeFormatter(number_format='0.0f',
                                        degree_symbol='')

    for i, a in enumerate(wrap_fields):

        if i == len(wrap_fields)-1:
            levels = cp_info['levelsdiff']
            cmap = cp_info['cmapdiff']
            norm = cp_info['normdiff']
        elif i == len(wrap_fields)-2:
            levels = cp_info['levelspctdiff']
            cmap = cp_info['cmappct']
            norm = cp_info['pctnorm']
        else:
            levels = cp_info['levels1']
            cmap = cp_info['cmap1']
            norm = cp_info['norm1']
        
        # Unstructured grid check
        if not unstructured:
            levs = np.unique(np.array(levels))
            if len(levs) < 2:
                img.append(ax[i].contourf(lons,lats,a,colors="w",transform=ccrs.PlateCarree(),transform_first=True))
                ax[i].text(0.4, 0.4, empty_message, transform=ax[i].transAxes, bbox=props)
            else:
                img.append(ax[i].contourf(lons, lats, a, levels=levels, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), transform_first=True, **cp_info['contourf_opt']))
            #End if
        else:
            #configure for polycollection plotting
            #TODO, would be nice to have levels set from the info, above
            ac = a.to_polycollection(projection=proj)
            img.append(ac)
            #ac.norm(norm)
            ac.set_cmap(cmap)
            ac.set_antialiased(False)
            ac.set_transform(proj)
            ac.set_clim(vmin=levels[0],vmax=levels[-1])
            ax[i].add_collection(ac)
            if i > 0:
                cbar = plt.colorbar(ac, ax=ax[i], orientation='vertical', 
                                    pad=0.05, shrink=0.8, **cp_info['colorbar_opt'])
                #TODO keep variable attributes on dataarrays
                #cbar.set_label(wrap_fields[i].attrs['units'])
        # End if unstructured grid

        #ax[i].set_title("AVG: {0:.3f}".format(area_avg[i]), loc='right', fontsize=11)
        ax[i].set_title(f"Mean: {area_avg[i].item():5.2f}\nMax: {wrap_fields[i].max().item():5.2f}\nMin: {wrap_fields[i].min().item():5.2f}", 
                     loc='right', fontsize=tiFontSize)

        # add contour lines <- Unused for now -JN
        # TODO: add an option to turn this on -BM
        #cs.append(ax[i].contour(lon2, lat2, fields[i], transform=ccrs.PlateCarree(), colors='k', linewidths=1))
        #ax[i].clabel(cs[i], cs[i].levels, inline=True, fontsize=tiFontSize-2, fmt='%1.1f')
        #ax[i].text( 10, -140, "CONTOUR FROM {} to {} by {}".format(min(cs[i].levels), max(cs[i].levels), cs[i].levels[1]-cs[i].levels[0]),
        #bbox=dict(facecolor='none', edgecolor='black'), fontsize=tiFontSize-2)

    st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=18)
    st.set_y(0.85)

    #Set plot titles
    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"
    ax[0].set_title(case_title, loc='left', fontsize=tiFontSize)

    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"
        ax[1].set_title(base_title, loc='left', fontsize=tiFontSize)

    """
    #Set stats: area_avg
    ax[0].set_title(f"Mean: {mdlfld.weighted(wgt).mean().item():5.2f}\nMax: {mdlfld.max():5.2f}\nMin: {mdlfld.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    ax[1].set_title(f"Mean: {obsfld.weighted(wgt).mean().item():5.2f}\nMax: {obsfld.max():5.2f}\nMin: {obsfld.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    ax[2].set_title(f"Mean: {pctld.weighted(wgt).mean().item():5.2f}\nMax: {pctld.max():5.2f}\nMin: {pctld.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    ax[3].set_title(f"Mean: {diffld.weighted(wgt).mean().item():5.2f}\nMax: {diffld.max():5.2f}\nMin: {diffld.min():5.2f}", loc='right',
                       fontsize=tiFontSize)
    """
    # set rmse title:
    ax[3].set_title(f"RMSE: {d_rmse:.3f}", fontsize=tiFontSize)
    ax[3].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=tiFontSize)
    ax[2].set_title("Test % Diff Baseline", loc='left', fontsize=tiFontSize,fontweight="bold")

    for a in ax:
        a.spines['geo'].set_linewidth(1.5) #cartopy's recommended method
        a.coastlines()
        a.set_xticks(np.linspace(-180, 120, 6), crs=ccrs.PlateCarree())
        a.set_yticks(np.linspace(-90, 90, 7), crs=ccrs.PlateCarree())
        a.tick_params('both', length=5, width=1.5, which='major')
        a.tick_params('both', length=5, width=1.5, which='minor')
        a.xaxis.set_major_formatter(lon_formatter)
        a.yaxis.set_major_formatter(lat_formatter)

    # __COLORBARS__
    cb_mean_ax = inset_axes(ax2,
                    width="5%",  # width = 5% of parent_bbox width
                    height="100%",  # height : 100%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0, 1, 1),
                    bbox_transform=ax2.transAxes,
                    borderpad=0,
                    )
    fig.colorbar(img[1], cax=cb_mean_ax, **cp_info['colorbar_opt'])
    

    cb_pct_ax = inset_axes(ax3,
                    width="5%",  # width = 5% of parent_bbox width
                    height="100%",  # height : 100%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0, 1, 1),
                    bbox_transform=ax3.transAxes,
                    borderpad=0,
                    )
    PCT_CB = fig.colorbar(img[2], cax=cb_pct_ax, **cp_info['colorbar_opt'])
    PCT_CB.ax.set_ylabel="%"

    cb_diff_ax = inset_axes(ax4,
                    width="5%",  # width = 5% of parent_bbox width
                    height="100%",  # height : 100%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0, 1, 1),
                    bbox_transform=ax4.transAxes,
                    borderpad=0,
                    )
    fig.colorbar(img[3], cax=cb_diff_ax, **cp_info['colorbar_opt'])

    # Write final figure to file
    fig.savefig(wks, bbox_inches='tight', dpi=300)

    #Close plots:
    plt.close()