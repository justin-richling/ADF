def make_polar_plot(wks, case_nickname, base_nickname,
                    case_climo_yrs, baseline_climo_yrs,
                    d1:xr.DataArray, d2:xr.DataArray, difference:Optional[xr.DataArray]=None,pctchange:Optional[xr.DataArray]=None,
                    domain:Optional[list]=None, hemisphere:Optional[str]=None, obs=False, **kwargs):

    """Make a stereographic polar plot for the given data and hemisphere.

    Parameters
    ----------
    wks : str or Path
        output file path
    case_nickname : str
        short case name for `d1`
    base_nickname : str
        short case name for `d2`
    case_climo_yrs : list
        years for case `d1`, used for annotation
    baseline_climo_yrs : list
        years for case `d2`, used for annotation
    d1, d2 : xr.DataArray
        input data, must contain dimensions `lat` and `lon`
    difference : xr.DataArray, optional
        data to use as the difference, otherwise `d2 - d1`
    pctchange : xr.DataArray, optional data to use as the percent change
    domain : list, optional
        the domain to plot, specified as
        [west_longitude, east_longitude, south_latitude, north_latitude]
        Defaults to pole to 45-degrees, all longitudes
    hemisphere : {'NH', 'SH'}, optional
        Hemsiphere to plot
    kwargs : dict, optional
        variable-dependent options for plots, See Notes.

    Notes
    -----
    - Uses contourf. No contour lines (yet).
    - kwargs is checked for:
        + `colormap`
        + `contour_levels`
        + `contour_levels_range`
        + `diff_contour_levels`
        + `diff_contour_range`
        + `diff_colormap`
        + `units`
    """
    if difference is None:
        dif = d2 - d1
    else:
        dif = difference
        
    if  pctchange is None:
        pct = (d2 - d1) / np.abs(d1) * 100.0
    else:
        pct = pctchange
        
    #check if pct has NaN's or Inf values and if so set them to 0 to prevent plotting errors
    pct = pct.where(np.isfinite(pct), np.nan)
    pct = pct.fillna(0.0)

    if hemisphere.upper() == "NH":
        proj = ccrs.NorthPolarStereo()
    elif hemisphere.upper() == "SH":
        proj = ccrs.SouthPolarStereo()
    else:
        raise AdfError(f'[make_polar_plot] hemisphere not specified, must be NH or SH; hemisphere set as {hemisphere}')

    if domain is None:
        if hemisphere.upper() == "NH":
            domain = [-180, 180, 45, 90]
        else:
            domain = [-180, 180, -90, -45]

    # statistics for annotation (these are scalars):
    d1_region_mean, d1_region_max, d1_region_min = domain_stats(d1, domain)
    d2_region_mean, d2_region_max, d2_region_min = domain_stats(d2, domain)
    dif_region_mean, dif_region_max, dif_region_min = domain_stats(dif, domain)
    pct_region_mean, pct_region_max, pct_region_min = domain_stats(pct, domain)

    #downsize to the specified region; makes plotting/rendering/saving much faster
    d1 = d1.sel(lat=slice(domain[2],domain[3]))
    d2 = d2.sel(lat=slice(domain[2],domain[3]))
    dif = dif.sel(lat=slice(domain[2],domain[3]))
    pct = pct.sel(lat=slice(domain[2],domain[3]))

    # add cyclic point to the data for better-looking plot
    d1_cyclic, lon_cyclic = add_cyclic_point(d1, coord=d1.lon)
    d2_cyclic, _ = add_cyclic_point(d2, coord=d2.lon)  # since we can take difference, assume same longitude coord.
    dif_cyclic, _ = add_cyclic_point(dif, coord=dif.lon)
    pct_cyclic, _ = add_cyclic_point(pct, coord=pct.lon)

    # -- deal with optional plotting arguments that might provide variable-dependent choices

    # determine levels & color normalization:
    minval    = np.min([np.min(d1), np.min(d2)])
    maxval    = np.max([np.max(d1), np.max(d2)])
    absmaxdif = np.max(np.abs(dif))
    absmaxpct = np.max(np.abs(pct))

    if 'colormap' in kwargs:
        cmap1 = kwargs['colormap']
    else:
        cmap1 = 'coolwarm'

    if 'contour_levels' in kwargs:
        levels1 = kwargs['contour_levels']
        norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
    elif 'contour_levels_range' in kwargs:
        assert len(kwargs['contour_levels_range']) == 3, "contour_levels_range must have exactly three entries: min, max, step"
        levels1 = np.arange(*kwargs['contour_levels_range'])
        norm1 = mpl.colors.Normalize(vmin=min(levels1), vmax=max(levels1))
    else:
        levels1 = np.linspace(minval, maxval, 12)
        norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)

    if ('colormap' not in kwargs) and ('contour_levels' not in kwargs):
        norm1, cmap1 = get_difference_colors(levels1)  # maybe these are better defaults if nothing else is known.

    if "diff_contour_levels" in kwargs:
        levelsdiff = kwargs["diff_contour_levels"]  # a list of explicit contour levels
    elif "diff_contour_range" in kwargs:
            assert len(kwargs['diff_contour_range']) == 3, "diff_contour_range must have exactly three entries: min, max, step"
            levelsdiff = np.arange(*kwargs['diff_contour_range'])
    else:
        # set levels for difference plot (with a symmetric color bar):
        levelsdiff = np.linspace(-1*absmaxdif, absmaxdif, 12)
    #End if
    
    if "pct_diff_contour_levels" in kwargs:
        levelspctdiff = kwargs["pct_diff_contour_levels"]  # a list of explicit contour levels
    elif "pct_diff_contour_range" in kwargs:
            assert len(kwargs['pct_diff_contour_range']) == 3, "pct_diff_contour_range must have exactly three entries: min, max, step"
            levelspctdiff = np.arange(*kwargs['pct_diff_contour_range'])
    else:
        levelspctdiff = [-100,-75,-50,-40,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,40,50,75,100]
    pctnorm = mpl.colors.BoundaryNorm(levelspctdiff,256)

    #NOTE: Sometimes the contour levels chosen in the defaults file
    #can result in the "contourf" software stack generating a
    #'TypologyException', which should manifest itself as a
    #"PredicateError", but due to bugs in the stack itself
    #will also sometimes raise an AttributeError.

    #To prevent this from happening, the polar max and min values
    #are calculated, and if the default contour values are significantly
    #larger then the min-max values, then the min-max values are used instead:
    #-------------------------------
    if max(levels1) > 10*maxval:
        levels1 = np.linspace(minval, maxval, 12)
        norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    elif minval < 0 and min(levels1) < 10*minval:
        levels1 = np.linspace(minval, maxval, 12)
        norm1 = mpl.colors.Normalize(vmin=minval, vmax=maxval)
    #End if

    if max(np.abs(levelsdiff)) > 10*absmaxdif:
        levelsdiff = np.linspace(-1*absmaxdif, absmaxdif, 12)
    
    
    #End if
    #-------------------------------

    # Difference options -- Check in kwargs for colormap and levels
    if "diff_colormap" in kwargs:
        cmapdiff = kwargs["diff_colormap"]
        dnorm, _ = get_difference_colors(levelsdiff)  # color map output ignored
    else:
        dnorm, cmapdiff = get_difference_colors(levelsdiff)  
        
    # Pct Difference options -- Check in kwargs for colormap and levels
    if "pct_diff_colormap" in kwargs:
        cmappct = kwargs["pct_diff_colormap"]        
    else:
        cmappct = "PuOr_r"
    #End if

    # -- end options

    lons, lats = np.meshgrid(lon_cyclic, d1.lat)

    fig = plt.figure(figsize=(10,10))
    gs = mpl.gridspec.GridSpec(2, 4, wspace=0.9)

    ax1 = plt.subplot(gs[0, :2], projection=proj)
    ax2 = plt.subplot(gs[0, 2:], projection=proj)
    ax3 = plt.subplot(gs[1, :2], projection=proj)
    ax4 = plt.subplot(gs[1, 2:], projection=proj)

    levs = np.unique(np.array(levels1))
    levs_diff = np.unique(np.array(levelsdiff))
    levs_pctdiff = np.unique(np.array(levelspctdiff))

    if len(levs) < 2:
        img1 = ax1.contourf(lons, lats, d1_cyclic, transform=ccrs.PlateCarree(), colors="w", norm=norm1)
        ax1.text(0.4, 0.4, empty_message, transform=ax1.transAxes, bbox=props)

        img2 = ax2.contourf(lons, lats, d2_cyclic, transform=ccrs.PlateCarree(), colors="w", norm=norm1)
        ax2.text(0.4, 0.4, empty_message, transform=ax2.transAxes, bbox=props)
    else:
        img1 = ax1.contourf(lons, lats, d1_cyclic, transform=ccrs.PlateCarree(), cmap=cmap1, norm=norm1, levels=levels1)
        img2 = ax2.contourf(lons, lats, d2_cyclic, transform=ccrs.PlateCarree(), cmap=cmap1, norm=norm1, levels=levels1)

    if len(levs_pctdiff) < 2:
        img3 = ax3.contourf(lons, lats, pct_cyclic, transform=ccrs.PlateCarree(), colors="w", norm=pctnorm, transform_first=True)
        ax3.text(0.4, 0.4, empty_message, transform=ax3.transAxes, bbox=props)
    else:
        img3 = ax3.contourf(lons, lats, pct_cyclic, transform=ccrs.PlateCarree(), cmap=cmappct, norm=pctnorm, levels=levelspctdiff, transform_first=True)

    if len(levs_diff) < 2:
        img4 = ax4.contourf(lons, lats, dif_cyclic, transform=ccrs.PlateCarree(), colors="w", norm=dnorm)
        ax4.text(0.4, 0.4, empty_message, transform=ax4.transAxes, bbox=props)
    else:
        img4 = ax4.contourf(lons, lats, dif_cyclic, transform=ccrs.PlateCarree(), cmap=cmapdiff, norm=dnorm, levels=levelsdiff)
        
    #Set Main title for subplots:
    st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=18)
    st.set_y(0.95)

    #Set plot titles
    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"
    ax1.set_title(case_title, loc='left', fontsize=6) #fontsize=tiFontSize

    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
        ax2.set_title(base_title, loc='left', fontsize=6) #fontsize=tiFontSize
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"
        ax2.set_title(base_title, loc='left', fontsize=6)

    ax1.text(-0.2, -0.10, f"Mean: {d1_region_mean:5.2f}\nMax: {d1_region_max:5.2f}\nMin: {d1_region_min:5.2f}", transform=ax1.transAxes)

    ax2.text(-0.2, -0.10, f"Mean: {d2_region_mean:5.2f}\nMax: {d2_region_max:5.2f}\nMin: {d2_region_min:5.2f}", transform=ax2.transAxes)

    ax3.text(-0.2, -0.10, f"Mean: {pct_region_mean:5.2f}\nMax: {pct_region_max:5.2f}\nMin: {pct_region_min:5.2f}", transform=ax3.transAxes)
    ax3.set_title("Test % diff Baseline", loc='left', fontsize=8)

    ax4.text(-0.2, -0.10, f"Mean: {dif_region_mean:5.2f}\nMax: {dif_region_max:5.2f}\nMin: {dif_region_min:5.2f}", transform=ax4.transAxes)
    ax4.set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=8)

    if "units" in kwargs:
        ax2.set_ylabel(kwargs["units"])
        ax4.set_ylabel(kwargs["units"])
    else:
        ax2.set_ylabel(f"{d1.units}")
        ax4.set_ylabel(f"{d1.units}")

    [a.set_extent(domain, ccrs.PlateCarree()) for a in [ax1, ax2, ax3, ax4]]
    [a.coastlines() for a in [ax1, ax2, ax3, ax4]]

    # __Follow the cartopy gallery example to make circular__:
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpl.path.Path(verts * radius + center)
    [a.set_boundary(circle, transform=a.transAxes) for a in [ax1, ax2, ax3, ax4]]

    # __COLORBARS__
    cb_mean_ax = inset_axes(ax2,
                    width="5%",  # width = 5% of parent_bbox width
                    height="90%",  # height : 90%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0.05, 1, 1),
                    bbox_transform=ax2.transAxes,
                    borderpad=0,
                    )
    fig.colorbar(img1, cax=cb_mean_ax)
    
    cb_pct_ax = inset_axes(ax3,
                    width="5%",  # width = 5% of parent_bbox width
                    height="90%",  # height : 90%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0.05, 1, 1),
                    bbox_transform=ax3.transAxes,
                    borderpad=0,
                    )  

    cb_diff_ax = inset_axes(ax4,
                    width="5%",  # width = 5% of parent_bbox width
                    height="90%",  # height : 90%
                    loc='lower left',
                    bbox_to_anchor=(1.05, 0.05, 1, 1),
                    bbox_transform=ax4.transAxes,
                    borderpad=0,
                    )      
                    
    fig.colorbar(img3, cax=cb_pct_ax)
    
    fig.colorbar(img4, cax=cb_diff_ax)

    # Save files
    fig.savefig(wks, bbox_inches='tight', dpi=300)

    # Close figures to avoid memory issues:
    plt.close(fig)