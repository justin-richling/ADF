def plot_unstructured_map_and_save(wks, case_nickname, base_nickname,
                                   case_climo_yrs, baseline_climo_yrs,
                                   mdlfld, obsfld, diffld, pctld, wgt,
                                   obs=False, projection='global',**kwargs):

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
    mdlfld : uxarray.DataArray
        input data for case, needs units and long name attrubutes
    obsfld : uxarray.DataArray
        input data for base case, needs units and long name attrubutes 
    diffld : uxarray.DataArray
        input difference data, needs units and long name attrubutes
    pctld : uxarray.DataArray
        input percent difference data, needs units and long name attrubutes
    wgt : uxarray.DataArray
        weights assumed to be (area*landfrac)/(area*landfrac).sum()
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
    
    # prepare info for plotting
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
        
    if 'tiFontSize' in kwargs:
        tiFontSize = kwargs.pop('tiFontSize')
    else:
        tiFontSize = 8
        
    #generate a dictionary of contour plot settings:
    cp_info = prep_contour_plot(mdlfld, obsfld, diffld, pctld, **kwargs)
    
    if projection == 'global':
        transform = ccrs.PlateCarree()
        proj = ccrs.PlateCarree()
        figsize= (14, 7)
    elif projection == 'arctic':
        transform = ccrs.NorthPolarStereo()
        proj = ccrs.NorthPolarStereo()
        figsize = (8, 8)
        
    #nice formatting for tick labels
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    lon_formatter = LongitudeFormatter(number_format='0.0f',
                                    degree_symbol='',
                                    dateline_direction_label=False)
    lat_formatter = LatitudeFormatter(number_format='0.0f',
                                  degree_symbol='')

    # create figure object
    fig, axs = plt.subplots(2,2,
        figsize=figsize,
        facecolor="w",
        constrained_layout=True,
        subplot_kw=dict(projection=proj),
        **cp_info['subplots_opt']
    )
    axs=axs.flatten()
    
    # Loop over data arrays to make plots
    for i, a in enumerate(wrap_fields):
        if i == len(wrap_fields)-2:
            levels = cp_info['levelsdiff']
            cmap = cp_info['cmapdiff']
            norm = cp_info['normdiff']
        elif i == len(wrap_fields)-1:
            levels = cp_info['levelspctdiff']
            cmap = cp_info['cmappct']
            norm = cp_info['pctnorm']
        else:
            levels = cp_info['levels1']
            cmap = cp_info['cmap1']
            norm = cp_info['norm1']
    
        levs = np.unique(np.array(levels))
    
        #configure for polycollection plotting
        #TODO, would be nice to have levels set from the info, above
        ac = a.to_polycollection(projection=proj)
        #ac.norm(norm)
        ac.set_cmap(cmap)
        ac.set_antialiased(False)
        ac.set_transform(transform)
        ac.set_clim(vmin=levels[0],vmax=levels[-1])
        axs[i].add_collection(ac)
        if i > 0:
            cbar = plt.colorbar(ac, ax=axs[i], orientation='vertical', 
                                pad=0.05, shrink=0.8, **cp_info['colorbar_opt'])
            #TODO keep variable attributes on dataarrays
            #cbar.set_label(wrap_fields[i].attrs['units'])
        #Set stats: area_avg
        axs[i].set_title(f"Mean: {area_avg[i].item():5.2f}\nMax: {wrap_fields[i].max().item():5.2f}\nMin: {wrap_fields[i].min().item():5.2f}", 
                     loc='right', fontsize=tiFontSize)
   
    # Custom setting for each subplot
    for a in axs:
        a.coastlines()
        if projection=='global':
            a.set_global()
            a.spines['geo'].set_linewidth(1.5) #cartopy's recommended method
            a.set_xticks(np.linspace(-180, 120, 6), crs=proj)
            a.set_yticks(np.linspace(-90, 90, 7), crs=proj)
            a.tick_params('both', length=5, width=1.5, which='major')
            a.tick_params('both', length=5, width=1.5, which='minor')
            a.xaxis.set_major_formatter(lon_formatter)
            a.yaxis.set_major_formatter(lat_formatter)
        elif projection == 'arctic':
            a.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
            # __Follow the cartopy gallery example to make circular__:
            # Compute a circle in axes coordinates, which we can use as a boundary
            # for the map. We can pan/zoom as much as we like - the boundary will be
            # permanently circular.
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpl.path.Path(verts * radius + center)
            a.set_boundary(circle, transform=a.transAxes)
            a.gridlines(draw_labels=False, crs=ccrs.PlateCarree(), 
                        lw=1, color="gray",y_inline=True, 
                        xlocs=range(-180,180,90), ylocs=range(0,90,10))
    
    st = fig.suptitle(wks.stem[:-5].replace("_"," - "), fontsize=18)
    st.set_y(0.85)

    #Set plot titles
    case_title = "$\mathbf{Test}:$"+f"{case_nickname}\nyears: {case_climo_yrs[0]}-{case_climo_yrs[-1]}"
    axs[0].set_title(case_title, loc='left', fontsize=tiFontSize)
    if obs:
        obs_var = kwargs["obs_var_name"]
        obs_title = kwargs["obs_file"][:-3]
        base_title = "$\mathbf{Baseline}:$"+obs_title+"\n"+"$\mathbf{Variable}:$"+f"{obs_var}"
        axs[1].set_title(base_title, loc='left', fontsize=tiFontSize)
    else:
        base_title = "$\mathbf{Baseline}:$"+f"{base_nickname}\nyears: {baseline_climo_yrs[0]}-{baseline_climo_yrs[-1]}"
        axs[1].set_title(base_title, loc='left', fontsize=tiFontSize)
    axs[2].set_title("$\mathbf{Test} - \mathbf{Baseline}$", loc='left', fontsize=tiFontSize)
    axs[2].set_title(f"RMSE: {d_rmse:.3f}", fontsize=tiFontSize)
    axs[3].set_title("Test % Diff Baseline", loc='left', fontsize=tiFontSize,fontweight="bold")
        
    fig.savefig(wks, bbox_inches='tight', dpi=300)
    
    #Close plots:
    plt.close()