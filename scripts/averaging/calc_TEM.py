import xarray as xr
import numpy as np
from scipy import integrate
from numpy import ma
from datetime import date
import matplotlib.pyplot as plt
from pathlib import Path
from glob import glob


def calc_TEM(adf):
    overwrite_output = True

    #Special ADF variables
    #CAM simulation variables (these quantities are always lists):
    case_names    = adf.get_cam_info("cam_case_name", required=True)    

    #Use test case settings, which are already lists:
    case_names    = adf.get_cam_info("cam_case_name", required=True)
    cam_hist_locs = adf.get_cam_info("cam_hist_loc", required=True)

    #New TEM netCDF file save location
    output_loc = adf.get_basic_info("tem_loc")
    
    #If path not specified, skip TEM calculation?
    if output_loc is None:
        return
    else:
        #Notify user that script has started:
        print("\n  Generating CAM TEM diagnostics files...")
    
    start_years   = adf.climo_yrs["syears"]
    end_years     = adf.climo_yrs["eyears"]

    #Set default to h4
    #TODO: Read this history file number from the yaml file?
    hist_num = "h4"

    #Loop over cases:
    for case_idx, case_name in enumerate(case_names):

        print(f"\t Processing TEM diagnostics for case '{case_name}' :")

        #Extract start and end year values:
        start_year = start_years[case_idx]
        end_year   = end_years[case_idx]

        #Create path object for the CAM history file(s) location:
        starting_location = Path(cam_hist_locs[case_idx])

        #Check that path actually exists:
        if not starting_location.is_dir():
            emsg = "Provided 'cam_hist_loc' directory '{starting_location}' not found."
            emsg += " Script is ending here."

            adf.end_diag_fail(emsg)
        #End if

        #Check if history files actually exist. If not then kill script:
        hist_str = '*.cam.'+hist_num
        if not list(starting_location.glob(hist_str+'.*.nc')):
            emsg = f"No CAM history {hist_str} files found in '{starting_location}'."
            emsg += " Script is ending here."
            adf.end_diag_fail(emsg)
        #End if

        #Check if TEM file already exists and over-writing is allowed:
        if Path(output_loc).is_file() and overwrite_output:
            #If so, then delete current file:
            output_loc.unlink()
        #End if

        # open input files
        shist = glob(f"{starting_location}/*h4.{start_year}*.nc")
        ehist = glob(f"{starting_location}/*h4.{end_year}*.nc")
        hist_files = sorted(shist + ehist)

        ds = xr.open_mfdataset(hist_files)

        print("ds.time.values",ds.time.values,"\n")

        #iterate over the times in a dataset
        for idx,_ in enumerate(ds.time.values):
            if idx == 0:
                dstem0 = calc_tem(ds.squeeze().isel(time=idx))
            else:
                dstem = calc_tem(ds.squeeze().isel(time=idx))
                dstem0 = xr.concat([dstem0, dstem],'time')
            #End if
        #End if    

        #Update the attributes
        dstem0.attrs = ds.attrs
        dstem0.attrs['created'] = str(date.today())
        dstem0['lev']=ds['lev']

        output_loc_idx = Path(output_loc) / case_name
        #Check if re-gridded directory exists, and if not, then create it:
        if not output_loc_idx.is_dir():
            print(f"    {output_loc_idx} not found, making new directory")
            output_loc_idx.mkdir(parents=True)
        #End if

        # write output to a netcdf file
        dstem0.to_netcdf(output_loc_idx / f'{case_name}.TEMdiag.nc', 
                            unlimited_dims='time', 
                            mode = 'w' )

    #Notify user that script has ended:
    print("  ...TEM diagnostics have been calculated successfully.")




def calc_tem(ds):
    """
    # calc_tem() function to calculate TEM diagnostics on CAM/WACCM output
    # This assumes the data have already been organized into zonal mean fluxes
    # Uzm, THzm, VTHzm, Vzm, UVzm, UWzm, Wzm
    # note that caculations are performed on model interface levels, which is ok
    # in the stratosphere but not in the troposphere.  If interested in tropospheric
    # TEM diagnostics, make sure input fields have been interpolated to true pressure levels.

    # The code follows the 'TEM recipe' from Appendix A of Gerber, E. P. and Manzini, E.:
    # The Dynamics and Variability Model Intercomparison Project (DynVarMIP) for CMIP6:
    # assessing the stratosphere–troposphere system, Geosci. Model Dev., 9, 3413–3425, 
    # https://doi.org/10.5194/gmd-9-3413-2016, 2016 and the corrigendum.

    # pdf available here: https://gmd.copernicus.org/articles/9/3413/2016/gmd-9-3413-2016.pdf, 
    # https://gmd.copernicus.org/articles/9/3413/2016/gmd-9-3413-2016-corrigendum.pdf

    # Output from post-processing function

    # Table A1. Momentum budget variable list (2-D monthly / daily zonal means, YZT).

    # Name      Long name [unit]

    # epfy      northward component of the Eliassen–Palm flux [m3 s−2]
    # epfz      upward component of the Eliassen–Palm flux [m3 s−2]
    # vtem      Transformed Eulerian mean northward wind [m s−1] 
    # wtem      Transformed Eulerian mean upward wind [m s−1]
    # psitem    Transformed Eulerian mean mass stream function [kg s−1]
    # utendepfd tendency of eastward wind due to Eliassen–Palm flux divergence [m s−2]
    # utendvtem tendency of eastward wind due to TEM northward wind advection and the Coriolis term [m s−2] 
    # utendwtem tendency of eastward wind due to TEM upward wind advection [m s−2]

    # this utility based on python code developed by Isla Simpson 25 Feb 2021
    # initial coding of stand alone function by Dan Marsh 16 Dec 2022

    # NOTE: function expects an xarray dataset with dataarrays of dimension (nlev,nlat)
    # to process more than one timestep interate over time. See calcTEM.ipynb notebook
    # for an example of processed a file or files with more than one timestep.
    """

    # constants for TEM calculations
    p0 = 101325. 
    a = 6.371e6 
    om = 7.29212e-5
    H = 7000.
    g0 = 9.80665

    nlat = ds['zalat'].size
    nlev = ds['lev'].size

    latrad = np.radians(ds.zalat)
    coslat = np.cos(latrad)
    coslat2d = np.tile(coslat,(nlev,1))
    
    pre = ds['lev']*100. # pressure levels in Pascals
    f = 2.*om*np.sin(latrad[:])
    f2d = np.tile(f,(nlev,1))
    
    # change missing values to NaNs
    uzm = ds['Uzm']
    uzm.values = ma.masked_greater_equal(uzm, 1e33)
    vzm = ds['Vzm']
    vzm.values = ma.masked_greater_equal(vzm, 1e33)
    wzm = ds['Wzm']
    wzm.values = ma.masked_greater_equal(wzm, 1e33)
    thzm = ds['THzm']
    thzm.values = ma.masked_greater_equal(thzm, 1e33)

    uvzm = ds['UVzm']
    uvzm.values = ma.masked_greater_equal(uvzm, 1e33)
    uwzm = ds['UWzm']
    uwzm.values = ma.masked_greater_equal(uwzm, 1e33)
    vthzm = ds['VTHzm']
    vthzm.values = ma.masked_greater_equal(vthzm, 1e33)
    
    # convert w terms from m/s to Pa/s
    wzm  = -1.*wzm*pre/H
    uwzm = -1.*uwzm*pre/H

    # compute the latitudinal gradient of U
    dudphi = (1./a)*np.gradient(uzm*coslat2d, 
                                latrad, 
                                axis=1)
    
    # compute the vertical gradient of theta and u
    dthdp = np.gradient(thzm, 
                        pre, 
                        axis=0)
    
    dudp = np.gradient(uzm,
                       pre,
                       axis=0)

    # compute eddy streamfunction and its vertical gradient
    psieddy = vthzm/dthdp
    dpsidp = np.gradient(psieddy,
                         pre,
                         axis=0)

    # (1/acos(phii))**d(psi*cosphi/dphi) for getting w*
    dpsidy = (1./(a*coslat2d)) \
           * np.gradient(psieddy*coslat2d,
                         latrad, 
                         axis=1)

    # TEM vertical velocity (Eq A7 of dynvarmip)
    wtem = wzm+dpsidy    
    
    # utendwtem (Eq A10 of dynvarmip)
    utendwtem = -1.*wtem*dudp

    # vtem (Eq A6 of dynvarmip)
    vtem = vzm-dpsidp
    
    # utendvtem (Eq A9 of dynvarmip)
    utendvtem = vtem*(f2d - dudphi) 

    # calculate E-P fluxes
    epfy = a*coslat2d*(dudp*psieddy - uvzm) # Eq A2
    epfz = a*coslat2d*((f2d-dudphi)*psieddy - uwzm) # Eq A3

    # calculate E-P flux divergence and zonal wind tendency 
    # due to resolved waves (Eq A5)
    depfydphi = (1./(a*coslat2d)) \
              * np.gradient(epfy*coslat2d,
                            latrad, 
                            axis=1)
        
    depfzdp = np.gradient(epfz,
                          pre,
                          axis=0)
    
    utendepfd = (depfydphi + depfzdp)/(a*coslat2d)
    utendepfd = xr.DataArray(utendepfd, coords = ds.Uzm.coords, name='utendepfd')
                             
    # TEM stream function, Eq A8
    topvzm = np.zeros([1,nlat])
    vzmwithzero = np.concatenate((topvzm, vzm), axis=0)
    prewithzero = np.concatenate((np.zeros([1]), pre))
    intv = integrate.cumtrapz(vzmwithzero,prewithzero,axis=0)
    psitem = (2*np.pi*a*coslat2d/g0)*(intv - psieddy)
   
    # final scaling of E-P fluxes and divergence to transform to log-pressure
    epfy = epfy*pre/p0      # A13
    epfz = -1.*(H/p0)*epfz  # A14
    wtem = -1.*(H/pre)*wtem # A16

    # 
    # add long name and unit attributes to TEM diagnostics
    epfy.attrs['long_name'] = 'northward component of E-P flux'
    epfy.attrs['units'] = 'm3/s2'
    
    epfz.attrs['long_name'] = 'upward component of E-P flux'
    epfz.attrs['units'] = 'm3/s2'

    vtem.attrs['long_name'] = 'Transformed Eulerian mean northward wind'
    vtem.attrs['units'] = 'm/s'
    
    wtem.attrs['long_name'] = 'Transformed Eulerian mean upward wind'
    wtem.attrs['units'] = 'm/s'
    
    psitem.attrs['long_name'] = 'Transformed Eulerian mean mass stream function'
    psitem.attrs['units'] = 'kg/s'
    
    utendepfd.attrs['long_name'] = 'tendency of eastward wind due to Eliassen-Palm flux divergence'
    utendepfd.attrs['units'] = 'm/s2'
    
    utendvtem.attrs['long_name'] = 'tendency of eastward wind due to TEM northward wind advection and the coriolis term'
    utendvtem.attrs['units'] = 'm/s2'
    
    utendwtem.attrs['long_name'] = 'tendency of eastward wind due to TEM upward wind advection'
    utendwtem.attrs['units'] = 'm/s2'
 
    epfy.values = np.float32(epfy.values)
    epfz.values = np.float32(epfz.values)
    wtem.values = np.float32(wtem.values)
    psitem.values = np.float32(psitem.values)
    utendepfd.values = np.float32(utendepfd.values)
    utendvtem.values = np.float32(utendvtem.values)
    utendwtem.values = np.float32(utendwtem.values)

    dstem = xr.Dataset(data_vars=dict(date = ds.date,
                                      datesec = ds.datesec,
                                      time_bnds = ds.time_bnds,
                                      uzm = uzm,
                                      vzm = vzm, 
                                      epfy = epfy,
                                      epfz = epfz,
                                      vtem = vtem,
                                      wtem = wtem,
                                      psitem = psitem,
                                      utendepfd = utendepfd,
                                      utendvtem = utendvtem,
                                      utendwtem = utendwtem
                                      )) 
    
    return dstem