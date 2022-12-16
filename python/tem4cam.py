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


import xarray as xr
import numpy as np
from scipy import integrate
from numpy import ma

def calc_tem(ds):

    # constants for TEM calculations
    p0 = 101325. 
    a = 6.371e6 
    om = 7.29212e-5
    H = 7000.
    g0 = 9.80665

    nlat = ds['lat'].size
    nlev = ds['ilev'].size

    latrad = np.radians(ds.lat)
    coslat = np.cos(latrad)
    coslat2d = np.tile(coslat,(nlev,1))
    
    pre = ds['ilev']*100. # pressure levels in Pascals
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

    # add long name and unit attributes to TEM diagnostics
    epfy.attrs['long_name'] = 'northward component of E-P flux'
    epfy.attrs['units'] = 'm3/s2'
    
    epfz.attrs['long_name'] = 'upward component of E-P flux'
    epfz.attrs['units'] = 'm2/s2'

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