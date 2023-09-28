import numpy as np
import xarray as xr
import pandas as pd
import gzip
import struct
import cartopy.crs as ccrs
import re




def set_BOEASE2():
    '''
    Set up the EASE2 grid for Brighter Ocean 
    '''    
    NXBO = 263
    DXBO = 25e3
    XLIM_BO = NXBO*DXBO/2 #3287.5e3
    YLIM_BO = XLIM_BO #3287.5e3

    return set_ease(version=2,XLIM=XLIM_BO,YLIM=YLIM_BO)



def set_ease(
        XLIM = 5000e3,
        YLIM = 5000e3,
        version = 2,
        res = 25,
):
    '''
    Set EASE or EASE2 grid for the Arctic region
    '''
    lon0 = 0
    lat0 = 90
    pcproj = ccrs.PlateCarree()
    if version==2:
        easeprj = ccrs.LambertAzimuthalEqualArea(central_longitude=lon0,central_latitude=lat0)                      
        DX = res*1e3
        DY = res*1e3
    elif version==1:
        easeprj = ccrs.epsg('3408')
        if res==12.5:
            XLIM=9030574.08
            YLIM=9030574.08
            DX=12533.76
            DY=12533.76


    nx   = int(XLIM*2/DX)
    ny   = int(XLIM*2/DY)

    xgrid = np.linspace(-XLIM,XLIM,nx+1)
    ygrid = np.linspace(-YLIM,YLIM,ny+1)
        
    xc = (xgrid[1:]+xgrid[:-1])/2
    yc = (ygrid[1:]+ygrid[:-1])/2
    xcs,ycs = np.meshgrid(xc,yc)
    lonc, latc = ptproj(easeprj,pcproj, xcs,ycs)
     
    mds = xr.Dataset()
    mds['X'] = xc
    mds['Y'] = yc
    mds['lat'] = (['Y','X'],latc)
    mds['lon'] = (['Y','X'],lonc)
    

    mds.attrs['XLIM'  ] = XLIM
    mds.attrs['YLIM'  ] = YLIM
    mds.attrs['dx'    ] = DX
    mds.attrs['dy'    ] = DY
    mds.attrs['lat0'  ] = lat0
    mds.attrs['lon0'  ] = lon0
    mds.attrs['grid'  ] = f'EASE {version}'
    mds.attrs['res'] = f'{res} km'
    return mds, easeprj


def ptproj(prj_in,prj_out,x_in,y_in):
    '''
    Convert the coordinates of input projection to the output projection
    
    Parameters:
    ----------
    prj_in : the map projection of the input coordinates
    prj_out: the map projection of the output coordinates
    x_in   : x coordinates in prj_in
    y_in   : y coordinates in prj_in
    #list   : A list of points or a single point, boolean, default True.
    
    Returns:
    -------
    x_out: x coordinates in prj_out
    y_out: y coordinates in prj_out
    '''
    if isinstance(x_in,(list,pd.core.series.Series,np.ndarray)):
        pts = prj_out.transform_points(prj_in,x_in,y_in)
        x_out = pts[...,0]; y_out = pts[...,1]
    else:
        x_out, y_out = prj_out.transform_point(x_in, y_in, prj_in)
        
    return x_out,y_out