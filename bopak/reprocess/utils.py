import numpy as np
import xarray as xr
import pandas as pd
import gzip
import struct
import cartopy.crs as ccrs

def set_ease(
        DX = 25e3,
        DY = 25e3,
        XLIM = 5000e3,
        YLIM = 5000e3,
):
    '''
    '''
    lon0 = 0
    lat0 = 90
    pcproj = ccrs.PlateCarree()
    ease2prj = ccrs.LambertAzimuthalEqualArea(central_longitude=lon0,central_latitude=lat0)

    nx   = int(XLIM*2/DX)
    ny   = int(XLIM*2/DY)
    print(ny,nx)
    xgrid = np.linspace(-XLIM,XLIM,nx+1)
    ygrid = np.linspace(-YLIM,YLIM,ny+1)
        
    xc = (xgrid[1:]+xgrid[:-1])/2
    yc = (ygrid[1:]+ygrid[:-1])/2
    xcs,ycs = np.meshgrid(xc,yc)
    lonc, latc = ptproj(ease2prj,pcproj, xcs,ycs)
     
    mds = xr.Dataset()
    mds['X'] = xc
    mds['Y'] = yc
    mds['lat'] = (['Y','X'],latc)
    mds['lon'] = (['Y','X'],lonc)
    

    mds.attrs['XLIM'  ] = XLIM
    mds.attrs['YLIM'  ] = YLIM
    mds.attrs['dx'    ] = DX
    mds.attrs['dy'    ] = DY
    mds.attrs['grid'  ] = 'EASE2.0'
    mds.attrs['lat0'  ] = lat0
    mds.attrs['lon0'  ] = lon0
    return mds, ease2prj


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