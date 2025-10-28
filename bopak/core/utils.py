import numpy as np
import xarray as xr

def valid_domain(yr_in,grid_fn):
    '''
    Parameters:
    ----------
    yr_in: int, year of domain
    grid_fn: the path to the grid file on EASE2 grid

    Returns:
    -------
    valid_range: xr.DataArray, the mask of valid range to exclude land and polehole
    '''
    with xr.open_dataset(grid_fn) as dsgrid:
            
        ocnmask = dsgrid['landmask']==0
        if yr_in<1988:
            polehole = dsgrid['p1mask']>0
        elif yr_in<2007:
            polehole = dsgrid['p2mask']>0
        else:
            polehole = dsgrid['p3mask']>0
       
        valid_range = ocnmask&(~polehole)

    return valid_range