import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import xesmf
import warnings

from pathlib import Path

from .utils import set_BOEASE2, ptproj

warnings.filterwarnings("ignore")


class pond:
    def __init__(self,fn):
        '''
        '''
        assert Path(fn).exists(), f'{fn} does not exist. Please enter the full path to the file of the MODIS_melt poind fraction data.'
        self.fn = Path(fn)
        self.load_data()
        return
    
    
    def load_data(self,):
        '''
        '''
        with xr.open_dataset(self.fn) as ds:
            ds.load()
            self.data = ds
        return
    
    
    def set_regridder(
            self,mds, 
            regridder_fn=None,
            method = 'bilinear',
            reuse_weight = False,    ):
        '''
        '''
        if regridder_fn is None: 
            regridder_fn = Path('./.regridder/MPF2EASE2.nc')
            if not regridder_fn.parent.exists(): regridder_fn.parent.mkdir(parents=True,exist_ok=True)
        self.regridder_fn = regridder_fn

        dsin = xr.Dataset(coords=dict(
            lat=(['y','x'],self.data.lat.values),
            lon=(['y','x'],self.data.lon.values) ))
        
        dsout = xr.Dataset(coords=dict(
            lat=(['Y','X'],mds.lat.values),
            lon=(['Y','X'],mds.lon.values)))
        
        assert Path(regridder_fn).exists() or not reuse_weight, \
            f'Cannot reuse regridder weight without regridder file: {regridder_fn}'


        self.regridder = xesmf.Regridder(
            dsin,dsout,method,reuse_weights=reuse_weight,filename=regridder_fn,periodic=False,)

        return 
    
    
    def regrid(self):
        '''
        '''
        assert hasattr(self,'regridder'), 'Please set up regridder first.'
        
        self.data_ease2 = self.regridder(self.data,skipna=True)
         
        return 