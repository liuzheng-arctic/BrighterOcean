import numpy as np
import pandas as pd
import xarray as xr
import gzip
import struct
import re
import xesmf
from pathlib import Path


PIOMAS_DIMS = (120,360)


class piomas:
    def __init__(self,fn,grid_fn = None):
        '''
        Parameters:
        ----------
        fn: str, full path to the PIOMAS gzipped data
        grid_fn: str, full path to the PIOMAS grid data
        '''
        assert Path(fn).exists(), f'{fn} does not exist. Please enter the full path to the file of the PIOMAS daily ice thickness data.'
        self.fn = Path(fn)
        self.dims = PIOMAS_DIMS
        if grid_fn is not None:
            self.set_grid(grid_fn=grid_fn)
        # if grid_fn is not None:
        #     self.grid_fn = grid_fn
        #     lon, lat = self.read_grid(grid_fn)
        #     self.add_grid(lon,lat)
        return
    
    def load_data(self,grid_fn=None,lat=None,lon=None):
        '''
        '''
        assert grid_fn is not None or (lat is not None and lon is not None), 'Please provide grid filename or the lat/lon arrays.'
        
        if grid_fn is not None:
            self.set_grid(grid_fn=grid_fn)
        else:
            self.set_grid(lon=lon,lat=lat)

        data3d = self.read_hiday(self.fn)
        self.data['ice thickness'] = (['time','Y','X'],data3d)
        return
    
    
    def set_grid(self,lon=None, lat=None,grid_fn=None):
        '''
        add lon, lat of grid to instance

        Parameters:
        ----------
        lon: np.ndarray
        '''
        if grid_fn is not None:
            self.grid_fn = grid_fn
            lon,lat = self.read_grid(grid_fn)
        self.lon = lon
        self.lat = lat
        self.time = self.daily_time()
        ds = xr.Dataset()
        ds['time'] = self.time
        ds['lon'] = (['Y','X'],lon)
        ds['lat'] = (['Y','X'],lat)
        self.data = ds
        return 
    
    def set_regridder(
            self,mds, 
            regridder_fn=None,
            reuse_weight = False,       
            method='bilinear',     
            ):
        '''
        '''
        assert hasattr(self,'data'), 'Please setup grid for PIOMAS data first.'

        if regridder_fn is None: 
            regridder_fn = Path('./.regridder/PIOMAS2EASE2.nc')
            if not regridder_fn.parent.exists(): regridder_fn.parent.mkdir(parents=True,exist_ok=True)
        self.regridder_fn = regridder_fn

        dsin = xr.Dataset(coords=dict(lat=(['Y','X'],self.data.lat.values),lon=(['Y','X'],self.data.lon.values)))
        dsout = xr.Dataset(coords=dict(lat=(['Y','X'],mds.lat.values),
                                    lon=(['Y','X'],mds.lon.values)))
        self.regridder = xesmf.Regridder(dsin,dsout,method=method,periodic=True,
                                reuse_weights=reuse_weight,filename=self.regridder_fn)
        return 
    
    
    def regrid(self,):
        '''
        '''
        assert hasattr(self,'regridder'), 'Please setup regridder first.'
        self.data_ease2 = self.regridder(self.data)
        return
    
    def save_regrid(self,outfn):
        '''
        '''
        assert hasattr(self,'data_ease2'), 'Please regrid data first.'
        fcomp = dict(zlib=True, complevel=5, dtype='float32')
        encoding = {var: fcomp for var in self.data_ease2.data_vars}
        self.data_ease2.to_netcdf(outfn,encoding=encoding)
        return
    

    @staticmethod
    def read_hiday(fn):
        '''
        Read gzipped PIOMAS daily sea ice thickness data.

        Parameters:
        ----------
        fn: str, full path to the PIOMAS gzipped data

        Return:
        ------
        data3d: np.ndarray, sea ice thickness data with dimension of [day,Y,X]

        '''
        
        with gzip.open(fn, mode='rb') as file:

            fileContent = file.read()
            data = struct.unpack("f" * (len(fileContent)// 4), fileContent)
            ny,nx = PIOMAS_DIMS
            nday = int(len(data)/ny/nx)
            data3d = np.array(data).reshape((nday,ny,nx))
        return data3d
    

    
    @staticmethod
    def read_grid(grid_fn):
        '''
        Read PIOMAS grid data from file

        '''
        grid_data = np.loadtxt(grid_fn)
        lon, lat = grid_data.ravel().reshape((2,)+PIOMAS_DIMS)
        return lon, lat
    
    @property
    def year(self):
        return self._year()
    
    @property
    def year(self):
        return self._year()
    
    def daily_time(self):  
        '''
        For leap years, the last day is dropped in PIOMAS for similicity. 
        '''
        isleap = pd.to_datetime(f'{self.year}0101',format='%Y%m%d').is_leap_year
        end_date = '1230' if isleap else '1231'
        return pd.date_range(self.year+'0101',self.year+end_date,freq='D')
    
    def _year(self):
        rx = re.compile('H(\d{4})')
        return rx.findall(self.fn.name).pop()
    #def regrid(self,)
