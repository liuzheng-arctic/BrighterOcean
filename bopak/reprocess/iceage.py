import numpy as np
import pandas as pd
import xarray as xr
import re
import xesmf
import cartopy.crs as ccrs
from pathlib import Path
import cftime

from .utils import set_ease, set_BOEASE2, ptproj

class iceage:
    def __init__(self,fn):
        '''
        '''
        assert Path(fn).exists(), f'{fn} does not exist. Please enter the full path to the file of the PIOMAS daily ice thickness data.'
        self.fn = Path(fn)
        # if grid_fn is not None:
        #      self.set_grid(grid_fn=grid_fn)
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
            regridder_fn = Path('./.regridder/ICEAGE2EASE2.nc')
            if not regridder_fn.parent.exists(): regridder_fn.parent.mkdir(parents=True,exist_ok=True)
        self.regridder_fn = regridder_fn

        dsin = xr.Dataset(coords=dict(
            lat=(['y','x'],self.data.latitude.values),
            lon=(['y','x'],self.data.longitude.values) ))
        
        dsout = xr.Dataset(coords=dict(
            lat=(['Y','X'],mds.lat.values),
            lon=(['Y','X'],mds.lon.values)))
        
        assert Path(regridder_fn).exists() or not reuse_weight, \
            f'Cannot reuse regridder weight without regridder file: {regridder_fn}'


        self.regridder = xesmf.Regridder(
            dsin,dsout,method,reuse_weights=reuse_weight,filename=regridder_fn)

        return 
    
    
    
    def regrid(self,mds=None,vn_mask='age_of_sea_ice'):
        '''
        '''
        assert hasattr(self,'regridder'), 'Please set up regridder first.'
        if mds is None: mds, _ = set_BOEASE2()
        #if sicprj is None: mds_sic, sicprj = self.set_proj()

        
        ds_mask = self._add_mask(vn_mask=vn_mask)
        #tdata = self.data.where(self.data['age_of_sea_ice']<20)
        ds_regrid = self.regridder(ds_mask,skipna=True)
        ds_regrid['Y'] = mds['Y']
        ds_regrid['X'] = mds['X']
        #ds_grid['age_of_sea_ice'].round()
        self.data_ease2 = ds_regrid
        self._apply_mask(vn_mask=vn_mask)
        self.data_ease2 = np.ceil(self.data_ease2)
         
        return #ds_regrid
    

    def _add_mask(self,vn_mask='age_of_sea_ice'):
        '''
        '''
        # --- set land/lake/coast to nan to avoid contamination
        vns = [vn_mask]
        it = 0
        vmask = xr.where( (self.data[vn_mask][it]>=0)&(self.data[vn_mask][it]<20), 1, 0 )
        ds_mask = self.data[vns]
        ds_mask['nMask'] = vmask
        return ds_mask
    
    def _apply_mask(self,vn_mask='age_of_sea_ice'):
        '''
        '''
        tdat = self.data_ease2[vn_mask]
        #self.data_ease2[vn_mask] = xr.where(self.data_ease2['nMask']==1, tdat, np.nan).transpose('time','Y','X')
        self.data_ease2[vn_mask] = tdat.where(self.data_ease2['nMask']==1)
        return
    
    def interp(self,t_stamp,
               grid='org',
               format='%Y-%m-%d',
               vn = 'age_of_sea_ice',
               method = 'linear',
               round = True,
               offset = 3,
               ):
        '''
        Three ways:
        (1) method='linear' and round=True to interpolate over time
        (2) choose closest time stamp, set round=False to avoid extra work, use offset=0
        (3) Similar to (2) but use offset=3. This may make more sense because 
            the time stamps in the IA dataset are not centered in the week.  
        Both (2) and (3) use a signle value for a week. 

        '''
        assert grid in ['org','regrid'], 'Unknow grid type.'
        assert type(t_stamp) in [str,pd.Timestamp]
        if grid == 'org':
            ds = self.data[vn]
        elif grid == 'regrid':
            ds = self.data_ease2[vn]
        tpd = pd.to_datetime(t_stamp,format=format)
        t_ctr = (tpd+pd.Timedelta(offset,'D')).strftime(format=format)
        t_intp = cftime.DatetimeJulian.strptime(t_ctr,format,calendar='julian')
        ods = ds.interp(time=t_intp,method=method)
        if round: ods = ods.round()
        return ods

    
    @staticmethod
    def set_proj():
        '''
        '''
        EASE_VER = 1
        RES = 12.5
        return  set_ease(version=EASE_VER,res=RES)
    
    
      
    @property
    def year(self):
        return self._year()
    
    def daily_time(self):  
        return pd.date_range(self.year+'0101',self.year+'1231',freq='D')
    
    def _year(self):
        rx = re.compile('(\d{8})')
        return rx.findall(self.fn.name).pop()[:4]