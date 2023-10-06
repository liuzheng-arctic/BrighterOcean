import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import xesmf

import re

from pathlib import Path
from .utils import ptproj, set_BOEASE2

class sic:
    def __init__(self,fn,grid_fn = None):
        assert Path(fn).exists(), f'{fn} does not exist. Please enter the full path to the file of the PIOMAS daily ice thickness data.'
        self.fn = Path(fn)
        if grid_fn is not None:
             self.set_grid(grid_fn=grid_fn)
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
            reuse_weight = False,    ):
        '''
        '''
        assert hasattr(self,'grid'), 'Please load the grid data from the anscillary grid first'
        if regridder_fn is None: 
            regridder_fn = Path('./.regridder/CDRSIC2EASE2.nc')
            if not regridder_fn.parent.exists(): regridder_fn.parent.mkdir(parents=True,exist_ok=True)
        self.regridder_fn = regridder_fn

        dsin = xr.Dataset(coords=dict(
            lat=(['Y','X'],self.grid.latitude.values),
            lon=(['Y','X'],self.grid.longitude.values) ))
        
        dsout = xr.Dataset(coords=dict(
            lat=(['Y','X'],mds.lat.values),
            lon=(['Y','X'],mds.lon.values)))
        
        assert Path(regridder_fn).exists() or not reuse_weight, \
            f'Cannot reuse regridder weight without regridder file: {regridder_fn}'


        self.regridder = xesmf.Regridder(
            dsin,dsout,'bilinear',reuse_weights=reuse_weight,filename=regridder_fn)

        return 
    
    
    def regrid(self,mds=None,sicprj=None,vn_mask='cdr_seaice_conc',add_mask=True):
        '''
        '''
        assert hasattr(self,'regridder'), 'Please set up regridder first.'
        if mds is None: mds, _ = set_BOEASE2()
        if sicprj is None: mds_sic, sicprj = self.set_proj(self.grid)

        ## --- for testing contamination 
        # data_mask = self.data
        # if add_mask: data_mask = self.mask2nan(vn_mask=vn_mask)
        # --- add mask to filter nan and out-of-bound pixels
        #     need skipna to avoid spotty missing values due to reprojection
        #     na_thres 0 seems to work but keep at 1 to be safe
        #     Refer to Masking page on xESMF doc for more details.
        ds_mask = self._add_mask(vn_mask=vn_mask)
        ds_regrid = self.regridder(ds_mask,skipna=True,na_thres=1)
        # --- interpolated lat/lon is always trash, esp. lon near dateline
        #     use the lat/lon set by the regridder instead
        for vn in ['Latitude','Longitude']:
            if vn in ds_regrid.data_vars:
                ds_regrid = ds_regrid.drop([vn])
        self.data_ease2 = ds_regrid
        if add_mask: self._apply_mask(vn_mask=vn_mask)
        # --- the following for out-of-bound reprojection mask 
        #     seems obsolete with explicit masking. 
        # pcproj = ccrs.PlateCarree()
        # xbo,ybo = ptproj(pcproj,sicprj,mds.lon.values,mds.lat.values)
        # self.data_ease2 = ds_regrid.where(
        #         (xbo<self.grid.x.values.max())&(xbo>self.grid.x.values.min())&
        #         (ybo<self.grid.y.values.max())&(ybo>self.grid.y.values.min())
        #         )
        self.data_ease2['X'] = mds.X 
        self.data_ease2['Y'] = mds.Y

        return #self.data_ease2
    
    def _add_mask(self,vn_mask='cdr_seaice_conc'):
        '''
        '''
        # --- set land/lake/coast to nan to avoid contamination
        vns = [vn for vn in self.data.data_vars if vn!='projection']
        vmask = xr.where( self.data[vn_mask]<2, 1, 0 )
        ds_mask = self.data[vns]
        ds_mask['nMask'] = vmask
        return ds_mask
    
    def _apply_mask(self,vn_mask='cdr_seaice_conc'):
        '''
        '''
        tdat = self.data_ease2[vn_mask]
        self.data_ease2[vn_mask] = xr.where(self.data_ease2['nMask']==1, tdat, np.nan)
        return
    
    def set_grid(self,grid_fn=None,ds_grid=None):
        '''
        '''
        assert grid_fn is not None or ds_grid is not None, 'Please provide grid file name or xr.Dataset of the grid data.'
        if grid_fn is not None:
            self.grid = self.load_grid(grid_fn)
            self.grid_fn = grid_fn
        else:
            self.grid = ds_grid
        return 

    @staticmethod
    def load_grid(grid_fn):
        '''
        '''
        assert Path(grid_fn).exists(), f'{grid_fn} does not exist.'
        with xr.open_dataset(grid_fn) as ds_grid:
            ds_grid.load()
        return ds_grid
    


    @staticmethod
    def set_proj(ds_grid=None):
        '''
        '''
        # using values from CDR_SIC data and cross checked with the anscillary data for NH
        # This may be changed in future data release. NEED TO CHECK BEFORE USE!!
        if ds_grid is None:
            # hard coded for current version
            sic_globe = ccrs.Globe(semimajor_axis=6378273.0,
                                semiminor_axis=6356889.449)
            sicprj = ccrs.NorthPolarStereo(
                central_longitude=-45, true_scale_latitude=70,
                globe=sic_globe)    

            XLOW = -3837500.
            XHIGH = 3737500.
            YLOW = -5337500.
            YHIGH = 5837500.
            DX = 25e3
            DY = 25e3
        else:
            # if the next version uses the anscillary grid file of this format
            pcproj = ccrs.PlateCarree()
            sic_globe = ccrs.Globe(semimajor_axis=ds_grid.crs.attrs['semimajor_axis'],
                                semiminor_axis=ds_grid.crs.attrs['semiminor_axis'])
            sicprj = ccrs.NorthPolarStereo(
                central_longitude=ds_grid.crs.attrs['longitude_of_projection_origin'],
                true_scale_latitude=ds_grid.crs.attrs['standard_parallel'],
                globe=sic_globe)
            XLOW = ds_grid.x.values.min()
            YLOW = ds_grid.y.values.min()
            XHIGH = ds_grid.x.values.max()
            YHIGH = ds_grid.y.values.max()
            DX = abs(np.diff(ds_grid.x.values).mean())
            DY = abs(np.diff(ds_grid.y.values).mean())

        xc = np.arange(XLOW,XHIGH+DX,DX)
        yc = np.arange(YLOW,YHIGH+DY,DY)
        xcs,ycs = np.meshgrid(xc,yc)
        lonc, latc = ptproj(sicprj,pcproj, xcs,ycs)
     
        mds = xr.Dataset()
        mds['X'] = xc
        mds['Y'] = yc
        mds['lat'] = (['Y','X'],latc)
        mds['lon'] = (['Y','X'],lonc)        

        mds.attrs['XLOW'  ] = XLOW
        mds.attrs['YLOW'  ] = YLOW
        mds.attrs['XHIGH'  ] = XHIGH
        mds.attrs['YHIGH'  ] = YHIGH
        mds.attrs['dx'    ] = DX
        mds.attrs['dy'    ] = DY
        mds.attrs['lat_true'] = ds_grid.crs.attrs['standard_parallel']
        mds.attrs['lon0'  ] = ds_grid.crs.attrs['longitude_of_projection_origin']
        mds.attrs['grid'  ] = ds_grid.crs.attrs['grid_mapping_name']
      
        return mds, sicprj
    

      
    @property
    def time(self):
        return self._time()
    
    
    
    def _time(self):
        rx = re.compile('(\d{8})')
        return pd.to_datetime(rx.findall(self.fn.name).pop())