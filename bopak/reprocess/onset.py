import numpy as np
import pandas as pd
import xarray as xr
from .sic import sic

from .utils import ptproj, set_BOEASE2

class onset(sic):


    def load_data(self,ds_grid=None):
        '''
        '''
        with xr.open_dataset(self.fn) as ds:
            ds.load()
            if ds_grid is None and hasattr(self,'grid'): ds_grid = self.grid
            if ds_grid is not None:
                ds['x'] = ds_grid['x'].values
                ds['y'] = ds_grid['y'].values
            #ds = ds.sortby('y')
            self.data = ds#.sortby('y')
        self.sic_vns = ['Melt','Earlymelt','Freeze','Earlyfreeze']
        self.source = 'onset'
        return    
    def _add_mask(self,vn_mask=None,mask_value=2):
        '''
        '''
        # --- set land/lake/coast to nan to avoid contamination
        ds_mask = self.data
        for vn in ['Melt','Earlymelt','Freeze','Earlyfreeze']:
            ds_mask[f'nMask_{vn}'] = xr.where( self.data[vn]>0, 1, 0 )
        return ds_mask
    
    def _apply_mask(self,vn_mask='cdr_seaice_conc'):
        '''
        '''
        for vn in ['Melt','Earlymelt','Freeze','Earlyfreeze']:
            tdat = self.data_ease2[vn]
            self.data_ease2[vn] = xr.where(self.data_ease2[f'nMask_{vn}']==1, tdat, np.nan)
        return
    
    def regrid(self,mds=None,proj=None,vn_mask=None,add_mask=True):
        '''
        Due to the complicated changes to regrid in sic to accommodate bt/nt data, 
        this function has to be a local version different from sic. 
        '''
        #super().regrid(mds=mds,sicprj=proj,vn_mask=vn_mask,add_mask=add_mask,mask_value=mask_value)
        
        assert hasattr(self,'regridder'), 'Please set up regridder first.'
        if mds is None: mds, _ = set_BOEASE2()
        if proj is None: mds_sic, sicprj = self.set_proj(self.grid)

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
        for vn in ['Melt','Earlymelt','Freeze','Earlyfreeze']:
            self.data_ease2[vn] = self.data_ease2[vn].round()
        return #* self.data_ease2

    def _time(self):
        return pd.to_datetime(self.fn.name[:4])