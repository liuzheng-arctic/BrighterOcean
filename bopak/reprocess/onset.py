import numpy as np
import pandas as pd
import xarray as xr
from .sic import sic

class onset(sic):
    
    def _add_mask(self,vn_mask=None):
        '''
        '''
        # --- set land/lake/coast to nan to avoid contamination
        ds_mask = self.data
        for vn in ['Melt','Earlymelt','Freeze','Earlyfreeze']:
            ds_mask[f'nMask_{vn}'] = xr.where( self.data[vn]<0, 0, 1 )
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
        '''
        super().regrid(mds=mds,sicprj=proj,vn_mask=vn_mask,add_mask=add_mask)
        return #self.data_ease2

    def _time(self):
        return pd.to_datetime(self.fn.name[:4])