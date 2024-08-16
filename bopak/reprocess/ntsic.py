import numpy as np
import pandas as pd
import xarray as xr
from .sic import sic

class ntsic(sic):
    
    
    
    def load_data(self,):
        '''
        '''
        with xr.open_dataset(self.fn) as ds:
            ds.load()
            vn_sics =  [vn for vn in ds if vn!='crs']
            assert len(vn_sics)==1,f'More than one/No sea ice data variable in {self.fn}'
            self.data = ds
            self.vn_sic = vn_sics[0]
            self.sic_vns = [self.vn_sic]
        return
    
    def regrid(self,mds=None,proj=None,vn_mask=None,add_mask=True,mask_value=1.01):
        '''
        '''
        if vn_mask is None: vn_mask = self.vn_sic
        super().regrid(mds=mds,sicprj=proj,vn_mask=vn_mask,add_mask=add_mask,mask_value=mask_value)
        self.data_ease2 = self.data_ease2.rename({self.vn_sic:'ntsic'})
        self.data_ease2.attrs['vn_sic'] = self.vn_sic
        return
    

