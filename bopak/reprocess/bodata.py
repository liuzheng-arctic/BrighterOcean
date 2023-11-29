import numpy as np
import pandas as pd
import xarray as xr

import timeit
import warnings
from pathlib import Path

# from bopak.reprocess import albedo, onset, sic, iceage
# from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from .utils import set_BOEASE2
from .database import bodb
from .albedo import OCN_ALB

KAPPA0 = 1 # [1/m]

class bodata_annual:
    def __init__(
        self,year,mds=None,boprj=None,
        albdir=None, onsetdir=None,iadir=None,sicdir=None,
        eradir=None,piodir=None,combine_db=None,
        ):
        '''
        '''

        self.year = year
        self.set_grid(mds=mds,boprj=boprj)
        self.vns = ['albedo','sic','onset','iceage','piomas','era5']
        if combine_db is None:
            combine_db = bodb(
                albdir=albdir,sicdir=sicdir,onsetdir=onsetdir,
                iadir=iadir, eradir=eradir,piodir=piodir)
            self.db = combine_db.query(year)
        else:
            assert isinstance(combine_db, bodb), f'combine_db must be an instance of class bodb defined in bopak.reprocess.database .'
            self.db = combine_db.query(year)

        return
    
    def process(self,):
        '''
        '''     
        for vn in self.vns:
            assert vn in self.db.keys(), f'Path or database for dataset {vn} is not provided.'
        self._load_data()

        self.data['h_nIC_MYI'] = self.data['insol']*(1-self.data['alb_MYI'])
        self.data['h_ice_MYI'] = self.data['h_nIC_MYI']*self.data['sic']
        self.data['h_ocn_SIC'] = self.data['insol']*(1-OCN_ALB)*(1-self.data['sic'])
        self.data['h_all_MYI'] = self.data['h_ice_MYI'] + self.data['h_ocn_SIC']

        self.data['h_nIC_FYI'] = self.data['insol']*(1-self.data['alb_FYI'])
        self.data['h_ice_FYI'] = self.data['h_nIC_FYI']*self.data['sic']
        self.data['h_all_FYI'] = self.data['h_ice_FYI'] + self.data['h_ocn_SIC']
        
        self.data['h_nIC_AGE'] = self.data['insol']*(1-self.data['alb_AGE'])
        self.data['h_ice_AGE'] = self.data['h_nIC_AGE']*self.data['sic']
        self.data['h_all_AGE'] = self.data['h_ice_AGE'] + self.data['h_ocn_SIC']

        ext0 = np.exp(-KAPPA0*self.data['sith'])
        self.data['hthr_ice_MYI'] = self.data['h_ice_MYI'] * ext0
        self.data['hthr_ice_FYI'] = self.data['h_ice_FYI'] * ext0

        
        return
    
    def _load_data(self):
        '''
        '''
        nt = len(self.time)
        ny,nx = self.data.lat.shape

        prd = 'sic'; vn = 'cdr_seaice_conc'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['sic'] = (['time','Y','X'], ds[vn].values)

        prd = 'iceage'; vn = 'age_of_sea_ice'
        self.data['iceage'] = (['time','Y','X'], np.zeros((nt,ny,nx)))
        with xr.open_dataset(self.db[prd]) as ds:
            for i_week in range(52):
                i0 = i_week*7
                i1 = i0 + 7
                if i_week==51: i1 = nt
                self.data['iceage'][i0:i1] =  ds[vn].values[i_week]

        prd = 'era5'; vn = 'ssrd'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['insol'] = (['time','Y','X'], ds[vn].values/86400.)

        prd = 'onset'
        with xr.open_dataset(self.db[prd]) as ds:   
                self.data['emelt'] = (['Y','X'], ds['Earlymelt'].values)
                self.data['fmelt'] = (['Y','X'], ds['Melt'].values)
                self.data['efreeze'] = (['Y','X'], ds['Earlyfreeze'].values)
                self.data['ffreeze'] = (['Y','X'], ds['Freeze'].values)

        prd = 'albedo'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['alb_FYI'] = (['time','Y','X'], ds['first_year'].values)
            self.data['alb_MYI'] = (['time','Y','X'], ds['multi_year'].values)
            self.data['alb_AGE'] = (['time','Y','X'], ds['actual_age'].values)

        prd = 'piomas'; vn = 'ice thickness'
        with xr.open_dataset(self.db[prd]) as ds:
            print(prd,vn,ds[vn].shape)
            self.data['sith'] = (['time','Y','X'], np.zeros((nt,ny,nx)) )
            self.data['sith'][:-1] = ds[vn].values
            self.data['sith'][-1] = ds[vn].values[-1]
        return

    def set_grid(self,mds=None,boprj=None):
        '''
        '''
        if mds is None or boprj is None:
            mds, boprj = set_BOEASE2()
        self.mds = mds
        self.boprj = boprj

        ds = xr.Dataset()
        ds['X'] = mds['X']
        ds['Y'] = mds['Y']
        ds['lat'] = (['Y','X'],mds['lat'].values)
        ds['lon'] = (['Y','X'],mds['lon'].values)
        self.time = self.daily_time()
        ds['time'] = self.time
        self.data = ds
        return


    def daily_time(self):  
        yrstr = str(self.year)
        return pd.date_range(yrstr+'0101',yrstr+'1231',freq='D')