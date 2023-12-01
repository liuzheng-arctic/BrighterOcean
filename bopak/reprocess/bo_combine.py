import numpy as np
import pandas as pd
import xarray as xr
import json

import timeit
import warnings
from pathlib import Path

# from bopak.reprocess import albedo, onset, sic, iceage
# from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from .utils import set_BOEASE2
from .database import combined_db
from .albedo import OCN_ALB

KAPPA0 = 1 # [1/m]
KA_ICE = 1
KA_POND = 0.75
MPF_DAYS = 128

class combined_data:
    def __init__(
        self,year,mds=None,boprj=None,
        albdir=None, onsetdir=None,iadir=None,sicdir=None,
        eradir=None,piodir=None,ponddir=None,db_in=None,
        ):
        '''
        '''

        self.year = year
        self.set_grid(mds=mds,boprj=boprj)
        self.vns = ['albedo','sic','onset','iceage','piomas','era5']
        if db_in is None:
            db_in = combined_db(
                albdir=albdir,sicdir=sicdir,onsetdir=onsetdir,
                iadir=iadir, eradir=eradir,piodir=piodir,ponddir=ponddir)
            self.db = db_in.query(year)
        else:
            assert isinstance(db_in, combined_db), f'db_in must be an instance of class bodb defined in bopak.reprocess.database .'
            self.db = db_in.query(year)

        return
    
    def process(self,):
        '''
        '''     
        for vn in self.vns:
            assert vn in self.db.keys(), f'Path or database for dataset {vn} is not provided.'
        self.load_data()

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

        for vn in ['h_nIC_MYI','h_nIC_FYI','h_nIC_AGE','h_ocn_SIC',
                   'h_ice_MYI','h_ice_FYI','h_ice_AGE']:
            self.data[f'accu_{vn}'] = self.cumsum(vn)*86400./1e6 # W/m2 to J/m2 to MJ/m2
        self.data['accu_h_all_MYI'] = self.data['accu_h_ice_MYI'] + self.data['accu_h_ocn_SIC'] 
        self.data['accu_h_all_FYI'] = self.data['accu_h_ice_FYI'] + self.data['accu_h_ocn_SIC'] 
        self.data['accu_h_all_AGE'] = self.data['accu_h_ice_AGE'] + self.data['accu_h_ocn_SIC'] 

        ext0 = np.exp(-KAPPA0*self.data['sith'])
        self.data['hthr_ice_MYI_KAPPA0'] = self.data['h_ice_MYI'] * ext0
        self.data['hthr_ice_FYI_KAPPA0'] = self.data['h_ice_FYI'] * ext0
        self.data['hthr_ice_AGE_KAPPA0'] = self.data['h_ice_AGE'] * ext0
        for vn in ['hthr_ice_MYI_KAPPA0','hthr_ice_FYI_KAPPA0','hthr_ice_AGE_KAPPA0']:
            self.data[f'accu_{vn}'] = self.cumsum(vn)*86400./1e6 

        # ----- For ponded ice ----
        # Kappa = (1-MeltPondFraction)*KA_ICE + MeltPondFraction*KA_POND
        #       = KA_ICE + MeltPondFraction*(KA_POND-KA_ICE)
        if 'pond' in self.db.keys():
            ext_pond = np.exp(-self.data['kappa']*self.data['sith'])
            self.data['hthr_ice_MYI_pond'] = self.data['h_ice_MYI'] * ext_pond
            self.data['hthr_ice_FYI_pond'] = self.data['h_ice_FYI'] * ext_pond
            self.data['hthr_ice_AGE_pond'] = self.data['h_ice_AGE'] * ext_pond
            for vn in ['hthr_ice_MYI_pond','hthr_ice_FYI_pond','hthr_ice_AGE_pond']:
                self.data[f'accu_{vn}'] = self.cumsum(vn)*86400./1e6         
        
            self.data['kappa'] = self.data['kappa'].where(self.data['land']==0,-.999)
            self.data['kappa_org'] = self.data['kappa_org'].where(self.data['land']==0,-.999)

        self._add_column_number()
        self._add_fill_values()
        return #mpf_mask

    def to_file(self,outroot,ftype='txt',freq='D'):
        '''
        '''

        assert ftype in ['txt','nc','csv', 'netcdf'], 'ftype has to be txt, csv, or nc.'

        txt_year = self.time[0].strftime(format='%Y')
        outdir = Path(outroot)/ftype
        if freq=='D': outdir = outdir/txt_year
        if not outdir.exists(): outdir.mkdir(parents=True,exist_ok=True)
        
        with open(outroot/f'BO_var_id_{txt_year}.json','w') as fid:
            json.dump(self.vdict,fid)

        if ftype in ['txt','csv']:
            assert freq=='D', 'Output in txt or csv has to be daily.'
            odf = pd.DataFrame()
            vdict = self.vdict
            vn2ds = ['lat','lon','emelt','fmelt','efreeze','ffreeze']
            ncols = len(vdict.keys())
            nt = len(self.time)
            for it in range(nt):
                txt_date = self.time[it].strftime(format='%Y-%m-%d')
                ofn = outdir/f'BO_EASE2_gridded_25km_{txt_date}.{ftype}'
                if self.time[it].day==1: print(it,ofn)
                for icol in range(ncols):
                    vn = vdict[icol+2]
                    if vn in vn2ds:
                        dat = self.data[vn].values.ravel()
                    else:
                        dat = self.data[vn].values[it].ravel()
                    odf[vn] = dat
                odf.index += 1
                odat = np.column_stack([odf.index,odf.values])
                np.savetxt(ofn,odat,fmt='%6d\t'+'%10.4f '*ncols,delimiter='\t')
        else:
            fcomp = dict(zlib=True, complevel=5, dtype='float32')
            encoding = {var: fcomp for var in self.data.data_vars}
            if freq=='D':
                nt = len(self.time)
                for it in range(nt):
                    ods = self.data.isel(time=it)
                    txt_date = self.time[it].strftime(format='%Y-%m-%d')
                    ofn = outdir/f'BO_EASE2_gridded_25km_{txt_date}.{ftype}'
                    if self.time[it].day==1: print(it,ofn)
                    ods.to_netcdf(ofn,encoding=encoding)
            else: 
                ods = self.data
                ofn = outdir/f'BO_EASE2_gridded_25km_{txt_year}.{ftype}'
                print(txt_year,ofn)
                self.data.to_netcdf(ofn,encoding=encoding)

        return

    def _daily_arr(self,it):
        '''
        '''
        odf = pd.DataFrame()
        vdict = {}
        # --- sort variables ----
        for vn in self.data.data_vars:
            if 'ColNo' in self.data[vn].attrs.keys():
                colno = self.data[vn].attrs['ColNo']
                vdict.update({colno:vn})
        vn2ds = ['lat','lon','emelt','fmelt','efreeze','ffreeze']
        ncols = len(vdict.keys())
        for icol in range(ncols):
            vn = vdict[icol+2]
            if vn in vn2ds:
                dat = tbda.data[vn].values.ravel()
            else:
                dat = tbda.data[vn].values[it].ravel()
            odf[vn] = dat
        odf.index += 1
        odat = np.column_stack([odf.index,odf.values])
        return odat
    
    def cumsum(self,vn):
        '''
        The direct cumsum might be problematic. So leave this interface for future changes.
        '''
        dat = self.data[vn]
        #dmask = (~dat.isnull())
        dat_sum = dat.cumsum(dim='time')
        #dmask_sum = dmask.cumsum(dim='time')
        #nday = xr.DataArray(self.time.dayofyear.values, coords=self.data.time.coords)
        return dat_sum
    
    def load_data(self):
        '''
        '''
        nt = len(self.time)
        ny,nx = self.data.lat.shape

        prd = 'era5'; vn = 'ssrd'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['insol'] = (['time','Y','X'], ds[vn].values/86400.)
            self.data['insol'].attrs['unit'] = 'W/m2'

        prd = 'sic'; vn = 'cdr_seaice_conc'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['sic'] = (['time','Y','X'], ds[vn].values)

        prd = 'onset'
        with xr.open_dataset(self.db[prd]) as ds:   
                self.data['emelt'] = (['Y','X'], ds['Earlymelt'].values)
                self.data['fmelt'] = (['Y','X'], ds['Melt'].values)
                self.data['efreeze'] = (['Y','X'], ds['Earlyfreeze'].values)
                self.data['ffreeze'] = (['Y','X'], ds['Freeze'].values)

        prd = 'iceage'; vn = 'age_of_sea_ice'
        self.data['iceage'] = (['time','Y','X'], np.zeros((nt,ny,nx)))
        with xr.open_dataset(self.db[prd]) as ds:
            for i_week in range(52):
                i0 = i_week*7
                i1 = i0 + 7
                if i_week==51: i1 = nt
                self.data['iceage'][i0:i1] =  ds[vn].values[i_week]

        prd = 'albedo'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['alb_FYI'] = (['time','Y','X'], ds['first_year'].values)
            self.data['alb_MYI'] = (['time','Y','X'], ds['multi_year'].values)
            self.data['alb_AGE'] = (['time','Y','X'], ds['actual_age'].values)

        prd = 'piomas'; vn = 'ice thickness'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['sith'] = (['time','Y','X'], np.zeros((nt,ny,nx)) )
            if self.time[0].is_leap_year:
                self.data['sith'][:-1] = ds[vn].values
                self.data['sith'][-1] = ds[vn].values[-1]
            else:
                self.data['sith'][:] = ds[vn].values

        prd = 'pond'; vn ='mpf'
        if prd in self.db.keys():
            with xr.open_dataset(self.db[prd]) as ds:
                self.data['mpf'] = (['time','Y','X'], np.zeros((nt,ny,nx)) )
                if 'land' not in self.data.data_vars:
                    self.data['land'] = (['Y','X'],ds['land'][0].values) 
                doy = ds.time.dt.dayofyear.values
                for it, tdoy in enumerate(doy):
                    i0 = doy[it]-1
                    i1 = doy[it]+7
                    self.data['mpf'][i0:i1] = ds[vn].values[it]/100

            self.data['kappa_org'] = KA_ICE + self.data['mpf']*(KA_POND-KA_ICE)
            self.data['mpf_fill'] = self._fill_kappa('mpf')
            self.data['kappa'] = self._fill_kappa('kappa_org')
                    
        return #ds


    @property
    def vdict(self):
        return self._vdict()

    def _vdict(self):
        '''
        '''
        vdict = {}
        # --- sort variables ----
        for vn in self.data.data_vars:
            if 'ColNo' in self.data[vn].attrs.keys():
                colno = self.data[vn].attrs['ColNo']
                vdict.update({colno:vn})
        return vdict
    
    def _add_fill_values(self):
        '''
        '''
        for vn in self.data.data_vars:
            if vn in ['emelt','fmelt','efreeze','ffreeze','iceage']:
                #self.data[vn] = self.data[vn].where(~self.data[vn].isnull(),-99)
                self.data[vn] = self.data[vn].fillna(-99)
                self.data[vn].attrs['_FillValue'] = -99
            elif 'alb' in vn or vn in ['mpf','mpf_fill','kappa','kappa_fill','sic']:
                #self.data[vn] = self.data[vn].where(self.data[vn].isnull(),-.999)
                self.data[vn] = self.data[vn].fillna(-.999)
                self.data[vn].attrs['_FillValue'] = -.999
            else:
                #self.data[vn] = self.data[vn].where(self.data[vn].isnull(),-999.99)
                self.data[vn] = self.data[vn].fillna(-999.99)
                self.data[vn].attrs['_FillValue'] = -999.99
        
        return
    
    def _add_column_number(self):
        '''
        '''
        self.data['lat'].attrs['ColNo'] = 2
        self.data['lon'].attrs['ColNo'] = 3
        self.data['insol'].attrs['ColNo'] = 4
        self.data['sic'].attrs['ColNo'] = 5
        self.data['emelt'].attrs['ColNo'] = 6
        self.data['fmelt'].attrs['ColNo'] = 7
        self.data['efreeze'].attrs['ColNo'] = 8
        self.data['ffreeze'].attrs['ColNo'] = 9
        self.data['iceage'].attrs['ColNo'] = 10
        self.data['alb_MYI'].attrs['ColNo'] = 11
        self.data['alb_FYI'].attrs['ColNo'] = 12
        self.data['alb_AGE'].attrs['ColNo'] = 13
        self.data['h_nIC_MYI'].attrs['ColNo'] = 14
        self.data['h_ice_MYI'].attrs['ColNo'] = 15
        self.data['h_ocn_SIC'].attrs['ColNo'] = 16
        self.data['h_all_MYI'].attrs['ColNo'] = 17
        self.data['accu_h_nIC_MYI'].attrs['ColNo'] = 18
        self.data['accu_h_ice_MYI'].attrs['ColNo'] = 19
        self.data['accu_h_ocn_SIC'].attrs['ColNo'] = 20
        self.data['accu_h_all_MYI'].attrs['ColNo'] = 21
        self.data['h_nIC_FYI'].attrs['ColNo'] = 22
        self.data['h_ice_FYI'].attrs['ColNo'] = 23
        self.data['h_all_FYI'].attrs['ColNo'] = 24
        self.data['accu_h_nIC_FYI'].attrs['ColNo'] = 25
        self.data['accu_h_ice_FYI'].attrs['ColNo'] = 26
        self.data['accu_h_all_FYI'].attrs['ColNo'] = 27
        self.data['h_nIC_AGE'].attrs['ColNo'] = 28
        self.data['h_ice_AGE'].attrs['ColNo'] = 29
        self.data['h_all_AGE'].attrs['ColNo'] = 30
        self.data['accu_h_nIC_AGE'].attrs['ColNo'] = 31
        self.data['accu_h_ice_AGE'].attrs['ColNo'] = 32
        self.data['accu_h_all_AGE'].attrs['ColNo'] = 33
        self.data['sith'].attrs['ColNo'] = 34
        self.data['hthr_ice_MYI_KAPPA0'].attrs['ColNo'] = 35
        self.data['hthr_ice_FYI_KAPPA0'].attrs['ColNo'] = 36
        self.data['hthr_ice_AGE_KAPPA0'].attrs['ColNo'] = 37
        self.data['accu_hthr_ice_MYI_KAPPA0'].attrs['ColNo'] = 38
        self.data['accu_hthr_ice_FYI_KAPPA0'].attrs['ColNo'] = 39
        self.data['accu_hthr_ice_AGE_KAPPA0'].attrs['ColNo'] = 40
        if 'pond' in self.db.keys():
            self.data['mpf'].attrs['ColNo'] = 41
            self.data['kappa'].attrs['ColNo'] = 42
            self.data['hthr_ice_MYI_pond'].attrs['ColNo'] = 43
            self.data['hthr_ice_FYI_pond'].attrs['ColNo'] = 44
            self.data['hthr_ice_AGE_pond'].attrs['ColNo'] = 45
            self.data['accu_hthr_ice_MYI_pond'].attrs['ColNo'] = 46
            self.data['accu_hthr_ice_FYI_pond'].attrs['ColNo'] = 47
            self.data['accu_hthr_ice_AGE_pond'].attrs['ColNo'] = 48


        return
    
    def _fill_kappa(self,vn):
        '''
        '''
        dat_fill3d = xr.apply_ufunc(
            self._fill_kappa_1d,
            self.data[vn],
            input_core_dims= [['time']], 
            output_core_dims= [['time']], 
            exclude_dims=set(('time',)),
            vectorize=True)
        return dat_fill3d.transpose('time',...)

    def _fill_kappa_1d(self,x):
        '''
        Use previous valid value to fill.
        This works because before day 129, the fill values are 0 for mpf and 1 for kappa. 
        '''
        y = x*1.
        idx = np.where(np.isnan(y))[0]
        for i in idx:
            y[i] = y[i-1]
        return y

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