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
from .albedo import OCN_ALB, POND_ALB

KA_ICE = 1
KA_POND = 0.75
MPF_DAYS = 128
NMISS0_SIC = 39135

Y_0 = 0.58
KAPPA0 = 1 # [1/m]
KAPPA1 = 1
KAPPA2 = 2

KAPPA_DEFT = 1
KAPPA_SNOW = 3 # Use 3 to allow for biology under snow covered ice.
KAPPA_INF  = 10000
KAPPA_POND = 0.7

class combined_data:
    def __init__(
        self,year,mds=None,boprj=None,
        albdir=None, onsetdir=None,iadir=None,sicdir=None,
        sic_intp_fn=None, #sic_miss_fn=None,
        eradir=None,piodir=None,ponddir=None,db_in=None,
        ):
        '''
        '''

        self.year = year
        self.sic_intp_fn = sic_intp_fn
        #self.sic_miss_fn = sic_miss_fn
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
    
    def process(self,interp_sic=True,MYI_out=True,FYI_out=True,full_out=True,skipna=True):
        '''
        '''     
        for vn in self.vns:
            print(vn)
            assert vn in self.db.keys(), f'Path or database for dataset {vn} is not provided.'
        self.load_data_all(interp_sic=interp_sic)

        self.data.attrs['OCN_ALB'] = OCN_ALB
        self.data.attrs['KAPPA_DEFT'] = KAPPA_DEFT
        self.data.attrs['KAPPA_SNOW'] = KAPPA_SNOW
        self.data.attrs['KAPPA_POND'] = KAPPA_POND

        sith = self.data['sith']/self.data['sic']    

        
        self.data['h_ocn_SIC'] = self.data['insol']*(1-OCN_ALB)*(1-self.data['sic'])

        #self.data['hrfl_ocn_SIC'] = self.data['insol']*(1.-self.data['sic'])*OCN_ALB

        alb_type_list = ['AGE']
        if MYI_out: alb_type_list.append('MYI')
        if FYI_out: alb_type_list.append('FYI')
        if full_out: alb_tpe_list = ['AGE','MYI','FYI']

        for alb_type in alb_type_list:

            KAPPA_ALB = xr.where(self.data[f'alb_{alb_type}']>.7,KAPPA_SNOW,KAPPA_DEFT)
            KAPPA_ALB = xr.where(self.data[f'alb_{alb_type}']<.5,KAPPA_POND,KAPPA_ALB)
            ext_alb = np.exp(-KAPPA_ALB*sith)
            ext0 = np.exp(-KAPPA0*sith)

            
            self.data[f'h_nIC_{alb_type}'] = self.data['insol']*(1-self.data[f'alb_{alb_type}'])
            self.data[f'h_ice_{alb_type}'] = self.data[f'h_nIC_{alb_type}']*self.data['sic'] 

            self.data[f'hthr_ice_{alb_type}_KAPPA0'] = Y_0 * self.data[f'h_ice_{alb_type}'] * ext0
            self.data[f'hthr_ice_{alb_type}_KAPPA_ALB'] = Y_0 * self.data[f'h_ice_{alb_type}'] * ext_alb
                

            # self.data[f'hrfl_ice_{alb_type}'] = self.data['insol']*self.data['sic']*self.data[f'alb_{alb_type}']
            # self.data[f'habs_ice_top_{alb_type}'] = (1.-Y_0) * self.data[f'h_ice_{alb_type}'] 
            # self.data[f'habs_ice_dep_{alb_type}'] = Y_0 * self.data[f'h_ice_{alb_type}'] * (1.-ext_alb)
            # for vn in [f'hrfl_ice_{alb_type}',f'habs_ice_top_{alb_type}',\
            #            f'habs_ice_dep_{alb_type}',f'hthr_ice_{alb_type}']:
            #     self.data[f'accu_{vn}_KAPPA_ALB'] = self.cumsum(vn)*86400./1e6
            

        # for vn in ['h_nIC_MYI','h_nIC_FYI','h_nIC_AGE','h_ocn_SIC',
        #         'h_ice_MYI','h_ice_FYI','h_ice_AGE']:
        #     self.data[f'accu_{vn}'] = self.cumsum(vn)*86400./1e6 # W/m2 to J/m2 to MJ/m2
        # self.data['h_all_AGE'] = self.data['h_ice_AGE'] + self.data['h_ocn_SIC']
        # self.data['h_all_MYI'] = self.data['h_ice_MYI'] + self.data['h_ocn_SIC']
        # self.data['h_all_FYI'] = self.data['h_ice_FYI'] + self.data['h_ocn_SIC']  
        # self.data['h_ice_MYI'] = self.data['h_nIC_MYI']*self.data['sic']
        # self.data['h_ice_FYI'] = self.data['h_nIC_FYI']*self.data['sic']
        # self.data['h_nIC_MYI'] = self.data['insol']*(1-self.data['alb_MYI'])
        # self.data['h_nIC_FYI'] = self.data['insol']*(1-self.data['alb_FYI'])
        # self.data['accu_h_all_MYI'] = self.data['accu_h_ice_MYI'] + self.data['accu_h_ocn_SIC'] 
        # self.data['accu_h_all_FYI'] = self.data['accu_h_ice_FYI'] + self.data['accu_h_ocn_SIC'] 
        # self.data['accu_h_all_AGE'] = self.data['accu_h_ice_AGE'] + self.data['accu_h_ocn_SIC'] 
        

        # self.data['hthr_ice_MYI_KAPPA0'] = Y_0 * self.data['h_ice_MYI'] * ext0
        # self.data['hthr_ice_FYI_KAPPA0'] = Y_0 * self.data['h_ice_FYI'] * ext0


        # ----- For ponded ice ----
        # --> NOTE:: TO BE UPDATED!!
        if 'pond' in self.db.keys():
            ext_pond = np.exp(-self.data['kappa']*sith)
            self.data['hthr_ice_MYI_pond'] = Y_0 * self.data['h_ice_MYI'] * ext_pond
            self.data['hthr_ice_FYI_pond'] = Y_0 * self.data['h_ice_FYI'] * ext_pond
            self.data['hthr_ice_AGE_pond'] = Y_0 * self.data['h_ice_AGE'] * ext_pond
            for vn in ['hthr_ice_MYI_pond','hthr_ice_FYI_pond','hthr_ice_AGE_pond']:
                self.data[f'accu_{vn}'] = self.cumsum(vn)*86400./1e6         
        
            self.data['kappa'] = self.data['kappa'].where(self.data['land']==0,-.999)
            self.data['kappa_org'] = self.data['kappa_org'].where(self.data['land']==0,-.999)

        if full_out:
            # --- These variables can be reproduced just by cumsum or simple addition. No linearity involved. 
            # --- add all heat input to ocean/ice
            # --- add accumulation.  
            self.data['accu_h_ocn_SIC']  = self.cumsum('h_ocn_SIC')*86400./1e6
            for alb_type in alb_type_list:                           
                self.data[f'h_all_{alb_type}'] = self.data[f'h_ice_{alb_type}'] + self.data['h_ocn_SIC'] 
                for vn in [f'h_nIC_{alb_type}',f'h_ice_{alb_type}']:
                    self.data[f'accu_{vn}'] = self.cumsum(vn)*86400./1e6 # W/m2 to J/m2 to MJ/m2  
                self.data[f'accu_h_all_{alb_type}'] = self.data[f'accu_h_ice_{alb_type}'] + self.data['accu_h_ocn_SIC'] 
                for vn in [f'hthr_ice_{alb_type}_KAPPA_ALB', f'hthr_ice_{alb_type}_KAPPA0']:
                    self.data[f'accu_{vn}'] = self.cumsum(vn)*86400./1e6
            self._add_column_number()
            self._add_data_description()

        #self._add_fill_values()
        return #mpf_mask

    def to_file(self,outroot,ftype='nc',freq='Y'):
        '''
        '''

        assert ftype in ['txt','nc','csv', 'netcdf'], 'ftype has to be txt, csv, or nc.'

        txt_year = self.time[0].strftime(format='%Y')
        outdir = Path(outroot)/ftype
        if freq=='D': outdir = outdir/txt_year
        if not outdir.exists(): outdir.mkdir(parents=True,exist_ok=True)
        
        with open(outroot/f'BO_var_id_{txt_year}.json','w') as fid:
            json.dump(self.vdict,fid)

        ds_fill = self._add_fill_values()
        if ftype in ['txt','csv']:
            assert freq=='D', 'Output in txt or csv has to be daily.'
            vdict = self.vdict
            vn2ds = ['lat','lon','emelt','fmelt','efreeze','ffreeze']
            ncols = len(vdict.keys())
            nt = len(self.time)
            for it in range(nt):
                odf = pd.DataFrame()
                txt_date = self.time[it].strftime(format='%Y-%m-%d')
                ofn = outdir/f'BO_EASE2_gridded_25km_{txt_date}.{ftype}'
                if self.time[it].day==1: print(it,ofn)
                for icol in range(ncols):
                    vn = vdict[icol+2]
                    if vn in vn2ds:
                        dat = ds_fill[vn].values.ravel()
                    else:
                        dat = ds_fill[vn].values[it].ravel()
                    odf[vn] = dat
                odf.index += 1
                odat = np.column_stack([odf.index,odf.values])
                np.savetxt(ofn,odat,fmt='%6d\t'+'%10.4f '*ncols,delimiter='\t')
        else:
            fcomp = dict(zlib=True, complevel=5, dtype='float32')
            encoding = {var: fcomp for var in self.data.data_vars}
            for vn in self.data.data_vars:
                encoding[vn].update({'_FillValue':ds_fill[vn].attrs['_FillValue']})
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
                ods.to_netcdf(ofn,encoding=encoding)

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
                dat = self.data[vn].values.ravel()
            else:
                dat = self.data[vn].values[it].ravel()
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
    
    def load_frw_data(self,interp_sic=True,sic_vn='nsidc_bt_seaice_conc',sic_intp_fn=None):
        '''
        Load and compute FRW data. 

        Parameters:
        ----------
        sic_vn: str, default nsidc_bt_seaice_conc. The sic variable name in the CDR data. 
        interp_sic: bool, default True. The flag to interpolate missing sea ice data. 
        sic_intp_fn: path to the data file to interpolate the sic data, default None. 

        '''

        prd = 'era5'; vn = 'ssrd'
        with xr.open_dataset(self.db[prd]) as ds:
            self.data['insol'] = (['time','Y','X'], ds[vn].values/86400.)
            self.data['insol'].attrs['unit'] = 'W/m2'

        prd = 'sic'; vn = sic_vn 
        sic_vns = ['cdr_seaice_conc','nsidc_bt_seaice_conc','nsidc_nt_seaice_conc']       
        with xr.open_dataset(self.db[prd]) as ds:
            print(vn)
            self.data['sic'] = (['time','Y','X'], ds[vn].values)
            self.data[vn] = (['time','Y','X'], ds[vn].values)
            for vn in sic_vns:
                if vn!=sic_vn: self.data[vn] = (['time','Y','X'], ds[sic_vn].values)
        assert not( interp_sic and self.sic_intp_fn is None ), f'Set interp_sic to False or provide an sic_intp_fn.'
        if interp_sic:
            self._interpolate_sic(sic_intp_fn=sic_intp_fn)


                    
        return #ds


    def load_data(
            self,prd,
            prddir=None,prdvn=None,
            interp_sic=True,
            sic_intp_fn=None,
            ):

            '''
            Load a specific dataset. 

            Parameters:
            ----------
            prd: product name, including era5, sic, onset, iceage, piomas, albedo, and pond.
            prddir: the path to the data product, default None. 
            prdvn: the name of the variable to load. For sic, this can be cdr_seaice_conc, 
                   nsidc_bt_seaice_conc, nsidc_nt_seaice_conc, with default using BT. 
            interp_sic: bool, default True. The flag to interpolate missing sea ice data. 
            sic_intp_fn: path to the data file to interpolate the sic data, default None. 
            '''
        
            assert prd in self.db.keys() or prddir is not None, \
            f'Please provide directory for {prd} data on EASE2 grid.'
            if prddir is None: prddir = self.db[prd]
            
            nt = len(self.time)
            ny,nx = self.data.lat.shape

            if prd=='era5':
                if prdvn is None: prdvn = 'ssrd'
                with xr.open_dataset(prddir) as ds:
                    self.data['insol'] = (['time','Y','X'], ds[prdvn].values/86400.)
                    self.data['insol'].attrs['unit'] = 'W/m2'
            
            if prd=='sic':
                if prdvn is None: prdvn = 'nsidc_bt_seaice_conc' #prdvn = 'cdr_seaice_conc'
                with xr.open_dataset(prddir) as ds:
                    self.data['sic'] = (['time','Y','X'], ds[prdvn].values)
                    self.data['sic'].attrs['long name'] = prdvn
                if interp_sic:
                    self._interpolate_sic(sic_intp_fn=sic_intp_fn)   

            if prd == 'onset':
                with xr.open_dataset(prddir) as ds:   
                        self.data['emelt'] = (['Y','X'], ds['Earlymelt'].values)
                        self.data['fmelt'] = (['Y','X'], ds['Melt'].values)
                        self.data['efreeze'] = (['Y','X'], ds['Earlyfreeze'].values)
                        self.data['ffreeze'] = (['Y','X'], ds['Freeze'].values)
            
            if prd == 'iceage':
                if prdvn is None: prdvn = 'age_of_sea_ice'
                self.data['iceage'] = (['time','Y','X'], np.zeros((nt,ny,nx)))
                with xr.open_dataset(prddir) as ds:
                    for i_week in range(52):
                        i0 = i_week*7
                        i1 = i0 + 7
                        if i_week==51: i1 = nt
                        self.data['iceage'][i0:i1] =  ds[prdvn].values[i_week]

            if prd == 'albedo':
                with xr.open_dataset(prddir) as ds:
                    self.data['alb_FYI'] = (['time','Y','X'], ds['first_year'].values)
                    self.data['alb_MYI'] = (['time','Y','X'], ds['multi_year'].values)
                    self.data['alb_AGE'] = (['time','Y','X'], ds['actual_age'].values)

            if prd == 'piomas':
                if prdvn is None: prdvn = 'ice thickness'
                with xr.open_dataset(prddir) as ds:
                    self.data['sith'] = (['time','Y','X'], np.zeros((nt,ny,nx)) )
                    if self.time[0].is_leap_year:
                        self.data['sith'][:-1] = ds[prdvn].values
                        self.data['sith'][-1] = ds[prdvn].values[-1]
                    else:
                        self.data['sith'][:] = ds[prdvn].values
                self.data['sith'] = xr.where(self.data['sic']>0.15,self.data['sith']/self.data['sic'],0)

            if prd == 'pond':
                if prdvn is None: prdvn ='mpf'
                with xr.open_dataset(prddir) as ds:
                    self.data['mpf'] = (['time','Y','X'], np.zeros((nt,ny,nx)) )
                    if 'land' not in self.data.data_vars:
                        self.data['land'] = (['Y','X'],ds['land'][0].values) 
                    doy = ds.time.dt.dayofyear.values
                    for it, tdoy in enumerate(doy):
                        i0 = doy[it]-1
                        i1 = doy[it]+7
                        self.data['mpf'][i0:i1] = ds[prdvn].values[it]/100

                self.data['kappa_org'] = KA_ICE + self.data['mpf']*(KA_POND-KA_ICE)
                self.data['mpf_fill'] = self._fill_kappa('mpf')
                self.data['kappa'] = self._fill_kappa('kappa_org')

            return
    
    def load_data_all(
            self,
            albdir=None, onsetdir=None,iadir=None,sicdir=None,
            eradir=None,piodir=None,ponddir=None,
            sic_vn=None,
            interp_sic=True, sic_intp_fn=None):
        '''
        Load all dataset needed for BO.

        Parameters:
        ----------
        albdir: default None, the directory of the albedo data. 
        onsetdir: default None, the directory of the melt onset data. 
        iadir: default None, the directory of the ice age data. 
        sicdir: default None, the directory of the sic data. 
        eradir: default None, the directory of the era5 swdn data. 
        piodir: default None, the directory of the Piomas data. 
        ponddir: default None, the directory of the melt pond data. 
        interp_sic: bool, default True. The flag to interpolate missing data. 
        sic_intp_fn: default None, the path to the file to interpolate sic data. 
        '''

        prd_dir_dict = {
            'era5':eradir, 'sic':sicdir, 'onset': onsetdir,
            'iceage':iadir, 'albedo':albdir, 'piomas':piodir,
            'pond': ponddir,
        }
        prd_vn_dict = { x:None for x in prd_dir_dict.keys() }
        prd_vn_dict['sic'] = sic_vn

        for prd in prd_dir_dict.keys():
            if prd in self.db.keys() or prd_dir_dict[prd] is not None:
                prddir = prd_dir_dict[prd] 
                prdvn = prd_vn_dict[prd]
                if prd_dir_dict[prd] is None: prddir = self.db[prd]
                self.load_data(prd, prddir=prddir,prdvn=prdvn,
                               interp_sic=interp_sic,sic_intp_fn=sic_intp_fn)
            else:
                print(f'*** No data directory provided for {prd} data on EASE2 grid.')
                    
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
    
    def _interpolate_sic(self,sic_intp_fn=None):
        '''
        '''
        if sic_intp_fn is None: sic_intp_fn = self.sic_intp_fn
        assert sic_intp_fn is not None, 'Please provide the sic interpolation file.'
        assert Path(self.sic_intp_fn).exists(), f'{self.sic_intp_fn} does not exist.'
        assert 'sic' in self.data.data_vars, 'Please load SIC data first.'
        with xr.open_dataset(sic_intp_fn) as sicds:
            #print(sic_intp_fn)
            sicds.load()
        # --- make sure time stamp matches.
        t_in_sicds = self.data.time.isin(sicds.time)
        t_in_year = sicds.time.dt.year==self.year
        assert (self.data['time'][t_in_sicds]==sicds['time'][t_in_year]).all(), \
            f'Time stamps in SIC data does not match interpolation data for year {self.year}'
        self.data['sic'].loc[dict(time=t_in_sicds)] = sicds.sel(time=t_in_year)['sic']
        return
    
    def _check_sic_nan(self):
        '''
        Returns:
        -------
        missing: bool, if there are missing data in SIC
        df: pd.DataFrame, a list of dates and number of missing cell in SIC data. 
        '''
        assert 'sic' in self.data.data_vars, 'Please load SIC data first.'
        n_nans = self.data['sic'].isnull().sum(dim=['Y','X']).values 
        f_miss = n_nans!=NMISS0_SIC
        missing = f_miss.any()
        df = pd.DataFrame(dict(time=self.time[f_miss],n_miss=n_nans[f_miss]-NMISS0_SIC))
        return missing, df

    
    def _add_fill_values(self):
        '''
        '''
        ds_fill = self.data.copy(deep=True)
        for vn in self.data.data_vars:
            if vn in ['emelt','fmelt','efreeze','ffreeze','iceage']:
                #self.data[vn] = self.data[vn].where(~self.data[vn].isnull(),-99)
                ds_fill[vn] = ds_fill[vn].fillna(-99)
                ds_fill[vn].attrs['_FillValue'] = -99
                # self.data[vn].attrs['_FillValue'] = -99
            elif 'alb' in vn or vn in ['mpf','mpf_fill','kappa','kappa_fill','sic']:
                #self.data[vn] = self.data[vn].where(self.data[vn].isnull(),-.999)
                ds_fill[vn] = ds_fill[vn].fillna(-.999)
                ds_fill[vn].attrs['_FillValue'] = -.999
                # self.data[vn].attrs['_FillValue'] = -.999
            else:
                #self.data[vn] = self.data[vn].where(self.data[vn].isnull(),-999.99)
                ds_fill[vn] = ds_fill[vn].fillna(-999.99)
                ds_fill[vn].attrs['_FillValue'] = -999.99
                # self.data[vn].attrs['_FillValue'] = -999.99
        
        return ds_fill
    
    def _add_data_description(self):

        self.data['insol'].attrs['description'] = \
        'Surface downwelling radation [W/m2]'
        self.data['sic'].attrs['description'] = \
        'Sea ice concentration'

        self.data['emelt'].attrs['description'] = \
        'Early melt onset date'
        self.data['fmelt'].attrs['description'] = \
        'Full melt onset date'
        self.data['efreeze'].attrs['description'] = \
        'Early freeze onset date'
        self.data['ffreeze'].attrs['description'] = \
        'Full freeze onset date'
        self.data['iceage'].attrs['description'] = \
        'Sea ice age'

        self.data['alb_MYI'].attrs['description'] = \
        'Albedo of ice using the multi-year ice scheme'
        self.data['alb_FYI'].attrs['description'] = \
        'Albedo of ice using the first-year ice scheme'
        self.data['alb_AGE'].attrs['description'] = \
        'Albedo of ice using actual ice age'

        self.data['h_nIC_MYI'].attrs['description'] = \
        'heat input to ice, no IC dependence, using multi-year ice albedo'
        self.data['h_ice_MYI'].attrs['description'] = \
        'heat input to ice portion of the cell, using multi-year ice albedo'
        self.data['h_ocn_SIC'].attrs['description'] = \
        'heat input to open water portion of cell, using multi-year ice albedo'
        self.data['h_all_MYI'].attrs['description'] = \
        'total heat input to cell, using multi-year ice albedo'
        self.data['accu_h_nIC_MYI'].attrs['description'] = \
        'accumulated heat input to ice, no IC dependence, using multi-year ice albedo'
        self.data['accu_h_ice_MYI'].attrs['description'] = \
        'accumulated heat input to ice portion of the cell, using multi-year ice albedo'
        self.data['accu_h_ocn_SIC'].attrs['description'] = \
        'accumulated heat input to open water portion of cell, using multi-year ice albedo'
        self.data['accu_h_all_MYI'].attrs['description'] = \
        'accumulated total heat input to cell, using multi-year ice albedo'

        self.data['h_nIC_FYI'].attrs['description'] = \
        'heat input to ice, no IC dependence, using first-year ice albedo'
        self.data['h_ice_FYI'].attrs['description'] = \
        'heat input to ice portion of the cell, using first-year ice albedo'
        self.data['h_all_FYI'].attrs['description'] = \
        'total heat input to cell, using first-year ice albedo'
        self.data['accu_h_nIC_FYI'].attrs['description'] = \
        'accumulated heat input to ice, no IC dependence, using first-year ice albedo'
        self.data['accu_h_ice_FYI'].attrs['description'] = \
        'accumulated heat input to ice portion of the cell, using first-year ice albedo'
        self.data['accu_h_all_FYI'].attrs['description'] = \
        'accumulated total heat input to cell, using first-year ice albedo'

        self.data['h_nIC_AGE'].attrs['description'] = \
        'heat input to ice, no IC dependence, using ice albedo by ice age'
        self.data['h_ice_AGE'].attrs['description'] = \
        'heat input to ice portion of the cell, using ice albedo by ice age'
        self.data['h_all_AGE'].attrs['description'] = \
        'total heat input to cell, using ice albedo by ice age'
        self.data['accu_h_nIC_AGE'].attrs['description'] = \
        'accumulated heat input to ice, no IC dependence, using ice albedo by ice age'
        self.data['accu_h_ice_AGE'].attrs['description'] = \
        'accumulated heat input to ice portion of the cell, using ice albedo by ice age'
        self.data['accu_h_all_AGE'].attrs['description'] = \
        'accumulated total heat input to cell, using ice albedo by ice age'

        self.data['sith'].attrs['description'] = \
        'Sea ice thickness from PIOMAS'


        self.data['hthr_ice_MYI_KAPPA0'].attrs['description'] = \
        'heat through ice portion of the cell using multi-year albedo and constant kappa of 1'
        self.data['hthr_ice_FYI_KAPPA0'].attrs['description'] = \
        'heat through ice portion of the cell using first-year albedo and constant kappa of 1'
        self.data['hthr_ice_AGE_KAPPA0'].attrs['description'] = \
        'heat through ice portion of the cell using ice albedo by ice age and constant kappa of 1'
        self.data['accu_hthr_ice_MYI_KAPPA0'].attrs['description'] = \
        'accumulated heat through ice portion of the cell using multi-year albedo and constant kappa of 1'
        self.data['accu_hthr_ice_FYI_KAPPA0'].attrs['description'] = \
        'accumulated heat through ice portion of the cell using first-year albedo and constant kappa of 1'
        self.data['accu_hthr_ice_AGE_KAPPA0'].attrs['description'] = \
        'accumulated heat through ice portion of the cell using ice albedo by ice age and constant kappa of 1'

        
        self.data['hthr_ice_MYI_KAPPA_ALB'].attrs['description'] = \
        'heat through ice portion of the cell using multi-year albedo and kappa depending on ice albedo'
        self.data['hthr_ice_FYI_KAPPA_ALB'].attrs['description'] = \
        'heat through ice portion of the cell using first-year albedo and kappa depending on ice albedo'
        self.data['hthr_ice_AGE_KAPPA_ALB'].attrs['description'] = \
        'heat through ice portion of the cell using ice albedo by ice age and kappa depending on ice albedo'
        self.data['accu_hthr_ice_MYI_KAPPA_ALB'].attrs['description'] = \
        'accumulated heat through ice portion of the cell using multi-year albedo and kappa depending on ice albedo'
        self.data['accu_hthr_ice_FYI_KAPPA_ALB'].attrs['description'] = \
        'accumulated heat through ice portion of the cell using first-year albedo and kappa depending on ice albedo'
        self.data['accu_hthr_ice_AGE_KAPPA_ALB'].attrs['description'] = \
        'accumulated heat through ice portion of the cell using ice albedo by ice age and kappa depending on ice albedo'


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