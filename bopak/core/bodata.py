import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import json
import re
from .utils import valid_domain
from ..reprocess import set_BOEASE2
from ..reprocess.albedo import OCN_ALB
from ..reprocess.bo_combine import Y_0

class bodata_daily:
    def __init__(self,fn,ftype='txt',fill_value=True,mds=None,vdict=None,vdictfn=None):
        '''
        '''
        self.fn = Path(fn)
        self.ftype = ftype
        self.time = pd.to_datetime(fn.stem.split('_')[-1])

        vdictfn0 = Path(__file__).resolve().parent/'BO_var_id.json'
        if vdict is not None:
            self.vdict = vdict
        elif vdictfn is not None:
            with open(vdictfn) as fid:
                self.vdict = json.load(fid)
        elif vdictfn0.exists():
            with open(vdictfn0) as fid:
                self.vdict = json.load(fid)


        if fn.exists(): 
            assert hasattr(self,'vdict'), 'Need variable mapping info, provide vdict or vdictfn.'
            self.load_data(fill_value=fill_value,mds=mds)
        return

    def load_data(self,fill_value=True,mds=None):
        '''
        '''
        vdict = self.vdict
        if mds is None: mds, boprj = set_BOEASE2()
        if self.ftype=='txt':
            dat = np.loadtxt(self.fn)
            ds = xr.Dataset()
            ds['lat'] = mds.lat
            ds['lon'] = mds.lon
            shp2d = ds.lat.shape
            for vkey in vdict.keys():
                vn = vdict[vkey]      
                icol = int(vkey)-1
                if icol<len(dat.T):
                    tdat = dat.T[icol].reshape(shp2d)
                    if fill_value:        
                        if vn in ['emelt','fmelt','efreeze','ffreeze','iceage']:  
                            tdat[tdat==-99] = np.nan
                        elif 'alb' in vn or vn in ['mpf','mpf_fill','kappa','kappa_fill','sic']:
                            tdat[tdat==-.999] = np.nan
                        else:
                            tdat[tdat==-999.99] = np.nan
                    ds[vn] = (['Y','X'],tdat) 
            ds = ds.expand_dims(time=[self.time])
            self.data = ds
        return

class bodata3d:

    def __init__(self,fn=None,ds=None,vns=None):
        '''
        '''
        if fn is not None:
            self.fn = Path(fn)
            self.year = int(self.fn.stem.split('_')[-1])
        # if ds is None:
        #     self.load_data(fn,vns=vns)
        # else:
        #     self.ds = ds    
        return 
    
    def load_data(self,fn,vns=None):
        '''
        '''
        print(fn)
        with xr.open_dataset(fn) as ds:
            if vns is None: vns = list(ds.data_vars)
            self.data = ds[vns]#.load()
        return
    
    def add_cumsum(self,vns=None,flux=True):
        '''
        '''
        if vns is None: vns = self.data.data_vars
        for vn in vns:
            self.data[f'accu_{vn}'] = self.data[vn].cumsum(dim='time')
        return 
    

    
    @staticmethod
    def annual_solar_partition(fn,yr_in,gridfn=None,vns=None, full_data=False,skipna=True):
        '''
        '''

        ods = xr.Dataset()
        with xr.open_dataset(fn) as ds:
            for vn in ['h_ocn_SIC','h_ice_AGE','hthr_ice_AGE_KAPPA_ALB']:
                assert vn in ds.data_vars, print(f'{vn} is not found in data: {fn}')
            
            # if 'OCN_ALB' in ds.attrs.keys(): OCN_ALB = ds.attrs['OCN_ALB']
            # if 'Y_0' in ds.attrs.keys(): OCN_ALB = ds.attrs['Y_0']


            ods.attrs = ds.attrs
            ods['sic'] = ds['sic'].mean(dim='time',skipna=skipna)
            ods['sith'] = ds['sith'].mean(dim='time',skipna=skipna)
            ods['sith_norm'] = (ds['sith']/ds['sic']).mean(dim='time',skipna=skipna)
            ods['alb_AGE'] = ds['alb_AGE'].mean(dim='time',skipna=skipna)
            ods['insol'] = ds['insol'].mean(dim='time',skipna=skipna)
            ods['fmelt'] = ds['fmelt']
            ods['ffreeze'] = ds['ffreeze']
            ods['emelt'] = ds['emelt']
            ods['efreeze'] = ds['efreeze']
            ods['lat'] = ds['lat']
            ods['lon'] = ds['lon']
            
            ds['hrfl_ocn_SIC'] = ds['insol']*(1.-ds['sic'])*OCN_ALB
            ds['habs_ocn_SIC'] = ds['insol']*(1.-ds['sic'])*(1.-OCN_ALB)
            ds['hall_ice_SIC'] = ds['insol']*ds['sic']
            alb_type_list = ['AGE']
            if full_data: alb_type_list = ['AGE','MYI','FYI']
            for alb_type in alb_type_list:
                ds[f'hrfl_ice_{alb_type}'] = ds['insol'] * ds['sic'] * ds[f'alb_{alb_type}']
                ds[f'habs_ice_top_{alb_type}'] = (1.-Y_0) * ds[f'h_ice_{alb_type}'] 
                ds[f'habs_ice_dep_{alb_type}_KAPPA_ALB'] = Y_0 * ds[f'h_ice_{alb_type}'] - xr.where(ds['sic']>0,ds[f'hthr_ice_{alb_type}_KAPPA_ALB'],0)
                ds[f'habs_ice_dep_{alb_type}_KAPPA0'] = Y_0 * ds[f'h_ice_{alb_type}'] - xr.where(ds['sic']>0,ds[f'hthr_ice_{alb_type}_KAPPA0'],0)
                for vn in [
                    f'h_nIC_{alb_type}',f'h_ice_{alb_type}',
                    f'hthr_ice_{alb_type}_KAPPA_ALB', f'hthr_ice_{alb_type}_KAPPA0',
                    f'hrfl_ice_{alb_type}',f'habs_ice_top_{alb_type}',
                    f'habs_ice_dep_{alb_type}_KAPPA_ALB',f'habs_ice_dep_{alb_type}_KAPPA0',
                    ]: 
                    ods[f'accu_{vn}'] = ds[vn].sum(dim='time',skipna=skipna)*86400./1e6 # W/m2 to J/m2 to MJ/m2  
            ods['accu_hall_ice_SIC'] = ds['hall_ice_SIC'].sum(dim='time',skipna=skipna)*86400./1e6 # W/m2 to J/m2 to MJ/m2  
            ods['accu_h_ocn_SIC'] = ds['h_ocn_SIC'].sum(dim='time',skipna=skipna)*86400./1e6 # W/m2 to J/m2 to MJ/m2  
            ods['accu_habs_ocn_SIC'] = ds['habs_ocn_SIC'].sum(dim='time',skipna=skipna)*86400./1e6 # W/m2 to J/m2 to MJ/m2  
            ods['accu_hrfl_ocn_SIC'] = ds['hrfl_ocn_SIC'].sum(dim='time',skipna=skipna)*86400./1e6 # W/m2 to J/m2 to MJ/m2  
            ods['accu_insol'] = ds['insol'].sum(dim='time',skipna=skipna)*86400./1e6 # W/m2 to J/m2 to MJ/m2  

            if gridfn is not None:
                valid_range = valid_domain(yr_in,gridfn)
                ods = ods.where(valid_range)
            ods = ods.assign_coords(time=ds.time[0].dt.year.values)
            ods = ods.expand_dims("time")
        return ods
    
    @staticmethod
    def database(boroot,ftype='txt',freq='D'):
        '''
        '''
        if freq=='D':
            fns_all = sorted( boroot.glob(f'**/*.{ftype}') )
            rx = re.compile('\d{4}-\d{2}-\d{2}')
        else:
            fns_all = sorted( boroot.glob(f'*.nc') )
            rx = re.compile('\d{4}')
        fns = []
        dates = []
        for fn in fns_all:
            if rx.search(fn.stem):
                fns.append(fn)
                dates.append(rx.findall(fn.stem)[-1])
        odf = pd.DataFrame(dict(time=pd.to_datetime(dates),fn=fns))

        return odf
        


def bodb(boroot,ftype='txt',freq='D'):
    '''
    '''
    if freq=='D':
        fns_all = sorted( boroot.glob(f'**/*.{ftype}') )
        rx = re.compile('\d{4}-\d{2}-\d{2}')
    else:
        fns_all = sorted( boroot.glob(f'*.nc') )
        rx = re.compile('\d{4}')
    fns = []
    dates = []
    for fn in fns_all:
        if rx.search(fn.stem):
            fns.append(fn)
            dates.append(rx.findall(fn.stem)[-1])
    odf = pd.DataFrame(dict(time=pd.to_datetime(dates),fn=fns))

    return odf