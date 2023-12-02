import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import json
import re
from .reprocess import set_BOEASE2

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
    def __init__(self,fn,ds=None):
        self.fn = Path(fn)
        self.year = int(self.fn.stem.split('_')[-1])
        if ds is None:
            self.load_data(fn)
        else:
            self.ds = ds
    
        return 
    def load_data(self,fn):
        '''
        '''
        with xr.open_dataset(fn) as ds:
            self.data = ds
        return


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