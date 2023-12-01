import numpy as np
import pandas as pd
import re
from pathlib import Path


def albedo_db(ALB_ROOT):
    '''
    '''
    fns = sorted( ALB_ROOT.glob('*.nc') )
    rx = re.compile('(\d{4})')
    yrs = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(yrs), fn=fns))

def iceage_db(ICEAGE_ROOT):
    '''
    '''
    fns = sorted(Path(ICEAGE_ROOT).glob('*.nc'))
    rx = re.compile('(\d{8})')
    yrs = [rx.findall(x.name).pop()[:4] for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(yrs),fn=fns))

def onset_db_bin(ONSET_ROOT):
    '''
    '''
    melt_fns = sorted( Path(ONSET_ROOT).glob('*731smelt.*') )
    freeze_fns = sorted( Path(ONSET_ROOT).glob('*731sfreeze.*') )
    return melt_fns, freeze_fns



def onset_db(ONSET_ROOT):
    '''
    '''
    fns = sorted( Path(ONSET_ROOT).glob('*731smelt*.nc') )
    yrs = [x.name[:4] for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(yrs),fn=fns))

def piomas_db_ease2(PIOMAS_ROOT):
    '''
    '''

    fns = sorted( Path(PIOMAS_ROOT).glob('hiday*.nc') )
    rx = re.compile('H(\d{4})')
    yrs = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(yrs),fn=fns))
    
def piomas_db(PIOMAS_ROOT):
    '''
    '''

    fns = sorted( (Path(PIOMAS_ROOT)/'hiday').glob('*.gz') )
    rx = re.compile('H(\d{4})')
    yrs = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(yrs),fn=fns))


def sic_db(SIC_ROOT):
    '''
    '''
    fns = sorted( Path(SIC_ROOT).glob('[0-9]'*4+'/*.nc') )
    rx = re.compile('(\d{8})')
    tt = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(tt),fn=fns))

def sic_db_ease2(SIC_ROOT):
    '''
    '''
    fns = sorted( Path(SIC_ROOT).glob('*'+'[0-9]'*4+'*.nc') )
    rx = re.compile('(\d{4})')
    tt = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(tt),fn=fns))

def era5_db(ERA_ROOT):
    '''
    '''
    
    fns = sorted( Path(ERA_ROOT).glob('ERA*'+'[0-9]'*4+'*.nc') )
    rx = re.compile('(\d{4})')
    tt = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(tt),fn=fns))

def pond_db(MP_ROOT):
    '''
    '''
    fns = sorted( Path(MP_ROOT).glob('MODIS*'+'[0-9]'*7+'*.nc') )
    rx = re.compile('(\d{7})')
    tt = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(tt,format='%Y%j'),fn=fns))

def pond_db_ease2(MP_ROOT):
    '''
    '''
    fns = sorted( Path(MP_ROOT).glob('MODIS*'+'[0-9]'*4+'*.nc') )
    rx = re.compile('(\d{4})')
    tt = [rx.findall(x.name).pop() for x in fns]
    return pd.DataFrame(dict(time=pd.to_datetime(tt),fn=fns))


class combined_db:
    def __init__(
        self,
        albdir=None,
        sicdir=None,
        piodir=None,
        eradir=None,
        onsetdir=None,
        iadir=None,
        ponddir=None
    ):
        dbdict = {}
        if albdir is not None:
            dbdf = albedo_db(albdir)
            assert len(dbdf)>0, f'Cannot find matching files under {albdir}'
            dbdict['albedo'] = dbdf
        if sicdir is not None:
            dbdf = sic_db_ease2(sicdir)
            assert len(dbdf)>0, f'Cannot find matching files under {sicdir}'
            dbdict['sic'] = dbdf
        if piodir is not None:
            dbdf = piomas_db_ease2(piodir)
            assert len(dbdf)>0, f'Cannot find matching files under {piodir}'
            dbdict['piomas'] = dbdf
        if eradir is not None:
            dbdf = era5_db(eradir)
            assert len(dbdf)>0, f'Cannot find matching files under {eradir}'
            dbdict['era5'] = dbdf
        if onsetdir is not None:
            dbdf = onset_db(onsetdir)
            assert len(dbdf)>0, f'Cannot find matching files under {onsetdir}'
            dbdict['onset'] = dbdf
        if iadir is not None:
            dbdf = iceage_db(iadir)
            assert len(dbdf)>0, f'Cannot find matching files under {iadir}'
            dbdict['iceage'] = dbdf
        if ponddir is not None:
            dbdf = pond_db_ease2(ponddir)
            assert len(dbdf)>0, f'Cannot find matching files under {ponddir}'
            dbdict['pond'] = dbdf
        self.dbdict = dbdict
        return

    def query(self,year):
        '''
        '''
        dbdict = {}
        for vn in self.dbdict.keys():
            tdf = self.dbdict[vn]
            dfyr = tdf[tdf.time.dt.year==year]
            if len(dfyr)==0:
                print(f'Warning: no data found in year {year} for {vn}.')
            else:
                dbdict[vn] = dfyr.fn.values[0]
        return dbdict




