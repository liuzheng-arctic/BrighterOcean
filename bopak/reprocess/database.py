import numpy as np
import pandas as pd
import re
from pathlib import Path


def albedo_db(ALB_ROOT):
    '''
    '''
    fns = sorted( ALB_ROOT.glob('*.nc') )
    rx = re.compile('(\d{4})')
    yrs = map(int,[rx.findall(x.name).pop() for x in fns])
    return pd.DataFrame(dict(year=yrs, fn=fns))

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