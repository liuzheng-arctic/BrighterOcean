import numpy as np
import pandas as pd
import re
from pathlib import Path


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