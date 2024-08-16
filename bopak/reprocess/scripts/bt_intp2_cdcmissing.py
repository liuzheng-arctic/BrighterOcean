import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
from pathlib import Path
from datetime import datetime
import timeit
import re
import warnings
warnings.filterwarnings("ignore")



from bopak.reprocess import albedo, onset, sic, iceage
from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from bopak.reprocess import set_BOEASE2



BTROOT = Path('/data/BO/EASE2/BTSIC/')
df_cdc = pd.read_csv('cdc_sic_missing.csv',parse_dates=['time'])
rx = re.compile('(\d{8})')
fns_day = sorted( BTROOT.glob('*_N25km_'+'[0-9]'*8+'_v4.0_BOEASE2.nc') )

tt = [rx.findall(x.stem).pop() for x in fns_day]
bt_db = pd.DataFrame(dict(time=pd.to_datetime(tt),fn=fns_day))


with xr.open_dataset(bt_db.iloc[0].fn) as ds:
    
ods = xr.Dataset()
ods['time'] = df_cdc.time




print('here')

