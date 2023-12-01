import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import datetime
import timeit

import warnings
warnings.filterwarnings("ignore")

from bopak.reprocess import pond
from bopak.reprocess import pond_db
from bopak.reprocess import set_BOEASE2

SYEAR = 2000
EYEAR = 2011

outdir = Path('/data/BO/EASE2/pond')
              
MPF_ROOT = Path('/data/BO/MODIS_MeltPonds/')
pond_df = pond_db(MPF_ROOT)

fcomp = dict(zlib=True, complevel=5, dtype='float32')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)


# --- setup map for EASE2
mds, boprj = set_BOEASE2()

reuse_weight = False

txt_now   = datetime.now().strftime(format='%Y-%m-%dT%H:%M:%S')
print('========================================= >>>')
print(f'Time now: {txt_now}')
t_start = timeit.default_timer()
t_start0 = t_start
for iyr, tyr in enumerate(range(SYEAR, EYEAR+1)):
    yr_txt = str(tyr)

    print(f'Processing year: {yr_txt}')
    df_year = pond_df[pond_df.time.dt.year==tyr]
    nday = len(df_year)
    ods = xr.Dataset()
    print('-----------------------------')
    for iday in range(nday):
        fn = df_year.iloc[iday].fn
        print(fn)
        dat = pond(fn)
        dat.set_regridder(mds,reuse_weight = reuse_weight)
        reuse_weight = True
        dat.regrid()
        if iday==0:
            ods= dat.data_ease2
        else:
            ods = xr.concat([ods,dat.data_ease2],dim='time')

    print('-----------------------------')
    encoding = {var: fcomp for var in ods.data_vars}
    outfn = outdir/'_'.join(fn.stem.split('__')[:-2]+[f'EASE2_25km_{yr_txt}.nc'])
    ods.to_netcdf(outfn,encoding=encoding)
print('Done!')