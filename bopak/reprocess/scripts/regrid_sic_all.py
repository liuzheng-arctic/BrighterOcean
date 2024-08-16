#%%
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
from pathlib import Path
from datetime import datetime
import timeit

import warnings
warnings.filterwarnings("ignore")

from bopak.reprocess import albedo, onset, sic, iceage
from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from bopak.reprocess import set_BOEASE2


# %%
SYEAR = 1984
EYEAR = 2022
# SYEAR = 1984
# EYEAR = 2022
# SYEAR = 1979
# EYEAR = 1983


outdir = Path('/data/BO/EASE2/SIC/cdr_all')


SIC_ROOT = Path('/data/BO/CDR_SIC/')
grid_fn = '/data/BO/CDR_SIC/G02202-cdr-ancillary-nh.nc'
# %%

# --- set up database
sic_df = sic_db(SIC_ROOT)
# %%
fcomp = dict(zlib=True, complevel=5, dtype='float32')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)
# %%

# --- setup map for EASE2
mds, boprj = set_BOEASE2()
# --- grid/projection info for sic and onset
sic_grid = sic.load_grid(grid_fn)
mds_sic, sicprj = sic.set_proj(sic_grid)

reuse_weight_sic = False


# --- not used. Not continous evolution by design
#     init info not used. 
# load_init = False
# init_fn = None
# alb_init = None

txt_now   = datetime.now().strftime(format='%Y-%m-%dT%H:%M:%S')
print('========================================= >>>')
print(f'Time now: {txt_now}')
t_start = timeit.default_timer()
t_start0 = t_start
for iyr, tyr in enumerate(range(SYEAR, EYEAR+1)):

    yr_txt = str(tyr)

    print(f'Processing year: {yr_txt}')
    # t_daily_yr = pd.date_range(yr_txt+'0101',yr_txt+'1231',freq='D')
    # nday = len(t_daily_yr)
    df_year = sic_df[sic_df.time.dt.year==tyr]
    nday = len(df_year)
    ods = xr.Dataset()
    #ods['time'] = df_year.time
    
    for iday in range(nday):
        fn = df_year.iloc[iday].fn
        if df_year.iloc[iday].time.day==1: 
            print('-----------------------------')
            print(fn)
            print('-----------------------------')
        dat_sic = sic(fn)
        dat_sic.set_grid(ds_grid=sic_grid)
        dat_sic.set_regridder(mds,reuse_weight = reuse_weight_sic)
        reuse_weight_sic = True
        dat_sic.regrid(mds=mds,sicprj=sicprj)
        if iday==0:
            ods= dat_sic.data_ease2
        else:
            ods = xr.concat([ods,dat_sic.data_ease2],dim='tdim')


    encoding = {var: fcomp for var in ods.data_vars}
    outfn = outdir/f'seaice_conc_daily_nh_{yr_txt}_f17_v04r00.nc'
    ods.to_netcdf(outfn,encoding=encoding)

    # t_end = timeit.default_timer()
    # print(f'Saved to {outfn} .')
    # print(f'Year {yr_txt}: finished in {t_end-t_start}s')
    # t_start = t_end

t_end = timeit.default_timer()
print(f'{SYEAR:4d} to {EYEAR:4d}: finished in {(t_end-t_start0)/60:.1f} minutes.')
txt_now   = datetime.now().strftime(format='%Y-%m-%dT%H:%M:%S')
print(f'Time now: {txt_now}')
print('<<< =========================================')

# %%
