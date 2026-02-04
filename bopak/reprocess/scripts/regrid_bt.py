#%%
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

from bopak.reprocess import albedo, onset, sic, iceage,btsic
from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from bopak.reprocess import set_BOEASE2

#%%


SYEAR = 2023
EYEAR = 2024


BTROOT = Path('/data/BO/BTSIC')
outdir = Path('/data/BO/EASE2/BTSIC/20234')
grid_fn = '/data/BO/CDR_SIC/G02202-cdr-ancillary-nh.nc'

rx = re.compile('(\d{8})')
fns_day = sorted( BTROOT.glob('*_N25km_'+'[0-9]'*8+'_v4.0.nc') )

tt = [rx.findall(x.stem).pop() for x in fns_day]
bt_db = pd.DataFrame(dict(time=pd.to_datetime(tt),fn=fns_day))


mds, boprj = set_BOEASE2()
# --- grid/projection info for sic and onset
sic_grid = sic.load_grid(grid_fn)
mds_sic, sicprj = sic.set_proj(sic_grid)
reuse_weight_sic = False

fcomp = dict(zlib=True, complevel=5, dtype='float32')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)


#%%
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
    df_year = bt_db[bt_db.time.dt.year==tyr]
    nday = len(df_year)
    ods = xr.Dataset()
    #ods['time'] = df_year.time
    dslist = []
    #nday = 50
    for iday in range(nday):
        fn = df_year.iloc[iday].fn
        if df_year.iloc[iday].time.day==1: 
            print('-----------------------------')
            print(fn)
            print('-----------------------------')
        dat_sic = btsic(fn)
        dat_sic.set_grid(ds_grid=sic_grid)
        dat_sic.set_regridder(mds,reuse_weight = reuse_weight_sic)
        reuse_weight_sic = True
        dat_sic.regrid(mds=mds,proj=sicprj)
        dslist.append(dat_sic.data_ease2)
    ods = xr.concat(dslist,dim='time')


    encoding = {var: fcomp for var in ods.data_vars}
    outfn = outdir/f'NSIDC0079_SEAICE_PS_N25km_{yr_txt}_v4.0_BOEASE2.nc'
    #outfn = outdir/f'seaice_conc_daily_nh_{yr_txt}_f17_v04r00.nc'
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


#%%

# ndays = len(t_all)
# bt_db_regrid = bt_db#[bt_db.time.isin(t_all)]
# ndays = len(bt_db_regrid)
# reuse_weight_sic = True


# fcomp = dict(zlib=True, complevel=5, dtype='float32')
# if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)
# #ndays = 2

# #%%
# for iday in range(ndays):
#     fn = bt_db_regrid.iloc[iday].fn
    
#     print(iday,fn)
#     # dat_sic = sic(fn)
#     # vn_sics =  [vn for vn in dat_sic.data.data_vars if vn!='crs']
#     # assert len(vn_sics)==1,f'More than one sea ice data variable in {fn}'
#     # vn_sic = vn_sics[0]
        
#     # dat_sic.set_grid(ds_grid=sic_grid)
#     # dat_sic.set_regridder(mds,reuse_weight = reuse_weight_sic)
#     # reuse_weight_sic = True
#     # dat_sic.regrid(mds=mds,sicprj=sicprj,vn_mask=vn_sic,mask_value=1.1)

#     dat_bt = btsic(fn)
#     dat_bt.set_grid(ds_grid=sic_grid)
#     dat_bt.set_regridder(mds,reuse_weight = reuse_weight_sic)
#     reuse_weight_sic = True
#     dat_bt.regrid(mds=mds,proj=sicprj)


#     ods= dat_bt.data_ease2
#     encoding = {var: fcomp for var in ods.data_vars}
#     outfn = outdir/f'{fn.stem}_BOEASE2.nc'
#     ods.to_netcdf(outfn,encoding=encoding)