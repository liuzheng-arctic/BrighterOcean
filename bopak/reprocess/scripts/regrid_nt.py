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

from bopak.reprocess import albedo, onset, sic, iceage,ntsic
from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from bopak.reprocess import set_BOEASE2


BTROOT = Path('/data/BO/NTSIC')
outdir = Path('/data/BO/EASE2/NTSIC/')
grid_fn = '/data/BO/CDR_SIC/G02202-cdr-ancillary-nh.nc'

df_cdr_nt = pd.read_csv('cdr_nt_sic_missing.csv',parse_dates=['time'])
rx = re.compile('(\d{8})')
fns_day = sorted( BTROOT.glob('*_N25km_'+'[0-9]'*8+'_v2.0.nc') )

tt = [rx.findall(x.stem).pop() for x in fns_day]
nt_db = pd.DataFrame(dict(time=pd.to_datetime(tt),fn=fns_day))


mds, boprj = set_BOEASE2()
# --- grid/projection info for sic and onset
sic_grid = sic.load_grid(grid_fn)
mds_sic, sicprj = sic.set_proj(sic_grid)

def check_nt(t0):
    fn = nt_db[nt_db.time==t0].iloc[0].fn
    with xr.open_dataset(fn) as ds:
        ds.load()

    fnd = False
    for vn in ds.data_vars:
        if 'ICECON' in vn:
            fnd = True
    return fnd

# t_all = []
# for i in range(len(df_cdr_nt)):
#     ti = df_cdr_nt.iloc[i].time
#     t_in_nt = check_nt(ti)
#     if t_in_nt:
#         if ti not in t_all: t_all.append(ti)
#     else:
#         t0 = ti - pd.Timedelta(1,'d')
#         t1 = ti + pd.Timedelta(1,'d')
#         if check_nt(t0) and check_nt(t1):
#             if t0 not in t_all: t_all.append(t0)
#             if t1 not in t_all: t_all.append(t1)


t_all = []
for i in range(len(df_cdr_nt)):
    ti = df_cdr_nt.iloc[i].time
    t0 = ti - pd.Timedelta(1,'d')
    t1 = ti + pd.Timedelta(1,'d')
    t_in_nt = check_nt(ti)
    t0_in_nt = check_nt(t0)
    t1_in_nt = check_nt(t1)
    if ti not in t_all and t_in_nt: t_all.append(ti)
    if t0 not in t_all and t0_in_nt: t_all.append(t0)
    if t1 not in t_all and t1_in_nt: t_all.append(t1)

ndays = len(t_all)
nt_db_regrid = nt_db[nt_db.time.isin(t_all)]
reuse_weight_sic = False
regridder_fn = Path('./.regridder/CDRNT_SIC2EASE2.nc')

fcomp = dict(zlib=True, complevel=5, dtype='float32')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)
#ndays = 2
for iday in range(ndays):
    fn = nt_db_regrid.iloc[iday].fn
    
    print(iday,fn)
    # dat_sic = sic(fn)
    # vn_sics =  [vn for vn in dat_sic.data.data_vars if vn!='crs']
    # assert len(vn_sics)==1,f'More than one sea ice data variable in {fn}'
    # vn_sic = vn_sics[0]
        
    # dat_sic.set_grid(ds_grid=sic_grid)
    # dat_sic.set_regridder(mds,reuse_weight = reuse_weight_sic)
    # reuse_weight_sic = True
    # dat_sic.regrid(mds=mds,sicprj=sicprj,vn_mask=vn_sic,mask_value=1.1)

    dat_nt = ntsic(fn)
    dat_nt.set_grid(ds_grid=sic_grid)
    dat_nt.set_regridder(mds,reuse_weight = reuse_weight_sic,regridder_fn=regridder_fn)
    reuse_weight_sic = True
    dat_nt.regrid(mds=mds,proj=sicprj)


    ods= dat_nt.data_ease2
    encoding = {var: fcomp for var in ods.data_vars}
    outfn = outdir/f'{fn.stem}_BOEASE2.nc'
    ods.to_netcdf(outfn,encoding=encoding)