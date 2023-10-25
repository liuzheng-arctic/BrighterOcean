#%%
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
from pathlib import Path
from datetime import datetime
import timeit


from bopak.reprocess import albedo, onset, sic, iceage
from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from bopak.reprocess import set_BOEASE2


# %%
SYEAR = 1984
EYEAR = 2022
SYEAR = 2022
EYEAR = 2022
alb_version = 'V0.04'


outdir = Path('/data/BO/EASE2/albedo/{alb_version}')
ONSET_OUT = Path(f'/data/BO/EASE2/onset')
IA_OUT = Path('/data/BO/EASE2/ICEAGE/V4')

SIC_ROOT = Path('/data/BO/CDR_SIC/')
ONSETNC_ROOT = Path('/data/BO/MeltFreezeOnset_NetCDF/')
IA_ROOT = Path('/data/BO/ICEAGE/V4')
grid_fn = '/data/BO/CDR_SIC/G02202-cdr-ancillary-nh.nc'
# %%

# --- set up database
onset_df = onset_db(ONSETNC_ROOT)
sic_df = sic_db(SIC_ROOT)
ia_df = iceage_db(IA_ROOT)
# %%
fcomp = dict(zlib=True, complevel=5, dtype='float32')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)
if not ONSET_OUT.exists(): ONSET_OUT.mkdir(parents=True, exist_ok=True)
if not IA_OUT.exists(): IA_OUT.mkdir(parents=True, exist_ok=True)
# %%

# --- setup map for EASE2
mds, boprj = set_BOEASE2()
# --- grid/projection info for sic and onset
sic_grid = sic.load_grid(grid_fn)
mds_sic, sicprj = sic.set_proj(sic_grid)

reuse_weight_onset = False
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
    t_daily_yr = pd.date_range(yr_txt+'0101',yr_txt+'1231',freq='D')
    nday = len(t_daily_yr)
    #nday = 
    
    fn_onset = onset_df[onset_df.time.dt.year==tyr].iloc[0].fn
    fn_ia = ia_df[ia_df.time.dt.year==tyr].iloc[0].fn

    dat_ia = iceage(fn_ia)
    dat_onset = onset(fn_onset)
    dat_onset.set_grid(ds_grid=sic_grid)

    # --- Setup regridder. Generate the regridder file in the first time, 
    #     reuse in the following years.
    dat_onset.set_regridder(mds,reuse_weight = reuse_weight_onset)
    dat_ia.set_regridder(mds,method='bilinear',reuse_weight = reuse_weight_onset)
    reuse_weight_onset = True
    reuse_weight_ia = True
    dat_onset.regrid(mds=mds,proj=sicprj)
    dat_ia.regrid(mds=mds)

    fn_onset_out = ONSET_OUT/fn_onset.name
    fn_ia_out = IA_OUT/fn_ia.name
    encoding_onset = {var: fcomp for var in dat_onset.data_ease2.data_vars}
    encoding_ia = {var: fcomp for var in dat_ia.data_ease2.data_vars}
    dat_onset.data_ease2.to_netcdf(fn_onset_out,encoding=encoding_onset)
    dat_ia.data_ease2.to_netcdf(fn_ia_out,encoding=encoding_ia)
 

    # --- load init for albedo from previous year
    #     Not really useful because it is not continuous by design. 
    # if load_init:
    #     if init_fn is None:
    #         init_fn = outdir/f'albedo_BO_{tyr-1:4d}_{alb_version}.nc'
    #     assert init_fn.exists(), f'Need a valid init_fn to load_init for albedo. {init_fn} not found.'
    #     alb_init = albedo.load_init(init_fn)

    # --- setup albedo for processing
    alb = albedo(yr_txt, onset_date=dat_onset.data_ease2,mds=mds, boprj=boprj)

    for iday in range(nday):

        t_day = t_daily_yr[iday]

        iceage_iday = dat_ia.interp(t_day,grid='EASE2')
        
        tdat_fy = alb.first_year(t_day)
        tdat_my = alb.multi_year(t_day)
        tdat_ia = alb.actual_age(t_day, iceage_iday)
        
        if t_day.day==1:
            print(iday, t_day)
    
    # --- Not used: for continous evolution. 
    # if iday==nday-1:
    #     alb_init = alb.alb_fnl

    encoding = {var: fcomp for var in alb.data.data_vars}
    outfn = outdir/f'albedo_BO_{yr_txt}_{alb_version}.nc'
    alb.data.to_netcdf(outfn,encoding=encoding)

    t_end = timeit.default_timer()
    print(f'Saved to {outfn} .')
    print(f'Year {yr_txt}: finished in {t_end-t_start}s')
    t_start = t_end

t_end = timeit.default_timer()
print(f'{SYEAR:4d} to {EYEAR:4d}: finished in {(t_end-t_start0)/60:.1f} minutes.')
txt_now   = datetime.now().strftime(format='%Y-%m-%dT%H:%M:%S')
print(f'Time now: {txt_now}')
print('<<< =========================================')

# %%
