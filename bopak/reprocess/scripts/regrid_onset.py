#%%
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
from pathlib import Path

import warnings
warnings.filterwarnings("ignore")

from bopak.reprocess import albedo, onset, sic, iceage,btsic
from bopak.reprocess import onset_db,sic_db,ptproj, iceage_db
from bopak.reprocess import set_BOEASE2

#%%
SYEAR = 2023
EYEAR = 2024
outdir = Path('/data/BO/EASE2/onset')
ONSET_ROOT = Path('/data/BO/MeltFreezeOnset/')
ONSETNC_ROOT = Path('/data/BO/MeltFreezeOnset_NetCDF/')
grid_fn = '/data/BO/CDR_SIC/G02202-cdr-ancillary-nh.nc'
dbin = onset_db(ONSETNC_ROOT)
dbrun = dbin[(dbin.time.dt.year>=SYEAR)&(dbin.time.dt.year>=SYEAR)]


mds, boprj = set_BOEASE2()
# --- use the sic grid to set up projection because onset use the same grid
ds_grid = sic.load_grid(grid_fn)
mds_sic, sicprj = sic.set_proj(ds_grid)
reuse_weight_onset = False
fcomp = dict(zlib=True, complevel=5, dtype='float32')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)
# %%

nyears = len(dbrun)
#nyears = 1

for iyr in range(nyears):
    tyr = dbrun.iloc[iyr].time.year
    onset_fn = dbrun.iloc[iyr].fn
    dat_onset = onset(onset_fn,ds_grid=ds_grid)
    dat_onset.set_regridder(mds,reuse_weight = reuse_weight_onset)
    reuse_weight_onset = True
    dat_onset.regrid()
    encoding = {var: fcomp for var in dat_onset.data_ease2.data_vars}
    outfn = outdir/onset_fn.name
    dat_onset.data_ease2.to_netcdf(outfn,encoding=encoding)
