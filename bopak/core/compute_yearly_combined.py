#%%
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
from pathlib import Path

from bopak import set_BOEASE2,  bodata3d, bodb
from bopak.reprocess.utils import loc_ease2

mds,boprj = set_BOEASE2()
pcproj = ccrs.PlateCarree()
XLIM = mds.XLIM
YLIM = mds.YLIM
xx,yy = np.meshgrid(mds.X.values, mds.Y.values)


gridfn_ease2 = '/data/BO/EASE2/regridded_G02202-cdr-ancillary-nh.nc'
BOROOT =  Path('/data/BO/EASE2/combined')
CASEID = 'kappa_alb_snow3_y058_part_albV006_T001'
bo_in_dir = BOROOT/CASEID/'nc'
#BOROOT = Path('/data/BO/EASE2/combined/kappa_alb_snow3_y058_part_albV006_T001/nc')
datdb = bodb(bo_in_dir,ftype='nc',freq='Y')



outdir = BOROOT/f'yearly/{CASEID}'
if not outdir.exists(): outdir.mkdir(parents=True,exist_ok=True)
#outdir = BOROOT.parent.parent/'yearly/kappa_alb_snow3_y058_refl_skipnaF'
ofn = outdir/'BO_yearly.nc'
ofn = outdir/'BO_yearly_2yr_example.nc'
#%%
dslist = []
#tdb = datdb[datdb.time.dt.year.isin([1998,2020])]
#tdb = datdb[datdb.time.dt.year.isin([1987])]
tdb = datdb
for i in range(len(tdb)):
#for i in range(2):
    tbo = bodata3d(tdb.iloc[i].fn)
    ods = bodata3d.annual_solar_partition(tbo.fn,tbo.year,gridfn=gridfn_ease2)
    dslist.append(ods)
    print(tdb.iloc[i].fn)
ovns = [vn for vn in dslist[0].data_vars if vn not in ['lat','lon']]
dslist_nocoord = [x[ovns] for x in dslist]
dsyearly = xr.concat(dslist_nocoord ,dim='time')
fcomp = dict(zlib=True, complevel=5, dtype='float32')
encoding = {var: fcomp for var in dsyearly.data_vars}
dsyearly['lat'] = dslist[0]['lat']#[0]
dsyearly['lon'] = dslist[0]['lon']#[0]
dsyearly.to_netcdf(ofn,encoding=encoding)
