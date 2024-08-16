import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import xesmf
import timeit
from pathlib import Path


from bopak.reprocess import set_BOEASE2, ptproj


SBND = 35
fERAroot = Path('/data/BO/ERA5/daily')
outdir = Path('/data/BO/EASE2/ERA5/')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)
fns = sorted( fERAroot.glob('*.nc') )

# --- setting up regridder and projection 
mds, boprj = set_BOEASE2()
pcproj = ccrs.PlateCarree()

with xr.open_dataset(fns[0]) as ds:
    nds = ds.where(ds.latitude>=SBND,drop=True)

lon,lat = np.meshgrid(nds.longitude.values,nds.latitude.values)
dsin = xr.Dataset(coords=dict(
    lat=(['y','x'],lat), lon=(['y','x'],lon) ))
dsout = xr.Dataset(coords=dict(
    lat=(['Y','X'],mds.lat.values),
    lon=(['Y','X'],mds.lon.values)))

reuse_weight = False
method = 'bilinear'
regridder_fn = Path('./.regridder/ERA5daily.nc')
regridder = xesmf.Regridder(
            dsin,dsout,method,reuse_weights=reuse_weight,filename=regridder_fn)

reuse_weight = True
fcomp = dict(zlib=True, complevel=5, dtype='float32')
t_start = timeit.default_timer()
#for ifn,fn in enumerate(fns):
for ifn,fn in enumerate([fns[32]]):
    
    outfn = outdir/fn.name
    with xr.open_dataset(fn) as ds:
        nds = ds.where(ds.latitude>=SBND,drop=True)
    ods = regridder(nds)
    encoding = {var: fcomp for var in ods.data_vars}
    ods.to_netcdf(outfn,encoding=encoding)
    delta_time = timeit.default_timer() - t_start
    t_start = timeit.default_timer()
    print(f'Time for {ifn:2d} {fn.stem}: {delta_time}s.')
