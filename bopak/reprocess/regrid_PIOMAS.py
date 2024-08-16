import numpy as np
import xarray as xr
import pandas as np
import cartopy.crs as ccrs
import xesmf
import timeit
from pathlib import Path 



from bopak.reprocess import set_BOEASE2, ptproj, piomas_db, piomas

PIOMAS_ROOT = Path('/data/PIOMAS')
outdir = Path('/data/BO/PIOMAS')
if not outdir.exists(): outdir.mkdir(parents=True, exist_ok=True)

dbdf = piomas_db(PIOMAS_ROOT)
fns = dbdf.fn.values
grid_fn = PIOMAS_ROOT/'grid.dat'

# --- setting up regridder and projection 
mds, boprj = set_BOEASE2()
pcproj = ccrs.PlateCarree()

# pdat = piomas(fns[0])
# lon,lat = piomas.read_grid(grid_fn)
# pdat.load_data(lon=lon,lat=lat)
# pdat.set_regridder(mds, reuse_weight=False)

reuse_weight = False
method = 'bilinear'
regridder_fn = Path('./.regridder/PIOMAS2EASE2.nc')

fcomp = dict(zlib=True, complevel=5, dtype='float32')
t_start = timeit.default_timer()
for ifn,fn in enumerate(fns):
    
    outfn = outdir/f'{fn.stem}.nc'
    pdat = piomas(fn)
    if ifn==0: lon,lat = piomas.read_grid(grid_fn)
    pdat.load_data(lon=lon,lat=lat)
    pdat.set_regridder(mds, reuse_weight=reuse_weight,method=method, regridder_fn=regridder_fn)
    if ifn==0: reuse_weight = True
    pdat.regrid()
    pdat.save_regrid(outfn=outfn)

    delta_time = timeit.default_timer() - t_start
    t_start = timeit.default_timer()
    print(f'Time for {ifn:2d} {fn.stem}: {delta_time}s.')