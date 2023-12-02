import numpy as np 
import xarray as xr 
import pandas as pd 
from pathlib import Path 
import timeit


fERAroot = Path('/data/BO/ERA5/hourly')
outroot = Path('/data/BO/ERA5/daily')
fns = sorted( fERAroot.glob('*.nc') )

fcomp = dict(zlib=True, complevel=5, dtype='float32')
t_start = timeit.default_timer()
for ifn,fn in enumerate(fns):
    
    outfn = outroot/fn.name
    with xr.open_dataset(fn) as ds:
        ds_day = ds.groupby('time.dayofyear').sum(dim='time')
    encoding = {var: fcomp for var in ds_day.data_vars}
    ds_day.to_netcdf(outfn,encoding=encoding)
    delta_time = timeit.default_timer() - t_start
    t_start = timeit.default_timer()
    print(f'Time for {ifn:2d} {fn.stem}: {delta_time}s.')
    


