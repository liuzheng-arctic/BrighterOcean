import numpy as np
import xarray as xr 
import pandas as pd
from pathlib import Path

from bopak.reprocess import combined_data as cbd



BOROOT = Path('/data/BO/EASE2')


ftype = 'nc'
freq = 'Y'
# ftype = 'txt'
# freq = 'D'
outroot = BOROOT/'combined' 
SYEAR = 2011
EYEAR = 2011
# SYEAR = 2005
# EYEAR = 2005


albdir = BOROOT/'albedo/V0.06'
btdir = BOROOT/'BTSIC'
sicdir = BOROOT/'SIC'
onsetdir = BOROOT/'onset'
iadir = BOROOT/'ICEAGE/V4'
eradir = BOROOT/'ERA5/'
piodir = BOROOT/'PIOMAS'
ponddir = BOROOT/'pond'

for tyr in range(SYEAR, EYEAR+1):
    tdat = cbd(
        tyr,
        albdir=albdir,
        sicdir=sicdir,
        btdir=btdir,
        onsetdir=onsetdir,
        iadir=iadir, 
        eradir=eradir,
        piodir=piodir,
        ponddir=ponddir)
    tdat.process()
    tdat.to_file(outroot,ftype=ftype,freq=freq)