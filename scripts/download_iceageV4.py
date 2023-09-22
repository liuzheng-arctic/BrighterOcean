import os
from pathlib import Path


SYEAR = 1985
EYEAR = 2022
outroot = Path('/data/ICEAGE/V4')
urlroot = 'https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0611_seaice_age_v4/data/'
if not outroot.exists(): outroot.mkdir(parents=True,exist_ok=True)

wget_pre = f'wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off -nd --directory-prefix={outroot}'

for tyr in range(SYEAR,EYEAR+1):
    fn = f'iceage_nh_12.5km_{tyr}0101_{tyr}1231_v4.1.nc'
    url = urlroot+fn
    os.system(f'{wget_pre} {url}')
