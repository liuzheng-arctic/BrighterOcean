import requests
import shutil
from pathlib import Path


SYEAR = 1979
EYEAR = 2022
outroot = Path('/data/PIOMAS/hiday')
urlroot = 'https://pscfiles.apl.washington.edu/zhang/PIOMAS/data/v2.1/hiday/'
if not outroot.exists(): outroot.mkdir(parents=True,exist_ok=True)

for tyr in range(SYEAR,EYEAR+1):
    fn = f'hiday.H{tyr:4d}.gz'
    url = urlroot+fn
    print(url)
    r = requests.get(url, auth=('anonymous', 'anonymous'), verify=False,stream=True)
    with open(outroot/fn, 'wb') as fid:
        shutil.copyfileobj(r.raw, fid)
