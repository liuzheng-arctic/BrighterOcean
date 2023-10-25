import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from cartopy.feature import OCEAN 
from cartopy.feature import LAND

from .database import iceage_db, onset_db, albedo_db, sic_db_ease2
from .utils import set_BOEASE2, loc_ease2
from .sic import sic
from .onset import onset
from .iceage import iceage

class test_year:

    def __init__(self,
                 year,
                 ia_root = None,
                 onset_root = None,
                 alb_root = None,
                 sic_root = None,
                 ):
        '''
        '''
        self.year = year
        self.ia_root = ia_root
        self.onset_root = onset_root
        self.alb_root = alb_root
        self.sic_root = sic_root
        # if ia_root is not None: self.load_iceage(ia_root)
        # if onset_root is not None: self.load_onset(onset_root)
        # if alb_root is not None: self.load_albedo(alb_root)
        return
    

    def load_data(self,):
        
        '''
        '''
        tyr = self.year
        onset_df = onset_db(self.onset_root)
        ia_df = iceage_db(self.ia_root)
        alb_df = albedo_db(self.alb_root)

        if self.sic_root is not None:
            sic_df = sic_db_ease2(self.sic_root)
            fn_sic = sic_df[sic_df.time.dt.year==tyr].iloc[0].fn
            with xr.open_dataset(fn_sic) as ds:
                ds.load()
            self.sic = ds

        mds, boprj = set_BOEASE2()
        self.mds = mds
        self.boprj = boprj

        fn_onset = onset_df[onset_df.time.dt.year==tyr].iloc[0].fn
        fn_ia = ia_df[ia_df.time.dt.year==tyr].iloc[0].fn
        fn_alb = alb_df[alb_df.year==tyr].iloc[0].fn
        print(fn_onset,fn_ia,fn_alb)

        with xr.open_dataset(fn_onset) as ds:
            ds.load()
        self.onset = ds

        with xr.open_dataset(fn_ia) as ds:
            ds.load()
        self.iceage = ds

        with xr.open_dataset(fn_alb) as ds:
            ds.load()
        self.alb = ds

        with xr.open_dataset(fn_alb) as ds:
            ds.load()
        self.alb = ds


        return
    

    def plot_onset(self,ftsz=16):
        '''
        '''
        tyr = self.year
        mds = self.mds
        boprj = self.boprj
        XLIM = mds.XLIM
        YLIM = mds.YLIM
        
        xx,yy = np.meshgrid(mds.X,mds.Y)
        vns = ['Earlymelt','Melt','Earlyfreeze','Freeze']
        mlvl = levels=np.linspace(80,210,27)
        flvl = levels=np.linspace(210,360,31)
        fig, axs = plt.subplots(2,2,figsize=(11,10),subplot_kw=dict(projection=self.boprj),sharex=True,sharey=True)
        chdls = []
        for iplt in range(4):
            irow = iplt//2
            icol = iplt%2
            vn = vns[iplt]
            ax = axs[irow,icol]
            ax.set_extent([-XLIM,XLIM,-YLIM,YLIM],crs=boprj)
            ax.coastlines(color='.5')
            ax.gridlines()
            lvls = mlvl if iplt<2 else flvl
            chdl = ax.contourf(xx,yy,self.onset[vn],cmap='RdBu_r',levels=lvls,extend='both')
            chdls.append(chdl)
            ax.add_feature(OCEAN,color='c',alpha=.5)
            ax.add_feature(LAND,color='.6')
            ax.annotate(vn,(-XLIM*.9,YLIM*.85),fontsize=ftsz)
        fig.subplots_adjust(wspace=0.01,hspace=.05,top=.95,bottom=.05,left=.05,right=.90)
        fig.suptitle(f'Year: {tyr}')
        cax = fig.add_axes([.92,.52,.03,.43])
        fig.colorbar(chdls[0],ax=ax,cax=cax)
        cax = fig.add_axes([.92,.05,.03,.43])
        fig.colorbar(chdls[-1],ax=ax,cax=cax)



        return fig, axs
    

    def plot_accordion(
            self,lat,lon,
            ftsz=16,
            lw=2, 
            fig = None, 
            ax = None,
            add_xlabel=True,
            xlim=None,
            #linestyle='-',
            marker=None,
            ms=None,
            add_sic=False,
            ):
        '''
        '''

        tyr = self.year
        iy,ix = loc_ease2(lon,lat)

        tt = self.alb.time
        alb0 = self.alb.isel(X=ix,Y=iy)
        onset0 = self.onset.isel(X=ix,Y=iy)

        # --- expand iceage to daily mask for multi-year
        ia0 = self.iceage.isel(X=ix,Y=iy)['age_of_sea_ice'].values
        ia7 = np.repeat(ia0,7)
        nn = 2 if tt.dt.is_leap_year.any() else 1
        for i in range(nn):
            ia7 = np.append(ia7,ia0[-1])
        ia_mask = xr.where(ia7>1,1,0)

        fy = alb0['first_year']
        my = alb0['multi_year']
        aa = alb0['actual_age']

        em_jday = onset0['Earlymelt'].values.astype(int)
        md_jday = onset0['Melt'].values.astype(int)
        fr_jday = onset0['Freeze'].values.astype(int)
        ef_jday = onset0['Earlyfreeze'].values.astype(int)

        assert ~np.isnan(onset0['Earlymelt'].values) and \
            ~np.isnan(onset0['Melt'].values) and \
            ~np.isnan(onset0['Earlyfreeze'].values) and\
            ~np.isnan(onset0['Melt'].values), f'No onset date for this location: lat={lat}, lon={lon}.'
            

        em = pd.to_datetime(em_jday+tyr*1000,format='%Y%j')
        md = pd.to_datetime(md_jday+tyr*1000,format='%Y%j')
        fr = pd.to_datetime(fr_jday+tyr*1000,format='%Y%j')
        ef = pd.to_datetime(ef_jday+tyr*1000,format='%Y%j')

        em_str = em.strftime('%Y-%m-%d')
        md_str = md.strftime('%Y-%m-%d')
        fr_str = fr.strftime('%Y-%m-%d')
        ef_str = ef.strftime('%Y-%m-%d')

        if fig is None or ax is None: 
            fig, ax = plt.subplots(1,1,figsize=(10,5))

        ax.plot(tt,fy,'g',tt,my,'k',tt,aa,'b',lw=2,marker=marker,ms=ms)
        ax.legend(['first-year','multi-year','actual-age'],loc=3,fontsize=ftsz)
        ax.axvline(x=em,ymax=.9,linestyle=':',color='r',alpha=.7)
        ax.axvline(x=md,ymax=.9,linestyle='-.',color='r',alpha=.7)
        ax.axvline(x=ef,ymax=.3,linestyle=':',color='k',alpha=.7)
        ax.axvline(x=fr,ymax=.3,linestyle='-.',color='k',alpha=.7)
        ax.fill_between(tt,ia_mask,color='.5',alpha=.3)

        x0 = pd.to_datetime(str(tyr)+'0101')
        if xlim is not None:
            xlim_pd = pd.to_datetime(xlim)
            ax.set_xlim(xlim_pd)
            x0 = xlim_pd[0]

        #ax.set_ylim([0,.9])
        ax.set_ylim([0,1.01])
        if self.sic_root is not None and add_sic:
            sic0 = self.sic['cdr_seaice_conc'].isel(X=ix,Y=iy)
            ax.plot(tt,sic0,'m-',lw=lw)
        ax.grid(True)
        ax.annotate(f'Early Melt: {em_str}: {em_jday}',(x0,.8))
        ax.annotate(f'Melt: {md_str}: {md_jday}',(x0,.7))
        ax.annotate(f'Early freeze: {ef_str}: {ef_jday}',(x0,.6))
        ax.annotate(f'Freeze: {fr_str}: {fr_jday}',(x0,.5))
        ax.annotate(f'Gray shading: ice age>1',(x0,.4))
        ax.set_title(f'{tyr}: lat={lat}, lon={lon}',fontsize=ftsz)
        ax.set_ylabel('albedo',fontsize=ftsz-2)
        if add_xlabel: ax.set_xlabel('date',fontsize=ftsz-2)
        return fig, ax


    def compare_accordion(self, loc0, loc1,
            xlim=None,
            #linestyle='-',
            marker=None,
            ms=None,):
        '''
        '''
        fig, axs = plt.subplots(2,1,figsize=(10,11),sharex=False)
        
        fig, ax = self.plot_accordion(*loc0, fig=fig, ax=axs[0],add_xlabel=False,xlim=xlim,marker=marker,ms=ms)
        fig, ax = self.plot_accordion(*loc1, fig=fig, ax=axs[1],xlim=xlim,marker=marker,ms=ms)
        fig.subplots_adjust()
        return fig, axs

    
 