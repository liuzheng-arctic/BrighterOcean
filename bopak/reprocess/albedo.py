import numpy as np
import pandas as pd
import xarray as xr

from pathlib import Path

from .utils import set_BOEASE2

SNOW_ALB = 0.85
OCN_ALB = 0.07
POND_ALB = 0.46
MELT_MIN = 75
MELT_MAX = 210
FREEZE_MIN = 210
FREEZE_MAX = 410
SIC_MIN = 0.15
T_ALL_STG_FIRST = 21 + np.ceil( (0.54-0.2)/.0083 )
T_ALL_STG_MULTI = 21 + 1 + 1 + np.ceil( (0.56-0.2)/.0029 )
T_FREEZE_MULTI = np.ceil( (SNOW_ALB-POND_ALB)/0.026 )

class albedo:

    def __init__(
            self,
            year,
            onset_date,  
            alb_init = None,
            mds = None,
            boprj = None,
    ):
        '''
        '''
        self.year = year
        self.onset_date = onset_date
        self.alb_init = alb_init
        self.set_grid(mds=mds,boprj=boprj)
        return 
    
    def set_grid(self,mds=None,boprj=None):
        '''
        '''
        if mds is  None or boprj is None:
            mds, boprj = set_BOEASE2()
        tt = pd.date_range(self.year+'0101',self.year+'1231',freq='D')
        ds = xr.Dataset()
        ds['X'] = mds.X 
        ds['Y'] = mds.Y
        ds['time'] = tt
        ds = ds.assign_coords(dict(lat=mds.lat,lon=mds.lon))

        ds['first_year'] = xr.DataArray(np.nan,coords=dict(X=ds.X,Y=ds.Y,time=ds.time),dims=('time','Y','X'))
        ds['multi_year'] = xr.DataArray(np.nan,coords=dict(X=ds.X,Y=ds.Y,time=ds.time),dims=('time','Y','X'))
        ds['actual_age'] = xr.DataArray(np.nan,coords=dict(X=ds.X,Y=ds.Y,time=ds.time),dims=('time','Y','X'))

        ds_fnl = xr.Dataset()
        ds_fnl = ds_fnl.assign_coords(dict(lat=mds.lat,lon=mds.lon))
        for vn in ds.data_vars:
            if self.alb_init is not None:
                ds[vn][0] = self.alb_init[vn]
            ds_fnl[vn] = xr.DataArray(np.nan,coords=dict(X=ds.X,Y=ds.Y),dims=('Y','X'))

        self.data = ds
        # --- now, use alb_fnl to track albedo on day efreeze-1 for multi-year and on freeze-1 for first-year
        self.alb_fnl = ds_fnl
        return
    
    def first_year(self, 
                  t_stamp, 
                  melt_offset = 0,
                  freeze_offset = 0,
                  format='%Y-%m-%d'):
        '''
        '''
        t0 = pd.to_datetime(t_stamp,format=format)
        doy = t0.day_of_year
        it = doy - 1
        tdat = self.data['first_year'][it]
        #return

        melt = self.onset_date['Melt'] + melt_offset
        freeze = self.onset_date['Freeze'] + freeze_offset
        # melt = self.onset_date['Earlymelt']
        # freeze = self.onset_date['Earlyfreeze']

        valid_melt = (melt>=MELT_MIN) & (melt<=MELT_MAX)
        valid_freeze = (freeze>=FREEZE_MIN) & (melt<=FREEZE_MAX)
        valid_onset = valid_melt & valid_freeze

        #stg0 = (doy<melt) & (valid_melt)
        stg0 = (doy<melt) & (valid_onset)
        stg1 = (doy>=melt) & (doy<melt+7) & (doy<freeze) & valid_onset
        stg2 = (doy>=melt+7) & (doy<melt+14) & (doy<freeze) & valid_onset
        stg3 = (doy>=melt+14) & (doy<melt+21) & (doy<freeze) & valid_onset
        stg4 = (doy>=melt+21) & (doy<melt+T_ALL_STG_FIRST) & (doy<freeze) & valid_onset
        stg4_end = (doy==melt+T_ALL_STG_FIRST) & (doy<freeze) & valid_onset
        stg5 = (doy>=melt+T_ALL_STG_FIRST+1) & (doy<freeze) & valid_onset
        stg6 = (doy>=freeze) & valid_onset

        
        stg_em = (doy==freeze-1) & valid_onset

        R1 = (0.6-SNOW_ALB)/7.
        R2 = (0.32-.6)/7.
        R3 = (0.54-0.32)/7.
        R4 = -0.0083
        R6 = 0.0082

        #return tdat
        tdat = xr.where( stg0, SNOW_ALB, tdat )
        #return tdat
        # --- Stage 1: 0.85 --> 0.6
        tdat = xr.where( stg1, SNOW_ALB + (doy-melt)*R1, tdat )
        # --- Stage 2: 0.6 --> 0.32
        tdat = xr.where( stg2, 0.6 + (doy-melt-7)*R2, tdat )
        # --- Stage 3: 0.32 --> 0.54
        tdat = xr.where( stg3, 0.32 + (doy-melt-14)*R3, tdat )
        # --- Stage 4: 0.054 --> 0.1997, which is set to 0.2
        tdat = xr.where( stg4, 0.54 + (doy-melt-21)*R4, tdat )
        tdat = xr.where( stg4_end, 0.2, tdat )
        # --- Stage 5: 0.2 to 0.07
        tdat = xr.where( stg5, 0.07, tdat )

        
        tdat_em = self.alb_fnl['first_year']
        tdat_em = xr.where( stg_em, tdat, tdat_em )
        self.alb_fnl['first_year'] = tdat_em

        # --- Stage 6: freezing, ramp up to 0.85
        #tdat = xr.where( stg6, OCN_ALB + (doy-freeze)*R6 , tdat )
        tdat = xr.where( stg6 & (tdat_em<=OCN_ALB), OCN_ALB + (doy-freeze)*R6 , tdat )
        tdat = xr.where( stg6 & (tdat_em>OCN_ALB), tdat_em + (doy-freeze)*R6 , tdat )
        # tdat = xr.where( stg6 & (freeze>=melt+T_ALL_STG_FIRST+1), OCN_ALB + (doy-freeze)*R6 , tdat )
        # tdat = xr.where( stg6 & (freeze<melt+T_ALL_STG_FIRST+1) , 0.54 + (freeze-melt-21)*R4 + (doy-freeze)*R6 , tdat )
        
        tdat = xr.where( stg6 & (tdat>SNOW_ALB) , SNOW_ALB, tdat )

        # --- assign tdat back to data
        self.data['first_year'][it] = tdat
        #self.data['first_year'][it] = xr.where( sic>=SIC_MIN, tdat, OCN_ALB )

        if doy==len(self.data.time):
            self.alb_fnl['first_year'] = self.data['first_year'][it]

        return #tdat
     
    def multi_year(self, 
                  t_stamp, 
                  melt_offset = 0,
                  freeze_offset = 0,
                  format='%Y-%m-%d'):
        '''
        '''
        t0 = pd.to_datetime(t_stamp,format=format)
        doy = t0.day_of_year
        it = doy - 1
        tdat = self.data['multi_year'][it]
        # tdat0 = self.data['multi_year'][0].values*1.
        # if it>1:
        #     tdat0 = self.data['multi_year'][it-1].values*1.

        melt = self.onset_date['Melt'] + melt_offset
        freeze = self.onset_date['Freeze'] + freeze_offset

        valid_melt = (melt>=MELT_MIN) & (melt<=MELT_MAX)
        valid_freeze = (freeze>=FREEZE_MIN) & (melt<=FREEZE_MAX)
        valid_onset = valid_melt & valid_freeze 

        #stg0 = (doy<emelt) & valid_onset
        #stg1 = (doy>=emelt) &(doy<melt) & valid_onset
        stg1 = (doy<melt) & valid_onset
        stg2 = (doy>=melt) & (doy<=melt+15) & (doy<freeze) & valid_onset
        stg3 = (doy>melt+15) & (doy<=melt+22) & (doy<freeze) & valid_onset
        # stg4 = (doy>melt+22) & (doy<melt+T_ALL_STG_MULTI) & (doy<efreeze) & valid_onset
        # stg4_end = (doy>=melt+T_ALL_STG_MULTI) & (doy<efreeze) & valid_onset
        #stg4 = (doy>melt+22) & (doy<melt+T_ALL_STG_MULTI) & (doy<=efreeze)&  valid_onset
        stg4 = (doy>melt+22) & (doy<melt+T_ALL_STG_MULTI) & (doy<freeze)&  valid_onset
        stg4_end = (doy>=melt+T_ALL_STG_MULTI) & (doy<freeze)&  valid_onset
        #stg5 = (doy>=efreeze) & (doy<freeze) & valid_onset
        stg6 = (doy>=freeze) & valid_onset

        stg_em = (doy==freeze-1) & valid_onset

        R2 = (0.71-SNOW_ALB)/15.
        R3 = (0.5-0.70)/6
        R4 = -0.0029
        R6 = 0.026

        #tdat = xr.where( stg0, SNOW_ALB, tdat )
        # --- Stage 1
        tdat = xr.where( stg1, SNOW_ALB, tdat )
        # --- Stage 2
        tdat = xr.where( stg2, SNOW_ALB + (doy-melt)*R2, tdat )
        # --- Stage 3
        tdat = xr.where( stg3, 0.70 + (doy-melt-15-1)*R3, tdat )
        #tdat = xr.where( stg3_end, 0.56, tdat )
        # --- Stage 4: 
        tdat = xr.where( stg4, 0.56 + (doy-melt-22-1)*R4, tdat )
        tdat = xr.where( stg4&(tdat<.2) , 0.2, tdat   )
        # --- Stage 4_end: 
        tdat = xr.where( stg4_end, 0.2, tdat  )

        tdat_em = self.alb_fnl['multi_year']
        tdat_em = xr.where( stg_em, tdat, tdat_em )
        self.alb_fnl['multi_year'] = tdat_em
        # --- Stage 5: 
        # tdat = xr.where( stg5&((tdat_em>=0.46)|tdat.isnull()), tdat_em, tdat )
        # tdat = xr.where( stg5&((tdat_em<0.46)|tdat.isnull()), 0.46, tdat )

        # --- Stage 6: 
        tdat = xr.where( stg6&(tdat_em<=.46) , 0.46 + (doy-freeze)*R6 , tdat )   
        tdat = xr.where( stg6&(tdat_em>.46) , tdat_em + (doy-freeze)*R6 , tdat )          
        tdat = xr.where( tdat>SNOW_ALB , SNOW_ALB, tdat )

        # --- assign tdat back to data
        self.data['multi_year'][it] = tdat

        if doy==len(self.data.time):
            self.alb_fnl['multi_year'] = self.data['multi_year'][it]
            
        return #tdat   


    
    def multi_year_all_onset(self, 
                  t_stamp, 
                  format='%Y-%m-%d'):
        '''
        '''
        t0 = pd.to_datetime(t_stamp,format=format)
        doy = t0.day_of_year
        it = doy - 1
        tdat = self.data['multi_year'][it]
        # tdat0 = self.data['multi_year'][0].values*1.
        # if it>1:
        #     tdat0 = self.data['multi_year'][it-1].values*1.

        emelt = self.onset_date['Earlymelt']
        melt = self.onset_date['Melt']
        freeze = self.onset_date['Freeze']
        efreeze = self.onset_date['Earlyfreeze']

        valid_melt = (melt>=MELT_MIN) & (melt<=MELT_MAX)
        valid_freeze = (freeze>=FREEZE_MIN) & (melt<=FREEZE_MAX)
        valid_emelt = (emelt>=MELT_MIN) & (emelt<=MELT_MAX)
        valid_efreeze = (efreeze>=FREEZE_MIN) & (emelt<=FREEZE_MAX)
        valid_onset = valid_melt & valid_emelt & valid_freeze & valid_efreeze

        stg0 = (doy<emelt) & valid_onset
        stg1 = (doy>=emelt) &(doy<melt) & valid_onset
        stg2 = (doy>=melt) & (doy<=melt+15) & (doy<efreeze) & valid_onset
        stg3 = (doy>melt+15) & (doy<=melt+22) & (doy<efreeze) & valid_onset
        # stg4 = (doy>melt+22) & (doy<melt+T_ALL_STG_MULTI) & (doy<efreeze) & valid_onset
        # stg4_end = (doy>=melt+T_ALL_STG_MULTI) & (doy<efreeze) & valid_onset
        stg4 = (doy>melt+22) & (doy<melt+T_ALL_STG_MULTI) & (doy<=efreeze)&  valid_onset
        stg5 = (doy>=efreeze) & (doy<freeze) & valid_onset
        stg6 = (doy>=freeze) & valid_onset

        stg_em = (doy==efreeze-1) & valid_onset

        R2 = (0.71-0.81)/15.
        R3 = (0.5-0.70)/6
        R4 = -0.0029
        R6 = 0.026

        tdat = xr.where( stg0, SNOW_ALB, tdat )
        # --- Stage 1
        tdat = xr.where( stg1, 0.81, tdat )
        # --- Stage 2
        tdat = xr.where( stg2, 0.81 + (doy-melt)*R2, tdat )
        # --- Stage 3
        tdat = xr.where( stg3, 0.70 + (doy-melt-15-1)*R3, tdat )
        #tdat = xr.where( stg3_end, 0.56, tdat )
        # --- Stage 4: 
        tdat = xr.where( stg4, 0.56 + (doy-melt-22-1)*R4, tdat )
        tdat = xr.where( stg4&(tdat<.2) , 0.2, tdat   )

        tdat_em = self.alb_fnl['multi_year']
        tdat_em = xr.where( stg_em, tdat, tdat_em )
        self.alb_fnl['multi_year'] = tdat_em
        # --- Stage 5: 
        tdat = xr.where( stg5&((tdat_em>=0.46)|tdat.isnull()), tdat_em, tdat )
        tdat = xr.where( stg5&((tdat_em<0.46)|tdat.isnull()), 0.46, tdat )

        #tdat = xr.where( stg5&(efreeze>=melt+T_ALL_STG_MULTI), 0.46, tdat )
        #tdat = xr.where( stg5&(efreeze<melt+T_ALL_STG_MULTI), 0.56 + (efreeze-melt-22-1)*R4, tdat )
        #tdat = tdat( stg5&(tdat>=0.46), 0.46 )
        # --- Stage 6: 
        tdat = xr.where( stg6&(tdat_em<=.46) , 0.46 + (doy-freeze)*R6 , tdat )   
        tdat = xr.where( stg6&(tdat_em>.46) , tdat_em + (doy-freeze)*R6 , tdat )    
        # tdat = xr.where( stg6 & (freeze>=melt+T_ALL_STG_MULTI), 0.46 + (doy-freeze)*R6 , tdat )   
        # tdat = xr.where( stg6 & (freeze<melt+T_ALL_STG_MULTI), 0.56 + (efreeze-melt-22-1)*R4 + (doy-freeze)*R6 , tdat )        
        tdat = xr.where( tdat>SNOW_ALB , SNOW_ALB, tdat )

        # --- assign tdat back to data
        self.data['multi_year'][it] = tdat

        if doy==len(self.data.time):
            self.alb_fnl['multi_year'] = self.data['multi_year'][it]
            
        return #tdat
    

    def actual_age(self,
                   t_stamp,
                   ice_age,
                   format='%Y-%m-%d'):
        '''
        '''
        t0 = pd.to_datetime(t_stamp,format=format)
        doy = t0.day_of_year
        it = doy - 1
        
        
        fy_it = self.data['first_year'][it]
        my_it = self.data['multi_year'][it]

        ia_it = xr.where( ice_age<=1 , fy_it, my_it )
        # --- filter by ice age data. Remove pixels without ice age too. 
        ia_it = ia_it.where( (ice_age>=0)&(ice_age<20) )
        self.data['actual_age'][it] = ia_it

        if doy==len(self.data.time):
            self.alb_fnl['actual_age'] = self.data['actual_age'][it]

        return #ia_it
    

    @staticmethod
    def load_init(fn):
        '''
        '''
        with xr.open_dataset(fn) as ds:
            ds.load()
        ds_init = xr.Dataset()
        for vn in ds.data_vars:
            ds_init[vn] = ds[vn][-1]
        return ds_init