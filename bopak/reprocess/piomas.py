import numpy as np
import gzip
import struct
import xesmf

PIOMAS_DIMS = (120,360)

class piomas:
    def __init__(self,fn,grid_fn = None):
        '''
        Parameters:
        ----------
        fn: str, full path to the PIOMAS gzipped data
        grid_fn: str, full path to the PIOMAS grid data
        '''
        self.fn = fn
        self.dims = PIOMAS_DIMS
        # if grid_fn is not None:
        #     self.grid_fn = grid_fn
        #     lon, lat = self.read_grid(grid_fn)
        #     self.add_grid(lon,lat)
        self.data = self.read_hiday(self.fn)
        return
    
    #def regrid(self,)
    
    @staticmethod
    def read_hiday(fn):
        '''
        Read gzipped PIOMAS daily sea ice thickness data.

        Parameters:
        ----------
        fn: str, full path to the PIOMAS gzipped data

        Return:
        ------
        data3d: np.ndarray, sea ice thickness data with dimension of [day,Y,X]

        '''
        
        with gzip.open(fn, mode='rb') as file:

            fileContent = file.read()
            data = struct.unpack("f" * (len(fileContent)// 4), fileContent)
            ny,nx = PIOMAS_DIMS
            nday = int(len(data)/ny/nx)
            data3d = np.array(data).reshape((nday,ny,nx))
        return data3d
    

    # def add_grid(lon, lat):
    #     '''
    #     add lon, lat of grid to instance

    #     Parameters:
    #     ----------
    #     lon: np.ndarray
    #     '''
    #     self.lon = lon
    #     self.lat = lat
    #     return 
    
    @staticmethod
    def read_grid(grid_fn):
        '''
        Read PIOMAS grid data from file

        '''
        grid_data = np.loadtxt(grid_fn)
        lon, lat = grid_data.ravel().reshape((2,)+PIOMAS_DIMS)
        return lon, lat
