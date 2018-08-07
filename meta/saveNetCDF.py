# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Save Output Files into NetCDF files
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.06
Last Update     : 2018.08.06
Contributor     :
Description     : This module aims to save the output fields into the netCDF files.
                  It outputs netCDF 4 files. The data has been compressed.
                  
Return Values   : netCDF files
"""

import sys
import os
import numpy as np
from netCDF4 import Dataset

class savenc:
    def __init__(self):
        """
        Save the output fields into netCDF files.
        """
        print ("Save output fields as netCDF4 files.")
        
    def ncCorrect(self, uc, vc, year, lat, lon, path):
        """
        Save the baratropic corrected wind into netCDF files.
        param uc:
        param vc:
        
        return: netCDF4 files containing uc and vc on the native grid.
        rtype: netCDF4
        """
        logging.info("Start creating netcdf file for baratropic corrected wind fields at each grid point.")
        data_wrap = Dataset(os.path.join(path, 'model_subdaily_075_{}_{}_uvc_point.nc'.format(year[0],year[-1])),
                            'w',format = 'NETCDF4')
        # create dimensions for netcdf data
        year_wrap_dim = data_wrap.createDimension('year',len(year))
        month_wrap_dim = data_wrap.createDimension('month',12)
        lat_wrap_dim = data_wrap.createDimension('latitude',len(lat))
        lon_wrap_dim = data_wrap.createDimension('longitude',len(lon))
        # create 1-dimension variables
        year_warp_var = data_wrap.createVariable('year',np.int32,('year',))
        month_warp_var = data_wrap.createVariable('month',np.int32,('month',))
        lat_warp_var = data_wrap.createVariable('latitude',np.float32,('latitude',))
        lon_warp_var = data_wrap.createVariable('longitude',np.float32,('longitude',))
        # create 4-dimension variables
        uc_warp_var = data_wrap.createVariable('uc',np.float32,('year','month','latitude','longitude'))
        vc_warp_var = data_wrap.createVariable('vc',np.float32,('year','month','latitude','longitude'))
        # global attributes
        data_wrap.description = 'Monthly mean baratropic corrected wind fields.'
        # variable attributes
        lat_warp_var.units = 'degree_north'
        lon_warp_var.units = 'degree_east'
        uc_warp_var.units = 'm/s'
        vc_warp_var.units = 'm/s'
        
        uc_warp_var.long_name = 'zonal barotropic correction wind'
        vc_warp_var.long_name = 'meridional barotropic correction wind'
        # writing data
        year_warp_var[:] = year
        month_warp_var[:] = np.arange(1,13,1)
        lat_warp_var[:] = lat
        lon_warp_var[:] = lon
        uc_warp_var[:] = uc
        vc_warp_var[:] = vc
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")
        
        
        