# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Save Output Files into NetCDF files
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.09
Last Update     : 2018.08.15
Contributor     :
Description     : This module aims to integrate all the single output fields into 
                  a single netCDF file. Simple postprocessing, such as zonal integral,
                  ocean basin division, will be carried out. It outputs netCDF 4 files.
                  The data has been compressed.
                  
Return Values   : netCDF files
"""

import sys
import os
import numpy as np
import logging
from netCDF4 import Dataset

class assembly:
    def __init__(self, year_start, year_end, in_path, out_path):
        """
        Assemble all the output netCDF files into a single netCDF file.
        param year_start: the starting time for the given files
        param year_end: the ending time for the given files
        param in_path: path of the input files
        param out_path: path of the output files
        
        """
        print ("Assemble all the output netCDF files into a single netCDF file.")
        self.year_start = year_start
        self.year_end = year_end
        self.in_path = in_path
        self.out_path = out_path
        
        
    def ncAMET(self, name='ERAI'):
        """
        Put all the single netCDF files into one file for the sake of postprocessing.
        """
        logging.info("Start wrap-up all the netcdf files of AMET and calculate the zonal mean of each variable.")
        if name == 'ERAI':
            alias = 'era'
        elif name == 'MERRA2':
            alias = 'merra'
        elif name == 'JRA55':
            alias = 'jra'
        else:
            raise IOError("This dataset is not supported in this module.")
        # create time dimensions
        year = np.arange(self.year_start, self.year_end+1, 1)
        month = np.arange(1, 13, 1)
        # load dimensions from input files

    def ncOMET(self, tmaskpac, tmaskatl, name='ORAS4'):
        """
        Put all the single netCDF files into one file for the sake of postprocessing.
        The 3D and 4D variables are compressed.
        param tmaskpac: land-sea mask for the Pacific Ocean  [jj, ji]
        param tmaskatl: land-sea mask for the Atlantic Ocean  [jj, ji]
        """
        logging.info("Start wrap-up all the netcdf files of OMET and calculate the zonal mean of each variable.")
        if name == 'ORAS4':
            alias = 'oras'
        elif name == 'GLOTYS2V3':
            alias = 'glorys'
        elif name == 'SODA3':
            alias = 'soda'
        else:
            raise IOError("This dataset is not supported in this module.")           
        # create time dimensions
        year = np.arange(self.year_start, self.year_end+1, 1)
        month = np.arange(1, 13, 1)
        # load dimensions from input files
        data_example = Dataset(os.path.join(self.in_path,
                               '{}_model_monthly_{}_omet_point.nc'.format(alias,self.year_start)))
        depth = data_example.variables['depth'][:]
        latitude_aux = data_example.variables['latitude_aux'][:]
        gphiv = data_example.variables['gphiv'][:]
        glamv = data_example.variables['glamv'][:]
        # dimension size
        t = len(year)
        z = len(depth)
        jj, ji = gphiv.shape
        # create arrays to store all the data and postprocess
        E = np.zeros((t, len(month), jj, ji), dtype=float)
        E_100 = np.zeros((t, len(month), jj, ji), dtype=float)
        E_300 = np.zeros((t, len(month), jj, ji), dtype=float)
        E_700 = np.zeros((t, len(month), jj, ji), dtype=float)
        E_2000 = np.zeros((t, len(month), jj, ji), dtype=float)
        E_pac_vert = np.zeros((t, len(month), z, jj), dtype=float)
        E_atl_vert = np.zeros((t, len(month), z, jj), dtype=float)
        E_vert = np.zeros((t, len(month), z, jj), dtype=float)
        
        E_int = np.zeros((t, len(month), jj), dtype=float)
        E_100_int = np.zeros((t, len(month), jj), dtype=float)
        E_300_int = np.zeros((t, len(month), jj), dtype=float)
        E_700_int = np.zeros((t, len(month), jj), dtype=float)
        E_2000_int = np.zeros((t, len(month), jj), dtype=float)    
        
        E_pac_int = np.zeros((t, len(month), jj), dtype=float)
        E_pac_100_int = np.zeros((t, len(month), jj), dtype=float)
        E_pac_300_int = np.zeros((t, len(month), jj), dtype=float)
        E_pac_700_int = np.zeros((t, len(month), jj), dtype=float)
        E_pac_2000_int = np.zeros((t, len(month), jj), dtype=float)

        E_atl_int = np.zeros((t, len(month), jj), dtype=float)
        E_atl_100_int = np.zeros((t, len(month), jj), dtype=float)
        E_atl_300_int = np.zeros((t, len(month), jj), dtype=float)
        E_atl_700_int = np.zeros((t, len(month), jj), dtype=float)
        E_atl_2000_int = np.zeros((t, len(month), jj), dtype=float)
        # increase the dimension of sub-basin masks for the sake of efficiency
        tmaskpac_3D = np.repeat(tmaskpac[np.newaxis,:,:],len(month),0)
        tmaskatl_3D = np.repeat(tmaskatl[np.newaxis,:,:],len(month),0)
        for i in year:
            data_key = Dataset(os.path.join(self.out_path,
                               '{}_model_monthly_{}_omet_point.nc'.format(alias,i)))
            E[i-self.year_start,:,:,:] = data_key.variables['E_total'][:]
            E_100[i-self.year_start,:,:,:] = data_key.variables['E_100'][:]
            E_300[i-self.year_start,:,:,:] = data_key.variables['E_300'][:]
            E_700[i-self.year_start,:,:,:] = data_key.variables['E_700'][:]
            E_2000[i-self.year_start,:,:,:] = data_key.variables['E_2000'][:]
            E_pac_vert[i-self.year_start,:,:,:] = data_key.variables['E_pac_vert'][:]
            E_atl_vert[i-self.year_start,:,:,:] = data_key.variables['E_atl_vert'][:]
            E_vert[i-self.year_start,:,:,:] = data_key.variables['E_vert'][:]
            
            E_int[i-self.year_start,:,:] = np.sum(E[i-self.year_start,:,:,:], 2)
            E_100_int[i-self.year_start,:,:] = np.sum(E_100[i-self.year_start,:,:,:], 2)
            E_300_int[i-self.year_start,:,:] = np.sum(E_300[i-self.year_start,:,:,:], 2)
            E_700_int[i-self.year_start,:,:] = np.sum(E_700[i-self.year_start,:,:,:], 2)
            E_2000_int[i-self.year_start,:,:] = np.sum(E_2000[i-self.year_start,:,:,:], 2)
            
            E_pac_int[i-self.year_start,:,:] = np.sum(E[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            E_pac_100_int[i-self.year_start,:,:] = np.sum(E_100[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            E_pac_300_int[i-self.year_start,:,:] = np.sum(E_300[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            E_pac_700_int[i-self.year_start,:,:] = np.sum(E_700[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            E_pac_2000_int[i-self.year_start,:,:] = np.sum(E_2000[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            
            E_atl_int[i-self.year_start,:,:] = np.sum(E[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            E_atl_100_int[i-self.year_start,:,:] = np.sum(E_100[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            E_atl_300_int[i-self.year_start,:,:] = np.sum(E_300[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            E_atl_700_int[i-self.year_start,:,:] = np.sum(E_700[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            E_atl_2000_int[i-self.year_start,:,:] = np.sum(E_2000[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
        # create netCDF files
        data_wrap = Dataset(os.path.join(self.out_path,
                            '{}_model_monthly_{}_{}_omet.nc'.format(alias, self.year_start, self.year_end)),
                            'w',format = 'NETCDF4')
        # create dimensions for netcdf data
        year_wrap_dim = data_wrap.createDimension('year',t)
        month_wrap_dim = data_wrap.createDimension('month',12)
        depth_wrap_dim = data_wrap.createDimension('depth',z)
        lat_wrap_dim = data_wrap.createDimension('jj',jj)
        lon_wrap_dim = data_wrap.createDimension('ji',ji)
        # create 1-dimension variables
        year_wrap_var = data_wrap.createVariable('year',np.int32,('year',))
        month_wrap_var = data_wrap.createVariable('month',np.int32,('month',))
        depth_wrap_var = data_wrap.createVariable('depth',np.float32,('depth',))
        lat_wrap_var = data_wrap.createVariable('latitude_aux',np.float32,('jj',))
        # create 2-dimension variables
        gphiv_wrap_var = data_wrap.createVariable('gphiv',np.float32,('jj','ji'))
        glamv_wrap_var = data_wrap.createVariable('glamv',np.float32,('jj','ji'))
        tmaskpac_wrap_var = data_wrap.createVariable('tmaskpac',np.int32,('jj','ji'))
        tmaskatl_wrap_var = data_wrap.createVariable('tmaskatl',np.int32,('jj','ji'))
        # create 3-dimension variables
        E_int_wrap_var = data_wrap.createVariable('E_int',np.float32,('year','month','jj'),zlib=True)
        E_100_int_wrap_var = data_wrap.createVariable('E_100_int',np.float32,('year','month','jj'),zlib=True)
        E_300_int_wrap_var = data_wrap.createVariable('E_300_int',np.float32,('year','month','jj'),zlib=True)
        E_700_int_wrap_var = data_wrap.createVariable('E_700_int',np.float32,('year','month','jj'),zlib=True)
        E_2000_int_wrap_var = data_wrap.createVariable('E_2000_int',np.float32,('year','month','jj'),zlib=True)

        E_pac_int_wrap_var = data_wrap.createVariable('E_pac_int',np.float32,('year','month','jj'),zlib=True)
        E_pac_100_int_wrap_var = data_wrap.createVariable('E_pac_100_int',np.float32,('year','month','jj'),zlib=True)
        E_pac_300_int_wrap_var = data_wrap.createVariable('E_pac_300_int',np.float32,('year','month','jj'),zlib=True)
        E_pac_700_int_wrap_var = data_wrap.createVariable('E_pac_700_int',np.float32,('year','month','jj'),zlib=True)
        E_pac_2000_int_wrap_var = data_wrap.createVariable('E_pac_2000_int',np.float32,('year','month','jj'),zlib=True)

        E_atl_int_wrap_var = data_wrap.createVariable('E_atl_int',np.float32,('year','month','jj'),zlib=True)
        E_atl_100_int_wrap_var = data_wrap.createVariable('E_atl_100_int',np.float32,('year','month','jj'),zlib=True)
        E_atl_300_int_wrap_var = data_wrap.createVariable('E_atl_300_int',np.float32,('year','month','jj'),zlib=True)
        E_atl_700_int_wrap_var = data_wrap.createVariable('E_atl_700_int',np.float32,('year','month','jj'),zlib=True)
        E_atl_2000_int_wrap_var = data_wrap.createVariable('E_atl_2000_int',np.float32,('year','month','jj'),zlib=True)    
        # create 4-dimension variables
        E_0_wrap_var = data_wrap.createVariable('E_total',np.float32,('year','month','jj','ji'),zlib=True)
        E_100_wrap_var = data_wrap.createVariable('E_100',np.float32,('year','month','jj','ji'),zlib=True)
        E_300_wrap_var = data_wrap.createVariable('E_300',np.float32,('year','month','jj','ji'),zlib=True)
        E_700_wrap_var = data_wrap.createVariable('E_700',np.float32,('year','month','jj','ji'),zlib=True)
        E_2000_wrap_var = data_wrap.createVariable('E_2000',np.float32,('year','month','jj','ji'),zlib=True)
        E_vert_wrap_var = data_wrap.createVariable('E_vert',np.float32,('year','month','depth','jj'),zlib=True)
        E_pac_vert_wrap_var = data_wrap.createVariable('E_pac_vert',np.float32,('year','month','depth','jj'),zlib=True)
        E_atl_vert_wrap_var = data_wrap.createVariable('E_atl_vert',np.float32,('year','month','depth','jj'),zlib=True)        
        # global attributes
        data_wrap.description = 'Monthly mean meridional energy transport fields.'
        # variable attributes
        E_0_wrap_var.units = 'Tera Watt'
        E_100_wrap_var.units = 'Tera Watt'
        E_300_wrap_var.units = 'Tera Watt'
        E_700_wrap_var.units = 'Tera Watt'
        E_2000_wrap_var.units = 'Tera Watt'
        E_vert_wrap_var.units = 'Tera Watt'
        E_pac_vert_wrap_var.units = 'Tera Watt'
        E_atl_vert_wrap_var.units = 'Tera Watt'
        
        E_int_wrap_var.units = 'Tera Watt'
        E_100_int_wrap_var.units = 'Tera Watt'
        E_300_int_wrap_var.units = 'Tera Watt'
        E_700_int_wrap_var.units = 'Tera Watt'
        E_2000_int_wrap_var.units = 'Tera Watt'
        
        E_pac_int_wrap_var.units = 'Tera Watt'
        E_pac_100_int_wrap_var.units = 'Tera Watt'
        E_pac_300_int_wrap_var.units = 'Tera Watt'
        E_pac_700_int_wrap_var.units = 'Tera Watt'
        E_pac_2000_int_wrap_var.units = 'Tera Watt'

        E_atl_int_wrap_var.units = 'Tera Watt'
        E_atl_100_int_wrap_var.units = 'Tera Watt'
        E_atl_300_int_wrap_var.units = 'Tera Watt'
        E_atl_700_int_wrap_var.units = 'Tera Watt'
        E_atl_2000_int_wrap_var.units = 'Tera Watt'
        
        tmaskpac_wrap_var.long_name = 'land sea mask for the Pacific basin'
        tmaskatl_wrap_var.long_name = 'land sea mask for the Atlantic basin'
        
        E_0_wrap_var.long_name = 'total meridional energy transport over the entire column'
        E_100_wrap_var.long_name = 'total meridional energy transport from surface to 100m'
        E_300_wrap_var.long_name = 'total meridional energy transport from surface to 300m'
        E_700_wrap_var.long_name = 'total meridional energy transport from surface to 700m'
        E_2000_wrap_var.long_name = 'total meridional energy transport from surface to 2000m'
        E_vert_wrap_var.long_name = 'zonal integral of meridional energy transport for the globe'
        E_pac_vert_wrap_var.long_name = 'zonal integral of meridional energy transport in the Pacific Ocean'
        E_pac_vert_wrap_var.long_name = 'zonal integral of meridional energy transport in the Pacific Ocean'
        
        E_int_wrap_var.long_name = 'zonal integral of total meridional energy transport over the entire column'
        E_100_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 100m'
        E_300_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 300m'
        E_700_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 700m'
        E_2000_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 2000m'
        
        E_pac_int_wrap_var.long_name = 'zonal integral of total meridional energy transport over the entire column in the Pacific'
        E_pac_100_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 100m in the Pacific'
        E_pac_300_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 300m in the Pacific'
        E_pac_700_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 700m in the Pacific'
        E_pac_2000_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 2000m in the Pacific'
        
        E_atl_int_wrap_var.long_name = 'zonal integral of total meridional energy transport over the entire column in the Atlantic'
        E_atl_100_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 100m in the Atlantic'
        E_atl_300_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 300m in the Atlantic'
        E_atl_700_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 700m in the Atlantic'
        E_atl_2000_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 2000m in the Atlantic'
        # writing data
        year_wrap_var[:] = year
        month_wrap_var[:] = np.arange(1,13,1)
        depth_wrap_var[:] = depth
        lat_wrap_var[:] = latitude_aux
        gphiv_wrap_var[:] = gphiv
        glamv_wrap_var[:] = glamv
        tmaskpac_wrap_var[:] = tmaskpac
        tmaskatl_wrap_var[:] = tmaskatl
        
        E_0_wrap_var[:] = E
        E_100_wrap_var[:] = E_100
        E_300_wrap_var[:] = E_300
        E_700_wrap_var[:] = E_700
        E_2000_wrap_var[:] = E_2000
        E_vert_wrap_var[:] = E_vert
        E_pac_vert_wrap_var[:] = E_pac_vert
        E_atl_vert_wrap_var[:] = E_atl_vert

        E_int_wrap_var[:] = E_int
        E_100_int_wrap_var[:] = E_100_int
        E_300_int_wrap_var[:] = E_300_int
        E_700_int_wrap_var[:] = E_700_int
        E_2000_int_wrap_var[:] = E_2000_int

        E_pac_int_wrap_var[:] = E_pac_int
        E_pac_100_int_wrap_var[:] = E_pac_100_int
        E_pac_300_int_wrap_var[:] = E_pac_300_int
        E_pac_700_int_wrap_var[:] = E_pac_700_int
        E_pac_2000_int_wrap_var[:] = E_pac_2000_int

        E_atl_int_wrap_var[:] = E_atl_int
        E_atl_100_int_wrap_var[:] = E_atl_100_int
        E_atl_300_int_wrap_var[:] = E_atl_300_int
        E_atl_700_int_wrap_var[:] = E_atl_700_int
        E_atl_2000_int_wrap_var[:] = E_atl_2000_int
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")        
        
        
    def ohc(self):
        """
        Put all the single netCDF files into one file for the sake of postprocessing.
        """
        logging.info("Start wrap-up all the netcdf files of OHC and calculate the zonal mean of each variable.")
        