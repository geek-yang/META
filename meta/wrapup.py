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
        
        
    def amet(self, name='ERAI'):
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

    def omet(self, name='ORAS4'):
        """
        Put all the single netCDF files into one file for the sake of postprocessing.
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
        data_example = Dataset(os.path.join(self.in_path, '{}_model_monthly_{}_omet_point.nc'.format(alias,self.year_start)
        depth = data_example.variables['depth'][:]
        latitude_aux = data_example.variables['latitude_aux'][:]
        
        # dimension size
        t = len(year)
        z = len(depth)
        jj, ji = .shape
        # create netCDF files
        data_wrap = Dataset(os.path.join(self.out_path, '{}_model_monthly_{}_{}_omet_point.nc'.format(alias, self.year_start, self.year_end)),
                            'w',format = 'NETCDF4')
        # create dimensions for netcdf data
        year_wrap_dim = data_wrap.createDimension('year',)
        month_wrap_dim = data_wrap.createDimension('month',12)
        depth_wrap_dim = data_wrap.createDimension('depth',z)
        lat_wrap_dim = data_wrap.createDimension('jj',jj)
        lon_wrap_dim = data_wrap.createDimension('ji',ji)
        
    def ohc(self):
        """
        Put all the single netCDF files into one file for the sake of postprocessing.
        """
        logging.info("Start wrap-up all the netcdf files of OHC and calculate the zonal mean of each variable.")
        