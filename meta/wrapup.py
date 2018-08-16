# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Save Output Files into NetCDF files
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.09
Last Update     : 2018.08.16
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
        param name: name of the reanalysis products. There are options
        - ERAI (default)
        - MERRA2
        - JRA55
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
        data_example = Dataset(os.path.join(self.in_path,
                               '{}_model_subdaily_{}_amet_point.nc'.format(alias,self.year_start)))
        level = data_example.variables['level'][:]
        lat = data_example.variables['latitude'][:]
        lon = data_example.variables['longitude'][:]
        # create space for the output
        # the number at the end of each name indicates the integral
        # from surface to a certain height (hPa)
        # the results will be saved per year to save memory
        # AMET in the entire column
        E_0 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        cpT_0 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        Lvq_0 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        gz_0 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        uv2_0 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # AMET upto 200hPA
        E_200 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        cpT_200 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        Lvq_200 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        gz_200 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        uv2_200 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # AMET upto 500hPA
        E_500 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        cpT_500 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        Lvq_500 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        gz_500 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        uv2_500 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # AMET upto 850hPA
        E_850 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        cpT_850 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        Lvq_850 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        gz_850 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        uv2_850 = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # AMET vertical profile following the sigma level
        E_vert = np.zeros((len(year),len(month),len(level),len(lat)), dtype=float)
        cpT_vert = np.zeros((len(year),len(month),len(level),len(lat)), dtype=float)
        Lvq_vert = np.zeros((len(year),len(month),len(level),len(lat)), dtype=float)
        gz_vert = np.zeros((len(year),len(month),len(level),len(lat)), dtype=float)
        uv2_vert = np.zeros((len(year),len(month),len(level),len(lat)), dtype=float)
        
        E_0_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        cpT_0_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        Lvq_0_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        gz_0_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        uv2_0_int = np.zeros((len(year),len(month),len(lat)), dtype=float)

        E_200_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        cpT_200_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        Lvq_200_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        gz_200_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        uv2_200_int = np.zeros((len(year),len(month),len(lat)), dtype=float)

        E_500_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        cpT_500_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        Lvq_500_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        gz_500_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        uv2_500_int = np.zeros((len(year),len(month),len(lat)), dtype=float)

        E_850_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        cpT_850_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        Lvq_850_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        gz_850_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        uv2_850_int = np.zeros((len(year),len(month),len(lat)), dtype=float)
        for i in year:
            data_key = Dataset(os.path.join(self.in_path,
                               '{}_model_subdaily_{}_amet_point.nc'.format(alias, i)))
            E_0[i-self.year_start,:,:,:] = data_key.variables['E_total'][:]
            cpT_0[i-self.year_start,:,:,:] = data_key.variables['cpT_total'][:]
            Lvq_0[i-self.year_start,:,:,:] = data_key.variables['Lvq_total'][:]
            gz_0[i-self.year_start,:,:,:] = data_key.variables['gz_total'][:]
            uv2_0[i-self.year_start,:,:,:] = data_key.variables['uv2_total'][:]

            E_200[i-self.year_start,:,:,:] = data_key.variables['E_200'][:]
            cpT_200[i-self.year_start,:,:,:] = data_key.variables['cpT_200'][:]
            Lvq_200[i-self.year_start,:,:,:] = data_key.variables['Lvq_200'][:]
            gz_200[i-self.year_start,:,:,:] = data_key.variables['gz_200'][:]
            uv2_200[i-self.year_start,:,:,:] = data_key.variables['uv2_200'][:]

            E_500[i-self.year_start,:,:,:] = data_key.variables['E_500'][:]
            cpT_500[i-self.year_start,:,:,:] = data_key.variables['cpT_500'][:]
            Lvq_500[i-self.year_start,:,:,:] = data_key.variables['Lvq_500'][:]
            gz_500[i-self.year_start,:,:,:] = data_key.variables['gz_500'][:]
            uv2_500[i-self.year_start,:,:,:] = data_key.variables['uv2_500'][:]

            E_850[i-self.year_start,:,:,:] = data_key.variables['E_850'][:]
            cpT_850[i-self.year_start,:,:,:] = data_key.variables['cpT_850'][:]
            Lvq_850[i-self.year_start,:,:,:] = data_key.variables['Lvq_850'][:]
            gz_850[i-self.year_start,:,:,:] = data_key.variables['gz_850'][:]
            uv2_850[i-self.year_start,:,:,:] = data_key.variables['uv2_850'][:]
            
            E_vert[i-self.year_start,:,:,:] = data_key.variables['E_vert'][:]
            cpT_vert[i-self.year_start,:,:,:] = data_key.variables['cpT_vert'][:]
            Lvq_vert[i-self.year_start,:,:,:] = data_key.variables['Lvq_vert'][:]
            gz_vert[i-self.year_start,:,:,:] = data_key.variables['gz_vert'][:]
            uv2_vert[i-self.year_start,:,:,:] = data_key.variables['uv2_vert'][:]
            
            E_0_int[i-self.year_start,:,:] = np.sum(E_0[i-self.year_start,:,:,:], 2)
            cpT_0_int[i-self.year_start,:,:] = np.sum(cpT_0[i-self.year_start,:,:,:], 2)
            Lvq_0_int[i-self.year_start,:,:] = np.sum(Lvq_0[i-self.year_start,:,:,:], 2)
            gz_0_int[i-self.year_start,:,:] = np.sum(gz_0[i-self.year_start,:,:,:], 2)
            uv2_0_int[i-self.year_start,:,:] = np.sum(uv2_0[i-self.year_start,:,:,:], 2)

            E_200_int[i-self.year_start,:,:] = np.sum(E_200[i-self.year_start,:,:,:], 2)
            cpT_200_int[i-self.year_start,:,:] = np.sum(cpT_200[i-self.year_start,:,:,:], 2)
            Lvq_200_int[i-self.year_start,:,:] = np.sum(Lvq_200[i-self.year_start,:,:,:], 2)
            gz_200_int[i-self.year_start,:,:] = np.sum(gz_200[i-self.year_start,:,:,:], 2)
            uv2_200_int[i-self.year_start,:,:] = np.sum(uv2_200[i-self.year_start,:,:,:], 2)

            E_500_int[i-self.year_start,:,:] = np.sum(E_500[i-self.year_start,:,:,:], 2)
            cpT_500_int[i-self.year_start,:,:] = np.sum(cpT_500[i-self.year_start,:,:,:], 2)
            Lvq_500_int[i-self.year_start,:,:] = np.sum(Lvq_500[i-self.year_start,:,:,:], 2)
            gz_500_int[i-self.year_start,:,:] = np.sum(gz_500[i-self.year_start,:,:,:], 2)
            uv2_500_int[i-self.year_start,:,:] = np.sum(uv2_500[i-self.year_start,:,:,:], 2)

            E_850_int[i-self.year_start,:,:] = np.sum(E_850[i-self.year_start,:,:,:], 2)
            cpT_850_int[i-self.year_start,:,:] = np.sum(cpT_850[i-self.year_start,:,:,:], 2)
            Lvq_850_int[i-self.year_start,:,:] = np.sum(Lvq_850[i-self.year_start,:,:,:], 2)
            gz_850_int[i-self.year_start,:,:] = np.sum(gz_850[i-self.year_start,:,:,:], 2)
            uv2_850_int[i-self.year_start,:,:] = np.sum(uv2_850[i-self.year_start,:,:,:], 2)            
        # create netCDF files
        data_wrap = Dataset(os.path.join(self.out_path,
                            '{}_model_subdaily_{}_{}_amet.nc'.format(alias, self.year_start, self.year_end)),
                            'w',format = 'NETCDF4')
        # create dimensions for netcdf data
        year_wrap_dim = data_wrap.createDimension('year',len(year))
        month_wrap_dim = data_wrap.createDimension('month',12)
        level_wrap_dim = data_wrap.createDimension('level',len(level))
        lat_wrap_dim = data_wrap.createDimension('latitude',len(lat))
        lon_wrap_dim = data_wrap.createDimension('longitude',len(lon))
        # create 1-dimension variables
        year_wrap_var = data_wrap.createVariable('year',np.int32,('year',))
        month_wrap_var = data_wrap.createVariable('month',np.int32,('month',))
        level_wrap_var = data_wrap.createVariable('level',np.float32,('level',))
        lat_wrap_var = data_wrap.createVariable('latitude',np.float32,('latitude',))
        lon_wrap_var = data_wrap.createVariable('longitude',np.float32,('longitude',))
        # create 3-dimension variables
        E_0_wrap_var = data_wrap.createVariable('E_total',np.float32,('year','month','latitude','longitude'),zlib=True)
        cpT_0_wrap_var = data_wrap.createVariable('cpT_total',np.float32,('year','month','latitude','longitude'),zlib=True)
        Lvq_0_wrap_var = data_wrap.createVariable('Lvq_total',np.float32,('year','month','latitude','longitude'),zlib=True)
        gz_0_wrap_var = data_wrap.createVariable('gz_total',np.float32,('year','month','latitude','longitude'),zlib=True)
        uv2_0_wrap_var = data_wrap.createVariable('uv2_total',np.float32,('year','month','latitude','longitude'),zlib=True)
        E_200_wrap_var = data_wrap.createVariable('E_200',np.float32,('year','month','latitude','longitude'),zlib=True)
        cpT_200_wrap_var = data_wrap.createVariable('cpT_200',np.float32,('year','month','latitude','longitude'),zlib=True)
        Lvq_200_wrap_var = data_wrap.createVariable('Lvq_200',np.float32,('year','month','latitude','longitude'),zlib=True)
        gz_200_wrap_var = data_wrap.createVariable('gz_200',np.float32,('year','month','latitude','longitude'),zlib=True)
        uv2_200_wrap_var = data_wrap.createVariable('uv2_200',np.float32,('year','month','latitude','longitude'),zlib=True)
        E_500_wrap_var = data_wrap.createVariable('E_500',np.float32,('year','month','latitude','longitude'),zlib=True)
        cpT_500_wrap_var = data_wrap.createVariable('cpT_500',np.float32,('year','month','latitude','longitude'),zlib=True)
        Lvq_500_wrap_var = data_wrap.createVariable('Lvq_500',np.float32,('year','month','latitude','longitude'),zlib=True)
        gz_500_wrap_var = data_wrap.createVariable('gz_500',np.float32,('year','month','latitude','longitude'),zlib=True)
        uv2_500_wrap_var = data_wrap.createVariable('uv2_500',np.float32,('year','month','latitude','longitude'),zlib=True)
        E_850_wrap_var = data_wrap.createVariable('E_850',np.float32,('year','month','latitude','longitude'),zlib=True)
        cpT_850_wrap_var = data_wrap.createVariable('cpT_850',np.float32,('year','month','latitude','longitude'),zlib=True)
        Lvq_850_wrap_var = data_wrap.createVariable('Lvq_850',np.float32,('year','month','latitude','longitude'),zlib=True)
        gz_850_wrap_var = data_wrap.createVariable('gz_850',np.float32,('year','month','latitude','longitude'),zlib=True)
        uv2_850_wrap_var = data_wrap.createVariable('uv2_850',np.float32,('year','month','latitude','longitude'),zlib=True)
        E_vert_wrap_var = data_wrap.createVariable('E_vert',np.float32,('year','month','level','latitude'),zlib=True)
        cpT_vert_wrap_var = data_wrap.createVariable('cpT_vert',np.float32,('year','month','level','latitude'),zlib=True)
        Lvq_vert_wrap_var = data_wrap.createVariable('Lvq_vert',np.float32,('year','month','level','latitude'),zlib=True)
        gz_vert_wrap_var = data_wrap.createVariable('gz_vert',np.float32,('year','month','level','latitude'),zlib=True)
        uv2_vert_wrap_var = data_wrap.createVariable('uv2_vert',np.float32,('year','month','level','latitude'),zlib=True)
        
        E_0_int_wrap_var = data_wrap.createVariable('E_total_int',np.float32,('year','month','latitude'),zlib=True)
        cpT_0_int_wrap_var = data_wrap.createVariable('cpT_total_int',np.float32,('year','month','latitude'),zlib=True)
        Lvq_0_int_wrap_var = data_wrap.createVariable('Lvq_total_int',np.float32,('year','month','latitude'),zlib=True)
        gz_0_int_wrap_var = data_wrap.createVariable('gz_total_int',np.float32,('year','month','latitude'),zlib=True)
        uv2_0_int_wrap_var = data_wrap.createVariable('uv2_total_int',np.float32,('year','month','latitude'),zlib=True)
        E_200_int_wrap_var = data_wrap.createVariable('E_200_int',np.float32,('year','month','latitude'),zlib=True)
        cpT_200_int_wrap_var = data_wrap.createVariable('cpT_200_int',np.float32,('year','month','latitude'),zlib=True)
        Lvq_200_int_wrap_var = data_wrap.createVariable('Lvq_200_int',np.float32,('year','month','latitude'),zlib=True)
        gz_200_int_wrap_var = data_wrap.createVariable('gz_200_int',np.float32,('year','month','latitude'),zlib=True)
        uv2_200_int_wrap_var = data_wrap.createVariable('uv2_200_int',np.float32,('year','month','latitude'),zlib=True)
        E_500_int_wrap_var = data_wrap.createVariable('E_500_int',np.float32,('year','month','latitude'),zlib=True)
        cpT_500_int_wrap_var = data_wrap.createVariable('cpT_500_int',np.float32,('year','month','latitude'),zlib=True)
        Lvq_500_int_wrap_var = data_wrap.createVariable('Lvq_500_int',np.float32,('year','month','latitude'),zlib=True)
        gz_500_int_wrap_var = data_wrap.createVariable('gz_500_int',np.float32,('year','month','latitude'),zlib=True)
        uv2_500_int_wrap_var = data_wrap.createVariable('uv2_500_int',np.float32,('year','month','latitude'),zlib=True)
        E_850_int_wrap_var = data_wrap.createVariable('E_850_int',np.float32,('year','month','latitude'),zlib=True)
        cpT_850_int_wrap_var = data_wrap.createVariable('cpT_850_int',np.float32,('year','month','latitude'),zlib=True)
        Lvq_850_int_wrap_var = data_wrap.createVariable('Lvq_850_int',np.float32,('year','month','latitude'),zlib=True)
        gz_850_int_wrap_var = data_wrap.createVariable('gz_850_int',np.float32,('year','month','latitude'),zlib=True)
        uv2_850_int_wrap_var = data_wrap.createVariable('uv2_850_int',np.float32,('year','month','latitude'),zlib=True)        
        # global attributes
        data_wrap.description = 'Monthly mean meridional energy transport fields.'        
        # variable attributes
        E_0_wrap_var.units = 'Tera Watt'
        cpT_0_wrap_var.units = 'Tera Watt'
        Lvq_0_wrap_var.units = 'Tera Watt'
        gz_0_wrap_var.units = 'Tera Watt'
        uv2_0_wrap_var.units = 'Tera Watt'
        E_200_wrap_var.units = 'Tera Watt'
        cpT_200_wrap_var.units = 'Tera Watt'
        Lvq_200_wrap_var.units = 'Tera Watt'
        gz_200_wrap_var.units = 'Tera Watt'
        uv2_200_wrap_var.units = 'Tera Watt'
        E_500_wrap_var.units = 'Tera Watt'
        cpT_500_wrap_var.units = 'Tera Watt'
        Lvq_500_wrap_var.units = 'Tera Watt'
        gz_500_wrap_var.units = 'Tera Watt'
        uv2_500_wrap_var.units = 'Tera Watt'
        E_850_wrap_var.units = 'Tera Watt'
        cpT_850_wrap_var.units = 'Tera Watt'
        Lvq_850_wrap_var.units = 'Tera Watt'
        gz_850_wrap_var.units = 'Tera Watt'
        uv2_850_wrap_var.units = 'Tera Watt'
        E_vert_wrap_var.units = 'Tera Watt'
        cpT_vert_wrap_var.units = 'Tera Watt'
        Lvq_vert_wrap_var.units = 'Tera Watt'
        gz_vert_wrap_var.units = 'Tera Watt'
        uv2_vert_wrap_var.units = 'Tera Watt'

        E_0_int_wrap_var.units = 'Tera Watt'
        cpT_0_int_wrap_var.units = 'Tera Watt'
        Lvq_0_int_wrap_var.units = 'Tera Watt'
        gz_0_int_wrap_var.units = 'Tera Watt'
        uv2_0_int_wrap_var.units = 'Tera Watt'
        E_200_int_wrap_var.units = 'Tera Watt'
        cpT_200_int_wrap_var.units = 'Tera Watt'
        Lvq_200_int_wrap_var.units = 'Tera Watt'
        gz_200_int_wrap_var.units = 'Tera Watt'
        uv2_200_int_wrap_var.units = 'Tera Watt'
        E_500_int_wrap_var.units = 'Tera Watt'
        cpT_500_int_wrap_var.units = 'Tera Watt'
        Lvq_500_int_wrap_var.units = 'Tera Watt'
        gz_500_int_wrap_var.units = 'Tera Watt'
        uv2_500_int_wrap_var.units = 'Tera Watt'
        E_850_int_wrap_var.units = 'Tera Watt'
        cpT_850_int_wrap_var.units = 'Tera Watt'
        Lvq_850_int_wrap_var.units = 'Tera Watt'
        gz_850_int_wrap_var.units = 'Tera Watt'
        uv2_850_int_wrap_var.units = 'Tera Watt'        
        
        E_0_wrap_var.long_name = 'total meridional energy transport over the entire column'
        cpT_0_wrap_var.long_name = 'meridional internal energy transport over the entire column'
        Lvq_0_wrap_var.long_name = 'meridional latent heat transport over the entire column'
        gz_0_wrap_var.long_name = 'meridional geopotential energy transport over the entire column'
        uv2_0_wrap_var.long_name = 'meridional kinetic energy transport over the entire column'
        E_200_wrap_var.long_name = 'total meridional energy transport from surface to 200hPa'
        cpT_200_wrap_var.long_name = 'meridional internal energy transport from surface 200hPa'
        Lvq_200_wrap_var.long_name = 'meridional latent heat transport from surface 200hPa'
        gz_200_wrap_var.long_name = 'meridional geopotential energy transport from surface 200hPa'
        uv2_200_wrap_var.long_name = 'meridional kinetic energy transport from surface 200hPa'
        E_500_wrap_var.long_name = 'total meridional energy transport from surface to 500hPa'
        cpT_500_wrap_var.long_name = 'meridional internal energy transport from surface 500hPa'
        Lvq_500_wrap_var.long_name = 'meridional latent heat transport from surface 500hPa'
        gz_500_wrap_var.long_name = 'meridional geopotential energy transport from surface 500hPa'
        uv2_500_wrap_var.long_name = 'meridional kinetic energy transport from surface 500hPa'
        E_850_wrap_var.long_name = 'total meridional energy transport from surface to 850hPa'
        cpT_850_wrap_var.long_name = 'meridional internal energy transport from surface 850hPa'
        Lvq_850_wrap_var.long_name = 'meridional latent heat transport from surface 850hPa'
        gz_850_wrap_var.long_name = 'meridional geopotential energy transport from surface 850hPa'
        uv2_850_wrap_var.long_name = 'meridional kinetic energy transport from surface 850hPa'
        E_vert_wrap_var.long_name = 'vertical profile of total meridional energy transport'
        cpT_vert_wrap_var.long_name = 'vertical profile of internal energy transport'
        Lvq_vert_wrap_var.long_name = 'vertical profile of latent heat transport'
        gz_vert_wrap_var.long_name = 'vertical profile of geopotential energy transport'
        uv2_vert_wrap_var.long_name = 'vertical profile of kinetic energy transport'

        E_0_int_wrap_var.long_name = 'zonal integral of total meridional energy transport over the entire column'
        cpT_0_int_wrap_var.long_name = 'zonal integral of meridional internal energy transport over the entire column'
        Lvq_0_int_wrap_var.long_name = 'zonal integral of meridional latent heat transport over the entire column'
        gz_0_int_wrap_var.long_name = 'zonal integral of meridional geopotential energy transport over the entire column'
        uv2_0_int_wrap_var.long_name = 'zonal integral of meridional kinetic energy transport over the entire column'
        E_200_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 200hPa'
        cpT_200_int_wrap_var.long_name = 'zonal integral of meridional internal energy transport from surface 200hPa'
        Lvq_200_int_wrap_var.long_name = 'zonal integral of meridional latent heat transport from surface 200hPa'
        gz_200_int_wrap_var.long_name = 'zonal integral of meridional geopotential energy transport from surface 200hPa'
        uv2_200_int_wrap_var.long_name = 'zonal integral of meridional kinetic energy transport from surface 200hPa'
        E_500_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 500hPa'
        cpT_500_int_wrap_var.long_name = 'zonal integral of meridional internal energy transport from surface 500hPa'
        Lvq_500_int_wrap_var.long_name = 'zonal integral of meridional latent heat transport from surface 500hPa'
        gz_500_int_wrap_var.long_name = 'zonal integral of meridional geopotential energy transport from surface 500hPa'
        uv2_500_int_wrap_var.long_name = 'zonal integral of meridional kinetic energy transport from surface 500hPa'
        E_850_int_wrap_var.long_name = 'zonal integral of total meridional energy transport from surface to 850hPa'
        cpT_850_int_wrap_var.long_name = 'zonal integral of meridional internal energy transport from surface 850hPa'
        Lvq_850_int_wrap_var.long_name = 'zonal integral of meridional latent heat transport from surface 850hPa'
        gz_850_int_wrap_var.long_name = 'zonal integral of meridional geopotential energy transport from surface 850hPa'
        uv2_850_int_wrap_var.long_name = 'zonal integral of meridional kinetic energy transport from surface 850hPa'
        # writing data
        year_wrap_var[:] = year
        month_wrap_var[:] = np.arange(1,13,1)
        level_wrap_var[:] = level
        lat_wrap_var[:] = lat
        lon_wrap_var[:] = lon
        
        E_0_wrap_var[:] = E_0
        cpT_0_wrap_var[:] = cpT_0
        Lvq_0_wrap_var[:] = Lvq_0
        gz_0_wrap_var[:] = gz_0
        uv2_0_wrap_var[:] = uv2_0
        E_200_wrap_var[:] = E_200
        cpT_200_wrap_var[:] = cpT_200
        Lvq_200_wrap_var[:] = Lvq_200
        gz_200_wrap_var[:] = gz_200
        uv2_200_wrap_var[:] = uv2_200
        E_500_wrap_var[:] = E_500
        cpT_500_wrap_var[:] = cpT_500
        Lvq_500_wrap_var[:] = Lvq_500
        gz_500_wrap_var[:] = gz_500
        uv2_500_wrap_var[:] = uv2_500
        E_850_wrap_var[:] = E_850
        cpT_850_wrap_var[:] = cpT_850
        Lvq_850_wrap_var[:] = Lvq_850
        gz_850_wrap_var[:] = gz_850
        uv2_850_wrap_var[:] = uv2_850
        E_vert_wrap_var[:] = E_vert
        cpT_vert_wrap_var[:] = cpT_vert
        Lvq_vert_wrap_var[:] = Lvq_vert
        gz_vert_wrap_var[:] = gz_vert
        uv2_vert_wrap_var[:] = uv2_vert

        E_0_int_wrap_var[:] = E_0_int
        cpT_0_int_wrap_var[:] = cpT_0_int
        Lvq_0_int_wrap_var[:] = Lvq_0_int
        gz_0_int_wrap_var[:] = gz_0_int
        uv2_0_int_wrap_var[:] = uv2_0_int
        E_200_int_wrap_var[:] = E_200_int
        cpT_200_int_wrap_var[:] = cpT_200_int
        Lvq_200_int_wrap_var[:] = Lvq_200_int
        gz_200_int_wrap_var[:] = gz_200_int
        uv2_200_int_wrap_var[:] = uv2_200_int
        E_500_int_wrap_var[:] = E_500_int
        cpT_500_int_wrap_var[:] = cpT_500_int
        Lvq_500_int_wrap_var[:] = Lvq_500_int
        gz_500_int_wrap_var[:] = gz_500_int
        uv2_500_int_wrap_var[:] = uv2_500_int
        E_850_int_wrap_var[:] = E_850_int
        cpT_850_int_wrap_var[:] = cpT_850_int
        Lvq_850_int_wrap_var[:] = Lvq_850_int
        gz_850_int_wrap_var[:] = gz_850_int
        uv2_850_int_wrap_var[:] = uv2_850_int
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")        

    def ncOMET(self, tmaskpac, tmaskatl, name='ORAS4', temporal='monthly'):
        """
        Put all the single netCDF files into one file for the sake of postprocessing.
        The 3D and 4D variables are compressed.
        param tmaskpac: land-sea mask for the Pacific Ocean  [jj, ji]
        param tmaskatl: land-sea mask for the Atlantic Ocean  [jj, ji]
        param name: name of the reanalysis products. There are options
        - ORAS4 (default)
        - GLORYS2V3
        - SODA3
        param temporal: the temporal resolution that the original calculation is based on. Two options below
        - monthly (default)
        - 5daily 
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
                               '{}_model_{}_{}_omet_point.nc'.format(alias, temporal, self.year_start)))
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
            data_key = Dataset(os.path.join(self.in_path,
                               '{}_model_{}_{}_omet_point.nc'.format(alias, temporal, i)))
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
                            '{}_model_{}_{}_{}_omet.nc'.format(alias, temporal, self.year_start, self.year_end)),
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
        
        
    def ncOHC(self, tmaskpac, tmaskatl, name='ORAS4', temporal='monthly'):
        """
        Put all the single netCDF files into one file for the sake of postprocessing.
        The 3D and 4D variables are compressed.
        param tmaskpac: land-sea mask for the Pacific Ocean  [jj, ji]
        param tmaskatl: land-sea mask for the Atlantic Ocean  [jj, ji]        
        param name: name of the reanalysis products. There are options
        - ORAS4 (default)
        - GLORYS2V3
        - SODA3
        param temporal: the temporal resolution that the original calculation is based on. Two options below
        - monthly (default)
        - 5daily 
        """
        logging.info("Start wrap-up all the netcdf files of OHC and calculate the zonal mean of each variable.")
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
                               '{}_model_{}_{}_ohc_point.nc'.format(alias, temporal, self.year_start)))
        depth = data_example.variables['depth'][:]
        latitude_aux = data_example.variables['latitude_aux'][:]
        gphit = data_example.variables['gphit'][:]
        glamt = data_example.variables['glamt'][:]
        # dimension size
        t = len(year)
        z = len(depth)
        jj, ji = gphit.shape
        # create arrays to store all the data and postprocess
        OHC = np.zeros((t, len(month), jj, ji), dtype=float)
        OHC_100 = np.zeros((t, len(month), jj, ji), dtype=float)
        OHC_300 = np.zeros((t, len(month), jj, ji), dtype=float)
        OHC_700 = np.zeros((t, len(month), jj, ji), dtype=float)
        OHC_2000 = np.zeros((t, len(month), jj, ji), dtype=float)
        OHC_pac_vert = np.zeros((t, len(month), z, jj), dtype=float)
        OHC_atl_vert = np.zeros((t, len(month), z, jj), dtype=float)
        OHC_vert = np.zeros((t, len(month), z, jj), dtype=float)
        
        OHC_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_100_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_300_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_700_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_2000_int = np.zeros((t, len(month), jj), dtype=float)    
        
        OHC_pac_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_pac_100_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_pac_300_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_pac_700_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_pac_2000_int = np.zeros((t, len(month), jj), dtype=float)

        OHC_atl_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_atl_100_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_atl_300_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_atl_700_int = np.zeros((t, len(month), jj), dtype=float)
        OHC_atl_2000_int = np.zeros((t, len(month), jj), dtype=float)
        # increase the dimension of sub-basin masks for the sake of efficiency
        tmaskpac_3D = np.repeat(tmaskpac[np.newaxis,:,:],len(month),0)
        tmaskatl_3D = np.repeat(tmaskatl[np.newaxis,:,:],len(month),0)
        for i in year:
            data_key = Dataset(os.path.join(self.in_path,
                               '{}_model_{}_{}_ohc_point.nc'.format(alias, temporal, i)))
            OHC[i-self.year_start,:,:,:] = data_key.variables['OHC_total'][:]
            OHC_100[i-self.year_start,:,:,:] = data_key.variables['OHC_100'][:]
            OHC_300[i-self.year_start,:,:,:] = data_key.variables['OHC_300'][:]
            OHC_700[i-self.year_start,:,:,:] = data_key.variables['OHC_700'][:]
            OHC_2000[i-self.year_start,:,:,:] = data_key.variables['OHC_2000'][:]
            OHC_pac_vert[i-self.year_start,:,:,:] = data_key.variables['OHC_pac_vert'][:]
            OHC_atl_vert[i-self.year_start,:,:,:] = data_key.variables['OHC_atl_vert'][:]
            OHC_vert[i-self.year_start,:,:,:] = data_key.variables['OHC_vert'][:]
            
            OHC_int[i-self.year_start,:,:] = np.sum(OHC[i-self.year_start,:,:,:], 2)
            OHC_100_int[i-self.year_start,:,:] = np.sum(OHC_100[i-self.year_start,:,:,:], 2)
            OHC_300_int[i-self.year_start,:,:] = np.sum(OHC_300[i-self.year_start,:,:,:], 2)
            OHC_700_int[i-self.year_start,:,:] = np.sum(OHC_700[i-self.year_start,:,:,:], 2)
            OHC_2000_int[i-self.year_start,:,:] = np.sum(OHC_2000[i-self.year_start,:,:,:], 2)
            
            OHC_pac_int[i-self.year_start,:,:] = np.sum(OHC[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            OHC_pac_100_int[i-self.year_start,:,:] = np.sum(OHC_100[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            OHC_pac_300_int[i-self.year_start,:,:] = np.sum(OHC_300[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            OHC_pac_700_int[i-self.year_start,:,:] = np.sum(OHC_700[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            OHC_pac_2000_int[i-self.year_start,:,:] = np.sum(OHC_2000[i-self.year_start,:,:,:] * tmaskpac_3D, 2)
            
            OHC_atl_int[i-self.year_start,:,:] = np.sum(OHC[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            OHC_atl_100_int[i-self.year_start,:,:] = np.sum(OHC_100[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            OHC_atl_300_int[i-self.year_start,:,:] = np.sum(OHC_300[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            OHC_atl_700_int[i-self.year_start,:,:] = np.sum(OHC_700[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
            OHC_atl_2000_int[i-self.year_start,:,:] = np.sum(OHC_2000[i-self.year_start,:,:,:] * tmaskatl_3D, 2)
        # create netCDF files
        data_wrap = Dataset(os.path.join(self.out_path,
                            '{}_model_{}_{}_{}_ohc.nc'.format(alias, temporal, self.year_start, self.year_end)),
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
        gphit_wrap_var = data_wrap.createVariable('gphit',np.float32,('jj','ji'))
        glamt_wrap_var = data_wrap.createVariable('glamt',np.float32,('jj','ji'))
        tmaskpac_wrap_var = data_wrap.createVariable('tmaskpac',np.int32,('jj','ji'))
        tmaskatl_wrap_var = data_wrap.createVariable('tmaskatl',np.int32,('jj','ji'))
        # create 3-dimension variables
        OHC_int_wrap_var = data_wrap.createVariable('OHC_int',np.float32,('year','month','jj'),zlib=True)
        OHC_100_int_wrap_var = data_wrap.createVariable('OHC_100_int',np.float32,('year','month','jj'),zlib=True)
        OHC_300_int_wrap_var = data_wrap.createVariable('OHC_300_int',np.float32,('year','month','jj'),zlib=True)
        OHC_700_int_wrap_var = data_wrap.createVariable('OHC_700_int',np.float32,('year','month','jj'),zlib=True)
        OHC_2000_int_wrap_var = data_wrap.createVariable('OHC_2000_int',np.float32,('year','month','jj'),zlib=True)

        OHC_pac_int_wrap_var = data_wrap.createVariable('OHC_pac_int',np.float32,('year','month','jj'),zlib=True)
        OHC_pac_100_int_wrap_var = data_wrap.createVariable('OHC_pac_100_int',np.float32,('year','month','jj'),zlib=True)
        OHC_pac_300_int_wrap_var = data_wrap.createVariable('OHC_pac_300_int',np.float32,('year','month','jj'),zlib=True)
        OHC_pac_700_int_wrap_var = data_wrap.createVariable('OHC_pac_700_int',np.float32,('year','month','jj'),zlib=True)
        OHC_pac_2000_int_wrap_var = data_wrap.createVariable('OHC_pac_2000_int',np.float32,('year','month','jj'),zlib=True)

        OHC_atl_int_wrap_var = data_wrap.createVariable('OHC_atl_int',np.float32,('year','month','jj'),zlib=True)
        OHC_atl_100_int_wrap_var = data_wrap.createVariable('OHC_atl_100_int',np.float32,('year','month','jj'),zlib=True)
        OHC_atl_300_int_wrap_var = data_wrap.createVariable('OHC_atl_300_int',np.float32,('year','month','jj'),zlib=True)
        OHC_atl_700_int_wrap_var = data_wrap.createVariable('OHC_atl_700_int',np.float32,('year','month','jj'),zlib=True)
        OHC_atl_2000_int_wrap_var = data_wrap.createVariable('OHC_atl_2000_int',np.float32,('year','month','jj'),zlib=True)    
        # create 4-dimension variables
        OHC_0_wrap_var = data_wrap.createVariable('OHC_total',np.float32,('year','month','jj','ji'),zlib=True)
        OHC_100_wrap_var = data_wrap.createVariable('OHC_100',np.float32,('year','month','jj','ji'),zlib=True)
        OHC_300_wrap_var = data_wrap.createVariable('OHC_300',np.float32,('year','month','jj','ji'),zlib=True)
        OHC_700_wrap_var = data_wrap.createVariable('OHC_700',np.float32,('year','month','jj','ji'),zlib=True)
        OHC_2000_wrap_var = data_wrap.createVariable('OHC_2000',np.float32,('year','month','jj','ji'),zlib=True)
        OHC_vert_wrap_var = data_wrap.createVariable('OHC_vert',np.float32,('year','month','depth','jj'),zlib=True)
        OHC_pac_vert_wrap_var = data_wrap.createVariable('OHC_pac_vert',np.float32,('year','month','depth','jj'),zlib=True)
        OHC_atl_vert_wrap_var = data_wrap.createVariable('OHC_atl_vert',np.float32,('year','month','depth','jj'),zlib=True)        
        # global attributes
        data_wrap.description = 'Monthly mean ocean heat content fields.'
        # variable attributes
        OHC_0_wrap_var.units = 'Tera Joule'
        OHC_100_wrap_var.units = 'Tera Joule'
        OHC_300_wrap_var.units = 'Tera Joule'
        OHC_700_wrap_var.units = 'Tera Joule'
        OHC_2000_wrap_var.units = 'Tera Joule'
        OHC_vert_wrap_var.units = 'Tera Joule'
        OHC_pac_vert_wrap_var.units = 'Tera Joule'
        OHC_atl_vert_wrap_var.units = 'Tera Joule'
        
        OHC_int_wrap_var.units = 'Tera Joule'
        OHC_100_int_wrap_var.units = 'Tera Joule'
        OHC_300_int_wrap_var.units = 'Tera Joule'
        OHC_700_int_wrap_var.units = 'Tera Joule'
        OHC_2000_int_wrap_var.units = 'Tera Joule'
        
        OHC_pac_int_wrap_var.units = 'Tera Joule'
        OHC_pac_100_int_wrap_var.units = 'Tera Joule'
        OHC_pac_300_int_wrap_var.units = 'Tera Joule'
        OHC_pac_700_int_wrap_var.units = 'Tera Joule'
        OHC_pac_2000_int_wrap_var.units = 'Tera Joule'

        OHC_atl_int_wrap_var.units = 'Tera Joule'
        OHC_atl_100_int_wrap_var.units = 'Tera Joule'
        OHC_atl_300_int_wrap_var.units = 'Tera Joule'
        OHC_atl_700_int_wrap_var.units = 'Tera Joule'
        OHC_atl_2000_int_wrap_var.units = 'Tera Joule'
        
        tmaskpac_wrap_var.long_name = 'land sea mask for the Pacific basin'
        tmaskatl_wrap_var.long_name = 'land sea mask for the Atlantic basin'
        
        OHC_0_wrap_var.long_name = 'total ocean heat content over the entire depth'
        OHC_100_wrap_var.long_name = 'total ocean heat content from surface to 100m'
        OHC_300_wrap_var.long_name = 'total ocean heat content from surface to 300m'
        OHC_700_wrap_var.long_name = 'total ocean heat content from surface to 700m'
        OHC_2000_wrap_var.long_name = 'total ocean heat content from surface to 2000m'
        OHC_vert_wrap_var.long_name = 'zonal integral of ocean heat content for the globe'
        OHC_pac_vert_wrap_var.long_name = 'zonal integral of ocean heat content in the Pacific Ocean'
        OHC_pac_vert_wrap_var.long_name = 'zonal integral of ocean heat content in the Pacific Ocean'
        
        OHC_int_wrap_var.long_name = 'zonal integral of total ocean heat content over the entire depth'
        OHC_100_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 100m'
        OHC_300_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 300m'
        OHC_700_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 700m'
        OHC_2000_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 2000m'
        
        OHC_pac_int_wrap_var.long_name = 'zonal integral of total ocean heat content over the entire depth in the Pacific'
        OHC_pac_100_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 100m in the Pacific'
        OHC_pac_300_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 300m in the Pacific'
        OHC_pac_700_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 700m in the Pacific'
        OHC_pac_2000_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 2000m in the Pacific'
        
        OHC_atl_int_wrap_var.long_name = 'zonal integral of total ocean heat content over the entire depth in the Atlantic'
        OHC_atl_100_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 100m in the Atlantic'
        OHC_atl_300_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 300m in the Atlantic'
        OHC_atl_700_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 700m in the Atlantic'
        OHC_atl_2000_int_wrap_var.long_name = 'zonal integral of total ocean heat content from surface to 2000m in the Atlantic'
        # writing data
        year_wrap_var[:] = year
        month_wrap_var[:] = np.arange(1,13,1)
        depth_wrap_var[:] = depth
        lat_wrap_var[:] = latitude_aux
        gphit_wrap_var[:] = gphit
        glamt_wrap_var[:] = glamt
        tmaskpac_wrap_var[:] = tmaskpac
        tmaskatl_wrap_var[:] = tmaskatl
        
        OHC_0_wrap_var[:] = OHC
        OHC_100_wrap_var[:] = OHC_100
        OHC_300_wrap_var[:] = OHC_300
        OHC_700_wrap_var[:] = OHC_700
        OHC_2000_wrap_var[:] = OHC_2000
        OHC_vert_wrap_var[:] = OHC_vert
        OHC_pac_vert_wrap_var[:] = OHC_pac_vert
        OHC_atl_vert_wrap_var[:] = OHC_atl_vert

        OHC_int_wrap_var[:] = OHC_int
        OHC_100_int_wrap_var[:] = OHC_100_int
        OHC_300_int_wrap_var[:] = OHC_300_int
        OHC_700_int_wrap_var[:] = OHC_700_int
        OHC_2000_int_wrap_var[:] = OHC_2000_int

        OHC_pac_int_wrap_var[:] = OHC_pac_int
        OHC_pac_100_int_wrap_var[:] = OHC_pac_100_int
        OHC_pac_300_int_wrap_var[:] = OHC_pac_300_int
        OHC_pac_700_int_wrap_var[:] = OHC_pac_700_int
        OHC_pac_2000_int_wrap_var[:] = OHC_pac_2000_int

        OHC_atl_int_wrap_var[:] = OHC_atl_int
        OHC_atl_100_int_wrap_var[:] = OHC_atl_100_int
        OHC_atl_300_int_wrap_var[:] = OHC_atl_300_int
        OHC_atl_700_int_wrap_var[:] = OHC_atl_700_int
        OHC_atl_2000_int_wrap_var[:] = OHC_atl_2000_int
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")
        