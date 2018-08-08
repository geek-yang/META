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
import logging
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
        param uc: baratropic corrected zonal winds
        param vc: baratropic corrected meridional winds
        
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
    
    def ncAMET(self, E_0, cpT_0, Lvq_0, gz_0, uv2_0,
               E_200, cpT_200, Lvq_200, gz_200, uv2_200,
               E_500, cpT_500, Lvq_500, gz_500, uv2_500,
               E_850, cpT_850, Lvq_850, gz_850, uv2_850,
               E_vert, cpT_vert, Lvq_vert, gz_vert, uv2_vert,
               year, level, lat, lon, path):
        """
        Save the AMET and its components upto certain levels into netCDF files.
        param E_0: total meridional energy transport over the entire column
        param cpT_0: internal energy transport over the entire column
        param Lvq_0: latent heat transport over the entire column
        param gz_0: geopotential energy transport over the entire column
        param uv2_0: kinetic energy transport over the entire column
        param E_200: total meridional energy transport upto 200 hPa
        param cpT_200: internal energy transport upto 200 hPa
        param Lvq_200: latent heat transport over upto 200 hPa
        param gz_200: geopotential energy transport upto 200 hPa
        param uv2_200: kinetic energy transport upto 200 hPa
        param E_500: total meridional energy transport upto 500 hPa
        param cpT_500: internal energy transport upto 500 hPa
        param Lvq_500: latent heat transport over upto 500 hPa
        param gz_500: geopotential energy transport upto 500 hPa
        param uv2_500: kinetic energy transport upto 500 hPa     
        param E_850: total meridional energy transport upto 850 hPa
        param cpT_850: internal energy transport upto 850 hPa
        param Lvq_850: latent heat transport over upto 850 hPa
        param gz_850: geopotential energy transport upto 850 hPa
        param uv2_850: kinetic energy transport upto 850 hPa
        param E_vert: vertical profile of total meridional energy transport
        param cpT_vert: vertical profile of internal energy transport
        param Lvq_vert: vertical profile of latent heat transport
        param gz_vert: vertical profile of geopotential energy transport
        param uv2_vert: vertical profile of kinetic energy transport
        
        return: netCDF4 files containing AMET and its components on the native grid.
        rtype: netCDF4
        """   
        logging.info("Start creating netcdf file for AMET and its components at each grid point.")
        data_wrap = Dataset(os.path.join(path, 'model_subdaily_075_{}_amet_point.nc'.format(year)),
                            'w',format = 'NETCDF4')
        # create dimensions for netcdf data
        month_wrap_dim = data_wrap.createDimension('month',12)
        level_wrap_dim = data_wrap.createDimension('level',len(level))
        lat_wrap_dim = data_wrap.createDimension('latitude',len(lat))
        lon_wrap_dim = data_wrap.createDimension('longitude',len(lon))
        # create 1-dimension variables
        month_warp_var = data_wrap.createVariable('month',np.int32,('month',))
        level_warp_var = data_wrap.createVariable('level',np.float32,('level',))
        lat_warp_var = data_wrap.createVariable('latitude',np.float32,('latitude',))
        lon_warp_var = data_wrap.createVariable('longitude',np.float32,('longitude',))
        # create 4-dimension variables
        E_0_warp_var = data_wrap.createVariable('E_total',np.float32,('month','latitude','longitude'))
        cpT_0_warp_var = data_wrap.createVariable('cpT_total',np.float32,('month','latitude','longitude'))
        Lvq_0_warp_var = data_wrap.createVariable('Lvq_total',np.float32,('month','latitude','longitude'))
        gz_0_warp_var = data_wrap.createVariable('gz_total',np.float32,('month','latitude','longitude'))
        uv2_0_warp_var = data_wrap.createVariable('uv2_total',np.float32,('month','latitude','longitude'))
        E_200_warp_var = data_wrap.createVariable('E_200',np.float32,('month','latitude','longitude'))
        cpT_200_warp_var = data_wrap.createVariable('cpT_200',np.float32,('month','latitude','longitude'))
        Lvq_200_warp_var = data_wrap.createVariable('Lvq_200',np.float32,('month','latitude','longitude'))
        gz_200_warp_var = data_wrap.createVariable('gz_200',np.float32,('month','latitude','longitude'))
        uv2_200_warp_var = data_wrap.createVariable('uv2_200',np.float32,('month','latitude','longitude'))
        E_500_warp_var = data_wrap.createVariable('E_500',np.float32,('month','latitude','longitude'))
        cpT_500_warp_var = data_wrap.createVariable('cpT_500',np.float32,('month','latitude','longitude'))
        Lvq_500_warp_var = data_wrap.createVariable('Lvq_500',np.float32,('month','latitude','longitude'))
        gz_500_warp_var = data_wrap.createVariable('gz_500',np.float32,('month','latitude','longitude'))
        uv2_500_warp_var = data_wrap.createVariable('uv2_500',np.float32,('month','latitude','longitude'))
        E_850_warp_var = data_wrap.createVariable('E_850',np.float32,('month','latitude','longitude'))
        cpT_850_warp_var = data_wrap.createVariable('cpT_850',np.float32,('month','latitude','longitude'))
        Lvq_850_warp_var = data_wrap.createVariable('Lvq_850',np.float32,('month','latitude','longitude'))
        gz_850_warp_var = data_wrap.createVariable('gz_850',np.float32,('month','latitude','longitude'))
        uv2_850_warp_var = data_wrap.createVariable('uv2_850',np.float32,('month','latitude','longitude'))
        E_vert_warp_var = data_wrap.createVariable('E_vert',np.float32,('month','level','latitude'))
        cpT_vert_warp_var = data_wrap.createVariable('cpT_vert',np.float32,('month','level','latitude'))
        Lvq_vert_warp_var = data_wrap.createVariable('Lvq_vert',np.float32,('month','level','latitude'))
        gz_vert_warp_var = data_wrap.createVariable('gz_vert',np.float32,('month','level','latitude'))
        uv2_vert_warp_var = data_wrap.createVariable('uv2_vert',np.float32,('month','level','latitude'))
        # global attributes
        data_wrap.description = 'Monthly mean meridional energy transport fields.'
        # variable attributes
        E_0_warp_var.units = 'Tera Watt'
        cpT_0_warp_var.units = 'Tera Watt'
        Lvq_0_warp_var.units = 'Tera Watt'
        gz_0_warp_var.units = 'Tera Watt'
        uv2_0_warp_var.units = 'Tera Watt'
        E_200_warp_var.units = 'Tera Watt'
        cpT_200_warp_var.units = 'Tera Watt'
        Lvq_200_warp_var.units = 'Tera Watt'
        gz_200_warp_var.units = 'Tera Watt'
        uv2_200_warp_var.units = 'Tera Watt'
        E_500_warp_var.units = 'Tera Watt'
        cpT_500_warp_var.units = 'Tera Watt'
        Lvq_500_warp_var.units = 'Tera Watt'
        gz_500_warp_var.units = 'Tera Watt'
        uv2_500_warp_var.units = 'Tera Watt'
        E_850_warp_var.units = 'Tera Watt'
        cpT_850_warp_var.units = 'Tera Watt'
        Lvq_850_warp_var.units = 'Tera Watt'
        gz_850_warp_var.units = 'Tera Watt'
        uv2_850_warp_var.units = 'Tera Watt'
        E_vert_warp_var.units = 'Tera Watt'
        cpT_vert_warp_var.units = 'Tera Watt'
        Lvq_vert_warp_var.units = 'Tera Watt'
        gz_vert_warp_var.units = 'Tera Watt'
        uv2_vert_warp_var.units = 'Tera Watt'

        E_0_warp_var.long_name = 'total meridional energy transport over the entire column'
        cpT_0_warp_var.long_name = 'meridional internal energy transport over the entire column'
        Lvq_0_warp_var.long_name = 'meridional latent heat transport over the entire column'
        gz_0_warp_var.long_name = 'meridional geopotential energy transport over the entire column'
        uv2_0_warp_var.long_name = 'meridional kinetic energy transport over the entire column'
        E_200_warp_var.long_name = 'total meridional energy transport from surface to 200hPa'
        cpT_200_warp_var.long_name = 'meridional internal energy transport from surface 200hPa'
        Lvq_200_warp_var.long_name = 'meridional latent heat transport from surface 200hPa'
        gz_200_warp_var.long_name = 'meridional geopotential energy transport from surface 200hPa'
        uv2_200_warp_var.long_name = 'meridional kinetic energy transport from surface 200hPa'
        E_500_warp_var.long_name = 'total meridional energy transport from surface to 500hPa'
        cpT_500_warp_var.long_name = 'meridional internal energy transport from surface 500hPa'
        Lvq_500_warp_var.long_name = 'meridional latent heat transport from surface 500hPa'
        gz_500_warp_var.long_name = 'meridional geopotential energy transport from surface 500hPa'
        uv2_500_warp_var.long_name = 'meridional kinetic energy transport from surface 500hPa'
        E_850_warp_var.long_name = 'total meridional energy transport from surface to 850hPa'
        cpT_850_warp_var.long_name = 'meridional internal energy transport from surface 850hPa'
        Lvq_850_warp_var.long_name = 'meridional latent heat transport from surface 850hPa'
        gz_850_warp_var.long_name = 'meridional geopotential energy transport from surface 850hPa'
        uv2_850_warp_var.long_name = 'meridional kinetic energy transport from surface 850hPa'
        E_vert_warp_var.long_name = 'vertical profile of total meridional energy transport'
        cpT_vert_warp_var.long_name = 'vertical profile of internal energy transport'
        Lvq_vert_warp_var.long_name = 'vertical profile of latent heat transport'
        gz_vert_warp_var.long_name = 'vertical profile of geopotential energy transport'
        uv2_vert_warp_var.long_name = 'vertical profile of kinetic energy transport'
        # writing data
        month_warp_var[:] = np.arange(1,13,1)
        level_warp_var[:] = level
        lat_warp_var[:] = lat
        lon_warp_var[:] = lon
        
        E_0_warp_var[:] = E_0
        cpT_0_warp_var[:] = cpT_0
        Lvq_0_warp_var[:] = Lvq_0
        gz_0_warp_var[:] = gz_0
        uv2_0_warp_var[:] = uv2_0
        E_200_warp_var[:] = E_200
        cpT_200_warp_var[:] = cpT_200
        Lvq_200_warp_var[:] = Lvq_200
        gz_200_warp_var[:] = gz_200
        uv2_200_warp_var[:] = uv2_200
        E_500_warp_var[:] = E_500
        cpT_500_warp_var[:] = cpT_500
        Lvq_500_warp_var[:] = Lvq_500
        gz_500_warp_var[:] = gz_500
        uv2_500_warp_var[:] = uv2_500
        E_850_warp_var[:] = E_850
        cpT_850_warp_var[:] = cpT_850
        Lvq_850_warp_var[:] = Lvq_850
        gz_850_warp_var[:] = gz_850
        uv2_850_warp_var[:] = uv2_850
        E_vert_warp_var[:] = E_vert
        cpT_vert_warp_var[:] = cpT_vert
        Lvq_vert_warp_var[:] = Lvq_vert
        gz_vert_warp_var[:] = gz_vert
        uv2_vert_warp_var[:] = uv2_vert
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")
    
    def ncEddy(self):
        """
        Save the eddy components of AMET and its components upto certain levels into netCDF files.
        """
        logging.info("Create netcdf files successfully!!")
        