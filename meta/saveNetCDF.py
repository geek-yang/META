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

    def ncCorrect(self, uc, vc, year, lat, lon, path, name='ERAI'):
        """
        Save the baratropic corrected wind into netCDF files.
        param uc: baratropic corrected zonal winds
        param vc: baratropic corrected meridional winds

        return: netCDF4 files containing uc and vc on the native grid.
        rtype: netCDF4
        """
        logging.info("Start creating netcdf file for baratropic corrected wind fields at each grid point.")
        if name == 'ERAI':
            data_wrap = Dataset(os.path.join(path, 'era_model_subdaily_{}_{}_uvc_point.nc'.format(year[0],year[-1])),
                                'w',format = 'NETCDF4')
        elif name == 'MERRA2':
            data_wrap = Dataset(os.path.join(path, 'merra_model_subdaily_{}_{}_uvc_point.nc'.format(year[0],year[-1])),
                                'w',format = 'NETCDF4')
        elif name == 'JRA55':
            data_wrap = Dataset(os.path.join(path, 'jra_model_subdaily_{}_{}_uvc_point.nc'.format(year[0],year[-1])),
                                'w',format = 'NETCDF4')
        else:
            raise IOError("This dataset is not supported in this module.")
        # create dimensions for netcdf data
        year_wrap_dim = data_wrap.createDimension('year',len(year))
        month_wrap_dim = data_wrap.createDimension('month',12)
        lat_wrap_dim = data_wrap.createDimension('latitude',len(lat))
        lon_wrap_dim = data_wrap.createDimension('longitude',len(lon))
        # create 1-dimension variables
        year_wrap_var = data_wrap.createVariable('year',np.int32,('year',))
        month_wrap_var = data_wrap.createVariable('month',np.int32,('month',))
        lat_wrap_var = data_wrap.createVariable('latitude',np.float32,('latitude',))
        lon_wrap_var = data_wrap.createVariable('longitude',np.float32,('longitude',))
        # create 4-dimension variables
        uc_wrap_var = data_wrap.createVariable('uc',np.float32,('year','month','latitude','longitude'))
        vc_wrap_var = data_wrap.createVariable('vc',np.float32,('year','month','latitude','longitude'))
        # global attributes
        data_wrap.description = 'Monthly mean baratropic corrected wind fields.'
        # variable attributes
        lat_wrap_var.units = 'degree_north'
        lon_wrap_var.units = 'degree_east'
        uc_wrap_var.units = 'm/s'
        vc_wrap_var.units = 'm/s'

        uc_wrap_var.long_name = 'zonal barotropic correction wind'
        vc_wrap_var.long_name = 'meridional barotropic correction wind'
        # writing data
        year_wrap_var[:] = year
        month_wrap_var[:] = np.arange(1,13,1)
        lat_wrap_var[:] = lat
        lon_wrap_var[:] = lon
        uc_wrap_var[:] = uc
        vc_wrap_var[:] = vc
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")

    def ncInter(self,):
        """
        Save the intermediate fields for the computation of mass correction with
        NCL on spherical harmonics.
        """
        

    def ncAMET(self, E_0, cpT_0, Lvq_0, gz_0, uv2_0,
               E_200, cpT_200, Lvq_200, gz_200, uv2_200,
               E_500, cpT_500, Lvq_500, gz_500, uv2_500,
               E_850, cpT_850, Lvq_850, gz_850, uv2_850,
               E_vert, cpT_vert, Lvq_vert, gz_vert, uv2_vert,
               year, level, lat, lon, path, name='ERAI'):
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
        if name == 'ERAI':
            data_wrap = Dataset(os.path.join(path, 'era_model_subdaily_{}_amet_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        elif name == 'MERRA2':
            data_wrap = Dataset(os.path.join(path, 'merra_model_subdaily_{}_amet_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        elif name == 'JRA55':
            data_wrap = Dataset(os.path.join(path, 'jra_model_subdaily_{}_amet_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        else:
            raise IOError("This dataset is not supported in this module.")
        # create dimensions for netcdf data
        month_wrap_dim = data_wrap.createDimension('month',12)
        level_wrap_dim = data_wrap.createDimension('level',len(level))
        lat_wrap_dim = data_wrap.createDimension('latitude',len(lat))
        lon_wrap_dim = data_wrap.createDimension('longitude',len(lon))
        # create 1-dimension variables
        month_wrap_var = data_wrap.createVariable('month',np.int32,('month',))
        level_wrap_var = data_wrap.createVariable('level',np.float32,('level',))
        lat_wrap_var = data_wrap.createVariable('latitude',np.float32,('latitude',))
        lon_wrap_var = data_wrap.createVariable('longitude',np.float32,('longitude',))
        # create 3-dimension variables
        E_0_wrap_var = data_wrap.createVariable('E_total',np.float32,('month','latitude','longitude'))
        cpT_0_wrap_var = data_wrap.createVariable('cpT_total',np.float32,('month','latitude','longitude'))
        Lvq_0_wrap_var = data_wrap.createVariable('Lvq_total',np.float32,('month','latitude','longitude'))
        gz_0_wrap_var = data_wrap.createVariable('gz_total',np.float32,('month','latitude','longitude'))
        uv2_0_wrap_var = data_wrap.createVariable('uv2_total',np.float32,('month','latitude','longitude'))
        E_200_wrap_var = data_wrap.createVariable('E_200',np.float32,('month','latitude','longitude'))
        cpT_200_wrap_var = data_wrap.createVariable('cpT_200',np.float32,('month','latitude','longitude'))
        Lvq_200_wrap_var = data_wrap.createVariable('Lvq_200',np.float32,('month','latitude','longitude'))
        gz_200_wrap_var = data_wrap.createVariable('gz_200',np.float32,('month','latitude','longitude'))
        uv2_200_wrap_var = data_wrap.createVariable('uv2_200',np.float32,('month','latitude','longitude'))
        E_500_wrap_var = data_wrap.createVariable('E_500',np.float32,('month','latitude','longitude'))
        cpT_500_wrap_var = data_wrap.createVariable('cpT_500',np.float32,('month','latitude','longitude'))
        Lvq_500_wrap_var = data_wrap.createVariable('Lvq_500',np.float32,('month','latitude','longitude'))
        gz_500_wrap_var = data_wrap.createVariable('gz_500',np.float32,('month','latitude','longitude'))
        uv2_500_wrap_var = data_wrap.createVariable('uv2_500',np.float32,('month','latitude','longitude'))
        E_850_wrap_var = data_wrap.createVariable('E_850',np.float32,('month','latitude','longitude'))
        cpT_850_wrap_var = data_wrap.createVariable('cpT_850',np.float32,('month','latitude','longitude'))
        Lvq_850_wrap_var = data_wrap.createVariable('Lvq_850',np.float32,('month','latitude','longitude'))
        gz_850_wrap_var = data_wrap.createVariable('gz_850',np.float32,('month','latitude','longitude'))
        uv2_850_wrap_var = data_wrap.createVariable('uv2_850',np.float32,('month','latitude','longitude'))
        E_vert_wrap_var = data_wrap.createVariable('E_vert',np.float32,('month','level','latitude'))
        cpT_vert_wrap_var = data_wrap.createVariable('cpT_vert',np.float32,('month','level','latitude'))
        Lvq_vert_wrap_var = data_wrap.createVariable('Lvq_vert',np.float32,('month','level','latitude'))
        gz_vert_wrap_var = data_wrap.createVariable('gz_vert',np.float32,('month','level','latitude'))
        uv2_vert_wrap_var = data_wrap.createVariable('uv2_vert',np.float32,('month','level','latitude'))
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
        # writing data
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
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")

    def ncEddyamet(self):
        """
        Save the eddy components of AMET and its components upto certain levels into netCDF files.
        """
        logging.info("Create netcdf files successfully!!")

    def ncOMET(self, E, E_100, E_300, E_700, E_2000, E_vert, E_pac_vert, E_atl_vert,
               year, gphiv, glamv, nav_lev, z, jj, ji, path, name='ORAS4'):
        """
        Save the OMET upto certain levels into netCDF files.
        param E: total meridional energy transport over the entire depth
        param E_100: total meridional energy transport from sea surface to 100m
        param E_300: total meridional energy transport from sea surface to 300m
        param E_700: total meridional energy transport from sea surface to 700m
        param E_2000: total meridional energy transport from sea surface to 2000m
        param E_vert: vertical profile of total meridional energy transport
        param E_pac_vert: vertical profile of total meridional energy transport in the Pacific Ocean
        param E_atl_vert: vertical profile of total meridional energy transport in the Atlantic Ocean

        return: netCDF4 files containing OMET on the native grid.
        rtype: netCDF4
        """
        logging.info("Start creating netcdf file for OMET at each grid point.")
        if name == 'ORAS4':
            data_wrap = Dataset(os.path.join(path, 'oras_model_monthly_{}_omet_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        elif name == 'GLOTYS2V3':
            data_wrap = Dataset(os.path.join(path, 'glorys_model_monthly_{}_omet_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        elif name == 'SODA3':
            data_wrap = Dataset(os.path.join(path, 'soda_model_5daily_{}_omet_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        else:
            raise IOError("This dataset is not supported in this module.")
         # create dimensions for netcdf data
        month_wrap_dim = data_wrap.createDimension('month',12)
        depth_wrap_dim = data_wrap.createDimension('depth',z)
        lat_wrap_dim = data_wrap.createDimension('jj',jj)
        lon_wrap_dim = data_wrap.createDimension('ji',ji)
        # create 1-dimension variables
        month_wrap_var = data_wrap.createVariable('month',np.int32,('month',))
        depth_wrap_var = data_wrap.createVariable('depth',np.float32,('depth',))
        lat_wrap_var = data_wrap.createVariable('latitude_aux',np.float32,('jj',))
        # create 2-dimension variables
        gphiv_wrap_var = data_wrap.createVariable('gphiv',np.float32,('jj','ji'))
        glamv_wrap_var = data_wrap.createVariable('glamv',np.float32,('jj','ji'))
        # create 3-dimension variables
        E_0_wrap_var = data_wrap.createVariable('E_total',np.float32,('month','jj','ji'))
        E_100_wrap_var = data_wrap.createVariable('E_100',np.float32,('month','jj','ji'))
        E_300_wrap_var = data_wrap.createVariable('E_300',np.float32,('month','jj','ji'))
        E_700_wrap_var = data_wrap.createVariable('E_700',np.float32,('month','jj','ji'))
        E_2000_wrap_var = data_wrap.createVariable('E_2000',np.float32,('month','jj','ji'))
        E_vert_wrap_var = data_wrap.createVariable('E_vert',np.float32,('month','depth','jj'))
        E_pac_vert_wrap_var = data_wrap.createVariable('E_pac_vert',np.float32,('month','depth','jj'))
        E_atl_vert_wrap_var = data_wrap.createVariable('E_atl_vert',np.float32,('month','depth','jj'))
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

        E_0_wrap_var.long_name = 'total meridional energy transport over the entire column'
        E_100_wrap_var.long_name = 'total meridional energy transport from surface to 100m'
        E_300_wrap_var.long_name = 'total meridional energy transport from surface to 300m'
        E_700_wrap_var.long_name = 'total meridional energy transport from surface to 700m'
        E_2000_wrap_var.long_name = 'total meridional energy transport from surface to 2000m'
        E_vert_wrap_var.long_name = 'zonal integral of meridional energy transport for the globe'
        E_pac_vert_wrap_var.long_name = 'zonal integral of meridional energy transport in the Pacific Ocean'
        E_pac_vert_wrap_var.long_name = 'zonal integral of meridional energy transport in the Pacific Ocean'
        # writing data
        month_wrap_var[:] = np.arange(1,13,1)
        depth_wrap_var[:] = nav_lev
        lat_wrap_var[:] = gphiv[:,96]
        gphiv_wrap_var[:] = gphiv
        glamv_wrap_var[:] = glamv

        E_0_wrap_var[:] = E
        E_100_wrap_var[:] = E_100
        E_300_wrap_var[:] = E_300
        E_700_wrap_var[:] = E_700
        E_2000_wrap_var[:] = E_2000
        E_vert_wrap_var[:] = E_vert
        E_pac_vert_wrap_var[:] = E_pac_vert
        E_atl_vert_wrap_var[:] = E_atl_vert
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")

    def ncOHC(self, OHC, OHC_100, OHC_300, OHC_700, OHC_2000, OHC_vert, OHC_pac_vert, OHC_atl_vert,
               year, nav_lat, nav_lon, nav_lev, gphiv, z, jj, ji, path, name='ORAS4'):
        """
        Save the OHC upto certain levels into netCDF files.
        param OHC: total ocean heat content over the entire depth
        param OHC_100: total ocean heat content from surface to 100m
        param OHC_300: total ocean heat content from surface to 300m
        param OHC_700: total ocean heat content from surface to 700m
        param OHC_2000: total ocean heat content from surface to 2000m
        param OHC_vert: Vertical profile of ocean heat content
        param OHC_pac_vert: Vertical profile of ocean heat content in the Pacific Ocean
        param OHC_atl_vert: Vertical profile of ocean heat content in the Atlantic Ocean

        return: netCDF4 files containing OHC on the native grid.
        rtype: netCDF4
        """
        logging.info("Start creating netcdf file for OHC at each grid point.")
        if name == 'ORAS4':
            data_wrap = Dataset(os.path.join(path, 'oras_model_monthly_{}_ohc_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        elif name == 'GLOTYS2V3':
            data_wrap = Dataset(os.path.join(path, 'glorys_model_monthly_{}_ohc_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        elif name == 'SODA3':
            data_wrap = Dataset(os.path.join(path, 'soda_model_5daily_{}_ohc_point.nc'.format(year)),
                                'w',format = 'NETCDF4')
        else:
            raise IOError("This dataset is not supported in this module.")
        # create dimensions for netcdf data
        month_wrap_dim = data_wrap.createDimension('month',12)
        depth_wrap_dim = data_wrap.createDimension('depth',z)
        lat_wrap_dim = data_wrap.createDimension('jj',jj)
        lon_wrap_dim = data_wrap.createDimension('ji',ji)
        # create 1-dimension variables
        month_wrap_var = data_wrap.createVariable('month',np.int32,('month',))
        depth_wrap_var = data_wrap.createVariable('depth',np.float32,('depth',))
        lat_wrap_var = data_wrap.createVariable('latitude_aux',np.float32,('jj',))
        # create 2-dimension variables
        gphit_wrap_var = data_wrap.createVariable('gphit',np.float32,('jj','ji'))
        glamt_wrap_var = data_wrap.createVariable('glamt',np.float32,('jj','ji'))
        # create 3-dimension variables
        OHC_0_wrap_var = data_wrap.createVariable('OHC_total',np.float32,('month','jj','ji'))
        OHC_100_wrap_var = data_wrap.createVariable('OHC_100',np.float32,('month','jj','ji'))
        OHC_300_wrap_var = data_wrap.createVariable('OHC_300',np.float32,('month','jj','ji'))
        OHC_700_wrap_var = data_wrap.createVariable('OHC_700',np.float32,('month','jj','ji'))
        OHC_2000_wrap_var = data_wrap.createVariable('OHC_2000',np.float32,('month','jj','ji'))
        OHC_vert_wrap_var = data_wrap.createVariable('OHC_vert',np.float32,('month','depth','jj'))
        OHC_pac_vert_wrap_var = data_wrap.createVariable('OHC_pac_vert',np.float32,('month','depth','jj'))
        OHC_atl_vert_wrap_var = data_wrap.createVariable('OHC_atl_vert',np.float32,('month','depth','jj'))
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

        OHC_0_wrap_var.long_name = 'total ocean heat content over the entire depth'
        OHC_100_wrap_var.long_name = 'total ocean heat content from surface to 100m'
        OHC_300_wrap_var.long_name = 'total ocean heat content from surface to 300m'
        OHC_700_wrap_var.long_name = 'total ocean heat content from surface to 700m'
        OHC_2000_wrap_var.long_name = 'total ocean heat content from surface to 2000m'
        OHC_vert_wrap_var.long_name = 'zonal integral of ocean heat content for the globe'
        OHC_pac_vert_wrap_var.long_name = 'zonal integral of ocean heat content in the Pacific Ocean'
        OHC_pac_vert_wrap_var.long_name = 'zonal integral of ocean heat content in the Pacific Ocean'
        # writing data
        month_wrap_var[:] = np.arange(1,13,1)
        depth_wrap_var[:] = nav_lev
        lat_wrap_var[:] = gphiv[:,96]
        gphit_wrap_var[:] = nav_lat
        glamt_wrap_var[:] = nav_lon

        OHC_0_wrap_var[:] = OHC
        OHC_100_wrap_var[:] = OHC_100
        OHC_300_wrap_var[:] = OHC_300
        OHC_700_wrap_var[:] = OHC_700
        OHC_2000_wrap_var[:] = OHC_2000
        OHC_vert_wrap_var[:] = OHC_vert
        OHC_pac_vert_wrap_var[:] = OHC_pac_vert
        OHC_atl_vert_wrap_var[:] = OHC_atl_vert
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")

    def ncEddyomet(self, E_eddy_steady_mean, E_100_eddy_steady_mean, E_300_eddy_steady_mean,
                   E_700_eddy_steady_mean, E_2000_eddy_steady_mean, E_pac_eddy_steady_mean,
                   E_pac_100_eddy_steady_mean, E_pac_300_eddy_steady_mean, E_pac_700_eddy_steady_mean,
                   E_pac_2000_eddy_steady_mean, E_atl_eddy_steady_mean, E_atl_100_eddy_steady_mean,
                   E_atl_300_eddy_steady_mean, E_atl_700_eddy_steady_mean, E_atl_2000_eddy_steady_mean,
                   E_eddy_stationary_mean, E_100_eddy_stationary_mean, E_300_eddy_stationary_mean,
                   E_700_eddy_stationary_mean, E_2000_eddy_stationary_mean, E_vert_eddy_stationary_mean,
                   E_pac_eddy_stationary_mean, E_pac_100_eddy_stationary_mean, E_pac_300_eddy_stationary_mean,
                   E_pac_700_eddy_stationary_mean, E_pac_2000_eddy_stationary_mean,
                   E_pac_vert_eddy_stationary_mean, E_atl_eddy_stationary_mean, E_atl_100_eddy_stationary_mean,
                   E_atl_300_eddy_stationary_mean, E_atl_700_eddy_stationary_mean,
                   E_atl_2000_eddy_stationary_mean, E_atl_vert_eddy_stationary_mean,
                   year, gphiv, glamv, nav_lev, z, jj, ji, path, name='ORAS4'):
        """
        Save the eddy components of OMET and its components upto certain levels into netCDF files.
        param E_eddy_steady_mean: energy transport by steady mean flow over the entire depth
        param E_100_eddy_steady_mean: energy transport by steady mean flow from surface to 100m

        param name: name of the reanalysis products. There are options
        - ORAS4 (default)
        - GLORYS2V3
        - SODA3
        param temporal: the temporal resolution that the original calculation is based on. Two options below
        - monthly (default)
        - 5daily
        """
        logging.info("Start creating netcdf file for eddy components at each grid point.")
        if name == 'ORAS4':
            alias = 'oras'
        elif name == 'GLOTYS2V3':
            alias = 'glorys'
        elif name == 'SODA3':
            alias = 'soda'
        else:
            raise IOError("This dataset is not supported in this module.")

        data_wrap = Dataset(os.path.join(path,'{}_model_monthly_{}_E_eddy_point.nc'.format(alias, year)),
                            'w',format = 'NETCDF4')
        # create dimensions for netcdf data
        month_wrap_dim = data_wrap.createDimension('month',12)
        depth_wrap_dim = data_wrap.createDimension('depth',z)
        lat_wrap_dim = data_wrap.createDimension('jj',jj)
        lon_wrap_dim = data_wrap.createDimension('ji',ji)
        # create 1-dimension variables
        month_wrap_var = data_wrap.createVariable('month',np.int32,('month',))
        depth_wrap_var = data_wrap.createVariable('depth',np.float32,('depth',))
        lat_wrap_var = data_wrap.createVariable('latitude_aux',np.float32,('jj',))
        # create 2-dimension variables
        gphiv_wrap_var = data_wrap.createVariable('gphiv',np.float32,('jj','ji'))
        glamv_wrap_var = data_wrap.createVariable('glamv',np.float32,('jj','ji'))
        # create 2-dimension variables
        E_0_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_eddy_steady_mean',np.float32,('month','jj'))
        E_100_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_100_eddy_steady_mean',np.float32,('month','jj'))
        E_300_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_300_eddy_steady_mean',np.float32,('month','jj'))
        E_700_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_700_eddy_steady_mean',np.float32,('month','jj'))
        E_2000_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_2000_eddy_steady_mean',np.float32,('month','jj'))
        E_pac_0_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_pac_eddy_steady_mean',np.float32,('month','jj'))
        E_pac_100_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_pac_100_eddy_steady_mean',np.float32,('month','jj'))
        E_pac_300_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_pac_300_eddy_steady_mean',np.float32,('month','jj'))
        E_pac_700_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_pac_700_eddy_steady_mean',np.float32,('month','jj'))
        E_pac_2000_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_pac_2000_eddy_steady_mean',np.float32,('month','jj'))
        E_atl_0_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_atl_eddy_steady_mean',np.float32,('month','jj'))
        E_atl_100_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_atl_100_eddy_steady_mean',np.float32,('month','jj'))
        E_atl_300_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_atl_300_eddy_steady_mean',np.float32,('month','jj'))
        E_atl_700_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_atl_700_eddy_steady_mean',np.float32,('month','jj'))
        E_atl_2000_eddy_steady_mean_wrap_var = data_wrap.createVariable('E_atl_2000_eddy_steady_mean',np.float32,('month','jj'))
        # create 3-dimension variables
        E_0_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_100_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_100_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_300_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_300_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_700_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_700_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_2000_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_2000_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_vert_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_vert_eddy_stationary_mean',np.float32,('month','depth','jj'))
        E_pac_0_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_pac_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_pac_100_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_pac_100_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_pac_300_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_pac_300_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_pac_700_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_pac_700_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_pac_2000_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_pac_2000_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_pac_vert_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_pac_vert_eddy_stationary_mean',np.float32,('month','depth','jj'))
        E_atl_0_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_atl_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_atl_100_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_atl_100_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_atl_300_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_atl_300_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_atl_700_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_atl_700_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_atl_2000_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_atl_2000_eddy_stationary_mean',np.float32,('month','jj','ji'))
        E_atl_vert_eddy_stationary_mean_wrap_var = data_wrap.createVariable('E_atl_vert_eddy_stationary_mean',np.float32,('month','depth','jj'))
        # global attributes
        data_wrap.description = 'Monthly mean eddy components of meridional energy transport fields.'
        # variable attributes
        E_0_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_100_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_300_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_700_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_2000_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_pac_0_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_pac_100_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_pac_300_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_pac_700_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_pac_2000_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_atl_0_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_atl_100_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_atl_300_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_atl_700_eddy_steady_mean_wrap_var.units = 'Tera Joule'
        E_atl_2000_eddy_steady_mean_wrap_var.units = 'Tera Joule'

        E_0_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_100_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_300_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_700_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_2000_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_vert_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_pac_0_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_pac_100_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_pac_300_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_pac_700_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_pac_2000_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_pac_vert_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_atl_0_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_atl_100_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_atl_300_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_atl_700_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_atl_2000_eddy_stationary_mean_wrap_var.units = 'Tera Joule'
        E_atl_vert_eddy_stationary_mean_wrap_var.units = 'Tera Joule'

        E_0_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow over the entire depth'
        E_100_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 100m'
        E_300_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 300m'
        E_700_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 700m'
        E_2000_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 2000m'
        E_pac_0_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow over the entire depth in the Pacific'
        E_pac_100_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 100m in the Pacific'
        E_pac_300_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 300m in the Pacific'
        E_pac_700_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 700m in the Pacific'
        E_pac_2000_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 2000m in the Pacific'
        E_atl_0_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow over the entire depth in the Atlantic'
        E_atl_100_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 100m in the Atlantic'
        E_atl_300_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 300m in the Atlantic'
        E_atl_700_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 700m in the Atlantic'
        E_atl_2000_eddy_steady_mean_wrap_var.long_name = 'Energy transport by steady mean flow from surface to 2000m in the Atlantic'

        E_0_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow over the entire depth'
        E_100_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 100m'
        E_300_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 300m'
        E_700_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 700m'
        E_2000_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 2000m'
        E_vert_eddy_stationary_mean_wrap_var.long_name = 'Vertical profile of energy transport by stationary mean flow'
        E_pac_0_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow over the entire depth in the Pacific'
        E_pac_100_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 100m in the Pacific'
        E_pac_300_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 300m in the Pacific'
        E_pac_700_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 700m in the Pacific'
        E_pac_2000_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 2000m in the Pacific'
        E_pac_vert_eddy_stationary_mean_wrap_var.long_name = 'Vertical profile of energy transport by stationary mean flow in the Pacific'
        E_atl_0_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow over the entire depth in the Atlantic'
        E_atl_100_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 100m in the Atlantic'
        E_atl_300_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 300m in the Atlantic'
        E_atl_700_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 700m in the Atlantic'
        E_atl_2000_eddy_stationary_mean_wrap_var.long_name = 'Energy transport by stationary mean flow from surface to 2000m in the Atlantic'
        E_atl_vert_eddy_stationary_mean_wrap_var.long_name = 'Vertical profile of energy transport by stationary mean flow in the Atlantic'
        # writing data
        month_wrap_var[:] = np.arange(1,13,1)
        depth_wrap_var[:] = nav_lev
        lat_wrap_var[:] = gphiv[:,96]
        gphiv_wrap_var[:] = gphiv
        glamv_wrap_var[:] = glamv

        E_0_eddy_steady_mean_wrap_var[:] = E_eddy_steady_mean
        E_100_eddy_steady_mean_wrap_var[:] = E_100_eddy_steady_mean
        E_300_eddy_steady_mean_wrap_var[:] = E_300_eddy_steady_mean
        E_700_eddy_steady_mean_wrap_var[:] = E_700_eddy_steady_mean
        E_2000_eddy_steady_mean_wrap_var[:] = E_2000_eddy_steady_mean
        E_pac_0_eddy_steady_mean_wrap_var[:] = E_pac_eddy_steady_mean
        E_pac_100_eddy_steady_mean_wrap_var[:] = E_pac_100_eddy_steady_mean
        E_pac_300_eddy_steady_mean_wrap_var[:] = E_pac_300_eddy_steady_mean
        E_pac_700_eddy_steady_mean_wrap_var[:] = E_pac_700_eddy_steady_mean
        E_pac_2000_eddy_steady_mean_wrap_var[:] = E_pac_2000_eddy_steady_mean
        E_atl_0_eddy_steady_mean_wrap_var[:] = E_atl_eddy_steady_mean
        E_atl_100_eddy_steady_mean_wrap_var[:] = E_atl_100_eddy_steady_mean
        E_atl_300_eddy_steady_mean_wrap_var[:] = E_atl_300_eddy_steady_mean
        E_atl_700_eddy_steady_mean_wrap_var[:] = E_atl_700_eddy_steady_mean
        E_atl_2000_eddy_steady_mean_wrap_var[:] = E_atl_2000_eddy_steady_mean

        E_0_eddy_stationary_mean_wrap_var[:] = E_eddy_stationary_mean
        E_100_eddy_stationary_mean_wrap_var[:] = E_100_eddy_stationary_mean
        E_300_eddy_stationary_mean_wrap_var[:] = E_300_eddy_stationary_mean
        E_700_eddy_stationary_mean_wrap_var[:] = E_700_eddy_stationary_mean
        E_2000_eddy_stationary_mean_wrap_var[:] = E_2000_eddy_stationary_mean
        E_vert_eddy_stationary_mean_wrap_var[:] = E_vert_eddy_stationary_mean
        E_pac_0_eddy_stationary_mean_wrap_var[:] = E_pac_eddy_stationary_mean
        E_pac_100_eddy_stationary_mean_wrap_var[:] = E_pac_100_eddy_stationary_mean
        E_pac_300_eddy_stationary_mean_wrap_var[:] = E_pac_300_eddy_stationary_mean
        E_pac_700_eddy_stationary_mean_wrap_var[:] = E_pac_700_eddy_stationary_mean
        E_pac_2000_eddy_stationary_mean_wrap_var[:] = E_pac_2000_eddy_stationary_mean
        E_pac_vert_eddy_stationary_mean_wrap_var[:] = E_pac_vert_eddy_stationary_mean
        E_atl_0_eddy_stationary_mean_wrap_var[:] = E_atl_eddy_stationary_mean
        E_atl_100_eddy_stationary_mean_wrap_var[:] = E_atl_100_eddy_stationary_mean
        E_atl_300_eddy_stationary_mean_wrap_var[:] = E_atl_300_eddy_stationary_mean
        E_atl_700_eddy_stationary_mean_wrap_var[:] = E_atl_700_eddy_stationary_mean
        E_atl_2000_eddy_stationary_mean_wrap_var[:] = E_atl_2000_eddy_stationary_mean
        E_atl_vert_eddy_stationary_mean_wrap_var[:] = E_atl_vert_eddy_stationary_mean
        # close the file
        data_wrap.close()
        logging.info("Create netcdf files successfully!!")
