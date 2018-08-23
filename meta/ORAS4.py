# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from ORAS4
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.09
Last Update     : 2018.08.09
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of ECMWF. It provides an
                  entrance for the following computation includes the mass budget
                  correction, quantification of meridional energy transport, decomposition
                  of eddies.
                  
                  ORAS4 is a state-of-the-art ocean reanalysis product produced by ECMWF.
                  It spans from 1958 to 2016. Natively it is generated on ORCA1 grid and 42 
                  vertical levels.
                  
                  The processing unit is monthly data, for the sake of memory saving.
                  
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is 
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /ORAS4
                      /theta
                          /thetao_oras4_1m_1958_grid_T.nc
                          /thetao_oras4_1m_1959_grid_T.nc
                          ...
                          /thetao_oras4_1m_2016_grid_T.nc
                      /v
                          /vo_oras4_1m_1958_grid_V.nc
                          /vo_oras4_1m_1959_grid_V.nc
                          ...
                          /vo_oras4_1m_2016_grid_V.nc
                      /u
                      /s
                          
                  The name rule follows the default name provided by the uni-Hamburg FTP
                  when downloading the files.
                  
                  T /v fields from ORAS4
                  masked array? [yes]
                  containing mask info? [no]
                  mask from ORAS4
                  Boolean? [no] -> sea-1 land-0  
"""

import sys
import os
import numpy as np
import logging
from netCDF4 import Dataset
import meta.omet
#import matplotlib
import meta.saveNetCDF

class oras:
    def __init__(self, path, out_path):
        """
        Initialize the extraction of fields from ORAS4.
        
        The data is on original ORCA1 grid. As the interpolation can introduce
        large errors to the computation of energy transport, we will follow the
        model level. The determination of depth level is based on the standard
        depth info given by the bathemetry 
        
        param path: the root path of the input fields
        param out_path: the location of output files
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        param z_100: index of depth level upto 100m
        param z_300: index of depth level upto 300m
        param z_700: index of depth level upto 700m
        param z_2000: index of depth level upto 2000m
        """
        self.path = path
        self.out_path = out_path
        # number of levels for certain depth - from sea surface to certain depth
        # following the index rule of python
        # due to the vertical resolution, we choose the nearest layer but no deeper
        self.z_100 = 9
        self.z_300 = 18
        self.z_700 = 23
        self.z_2000 = 29
    
    @staticmethod
    def setConstants():
        """
        Define constants used in the calculations. The constants include:
        const g: gravititional acceleration     [m / s2]
        const R: radius of the earth            [m]
        const cp: heat capacity of sea water    [J/(Kg*C)]
        const rho: sea water density            [Kg/m3]
        """
        # define the constant:
        constant ={'g' : 9.80616,      # gravititional acceleration [m / s2]
                   'R' : 6371009,      # radius of the earth [m]
                   'cp': 3987,         # heat capacity of sea water [J/(Kg*C)]
                   'rho': 1027,        # sea water density [Kg/m3]
                   }
        return constant
    
    def omet(self, year_start, year_end, path_mask, path_submask):
        """
        Quantify Meridional Energy Transport.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param path_mask: location of the file containing land sea mask and coordinate info
        param path_submask: land-sea mask of the sub-basin
        
        return: arrays containing OMET upto differnt depth levels
        rtype: netCDF4
        """
        # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.out_path,'history_omet.log') ,
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        year = np.arange(year_start, year_end+1, 1)
        month = np.arange(1, 13, 1)
        #month_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12' ]
        # extract information from land-sea mask file
        nav_lat, nav_lon, nav_lev, tmask, vmask, e1t, e2t, e1v, e2v,\
        gphiv, glamv, mbathy, e3t_0, e3t_ps = self.mask(path_mask)
        logging.info('Successfully retrieving the ORCA1 coordinate and mask info!')
        # obtain the land-sea mask of sub-ocean basins
        subbasin_key = Dataset(path_submask)
        tmaskpac = subbasin_key.variables['tmaskpac_ORCA1'][:]
        tmaskatl = subbasin_key.variables['tmaskatl_ORCA1'][:]
        # obtain dimensions
        jj, ji = nav_lat.shape
        z = len(nav_lev)
        # construct partial cell depth matrix
        # the size of partial cell is given by e3t_ps
        # for the sake of simplicity of the code, just calculate the difference between e3t_0 and e3t_ps
        # then minus this adjustment when calculate the OMET at each layer with mask
        # Attention! Since python start with 0, the partial cell info given in mbathy should incoporate with this
        e3t_adjust = np.zeros((z,jj,ji),dtype = float)
        for i in np.arange(1,z,1): # start from 1
            for j in np.arange(jj):
                for k in np.arange(ji):
                    if i == mbathy[j,k]:
                        e3t_adjust[i-1,j,k] = e3t_0[i-1] - e3t_ps[j,k] # python start with 0, so i-1
        # create arrays to store the output data
        E = np.zeros((len(month), jj, ji), dtype=float)
        E_100 = np.zeros((len(month), jj, ji), dtype=float)
        E_300 = np.zeros((len(month), jj, ji), dtype=float)
        E_700 = np.zeros((len(month), jj, ji), dtype=float)
        E_2000 = np.zeros((len(month), jj, ji), dtype=float)
        E_pac_vert = np.zeros((len(month), z, jj), dtype=float)
        E_atl_vert = np.zeros((len(month), z, jj), dtype=float)
        E_vert = np.zeros((len(month), z, jj), dtype=float)
        # loop for the computation of OMET
        for i in year:
            logging.info("Start retrieving variables for {}(y)".format(i))
            datapath_theta = os.path.join(self.path, 'theta', 'thetao_oras4_1m_{}_grid_T.nc'.format(i))
            datapath_v = os.path.join(self.path, 'v', 'vo_oras4_1m_{}_grid_V.nc'.format(i))
            theta_key = Dataset(datapath_theta)
            v_key = Dataset(datapath_v)
            theta = theta_key.variables['thetao'][:] # the unit of theta is Celsius!
            theta_nomask = theta.data
            theta_nomask[theta.mask == True] = 0
            v = v_key.variables['vo'][:]
            v_nomask = v.data
            v_nomask[v.mask == True] = 0
            logging.info("Extracting variables successfully!")
            # calculate the meridional velocity at T grid
            T_vgrid = np.zeros((len(month),z,jj,ji),dtype=float)
            for c in np.arange(jj):
                if c == jj-1:
                    T_vgrid[:,:,c,:] = theta_nomask[:,:,c,:]
                else:
                    T_vgrid[:,:,c,:] = (theta_nomask[:,:,c,:] + theta_nomask[:,:,c+1,:])/2
            # as the module omet is designed to process data per month, we add a loop for month
            for j in month:
                OMET = meta.omet.met()
                E[j-1,:,:], E_100[j-1,:,:], E_300[j-1,:,:], E_700[j-1,:,:], \
                E_2000[j-1,:,:], E_vert[j-1,:,:], E_pac_vert[j-1,:,:],\
                E_atl_vert[j-1,:,:] = OMET.calc_met(T_vgrid[j-1,:,:,:], v_nomask[j-1,:,:,:],
                                                    vmask, e1v, e3t_0, e3t_adjust, tmaskpac, 
                                                    tmaskatl, z, jj, ji, self.z_100,self.z_300,
                                                    self.z_700, self.z_2000)
            # save output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncOMET(E, E_100, E_300, E_700, E_2000, E_vert, E_pac_vert,
                           E_atl_vert, i, gphiv, glamv, nav_lev, z, jj, ji, self.out_path)
        logging.info('The entire pipeline is complete!')

    def ohc(self, year_start, year_end, path_mask, path_submask):
        """
        Quantify ocean heat content upto certain depths.
        param 
        
        return:
        rtype: netCDF4
        """
        # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.out_path,'history_ohc.log') ,
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        year = np.arange(year_start, year_end+1, 1)
        month = np.arange(1, 13, 1)
        #month_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12' ]
        # extract information from land-sea mask file
        nav_lat, nav_lon, nav_lev, tmask, vmask, e1t, e2t, e1v, e2v,\
        gphiv, glamv, mbathy, e3t_0, e3t_ps = self.mask(path_mask)
        logging.info('Successfully retrieving the ORCA1 coordinate and mask info!')
        # obtain the land-sea mask of sub-ocean basins
        subbasin_key = Dataset(path_submask)
        tmaskpac = subbasin_key.variables['tmaskpac_ORCA1'][:]
        tmaskatl = subbasin_key.variables['tmaskatl_ORCA1'][:]
        # obtain dimensions
        jj, ji = nav_lat.shape
        z = len(nav_lev)
        # construct partial cell depth matrix
        # the size of partial cell is given by e3t_ps
        # for the sake of simplicity of the code, just calculate the difference between e3t_0 and e3t_ps
        # then minus this adjustment when calculate the OMET at each layer with mask
        # Attention! Since python start with 0, the partial cell info given in mbathy should incoporate with this
        e3t_adjust = np.zeros((z,jj,ji),dtype = float)
        for i in np.arange(1,z,1): # start from 1
            for j in np.arange(jj):
                for k in np.arange(ji):
                    if i == mbathy[j,k]:
                        e3t_adjust[i-1,j,k] = e3t_0[i-1] - e3t_ps[j,k] # python start with 0, so i-1
        # create arrays to store the output data
        OHC = np.zeros((len(month), jj, ji), dtype=float)
        OHC_100 = np.zeros((len(month), jj, ji), dtype=float)
        OHC_300 = np.zeros((len(month), jj, ji), dtype=float)
        OHC_700 = np.zeros((len(month), jj, ji), dtype=float)
        OHC_2000 = np.zeros((len(month), jj, ji), dtype=float)
        OHC_pac_vert = np.zeros((len(month), z, jj), dtype=float)
        OHC_atl_vert = np.zeros((len(month), z, jj), dtype=float)
        OHC_vert = np.zeros((len(month), z, jj), dtype=float)
        # loop for the computation of OHC
        for i in year:
            logging.info("Start retrieving variables for {}(y)".format(i))
            datapath_theta = os.path.join(self.path, 'theta', 'thetao_oras4_1m_{}_grid_T.nc'.format(i))
            theta_key = Dataset(datapath_theta)
            theta = theta_key.variables['thetao'][:] # the unit of theta is Celsius!
            theta_nomask = theta.data
            theta_nomask[theta.mask == True] = 0
            logging.info("Extracting variables successfully!")
            # as the module omet is designed to process data per month, we add a loop for month
            for j in month:
                OMET = meta.omet.met()
                OHC[j-1,:,:], OHC_100[j-1,:,:], OHC_300[j-1,:,:], OHC_700[j-1,:,:],\
                OHC_2000[j-1,:,:], OHC_vert[j-1,:,:], OHC_pac_vert[j-1,:,:],\
                OHC_atl_vert[j-1,:,:] = OMET.calc_ohc(theta_nomask[j-1,:,:,:], tmask, e1t, e2t,
                                                      e3t_0, e3t_adjust, tmaskpac, tmaskatl, z,
                                                      jj, ji, self.z_100,self.z_300, self.z_700,
                                                      self.z_2000)
            # save output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncOHC(OHC, OHC_100, OHC_300, OHC_700, OHC_2000, OHC_vert, OHC_pac_vert,
                           OHC_atl_vert, i, nav_lat, nav_lon, nav_lev, gphiv, z, jj, ji, self.out_path)                
        logging.info('The entire pipeline is complete!')
        
    def psi(self, year_start, year_end, path_mask, path_submask):
        """
        Quantify the mass transport upto certain depths.
        param path_mask: location of the file containing land sea mask and coordinate info
        
        return:
        rtype: netCDF4        
        """
        # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.out_path,'history_psi.log') ,
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        year = np.arange(year_start, year_end+1, 1)
        month = np.arange(1, 13, 1)
        #month_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12' ]
        # extract information from land-sea mask file
        nav_lat, nav_lon, nav_lev, tmask, vmask, e1t, e2t, e1v, e2v,\
        gphiv, glamv, mbathy, e3t_0, e3t_ps = self.mask(path_mask)
        logging.info('Successfully retrieving the ORCA1 coordinate and mask info!')
        # obtain the land-sea mask of sub-ocean basins
        subbasin_key = Dataset(path_submask)
        tmaskpac = subbasin_key.variables['tmaskpac_ORCA1'][:]
        tmaskatl = subbasin_key.variables['tmaskatl_ORCA1'][:]
        # obtain dimensions
        jj, ji = nav_lat.shape
        z = len(nav_lev)
        # construct partial cell depth matrix
        # the size of partial cell is given by e3t_ps
        # for the sake of simplicity of the code, just calculate the difference between e3t_0 and e3t_ps
        # then minus this adjustment when calculate the OMET at each layer with mask
        # Attention! Since python start with 0, the partial cell info given in mbathy should incoporate with this
        e3t_adjust = np.zeros((z,jj,ji),dtype = float)
        for i in np.arange(1,z,1): # start from 1
            for j in np.arange(jj):
                for k in np.arange(ji):
                    if i == mbathy[j,k]:
                        e3t_adjust[i-1,j,k] = e3t_0[i-1] - e3t_ps[j,k] # python start with 0, so i-1
        logging.info('The entire pipeline is complete!')
        
    def mask(self, path_mask):
        """
        Extract information from land-sea mask.
        param path_mask: location of the file containing land sea mask and coordinate info
        
        return: numpy arrays with coordinate and land-sea mask info.
        rtype: numpy arrays
        """
        logging.info('Start retrieving the ORCA1 coordinate and mask info.')
        mesh_mask_key = Dataset(path_mask)
        # lat-lon-depth coordinate info
        nav_lat = mesh_mask_key.variables['nav_lat'][:]
        nav_lon = mesh_mask_key.variables['nav_lon'][:]
        nav_lev = mesh_mask_key.variables['nav_lev'][:]
        # lat-lon coordinate of V grid
        gphiv = mesh_mask_key.variables['gphiv'][0,:,:] # lat from -78 to -89
        glamv = mesh_mask_key.variables['glamv'][0,:,:] # lon from -179 to 179
        # land-sea mask
        tmask = mesh_mask_key.variables['tmask'][0,:,:,:]
        #umask = mesh_mask_key.variables['umask'][0,:,:,:]
        vmask = mesh_mask_key.variables['vmask'][0,:,:,:]
        # grid spacing scale factors (zonal)
        e1t = mesh_mask_key.variables['e1t'][0,:,:]
        e2t = mesh_mask_key.variables['e2t'][0,:,:]
        #e1u = mesh_mask_key.variables['e1u'][0,:,:]
        #e2u = mesh_mask_key.variables['e2u'][0,:,:]
        e1v = mesh_mask_key.variables['e1v'][0,:,:]
        e2v = mesh_mask_key.variables['e2v'][0,:,:]
        # take the bathymetry
        mbathy = mesh_mask_key.variables['mbathy'][0,:,:]
        # depth of each layer
        e3t_0 = mesh_mask_key.variables['e3t_0'][0,:]
        e3t_ps = mesh_mask_key.variables['e3t_ps'][0,:,:]
        
        return nav_lat, nav_lon, nav_lev, tmask, vmask, e1t, e2t,\
               e1v, e2v, gphiv, glamv, mbathy, e3t_0, e3t_ps
    
    def eddies(self, year_start, year_end, path_mask, path_submask):
        """
        Decompose the eddy contributions for the Meridional Energy Transport.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param path_mask: location of the file containing land sea mask and coordinate info
        param path_submask: land-sea mask of the sub-basin
        
        return: arrays containing OMET upto differnt depth levels
        rtype: netCDF4
        """
        # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.out_path,'history_eddy.log') ,
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        year = np.arange(year_start, year_end+1, 1)
        month = np.arange(1, 13, 1)
        #month_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12' ]
        # extract information from land-sea mask file
        nav_lat, nav_lon, nav_lev, tmask, vmask, e1t, e2t, e1v, e2v,\
        gphiv, glamv, mbathy, e3t_0, e3t_ps = self.mask(path_mask)
        logging.info('Successfully retrieving the ORCA1 coordinate and mask info!')
        # obtain the land-sea mask of sub-ocean basins
        subbasin_key = Dataset(path_submask)
        tmaskpac = subbasin_key.variables['tmaskpac_ORCA1'][:]
        tmaskatl = subbasin_key.variables['tmaskatl_ORCA1'][:]
        # obtain dimensions
        jj, ji = nav_lat.shape
        z = len(nav_lev)
        # construct partial cell depth matrix
        # the size of partial cell is given by e3t_ps
        # for the sake of simplicity of the code, just calculate the difference between e3t_0 and e3t_ps
        # then minus this adjustment when calculate the OMET at each layer with mask
        # Attention! Since python start with 0, the partial cell info given in mbathy should incoporate with this
        e3t_adjust = np.zeros((z,jj,ji),dtype = float)
        for i in np.arange(1,z,1): # start from 1
            for j in np.arange(jj):
                for k in np.arange(ji):
                    if i == mbathy[j,k]:
                        e3t_adjust[i-1,j,k] = e3t_0[i-1] - e3t_ps[j,k] # python start with 0, so i-1
        # create arrays to store the output data
        # steady mean
        E_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_100_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_300_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_700_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_2000_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_pac_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_pac_100_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_pac_300_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_pac_700_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_pac_2000_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_atl_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_atl_100_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_atl_300_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_atl_700_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        E_atl_2000_eddy_steady_mean = np.zeros((len(month), jj), dtype=float)
        # stationary mean
        E_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_100_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_300_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_700_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_2000_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_vert_eddy_stationary_mean = np.zeros((len(month), z, jj), dtype=float)
        E_pac_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_pac_100_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_pac_300_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_pac_700_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_pac_2000_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)        
        E_pac_vert_eddy_stationary_mean = np.zeros((len(month), z, jj), dtype=float)
        E_atl_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_atl_100_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_atl_300_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_atl_700_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)
        E_atl_2000_eddy_stationary_mean = np.zeros((len(month), jj, ji), dtype=float)        
        E_atl_vert_eddy_stationary_mean = np.zeros((len(month), z, jj), dtype=float)
        # loop for the computation of OMET
        for i in year:
            logging.info("Start retrieving variables for {}(y)".format(i))
            datapath_theta = os.path.join(self.path, 'theta', 'thetao_oras4_1m_{}_grid_T.nc'.format(i))
            datapath_v = os.path.join(self.path, 'v', 'vo_oras4_1m_{}_grid_V.nc'.format(i))
            theta_key = Dataset(datapath_theta)
            v_key = Dataset(datapath_v)
            theta = theta_key.variables['thetao'][:] # the unit of theta is Celsius!
            theta_nomask = theta.data
            theta_nomask[theta.mask == True] = 0
            v = v_key.variables['vo'][:]
            v_nomask = v.data
            v_nomask[v.mask == True] = 0
            logging.info("Extracting variables successfully!")
            # calculate the meridional velocity at T grid
            T_vgrid = np.zeros((len(month),z,jj,ji),dtype=float)
            for c in np.arange(jj):
                if c == jj-1:
                    T_vgrid[:,:,c,:] = theta_nomask[:,:,c,:]
                else:
                    T_vgrid[:,:,c,:] = (theta_nomask[:,:,c,:] + theta_nomask[:,:,c+1,:])/2
            # as the module omet is designed to process data per month, we add a loop for month
            for j in month:
                OMET = meta.omet.eddy()
                E_eddy_steady_mean[j-1,:], E_100_eddy_steady_mean[j-1,:], E_300_eddy_steady_mean[j-1,:],\
                E_700_eddy_steady_mean[j-1,:], E_2000_eddy_steady_mean[j-1,:], E_pac_eddy_steady_mean[j-1,:],\
                E_pac_100_eddy_steady_mean[j-1,:], E_pac_300_eddy_steady_mean[j-1,:], E_pac_700_eddy_steady_mean[j-1,:],\
                E_pac_2000_eddy_steady_mean[j-1,:], E_atl_eddy_steady_mean[j-1,:], E_atl_100_eddy_steady_mean[j-1,:],\
                E_atl_300_eddy_steady_mean[j-1,:], E_atl_700_eddy_steady_mean[j-1,:], E_atl_2000_eddy_steady_mean[j-1,:],\
                E_eddy_stationary_mean[j-1,:,:], E_100_eddy_stationary_mean[j-1,:,:], E_300_eddy_stationary_mean[j-1,:,:],\
                E_700_eddy_stationary_mean[j-1,:,:], E_2000_eddy_stationary_mean[j-1,:,:], E_vert_eddy_stationary_mean[j-1,:,:],\
                E_pac_eddy_stationary_mean[j-1,:,:], E_pac_100_eddy_stationary_mean[j-1,:,:],\
                E_pac_300_eddy_stationary_mean[j-1,:,:], E_pac_700_eddy_stationary_mean[j-1,:,:],\
                E_pac_2000_eddy_stationary_mean[j-1,:,:], E_pac_vert_eddy_stationary_mean[j-1,:,:],\
                E_atl_eddy_stationary_mean[j-1,:,:], E_atl_100_eddy_stationary_mean[j-1,:,:],\
                E_atl_300_eddy_stationary_mean[j-1,:,:], E_atl_700_eddy_stationary_mean[j-1,:,:],\
                E_atl_2000_eddy_stationary_mean[j-1,:,:], E_atl_vert_eddy_stationary_mean[j-1,:,:]\
                = OMET.calc_eddy(T_vgrid[j-1,:,:,:], v_nomask[j-1,:,:,:], vmask, e1v, e3t_0,
                                e3t_adjust, tmaskpac, tmaskatl, z, jj, ji, self.z_100,self.z_300,
                                self.z_700, self.z_2000)
            # save output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncEddyomet(E_eddy_steady_mean, E_100_eddy_steady_mean, E_300_eddy_steady_mean,
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
                               i, gphiv, glamv, nav_lev, z, jj, ji, self.out_path) 
        logging.info('The entire pipeline is complete!')
                