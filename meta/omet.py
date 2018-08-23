# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Calculate Meridional Energy Transport in the Ocean with Reanalysis
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.03
Last Update     : 2018.08.08
Contributor     :
Description     : This module provides a method to perform the computation
                  of meridional energy transport in the ocean. Moreover,
                  the decomposition of standing eddy and transient eddy is also
                  included. It works with ocean reanalysis datasets. The landsea
                  mask are required.
Return Values   : netCDF files
Caveat!         :
"""

import numpy as np
from netCDF4 import Dataset
import os
import platform
import sys
import logging

def setConstants():
    """
    Define constants used in the calculations. The constants include:
    const g: gravititional acceleration     [m / s2]
    const R: radius of the earth            [m]
    const cp: heat capacity of sea water    [J/(Kg*C)]
    const rho: sea water density            [Kg/m3]

    returns: dictionary with constants for g, R. cp, rho
    rtype: dict
    """
    # define the constant:
    constant = {'g' : 9.80616,      # gravititional acceleration [m / s2]
                'R' : 6371009,      # radius of the earth [m]
                'cp': 3987,         # heat capacity of sea water [J/(Kg*C)]
                'rho': 1027,        # sea water density [Kg/m3]
                }
        
    return constant

class met:
    def __init__(self):
        """
        Quantify the meridional energy transport upto certain depth levels.
        To avoid the error due to the filled values for mask, the input data
        should not contain mask and the filled values should be 0.
        """
        print ("Start quantifying the meridional energy transport / ocean heat content in the ocean.")

    def calc_met(self, T, v, mask, e1v, e3t_0, e3t_adjust, tmaskpac, tmaskatl, 
                 z, jj, ji, z_100, z_300, z_700, z_2000, time_scale='monthly'):
        """
        Calculate the meridional energy transport in the ocean.
        param T: interpolation of potential temperature on the vector grid (v grid)    [z, jj, ji]
        param v: meridional current velocity                   [z, jj, ji]
        param mask: land-sea mask on vector grid               [z, jj, ji]
        param e1v: width of the grid box in zonal direction    [jj, ji]
        param e3t_0: depth of each grid box                    [z]
        param e3t_adjust: adjustment of grid box depth based on the definition of partial grid  [z, jj, ji]
        param tmaskpac: land-sea mask for the Pacific Ocean  [jj, ji]
        param tmaskatl: land-sea mask for the Atlantic Ocean  [jj, ji]
        param z: number of levels    [z]
        param jj: number of grids in latitudinal axis
        param ji: number of grids in longitudinal axis
        param z_100: index of depth level upto 100m
        param z_300: index of depth level upto 300m
        param z_700: index of depth level upto 700m
        param z_2000: index of depth level upto 2000m
        param time_scale: time scale of the input data (affect dimensions), two options are available
        - 'monthly' (default) the input data is monthly mean
        - 'submonthly' the input data is weekly / 5 daily / daily
        
        return: arrays of OMET upto certain depth
        rtype: numpy arrays
        """
        # generate the contants
        constant = setConstants()
        if time_scale == 'monthly':
            # calculate heat flux at each grid point
            E_flux = np.zeros((z,jj,ji),dtype=float)
            E_flux_pac = np.zeros((z,jj,ji),dtype=float)
            E_flux_atl = np.zeros((z,jj,ji),dtype=float)
            for i in np.arange(z):
                E_flux[i,:,:] = constant['rho'] * constant['cp'] * v[i,:,:] *\
                                T[i,:,:] * e1v * e3t_0[i] * mask[i,:,:] -\
                                constant['rho'] * constant['cp'] * v[i,:,:] *\
                                T[i,:,:] * e1v * e3t_adjust[i,:,:] * mask[i,:,:]
                E_flux_pac[i,:,:] = E_flux[i,:,:] * tmaskpac
                E_flux_atl[i,:,:] = E_flux[i,:,:] * tmaskatl
            # take the vertical integral upto certain layer
            # change the unit to Tera Watt
            E_flux_0 = np.sum(E_flux,0) /1e+12
            E_flux_100 = np.sum(E_flux[:z_100+1,:,:],0) /1e+12
            E_flux_300 = np.sum(E_flux[:z_300+1,:,:],0) /1e+12
            E_flux_700 = np.sum(E_flux[:z_700+1,:,:],0) /1e+12
            E_flux_2000 = np.sum(E_flux[:z_2000+1,:,:],0) /1e+12
            E_vert = np.sum(E_flux,2) /1e+12
            E_pac_vert = np.sum(E_flux_pac,2) /1e+12
            E_atl_vert = np.sum(E_flux_atl,2) /1e+12
            logging.info("The calculation of meridional energy flux is finished!")
        elif time_scale == 'submonthly':
            print ('This function will be added soon.')
        else:
            IOError("The time scale of the input fields are not supported by this module.")
        
        return E_flux_0, E_flux_100, E_flux_300, E_flux_700, E_flux_2000, E_vert, E_pac_vert, \
               E_atl_vert

    def calc_ohc(self, T, mask, e1t, e2t, e3t_0, e3t_adjust, tmaskpac, tmaskatl, 
            z, jj, ji, z_100, z_300, z_700, z_2000, time_scale='monthly'):
        """
        Quantify ocean heat content upto certain depths.
        param T: interpolation of potential temperature on the vector grid (v grid)    [z, jj, ji]
        param mask: land-sea mask on vector grid               [z, jj, ji]
        param e1t: width of the grid box in zonal direction on T grid    [jj, ji]
        param e2t: height of the grid box in zonal direction on T grid   [jj, ji]
        param e3t_0: depth of each grid box                    [z]
        param e3t_adjust: adjustment of grid box depth based on the definition of partial grid  [z, jj, ji]
        param z_100: index of depth level upto 100m
        param z_300: index of depth level upto 300m
        param z_700: index of depth level upto 700m
        param z_2000: index of depth level upto 2000m
        param time_scale: time scale of the input data (affect dimensions), two options are available
        - 'monthly' (default) the input data is monthly mean
        - 'submonthly' the input data is weekly / 5 daily / daily
        
        return: arrays of OHC upto certain depth
        rtype: numpy arrays 
        """
        # generate the contants
        constant = setConstants()
        if time_scale == 'monthly':
            # calculate heat flux at each grid point
            OHC_block = np.zeros((z,jj,ji),dtype=float)
            OHC_block_pac = np.zeros((z,jj,ji),dtype=float)
            OHC_block_atl = np.zeros((z,jj,ji),dtype=float)
            for i in np.arange(z):
                OHC_block[i,:,:] = constant['rho'] * constant['cp'] * T[i,:,:] *\
                                e1t * e2t * e3t_0[i] * mask[i,:,:] -\
                                constant['rho'] * constant['cp'] * T[i,:,:] *\
                                e1t * e2t * e3t_adjust[i,:,:] * mask[i,:,:]
                OHC_block_pac[i,:,:] = OHC_block[i,:,:] * tmaskpac
                OHC_block_atl[i,:,:] = OHC_block[i,:,:] * tmaskatl
            # take the vertical integral upto certain layer
            # change the unit to Tera Joule
            OHC_block_0 = np.sum(OHC_block,0) /1e+12
            OHC_block_100 = np.sum(OHC_block[:z_100+1,:,:],0) /1e+12
            OHC_block_300 = np.sum(OHC_block[:z_300+1,:,:],0) /1e+12
            OHC_block_700 = np.sum(OHC_block[:z_700+1,:,:],0) /1e+12
            OHC_block_2000 = np.sum(OHC_block[:z_2000+1,:,:],0) /1e+12
            OHC_vert = np.sum(OHC_block,2) /1e+12
            OHC_pac_vert = np.sum(OHC_block_pac,2) /1e+12
            OHC_atl_vert = np.sum(OHC_block_atl,2) /1e+12
            logging.info("The calculation of OHC is finished!")
        elif time_scale == 'submonthly':
            print ('This function will be added soon.')
        else:
            IOError("The time scale of the input fields are not supported by this module.")
        
        return OHC_block_0, OHC_block_100, OHC_block_300, OHC_block_700, OHC_block_2000,\
               OHC_vert, OHC_pac_vert, OHC_atl_vert
        
class eddy:
    def __init__(self):
        """
        Decompose the eddy components of AMET. The results include mean transport,
        standing eddy transport. For the ocean the data is always monthly mean, so
        it is not possible to calculate transient eddy transport directly. The transient
        eddy can be taken as the residual from overall energy transport.
        """
        print ("Start decomposing the meridional energy transport into eddy components.")

    def calc_eddy(self, T, v, mask, e1v, e3t_0, e3t_adjust, tmaskpac, tmaskatl, 
                  z, jj, ji, z_100, z_300, z_700, z_2000, time_scale='monthly'):
        """
        Decompose the eddy components of AMET.
        param T: interpolation of potential temperature on the vector grid (v grid)    [z, jj, ji]
        param v: meridional current velocity                   [z, jj, ji]
        param mask: land-sea mask on vector grid               [z, jj, ji]
        param e1v: width of the grid box in zonal direction    [jj, ji]
        param e3t_0: depth of each grid box                    [z]
        param e3t_adjust: adjustment of grid box depth based on the definition of partial grid  [z, jj, ji]
        param tmaskpac: land-sea mask for the Pacific Ocean  [jj, ji]
        param tmaskatl: land-sea mask for the Atlantic Ocean  [jj, ji]
        param z: number of levels    [z]
        param jj: number of grids in latitudinal axis
        param ji: number of grids in longitudinal axis
        param z_100: index of depth level upto 100m
        param z_300: index of depth level upto 300m
        param z_700: index of depth level upto 700m
        param z_2000: index of depth level upto 2000m
        param time_scale: time scale of the input data (affect dimensions), two options are available
        - 'monthly' (default) the input data is monthly mean
        - 'submonthly' the input data is weekly / 5 daily / daily
        
        return: arrays of OMET upto certain depth
        rtype: numpy arrays
        """
        logging.info("The calculation of standing eddy is finished!")
        # generate the contants
        constant = setConstants()
        if time_scale == 'monthly':
            # create arrays to save intermedium variables
            E_stationary_mean = np.zeros((z,jj,ji),dtype=float)
            E_stationary_mean_pac = np.zeros((z,jj,ji),dtype=float)
            E_stationary_mean_atl = np.zeros((z,jj,ji),dtype=float)
            E_steady_mean = np.zeros((z,jj),dtype=float)
            E_steady_mean_pac = np.zeros((z,jj),dtype=float)
            E_steady_mean_atl = np.zeros((z,jj),dtype=float)
            for i in np.arange(z):
                # take zonal mean for the entire globe and Pacific & Atlantic basins
                # ! the mean should be taken with the weight of gird box width
                T_zonal_mean = np.sum(T[i,:,:] * mask[i,:,:] * e1v, 1) / \
                               np.sum(mask[i,:,:] * e1v, 1)
                v_zonal_mean = np.sum(v[i,:,:] * mask[i,:,:] * e1v, 1) / \
                               np.sum(mask[i,:,:] * e1v, 1)
                T_zonal_mean_pac = np.sum(T[i,:,:] * tmaskpac * mask[i,:,:] * e1v, 1) / \
                                   np.sum(tmaskpac * mask[i,:,:] * e1v, 1)
                v_zonal_mean_pac = np.sum(v[i,:,:] * tmaskpac * mask[i,:,:] * e1v, 1) / \
                                   np.sum(tmaskpac * mask[i,:,:] * e1v, 1)
                T_zonal_mean_atl = np.sum(T[i,:,:] * tmaskatl * mask[i,:,:] * e1v, 1) / \
                                   np.sum(tmaskatl * mask[i,:,:] * e1v, 1)
                v_zonal_mean_atl = np.sum(v[i,:,:] * tmaskatl * mask[i,:,:] * e1v, 1) / \
                                   np.sum(tmaskatl * mask[i,:,:] * e1v, 1)
                ##############  steady mean eddy ###############
                E_steady_mean[i,:] = constant['rho'] * constant['cp'] * T_zonal_mean * \
                                     v_zonal_mean * (np.sum(e1v,1) / np.sum(mask[i,:,:],1)) * e3t_0[i]
                E_steady_mean_pac[i,:] = constant['rho'] * constant['cp'] * T_zonal_mean_pac * \
                                         v_zonal_mean_pac * e3t_0[i] *\
                                         (np.sum(e1v,1) / np.sum(mask[i,:,:] * tmaskpac,1))
                E_steady_mean_atl[i,:] = constant['rho'] * constant['cp'] * T_zonal_mean_atl * \
                                         v_zonal_mean_atl * e3t_0[i] *\
                                         (np.sum(e1v,1) / np.sum(mask[i,:,:] * tmaskatl,1))                                                      
                ##############  stationary mean eddy ###############
                T_zonal_mean_enlarge = np.repeat(T_zonal_mean_pac[:,np.newaxis], ji, 1)
                v_zonal_mean_enlarge = np.repeat(v_zonal_mean_pac[:,np.newaxis], ji, 1)
                # globe
                E_stationary_mean[i,:,:] = constant['rho'] * constant['cp'] * (v[i,:,:] -
                                           v_zonal_mean_enlarge) * (T[i,:,:] -
                                           T_zonal_mean_enlarge) * e1v * e3t_0[i] *\
                                           mask[i,:,:] - constant['rho'] * constant['cp'] *\
                                           (v[i,:,:] - v_zonal_mean_enlarge) * (T[i,:,:] -
                                           T_zonal_mean_enlarge) * e1v * e3t_adjust[i,:,:] *\
                                           mask[i,:,:]
                # pacific
                T_zonal_mean_pac_enlarge = np.repeat(T_zonal_mean_pac[:,np.newaxis], ji, 1)
                v_zonal_mean_pac_enlarge = np.repeat(v_zonal_mean_pac[:,np.newaxis], ji, 1)
                E_stationary_mean_pac[i,:,:] = constant['rho'] * constant['cp'] * (v[i,:,:] -
                                           v_zonal_mean_pac_enlarge) * (T[i,:,:] -
                                           T_zonal_mean_pac_enlarge) * e1v * e3t_0[i] *\
                                           mask[i,:,:] * tmaskpac - constant['rho'] * constant['cp'] *\
                                           (v[i,:,:] - v_zonal_mean_pac_enlarge) * (T[i,:,:] -
                                           T_zonal_mean_pac_enlarge) * e1v * e3t_adjust[i,:,:] *\
                                           mask[i,:,:] * tmaskpac
                # atlantic
                T_zonal_mean_atl_enlarge = np.repeat(T_zonal_mean_pac[:,np.newaxis], ji, 1)
                v_zonal_mean_atl_enlarge = np.repeat(v_zonal_mean_pac[:,np.newaxis], ji, 1)
                E_stationary_mean_atl[i,:,:] = constant['rho'] * constant['cp'] * (v[i,:,:] -
                                           v_zonal_mean_atl_enlarge) * (T[i,:,:] -
                                           T_zonal_mean_atl_enlarge) * e1v * e3t_0[i] *\
                                           mask[i,:,:] * tmaskatl - constant['rho'] * constant['cp'] *\
                                           (v[i,:,:] - v_zonal_mean_atl_enlarge) * (T[i,:,:] -
                                           T_zonal_mean_atl_enlarge) * e1v * e3t_adjust[i,:,:] *\
                                           mask[i,:,:] * tmaskatl
            # take the vertical integral upto certain layer
            ##############  steady mean eddy ###############
            # change the unit to Tera Watt
            E_steady_mean_0 = np.sum(E_steady_mean,0) / 1e+12
            E_steady_mean_100 = np.sum(E_steady_mean[:z_100+1,:],0) / 1e+12
            E_steady_mean_300 = np.sum(E_steady_mean[:z_300+1,:],0) / 1e+12
            E_steady_mean_700 = np.sum(E_steady_mean[:z_700+1,:],0) / 1e+12
            E_steady_mean_2000 = np.sum(E_steady_mean[:z_2000+1,:],0) / 1e+12
            E_steady_mean_pac_0 = np.sum(E_steady_mean_pac,0) / 1e+12
            E_steady_mean_pac_100 = np.sum(E_steady_mean_pac[:z_100+1,:],0) / 1e+12
            E_steady_mean_pac_300 = np.sum(E_steady_mean_pac[:z_300+1,:],0) / 1e+12
            E_steady_mean_pac_700 = np.sum(E_steady_mean_pac[:z_700+1,:],0) / 1e+12
            E_steady_mean_pac_2000 = np.sum(E_steady_mean_pac[:z_2000+1,:],0) / 1e+12
            E_steady_mean_atl_0 = np.sum(E_steady_mean_atl,0) / 1e+12
            E_steady_mean_atl_100 = np.sum(E_steady_mean_atl[:z_100+1,:],0) / 1e+12
            E_steady_mean_atl_300 = np.sum(E_steady_mean_atl[:z_300+1,:],0) / 1e+12
            E_steady_mean_atl_700 = np.sum(E_steady_mean_atl[:z_700+1,:],0) / 1e+12
            E_steady_mean_atl_2000 = np.sum(E_steady_mean_atl[:z_2000+1,:],0) / 1e+12
            ##############  stationary mean eddy ###############
            # change the unit to Tera Watt
            E_stationary_mean_0 = np.sum(E_stationary_mean,0) / 1e+12
            E_stationary_mean_100 = np.sum(E_stationary_mean[:z_100+1,:,:],0) / 1e+12
            E_stationary_mean_300 = np.sum(E_stationary_mean[:z_300+1,:,:],0) / 1e+12
            E_stationary_mean_700 = np.sum(E_stationary_mean[:z_700+1,:,:],0) / 1e+12
            E_stationary_mean_2000 = np.sum(E_stationary_mean[:z_2000+1,:,:],0) / 1e+12
            E_stationary_mean_vert = np.sum(E_stationary_mean,2) / 1e+12
            E_stationary_mean_pac_0 = np.sum(E_stationary_mean_pac,0) / 1e+12
            E_stationary_mean_pac_100 = np.sum(E_stationary_mean_pac[:z_100+1,:,:],0) / 1e+12
            E_stationary_mean_pac_300 = np.sum(E_stationary_mean_pac[:z_300+1,:,:],0) / 1e+12
            E_stationary_mean_pac_700 = np.sum(E_stationary_mean_pac[:z_700+1,:,:],0) / 1e+12
            E_stationary_mean_pac_2000 = np.sum(E_stationary_mean_pac[:z_2000+1,:,:],0) / 1e+12
            E_stationary_mean_pac_vert = np.sum(E_stationary_mean_pac,2) / 1e+12
            E_stationary_mean_atl_0 = np.sum(E_stationary_mean_atl,0) / 1e+12
            E_stationary_mean_atl_100 = np.sum(E_stationary_mean_atl[:z_100+1,:,:],0) / 1e+12
            E_stationary_mean_atl_300 = np.sum(E_stationary_mean_atl[:z_300+1,:,:],0) / 1e+12
            E_stationary_mean_atl_700 = np.sum(E_stationary_mean_atl[:z_700+1,:,:],0) / 1e+12
            E_stationary_mean_atl_2000 = np.sum(E_stationary_mean_atl[:z_2000+1,:,:],0) / 1e+12
            E_stationary_mean_atl_vert = np.sum(E_stationary_mean_atl,2) / 1e+12
            logging.info("The calculation of meridional energy flux is finished!")
        elif time_scale == 'submonthly':
            print ('This function will be added soon.')
        else:
            IOError("The time scale of the input fields are not supported by this module.")
            
        return E_steady_mean_0, E_steady_mean_100, E_steady_mean_300, E_steady_mean_700,\
               E_steady_mean_2000, E_steady_mean_pac_0, E_steady_mean_pac_100, E_steady_mean_pac_300,\
               E_steady_mean_pac_700, E_steady_mean_pac_2000, E_steady_mean_atl_0, E_steady_mean_atl_100, \
               E_steady_mean_atl_300, E_steady_mean_atl_700, E_steady_mean_atl_2000, E_stationary_mean_0, \
               E_stationary_mean_100, E_stationary_mean_300, E_stationary_mean_700, E_stationary_mean_2000,\
               E_stationary_mean_vert, E_stationary_mean_pac_0, E_stationary_mean_pac_100, \
               E_stationary_mean_pac_300, E_stationary_mean_pac_700, E_stationary_mean_pac_2000,\
               E_stationary_mean_pac_vert, E_stationary_mean_atl_0, E_stationary_mean_atl_100,\
               E_stationary_mean_atl_300, E_stationary_mean_atl_700, E_stationary_mean_atl_2000,\
               E_stationary_mean_atl_vert
    