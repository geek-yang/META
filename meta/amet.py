# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Calculate Meridional Energy Transport in the Atmosphere with Reanalysis
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.01
Last Update     : 2018.08.08
Contributor     :
Description     : This module provides a method to perform the computation
                  of meridional energy transport in the atmosphere. Moreover,
                  the decomposition of standing eddy and transient eddy is also
                  included. It works with atmosphere reanalysis datasets. In
                  order to obtain meaningful values, the mass budget correction
                  must be carried out before hand. So, this module requires input
                  of baratropic wind correction terms from the mass budget
                  conservation.
Return Values   : netCDF files
Caveat!         :
"""

import numpy as np
from netCDF4 import Dataset
import os
import platform
import sys
import logging
#import matplotlib
# generate images without having a window appear
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

def setConstants():
    '''
    Define constants used in the calculations. The constants include:
    const g: gravititional acceleration [m / s2]
    const R: radius of the earth [m]
    const cp: heat capacity of air [J/(Kg*K)]
    const Lv: Latent heat of vaporization [J/Kg]
    const R_dry: gas constant of dry air [J/(kg*K)]
    const R_vap: gas constant for water vapour [J/(kg*K)]

    returns: dictionary with constants for g, R. cp, Lv, R_dry and R_vap
    rtype: dict
    '''
    # define the constant:
    constant = {'g' : 9.80616,      # gravititional acceleration [m / s2]
                'R' : 6371009,      # radius of the earth [m]
                'cp': 1004.64,      # heat capacity of air [J/(Kg*K)]
                'Lv': 2264670,      # Latent heat of vaporization [J/Kg]
                'R_dry' : 286.9,    # gas constant of dry air [J/(kg*K)]
                'R_vap' : 461.5,    # gas constant for water vapour [J/(kg*K)]
                }

    return constant

class met:
    def __init__(self):
        """
        Quantify the meridional energy transport and its components upto
        certain pressure levels.
        """
        print ("Start quantifying the meridional energy transport in the atmosphere.")

    def calc_met(self, T, q, sp, u, v, gz, A, B, t, h, y, x,
                 lat, lat_unit, vc):
        """
        Calculate the meridional energy transport and its components in the atmosphere.
        All the input files should contain the fields for the entire month.
        Caveat! The integral is taken from TOA to the Surface!
        param q: Specific Humidity [kg/kg]
        param sp: Surface Pressure [Pa]
        param u: Zonal Wind [m/s]
        param v: Meridional Wind [m/s]
        param gz: Geopotential [m2/s2]
        param A: Constant A for Defining Sigma Level
        param B: Constant B for Defining Sigma Level
        param t: time dimension of input fields
        param h: level dimension of input fields
        param y: latitude dimension of input fields
        param x: longitude dimension of input fields
        param lat: latitude
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        param vc: Baratropic Corrected Meridional Winds

        returns: AMET and its components (internal, latent, geopotential, kinetic energy)
                 upto different heights.
        rtype: numpy array
        """
        # generate the contants
        constant = setConstants()
        # calculate dp
        dp_level = np.zeros((t, h, y, x),dtype = float)
        for i in np.arange(h):
            dp_level[:,i,:,:] = (A[i+1] + B[i+1] * sp) - (A[i] + B[i] * sp)
        # calculate each component of total energy and take the vertical integral
        # Internal Energy cpT
        internal_energy = constant['cp'] * (1-q) * T * dp_level / constant['g']
        del T
        internal_flux_int = np.mean(np.sum(internal_energy * v,1),0)
        internal_flux_int_correct = np.mean(np.sum(internal_energy,1),0) * vc
        del internal_energy
        logging.info("The calculation of internal energy flux is finished!")
        # Latent heat Lvq
        latent_energy = constant['Lv'] * q * dp_level / constant['g']
        del q
        latent_flux_int = np.mean(np.sum(latent_energy * v,1),0)
        latent_flux_int_correct = np.mean(np.sum(latent_energy,1),0) * vc
        del latent_energy
        logging.info("The calculation of latent heat flux is finished!")
        # Geopotential Energy gz
        geopotential_energy = gz * dp_level / constant['g']
        geopotential_flux_int = np.mean(np.sum(geopotential_energy * v,1),0)
        del gz
        geopotential_flux_int_correct = np.mean(np.sum(geopotential_energy,1),0) * vc
        del geopotential_energy
        logging.info("The calculation of geopotential energy flux is finished!")
        # Kinetic Energy u2+v2
        kinetic_energy = 0.5 * (u**2 + v**2) * dp_level / constant['g']
        del u
        kinetic_flux_int = np.mean(np.sum(kinetic_energy * v,1),0)
        kinetic_flux_int_correct = np.mean(np.sum(kinetic_energy,1),0) * vc
        del kinetic_energy, v
        logging.info("The calculation of kinetic energy flux is finished!")
        # the earth is taken as a perfect sphere, instead of a ellopsoid
        dx = 2 * np.pi * constant['R'] * np.cos(2 * np.pi * lat / 360) / x
        # plugin the weight of grid box width and apply the correction
        #dx[0] = 0
        # create arrays for each energy transport components after correction
        E_internal = np.zeros((y, x),dtype=float)
        E_latent = np.zeros((y, x),dtype=float)
        E_geopotential = np.zeros((y, x),dtype=float)
        E_kinetic = np.zeros((y, x),dtype=float)
        E_total = np.zeros((y, x),dtype=float)

        E_internal_correct = np.zeros((y, x),dtype=float)
        E_latent_correct = np.zeros((y, x),dtype=float)
        E_geopotential_correct = np.zeros((y, x),dtype=float)
        E_kinetic_correct = np.zeros((y, x),dtype=float)
        E_total_correct = np.zeros((y, x),dtype=float)
        # apply correction and weight of grid width
        # also change the unit from Watt to Tera Watt
        for i in np.arange(y):
            E_internal[i,:] = (internal_flux_int[i,:] - internal_flux_int_correct[i,:]) * dx[i]/1e+12
            E_latent[i,:] = (latent_flux_int[i,:] - latent_flux_int_correct[i,:]) * dx[i]/1e+12
            E_geopotential[i,:] = (geopotential_flux_int[i,:] - geopotential_flux_int_correct[i,:]) * dx[i]/1e+12
            E_kinetic[i,:] = (kinetic_flux_int[i,:] - kinetic_flux_int_correct[i,:]) * dx[i]/1e+12

            E_internal_correct[i,:] = (internal_flux_int_correct[i,:]) * dx[i]/1e+12
            E_latent_correct[i,:] = (latent_flux_int_correct[i,:]) * dx[i]/1e+12
            E_geopotential_correct[i,:] = (geopotential_flux_int_correct[i,:]) * dx[i]/1e+12
            E_kinetic_correct[i,:] = (kinetic_flux_int_correct[i,:]) * dx[i]/1e+12

        E_total = E_internal + E_latent + E_geopotential + E_kinetic
        E_total_correct = E_internal_correct + E_latent_correct + E_geopotential_correct + E_kinetic_correct

        logging.info("Weight the energy transport by the width of grid box and apply correction!")
        logging.info("Computation of meridional energy transport on model level is finished!")

        return E_total, E_internal, E_latent, E_geopotential, E_kinetic, \
               E_total_correct, E_internal_correct, E_latent_correct, E_geopotential_correct, E_kinetic_correct

class eddy:
    def __init__(self):
        """
        Decompose the eddy components of AMET. The results include mean transport,
        standing eddy transport and transient eddy transport.
        """
        print ("Start decomposing the meridional energy transport into eddy components.")


    def calc_mean(self):
        """
        Calculate the temporal and spatial mean of certain variables.
        """
        logging.info("The calculation of temporal mean is finished!")

    def calc_eddy(self):
        """
        Decompose the eddy components of AMET.
        """
        logging.info("The calculation of standing eddy is finished!")
