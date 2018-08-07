# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Calculate Meridional Energy Transport in the Atmosphere with Reanalysis
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.01
Last Update     : 2018.08.07
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
                 lat, lat_unit, vc, p_200, p_500, p_850):
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
        
        returns: AMET and its components (internal, latent, geopotential, kinetic energy)
                 upto different heights.
        rtype: numpy array       
        """
        dp_level = np.zeros((t, h, y, x),dtype = float)
        for i in index_level:
            dp_level[:,i,:,:] = (A[i+1] + B[i+1] * sp) - (A[i] + B[i] * sp)
        # calculate each component of total energy and take the vertical integral 
        # include mass correction component
        # increase the dimension of baratropic corrected winds to speed up
        # the calculation
        vc_4D = np.repeat(vc[:,np.newaxis,:,:],h,1)
        # Internal Energy cpT
        internal_flux = constant['cp'] * v * T * dp_level / constant['g']
        internal_flux_0 = np.mean(np.sum(internal_flux,1),0)
        internal_flux_200 = np.mean(np.sum(internal_flux[:,p_200:,:,:],1),0)
        internal_flux_500 = np.mean(np.sum(internal_flux[:,p_500:,:,:],1),0)
        internal_flux_850 = np.mean(np.sum(internal_flux[:,p_850:,:,:],1),0)
        del internal_flux
        internal_flux_correct = vc_4D * constant['cp'] * T * dp_level / constant['g']
        del T
        # Latent heat Lvq
        
        return
        
class eddy:
    def __init__(self):
        """
        Decompose the eddy components of AMET. The results include mean transport,
        standing eddy transport and transient eddy transport.
        """
        print ("Start decomposing the meridional energy transport into eddy components.")
    
    def calc_eddy(self):
        """
        
        """
    
    def calc_(self):
        
        