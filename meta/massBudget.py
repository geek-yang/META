# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Mass Budget Correction for Atmosphere Reanalysis Products
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.01
Last Update     : 2018.08.07
Contributor     :
Description     : This module provides a method to perform mass budget
                  correction for further quantification of meridional
                  energy transport in the atmosphere. The outputs are
                  baratropic wind correction terms for the mass budget
                  conservation.
                  
                  The mass budget correction is based on the method provided
                  by Trenberth (1991). It assumes the mass imbalance is due
                  to the error in barotropic winds. The intermediate variables
                  are named following the definitions given by the paper. For
                  more details about application of this method, please visit
                  the the NCAR webpage:
                  http://www.cgd.ucar.edu/cas/catalog/newbudgets/index.html
Return Values   : netCDF files
Caveat!         : The program is designed to run on large cluster. It has a
                  memory requirement for 64GB.
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

##########################################################################
###########################   Units vacabulory   #########################

# cpT:  [J / kg K] * [K]     = [J / kg]
# Lvq:  [J / kg] * [kg / kg] = [J / kg]
# gz in [m2 / s2] = [ kg m2 / kg s2 ] = [J / kg]

# multiply by v: [J / kg] * [m / s] => [J m / kg s]
# sum over longitudes [J m / kg s] * [ m ] = [J m2 / kg s]

# integrate over pressure: dp: [Pa] = [N m-2] = [kg m2 s-2 m-2] = [kg s-2]
# [J m2 / kg s] * [Pa] = [J m2 / kg s] * [kg / s2] = [J m2 / s3]
# and factor 1/g: [J m2 / s3] * [s2 /m2] = [J / s] = [Wat]
##########################################################################  
class correction:
    def __init__(self):
        print ("Start the mass budget correction.")
    
    @staticmethod
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
        
    def massCorrect(self, q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                    t, h, y, x, lat, lat_unit):
        """
        Perform mass budget correction. It is based on the hypothesis that the
        mass imbalance mainly comes from the baratropic winds.
        All the input files should contain the fields for the entire month.
        param q: Specific Humidity [kg/kg]
        param sp: Surface Pressure [Pa]
        param u: Zonal Wind [m/s]
        param v: Meridional Wind [m/s]
        param A: Constant A for Defining Sigma Level
        param B: Constant B for Defining Sigma Level
        param t: time dimension of input fields
        param h: level dimension of input fields
        param y: latitude dimension of input fields
        param x: longitude dimension of input fields
        param lat: latitude
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        
        returns: barotropic corrected wind components in meridional (vc)
                 zonal (uc) direction
        rtype: numpy array        
        """
        # create constants
        constant = self.setConstants()
        # calculate the delta pressure for the tendency terms
        dp_level_start = np.zeros((h, y, x),dtype = float) # start of the current month
        dp_level_end = np.zeros((h, y, x),dtype = float) # end of the current month
        dp_level_last = np.zeros((h, y, x),dtype = float) # last day of the last month
        dp_level_next = np.zeros((h, y, x),dtype = float) # first day of the next month
        # calculate the index of pressure levels
        index_level = np.arange(h)
        # use matrix A and B to calculate dp based on half pressure level
        for i in index_level:
            dp_level_start[i,:,:] = (A[i+1] + B[i+1] * sp[0,:,:]) - \
                                    (A[i] + B[i] * sp[0,:,:])
            dp_level_end[i,:,:] = (A[i+1] + B[i+1] * sp[-1,:,:]) - \
                                  (A[i] + B[i] * sp[-1,:,:])
            dp_level_last[i,:,:] = (A[i+1] + B[i+1] * sp_last) - \
                                   (A[i] + B[i] * sp_last)
            dp_level_next[i,:,:] = (A[i+1] + B[i+1] * sp_next) - \
                                   (A[i] + B[i] * sp_next)
        # calculte the precipitable water tendency and take the vertical integral
        moisture_start = np.sum((q[0,:,:,:] * dp_level_start), 0) # start of the current month
        moisture_end = np.sum((q[-1,:,:,:] * dp_level_end), 0) # end of the current month
        moisture_last = np.sum((q_last * dp_level_last), 0) # last day of the last month
        moisture_next = np.sum((q_next * dp_level_next), 0) # first day of the next month        
        # compute the moisture tendency (one day has 86400s)
        moisture_tendency = ((moisture_end + moisture_next) / 2 -
                             (moisture_last + moisture_start) / 2) / (30*86400) / constant['g']
        logging.info("The calculation of precipitable water tendency is finished!")
        # calculate the delta pressure for the current month
        sp_mean = np.mean(sp,0)
        dp_level = np.zeros((t, h, y, x),dtype = float)
        for i in index_level:
            dp_level[:,i,:,:] = (A[i+1] + B[i+1] * sp) - (A[i] + B[i] * sp)
        # calculte the mean moisture flux for a certain month and take the vertical integral
        moisture_flux_u_int = np.sum((u * q * dp_level / constant['g']),1)
        moisture_flux_v_int = np.sum((v * q * dp_level / constant['g']),1)
        # calculate zonal & meridional grid size on earth
        # the earth is taken as a perfect sphere, instead of a ellopsoid
        dx = 2 * np.pi * constant['R'] * np.cos(2 * np.pi * lat / 360) / x
        dy = np.pi * constant['R'] / lat_unit
        # calculate the divergence of moisture flux
        div_moisture_flux_u = np.zeros((t, y, x),dtype = float)
        div_moisture_flux_v = np.zeros((t, y, x),dtype = float)
        # Pay attention to the coordinate and sign!
        # zonal moisture flux divergence
        for i in np.arange(y):
            for j in np.arange(x):
                # the longitude could be from 0 to 360 or -180 to 180, but the index remains the same
                if j == 0:
                    div_moisture_flux_u[:,i,j] = (moisture_flux_u_int[:,i,j+1] - moisture_flux_u_int[:,i,-1]) / (2 * dx[i])
                elif j == (x-1) :
                    div_moisture_flux_u[:,i,j] = (moisture_flux_u_int[:,i,0] - moisture_flux_u_int[:,i,j-1]) / (2 * dx[i])
                else:
                    div_moisture_flux_u[:,i,j] = (moisture_flux_u_int[:,i,j+1] - moisture_flux_u_int[:,i,j-1]) / (2 * dx[i])
        # meridional moisture flux divergence
        for i in np.arange(y):
            if i == 0:
                div_moisture_flux_v[:,i,:] = -(moisture_flux_v_int[:,i+1,:] - moisture_flux_v_int[:,i,:]) / (2 * dy)
            elif i == (y-1):
                div_moisture_flux_v[:,i,:] = -(moisture_flux_v_int[:,i,:] - moisture_flux_v_int[:,i-1,:]) / (2 * dy)
            else:
                div_moisture_flux_v[:,i,:] = -(moisture_flux_v_int[:,i+1,:] - moisture_flux_v_int[:,i-1,:]) / (2 * dy)     
        logging.info("The calculation of divergent verically integrated moisture flux is finished!")
        # calculate evaporation minus precipitation
        E_P = np.zeros((y, x),dtype = float)
        E_P = moisture_tendency + np.mean(div_moisture_flux_u,0) + np.mean(div_moisture_flux_v,0)
        logging.info("Computation of E-P on each grid point is finished!")
        sp_tendency = ((sp[-1,:,:] + sp_next) / 2 - (sp_last + sp[0,:,:]) / 2 ) / (30 * 86400)
        logging.info("The calculation of surface pressure tendency is finished!")
        # calculte the mean mass flux for a certain month and take the vertical integral
        mass_flux_u_int = np.sum((u * dp_level / constant['g']),1)
        mass_flux_v_int = np.sum((v * dp_level / constant['g']),1)
        # calculate the divergence of moisture flux
        div_mass_flux_u = np.zeros((t,y,x),dtype = float)
        div_mass_flux_v = np.zeros((t,y,x),dtype = float)
        # zonal mass flux divergence
        for i in np.arange(y):
            for j in np.arange(x):
                # the longitude could be from 0 to 360 or -180 to 180, but the index remains the same
                if j == 0:
                    div_mass_flux_u[:,i,j] = (mass_flux_u_int[:,i,j+1] - mass_flux_u_int[:,i,-1]) / (2 * dx[i])
                elif j == (x-1) :
                    div_mass_flux_u[:,i,j] = (mass_flux_u_int[:,i,0] - mass_flux_u_int[:,i,j-1]) / (2 * dx[i])
                else:
                    div_mass_flux_u[:,i,j] = (mass_flux_u_int[:,i,j+1] - mass_flux_u_int[:,i,j-1]) / (2 * dx[i])                    
        # meridional mass flux divergence
        for i in np.arange(y):
            if i == 0:
                div_mass_flux_v[:,i,:] = -(mass_flux_v_int[:,i+1,:] - mass_flux_v_int[:,i,:]) / (2 * dy)
            elif i == (y-1):
                div_mass_flux_v[:,i,:] = -(mass_flux_v_int[:,i,:] - mass_flux_v_int[:,i-1,:])/ (2 * dy)
            else:
                div_mass_flux_v[:,i,:] = -(mass_flux_v_int[:,i+1,:] - mass_flux_v_int[:,i-1,:]) / (2 * dy)        
        mass_residual = np.zeros((y,x),dtype = float)
        mass_residual = sp_tendency + constant['g'] * (np.mean(div_mass_flux_u,0) +\
                        np.mean(div_mass_flux_v,0)) - constant['g'] * E_P
        logging.info("Computation of mass residual on each grid point is finished!")
        # calculate precipitable water
        precipitable_water = q * dp_level / constant['g']
        # take the vertical integral
        precipitable_water_int = np.mean(np.sum(precipitable_water,1),0)
        # calculate barotropic correction wind
        uc = np.zeros((y,x),dtype = float)
        vc = np.zeros((y,x),dtype = float)
        vc = mass_residual * dy / (sp_mean - constant['g'] * precipitable_water_int)
        vc[0,:] = 0 # Modification at polar points
        for i in np.arange(y):
            uc[i,:] = mass_residual[i,:] * dx[i] / (sp_mean[i,:] -
                      constant['g'] * precipitable_water_int[i,:])
        logging.info("Computation of barotropic correction wind on each grid point is finished!")
        
        return uc, vc