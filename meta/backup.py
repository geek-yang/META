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
    def __init__(self, q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                 t, h, y, x, lat, lat_unit):
        """
        Get input fields for the mass budget correction. All the input files
        should contain the fields for the entire month.
        param T: Absolute temperature [K]
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
        self.q = q
        self.sp = sp
        self.u = u
        self.v = v
        self.q_last = q_last
        self.q_next = q_next
        self.sp_last = sp_last
        self.sp_next = sp_next
        self.A = A
        self.B = B
        self.t = t
        self.h = h
        self.y = y
        self.x = x
        self.lat = lat
        self.lat_unit = lat_unit
        # save memory
        del q, u, v
        # perform the calculation of mass budget correction 
        uc, vc = self.massCorrect()
        return uc, vc
    
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
        
    def massCorrect(self):
        """
        Perform mass budget correction. It is based on the hypothesis that the
        mass imbalance mainly comes from the baratropic winds.
        
        """
        # create constants
        constant = self.setConstants()
        # calculate the delta pressure for the tendency terms
        dp_level_start = np.zeros((self.h, self.y, self.x),dtype = float) # start of the current month
        dp_level_end = np.zeros((self.h, self.y, self.x),dtype = float) # end of the current month
        dp_level_last = np.zeros((self.h, self.y, self.x),dtype = float) # last day of the last month
        dp_level_next = np.zeros((self.h, self.y, self.x),dtype = float) # first day of the next month
        # calculate the index of pressure levels
        index_level = np.arange(self.h)
        # use matrix A and B to calculate dp based on half pressure level
        for i in index_level:
            dp_level_start[i,:,:] = (self.A[i+1] + self.B[i+1] * self.sp[0,:,:]) - \
                                    (self.A[i] + self.B[i] * self.sp[0,:,:])
            dp_level_end[i,:,:] = (self.A[i+1] + self.B[i+1] * self.sp[-1,:,:]) - \
                                  (self.A[i] + self.B[i] * self.sp[-1,:,:])
            dp_level_last[i,:,:] = (self.A[i+1] + self.B[i+1] * self.sp_last) - \
                                   (self.A[i] + self.B[i] * self.sp_last)
            dp_level_next[i,:,:] = (self.A[i+1] + self.B[i+1] * self.sp_next) - \
                                   (self.A[i] + self.B[i] * self.sp_next)
        # calculte the precipitable water tendency and take the vertical integral
        moisture_start = np.sum((self.q[0,:,:,:] * dp_level_start), 0) # start of the current month
        moisture_end = np.sum((self.q[-1,:,:,:] * dp_level_end), 0) # end of the current month
        moisture_last = np.sum((self.q_last * dp_level_last), 0) # last day of the last month
        moisture_next = np.sum((self.q_next * dp_level_next), 0) # first day of the next month        
        # compute the moisture tendency (one day has 86400s)
        moisture_tendency = ((moisture_end + moisture_next) / 2 -
                             (moisture_last + moisture_start) / 2) / (30*86400) / constant['g']
        logging.info("The calculation of precipitable water tendency is finished!")
        # calculate the delta pressure for the current month
        sp_mean = np.mean(self.sp,0)
        dp_level = np.zeros((self.t,self.h,self.y,self.x),dtype = float)
        for i in index_level:
            dp_level[:,i,:,:] = (self.A[i+1] + self.B[i+1] * self.sp) - (self.A[i] + self.B[i] * self.sp)
        # calculte the mean moisture flux for a certain month
        moisture_flux_u = self.u * self.q * dp_level / constant['g']
        moisture_flux_v = self.v * self.q * dp_level / constant['g']
        # take the vertical integral
        moisture_flux_u_int = np.sum(moisture_flux_u,1)
        moisture_flux_v_int = np.sum(moisture_flux_v,1)
        # delete intermedium variables to save memory
        del moisture_flux_u, moisture_flux_v
        # calculate zonal & meridional grid size on earth
        # the earth is taken as a perfect sphere, instead of a ellopsoid
        dx = 2 * np.pi * constant['R'] * np.cos(2 * np.pi * self.lat / 360) / self.x
        dy = np.pi * constant['R'] / self.lat_unit
        # calculate the divergence of moisture flux
        div_moisture_flux_u = np.zeros((self.t,self.y,self.x),dtype = float)
        div_moisture_flux_v = np.zeros((self.t,self.y,self.x),dtype = float)
        # Pay attention to the coordinate and sign!
        # zonal moisture flux divergence
        for i in np.arange(self.y):
            for j in np.arange(self.x):
                # the longitude could be from 0 to 360 or -180 to 180, but the index remains the same
                if j == 0:
                    div_moisture_flux_u[:,i,j] = (moisture_flux_u_int[:,i,j+1] - moisture_flux_u_int[:,i,-1]) / (2 * dx[i])
                elif j == (self.x-1) :
                    div_moisture_flux_u[:,i,j] = (moisture_flux_u_int[:,i,0] - moisture_flux_u_int[:,i,j-1]) / (2 * dx[i])
                else:
                    div_moisture_flux_u[:,i,j] = (moisture_flux_u_int[:,i,j+1] - moisture_flux_u_int[:,i,j-1]) / (2 * dx[i])
        # meridional moisture flux divergence
        for i in np.arange(self.y):
            if i == 0:
                div_moisture_flux_v[:,i,:] = -(moisture_flux_v_int[:,i+1,:] - moisture_flux_v_int[:,i,:]) / (2 * dy)
            elif i == (self.x-1):
                div_moisture_flux_v[:,i,:] = -(moisture_flux_v_int[:,i,:] - moisture_flux_v_int[:,i-1,:]) / (2 * dy)
            else:
                div_moisture_flux_v[:,i,:] = -(moisture_flux_v_int[:,i+1,:] - moisture_flux_v_int[:,i-1,:]) / (2 * dy)     
        logging.info("The calculation of divergent verically integrated moisture flux is finished!")
        # calculate evaporation minus precipitation
        E_P = np.zeros((self.y, self.x),dtype = float)
        E_P = moisture_tendency + np.mean(div_moisture_flux_u,0) + np.mean(div_moisture_flux_v,0)
        logging.info("Computation of E-P on each grid point is finished!")
        sp_tendency = ((self.sp[-1,:,:] + self.sp_next) / 2 - (self.sp_last + self.sp[0,:,:]) / 2 ) / (30 * 86400)
        logging.info("The calculation of surface pressure tendency is finished!")
        # calculte the mean mass flux for a certain month
        mass_flux_u = self.u * dp_level / constant['g']
        mass_flux_v = self.v * dp_level / constant['g']
        # take the vertical integral
        mass_flux_u_int = np.sum(mass_flux_u,1)
        mass_flux_v_int = np.sum(mass_flux_v,1)
        # calculate the divergence of moisture flux
        div_mass_flux_u = np.zeros((self.t,self.y,self.x),dtype = float)
        div_mass_flux_v = np.zeros((self.t,self.y,self.x),dtype = float)
        # zonal mass flux divergence
        for i in np.arange(self.y):
            for j in np.arange(self.x):
                # the longitude could be from 0 to 360 or -180 to 180, but the index remains the same
                if j == 0:
                    div_mass_flux_u[:,i,j] = (mass_flux_u_int[:,i,j+1] - mass_flux_u_int[:,i,-1]) / (2 * dx[i])
                elif j == (len(longitude)-1) :
                    div_mass_flux_u[:,i,j] = (mass_flux_u_int[:,i,0] - mass_flux_u_int[:,i,j-1]) / (2 * dx[i])
                else:
                    div_mass_flux_u[:,i,j] = (mass_flux_u_int[:,i,j+1] - mass_flux_u_int[:,i,j-1]) / (2 * dx[i])                    
        # meridional mass flux divergence
        for i in np.arange(self.y):
            if i == 0:
                div_mass_flux_v[:,i,:] = -(mass_flux_v_int[:,i+1,:] - mass_flux_v_int[:,i,:]) / (2 * dy)
            elif i == (self.y-1):
                div_mass_flux_v[:,i,:] = -(mass_flux_v_int[:,i,:] - mass_flux_v_int[:,i-1,:])/ (2 * dy)
            else:
                div_mass_flux_v[:,i,:] = -(mass_flux_v_int[:,i+1,:] - mass_flux_v_int[:,i-1,:]) / (2 * dy)        
        mass_residual = np.zeros((self.y,self.x),dtype = float)
        mass_residual = sp_tendency + constant['g'] * (np.mean(div_mass_flux_u,0) +\
                        np.mean(div_mass_flux_v,0)) - constant['g'] * E_P
        # delete intermedium variables to save memory
        del mass_flux_u, mass_flux_v
        logging.info("Computation of mass residual on each grid point is finished!")
        # calculate precipitable water
        precipitable_water = self.q * dp_level / constant['g']
        # take the vertical integral
        precipitable_water_int = np.mean(np.sum(precipitable_water,1),0)
        # calculate barotropic correction wind
        uc = np.zeros((self.y,self.x),dtype = float)
        vc = np.zeros((self.y,self.x),dtype = float)
        vc = mass_residual * dy / (sp_mean - constant['g'] * precipitable_water_int)
        vc[0,:] = 0 # Modification at polar points
        for i in np.arange(len(latitude)):
            uc[i,:] = mass_residual[i,:] * dx[i] / (sp_mean[i,:] -
                      constant['g'] * precipitable_water_int[i,:])
        logging.info("Computation of barotropic correction wind on each grid point is finished!")
        
        return uc, vc
