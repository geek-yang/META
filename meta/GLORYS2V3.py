# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from GLORYS2V3
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.09
Last Update     : 2018.08.09
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of ECMWF. It provides an
                  entrance for the following computation includes the mass budget
                  correction, quantification of meridional energy transport, decomposition
                  of eddies.
                  
                  GLORYS2V3 is a state-of-the-art ocean reanalysis product produced by Mercator
                  Ocean in France. It spans from 1993 to 2014. Natively it is generated on 
                  ORCA025 grid and 75 vertical levels.
                  
                  The processing unit is monthly data, for the sake of memory saving.
                  
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is 
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /GLORYS2V3
                      /T
                          /GLORYS2V3_ORCA025_19930115_R20130808_gridT.nc
                          /GLORYS2V3_ORCA025_19930215_R20130808_gridT.nc
                          ...
                          /GLORYS2V3_ORCA025_20141215_R20151218_gridT.nc
                      /UV
                          /GLORYS2V3_ORCA025_19930115_R20130808_gridUV.nc
                          /GLORYS2V3_ORCA025_19930215_R20130808_gridUV.nc
                          ...
                          /GLORYS2V3_ORCA025_20141215_R20151218_gridUV.nc
                          
                  The name rule follows the default name provided by the Mercator Ocean FTP
                  when downloading the files.
                  
                  T /v fields from GLORYS2V3
                  masked array? [yes]
                  containing mask info? [yes]
                  mask from GLORYS2V3
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

class glorys:
    def __init__(self, path, out_path):
        """
        Initialize the extraction of fields from GLORYS2V3.
        
        The data is on original ORCA025 grid. As the interpolation can introduce
        large errors to the computation of energy transport, we will follow the
        model level. The determination of depth level is based on the standard
        depth info given by the bathemetry 
        
        param path: the root path of the input fields
        param out_path: the location of output files
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        param h_100: index of depth level upto 100m
        param p_300: index of depth level upto 300m
        param p_700: index of depth level upto 700m
        param p_2000: index of depth level upto 2000m
        """
        self.path = path
        self.out_path = out_path
        # number of levels for certain pressure
        self.z_100 = 29
        self.z_300 = 38
        self.z_700 = 48
        self.z_2000 = 48
    
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
    
    def omet(self, year_start, year_end, path_mask)
        """
        Quantify Meridional Energy Transport.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param path_mask: location of the land sea mask
        
        return: arrays containing OMET and its components upto differnt pressure levels
        rtype: netCDF4
        """
            theta = theta_key.variables['thetao'][:] # the unit of theta is Celsius!
            theta_nomask = theta.data
            theta_nomask[theta.mask == True] = 0
            v = v_key.variables['vo'][:]
            v_nomask = v.data
            v_nomask[v.mask == True] = 0
        
    def ohc(self, year_start, year_end, path_mask):
        
    def psi(self, year_start, year_end, path_mask):
    
    def mask(self, path_mask):