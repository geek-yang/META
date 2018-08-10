# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from SODA3
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.10
Last Update     : 2018.08.10
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of  
                  It provides anentrance for the following computation includes the mass 
                  budget correction, quantification of meridional energy transport, 
                  decomposition of eddies.
                  
                  SODA3 is a state-of-the-art ocean reanalysis product produced by
                  University of Maryland. It spans from 1980 to 2015. Natively it is 
                  generated on MOM5 grid and 50 vertical levels.
                  
                  The processing unit is monthly data, for the sake of memory saving.
                  
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is 
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /ORAS4
                      /soda1980
                          /soda3.4.1_5dy_ocean_or_1980_01_03.nc
                          /soda3.4.1_5dy_ocean_or_1980_01_08.nc
                          ...
                      /soda1981
                          /soda3.4.1_5dy_ocean_or_1981_01_02.nc
                          /soda3.4.1_5dy_ocean_or_1981_01_07.nc
                      ...
                      /soda2015

                          
                  The name rule follows the default name provided by the Univ of Maryland FTP
                  when downloading the files.
                  
                  T /v fields from SODA3
                  masked array? [yes]
                  containing mask info? [yes]
                  mask from SODA3
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

class soda:
    def __init__(self, path, out_path):
        """
        Initialize the extraction of fields from SODA3.
        
        The data is on original MOM55 grid. As the interpolation can introduce
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
        param path_mask: location of the file containing land sea mask and coordinate info
        
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
        
            theta = theta_key.variables['thetao'][:] # the unit of theta is Celsius!
            theta_nomask = theta.data
            theta_nomask[theta.mask == True] = 0
            v = v_key.variables['vo'][:]
            v_nomask = v.data
            v_nomask[v.mask == True] = 0        