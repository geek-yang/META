# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from ERA-Interim
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.03
Last Update     : 2018.08.07
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of ECMWF. It provides an
                  entrance for the following computation includes the mass budget
                  correction, quantification of meridional energy transport, decomposition
                  of eddies.
                  
                  The processing unit is monthly data, for the sake of memory saving.
                  
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is 
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /ERAI
                      /era1979
                          /model_daily_075_1979_1_T_q.nc
                          /model_daily_075_1979_1_u_v.nc
                          /model_daily_075_1979_1
                          ...
                          /model_daily_075_1979_12_T_q.nc
                          ...
                          
                  Please name the files as shown in the folder tree after downloading from MARS.
                  It is recommended to combine two variables in a single file (e.g. T_q includse 
                  temperature and specific humidity), as it can save downloading time by reducing
                  the times of requests.
                  In addition, for general cases, the module can work with files containing only 1
                  field, too. This function will be added soon.
"""

import sys
import os
import numpy as np
import logging
from netCDF4 import Dataset
import massBudget
import amet
import matplotlib
import saveNetCDF
# generate images without having a window appear
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

class erai:
    def __init__(self, path, out_path):
        """
        Initialize the extraction of fields from ERA-Interim.
        param path: the root path of the
        param out_path: the location of output files
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        """
        self.path = path
        self.out_path = out_path
        # 0.75 deg per grid box latitudinally
        self.lat_unit = 240
    
    @staticmethod
    def defineSigmaLevels():
        """
        Definine sigma levels. For more information, please visit the website of ECMWF.
        Since there are 60 model levels, there are 61 half levels, so it is for A and B values.
        
        returns: tuple containing arrays with A and B values for the definition of
                 sigma levellist
        rtype: tuple
        """
        # The unit of A is Pa!
        # It is from TOA to the surface
        A = np.array([
            0.0000000000e+000, 2.0000000000e+001, 3.8425338745e+001, 6.3647796631e+001, 9.5636962891e+001,
            1.3448330688e+002, 1.8058435059e+002, 2.3477905273e+002, 2.9849584961e+002, 3.7397192383e+002,
            4.6461816406e+002, 5.7565112305e+002, 7.1321801758e+002, 8.8366040039e+002, 1.0948347168e+003,
            1.3564746094e+003, 1.6806403809e+003, 2.0822739258e+003, 2.5798886719e+003, 3.1964216309e+003,
            3.9602915039e+003, 4.9067070313e+003, 6.0180195313e+003, 7.3066328125e+003, 8.7650546875e+003,
            1.0376125000e+004, 1.2077445313e+004, 1.3775324219e+004, 1.5379804688e+004, 1.6819472656e+004,
            1.8045183594e+004, 1.9027695313e+004, 1.9755109375e+004, 2.0222203125e+004, 2.0429863281e+004,
            2.0384480469e+004, 2.0097402344e+004, 1.9584328125e+004, 1.8864750000e+004, 1.7961359375e+004,
            1.6899468750e+004, 1.5706449219e+004, 1.4411125000e+004, 1.3043218750e+004, 1.1632757813e+004,
            1.0209500000e+004, 8.8023554688e+003, 7.4388046875e+003, 6.1443164063e+003, 4.9417773438e+003,
            3.8509133301e+003, 2.8876965332e+003, 2.0637797852e+003, 1.3859125977e+003, 8.5536181641e+002,
            4.6733349609e+002, 2.1039389038e+002, 6.5889236450e+001, 7.3677425385e+000, 0.0000000000e+000,
            0.0000000000e+000,],dtype=float)
        B = np.array([
            0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000,
            0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000,
            0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000,
            0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000,
            0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 0.0000000000e+000, 7.5823496445e-005,
            4.6139489859e-004, 1.8151560798e-003, 5.0811171532e-003, 1.1142909527e-002, 2.0677875727e-002,
            3.4121163189e-002, 5.1690407097e-002, 7.3533833027e-002, 9.9674701691e-002, 1.3002252579e-001,
            1.6438430548e-001, 2.0247590542e-001, 2.4393314123e-001, 2.8832298517e-001, 3.3515489101e-001,
            3.8389211893e-001, 4.3396294117e-001, 4.8477154970e-001, 5.3570991755e-001, 5.8616840839e-001,
            6.3554745913e-001, 6.8326860666e-001, 7.2878581285e-001, 7.7159661055e-001, 8.1125342846e-001,
            8.4737491608e-001, 8.7965691090e-001, 9.0788388252e-001, 9.3194031715e-001, 9.5182150602e-001,
            9.6764522791e-001, 9.7966271639e-001, 9.8827010393e-001, 9.9401944876e-001, 9.9763011932e-001,
            1.0000000000e+000,],dtype=float)
        
        return (A, B)
    
    def massCorrect(self, year_start, year_end, example, fields=2):
        """
        Mass budget correction.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param fields: number of fields contained in one file, two options available
        - 1 each file contains only 1 field
        - 2 (default) each file contains 2 fields
        param example: an example input file for loading dimensions (lat, lon)
        """
        # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.out_path,'history_massBudget.log') ,
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        year = np.arange(year_start, year_end+1, 1)
        month = np.arange(1, 13, 1)
        #month_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12' ]
        # define sigma level
        A, B = self.defineSigmaLevels()
        # use example input file to load the basic dimensions information
        example_key = Dataset(example)
        #time = example_key['time'][:]
        lat = example_key['latitude'][:]
        lon = example_key['longitude'][:]
        level = example_key['level'][:]
        # create space for the output
        uc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        vc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # loop for the computation of divergent corrected winds
        if fields == 2:
            for i in year:
                for j in month:
                    logging.info("Start retrieving variables for {}(y)-{}(m)".format(i, j))
                    datapath_T_q = os.path.join(self.path,'era{}'.format(i),
                                                'model_daily_075_{}_{}_T_q.nc'.format(i, j))
                    datapath_u_v = os.path.join(self.path,'era{}'.format(i),
                                                'model_daily_075_{}_{}_u_v.nc'.format(i, j))
                    datapath_z_lnsp = os.path.join(self.path,'era{}'.format(i),
                                                'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j))
                    # extract fields for the calculation of tendency terms
                    if j == 1:
                        datapath_q_last = os.path.join(self.path,'era{}'.format(i-1),
                                                       'model_daily_075_{}_{}_T_q.nc'.format(i-1, 12))
                        datapath_q_next = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{}'.format(i-1),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i-1, 12))
                        datapath_lnsp_next = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j+1))
                        if i == year_start:
                            datapath_q_last = datapath_T_q
                            datapath_lnsp_last = datapath_z_lnsp
                    elif j == 12:
                        datapath_q_last = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{}'.format(i+1),
                                                       'model_daily_075_{}_{}_T_q.nc'.format(i+1, 1))
                        datapath_lnsp_last = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{}'.format(i+1),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i+1, 1))
                        if i == year_end:
                            datapath_q_next = datapath_T_q
                            datapath_lnsp_next = datapath_z_lnsp
                    else:
                        datapath_q_last = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j+1))
                    # get all the variables for the mass budget correction
                    T_q_key = Dataset(datapath_T_q)
                    u_v_key = Dataset(datapath_u_v)
                    z_lnsp_key = Dataset(datapath_z_lnsp)
                    # get the variable keys for the calculation of tendency during mass budget correction
                    q_last_key = Dataset(datapath_q_last)
                    q_next_key = Dataset(datapath_q_next)
                    lnsp_last_key = Dataset(datapath_lnsp_last)
                    lnsp_next_key = Dataset(datapath_lnsp_next)
                    logging.info("Get the key of all the required variables for {}(y)-{}(m)".format(i, j))
                    # extract variables
                    q = T_q_key.variables['q'][:]
                    lnsp = z_lnsp_key.variables['lnsp'][:]
                    u = u_v_key.variables['u'][:]
                    v = u_v_key.variables['v'][:]
                    # get time dimension
                    time = T_q_key.variables['time'][:]
                    # extract variables for the calculation of tendency
                    q_last = q_last_key.variables['q'][-1,:,:,:]
                    q_next = q_next_key.variables['q'][0,:,:,:]
                    lnsp_last = lnsp_last_key.variables['lnsp'][-1,:,:]
                    lnsp_next = lnsp_next_key.variables['lnsp'][0,:,:]
                    # calculate sp
                    sp = np.exp(lnsp)
                    sp_last = np.exp(lnsp_last)
                    sp_next = np.exp(lnsp_next)
                    del lnsp, lnsp_last, lnsp_next
                    logging.info("Extract all the required variables for {}(y)-{}(m) successfully!".format(i, j))
                    # start the mass correction
                    uc, vc = massBudget.correction(q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                                                   len(time), len(level), len(lat), len(lon), lat, self.lat_unit)
                    # save the output to the data pool
                    uc_pool[i-year_start,j-1,:,:] = uc
                    vc_pool[i-year_start,j-1,:,:] = vc
            # export output as netCDF files
            saveNetCDF.savenc.ncCorrect(uc_pool, vc_pool, year, lat, lon, self.out_path)
        elif fields == 1:
            print ("This function will be added soon")
        else:
            IOError("Please follow the naming rule as described in the documentation!")
    
    def amet(self, year_start, year_end):
        """
        Quantify Meridional Energy Transport.
        
        """
         # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.out_path,'history_amet.log') ,
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        year = np.arange(year_start, year_end+1, 1)
        month = np.arange(1, 13, 1)
        #month_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12' ]
        # define sigma level
        A, B = self.defineSigmaLevels()
        # use example input file to load the basic dimensions information
        example_key = Dataset(example)
        time = example_key['time'][:]
        lat = example_key['latitude'][:]
        lon = example_key['longitude'][:]
        level = example_key['level'][:]
        # create space for the output
        E = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        cpT = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        Lvq = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        gz = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        
    def eddies(self, year_start, year_end):
        """
        Decompose eddy components for the AMET.
        """
          # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.out_path,'history_eddies.log') ,
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        year = np.arange(year_start, year_end+1, 1)
        month = np.arange(1, 13, 1)
        #month_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12' ]
        # define sigma level
        A, B = self.defineSigmaLevels()
        # use example input file to load the basic dimensions information
        example_key = Dataset(example)
        time = example_key['time'][:]
        lat = example_key['latitude'][:]
        lon = example_key['longitude'][:]
        level = example_key['level'][:]   
