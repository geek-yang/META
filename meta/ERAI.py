# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from ERA-Interim
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.03
Last Update     : 2019.06.11
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of ECMWF. It provides an
                  entrance for the following computation includes the mass budget
                  correction, quantification of meridional energy transport, decomposition
                  of eddies.

                  ERA-Interim is a state-of-the-art atmosphere reanalysis product produced
                  by ECMWF. It spans from 1979 to 2017. Natively it is generated on a hybrid
                  sigma grid with a horizontal resolution of 0.75 x 0.75 deg and 60 vertical
                  levels.

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
                  It is recommended to combine two surface variables in a single file (e.g.
                  model_daily_075_{}_{}_z_lnsp.nc includes surface pressure and surface geopotential)
                  and combine four 3D variables in a single file (e.g. model_daily_075_{}_{}_T_q_u_v.nc
                  includes air temperature, specific humidity, zonal and meridional wind), as it can save
                  downloading time by reducing the times of requests. It can also work with files containing
                  two variables (e.g. model_daily_075_{}_{}_T_q includes temperature and specific humidity).
"""

##########################################################################
###########################   Units vacabulory   #########################
# cpT:  [J / kg K] * [K]     = [J / kg]
# Lvq:  [J / kg] * [kg / kg] = [J / kg]
# gz is [m2 / s2] = [ kg m2 / kg s2 ] = [J / kg]

# multiply by v: [J / kg] * [m / s] => [J m / kg s]
# sum over longitudes [J m / kg s] * [ m ] = [J m2 / kg s]

# integrate over pressure: dp: [Pa] = [N m-2] = [kg m2 s-2 m-2] = [kg s-2]
# [J m2 / kg s] * [Pa] = [J m2 / kg s] * [kg / s2] = [J m2 / s3]
# and factor 1/g: [J m2 / s3] * [s2 /m2] = [J / s] = [Wat]
##########################################################################

import sys
import os
import numpy as np
import logging
from netCDF4 import Dataset
import meta.massBudget
import meta.amet
#import matplotlib
import meta.saveNetCDF
# generate images without having a window appear
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

class erai:
    def __init__(self, path, out_path, package_path):
        """
        Initialize the extraction of fields from ERA-Interim.

        The data is on hybrid sigma levels. As the interpolation can introduce
        large errors to the computation of energy transport, we will follow the
        model level. The determination of reference levels is based on the estimation of
        pressure on each level with standard surface pressure 1013.25 hPa
        (see ERA-Interim archive by ECMWF.). For the actual calculation, the varying
        surface pressure should be taken into account.

        param path: the root path of the input fields
        param out_path: the location of output files
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        """
        self.path = path
        self.out_path = out_path
        # 0.75 deg per grid box latitudinally
        self.lat_unit = 240
        self.package_path = package_path

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

    def massCorrect(self, year_start, year_end, example, fields=1, method='SH'):
        """
        Mass budget correction.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param example: an example input file for loading dimensions (lat, lon)
        param fields: number of fields contained in one file, two options available
        - 1 (default) two seperate files with T,q,u,v on multiple sigma levels and lnsp,z on surface
        - 2 three seperate files, T,q and u,v and lnsp,z
        param method: numerical methods for mass correction
        - FD Calculate divergence/inverse Laplacian/gradient through finite difference
        - SH (default) Calculate divergence/inverse Laplacian/gradient through spherical harmonics
        Caveat! In order to make use of NCL, all the input fields should have ascending lat.
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
        lat = example_key.variables['latitude'][::-1]
        lon = example_key.variables['longitude'][:]
        level = example_key.variables['level'][:]
        # create space for the output
        uc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        vc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # loop for the computation of divergent corrected winds
        if fields == 1:
            for i in year:
                for j in month:
                    logging.info("Start retrieving variables for {0}(y)-{1}(m)".format(i, j))
                    datapath_T_q_u_v = os.path.join(self.path,'era{0}'.format(i),
                                                    'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j))
                    datapath_z_lnsp = os.path.join(self.path,'era{0}'.format(i),
                                                   'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j))
                    # extract fields for the calculation of tendency terms
                    if j == 1:
                        datapath_q_last = os.path.join(self.path,'era{0}'.format(i-1),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i-1, 12))
                        datapath_q_next = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{0}'.format(i-1),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i-1, 12))
                        datapath_lnsp_next = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j+1))
                        if i == year_start:
                            datapath_q_last = datapath_T_q_u_v
                            datapath_lnsp_last = datapath_z_lnsp
                    elif j == 12:
                        datapath_q_last = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{0}'.format(i+1),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i+1, 1))
                        datapath_lnsp_last = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{0}'.format(i+1),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i+1, 1))
                        if i == year_end:
                            datapath_q_next = datapath_T_q_u_v
                            datapath_lnsp_next = datapath_z_lnsp
                    else:
                        datapath_q_last = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j+1))
                    # get all the variables for the mass budget correction
                    T_q_u_v_key = Dataset(datapath_T_q_u_v)
                    z_lnsp_key = Dataset(datapath_z_lnsp)
                    # get the variable keys for the calculation of tendency during mass budget correction
                    q_last_key = Dataset(datapath_q_last)
                    q_next_key = Dataset(datapath_q_next)
                    lnsp_last_key = Dataset(datapath_lnsp_last)
                    lnsp_next_key = Dataset(datapath_lnsp_next)
                    logging.info("Get the key of all the required variables for {0}(y)-{1}(m)".format(i, j))
                    # extract variables
                    q = T_q_u_v_key.variables['q'][:,:,::-1,:]
                    lnsp = z_lnsp_key.variables['lnsp'][:,::-1,:]
                    u = T_q_u_v_key.variables['u'][:,:,::-1,:]
                    v = T_q_u_v_key.variables['v'][:,:,::-1,:]
                    # get time dimension
                    time = T_q_u_v_key.variables['time'][:]
                    # extract variables for the calculation of tendency
                    q_last = q_last_key.variables['q'][-1,:,::-1,:]
                    q_next = q_next_key.variables['q'][0,:,::-1,:]
                    lnsp_last = lnsp_last_key.variables['lnsp'][-1,::-1,:]
                    lnsp_next = lnsp_next_key.variables['lnsp'][0,::-1,:]
                    # calculate sp
                    sp = np.exp(lnsp)
                    sp_last = np.exp(lnsp_last)
                    sp_next = np.exp(lnsp_next)
                    del lnsp, lnsp_last, lnsp_next
                    logging.info("Extract all the required variables for {0}(y)-{1}(m) successfully!".format(i, j))
                    if method == 'SH':
                        # start the mass correction
                        SinkSource = meta.massBudget.correction_SH()
                        SinkSource.massCorrect(q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                                                        len(time), len(level), len(lat), len(lon), lat, lon,
                                                        self.lat_unit, self.out_path, self.package_path)
                        del u, v, q # save memory
                        # call bash to execute ncl script via subprocess
                        uc, vc = SinkSource.massCorrect(self.out_path, self.package_path)
                    elif method == 'FD':
                        # start the mass correction
                        SinkSource = meta.massBudget.correction_FD()
                        uc, vc = SinkSource.massCorrect(q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                                                    len(time), len(level), len(lat), len(lon), lat, self.lat_unit)
                    else:
                        IOError("Please choose the methods listed in the documentation!")
                    # save the output to the data pool
                    uc_pool[i-year_start,j-1,:,:] = uc
                    vc_pool[i-year_start,j-1,:,:] = vc
            # export output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncCorrect(uc_pool, vc_pool, year, lat, lon, self.out_path)
        elif fields == 2:
            for i in year:
                for j in month:
                    logging.info("Start retrieving variables for {0}(y)-{1}(m)".format(i, j))
                    datapath_T_q = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_T_q.nc'.format(i, j))
                    datapath_u_v = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_u_v.nc'.format(i, j))
                    datapath_z_lnsp = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j))
                    # extract fields for the calculation of tendency terms
                    if j == 1:
                        datapath_q_last = os.path.join(self.path,'era{0}'.format(i-1),
                                                       'model_daily_075_{0}_{1}_T_q.nc'.format(i-1, 12))
                        datapath_q_next = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{0}'.format(i-1),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i-1, 12))
                        datapath_lnsp_next = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j+1))
                        if i == year_start:
                            datapath_q_last = datapath_T_q
                            datapath_lnsp_last = datapath_z_lnsp
                    elif j == 12:
                        datapath_q_last = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{0}'.format(i+1),
                                                       'model_daily_075_{0}_{1}_T_q.nc'.format(i+1, 1))
                        datapath_lnsp_last = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{0}'.format(i+1),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i+1, 1))
                        if i == year_end:
                            datapath_q_next = datapath_T_q
                            datapath_lnsp_next = datapath_z_lnsp
                    else:
                        datapath_q_last = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{0}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{0}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j+1))
                    # get all the variables for the mass budget correction
                    T_q_key = Dataset(datapath_T_q)
                    u_v_key = Dataset(datapath_u_v)
                    z_lnsp_key = Dataset(datapath_z_lnsp)
                    # get the variable keys for the calculation of tendency during mass budget correction
                    q_last_key = Dataset(datapath_q_last)
                    q_next_key = Dataset(datapath_q_next)
                    lnsp_last_key = Dataset(datapath_lnsp_last)
                    lnsp_next_key = Dataset(datapath_lnsp_next)
                    logging.info("Get the key of all the required variables for {0}(y)-{1}(m)".format(i, j))
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
                    logging.info("Extract all the required variables for {0}(y)-{1}(m) successfully!".format(i, j))
                    # choose the methods for mass correction
                    if method == 'SH':
                        # start the mass correction
                        SinkSource = meta.massBudget.correction_SH()
                        SinkSource.massCorrect(q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                                               len(time), len(level), len(lat), len(lon), lat, lon,
                                               self.lat_unit, self.out_path, self.package_path)
                        del u, v, q # save memory
                        # call bash to execute ncl script via subprocess
                        uc, vc = SinkSource.massCorrect(self.out_path, self.package_path)
                    elif method == 'FD':
                        # start the mass correction
                        SinkSource = meta.massBudget.correction_FD()
                        uc, vc = SinkSource.massCorrect(q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                                                        len(time), len(level), len(lat), len(lon), lat, self.lat_unit)
                    else:
                        IOError("Please choose the methods listed in the documentation!")
                    # save the output to the data pool
                    uc_pool[i-year_start,j-1,:,:] = uc
                    vc_pool[i-year_start,j-1,:,:] = vc
            # export output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncCorrect(uc_pool, vc_pool, year, lat, lon, self.out_path)
        else:
            IOError("Please follow the naming rule as described in the documentation!")

    def amet(self, year_start, year_end, path_uvc, fields=1):
        """
        Quantify Meridional Energy Transport.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param path_uvc: location of the baratropic corrected winds
        param fields: number of fields contained in one file, two options available
        - 1 (default) two seperate files with T,q,u,v on multiple sigma levels and lnsp,z on surface
        - 2 three seperate files, T,q and u,v and lnsp,z
        param example: an example input file for loading dimensions (level)

        return: arrays containing AMET and its components upto differnt pressure levels
        rtype: netCDF4
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
        uvc_key = Dataset(path_uvc)
        lat = uvc_key.variables['latitude'][:]
        lon = uvc_key.variables['longitude'][:]
        #uc = uvc_key['uc'][:]
        vc = uvc_key.variables['vc'][:]
        # calculate the reference levels based on A & B and standard surface pressure
        half_level = A + B * 101325
        level = (half_level[1:] + half_level[:-1]) / 2
        # create space for the output
        # the number at the end of each name indicates the integral
        # from surface to a certain height (hPa)
        # the results will be saved per year to save memory
        # AMET in the entire column
        E = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        cpT = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        Lvq = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        gz = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        uv2 = np.zeros((len(month),len(lat),len(lon)), dtype=float)

        E_c = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        cpT_c = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        Lvq_c = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        gz_c = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        uv2_c = np.zeros((len(month),len(lat),len(lon)), dtype=float)
        # loop for the computation of AMET
        if fields == 1:
            for i in year:
                for j in month:
                    logging.info("Start retrieving variables for {0}(y)-{1}(m)".format(i, j))
                    datapath_T_q_u_v = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j))
                    datapath_z_lnsp = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j))
                    # get all the variables for the mass budget correction
                    T_q_u_v_key = Dataset(datapath_T_q_u_v)
                    z_lnsp_key = Dataset(datapath_z_lnsp)
                    logging.info("Get the keys of all the required variables for {0}(y)-{1}(m)".format(i, j))
                    # extract variables
                    T = T_q_u_v_key.variables['t'][:,:,::-1,:]
                    q = T_q_u_v_key.variables['q'][:,:,::-1,:]
                    lnsp = z_lnsp_key.variables['lnsp'][:,::-1,:]
                    z = z_lnsp_key.variables['z'][:,::-1,:]
                    u = T_q_u_v_key.variables['u'][:,:,::-1,:]
                    v = T_q_u_v_key.variables['v'][:,:,::-1,:]
                    # get time dimension
                    time = T_q_u_v_key.variables['time'][:]
                    #level = T_q_key.variables['level'][:]
                    # calculate sp
                    sp = np.exp(lnsp)
                    # calculate geopotential
                    print ('Calculate geopotential height on each model level.')
                    z_model = self.calc_gz(T, q, sp, z, A, B, len(time),
                                           len(level), len(lat), len(lon))
                    logging.info("Extracting variables successfully!")
                    AMET = meta.amet.met()
                    E[j-1,:,:], cpT[j-1,:,:], Lvq[j-1,:,:], gz[j-1,:,:],\
                    uv2[j-1,:,:], E_c[j-1,:,:], cpT_c[j-1,:,:], Lvq_c[j-1,:,:],\
                    gz_c[j-1,:,:], uv2_c[j-1,:,:] = AMET.calc_met(T, q, sp, u, v, z_model[:,:,::-1,:],
                                                                  A, B, len(time), len(level),
                                                                  len(lat), len(lon), lat,
                                                                  self.lat_unit, vc[i-year_start,j-1,:,:])
                # save output as netCDF files
                packing = meta.saveNetCDF.savenc()
                packing.ncAMET(E, cpT, Lvq, gz, uv2,
                               E_c, cpT_c, Lvq_c, gz_c, uv2_c,
                               i, level, lat, lon, self.out_path)
        elif fields == 2:
            for i in year:
                for j in month:
                    logging.info("Start retrieving variables for {0}(y)-{1}(m)".format(i, j))
                    datapath_T_q = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_T_q.nc'.format(i, j))
                    datapath_u_v = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_u_v.nc'.format(i, j))
                    datapath_z_lnsp = os.path.join(self.path,'era{0}'.format(i),
                                                'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j))
                    # get all the variables for the mass budget correction
                    T_q_key = Dataset(datapath_T_q)
                    u_v_key = Dataset(datapath_u_v)
                    z_lnsp_key = Dataset(datapath_z_lnsp)
                    logging.info("Get the keys of all the required variables for {0}(y)-{1}(m)".format(i, j))
                    # extract variables
                    T = T_q_key.variables['t'][:]
                    q = T_q_key.variables['q'][:]
                    lnsp = z_lnsp_key.variables['lnsp'][:]
                    z = z_lnsp_key.variables['z'][:]
                    u = u_v_key.variables['u'][:]
                    v = u_v_key.variables['v'][:]
                    # get time dimension
                    time = T_q_key.variables['time'][:]
                    #level = T_q_key.variables['level'][:]
                    # calculate sp
                    sp = np.exp(lnsp)
                    # calculate geopotential
                    print ('Calculate geopotential on each model level.')
                    z_model = self.calc_gz(T, q, sp, z, A, B, len(time),
                                      len(level), len(lat), len(lon))
                    logging.info("Extracting variables successfully!")
                    AMET = meta.amet.met()
                    E[j-1,:,:], cpT[j-1,:,:], Lvq[j-1,:,:], gz[j-1,:,:],\
                    uv2[j-1,:,:], E_c[j-1,:,:], cpT_c[j-1,:,:], Lvq_c[j-1,:,:],\
                    gz_c[j-1,:,:], uv2_c[j-1,:,:] = AMET.calc_met(T, q, sp, u, v, z_model[:,:,::-1,:],
                                                                  A, B, len(time), len(level),
                                                                  len(lat), len(lon), lat,
                                                                  self.lat_unit, vc[i-year_start,j-1,:,:])
                # save output as netCDF files
                packing = meta.saveNetCDF.savenc()
                packing.ncAMET(E, cpT, Lvq, gz, uv2,
                               E_c, cpT_c, Lvq_c, gz_c, uv2_c,
                               i, level, lat, lon, self.out_path)
        else:
            IOError("Please follow the naming rule as described in the documentation!")

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
        time = example_key.variables['time'][:]
        lat = example_key.variables['latitude'][:]
        lon = example_key.variables['longitude'][:]
        level = example_key.variables['level'][:]

    def calc_gz(self, T, q, sp, z, A, B, t, h, y, x):
        """
        Calculate geopotential on sigma levels.
        The method is given in ECMWF IFS 9220, from section 2.20 - 2.23.

        param T: Absolute Temperature  [K]
        param q: Specific Humidity     [kg/kg]
        param sp: Surface Pressure     [Pa]
        param z: Surface Geopotential  [m2/s2]
        param A: Constant A for Defining Sigma Level [Pa]
        param B: Constant B for Defining Sigma Level
        param t: time dimension of input fields
        param h: level dimension of input fields
        param y: latitude dimension of input fields
        param x: longitude dimension of input fields

        return: An array of geopotential.
        rtype: numpy array
        """
        logging.info("Start the computation of geopotential on model level.")
        # call the function to generate contants
        constant = self.setConstants()
        # define the half level pressure matrix
        p_half_plus = np.zeros((t, h, y, x),dtype = float)
        p_half_minus = np.zeros((t, h, y, x),dtype = float)
        # calculate the index of pressure levels
        index_level = np.arange(h)
        # calculate the pressure at each half level
        for i in index_level:
            p_half_plus[:,i,:,:] = A[i+1] + B[i+1] * sp
            p_half_minus[:,i,:,:] = A[i] + B[i] * sp
        # calculate full pressure level
        #level_full = (p_half_plus + p_half_minus) / 2
        # compute the moist temperature (virtual temperature)
        Tv = T * (1 + (constant['R_vap'] / constant['R_dry'] - 1) * q)
        # initialize the first half level geopotential
        gz_half = np.zeros((t, y, x),dtype =float)
        # initialize the full level geopotential
        gz = np.zeros((t, h, y, x),dtype = float)
        # Calculate the geopotential at each level
        # The integral should be taken from surface level to the TOA
        for i in index_level:
            # reverse the index
            i_inverse = h - 1 - i
            # the ln(p_plus/p_minus) is calculated, alpha is defined
            # an exception lies in the TOA
            # see equation 2.23 in ECMWF IFS 9220
            if i_inverse == 0:
                ln_p = np.log(p_half_plus[:,i_inverse,:,:]/10)
                alpha = np.log(2)
            else:
                ln_p = np.log(p_half_plus[:,i_inverse,:,:]/p_half_minus[:,i_inverse,:,:])
                delta_p = p_half_plus[:,i_inverse,:,:] - p_half_minus[:,i_inverse,:,:]
                alpha = 1 - p_half_minus[:,i_inverse,:,:] / delta_p * ln_p
            # calculate the geopotential of the full level (exclude surface geopotential)
            # see equation 2.22 in ECMWF IFS 9220
            gz_full = gz_half + alpha * constant['R_dry'] * Tv[:,i_inverse,:,:]
            # add surface geopotential to the full level
            # see equation 2.21 in ECMWF IFS 9220
            gz[:,i_inverse,:,:] = z + gz_full
            # renew the half level geopotential for next loop step (from p_half_minus level to p_half_plus level)
            # see equation 2.20 in ECMWF IFS 9220
            gz_half = gz_half + ln_p * constant['R_dry'] * Tv[:,i_inverse,:,:]
        logging.info("Computation of geopotential on model level is finished!")

        return gz
