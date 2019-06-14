# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from MERRA2
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.03
Last Update     : 2019.06.12
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of NASA. It provides an
                  entrance for the following computation includes the mass budget
                  correction, quantification of meridional energy transport, decomposition
                  of eddies.

                  MERRA2 is a state-of-the-art atmosphere reanalysis product produced
                  by NASA. It spans from 1980 to 2017. Natively it is generated on a hybrid
                  sigma grid with 72 vertical levels.

                  The processing unit is monthly data, for the sake of memory saving.
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /MERRA2
                      /merra1979
                          /MERRA2_200.inst3_3d_asm_Nv.19790101.nc4.nc
                          /MERRA2_200.inst3_3d_asm_Nv.19790102.nc4.nc
                          ...

                  Please use the default names after downloading from NASA.
                  The files are in netCDF4 format. By default, the latitude
                  ascends.
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

class merra2:
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
        self.lat_unit = 360
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
        # A and B values for the definition of sigma levelist
        # Since there are 72 model levels, there are 73 half levels, so it is for A and B values
        # the unit of A is hPa!!!!!!!!!!!!
        # from surface to TOA
        A = np.array([
            0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01,
            3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
            7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, 1.091817e+02, 1.189586e+02,
            1.286959e+02, 1.429100e+02, 1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
            2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02,
            2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
            7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01, 4.017541e+01, 3.381001e+01,
            2.836781e+01, 2.373041e+01, 1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,
            9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00, 4.076571e+00, 3.276431e+00,
            2.620211e+00, 2.084970e+00, 1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
            6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01, 2.113490e-01, 1.594950e-01,
            1.197030e-01, 8.934502e-02, 6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
            1.000000e-02,],dtype=float)
            # reverse A
        A = A[::-1]
            # the unit of B is 1!!!!!!!!!!!!
            # from surfac eto TOA
        B = np.array([
            1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
            8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
            7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01,
            5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
            2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
            6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00,],dtype=float)
            # reverse B
        B = B[::-1]

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

    def massCorrect(self, year_start, year_end, example, method='SH'):
        """
        Mass budget correction.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param example: an example input file for loading dimensions (lat, lon)
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
        lat = example_key.variables['latitude'][:]
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
                    datapath_T_q_u_v = os.path.join(self.path,'era{}'.format(i),
                                                    'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j))
                    datapath_z_lnsp = os.path.join(self.path,'era{}'.format(i),
                                                   'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j))
                    # extract fields for the calculation of tendency terms
                    if j == 1:
                        datapath_q_last = os.path.join(self.path,'era{}'.format(i-1),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i-1, 12))
                        datapath_q_next = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{0}_{1}_T_q_u_v.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{}'.format(i-1),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i-1, 12))
                        datapath_lnsp_next = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{0}_{1}_z_lnsp.nc'.format(i, j+1))
                        if i == year_start:
                            datapath_q_last = datapath_T_q_u_v
                            datapath_lnsp_last = datapath_z_lnsp
                    elif j == 12:
                        datapath_q_last = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q_u_v.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{}'.format(i+1),
                                                       'model_daily_075_{}_{}_T_q_u_v.nc'.format(i+1, 1))
                        datapath_lnsp_last = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{}'.format(i+1),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i+1, 1))
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
        else:
            IOError("Please follow the naming rule as described in the documentation!")

    def amet(self, year_start, year_end, path_uvc):
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
