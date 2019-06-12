# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from JRA55
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2019.06.11
Last Update     : 2019.06.12
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of NCAR/UCAR Research
                  Data Archive. It provides an entrance for the following computation
                  includes the mass budget correction, quantification of meridional
                  energy transport, decomposition of eddies.

                  JRA55 is a state-of-the-art atmosphere reanalysis product produced
                  by JMA (Japan). It spans from 1979 to 2015. Natively it is generated on a hybrid
                  sigma grid with a horizontal resolution of 320 (lat) x 640 (lon) and 60 vertical
                  levels.

                  The processing unit is monthly data, for the sake of memory saving.

Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /JRA55
                      /jra1979
                          /anl_mdl.007_hgt.reg_tl319.1979010100_1979011018
                          /anl_mdl.007_hgt.reg_tl319.1979011100_1979012018
                          /anl_mdl.007_hgt.reg_tl319.1979012100_1979013118
                          /anl_mdl.007_hgt.reg_tl319.1979020100_1979021018
                          ...
                          /anl_mdl.011_tmp.reg_tl319.1979010100_1979011018
                          ...
                          /anl_mdl.033_ugrd.reg_tl319.1979010100_1979011018
                          ...
                          /anl_mdl.034_vgrd.reg_tl319.1979010100_1979011018
                          ...
                          /anl_mdl.051_spfh.reg_tl319.1979010100_1979011018
                          ...

                  Please use the default names after downloading from NCAR/UCAR Research
                  Data Archive. The files are in GRIB format.
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
import pygrib

class jra55:
    def __init__(self, path, out_path):
        """
        Initialize the extraction of fields from ERA-Interim.

        The data is on hybrid sigma levels. As the interpolation can introduce
        large errors to the computation of energy transport, we will follow the
        model level. The determination of levels is based on the estimation of
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
        self.lat_unit = 319
        #self.package_path = package_path

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
        # Since there are 60 model levels, there are 61 half levels, so it is for A and B values
        # the unit of A is Pa!!!!!!!!!!!!
        # from surface to TOA
        A = np.array([0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                      0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 1.33051011E+02, 3.64904149E+02,
                      6.34602716E+02, 9.59797167E+02, 1.34768004E+03, 1.79090740E+03, 2.29484169E+03,
                      2.84748478E+03, 3.46887149E+03, 4.16295646E+03, 4.89188083E+03, 5.67182424E+03,
                      6.47671300E+03, 7.29746989E+03, 8.12215979E+03, 8.91408220E+03, 9.65618191E+03,
                      1.03294362E+04, 1.09126384E+04, 1.13696478E+04, 1.16953716E+04, 1.18612531E+04,
                      1.18554343E+04, 1.16633554E+04, 1.12854041E+04, 1.07299494E+04, 1.00146151E+04,
                      9.16724704E+03, 8.22624491E+03, 7.20156898E+03, 6.08867301E+03, 4.95000000E+03,
                      4.00000000E+03, 3.23000000E+03, 2.61000000E+03, 2.10500000E+03, 1.70000000E+03,
                      1.37000000E+03, 1.10500000E+03, 8.93000000E+02, 7.20000000E+02, 5.81000000E+02,
                      4.69000000E+02, 3.77000000E+02, 3.01000000E+02, 2.37000000E+02, 1.82000000E+02,
                      1.36000000E+02, 9.70000000E+01, 6.50000000E+01, 3.90000000E+01, 2.00000000E+01,
                      0.00000000E+00],dtype=float)
        B = np.array([1.00000000E+00, 9.97000000E-01, 9.94000000E-01, 9.89000000E-01, 9.82000000E-01,
                      9.72000000E-01, 9.60000000E-01, 9.46000000E-01, 9.26669490E-01, 9.04350959E-01,
                      8.79653973E-01, 8.51402028E-01, 8.19523200E-01, 7.85090926E-01, 7.48051583E-01,
                      7.09525152E-01, 6.68311285E-01, 6.24370435E-01, 5.80081192E-01, 5.34281758E-01,
                      4.88232870E-01, 4.42025301E-01, 3.95778402E-01, 3.50859178E-01, 3.07438181E-01,
                      2.65705638E-01, 2.25873616E-01, 1.89303522E-01, 1.55046284E-01, 1.24387469E-01,
                      9.64456568E-02, 7.23664463E-02, 5.21459594E-02, 3.57005059E-02, 2.28538495E-02,
                      1.33275296E-02, 6.73755092E-03, 2.48431020E-03, 1.13269915E-04, 0.00000000E+00,
                      0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                      0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                      0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                      0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                      0.00000000E+00],dtype=float)

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
        namelist_month = ['01','02','03','04','05','06','07','08','09','10','11','12']
        # define sigma level
        A, B = self.defineSigmaLevels()
        # use example input file to load the basic dimensions information
        example_grbs = pygrib.open(example)
        example_key = example_grbs.message(1)
        lats, lons = benchmark_key.latlons()
        #time = example_key['time'][:]
        lat = lats[:,0]
        lon = lons[0,:]
        level = np.arange(60)
        example_grbs.close()
        # create space for the output
        uc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        vc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # loop for the computation of divergent corrected winds
        for i in year:
            # set the message counter for the extraction of surface field
            counter_surface = 0
            for j in month:
                logging.info("Start retrieving variables for {0}(y)-{1}(m)".format(i, j))
                # determine how many days are there in a month
                # for the first 10 days
                key_10d_ugrd = pygrib.open(os.path.join(self.path,'jra{0}'.format(i),
                                           'anl_mdl.033_ugrd.reg_tl319.{0}{1}0100_{2}{3}1018'.format(i,namelist_month[j-1],i,namelist_month[j-1])))
                key_10d_vgrd = pygrib.open(os.path.join(self.path,'jra{0}'.format(i),
                                           'anl_mdl.033_vgrd.reg_tl319.{0}{1}0100_{2}{3}1018'.format(i,namelist_month[j-1],i,namelist_month[j-1])))
                key_10d_spfh = pygrib.open(os.path.join(self.path,'jra{0}'.format(i),
                                           'anl_mdl.033_spfh.reg_tl319.{0}{1}0100_{2}{3}1018'.format(i,namelist_month[j-1],i,namelist_month[j-1])))



    key_10d_ugrd = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.033_ugrd.reg_tl319.%d%s0100_%d%s1018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    key_10d_vgrd = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.034_vgrd.reg_tl319.%d%s0100_%d%s1018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    key_10d_spfh = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.051_spfh.reg_tl319.%d%s0100_%d%s1018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    # for the second 10 days
    key_20d_hgt = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.007_hgt.reg_tl319.%d%s1100_%d%s2018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    key_20d_tmp = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.011_tmp.reg_tl319.%d%s1100_%d%s2018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    key_20d_ugrd = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.033_ugrd.reg_tl319.%d%s1100_%d%s2018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    key_20d_vgrd = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.034_vgrd.reg_tl319.%d%s1100_%d%s2018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    key_20d_spfh = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.051_spfh.reg_tl319.%d%s1100_%d%s2018' %(year,namelist_month[month-1],year,namelist_month[month-1]))
    # for the rest of days
    if month in long_month_list:
        last_day = 31
    elif month == 2:
        if year in leap_year_list:
            last_day = 29
        else:
            last_day = 28
    else:
        last_day = 30
    # deal with the changing last day of each month
    key_30d_hgt = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.007_hgt.reg_tl319.%d%s2100_%d%s%d18' %(year,namelist_month[month-1],year,namelist_month[month-1],last_day))
    key_30d_tmp = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.011_tmp.reg_tl319.%d%s2100_%d%s%d18' %(year,namelist_month[month-1],year,namelist_month[month-1],last_day))
    key_30d_ugrd = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.033_ugrd.reg_tl319.%d%s2100_%d%s%d18' %(year,namelist_month[month-1],year,namelist_month[month-1],last_day))
    key_30d_vgrd = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.034_vgrd.reg_tl319.%d%s2100_%d%s%d18' %(year,namelist_month[month-1],year,namelist_month[month-1],last_day))
    key_30d_spfh = pygrib.open(datapath + os.sep + 'jra%d' % (year) + os.sep + 'anl_mdl.051_spfh.reg_tl319.%d%s2100_%d%s%d18' %(year,namelist_month[month-1],year,namelist_month[month-1],last_day))
    print "Retrieving datasets successfully and return the variable key!"

                    # extract fields for the calculation of tendency terms
                    if j == 1:
                        datapath_q_last = os.path.join(self.path,'era{}'.format(i-1),
                                                       'model_daily_075_{}_{}_T_q_u_v.nc'.format(i-1, 12))
                        datapath_q_next = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q_u_v.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{}'.format(i-1),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i-1, 12))
                        datapath_lnsp_next = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j+1))
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
                        datapath_q_last = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q_u_v.nc'.format(i, j-1))
                        datapath_q_next = os.path.join(self.path,'era{}'.format(i),
                                                       'model_daily_075_{}_{}_T_q_u_v.nc'.format(i, j+1))
                        datapath_lnsp_last = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j-1))
                        datapath_lnsp_next = os.path.join(self.path,'era{}'.format(i),
                                                          'model_daily_075_{}_{}_z_lnsp.nc'.format(i, j+1))
                    # get all the variables for the mass budget correction
                    T_q_u_v_key = Dataset(datapath_T_q_u_v)
                    z_lnsp_key = Dataset(datapath_z_lnsp)
                    # get the variable keys for the calculation of tendency during mass budget correction
                    q_last_key = Dataset(datapath_q_last)
                    q_next_key = Dataset(datapath_q_next)
                    lnsp_last_key = Dataset(datapath_lnsp_last)
                    lnsp_next_key = Dataset(datapath_lnsp_next)
                    logging.info("Get the key of all the required variables for {}(y)-{}(m)".format(i, j))
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
                    logging.info("Extract all the required variables for {}(y)-{}(m) successfully!".format(i, j))
                    if method == 'SH':
                        # start the mass correction
                        SinkSource = meta.massBudget.correction_SH()
                        uc, vc = SinkSource.massCorrect(q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                                                        len(time), len(level), len(lat), len(lon), lat, lon,
                                                        self.lat_unit, self.out_path)
                    elif method == 'FD':
                        # start the mass correction
                        SinkSource = meta.massBudget.correction_FD()
                        uc, vc = SinkSource.massCorrect(q, sp, u, v, q_last, q_next, sp_last, sp_next, A, B,
                                                    len(time), len(level), len(lat), len(lon), lat, self.lat_unit)
                    else:
                        IOError("Please follow the naming rule as described in the documentation!")
                    # save the output to the data pool
                    uc_pool[i-year_start,j-1,:,:] = uc
                    vc_pool[i-year_start,j-1,:,:] = vc
            # export output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncCorrect(uc_pool, vc_pool, year, lat, lon, self.out_path)


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
        lat = uvc_key['latitude'][:]
        lon = uvc_key['longitude'][:]
        #uc = uvc_key['uc'][:]
        vc = uvc_key['vc'][:]
        # calculate the levels based on A & B and standard surface pressure
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
        # loop for the computation of AMET
        for i in year:
            for j in month:
                logging.info("Start retrieving variables for {0}(y)-{1}(m)".format(i, j))
                # determine how many days are there in a month
                # for the first 10 days
                key_10d_hgt = pygrib.open(os.path.join(self.path,'jra{0}'.format(i),
                                           'anl_mdl.033_hgt.reg_tl319.{0}{1}0100_{2}{3}1018'.format(i,namelist_month[j-1],i,namelist_month[j-1])))
                key_10d_tmp = pygrib.open(os.path.join(self.path,'jra{0}'.format(i),
                                           'anl_mdl.033_tmp.reg_tl319.{0}{1}0100_{2}{3}1018'.format(i,namelist_month[j-1],i,namelist_month[j-1])))
                key_10d_vgrd = pygrib.open(os.path.join(self.path,'jra{0}'.format(i),
                                           'anl_mdl.033_vgrd.reg_tl319.{0}{1}0100_{2}{3}1018'.format(i,namelist_month[j-1],i,namelist_month[j-1])))
                key_10d_spfh = pygrib.open(os.path.join(self.path,'jra{0}'.format(i),
                                           'anl_mdl.033_spfh.reg_tl319.{0}{1}0100_{2}{3}1018'.format(i,namelist_month[j-1],i,namelist_month[j-1])))

                    # get all the variables for the mass budget correction
                    T_q_u_v_key = Dataset(datapath_T_q_u_v)
                    z_lnsp_key = Dataset(datapath_z_lnsp)
                    logging.info("Get the keys of all the required variables for {}(y)-{}(m)".format(i, j))
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
                    print ('Calculate geopotential on each model level.')
                    gz = self.calc_gz(T, q, sp, z, A, B, len(time),
                                      len(level), len(lat), len(lon))
                    logging.info("Extracting variables successfully!")
                    AMET = meta.amet.met()
                    E[j-1,:,:], cpT[j-1,:,:], Lvq[j-1,:,:], gz[j-1,:,:],\
                    uv2[j-1,:,:], E_c[j-1,:,:], cpT_c[j-1,:,:], Lvq_c[j-1,:,:],
                    gz_c[j-1,:,:], uv2_c[j-1,:,:] = AMET.calc_met(T, q, sp, u, v, gz[:,::-1,:],
                                                                  A, B, len(time), len(level),
                                                                  len(lat), len(lon), lat,
                                                                  self.lat_unit, vc[i-year_start,j-1,:,:])
                # save output as netCDF files
                packing = meta.saveNetCDF.savenc()
                packing.ncAMET(E, cpT, Lvq, gz, uv2,
                               E_c, cpT_c, Lvq_c, gz_c, uv2_c,
                               i, level, lat, lon, self.out_path)


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
