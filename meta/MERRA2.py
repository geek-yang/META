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
                      /merra1980
                          /MERRA2_200.inst3_3d_asm_Nv.19800101.nc4.nc
                          /MERRA2_200.inst3_3d_asm_Nv.19800102.nc4.nc
                          ...

                  Please use the default names after downloading from NASA.
                  The files are in netCDF4 format. By default, the latitude
                  ascends.

                  The native grid of MERRA2 is different from normal reanalyses.
                  It is not natively generated on a reduced Gaussian Grid. All the
                  variables are computed on a cubed-sphere grid using GEOS-5 model.
                  As a result, it is not suitable to bring the point data back to
                  spherical harmonics for the calculations of divergence and gradients.
                  More information about its grid set-up is given in the official
                  documentation:
                  https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
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
        A = A[::-1] * 100 # change unit to Pa
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
        # date and time arrangement
        # namelist of month and days for file manipulation
        namelist_month = ['01','02','03','04','05','06','07','08','09','10','11','12']
        namelist_day = ['01','02','03','04','05','06','07','08','09','10',
                        '11','12','13','14','15','16','17','18','19','20',
                        '21','22','23','24','25','26','27','28','29','30',
                        '31']
        # index of months
        index_days_long = np.arange(31)
        index_days_short = np.arange(30)
        index_days_Feb_short = np.arange(28)
        index_days_Feb_long = np.arange(29)
        long_month_list = np.array([1,3,5,7,8,10,12])
        leap_year_list = np.array([1976,1980,1984,1988,1992,1996,2000,2004,2008,2012,2016,2020])
        # define sigma level
        A, B = self.defineSigmaLevels()
        # use example input file to load the basic dimensions information
        example_key = Dataset(example)
        #time = example_key['time'][:]
        lat = example_key.variables['lat'][:]
        lon = example_key.variables['lon'][:]
        level = example_key.variables['lev'][:]
        # create space for the output
        uc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        vc_pool = np.zeros((len(year),len(month),len(lat),len(lon)), dtype=float)
        # loop for the computation of divergent corrected winds
        # loop for calculation
        for i in year:
            ###################################################################
            ######                   begin the month loop                ######
            ###################################################################
            for j in month:
                # determine how many days are there in a month
                if j in long_month_list:
                    days = index_days_long
                elif j == 2:
                    if i in leap_year_list:
                        days = index_days_Feb_long
                    else:
                        days = index_days_Feb_short
                else:
                    days = index_days_short
                last_day = len(days)
                # create space for mass budget correction variables
                sp_mean_pool = np.zeros((len(days),len(lat),len(lon)), dtype=float)
                moisture_flux_u_int_pool = np.zeros((len(days)*8,len(lat),len(lon)), dtype=float)
                moisture_flux_v_int_pool = np.zeros((len(days)*8,len(lat),len(lon)), dtype=float)
                mass_flux_u_int_pool = np.zeros((len(days)*8,len(lat),len(lon)), dtype=float)
                mass_flux_v_int_pool = np.zeros((len(days)*8,len(lat),len(lon)), dtype=float)
                precipitable_water_int_pool = np.zeros((len(days),len(lat),len(lon)), dtype=float)
                ###################################################################
                ######                    begin the day loop                 ######
                ###################################################################
                for k in days:
                    logging.info("Start retrieving variables T,q,u,v,sp,z for from {0} (y) - {1} (m) - {2} (d) ".format(i,namelist_month[j-1],namelist_day[k]))
                    if year < 1992:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_100.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    elif year < 2001:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_200.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    elif year < 2011:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_300.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    else:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_400.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    # get the variable keys
                    var_key = Dataset(datapath_var)
                    # The shape of each variable is (8,72,361,576)
                    logging.info("Retrieving variables for from {0} (y) - {1} (m) - {2} (d) successfully!".format(i,namelist_month[j-1],namelist_day[k]))
                    # get the variable keys
                    var_key = Dataset(datapath_var)
                    # get the variabels
                    # The shape of each variable is (8,72,361,576)
                    q = var_key.variables['QV'][:]
                    sp = var_key.variables['PS'][:] #(8,361,576)
                    u = var_key.variables['U'][:]
                    v = var_key.variables['V'][:]
                    # get the basic shape
                    tt, hh, yy, xx = q.shape
                    # for the computation of tendency terms in mass correction
                    # get the value of the first day
                    if k == 0:
                        q_start = q[0,:,:,:]
                        sp_start = sp[0,:,:]
                    elif k == days[-1]:
                        q_end = q[-1,:,:,:]
                        sp_end = sp[-1,:,:]
                    # compute flux terms
                    if method == 'SH':
                        Flux = meta.massBudget.correction_SH()
                        sp_mean_pool[k,:,:], moisture_flux_u_int_pool[k*8:k*8+8,:,:],\
                        moisture_flux_v_int_pool[k*8:k*8+8,:,:], mass_flux_u_int_pool[k*8:k*8+8,:,:],\
                        mass_flux_v_int_pool[k*8:k*8+8,:,:],\
                        precipitable_water_int_pool[k,:,:] = Flux.massFlux(q, sp, u, v, A, B, tt, hh, yy, xx)
                    elif method == 'FD':
                        print("still under construction.")
                    else:
                        IOError("Please choose the methods listed in the documentation!")
                ###################################################################
                ######                     end of day loop                   ######
                ###################################################################
                # take the mean for mean terms
                sp_mean = np.mean(sp_mean_pool,0)
                precipitable_water = np.mean(precipitable_water_int_pool,0)
                # regarding the tendency term
                if j == 1:
                    # last month
                    try:
                        if year < 1992:
                            datapath_last = os.path.join(self.path, 'merra{0}'.format(i-1), 'MERRA2_100.inst3_3d_asm_Nv.{0}1231.nc4.nc'.format(i-1))
                        elif year < 2001:
                            datapath_last = os.path.join(self.path, 'merra{0}'.format(i-1), 'MERRA2_200.inst3_3d_asm_Nv.{0}1231.nc4.nc'.format(i-1))
                        elif year < 2011:
                            datapath_last = os.path.join(self.path, 'merra{0}'.format(i-1), 'MERRA2_300.inst3_3d_asm_Nv.{0}1231.nc4.nc'.format(i-1))
                        else:
                            datapath_last = os.path.join(self.path, 'merra{0}'.format(i-1), 'MERRA2_400.inst3_3d_asm_Nv.{0}1231.nc4.nc'.format(i-1))
                        # get the variable key
                        var_last = Dataset(datapath_last)
                        sp_last = var_last.variables['PS'][-1,:,:]
                        q_last = var_last.variables['QV'][-1,:,:,:]
                    except:
                        sp_last = sp_start
                        q_last = q_start
                    # next month
                    if year < 1992:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_100.inst3_3d_asm_Nv.{0}0201.nc4.nc'.format(i))
                    elif year < 2001:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_200.inst3_3d_asm_Nv.{0}0201.nc4.nc'.format(i))
                    elif year < 2011:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_300.inst3_3d_asm_Nv.{0}0201.nc4.nc'.format(i))
                    else:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_400.inst3_3d_asm_Nv.{0}0201.nc4.nc'.format(i))
                    # get the variable key
                    var_next = Dataset(datapath_next)
                    sp_next = var_next.variables['PS'][0,:,:]
                    q_next = var_next.variables['QV'][0,:,:,:]
                elif j == 12:
                    # last month
                    if year < 1992:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_100.inst3_3d_asm_Nv.{0}1130.nc4.nc'.format(i))
                    elif year < 2001:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_200.inst3_3d_asm_Nv.{0}1130.nc4.nc'.format(i))
                    elif year < 2011:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_300.inst3_3d_asm_Nv.{0}1130.nc4.nc'.format(i))
                    else:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_400.inst3_3d_asm_Nv.{0}1130.nc4.nc'.format(i))
                    # get the variable key
                    var_last = Dataset(datapath_last)
                    sp_last = var_last.variables['PS'][-1,:,:]
                    q_last = var_last.variables['QV'][-1,:,:,:]
                    # next month
                    try:
                        if year < 1992:
                            datapath_next = os.path.join(self.path, 'merra{0}'.format(i+1), 'MERRA2_100.inst3_3d_asm_Nv.{0}0131.nc4.nc'.format(i+1))
                        elif year < 2001:
                            datapath_next = os.path.join(self.path, 'merra{0}'.format(i+1), 'MERRA2_200.inst3_3d_asm_Nv.{0}0131.nc4.nc'.format(i+1))
                        elif year < 2011:
                            datapath_next = os.path.join(self.path, 'merra{0}'.format(i+1), 'MERRA2_300.inst3_3d_asm_Nv.{0}0131.nc4.nc'.format(i+1))
                        else:
                            datapath_next = os.path.join(self.path, 'merra{0}'.format(i+1), 'MERRA2_400.inst3_3d_asm_Nv.{0}0131.nc4.nc'.format(i+1))
                        # get the variable key
                        var_next = Dataset(datapath_next)
                        sp_next = var_next.variables['PS'][0,:,:]
                        q_next = var_next.variables['QV'][0,:,:,:]
                    except:
                        sp_next = sp_end
                        q_next = q_end
                else:
                    # check the number of days in last month
                    j_last = j-1
                    j_next = j+1
                    if j_last in long_month_list:
                        last_month_last_day = 31
                    elif j_last == 2:
                        if i in leap_year_list:
                            last_month_last_day = 29
                        else:
                            last_month_last_day = 28
                    else:
                        last_month_last_day = 30
                    # last month
                    if year < 1992:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_100.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j_last-1],last_month_last_day))
                    elif year < 2001:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_200.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j_last-1],last_month_last_day))
                    elif year < 2011:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_300.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j_last-1],last_month_last_day))
                    else:
                        datapath_last = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_400.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j_last-1],last_month_last_day))
                    # get the variable key
                    var_last = Dataset(datapath_last)
                    sp_last = var_last.variables['PS'][-1,:,:]
                    q_last = var_last.variables['QV'][-1,:,:,:]
                    # next month
                    if year < 1992:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_100.inst3_3d_asm_Nv.{0}{1}01.nc4.nc'.format(i,namelist_month[j_next-1]))
                    elif year < 2001:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_200.inst3_3d_asm_Nv.{0}{1}01.nc4.nc'.format(i,namelist_month[j_next-1]))
                    elif year < 2011:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_300.inst3_3d_asm_Nv.{0}{1}01.nc4.nc'.format(i,namelist_month[j_next-1]))
                    else:
                        datapath_next = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_400.inst3_3d_asm_Nv.{0}{1}01.nc4.nc'.format(i,namelist_month[j_next-1]))
                    # get the variable key
                    var_next = Dataset(datapath_next)
                    sp_next = var_next.variables['PS'][0,:,:]
                    q_next = var_next.variables['QV'][0,:,:,:]
                ################################################################
                ######     final step for the mass budget correction      ######
                ################################################################
                # take basic shape
                hh, yy, xx = q_last.shape
                if method == 'SH':
                    SinkSource = meta.massBudget.correction_SH()
                    moisture_tendency, sp_tendency = SinkSource.massTendency(q_last, q_next, sp_last, sp_next, q_start,
                                                                             q_end, sp_start, sp_end, A, B, last_day,
                                                                             hh, yy, xx)
                    SinkSource.internc(sp_mean, moisture_flux_u_int_pool, moisture_flux_v_int_pool,
                                       mass_flux_u_int_pool, mass_flux_v_int_pool, precipitable_water,
                                       moisture_tendency, sp_tendency, lat, lon, last_day, self.lat_unit,
                                       self.out_path)
                    # call bash to execute ncl script via subprocess
                    uc, vc = SinkSource.massCorrect(self.out_path, self.package_path, method='FG')
                elif method == 'FD':
                    print("still under construction.")
                else:
                    IOError("Please choose the methods listed in the documentation!")
                # save the output to the data pool
                uc_pool[i-year_start,j-1,:,:] = uc
                vc_pool[i-year_start,j-1,:,:] = vc
            ###################################################################
            ######                  end of the month loop                ######
            ###################################################################
            # export output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncCorrect(uc_pool, vc_pool, year, lat, lon, self.out_path, name='MERRA2')

    def amet_memoryWise(self, year_start, year_end, path_uvc):
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
        # date and time arrangement
        # namelist of month and days for file manipulation
        namelist_month = ['01','02','03','04','05','06','07','08','09','10','11','12']
        namelist_day = ['01','02','03','04','05','06','07','08','09','10',
                        '11','12','13','14','15','16','17','18','19','20',
                        '21','22','23','24','25','26','27','28','29','30',
                        '31']
        # index of months
        index_days_long = np.arange(31)
        index_days_short = np.arange(30)
        index_days_Feb_short = np.arange(28)
        index_days_Feb_long = np.arange(29)
        long_month_list = np.array([1,3,5,7,8,10,12])
        leap_year_list = np.array([1976,1980,1984,1988,1992,1996,2000,2004,2008,2012,2016,2020])
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
        for i in year:
            ###################################################################
            ######                   begin the month loop                ######
            ###################################################################
            for j in month:
                # determine how many days are there in a month
                if j in long_month_list:
                    days = index_days_long
                elif j == 2:
                    if i in leap_year_list:
                        days = index_days_Feb_long
                    else:
                        days = index_days_Feb_short
                else:
                    days = index_days_short
                last_day = len(days)
                # create space for the output
                E_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                cpT_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                Lvq_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                gz_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                uv2_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)

                E_c_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                cpT_c_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                Lvq_c_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                gz_c_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                uv2_c_day = np.zeros((last_day,len(lat),len(lon)), dtype=float)
                ###################################################################
                ######                    begin the day loop                 ######
                ###################################################################
                for k in days:
                    logging.info("Start retrieving variables T,q,u,v,sp,z for from {0} (y) - {1} (m) - {2} (d) ".format(i,namelist_month[j-1],namelist_day[k]))
                    if year < 1992:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_100.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    elif year < 2001:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_200.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    elif year < 2011:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_300.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    else:
                        datapath_var = os.path.join(self.path, 'merra{0}'.format(i), 'MERRA2_400.inst3_3d_asm_Nv.{0}{1}{2}.nc4.nc'.format(i,namelist_month[j-1],namelist_day[k]))
                    # get the variable keys
                    var_key = Dataset(datapath_var)
                    # The shape of each variable is (8,72,361,576)
                    logging.info("Retrieving variables for from {0} (y) - {1} (m) - {2} (d) successfully!".format(i,namelist_month[j-1],namelist_day[k]))
                    # get the variable keys
                    var_key = Dataset(datapath_var)
                    # compute gz
                    z_model = self.calc_gz(var_key)
                    print ('Calculate geopotential height on each model level.')
                    # get the variabels
                    # The shape of each variable is (8,72,361,576)
                    T = var_key.variables['T'][:]
                    q = var_key.variables['QV'][:]
                    sp = var_key.variables['PS'][:] #(8,361,576)
                    u = var_key.variables['U'][:]
                    v = var_key.variables['V'][:]
                    logging.info("Extracting variables successfully!")
                    # get the basic shape
                    tt, hh, yy, xx = q.shape
                    AMET = meta.amet.met()
                    E_day[j-1,:,:], cpT_day[j-1,:,:], Lvq_day[j-1,:,:], gz_day[j-1,:,:],\
                    uv2_day[j-1,:,:], E_c_day[j-1,:,:], cpT_c_day[j-1,:,:], Lvq_c_day[j-1,:,:],\
                    gz_c_day[j-1,:,:], uv2_c_day[j-1,:,:] = AMET.calc_met(T, q, sp, u, v, z_model[:,:,::-1,:],
                                                                          A, B, tt, hh, len(lat), len(lon), lat,
                                                                          self.lat_unit, vc[i-year_start,j-1,:,:])
                ###################################################################
                ######                     end of day loop                   ######
                ###################################################################
                E[j-1,:,:] = np.mean(E_day,0)
                cpT[j-1,:,:] = np.mean(cpT_day,0)
                Lvq[j-1,:,:] = np.mean(Lvq_day,0)
                gz[j-1,:,:] = np.mean(gz_day,0)
                uv2[j-1,:,:] = np.mean(uv2_day,0)
                E_c[j-1,:,:] = np.mean(E_c_day,0)
                cpT_c[j-1,:,:] = np.mean(cpT_c_day,0)
                Lvq_c[j-1,:,:] = np.mean(Lvq_c_day,0)
                gz_c[j-1,:,:] = np.mean(gz_c_day,0)
                uv2_c[j-1,:,:] = np.mean(uv2_c_day,0)
            ###################################################################
            ######                  end of the month loop                ######
            ###################################################################
            # save output as netCDF files
            packing = meta.saveNetCDF.savenc()
            packing.ncAMET(E, cpT, Lvq, gz, uv2, E_c, cpT_c,
                           Lvq_c, gz_c, uv2_c, i, level,
                           lat, lon, self.out_path, name='MERRA2')
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

    def calc_gz(self, var_key):
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
        # define sigma level
        A, B = self.defineSigmaLevels()
        # extract variables
        T = var_key.variables['T'][:]
        q = var_key.variables['QV'][:]
        sp = var_key.variables['PS'][:]
        z = var_key.variables['PHIS'][:]
        # basic dimension
        t, h, y, x = T.shape
        # define the half level pressure matrix
        p_half_plus = np.zeros((t, h, y, x),dtype = float)
        p_half_minus = np.zeros((t, h, y, x),dtype = float)
        # calculate the index of pressure levels
        index_level = np.arange(h)
        # calculate the pressure at each half level
        for i in index_level:
            p_half_plus[:,i,:,:] = A[i+1] + B[i+1] * sp # A is Pa already
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
