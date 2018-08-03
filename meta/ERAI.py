# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from ERA-Interim
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.03
Last Update     : 2018.08.03
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
import numpy as np

class erai:
    def __init__(self, path):
        """
        Initialize the extraction of fields from ERA-Interim.
        param path: the root path of the 
        """
        self.path = path
    
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
    
