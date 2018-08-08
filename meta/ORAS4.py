# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from ORAS4
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.09
Last Update     : 2018.08.09
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of ECMWF. It provides an
                  entrance for the following computation includes the mass budget
                  correction, quantification of meridional energy transport, decomposition
                  of eddies.
                  
                  ERA-Interim is a state-of-the-art ocean reanalysis product produced by ECMWF.
                  It spans from 1958 to 2016. Natively it is generated on ORCA1 grid and 40 vertical
                  levels.
                  
                  The processing unit is monthly data, for the sake of memory saving.
                  
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is 
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /ERAI
                      /oras1958
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