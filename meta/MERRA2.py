# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from MERRA2
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.03
Last Update     : 2018.08.03
Contributor     :
Description     : This module aims to load fields from the standard netCDF files
                  downloaded directly from online data system of NASA. It provides an
                  entrance for the following computation includes the mass budget
                  correction, quantification of meridional energy transport, decomposition
                  of eddies.
                  
                  The processing unit is monthly data, for the sake of memory saving.
Return Values   : netCDF files
Caveat!         :
"""
