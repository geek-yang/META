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

                  MERRA2 is a state-of-the-art atmosphere reanalysis product produced
                  by NASA. It spans from 1980 to 2017. Natively it is generated on a hybrid
                  sigma grid with a horizontal resolution of 0.75 x 0.75 deg and 72 vertical
                  levels.

                  The processing unit is monthly data, for the sake of memory saving.
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:
                  /JRA55
                      /merra1979
                          /MERRA2_200.inst3_3d_asm_Nv.19790101.nc4.nc
                          /MERRA2_200.inst3_3d_asm_Nv.19790102.nc4.nc
                          ...

                  Please use the default names after downloading from NASA.
                  The files are in netCDF4 format.
"""
