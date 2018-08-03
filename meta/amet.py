# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Calculate Meridional Energy Transport in the Atmosphere with Reanalysis
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.01
Last Update     : 2018.08.01
Contributor     :
Description     : This module provides a method to perform the computation
                  of meridional energy transport in the atmosphere. Moreover,
                  the decomposition of standing eddy and transient eddy is also
                  included. It works with atmosphere reanalysis datasets. In
                  order to obtain meaningful values, the mass budget correction
                  must be carried out before hand. So, this module requires input
                  of baratropic wind correction terms from the mass budget
                  conservation.
Return Values   : netCDF files
Caveat!         :
"""