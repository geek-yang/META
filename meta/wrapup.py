# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Save Output Files into NetCDF files
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.09
Last Update     : 2018.08.09
Contributor     :
Description     : This module aims to integrate all the single output fields into 
                  a single netCDF file. Simple postprocessing, such as zonal integral,
                  ocean basin division, will be carried out. It outputs netCDF 4 files.
                  The data has been compressed.
                  
Return Values   : netCDF files
"""

import sys
import os
import numpy as np
import logging
from netCDF4 import Dataset

class assembly:
    def __init__(self):
        """
        Assemble all the output netCDF files into a single netCDF file.
        """
        print ("Assemble all the output netCDF files into a single netCDF file.")
        
    def 