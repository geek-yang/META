# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Calculate Meridional Energy Transport in the Atmosphere with Reanalysis
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.01
Last Update     : 2018.08.07
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

import numpy as np
from netCDF4 import Dataset
import os
import platform
import sys
import logging
#import matplotlib
# generate images without having a window appear
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

class met:
    def __init__(self):
        """
        Quantify the meridional energy transport.
        param T:
        """
        self.T = T
        
class eddy:
    def __init__(self):
        """
        Decompose the eddy components of AMET. The results include mean transport,
        standing eddy transport and transient eddy transport.
        param T:
        """
        self.T = T