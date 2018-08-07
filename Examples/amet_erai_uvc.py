#!/usr/bin/env python

import numpy as np
import sys
sys.path.append("../")
import scipy
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import os
import platform
#import statistics
import meta.ERAI

print (platform.architecture())
print (os.path)

if __name__=="__main__":
    # sample
    ################################   Input zone  ######################################
    # specify data path
    datapath_ERAI = '/home/yang/workbench/test'
    output_path = '/home/yang/workbench/test'
    example = '/home/yang/workbench/test/era1991/model_daily_075_1991_4_T_q.nc'
    #####################################################################################
    print ('*********************** call functions *************************')
    instance = ERAI.erai(datapath_ERAI,output_path)
    instance.massCorrect(1991,1992,example)
