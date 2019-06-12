#!/usr/bin/env python

import numpy as np
import sys
sys.path.append("/root/runtime/META")
import scipy
from netCDF4 import Dataset
#import matplotlib
#import matplotlib.pyplot as plt
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
    #datapath_ERAI = '/home/yang/workbench/test'
    datapath_ERAI = "/project/Reanalysis/ERA_Interim/Subdaily/Model"
    #output_path = "/home/yang/workbench/test"
    output_path = "/project/Reanalysis/ERA_Interim/Subdaily/Model/HPC_Output/meta_out/SphericalHarm"
    package_path = "/root/runtime/META"
    #example = '/home/yang/workbench/test/era1991/model_daily_075_1991_4_T_q.nc"
    example = "/project/Reanalysis/ERA_Interim/Subdaily/Model/era1991/model_daily_075_1991_4_T_q_u_v.nc"
    #uvc_path = "/home/yang/workbench/test/model_subdaily_075_1991_1991_uvc_point.nc'
    uvc_path = "/project/Reanalysis/ERA_Interim/Subdaily/Model/HPC_Output/meta_out/SphericalHarm/model_daily_075_1979_2016_E_point.nc"
    #####################################################################################
    print ('*********************** call functions *************************')
    instance = meta.ERAI.erai(datapath_ERAI, output_path, package_path)
    #instance.massCorrect(2005,2005,example)
    instance.amet(2005,2005,uvc_path)
