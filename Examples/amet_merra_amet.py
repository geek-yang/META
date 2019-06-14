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
import meta.MERRA2

print (platform.architecture())
print (os.path)

if __name__=="__main__":
    # sample
    ################################   Input zone  ######################################
    # specify data path
    datapath_MERRA2 = "/projects/0/blueactn/reanalysis/MERRA2/subdaily"
    #output_path = "/home/yang/workbench/test"
    output_path = "/projects/0/blueactn/reanalysis/MERRA2/output"
    package_path = "/root/runtime/META"
    #example = '/home/yang/workbench/test/era1991/model_daily_075_1991_4_T_q.nc"
    example = "/projects/0/blueactn/reanalysis/MERRA2/subdaily/merra1980/MERRA2_100.inst3_3d_asm_Nv.19801221.nc4.nc"
    #uvc_path = "/home/yang/workbench/test/model_subdaily_075_1991_1991_uvc_point.nc'
    uvc_path = "/project/Reanalysis/ERA_Interim/Subdaily/Model/HPC_Output/meta_out/SphericalHarm/era_model_subdaily_2005_2005_uvc_point.nc"
    #####################################################################################
    print ('*********************** call functions *************************')
    instance = meta.MERRA2.merra(datapath_MERRA2, output_path, package_path)
    #instance.massCorrect(2005,2005,example)
    instance.amet(2005,2005,uvc_path)
