#!/usr/bin/env python

import numpy as np
import sys
sys.path.append("/home/lwc16308/META")
import scipy
from netCDF4 import Dataset
#import matplotlib
#import matplotlib.pyplot as plt
import os
import platform
#import statistics
import meta.JRA55 # for 1979-2013
import meta.JRA55_newLayout # for 2014 onward

print (platform.architecture())
print (os.path)

if __name__=="__main__":
    # sample
    ################################   Input zone  ######################################
    # specify data path
    datapath_JRA55 = "/projects/0/blueactn/reanalysis/JRA55/subdaily"
    output_path = "/home/lwc16308/reanalysis/JRA55/output/SH"
    package_path = "/home/lwc16308/META"
    example = "/projects/0/blueactn/reanalysis/JRA55/subdaily/jra2000/anl_mdl.011_tmp.reg_tl319.2000010100_2000011018"
    uvc_path = "/home/lwc16308/reanalysis/JRA55/output/SH/jra_model_subdaily_2005_2005_uvc_point.nc"
    #####################################################################################
    print ('*********************** call functions *************************')
    #instance = meta.JRA55.jra55(datapath_JRA55, output_path, package_path)
    instance = meta.JRA55_newLayout.jra55(datapath_JRA55, output_path, package_path)
    #instance.amet(1979,2016,uvc_path)
    instance.amet_memoryWise(2005,2005,uvc_path)
