#!/usr/bin/env python

import numpy as np
import sys
sys.path.append("/home/yang/NLeSC/Computation_Modeling/Bjerknes/Scripts/META")
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
    #datapath_JRA55 = '/home/yang/workbench/test'
    datapath_JRA55 = "/projects/0/blueactn/reanalysis/JRA55/subdaily"
    output_path = "/home/lwc16308/reanalysis/JRA55/output"
    package_path = "/home/yang/NLeSC/Computation_Modeling/Bjerknes/Scripts/META"
    example = "/projects/0/blueactn/reanalysis/JRA55/subdaily/jra2000/anl_mdl.011_tmp.reg_tl319.2000010100_2000011018"
    uvc_path = "/project/Reanalysis/ERA_Interim/Subdaily/Model/HPC_Output/model_daily_075_1979_2016_E_point.nc"
    #####################################################################################
    print ('*********************** call functions *************************')
    instance = meta.JRA55.jra55(datapath_JRA55, output_path, package_path)
    #instance.massCorrect(1991,1991,example)
    instance.amet(1979,2016,uvc_path)
