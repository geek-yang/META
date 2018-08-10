#!/usr/bin/env python

import numpy as np
import sys
sys.path.append("../")
import scipy
from netCDF4 import Dataset
#import matplotlib
#import matplotlib.pyplot as plt
import os
import platform
#import statistics
import meta.ORAS4

print (platform.architecture())
print (os.path)

if __name__=="__main__":
    # sample
    ################################   Input zone  ######################################
    # specify data path
    #datapath_ORAS4 = '/mnt/Associate/DataBase/ORAS/ORAS4/Monthly/model'
    datapath_ORAS4 = '/project/Reanalysis/ORAS4/Monthly/Model'
    #output_path = '/home/yang/workbench/test'
    output_path = '/project/Reanalysis/ORAS4/Monthly/Model/output/meta_out/met'
    #mask_path = '/home/yang/workbench/Core_Database_AMET_OMET_reanalysis/ORAS4/mesh_mask.nc'
    mask_path = '/project/Reanalysis/ORAS4/Monthly/Model/mesh_mask.nc'
    #mask_subbasin_path = '/home/yang/workbench/Core_Database_AMET_OMET_reanalysis/land_sea_mask/mask_subbasin_ORCA_MOM.nc'
    mask_subbasin_path = '/project/Reanalysis/mask_subbasin_ORCA_MOM.nc'
    #####################################################################################
    print ('*********************** call functions *************************')
    instance = meta.ORAS4.oras(datapath_ORAS4,output_path)
    instance.omet(1958,2017,mask_path, mask_subbasin_path)