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
#import meta.ORAS4
import meta.wrapup

print (platform.architecture())
print (os.path)

if __name__=="__main__":
    # sample
    ################################   Input zone  ######################################
    # specify data path
    #input_path = '/project/Reanalysis/ORAS4/Monthly/Model/output/meta_out/met'
    #input_path = '/project/Reanalysis/ORAS4/Monthly/Model/output/meta_out/ohc'
    input_path = '/project/Reanalysis/ERA_Interim/Subdaily/Model/HPC_Output/meta_out/met'
    #output_path = '/project/Reanalysis/ORAS4/Monthly/Model/output/meta_out'
    output_path = '/project/Reanalysis/ERA_Interim/Subdaily/Model/HPC_Output/meta_out'
    #mask_path = '/project/Reanalysis/ORAS4/Monthly/Model/mesh_mask.nc'
    #mask_subbasin_path = '/project/Reanalysis/mask_subbasin_ORCA_MOM.nc'
    #####################################################################################
    print ('*********************** call functions *************************')
    #oras_instance = meta.wrapup.assembly(1958, 2017, input_path, output_path)
    erai_instance = meta.wrapup.assembly(1979, 2017, input_path, output_path)
    # load subbasin mask
    #subbasin_key = Dataset(mask_subbasin_path)
    #tmaskpac = subbasin_key.variables['tmaskpac_ORCA1'][:]
    #tmaskatl = subbasin_key.variables['tmaskatl_ORCA1'][:]
    #oras_instance.ncOMET(tmaskpac, tmaskatl)
    #oras_instance.ncOHC(tmaskpac, tmaskatl)
    erai_instance.ncAMET()
    