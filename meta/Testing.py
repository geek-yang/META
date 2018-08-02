#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib
import os
import statistics

if __name__=="__main__":
    # sample
    datapath_ERAI_fields = '/home/yang/workbench/Core_Database_AMET_OMET_reanalysis/ERAI/regression'
    dataset_ERAI_fields_SIC_SST_SLP = Dataset(datapath_ERAI_fields + os.sep + 'surface_ERAI_monthly_regress_1979_2016.nc')
    SST_ERAI_series = dataset_ERAI_fields_SIC_SST_SLP.variables['sst'][:]
    print (SST_ERAI_series.shape)
    instance = statistics.operator(SST_ERAI_series)
    print (type(instance))
    SST_ERAI_series_white = instance.anomaly(mode=1) # input doesn't have dimension [month]
    print (SST_ERAI_series_white.shape)
