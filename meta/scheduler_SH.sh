#!/bin/bash
#Copyright Netherlands eScience Center
#Function        : Scheduler for Mass Budget Correction Module
#Author          : Yang Liu (y.liu@esciencecenter.nl)
#First Built     : 2019.06.11
#Last Update     : 2019.06.11
#Contributor     :
#Description     : This module works together with the python meta.massBudget module.
#                  It acts as a scheduler to call NCL in meta.massBudget.correction_SH
#                  methods.

#Return Values   : netCDF files
#Caveat!         : The program is designed to run on large cluster. It has a
#                  memory requirement for 64GB.
# parse argument from python
out_path=$1
#year=$2
#month=$3
# setup the environment for ncl script
export path=${out_path}
# call ncl function
ncl calc_SH.ncl
# delete the temporary file
file_path="${out_path}mass_correct_temp_flux.nc"
rm ${file_path}
