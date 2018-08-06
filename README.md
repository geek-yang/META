# META :globe_with_meridians:
**Meridional Energy Transport Analyzer**, in short as **META**, is a python library for the calculation of meridional energy transport and further analysis based on climate data. It is able to deal with the calculations of MET in both the atmosphere (AMET) and ocean (OMET). In addition, it provides diagnostic modules to perform statistical operation on MET.<br/>

Currently, MET only works with several reanalysis datasets shown below for the quantification of MET. It will embrace more datasets soon in the next update, especially with regards to the outputs from numerical models (e.g. EC-earth). <br />

## Function
META can serve as a calculator for the following tasks: <br>
* Quantification of MET in the atmosphere :cloud: <br>
   mass budget correction <br>
   AMET <br>
   eddy decomposition <br>
* Quantification of MET in the ocean :ocean: <br>
   OMET <br>
   eddy decomposition <br>
   OHC <br>
   MOC <br>
* Statistical operations: :computer: <br>
   detrend (polynomial fit) <br>
   linear regression (time series, time series with spatial distribution, lead / lag) <br>
   visualization <br>


## Reanalysis
This library is designed to work with the atmospheric and oceanic reanalysis products in the following list: <br>
* ERA-Interim     [ECMWF] <br>
* MERRA2          [NASA]  <br>
* JRA55           [JMA]  <br>
* ORAS4           [ECMWF] <br>
* GLORYS2V3       [Mercator Ocean] <br>
* SODA3           [Univ. Maryland] <br>

## Dependency
META is tested on python 3 and has the following dependencies:
* numpy
* matplotlib
* netcdf
* scipy
* cartopy
* iris


