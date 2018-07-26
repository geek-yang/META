# META
<b>Meridional Energy Transport Analyzer<b>, in short as <i> META <i>, is a python library for the calculation of meridional energy transport and further analysis based on climate data. It is able to deal with the calculations of MET in both the atmosphere (AMET) and ocean (OMET). In addition, it provides diagnostic modules to perform statistical operation on MET.<br/>

Currently, MET only works with several reanalysis datasets shown below for the quantification of MET. It will embrace more datasets soon in the next update, especially with regards to the outputs from numerical models (e.g. EC-earth). <br />

## Function
META can serve as a calculator for the following tasks:
* Quantification of MET in the atmosphere
** calculate AMET
** 

## Reanalysis
This library is designed to work with the atmospheric and oceanic reanalysis products in the following list: <br>
* ERA-Interim          .[ECMWF]. <br>
* MERRA2               .[NASA].  <br>
* JRA55                .[JMA].   <br>
* ORAS4                .[ECMWF]. <br>

## Dependency
META is tested on python 3 and has the following dependencies:
* numpy
* matplotlib
* netcdf
* scipy
* cartopy
* iris


