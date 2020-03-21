[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3609505.svg)](https://doi.org/10.5281/zenodo.3609505)
# META :globe_with_meridians:
**Meridional Energy Transport Analyzer**, in short as **META**, is a python library designed for the calculation of meridional energy transport (MET) and extended analysis with climatological data. It is able to deal with the calculations of MET in both the atmosphere (AMET) and ocean (OMET). In addition, it provides diagnostic modules to perform statistical operations on MET.<br/>

Initially, MET was developed to work with 6 state-of-the-art reanalysis datasets. It can be employed to work different meteorological data sets, for instance, the outputs from numerical models (e.g. EC-earth). <br />

## Function :dart:
META can serve as a calculator for the following tasks: <br>
* Quantification of MET in the atmosphere :cloud: <br>
   mass budget correction <br>
   AMET quantificaton <br>
* Quantification of MET in the ocean :ocean: <br>
   OMET <br>
   eddy decomposition <br>
   ocean heat content (OHC) <br>
* Statistical operations and visualization: :computer: <br>
   detrend signals (polynomial fit) <br>
   linear regression (time series, time series with spatial distribution, lead / lag) <br>
   visualization <br>


## Reanalysis :satellite:
This library is designed to work with multiple state-of-the-art atmospheric and oceanic reanalysis products in the following list: <br>
* ERA-Interim     [ECMWF] <br>
* MERRA2          [NASA]  <br>
* JRA55           [JMA]  <br>
* ORAS4           [ECMWF] <br>
* GLORYS2V3       [Mercator Ocean] <br>
* SODA3           [Univ. Maryland] <br>

## Dependency :books:
META is tested on python 2.6, 2.7 and 3.6 and has the following dependencies:
* numpy
* matplotlib
* netCDF4
* scipy
* cartopy
* iris
* pygrib

It also requires NCL for the barotropic mass correction as the computation may take place via Spherical Harmonics.

## Modules :floppy_disk:
Here is a brief introduction of all the modules included in this package:
* amet: quantify atmospheric meridional energy transport
* omet: quantify oceanic meridional energy transport and ocean heat content
* massBudget: mass budget correction on structured grid / Gaussian grid
* statistics: diagnostics tools
* saveNetCDF: netcdf file wrapper
* visualizer: plotting tools

## Get start :clapper:
All the modules have been tested with six reanalysis data sets. Operations within given reanalysis data sets, including mass budget correction and heat transport quantification, are illustrated using Jupyter Notebook in the folder "Examples". <br> 

In order to use the existing workflow, the input files should be organzied in a certain structure with certain file names (the file names can be customized in each script named after reanalysis product, e.g. "MERRA.py"). Since different data sets have their own naming convention and saved structure, the file structure and file names are listed in the beginning of each script. Please check the code of your target reanalysis product.<br>

For more information about how to use/customize each module, please check the comments in the code at the beginning of each function. <br>

## Cite our work :heart:
If you use this package, please cite it via: <br>
```latex
@Article{esd-11-77-2020,
AUTHOR = {Liu, Y. and Attema, J. and Moat, B. and Hazeleger, W.},
TITLE = {Synthesis and evaluation of historical meridional heat transport from midlatitudes towards the Arctic},
JOURNAL = {Earth System Dynamics},
VOLUME = {11},
YEAR = {2020},
NUMBER = {1},
PAGES = {77--96},
URL = {https://www.earth-syst-dynam.net/11/77/2020/},
DOI = {10.5194/esd-11-77-2020}
}
```
