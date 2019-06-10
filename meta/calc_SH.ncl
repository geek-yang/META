; Copyright Netherlands eScience Center
; Function        :  Computation operatior via spherical harmonics
; Author          : Yang Liu (y.liu@esciencecenter.nl)
; First Built     : 2019.06.10
; Last Update     : 2019.06.10
; Description     : This module uses NCL to calculate divergence, inverse Laplacian,
;                   and gradients on spectral domain with spherical harmonics.
;
;                   The mass budget correction is based on the method provided
;                   by Trenberth (1991). More information is available through:
;                   http://www.cgd.ucar.edu/cas/catalog/newbudgets/index.html
; Return Values   : netCDF files
; Caveat!         : All the input fields must have asending latitude
; ##################################################################
; ############      collect fields for mass correction  ############
; ##################################################################
; get environment variable from bash
dir = getenv("path")
; get file name
filename = "mass_correct_tend_flux.nc"
; load file
f = addfile(dir + filename,"r")
; extract each variable
moisture_flux_u = f-> moisture_flux_u_int
moisture_flux_v = f-> moisture_flux_v_int
time = f->time
lat = f->latitude
lon = f->longitude
; ##################################################################
; ############           calculate mass residual        ############
; ##################################################################
div_moisture_flux = uv2dvF(moisture_flux_u,moisture_flux_v)
