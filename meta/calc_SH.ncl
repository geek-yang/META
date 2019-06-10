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
; ############               define constant            ############
; ##################################################################
; gravititional acceleration [m / s2]
g = 9.80616
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
; calculate divergence of moisture flux
div_moisture_flux = uv2dvF(moisture_flux_u,moisture_flux_v)
; calculate evaporation - precipitation
E_P = moisture_tendency + dim_avg_n(div_moisture_flux,0)
; calculate divergence of mass flux
div_mass_flux = uv2dvF(mass_flux_u, mass_flux_v)
; calculate mass residual
mass_residual = sp_tendency + g * dim_avg_n(div_mass_flux,0) - g * E_P
; ##################################################################
; ############    calculate corrected barotropic wind   ############
; ##################################################################
; inverse Laplacian
X = ilapsF(mass_residual,0) ; solving the Poisson equation (k=0)
; create arrays to hold output
grad_uc = X
grad_vc = X
; compute gradient arrays
gradsf(X, grad_uc, grad_vc)
; calculate barotropic winds
uc = grad_uc / (sp_mean - g * precipitable_water_int)
vc = grad_vc / (sp_mean - g * precipitable_water_int)
; save netCDF file
https://www.ncl.ucar.edu/Applications/method_2.shtml
https://www.ncl.ucar.edu/Applications/write_netcdf.shtml
