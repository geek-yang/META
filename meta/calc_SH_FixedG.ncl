; Copyright Netherlands eScience Center
; Function        :  Computation operatior via spherical harmonics
; Author          : Yang Liu (y.liu@esciencecenter.nl)
; First Built     : 2019.06.10
; Last Update     : 2019.06.11
; Description     : This module uses NCL to calculate divergence, inverse Laplacian,
;                   and gradients on spectral domain with spherical harmonics.
;
;                   The input fields should be on the ordinary lat/lon grid.
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
filename = "mass_correct_temp_flux.nc"
; load file
f = addfile(dir + filename,"r")
; extract each variable
moisture_flux_u = f->moisture_flux_u
moisture_flux_v = f->moisture_flux_v
mass_flux_u = f->mass_flux_u
mass_flux_v = f->mass_flux_v
sp_mean = f->sp_mean
moisture_tendency = f->moisture_tendency
sp_tendency = f->sp_tendency
precipitable_water = f->precipitable_water
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
uc = grad_uc / (sp_mean - g * precipitable_water)
vc = grad_vc / (sp_mean - g * precipitable_water)
; ##################################################################
; ############            create netCDF files           ############
; ##################################################################
; save netCDF file with NCL
;https://www.ncl.ucar.edu/Applications/method_2.shtml
;https://www.ncl.ucar.edu/Applications/write_netcdf.shtml
;===================================================================
; Assume variables exist and that they have associated meta data:
; (a) coordinate variables time, lev, lat, lon and (b) attributes
;===================================================================
nlat  = dimsizes(lat)
nlon  = dimsizes(lon)
diro = dir                       ; Output directory
filo = "temp_uvc.nc"             ; Output file
system("/bin/rm -f " + diro + filo) ; remove if exists
fout  = addfile (diro + filo, "c")  ; open output file
;===================================================================
; explicitly declare file definition mode. Improve efficiency.
;===================================================================
setfileoption(fout,"DefineMode",True)
;===================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;===================================================================
dimNames = (/"lat", "lon"/)
dimSizes = (/nlat,  nlon/)
dimUnlim = (/False, False/)
filedimdef(fout,dimNames,dimSizes,dimUnlim)
;===================================================================
; predefine the the dimensionality of the variables to be written out
;===================================================================
; Here we are using NCL functions to facilitate defining
; each variable's dimension name(s) and type.
; The following could be replaced with explicit, user defined dimension
; names different from those associated with the variable in memory.
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via:
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/))
;===================================================================
filevardef(fout, "lat"  ,typeof(lat), (/"lat"/))
filevardef(fout, "lon"  ,typeof(lon), (/"lon"/))
filevardef(fout, "uc"   ,typeof(uc) , (/"lat","lon"/))
filevardef(fout, "vc"   ,typeof(vc) , (/"lat","lon"/))
;===================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;====================================================================
;lat@long_name = "Latitude"
;lon@long_name = "Longitude"
;uc@long_name = "zonal correct barotropic wind"
;vc@long_name = "meridional correct barotropic wind"
;lat@unit = "degree north"
;lon@unit = "degree east"
;@unit = "m/s"
;vc@unit = "m/s"
;filevarattdef(fout,"lat",lat)   ; copy lat attributes
;filevarattdef(fout,"lon",lon)   ; copy lon attributes
;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
setfileoption(fout,"DefineMode",False)
;===================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;====================================================================
fout->lat    = (/lat/)
fout->lon    = (/lon/)
fout->uc     = (/uc/)
fout->vc     = (/vc/)
print("Successfully create netCDF file for barotropic correction winds!")
