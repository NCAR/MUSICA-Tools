;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("create_contvolc_files.ncl")
procedure create_contvolc_files(resolution,dstPath,grid_file)

;
begin
;
debug = False
;
;  Provide a file containing SE grid information here:
;   The File Must provide the lat/lon/area associated with each SE gridpoint.
;   Set the ShortName value to identify the grid in the names of created files.
;
ShortName = resolution

out_dir = dstPath 
;

; scaling of the original Andres and Kasgnoc (JGR, 1998)
; This takes into account non-SO2 emissions, in addition to
; a considered underestimate of the emissions in this paper
;
scaling = 2.6325
;
; read original file
;
in = addfile("/glade/work/tilmes/Volcanoes/Input/volc_continuous.nc","r")
volc_alt  = in->elevation
volc_emis = in->so2
volc_lat  = in->lat
volc_lon  = in->lon
delete(in)
;
; change volc_lon to be 0-360
;
volc_lon = where(volc_lon.lt.0.,volc_lon+360.,volc_lon)
;
; change units of volc_emis from Mg(SO2)/day to molec./s
;
; Avogadro number is 6.02214e23 molecules/g
; molecular weight of SO2 is 64
; 1 day = 86400. seconds
;
volc_emis = volc_emis * 1.e6 * 6.02214e23/64./86400. * scaling
;
; read lat-lon-area of CESM SE grid: 
;
in   = addfile(grid_file,"r")
se_lat  = in->lat
se_lon  = in->lon
se_area = in->area
delete(in)
ncol = dimsizes(se_lat)
;
; define Earth radius in cm for horizontal grid area
;
r_earth = 6.37122e6 * 1.e2       ; from models/csm_share/shr/shr_const_mod.F90, expressed in cm
;
; number of altitude layers (500 m each)
;
nalt = 20
z_mid = ispan(0,nalt-1,1) * 500. + 250.
z_int = ispan(0,nalt  ,1) * 500.
;
if ( debug ) then
  print(z_mid)
  print(z_int)
end if
;
; define time array (monthly time slices)
;
ntim = 9 * 12
;
date = new( ntim, integer)
time = new( ntim, double)
units = "days since 850-01-01 00:00"
;
; create global attributes for the output files    (*** SIMONE: NEED TO FIX THESE ATTRIBUTES UP ***)
;
fAtt                  = True            ; assign file attributes
fAtt@cesm_contact     = "Simone Tilmes (tilmes@ucar.edu), Jean-Francois Lamarque (lamar@ucar.edu)"
fAtt@source_file      = "http://aerocom.met.no/download/emissions/AEROCOM_HC/volc/Old/cont-expl-splitted/volc_continuous.nc"
fAtt@data_script      = "/glade/u/home/tilmes/ncl/regrid/Volcanoes/create_files.ncl"
creation_date = systemfunc("date +%Y%m%d")
fAtt@creation_date    = creation_date
;
cnt = 0
do nyr=0,8
;
  if ( nyr.eq.0 ) then
    year = 850
  end if
  if ( nyr.eq.1 ) then
    year = 1850
  end if
  if ( nyr.eq.2 ) then
    year = 1910
  end if
  if ( nyr.eq.3 ) then
    year = 1949
  end if
  if ( nyr.eq.4 ) then
    year = 1960
  end if
  if ( nyr.eq.5 ) then
    year = 2000
  end if
  if ( nyr.eq.6 ) then
    year = 2005
  end if
  if ( nyr.eq.7 ) then
    year = 2010
  end if
  if ( nyr.eq.8 ) then
    year = 5000
  end if
;
  do nmo=1,12
;
    date(cnt) = year * 10000 + nmo * 100 + 15
    time(cnt) = cd_inv_calendar(year,nmo,15,0,0,0,units,0)
    cnt = cnt + 1
;
  end do
end do    
;
if (debug ) then
  print(date)
  print(time)
end if
;
; define 4-d array
;
volc_4d = new( (/ntim,nalt,ncol/), double)
volc_4d = 0.
ps = new( (/ntim,ncol/), double)
ps = 101325.
;
; get the set of SE indicies for the volcanos
;
volc_icol = new( dimsizes(volc_alt), integer)
do n=0,dimsizes(volc_alt)-1
  volc_icol(n) = 0
  distmin = gc_latlon(volc_lat(n),volc_lon(n),se_lat(0),se_lon(0),2,4)
  do ii=0,(ncol-1)
    distval = gc_latlon(volc_lat(n),volc_lon(n),se_lat(ii),se_lon(ii),2,4)
    if(distval.lt.distmin) then
      distmin = distval
      volc_icol(n) = ii
    end if
  end do
  print(" n="+n+" "+volc_lat(n)+" "+se_lat(volc_icol(n))+" "+volc_lon(n)+" "+se_lon(volc_icol(n)))
end do
;
; loop over volcanoes
;
ind_alt_top = ind_nearest_coord(volc_alt      ,z_mid,0)
ind_alt_bot = ind_nearest_coord(volc_alt*2./3.,z_mid,0)
do n=0,dimsizes(volc_alt)-1
;
  area = r_earth * r_earth * se_area(volc_icol(n)) ; in cm^2
  vol  = area * 500. * 100.                                             ; in cm^3, vertical spacing is 500 m
;
  nlev = ind_alt_top(n) - ind_alt_bot(n) + 1
;
  if ( debug ) then
    print(n)
    print(volc_icol(n))
    print(nlev)
    print(ind_alt_top(n))
    print(ind_alt_bot(n))
  end if
;
; spread emissions evenly over altitude range
;
  emis = volc_emis(n)/nlev
  do k=0,nlev-1
    kk = ind_alt_bot(n) + k
    volc_4d(:,kk,volc_icol(n)) = volc_4d(:,kk,volc_icol(n)) + emis/vol ; molecules/cm3/s
  end do
;
end do
;
; define pressure grid
;
scale_height = 7600.
nlev = dimsizes(z_mid)
p_mid = 1013.25 * exp(-z_mid(::-1)/scale_height)
p_int = 1013.25 * exp(-z_int(::-1)/scale_height)
;
if ( debug ) then
  print(p_mid)
end if
;
p_mid@units = "hPa"
p_int@units = "hPa"
;
; define volc_4d attributes and coordinates
;
volc_4d!0 = "time"
volc_4d!1 = "lev"
volc_4d!2 = "ncol"
;
volc_4d&time = time
volc_4d&lev  = p_mid
;
volc_4d@units = "molec/cm3/s"
;
ps!0 = "time"
ps!1 = "ncol"
;
ps&time = time
;
ps@units     = "Pa"
ps@long_name = "Surface pressure"
;
; create SO2 output file
; Assumption: 2.5% directly goes to SO4
;
;
file_name = out_dir+"/emissions-cmip6_SO2_contvolcano_vertical_850-5000_"+ShortName+"_c" + creation_date + ".nc"
system("/bin/rm -f " + file_name)   ; remove any pre-existing file
ncdf = addfile(file_name,"c")  ; open output netCDF file
;
; add global attributes to the output file
;
fAtt@data_summary     = "Vertically distributed SO2 emissions from continuously outgassing volcanoes, using volc_continuous.nc from T. Diehl (JRC/Ispra)"
fAtt@molecular_weight = "64."
fileattdef( ncdf, fAtt )
;
; make time an UNLIMITED dimension
;
filedimdef(ncdf,"time",-1,True) 
;
var_out = volc_4d
var_out = 0.975 * volc_4d(:,::-1,:)
var_out@long_name = "SO2 emissions from continuously outgassing volcanoes"
ncdf->contvolcano = var_out
;
date!0 = "time"
date&time = time
date@long_name = "YYYYMMDD"
ncdf->date = date
ncdf->PS   = ps
ncdf->lat  = se_lat
ncdf->lon  = se_lon
ncdf->area = se_area
;
command = "ncks -h -A /glade/work/tilmes/Volcanoes/Output/p.nc " + file_name
system(command)
;
delete(ncdf)
;
;
print("done with SO2")
;
; create SO4 output file
; Assumption: 2.5% directly goes to SO4
;
file_name = out_dir+"/emissions-cmip6_SO4_contvolcano_vertical_850-5000_"+ShortName+"_c" + creation_date + ".nc"
system("/bin/rm -f " + file_name)   ; remove any pre-existing file
ncdf = addfile(file_name,"c")  ; open output netCDF file
;
; add global attributes to the output file
;
fAtt@data_summary     = "Vertically distributed SO4 emissions from continuously outgassing volcanoes, using volc_continuous.nc from T. Diehl (JRC/Ispra)"
fAtt@molecular_weight = "96."
fileattdef( ncdf, fAtt )
;
; make time an UNLIMITED dimension
;
filedimdef(ncdf,"time",-1,True)
;
var_out = volc_4d
var_out = 0.025 * volc_4d(:,::-1,:)
var_out@long_name = "SO4 emissions from continuously outgassing volcanoes"
ncdf->contvolcano = var_out
ncdf->date        = date
ncdf->PS          = ps
ncdf->lat         = se_lat
ncdf->lon         = se_lon
ncdf->area        = se_area
;
delete(ncdf)
;
command = "ncks -h -A /glade/work/tilmes/Volcanoes/Output/p.nc " + file_name
system(command)
;
print("done with SO4")
;
; create so4_a1 and so4_a2 output files
; Assumption: 2.5% directly goes to SO4, split evenly between modes 1 and 2
;
do n=0,1
;
  if (n.eq.0) then
    name = "a1"
  end if
  if (n.eq.1) then
    name = "a2"
  end if
;
  file_name = out_dir+"emissions-cmip6_so4_" + name + "_contvolcano_vertical_850-5000_"+ShortName+"_c" + creation_date + ".nc"
;
  system("/bin/rm -f " + file_name)   ; remove any pre-existing file
  ncdf = addfile(file_name,"c")  ; open output netCDF file
;
; add global attributes to the output file
;
  fAtt@data_summary     = "Vertically distributed NH4SO4 emissions in MAM mode " + name + " from continuously outgassing volcanoes, using volc_continuous.nc from T. Diehl (JRC/Ispra)"
  fAtt@molecular_weight = "114."
  fileattdef( ncdf, fAtt )
;
; make time an UNLIMITED dimension
;
  filedimdef(ncdf,"time",-1,True)
;
  var_out = volc_4d
  var_out = 0.025 * 0.5 * volc_4d(:,::-1,:)
  var_out@long_name = "so4" + name + " emissions from continuously outgassing volcanoes"
  ncdf->contvolcano = var_out
  ncdf->date        = date
  ncdf->PS          = ps
  ncdf->lat         = se_lat
  ncdf->lon         = se_lon
  ncdf->area        = se_area
;
  delete(ncdf)
;
 command = "ncks -h -A /glade/work/tilmes/Volcanoes/Output/p.nc " + file_name
  system(command)
;
end do
;
print("done with so4_a*")
;
; create num_a1 and num_a2 files
;
do n=0,1
;
  if (n.eq.0) then
    name = "a1"
  end if
  if (n.eq.1) then
    name = "a2"
  end if
;
  file_name = out_dir+"emissions-cmip6_num_" + name + "_so4_contvolcano_vertical_850-5000_"+ShortName+"_c" + creation_date + ".nc"
;
  system("/bin/rm -f " + file_name)   ; remove any pre-existing file
  ncdf = addfile(file_name,"c")  ; open output netCDF file
;
; add global attributes to the output file
;
  fAtt@data_summary     = "Vertically distributed number emissions in MAM mode " + name + " from continuously outgassing volcanoes, using volc_continuous.nc from T. Diehl (JRC/Ispra)"
  fAtt@molecular_weight = "1."
  fileattdef( ncdf, fAtt )
;
; make time an UNLIMITED dimension
;
  filedimdef(ncdf,"time",-1,True)
;
; values are from the supplement of Liu at al., GMD, MAM3 paper
;
  diam = 0.134e-6 ; m
  rho  = 1770.    ; kg/m3
  mw   = 113.     ; g/mole
  mass_particle = rho *(acos(-1.)/6.) * (diam)^3  ; mass per particle (kg/particle)
;
  var_out = volc_4d
  var_out = 0.025 * 0.5 * volc_4d(:,::-1,:)  * mw/mass_particle
  var_out@long_name = "num" + name + " so4 emissions from continuously outgassing volcanoes"
  ncdf->contvolcano = var_out
  ncdf->date       = date
  ncdf->PS         = ps
  ncdf->lat        = se_lat
  ncdf->lon        = se_lon
  ncdf->area       = se_area
;
  delete(ncdf)
;
  command = "ncks -h -A /glade/work/tilmes/Volcanoes/Output/p.nc " + file_name
  system(command)
;
end do
;
print("done with num_a*")
;
end
