;
; read original CEDS/CMIP6 SSP emissions files
; regrid to 1 deg and convert from kg/m2/s to molecules/cm2/s
; sum all anthro sectors to one variable
;  (requires gnu compiler, not intel, on cheyenne/geyser)
;

begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; USER CHANGES START
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Location of files created with Regrid_fv2se_cams_anthro.ncl
dir_inp = "/glade/scratch/$USER/emis/cams_ne30np4/"
;dir_inp = "/glade/scratch/$USER/emis/cams_ne0np4_Asia03_ne30x8/"
;
version= "v4.2_R2.1"

; model resolution (as written in filenames)
;resolution = "ne0CONUSne30x8"
resolution = "ne30np4"
;resolution = "ne0np4_Asia03_ne30x8"
;resolution = "ne30pg3"
;
; Directory for new files
dir_out = "/glade/scratch/$USER/emis/cesm_ne30np4/"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; USER CHANGES END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; constants
;
  PI  = get_pi("float" )
  avog = 6.022e23
  earth_rad  = 6.371e8 ; in cm
  earth_area = 4. * acos(-1.) * earth_rad * earth_rad
; define creation date
;
      creation_date = systemfunc ("date +%Y%m%d")
;
; make sure output dir exists
;
  system("mkdir -p " + dir_out)
;
; loop over species
;
   do nspc=0,28
;
    if ( nspc.eq.0 ) then
      spc_name = "BC"
      spc_name1 = "bc"
      spc_name2 = "bc_a4"
      mw = 12.
      sf = 1.0
    end if
    if ( nspc.eq.1 ) then
      spc_name = "CO"
      spc_name1 = "co"
      spc_name2 = "CO"
      mw = 28.
      sf = 1.0
    end if
    if ( nspc.eq.2 ) then
      spc_name = "NH3"
      spc_name1 = "nh3"
      spc_name2 = "NH3"
      mw = 17.
      sf = 1.0
    end if
    if ( nspc.eq.3 ) then
      spc_name = "NOx"
      spc_name1 = "nox"
      spc_name2 = "NOx"
      mw = 30.
      ;sf = 1.0
      ; different factor required because of mistake in original emissions
      sf = 46./30. 
    end if
    if ( nspc.eq.4 ) then
      spc_name = "OC"
      spc_name1 = "oc"
      spc_name2 = "pom_a4"
       sf = 1.4
      mw = 12.
    end if
    if ( nspc.eq.5 ) then
      spc_name = "SO2"
      spc_name1 = "so2"
      spc_name2 = "SO2"
      mw = 64.
      sf = 1.0
    end if
    if ( nspc.eq.6 ) then
     spc_name = "C2H6"
     spc_name1 = (/"ethane"/)
     mw = 30.
     sf = 1.0
    end if
    if ( nspc.eq.7 ) then
     spc_name = "C3H8"
     spc_name1 = (/"propane"/)
     mw = 44.
     sf = 1.0
   end if
   if ( nspc.eq.8 ) then
     spc_name = "BIGALK"
      spc_name1 = (/"butanes","pentanes","hexanes","esters","ethers"/)
      mw = 72.
      mw_spec_set = (/58., 72., 86., 184., 81./) 
      sf = 1.0
   end if
    if ( nspc.eq.9 ) then
     spc_name = "C2H4"
      spc_name1= (/"ethene"/)
      mw = 28.
      sf = 1.0
   end if
    if ( nspc.eq.10 ) then
     spc_name = "C3H6"
      spc_name1 = (/"propene"/)
      mw = 42.
      sf = 1.0
   end if
    if ( nspc.eq.11 ) then
     spc_name = "C2H2"
      spc_name1 = (/"ethyne"/)
      mw = 26.
      sf = 1.0
   end if
    if ( nspc.eq.12) then
     spc_name = "BIGENE"
      ;spc_name1 = (/"olef"/)
      spc_name1 = (/"other-alkenes-and-alkynes"/)
      mw = 56.
      sf = 1.0
   end if
   if ( nspc.eq.13 ) then
     spc_name = "BENZENE"
      spc_name1= (/"benzene"/)
      mw = 78.
      sf = 1.0
   end if
    if ( nspc.eq.14 ) then
     spc_name = "TOLUENE"
      spc_name1 = (/"toluene"/)
      mw = 92.
      sf = 1.0
   end if
    if ( nspc.eq.15 ) then
     spc_name = "XYLENES"
      spc_name1 = (/"xylene","trimethylbenzene","other-aromatics"/)
      mw = 106.
      mw_spec_set = (/106., 120., 126./) 
      sf = 1.0
   end if
    if ( nspc.eq.16 ) then
     spc_name = "CH2O"
      spc_name1 = (/"formaldehyde"/)
      mw = 30.
      sf = 1.0
   end if
    if ( nspc.eq.17 ) then
     spc_name = "CH3CHO"
      spc_name1 = (/"other-aldehydes"/)
      mw = 44.
      sf = 1.
   end if
  if ( nspc.eq.18 ) then
     spc_name = "CH3OH"
      spc_name1 = (/"alcohols"/)
      sf = 0.15
      mw = 32.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.19 ) then
     spc_name = "C2H5OH"
      spc_name1 = (/"alcohols"/)
      sf = 0.85
      mw = 46.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.20 ) then
     spc_name = "CH3COCH3"
      spc_name1 = (/"ketones"/)
      sf = 0.2
      mw = 58.
   end if
    if ( nspc.eq.21 ) then
     spc_name = "MEK"
      spc_name1 = (/"ketones"/)
      sf = 0.8
      mw = 72.
   end if
    if ( nspc.eq.22 ) then
     spc_name = "HCOOH"
      spc_name1 = (/"acids"/)
      sf = 0.5
      mw = 46.
   end if
  if ( nspc.eq.23 ) then
     spc_name = "CH3COOH"
      spc_name1 = (/"acids"/)
      sf = 0.5
      mw = 60.
   end if
    if ( nspc.eq.24 ) then
     spc_name = "ISOP"
      spc_name1 = (/"isoprene"/)
      sf = 1.
      mw = 68.
   end if
    if ( nspc.eq.25 ) then
     spc_name = "MTERP"
      spc_name1  =(/"monoterpenes"/)
      sf = 1.
      mw = 136.
   end if
    if ( nspc.eq.26 ) then
     spc_name = "HCN"
      spc_name1 = "co"
      sf = 0.003
      mw = 27.
    end if
    if ( nspc.eq.27 ) then
     spc_name = "CH3CN"
     spc_name1 = "co"
      sf = 0.002
      mw = 41.
    end if
    if ( nspc.eq.28) then
     spc_name = "SVOC"
     spc_name1 = "SVOC"
     sf = 1.
     mw = 310.
    end if

    nspec1 = dimsizes(spc_name1)
    str1 = str_join(spc_name1,"+")
    if sf .lt. 1 then
       str2 = sprintf("%4.2e", sf)
    else
      str2= " "
    end if
    hist_spec = "CAMS species: "+str2+str1
    print(hist_spec)

if nspc .le. 27 then

  ced_cnt = 0
  do ced  = 0, nspec1-1
    files = systemfunc("ls " + dir_inp + resolution +"_CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_"+spc_name1(ced) + "_v4.2-R2.1.nc")
    exists = isfilepresent(files)
   if(exists) then
    ced_cnt = ced_cnt+1
    print(spc_name1(ced))
   end if
  end do
  if  ced_cnt.eq.nspec1 then
    file_available = 1
  else 
    file_available = 0
    if nspec1.gt.1 then
     delete(mw_spec_set)
    end if
  end if
  print("file_available "+file_available)

if file_available.eq.1 then
  do ced  = 0, nspec1-1
;
; get anthro files
    ;files = systemfunc("ls " + dir_inp + "CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_"+spc_name1(ced) + version+".nc")
    ;files = systemfunc("ls " + dir_inp + "CANTv3.1_Glb_0.1x0.1_anthro_"+spc_name1(ced) + version+".nc")
    files = systemfunc("ls " + dir_inp + resolution +"_CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_"+spc_name1(ced) + "_v4.2-R2.1.nc")
    nfils = dimsizes(files)

; should be just 1 file per species
    nfil=0
;
      file_inp = files(nfil)
      in = addfile(file_inp,"r")
   if ced.eq.0 then
      lat_in   = in->lat
      lon_in   = in->lon
      ncol = dimsizes(lat_in)
      area = in->area
      time      = in->time
      ntim = dimsizes(time)
      time0 = ispan(0,ntim-1,1)
      printVarSummary(time)
      time0 = (/0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366, 397, 425,   \
    456, 486, 517, 547, 578, 609, 639, 670, 700, 731, 762, 790, 821, 851,   \
    882, 912, 943, 974, 1004, 1035, 1065, 1096, 1127, 1155, 1186, 1216, 1247,  \
    1277, 1308, 1339, 1369, 1400, 1430, 1461, 1492, 1521, 1552, 1582, 1613,   \
    1643, 1674, 1705, 1735, 1766, 1796, 1827, 1858, 1886, 1917, 1947, 1978,   \
    2008, 2039, 2070, 2100, 2131, 2161, 2192, 2223, 2251, 2282, 2312, 2343,   \
    2373, 2404, 2435, 2465, 2496, 2526, 2557, 2588, 2616, 2647, 2677, 2708,   \
    2738, 2769, 2800, 2830, 2861, 2891, 2922, 2953, 2982, 3013, 3043, 3074,   \
    3104, 3135, 3166, 3196, 3227, 3257, 3288, 3319, 3347, 3378, 3408, 3439,   \
    3469, 3500, 3531, 3561, 3592, 3622, 3653, 3684, 3712, 3743, 3773, 3804,   \
    3834, 3865, 3896, 3926, 3957, 3987, 4018, 4049, 4077, 4108, 4138, 4169,   \
    4199, 4230, 4261, 4291, 4322, 4352, 4383, 4414, 4443, 4474, 4504, 4535,   \
    4565, 4596, 4627, 4657, 4688, 4718/)
      
      delete(time)
      time=time0
      time!0="time"
      time@units="days since 2008-01-15 00:00:00"
      time@calendar="gregorian"
     printVarSummary(time)
;
; define date array
;
     
      utc_date = cd_calendar(time, 0)
      date = doubletointeger(time)
      date@units = "YYYYMMDD"
      date@long_name = "date"
      date!0 = "time"
      date&time = time
      date = tointeger(utc_date(:,0)) * 10000 + tointeger(utc_date(:,1)) * 100 + tointeger(utc_date(:,2))
      print(date)
  end if
;

;
      print("Processing " + file_inp)
;
; find emissions variable name
;
      vnames = getfilevarnames(in)
      nnames = dimsizes(vnames)
print(vnames)

; get sum of all 
   var1 = in->ind
   mw_spec = stringtofloat(var1@molecular_weight)
   print(mw_spec)
   if mw_spec.ne.mw then
     print("mw_spec different from mw")
   end if
   mw_s = mw
   if nspec1 .gt. 1 then
     mw_s = mw_spec_set(ced)
   end if

    print("read ene")
      sect_energy = in->ene 
    print("read ind")
      sect_industry = in->ind 
      sect_industry2 = in->fef 
      sect_solvents= in->slv
      sect_transport = in->tro 
      sect_transport2 = in->tnr 
      sect_resident = in->res 
      sect_agriculture= in->awb
      sect_agriculture2= in->ags
      sect_agriculture3= in->agl
      sect_waste = in->swd 
    print("read shp")
      sect_shipping = in->shp 
    
    printVarSummary(sect_energy)

   if ced.eq.0 then
      sect_other1=sect_agriculture
      sect_other1=sf*(sect_agriculture+sect_agriculture2+sect_agriculture3+sect_solvents+sect_waste)
      sect_energy1=sf*sect_energy
      sect_industry1=sf*(sect_industry+sect_industry2)
      sect_transport1=sf*(sect_transport+sect_transport2)
      sect_resident1=sf*sect_resident
      sect_shipping1=sf*sect_shipping
      
      sect_other1= sect_other1* 1.e3 * avog/mw_s* 1.e-4
      sect_energy1 = sect_energy1 * 1.e3 * avog/mw_s* 1.e-4
      sect_industry1 = sect_industry1* 1.e3 * avog/mw_s* 1.e-4
      sect_transport1 = sect_transport1 * 1.e3 * avog/mw_s* 1.e-4
      sect_resident1 = sect_resident1 * 1.e3 * avog/mw_s* 1.e-4
      sect_shipping1 = sect_shipping1 * 1.e3 * avog/mw_s* 1.e-4

   else
      ; apply unit conversion for each species in the list separately assuming mw in file is correct
      ; since not correct, read correct mw:

      ; mw_s = mw_spec
      mw_s = mw_spec_set(ced)

      sect_other1_ced=sect_agriculture
      sect_other1_ced=sf*(sect_agriculture+sect_agriculture2+sect_agriculture3+sect_solvents+sect_waste)
      sect_energy1_ced=sect_energy
      sect_energy1_ced=sf*sect_energy
      sect_industry1_ced=sect_industry
      sect_industry1_ced=sf*(sect_industry+sect_industry2)
      sect_transport1_ced=sect_transport
      sect_transport1_ced=sf*(sect_transport+sect_transport2)
      sect_resident1_ced=sect_resident
      sect_resident1_ced=sf*sect_resident
      sect_shipping1_ced=sect_shipping
      sect_shipping1_ced=sf*sect_shipping
      
      sect_other1_ced= sect_other1_ced* 1.e3 * avog/mw_s* 1.e-4
      sect_energy1_ced = sect_energy1_ced * 1.e3 * avog/mw_s* 1.e-4
      sect_industry1_ced = sect_industry1_ced* 1.e3 * avog/mw_s* 1.e-4
      sect_transport1_ced = sect_transport1_ced * 1.e3 * avog/mw_s* 1.e-4
      sect_resident1_ced = sect_resident1_ced * 1.e3 * avog/mw_s* 1.e-4
      sect_shipping1_ced = sect_shipping1_ced * 1.e3 * avog/mw_s* 1.e-4


      sect_other1=sect_other1 + sect_other1_ced 
      sect_energy1=sect_energy1 +  sect_energy1_ced
      sect_industry1=sect_industry1 + sect_industry1_ced
      sect_transport1=sect_transport1 + sect_transport1_ced
      sect_resident1=sect_resident1 + sect_resident1_ced
      sect_shipping1=sect_shipping1 + sect_shipping1_ced
   end if  
     
 end do ;ced

      other_out1 = sect_other1
      transport_out1 = sect_transport1
      energy_out1 = sect_energy1
      industry_out1 = sect_industry1
      resident_out1 = sect_resident1
      shipping_out1 = sect_shipping1


;     sector = (/"other","energy","industry","transport","resident","shipping"/)

      other_out  = new((/ntim,ncol/),float)
      other_out!0 = "time"
      other_out!1 = "ncol"
      other_out&time = time
      other_out = (/doubletofloat(other_out1)/)
      other_out@units = "molecules/cm2/s"
      other_out@molecular_weight=mw
      other_out@molecular_weight_units = "g/mole"
      other_out@short_name=other_out1@short_name
      other_out@sector="agriculture,solvents,waste"
      printVarSummary(other_out)


      energy_out  = new((/ntim,ncol/),float)
      energy_out!0 = "time"
      energy_out!1 = "ncol"
      energy_out&time = time
      energy_out@units = "molecules/cm2/s"
      energy_out@molecular_weight=mw
      energy_out@molecular_weight_units = "g/mole"
      energy_out@short_name=energy_out1@short_name
      energy_out@sector="energy"
      energy_out = (/doubletofloat(energy_out1)/)

      industry_out  = new((/ntim,ncol/),float)
      industry_out!0 = "time"
      industry_out!1 = "ncol"
      industry_out&time = time
      industry_out@units = "molecules/cm2/s"
      industry_out@molecular_weight=mw
      industry_out@molecular_weight_units = "g/mole"
      industry_out@short_name=industry_out1@short_name
      industry_out@sector="industry"
      industry_out = (/doubletofloat(industry_out1)/)
     
      resident_out  = new((/ntim,ncol/),float)
      resident_out!0 = "time"
      resident_out!1 = "ncol"
      resident_out&time = time
      resident_out@units = "molecules/cm2/s"
      resident_out@molecular_weight=mw
      resident_out@molecular_weight_units = "g/mole"
      resident_out@short_name=resident_out1@short_name
      resident_out@sector="resident"
      resident_out = (/doubletofloat(resident_out1)/)


      transport_out  = new((/ntim,ncol/),float)
      transport_out!0 = "time"
      transport_out!1 = "ncol"
      transport_out&time = time
      transport_out@units = "molecules/cm2/s"
      transport_out@molecular_weight=mw
      transport_out@molecular_weight_units = "g/mole"
      transport_out@short_name=transport_out1@short_name
      transport_out@sector="transport"
      transport_out = (/doubletofloat(transport_out1)/)

      shipping_out  = new((/ntim,ncol/),float)
      shipping_out!0 = "time"
      shipping_out!1 = "ncol"
      shipping_out&time = time
      shipping_out@units = "molecules/cm2/s"
      shipping_out@molecular_weight=mw
      shipping_out@molecular_weight_units = "g/mole"
      shipping_out@short_name=shipping_out1@short_name
      shipping_out@sector="shipping"
      shipping_out = (/doubletofloat(shipping_out1)/)



      

      ; output to netcdf file
;
if nspc .le. 5  then
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name2 + "_"+version+ "_c" + tostring(creation_date) + ".nc"
else
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name + "_"+version+ "_c" + tostring(creation_date) + ".nc"
end if
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

printVarSummary(other_out)
print(min(other_out))
print(max(other_out))
print(min(energy_out))
print(max(energy_out))
      ncdf->emiss_other= other_out
      ncdf->emiss_ene= energy_out
      ncdf->emiss_ind= industry_out
      ncdf->emiss_res= resident_out
      ncdf->emiss_ship= shipping_out
      ncdf->emiss_tran= transport_out
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;

      
    if spc_name .eq. "BC" .or. spc_name .eq. "OC" then
        PI  = get_pi("float" )
   ; calculate number emissions
        diam = 0.134e-6
        if spc_name .eq. "BC" then
           rho = 1700.
        else
          rho = 1000.
        end if
        mw = 12.
        mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
        other_out = other_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        energy_out = energy_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        industry_out = industry_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        resident_out = resident_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        shipping_out = shipping_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        transport_out = transport_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
   
      other_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      energy_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      industry_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      shipping_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      transport_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"
      resident_out@units = "(particles/cm2/s)(molecules/mole)(g/kg)"

        if spc_name .eq. "BC" then
         other_out@long_name =  "number of BC Emissions"
         energy_out@long_name =  "number of BC Emissions"
         industry_out@long_name =  "number of BC Emissions"
         shipping_out@long_name =  "number of BC Emissions"
         transport_out@long_name =  "number of BC Emissions"
         resident_out@long_name =  "number of BC Emissions"
        end if
        if spc_name .eq. "OC" then
         other_out@long_name =  "number of (1.4*OC) POM Emissions"
         energy_out@long_name =  "number of (1.4*OC) POM Emissions"
         industry_out@long_name =  "number of (1.4*OC) POM Emissions"
         shipping_out@long_name =  "number of (1.4*OC) POM Emissions"
         transport_out@long_name =  "number of (1.4*OC) POM Emissions"
         resident_out@long_name =  "number of (1.4*OC) POM Emissions"
        end if

; output to netcdf file
;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_num_"+spc_name2 + "_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of number "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->num_emiss_other= other_out
      ncdf->num_emiss_ene= energy_out
      ncdf->num_emiss_ind= industry_out
      ncdf->num_emiss_res= resident_out
      ncdf->num_emiss_ship= shipping_out
      ncdf->num_emiss_tran= transport_out
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)
;;;;;;;;;;;

   end if ; bc and pom 


     if spc_name .eq. "SO2" then
      mw_so2 = 64.
      mw_so4_bam = 96.
      mw_so4_mam = 115.
      other_out_so2 = other_out
      other_out_so4 = other_out
      other_out_so2 = 0.975*other_out_so2
      other_out_so4 = 0.025*other_out_so4
      energy_out_so2 = energy_out
      energy_out_so4 = energy_out
      energy_out_so2 = 0.975*energy_out_so2
      energy_out_so4 = 0.025*energy_out_so4
      industry_out_so2 = industry_out
      industry_out_so4 = industry_out
      industry_out_so2 = 0.975*industry_out_so2
      industry_out_so4 = 0.025*industry_out_so4
      resident_out_so2 = resident_out 
      resident_out_so4 = resident_out 
      resident_out_so2 = 0.975*resident_out_so2
      resident_out_so4 = 0.025*resident_out_so4
      shipping_out_so2 = shipping_out
      shipping_out_so4 = shipping_out
      shipping_out_so2 = 0.975*shipping_out_so2
      shipping_out_so4 = 0.025*shipping_out_so4
      transport_out_so2 = transport_out
      transport_out_so4 = transport_out
      transport_out_so2 = 0.975*transport_out_so2
      transport_out_so4 = 0.025*transport_out_so4

      other_out_so2@long_name="0.975*SO2 Other Emissions" ; surface
      other_out_so4@long_name="0.025*SO2 Other Emissions" ; surface
      industry_out_so2@long_name="0.975*SO2 Industry Emissions" ; surface
      industry_out_so4@long_name="0.025*SO2 Industry Emissions" ; surface
      energy_out_so2@long_name="0.975*SO2 Energy Emissions" ; surface
      energy_out_so4@long_name="0.025*SO2 Energy Emissions" ; surface
      resident_out_so2@long_name="0.975*SO2 Resident Emissions" ; surface
      resident_out_so4@long_name="0.025*SO2 Resident Emissions" ; surface
      shipping_out_so2@long_name="0.975*SO2 Shipping Emissions" ; surface
      shipping_out_so4@long_name="0.025*SO2 Shipping Emissions" ; surface
      transport_out_so2@long_name="0.975*SO2 Transport Emissions" ; surface
      transport_out_so4@long_name="0.025*SO2 Transport Emissions" ; surface
      
       ; write so4_a1: sulfate accumulation mode for energy+industrial
       ; vertically distributed: 100-300m
       nalt = 8
       dz = 0.05 ;km
       dalt = 2.e4 ;cm = 200m
       altitude = fspan(0,nalt-1,nalt) * dz + 0.5*dz
       altitude_int = fspan(0,nalt,nalt+1) * dz
       altitude!0="altitude"
       altitude_int!0="altitude_int"
       altitude@units="km"
       altitude@long_name="Altitude"
       altitude_int@units="km"
       altitude_int@long_name="Altitude interfaces"

       energy_out_so4_vert = new((/ntim,nalt,ncol/),"float")
       industry_out_so4_vert = new((/ntim,nalt,ncol/),"float")
       energy_out_so4_vert(:,:,:) = 0.
       industry_out_so4_vert(:,:,:) = 0.
       energy_out_so4_vol = energy_out_so4
       energy_out_so4_vol = energy_out_so4/dalt
       industry_out_so4_vol = industry_out_so4
       industry_out_so4_vol = industry_out_so4/dalt
       do ialt=2,5
         energy_out_so4_vert(:,ialt,:) = energy_out_so4_vol
         industry_out_so4_vert(:,ialt,:) = industry_out_so4_vol
       end do
       energy_out_so4_vert!0="time"
       energy_out_so4_vert!1="altitude"
       energy_out_so4_vert!2="ncol"
       energy_out_so4_vert@units="kg m-3 s-1"
       energy_out_so4_vert@long_name="0.025*SO2 Anthropogenic Emissions from Energy"
       energy_out_so4_vert&altitude=altitude

       industry_out_so4_vert!0="time"
       industry_out_so4_vert!1="altitude"
       industry_out_so4_vert!2="ncol"
       industry_out_so4_vert@units="kg m-3 s-1"
       industry_out_so4_vert@long_name="0.025*SO2 Anthropogenic Emissions from Industry"
       industry_out_so4_vert&altitude=altitude

printVarSummary(energy_out_so4_vert)
printVarSummary(industry_out_so4_vert)

; spec_name_sec = (/"SO2_anthro-ag-ship-res_surface","SO2_anthro-ene_surface","so4_a1_anthro-ene_vertical","so4_a1_anthro-ag-ship_surface","so4_a2_anthro-res_surface","num_so4_a1_anthro-ene_vertical","num_so4_a1_anthro-ag-ship_surface","num_so4_a2_anthro-res_surface"/)
      printVarSummary(time)
      

  
   ; calculate number emission for so4 
   ;  sector = (/"other","energy","industry","transport","resident","shipping"/)
     ;varname = 'emiss_ag_sol_was'
        diam = 0.134e-6
        rho = 1770.
        mw2 = 115.
        mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
        num_other_out_so4 = other_out_so4 
        num_other_out_so4= num_other_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_other_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_other_out_so4@long_name = "particle number emissions of SO2*0.025 "

    ; if sec.eq.1 .or. sec.eq.2 .or.  sec.eq.5 then
    ; varname = 'emiss_shipping'
    ;varname = 'emiss_ene_ind'
        diam = 0.261e-6
        rho = 1770.
        mw2 = 115.
        num_energy_out_so4_vert = energy_out_so4_vert
        num_energy_out_so4_vert= num_energy_out_so4_vert*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_energy_out_so4_vert@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_energy_out_so4_vert@long_name = "particle number emissions of SO2*0.025 "
        num_industry_out_so4_vert = industry_out_so4_vert
        num_industry_out_so4_vert= num_industry_out_so4_vert*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_industry_out_so4_vert@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_industry_out_so4_vert@long_name = "particle number emissions of SO2*0.025 "
        num_shipping_out_so4 = shipping_out_so4 
        num_shipping_out_so4= num_shipping_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_shipping_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_shipping_out_so4@long_name = "particle number emissions of SO2*0.025 "
     ;end if
     ; if sec.eq.3 .or. sec.eq.4 then
     ;varname = 'emiss_res_tran'
        diam = 0.0504e-6
        rho = 1770.
        mw2 = 115.
        num_resident_out_so4= resident_out_so4
        num_resident_out_so4= num_resident_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_resident_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_resident_out_so4@long_name = "particle number emissions of SO2*0.025 "
        num_transport_out_so4= transport_out_so4
        num_transport_out_so4= num_transport_out_so4*mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
        num_transport_out_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
        num_transport_out_so4@long_name = "particle number emissions of SO2*0.025 "
     ; end if


      ; output to netcdf file
; spec_name_sec = (/"SO2_anthro-ag-ship-res_surface","SO2_anthro-ene_surface","so4_a1_anthro-ene_vertical","so4_a1_anthro-ag-ship_surface","so4_a2_anthro-res_surface","num_so4_a1_anthro-ene_vertical","num_so4_a1_anthro-ag-ship_surface","num_so4_a2_anthro-res_surface"/)
;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_SO2_anthro_ag-ship-res-surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so2
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_other= other_out_so2
      ncdf->emiss_res= resident_out_so2
      ncdf->emiss_tran= transport_out_so2
      ncdf->emiss_ship= shipping_out_so2
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_SO2_anthro_ene-surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so2
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_ene= energy_out_so2
      ncdf->emiss_ind= industry_out_so2
      ncdf->date         = date
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_so4_a1_anthro-ag-ship_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_other= other_out_so4
      ncdf->emiss_ship= shipping_out_so4
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a1_anthro-ag-ship_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_other= num_other_out_so4
      ncdf->emiss_ship= num_shipping_out_so4
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_so4_a2_anthro-res_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_tran= transport_out_so4
      ncdf->emiss_res= resident_out_so4
      ncdf->date         = date
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a2_anthro-res_surface_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_tran= num_transport_out_so4
      ncdf->emiss_res= num_resident_out_so4
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;



           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_so4_a1_anthro-ene_vertical_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_ene= energy_out_so4_vert
      ncdf->emiss_ind= industry_out_so4_vert
      ncdf->altitude = altitude
      ncdf->altitude_int = altitude_int
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;

           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a1_anthro-ene_vertical_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)

      ncdf->emiss_ene= energy_out_so4_vert
      ncdf->emiss_ind= industry_out_so4_vert
      ncdf->altitude = altitude
      ncdf->altitude_int = altitude_int
      ncdf->date         = date
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
      ncdf->time         = (/time/)
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

;;;;;;;;;;;

           file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_num_so4_a1_anthro-ene_vertical_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_2020_sector.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      file_att@molecular_weight = mw_so4_mam
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)

      ncdf->emiss_ene= num_energy_out_so4_vert
      ncdf->emiss_ind= num_industry_out_so4_vert
      ncdf->date         = date
      ncdf->time         = (/time/)
      ncdf->lat         = lat_in
      ncdf->lon         = lon_in
      ncdf->area=area
;     ncdf->time_bnds    = time_bnds
      ncdf->altitude = altitude
      ncdf->altitude_int = altitude_int
      delete(ncdf)

;;;;;;;;;;;


   end if  ; SO2
;

;



;
; clean up
;
print("clean up")
      ;delete(wk)
      delete(in)
      delete(file_inp)
      delete(time)
      ;delete(time_bnds)
      delete(utc_date)
      delete(date)
;    delete(var_name)
    delete(mw)
    delete(spc_name)
if nspec1 .gt. 1 then
    delete(mw_spec_set)
end if
; 
   ; end loop over sectors 
; end loop over species
;
 else 
  print("file does not exist")
 end if ; file exists

else 

; npec_name = "SVOC"
 if spc_name.eq. "SVOC"
   specs_in = "pom_a4"
   hist_spec = "SVOC=0.6*pom_a4"
   fact = 0.6
  end if
 ; read required files
   cnt_s = 0
   nsp = dimsizes(specs_in)
   do sp = 0,nsp-1
    files = systemfunc("ls " + dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+specs_in(sp)+"_"+version+ "_c*.nc")
    print(files)
    files1=files(0)
    in = addfile(files1,"r")

    if cnt_s.eq.0 then
     time1 = in->time
      date = in->date
      lat = in->lat
      lon = in->lon
      area = in->area
      emiss_o = in->emiss_other
      emiss_e = in->emiss_ene
      emiss_i = in->emiss_ind
      emiss_t = in->emiss_tran
      emiss_r = in->emiss_res
      emiss_s = in->emiss_ship

      mw_hc =  emiss_o@molecular_weight
     if spc_name.eq."SOAGx1.5" then
      emiss_o = yield(sp)*emiss_o*mw_hc/mw
      emiss_e = yield(sp)*emiss_e*mw_hc/mw
      emiss_i = yield(sp)*emiss_i*mw_hc/mw
      emiss_t = yield(sp)*emiss_t*mw_hc/mw
      emiss_r = yield(sp)*emiss_r*mw_hc/mw
      emiss_s = yield(sp)*emiss_s*mw_hc/mw
     else
      emiss_o = emiss_o *mw_hc/mw
      emiss_e = emiss_e *mw_hc/mw
      emiss_i = emiss_i *mw_hc/mw
      emiss_t = emiss_t *mw_hc/mw
      emiss_r = emiss_r *mw_hc/mw
      emiss_s = emiss_s *mw_hc/mw
     end if
    else
     emiss_o_hc1 = in->emiss_other
      emiss_e_hc1 = in->emiss_ene
      emiss_i_hc1 = in->emiss_ind
      emiss_t_hc1 = in->emiss_tran
      emiss_r_hc1 = in->emiss_res
      emiss_s_hc1 = in->emiss_ship
      mw_hc =  emiss_o_hc1@molecular_weight
     if spc_name.eq."SOAGx1.5" then
       emiss_o = emiss_o + yield(sp)*emiss_o_hc1*mw_hc/mw
       emiss_e = emiss_e + yield(sp)*emiss_e_hc1*mw_hc/mw
       emiss_i = emiss_i + yield(sp)*emiss_i_hc1*mw_hc/mw
       emiss_t = emiss_t + yield(sp)*emiss_t_hc1*mw_hc/mw
       emiss_r = emiss_r + yield(sp)*emiss_r_hc1*mw_hc/mw
       emiss_s = emiss_s + yield(sp)*emiss_s_hc1*mw_hc/mw
     else
       emiss_o = emiss_o + emiss_o_hc1*mw_hc/mw
       emiss_e = emiss_e + emiss_e_hc1*mw_hc/mw
       emiss_i = emiss_i + emiss_i_hc1*mw_hc/mw
       emiss_t = emiss_t + emiss_t_hc1*mw_hc/mw
       emiss_r = emiss_r + emiss_r_hc1*mw_hc/mw
       emiss_s = emiss_s + emiss_s_hc1*mw_hc/mw
     end if
      delete(emiss_o_hc1)
      delete(emiss_e_hc1)
      delete(emiss_i_hc1)
      delete(emiss_t_hc1)
      delete(emiss_r_hc1)
      delete(emiss_s_hc1)
    end if
    cnt_s = cnt_s+1
   end do  ; sp

    ;IVOC = 0.2*(HCs)
    ;SVOC = 0.6*pom_a4
      emiss_o = fact*emiss_o
      emiss_e = fact*emiss_e
      emiss_i = fact*emiss_i
      emiss_t = fact*emiss_t
      emiss_r = fact*emiss_r
      emiss_s = fact*emiss_s
 

      emiss_other = emiss_o
      emiss_other@history= hist_spec
      emiss_other@long_name= spc_name+" CAMS emissions"
      emiss_ene = emiss_e
      emiss_ene@history= hist_spec
      emiss_ene@long_name= spc_name+" CAMS emissions"
      emiss_ind = emiss_i
      emiss_ind@history= hist_spec
      emiss_ind@long_name= spc_name+" CAMS emissions"
      emiss_tran = emiss_t
      emiss_tran@history= hist_spec
      emiss_tran@long_name= spc_name+" CAMS emissions"
      emiss_res = emiss_r
      emiss_res@history= hist_spec
      emiss_res@long_name= spc_name+" CAMS emissions"
      emiss_ship  = emiss_s
      emiss_ship@history= hist_spec
      emiss_ship@long_name= spc_name+" CAMS emissions"

      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name + "_"+version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")

      file_att               = True
      file_att@molecular_weight = mw
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@creation_date = creation_date
      file_att@author        = "Simone Tilmes (tilmes@ucar.edu)"
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_se.ncl"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
      ncdf->emiss_other= emiss_other
      ncdf->emiss_ene= emiss_ene
      ncdf->emiss_ind= emiss_ind
      ncdf->emiss_res= emiss_res
      ncdf->emiss_ship= emiss_ship
      ncdf->emiss_tran= emiss_tran

      ncdf->date         = date
      ncdf->time         = time1
      ncdf->lat          = lat
      ncdf->lon          = lon
      ncdf->area = area
      delete(ncdf)
      delete(time1)
      delete(date)
      delete(specs_in)

  if spc_name.eq."SOAGx1.5" then
   delete(yield)
  end if
 ; end do; SOAG, IVOC, SVOC
 
 end if ;  spc < 27
  delete(spc_name1)
 end do  ;  species
;
end
