; rename_cams_anthro_se.ncl

begin

; Create CESM-ready emissions from SE grid CAMS emissions
;
;   transform regridded CAMS anthropogenic emissions into CESM species
;   run after running Regrid_fv2se_cams_anthro.ncl
;
; usage: ncl rename_cams_anthro_se.ncl

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; USER CHANGES START
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Location of files created with Regrid_fv2se_cams_anthro.ncl
dir_inp = "/glade/scratch/$USER/emis/cams_ne30np4/"
;
; version of CAMS files (in filename)
version= "_v3.1"
;
; model resolution (as written in filenames)
;resolution = "ne0CONUSne30x8"
resolution = "ne30np4"
;resolution = "ne30pg3"
;
; Directory for new files
dir_out = "/glade/scratch/$USER/emis/cesm_ne30np4/"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; USER CHANGES END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print(dir_out)
print(dir_inp)
;
; constants
;
  PI  = get_pi("float" )
  avog = 6.022e23
  earth_rad  = 6.371e8 ; in cm
  earth_area = 4. * acos(-1.) * earth_rad * earth_rad
;
; define creation date
;
  creation_date = systemfunc ("date +%Y%m%d")
;
; make sure output dir exists
;
  system("mkdir -p " + dir_out)
;
; loop over species for CESM (MOZART-TS1)
;
;  do nspc=0,33
  do nspc=1,1
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
      mw = 28.
      sf = 1.0
    end if
    if ( nspc.eq.2 ) then
      spc_name = "NH3"
      spc_name1 = "nh3"
      mw = 17.
      sf = 1.0
    end if
    if ( nspc.eq.3 ) then
      spc_name = "NOx"
      spc_name1 = "nox"
      mw = 46.
      sf = 1.0
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
      spc_name1 = (/"olef"/)
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
      spc_name1 = (/"xylenes","trimethylbenzenes","arom"/)
      mw = 106.
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
      spc_name1 = (/"aldehydes"/)
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
      spc_name1 = (/"isoprenes"/)
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
     spc_name = "IVOC"
     spc_name1 = "IVOC"
     sf = 1.
     mw = 184.
    end if
    if ( nspc.eq.29) then
     spc_name = "SVOC"
     spc_name1 = "SVOC"
     sf = 1.
     mw = 310.
    end if
    if ( nspc.eq.30) then
     spc_name = "SOAGx1.5"
     spc_name1 = "SOAGx1.5"
     sf = 1.5
     mw = 12.
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

; 
; Combine multiple CAMS species as needed to create CESM species
;
 if nspc .le. 27 then

  do ced  = 0, nspec1-1
;
; Get each CAMS file
;
    files = systemfunc("ls " + dir_inp + resolution + "_CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_"+spc_name1(ced) + version+".nc")
    nfils = dimsizes(files)
    print(files)
    ; should be just 1 file per CAMS species
    nfil=0
    file_inp = files(nfil)
    in = addfile(file_inp,"r")
;
; read coordinates
;
    lat_out   = in->lat
    lon_out   = in->lon
    time      = in->time
    area = in->area
    ntim = dimsizes(time)
;
    print("Processing " + file_inp)

   ; Use only 'sum' variable, not individual sectors
   var1 = in->sum
   fillv = var1@_FillValue
   var1=where(var1.eq.fillv,0.,var1)
   delete(fillv)
   emiss_inp = var1
   mw_spec = var1@molecular_weight
   printVarSummary(emiss_inp)
;
; convert from kg/m2/s to molec/cm2/s
;
;                            kg to g         m2 to cm2
;                               |                |
      emiss_inp = emiss_inp * 1.e3 * avog/mw_spec * 1.e-4

 if ced.eq.0 then
      emiss_out1 = emiss_inp
      emiss_out1 = sf * emiss_inp
     else
      emiss_out1 = emiss_out1 + sf * emiss_inp
     end if
  end do ;ced
;
; regrid to CESM lat-lon grid (not done)
;
;     emiss_out = area_conserve_remap_Wrap(lon_inp,lat_inp,emiss_out1,lon_out,lat_out,False)
;     delete(emiss_out1)
;     delete(emiss_out@units)
; no regridding
    emiss_out = emiss_out1 
    emiss_out@units = "molecules/cm2/s"
    emiss_out@molecular_weight = mw_spec
    emiss_out@molecular_weight_units = "g/mole"
    delete(emiss_out1)

;
; compute global average (in molec/cm2/s)
; annual average in Tg/year
;
     area_emis = conform(emiss_out,area,1)
     emiss_out_gl = emiss_out*area_emis
     emiss_out_gl1 = dim_sum_n_Wrap(emiss_out_gl,1) / dim_sum_n_Wrap(area_emis,1)

      dims = dimsizes(emiss_out_gl1)
      nyr = dims(0)/12
      emiss_out_gl_ann = emiss_out_gl1(0:nyr-1)
      do n=0,nyr-1
        n1 = n * 12
        n2 = n1 + 11
        emiss_out_gl_ann(n) = dim_avg_Wrap(emiss_out_gl1(n1:n2))
      end do
      emiss_out_gl_ann = emiss_out_gl_ann * earth_area * 86400. * 365 * mw * 1.e-12/avog
;
     print(emiss_out_gl_ann(0))
     delete(emiss_out_gl)
     delete(emiss_out_gl1)
     delete(emiss_out_gl_ann)
     delete(area_emis)


   if spc_name .eq. "SO2" then
      sect_agriculture  = in->agr
      sect_energy = in->ene 
      sect_industry = in->ind 
      sect_industry2 = in->fef 
      sect_solvents= in->slv
      sect_transport = in->tro 
      sect_transport2 = in->tnr 
      sect_resident = in->res 
      sect_agriculture2= in->mma
      sect_waste = in->swd 
      sect_shipping = in->shp 
      fillv = emiss_out@_FillValue

      sect_agriculture=where(sect_agriculture.eq.fillv,0.,sect_agriculture)
      sect_energy=where(sect_energy.eq.fillv,0.,sect_energy)
      sect_industry=where(sect_industry.eq.fillv,0.,sect_industry)
      sect_industry2=where(sect_industry2.eq.fillv,0.,sect_industry2)
      sect_solvents=where(sect_solvents.eq.fillv,0.,sect_solvents)
      sect_transport=where(sect_transport.eq.fillv,0.,sect_transport)
      sect_transport2=where(sect_transport2.eq.fillv,0.,sect_transport2)
      sect_resident=where(sect_resident.eq.fillv,0.,sect_resident)
      sect_agriculture2=where(sect_agriculture2.eq.fillv,0.,sect_agriculture2)
      sect_waste=where(sect_waste.eq.fillv,0.,sect_waste)
      sect_shipping=where(sect_shipping.eq.fillv,0.,sect_shipping)


      sect_industry=sect_industry+sect_industry2
      sect_transport=sect_transport+sect_transport2
      sect_agriculture=sect_agriculture+sect_agriculture2
     
      sect_agriculture = sect_agriculture * 1.e3 * avog/mw_spec * 1.e-4
      sect_energy = sect_energy * 1.e3 * avog/mw_spec * 1.e-4
      sect_industry = sect_industry * 1.e3 * avog/mw_spec * 1.e-4
      sect_transport = sect_transport * 1.e3 * avog/mw_spec * 1.e-4
      sect_resident = sect_resident * 1.e3 * avog/mw_spec * 1.e-4
      sect_solvents = sect_solvents * 1.e3 * avog/mw_spec * 1.e-4
      sect_waste = sect_waste * 1.e3 * avog/mw_spec * 1.e-4
      sect_shipping = sect_shipping * 1.e3 * avog/mw_spec * 1.e-4
      mw_so2 = 64.
      mw_so4_bam = 96.
      mw_so4_mam = 115.
      ag_sol_was_so2_out = sect_agriculture
      ag_sol_was_so4_out = sect_agriculture
      ag_sol_was_so2_out = (sect_agriculture+sect_solvents+sect_waste)*0.975
      ag_sol_was_so4_out = (sect_agriculture+sect_solvents+sect_waste)*0.025
      ag_sol_was_so4_out@units="molecules/cm2/s"
      ship_so2_out = sect_shipping
      ship_so4_out = sect_shipping
      ship_so2_out = sect_shipping*0.975
      ship_so4_out = sect_shipping*0.025
      ship_so4_out@units="molecules/cm2/s"
      res_tran_so2_out = sect_resident
      res_tran_so4_out = sect_resident
      res_tran_so2_out = (sect_resident+sect_transport)*0.975
      res_tran_so4_out = (sect_resident+sect_transport)*0.025
      res_tran_so4_out@units="molecules/cm2/s"
      ene_ind_so2_surf_out = sect_energy
      ene_ind_so2_surf_out = sect_energy+sect_industry
printVarSummary(ene_ind_so2_surf_out)
printVarSummary(res_tran_so2_out)
printVarSummary(ship_so2_out)

  ;   ag_sol_was_so2_out = area_conserve_remap_Wrap(lon_inp,lat_inp,ag_sol_was_so2,lon_out,lat_out,False)
  ;   ag_sol_was_so4_out = area_conserve_remap_Wrap(lon_inp,lat_inp,ag_sol_was_so4,lon_out,lat_out,False)
  ;   ship_so2_out = area_conserve_remap_Wrap(lon_inp,lat_inp,ship_so2,lon_out,lat_out,False)
  ;   ship_so4_out = area_conserve_remap_Wrap(lon_inp,lat_inp,ship_so4,lon_out,lat_out,False)
  ;   res_tran_so2_out = area_conserve_remap_Wrap(lon_inp,lat_inp,res_tran_so2,lon_out,lat_out,False)
  ;   res_tran_so4_out = area_conserve_remap_Wrap(lon_inp,lat_inp,res_tran_so4,lon_out,lat_out,False)
  ;   ene_ind_so2_surf_out = area_conserve_remap_Wrap(lon_inp,lat_inp,ene_ind_so2_surf,lon_out,lat_out,False)

      ; write so4_a1: sulfate accumulation mode for energy+industrial
  ; vertically distributed: 100-300m
     ncol = dimsizes(lat_out)
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
     ene_ind_so2_vert = new((/ntim,nalt,ncol/),"float")
     ene_ind_so2_vert(:,:,:)=0.
     ene_ind_so2_surf_vol = ene_ind_so2_surf_out/dalt
     ene_ind_so2_surf_vol=where(ene_ind_so2_surf_out.eq.0.,0.,ene_ind_so2_surf_vol)
     do ialt=3,6
       ene_ind_so2_vert(:,ialt,:) = ene_ind_so2_surf_vol
     end do
     ene_ind_so2_vert!0="time"
     ene_ind_so2_vert!1="altitude"
     ene_ind_so2_vert!2="ncol"
     ene_ind_so2_vert@units="molecules/cm3/s"
     ene_ind_so2_vert@long_name="SO2 Anthropogenic Emissions from Industry"
     ene_ind_so2_vert&time=ene_ind_so2_surf_out&time
     ene_ind_so2_vert&altitude=altitude

     ene_ind_so2_surf_out!0="time"
     ene_ind_so2_surf_out!1="ncol"
     ene_ind_so2_surf_out@units="molecules/cm2/s"
     ene_ind_so2_surf_out@long_name="SO2 Anthropogenic Emissions from Industry"
     ene_ind_so2_surf_out&time=emiss_out&time
     ene_ind_so4_surf_out = ene_ind_so2_surf_out
     ene_ind_so4_surf_out@long_name="SO4 Anthropogenic Emissions from Industry"

     ene_ind_so4_vert = ene_ind_so2_vert
     ene_ind_so4_vert@long_name="SO4 Anthropogenic Emissions from Industry"

      ene_ind_so2_surf_out= 0.975*(ene_ind_so2_surf_out)
      ene_ind_so4_surf_out = 0.025*(ene_ind_so2_surf_out)
      ene_ind_so2_vert= 0.975*(ene_ind_so2_vert)
      ene_ind_so4_vert = 0.025*(ene_ind_so2_vert)
printVarSummary(ene_ind_so4_surf_out)
printVarSummary(ene_ind_so4_vert)
  
   ; calculate number emissions
  ;varname = 'emiss_ag_sol_was'
      diam = 0.134e-6
      rho = 1770.
      mw = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_ag_sol_was_so4 = ag_sol_was_so4_out
      num_ag_sol_was_so4 = ag_sol_was_so4_out *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ag_sol_was_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
   ; varname = 'emiss_shipping'
      diam = 0.261e-6
      rho = 1770.
      mw = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_ship_so4 = ship_so4_out
      num_ship_so4 = ship_so4_out *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ship_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
  ;varname = 'emiss_res_tran'
      diam = 0.0504e-6
      rho = 1770.
      mw = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_res_tran_so4 = res_tran_so4_out
      num_res_tran_so4 = res_tran_so4_out *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_res_tran_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
  ;varname = 'emiss_ene_ind'
      diam = 0.261e-6
      rho = 1770.
      mw = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_ene_ind_so4_surf = ene_ind_so4_surf_out
      num_ene_ind_so4_vert = ene_ind_so4_vert
      num_ene_ind_so4_surf = ene_ind_so4_surf_out *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ene_ind_so4_vert = ene_ind_so4_vert *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ene_ind_so4_vert@units="(particles/cm2/s)(molecules/mole)(g/kg)"
end if

;
   if spc_name .eq. "BC" .or. spc_name .eq. "OC" then
   sec = 0
   PI  = get_pi("float" )
   ; calculate number emissions
      diam = 0.134e-6
     if spc_name .eq. "BC" then
      rho = 1700.
     else
      rho = 1000.
     end if
      emis_num = emiss_out
      mw = 12.
      sector = (/"Anthropogenic","Biomass Burning","Aviation"/)
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      emis_num = emiss_out*mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      emis_num@units="(particles/cm2/s)(molecules/mole)(g/kg)"
      emis_num@long_name = "number of BC "+sector(sec)+" Emissions"
      if spc_name .eq. "OC" then
       emiss_out@long_name = "(1.4*BC) POM "+sector(sec)+" Emissions"
       emis_num@long_name = "number of (1.4*BC) POM "+sector(sec)+" Emissions"
      end if
   end if

        emiss_out@long_name= spc_name+" CAMS Global anthropogenic emissions"
     if nspc.gt.5 then
        emiss_out@history= hist_spec
      end if

;
; extract time range string from input file name
;
  ;   length = strlen(file_inp)
  ;   char   = stringtochar(file_inp)
  ;   time_range = char(length-16:length-4)
  ;   delete(char)
;print(time_range)
;
;
; output to netcdf file
;
;      file_out = dir_out + "emissions-cams-"+version+"_" + spc_name + "_anthro_surface_mol_" + time_range + "_" + resolution + "_c" + tostring(creation_date) + ".nc"
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name + version+ "_c" + tostring(creation_date) + ".nc"
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
      file_att@data_creator= $USER  (email_user) 
      file_att@cesm_contact= $USER  (email_user)               
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_se.ncl"
      file_att@data_summary = "CAMS Emissions have been manipulated for use in CESM2" 
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)" 
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
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
;
      ncdf->emiss_anthro = emiss_out
      ncdf->date         = date
      ncdf->time         = time
      ncdf->lat          = lat_out
      ncdf->lon          = lon_out
      ncdf->area         = area
;     ncdf->time_bnds    = time_bnds
      delete(ncdf)

   if spc_name .eq. "BC" .or. spc_name .eq. "OC" then
  ;
; output to netcdf file
;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name2 + version+ "_c" + tostring(creation_date) + ".nc"
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
      file_att@data_creator= $USER  (email_user)
      file_att@cesm_contact= $USER  (email_user)
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_se.ncl"
      file_att@data_summary = "CAMS Emissions have been manipulated for use in CESM2"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)" 
      file_att@data_title = "Number of "+sector(sec)+" emissions of "+spc_name
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
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
;
      ncdf->emiss_anthro= emiss_out
      ncdf->date         = date
      ncdf->time         = time
      ncdf->lat          = lat_out
      ncdf->lon          = lon_out
      ncdf->area = area
      ;ncdf->time_bnds    = time_bnds
;
; clean up
      delete(ncdf)

   ; output to netcdf num file
;
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_num_"+spc_name2 + version+ "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_creator= $USER  (email_user)
      file_att@cesm_contact= $USER  (email_user)
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_se.ncl"
      file_att@data_summary = "CAMS Emissions have been manipulated for use in CESM2"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)" 
      file_att@data_title = "Number of "+sector(sec)+" emissions of "+spc_name2
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
;
      ncdf->emiss_anthro= emis_num
      ncdf->date         = date
      ncdf->time         = time
      ncdf->lat          = lat_out
      ncdf->lon          = lon_out
      ncdf->area = area
      ;ncdf->time_bnds    = time_bnds
;
; clean up
      delete(ncdf)
      delete(emis_num)

  end if ; bc and pom


  ; for SO2 add sector specific output
  if spc_name .eq. "SO2" then
 ;  print out additional files
  spec_name_sec = (/"SO2_anthro-ag-ship-res_surface","SO2_anthro-ene_surface","so4_a1_anthro-ene_vertical","so4_a1_anthro-ag-ship_surface","so4_a2_anthro-res_surface","num_so4_a1_anthro-ene_vertical","num_so4_a1_anthro-ag-ship_surface","num_so4_a2_anthro-res_surface"/)

;  s1var = (/2,5/)
;  do s2 = 0,1
;   s1=s1var(s2)
   do s1 = 0,7
     if s1.eq.0 then
       emiss_out1 =  ag_sol_was_so2_out
       emiss_out2 =  res_tran_so2_out
       emiss_out3 =  ship_so2_out
       emiss_out1@long_name = "SO2*0.975 anthro agriculture+waste emissions"
       emiss_out2@long_name = "SO2*0.975 anthro residential+transportation emissions"
       emiss_out3@long_name = "SO2*0.975 anthro shipping emissions"
     end if
     if s1.eq.1 then
       emiss_out1 = ene_ind_so2_surf_out
       emiss_out1@long_name = "SO2*0.975 anthro energy+industrial emissions"
     end if
     if s1.eq.2 then
       emiss_out1 = ene_ind_so4_vert
       emiss_out1@long_name = "SO2*0.025 anthro energy+industry emissions"
     end if
     if s1.eq.3 then
       emiss_out1 =  ag_sol_was_so4_out
       emiss_out2 =  ship_so4_out
       emiss_out1@long_name = "SO2*0.025 anthro agriculture+waste emissions"
       emiss_out2@long_name = "SO2*0.025 anthro shipping emissions"
     end if
     if s1.eq.4 then
       emiss_out1 =  res_tran_so4_out
       emiss_out1@long_name = "SO2*0.975 anthro residential+transportation emissions"
     end if
     if s1.eq.5 then
       emiss_out1 = num_ene_ind_so4_vert
       emiss_out1@long_name = "particle number emissions of SO2*0.025 anthro energy+industry"
     end if
     if s1.eq.6 then
       emiss_out1 =  num_ag_sol_was_so4
       emiss_out2 =  num_ship_so4
       emiss_out1@long_name = "particle number emissions of SO2*0.025 anthro agriculture+waste"
       emiss_out2@long_name = "particle number emissions of SO2*0.025 anthro shipping"
     end if
     if s1.eq.7 then
       emiss_out1 =  num_res_tran_so4
       emiss_out1@long_name = "particle number emissions of SO2*0.975 anthro residential+transportation"
     end if

print(spec_name_sec(s1))
      spec_name11 = spec_name_sec(s1)
      file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_"+spec_name_sec(s1)+ version+ "_c" + tostring(creation_date) + ".nc"
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
      file_att@data_title = "Emissions of "+spec_name11+" for CESM"
      file_att@data_creator= $USER  (email_user)
      file_att@cesm_contact= $USER  (email_user)
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_se.ncl"
      file_att@data_summary = "CAMS Emissions have been manipulated for use in CESM2"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)" 
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
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
;
      ncdf->emiss= emiss_out1
     if s1.eq.0 .or. s1.eq.3 .or. s1.eq.6 then
      ncdf->emiss2= emiss_out2
      delete(emiss_out2)
     end if
     if s1.eq.0 then
      ncdf->emiss3 = emiss_out3
      delete(emiss_out3)
     end if
      ncdf->date         = date
      ncdf->time         = time
      ncdf->lat          = lat_out
      ncdf->lon          = lon_out
      ncdf->area = area
      if s1.eq.2 .or. s1.eq.5 then
       ncdf->altitude = altitude
       ncdf->altitude_int = altitude_int
      end if
      ;ncdf->time_bnds    = time_bnds
      delete(ncdf)
      delete(emiss_out1)
  end do
;
  end if
;
; clean up
;
      ;delete(wk)
      delete(emiss_inp)
      delete(emiss_out)
      delete(in)
      delete(file_inp)
      delete(time)
      ;delete(time_bnds)
      delete(utc_date)
      delete(date)
;    delete(var_name)
    delete(mw)
    delete(spc_name)
    delete(spc_name1)
;
;  end loop over species
;
  end if ; spec = 1-27

  if nspc.gt.27 then
do sec=0,0
 if spc_name.eq. "IVOC"
   specs_in = (/"C3H6", "C3H8", "C2H6", "C2H4", "BIGENE", "BIGALK", "CH3COCH3", "MEK", "CH3CHO", "CH2O", "BENZENE", "TOLUENE", "XYLENES"/)
   hist_spec = "IVOC=0.2*("+str_join(specs_in,"+")+")"
   fact = 0.2
 end if
  if spc_name.eq. "SVOC"
   specs_in = "pom_a4"
   hist_spec = "SVOC=0.6*pom_a4"
   fact = 0.6
  end if
 if spc_name.eq. "SOAGx1.5"
; if sec.eq.0 then
;  specs_in  = (/"BIGALK","BIGENE","TOLUENE","BENZENE","XYLENES"/)
;  yield  = (/0.05, 0.05, 0.15, 0.15, 0.15/)
;  hist_spec = "SOAG=1.5* yields(0.05, 0.05, 0.15, 0.15, 0.15) from ("+str_join(specs_in,"+")+")"
; end if
; if sec.eq.1 then
   specs_in  = (/"BIGALK","BIGENE","TOLUENE","BENZENE","XYLENES","ISOP","MTERP"/)
   yield  = (/0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25/)
   hist_spec = "SOAG=1.5* yields(0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25) from ("+str_join(specs_in,"+")+")"
; end if
   fact = 1.5
 end if
 ; read required files
   cnt_s = 0
   nsp = dimsizes(specs_in)
   do sp = 0,nsp-1
    files = systemfunc("ls " + dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+specs_in(sp)+ version+ "_c*.nc")
     print(files)
    files1=files(0)
    in = addfile(files1,"r")
    if cnt_s.eq.0 then
      time1 = in->time
      date = in->date
      area = in->area
      lat = in->lat
      lon = in->lon
      wk = in->emiss_anthro
      printVarSummary(wk)
      mw_hc =  wk@molecular_weight
     if spc_name.eq."SOAGx1.5" then
      wk = yield(sp)*wk*mw_hc/mw
     else
      wk = wk*mw_hc/mw
     end if
    else
      emis_hc1 = in->emiss_anthro
      printVarSummary(emis_hc1)
      mw_hc =  emis_hc1@molecular_weight
      print(mw_hc)
      print(mw)
      printVarSummary(wk)
     if spc_name.eq."SOAGx1.5" then
      wk = wk + yield(sp)*emis_hc1*mw_hc/mw
     else
      wk = wk + emis_hc1*mw_hc/mw
     end if
      delete(emis_hc1)
    end if
    ;IVOC = 0.2*(HCs)
    ;SVOC = 0.6*pom_a4
     wk = fact*wk
    cnt_s = cnt_s+1
   end do
  emiss_out = wk
  emiss_out@history= hist_spec
  file_out = dir_out + "CAMS-GLOB-ANT_Glb_"+ resolution + "_anthro_surface_"+spc_name + version+ "_c" + tostring(creation_date) + ".nc"
  system("/bin/rm -f " + file_out) ; remove any pre-existing file
  ncdf = addfile(file_out,"c")
  printVarSummary(emiss_out)
; emiss_out@units="molecules/cm2/s"
  emiss_out@long_name= spc_name+" CAMS emissions"
  system("/bin/rm -f " + file_out) ; remove any pre-existing file
  ncdf = addfile(file_out,"c")

      file_att               = True
      file_att@molecular_weight = mw
      file_att@data_title = "Emissions of "+spc_name+" for CESM"
      file_att@creation_date = creation_date
      file_att@data_creator= $USER  (email_user)
      file_att@cesm_contact= $USER  (email_user)
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_cams/regrid_cams_anthro_se.ncl"
      file_att@data_summary = "CAMS Emissions have been manipulated for use in CESM2"
      file_att@data_source_url = "https://eccad3.sedoo.fr" ;
      file_att@data_reference = "Granier, C., et al, GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)" 
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
  ncdf->emiss= emiss_out
  ncdf->date         = date
  ncdf->time         = time1
  ncdf->lat          = lat
  ncdf->lon          = lon
  ncdf->area = area

  delete(ncdf)
  delete(time1)
  delete(date)
  delete(emiss_out)
  delete(wk)
  delete(specs_in)

 if spc_name.eq."SOAGx1.5" then
  delete(yield)
 end if
 end do
 end if  ; SOAG, IVOC, SVOC
;
  end do   
;
end
