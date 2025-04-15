; This scripts conservely regrids historical CMIP6 emissions from fv09 to other grid resolutions
; 
; Different emission sectors will be regridded using different programs and all are required for running CESM.
;
; Simone Tilmes, May 21st 2020, tilmes@ucar.edu
;

load "Regrid_fv2se_cmip6_hist_anthro.ncl"
load "Regrid_fv2se_cmip6_hist_bb.ncl"
load "Regrid_fv2se_cmip6_hist_voc.ncl"
load "Regrid_fv2se_cmip6_other.ncl"
load "Regrid_fv2se_cmip6_hist_air.ncl"
load "create_contvolc_files.ncl"

; these should not be changed
  ingrid   = "f09"
  interp   = "conserve"

; don't change srcPath 
  srcPath_cmip = "/glade/p/acom/acom-climate/cmip6inputs/emis/hist/"
  srcPath_other = "/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015/"
  srcPath_e90 = "/glade/p/acom/acom-climate/cmip6inputs/emis/hist/"
  srcPath_soag = "/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp/"

;;;;;;;;;;;;;;;;;;;;;;
; USER changes start 
; change <user> to your user name 
;;;;;;;;;;;;;;;;;;;;;;
; define/change grid information, destination location for the regridded files, and specifc period
; options for outgrid 

; PLEASE  add your own directory if you regrid to a new grid
  grid_dir ="/glade/work/tilmes/refgrids/"
  wght_dir ="/glade/work/tilmes/regrid/Grid_Weights/"
;
; these have to be adjusted to the grid you are going to interpolate to:
  resolution= "ne30np4"
  ;resolution = "ne30pg3"
  ;resolution = "conus_30_x8"
  ;resolution = "C96"
  dstPath = "/glade/scratch/<user>/emissions/historical_"+resolution"+/"
  dstPath_other = "/glade/scratch/<user>/cmip6inputs/historical_"+resolution+"/"

  system("mkdir "+dstPath)
  system("mkdir "+dstPath_other)
  period = "all"
  ;period = "200001-201412"

; add your resolution files here:
  outgrid = resolution
  if(outgrid.eq."ne30np4") then
    dstGridName = grid_dir+"ne30np4_091226_pentagons.nc"
    grid_file_volc = "/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/se/cami-mam3_0000-01-01_ne30np4_L30_c100618_conserve.nc"
  end if
  if(outgrid.eq."ne30pg3") then
    dstGridName = grid_dir+"ne30pg3.nc"
    grid_file_volc = "/glade/p/acom/acom-climate/cmip6inputs/historical_ne30pg3/emissions-cmip6_C2H4_other_surface_1750_2015_ne30_ne30_c20191022.nc"  ; CSLAM
  end if
  if(outgrid.eq."ne120") then
    dstGridName = grid_dir+"ne120np4_pentagons_100310.nc"
  end if
  if(outgrid.eq."conus_30_x8") then
    dstGridName = grid_dir+"conus_30_x8.g_scrip.nc"
    grid_file_volc = "/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/se/fchist_ne0conus30x8_L32_2012-11-01_c190906.nc"  ; se Conus
  end if
  if(outgrid.eq."C96") then
    dstGridName = grid_dir+"C96_SCRIP_desc.181018.nc"
    grid_file_volc = "/glade/p/acom/acom-climate/cmip6inputs/historical_C96/emissions-cmip6_CO_other_surface_1750_2015_C96_c20191202.nc" ; C96
  end if

;;;;;;;;;;;;;;;;;;;;;;
; USER changes done
;;;;;;;;;;;;;;;;;;;;;;


  srcGridName = wght_dir+"SCRIP_"+ingrid+".nc"

  ; Open the first input file and read in lat/lon values
  ;  to create ESMF grid info.
  ;---------------------------------------------------------
  srcGridName = wght_dir+"SCRIP_"+ingrid+".nc"

  if(isfilepresent(srcGridName)) then
    print(" Removing old ESMF weights:"+srcGridName)
    system("rm -f "+srcGridName)
  end if
  print(" Creating SCRIP file:"+srcGridName)


 ; Remove existing ESMF weights, create new values
  ;---------------------------------------------------------
  wgtGridName = wght_dir+"ESMFmap_"+ingrid+"_"+outgrid+"_"+interp+".nc"

  if(isfilepresent(wgtGridName)) then
;   print(" Use existing ESMF weights:"+wgtGridName)
; else
     print(" Removing old ESMF weights:"+wgtGridName)
     system("rm -f "+wgtGridName)
   end if
  print(" Creating ESMF weights:"+wgtGridName)
  print (" ")


;;;;;;;;;;;;;;;;;;;;;;
; you can comment out those files that you don't want to run 
;;;;;;;;;;;;;;;;;;;;;


  print("Regrid_fv2se_cmip6_hist_anthro")
  Regrid_fv2se_cmip6_hist_anthro(resolution,dstPath,period,dstGridName,srcGridName,wgtGridName)
  print("Regrid_fv2se_cmip6_hist_bb")
  Regrid_fv2se_cmip6_hist_bb(resolution,dstPath,period,dstGridName,srcGridName,wgtGridName)
  print("Regrid_fv2se_cmip6_hist_voc")
  Regrid_fv2se_cmip6_hist_voc(resolution,dstPath,period,dstGridName,srcGridName,wgtGridName)
  print("Regrid contvolc")
  create_contvolc_files(resolution,dstPath_other,grid_file_volc)
  print("Regrid_fv2se_cmip6_other")
  vari = "other"
  Regrid_fv2se_cmip6_other(resolution,srcPath_other,dstPath_other,period,dstGridName,srcGridName,wgtGridName,vari)
  print("Regrid_fv2se_cmip6_SOAG bio")
  vari = "SOAGx1.5_biogenic"
  Regrid_fv2se_cmip6_other(resolution,srcPath_soag,dstPath,period,dstGridName,srcGridName,wgtGridName,vari)
;; only if desired
; vari = "E90"
; Regrid_fv2se_cmip6_other(resolution,srcPath_other,dstPath_other,period,dstGridName,srcGridName,wgtGridName,vari)
; Regrid_fv2se_cmip6_hist_air(resolution,dstPath,period,dstGridName,srcGridName,wgtGridName)

