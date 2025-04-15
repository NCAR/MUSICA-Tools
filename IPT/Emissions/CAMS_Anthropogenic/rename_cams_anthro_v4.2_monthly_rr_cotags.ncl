;
; read original CEDS/CMIP6 SSP emissions files
; regrid to 1 deg and convert from kg/m2/s to molecules/cm2/s
; sum all anthro sectors to one variable
;  (requires gnu compiler, not intel, on cheyenne/geyser)
;

begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Location of files created with Regrid_fv2se_cams_anthro.ncl
dir_inp = "/glade/p/acom/MUSICA/emissions/cams/CAMS-GLOB-ANTv4.2/ne30np4/"
file_co = dir_inp+"CAMS-GLOB-ANT_2000-2020_ne30np4_anthro_surface_CO_v4.2_c20210220.nc"
;
version= "v4.2"

; model resolution (as written in filenames)
;resolution = "ne0CONUSne30x8"
resolution = "ne30np4"
;resolution = "ne0np4_Asia03_ne30x8"
;resolution = "ne30pg3"
;
; Directory for new files
;dir_out = "/glade/scratch/$USER/emis/cesm_ne30np4/"
;dir_out = "/glade/p/acom/MUSICA/emissions/cams/CAMS-GLOB-ANTv4.2/ne30np4/"
;dir_out = "/glade/p/acom/MUSICA/emissions/cams/CAMS-GLOB-ANTv4.2/ne0np4.Asia03.ne30x8/"
dir_out = "/glade/p/acom/MUSICA/emissions/cams/CAMS-GLOB-ANTv4.2/ne30np4/"
;
creation_date = systemfunc ("date +%Y%m%d")

regions = (/"NAmerica","Europe","NAfrica","WAsia","EAsia","SAmerica","SAfrica","Australia"/)
lat_min = (/17.,30.,0.,0.,0.,-55.,-55.,-55./)
lat_max = (/89.9,89.9,30.,89.9,89.9,17.,0.,0./)
lon_min1 = (/181.,335.,335.,60.,100.,181.,335.,75./)
lon_max1 = (/335.,360.,360.,100.,180.,335.,360.,180./)
lon_min2 = (/181.,0.,    0.,60.,100.,181.,0.,75./)
lon_max2 = (/335.,60.,  60.,100.,180.,335.,75.,180./)
nre = dimsizes(regions)
re_num = sprinti("%0.2i", ispan(1,nre,1))
print(re_num)

do r=0,nre-1 
 file_cotag = dir_inp + "CAMS-GLOB-ANT_2000-2020_ne30np4_anthro_surface_CO"+re_num(r)+"_"+version+"_"+creation_date+".nc" 
print(file_cotag)
 command = "cp "+file_co+" "+file_cotag
 print((/command/))
 system(command)
 fc = addfile(file_cotag,"w")
 lat = fc->lat 
 lon = fc->lon 
 emiss_other=fc->emiss_other
 emiss_ene=fc->emiss_ene
 emiss_ind=fc->emiss_ind
 emiss_res=fc->emiss_res
 emiss_tran=fc->emiss_tran
 emiss_ship=fc->emiss_ship

 emiss_other1=emiss_other
 emiss_ene1=emiss_ene
 emiss_ind1=emiss_ind
 emiss_res1=emiss_res
 emiss_tran1=emiss_tran
 emiss_ship1=emiss_ship

 emiss_other1=0.              
 emiss_ene1=0.              
 emiss_ind1=0.              
 emiss_res1=0.              
 emiss_tran1=0.              
 emiss_ship1=0.              
 
 ind_co = ind((lat.ge.lat_min(r) .and. lat.le.lat_max(r) .and. lon.ge.lon_min1(r) .and. lon.le.lon_max1(r)) .or. (lat.ge.lat_min(r) .and. lat.le.lat_max(r) .and. lon.ge.lon_min2(r) .and. lon.le.lon_max2(r)))
 print(max(lat(ind_co)))
 print(min(lat(ind_co)))
 print(max(lon(ind_co)))
 print(min(lon(ind_co)))

 emiss_other1(:,ind_co) = emiss_other(:,ind_co) 
 emiss_ene1(:,ind_co) = emiss_ene(:,ind_co) 
 emiss_ind1(:,ind_co) = emiss_ind(:,ind_co) 
 emiss_res1(:,ind_co) = emiss_res(:,ind_co) 
 emiss_tran1(:,ind_co) = emiss_tran(:,ind_co) 
 emiss_ship1(:,ind_co) = emiss_ship(:,ind_co) 
 
 fc->emiss_other=emiss_other1
 fc->emiss_ene=emiss_ene1
 fc->emiss_ind=emiss_ind1
 fc->emiss_res=emiss_res1
 fc->emiss_tran=emiss_tran1
 fc->emiss_ship=emiss_ship1

 delete(ind_co)

end do

;
end
