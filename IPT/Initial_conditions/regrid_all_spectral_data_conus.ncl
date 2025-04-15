load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
function open_netcdf_for_regridding(sfile,dfile,interp_method,src_fname)
local regrid_file, regrid_fname, regrid_dir

; PLEASE FIND 'USER CHANGES' to adjust the program for your case
; replace any <your_...> with the appropriate directory or file name
;
begin
; USER CHANGES

;   add your destination directory and file name where your regrid file will be written to  
    regrid_dir   = ""<your_VRM_Files>/<your_grid>/inic/"
    regrid_fname = regrid_dir + src_fname+"_ne0CONUSne30x8_mt12.nc"

; END USER CHANGES        

    setfileoption("nc","Format","NetCDF4Classic")         

    system("/bin/rm -f "+regrid_fname)        ; delete any pre-existing file
    rfile = addfile(regrid_fname, "c")  ; open for writing

;---Predefine mode
    setfileoption(rfile,"DefineMode",True)

;---Create global attributes of the netCDF file containing regridded data
    global    = True
    global@separator2= "---------------------------------------"
    copy_VarAtts(sfile, global)
    global@separator1= "------- SOURCE FILE ATTRIBUTES --------"

    if (isatt(sfile,"title")) then
        global@TITLE = "REMAPPED: "+sfile@title
    end if

    global@remap         = "NCL: ESMF_regrid_with_weights"
    global@map_method    = interp_method
    global@creation_date = systemfunc("date")
    
    fileattdef( rfile, global )        ; copy file attributes to netCDF file

;---predefine the coordinate dimension names and their sizes
    src_dnames = getvardims(sfile)     ; dimension names
    dst_dnames = getvardims(dfile)
    src_dsizes = getfiledimsizes(sfile)      ; dimension sizes
    dst_dsizes = getfiledimsizes(dfile)
    itime    = ind(src_dnames.eq."time" )      ; time dimension
    ilev     = ind(src_dnames.eq."lev" )       ; level dimension 
    icell    = ind(dst_dnames.eq."ncol" )      ; cell dimension
    ncell    = dst_dsizes(icell)               ; get # of cells from destination grid
    nlev     = src_dsizes(ilev)
    dimNames = (/"time",   "lev", "ncol" /)
    dimSizes = (/ -1   ,    nlev,  ncell /)
    dimUnlim = (/ True ,    False, False /)
    filedimdef(rfile,dimNames,dimSizes,dimUnlim)

;--- The following are explicitly added because they are 'special'
   dst_lat = dfile->lat
   dst_lon = dfile->lon

   filevardef(rfile, "lat", typeof(dst_lat), getvardims(dst_lat))
   filevarattdef(rfile,"lat" ,dst_lat)

   filevardef(rfile, "lon", typeof(dst_lon), getvardims(dst_lon))
   filevarattdef(rfile,"lon" ,dst_lon)           
   return(rfile)
end

;----------------------------------------------------------------------
; Get list of variables to regrid
;----------------------------------------------------------------------
function get_file_vars_to_regrid(sfile,rfile)
local nv, nvars, vdimnames, vrank, vnames, vtypes, valid_types, \
src_time_name, src_cell_name, imsg, inmsg
begin
  vnames = getfilevarnames(sfile)
  vtypes = getfilevartypes(sfile,vnames)
  nvars  = dimsizes(vnames)
  vnames_to_regrid = new(nvars,string)
  valid_types = (/"float","double"/)

;---Specify a list of 1D variables on the *source SE file* that should NOT be copied
  var_src_exclude = (/"area"/)

;---Get dimension names for the cell and time dimensions
  src_cell_name = getfilevardims(sfile,"lat")
  src_time_name = getfilevardims(sfile,"time")

;---Loop through each variable on the source 
  do nv=0,nvars-1
    vdimnames := getfilevardims(sfile,vnames(nv))
    vdimsizes := getfilevardimsizes(sfile,vnames(nv))
    vrank     = dimsizes(vdimnames)
    if(.not.any(vnames(nv).eq.var_src_exclude).and.any(vtypes(nv).eq.valid_types).and.\
       (vdimnames(0).eq.src_time_name.and.vdimnames(vrank-1).eq.src_cell_name)) then
      vnames_to_regrid(nv) = vnames(nv)    ; store so we can return it
      filevardef(rfile, vnames(nv), vtypes(nv), vdimnames)
      dumatts = new( 1, vtypes(nv))
      varatts = getfilevaratts(sfile, vnames(nv))
      if (.not.ismissing(varatts(0))) then
          natts = dimsizes(varatts)
          do na=0,natts-1
             dumatts@$varatts(na)$ = sfile->$vnames(nv)$@$varatts(na)$
          end do
          filevarattdef(rfile, vnames(nv) , dumatts)
      end if
      delete(varatts)   ; may change size next iteration
      delete(dumatts)
    end if
  end do
  imsg  = ind(ismissing(vnames_to_regrid))
  inmsg = ind(.not.ismissing(vnames_to_regrid))
  print("Variables we will NOT regrid : " + str_join(vnames(imsg),","))
  print("Variables we WILL regrid     : " + str_join(vnames(inmsg),","))

;
; Turn off define mode because we know after this point, we will
; not be writing any more metadata.
;
  setfileoption(rfile,"DefineMode",False)

  return(vnames_to_regrid(inmsg))
end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
;---Source lat/lon mesh and variable to regrid
; USER CHANGES        
; Change to the starting date you need.
; existing ne30_ne30 are locally on cheyenne at: /glade/p/acom/MUSICA/init/ne30_ne30/ 
  src_dir   = "/glade/p/acom/MUSICA/init/ne30_ne30/"
 ;change YYYY-MM to the year and months you want to select
  src_fname = "f.e22.FCnudged.ne30_ne30_mg17.release-cesm2.2.0_spinup.2010_2020.001.cam.i.<restart_date>-01-00000.nc"
; END USER CHANGES        

;---Source lat/lon mesh
  src_file = addfile(src_dir+src_fname,"r")
  src_lat = src_file->lat    
  src_lon = src_file->lon
  src_dims = dimsizes(src_lat)    ; 48602 cells


;---Destination lat/lon mesh
; USER CHANGES: already changed for CONUS      
  dst_dir   = "/glade/p/acom/MUSICA/init/ne0CONUSne30x8/atm/"
  dst_fname = "f.e22.FCcotagsNudged.ne0CONUSne30x8.cesm220.2012-01.cam.i.2012-02-01-00000.nc"
  dfile = addfile(dst_dir + dst_fname,"r")
; END USER CHANGES        

  dst_lat = dfile->lat
  dst_lon = dfile->lon
  dst_dims = dimsizes(dst_lat)    ; 174098 cells
  dst_area = dfile->area
  dst_lev = dfile->lev
  hyam = dfile->hyam
  hybm = dfile->hybm
  hyai = dfile->hyai
  hybi = dfile->hybi
  ;P0 = dfile->P0
  ilev = dfile->ilev
  time = dfile->time
  date = dfile->date
  datesec = dfile->datesec
  time_bnds= dfile->time_bnds
  ndbase= dfile->ndbase
  nsbase= dfile->nsbase
  nbdate= dfile->nbdate
  nbsec= dfile->nbsec
  mdt= dfile->mdt
  ndcur= dfile->ndcur
  nscur= dfile->nscur
  nsteph= dfile->nsteph
  date_written= dfile->date_written
  time_written= dfile->time_written

;make sure area does not contain any _FillValues
if(any(ismissing(dst_area))) then
    print("Your area contains some missing values. Check your file.")
    exit()
end if
;remove the FillValue on area as netcdf cannot write the file with it.
 delete(dst_area@_FillValue)

;---Choose regridding method and weights file name
  interp_method = "bilinear"
  weights_fname = "spectral_" + src_dims + "_cells_to_spectral_"+ dst_dims + "_cells_" + interp_method + ".nc"

;---Open NetCDF file for writing regridded variables and define some metadata on the file
  regrid_file = open_netcdf_for_regridding(src_file,dfile,interp_method,src_fname)

;---Get variables we want to regrid and define them on NetCDF file
  vars_to_regrid = get_file_vars_to_regrid(src_file,regrid_file)

;---Debug prints
  print("======================================================================")
  print("Source mesh                 : " + src_dims + " cells")
  print("Destination mesh            : " + dst_dims + " cells")
  print("Interpolation method        : " + interp_method)
  print("Weights file name           : " + weights_fname)
  print("Variables being regridded   : " + str_join(vars_to_regrid,","))
  print("======================================================================")

;---Set up regridding options
  Opt                   = True

  Opt@InterpMethod      = interp_method

  Opt@SrcGridLat        = src_lat           ; source grid
  Opt@SrcGridLon        = src_lon
  Opt@SrcGridType       = "unstructured"

  Opt@DstGridLat        = dst_lat           ; destination grid
  Opt@DstGridLon        = dst_lon
  Opt@DstGridType       = "unstructured"

  Opt@WgtFileName    = weights_fname 

  Opt@ForceOverwrite    = True
; Opt@PrintTimings      = True
; Opt@Debug             = True     ; This will spit out a lot of debug information!

; Set these to True only if you already have the weights file.
; Otherwise, the first time in the loop, the weights file will
; get created for you and then used for subsequent variables.
;
; Opt@SkipSrcGrid = True
; Opt@SkipWgtGen  = True
; Opt@SkipDstGrid = True
;

  regrid_file->lat=dst_lat
  regrid_file->lon=dst_lon
  regrid_file->lev=dst_lev
  regrid_file->hyam=hyam
  regrid_file->hybm=hybm
  regrid_file->hyai=hyai
  regrid_file->hybi=hybi
  ;regrid_file->P0=P0
  regrid_file->ilev=ilev
  regrid_file->time=time
  regrid_file->date=date
  regrid_file->datesec=datesec
  regrid_file->time_bnds=time_bnds
  regrid_file->ndbase=ndbase
  regrid_file->nsbase=nsbase
  regrid_file->nbdate=nbdate
  regrid_file->nbsec=nbsec
  regrid_file->mdt=mdt
  regrid_file->ndcur=ndcur
  regrid_file->nscur=nscur
  regrid_file->nsteph=nsteph
  regrid_file->date_written=date_written
  regrid_file->time_written=time_written
  regrid_file->area=dst_area

;---Loop through all the variables to regrid, regrid them, and write to the open NetCDF file.
  nvars = dimsizes(vars_to_regrid)
  do nv=0,nvars-1
    print("======================================================================")
    print("Regridding '" + vars_to_regrid(nv))
    var := src_file->$vars_to_regrid(nv)$
    printVarSummary(var)
    printMinMax(var,0)

    var_regrid := ESMF_regrid(var,Opt)               ; Do the regridding
    regrid_file->$vars_to_regrid(nv)$ = (/var_regrid/)   ; Write only the values to the file since we've
                                                         ; already written attributes and coordinate arrays
    printVarSummary(var_regrid)
    printMinMax(var_regrid,0)
;
; If the weights file has already been generated, then set
; these special Skip resources for much faster regridding
; on the rest of the variables. The regridding will be done 
; using the weights file.
;
    if(nv.eq.0) then
      Opt@SkipSrcGrid = True
      Opt@SkipWgtGen  = True
      Opt@SkipDstGrid = True
    end if
  end do

end
