load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ; Interpolate FV grid to SE grid 
  ;
  ; ASSUMPTIONS: for regridding CAMS 0.1 deg emissions files
  ;    CAMS emissions files have multiple variables: 
  ;         for each emissions sector AND the sum of all sectors
  ;
  ;    *** ESMF problems in NCL:
  ;              (1) The 'current' version of ESMF_regrid_with_weight for 'conservative'
  ;                  interpolation does not work properly for arrays with more dimensions
  ;                  than (lat,lon) --> (ncol)
  ;              (2) The areas (corners) from 'rectilinear_to_SCRIP' are crap for wrap 
  ;                  arround longitudes. 
  ;          ** Explicitly add wrap around points for (fv09) so that the [0.,360.] 
  ;             region is completely covered with properly sized areas. 
  ;          ** Use an older version of the ESMF code (before someone broke it) that will 
  ;             work on arrays of any size. 
  ;          ** Explicilty load Lat,Lon data with wrap arounds and interpolate at each 
  ;             level when mapping to SE. 
  ;             
  ; USAGE: ncl Regrid_fv2se_cams_anthro.ncl
  ;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  USER CHANGES START
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Directory to read grid defintions
  grid_dir ="/glade/work/tilmes/refgrids/"

  ; Directory to write new grid weights (change to your directory)
  wght_dir ="/glade/scratch/$USER/emis/Grid_Weights/"

  ; Type of interpolation
  interp   = "conserve"

  ; Resolution of original file (for weight filename)
  ingrid   = "0.1x0.1"

  ; Resolution of new grid
  outgrid   = "ne30np4"
  ;outgrid   = "ne30pg3"
  ;outgrid  = "conus_30_x8"

  ; Directory of original emissions files - all files will be processed
  srcPath = "/glade/p/acom/acom-climate/tilmes/emis/download/"

  ; Directory for new files
  dstPath = "/glade/scratch/$USER/emis/cams_ne30np4/"

  ; Set output grid file and set coordinate values
  ;   SE needs to use info from pentagons file.
  ;-------------------------------------------------
  if(outgrid.eq."ne30np4") then
    dstGridName = grid_dir+"ne30np4_091226_pentagons.nc"
  end if
  if(outgrid.eq."ne120") then
    dstGridName = grid_dir+"ne120np4_pentagons_100310.nc"
  end if
  if(outgrid.eq."conus_30_x8") then
    dstGridName = grid_dir+"conus_30_x8.g_scrip.nc"
  end if
  if(outgrid.eq."ne30pg3") then
    dstGridName = grid_dir+"ne30pg3.nc"
  end if
  print(" Using SE pentagons file:"+dstGridName)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  USER CHANGES END 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Set some constants for readability
  ;------------------------------------
  COPY    = 0
  DISCARD = 1
  INTERP  = 2
  SLICE   = 3
  NONE    = 4

  ; Get the list of original CAMS 0.1deg filenames to process
  ;------------------------------------------
  syscmd = "(cd "+srcPath+" ; ls -1 CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_*.nc)"
  Flist = systemfunc(syscmd)
  Fnum  = dimsizes(Flist)
  print(" ")
  print(" Number of files to process="+Fnum)
  print(" ")
  print(Flist)

  if(isfilepresent(dstGridName)) then
    dstGrid = addfile(dstGridName,"r")
    dstLon  = dstGrid->grid_center_lon
    dstLon!0          = "ncol"
    dstLon@units      = "degrees_east"
    dstLon@long_name  = "longitude"
    dstLat  = dstGrid->grid_center_lat
    dstLat!0          = "ncol"
    dstLat@units      = "degrees_north"
    dstLat@long_name  = "latgitude"
    dstArea = dstGrid->grid_area
    dstArea!0         = "ncol"
    dstArea@units     = "normalized"
    dstArea@long_name = "gridpoint area"
    dstNcol           = dimsizes(dstArea)
    delete(dstGrid)
  else
    print(" ERROR: outgrid for "+outgrid+" not found")
    exit
  end if

  ; Open the first input file and read in lat/lon values
  ;  to create ESMF grid info.
  ;---------------------------------------------------------
  srcGridName = wght_dir+"SCRIP_"+ingrid+"_cams.nc"

  if(isfilepresent(srcGridName)) then
    print(" Removing old ESMF weights:"+srcGridName)
    system("rm -f "+srcGridName)
  end if
  print(" Creating SCRIP file:"+srcGridName)

 do s=0,Fnum-1
  srcFileName = srcPath+Flist(s)
  if(isfilepresent(srcFileName)) then
    print(" Using FV lat/lon from file:"+srcFileName)
    srcFile = addfile(srcFileName,"r")
    srcLonX = srcFile->lon
    srcLat = srcFile->lat
    delete(srcFile)
    srcNlonX=dimsizes(srcLonX)
    srcNlat =dimsizes(srcLat)
    srcNlon = srcNlonX+3
    srcLon = new(srcNlon,typeof(srcLonX))
    srcLon(1:srcNlonX)=srcLonX(0:(srcNlonX-1))
    srcLon(0)=srcLonX(srcNlonX-1)-360.
    srcLon(srcNlonX+1)=srcLonX(0)+360.
    srcLon(srcNlonX+2)=srcLonX(1)+360.
    ESMFopt               = True
    ESMFopt@ForceOverwrite= True
    ESMFopt@Debug         = False ; True
    ESMFopt@PrintTimings  = False ; True
    ESMFopt@SrcESMF       = False
    ESMFopt@DstESMF       = False
    rectilinear_to_SCRIP(srcGridName,srcLat,srcLon,ESMFopt)
    delete(ESMFopt)
    wrkArray = new( (/ srcNlat, srcNlon /),"double")
  else
    print(" ERROR: Cannot find FV file:"+srcFileName)
    exit
  end if

  ; Remove existing ESMF weights, create new values
  ;---------------------------------------------------------
  wgtGridName = wght_dir+"ESMFmap_"+ingrid+"_"+outgrid+"_"+interp+"_cams.nc"

  if(isfilepresent(wgtGridName)) then
     print(" Removing old ESMF weights:"+wgtGridName)
     system("rm -f "+wgtGridName)
   end if
  print(" Creating ESMF weights:"+wgtGridName)
  print (" ")

  ESMFopt               = True
  ESMFopt@ForceOverwrite= True
  ESMFopt@SrcGridLat    = srcLat
  ESMFopt@SrcGridLon    = srcLon
  ESMFopt@DstGridLat    = dstLat
  ESMFopt@DstGridLon    = dstLon
  ESMFopt@InterpMethod  = interp
  ESMFopt@Debug         = False ; True
  ESMFopt@PrintTimings  = False ; True
  ESMFopt@SrcESMF       = False
  ESMFopt@DstESMF       = False
  ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtGridName,ESMFopt)

  ; Open Source file
  ;------------------
  srcFileName = srcPath+Flist(s)
  if(isfilepresent(srcFileName)) then
    srcFile = addfile(srcFileName,"r")
    print(" ")
    print(" Processing File: "+srcFileName)
    print(" ")
  else
    print(" ")
    print("ERROR: Source Data Not Found?? "+srcFileName)
    print(" ")
    exit
  end if
 
  ; Get attributes from source file, set output attributes
  ;-------------------------------------------------------
  File_in = Flist(s)
  dstAtts               = True
  dstAtts@srcAttEnd     = "---SOURCE-FILE-ATTRIBUTES-------------------"
  copy_VarAtts(srcFile,dstAtts)
  dstAtts@srcAttBegin   = "---SOURCE-FILE-ATTRIBUTES-------------------"
  dstAtts@interp_method =interp 
  dstAtts@source_file   = File_in
  dstAtts@interp        = " FV Source data interpolated to SE grid using ESMF "

  ; Get source dimension names and sizes
  ;-------------------------------------
  srcDname = getvardims     (srcFile)
  srcDsize = getfiledimsizes(srcFile)
   
  srcNlat    = srcDsize(ind(srcDname .eq. "lat"   ))
  srcNlonX   = srcDsize(ind(srcDname .eq. "lon"   ))
  srcNtime   = srcDsize(ind(srcDname .eq. "time"  ))

  ; Set output dimension names and sizes
  ;-------------------------------------
  dstDname = (/ "time" , "ncol"  /)
  dstDsize = (/   -1   , dstNcol  /)
  dstDulim = (/  True  , False    /)

  ; read in coordinate variables
  ;-------------------------------
  dstTime   = srcFile->time

  ; Get a list a variables in this file
  ;------------------------------------
  srcVarNames = getfilevarnames(srcFile)
  srcVarNum   = dimsizes(srcVarNames)

  ; Construct the output filename
  ;------------------------------------
  File_out = outgrid+"_"+Flist(s)

  ; Create the output file
  ;----------------------------
  dstFileName = dstPath+File_out
  if(isfilepresent(dstFileName)) then
    print("   Removing Existing File: "+dstFileName)
    syscmd = "rm "+dstFileName
    Stat=systemfunc(syscmd)
  end if
  print("   Creating: "+dstFileName)
  dstFile = addfile(dstFileName,"c")

  ; write out file attributes and dimensions
  ;------------------------------------------
  setfileoption(dstFile,"DefineMode",True)
  fileattdef(dstFile,dstAtts)
  filedimdef(dstFile,dstDname,dstDsize,dstDulim)

  ; Add Coordinate values
  ;-----------------------
  filevardef(dstFile,"time"  ,typeof(dstTime)  ,getvardims(dstTime)  )
  filevarattdef(dstFile,"time"  ,dstTime  )
      
  ; define SE grid variables
  ;-----------------------------
  filevardef(dstFile,"lon" ,"double",(/"ncol"/))
  filevardef(dstFile,"lat" ,"double",(/"ncol"/))
  filevardef(dstFile,"area","double",(/"ncol"/))

  ; Loop thru and define grid variables
  ;----------------------------------------
  do nv=0,(srcVarNum-1)

    ; skip over FV variables we don't want
    ;---------------------------------------
    if(str_lower(srcVarNames(nv)).eq."latitude"        ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."longitude"        ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."time"       ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."lat_bnds"   ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."lon_bnds"   ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."time_bnds"  ) then continue end if

    ; Get source info for this variable
    ;------------------------------------
    srcVarType  = getfilevartypes(srcFile,srcVarNames(nv))
    srcDimNames = getfilevardims (srcFile,srcVarNames(nv))
    srcNdims    = dimsizes(srcDimNames)

    ; get indices of lat,lon, time dimensions
    ;-----------------------------------------
    ind_lat  = -1
    ind_lon  = -1
    ind_time = -1
    do nd=0,(srcNdims-1)
      if((srcDimNames(nd).eq."lat").or.(srcDimNames(nd).eq."slat") .or.(srcDimNames(nd).eq."latitude")) then 
        ind_lat=nd 
      end if
      if((srcDimNames(nd).eq."lon").or.(srcDimNames(nd).eq."slon") .or.(srcDimNames(nd).eq."longitude")) then 
        ind_lon=nd 
      end if
      if(srcDimNames(nd).eq."time") then 
        ind_time=nd 
      end if
    end do ; nd=0,(srcNdims-1)

    ; set variable processing options
    ;---------------------------------
    VarProc = COPY
    if( (ind_lat.ne. -1).or.(ind_lon.ne. -1)) then
      VarProc = COPY ; DISCARD
      if( (ind_lat.ne. -1).and.(ind_lon.ne. -1)) then
        VarProc = INTERP
      end if
    end if
  
    ; define variables in output file
    ;-----------------------------------
    print(dstFile)
    print(srcVarNames(nv))
    print(srcVarType)
    print(srcDimNames)
    if(VarProc .eq. COPY) then
      filevardef(dstFile,srcVarNames(nv),srcVarType,srcDimNames)
    end if
    if(VarProc .eq. INTERP) then
      dstDimNames = new((srcNdims-1),typeof(srcDimNames))
      md=-1
      do nd=0,(srcNdims-1)
        if(srcDimNames(nd).ne."lat") then 
          md = md +1
          if(srcDimNames(nd).eq."lon") then 
            dstDimNames(md) = "ncol"
          else
            dstDimNames(md) = srcDimNames(nd)
          end if
        end if
      end do ; nd=0,(srcNdims-1)
      filevardef(dstFile,srcVarNames(nv),srcVarType,dstDimNames)
      delete(dstDimNames)
    end if

    ; keep this loop clean
    ;-----------------------
    delete(srcNdims   )
    delete(srcDimNames)
    delete(srcVarType )
  end do ;nv=0,(srcVarNum-1)

  ; end of NETCDF file definitions
  ;--------------------------------
  setfileoption(dstFile,"DefineMode",False)

  ; write out SE grid information 
  ;-------------------------------
  dstFile->time   = (/ dstTime   /)  
  dstFile->lon    =    dstLon  
  dstFile->lat    =    dstLat  
  dstFile->area   =    dstArea 

  ; Loop thru and process grid variables
  ;----------------------------------------
  do nv=0,(srcVarNum-1)

    ; skip over FV variables we don't want
    ;---------------------------------------
    if(str_lower(srcVarNames(nv)).eq."lat"   ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."lon"   ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."time"  ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."lat_bnds"   ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."lon_bnds"   ) then continue end if
    if(str_lower(srcVarNames(nv)).eq."time_bnds"  ) then continue end if

    ; Get source info for this variable
    ;------------------------------------
    srcVarType  = getfilevartypes(srcFile,srcVarNames(nv))
    srcDimNames = getfilevardims (srcFile,srcVarNames(nv))
    srcNdims    = dimsizes(srcDimNames)

    ; get indices of lat,lon, time dimensions
    ;-----------------------------------------
    ind_lat  = -1
    ind_lon  = -1
    ind_time = -1
    do nd=0,(srcNdims-1)
      if((srcDimNames(nd).eq."lat").or.(srcDimNames(nd).eq."slat") .or.(srcDimNames(nd).eq."latitude")) then 
        ind_lat=nd 
      end if
      if((srcDimNames(nd).eq."lon").or.(srcDimNames(nd).eq."slon") .or.(srcDimNames(nd).eq."longitude")) then 
        ind_lon=nd 
      end if
      if(srcDimNames(nd).eq."time") then 
        ind_time=nd 
      end if
    end do ; nd=0,(srcNdims-1)

    ; set variable processing options
    ;---------------------------------
    VarProc = COPY
    if( (ind_lat.ne. -1).or.(ind_lon.ne. -1)) then
      VarProc = COPY ; DISCARD
      if( (ind_lat.ne. -1).and.(ind_lon.ne. -1)) then
        VarProc = INTERP
      end if
    end if

    ; Process variable
    ;--------------------
    Var = srcFile->$srcVarNames(nv)$

    if(VarProc.eq.COPY) then
      dstFile->$srcVarNames(nv)$ = Var
    end if
    if(VarProc.eq.INTERP) then
     print(" INTERP VAR: "+srcVarNames(nv))
     Vinterp = new( (/ srcNtime, dstNcol /),"float")
     Vinterp!0 = "time"
     Vinterp!1 = "ncol"

      ; Loop over the number of times 
      ;------------------------------
      do nt=0,(srcNtime-1  )
        wrkArray(:,:) = wrkArray@_FillValue
        wrkArray(:,1:(srcNlon-3)) = todouble(Var(nt,:,0:(srcNlon-4)))
        wrkArray(:,  0          ) = wrkArray(:,(srcNlon-3))
        wrkArray(:, (srcNlon-2) ) = wrkArray(:,1)
        wrkArray(:, (srcNlon-1) ) = wrkArray(:,2)
        a = tofloat(ESMF_regrid_with_weights(wrkArray,wgtGridName,ESMFopt))
        Vinterp(nt,:) = tofloat(ESMF_regrid_with_weights(wrkArray,wgtGridName,ESMFopt))
      end do

      copy_VarAtts(Var,Vinterp)
      copy_VarCoords(dstFile->$srcVarNames(nv)$,Vinterp)
      dstFile->$srcVarNames(nv)$ = Vinterp
      delete(Vinterp)
    end if
;;  if(VarProc.eq.DISCARD) then
;;    print(" DISCARD VAR: "+srcVarNames(nv))
;;  end if
 
    ; keep this loop clean
    ;-----------------------
    delete(Var)
    delete(srcNdims   )
    delete(srcDimNames)
    delete(srcVarType )
  end do ;nv=0,(srcVarNum-1)

  ; Close Output file
  ;-------------------
  delete(dstFile)
  delete(dstTime)
 end do ; s (number of file)

end

