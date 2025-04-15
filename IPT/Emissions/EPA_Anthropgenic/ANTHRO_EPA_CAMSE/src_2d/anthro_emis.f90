
   program anthro_emis

   use anthro_types
   use utils
   use camse_utils, only : camse_init
   use area_mapper, only : xlong => lon, xlat => lat
   use mapper_types
   use netcdf_utils, only : ngatts
   use misc_definitions_module, only : PROJ_PS, PROJ_LATLON, PROJ_CASSINI, PROJ_CAMSE
   use mo_calendar, only : diffdat, addsec2dat, newdate
   use stack, only : stkFileRead, stkFileCleanup
   use data_file_utils, only : anthro_dir, data_yrs_offset
   use data_file_utils, only : src_lon_dim_name, src_lon_var_name
   use data_file_utils, only : src_lat_dim_name, src_lat_var_name
   use data_file_utils, only : data_file_init, get_src_time_ndx
   use data_file_utils, only : read_src_data, tinterp_src_data
   use EPA, only             : init_table, translate_date

   implicit none

!-----------------------------------------------------------------
!     namelist control variables
!-----------------------------------------------------------------

   integer                :: nemis
   integer                :: emissions_zdim_stag = 10
   integer                :: domains = 1
   integer                :: output_interval = 3600      ! seconds
   character(len=namsize) :: mdlFilenm
   character(len=linsize) :: model_dir
   character(len=linsize) :: src_file_prefix
   character(len=linsize) :: src_file_suffix
   character(len=linsize) :: sec_file_prefix
   character(len=linsize) :: sec_file_suffix
   character(len=linsize) :: stk_file_prefix
   character(len=linsize) :: stk_file_suffix
   character(len=linsize) :: stk_grp_file_prefix
   character(len=linsize) :: stk_grp_file_suffix
   character(len=linsize) :: sectorlist_flnm
   character(len=linsize) :: smk_merge_flnm
   character(len=linsize) :: cat_var_prefix
   character(len=linsize) :: cat_var_suffix
   character(len=32)      :: currDateStr
   character(len=3)       :: numa
   character(len=linsize) :: emis_map(linemax)
   character(len=namsize) :: sub_categories(maxsrc)
   character (LEN=64)     :: mdl_lon_dim_name = ' '
   character (LEN=64)     :: mdl_lon_var_name = ' '
   character (LEN=64)     :: mdl_lat_dim_name = ' '
   character (LEN=64)     :: mdl_lat_var_name = ' '
   character (LEN=19)     :: start_output_time = ' '
   character (LEN=19)     :: stop_output_time  = ' '
   character (LEN=19)     :: start_data_time = ' '
   character (LEN=19)     :: stop_data_time  = ' '
   character (LEN=16)     :: date_frmt  = ' '
   character (LEN=16)     :: mdlType  = 'WRF'
   logical                :: serial_output = .false.

   namelist /control/ anthro_dir, model_dir, emis_map, domains, src_names, &
                      src_file_prefix, src_file_suffix, cat_var_prefix, &
                      cat_var_suffix, emissions_zdim_stag, sub_categories, &
                      date_frmt, start_output_time, stop_output_time, output_interval, &
                      start_data_time, stop_data_time, serial_output, data_yrs_offset, diag_level, &
                      src_lon_dim_name, src_lon_var_name, src_lat_dim_name, src_lat_var_name, &
                      currDateStr, sec_file_prefix, sec_file_suffix, stk_file_prefix, stk_file_suffix,  &
                      stk_grp_file_prefix, stk_grp_file_suffix, &
                      mdlType, mdlFilenm, mdl_lon_dim_name, mdl_lon_var_name, mdl_lat_dim_name, mdl_lat_var_name

!-----------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------
   integer, parameter :: lower = 0
   integer, parameter :: upper = 1
   real, parameter    :: confac_gas = 3.6e12             !(kg/m^2/s -> mole/km^2/hr)
   real, parameter    :: confac_aer = 1.e9               !(kg/m^2/s -> ug/m^2/s)
   real, parameter    :: secphr     = 3600.              ! seconds per hour
   real, parameter    :: dz         = 50.                ! CAM ref height for emissions (meters)

   integer :: nc, ns, ns1, it
   integer :: timeNdx
   integer :: src
   INTEGER :: icnt
   INTEGER :: ids, ide, jds, jde
   INTEGER :: i, j, n, lev
   INTEGER :: ii, jj
   integer :: ncid
   integer :: cat_ndx
   integer :: map_proj
   integer :: ierr, astat, istat
   integer :: dimid, varid
   integer :: nlon_src, nlat_src
   integer :: domain
   integer :: WRFProjNdx, SrcProjNdx
   integer :: currDate, wrkDate
   integer :: xndx_src(2)
   integer :: yndx_src(2)
   integer, allocatable :: nNzSrcWater(:,:)
   integer, allocatable :: nNzSrcLand(:,:)
   integer, allocatable :: ix(:,:,:)                        ! index used by interpolation
   integer, allocatable :: jy(:,:,:)                        ! index used by interpolation

   real    :: scale_factor
   real    :: wrk_sum
   real    :: ds1, ds2
   real    :: xl, xu
   real    :: yl, yu, dy
   real    :: wrf_lon_min
   real    :: wrf_lon_max
   real    :: wrf_lat_min
   real    :: wrf_lat_max
   real    :: cen_lon
   real    :: cen_lat
   real    :: stand_lon
   real    :: truelat1
   real    :: truelat2
   real    :: loninc
   real    :: latinc
   real    :: knowni
   real    :: knownj
   real    :: dx
   real    :: dxRatio
   real    :: dburdenSrf, mburdenSrf
   real, allocatable :: totalSrcLand(:,:)
   real, allocatable :: totalSrcWater(:,:)
   real, allocatable :: wrk_emis(:,:,:)
   real, allocatable :: ax(:,:,:)                        ! weight coef. all domain
   real, allocatable :: by(:,:,:)                        ! weight coef. all domain

   character(len=19)   :: proj_name(0:3) = (/ 'LATLON             ', 'LAMBERT            ', &
                                              'POLAR STEREOGRAPHIC', 'MERCATOR           ' /)
   character(len=linsize) :: wrk_file_prefix
   character(len=linsize) :: wrk_file_suffix
   CHARACTER (LEN=192) :: filespec
   CHARACTER (LEN=132) :: varname
   CHARACTER (LEN=80)  :: message
   CHARACTER (LEN=80)  :: attribute
   CHARACTER (LEN=80)  :: units_attribute
   CHARACTER (LEN=80)  :: description_attribute
   CHARACTER (LEN=80)  :: stagger_attribute
   CHARACTER (LEN=80)  :: coor_attribute
   CHARACTER (LEN=80)  :: memord_attribute
   CHARACTER (LEN=80)  :: outpname
   CHARACTER (LEN=32)  :: wrkDateStr
   CHARACTER (LEN=19)  :: Times(1)
   CHARACTER (LEN=16)  :: wrkStr
   CHARACTER (LEN=3)   :: num

   logical :: Error
   logical :: has_area_map
   logical :: new_grid
   logical :: lexist
   logical :: is_EPA, is_EPA_SECTOR, is_EPA_STACK
   logical :: mdl_is_WRF, mdl_is_CAM, mdl_is_CAMFV, mdl_is_CAMSE
   logical, allocatable :: cat_active(:)
   logical, allocatable :: waterMask(:,:), landMask(:,:)

   type(dates)                :: start_output
   type(dates)                :: stop_output
   type(dates)                :: loop_date
   TYPE(model_grid_type)      :: model_grid
   type(data_file_type), allocatable  :: data_file(:)

!---------------------------------------------------------------------
!  include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

!---------------------------------------------------------------------
!  set namelist variables default values
!---------------------------------------------------------------------
   src_names(:) = ' '
   model_dir    = '.'
   anthro_dir = '.'
   src_file_prefix = ' '
   src_file_prefix = ' '
   sec_file_suffix = ' '
   sec_file_suffix = ' '
   stk_file_prefix = ' '
   stk_file_suffix = ' '
   stk_grp_file_prefix = 'stack_groups_'
   stk_grp_file_suffix = ' '
   currDateStr     = ' '
   cat_var_prefix  = 'emiss_'
   cat_var_suffix  = ' '
   mdlFilenm       = ' '
   sub_categories(:) = ' ' 
   emis_map(:)       = ' ' 
   mdl_lon_dim_name  = 'lon'
   mdl_lon_var_name  = 'lon'
   mdl_lat_dim_name  = 'lat'
   mdl_lat_var_name  = 'lat'
!-----------------------------------------------------------------
!     read control variables
!-----------------------------------------------------------------
   read(*,nml=control,iostat=istat)
   if( istat /= 0 ) then
     write(*,*) 'anthro_emis: failed to read namelist; error = ',istat
     stop
   end if
!-----------------------------------------------------------------
!     check namelist inputs
!-----------------------------------------------------------------
   if( domains < 1 ) then
     write(*,*) 'anthro_emis: domains must be >= 1'
     stop 'Namelist err'
   endif
   call upcase( trim(mdlType),wrkStr )
   if( trim(wrkStr) == 'WRF' .or. trim(wrkStr) == 'CAMFV' .or. &
       trim(wrkStr) == 'CAMSE' ) then
     mdlType = wrkStr
     model_grid%filename = mdlFilenm
   else
     write(*,*) 'anthro_emis: model type ',trim(mdlType)
     write(*,*) '             must be WRF, CAMFV, or CAMSE'
     stop 'Input parm error'
   endif
   mdl_is_WRF   = trim(mdlType) == 'WRF'
   mdl_is_CAMFV = trim(mdlType) == 'CAMFV'
   mdl_is_CAMSE = trim(mdlType) == 'CAMSE'
   mdl_is_CAM   = mdl_is_CAMFV .or. mdl_is_CAMSE
   model_grid%mdl_is_WRF = mdl_is_WRF
   model_grid%mdl_is_CAM = mdl_is_CAM
   model_grid%mdl_is_CAMFV = mdl_is_CAMFV
   model_grid%mdl_is_CAMSE = mdl_is_CAMSE
   do n = 1,domains
     if( mdl_is_WRF) then
       filespec = trim(model_dir) // '/wrfinput_d'
       write(filespec(len_trim(filespec)+1:),'(i2.2)') n
     else
       filespec = trim(model_dir) // '/' // trim(mdlFilenm)
     endif
     inquire( file=trim(filespec),exist=lexist )
     if( .not. lexist ) then
       write(*,*) 'anthro_emis: model input file'
       write(*,*) trim(filespec)
       write(*,*) 'anthro_emis: does not exist'
       stop 'File err'
     endif
   end do
   if( mdl_is_CAM ) then
     emissions_zdim_stag = 8
     serial_output       = .true.
   elseif( mdl_is_WRF .and. emissions_zdim_stag < 1 ) then
     write(*,*) 'anthro_emis: emissions zdim must >= 1'
     stop 'Namelist err'
   endif
   if( serial_output .and. output_interval < 0 ) then
     write(*,*) 'anthro_emis: output_interval must be >= 0'
     stop 'Input parm error'
   endif
   if( currDateStr == ' ' ) then
     write(*,*) 'anthro_emis: current date string(currDateStr) not set'
     stop 'Input parm error'
   endif

!-----------------------------------------------------------------
!     set species maps
!-----------------------------------------------------------------
   call mapper( nemis, emis_map, sub_categories, data_file )

   if( any( data_file(:)%is_EPA_STACK ) ) then
     if( trim(stk_grp_file_suffix) == ' ' ) then
       stk_grp_file_suffix = trim(stk_file_suffix)
     endif
   endif
!-----------------------------------------------------------------
!     initialize EPA date translation
!-----------------------------------------------------------------
   if( any( data_file(:)%is_EPA ) ) then
     smk_merge_flnm  = trim(anthro_dir) // '/' // 'smk_merge_dates_' // currDateStr(1:4) // '.txt'
     sectorlist_flnm = trim(anthro_dir) // '/' // 'sectorlist_' // currDateStr(1:4) // 'fd_nata'
     Error = init_table(  sectorlist_flnm, smk_merge_flnm )
     if( Error ) then
       write(*,*) 'anthro_emis: EPA date error in init_table'
       stop
     endif
   endif

!-----------------------------------------------------------------
!     check that initial input files exist
!-----------------------------------------------------------------
src_loop_a: &
   do n = 1,nsrc_names
     if( .not. data_file(n)%is_EPA ) then
       wrkDateStr = currDateStr
       filespec = trim(anthro_dir) // '/' // trim(src_file_prefix) &
                                   // trim(src_names(n)) // '_' // trim(wrkDateStr) &
                                   // trim(src_file_suffix)
     else
       read(currDateStr,*) currDate
       wrkDate = translate_date( currDate, trim(src_names(n)) )
       write(wrkDateStr,'(i8)') wrkDate
       if( data_file(n)%is_EPA_SECTOR ) then
         filespec = trim(anthro_dir) // '/' // trim(src_names(n)) // '/' // trim(sec_file_prefix) &
                                     // trim(src_names(n)) // '_' // trim(wrkDateStr) &
                                     // trim(sec_file_suffix)
       else
         filespec = trim(anthro_dir) // '/' // trim(src_names(n)) // '/'// trim(stk_file_prefix) &
                                     // trim(src_names(n)) // '_' // trim(wrkDateStr) &
                                     // trim(stk_file_suffix)
       endif
     endif
     
     inquire( file=trim(filespec),exist=lexist )
     if( .not. lexist ) then
       write(*,*) 'anthro_emis: anthro source file'
       write(*,*) trim(filespec)
       write(*,*) 'anthro_emis: does not exist'
       stop 'File err'
     endif

     if( data_file(n)%is_EPA_STACK ) then
       filespec = trim(anthro_dir) // '/' // trim(src_names(n)) // '/' // trim(stk_grp_file_prefix) &
                                   // trim(src_names(n)) // trim(stk_grp_file_suffix)
       inquire( file=trim(filespec),exist=lexist )
       if( .not. lexist ) then
         write(*,*) 'anthro_emis: anthro stack group file'
         write(*,*) trim(filespec)
         write(*,*) 'anthro_emis: does not exist'
         stop 'File err'
       endif
     endif
     data_file(n)%src_name    = src_names(n)
     data_file(n)%flnmDateStr = wrkDateStr
     data_file(n)%flnmDate    = wrkDate
   end do src_loop_a

   write(*,*) 'main: nemis = ',nemis
   do n = 1,nemis
     write(*,*) '================================================='
     write(*,*) 'anthro_map(',n,'):'
     write(*,*) 'src species count = ',anthro_map(n)%src_cnt
     write(*,*) 'emis species name = ',anthro_map(n)%emis_name
     write(*,*) anthro_map(n)%src_wght(:anthro_map(n)%src_cnt)
     write(*,*) anthro_map(n)%src_var(:anthro_map(n)%src_cnt)
     write(*,*) 'cat wghts'
     do ns = 1,anthro_map(n)%src_cnt
       write(*,*) anthro_map(n)%cat_wght(:,ns)
     end do
   end do
   write(*,*) ' '
   write(*,*) 'src active'
   write(*,*) src_active(:nsrc_names)
   write(*,*) 'active src = ',count(src_active(:))
   write(*,*) '================================================='

   allocate( cat_active(n_sub_cats),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'anthro_emis: failed to allocate cat_active array'
     stop 'Alloc err'
   endif

   if( any( data_file(:)%is_EPA_SECTOR .and. mdl_is_WRF ) ) then
       if( allocated( nNzSrcWater ) ) then
        deallocate( nNzSrcWater )
       endif
       allocate( nNzSrcWater(nsrc_names,n_sub_cats),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'anthro_emis: failed to allocate nNzSrcWater array'
         stop 'Alloc err'
       endif
       if( allocated( nNzSrcLand ) ) then
        deallocate( nNzSrcLand )
       endif
       allocate( nNzSrcLand(nsrc_names,n_sub_cats),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'anthro_emis: failed to allocate nNzSrcLand array'
         stop 'Alloc err'
       endif
       if( allocated( totalSrcLand ) ) then
        deallocate( totalSrcLand )
       endif
       allocate( totalSrcLand(nsrc_names,n_sub_cats),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'anthro_emis: failed to allocate totalSrcLand array'
         stop 'Alloc err'
       endif
       if( allocated( totalSrcWater ) ) then
        deallocate( totalSrcWater )
       endif
       allocate( totalSrcWater(nsrc_names,n_sub_cats),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'anthro_emis: failed to allocate totalSrcWater array'
         stop 'Alloc err'
       endif
   endif

domain_loop : &
   do domain = 1,domains
     timeNdx = 0
     call model_init
!-----------------------------------------------------------------
!     setup data, output times
!-----------------------------------------------------------------
first_domain : &
     if( domain == 1 ) then
       if( start_output_time == ' ' ) then
         start_output_time = Times(1)
       endif
       call wrf2mz_time( start_output_time, start_output%date, start_output%secs )
       if( .not. serial_output ) then
         stop_output_time = start_output_time
       elseif( stop_output_time == ' ' ) then
         stop_output_time = start_output_time
       endif
       call wrf2mz_time( stop_output_time, stop_output%date, stop_output%secs )
       if( diffdat( start_output%date, start_output%secs, stop_output%date, stop_output%secs ) < 0. ) then
         write(*,*) 'anthro_emis: start output time > stop output time'
         stop 'Input parameter error'
       endif
       do src = 1,nsrc_names
         cat_active(:) = .false.
         if( src_active(src) ) then
           do ns = 1,nemis
             do ns1 = 1,anthro_map(ns)%src_cnt
               if( trim(anthro_map(ns)%src_var(ns1)) ==  trim(src_names(src)) ) then
                 cat_active(:) = cat_active(:) .or. anthro_map(ns)%cat_wght(:,ns1) /= 0.
               endif
             end do
           end do
           allocate( data_file(src)%cat_active(n_sub_cats),stat=astat )
           if( astat /= 0 ) then
             write(*,*) 'anthro_emis: allocate for data_file%cat_active failed; error = ', astat
            stop 'Alloc err'
           endif
           data_file(src)%cat_active(:) = cat_active(:)
           write(*,*) 'will use sub cats: ',cat_active(:)
         endif
       end do
     endif first_domain

     if( .not. allocated( wrk_emis ) ) then
       allocate( wrk_emis(ide,jde,emissions_zdim_stag),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'anthro_emis: allocate for wrk_emis failed; error = ', astat
        stop 'Alloc err'
       endif
     endif
     if( any( data_file(:)%is_EPA_SECTOR ) .and. mdl_is_WRF ) then
       if( .not. allocated( waterMask ) ) then
         allocate( waterMask(ide,jde),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'anthro_emis: failed to allocate totalSrcWater array'
           stop 'Alloc err'
         endif
       endif
       if( .not. allocated( landMask ) ) then
         allocate( landMask(ide,jde),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'anthro_emis: failed to allocate totalSrcWater array'
           stop 'Alloc err'
         endif
       endif
     endif
!-----------------------------------------------------------------
!     initialize data file type
!-----------------------------------------------------------------
src_loop_b: &
     do src = 1,nsrc_names
       if( src_active(src) ) then
         if( domain == 1 ) then
           data_file(src)%molecw   = molecw(src)
         endif
         if( .not. data_file(src)%is_EPA ) then
           data_file(src)%filename = trim(src_file_prefix) // trim(src_names(src)) // '_' &
                                     // trim(data_file(src)%flnmDateStr) // trim(src_file_suffix)
         elseif( data_file(src)%is_EPA_SECTOR ) then
           data_file(src)%filename = trim(sec_file_prefix) // trim(src_names(src)) // '_' &
                                     // trim(data_file(src)%flnmDateStr) // trim(sec_file_suffix)
         elseif( data_file(src)%is_EPA_STACK ) then
           data_file(src)%filename = trim(stk_file_prefix) // trim(src_names(src)) // '_' &
                                     // trim(data_file(src)%flnmDateStr) // trim(stk_file_suffix)
         endif
         data_file(src)%filespec = trim(anthro_dir) // '/' // trim(src_names(src)) // '/' // trim(data_file(src)%filename)
         if( data_file(src)%is_EPA_STACK ) then
           data_file(src)%stack%filename = trim(stk_grp_file_prefix) // trim(src_names(src)) // trim(stk_grp_file_suffix)
           data_file(src)%stack%filespec = &
             trim(anthro_dir) // '/' // trim(src_names(src)) // '/' // trim(data_file(src)%stack%filename)
         endif
         call data_file_init( data_file(src), start_output, cat_var_prefix, cat_var_suffix, &
                              domain, dx, ids, ide, jds, jde, &
                              currDate, model_grid, emissions_zdim_stag )
       endif
     end do src_loop_b
!-----------------------------------------------------------------
!     allocate emission array
!-----------------------------------------------------------------
     do ns = 1,nemis
       anthro_map(ns)%is_2d(1:nsrc_names) = .true.
       anthro_map(ns)%is_3d(1:nsrc_names) = .false.
       if( mdl_is_WRF .or. mdl_is_CAM ) then
         do src = 1,anthro_map(ns)%src_cnt
           ns1 = get_src_ndx( trim(anthro_map(ns)%src_var(src)) )
           if( src_active(ns1) ) then
             anthro_map(ns)%is_3d(ns1) = data_file(ns1)%is_EPA_STACK
             anthro_map(ns)%is_2d(ns1) = .not. anthro_map(ns)%is_3d(ns1)
           endif
         enddo
       endif
       if( any( anthro_map(ns)%is_2d(1:anthro_map(ns)%src_cnt) ) ) then
         if( allocated( anthro_map(ns)%emissionSrf ) ) then
           deallocate( anthro_map(ns)%emissionSrf )
         endif
         allocate( anthro_map(ns)%emissionSrf(ide,jde,emissions_zdim_stag),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'anthro_emis: failed to allocate emissionSrf array'
           stop 'Alloc err'
         endif
       endif
       if( any( anthro_map(ns)%is_3d(1:anthro_map(ns)%src_cnt) ) ) then
         if( allocated( anthro_map(ns)%emissionLay ) ) then
           deallocate( anthro_map(ns)%emissionLay )
         endif
         allocate( anthro_map(ns)%emissionLay(ide,jde,emissions_zdim_stag),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'anthro_emis: failed to allocate emissionLay array'
           stop 'Alloc err'
         endif
       endif
     end do

     loop_date = start_output
time_loop : &
     do
       do ns = 1,nemis
         if( allocated(anthro_map(ns)%emissionSrf) ) then
           anthro_map(ns)%emissionSrf(:,:,:) = 0.
         endif
         if( allocated(anthro_map(ns)%emissionLay) ) then
           anthro_map(ns)%emissionLay(:,:,:) = 0.
         endif
       end do
!      if( data_file(src)%is_EPA_SECTOR ) then
       if( any( data_file(:)%is_EPA_SECTOR .and. mdl_is_WRF ) ) then
         nNzSrcWater(:,:) = 0
         nNzSrcLand(:,:)  = 0
         totalSrcWater(:,:) = 0.
         totalSrcLand(:,:)  = 0.
       endif
       dburdenSrf = 0.
       mburdenSrf = 0.
src_loop_c : &
       do src = 1,nsrc_names
use_src: if( src_active(src) ) then
           write(*,*) ' '
           write(*,*) 'will use source file for ',trim(src_names(src))
           if( .not. data_file(src)%is_EPA ) then
             wrk_file_prefix = src_file_prefix ; wrk_file_suffix = src_file_suffix
           elseif( data_file(src)%is_EPA_SECTOR ) then
             wrk_file_prefix = sec_file_prefix ; wrk_file_suffix = sec_file_suffix
           elseif( data_file(src)%is_EPA_STACK ) then
             wrk_file_prefix = stk_file_prefix ; wrk_file_suffix = stk_file_suffix
           endif
           call get_src_time_ndx( data_file(src), loop_date, currDate, src_names(src), &
                                  wrk_file_prefix, wrk_file_suffix )

           cat_ndx = 0
cat_loop: do nc = 1,n_sub_cats
use_cat :    if( data_file(src)%cat_active(nc) ) then
               varname = trim(cat_var_prefix) // trim(sub_cats(nc)) // trim(cat_var_suffix)
               cat_ndx = cat_ndx + 1
               if( data_file(src)%read_lo_tndx ) then
                 call read_src_data( data_file(src), varname, data_file(src)%lo_buf_ndx, &
                                     data_file(src)%lo_tndx, data_file(src)%ncid_lo, cat_ndx )
               endif
               if( data_file(src)%t_interp ) then
                 if( data_file(src)%read_hi_tndx ) then
                   call read_src_data( data_file(src), varname, data_file(src)%hi_buf_ndx, &
                                       data_file(src)%hi_tndx, data_file(src)%ncid_hi, cat_ndx )
                 endif
               endif
               call tinterp_src_data( data_file(src), cat_ndx )
               call map_src_emissions( data_file(src), grid_specs(data_file(src)%grid_ndx) )
               do ns = 1,nemis
                 do ns1 = 1,anthro_map(ns)%src_cnt
                   if( trim(anthro_map(ns)%src_var(ns1)) ==  trim(src_names(src)) ) then
                     if( anthro_map(ns)%cat_wght(nc,ns1) /= 0.) then
                       if( anthro_map(ns)%is_gas ) then
                         scale_factor = anthro_map(ns)%src_wght(ns1) &
                                        * anthro_map(ns)%cat_wght(nc,ns1) * data_file(src)%con_fac(1) &
                                        / data_file(src)%molecw
                       else
                         scale_factor = anthro_map(ns)%src_wght(ns1) &
                                        * anthro_map(ns)%cat_wght(nc,ns1) * data_file(src)%con_fac(2)
                       endif
!                      if( .not. data_file(src)%is_EPA_STACK ) then
                       if( anthro_map(ns)%is_2d(ns1) ) then
                         anthro_map(ns)%emissionSrf(:,:,1) = anthro_map(ns)%emissionSrf(:,:,1) &
                                                           + scale_factor*wrk_emis(:,:,1)
                         if( nc == 1 .and. ns == 1 ) then
                           mburdenSrf = mburdenSrf + sum( wrk_emis(:,:,1) ) 
                           dburdenSrf = dburdenSrf + sum( data_file(src)%src_data(:,:) )
                           write(*,*) ' '
                           write(*,*) 'anthro_emis: datburdenSrf = ',dburdenSrf*3.1536e7*28.e-12,' Tg/yr'
                           write(*,*) 'anthro_emis: mdlburdenSrf = ',mburdenSrf*3.1536e7*28.e-12,' Tg/yr'
                         endif
                       elseif( anthro_map(ns)%is_3d(ns1) ) then
                         do lev = 1,emissions_zdim_stag
                           anthro_map(ns)%emissionLay(:,:,lev) = anthro_map(ns)%emissionLay(:,:,lev) &
                                                               + scale_factor*wrk_emis(:,:,lev)
                         enddo
                       endif
!                      else
!                        do lev = 1,emissions_zdim_stag
!                          anthro_map(ns)%emissionSrf(:,:,lev) = anthro_map(ns)%emissionSrf(:,:,lev) &
!                                                              + scale_factor*wrk_emis(:,:,lev)
!                        enddo
!                      endif
                       if( data_file(src)%is_EPA_SECTOR .and. mdl_is_WRF ) then
                         dxRatio = model_grid%dx/data_file(src)%dx
                         dxRatio = dxRatio*dxRatio
                         nNzSrcWater(src,nc) = count( wrk_emis(:,:,1) > 0. .and.  model_grid%landmask(:,:) == 0. )
                         nNzSrcLand(src,nc)  = count( wrk_emis(:,:,1) > 0. .and.  model_grid%landmask(:,:) == 1. )
                         totalSrcWater(src,nc) = secphr*dxRatio*sum( wrk_emis(:,:,1),mask=model_grid%landmask(:,:) == 0. )
                         totalSrcLand(src,nc)  = secphr*dxRatio*sum( wrk_emis(:,:,1),mask=model_grid%landmask(:,:) == 1. )
                         waterMask(:,:) = wrk_emis(:,:,1) > 0. .and.  model_grid%landmask(:,:) == 0.
                         landMask(:,:)  = wrk_emis(:,:,1) > 0. .and.  model_grid%landmask(:,:) == 1.
                         write(*,*) 'anthro_emis: there are ',count(wrk_emis(:,:,1) > 0.),' non-zero emissions'
                       endif
                     endif
                   endif
                 end do
               end do
             endif use_cat
           end do cat_loop

           if( diag_level > 300 ) then
             do ns = 1,nemis
             do ns1 = 1,anthro_map(ns)%src_cnt
               if( trim(anthro_map(ns)%src_var(ns1)) ==  trim(src_names(src)) ) then
                 do nc = 1,n_sub_cats
                   if( cat_active(nc) .and. anthro_map(ns)%cat_wght(nc,ns1) /= 0.) then
                     write(*,*) 'will use sub cats ',trim(sub_cats(nc)),' with wght = ',anthro_map(ns)%cat_wght(nc,ns1)
                   endif
                 end do
               endif
             end do
             end do
           endif
         endif use_src
       end do src_loop_c
!---------------------------------------------------------------------
!   write emission file
!---------------------------------------------------------------------
       if( serial_output ) then
         if( mdl_is_WRF ) then
           call write_wrf_emis( 0 )
         elseif( mdl_is_CAM ) then
           call write_cam_emis( diffdat( start_output%date, start_output%secs, &
                                         loop_date%date, loop_date%secs ) == 0. )
         endif
         call addsec2dat( output_interval, loop_date%date, loop_date%secs )
         if( diffdat( loop_date%date, loop_date%secs, stop_output%date, stop_output%secs ) < 0. ) then
           exit time_loop
         elseif( loop_date%date > currDate .and. any( data_file(:)%is_EPA ) ) then
           currDate = newdate( currDate, 1 )
         endif
       else
         call write_wrf_emis( 1 )
         call write_wrf_emis( 2 )
         exit time_loop
       endif
     end do time_loop
     
!---------------------------------------------------------------------
!   cleanup for next domain
!---------------------------------------------------------------------
     do src = 1,nsrc_names
       if( data_file(src)%ncid_lo /= 0 ) then
         ierr = nf_close( data_file(src)%ncid_lo )
       endif
       if( data_file(src)%ncid_hi /= 0 .and. &
           data_file(src)%ncid_hi /= data_file(src)%ncid_lo ) then
         ierr = nf_close( data_file(src)%ncid_hi )
       endif
       if( allocated( data_file(src)%emis ) ) then
         deallocate( data_file(src)%emis )
       endif
       if( allocated( data_file(src)%src_data ) ) then
         deallocate( data_file(src)%src_data )
       endif
       if( data_file(src)%is_EPA_STACK ) then
         call stkFileCleanup( data_file(src)%stack )
       endif
       if( data_file(src)%grid_ndx > 0 ) then
         call cleanup_grid( grid_specs(data_file(src)%grid_ndx), ide, jde )
       endif
     end do
     grid_cnt = 0
     do ns = 1,nemis
       if( allocated(anthro_map(ns)%emissionSrf) ) then
         deallocate( anthro_map(ns)%emissionSrf )
       endif
       if( allocated(anthro_map(ns)%emissionLay) ) then
         deallocate( anthro_map(ns)%emissionLay )
       endif
     end do
     deallocate( wrk_emis )
     if( mdl_is_WRF .and. domains > 1 ) then
       read(currDateStr,*) currDate
     endif
   end do domain_loop

   do ns = 1,nemis
     if( allocated( anthro_map(ns)%cat_wght ) ) then
       deallocate( anthro_map(ns)%cat_wght )
     endif
   end do
   if( allocated( anthro_map ) ) then
     deallocate( anthro_map )
   endif
   if( allocated( data_file ) ) then
     deallocate( data_file )
   endif

   write(*,*) ' '
   write(*,*) '----------------------------------'
   write(*,*) 'anthro_emis completed successfully'
   write(*,*) '----------------------------------'

   CONTAINS

   subroutine model_init
!---------------------------------------------------------------------
!   read model file
!---------------------------------------------------------------------

   use area_mapper, only : proj_init
   use netcdf_utils, only : get_glb_atts, handle_ncerr
   use constants_module, only : rad_per_deg, earth_radius_m, g

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer :: j, k
   integer :: astat 
   integer :: maxlonndx(1)
   integer, allocatable :: isltyp(:,:) 
   real, allocatable :: ph(:,:,:)
   real, allocatable :: phb(:,:,:)
   real, allocatable :: wrk(:,:)
   CHARACTER (LEN=132) :: inpname
   character(len=64) :: mdl_vrt_dim_name, mdl_vrt_var_name

   if( mdl_is_WRF ) then
     inpname = 'wrfinput_d'
     write(inpname(len_trim(inpname)+1:),'(i2.2)') domain
     mdl_lon_dim_name = 'west_east'
     mdl_lat_dim_name = 'south_north'
     mdl_vrt_dim_name = 'bottom_top'
   else
     inpname = model_grid%filename
   endif
   filespec = trim( model_dir ) // '/' // trim( inpname )
!---------------------------------------------------------------------
!   open model input file
!---------------------------------------------------------------------
   message = 'model_init: Failed to open ' // trim(inpname)
   call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
   if( mdl_is_CAMSE ) then
     call camse_init( ncid, model_grid )
     ide = model_grid%nPolygons ; jde = 1
!---------------------------------------------------------------------
!   get model dimesions
!---------------------------------------------------------------------
   else
     message = 'Failed to get longitude dimension id'
     call handle_ncerr( nf_inq_dimid( ncid, trim(mdl_lon_dim_name), dimid ), message )
     message = 'Failed to get longitude dimension'
     call handle_ncerr( nf_inq_dimlen( ncid, dimid, ide ), message )
     message = 'Failed to get latitude dimension id'
     call handle_ncerr( nf_inq_dimid( ncid, trim(mdl_lat_dim_name), dimid ), message )
     message = 'Failed to get latitude dimension'
     call handle_ncerr( nf_inq_dimlen( ncid, dimid, jde ), message )
   endif

   model_grid%nlons = ide
   model_grid%nlats = jde

   if( mdl_is_WRF ) then
     message = 'Failed to get vertical dimension id'
     call handle_ncerr( nf_inq_dimid( ncid, trim(mdl_vrt_dim_name), dimid ), message )
     message = 'Failed to get vertical dimension'
     call handle_ncerr( nf_inq_dimlen( ncid, dimid, model_grid%nlevs ), message )
!---------------------------------------------------------------------
!   get model map projection variables
!---------------------------------------------------------------------
     message = 'Failed to get map_proj'
     call handle_ncerr( nf_get_att_int( ncid, nf_global, 'MAP_PROJ', map_proj ), message )
   elseif( mdl_is_CAM ) then
!---------------------------------------------------------------------
!   CAM "reference" vertical levels are 50 m deep
!---------------------------------------------------------------------
     model_grid%nlevs = 8
     allocate( model_grid%z_at_w(1,model_grid%nlevs+1,1),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'model_init: failed to allocate z_at_w; error = ',astat
       stop 'Alloc error'
     endif
     model_grid%z_at_w(1,1:model_grid%nlevs+1,1) = (/ (dz*real(k),k=0,model_grid%nlevs) /)
     if( mdl_is_CAMFV ) then
       map_proj = PROJ_LATLON
       allocate( model_grid%mdl_cell_area(jde,1),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'model_init: failed to allocate mdl_cell_area; error = ',astat
         stop 'Alloc error'
       endif
     elseif( mdl_is_CAMSE ) then
       map_proj = PROJ_CAMSE
     endif
   endif
   if( .not. mdl_is_CAMSE ) then
     write(*,*) ' '
     write(*,*) 'model_init: MAP_PROJ is ',trim(proj_name(map_proj))
     write(*,*) ' '
   endif
proj_is_latlon: &
   if( map_proj == PROJ_LATLON .or. map_proj == PROJ_CASSINI ) then
     message = 'Failed to get ' // trim(mdl_lon_var_name) // ' variable id'
     call handle_ncerr( nf_inq_varid( ncid, trim(mdl_lon_var_name), varid ), message )
     if( mdl_is_WRF ) then
       allocate( wrk(ide,jde),model_grid%mdl_lons(ide,jde),stat=astat )
     elseif( mdl_is_CAMFV ) then
       allocate( wrk(ide,1),model_grid%mdl_lons(ide,1),stat=astat )
     endif
     if( astat /= 0 ) then
       write(*,*) 'model_init: failed to allocate model_grid mdl_lons variable; error = ',astat
       stop 'Alloc error'
     endif
     message = 'Failed to read ' // trim(mdl_lon_var_name) // ' variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     model_grid%mdl_lons(:,:) = wrk(:,:)
     loninc = abs(wrk(2,1) - wrk(1,1))
!---------------------------------------------------------------------
!   modify incoming lon grid for CAM model
!---------------------------------------------------------------------
     if( mdl_is_CAMFV ) then
       where( wrk(:,1) > 180. )
         wrk(:,1) = wrk(:,1) - 360.
       endwhere
       maxlonndx(:) = maxloc( wrk(:,1) )
       wrk(1,1) = wrk(maxlonndx(1)+1,1)
       model_grid%lon_shift = maxlonndx(1)
     endif
     cen_lon = wrk(1,1)
     if( mdl_is_CAMFV ) then
       deallocate( wrk )
       allocate( wrk(1,jde),model_grid%mdl_lats(1,jde),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'model_init: failed to allocate model_grid mdl_lats variable; error = ',astat
         stop 'Alloc error'
       endif
     endif
     message = 'Failed to get ' // trim(mdl_lat_dim_name) // ' variable id'
     call handle_ncerr( nf_inq_varid( ncid, trim(mdl_lat_dim_name), varid ), message )
     message = 'Failed to read ' // trim(mdl_lat_var_name) // ' variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     model_grid%mdl_lats(:,:) = wrk(:,:)
     cen_lat = wrk(1,1)
     latinc = abs(wrk(1,2) - cen_lat)
     knowni = 1.
     knownj = 1.
     deallocate( wrk )
     map_proj = PROJ_LATLON
!---------------------------------------------------------------------
!   setup model cell areas for CAMFV
!---------------------------------------------------------------------
     dx = 1.e4 * earth_radius_m * earth_radius_m * loninc * rad_per_deg
     if( mdl_is_CAMFV ) then
       model_grid%mdl_cell_area(1,1) = 1.e20 ; model_grid%mdl_cell_area(jde,1) = 1.e20
       do j = 2,model_grid%nlats-1
         model_grid%mdl_cell_area(j,1) = &
             dx * (sin((model_grid%mdl_lats(1,j) + latinc)*rad_per_deg) &
                       - sin((model_grid%mdl_lats(1,j) - latinc)*rad_per_deg))
       enddo
     dx = earth_radius_m * latinc * rad_per_deg
     endif
   elseif( .not. mdl_is_CAMSE ) then  proj_is_latlon
     message = 'model_init: Failed to get cen_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LON', cen_lon ), message )
     write(*,*) 'model_init: CEN_LON = ',cen_lon
     message = 'model_init: Failed to get cen_lat'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LAT', cen_lat ), message )
     write(*,*) 'model_init: CEN_LAT = ',cen_lat
     message = 'model_init: Failed to get stand_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'STAND_LON', stand_lon ), message )
     write(*,*) 'model_init: STAND_LON = ',stand_lon
     message = 'model_init: Failed to get truelat1'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT1', truelat1 ), message )
     write(*,*) 'model_init: TRUELAT1 = ',truelat1
     message = 'model_init: Failed to get truelat2'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT2', truelat2 ), message )
     write(*,*) 'model_init: TRUELAT2 = ',truelat2
     message = 'model_init: Failed to get dx'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'DX', dx ), message )
     write(*,*) 'model_init: DX = ',dx
   elseif( mdl_is_CAMSE ) then  proj_is_latlon
   endif proj_is_latlon
   if( .not. mdl_is_CAMSE ) then
     model_grid%dx = dx
!---------------------------------------------------------------------
!   initialize map projection for model
!---------------------------------------------------------------------
     call proj_init( map_proj, cen_lon, cen_lat, truelat1, truelat2, &
                     stand_lon, loninc, latinc, knowni, knownj, &
                     dx, ide, jde )
   endif

   ids = 1
   jds = 1

is_mdl_WRF: &
   if( mdl_is_WRF ) then
     message = 'Failed to get Times id'
     call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), message )
     message = 'Failed to read Times'
     call handle_ncerr( nf_get_var_text( ncid, varid, Times ), message )

     write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     write(*,*) 'model_init: time = ',trim(Times(1))
     write(*,*) 'model_init: grid dimensions'
     write(*,*) 'model_init: ids,ide,jds,jde'
     write(*,'(4i6)') ids,ide,jds,jde
     write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!---------------------------------------------------------------------
!   get wrfinput_<domain> global attributes
!---------------------------------------------------------------------
     call get_glb_atts( ncid )
!---------------------------------------------------------------------
!   if any data file is epa then allocate and set landmask
!---------------------------------------------------------------------
is_epa_type: &
    if( any( data_file(:)%is_EPA ) ) then
      if( allocated( model_grid%landmask ) ) then
        deallocate( model_grid%landmask )
      endif
      allocate( model_grid%landmask(ide,jde),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'model_init: failed to allocate landmask; error = ',astat
        stop
      endif
      if( allocated( isltyp) ) then
        deallocate( isltyp )
      endif
      allocate( isltyp(ide,jde),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'model_init: failed to allocate isltyp; error = ',astat
        stop
      endif
      message = 'Failed to get ISLTYP var id'
      call handle_ncerr( nf_inq_varid( ncid, 'ISLTYP', varid ), message )
      message = 'Failed to read ISLTYP variable'
      call handle_ncerr( nf_get_var_int( ncid, varid, isltyp ), message )
      do j = 1,jde
        where( isltyp(:,j) == 14 )
          model_grid%landmask(:,j) = 0.
        elsewhere
          model_grid%landmask(:,j) = 1.
        endwhere
      enddo
      deallocate( isltyp )
      
      write(*,*) ' '
      write(*,*) 'model_init: there are ',count(model_grid%landmask(:,:) == 0.),' water cells'
      write(*,*) ' '

!---------------------------------------------------------------------
!   if any data file is epa-stack then allocate and set vertical grid
!---------------------------------------------------------------------
is_stack: &
      if( any( data_file(:)%is_EPA_STACK ) ) then
        if( allocated( model_grid%z_at_w ) ) then
          deallocate( model_grid%z_at_w )
        endif
        allocate( model_grid%z_at_w(ide,model_grid%nlevs+1,jde),stat=astat )
        if( astat /= 0 ) then
          write(*,*) 'model_init: failed to allocate z_at_w; error = ',astat
          stop
        endif
        if( allocated( ph ) ) then
          deallocate( ph )
        endif
        allocate( ph(ide,jde,model_grid%nlevs+1),stat=astat )
        if( astat /= 0 ) then
          write(*,*) 'model_init: failed to allocate ph; error = ',astat
          stop
        endif
        message = 'Failed to get PH var id'
        call handle_ncerr( nf_inq_varid( ncid, 'PH', varid ), message )
        message = 'Failed to read PH variable'
        call handle_ncerr( nf_get_var_real( ncid, varid, ph ), message )
        if( allocated( phb ) ) then
          deallocate( phb )
        endif
        allocate( phb(ide,jde,model_grid%nlevs+1),stat=astat )
        if( astat /= 0 ) then
          write(*,*) 'model_init: failed to allocate phb; error = ',astat
          stop
        endif
        message = 'Failed to get PHB var id'
        call handle_ncerr( nf_inq_varid( ncid, 'PHB', varid ), message )
        message = 'Failed to read PHB variable'
        call handle_ncerr( nf_get_var_real( ncid, varid, phb ), message )
        do k = 1,model_grid%nlevs+1
          do j = jds,jde
            model_grid%z_at_w(ids:ide,k,j) = (ph(ids:ide,j,k) + phb(ids:ide,j,k))/g
          enddo
        enddo
        do k = 2,model_grid%nlevs+1
          do j = jds,jde
            model_grid%z_at_w(ids:ide,k,j) = model_grid%z_at_w(ids:ide,k,j) &
                                           - model_grid%z_at_w(ids:ide,1,j)
          enddo
        enddo
        do j = jds,jde
          model_grid%z_at_w(ids:ide,1,j) = 0.
        enddo
        deallocate( ph, phb )
      endif is_stack
    endif is_epa_type
    endif is_mdl_WRF
!---------------------------------------------------------------------
!   close model file
!---------------------------------------------------------------------
   message = 'model_init: Failed to close ' // trim(inpname)
   call handle_ncerr( nf_close( ncid ), message )       

   end subroutine model_init

   subroutine map_src_emissions( data_file, grid )
!---------------------------------------------------------------------
!   map src dataset to model grid
!---------------------------------------------------------------------

   use area_mapper, only : area_interp
   use constants_module, only : rad_per_deg, earth_radius_m
   use mapper_types

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
    type(data_file_type), intent(in) :: data_file
    type(grid_type), intent(in)      :: grid

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
    integer :: il, iu, jl, ju, k, n, stk
    integer :: status
    real    :: data_dx
    real    :: wrf_lon, wrf_lat
    real, allocatable :: raw_data(:,:)
    logical :: debug = .false.

not_STACK: &
    if( .not. data_file%is_EPA_STACK ) then
is_area_map : &
      if( grid%has_area_map ) then
!---------------------------------------------------------------------
!   area conserving interpolation
!---------------------------------------------------------------------
        call area_interp( data_file, grid%model_area_type, model_grid, wrk_emis(:,:,1), diag_level )
      else is_area_map
!---------------------------------------------------------------------
!   bilinear interpolation
!---------------------------------------------------------------------
!   form src coordinate limits
!---------------------------------------------------------------------
      if( diag_level > 100 ) then
        wrf_lon_min = minval( xlong(ids:ide,jds:jde) )
        wrf_lon_max = maxval( xlong(ids:ide,jds:jde) )
        wrf_lat_min = minval( xlat(ids:ide,jds:jde) )
        wrf_lat_max = maxval( xlat(ids:ide,jds:jde) )
        write(*,*) ' '
        write(*,'('' map_src_emissions: model lon limits = '',1p,2g25.16)') wrf_lon_min,wrf_lon_max
        write(*,'('' map_src_emissions: model lat limits = '',1p,2g25.16)') wrf_lat_min,wrf_lat_max
        write(*,*) ' '
        write(*,'('' map_src_emissions: src lon limits = '',1p,2g25.16)') grid%lon(1),grid%lon(grid%nlons)
        write(*,'('' map_src_emissions: src lat limits = '',1p,2g25.16)') grid%lat(1),grid%lat(grid%nlats)
        write(*,*) ' '
      endif
!---------------------------------------------------------------------
!   form dataset index limits
!---------------------------------------------------------------------
      write(*,*) 'map_src_emissions: count of points <,> data min,max lat = ',count(grid%ix(:,:,0) == grid%nlons )
      xndx_src(1) = minval( grid%ix(:,:,lower),mask=grid%jy(:,:,lower)>0 )
      xndx_src(2) = maxval( grid%ix(:,:,lower),mask=grid%jy(:,:,lower)>0 ) + 1
      write(*,*) 'xndx_src = ',xndx_src(:)
      write(*,*) 'map_src_emissions: count of points < data min lat = ',count(grid%jy(:,:,0) == -1)
      write(*,*) 'map_src_emissions: count of points > data max lat = ',count(grid%jy(:,:,0) == -2)
      yndx_src(1) = minval( grid%jy(:,:,lower),mask=grid%jy(:,:,lower)>0 )
      yndx_src(2) = maxval( grid%jy(:,:,lower),mask=grid%jy(:,:,lower)>0 ) + 1
      write(*,*) 'yndx_src = ',yndx_src(:)

      if( debug ) then
      write(*,*) ' '
      write(*,*) 'map_src_emissions: bilinear interp diagnostics'
      write(*,*) 'map_src_emissions: ix'
      write(*,*) grid%ix(ids,jds,:)
      write(*,*) 'map_src_emissions: ax'
      write(*,*) grid%ax(ids,jds,:)
      write(*,*) 'map_src_emissions: src lons'
      write(*,*) grid%lon(grid%ix(ids,jds,0)),grid%lon(grid%ix(ids,jds,1))
      write(*,*) 'map_src_emissions: wrf lon = ',xlong(ids,jds)
      write(*,*) 'map_src_emissions: jy'
      write(*,*) grid%jy(ids,jds,:)
      write(*,*) 'map_src_emissions: by'
      write(*,*) grid%by(ids,jds,:)
      write(*,*) 'map_src_emissions: src lats'
      write(*,*) grid%lat(grid%jy(ids,jds,0)),grid%lat(grid%jy(ids,jds,1))
      write(*,*) 'map_src_emissions: wrf lat = ',xlat(ids,jds)
      write(*,*) ' '
      do j = jds,jde
        do i = ids,ide
          if( grid%ix(i,j,lower) == grid%nlons ) then
      write(*,*) 'map_src_emissions: bilinear interp diagnostics'
      write(*,*) 'map_src_emissions: ix'
      write(*,*) grid%ix(i,j,:)
      write(*,*) 'map_src_emissions: ax'
      write(*,*) grid%ax(i,j,:)
      write(*,*) 'map_src_emissions: src lons'
      write(*,*) grid%lon(grid%ix(i,j,0)),grid%lon(grid%ix(i,j,1))
      write(*,*) 'map_src_emissions: wrf lon = ',xlong(i,j)
      write(*,*) 'map_src_emissions: jy'
      write(*,*) grid%jy(i,j,:)
      write(*,*) 'map_src_emissions: by'
      write(*,*) grid%by(i,j,:)
      write(*,*) 'map_src_emissions: src lats'
      write(*,*) grid%lat(grid%jy(i,j,0)),grid%lat(grid%jy(i,j,1))
      write(*,*) 'map_src_emissions: wrf lat = ',xlat(i,j)
            stop 'diagnostics'
          endif
        end do
      end do
      stop 'diagnostics'
      endif

!---------------------------------------------------------------------
!   allocate and read dataset variable
!---------------------------------------------------------------------
       if( .not. allocated( raw_data ) ) then
         allocate( raw_data(xndx_src(1):xndx_src(2),yndx_src(1):yndx_src(2)),stat=status )
         if( status /= 0 ) then
           write(*,*) 'map_src_emissions: allocate for raw_data failed; error = ', status
          stop 'Alloc err'
         endif
       endif
!-------------------------------------------------------------
!  transfer shifted raw data to final working array
!-------------------------------------------------------------
       do j = yndx_src(1),yndx_src(2)
         raw_data(:,j) = data_file%src_data(xndx_src(1):xndx_src(2),j)
       end do

       if( diag_level > 100 ) then
         write(*,*)  'dataset size = ',size(data_file%src_data)
         write(*,*)  'dataset min,max values = ',minval(data_file%src_data(:,:)),maxval(data_file%src_data(:,:))
       endif

!---------------------------------------------------------------------
!   set model anthro emission
!---------------------------------------------------------------------
       do j = jds,jde
         do i = ids,ide
           jl = grid%jy(i,j,0)
           if( jl > 0 ) then
             il = grid%ix(i,j,0) ; iu = il + 1
             ju = jl + 1
             wrk_sum  = raw_data(il,jl)*grid%ax(i,j,upper)*grid%by(i,j,upper) &
                      + raw_data(il,ju)*grid%ax(i,j,upper)*grid%by(i,j,lower) &
                      + raw_data(iu,jl)*grid%ax(i,j,lower)*grid%by(i,j,upper) &
                      + raw_data(iu,ju)*grid%ax(i,j,lower)*grid%by(i,j,lower)
           else
             wrk_sum  = 0.
           endif
           wrk_emis(i,j,1) = wrk_sum
         end do
       end do
       if( allocated( raw_data ) ) then
         deallocate( raw_data )
       endif
    endif is_area_map
    else not_STACK
      wrk_emis(:,:,:) = 0.
      do stk = 1,data_file%stack%nStk
        if( data_file%stack%dataMask(stk) ) then
          il = data_file%stack%mdl_i(stk) ; jl = data_file%stack%mdl_j(stk)
          k = data_file%stack%mdl_k(stk)
          wrk_emis(il,jl,k) = wrk_emis(il,jl,k) + data_file%stack%src_data(stk)
        endif
      enddo
    endif not_STACK

    end subroutine map_src_emissions

    subroutine write_wrf_emis( nfile )
!---------------------------------------------------------------------
!	... write the WRF netcdf anthro emission file
!---------------------------------------------------------------------
    use netcdf_utils, only : set_glb_atts, handle_ncerr, ngatts
    use utils, only : src_names, sub_cats

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
      integer, intent(in) :: nfile

!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
      integer :: lon_id
      integer :: lat_id
      integer :: time_id
      integer :: zdim_id
      integer :: string_id
      integer :: dims(4)
      integer :: start_ndx(4)
      integer :: length(4)
      integer :: astat
      integer :: m, cat, nt, nt1, src
      real, allocatable :: wrk_emis(:,:,:)
      character(len=132) :: message, text
      character(len=10)  :: ctime
      character(len=8)   :: cdate
      character(len=10)  :: t_string(2)
      character(len=19)  :: wrf_time

      write(num(2:3),'(i2.2)') domain
      if( serial_output ) then
        call mz2wrf_time( wrf_time, loop_date%date, loop_date%secs )
        outpname = 'wrfchemi_d' // num(2:3) // '_' // wrf_time
      else
        if( nfile == 1 ) then
          outpname = 'wrfchemi_00z_d' // num(2:3)
        else
          outpname = 'wrfchemi_12z_d' // num(2:3)
        endif
      endif
!-----------------------------------------------------------------------
!     	... create netcdf anthro emission file and enter define mode
!-----------------------------------------------------------------------
      message = 'write_emis: Failed to create ' // trim( outpname )
      call handle_ncerr( nf_create( trim( outpname ), nf_clobber, ncid ), message )

!-----------------------------------------------------------------------
!     	... define the dimensions
!-----------------------------------------------------------------------
      call handle_ncerr( nf_def_dim( ncid, 'west_east', ide, lon_id ), &
                         'write_emis: Failed to define longitude dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'south_north', jde, lat_id ), &
                         'write_emis: Failed to define latitude dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'emissions_zdim_stag', emissions_zdim_stag, zdim_id ), &
                         'write_emis: Failed to define emissions_zdim_stag dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'DateStrLen', 19, string_id ), &
                         'write_emis: Failed to define DateStrLen dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'Time', nf_unlimited, time_id ), &
                         'write_emis: Failed to create Time dimension' )
!-----------------------------------------------------------------------
!     	... define the variables
!-----------------------------------------------------------------------
      dims(1:2) = (/ string_id, time_id /)
      call handle_ncerr( nf_def_var( ncid, 'Times', nf_char, 2, dims(1:2), varid ), &
                         'write_emis: Failed to define Times variable' )
      dims(1:2) = (/ lon_id, lat_id /)
      call handle_ncerr( nf_def_var( ncid, 'XLONG', nf_float, 2, dims(1:2), varid ), &
                         'write_emis: Failed to define XLONG variable' )
      call handle_ncerr( nf_def_var( ncid, 'XLAT', nf_float, 2, dims(1:2), varid ), &
                         'write_emis: Failed to define XLAT variable' )

      dims(:) = (/ lon_id, lat_id, zdim_id, time_id /)
      do m = 1,nemis
        message = 'write_emis: Failed to define ' // trim(anthro_map(m)%emis_name)
        call handle_ncerr( nf_def_var( ncid, 'E_'//trim(anthro_map(m)%emis_name), nf_float, 4, dims, varid ), &
                           trim(message) )
      end do

!-----------------------------------------------------------------------
!   ... define variable attributes
!-----------------------------------------------------------------------
      varname = 'XLONG'
      units_attribute       = 'degree east'
      description_attribute = 'LONGITUDE, WEST IS NEGATIVE'
      stagger_attribute     = ''
      memord_attribute      = 'XY '
      coor_attribute        = ' '
      call write_attributes

      varname = 'XLAT'
      units_attribute       = 'degree north'
      description_attribute = 'LATITUDE, SOUTH IS NEGATIVE'
      stagger_attribute     = ''
      memord_attribute      = 'XY '
      coor_attribute        = ' '
      call write_attributes

      do m = 1,nemis
        varname = 'E_' // trim(anthro_map(m)%emis_name)
        if( anthro_map(m)%is_gas ) then
          units_attribute       = 'mol km^-2 hr^-1'
        else
          units_attribute       = 'ug m^-2 s^-1'
        endif
        description_attribute = 'EMISSIONS'
        stagger_attribute     = 'Z'
        memord_attribute      = 'XYZ'
        call write_attributes
is_SECTOR: &
        if( any( data_file(:)%is_EPA_SECTOR ) ) then
          message = 'write_attributes: Failed to get ' // trim(varname) // ' variable id'
          call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), message )
          do cat = 1,n_sub_cats
            if( trim(sub_cats(cat)) == trim(anthro_map(m)%emis_name) ) then
              do src = 1,nsrc_names
                message = 'write_attributes: Failed to create ' // trim(varname) // ' attribute'
                text = 'nNzWater_' // trim(src_names(src))
                call handle_ncerr( nf_put_att_int( ncid, varid, trim(text), nf_int, 1, nNzSrcWater(src,cat) ), message )
                text = 'nNzLand_' // trim(src_names(src))
                call handle_ncerr( nf_put_att_int( ncid, varid, trim(text), nf_int, 1, nNzSrcLand(src,cat) ), message )
                text = 'water_mole_per_hr_' // trim(src_names(src))
                call handle_ncerr( nf_put_att_real( ncid, varid, trim(text), nf_real, 1, totalSrcWater(src,cat) ), message )
                text = 'land_mole_per_hr_' // trim(src_names(src))
                call handle_ncerr( nf_put_att_real( ncid, varid, trim(text), nf_real, 1, totalSrcLand(src,cat) ), message )
              enddo
              exit
            endif
          enddo
        endif is_SECTOR
      end do

!-----------------------------------------------------------------------
!   ... define global attributes
!-----------------------------------------------------------------------
      message = 'global_attributes: Failed to write title'
      text    = 'Anthropogenic emissions'
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Title', len_trim(text), trim(text) ), message )
      message = 'global_attributes: Failed to write History'
      call date_and_time( cdate, ctime )
      t_string(1) = cdate(1:4) // '-' // cdate(5:6) // '-' // cdate(7:8)
      t_string(2) = ctime(1:2) // ':' // ctime(3:4)
      text    = 'Created on ' // trim(t_string(1)) // ' at ' // trim(t_string(2))
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'History', len_trim(text), trim(text) ), message )
      message = 'global_attributes: Author to write Files'
      text    = 'anthro_emis'
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Author', len_trim(text), trim(text) ), message )
!-----------------------------------------------------------------------
!     	... leave define mode
!-----------------------------------------------------------------------
      call handle_ncerr( nf_enddef( ncid ), 'write_emis: Failed to leave define mode' )

      allocate( wrk_emis(ide,jde,emissions_zdim_stag),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'write_emis: failed to allocate wrk_emis; error = ',astat
        stop 'Alloc err'
      endif

!-----------------------------------------------------------------------
!     	... write the variables
!-----------------------------------------------------------------------
      start_ndx(1:2) = (/ 1,1 /)
      length(1:2)    = (/ 19, 1 /)
      call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), &
                         'write_emis: Failed to get Times variable id' )
      if( serial_output ) then
        call handle_ncerr( nf_put_vara_text( ncid, varid, start_ndx(:2), length(:2), wrf_time ), &
                           'write_emis: Failed to write Times variable' )
      else
        wrf_time = Times(1)
        do nt = 0,11
          start_ndx(2) = nt+1
          write(wrf_time(12:13),'(i2.2)') (nfile-1)*12 + nt
          call handle_ncerr( nf_put_vara_text( ncid, varid, start_ndx(:2), length(:2), wrf_time ), &
                             'write_emis: Failed to write Times variable' )
        end do
      endif

      call handle_ncerr( nf_inq_varid( ncid, 'XLONG', varid ), &
                         'write_emis: Failed to get xlong variable id' )
      call handle_ncerr( nf_put_var_real( ncid, varid, xlong ), &
                         'write_emis: Failed to write xlong variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'XLAT', varid ), &
                         'write_emis: Failed to get xlat variable id' )
      call handle_ncerr( nf_put_var_real( ncid, varid, xlat ), &
                         'write_emis: Failed to write xlat variable' )
      start_ndx(:) = 1
      length(:)    = (/ ide, jde, emissions_zdim_stag, 1 /)
      do m = 1,nemis
        message = 'write_emis: Failed to write ' // trim(anthro_map(m)%emis_name)
        call handle_ncerr( nf_inq_varid( ncid, 'E_'//trim(anthro_map(m)%emis_name), varid ), &
                           trim(message) )
        wrk_emis(:,:,:) = 0.
        if( allocated(anthro_map(m)%emissionSrf) ) then
          wrk_emis(:,:,1) = anthro_map(m)%emissionSrf(:,:,1)
        endif
        if( allocated(anthro_map(m)%emissionLay) ) then
          wrk_emis(:,:,:) = wrk_emis(:,:,:) + anthro_map(m)%emissionLay(:,:,:)
        endif
        if( serial_output ) then
          start_ndx(4) = 1
          call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, wrk_emis ), &
                             trim(message) )
        else
          do nt = 1,12
            start_ndx(4) = nt
            call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, wrk_emis ), &
                               trim(message) )
          end do
        endif
      end do
!---------------------------------------------------------------------
!   close model file
!---------------------------------------------------------------------
      message = 'Failed to close ' // trim(outpname)
      call handle_ncerr( nf_close( ncid ), message )       

      deallocate( wrk_emis )

    end subroutine write_wrf_emis

    subroutine write_cam_emis( is_initial_time )
!---------------------------------------------------------------------
!	... write the CAM netcdf anthro emission file
!---------------------------------------------------------------------
    use netcdf_utils, only : set_glb_atts, handle_ncerr, ngatts
    use utils, only : src_names, sub_cats
    use constants_module, only : pi, earth_radius_m

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
      logical, intent(in) :: is_initial_time
!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
      real, parameter :: m2cm = 100.

      integer :: lon_id, ncol_id
      integer :: lat_id
      integer :: time_id
      integer :: zdim_id
      integer :: zdim_int_id
      integer :: string_id
      integer :: astat
      integer :: k, m, n, nt, nt1, nVert
      integer :: cat, file, src
      integer :: nShift
      integer :: dims(4)
      integer :: start_ndx(4)
      integer :: length(4)
      real    :: daysSince
      real    :: sArea
      real    :: wrkLons(ide)
      real    :: wrk_emis(ide,1,emissions_zdim_stag)
      character(len=132) :: message, text
      character(len=10)  :: ctime
      character(len=8)   :: cdate
      character(len=8)   :: vertname(2) = (/ 'surface ','vertical' /)
      character(len=10)  :: t_string(2)
      character(len=19)  :: wrf_time

      nShift = -model_grid%lon_shift
      call date_and_time( date=cdate )
      timeNdx = timeNdx + 1

emis_loop: &
      do m = 1,nemis
file_loop: &
        do file = 1,2
          if( file == 1 ) then
            if( .not. any(anthro_map(m)%is_2d(1:nsrc_names)) ) then 
              cycle file_loop
            endif
          else
            if( .not. any(anthro_map(m)%is_3d(1:nsrc_names)) ) then 
              cycle file_loop
            endif
          endif
is_first_time: &
          if( is_initial_time ) then
            if( mdl_is_CAMFV ) then
              outpname = 'emissions-EPA_' // trim(anthro_map(m)%emis_name) // '_' // &
                       trim(vertname(file)) // '_anthro_c' // cdate // '.nc'
            elseif( mdl_is_CAMSE ) then
              outpname = 'emissions-EPA_' // trim(anthro_map(m)%emis_name) // '_' // &
                     trim(vertname(file)) // '_anthro_c' // cdate // '_conus_30x8.nc'
            endif
!-----------------------------------------------------------------------
!     	... create netcdf anthro emission file and enter define mode
!-----------------------------------------------------------------------
            message = 'write_cam_emis: Failed to create ' // trim( outpname )
            call handle_ncerr( nf_create( trim( outpname ), nf_clobber, ncid ), message )

!-----------------------------------------------------------------------
!     	... define the dimensions
!-----------------------------------------------------------------------
            if( mdl_is_CAMFV ) then
              call handle_ncerr( nf_def_dim( ncid, 'lon', model_grid%nlons, lon_id ), &
                                 'write_cam_emis: Failed to define longitude dimension' )
              call handle_ncerr( nf_def_dim( ncid, 'lat', model_grid%nlats, lat_id ), &
                                 'write_cam_emis: Failed to define latitude dimension' )
            elseif( mdl_is_CAMSE ) then
              call handle_ncerr( nf_def_dim( ncid, 'ncol', model_grid%nlons, ncol_id ), &
                                 'write_cam_emis: Failed to define ncol dimension' )
            endif
            if( file == 2 ) then
              call handle_ncerr( nf_def_dim( ncid, 'altitude', model_grid%nlevs, zdim_id ), &
                                 'write_cam_emis: Failed to define altitude dimension' )
              call handle_ncerr( nf_def_dim( ncid, 'altitude_int', model_grid%nlevs+1, zdim_int_id ), &
                                 'write_cam_emis: Failed to define altitude interface dimension' )
            endif
            call handle_ncerr( nf_def_dim( ncid, 'time', nf_unlimited, time_id ), &
                               'write_cam_emis: Failed to create Time dimension' )
!-----------------------------------------------------------------------
!     	... define the variables
!-----------------------------------------------------------------------
            dims(1) = time_id
            call handle_ncerr( nf_def_var( ncid, 'time', nf_float, 1, dims, varid ), &
                               'write_cam_emis: Failed to define time variable' )
            if( mdl_is_CAMFV ) then
              dims(1) = lon_id
              call handle_ncerr( nf_def_var( ncid, 'lon', nf_float, 1, dims, varid ), &
                                 'write_cam_emis: Failed to define lon variable' )
              dims(1) = lat_id
              call handle_ncerr( nf_def_var( ncid, 'lat', nf_float, 1, dims, varid ), &
                                 'write_cam_emis: Failed to define lat variable' )
            elseif( mdl_is_CAMSE ) then
              dims(1) = ncol_id
              call handle_ncerr( nf_def_var( ncid, 'lon', nf_float, 1, dims, varid ), &
                                 'write_cam_emis: Failed to define lon variable' )
              call handle_ncerr( nf_def_var( ncid, 'lat', nf_float, 1, dims, varid ), &
                                 'write_cam_emis: Failed to define lat variable' )
            endif

            if( file == 1 ) then
              if( mdl_is_CAMFV ) then
                n = 3
                dims(:n) = (/ lon_id, lat_id, time_id /)
              elseif( mdl_is_CAMSE ) then
                n = 2
                dims(:n) = (/ ncol_id, time_id /)
              endif
            elseif( file == 2 ) then
              dims(1) = zdim_id
              call handle_ncerr( nf_def_var( ncid, 'altitude', nf_float, 1, dims, varid ), &
                                 'write_cam_emis: Failed to define altitude variable' )
              dims(1) = zdim_int_id
              call handle_ncerr( nf_def_var( ncid, 'altitude_int', nf_float, 1, dims, varid ), &
                                 'write_cam_emis: Failed to define altitude_int variable' )
              if( mdl_is_CAMFV ) then
                n = 4
                dims(:n) = (/ lon_id, lat_id, zdim_id, time_id /)
              elseif( mdl_is_CAMSE ) then
                n = 3
                dims(:n) = (/ ncol_id, zdim_id, time_id /)
              endif
            endif

            message = 'write_cam_emis: Failed to define emiss_anthro variable'
            call handle_ncerr( nf_def_var( ncid, 'emiss_anthro', nf_float, n, dims, varid ), &
                               trim(message) )
!-----------------------------------------------------------------------
!   ... define variable attributes
!-----------------------------------------------------------------------
            if( mdl_is_CAMFV ) then
              varname = 'lon'
            elseif( mdl_is_CAMSE ) then
              varname = 'lon'
            endif
            units_attribute       = 'degrees_east'
            description_attribute = 'Longitude'
            call write_attributes

            if( mdl_is_CAMFV ) then
              varname = 'lat'
            elseif( mdl_is_CAMSE ) then
              varname = 'lat'
            endif
            units_attribute       = 'degrees_north'
            description_attribute = 'Latitude'
            call write_attributes

            varname = 'time'
            units_attribute       = 'Time'
            description_attribute = 'days since 1900-01-01 00:00:00'
            call write_attributes

            if( file == 2 ) then
              varname = 'altitude'
              units_attribute       = 'km'
              description_attribute = 'Altitude midpoints'
              call write_attributes

              varname = 'altitude_int'
              units_attribute       = 'km'
              description_attribute = 'Altitude interfaces'
              call write_attributes
            endif

            varname = 'emiss_anthro'
            if( file == 1 ) then
              units_attribute       = 'molecules/cm2/s'
            else
              units_attribute       = 'molecules/cm3/s'
            endif
            description_attribute = trim(anthro_map(m)%emis_name) // ' EPA anthropogenic emissions'
            call write_attributes
!-----------------------------------------------------------------------
!   ... define global attributes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     	... leave define mode
!-----------------------------------------------------------------------
            call handle_ncerr( nf_enddef( ncid ), 'write_cam_emis: Failed to leave define mode' )
          else is_first_time
            if( mdl_is_CAMFV ) then
              outpname = 'emissions-EPA_' // trim(anthro_map(m)%emis_name) // '_' // &
                         trim(vertname(file)) // '_anthro_c' // cdate // '.nc'
            elseif( mdl_is_CAMSE ) then
              outpname = 'emissions-EPA_' // trim(anthro_map(m)%emis_name) // '_' // &
                     trim(vertname(file)) // '_anthro_c' // cdate // '_conus_30x8.nc'
            endif
            message = 'write_cam_emis: Failed to open file: ' // trim(outpname)
            call handle_ncerr( nf_open( trim(outpname), nf_write, ncid ), message )       
!           call handle_ncerr( nf_inq_dimid( ncid, 'time', time_id ), message )
!           message = 'write_cam_emis: Failed to get time dimension'
!           call handle_ncerr( nf_inq_dimlen( ncid, time_id, timeNdx ), message )
          endif is_first_time

!-----------------------------------------------------------------------
!     	... write the variables
!-----------------------------------------------------------------------
          start_ndx(1) = timeNdx
          length(1)    = 1
          daysSince = diffdat( 19000101, 0, loop_date%date, loop_date%secs ) 
          call handle_ncerr( nf_inq_varid( ncid, 'time', varid ), &
                             'write_cam_emis: Failed to get time variable id' )
          call handle_ncerr( nf_put_vara_real( ncid, varid, (/timeNdx/), (/1/), daysSince ), &
                             'write_cam_emis: Failed to write time variable' )

          if( is_initial_time ) then
            if( mdl_is_CAMFV ) then
              call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), &
                                 'write_cam_emis: Failed to get lon variable id' )
              call handle_ncerr( nf_put_var_real( ncid, varid, model_grid%mdl_lons(:,1) ), &
                                 'write_cam_emis: Failed to write lon variable' )
              call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), &
                                 'write_cam_emis: Failed to get lat variable id' )
              call handle_ncerr( nf_put_var_real( ncid, varid, model_grid%mdl_lats(1,:) ), &
                                 'write_cam_emis: Failed to write lat variable' )
            elseif( mdl_is_CAMSE ) then
              call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), &
                                 'write_cam_emis: Failed to get lon variable id' )
              wrkLons(:) = model_grid%mdl_poly(:)%cntr_lon 
              where( model_grid%mdl_poly(:)%cntr_lon < 0. )
                 wrkLons(:) = model_grid%mdl_poly(:)%cntr_lon +360.
              endwhere
              call handle_ncerr( nf_put_var_real( ncid, varid, wrkLons ), &
                                 'write_cam_emis: Failed to write lon variable' )
              call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), &
                                 'write_cam_emis: Failed to get lat variable id' )
              call handle_ncerr( nf_put_var_real( ncid, varid, model_grid%mdl_poly(:)%cntr_lat ), &
                                 'write_cam_emis: Failed to write lat variable' )
            endif
            if( mdl_is_CAM .and. file == 2 ) then
              call handle_ncerr( nf_inq_varid( ncid, 'altitude_int', varid ), &
                                 'write_cam_emis: Failed to get altitude_int variable id' )
              call handle_ncerr( nf_put_var_real( ncid, varid, 1.e-3*model_grid%z_at_w(1,:,1) ), &
                                 'write_cam_emis: Failed to write altitude_int variable' )
              call handle_ncerr( nf_inq_varid( ncid, 'altitude', varid ), &
                                 'write_cam_emis: Failed to get altitude variable id' )
              call handle_ncerr( nf_put_var_real( ncid, varid, &
                                 (/ (1.e-3*dz*(real(k)+.5),k=0,model_grid%nlevs-1) /) ), &
                                 'write_cam_emis: Failed to write altitude variable' )
            endif
          endif
          start_ndx(:) = 1
          if( file == 1 ) then
            if( mdl_is_CAMFV ) then
              n = 3
              length(:n)    = (/ ide, jde, 1 /)
            elseif( mdl_is_CAMSE ) then
              n = 2
              length(:n)    = (/ model_grid%nPolygons, 1 /)
            endif
            nVert = 1
            wrk_emis(:,:,:) = anthro_map(m)%emissionSrf(:,:,:)
          else
            if( mdl_is_CAMFV ) then
              n = 4
              length(:n)    = (/ ide, jde, emissions_zdim_stag, 1 /)
            elseif( mdl_is_CAMSE ) then
              n = 3
              length(:n)    = (/ model_grid%nPolygons, emissions_zdim_stag, 1 /)
            endif
            nVert = emissions_zdim_stag
            wrk_emis(:,:,:) = anthro_map(m)%emissionLay(:,:,:)
          endif
          message = 'write_cam_emis: Failed to get ' // trim(anthro_map(m)%emis_name) // 'id'
          call handle_ncerr( nf_inq_varid( ncid, 'emiss_anthro', varid ), trim(message) )
!---------------------------------------------------------------------
!   divide by model cell area and lon shift back to CAM grid
!---------------------------------------------------------------------
          if( mdl_is_CAMFV ) then
            do k = 1,emissions_zdim_stag
              do j = 1,model_grid%nlats
                if( file == 1 ) then
                  wrk_emis(:,j,k) = wrk_emis(:,j,k)/model_grid%mdl_cell_area(j,1)
                else
                  wrk_emis(:,j,k) = wrk_emis(:,j,k)/(5.e3*model_grid%mdl_cell_area(j,1))
                endif
              enddo
            enddo
            wrk_emis(:,j,k) = cshift( wrk_emis(:,j,k),nShift )
          elseif( mdl_is_CAMSE ) then
            sArea = m2cm*m2cm*earth_radius_m*earth_radius_m/(4.*pi)
            do k = 1,emissions_zdim_stag
              if( file == 1 ) then
                wrk_emis(:,1,k) = wrk_emis(:,1,k)/(model_grid%mdl_poly(:)%area*sArea)
                if( trim(anthro_map(m)%emis_name) == 'CO' ) then
                  if( any(wrk_emis(:,1,1) < 0.) ) then
                    write(*,*) 'write_cam_emis: there are ',count(wrk_emis(:,1,1) < 0.),' CO srf emis values < 0'
                    stop 'DEBUG'
                  endif
                endif
              else
                wrk_emis(:,1,k) = wrk_emis(:,1,k)/(m2cm*dz*model_grid%mdl_poly(:)%area*sArea)
              endif
            enddo
          endif
          start_ndx(n) = timeNdx
          message = 'write_cam_emis: Failed to write ' // trim(anthro_map(m)%emis_name)
          call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(1:n), length(1:n), wrk_emis(:,:,1:nVert) ), &
                             trim(message) )
!---------------------------------------------------------------------
!   close model file
!---------------------------------------------------------------------
          message = 'write_cam_emis: Failed to close ' // trim(outpname)
          call handle_ncerr( nf_close( ncid ), message )       
        enddo file_loop
      enddo emis_loop

   end subroutine write_cam_emis

   subroutine write_attributes
!---------------------------------------------------------------------
!   write common variable attributes
!---------------------------------------------------------------------

      use netcdf_utils, only : handle_ncerr

      message = 'write_attributes: Failed to get ' // trim(varname) // ' variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), message )
      message = 'write_attributes: Failed to create ' // trim(varname) // ' attribute'
      if( mdl_is_WRF ) then
        call handle_ncerr( nf_put_att_text( ncid, varid, 'stagger', &
                                            len_trim(stagger_attribute), trim(stagger_attribute) ), message )
        if( coor_attribute /= ' ' ) then
          call handle_ncerr( nf_put_att_text( ncid, varid, 'coordinates', &
                                              len_trim(coor_attribute), trim(coor_attribute) ), message )
        endif
        ii = 104
        call handle_ncerr( nf_put_att_int( ncid, varid, 'FieldType', nf_int, 1, ii ), message )
        call handle_ncerr( nf_put_att_text( ncid, varid, 'MemoryOrder', 3, memord_attribute ), message )
        call handle_ncerr( nf_put_att_text( ncid, varid, 'description', &
                                            len_trim(description_attribute), trim(description_attribute) ), message )
      elseif( mdl_is_CAM ) then
        call handle_ncerr( nf_put_att_text( ncid, varid, 'longname', &
                                            len_trim(description_attribute), trim(description_attribute) ), message )
        if( trim(varname) == 'time' ) then
          call handle_ncerr( nf_put_att_text( ncid, varid, 'calender', 9, 'Gregorian' ), message )
        endif
      endif
      call handle_ncerr( nf_put_att_text( ncid, varid, 'units', &
                                          len_trim(units_attribute), trim(units_attribute) ), message )

   end subroutine write_attributes

   end program anthro_emis
