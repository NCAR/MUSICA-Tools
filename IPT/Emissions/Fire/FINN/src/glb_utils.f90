
   module glb_utils
!------------------------------------------------------------
!  global domain utilities
!------------------------------------------------------------

   use wrf_utils, only : PROJ_LATLON, PI, RAD_PER_DEG

   implicit none

   private
   public :: glb_file
   public :: write_glb_fire_file
   public :: glb_file_final

   real, parameter :: kg2g   = 1.e3              ! kilograms -> grams
   real, parameter :: m2cm   = 100.              ! meters -> centimeters
   real, parameter :: REARTH = 6.37122E6         ! (m)
   real, parameter :: RSQ    = m2cm*m2cm*REARTH*REARTH  ! cm^2
   real, parameter :: AVOG   = 6.022e23          ! molecules/mole
   real, parameter :: s_per_day = 86400.         ! seconds/day
   real, parameter :: const0 = AVOG/s_per_day
   real, parameter :: const1 = const0/RSQ
   real, parameter :: ZERO      = 0.
   real, parameter :: ONE80     = 180.
   real, parameter :: THREE60   = 360.
 
   integer           :: time_ndx = 1
   real, allocatable :: area(:)
   real, allocatable :: accum_emiss(:,:,:)
   character(len=132), allocatable :: outpname(:)

!------------------------------------------------------------------
!     include files
!------------------------------------------------------------------
   include 'netcdf.inc'

   contains

   subroutine glb_file( UtilAttrs, proj, n_spc, n_days, start_date, end_date, &
                        n_files, fire_directory, fire_filename, output_timing, n_mnths, mdlGrd )
!---------------------------------------------------------------------
!   setup global output file
!---------------------------------------------------------------------
   
   use wrf_utils,  only : proj_info
   use fire_types, only : model_grid_type
   use utils,      only : wrf2fire_type, diffdat, geth_newdate, MaxChars
   use attr_types, only : UtilAttributes

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
    integer, intent(in)             :: n_spc
    integer, intent(in)             :: n_days
    integer, intent(in)             :: n_files
    integer, intent(in)             :: n_mnths
    character(len=*), intent(in)    :: start_date
    character(len=*), intent(in)    :: end_date
    character(len=*), intent(in)    :: fire_directory
    character(len=*), intent(in)    :: output_timing
    character(len=*), intent(in)    :: fire_filename(n_files)
    type(proj_info), intent(inout)  :: proj
    type(model_grid_type), intent(in) :: mdlGrd
    type(UtilAttributes), intent(in)  :: UtilAttrs
!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer            :: m, n
   integer            :: astat
   integer            :: ncid
   integer            :: dimid
   integer            :: varid
   integer            :: ncol_id
   integer            :: lon_id
   integer            :: lat_id
   integer            :: time_id
   integer            :: date_id
   integer            :: nlon, nlat
   integer            :: yr,mnth,day
   integer            :: dims(3)
   integer, allocatable  :: dates(:)
   real               :: dlon, dlat
   real               :: dlon2, dlat2
   real               :: lonmin, latmin
   real               :: lonmax, latmax
   real               :: wrk
   real               :: wrkLons(mdlGrd%nPolygons)
   real, allocatable  :: lons(:)
   real, allocatable  :: lats(:)
   real, allocatable  :: times(:)
   character(len=3)   :: num(2)
   character(len=10)  :: sinceDate, wrk_date
   character(len=8)   :: sDate, eDate
   character(len=MaxChars) :: message
   logical            :: is_LR
   logical            :: is_CAMSE, isnt_CAMSE

   is_CAMSE = mdlGrd%mdl_is_CAMSE
   isnt_CAMSE = .not. is_CAMSE
!---------------------------------------------------------------------
!   initialize map projection
!---------------------------------------------------------------------
is_not_CAMSE: &
   if( isnt_CAMSE ) then
     is_LR = .true.
     latmin = -90.; lonmin = 0.
     latmax =  90.; lonmax = THREE60
     select case( trim(UtilAttrs%resol) )
       case( '0.5x0.5' ) 
         is_LR = .false.
         dlon = .5; dlat = .5
         lonmin = lonmin + .5*dlon
         lonmax = lonmax - .5*dlon
         latmin = latmin + .5*dlat
         latmax = latmax - .5*dlat
         nlon = int((lonmax - lonmin)/dlon) + 1
         nlat = int((latmax - latmin)/dlat) + 1
       case( '1x1' ) 
         is_LR = .false.
         dlon = 1.; dlat = 1.
         lonmin = lonmin + .5*dlon
         lonmax = lonmax - .5*dlon
         latmin = latmin + .5*dlat
         latmax = latmax - .5*dlat
         nlon = int((lonmax - lonmin)/dlon) + 1
         nlat = int((latmax - latmin)/dlat) + 1
       case( 'T42LR' )
         nlon = 128; nlat = 64
         dlon = THREE60/real(nlon)
         dlat = ONE80/real(nlat-1)
       case( 'T63LR' )
         nlon = 192; nlat = 96
         dlon = THREE60/real(nlon)
         dlat = ONE80/real(nlat-1)
       case( 'T85LR' )
         nlon = 256; nlat = 128
         dlon = THREE60/real(nlon)
         dlat = ONE80/real(nlat-1)
       case( 'T106LR' )
         nlon = 320; nlat = 160
         dlon = THREE60/real(nlon)
         dlat = ONE80/real(nlat-1)
       case( 'T170LR' )
         nlon = 512; nlat = 256
         dlon = THREE60/real(nlon)
         dlat = ONE80/real(nlat-1)
       case( '1.9x2.5' )
         nlon = 144; nlat = 96
         dlon = THREE60/real(nlon)
         dlat = ONE80/real(nlat-1)
       case( '0.9x1.25' )
         nlon = 288; nlat = 192
         dlon = THREE60/real(nlon)
         dlat = ONE80/real(nlat-1)
       case default
         write(*,*) 'glb_file: Requested horizontal resolution, ',trim(UtilAttrs%resol)
         write(*,*) '          is not presently supported'
         stop 'Namelist error'
     end select

!---------------------------------------------------------------------
!   initialize map projection to latlon
!---------------------------------------------------------------------
     proj%code   = PROJ_LATLON
     proj%latinc = dlat
     proj%loninc = dlon
     proj%knowni = 1.
     proj%knownj = 1.
     proj%lat1   = latmin
     proj%lon1   = lonmin
     proj%ide    = nlon
     proj%jde    = nlat
     proj%init   = .true.

!---------------------------------------------------------------------
!   allocate working arrays
!---------------------------------------------------------------------
     allocate( area(nlat),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'glb_file: failed to allocate area array; error = ',astat
       stop 'Alloc error'
     endif
     allocate( lons(nlon),lats(nlat),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'glb_file: failed to allocate lons,lats  arrays; error = ',astat
       stop 'Alloc error'
     endif
     if( output_timing == 'monthly' ) then
       allocate( accum_emiss(nlon,nlat,n_spc),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'glb_file: failed to allocate accum_emiss array; error = ',astat
         stop 'Alloc error'
       endif
       accum_emiss(:,:,:) = ZERO
     endif
   else is_not_CAMSE
     allocate( lons(mdlGrd%nPolygons),lats(mdlGrd%nPolygons),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'glb_file: failed to allocate lons,lats  arrays; error = ',astat
       stop 'Alloc error'
     endif
     if( output_timing == 'monthly' ) then
       allocate( accum_emiss(mdlGrd%nPolygons,n_spc,1),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'glb_file: failed to allocate accum_emiss array; error = ',astat
         stop 'Alloc error'
       endif
       accum_emiss(:,:,:) = ZERO
     endif
   endif is_not_CAMSE

!---------------------------------------------------------------------
!   allocate working arrays
!---------------------------------------------------------------------
   allocate( outpname(n_spc),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'glb_file: failed to allocate outpname array; error = ',astat
     stop 'Alloc error'
   endif
   
   if( output_timing == 'daily' ) then
     allocate( times(n_days),dates(n_days),stat=astat )
   else
     allocate( times(n_mnths),dates(n_mnths),stat=astat )
   endif
   if( astat /= 0 ) then
     write(*,*) 'glb_file: failed to allocate times,dates  arrays; error = ',astat
     stop 'Alloc error'
   endif

   if( isnt_CAMSE ) then
!---------------------------------------------------------------------
!   setup coordinate, time arrays
!---------------------------------------------------------------------
     do m = 1,nlon
       lons(m) = lonmin + real(m-1)*dlon
     end do
     wrk = 4.e4*pi*rearth*rearth*sin( .5*dlat*rad_per_deg )/real(nlon)
     do m = 1,nlat
       lats(m) = latmin + real(m-1)*dlat
     end do
     do m = 2,nlat-1
       area(m) = const0 / (wrk * cos( lats(m)*rad_per_deg ) )
     end do
     if( is_LR ) then
       area(1) = const0 / (wrk * cos( .25*dlat*rad_per_deg ) )
       area(nlat) = area(1)
     else
       area(1)    = const0 / (wrk * cos( lats(1)*rad_per_deg ) )
       area(nlat) = const0 / (wrk * cos( lats(nlat)*rad_per_deg ) )
     endif
     sinceDate = '1990-01-01'
   else
     sinceDate = '1750-01-01'
   endif
   if( output_timing == 'daily' ) then
     call geth_newdate( wrk_date, start_date, -24 )
     do m = 1,n_days
       call geth_newdate( wrk_date, wrk_date, 24 )
       read(wrk_date(1:4),*) yr
       read(wrk_date(6:7),*) mnth
       read(wrk_date(9:10),*) day
       dates(m) = 10000*yr + 100*mnth + day
       times(m) = diffdat( (/ sinceDate,wrk_date /) )
     end do
   else
     wrk_date = start_date
     read(wrk_date(1:4),*) yr
     read(wrk_date(6:7),*) mnth
     read(wrk_date(9:10),*) day
     dates(1) = 10000*yr + 100*mnth + day
     times(1) = diffdat( (/ sinceDate,wrk_date /) )
     do m = 2,n_mnths
       mnth = mnth + 1
       if( m > 12 ) then
         mnth = 1
         yr   = yr + 1
       endif
       dates(m) = 10000*yr + 100*mnth + day
       write(num(1),'(i3)') 100+mnth
       write(num(2),'(i3)') 100+day
       write(wrk_date(1:4),'(i4)') yr
       wrk_date(5:) = '-' // num(1)(2:3) // '-' // num(2)(2:3)
       times(m) = diffdat( (/ sinceDate,wrk_date /) )
     end do
   endif

species_loop : &
   do m = 1,n_spc
     sDate = start_date(1:4) // start_date(6:7) //start_date(9:10)
     eDate = end_date(1:4) // end_date(6:7) //end_date(9:10)
     outpname(m) = 'emissions-finnv' // trim( UtilAttrs%FinnVers ) // '_' // trim(wrf2fire_type(m)%wrf_name) &
                   // '_' // trim( UtilAttrs%EmisType ) // '_' // sDate // '-' // eDate // '_' // trim(UtilAttrs%resol) // '.nc'
!    outpname(m) = 'emissions_' // trim(wrf2fire_type(m)%wrf_name) // '_' // trim(resol) // '.nc'
!-----------------------------------------------------------------------
!  create netcdf bio emission file and enter define mode
!-----------------------------------------------------------------------
     message = 'glb_file: Failed to create ' // trim( outpname(m) )
     call handle_ncerr( nf_create( trim( outpname(m) ), nf_clobber, ncid ), message )

!-----------------------------------------------------------------------
!  define the dimensions
!-----------------------------------------------------------------------
     if( isnt_CAMSE ) then
       call handle_ncerr( nf_def_dim( ncid, 'lon', nlon, lon_id ), &
                          'glb_file: Failed to define longitude dimension' )
       call handle_ncerr( nf_def_dim( ncid, 'lat', nlat, lat_id ), &
                          'glb_file: Failed to define latitude dimension' )
     else
       call handle_ncerr( nf_def_dim( ncid, 'ncol', mdlGrd%nPolygons, ncol_id ), &
                          'glb_file: Failed to define ncol dimension' )
     endif
     if( output_timing == 'daily' ) then
       call handle_ncerr( nf_def_dim( ncid, 'time', nf_unlimited, time_id ), &
                          'write_cam_emis: Failed to create time dimension' )
     else
       call handle_ncerr( nf_def_dim( ncid, 'time', nf_unlimited, time_id ), &
                          'write_cam_emis: Failed to create time dimension' )
     endif

!-----------------------------------------------------------------------
!  define the variables, set attributes
!-----------------------------------------------------------------------
     if( isnt_CAMSE ) then
       dims(1) = lon_id
     else
       dims(1) = ncol_id
     endif
     call handle_ncerr( nf_def_var( ncid, 'lon', nf_float, 1, dims(1), varid ), &
                        'glb_file: Failed to define lon variable' )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'units', 12, 'degrees_east' ), &
                        'glb_file: Failed to create lon units attribute' )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'long_name', 9, 'Longitude' ), &
                        'glb_file: Failed to create lon long_name attribute' )
     if( isnt_CAMSE ) then
       dims(1) = lat_id
     else
       dims(1) = ncol_id
     endif
     call handle_ncerr( nf_def_var( ncid, 'lat', nf_float, 1, dims(1), varid ), &
                        'glb_file: Failed to define lat variable' )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'units', 13, 'degrees_north' ), &
                        'glb_file: Failed to create lat units attribute' )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'long_name', 8, 'Latitude' ), &
                        'glb_file: Failed to create lat long_name attribute' )
     dims(1) = time_id
     call handle_ncerr( nf_def_var( ncid, 'time', nf_float, 1, dims(1), varid ), &
                        'glb_file: Failed to define time variable' )
     message = 'days since ' // sinceDate // ' 00:00:00'
     call handle_ncerr( nf_put_att_text( ncid, varid, 'units', len_trim(message), trim(message) ), &
                        'glb_file: Failed to create time units attribute' )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'long_name', 4, 'Time' ), &
                        'glb_file: Failed to create time long_name attribute' )
     message = 'Gregorian'
     call handle_ncerr( nf_put_att_text( ncid, varid, 'calendar', len_trim(message), trim(message) ), &
                        'glb_file: Failed to create time calendar attribute' )

     call handle_ncerr( nf_def_var( ncid, 'date', nf_int, 1, dims(1), varid ), &
                        'glb_file: Failed to define date variable' )
     message = 'YYYYMMDD'
     call handle_ncerr( nf_put_att_text( ncid, varid, 'units', len_trim(message), trim(message) ), &
                        'glb_file: Failed to create date units attribute' )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'long_name', 4, 'Date' ), &
                        'glb_file: Failed to create date long_name attribute' )

     if( isnt_CAMSE ) then
       dims(1:3) = (/ lon_id, lat_id, time_id /)
       call handle_ncerr( nf_def_var( ncid, 'fire', nf_float, 3, dims, varid ), &
                          'glb_file: Failed to define fire variable' )
     else
       dims(1:2) = (/ ncol_id, time_id /)
       call handle_ncerr( nf_def_var( ncid, 'bb', nf_float, 2, dims, varid ), &
                          'glb_file: Failed to define fire variable' )
     endif
!---------------------------------------------------------------------
!  write the species units, long name, and full mapping
!---------------------------------------------------------------------
     message = trim( wrf2fire_type(m)%units )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'units', len_trim(message), trim(message) ), &
                        'glb_file: Failed to create fire units attribute' )
     message = trim( wrf2fire_type(m)%map )
     call handle_ncerr( nf_put_att_text( ncid, varid, 'map', len_trim(message), trim(message) ), &
                        'glb_file: Failed to create fire map attribute' )
     message = trim( wrf2fire_type(m)%wrf_name ) // ' fire emissions'
     call handle_ncerr( nf_put_att_text( ncid, varid, 'long_name', len_trim(message), trim(message) ), &
                        'glb_file: Failed to create fire long_name attribute' )
!---------------------------------------------------------------------
!  write the global attributes
!---------------------------------------------------------------------
     call global_attributes( ncid, UtilAttrs%FinnVers, UtilAttrs%Model, UtilAttrs%resol, proj, n_files, &
                             fire_directory, fire_filename, output_timing )
!-----------------------------------------------------------------------
!  leave define mode
!-----------------------------------------------------------------------
     call handle_ncerr( nf_enddef( ncid ), 'glb_file: Failed to leave define mode' )
!---------------------------------------------------------------------
!  write the coordinate, time arrays
!---------------------------------------------------------------------
     message = 'glb_utils: Failed to get lon id'
     call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), trim(message) )
     message = 'glb_utils: Failed to write lon variable'
     if( isnt_CAMSE ) then
       call handle_ncerr( nf_put_var_real( ncid, varid, lons ), trim(message) )
     else
       wrkLons(:) = mdlGrd%mdl_poly(:)%cntr_lon
       where( wrkLons(:) < 0 )
         wrkLons(:) = wrkLons(:) + THREE60
       endwhere
       call handle_ncerr( nf_put_var_real( ncid, varid, wrkLons ), trim(message) )
     endif
     message = 'glb_utils: Failed to get lat id'
     call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), trim(message) )
     message = 'glb_utils: Failed to write lat variable'
     if( isnt_CAMSE ) then
       call handle_ncerr( nf_put_var_real( ncid, varid, lats ), trim(message) )
     else
       call handle_ncerr( nf_put_var_real( ncid, varid, mdlGrd%mdl_poly(:)%cntr_lat ), trim(message) )
     endif

     message = 'glb_utils: Failed to get time id'
     call handle_ncerr( nf_inq_varid( ncid, 'time', varid ), trim(message) )
     message = 'glb_utils: Failed to write time variable'
     do n = 1,size(times)
       call handle_ncerr( nf_put_vara_real( ncid, varid, (/n/), 1, times(n) ), trim(message) )
     enddo
     message = 'glb_utils: Failed to get date id'
     call handle_ncerr( nf_inq_varid( ncid, 'date', varid ), trim(message) )
     message = 'glb_utils: Failed to write date variable'
     do n = 1,size(times)
       call handle_ncerr( nf_put_vara_int( ncid, varid, (/n/), 1, dates(n) ), trim(message) )
     enddo
!    call handle_ncerr( nf_put_var_int( ncid, varid, dates ), trim(message) )
!---------------------------------------------------------------------
!  close file
!---------------------------------------------------------------------
     message = 'Failed to close ' // trim(outpname(m))
     call handle_ncerr( nf_close( ncid ), message )       
   end do species_loop

   deallocate( lons, lats, times, dates )
   
   end subroutine glb_file

   subroutine write_glb_fire_file( file, n_spc, m1, m2, lon_ndx, &
                                   lat_ndx, n_fire_spc, fire_emissions, proj, output_timing, &
                                   wrk_date, mdlGrd )
!---------------------------------------------------------------------
!  write global fire emission files
!---------------------------------------------------------------------

    use utils,      only : wrf2fire_type, diag_level, geth_newdate
    use fire_types, only : model_grid_type
    use wrf_utils,  only : proj_info

!---------------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------------
    integer, intent(in) :: file
    integer, intent(in) :: n_spc
    integer, intent(in) :: n_fire_spc
    integer, intent(in) :: m1, m2
    integer, intent(in) :: lon_ndx(m1:m2)
    integer, intent(in) :: lat_ndx(m1:m2)
    real, intent(in)    :: fire_emissions(n_fire_spc,m1:m2)
    character(len=*), intent(in) :: output_timing
    character(len=*), intent(in) :: wrk_date
    type(proj_info), intent(in)  :: proj
    type(model_grid_type), intent(in)  :: mdlGrd

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
   integer :: k, m, n
   integer :: n1, n2
   integer :: astat
   integer :: ncid
   integer :: varid
   integer :: next_day
   integer :: start_ndx(3), length(3)
   real    :: wrk_sum
   real    :: avrg_factor
   real, allocatable :: wrk_emiss(:,:)
   character(len=10)  :: next_date
   character(len=128) :: message
   logical :: do_output
   logical :: is_CAMSE, isnt_CAMSE

   is_CAMSE = mdlGrd%mdl_is_CAMSE
   isnt_CAMSE = .not. is_CAMSE

   if( isnt_CAMSE ) then
     allocate( wrk_emiss(proj%ide,proj%jde),stat=astat )
   else
     allocate( wrk_emiss(mdlGrd%nPolygons,1),stat=astat )
   endif
   if( astat /= 0 ) then
     write(*,*) 'write_glb_fire_file: failed to allocate wrk_emiss array; error = ',astat
     stop 'Alloc error'
   endif

   do_output = output_timing == 'daily'
   if( .not. do_output ) then
     call geth_newdate( next_date, wrk_date, 24 )
     read(next_date(9:10),*) next_day
     do_output = next_day == 1
   endif

spc_loop : &
   do n = 1,n_spc
     wrk_emiss(:,:) = ZERO
time_loop : &
     do m = m1,m2
       n1 = lon_ndx(m)
       if( isnt_CAMSE ) then
         n2 = lat_ndx(m)
       else
         n2 = 1
       endif
       wrk_sum = ZERO
       do k = 1,wrf2fire_type(n)%fire_cnt
         wrk_sum = wrk_sum &
                 + fire_emissions(wrf2fire_type(n)%fire_ndx(k,file),m) &
                   *wrf2fire_type(n)%fire_wght(k)
       end do
       if( isnt_CAMSE ) then
         wrk_sum = wrk_sum * area(n2)
         if( wrf2fire_type(n)%is_aerosol ) then
           wrk_sum = kg2g * wrk_sum / wrf2fire_type(n)%m_wght
         endif
       else
         wrk_sum = wrk_sum * const1 / mdlGrd%mdl_poly(n1)%area
         if( wrf2fire_type(n)%is_aerosol ) then
           wrk_sum = kg2g * wrk_sum / wrf2fire_type(n)%m_wght
         endif
       endif
       wrk_emiss(n1,n2) = wrk_emiss(n1,n2) + wrk_sum
     end do time_loop
     if( output_timing == 'monthly' ) then
       if( isnt_CAMSE ) then
         accum_emiss(:,:,n) = accum_emiss(:,:,n) + wrk_emiss(:,:)
       else
         accum_emiss(:,n,1) = accum_emiss(:,n,1) + wrk_emiss(:,1)
       endif
     endif

     if( do_output ) then
!---------------------------------------------------------------------
!   open global output file
!---------------------------------------------------------------------
       message = 'write_glb_fire_file: Failed to open ' // trim(outpname(n))
       call handle_ncerr( nf_open( trim(outpname(n)), nf_write, ncid ), trim(message) )       

       if( isnt_CAMSE ) then
         start_ndx(:) = (/ 1,1,time_ndx /)
         length(:)    = (/ proj%ide,proj%jde,1 /)
       else
         start_ndx(1:2) = (/ 1,time_ndx /)
         length(1:2)    = (/ mdlGrd%nPolygons,1 /)
       endif
!---------------------------------------------------------------------
!   write global output file
!---------------------------------------------------------------------
       message = 'write_glb_fire_file: Failed to get fire variable id'
       if( isnt_CAMSE ) then
         call handle_ncerr( nf_inq_varid( ncid, 'fire', varid ), trim(message) )
       else
         call handle_ncerr( nf_inq_varid( ncid, 'bb', varid ), trim(message) )
       endif
       message = 'write_glb_fire_file: Failed to write fire variable'
       if( output_timing == 'daily' ) then
         if( isnt_CAMSE ) then
           call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, wrk_emiss ), trim(message) )
         else
           call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(1:2), length(1:2), wrk_emiss(:,1) ), trim(message) )
         endif
       else
         read(wrk_date(9:10),*) next_day
         avrg_factor = 1./real(next_day)
         if( isnt_CAMSE ) then
           accum_emiss(:,:,n) = accum_emiss(:,:,n) * avrg_factor
           call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, accum_emiss(:,:,n) ), trim(message) )
           accum_emiss(:,:,n) = ZERO
         else
           accum_emiss(:,n,1) = accum_emiss(:,n,1) * avrg_factor
           call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(1:2), length(1:2), accum_emiss(:,n,1) ), trim(message) )
           accum_emiss(:,n,1) = ZERO
         endif
       endif
!---------------------------------------------------------------------
!  close file
!---------------------------------------------------------------------
       message = 'Failed to close ' // trim(outpname(n))
       call handle_ncerr( nf_close( ncid ), message )       
     endif
   end do spc_loop

   if( do_output ) then
     time_ndx = time_ndx + 1
   endif

   deallocate( wrk_emiss )

   end subroutine write_glb_fire_file

   subroutine global_attributes( ncid, FinnVers, Model, resol, proj, n_files, &
                                 fire_directory, fire_filename, output_timing )
!---------------------------------------------------------------------
!   write file global attributes
!---------------------------------------------------------------------
    use wrf_utils, only : proj_info
!---------------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------------
   integer, intent(in)          :: ncid
   integer, intent(in)          :: n_files
   character(len=*), intent(in) :: FinnVers
   character(len=*), intent(in) :: Model
   character(len=*), intent(in) :: resol
   character(len=*), intent(in) :: fire_directory
   character(len=*), intent(in) :: output_timing
   character(len=*), intent(in) :: fire_filename(n_files)
   type(proj_info), intent(in)  :: proj
!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
   integer            :: n
   character(len=10)  :: ctime
   character(len=8)   :: cdate
   character(len=132) :: text
   character(len=132) :: message
   character(len=10)  :: t_string(2)

   message = 'global_attributes: Failed to write title'
!  text    = 'FINNv15 daily fire emissions'
   text    = 'FINNv' // trim(FinnVers) // ' daily fire emissions'
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'title', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Failed to write title'
   text    = trim(Model)
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'TargetModel', len_trim(text), trim(text) ), message )
   text    = trim(resol)
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'TargetGrid', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Failed to write History'
   text    = trim(output_timing)
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Frequency', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Failed to write Frequency'
   call date_and_time( cdate, ctime )
   t_string(1) = cdate(1:4) // '-' // cdate(5:6) // '-' // cdate(7:8)
   t_string(2) = ctime(1:2) // ':' // ctime(3:4)
   text    = 'Created on ' // trim(t_string(1)) // ' at ' // trim(t_string(2))
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'History', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Failed to write Directory'
   text    = trim( fire_directory )
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Directory', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Failed to write Files'
   text    = trim(fire_filename(1))
   if( n_files > 1 ) then
     text = trim(fire_filename(1)) // ','
     do n = 2,n_files
       if( n /= n_files ) then
         text(len_trim(text)+1:) = trim(fire_filename(n)) // ','
       else
         text(len_trim(text)+1:) = trim(fire_filename(n))
       endif
     end do
   endif
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Files', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Author to write Files'
   text    = 'fire_emis'
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Author', len_trim(text), trim(text) ), message )

   end subroutine global_attributes

   subroutine glb_file_final
!---------------------------------------------------------------------
!  module "destructor"
!---------------------------------------------------------------------

   if( allocated( area ) ) then
     deallocate( area )
   endif
   if( allocated( outpname ) ) then
     deallocate( outpname )
   endif
   if( allocated( accum_emiss ) ) then
     deallocate( accum_emiss )
   endif

   end subroutine glb_file_final

   subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!	... netcdf error handling routine
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ret
   character(len=*), intent(in) :: mes

   if( ret /= nf_noerr ) then
      write(*,*) 'handle_ncerr: ',trim(mes)
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

   end module glb_utils
