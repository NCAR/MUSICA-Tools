
   module data_file_utils

   use anthro_types, only : linsize
   use netcdf_utils, only : handle_ncerr

   implicit none

!---------------------------------------------------------------------
!  include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   private

   public :: data_file_init
   public :: get_src_time_ndx
   public :: read_src_data
   public :: tinterp_src_data
   public :: anthro_dir
   public :: data_yrs_offset
   public :: src_lon_dim_name, src_lon_var_name
   public :: src_lat_dim_name, src_lat_var_name
   public :: stack_dx

   integer                :: data_yrs_offset     = 0         ! years( + or - )
   real                   :: stack_dx            = 0.
   logical                :: is_EPA, is_EPA_STACK, is_EPA_SECTOR
   character(len=linsize) :: anthro_dir
   character(len=64)      :: src_lon_dim_name = 'lon'
   character(len=64)      :: src_lon_var_name = 'lon'
   character(len=64)      :: src_lat_dim_name = 'lat'
   character(len=64)      :: src_lat_var_name = 'lat'

   CONTAINS

   subroutine data_file_init( data_file, start_output, cat_var_prefix, cat_var_suffix, &
                              domain, dx, ids, ide, jds, jde, &
                              currDate, mdlGrid, emissions_zdim_stag )
!---------------------------------------------------------------------
!  initialize data file type
!---------------------------------------------------------------------

   use constants_module, only : rad_per_deg, earth_radius_m
   use mapper_types
   use area_mapper,      only : area_interp_init
   use area_mapper,      only : ijll, lm_2_xy, set_wghts, pnt_in_quad
   use area_mapper,      only : xlong => lon, xlat => lat
   use anthro_types,     only : data_file_type, model_grid_type, dates
   use utils,            only : diag_level
   use misc_definitions_module, only : PROJ_LC_EPA
   use stack, only : stkFileRead, stkFileInit

!---------------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ids, ide
   integer, intent(in) :: jds, jde
   integer, intent(in) :: domain
   integer, intent(inout) :: currDate
   integer, intent(in) :: emissions_zdim_stag
   real, intent(in)    :: dx
   character(len=*), intent(in)        :: cat_var_prefix
   character(len=*), intent(in)        :: cat_var_suffix
   type(data_file_type), intent(inout) :: data_file
   type(dates), intent(in) :: start_output
   type(model_grid_type), intent(in) :: mdlGrid

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
    real, parameter :: m2km = 1.e-3
    real, parameter :: km2m = 1.e3

    integer :: i, j
    integer :: ip1, jp1
    integer :: il, iu, jl, ju
    integer :: n, varid, dimid
    integer :: nlon_src, nlat_src
    integer :: lon_dim_id, lat_dim_id
    integer :: lon_var_id, lat_var_id
    integer :: ncid
    integer :: astat, ierr, status
    integer :: xndx_src(2)
    integer :: yndx_src(2)
    integer :: maxind(1)
    integer, allocatable :: var_dimids(:)
    real    :: roundoff, checksum
    real    :: d1, d2, ds1
    real    :: data_dx, epa_dx
    real    :: wrf_lon, wrf_lat
    real    :: wrf_lon_min, wrf_lat_min
    real    :: wrf_lon_max, wrf_lat_max
    real    :: xy(2)
    real    :: wghts(4)
    real    :: quad_x(4), quad_y(4)
    real    :: src_dom_x(4), src_dom_y(4)
    real, allocatable :: src_lons(:), src_lons_2d(:,:)
    real, allocatable :: src_lats(:), src_lats_2d(:,:)
    real, allocatable :: wrk_emis(:,:)
    character(len=132) :: varname
    character(len=132) :: message
    logical :: new_grid
    logical :: has_area_map
    logical :: has_lon_shift
    logical :: reorder_lons
    logical :: reorder_lats
    logical :: found
    logical :: is2d
    type(grid_type), pointer :: grid
    type(proj_info) :: proj

    roundoff = 100.*epsilon(roundoff)
    is_EPA = data_file%is_EPA
    is_EPA_SECTOR = data_file%is_EPA_SECTOR
    is_EPA_STACK  = data_file%is_EPA_STACK

    write(*,*) ' '
    write(*,*) 'data_file_init: Initializing type for src emission file ' // trim(data_file%filename)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   open src emission dataset file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    message = 'data_file_init: Failed to open ' // trim(data_file%filespec)
    call handle_ncerr( nf_open( trim(data_file%filespec), nf_noclobber, ncid ), message )       
    if( .not. is_EPA_STACK ) then
!---------------------------------------------------------------------
!   get src dataset dimesions
!---------------------------------------------------------------------
      message = 'data_file_init: Failed to get lon dimension id'
      call handle_ncerr( nf_inq_dimid( ncid, trim(src_lon_dim_name), lon_dim_id ), message )
      message = 'data_file_init: Failed to get lon dimension'
      call handle_ncerr( nf_inq_dimlen( ncid, lon_dim_id, nlon_src ), message )
      message = 'data_file_init: Failed to get lat dimension id'
      call handle_ncerr( nf_inq_dimid( ncid, trim(src_lat_dim_name), lat_dim_id ), message )
      message = 'data_file_init: Failed to get lat dimension'
      call handle_ncerr( nf_inq_dimlen( ncid, lat_dim_id, nlat_src ), message )
      write(*,*) 'data_file_init:  nlon_src, nlat_src = ',nlon_src,nlat_src
    endif
not_EPA: &
    if( .not. is_EPA ) then
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   check if lon, lat variables have two horizontal dimensions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      message = 'data_file_init: Failed to get lon variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim(src_lon_var_name), lon_var_id ), message )
      message = 'data_file_init: Failed to get lon variable dim count'
      call handle_ncerr( nf_inq_varndims( ncid, lon_var_id, n ), message )
      allocate( var_dimids(n),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'data_file_init: Failed to allocate var_dimds; error = ',astat
        stop 'allocate failed'
      endif
      message = 'data_file_init: Failed to get lon variable dim ids'
      call handle_ncerr( nf_inq_vardimid( ncid, lon_var_id, var_dimids ), message )
      is2d = any( var_dimids(:) == lon_dim_id ) .and. any( var_dimids(:) == lat_dim_id ) 
      deallocate( var_dimids )
    endif not_EPA
not_STACK: &
    if( .not. is_EPA_STACK ) then
!---------------------------------------------------------------------
!   allocate data longitudes
!---------------------------------------------------------------------
      if( allocated( src_lons_2d ) ) then
        deallocate( src_lons_2d )
      endif
      allocate( src_lons_2d(nlon_src,nlat_src),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'data_file_init: Failed to allocate src_lons_2d; error = ',astat
        stop 'allocate failed'
      endif
!---------------------------------------------------------------------
!   allocate data latitudes
!---------------------------------------------------------------------
      if( allocated( src_lats_2d ) ) then
        deallocate( src_lats_2d )
      endif
      allocate( src_lats_2d(nlon_src,nlat_src),stat=ierr )
      if( ierr /= 0 ) then
        write(*,*) 'data_file_init: Failed to allocate src_lats_2d; error = ',ierr
        stop 'allocate failed'
      endif
    endif not_STACK
!---------------------------------------------------------------------
!   read src longitude variable
!---------------------------------------------------------------------
not_EPA_a: &
   if( .not. is_EPA ) then
      message = 'data_file_init: Failed to read lon variable'
      if( is2d ) then
        call handle_ncerr( nf_get_var_real( ncid, lon_var_id, src_lons_2d ), message )
      else
        call handle_ncerr( nf_get_var_real( ncid, lon_var_id, src_lons_2d(1:nlon_src,1) ), message )
        do j = 2,nlat_src
          src_lons_2d(1:nlon_src,j) = src_lons_2d(1:nlon_src,1)
        enddo
      endif
!---------------------------------------------------------------------
!   read src latitude variable
!---------------------------------------------------------------------
      message = 'data_file_init: Failed to get lat variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim(src_lat_var_name), varid ), message )
      message = 'data_file_init: Failed to read lat variable'
      if( is2d ) then
        call handle_ncerr( nf_get_var_real( ncid, varid, src_lats_2d), message )
      else
        call handle_ncerr( nf_get_var_real( ncid, varid, src_lats_2d(1,1:nlat_src) ), message )
        do i = 2,nlon_src
          src_lats_2d(i,1:nlat_src) = src_lats_2d(1,1:nlat_src)
        enddo
      endif
   elseif( is_EPA_SECTOR ) then
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   if EPA data src, get global attributes needed for projection mapping
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      is2d = .true.
      proj%code = PROJ_LC_EPA
      proj%re_m = earth_radius_m
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'P_ALP', proj%truelat1 ), message )       
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'P_BET', proj%truelat2 ), message )       
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'P_GAM', proj%stdlon ), message )       
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'XCENT', proj%xcent ), message )       
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'YCENT', proj%ycent ), message )       
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'XORIG', proj%xOrigin ), message )       
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'YORIG', proj%yOrigin ), message )       
      call handle_ncerr( nf_get_att_real( ncid, nf_global, 'XCELL', proj%dx ), message )       
      proj%xOrigin = m2km * proj%xOrigin
      proj%yOrigin = m2km * proj%yOrigin
      proj%dx      = m2km * proj%dx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   set src lon,lat at cell center
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do j = 1,nlat_src
        do i = 1,nlon_src
          call ijll( (real(i)-.5)*proj%dx, (real(j)-.5)*proj%dx, proj, &
                      src_lats_2d(i,j), src_lons_2d(i,j) )
        enddo
      enddo
   elseif( is_EPA_STACK ) then
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   read stack group data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     call stkFileRead( data_file%stack )
   endif not_EPA_a

not_STACK_a: &
   if( .not. is_EPA_STACK ) then
!---------------------------------------------------------------------
!   check that lons are monotonically increasing and reorder if not
!---------------------------------------------------------------------
      reorder_lons = .false.
      if( all( src_lons_2d(2:nlon_src,1:nlat_src) < src_lons_2d(1:nlon_src-1,1:nlat_src) ) ) then
        reorder_lons = .true.
        do j = 1,nlat_src
          do i = 1,nlon_src/2
            d1 = src_lons_2d(i,j)
            iu = nlon_src - i + 1
            src_lons_2d(i,j)  = src_lons_2d(iu,j)
            src_lons_2d(iu,j) = d1
          enddo
        enddo
      endif
!---------------------------------------------------------------------
!   make sure src_lons are in [-180,180]
!---------------------------------------------------------------------
     has_lon_shift = .false.
     if( all( src_lons_2d(:,:) >= 0. ) ) then
       do j = 1,nlat_src
         if( count( src_lons_2d(:,j) >= 180. ) > 0 ) then
           d2 = 360. + src_lons_2d(1,j) - src_lons_2d(nlon_src,j)
           d1 = src_lons_2d(2,j) - src_lons_2d(1,j)
           if( abs( d2 - d1 ) > 1.e-3*d1 ) then
             write(*,*) 'data_file_init: can not map data longitudes to WRF grid'
             stop
           endif
           src_lons_2d(:,j) = cshift( src_lons_2d(:,j),nlon_src/2 ) 
           src_lons_2d(1:nlon_src/2,j) = src_lons_2d(1:nlon_src/2,j) - 360.
           has_lon_shift = .true.
         endif
       enddo
     endif
!---------------------------------------------------------------------
!   check that lats are monotinicity increasing and reorder if not
!---------------------------------------------------------------------
     reorder_lats = .false.
     do i = 1,nlon_src
       if( all( src_lats_2d(i,2:nlat_src) < src_lats_2d(i,1:nlat_src-1) ) ) then
         reorder_lats = .true.
         do j = 1,nlat_src/2
           d1 = src_lats_2d(i,j)
           ju = nlat_src - j + 1
           src_lats_2d(i,j)  = src_lats_2d(i,ju)
           src_lats_2d(i,ju) = d1
         end do
       endif
     enddo
   endif not_STACK_a

!---------------------------------------------------------------------
!   determine interpolation type; bilinear or area conserving
!---------------------------------------------------------------------
    if( .not. is_EPA ) then
      data_dx = earth_radius_m * (src_lats_2d(1,2) - src_lats_2d(1,1)) * rad_per_deg
!   elseif( is_EPA_SECTOR ) then
    else
      data_dx = proj%dx*km2m
    endif
    if( .not. is_EPA_STACK ) then
      has_area_map = data_dx <= dx
      data_file%dx = data_dx
      write(*,*) 'data_file_init: data_dx,mdl_dx,has_area_map = ',data_dx,dx,has_area_map
    else
      has_area_map = .false.
    endif

!-------------------------------------------------------------
!   check for match against prior datasets
!-------------------------------------------------------------
   if( grid_cnt >= grid_max ) then
     write(*,*) 'data_file_init: reached grid cache max: ',grid_max
     stop
   endif
   grid_ndx = 0
   new_grid = .true.
   if( .not. is_EPA_STACK ) then
     do n = 1,grid_cnt
       if( grid_specs(n)%nlons /= nlon_src .or. grid_specs(n)%nlats /= nlat_src ) then
         cycle
       endif
       if( any( grid_specs(n)%lon_2d(:,:) /= src_lons_2d(:,:) ) ) then
         cycle
       endif
       if( any( grid_specs(n)%lat_2d(:,:) /= src_lats_2d(:,:) ) ) then
         cycle
       endif
       grid_ndx = n
       new_grid = .false.
       exit
     end do
   else
     new_grid = .false.
   endif
!-------------------------------------------------------------
!   new data grid to cache
!-------------------------------------------------------------
has_new_grid : &
   if( new_grid ) then
     grid_cnt = grid_cnt + 1
     grid => grid_specs(grid_cnt)
     grid%nlons = nlon_src
     grid%nlats = nlat_src
     grid%has_area_map  = has_area_map
     grid%has_lon_shift = has_lon_shift
     grid%reorder_lons  = reorder_lons
     grid%reorder_lats  = reorder_lats
     grid_ndx = grid_cnt
     allocate( grid%lon_2d(nlon_src,nlat_src), &
               grid%xedge_2d(nlon_src+1,nlat_src+1), &
               grid%xedge_2d_ij(nlon_src+1,nlat_src+1),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'data_file_init: Failed to allocate grid_lons; error = ',ierr
       stop 'allocate failed'
     endif
     allocate( grid%lat_2d(nlon_src,nlat_src), &
               grid%yedge_2d(nlon_src+1,nlat_src+1), &
               grid%yedge_2d_ij(nlon_src+1,nlat_src+1),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'data_file_init: Failed to allocate grid_lats; error = ',ierr
       stop 'allocate failed'
     endif

     grid%lon_2d(:,:) = src_lons_2d(:,:)
     grid%lat_2d(:,:) = src_lats_2d(:,:)
     write(*,*) 'data_file_init: file ' // trim(data_file%filename),' is a new grid'

is_area_map : &
     if( has_area_map ) then
!---------------------------------------------------------------------
!   form src longitude, latitude edges
!---------------------------------------------------------------------
not_EPA_b: &
       if( .not. is_EPA ) then
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   "interior" data cell edges
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do j = 1,nlat_src-1
           jp1 = j + 1
           do i = 1,nlon_src-1
             ip1 = i + 1
             quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                            grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
             quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                            grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
             call lm_2_xy( quad_x, quad_y, (/ .5,.5 /), xy )
             grid%xedge_2d(ip1,jp1) = xy(1)
             grid%yedge_2d(ip1,jp1) = xy(2)
           enddo
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   "border" data cell edges
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   lower edge
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         j = 1 ; jp1 = j + 1
         do i = 1,nlon_src-1
           ip1 = i + 1
           quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                          grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
           quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                          grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
           call lm_2_xy( quad_x, quad_y, (/ .5,-.5 /), xy )
           grid%xedge_2d(ip1,j) = xy(1)
           grid%yedge_2d(ip1,j) = xy(2)
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   upper edge
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         j = nlat_src-1 ; jp1 = j + 1
         do i = 1,nlon_src-1
           ip1 = i + 1
           quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                          grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
           quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                          grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
           call lm_2_xy( quad_x, quad_y, (/ .5,1.5 /), xy )
           grid%xedge_2d(ip1,jp1+1) = xy(1)
           grid%yedge_2d(ip1,jp1+1) = xy(2)
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   left edge
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         i = 1 ; ip1 = i + 1
         do j = 1,nlat_src-1
           jp1 = j + 1
           quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                          grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
           quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                          grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
           call lm_2_xy( quad_x, quad_y, (/ -.5,.5 /), xy )
           grid%xedge_2d(i,jp1) = xy(1)
           grid%yedge_2d(i,jp1) = xy(2)
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   right edge
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         i = nlon_src-1 ; ip1 = i + 1
         do j = 1,nlat_src-1
           jp1 = j + 1
           quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                          grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
           quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                          grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
           call lm_2_xy( quad_x, quad_y, (/ 1.5,.5 /), xy )
           grid%xedge_2d(ip1+1,jp1) = xy(1)
           grid%yedge_2d(ip1+1,jp1) = xy(2)
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   "corner" data cell edges
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   lower, left corner
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         i = 1 ; ip1 = i + 1
         j = 1 ; jp1 = j + 1
         quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                        grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
         quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                        grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
         call lm_2_xy( quad_x, quad_y, (/ -.5,-.5 /), xy )
         grid%xedge_2d(i,j) = xy(1)
         grid%yedge_2d(i,j) = xy(2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   lower, right corner
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         i = nlon_src-1 ; ip1 = i + 1
         j = 1 ; jp1 = j + 1
         quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                        grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
         quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                        grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
         call lm_2_xy( quad_x, quad_y, (/ 1.5,-.5 /), xy )
         grid%xedge_2d(ip1+1,j) = xy(1)
         grid%yedge_2d(ip1+1,j) = xy(2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   upper, right corner
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         i = nlon_src-1 ; ip1 = i + 1
         j = nlat_src-1 ; jp1 = j + 1
         quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                        grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
         quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                        grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
         call lm_2_xy( quad_x, quad_y, (/ 1.5,1.5 /), xy )
         grid%xedge_2d(ip1+1,jp1+1) = xy(1)
         grid%yedge_2d(ip1+1,jp1+1) = xy(2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   upper, left corner
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         i = 1 ; ip1 = i + 1
         j = nlat_src-1 ; jp1 = j + 1
         quad_x(:) = (/ grid%lon_2d(i,j),grid%lon_2d(ip1,j), &
                        grid%lon_2d(ip1,jp1),grid%lon_2d(i,jp1) /)
         quad_y(:) = (/ grid%lat_2d(i,j),grid%lat_2d(ip1,j), &
                        grid%lat_2d(ip1,jp1),grid%lat_2d(i,jp1) /)
         call lm_2_xy( quad_x, quad_y, (/ -.5,1.5 /), xy )
         grid%xedge_2d(i,jp1+1) = xy(1)
         grid%yedge_2d(i,jp1+1) = xy(2)
       else not_EPA_b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   set src cell edges (lon,lat)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do j = 0,nlat_src
           jp1 = j + 1
           do i = 0,nlon_src
             ip1 = i + 1
             call ijll( real(i)*proj%dx, real(j)*proj%dx, proj, d1, d2 )
             grid%xedge_2d(ip1,jp1) = real(d2,kind=8)
             grid%yedge_2d(ip1,jp1) = real(d1,kind=8)
           enddo
         enddo
       endif not_EPA_b

       write(*,'(''data_file_init: xcen_src(1,2)  = '',1p,2g22.15)') .5*(src_lons_2d(1:2,1)+src_lons_2d(2:3,1))
       write(*,'(''data_file_init: xedge(1,2) = '',1p,2g22.15)') grid_specs(grid_cnt)%xedge_2d(1:2,1)
       write(*,'(''data_file_init: dx = '',1pg22.15)') int( 1./(src_lons_2d(2,1) - src_lons_2d(1,1)) )

       write(*,'(''data_file_init: nlon_src,nlat_src = '',i6,1x,i6)') nlon_src,nlat_src
       write(*,'(''data_file_init: ycen_src  = '',1p,2g22.15)') .5*(src_lats_2d(1,nlat_src-2:nlat_src-1)+src_lats_2d(1,nlat_src-1:nlat_src) )
       write(*,'(''data_file_init: yedge_src = '',1p,2g22.15)') grid_specs(grid_cnt)%yedge_2d(1,nlat_src-1:nlat_src)

       allocate( grid%model_area_type(ide,jde),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'data_file_init; failed to allocate model_area_type: error = ',astat
         stop 'Alloc error'
       endif
       grid%model_area_type(:,:)%has_data = .false.
       grid%model_area_type(:,:)%active_dcell_cnt = 0
       grid%model_area_type(:,:)%total_dcell_cnt  = 0
       grid%model_area_type(:,:)%interior_dcell_cnt = 0
       grid%model_area_type(:,:)%partial_dcell_cnt  = 0
!---------------------------------------------------------------------
!   setup area conserving interpolation
!---------------------------------------------------------------------
       if( .not. allocated( wrk_emis ) ) then
         allocate( wrk_emis(ide,jde),stat=status )
         if( status /= 0 ) then
           write(*,*) 'data_file_init: allocate for wrk_emis failed; error = ',ierr
           stop 'Alloc error'
         endif
       endif
       call area_interp_init( mdlGrid%mdl_is_CAMFV, grid, grid%model_area_type, diag_level )
     else is_area_map
       ierr = 0
       allocate( grid%ax(ids:ide,jds:jde,0:1), &
                 grid%by(ids:ide,jds:jde,0:1),stat=status )
       ierr = ierr + status
       allocate( grid%ix(ids:ide,jds:jde,0:1), &
                 grid%jy(ids:ide,jds:jde,0:1),stat=status )
       ierr = ierr + status
       if( ierr /= 0 ) then
         write(*,*) 'data_file_init: allocate failed for ax,by,ix,jy'
         stop 'Alloc_err'
       endif

mlat_loop: &
       do j = jds,jde
mlon_loop: &
         do i = ids,ide
           xy(:) = (/ xlong(i,j),xlat(i,j) /)
           found = .false.
slat_loop: do jl = 1,nlat_src-1
             jp1 = jl + 1
             do il = 1,nlon_src-1
               ip1 = il + 1
               quad_x(:) = (/ src_lons_2d(il,jl),src_lons_2d(ip1,jl), &
                              src_lons_2d(ip1,jp1),src_lons_2d(il,jp1) /)
               quad_y(:) = (/ src_lats_2d(il,jl),src_lats_2d(ip1,jl), &
                              src_lats_2d(ip1,jp1),src_lats_2d(il,jp1) /)
               found = pnt_in_quad( xy, quad_x, quad_y )
               if( found ) then
                 grid%ix(i,j,0) = il
                 grid%jy(i,j,0) = jl
                 call set_wghts( quad_x, quad_y, xy, wghts ) 
                 grid%ax(i,j,0:1) = wghts(2:1:-1)
                 grid%by(i,j,0:1) = wghts(4:3:-1)
                 checksum = (wghts(1)+wghts(2))*(wghts(3)+wghts(4))
                 if( abs( checksum - 1. ) > roundoff ) then
                   write(*,*) 'data_file_init: wght(',i,',',j,') = ',sum(wghts(:))
                 endif
                 exit slat_loop
               endif
             enddo
           enddo slat_loop
           if( .not. found ) then
             grid%jy(i,j,0) = 0
           endif
         enddo mlon_loop
       enddo mlat_loop
       write(*,*) 'data_file_init: # model cells outside data grid = ',count(grid%jy(:,:,0) == 0)
       write(*,*) 'data_file_init: min x data cell = ',minval(grid%ix(:,:,0),mask=grid%jy(:,:,0)>0)
       write(*,*) 'data_file_init: max x data cell = ',maxval(grid%ix(:,:,0),mask=grid%jy(:,:,0)>0)
       write(*,*) 'data_file_init: min y data cell = ',minval(grid%jy(:,:,0),mask=grid%jy(:,:,0)>0)
       write(*,*) 'data_file_init: max y data cell = ',maxval(grid%jy(:,:,0),mask=grid%jy(:,:,0)>0)
     endif is_area_map
   endif has_new_grid
 
   if( is_EPA_STACK ) then
     call stkFileInit( data_file%stack, mdlGrid, emissions_zdim_stag )
   endif

!---------------------------------------------------------------------
!   if not set, get src molecular weight
!---------------------------------------------------------------------
   if( data_file%molecw == 0. ) then
     if( .not. is_EPA ) then
       call get_molecwght( ncid, cat_var_prefix, cat_var_suffix, data_file%molecw, data_file%missing_value )
     else
       data_file%molecw = 1.
     endif
   endif

!---------------------------------------------------------------------
!   set category units
!---------------------------------------------------------------------
   if( domain == 1 ) then
     if( .not. is_EPA ) then
       call get_units( ncid, cat_var_prefix, cat_var_suffix, data_file%molecw, &
                       data_file%con_fac, m2km*mdlGrid%dx, mdlGrid )
     else
       if( is_EPA_SECTOR ) then
         epa_dx = m2km*data_dx  
       else
         epa_dx = m2km*mdlGrid%dx
       endif
       call get_units( ncid, cat_var_prefix, cat_var_suffix, data_file%molecw, &
                       data_file%con_fac, epa_dx, mdlGrid )
     endif
   endif

   data_file%grid_ndx = grid_ndx

!---------------------------------------------------------------------
!   initialize timing for data file
!---------------------------------------------------------------------
   data_file%lo_tndx = 0
   data_file%hi_tndx = 0
   data_file%lo_buf_ndx = 1
   data_file%hi_buf_ndx = 2
   data_file%ncid_lo    = 0
   data_file%ncid_hi    = 0
   data_file%in_gap     = .false.
   data_file%t_interp   = .false.

!---------------------------------------------------------------------
!   initialize data file timing
!---------------------------------------------------------------------
   call data_file_timing_init( ncid, start_output, currDate, data_file ) 

   if( data_file%lo_tndx > 0 ) then
     data_file%lo_tndx = 0
     data_file%ncid_lo = ncid
     data_file%ncid_hi = ncid
   endif

   n = count( data_file%cat_active(:) )
   if( .not. is_EPA_STACK ) then
     allocate( data_file%emis(nlon_src,nlat_src,2,n), &
               data_file%src_data(nlon_src,nlat_src),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'data_file_init: allocate for emis,src_data failed; error = ',astat
       stop 'Alloc error'
     endif
   else
     allocate( data_file%stack%emis(data_file%stack%nStk,2,n), &
               data_file%stack%src_data(data_file%stack%nStk),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'data_file_init: allocate for stack emis,src_data failed; error = ',astat
       stop 'Alloc error'
     endif
   endif

   end subroutine data_file_init

   subroutine data_file_timing_init( ncid, interp_time, currDate, data_file )
!---------------------------------------------------------------
!     initialize data file timing info
!---------------------------------------------------------------

   use anthro_types, only : data_file_type, dates
   use mo_calendar, only  : diffdat

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
   integer, intent(in)          :: ncid
   integer, intent(inout)       :: currDate
   type(data_file_type), intent(inout) :: data_file
   type(dates),          intent(in)    :: interp_time

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
   integer :: status, tstat
   integer :: dimid, varid
   integer :: i, j, k, n 

!---------------------------------------------------------------
!     read times
!---------------------------------------------------------------
   call read_src_times( ncid, data_file%date, data_file%secs, data_file%ntimes, data_file%is_EPA )
   if( data_file%is_EPA ) then
!    call get_flnm_date( data_file%filename, date=data_file%flnmDate )
!    data_file%date(1:data_file%ntimes) = data_file%flnmDate
     data_file%date(1:data_file%ntimes) = currDate
   endif

   data_file%lo_tndx = lotim( interp_time%date, interp_time%secs, data_file )
   if( data_file%lo_tndx <= 0 ) then
!---------------------------------------------------------------
!     check either prior or next data file
!---------------------------------------------------------------
     status = nf_close( ncid )
     if( data_file%lo_tndx == 0 ) then
       call next_flnm( data_file%filename, incr=.false. )
     else
       call next_flnm( data_file%filename, incr=.true. )
     endif
!---------------------------------------------------------------
!     open the prior or next mozart netCDF file
!---------------------------------------------------------------
     if( .not. data_file%is_EPA ) then
       data_file%filespec = trim(anthro_dir) // '/' // trim(data_file%filename)
     else
       data_file%filespec = trim(anthro_dir) // '/' // trim(data_file%src_name) // '/' // trim(data_file%filename)
     endif
     status = nf_open( trim(data_file%filespec), nf_nowrite, ncid )
     if( status /= nf_noerr ) then
       write(*,*) 'data_file_timing_init : Failed to open ',trim(data_file%filespec)
       call handle_error( status )
     end if
     write(*,*) 'data_file_timing_init: opened ',trim(data_file%filespec)
     call read_src_times( ncid, data_file%date, data_file%secs, data_file%ntimes, data_file%is_EPA )
     if( data_file%is_EPA ) then
       call get_flnm_date( data_file%filename, date=data_file%flnmDate )
       data_file%date(1:data_file%ntimes) = data_file%flnmDate
     endif
     data_file%lo_tndx = lotim( interp_time%date, interp_time%secs, data_file )
     if( data_file%lo_tndx == 0 ) then
       write(*,*) 'data_file_timing_init: time ',interp_time%date,' ',interp_time%secs
       write(*,*) '                 is before '
       write(*,*) '                 ',data_file%date(1),' ',data_file%secs(1)
       write(*,*) '                 the first time in the data file ',trim(data_file%filename)
       status = nf_close( ncid )
       stop
     end if
   end if

   end subroutine data_file_timing_init

   subroutine read_src_times( ncid, date, secs, ntimes, dfile_is_EPA )
!---------------------------------------------------------------
!     read times from current src input file
!---------------------------------------------------------------

      use mo_calendar, only : newdate

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in)  :: ncid
      integer, intent(out) :: ntimes
      integer, allocatable :: date(:)
      integer, allocatable :: secs(:)
      logical, intent(in)  :: dfile_is_EPA

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: m
      integer :: dimid, varid
      integer :: istat, status
      integer :: sl, su
      integer :: yr, mn, dy
      integer :: hr, min, sec
      real, allocatable :: times(:)
      character(len=132) :: units_text
      character(len=8)   :: time_name
      logical :: got_it

!---------------------------------------------------------------
!     get time dimension and allocate arrays
!---------------------------------------------------------------
      if( .not. dfile_is_EPA ) then
        time_name = 'time'
      else
        time_name = 'TSTEP'
      endif
      status = nf_inq_dimid( ncid, trim(time_name), dimid )
      if( status /= nf_noerr )  then
        if( .not. dfile_is_EPA ) then
          time_name = 'date'
          status = nf_inq_dimid( ncid, trim(time_name), dimid )
          if( status /= nf_noerr )  then
            write(*,*) 'read_src_times: can not find either time or date dimension'
            call handle_error( status )
          endif
        else
          write(*,*) 'read_src_times: can not find ',trim(time_name),' dimension'
          call handle_error( status )
        endif
      endif

      status = nf_inq_dimlen( ncid, dimid, ntimes )
      if( status /= nf_noerr ) call handle_error( status )
      if( is_EPA ) then
        ntimes = max( ntimes - 1,1 )
      endif

      if( allocated( date ) ) then
        deallocate( date, secs )
      end if
      allocate( date(ntimes), secs(ntimes), stat=status )
      if( status /= 0 ) then
        write(*,*) 'failed to allocate date, secs; error = ',status
        stop 'Alloc err'
      end if

!---------------------------------------------------------------
!     check for "standard" date, datesec variables
!---------------------------------------------------------------
      got_it = .false.
      status = nf_inq_varid( ncid, 'date', varid )
      if( status == nf_noerr ) then
        status = nf_inq_varid( ncid, 'datesec', m )
        if( status == nf_noerr ) then
          status = nf_get_var_int( ncid, varid, date )
          if( status /= nf_noerr )  call handle_error( status )
          status = nf_get_var_int( ncid, m, secs )
          if( status /= nf_noerr )  call handle_error( status )
          got_it = .true.
        endif
      endif
!---------------------------------------------------------------
!     try to get date, datesec from time variable
!---------------------------------------------------------------
      if( .not. got_it ) then
not_EPA: &
        if( .not. dfile_is_EPA ) then
          status = nf_inq_varid( ncid, trim(time_name), varid )
          if( status /= nf_noerr )  call handle_error( status )
          units_text = ' '
          status = nf_get_att_text( ncid, varid, 'units', units_text )
          if( status /= nf_noerr )  call handle_error( status )
          sl = scan( units_text, '0123456789' )
          if( sl > 0 ) then
            su = index( units_text, '-', back=.true. )
            if( units_text(sl+4:sl+4) == '-' .and. su == sl+7 ) then
              read(units_text(sl:sl+3),*,iostat=istat) yr
              if( istat == 0 ) then
                read(units_text(sl+5:sl+6),*,iostat=istat) mn
                if( istat == 0 ) then
                  read(units_text(sl+8:sl+9),*,iostat=istat) dy
                endif
              endif
              if( istat == 0 ) then
                allocate( times(ntimes),stat=istat )
                if( istat /= 0 ) then
                  write(*,*) 'read_src_times: failed to allocate times; error = ',istat
                  stop 'Alloc err'
                endif
                status = nf_get_var_real( ncid, varid, times )
                if( status /= nf_noerr )  call handle_error( status )
                date(:) = dy + 100*(mn + 100*yr)
                do m = 1,ntimes
                  date(m) = newdate( date(m), int( times(m) ) )
                  secs(m) = 86400*(times(m) - aint(times(m)))
                end do
                got_it = .true.
              endif
            endif
          endif
        else not_EPA
          status = nf_inq_varid( ncid, 'TFLAG', varid )
          if( status /= nf_noerr )  call handle_error( status )
          status = nf_get_vara_int( ncid, varid, (/ 1,1,1 /), (/ 1,1,ntimes /), date )
          if( status /= nf_noerr )  call handle_error( status )
          status = nf_get_vara_int( ncid, varid, (/ 2,1,1 /), (/ 1,1,ntimes /), secs )
          if( status /= nf_noerr )  call handle_error( status )
          do m = 1,ntimes 
            yr = (date(m)/1000)*10000
            yr = yr + 101
            dy = mod( date(m),1000 ) - 1
            date(m) = newdate( yr, dy )

            sec = mod( secs(m),100 )
            hr = secs(m)/10000
            min = (secs(m)-10000*hr)/100
            secs(m) = sec + 60*(min + 60*hr)
          enddo
          got_it = .true.
        endif not_EPA
      endif

      if( .not. got_it ) then
        write(*,*) 'read_src_times: failed to read date,secs'
        stop
      endif

      if( data_yrs_offset /= 0 ) then
        do m = 1,ntimes
          date(m) = date(m) + 10000*data_yrs_offset
        end do
      endif

      end subroutine read_src_times

      subroutine get_src_time_ndx( data_file, loop_time, currDate, src_name, &
                                   src_file_prefix, src_file_suffix )
!---------------------------------------------------------------
!     get src time index for loop_time
!---------------------------------------------------------------

      use anthro_types, only : data_file_type, dates
      use mo_calendar,  only : diffdat
      use EPA,          only : translate_date

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in)                  :: currDate
      type(data_file_type), intent(inout)  :: data_file
      type(dates), intent(in)              :: loop_time
      character(len=*), intent(in)         :: src_name
      character(len=*), intent(in)         :: src_file_prefix
      character(len=*), intent(in)         :: src_file_suffix

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, n
      integer :: status
      integer :: ncid
      integer :: wrkDate
      character(len=128) :: filenm
      character(len=32)  :: wrkDateStr
      logical :: found

      write(*,*) ' '
      write(*,*) 'get_src_time_ndx; src_dir,src_fn = ',trim(data_file%filespec)
      write(*,*) 'get_src_time_ndx; interp_date,datesec,ntimes = ',loop_time%date,loop_time%secs,data_file%ntimes

      n = lotim( loop_time%date, loop_time%secs, data_file )
      if( n > 0 ) then
!---------------------------------------------------------------
!     loop time in present src dataset
!---------------------------------------------------------------
        write(*,*) 'get_src_time_ndx; tndx = ',n
        if( .not. data_file%in_gap ) then
          data_file%read_lo_tndx = n /= data_file%lo_tndx
          if( data_file%read_lo_tndx ) then
            data_file%lo_tndx = n
          endif
          data_file%in_gap  = .false.
          if( data_file%t_interp ) then
            data_file%read_hi_tndx = (n + 1) /= data_file%hi_tndx
            if( data_file%read_hi_tndx ) then
              data_file%hi_tndx = n + 1
            endif
          endif
        else
!---------------------------------------------------------------
!     special handling if prior time was in "gap"
!---------------------------------------------------------------
          data_file%in_gap  = .false.
          data_file%lo_tndx = n
          data_file%read_lo_tndx = .true.
          if( data_file%t_interp ) then
            data_file%hi_tndx = n + 1
            data_file%read_hi_tndx = .true.
          endif
        endif
        if( data_file%ncid_hi /= data_file%ncid_lo ) then
          status = nf_close( data_file%ncid_lo )
          if( status /= 0 ) then
            write(*,*) 'get_src_time_ndx: failed to close ',trim(data_file%filespec),' ; error = ',status
            stop
          end if
          data_file%ncid_lo = data_file%ncid_hi
        endif
        data_file%ncid_hi = data_file%ncid_lo
      else if( n < 0 ) then
!---------------------------------------------------------------
!     loop time after present src dataset
!---------------------------------------------------------------
         data_file%lo_tndx = data_file%ntimes
         if( .not. data_file%is_EPA ) then
           call next_flnm( data_file%filename, incr=.true. )
         else
           wrkDate = translate_date( currDate, trim(src_name) )
           write(wrkDateStr,'(i8)') wrkDate
           data_file%filename = trim(src_file_prefix) &
                                // trim(src_name) // '_' // trim(wrkDateStr) &
                                // trim(src_file_suffix)
         endif
!---------------------------------------------------------------
!     open the input netCDF file
!---------------------------------------------------------------
         if( .not. data_file%is_EPA ) then
           data_file%filespec = trim(anthro_dir) // '/' // trim(data_file%filename)
         else
           data_file%filespec = trim(anthro_dir) // '/' // trim(data_file%src_name) // '/' // trim(data_file%filename)
         endif
         status = nf_open( trim(data_file%filespec), nf_nowrite, ncid )
         if( status /= nf_noerr ) call handle_error( status )
         write(*,*) 'get_src_time_ndx: opened ',trim(data_file%filespec)
         data_file%gap_date = data_file%date(data_file%ntimes)
         data_file%gap_secs = data_file%secs(data_file%ntimes)
         call read_src_times( ncid, data_file%date, data_file%secs, data_file%ntimes, data_file%is_EPA )
         if( data_file%is_EPA ) then
           data_file%date(1:data_file%ntimes) = currDate
         endif
         n = lotim( loop_time%date, loop_time%secs, data_file )
         if( n > 0 ) then
           write(*,*) 'get_src_time_ndx; src_tndx = ',n
           status = nf_close( data_file%ncid_lo )
           if( status /= 0 ) then
             write(*,*) 'get_src_time_ndx: failed to close ',trim(data_file%filespec),' ; error = ',status
             stop
           end if
           data_file%in_gap  = .false.
           data_file%ncid_lo = ncid
           data_file%ncid_hi = ncid
           data_file%lo_tndx = n
           data_file%read_lo_tndx = .true.
           if( data_file%t_interp ) then
             data_file%hi_tndx = n + 1
             data_file%read_hi_tndx = .true.
           endif
         else if( n == 0 ) then
           data_file%in_gap  = .true.
           data_file%hi_tndx = 1
           data_file%ncid_hi = ncid
           data_file%dels = data_file%dels/diffdat( data_file%gap_date, data_file%gap_secs, data_file%date(1), data_file%secs(1) )
           data_file%t_interp = .true.
           data_file%read_lo_tndx = .true.
           data_file%read_hi_tndx = .true.
         else
           write(*,*) 'get_src_time_ndx: failed to find ',loop_time%date,' : ',loop_time%secs
           write(*,*) '                  in file ',trim(data_file%filename)
           status = nf_close( ncid )
           stop
         end if
      else if( data_file%in_gap ) then
        data_file%dels = data_file%dels/diffdat( data_file%gap_date, data_file%gap_secs, data_file%date(1), data_file%secs(1) )
        data_file%read_lo_tndx = .false.
        data_file%read_hi_tndx = .false.
      end if

      end subroutine get_src_time_ndx

      subroutine read_src_data( data_file, varname, bndx, tndx, ncid, cat_ndx )
!---------------------------------------------------------------
!     read source emission category
!---------------------------------------------------------------

      use anthro_types, only : data_file_type, stack_type
      use mapper_types, only : grid_type, grid_specs

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: bndx              ! buffer index
      integer, intent(in) :: tndx              ! time index
      integer, intent(in) :: ncid              ! netcdf file index
      integer, intent(in) :: cat_ndx           ! category buffer index
      type(data_file_type), intent(inout) :: data_file
      character(len=*), intent(in) :: varname

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: gndx
      integer :: n
      integer :: ndim
      integer :: status
      integer :: varid
      integer :: nstt(4), ncnt(4)

      gndx = data_file%grid_ndx
      if( .not. data_file%is_EPA_STACK ) then
        if( .not. data_file%is_EPA ) then
          nstt(1:3) = (/ 1, 1, tndx /)
          ncnt(1:3) = (/ grid_specs(gndx)%nlons, grid_specs(gndx)%nlats, 1 /)
          ndim = 3
        else
          nstt(:) = (/ 1, 1, 1, tndx /)
          ncnt(:) = (/ grid_specs(gndx)%nlons, grid_specs(gndx)%nlats, 1, 1 /)
          ndim = 4
        endif
        data_file%emis(:,:,bndx,cat_ndx) = 0.
      else
        nstt(:) = (/ 1, 1, 1, tndx /)
        ncnt(:) = (/ 1, data_file%stack%nStk, 1, 1 /)
        ndim = 4
      endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     read data, replace missing data with zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      status = nf_inq_varid( ncid, trim(varname), varid )
      if( status /= nf_noerr )  then
        write(*,*) 'read_src_data: failed to get id of ',trim(varname)
        call handle_error( status )
      end if

      if( .not. data_file%is_EPA_STACK ) then
        status = nf_get_vara_real( ncid, varid, nstt(1:ndim), ncnt(1:ndim), data_file%emis(1,1,bndx,cat_ndx) )
        if( status /= nf_noerr ) then
          write(*,*) 'read_src_data: failed to read ',trim(varname)
          call handle_error( status )
        end if
        where( data_file%emis(:,:,bndx,cat_ndx) == data_file%missing_value )
          data_file%emis(:,:,bndx,cat_ndx) = 0.
        endwhere
      else
        status = nf_get_vara_real( ncid, varid, nstt(1:ndim), ncnt(1:ndim), data_file%stack%emis(1,bndx,cat_ndx) )
        if( status /= nf_noerr ) then
          write(*,*) 'read_src_data: failed to read ',trim(varname)
          call handle_error( status )
        end if
        where( data_file%stack%emis(:,bndx,cat_ndx) == data_file%missing_value )
          data_file%stack%emis(:,bndx,cat_ndx) = 0.
        endwhere
      end if

      end subroutine read_src_data

      subroutine tinterp_src_data( data_file, cat_ndx )
!---------------------------------------------------------------
!     time interpolation of source emission
!---------------------------------------------------------------

      use anthro_types, only : data_file_type, stack_type
      use mapper_types, only : grid_type, grid_specs

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: cat_ndx
      type(data_file_type), intent(inout) :: data_file

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, j, iu, ju
      integer :: lo_ndx, hi_ndx
      integer :: astat
      integer :: gndx
      real    :: d1
      real    :: dels                ! linear interp factor
      real    :: delsm1              ! linear interp factor
      real, allocatable  :: data_col(:)

      gndx = data_file%grid_ndx
      dels = data_file%dels
      lo_ndx = data_file%lo_buf_ndx
      if( data_file%t_interp ) then
        hi_ndx = data_file%hi_buf_ndx
        delsm1 = 1. - dels
        if( .not. data_file%is_EPA_STACK ) then
          do j = 1,grid_specs(gndx)%nlats
            data_file%src_data(:,j) = data_file%emis(:,j,lo_ndx,cat_ndx) * dels &
                                    + data_file%emis(:,j,hi_ndx,cat_ndx) * delsm1
          end do
        else
          data_file%stack%src_data(:) = data_file%stack%emis(:,lo_ndx,cat_ndx) * dels &
                                      + data_file%stack%emis(:,hi_ndx,cat_ndx) * delsm1
        endif
      else
        if( .not. data_file%is_EPA_STACK ) then
          do j = 1,grid_specs(gndx)%nlats
            data_file%src_data(:,j) = data_file%emis(:,j,lo_ndx,cat_ndx)
          end do
        else
          data_file%stack%src_data(:) = data_file%stack%emis(:,lo_ndx,cat_ndx)
        endif
      endif

not_STACK: &
      if( .not. data_file%is_EPA_STACK ) then
      if( grid_specs(gndx)%reorder_lons ) then
!-------------------------------------------------------------
!  longitude reorder the src data
!-------------------------------------------------------------
        do j = 1,grid_specs(gndx)%nlats
          do i = 1,grid_specs(gndx)%nlons/2
            d1 = data_file%src_data(i,j)
            iu = grid_specs(gndx)%nlons - i + 1
            data_file%src_data(i,j)  = data_file%src_data(iu,j)
            data_file%src_data(iu,j) = d1
          end do
        end do
      endif

      if( grid_specs(gndx)%has_lon_shift ) then
!-------------------------------------------------------------
!  longitude shift the src data
!-------------------------------------------------------------
        do j = 1,grid_specs(gndx)%nlats
          data_file%src_data(:,j) = cshift( data_file%src_data(:,j),grid_specs(gndx)%nlons/2 )
        end do
      endif

      if( grid_specs(gndx)%reorder_lats ) then
!-------------------------------------------------------------
!  latitude reorder the src data
!-------------------------------------------------------------
        allocate( data_col(grid_specs(gndx)%nlons),stat=astat )
        if( astat /= 0 ) then
          write(*,*) 'tinerpt_src_data: failed to allocate data_col: error = ',astat
          stop 'Alloc error'
        endif
        do j = 1,grid_specs(gndx)%nlats/2
          ju = grid_specs(gndx)%nlats - j + 1
          data_col(:) = data_file%src_data(:,j)
          data_file%src_data(:,j)  = data_file%src_data(:,ju)
          data_file%src_data(:,ju) = data_col(:)
        end do
        deallocate( data_col )
      endif
      endif not_STACK

      end subroutine tinterp_src_data

      subroutine next_flnm( filenm, incr )
!---------------------------------------------------------------
!     increment mozart filename
!---------------------------------------------------------------

      use mo_calendar, only : newdate

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      character(len=*), intent(inout) :: filenm
      logical, intent(in)             :: incr

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: il, iu, ios, nd
      integer :: file_number
      integer :: date
      logical :: found
      character(len=6) :: frmt

      found = .false.
not_EPA: &
      if( .not. is_EPA ) then
        iu = scan( trim(filenm), '0123456789', back=.true. )
        do il = iu-1,1,-1
          if( index( '0123456789', filenm(il:il) ) == 0 ) then
            found = .true.
            exit
         endif
        end do
        if( found ) il = il + 1
      else not_EPA
        call get_flnm_date( trim(filenm), il, iu )
        found = iu - il == 7
      endif not_EPA

      if( .not. found ) then
        write(*,*) 'next_filenm: ',trim(filenm)
        write(*,*) 'next_filenm: filename does not contain date string'
        write(*,*) '             or sequence number'
        stop
      endif

      nd = iu - il + 1
      if( nd > 9 ) then
        if( .not. is_EPA ) then
          nd = 9
          il = iu - 8
        else
          write(*,*) 'next_filenm: filename does not contain <yyyymmdd> date string'
          stop
        endif
      endif

      write(*,*) ' '
      write(*,*) 'next_flnm; trying to increment file ',trim(filenm)
      write(*,*) 'next_flnm; il, iu, nd = ',il,iu,nd

not_EPA_a: &
      if( .not. is_EPA ) then
        read(filenm(il:iu),*,iostat=ios) file_number
        if( ios /= 0 ) then
           write(*,*) 'next_filenm: failed to read ',filenm(il:iu),' ; error = ',ios
           stop
        endif
      write(*,*) 'next_flnm; file_number = ',file_number
      if( incr ) then
        file_number = file_number + 1
      else
        file_number = file_number - 1
      endif

      write(frmt,'(''(i'',i1.1,''.'',i1.1,'')'')') nd,nd
      if( ios /= 0 ) then
        write(*,*) 'next_filenm: failed to write frmt; error = ',ios
        stop
      endif
      write(filenm(il:iu),frmt) file_number
      if( ios /= 0 ) then
        write(*,*) 'next_filenm: failed to write new filename; error = ',ios
        stop
      endif
      else not_EPA_a
        read(filenm(il:iu),*,iostat=ios) date
        if( ios /= 0 ) then
          stop 'Rd err'
        endif
        date = newdate( date, 1 )
        write(filenm(il:iu),'(i8)',iostat=ios) date
        if( ios /= 0 ) then
          stop 'Wr err'
        endif
      endif not_EPA_a

      write(*,*) 'next_flnm; new file = ',trim(filenm)

      end subroutine next_flnm

      subroutine get_flnm_date( flnm, start, end, date )
!-----------------------------------------------------------------------
! 	... return the date string in filename, if found
!-----------------------------------------------------------------------

      use mo_calendar, only : newdate

      character(len=*), intent(in)    :: flnm
      integer, optional, intent(out)  :: start, end
      integer, optional, intent(out)  :: date

      integer :: il, iu
      integer :: ios, pos, slen
      logical :: found

      slen = len_trim(flnm)
      il   = 1
char_loop: &
      do
        found = .false.
        do pos = il,slen
          if( index( '0123456789', flnm(pos:pos) ) /= 0 ) then
            il = pos
            found = .true.
            exit
          endif
        end do
        if( .not. found ) then
          write(*,*) 'get_flnm_date: ',trim(flnm)
          stop 'get_flnm_date: could not find date'
        endif
        do iu = il+1,len_trim(flnm)
          if( index( '0123456789', flnm(iu:iu) ) == 0 ) then
            found = .true.
            exit
          endif
        end do
        if( .not. found ) then
          write(*,*) 'get_flnm_date: ',trim(flnm)
          stop 'get_flnm_date: could not find date'
        endif
        if( iu - il /= 8 ) then
          il = iu
          cycle char_loop
        endif
        iu = iu - 1
        if( present(date) ) then
          read(flnm(il:iu),*,iostat=ios) date
          if( ios /= 0 ) then
            write(*,*) 'get_flnm_date: date string = ',trim(flnm(il:iu))
            stop 'get_flnm_date: Rd err'
          endif
          if( present(start) .and. present(end) ) then
            start = il ; end = iu
          endif
          exit char_loop
        elseif( present(start) .and. present(end) ) then
          start = il ; end = iu
          exit char_loop
        endif
      end do char_loop

      end subroutine get_flnm_date

      integer function lotim( cdate, csec, data_file )
!-----------------------------------------------------------------------
! 	... return the index of the time sample that is the lower
!           bound of the interval that contains the input date.  if
!           (cdate,csec) is earlier than the first time sample then 0 is
!           returned.  if (cdate,csec) is later than the last time sample then
!           -index is returned.  if (cdate,csec) is equal to the date of a
!           dynamics time sample then that index is returned.
!-----------------------------------------------------------------------

      use mo_calendar,  only : diffdat
      use anthro_types, only : data_file_type

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: cdate    ! date in yyyymmdd
      integer, intent(in) :: csec     ! seconds relative to date
      type(data_file_type), intent(inout) :: data_file

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer :: n
      real    :: dtime

!-----------------------------------------------------------------------
!     	... find latest date that is earlier than or equal to (date,sec)
!-----------------------------------------------------------------------
      do n = 1,data_file%ntimes
        dtime = diffdat( cdate, csec, data_file%date(n), data_file%secs(n) )
        if( dtime > 0. ) then
          lotim = n - 1
          if( lotim > 0 ) then
            data_file%dels = dtime/diffdat( data_file%date(lotim), data_file%secs(lotim), data_file%date(n), data_file%secs(n) )
            data_file%t_interp = data_file%dels /= 1.
          else
            data_file%dels = dtime
          endif
          exit
        endif
      end do

      if( n > data_file%ntimes ) then
        if( dtime == 0. ) then
          lotim = data_file%ntimes
          data_file%dels = 1.
          data_file%t_interp = .false.
        else
          lotim = -data_file%ntimes
          data_file%dels = dtime
        endif
      endif

      end function lotim

   subroutine get_molecwght( ncid, cat_var_prefix, cat_var_suffix, molecw, missing_value )
!---------------------------------------------------------------------
!	... try to get the source species molecular wght from src file
!---------------------------------------------------------------------

   use utils, only : n_sub_cats, sub_cats

!---------------------------------------------------------------------
!	... dummy args
!---------------------------------------------------------------------
   integer, intent(in) :: ncid
   real, intent(in)    :: missing_value
   real, intent(inout) :: molecw
   character(len=*), intent(in) :: cat_var_prefix
   character(len=*), intent(in) :: cat_var_suffix

!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
   integer :: m
   integer :: ncstat
   integer :: varid
   integer :: molecw_int
   integer(2) :: molecw_int2
   integer :: att_type, att_len
   real(8) :: molecw_8
   character(len=64)  :: var_nam
   character(len=132) :: message
   logical :: got_it

   got_it = .false.
!---------------------------------------------------------------------
!	... first try to get the variable molecular_weight
!---------------------------------------------------------------------
   ncstat = nf_inq_varid( ncid, 'molecular_weight', varid )
   if( ncstat == nf_noerr ) then
     ncstat = nf_get_var_real( ncid, varid, molecw )
     if( ncstat /= nf_noerr ) then
       write(*,*) 'get_molecwght: failed to read molecular weight: error = ',ncstat
       stop 'Netcdf error'
     endif
     got_it = molecw /= missing_value .and. molecw > 0. .and. molecw < 300.
   endif

   if( .not. got_it ) then
!---------------------------------------------------------------------
!	... next try to get global attribute molecular_weight
!---------------------------------------------------------------------
     ncstat = nf_inq_att( ncid, nf_global, 'molecular_weight', att_type, att_len )
     if( ncstat == nf_noerr ) then
       if( att_len == 1 ) then
         select case( att_type )
           case( nf_short )
             ncstat = nf_get_att_int2( ncid, nf_global, 'molecular_weight', molecw_int2 )
           case( nf_int )
             ncstat = nf_get_att_int( ncid, nf_global, 'molecular_weight', molecw_int )
           case( nf_float )
             ncstat = nf_get_att_real( ncid, nf_global, 'molecular_weight', molecw )
           case( nf_double )
             ncstat = nf_get_att_double( ncid, nf_global, 'molecular_weight', molecw_8 )
           case default
             ncstat = nf_noerr + 1
         end select
         if( ncstat /= nf_noerr ) then
           write(*,*) 'get_molecwght: failed to read molecular weight: error = ',ncstat
           stop 'Netcdf error'
         endif
         select case( att_type )
           case( nf_short )
             molecw = real(molecw_int2)
           case( nf_int )
             molecw = real(molecw_int)
           case( nf_double )
             molecw = real(molecw_8)
         end select
         got_it = .true.
       endif
     endif
   endif

   if( .not. got_it ) then
!---------------------------------------------------------------------
!	... finally try category attribute
!---------------------------------------------------------------------
sub_cat_loop : &
     do m = 1,n_sub_cats
       var_nam = trim( cat_var_prefix) // trim( sub_cats(m) ) // trim( cat_var_suffix )
       ncstat = nf_inq_varid( ncid, trim(var_nam), varid )
       if( ncstat == nf_noerr ) then
         ncstat = nf_inq_att( ncid, varid, 'molecular_weight', att_type, att_len )
         if( ncstat == nf_noerr ) then
           if( att_len == 1 ) then
             select case( att_type )
               case( nf_short )
                 ncstat = nf_get_att_int2( ncid, varid, 'molecular_weight', molecw_int2 )
               case( nf_int )
                 ncstat = nf_get_att_int( ncid, varid, 'molecular_weight', molecw_int )
               case( nf_float )
                 ncstat = nf_get_att_real( ncid, varid, 'molecular_weight', molecw )
               case( nf_double )
                 ncstat = nf_get_att_double( ncid, varid, 'molecular_weight', molecw_8 )
               case default
                 ncstat = nf_noerr + 1
             end select
             if( ncstat == nf_noerr ) then
               select case( att_type )
                 case( nf_short )
                   molecw = real(molecw_int2)
                 case( nf_int )
                   molecw = real(molecw_int)
                 case( nf_double )
                   molecw = real(molecw_8)
               end select
               got_it = .true.
             endif
           endif
         endif
       endif
       if( got_it ) then
         exit sub_cat_loop
       endif
     end do sub_cat_loop
   endif

   got_it = molecw /= missing_value .and. molecw > 0. .and. molecw < 300.
   if( .not. got_it ) then
     write(*,*) 'get_molecwght: failed to read molecular weight'
     call handle_error( ncstat )
     stop
   endif

   end subroutine get_molecwght

   subroutine get_units( ncid, cat_var_prefix, cat_var_suffix, molecw, &
                         con_fac, dx, mdlGrid )
!---------------------------------------------------------------------
!	... try to get the source species molecular wght from src file
!---------------------------------------------------------------------

   use anthro_types, only : model_grid_type
   use utils, only : n_sub_cats, sub_cats
   use utils, only : upcase

!---------------------------------------------------------------------
!	... dummy args
!---------------------------------------------------------------------
   integer, intent(in) :: ncid
   real, intent(in)    :: molecw
   real, intent(in)    :: dx                             ! grid edge in km
   real, intent(inout) :: con_fac(2)
   character(len=*), intent(in) :: cat_var_prefix
   character(len=*), intent(in) :: cat_var_suffix
   type(model_grid_type), intent(in) :: mdlGrid

!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
   real, parameter :: avogadro = 6.02214129e23
   real, parameter :: secphr   = 3600.

   integer :: m
   integer :: ncstat
   integer :: varid
   integer :: att_type, att_len
   character(len=64)  :: var_nam
   character(len=132) :: units, wrk_string

sub_cat_loop : &
   do m = 1,n_sub_cats
     var_nam = trim( cat_var_prefix) // trim( sub_cats(m) ) // trim( cat_var_suffix )
     ncstat = nf_inq_varid( ncid, trim(var_nam), varid )
     if( ncstat == nf_noerr ) then
       ncstat = nf_inq_att( ncid, varid, 'units', att_type, att_len )
       if( ncstat == nf_noerr ) then
         if( att_type == nf_char .and. att_len /= 0 ) then
           units = ' '
           ncstat = nf_get_att_text( ncid, varid, 'units', units )
           if( ncstat == nf_noerr ) then
             call upcase( units(:len_trim(units)), wrk_string )
             if( mdlGrid%mdl_is_WRF ) then
               select case( trim(wrk_string) )
                 case( 'KG M-2 S-1','KG/M2/S','KG/M^2/S' )
                   con_fac(1) = 3.6e12 ; con_fac(2) = 1.e9
                 case( 'MOLECULES/CM2/S','MOLECULES/CM^2/S','MOLECULES CM-2 S-1' )
                   con_fac(1) = 3.6e13*molecw/avogadro
                   con_fac(2) = 1.e10*molecw/avogadro
                 case( 'MOLES/S' )
                   con_fac(2) = 1./(dx*dx) ; con_fac(1) = secphr*con_fac(2)
                   con_fac(2) = molecw * con_fac(2)
                 case( 'G/S' )
                   con_fac(2) = 1./(dx*dx)
                 case default
                   con_fac(1) = 3.6e12 ; con_fac(2) = 1.e9
               end select
             elseif( mdlGrid%mdl_is_CAM ) then
!---------------------------------------------------------------
!     CAM units are molecules/cm^2/s for gas and aerosol
!---------------------------------------------------------------
               select case( trim(wrk_string) )
                 case( 'MOLECULES/CM2/S','MOLECULES/CM^2/S','MOLECULES CM-2 S-1' )
                   con_fac(1:2) = 1.
                 case( 'MOLES/S', 'G/S' )
                   con_fac(1:2) = avogadro
                 case default
                   con_fac(1:2) = 1.
               end select
             endif
             exit sub_cat_loop
           endif
         endif
       endif
     endif
   enddo sub_cat_loop

   write(*,*) '===================================================='
   write(*,*) 'get_units: con_fac(1,2) = ',con_fac(:)
   write(*,*) '===================================================='

   end subroutine get_units

   subroutine handle_error( status )
!---------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments :
!---------------------------------------------------------------
   integer, intent(in) :: status

!---------------------------------------------------------------
!     print the error information from processing NETcdf file
!---------------------------------------------------------------
   write(*,*) nf_strerror( status )
   stop 'Netcdf error'

   end subroutine handle_error

   end module data_file_utils
