
   module camse_utils

   implicit none

   private
   public :: camse_init, camse_mapper

   real, parameter :: ZERO    = 0.
   real, parameter :: ZERO_8  = 0._8
   real, parameter :: ONEHALF = .5
   real, parameter :: ONEHALF_8 = .5_8
   real, parameter :: ONE     = 1.
   real, parameter :: ONE_8   = 1._8
   real, parameter :: FOUR    = 4.
   real, parameter :: FOUR_8  = 4._8
   real, parameter :: TEN     = 10.
   real, parameter :: TEN_8   = 10_8
   real, parameter :: NINETY  = 90.
   real, parameter :: ONE80   = 180.
   real, parameter :: ONE80_8 = 180._8
   real, parameter :: THREE60 = 360.
   real, parameter :: THREE60_8 = 360._8
   real, parameter :: ROUNDOFF = 10.*epsilon(ROUNDOFF)

   integer :: maxVtx
   real(8) :: PI, D2R, R2D

   include 'netcdf.inc'

   CONTAINS

   subroutine camse_init( ncid, mdl_grd )

   use fire_types, only : mdl_poly_type, model_grid_type
   use netcdf_utils

!-----------------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------------
   integer, intent(in) :: ncid         ! model netcdf file id
   type(model_grid_type), intent(inout) :: mdl_grd

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   integer :: polyNdx, vtxNdx
   integer :: astat
   integer :: dimid, varid
   integer :: nMatchLon, nMatchLat, nMatchVtx
   integer :: m, mp1
   integer :: nVtx, nPolygons
   integer :: CapNdx
   integer :: minmaxLat(2)
   integer :: Quadcnt(4)
   integer, allocatable :: Quad(:)
   real    :: MatchVtxLon, MatchVtxLat
   real, allocatable :: wrk(:,:)
   real, allocatable :: wrkLon(:), wrkLat(:), xprod(:)
   character(len=64) :: mess

!-----------------------------------------------------------------------
!  get number polygons
!-----------------------------------------------------------------------
   mess = 'camse_init: Failed to get grid_size dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'grid_size', dimid ), mess )
   mess = 'camse_init: Failed to get grid_size dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, mdl_grd%nPolygons ), mess )
   nPolygons = mdl_grd%nPolygons

!-----------------------------------------------------------------------
!  get max number vertices/polygon
!-----------------------------------------------------------------------
   mess = 'camse_init: Failed to get grid_corners dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'grid_corners', dimid ), mess )
   mess = 'camse_init: Failed to get grid_corners dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, mdl_grd%maxPolyvtx ), mess )
   maxVtx = mdl_grd%maxPolyvtx

!-----------------------------------------------------------------------
!  allocate mdl_poly_type
!-----------------------------------------------------------------------
   allocate( mdl_grd%mdl_poly(nPolygons),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'camse_init: Failed to allocate model polygons;  error = ',astat
     stop 'AllocErr'
   endif
!-----------------------------------------------------------------------
!  allocate wrk
!-----------------------------------------------------------------------
   allocate( wrk(nPolygons,1),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'camse_init: Failed to allocate wrk;  error = ',astat
     stop 'AllocErr'
   endif
!-----------------------------------------------------------------------
!  get and distribute polygon center longitudes (degrees)
!-----------------------------------------------------------------------
   mess = 'camse_init: Failed to get grid_center_lon variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'grid_center_lon', varid ), mess )
   mess = 'camse_init: Failed to get grid_center_lon'
   call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), mess )
   do polyNdx = 1,nPolygons
     mdl_grd%mdl_poly(polyNdx)%cntr_lon = wrk(polyNdx,1)
   enddo
!-----------------------------------------------------------------------
!  get and distribute polygon center latitudes (degrees)
!-----------------------------------------------------------------------
   mess = 'camse_init: Failed to get grid_center_lat variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'grid_center_lat', varid ), mess )
   mess = 'camse_init: Failed to get grid_center_lon'
   call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), mess )
   do polyNdx = 1,nPolygons
     mdl_grd%mdl_poly(polyNdx)%cntr_lat = wrk(polyNdx,1)
   enddo
!-----------------------------------------------------------------------
!  get and distribute polygon areas (square radians)
!-----------------------------------------------------------------------
   mess = 'camse_init: Failed to get grid_center_area variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'grid_area', varid ), mess )
   mess = 'camse_init: Failed to get grid_area'
   call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), mess )
   do polyNdx = 1,nPolygons
     mdl_grd%mdl_poly(polyNdx)%area = wrk(polyNdx,1)
   enddo

   write(*,*) ' '
   write(*,'(''camse_init: mdl surf area = '',1pg15.7)') sum(wrk(:,1))

!-----------------------------------------------------------------------
!  reallocate wrk
!-----------------------------------------------------------------------
   deallocate( wrk )
   allocate( wrk( maxVtx,nPolygons),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'camse_init: Failed to allocate wrk;  error = ',astat
     stop 'AllocErr'
   endif
   do polyNdx = 1,nPolygons
     allocate( mdl_grd%mdl_poly(polyNdx)%vtx_lon(maxVtx), &
               mdl_grd%mdl_poly(polyNdx)%vtx_lat(maxVtx),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'camse_init: Failed to allocate polygon vtx lon,lat;  error = ',astat
       stop 'AllocErr'
     endif
   enddo
!-----------------------------------------------------------------------
!  get and distribute polygon vertex longitudes (degrees)
!-----------------------------------------------------------------------
   mess = 'camse_init: Failed to get grid_vertex_lon variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'grid_corner_lon', varid ), mess )
   mess = 'camse_init: Failed to get grid_vertex_lon'
   call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), mess )
   do polyNdx = 1,nPolygons
     mdl_grd%mdl_poly(polyNdx)%vtx_lon(:) = wrk(:,polyNdx)
   enddo
!-----------------------------------------------------------------------
!  get and distribute polygon vertex latitudes (degrees)
!-----------------------------------------------------------------------
   mess = 'camse_init: Failed to get grid_vertex_lat variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'grid_corner_lat', varid ), mess )
   mess = 'camse_init: Failed to get grid_vertex_lat'
   call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), mess )
   do polyNdx = 1,nPolygons
     mdl_grd%mdl_poly(polyNdx)%vtx_lat(:) = wrk(:,polyNdx)
   enddo

   deallocate( wrk )

!-----------------------------------------------------------------------
!  delineate the "unique" model grid polygon vertices
!-----------------------------------------------------------------------
   do polyNdx = 1,nPolygons
     matchVtxLon = mdl_grd%mdl_poly(polyNdx)%vtx_lon(maxVtx)
     matchVtxLat = mdl_grd%mdl_poly(polyNdx)%vtx_lat(maxVtx)
     nMatchLon   = count( mdl_grd%mdl_poly(polyNdx)%vtx_lon(:) == matchVtxLon )
     nMatchLat   = count( mdl_grd%mdl_poly(polyNdx)%vtx_lat(:) == matchVtxLat )
     nMatchVtx   = min( nMatchLon,nMatchLat )
     mdl_grd%mdl_poly(polyNdx)%nVtx = maxVtx - nMatchVtx + 1
     nVtx        = mdl_grd%mdl_poly(polyNdx)%nVtx
     if( nVtx < 3 ) then
       write(*,'(''camse_utils: Polygon '',i6,'' has fewer than 3 vertices'')') polyNdx
       Stop 'Poly-ERR'
     endif
   enddo

   allocate( Quad(maxVtx),wrkLon(maxVtx),wrkLat(maxVtx),xprod(maxVtx),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'camse_init: Failed to allocate Quad,wrkLon;  error = ',astat
     stop 'AllocErr'
   endif

   CapNdx = 0
   mdl_grd%mdl_poly(1:nPolygons)%active = .true.
   mdl_grd%mdl_poly(1:nPolygons)%x180   = .false.
!-----------------------------------------------------------------------
!  mark polygons with edges that "cross" 180 longitude
!-----------------------------------------------------------------------
poly_loop: &
   do polyNdx = 1,nPolygons
     nVtx = mdl_grd%mdl_poly(polyNdx)%nVtx
     wrkLon(1:nVtx) = mdl_grd%mdl_poly(polyNdx)%vtx_lon(1:nVtx) - mdl_grd%mdl_poly(polyNdx)%cntr_lon
     wrkLat(1:nVtx) = mdl_grd%mdl_poly(polyNdx)%vtx_lat(1:nVtx) - mdl_grd%mdl_poly(polyNdx)%cntr_lat
     do m = 1,nVtx
       if( m /= nVtx ) then
         mp1 = m + 1
       else
         mp1 = 1
       endif
       xprod(m) = wrkLon(m)*wrkLat(mp1) - wrkLon(mp1)*wrkLat(m)
     enddo
     if( any(xprod(1:nVtx) < ZERO) ) then
       write(*,'(''camse_init: negative cross prod for polygon,neg,pos '',i6,2i3)') &
         polyndx,count(xprod(1:nVtx) < ZERO),count(xprod(1:nVtx) > ZERO)
       write(*,'(''camse_init: poly lons = '',1p10g15.7)') mdl_grd%mdl_poly(polyndx)%vtx_lon(1:nVtx)
     endif
!-----------------------------------------------------------------------
!  put all mdl polygon vertices and centers into range (-180,180)
!-----------------------------------------------------------------------
     if( mdl_grd%mdl_poly(polyNdx)%cntr_lon > ONE80 ) then
       mdl_grd%mdl_poly(polyNdx)%cntr_lon = &
          mdl_grd%mdl_poly(polyNdx)%cntr_lon - THREE60
     endif
     where( mdl_grd%mdl_poly(polyNdx)%vtx_lon(1:nVtx) > ONE80 )
       mdl_grd%mdl_poly(polyNdx)%vtx_lon(1:nVtx) = &
          mdl_grd%mdl_poly(polyNdx)%vtx_lon(1:nVtx) - THREE60
     endwhere
     wrkLon(1:nVtx) = mdl_grd%mdl_poly(polyNdx)%vtx_lon(1:nVtx)
!-----------------------------------------------------------------------
!  place vtx in "quadrant" by longitude
!-----------------------------------------------------------------------
     Quadcnt(:) = 0
vtx_loop: &
     do vtxNdx = 1,nVtx
       if( ZERO <= wrkLon(vtxNdx) .and. wrkLon(vtxNdx) < NINETY ) then
         Quad(vtxNdx) = 1
         Quadcnt(1) = Quadcnt(1) + 1
       elseif( NINETY <= wrkLon(vtxNdx) .and. wrkLon(vtxNdx) < ONE80 ) then
         Quad(vtxNdx) = 2
         Quadcnt(2) = Quadcnt(2) + 1
       elseif( -ONE80 <= wrkLon(vtxNdx) .and. wrkLon(vtxNdx) < -NINETY ) then
         Quad(vtxNdx) = 3
         Quadcnt(3) = Quadcnt(3) + 1
       elseif( -NINETY <= wrkLon(vtxNdx) .and. wrkLon(vtxNdx) < ZERO ) then
         Quad(vtxNdx) = 4
         Quadcnt(4) = Quadcnt(4) + 1
       endif
     enddo vtx_loop
     if( count( Quadcnt(:) /= 0 ) > 2 ) then
       mdl_grd%mdl_poly(polyNdx)%active = .false.
       CapNdx = CapNdx + 1
       if( CapNdx < 3 ) then
         mdl_grd%PolarCapNdx(CapNdx) = polyNdx
       else
         write(*,'(''camse_init: more than two polar cap model grid cells; terminating'')')
         stop 'Runtime-Err'
       endif
     elseif( minval(wrkLon(1:nvtx)) * maxval(wrkLon(1:nVtx)) < 0. ) then
       if( any(Quad(:nVtx) == 2) .and. any(Quad(:nVtx) == 3) ) then
         mdl_grd%mdl_poly(polyNdx)%x180 = .true.
       endif
     endif
   enddo poly_loop

   deallocate( Quad,wrkLon,wrkLat,xprod )

!-----------------------------------------------------------------------
!  diagnostics
!-----------------------------------------------------------------------
   write(*,*) ' '
   write(*,'(''camse_init: There are '',i6,'' total mdl grid polygons'')') nPolygons
   write(*,'(''camse_init: There are '',i2,'' polar cap polygons'')') &
     count( .not. mdl_grd%mdl_poly(:)%active )
   write(*,'(''camse_init: There are '',i6,'' total mdl polygons crossing 180 longitude'')') &
     count( mdl_grd%mdl_poly(:)%x180 )
   minmaxLat(1:1) = minloc( abs(mdl_grd%mdl_poly(:)%cntr_lat),mask=.not. mdl_grd%mdl_poly(:)%x180 )
   minmaxLat(2:2) = maxloc( abs(mdl_grd%mdl_poly(:)%cntr_lat),mask=.not. mdl_grd%mdl_poly(:)%x180 )
   write(*,*) ' '
   write(*,'(''camse_init: masked polygon ('',i6,'') cntr_lat min = '',1pg15.7)') &
     minmaxLat(1),mdl_grd%mdl_poly(minmaxLat(1))%cntr_lat
   write(*,*) ' '
   write(*,'(''camse_init: masked polygon ('',i6,'') cntr_lat max = '',1pg15.7)') &
     minmaxLat(2),mdl_grd%mdl_poly(minmaxLat(2))%cntr_lat

   PI = FOUR*atan(ONE)
   D2R = PI/ONE80_8
   R2D = ONE80_8/PI

   end subroutine camse_init

   subroutine camse_mapper( longitudes, lats, mdlGrd, lon_ndx )
!-----------------------------------------------------------------------
!  map fire lon,lat to camse grid cell
!-----------------------------------------------------------------------

   use fire_types, only : mdl_poly_type, model_grid_type

!-----------------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------------
   real, intent(in)  :: longitudes(:)
   real, intent(in)  :: lats(:)
   integer, intent(out) :: lon_ndx(:)
   type(model_grid_type), intent(inout) :: mdlGrd

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   integer :: inside
   integer :: polyNdx, nVtx
   integer :: fireNdx, nFires
   integer :: nLost, nNeg
   real    :: lons(size(longitudes))
   real    :: fireLon, fireLat
   real    :: wrkfireLon
   real    :: wrkmdlLons(mdlGrd%maxPolyVtx)
   real(8) :: TransLon
   real(8) :: pnt(2)
   real(8) :: mPoly_x(mdlGrd%maxPolyVtx), mPoly_y(mdlGrd%maxPolyVtx)
   logical :: found

   type enclosingBox
     real :: minLon, maxLon
     real :: minLat, maxLat
   end type enclosingBox
   type(enclosingBox) :: mdlBox

   integer :: pnt_in_poly

   nFires = size( lons )
  
   lons(:)    = longitudes(:)
   lon_ndx(:) = 0
   where( lons(:) > ONE80 )
     lons(:) = lons(:) - THREE60
   endwhere

   write(*,*) ' '
   write(*,*) 'camse_mapper: diagnostics'
   write(*,'(''camse_mapper: (min,max) longitude = ('',1pg15.7,'','',g15.7,'')'')') minval(lons(:)),maxval(lons(:))
   write(*,'(''camse_mapper: (min,max) latitude  = ('',1pg15.7,'','',g15.7,'')'')') minval(lats(:)),maxval(lats(:))

   nLost = 0 ; nNeg = 0
fire_loop: &
   do fireNdx = 1,nFires
     fireLon = lons(fireNdx) ; fireLat = lats(fireNdx)
     found = .false.
     if( fireLon < ZERO ) then
       nNeg = nNeg + 1
     endif
poly_loop: &
     do polyNdx = 1,mdlGrd%nPolygons
poly_mask: &
       if( mdlGrd%mdl_poly(polyNdx)%active ) then
         nVtx = mdlGrd%mdl_poly(polyNdx)%nVtx
         wrkmdlLons(1:nVtx) = mdlGrd%mdl_poly(polyNdx)%vtx_lon(1:nVtx)
         wrkfireLon = fireLon
         if( mdlGrd%mdl_poly(polyNdx)%x180 ) then
           where( wrkmdlLons(1:nVtx) < ZERO )
             wrkmdlLons(1:nVtx) = wrkmdlLons(1:nVtx) + THREE60
           endwhere 
           if( wrkfireLon < ZERO ) then
             wrkfireLon = fireLon + THREE60
           endif
         endif
         mdlBox%minLon = minval(wrkmdlLons(1:nVtx))
         mdlBox%maxLon = maxval(wrkmdlLons(1:nVtx))
         mdlBox%minLat = minval(mdlGrd%mdl_poly(polyNdx)%vtx_lat(1:nVtx))
         mdlBox%maxLat = maxval(mdlGrd%mdl_poly(polyNdx)%vtx_lat(1:nVtx))
in_mdl_box: &
         if( pnt_n_Box() ) then
           pnt(:) = (/ real(wrkfireLon,8),real(fireLat,8) /)
           mPoly_x(1:nVtx) = real( wrkmdlLons(1:nVtx),8)
           mPoly_y(1:nVtx) = real( mdlGrd%mdl_poly(polyNdx)%vtx_lat(1:nVtx),8)
           inside = pnt_in_poly( nVtx, pnt, mPoly_x, mPoly_y )
           if( inside == 1 ) then
             found = .true.
             lon_ndx(fireNdx) = polyNdx
             exit poly_loop
           endif
         endif in_mdl_box
       endif poly_mask
     enddo poly_loop
     if( .not. found ) then
       nLost = nLost + 1
       write(*,'(''camse_mapper: fire @ (lon,lat) = ('',1pg15.7,'','',g15.7,'') not found; firendx = '',i6)') fireLon,fireLat,fireNdx
     endif
   enddo fire_loop

   write(*,*) ' '
   write(*,'(''camse_mapper: '',i6,''/'',i6,'' fires with lon < 0 '')') nNeg,nFires
   write(*,'(''camse_mapper: '',i6,''/'',i6,'' fires not mapped'')') nLost,nFires

   CONTAINS

   function pnt_n_Box() result(inBox)

   logical :: inBox

   inBox = mdlBox%minLon <= wrkfireLon .and. wrkfireLon <= mdlBox%maxLon
   if( inBox ) then
     inBox = mdlBox%minLat <= fireLat .and. fireLat <= mdlBox%maxLat
   endif

   end function pnt_n_Box


   end subroutine camse_mapper

   end module camse_utils
