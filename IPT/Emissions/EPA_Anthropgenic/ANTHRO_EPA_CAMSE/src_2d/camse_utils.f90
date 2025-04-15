
   module camse_utils

   implicit none

   private
   public :: camse_init, area_interp_init_camse

   real, parameter :: ZERO    = 0.
   real, parameter :: ZERO_8  = 0._8
   real, parameter :: ONE     = 1.
   real, parameter :: ONE_8   = 1._8
   real, parameter :: FOUR    = 4.
   real, parameter :: FOUR_8  = 4._8
   real, parameter :: NINETY  = 90.
   real, parameter :: ONE80   = 180.
   real, parameter :: ONE80_8 = 180._8
   real, parameter :: THREE60 = 360.
   real, parameter :: ROUNDOFF = 10.*epsilon(ROUNDOFF)

   integer :: maxVtx
   real(8) :: PI, D2R, R2D

   include 'netcdf.inc'

   CONTAINS

   subroutine camse_init( ncid, mdl_grd )

   use anthro_types, only : mdl_poly_type, model_grid_type
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
   integer :: nVtx, nPolygons
   integer :: Quadcnt(4)
   integer, allocatable :: Quad(:)
   real    :: MatchVtxLon, MatchVtxLat
   real, allocatable :: wrk(:,:)
   real, allocatable :: wrkLon(:)
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
     if( mod(nVtx,2) /= 0 ) then
       write(*,*) 'camse_utils: odd number of polygon vertices'
       Stop 'Poly-ERR'
     endif
   enddo

   allocate( Quad(maxVtx),wrkLon(maxVtx),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'camse_init: Failed to allocate Quad,wrkLon;  error = ',astat
     stop 'AllocErr'
   endif

   mdl_grd%mdl_poly(1:nPolygons)%active = .true.
!-----------------------------------------------------------------------
!  remove polygons with edges that "cross" 180 longitude
!-----------------------------------------------------------------------
poly_loop: &
   do polyNdx = 1,nPolygons
     nVtx = mdl_grd%mdl_poly(polyNdx)%nVtx
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
     elseif( minval(wrkLon(1:nvtx)) * maxval(wrkLon(1:nVtx)) < 0. ) then
       if( any(Quad(:nVtx) == 2) .and. any(Quad(:nVtx) == 3) ) then
         mdl_grd%mdl_poly(polyNdx)%active = .false.
       endif
     endif
   enddo poly_loop

   deallocate( Quad,wrkLon )

   PI = FOUR*atan(ONE)
   D2R = PI/ONE80_8
   R2D = ONE80_8/PI

   end subroutine camse_init

   subroutine area_interp_init_camse( proj, mdlGrid, d2mMap, mdl_area_type, diag_level )
!-----------------------------------------------------------------------
!  initialize area interpolation for camse model
!-----------------------------------------------------------------------

   use anthro_types, only : mdl_poly_type, model_grid_type
   use mapper_types, only : grid_type, area_type, proj_info
   use area_mapper,  only : poly_area
   use constants_module, only : earth_radius_m, pi
!-----------------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------------
   integer,               intent(in)    :: diag_level
   type(model_grid_type), intent(inout) :: mdlGrid
   type(grid_type),       intent(in)    :: d2mMap
   type(area_type),       intent(inout) :: mdl_area_type(:)
   type(proj_info),       intent(in)    :: proj
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   type enclosingBox
     real :: minLon, maxLon
     real :: minLat, maxLat
   end type enclosingBox

   real, parameter :: msq2kmsq = 1.e-6

   integer :: astat
   integer :: polyNdx
   integer :: lonNdx, latNdx
   integer :: nVtx
   integer :: xcellCnt
   integer :: nData_n_Mdl, nMdl_n_Data
   integer :: ndVtx
   integer :: minmaxNdx(2)
   integer :: xLonNdx(d2mMap%nlons*d2mMap%nlats)
   integer :: xLatNdx(d2mMap%nlons*d2mMap%nlats)
   real    :: swghts
   real    :: wghts(mdlGrid%nPolygons)
   real    :: dataCellArea(d2mMap%nlons,d2mMap%nlats)
   real(8) :: xsectArea
   real(8) :: dtotArea, xtotArea
   real(8) :: mCellArea, dCellArea
   real(8) :: dCellAreaMin, dCellAreaMax
   real(8) :: mCellAreaMin, mCellAreaMax
   real(8) :: dPoly_x(4), dPoly_y(4)
   real(8) :: dPoly_lam(4), dPoly_phi(4)
   real(8) :: mPoly_x(10), mPoly_y(10)
   real(8) :: mPoly_lam(10), mPoly_phi(10)
   logical :: dcellInsidemcell, mcellInsidedcell
   logical :: Mask(mdlGrid%nPolygons)
   logical :: Mask1(mdlGrid%nPolygons)
   type(enclosingBox) :: mdlBox, dataBox

   real(8) :: polyintersectarea

!-----------------------------------------------------------------------
!  set enclosing "box" for overall data domain
!-----------------------------------------------------------------------
   dataBox%minLon = minval(d2mMap%xedge_2d(:,:) )
   dataBox%minLat = minval(d2mMap%yedge_2d(:,:) )
   dataBox%maxLon = maxval(d2mMap%xedge_2d(:,:) )
   dataBox%maxLat = maxval(d2mMap%yedge_2d(:,:) )

   write(*,*) ' '
   write(*,*) 'There are ',count(mdlGrid%mdl_poly(:)%active),' active polygons before box test'

!-----------------------------------------------------------------------
!  weed out mdl polygons that are completely outside data grid
!-----------------------------------------------------------------------
   do polyNdx = 1,mdlGrid%nPolygons
     if( mdlGrid%mdl_poly(polyNdx)%active ) then
       nVtx = mdlGrid%mdl_poly(polyNdx)%nVtx
!-----------------------------------------------------------------------
!  set mdl cell enclosing "box"
!-----------------------------------------------------------------------
       mdlBox%minLon = minval(mdlGrid%mdl_poly(polyNdx)%vtx_lon(1:nVtx) )
       mdlBox%minLat = minval(mdlGrid%mdl_poly(polyNdx)%vtx_lat(1:nVtx) )
       mdlBox%maxLon = maxval(mdlGrid%mdl_poly(polyNdx)%vtx_lon(1:nVtx) )
       mdlBox%maxLat = maxval(mdlGrid%mdl_poly(polyNdx)%vtx_lat(1:nVtx) )
       if( .not. Boxes_overlap() ) then
         mdlGrid%mdl_poly(polyNdx)%active = .false.
       else
         mdlGrid%mdl_poly(polyNdx)%active = .true.
       endif
     endif
   enddo

   write(*,*) ' '
   write(*,*) 'There are ',count(mdlGrid%mdl_poly(:)%active),' active polygons after  box test'

   Mask(:) = mdlGrid%mdl_poly(:)%active
   mdl_area_type(1:mdlGrid%nPolygons)%has_data = .false.
   mdl_area_type(1:mdlGrid%nPolygons)%active_dcell_cnt = 0
   mdl_area_type(1:mdlGrid%nPolygons)%interior_dcell_cnt = 0
   mdl_area_type(1:mdlGrid%nPolygons)%partial_dcell_cnt = 0

   Mask1(:) = Mask(:) .and. mdlGrid%mdl_poly(:)%nVtx == 4
   write(*,*) ' '
   write(*,*) 'There are ',count(Mask1(:)),' active quad polygons after  box test'
   write(*,*) 'Min quad area = ',msq2kmsq*earth_radius_m**2*minval(mdlGrid%mdl_poly(:)%area,mask=Mask1)
   write(*,*) 'Max quad area = ',msq2kmsq*earth_radius_m**2*maxval(mdlGrid%mdl_poly(:)%area,mask=Mask1)
   minmaxNdx(1:1) = minloc(mdlGrid%mdl_poly(:)%area,mask=Mask)
   minmaxNdx(2:2) = maxloc(mdlGrid%mdl_poly(:)%area,mask=Mask)
   write(*,*) ' '
   write(*,*) 'Min poly area,#vtx = ',msq2kmsq*earth_radius_m**2*mdlGrid%mdl_poly(minmaxNdx(1))%area, &
                                    mdlGrid%mdl_poly(minmaxNdx(1))%nVtx
   write(*,*) 'Max poly area,#vtx = ',msq2kmsq*earth_radius_m**2*mdlGrid%mdl_poly(minmaxNdx(2))%area, &
                                    mdlGrid%mdl_poly(minmaxNdx(2))%nVtx

   dtotArea = 0._8
   dCellAreaMin = 1.e-3_8; dcellAreaMax = 0._8
   ndVtx    = 4
!-----------------------------------------------------------------------
!  set dataset total area
!-----------------------------------------------------------------------
   do latNdx = 1,d2mMap%nlats
     do lonNdx = 1,d2mMap%nlons
       dPoly_x(1) = real(d2mMap%xedge_2d(lonNdx,latNdx),8)
       dPoly_x(2) = real(d2mMap%xedge_2d(lonNdx+1,latNdx),8)
       dPoly_x(3) = real(d2mMap%xedge_2d(lonNdx+1,latNdx+1),8)
       dPoly_x(4) = real(d2mMap%xedge_2d(lonNdx,latNdx+1),8)
       dPoly_y(1) = real(d2mMap%yedge_2d(lonNdx,latNdx),8)
       dPoly_y(2) = real(d2mMap%yedge_2d(lonNdx+1,latNdx),8)
       dPoly_y(3) = real(d2mMap%yedge_2d(lonNdx+1,latNdx+1),8)
       dPoly_y(4) = real(d2mMap%yedge_2d(lonNdx,latNdx+1),8)
       dPoly_lam(1:4) = dPoly_x(1:4)*D2R
       dPoly_phi(1:4) = sin(dPoly_y(1:4)*D2R)
       dCellArea = poly_area( 4, dPoly_lam, dPoly_phi )
       dtotArea = dtotArea + dCellArea
       dCellAreaMin = min( dCellAreaMin,dCellArea )
       dCellAreaMax = max( dCellAreaMax,dCellArea )
       dataCellArea(lonNdx,latNdx) = real( dCellArea,4 )
     enddo
   enddo

   write(*,*) ' '
   dCellAreaMin = msq2kmsq*earth_radius_m**2*dCellAreaMin
   dCellAreaMax = msq2kmsq*earth_radius_m**2*dCellAreaMax
   write(*,*) 'Min,Max dCell area = ',dCellAreaMin,dCellAreaMax

   xtotArea = 0._8
   mCellAreaMin = 1.e-3_8; mcellAreaMax = 0._8
   nData_n_Mdl = 0 ; nMdl_n_Data = 0
!-----------------------------------------------------------------------
!  find intersection between mdl polygon and data quads
!-----------------------------------------------------------------------
poly_loop: &
   do polyNdx = 1,mdlGrid%nPolygons
mdl_cell_is_active: &
     if( Mask(polyNdx) ) then
       nVtx = mdlGrid%mdl_poly(polyNdx)%nVtx
       mPoly_x(1:nVtx) = real(mdlGrid%mdl_poly(polyNdx)%vtx_lon(1:nVtx),8)
       mPoly_y(1:nVtx) = real(mdlGrid%mdl_poly(polyNdx)%vtx_lat(1:nVtx),8)
       mPoly_lam(1:nVtx) = mPoly_x(1:nVtx)*D2R
       mPoly_phi(1:nVtx) = sin(mPoly_y(1:nVtx)*D2R)
       mCellArea = poly_area( nVtx, mPoly_lam, mPoly_phi ) 
       mCellAreaMin = min( mCellAreaMin,mCellArea )
       mCellAreaMax = max( mCellAreaMax,mCellArea )
!-----------------------------------------------------------------------
!  "box" enclosing mdl cell
!-----------------------------------------------------------------------
       mdlBox%minLon = real(minval(mPoly_x(1:nVtx)),4)
       mdlBox%maxLon = real(maxval(mPoly_x(1:nVtx)),4)
       mdlBox%minLat = real(minval(mPoly_y(1:nVtx)),4)
       mdlBox%maxLat = real(maxval(mPoly_y(1:nVtx)),4)

       xcellCnt = 0
dlat_loop: &
       do latNdx = 1,d2mMap%nlats
dlon_loop: &
         do lonNdx = 1,d2mMap%nlons
           dPoly_x(1) = real(d2mMap%xedge_2d(lonNdx,latNdx),8)
           dPoly_x(2) = real(d2mMap%xedge_2d(lonNdx+1,latNdx),8)
           dPoly_x(3) = real(d2mMap%xedge_2d(lonNdx+1,latNdx+1),8)
           dPoly_x(4) = real(d2mMap%xedge_2d(lonNdx,latNdx+1),8)
           dPoly_y(1) = real(d2mMap%yedge_2d(lonNdx,latNdx),8)
           dPoly_y(2) = real(d2mMap%yedge_2d(lonNdx+1,latNdx),8)
           dPoly_y(3) = real(d2mMap%yedge_2d(lonNdx+1,latNdx+1),8)
           dPoly_y(4) = real(d2mMap%yedge_2d(lonNdx,latNdx+1),8)
!-----------------------------------------------------------------------
!  "box" enclosing data cell
!-----------------------------------------------------------------------
           dataBox%minLon = real(minval(dPoly_x(:)),4)
           dataBox%minLat = real(minval(dPoly_y(:)),4)
           dataBox%maxLon = real(maxval(dPoly_x(:)),4)
           dataBox%maxLat = real(maxval(dPoly_y(:)),4)
!-----------------------------------------------------------------------
!  set xsection parameters if mdl, data enclosing "boxes" overlap
!-----------------------------------------------------------------------
overlap:   if( Boxes_overlap() ) then
             dPoly_lam(1:4) = dPoly_x(1:4)*D2R
             dPoly_phi(1:4) = sin(dPoly_y(1:4)*D2R)
             xsectArea = polyintersectarea( nData_n_Mdl, nMdl_n_Data, nVtx, ndVtx, mPoly_lam, mPoly_phi, dPoly_lam, dPoly_phi )
             if( xsectArea /= ZERO_8 ) then
               xcellCnt = xcellCnt + 1
               wghts(xcellCnt) = real(xsectArea,4)/dataCellArea(lonNdx,latNdx)
               xtotArea = xtotArea + xsectArea
               xLonNdx(xcellCnt) = lonNdx ; xLatNdx(xcellCnt) = latNdx
               mdl_area_type(polyndx)%partial_dcell_cnt &
                 = mdl_area_type(polyndx)%partial_dcell_cnt + 1
             endif
           endif overlap
         enddo dlon_loop
       enddo dlat_loop
       mdl_area_type(polyNdx)%active_dcell_cnt = xcellCnt
       if( xcellCnt > 0 ) then
         allocate( mdl_area_type(polyNdx)%wght(xcellCnt), &
                   mdl_area_type(polyNdx)%dcell_lon_ndx(xcellCnt), &
                   mdl_area_type(polyNdx)%dcell_lat_ndx(xcellCnt),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'area_interp_init_camse: Failed to mdl_area_type array',astat
           Stop 'Alloc-ERR'
         endif
         mdl_area_type(polyNdx)%has_data = .true.
         mdl_area_type(polyNdx)%wght(1:xcellCnt) = wghts(1:xcellCnt)
         mdl_area_type(polyNdx)%dcell_lon_ndx(1:xcellCnt) = xLonNdx(1:xcellCnt)
         mdl_area_type(polyNdx)%dcell_lat_ndx(1:xcellCnt) = xLatNdx(1:xcellCnt)
       endif
     endif mdl_cell_is_active
   enddo poly_loop

   write(*,*) ' '
   mCellAreaMin = msq2kmsq*earth_radius_m**2*mCellAreaMin
   mCellAreaMax = msq2kmsq*earth_radius_m**2*mCellAreaMax
   write(*,*) 'Min,Max mCell area = ',mCellAreaMin,mCellAreaMax
   write(*,*) ' '
   write(*,*) 'area_interp_init_camse: nData_n_Mdl,nMdl_n_Data = ',nData_n_Mdl,nMdl_n_Data
   write(*,*) 'area_interp_init_camse: min dCell cnt = ', &
      minval( mdl_area_type(1:mdlgrid%nPolygons)%active_dcell_cnt, &
              mask=mdl_area_type(1:mdlgrid%nPolygons)%active_dcell_cnt > 0 )
   write(*,*) 'area_interp_init_camse: max dCell cnt = ', &
      maxval(mdl_area_type(1:mdlgrid%nPolygons)%active_dcell_cnt)
   write(*,*) 'area_interp_init_camse: dtotArea,xtotArea = ',dtotArea,xtotArea

   CONTAINS

   function Boxes_overlap() result(Ovrlap)
     
     logical :: Ovrlap

     logical :: NoOvrlap

     NoOvrlap = (dataBox%minLon >= mdlBox%maxLon) .or. (dataBox%maxLon <= mdlBox%minLon)
     if( .not. NoOvrlap ) then
       NoOvrlap = (dataBox%minLat >= mdlBox%maxLat) .or. (dataBox%maxLat <= mdlBox%minLat)
     endif
     Ovrlap = .not. NoOvrlap

   end function Boxes_overlap

   end subroutine area_interp_init_camse

   end module camse_utils
