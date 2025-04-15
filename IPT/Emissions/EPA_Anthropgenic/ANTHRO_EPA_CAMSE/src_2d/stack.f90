
   module stack

   use netcdf_utils, only : handle_ncerr

   implicit none

   real, parameter :: ONE     = 1.
   real, parameter :: FOUR    = 4.
   real, parameter :: ONE80   = 180.

   real :: PI, D2R

   include 'netcdf.inc'

   CONTAINS

   subroutine stkFileRead( stkFile )

   use anthro_types, only : stack_type

   TYPE(stack_type), intent(inout) :: stkFile

   integer :: ncid
   integer :: dimid, varid
   integer :: status
   character(len=132) :: varname
   character(len=132) :: message

   PI = FOUR*atan(ONE)
   D2R = PI/ONE80

   write(*,*) ' '
   write(*,*) 'stkFileRead: Reading stack file ' // trim(stkFile%filename)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   open stack dataset file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   message = 'stkFileRead: Failed to open ' // trim(stkFile%filespec)
   call handle_ncerr( nf_open( trim(stkFile%filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get stack dataset dimesions
!---------------------------------------------------------------------
   message = 'stkFileRead: Failed to get ROW dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'ROW', dimid ), message )
   message = 'stkFileRead: Failed to get ROW dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, stkFile%nStk ), message )
!---------------------------------------------------------------------
!   allocate stack longitudes
!---------------------------------------------------------------------
   if( allocated( stkFile%longitude ) ) then
     deallocate( stkFile%longitude )
   endif
   allocate( stkFile%longitude(stkFile%nStk),stat=status )
   if( status /= 0 ) then
     write(*,*) 'stkFileRead: Failed to allocate longitude; error = ',status
     stop 'allocate failed'
   endif
!---------------------------------------------------------------------
!   allocate stack latitudes
!---------------------------------------------------------------------
   if( allocated( stkFile%latitude ) ) then
     deallocate( stkFile%latitude )
   endif
   allocate( stkFile%latitude(stkFile%nStk),stat=status )
   if( status /= 0 ) then
     write(*,*) 'stkFileRead: Failed to allocate latitude; error = ',status
     stop 'allocate failed'
   endif
!---------------------------------------------------------------------
!   allocate stack height
!---------------------------------------------------------------------
   if( allocated( stkFile%stkHt ) ) then
     deallocate( stkFile%stkHt )
   endif
   allocate( stkFile%stkHt(stkFile%nStk),stat=status )
   if( status /= 0 ) then
     write(*,*) 'stkFileRead: Failed to allocate stkHt; error = ',status
     stop 'allocate failed'
   endif
!---------------------------------------------------------------------
!   allocate model i,j,k
!---------------------------------------------------------------------
   if( allocated( stkFile%mdl_i ) ) then
     deallocate( stkFile%mdl_i )
   endif
   allocate( stkFile%mdl_i(stkFile%nStk),stat=status )
   if( status /= 0 ) then
     write(*,*) 'stkFileRead: Failed to allocate mdl_i; error = ',status
     stop 'allocate failed'
   endif
   if( allocated( stkFile%mdl_j ) ) then
     deallocate( stkFile%mdl_j )
   endif
   allocate( stkFile%mdl_j(stkFile%nStk),stat=status )
   if( status /= 0 ) then
     write(*,*) 'stkFileRead: Failed to allocate mdl_j; error = ',status
     stop 'allocate failed'
   endif
   if( allocated( stkFile%mdl_k ) ) then
     deallocate( stkFile%mdl_k )
   endif
   allocate( stkFile%mdl_k(stkFile%nStk),stat=status )
   if( status /= 0 ) then
     write(*,*) 'stkFileRead: Failed to allocate mdl_k; error = ',status
     stop 'allocate failed'
   endif

   stkFile%mdl_i(:) = 0 ; stkFile%mdl_j(:) = 0 ; stkFile%mdl_k(:) = 0
!---------------------------------------------------------------------
!   allocate stack mask
!---------------------------------------------------------------------
   if( allocated( stkFile%dataMask ) ) then
     deallocate( stkFile%dataMask )
   endif
   allocate( stkFile%dataMask(stkFile%nStk),stat=status )
   if( status /= 0 ) then
     write(*,*) 'stkFileRead: Failed to allocate dataMask; error = ',status
     stop 'allocate failed'
   endif

   stkFile%dataMask(:) = .false.

!---------------------------------------------------------------------
!   read stack longitude variable
!---------------------------------------------------------------------
   message = 'stkFileRead: Failed to get longitude id'
   call handle_ncerr( nf_inq_varid( ncid, 'LONGITUDE', varid ), message )
   message = 'stkFileRead: Failed to read longitude variable'
   call handle_ncerr( nf_get_var_real( ncid, varid, stkFile%longitude ), message )
!---------------------------------------------------------------------
!   read stack latitude variable
!---------------------------------------------------------------------
   message = 'stkFileRead: Failed to get latitude id'
   call handle_ncerr( nf_inq_varid( ncid, 'LATITUDE', varid ), message )
   message = 'stkFileRead: Failed to read latitude variable'
   call handle_ncerr( nf_get_var_real( ncid, varid, stkFile%latitude ), message )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   read stack height variable
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   message = 'stkFileRead: Failed to get height id'
   call handle_ncerr( nf_inq_varid( ncid, 'STKHT', varid ), message )
   message = 'stkFileRead: Failed to read height variable'
   call handle_ncerr( nf_get_var_real( ncid, varid, stkFile%stkHt ), message )

   status = nf_close( ncid )

   end subroutine stkFileRead

   subroutine stkFileInit( stkFile, mdlGrid, emissions_zdim_stag )

   use anthro_types, only : stack_type, model_grid_type
   use area_mapper,  only : ll_2_ij

   integer, intent(in)   :: emissions_zdim_stag
   TYPE(stack_type    )  :: stkFile
   TYPE(model_grid_type) :: mdlGrid

   real, parameter :: lower = .5

   integer :: im, jm, k, ktop
   integer :: stk
   real    :: x, y
   real    :: upperX, upperY
   real    :: stkZ
   real    :: gridZ(mdlGrid%nlevs+1)

   write(*,*) ' '
   write(*,*) 'stkFileInit: Initializing type for stack file ' // trim(stkFile%filename)

   if( mdlGrid%mdl_is_CAMFV ) then
     upperX = real(mdlGrid%nlons) + lower
     upperY = real(mdlGrid%nlats) + lower
     call stkFileInit_camfv( stkFile )
   elseif( mdlGrid%mdl_is_CAMSE ) then
     call stkFileInit_camse( stkFile, mdlGrid )
   endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   stack in model vertical domain?
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stk_loop: &
   do stk = 1,stkFile%nstk
     if( stkFile%dataMask(stk) ) then
       ktop = mdlGrid%nlevs
       im = stkFile%mdl_i(stk) ; jm = stkFile%mdl_j(stk)
       if( mdlGrid%mdl_is_CAM ) then
         gridZ(:) = mdlGrid%z_at_w(1,:,1)
       else
         gridZ(:) = mdlGrid%z_at_w(im,:,jm)
       endif
       stkZ = stkFile%stkHt(stk)
       if( stkZ >= gridZ(1) .and. stkZ <= gridZ(ktop+1) ) then
         do k = 2,ktop+1
           if( stkZ <= gridZ(k) ) then
             stkFile%mdl_k(stk) = k - 1
             exit
           endif 
         enddo
         if( stkFile%mdl_k(stk) > emissions_zdim_stag ) then
           stkFile%dataMask(stk) = .false.
         endif
       else
         stkFile%dataMask(stk) = .false.
       endif
     endif
   enddo stk_loop

   if( mdlGrid%mdl_is_CAMSE ) then
     if( any( .not. stkFile%dataMask(:) ) ) then
       write(*,*) '**************************'
       write(*,*) 'stkFileInit: failed to locate all stacks on camse grid'
       write(*,*) 'stkFileInit: located ',count(stkFile%dataMask(:)),' stacks'
       write(*,*) 'stkFileInit: Total stacks = ',stkFile%nstk
       write(*,*) '**************************'
       Stop 'stkFileInit-FAILED'
     endif
   endif

   CONTAINS

   subroutine stkFileInit_camfv( stkFile )

   use anthro_types, only : stack_type

   type(stack_type), intent(inout) :: stkFile

stk_loop: &
   do stk = 1,stkFile%nstk
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   stack in model horizontal domain?
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     call ll_2_ij( stkFile%latitude(stk), stkFile%longitude(stk), x, y ) 
     if( (lower <= x .and. x < upperX) .and. (lower <= y .and. y < upperY) ) then
       stkFile%dataMask(stk) = .true.
       stkFile%mdl_i(stk) = nint(x)
       stkFile%mdl_j(stk) = nint(y)
     else
       stkFile%dataMask(stk) = .false.
     endif
   enddo stk_loop

   end subroutine stkFileInit_camfv

   subroutine stkFileInit_camse( stkFile, mdlGrd )

   use anthro_types, only : stack_type, model_grid_type
   use area_mapper,  only : pnt_in_triangle

   TYPE(stack_type), intent(inout)   :: stkFile
   TYPE(model_grid_type), intent(in) :: mdlGrd

   type enclosingBox
     real :: minLon, maxLon
     real :: minLat, maxLat
   end type enclosingBox

   integer :: nVtx, polyNdx
   integer :: nStk, stkNdx
   integer :: triL, triU
   real    :: x(4), y(4)
   real    :: xr(4), yr(4)
   logical :: InTriangle
   logical :: stkMask(stkFile%nstk)
   type(enclosingBox) :: mdlBox


   nStk = stkFile%nstk
   stkMask(:) = .false.

polyLoop: &
   do polyNdx = 1,mdlGrd%nPolygons
     if( all( stkMask(:)) ) then
       exit polyLoop
     endif
ActiveMdlPoly: &
     if( mdlGrd%mdl_poly(polyNdx)%active ) then
       nVtx = mdlGrd%mdl_poly(polyNdx)%nVtx
       mdlBox%minLon = minval(mdlGrd%mdl_poly(polyNdx)%vtx_lon(1:nVtx))
       mdlBox%maxLon = maxval(mdlGrd%mdl_poly(polyNdx)%vtx_lon(1:nVtx))
       mdlBox%minLat = minval(mdlGrd%mdl_poly(polyNdx)%vtx_lat(1:nVtx))
       mdlBox%maxLat = maxval(mdlGrd%mdl_poly(polyNdx)%vtx_lat(1:nVtx))
       xr(1) = mdlGrd%mdl_poly(polyNdx)%cntr_lon*D2R
       yr(1) = sin(mdlGrd%mdl_poly(polyNdx)%cntr_lat*D2R)
stkLoop: &
       do stkNdx = 1,nStk
ActiveStk: &
         if( .not. stkMask(stkNdx) ) then
           x(4) = stkFile%longitude(stkNdx)
           y(4) = stkFile%latitude(stkNdx)
           if( x(4) < mdlBox%minLon .or. x(4) > mdlBox%maxLon ) then
             cycle stkLoop
           elseif( y(4) < mdlBox%minLat .or. y(4) > mdlBox%maxLat ) then
             cycle stkLoop
           endif
           xr(4) = x(4)*D2R
           yr(4) = sin(y(4)*D2R)
TriLoop:   do TriL = 1,nVtx
             triU = mod( TriL,nVtx ) + 1 
             xr(2) = mdlGrd%mdl_poly(polyNdx)%vtx_lon(TriL)*D2R
             yr(2) = sin(mdlGrd%mdl_poly(polyNdx)%vtx_lat(TriL)*D2R)
             xr(3) = mdlGrd%mdl_poly(polyNdx)%vtx_lon(TriU)*D2R
             yr(3) = sin(mdlGrd%mdl_poly(polyNdx)%vtx_lat(TriU)*D2R)
             InTriangle = pnt_in_triangle( xr,yr )
             if( InTriangle ) then
               stkFile%mdl_i(stkNdx) = polyNdx
               stkFile%mdl_j(stkNdx) = 1
               stkMask(stkNdx) = .true.
               exit TriLoop
             endif
           enddo TriLoop
         endif ActiveStk
       enddo stkLoop
     endif ActiveMdlPoly
   enddo polyLoop

   stkFile%dataMask(:) = stkMask(:)

   end subroutine stkFileInit_camse

   end subroutine stkFileInit

   subroutine stkFileCleanup( stkFile )

   use anthro_types, only : stack_type

   TYPE(stack_type)  :: stkFile

   if( allocated( stkFile%mdl_i ) ) then
     deallocate( stkFile%mdl_i )
   endif
   if( allocated( stkFile%mdl_j ) ) then
     deallocate( stkFile%mdl_j )
   endif
   if( allocated( stkFile%mdl_k ) ) then
     deallocate( stkFile%mdl_k )
   endif
   if( allocated( stkFile%longitude ) ) then
     deallocate( stkFile%longitude )
   endif
   if( allocated( stkFile%latitude ) ) then
     deallocate( stkFile%latitude )
   endif
   if( allocated( stkFile%stkHt ) ) then
     deallocate( stkFile%stkHt )
   endif
   if( allocated( stkFile%emis ) ) then
     deallocate( stkFile%emis )
   endif
   if( allocated( stkFile%src_data ) ) then
     deallocate( stkFile%src_data )
   endif

   end subroutine stkFileCleanup

   end module stack
