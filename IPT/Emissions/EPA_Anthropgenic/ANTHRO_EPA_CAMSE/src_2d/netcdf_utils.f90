
   module netcdf_utils

   use anthro_types, only : glb_att

   implicit none

   integer :: ngatts
   type(glb_att), allocatable :: attrs(:)

   private
   public :: get_glb_atts, set_glb_atts, handle_ncerr
   public :: ngatts

!---------------------------------------------------------------------
!  include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   CONTAINS

   subroutine get_glb_atts( ncid )
!---------------------------------------------------------------------
!   read the global attributes
!---------------------------------------------------------------------

   implicit none

   integer, intent(in) :: ncid

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: astat
   integer              :: attr_len
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   get global attr count
!---------------------------------------------------------------------
   message = 'glb_attr: Failed to get glb attr count'
   call handle_ncerr( nf_inq_natts( ncid, ngatts ), message )       
!---------------------------------------------------------------------
!   allocate variables
!---------------------------------------------------------------------
   if( ngatts > 0 ) then
     if( allocated( attrs ) ) then
       call dealloc_glb_atts
     endif
     allocate( attrs(ngatts),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'glb_attr: failed to allocate type glb_att'
       stop 'Alloc err'
     endif
     attrs(:)%name = ' '
   endif
!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' name'
     call handle_ncerr( nf_inq_attname( ncid, nf_global, m, attr_name ), message )       
     attrs(m)%name = attr_name
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' type,len'
     call handle_ncerr( nf_inq_att( ncid, nf_global, trim(attr_name), attrs(m)%type, attr_len ), message )       
     attrs(m)%len = attr_len
     message = 'glb_attr: Failed to get ' // trim(attr_name)
     select case( attrs(m)%type )
       case( nf_byte )
         allocate( attrs(m)%attr_byte(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_byte'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int1( ncid, nf_global, trim(attr_name), attrs(m)%attr_byte ), message )       
       case( nf_char )
         attrs(m)%attr_char = ' '
         call handle_ncerr( nf_get_att_text( ncid, nf_global, trim(attr_name), attrs(m)%attr_char ), message )       
       case( nf_short )
         allocate( attrs(m)%attr_short(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_short'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int2( ncid, nf_global, trim(attr_name), attrs(m)%attr_short ), message )       
       case( nf_int )
         allocate( attrs(m)%attr_int(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_int'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int( ncid, nf_global, trim(attr_name), attrs(m)%attr_int ), message )       
       case( nf_float )
         allocate( attrs(m)%attr_real(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_real'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_real( ncid, nf_global, trim(attr_name), attrs(m)%attr_real ), message )       
       case( nf_double )
         allocate( attrs(m)%attr_dbl(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_dbl'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_double( ncid, nf_global, trim(attr_name), attrs(m)%attr_dbl ), message )       
     end select
   end do glb_attr_loop

   end subroutine get_glb_atts

   subroutine set_glb_atts( ncid )
!---------------------------------------------------------------------
!   set the global attributes
!---------------------------------------------------------------------

   implicit none

   integer, intent(in) :: ncid

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: ncd_err
   integer              :: attr_len
   integer              :: attr_xtype
   integer              :: slen
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     attr_name = trim(attrs(m)%name)
     if( trim(attr_name) == 'TITLE' .or. trim(attr_name) == 'START_DATE' .or. &
         trim(attr_name) == 'SIMULATION_START_DATE' ) then
       cycle
     endif
     slen      = len_trim(attr_name)
     write(message,*) 'set_glb_att: Failed to define glb att ',trim(attr_name)
     attr_len   = attrs(m)%len
     attr_xtype = attrs(m)%type
     select case( attrs(m)%type )
       case( nf_byte )
         ncd_err = nf_put_att_int1( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_byte )
       case( nf_char )
         ncd_err = nf_put_att_text( ncid, nf_global, attr_name(:slen), attr_len, attrs(m)%attr_char )
       case( nf_short )
         ncd_err = nf_put_att_int2( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_short )
       case( nf_int )
         ncd_err = nf_put_att_int( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_int )
       case( nf_float )
         ncd_err = nf_put_att_real( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_real )
       case( nf_double )
         ncd_err = nf_put_att_double( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_dbl )
     end select
     call handle_ncerr( ncd_err, message )       
   end do glb_attr_loop

   end subroutine set_glb_atts

   subroutine dealloc_glb_atts
!---------------------------------------------------------------------
!   deallocate variables
!---------------------------------------------------------------------

   integer :: m

   if( allocated( attrs ) ) then
     do m = 1,ngatts
       select case( attrs(m)%type )
         case( nf_byte )
           deallocate( attrs(m)%attr_byte )
         case( nf_short )
           deallocate( attrs(m)%attr_short )
         case( nf_int )
           deallocate( attrs(m)%attr_int )
         case( nf_float )
           deallocate( attrs(m)%attr_real )
         case( nf_double )
           deallocate( attrs(m)%attr_dbl )
       end select
     end do
     deallocate( attrs )
   endif

   end subroutine dealloc_glb_atts

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
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

   end module netcdf_utils
