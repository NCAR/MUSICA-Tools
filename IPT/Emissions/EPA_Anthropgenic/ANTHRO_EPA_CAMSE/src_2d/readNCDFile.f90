
   program readNCDFile

   implicit none

   integer :: ncid, astat, varid, dimid
   integer :: files, ncols, ntimes, nalt
   real, allocatable :: emis(:,:,:)
   character(len=256) :: outpname
   character(len=256) :: message

!---------------------------------------------------------------------
!  include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

!---------------------------------------------------------------------
!  check surface emissions first, then vertical emissions
!---------------------------------------------------------------------
   do files = 1,2
     if( files == 1 ) then
       outpname = 'emissions-EPA_CO_surface_anthro_c20190221_conus_30x8.nc'
     else
       outpname = 'emissions-EPA_CO_vertical_anthro_c20190221_conus_30x8.nc'
     endif
     write(*,*) ' '
     write(*,*) 'ReadNCDFile: Opening file ',trim(outpname)
     message = 'ReadNCDFile: Failed to open file: ' // trim(outpname)
     call handle_ncerr( nf_open( trim(outpname), nf_write, ncid ), message )       

     message = 'ReadNCDFile: Failed to get ncol dim id'
     call handle_ncerr( nf_inq_dimid( ncid, 'ncol', dimid ), message )       
     message = 'ReadNCDFile: Failed to get ncol size'
     call handle_ncerr( nf_inq_dimlen( ncid, dimid, ncols ), message )       

     write(*,*) ' '
     write(*,*) 'ReadNCDFile: ncols  = ',ncols

     message = 'ReadNCDFile: Failed to get time dim id'
     call handle_ncerr( nf_inq_dimid( ncid, 'time', dimid ), message )       
     message = 'ReadNCDFile: Failed to get time size'
     call handle_ncerr( nf_inq_dimlen( ncid, dimid, ntimes ), message )       

     write(*,*) 'ReadNCDFile: ntimes = ',ntimes

     if( files == 2 ) then
       message = 'ReadNCDFile: Failed to get altitude dim id'
       call handle_ncerr( nf_inq_dimid( ncid, 'altitude', dimid ), message )       
       message = 'ReadNCDFile: Failed to get altitude size'
       call handle_ncerr( nf_inq_dimlen( ncid, dimid, nalt ), message )       
     else
       nalt = 1
     endif

     write(*,*) 'ReadNCDFile: nalt   = ',nalt

     if( allocated( emis ) ) then
       deallocate( emis )
     endif
     allocate( emis(ncols,ntimes,nalt),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'ReadNCDFile: Failed to allocate emis; error = ',astat
       Stop 'Alloc-ERR'
     endif

     message = 'ReadNCDFile: Failed to get emiss-anthro var id'
     call handle_ncerr( nf_inq_varid( ncid, 'emiss-anthro', varid ), message )       
     message = 'ReadNCDFile: Failed to read emiss-anthro'
     call handle_ncerr( nf_get_var_real( ncid, varid, emis ), message )       

     write(*,*) ' '
     write(*,*) 'ReadNCDFile: < 0 count = ',count(emis(:,:,:) < 0.)
     write(*,*) 'ReadNCDFile: = 0 count = ',count(emis(:,:,:) == 0.)
     write(*,*) 'ReadNCDFile: > 0 count = ',count(emis(:,:,:) > 0.)

     message = 'ReadNCDFile: Failed to close file: ' // trim(outpname)
     call handle_ncerr( nf_close(ncid), message )       
   enddo
   
   CONTAINS

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

   end program readNCDFile
