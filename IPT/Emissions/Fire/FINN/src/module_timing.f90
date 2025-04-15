!WRF:DRIVER_LAYER:UTIL

MODULE module_timing

   INTEGER, PARAMETER, PRIVATE :: cnmax = 30
   INTEGER, PRIVATE :: cn
   REAL, PRIVATE    :: elapsed_seconds , elapsed_seconds_total = 0

   REAL(8)          :: epoch_seconds_hires(cnmax)

CONTAINS

   SUBROUTINE init_module_timing
!------------------------------------------------------------------
! Initialize the high-res timer.  This is optional, but will allow
! higher precision.  Read hires_timer.c for details.
!------------------------------------------------------------------
      call init_hires_timer()
      cn = 0
   END SUBROUTINE init_module_timing


   SUBROUTINE start_timing
!------------------------------------------------------------------
!  start the timer
!------------------------------------------------------------------

      IMPLICIT NONE

      cn = cn + 1
      IF ( cn > cnmax ) THEN
        write(*,*) 'module_timing: clock nesting error (too many nests)'
        RETURN
      ENDIF

      call hires_timer(epoch_seconds_hires(cn))

   END SUBROUTINE start_timing


   SUBROUTINE end_timing ( string )
!------------------------------------------------------------------
!  stop timer, report timing
!------------------------------------------------------------------
   
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)   :: string

      character(len=512) :: buf
      REAL(8)            :: now_hires

!------------------------------------------------------------------
!  check timer count
!------------------------------------------------------------------
      IF( cn < 1 ) THEN
        write(*,*) 'module_timing: clock nesting error, cn<1'
        stop 'Timing err'
      ELSE IF ( cn > cnmax ) THEN
        write(*,*) 'module_timing: clock nesting error, cn>cnmax'
        stop 'Timing err'
      ENDIF

      call hires_timer(now_hires)

      elapsed_seconds = REAL(now_hires-epoch_seconds_hires(cn),4)

      elapsed_seconds_total = elapsed_seconds_total + elapsed_seconds

      write(*,*) ' '
      write(*,'(''Timing for '',A,'': '',F10.5,'' seconds'')') &
        TRIM(string),elapsed_seconds

      cn = cn - 1

   END SUBROUTINE end_timing

   FUNCTION now_time() result(timef)
!------------------------------------------------------------------
! This is a simple subroutine that returns the current time in
! seconds since some arbitrary reference point.  This routine is
! meant to be used to accumulate timing information.
!------------------------------------------------------------------

     implicit none

     real(8) :: timef

     call hires_timer(timef)

   END FUNCTION now_time

END MODULE module_timing

