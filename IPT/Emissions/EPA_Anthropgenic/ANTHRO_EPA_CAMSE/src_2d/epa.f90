
   module EPA

   implicit none

   private
   public :: init_table, translate_date

   integer, private  :: ndateTypes
   integer, private  :: daysPerYear
   integer, private, allocatable :: dateTable(:,:)
   character(len=16), private, allocatable :: dateTypeTable(:,:)
   character(len=*), private, parameter :: dateHdr(7) = &
      (/ 'all     ', 'aveday_N', 'aveday_Y', 'mwdss_N ', 'mwdss_Y ', 'week_N  ', 'week_Y  ' /)

   CONTAINS

   function init_table(  srcFileName, dateFile ) result( err )
!-----------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------
     character(len=*), intent(in) :: srcFileName       ! full pathname of epa sector to date type mapping
     character(len=*), intent(in) :: dateFile          ! full pathname of epa date mapping
     logical             :: err
!-----------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------
     integer, parameter :: unitno = 33

     integer :: istat, n
     integer :: beg, end
     character(len=64 ) :: cdummy
     logical :: lexist
 
     err = .false.
!-----------------------------------------------------------------
!  check for file existance
!-----------------------------------------------------------------
     inquire( file=trim(srcFileName),exist=lexist )
     if( .not. lexist ) then
       write(*,*) 'init_table: EPA sector file'
       write(*,*) trim(srcFileName)
       write(*,*) 'init_table: does not exist'
       err = .true.
       return
     endif
     inquire( file=trim(dateFile),exist=lexist )
     if( .not. lexist ) then
       write(*,*) 'init_table: EPA date file'
       write(*,*) trim(dateFile)
       write(*,*) 'init_table: does not exist'
       err = .true.
       return
     endif
!-----------------------------------------------------------------
!  open and read file
!-----------------------------------------------------------------
   open(unit=unitno,file=trim(srcFileName),iostat=istat)
   if( istat /= 0 ) then
     write(*,*) 'init_table: error opening ',trim(srcFileName),' = ',istat
     err = .true.
     return
   endif
   n = 0
   do
     read(unitno,*,iostat=istat) cdummy
     if( istat /= 0 ) then
       exit
     endif
     n = n + 1 
   enddo
   n = n - 1
   if( n <= 0 ) then
     write(*,*) 'init_table: ',trim(srcFileName),' is empty'
     err = .true.
     return
   endif
   ndateTypes = n
   allocate( dateTypeTable(n,2) )
   rewind( unit=unitno )
   read(unitno,*,iostat=istat) cdummy
   do n = 1,ndateTypes
     read(unitno,*,iostat=istat) dateTypeTable(n,1),cdummy,cdummy,dateTypeTable(n,2)
     if( istat /= 0 ) then
       write(*,*) 'init_table: failed to read line ',n+1,' of ',trim(srcFileName)
       write(*,*) '            error = ',istat
       err = .true.
       return
     endif
   enddo

   close( unit=unitno )

   open(unit=unitno,file=trim(dateFile),iostat=istat)
   if( istat /= 0 ) then
     write(*,*) 'init_table: error opening ',trim(dateFile),' = ',istat
     err = .true.
     return
   endif
   n = 0
   do
     read(unitno,*,iostat=istat) cdummy
     if( istat /= 0 ) then
       exit
     endif
     n = n + 1 
   enddo
   n = n - 1
   if( n <= 0 ) then
     write(*,*) 'init_table: ',trim(dateFile),' is empty'
     err = .true.
     return
   endif
   daysPerYear = n
   allocate( dateTable(n,7) )
   rewind( unit=unitno )
   read(unitno,*,iostat=istat) cdummy
   do n = 1,daysPerYear
     read(unitno,*,iostat=istat) dateTable(n,1:7)
     if( istat /= 0 ) then
       write(*,*) 'init_table: failed to read line ',n,' of ',trim(dateFile)
       write(*,*) '            error = ',istat
       err = .true.
       return
     endif
   enddo

   close( unit=unitno )

   end function init_table

   function translate_date(  inDate, src ) result( outDate )
!-----------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------
     integer, intent(in) :: inDate                    ! incoming date in yyyymmdd form
     integer             :: outDate                   ! output date in yyyymmdd form
     character(len=*), intent(in) :: src              ! epa sector type
!-----------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------
   integer :: m, n
   logical :: found

   found = .false.
   do n = 1,ndateTypes
     if( trim(src) == trim(dateTypeTable(n,1)) ) then
       found = .true.
       exit
     endif
   enddo

   if( .not. found ) then
     write(*,*) 'translate_table: date type ',trim(src),' not in table; terminating'
     stop
   endif

   found = .false.
   do m = 1,7
     if( trim(dateHdr(m)) == trim(dateTypeTable(n,2)) ) then
       found = .true.
       exit
     endif
   enddo

   if( .not. found ) then
     write(*,*) 'translate_table: date type ',trim(src),' not in table; terminating'
     stop
   endif

   found = .false.
   do n = 1,daysPerYear
     if( inDate == dateTable(n,1) ) then
       found = .true.
       exit
     endif
   enddo

   if( .not. found ) then
     write(*,*) 'translate_table: date ',inDate,' not in table; terminating'
     stop
   endif

   outDate = dateTable(n,m)

   end function translate_date

   end module EPA
