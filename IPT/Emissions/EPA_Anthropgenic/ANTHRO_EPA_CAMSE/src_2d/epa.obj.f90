
   module EPA

   implicit none

   private

   type :: EPATiming
     integer, private, allocatable :: dateTable(:,:)
     character(len=16), private, allocatable :: dateTypeTable(:,:)
     integer, private  :: ndateTypes
     integer, private  :: daysPerYear
     CONTAINS
     private
     procedure, public :: init_table
     procedure, public :: translate_date
   end type EPATiming

   CONTAINS

   function init_table( this, srcFileName, dateFile ) result( err )
!-----------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------
     class(EPATiming)    :: this
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
   do
     read(unit=unitno,*,iostat=istat) cdummy
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
   this%ndateTypes = n
   allocate( this%dateTypeTable(n,2) )
   rewind( unit=unitno )
   read(unit=unitno,*,iostat=istat) cdummy
   do n = 1,this%ndateTypes
     read(unit=unitno,*,iostat=istat) this%dateTypeTable(n,1),cdummy,cdummy,this%dateTypeTable(n,2)
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
   do
     read(unit=unitno,*,iostat=istat) cdummy
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
   this%daysPerYear = n
   allocate( this%dateTable(n,7) )
   rewind( unit=unitno )
   read(unit=unitno,*,iostat=istat) cdummy
   do n = 1,this%daysPerYear
     read(unit=unitno,*,iostat=istat) this%dateTable(n,1:7)
     if( istat /= 0 ) then
       write(*,*) 'init_table: failed to read line ',n,' of ',trim(dateFile)
       write(*,*) '            error = ',istat
       err = .true.
       return
     endif
   enddo

   close( unit=unitno )

   end function init_table

   function translate_date( this, inDate, src ) result( outDate )
!-----------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------
     class(EPATiming)    :: this
     integer, intent(in) :: inDate                    ! incoming date in yyyymmdd form
     integer             :: outDate                   ! output date in yyyymmdd form
     character(len=*), intent(in) :: src              ! epa sector type

   end function translate_date

   end module EPA
