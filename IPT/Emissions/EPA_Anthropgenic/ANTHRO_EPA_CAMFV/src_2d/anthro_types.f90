
   module anthro_types

   implicit none

   integer, parameter :: linemax = 500
   integer, parameter :: linsize = 512
   integer, parameter :: namsize = 64
   integer, parameter :: maxsrc = 50

   type glb_att
     integer :: len
     integer :: type
     character(len=132)  :: name
     integer(1), pointer :: attr_byte(:)
     integer(2), pointer :: attr_short(:)
     integer, pointer    :: attr_int(:)
     real, pointer       :: attr_real(:)
     real(8), pointer    :: attr_dbl(:)
     character(len=256)  :: attr_char
   end type

   type anthro_map_type
     integer           :: src_cnt                         ! count of src species
     character(len=namsize) :: emis_name                  ! emission species name
     character(len=namsize) :: src_var(maxsrc)            ! src species names
     real              :: src_wght(maxsrc)                ! multiplier for each src species
     real, allocatable :: cat_wght(:,:)                   ! multiplier for sub cats for each src species
     real, allocatable :: emissionSrf(:,:,:)              ! surface emission
     real, allocatable :: emissionLay(:,:,:)              ! layer production
     logical           :: is_gas                          ! .t. => gas phase, .f. => aerosol
     logical           :: is_2d(maxsrc)                   ! .t. => has srf emissions
     logical           :: is_3d(maxsrc)                   ! .t. => has 3d emissions
   end type anthro_map_type

   type dates
     integer :: date
     integer :: secs
   end type dates

   type stack_type
     integer :: nStk
     integer, allocatable :: mdl_i(:)
     integer, allocatable :: mdl_j(:)
     integer, allocatable :: mdl_k(:)
     logical, allocatable :: dataMask(:)
     real, allocatable :: longitude(:)
     real, allocatable :: latitude(:)
     real, allocatable :: stkHt(:)
     real, allocatable :: emis(:,:,:)
     real, allocatable :: src_data(:)
     character(len=256) :: filespec
     character(len=128) :: filename
   end type stack_type

   type data_file_type
     integer :: flnmDate
     integer :: ntimes
     integer :: ncid_lo, ncid_hi
     integer :: lo_tndx, hi_tndx
     integer :: lo_buf_ndx, hi_buf_ndx
     integer :: gap_date, gap_secs
     integer :: grid_ndx
     integer, allocatable :: date(:)
     integer, allocatable :: secs(:)
     real    :: missing_value = 1.e36
     real    :: dels
     real    :: molecw
     real    :: dx
     real    :: con_fac(2)
     real, allocatable :: emis(:,:,:,:)
     real, allocatable :: src_data(:,:)
     character(len=namsize) :: src_name
     character(len=256) :: filespec
     character(len=128) :: filename
     character(len=32)  :: flnmDateStr
     logical :: read_lo_tndx
     logical :: read_hi_tndx
     logical :: in_gap
     logical :: t_interp
     logical :: active
     logical :: is_EPA
     logical :: is_EPA_SECTOR, is_EPA_STACK
     logical, allocatable :: cat_active(:)
     TYPE(stack_type) :: stack
   end type data_file_type

   TYPE model_grid_type
     integer :: nlons
     integer :: nlats
     integer :: nlevs
     integer :: lon_shift = 0
     real    :: dx
     real, allocatable :: mdl_lons(:,:)
     real, allocatable :: mdl_lats(:,:)
     real, allocatable :: mdl_cell_area(:,:)
     real, allocatable :: z_at_w(:,:,:)                    ! height above ground at interfaces (meters)
     real, allocatable :: landmask(:,:)                    ! land/water mask (1/0)
     character(len=namsize) :: filename                    ! name of input file for CAM
     character(len=5 ) :: mdlType                          ! 'WRF' || 'CAMFV' || 'CAMSE'
     logical :: mdl_is_WRF
     logical :: mdl_is_CAM, mdl_is_CAMFV, mdl_is_CAMSE
   END TYPE model_grid_type

   end module anthro_types
