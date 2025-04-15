
   module fire_types

   implicit none

   integer, parameter :: namsize = 256
   integer, parameter :: linsize = 512

   TYPE mdl_poly_type
     integer :: nVtx
     real    :: cntr_lon
     real    :: cntr_lat
     real    :: area
     real, allocatable :: vtx_lon(:)
     real, allocatable :: vtx_lat(:)
     logical :: active
     logical :: x180
   END TYPE mdl_poly_type

   TYPE model_grid_type
     integer :: nlons
     integer :: nlats
     integer :: nlevs
     integer :: maxPolyVtx
     integer :: nPolygons
     integer :: lon_shift = 0
     integer :: PolarCapNdx(2) = 0
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
     TYPE(mdl_poly_type), allocatable :: mdl_poly(:)
   END TYPE model_grid_type

   end module fire_types
