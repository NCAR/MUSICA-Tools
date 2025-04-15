
   module attr_types

   implicit none

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
   end type glb_att

   type dom_glb_att
     integer :: ngatts
     type(glb_att), pointer :: attrs(:)
   end type dom_glb_att

   type UtilAttributes
     character(len=32)      :: Model = 'WRF'
     character(len=32)      :: resol = ' '
     character(len=32)      :: EmisType = 'bb_surface'
     character(len=16)      :: FinnVers = '1.5'
   end type UtilAttributes

   end module attr_types
