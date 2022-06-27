module class_Units
   use m_config
   use Constants
   use mpi_f08
   implicit none

   type, public :: units
      real(8) :: length, energy, inv_length, temperature, mag_dipol
   end type units
contains
   function init_units(cfg, me) result(ret)
      implicit none
      type(units)            :: ret
      type(CFG_t), intent(in):: cfg
      integer   , intent(in) :: me

      ret%length      = get_unit_conv("length",      cfg, me, .True.)
      ret%energy      = get_unit_conv("energy",      cfg, me, .True.)
      ret%inv_length  = get_unit_conv("inv_length",  cfg, me, .True.)
      ret%temperature = get_unit_conv("temperature", cfg, me, .True.)
      ret%mag_dipol   = get_unit_conv("mag_dipol",   cfg, me, .True.)
   end function

   function get_unit_conv(field_name, cfg, me, bcast) result(factor)
      implicit none
      character(len=*), intent(in)        :: field_name
      integer   , intent(in)              :: me
      logical, optional                   :: bcast
      logical                             :: bcast_loc
      type(CFG_t)                         :: cfg
      real(8)                             :: factor
      integer                             :: ierr
      character(len=300)                  :: unit_name

      if(present(bcast)) then
         bcast_loc = bcast
      else
         bcast_loc = .True.
      endif

      if(me == root) then
         call CFG_get(cfg, "units%" // trim(field_name), unit_name)

         select case(trim(unit_name))
         case ("a0")
            factor = 1.0d0
         case ("eV")
            factor = 0.03674932d0
         case ("Hartree")
            factor = 1d0
         case ("a0^-1")
            factor = 1.0d0
         case ("K")
            factor = 1.0d0
         case ("mu_B")
            factor = 2d0
            !case ("atm_M")
            !factor = 1d0
            !case ("mu_0")
            !factor = 0.1591549430d0 ! 1/(2pi)
         case default
            write (*,*) "Unit unknown"
            write (*,*) "Unit: ", trim(unit_name)
            write (*,*) "Field name: ", trim(field_name)
            call MPI_Abort(MPI_COMM_WORLD, 13, ierr)
         end select
      endif

      if(bcast_loc) then
         call MPI_Bcast(factor, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr)
         if(ierr /= 0) then
            write(*,*) me, "Unit Bcast failed"
            stop
         endif
      endif

   end function get_unit_conv
end module class_Units

