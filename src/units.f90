module class_Units
   use m_config
   use Constants
   use mpi_f08
   use stdlib_kinds, only: sp,dp,xdp,int32
   implicit none

   type, public :: units
      real(dp) :: length, energy, inv_length, temperature, mag_dipol
   end type units
contains
   function init_units(cfg, me, comm) result(ret)
      implicit none
      type(units)                   :: ret
      type(CFG_t), intent(in)       :: cfg
      integer(int32)   , intent(in) :: me
      type(MPI_Comm)   , intent(in) :: comm
      ret%length      = get_unit_conv("length",      cfg, me, comm, .True.)
      ret%energy      = get_unit_conv("energy",      cfg, me, comm, .True.)
      ret%inv_length  = get_unit_conv("inv_length",  cfg, me, comm, .True.)
      ret%temperature = get_unit_conv("temperature", cfg, me, comm, .True.)
      ret%mag_dipol   = get_unit_conv("mag_dipol",   cfg, me, comm, .True.)
   end function

   function get_unit_conv(field_name, cfg, me, comm, bcast) result(factor)
      implicit none
      character(len=*), intent(in)  :: field_name
      integer(int32)  , intent(in)  :: me
      type(MPI_Comm)  , intent(in)  :: comm
      type(CFG_t)                   :: cfg
      logical, optional             :: bcast
      logical                       :: bcast_loc
      real(dp)                      :: factor
      integer(int32)                :: ierr
      character(len=300)            :: unit_name

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
            call MPI_Abort(comm, 13, ierr)
         end select
      endif

      if(bcast_loc) then
         call MPI_Bcast(factor, 1, MPI_REAL8, root, comm, ierr)
         if(ierr /= 0) then
            write(*,*) me, "Unit Bcast failed"
            stop
         endif
      endif

   end function get_unit_conv
end module class_Units

