!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module Class_mpi_wrapper
    use m_config
    use mpi_f08
    implicit none

contains
   subroutine judft_comm_split(comm, color, key, new_comm)
      use mpi_f08
      implicit none
      integer, intent(in)    :: color, key
      type(MPI_Comm), intent(in) :: comm
      type(MPI_Comm), intent(inout) :: new_comm
      integer                :: ierr, err_handler

      call MPI_Comm_Split(comm,color,key,new_comm,ierr)
      if(ierr /= 0) then
            write(*,*) "Can't split comm"
            stop
      endif
      call MPI_Comm_create_errhandler(judft_mpi_error_handler, err_handler, ierr)
      if(ierr /= 0) then
            write(*,*) "Can't create Error handler"
            stop
      endif
      call MPI_Comm_Set_Errhandler(new_comm, err_handler, ierr)
      if(ierr /= 0) then
            write(*,*) "Can't assign Error handler to new_comm"
            stop
      endif
   end subroutine judft_comm_split

   subroutine judft_mpi_error_handler(comm, error_code)
      use mpi_f08
      implicit none
      integer  :: error_code
      type(MPI_Comm) :: comm
      integer             :: str_len, ierr
      character(len=3000) :: error_str

      call MPI_ERROR_STRING(error_code, error_str, str_len, ierr)
      write (*,*) "MPI failed with Error_code = ", error_code,new_line("A"),error_str(1:str_len)
      !write(*,*) "MPI failed with Error_code = " // int2str(error_code) // new_line("A") // &
      !                  error_str(1:str_len)
       end subroutine judft_mpi_error_handler
end module Class_mpi_wrapper
    
