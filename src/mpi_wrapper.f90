!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module Class_mpi_wrapper
    use mpi
contains   
   subroutine judft_comm_split(comm, color, key, new_comm)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      integer, intent(in)    :: comm, color, key
      integer, intent(inout) :: new_comm
#ifdef CPP_MPI
      integer                :: ierr, err_handler

      call MPI_COMM_SPLIT(comm,color,key,new_comm,ierr)
      if(ierr /= 0) call judft_error("Can't split comm")

      call MPI_Comm_create_errhandler(judft_mpi_error_handler, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't create Error handler")

      call MPI_Comm_Set_Errhandler(new_comm, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to new_comm")
#endif
   end subroutine judft_comm_split

   subroutine judft_mpi_error_handler(comm, error_code)
    #ifdef CPP_MPI
          use mpi
    #endif
          use m_judft
          implicit none
          integer  :: comm, error_code
          integer             :: str_len, ierr
          character(len=3000) :: error_str
    
    #ifdef CPP_MPI
          call MPI_ERROR_STRING(error_code, error_str, str_len, ierr)
          call judft_error("MPI failed with Error_code = " // int2str(error_code) // new_line("A") // &
                           error_str(1:str_len))
    #endif
       end subroutine judft_mpi_error_handler
end module Class_mpi_wrapper
    
