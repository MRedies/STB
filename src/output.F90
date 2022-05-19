Module  output
   use mpi
   use m_npy
   use Class_helper
   implicit None
   integer(4), parameter  :: std_out =  6
   interface print_mtx
      module procedure print_mtx_real, print_mtx_cmplx, print_vec_real,&
         print_mtx_real_no_unit, print_mtx_cmplx_no_unit, &
         print_vec_real_no_unit, print_vec_int, &
         print_vec_int_no_unit, print_mtx_int_no_unit, &
         print_mtx_int, print_vec_cmplx, print_vec_cmplx_no_unit
   end interface
contains
   subroutine print_mtx_int_no_unit(mtx)
      Implicit None
      integer(4), dimension(:,:), intent(in)    :: mtx

      call print_mtx_int(std_out,  mtx)
   end subroutine

   Subroutine  print_mtx_int(p_unit, mtx)
      Implicit None
      integer(4), dimension(:,:), intent(in)    :: mtx
      integer(4), intent(in)                 :: p_unit
      integer(4)                             :: i,j

      do i = 1, size(mtx,1)
         do j =  1, size(mtx,2)
            if(j < size(mtx,2)) then
               write (p_unit, "(I6)", advance="no") mtx(i,j)
            endif
         enddo
         write(p_unit, "(I6)") mtx(i, size(mtx,2))
      enddo

   End Subroutine print_mtx_int

   subroutine print_vec_int_no_unit(vec)
      implicit none
      integer(4), dimension(:), intent(in)   :: vec

      call print_vec_int(std_out, vec)
   end subroutine print_vec_int_no_unit

   subroutine print_vec_int(p_unit, vec)
      implicit none
      integer(4), dimension(:), intent(in)   :: vec
      integer(4), intent(in)                 :: p_unit
      integer(4)                             :: i

      do i =  1,size(vec)
         write(p_unit,"(I6)") vec(i)
      enddo
   end subroutine print_vec_int

   subroutine print_vec_real_no_unit(vec)
      implicit none
      real(8), dimension(:), intent(in)      :: vec

      call print_vec_real(std_out,  vec)
   end subroutine print_vec_real_no_unit

   subroutine print_vec_real(p_unit, vec)
      implicit none
      real(8), dimension(:), intent(in)      :: vec
      integer(4), intent(in)                 :: p_unit
      integer(4)                             :: i

      do i =  1, size(vec)
         write(p_unit, "(ES10.3)") vec(i)
      enddo
   end subroutine print_vec_real

   subroutine print_vec_cmplx_no_unit(vec)
      implicit none
      complex(8), dimension(:), intent(in)      :: vec

      call print_vec_cmplx(std_out,  vec)
   end subroutine print_vec_cmplx_no_unit

   subroutine print_vec_cmplx(p_unit, vec)
      implicit none
      complex(8), dimension(:), intent(in)      :: vec
      integer(4), intent(in)                 :: p_unit
      integer(4)                             :: i

      do i =  1, size(vec)
         write(p_unit, "(A, ES10.3, A, ES10.3, A)")&
            "(", real(vec(i)), "+ i*", aimag(vec(i)), ")"
      enddo
   end subroutine print_vec_cmplx

   subroutine print_mtx_real_no_unit(mtx)
      Implicit None
      real(8), dimension(:,:), intent(in)    :: mtx

      call print_mtx_real(std_out,  mtx)
   end subroutine

   Subroutine  print_mtx_real(p_unit, mtx)
      Implicit None
      real(8), dimension(:,:), intent(in)    :: mtx
      integer(4), intent(in)                 :: p_unit
      integer(4)                             :: i,j

      do i = 1, size(mtx,1)
         do j =  1, size(mtx,2)
            if(j < size(mtx,2)) then
               write (p_unit, "(ES18.5)", advance="no") mtx(i,j)
            endif
         enddo
         write(p_unit, "(ES18.5)") mtx(i, size(mtx,2))
      enddo

   End Subroutine print_mtx_real

   subroutine print_mtx_cmplx_no_unit(mtx)
      Implicit None
      complex(8), dimension(:,:), intent(in)    :: mtx

      call print_mtx_cmplx(std_out,  mtx)
   end subroutine print_mtx_cmplx_no_unit

   Subroutine  print_mtx_cmplx(p_unit, mtx)
      Implicit None
      complex(8), dimension(:,:), intent(in)    :: mtx
      integer(4), intent(in)                     :: p_unit
      integer(4)                                 :: i,j
      character(len=4)                           :: i_str
      character(len=4)                           :: deli = "j), "

      do i = 1, size(mtx,1)
         do j =  1, size(mtx,2)
            if(aimag(mtx(i,j))>=0.0d0) then
               i_str =  " +i*"
            else
               i_str =  " -i*"
            endif

            if(j < size(mtx,2)) then
               write (p_unit, "(SP,a,ES10.3,a,ES10.3,a)", advance="no") &
                  "(",real(mtx(i,j)), " ", aimag(mtx(i,j)), deli
            endif
         enddo
         write(p_unit, "(SP,a,ES15.8,a,ES15.8,a)") &
            "(",real(mtx(i, size(mtx,2))), " ", &
            aimag(mtx(i, size(mtx,2))), deli
      enddo
   End Subroutine print_mtx_cmplx

   subroutine create_dir(folder)
      implicit none
      character(len=*) :: folder
      integer(4)       :: succ

      call execute_command_line("rm " // trim(folder) // "*.npy", exitstat=succ)
      if(succ /= 0) then
         write (*,*) "Could not clear dir"
      endif
      call execute_command_line("mkdir -p " // folder, exitstat=succ)
      if(succ /= 0) then
         write (*,*) "Could not reate dir through cmd"
      endif
   end subroutine create_dir

End Module  output
