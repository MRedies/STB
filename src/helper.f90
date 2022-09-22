module Class_helper
   use m_config
   use Constants
   use mpi_f08
   use stdlib_kinds, only: sp,dp,xdp,int32
   implicit none

   character(len=1), parameter :: c_esc = achar(27)
   character(len=2), parameter :: c_start = c_esc // '['
   character(len=1), parameter :: c_end = 'm'
   character(len=*), parameter :: c_black = '30'
   character(len=*), parameter :: c_red = '31'
   character(len=*), parameter :: c_green = '32'
   character(len=*), parameter :: c_yellow = '33'
   character(len=*), parameter :: c_blue = '34'
   character(len=*), parameter :: c_magenta = '35'
   character(len=*), parameter :: c_cyan = '36'
   character(len=*), parameter :: c_white = '37'
   character(len=*), parameter :: c_clear = c_start // '0' // c_end

contains
   function cnorm2(vec) result(norm)
      implicit none
      complex(dp), intent(in)   :: vec(:)
      complex(dp)               :: norm

      ! see fortran defition for complex dot_product
      norm = sqrt(dot_product(vec,vec))
   end function cnorm2

   pure function my_norm2(vec) result(norm)
      implicit none
      real(xdp), intent(in) :: vec(:)
      real(xdp)             :: norm

      norm = sqrt(dot_product(vec,vec))
   end function my_norm2

   function omp_matvec(A,x) result(b)
      implicit none
      complex(dp), intent(in)    :: A(:,:), x(:)
      complex(dp)                :: b(size(x))
      integer(int32)                :: i, j

!$OMP         parallel do private(j) shared(A,x,b)
      do i = 1,size(x)
         b(i) = 0d0
         do j = 1,size(x)
            b(i) =  b(i) +  A(i,j) * x(j)
         enddo
      enddo
   end function omp_matvec

   function matvec(A,x) result(b)
      implicit none
      complex(dp), intent(in)    :: A(:,:), x(:)
      complex(dp)                :: b(size(x))
      integer(int32)                :: i, j

      do i = 1,size(x)
         b(i) = 0d0
         do j = 1,size(x)
            b(i) =  b(i) +  A(i,j) * x(j)
         enddo
      enddo
   end function matvec

   subroutine calc_zheevd_size(vn_flag, H, eig_val, lwork, lrwork, liwork)
      implicit none
      character(len=1)          :: vn_flag
      complex(dp), intent(in)    :: H(:,:)
      real(dp), intent(in)       :: eig_val(:)
      integer, intent(out)      :: lwork, lrwork, liwork
      complex(dp)                :: lwork_tmp
      real(dp)                   :: lrwork_tmp
      integer(int32)                :: N, info, minus_1

      minus_1 = - 1
      N       = size(H, dim = 1)

      call zheevd(vn_flag, 'U', N, H, N, eig_val, &
                  lwork_tmp, minus_1, lrwork_tmp, minus_1, liwork, minus_1, info)

      lwork  = int(lwork_tmp, kind  = 8)
      lrwork = int(lrwork_tmp, kind = 8)
   end subroutine calc_zheevd_size

   subroutine my_section(me, nProcs, length, first, last)
      implicit none
      integer, intent(in)     :: me, nProcs, length
      integer, intent(out)    :: first, last
      integer(int32)                 :: small_chunk, nLarge, i, chunk_sz

      small_chunk =  length / nProcs
      nLarge =  mod(length, nProcs)
      first =  1

      ! MPI Procs start at 0
      do i =  0, me
         if(i <  nLarge) then
            chunk_sz =  small_chunk + 1
         else
            chunk_sz =  small_chunk
         endif

         if(me ==  i) then
            last =  first + chunk_sz - 1
         else
            first =  first +  chunk_sz
         endif
      enddo
   end subroutine my_section

   subroutine sections(nProcs, length, num_elem, offset)
      implicit none
      integer, intent(in)    :: nProcs
      integer, intent(in)    :: length
      integer, intent(out)   :: num_elem(nProcs), offset(nProcs)
      integer(int32)                :: small_chunk, nLarge, i

      small_chunk =  length / nProcs
      nLarge =  mod(length, nProcs)

      ! MPI Procs start at 0
      do i =  1, nProcs
         if(i-1 <  nLarge) then
            num_elem(i) = small_chunk + 1
         else
            num_elem(i) = small_chunk
         endif
      enddo

      offset(1) = 0

      do i = 2,nProcs
         offset(i) =  offset(i-1) + num_elem(i-1)
      enddo
   end subroutine

   subroutine linspace(start, halt, n, x)
      implicit none
      real(dp), intent(in)    :: start, halt
      integer, intent(in) :: n
      real(dp), allocatable   :: x(:)
      real(dp)                :: step, curr
      integer(int32)             :: i

      step =  (halt - start) /  (n-1)
      curr =  start

      if(allocated(x)) then
         if(size(x) /= n) then
            deallocate(x)
         endif
      endif
      if(.not. allocated(x)) then
         allocate(x(n))
      endif

      do i = 1,n
         x(i) =  curr
         curr =  curr +  step
      enddo
   end subroutine linspace

   pure function cross_prod(a,b) result(c)
      implicit none
      real(dp), intent(in)   :: a(3), b(3)
      real(dp)               :: c(3)

      c(1) =  a(2) * b(3) - a(3) * b(2)
      c(2) =  a(3) * b(1) - a(1) * b(3)
      c(3) =  a(1) * b(2) - a(2) * b(1)
   end function cross_prod

   subroutine check_ierr(ierr, me_in, info, msg)
      implicit none
      integer, intent(inout)  :: ierr(:)
      integer, optional       :: me_in
      integer(int32)                 :: me, holder
      integer(int32)                 :: i
      character(len=*), optional :: info, msg(:)

      if(present(me_in)) then
         me = me_in
      else
         call MPI_Comm_rank(MPI_COMM_WORLD, me, holder)
      endif

      do i = 1,size(ierr)
         if(ierr(i) /= 0) then
            if(present(info)) then
               write (*, "(A,I5,A,I3,A)")  "[", me, "]  error at : ", i, info
            else
               write (*, "(A, I5, A, I3)")  "[", me, "] error at : ", i
            endif
            if(present(msg)) then
               write (*,"(I5, A)") me, msg
            endif
            call error_msg("encountered error. Aborting!", abort=.True.)
         endif
      enddo
      ierr =  0
   end subroutine check_ierr

   recursive subroutine qargsort(data, idx, first_in, last_in)
      implicit none
      real(dp)                 :: data(:)
      integer, allocatable :: idx(:)
      integer, optional    :: first_in, last_in
      integer(int32)              :: first, last, i, p

      !setting up a nice interface
      if(.not. present(first_in)) then
         if(allocated(idx)) then
            if(size(data) /= size(idx)) then
               deallocate(idx)
            endif
         endif

         if(.not. allocated(idx)) allocate(idx(size(data)))

         forall(i=1:size(idx)) idx(i) = i
         first = 1
         last  = size(data)
      else
         first = first_in
         last  = last_in
      endif

      !actual algo
      if(first < last) then
         p = partition(data, idx, first, last)
         call qargsort(data, idx, first, p - 1)
         call qargsort(data, idx, p + 1, last)
      endif
   end subroutine qargsort

   function partition(data, idx, first, last) result(p)
      implicit none
      real(dp)                 :: data(:), pivot
      integer, allocatable :: idx(:)
      integer(int32)              :: first, last, p, i, j, tmp

      pivot = data(idx(last))
      i = first - 1
      do j = first, last-1
         if(data(idx(j)) <= pivot) then
            i = i + 1

            !swap
            tmp = idx(i)
            idx(i) = idx(j)
            idx(j) = tmp
         endif
      enddo
      tmp = idx(i+1)
      idx(i+1) = idx(last)
      idx(last) = tmp

      p = i + 1
   end function partition

   function find_list_idx(list, elem) result(idx)
      implicit none
      integer, intent(in)   :: list(:), elem
      integer(int32)               :: idx

      idx = 1
      do while(list(idx) /= elem)
         idx = idx + 1
         if(idx > size(list)) then
            idx =  -1
            return
         endif
      enddo
   end function find_list_idx

   function dag(M) result(dag_M)
      implicit none
      complex(dp), intent(in) :: M(:,:)
      complex(dp)             :: dag_M(size(M,2), size(M,1))
      dag_m = conjg(transpose(M))
   end function

   subroutine gaussian_noise(out_arr, sigma, mu)
      implicit none
      real(dp), intent(out)  :: out_arr(:)
      real(dp), intent(in)   :: sigma, mu
      real(dp)               :: u(size(out_arr)), fac, u_odd(2)
      integer(int32)            :: even_end, i

      call random_number(u)
      even_end = size(out_arr) - mod(size(out_arr),2)

      do i = 1,even_end,2
         fac = sqrt(-2d0 *  log(u(i)))
         out_arr(i)   = fac * cos(2d0 * PI * u(i+1))
         out_arr(i+1) = fac * sin(2d0 * PI * u(i+1))
      enddo

      if(mod(size(out_arr),2) ==  1)then
         call random_number(u_odd)
         fac  = sqrt(-2d0 *  log(u_odd(1)))
         out_arr(size(out_arr))   = fac * cos(2d0 * PI * u(2))
      endif

      out_arr =  sigma * out_arr + mu
   end subroutine gaussian_noise

   function date_time() result(res)
      implicit none
      character(10)  :: time
      character(8)   :: date
      character(23)  :: res

      call date_and_time(time=time, date=date)

      res(1:2) = time(1:2)
      res(3:3) = ":"
      res(4:5) = time(3:4)
      res(6:6) = ":"
      res(7:12) = time(5:10)
      res(13:13) = " "

      res(14:15) =  date(7:8)
      res(16:16) =  "."
      res(17:18) =  date(5:6)
      res(19:19) =  "."
      res(20:23) =  date(1:4)
   end function date_time

   subroutine time_msg(msg)
      implicit none
      character(len=*), intent(in)   :: msg
      integer(int32)                        :: me, info

      call MPI_Comm_rank(MPI_COMM_WORLD, me, info)

      write (*,*) me, date_time() // ": " // msg
   end subroutine time_msg

   subroutine error_msg(msg, p_color, abort)
      implicit none
      character(len=*), intent(in) :: msg
      character(len=*), optional   :: p_color
      character(len=2)             :: code
      integer(int32)                   :: me, holder, info
      logical, optional            :: abort

      if(present(p_color)) then
         code = p_color
      else
         code =  c_red
      endif

      call MPI_Comm_rank(MPI_COMM_WORLD, me, holder)

      if(me == root) then
         write (*,*) c_start // code // c_end // msg // c_clear

         ! abort if necessary
         if(present(abort)) then
            if(abort) then
               call MPI_Abort(MPI_COMM_WORLD, 0, info)
            endif
         endif
      else

         ! if a non-root process is called also abort (a little later)
         if(present(abort)) then
            if(abort) then
               call sleep(5)
               call MPI_Abort(MPI_COMM_WORLD, 0, info)
            endif
         endif
      endif
   end subroutine error_msg

   function mtx_norm(M) result(norm)
      implicit none
      real(xdp), intent(in)     :: M(:,:)
      real(dp)                 :: norm

      norm =  my_norm2(reshape(M, [size(M)]))
   end function
end module Class_helper
