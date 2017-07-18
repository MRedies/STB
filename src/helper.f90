module Class_helper
    use m_config
    use Constants
    use mpi
    implicit none

contains
    function cnorm2(vec) result(norm)
        implicit none
        complex(8), intent(in)   :: vec(:)
        real(8)                  :: norm
    
        ! see fortran defition for complex dot_product
        norm = sqrt(dot_product(vec,vec))
    end function cnorm2

    function my_norm2(vec) result(norm)
        implicit none
        real(8), intent(in) :: vec(:)
        real(8)             :: norm

        norm = sqrt(dot_product(vec,vec))
    end function my_norm2

    function omp_matvec(A,x) result(b)
        implicit none
        complex(8), intent(in)    :: A(:,:), x(:)
        complex(8)                :: b(size(x))
        integer(4)                :: i, j

        !$omp parallel do private(j) shared(A,x,b)
        do i = 1,size(x)
            b(i) = 0d0
            do j = 1,size(x)
                b(i) =  b(i) +  A(i,j) * x(j)
            enddo
        enddo
    end function omp_matvec
    
    function matvec(A,x) result(b)
        implicit none
        complex(8), intent(in)    :: A(:,:), x(:)
        complex(8)                :: b(size(x))
        integer(4)                :: i, j

        do i = 1,size(x)
            b(i) = 0d0
            do j = 1,size(x)
                b(i) =  b(i) +  A(i,j) * x(j)
            enddo
        enddo

        !b = 0d0
        !do j = 1,size(x)
            !do i = 1,size(x)
                !b(i) =  b(i) +  A(i,j) * x(j)
            !enddo
        !enddo
    end function matvec

    subroutine calc_zheevd_size(vn_flag, H, eig_val, lwork, lrwork, liwork)
        implicit none
        character(len=1)          :: vn_flag
        complex(8), intent(in)    :: H(:,:)
        real(8), intent(in)       :: eig_val(:)
        integer(4), intent(out)   :: lwork, lrwork, liwork
        complex(8)                :: lwork_tmp
        real(8)                   :: lrwork_tmp
        integer(4)                :: N, info

        N =  size(H, dim=1)
        call zheevd(vn_flag, 'U', N, H, N, eig_val, &
            lwork_tmp, -1, lrwork_tmp, - 1, liwork, -1, info)
        lwork =  int(lwork_tmp)
        lrwork =  int(lrwork_tmp)
    end subroutine calc_zheevd_size

    subroutine my_section(me, nProcs, length, first, last)
        implicit none
        integer(4), intent(in)     :: me, nProcs, length
        integer(4), intent(out)    :: first, last
        integer(4)                 :: small_chunk, nLarge, i, chunk_sz

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
        integer(4), intent(in)    :: nProcs, length
        integer(4), intent(out)   :: num_elem(nProcs), offset(nProcs)
        integer(4)                :: small_chunk, nLarge, i

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
        real(8), intent(in)    :: start, halt
        integer(4), intent(in) :: n
        real(8), allocatable   :: x(:)
        real(8)                :: step, curr 
        integer(4)             :: i

        step =  (halt - start) /  (n-1)
        curr =  start

        if(allocated(x) .and. size(x) /= n) then
            deallocate(x)
        endif
        if(.not. allocated(x)) then 
            allocate(x(n))
        endif

        do i = 1,n
            x(i) =  curr
            curr =  curr +  step
        enddo
    end subroutine linspace

    function cross_prod(a,b) result(c)
        implicit none
        real(8), intent(in)   :: a(3), b(3)
        real(8)               :: c(3)

        c(1) =  a(2) * b(3) - a(3) * b(2)
        c(2) =  a(3) * b(1) - a(1) * b(3)
        c(3) =  a(1) * b(2) - a(2) * b(1)
    end function cross_prod


    subroutine check_ierr(ierr, me, info)
        implicit none
        integer(4), intent(in)     :: ierr(:), me
        integer(4)                 :: error, i
        character(len=*), optional :: info

        do i = 1,size(ierr)
            if(ierr(i) /= 0) then
                if(present(info)) then
                    write (*, "(A, I3, A, I3)")  "[", me, "] Bcast error at :", i, info
                else
                    write (*, "(A, I3, A, I3)")  "[", me, "] Bcast error at :", i
                endif
                call MPI_Barrier(MPI_COMM_WORLD, error)
                call MPI_Abort(MPI_COMM_WORLD, 0, error)
            endif
        enddo
    end subroutine check_ierr

    recursive subroutine qargsort(data, idx, first_in, last_in)
        implicit none
        real(8)                 :: data(:)
        integer(4), allocatable :: idx(:)
        integer(4), optional    :: first_in, last_in
        integer(4)              :: first, last, i, p

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
        real(8)                 :: data(:), pivot
        integer(4), allocatable :: idx(:)
        integer(4)              :: first, last, p, i, j, tmp

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
        integer(4), intent(in)   :: list(:), elem
        integer(4)               :: idx

        idx = 1
        do while(list(idx) /= elem)
            idx = idx + 1
            if(idx > size(list)) then
                idx =  -1
                return
            endif
        enddo
    end function find_list_idx
end module Class_helper
