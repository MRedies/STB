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

    function get_unit_conv(field_name, cfg, me, bcast) result(factor)
        implicit none
        character(len=*), intent(in)        :: field_name 
        integer(4), intent(in)              :: me
        logical, optional                   :: bcast
        logical                             :: bcast_loc
        type(CFG_t)                         :: cfg
        real(8)                             :: factor
        integer(4)                          :: ierr
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
                factor =  1.0d0
            case ("eV")
                factor =  0.03674932d0
            case ("a0^-1") 
                factor = 1.0d0
            case ("K")
                factor = 1.0d0
            case default
                write (*,*) "Unit unknown"
                stop
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

    subroutine check_ierr(ierr, me)
        implicit none
        integer(4), intent(in)  :: ierr(:), me
        integer(4)              :: i

        do i = 1,size(ierr)
            if(ierr(i) /= 0) then
                write (*, "(A, I3, A, I3)") &
                    "[", me, "] Bcast error at :", i
                stop
            endif
        enddo
    end subroutine check_ierr

        

end module Class_helper
