module Class_math_helper
    implicit none

contains
    function CNORM2(vec) result(norm)
        implicit none
        complex(8), intent(in)   :: vec(:)
        real(8)                  :: norm
        integer(4)               :: i

        norm =  0d0 
        do i =  1,size(vec)
            norm =  norm +  vec(i) * conjg(vec(i))
        enddo
    end function CNORM2

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
            if(i <= nLarge) then
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
            if(i <= nLarge) then
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

end module Class_math_helper
