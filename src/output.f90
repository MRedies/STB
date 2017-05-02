Module  output  
    implicit None
    interface print_mtx
        module procedure print_mtx_real, print_mtx_cmplx, print_vec_real
    end interface
contains
    subroutine print_vec_real(p_unit, vec)
        implicit none
        real(8), dimension(:), intent(in)      :: vec
        integer(4), intent(in)                 :: p_unit
        integer(4)                             :: i

        do i =  1, size(vec)
            write(p_unit, "(ES18.10)") vec(i)
        enddo
    end subroutine print_vec_real

    Subroutine  print_mtx_real(p_unit, mtx)
        Implicit None
        real(8), dimension(:,:), intent(in)    :: mtx
        integer(4), intent(in)                 :: p_unit
        integer(4)                             :: i,j
        
        do i = 1, size(mtx,1)
            do j =  1, size(mtx,2)
                if(j < size(mtx,2)) then 
                    write (p_unit, "(ES18.10)", advance="no") mtx(i,j)
                endif
            enddo
            write(p_unit, "(ES18.10)") mtx(i, size(mtx,2))
        enddo

    End Subroutine print_mtx_real

    Subroutine  print_mtx_cmplx(p_unit, mtx)
        Implicit None
        complex(16), dimension(:,:), intent(in)    :: mtx
        integer(4), intent(in)                     :: p_unit
        integer(4)                                 :: i,j
        character(len=3)                           :: i_str

        do i = 1, size(mtx,1)
            do j =  1, size(mtx,2)
                if(aimag(mtx(i,j))>=0.0d0) then
                    i_str =  "+i*"
                else
                    i_str =  "-i*"
                endif

                if(j < size(mtx,2)) then 
                    write (p_unit, "(ES18.10,a,ES18.10)", advance="no") &
                        real(mtx(i,j)), i_str, aimag(mtx(i,j))
                endif
            enddo
            write(p_unit, "(ES28.10,a,ES18.10)") &
                real(mtx(i, size(mtx,2))), i_str, aimag(mtx(i, size(mtx,2)))
        enddo
    End Subroutine print_mtx_cmplx
End Module  output
