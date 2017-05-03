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
            write(p_unit, "(E10.3)") vec(i)
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
                    write (p_unit, "(E10.3)", advance="no") mtx(i,j)
                endif
            enddo
            write(p_unit, "(E10.3)") mtx(i, size(mtx,2))
        enddo

    End Subroutine print_mtx_real

    Subroutine  print_mtx_cmplx(p_unit, mtx)
        Implicit None
        complex(8), dimension(:,:), intent(in)    :: mtx
        integer(4), intent(in)                     :: p_unit
        integer(4)                                 :: i,j
        character(len=4)                           :: i_str
        character(len=3)                           :: deli = "j), "

        do i = 1, size(mtx,1)
            do j =  1, size(mtx,2)
                if(aimag(mtx(i,j))>=0.0d0) then
                    i_str =  " +i*"
                else
                    i_str =  " -i*"
                endif

                if(j < size(mtx,2)) then 
                    write (p_unit, "(SP,a,E10.3,a,E10.3,a)", advance="no") &
                        "(",real(mtx(i,j)), " ", aimag(mtx(i,j)), deli
                endif
            enddo
            write(p_unit, "(SP,a,E10.3,a,E10.3,a)") &
                "(",real(mtx(i, size(mtx,2))), " ", &
                aimag(mtx(i, size(mtx,2))), deli
        enddo
    End Subroutine print_mtx_cmplx
End Module  output
