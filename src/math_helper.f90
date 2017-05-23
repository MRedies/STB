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
end module Class_math_helper
