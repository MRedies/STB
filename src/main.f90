program STB
    use Class_unit_cell
    implicit none
    integer(4)      :: num, i
    type(unit_cell) :: UC
    type(atom), dimension(:), allocatable   :: atoms 

    !call UC%init(200)
    !atoms =  UC%get_atoms()
    !write (*,*) "Index    ", "connected"
    !do i =  1, UC%get_num_atoms()
        !write (*,*) i, atoms(i)%neigh
    !enddo
    call test_mkl()
contains
    Subroutine  test_mkl()
    Implicit None
    external DGESV
    external PRINT_MATRIX, PRINT_INT_VECTOR
    integer(4), parameter         :: N =  5, NRHS = 3
    integer(4), parameter         :: LDA =  N, LDB = N
    integer(4)                    :: info
    integer(4), dimension(N)      :: IPIV
    real(8), dimension(LDA, N)    :: A
    real(8), dimension(LDB, NRHS) :: B

    A =  reshape((/6.80,-2.11, 5.66, 5.97, 8.23, &
                  -6.05,-3.30, 5.36,-4.44, 1.08, &
                  -0.45, 2.58,-2.70, 0.27, 9.04, &
                   8.32, 2.71, 4.35,-7.17, 2.14, &
                   -9.67,-5.14,-7.26, 6.08,-6.87 /), shape(A))
    
    B =  reshape((/4.02, 6.19,-8.22,-7.57,-3.03, &
                  -1.56, 4.00,-8.67, 1.75, 2.86, &
                   9.81,-4.09,-4.57,-8.61, 8.99 /), shape(B))

    call DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)

    write (*,*) "Info", info

    End Subroutine test_mkl
end program STB

