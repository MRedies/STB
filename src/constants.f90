module Constants
    !< mathematical parameter
    real(8), parameter     :: PI = 3.14159265359d0

    !< complex parameter
    complex(8), parameter :: i_unit = (0d0, 1d0)
    complex(8), parameter :: c_0    = (0d0, 0d0)
    complex(8), parameter :: c_1    = (1d0, 0d0)
    complex(8), parameter :: c_i    = (0d0, 1d0)

    !< pauli matricies
    complex(8), parameter :: sigma_x(2,2) &
           =  reshape((/0,1,&
                        1,0/), (/2,2/))
    complex(8), parameter :: sigma_y(2,2) &
           =  reshape((/c_0,   i_unit,&
                      -i_unit, c_0/), (/2,2/))
    complex(8), parameter :: sigma_z(2,2) &
           =  reshape((/1, 0,&
                        0,-1/), (/2,2/))
    !< Physical constants
    real(8), parameter    :: boltzmann_const = 3.16681d-6 !E_h/K (hartrees per kelvin)

    !< MPI Parameter
    integer(4), parameter :: root = 0

    !<
    !real(8), parameter  :: small_imag = 1d-24
    real(8), parameter  :: eta_sq =  1d-16
    real(8), parameter  :: eta    =  1d-8
    real(8), parameter  :: pos_eps    = 1d-6

    !> some useful angles
    real(8), parameter               :: deg_30 =  30.0 * PI / 180.0
    real(8), parameter               :: deg_60 =  60.0 * PI / 180.0

    !> physical constants in atomic units hartree
    real(8), parameter               :: speed_of_light =  137.035999d0 !> 1/(fine structure const)


    !> Lx
    complex(8), parameter :: Lx(3,3) = transpose(reshape([c_0, c_0,   c_0,&
                                                          c_0, c_0, - c_i,&
                                                          c_0, c_i,   c_0], [3,3]))

    complex(8), parameter :: Ly(3,3) = transpose(reshape([  c_0, c_0, c_i,&
                                                            c_0, c_0, c_0,&
                                                          - c_i, c_0, c_0], [3,3]))

    complex(8), parameter :: Lz(3,3) = transpose(reshape([c_0, - c_i, c_0,&
                                                          c_i,   c_0, c_0,&
                                                          c_0,   c_0, c_0], [3,3]))

    complex(8), parameter :: LxpILy(3,3) = Lx + i_unit * Ly
    complex(8), parameter :: LxmILy(3,3) = Lx - i_unit * Ly

    complex(8), parameter :: sqrt_2   = (0.70710678118d0, 0d0)
    complex(8), parameter :: i_sqrt_2 = (0d0,             0.70710678118d0)

    complex(8), parameter :: BT_cmplx_to_real(3,3) &
                           = transpose(reshape([sqrt_2,     c_0, - sqrt_2,   &
                                                - i_sqrt_2, c_0, - i_sqrt_2, &
                                                c_0,        c_1, c_0],       [3, 3]))

    complex(8), parameter :: BT_real_to_cmplx(3,3) &
                           = transpose(reshape([sqrt_2,   i_sqrt_2, c_0,  &
                                                c_0,      c_0,      c_1,  &
                                                - sqrt_2, i_sqrt_2, c_0], [3, 3]))
end module Constants
