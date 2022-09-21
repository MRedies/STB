module Constants
   !< mathematical parameter
   real(dp), parameter     :: PI = 3.14159265359d0

   !< complex parameter
   complex(dp), parameter :: i_unit = (0d0, 1d0)
   complex(dp), parameter :: c_0    = (0d0, 0d0)
   complex(dp), parameter :: c_1    = (1d0, 0d0)
   complex(dp), parameter :: c_i    = (0d0, 1d0)

   !< pauli matricies
   complex(dp), parameter :: sigma_x(2,2) &
                            =  reshape((/0,1,&
                                         1,0/), (/2,2/))
   complex(dp), parameter :: sigma_y(2,2) &
                            =  reshape((/c_0,   i_unit,&
                                         -i_unit, c_0/), (/2,2/))
   complex(dp), parameter :: sigma_z(2,2) &
                            =  reshape((/1, 0,&
                                         0,-1/), (/2,2/))
   !< Physical constants
   real(dp), parameter    :: boltzmann_const = 3.16681d-6 !E_h/K (hartrees per kelvin)

   !< MPI Parameter
   integer(int64), parameter :: root = 0

   !<
   !real(dp), parameter  :: small_imag = 1d-24
   real(dp), parameter  :: eta    =  1d-3!1d0/3d0*1d-4!10 meV
   real(dp), parameter  :: eta_sq =  eta**2
   real(dp), parameter  :: pos_eps    = 1d-6

   !> some useful angles
   real(dp), parameter               :: deg_30 =  30.0 * PI / 180.0
   real(dp), parameter               :: deg_60 =  60.0 * PI / 180.0

   !> physical constants in atomic units hartree
   real(dp), parameter               :: speed_of_light =  137.035999d0 !> 1/(fine structure const)

   !> Lx
   complex(dp), parameter :: Lx(3,3) = transpose(reshape([c_0, c_0,   c_0,&
                                                         c_0, c_0, - c_i,&
                                                         c_0, c_i,   c_0], [3,3]))

   complex(dp), parameter :: Ly(3,3) = transpose(reshape([  c_0, c_0, c_i,&
                                                         c_0, c_0, c_0,&
                                                         - c_i, c_0, c_0], [3,3]))

   complex(dp), parameter :: Lz(3,3) = transpose(reshape([c_0, - c_i, c_0,&
                                                         c_i,   c_0, c_0,&
                                                         c_0,   c_0, c_0], [3,3]))

   complex(dp), parameter :: LxpILy(3,3) = Lx + i_unit * Ly
   complex(dp), parameter :: LxmILy(3,3) = Lx - i_unit * Ly

   complex(dp), parameter :: sqrt_2   = (0.70710678118d0, 0d0)
   complex(dp), parameter :: i_sqrt_2 = (0d0,             0.70710678118d0)

   complex(dp), parameter :: BT_cmplx_to_real(3,3) &
                            = transpose(reshape([sqrt_2,     c_0, - sqrt_2,   &
                                                 - i_sqrt_2, c_0, - i_sqrt_2, &
                                                 c_0,        c_1, c_0],       [3, 3]))

   complex(dp), parameter :: BT_real_to_cmplx(3,3) &
                            = transpose(reshape([sqrt_2,   i_sqrt_2, c_0,  &
                                                 c_0,      c_0,      c_1,  &
                                                 - sqrt_2, i_sqrt_2, c_0], [3, 3]))
end module Constants
