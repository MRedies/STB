module Constants
    !< mathematical parameter
    real(8), parameter     :: PI = 3.14159265359d0

    !< complex parameter
    complex(8), parameter :: i_unit = (0d0, 1d0)
    complex(8), parameter :: c_0    = (0d0, 0d0)
    complex(8), parameter :: c_1    = (1d0, 0d0)
    
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
    real(8), parameter  :: small_imag = 1d-24
    real(8), parameter  :: pos_eps    = 1d-6
end module Constants 
