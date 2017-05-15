module Class_atom
    implicit none

    type atom
        real(8)  :: m_phi   !< azimuthal spin angle \f$\phi\f$
        real(8)  :: m_theta !< polar spin angle \f$\theta\f$
        real(8), dimension(3) :: pos !< Position in RS in atomic units
        integer(4)  :: n_neigh !< number of neighbours
        real(8), dimension(:), allocatable     :: hopping !< hopping term for a given connection
        integer(4), dimension(:), allocatable  :: neigh_idx !< index of neighbour atom
        real(8), dimension(:,:), allocatable   :: neigh_conn !< real space connection to neighbour. 
        !< First index connection, second element of connection.

    contains
        procedure :: get_m_cart =>  get_m_cart 
    end type atom
contains
    function get_m_cart(this) result(coord)
        implicit none
        Class(atom), intent(in)   :: this
        real(8), dimension(3)     :: coord
        
        ! assume r =  1
        coord(1) = sin(this%m_theta) *  cos(this%m_phi)
        coord(2) = sin(this%m_theta) *  sin(this%m_phi)
        coord(3) = cos(this%m_theta)
    end function get_m_cart
    
    function init_ferro(p_pos) result(ret)
        implicit none
        type(atom)                           :: ret
        real(8), dimension(3), intent(in) :: p_pos

        ret%m_phi      = 0.0d0
        ret%m_theta    = 0.0d0
        ret%pos        = p_pos
    end function

end module 
