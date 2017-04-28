module Class_atom
    implicit none

    type atom
        real(8)                     :: m_phi, m_theta
        integer(4), dimension(2)    :: pos
        ! only three neighours to avoid double counting
        ! right, up-right and up-left
        integer(4), dimension(3)    :: neigh
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
        integer(4), dimension(2), intent(in) :: p_pos

        ret%m_phi      = 0.0d0
        ret%m_theta    = 0.0d0
        ret%pos        = p_pos
        ret%neigh      = (/0,0,0/)
    end function

end module 
