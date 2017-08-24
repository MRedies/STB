module Class_atom
    implicit none
   
    enum, bind(c)  !> A or B site in graphene
        enumerator :: A_site, B_site, no_site
    end enum 

    type atom
        real(8)               :: m_phi   !> azimuthal spin angle \f$\phi\f$
                                         !> see german wikipedia, not english
        real(8)               :: m_theta !> polar spin angle \f$\theta\f$
                                         !> see german wikipedia, not english                 
        real(8), dimension(3) :: pos     !> Position in RS in atomic units
        integer(4)            :: n_neigh !> number of neighbours
        integer               :: site_type !> A or B site
        
        integer(4), allocatable  :: neigh_idx(:)  !> index of neighbour atom
        integer(4), allocatable  :: snd_neigh_idx_clk(:)  !> index of next ot nearest neighbour atom clock-wise only
        real(8), allocatable     :: neigh_conn(:,:) !> real space connection to neighbour. 
        real(8), allocatable     :: snd_neigh_conn_clk(:,:) !> real space_connection to 2nd nearest neigh 
        !> First index connection, second element of connection.

    contains
        procedure :: get_m_cart => get_m_cart 
        procedure :: set_sphere => set_sphere
        procedure :: set_m_cart => set_m_cart
        procedure :: free_atm   => free_atm
    end type atom
contains
    subroutine free_atm(self)
        implicit none
        class(atom)         :: self
        
        if(allocated(self%neigh_idx)) deallocate(self%neigh_idx)
        if(allocated(self%neigh_conn)) deallocate(self%neigh_conn)
    end subroutine free_atm

    function get_m_cart(self) result(coord)
        implicit none
        Class(atom), intent(in)   :: self
        real(8), dimension(3)     :: coord
        
        ! assume r =  1
        coord(1) = sin(self%m_theta) *  cos(self%m_phi)
        coord(2) = sin(self%m_theta) *  sin(self%m_phi)
        coord(3) = cos(self%m_theta)
    end function get_m_cart
    
    function init_ferro_z(p_pos, site) result(ret)
        implicit none
        type(atom)                 :: ret
        real(8), intent(in)        :: p_pos(3)
        integer, optional          :: site

        if(present(site)) then
            ret%site_type = site
        else
            ret%site_type =  no_site
        endif

        ret%m_phi      = 0d0 
        ret%m_theta    = 0d0
        ret%pos        = p_pos
    end function init_ferro_z

    subroutine set_m_cart(self,x,y,z)
        implicit none
        class(atom)           :: self
        real(8), intent(in)   :: x,y,z

        self%m_theta = acos(z / sqrt(x*x + y*y + z*z))
        self%m_phi   = atan2(y,x)
    end subroutine set_m_cart

    subroutine set_sphere(self, phi, theta)
        implicit none
        class(atom)           :: self
        real(8), intent(in)   :: phi, theta

        self%m_phi   = phi
        self%m_theta = theta
    end subroutine set_sphere 
end module 
