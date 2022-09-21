module Class_atom
    use Class_helper
    use Constants
    use mpi_f08
    implicit none

    enum, bind(c)  !> A or B site in graphene
        enumerator :: A_site, B_site, no_site
    end enum

    type atom
        real(8)                  :: m_phi   !> azimuthal spin angle \f$\phi\f$
                                            !> see german wikipedia, not english
        real(8)                  :: m_theta !> polar spin angle \f$\theta\f$
                                            !> see german wikipedia, not english
        real(8), dimension(3)    :: pos     !> Position in RS in atomic units
        integer(8)               :: site_type !> A or B site

        integer   , allocatable  :: neigh_idx(:)  !> index of neighbour atom
        real(8), allocatable     :: neigh_conn(:,:) !> real space connection to neighbour.
        integer(4), allocatable  :: conn_type(:) !> type of connection
        !> First index connection, second element of connection.

        integer                  :: me, nProcs

    contains
        procedure :: get_m_cart      => get_m_cart
        procedure :: set_sphere      => set_sphere
        procedure :: set_m_cart      => set_m_cart
        procedure :: free_atm        => free_atm
        procedure :: compare_to_root => compare_to_root
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

    function init_ferro_z(p_pos, comm, site) result(self)
        implicit none
        type(atom)                 :: self
        real(8), intent(in)        :: p_pos(3)
        type(MPI_Comm), intent(in) :: comm
        integer(8), optional       :: site
        integer                    :: ierr(2)

        call MPI_Comm_size(comm, self%nProcs, ierr(1))
        call MPI_Comm_rank(comm, self%me, ierr(2))
        call check_ierr(ierr, self%me, "init_ferro_z call")

        if(present(site)) then
            self%site_type = site
        else
            self%site_type =  no_site
        endif

        self%m_phi      = 0d0
        self%m_theta    = 0d0
        self%pos        = p_pos
    end function init_ferro_z

    subroutine set_m_cart(self,x,y,z)
        implicit none
        class(atom)           :: self
        real(8), intent(in)   :: x,y,z

        if(abs(my_norm2([x,y,z]) -  1d0) > 1d-4) then
            call error_msg("Spin not normed", abort=.True.)
        endif
        self%m_theta = acos(z / sqrt(x*x + y*y + z*z))
        self%m_phi   = atan2(y,x)
    end subroutine set_m_cart

    subroutine set_sphere(self, phi, theta)
        implicit none
        class(atom)           :: self
        real(8), intent(in)   :: phi, theta

        self%m_phi   = phi
        self%m_theta = theta
        !write(*,*) "self%m_theta",self%m_theta
    end subroutine set_sphere

    function compare_to_root(self,comm) result(success)
        implicit none
        class(atom)                :: self
        real(8)                    :: tmp, tmp_p(3)
        integer                    :: ierr(10), tmp_i,s1,s2
        type(MPI_Comm), intent(in) :: comm
        integer, allocatable    :: tmp_ivec(:)
        integer(8)              :: tmp_i8
        integer(4), allocatable :: tmp_i4vec(:)
        real(8), allocatable    :: tmp_rmtx(:,:)
        logical                 :: success

        success = .True.

        ! compare angles
        if(self%me == root) tmp = self%m_phi
        call MPI_Bcast(tmp, 1, MPI_REAL8, root, comm, ierr(1))
        if(abs(tmp - self%m_phi) > 1d-12) then
            call error_msg("m_phi doesn't match", abort=.True.)
            success = .False.
        endif

        if(self%me == root) tmp = self%m_theta
        call MPI_Bcast(tmp, 1, MPI_REAL8, root, comm, ierr(2))
        if(abs(tmp - self%m_theta) > 1d-12) then
            call error_msg("m_theta doesn't match", abort=.True.)
            success = .False.
        endif

        ! compare positions
        if(self%me == root) tmp_p = self%pos
        call MPI_Bcast(tmp_p, 3, MPI_REAL8, root, comm, ierr(3))
        if(my_norm2(tmp_p - self%pos) > 1d-11) then
            call error_msg("pos doesn't match", abort=.True.)
            success = .False.
        endif

        ! compare site_types
        if(self%me == root) tmp_i8 = self%site_type
        call MPI_Bcast(tmp_i8, 1, MPI_INTEGER8, root, comm, ierr(4))
        if(tmp_i8 /= self%site_type) then
            call error_msg("site_type doesn't match", abort=.True.)
            success = .False.
        endif

        ! compare neighbours
        if(self%me == root) tmp_i = size(self%neigh_idx)
        call MPI_Bcast(tmp_i, 1, MPI_INTEGER, root, comm, ierr(5))
        if(tmp_i /= size(self%neigh_idx)) then
            call error_msg("size(neigh_idx) doesn't match", abort=.True.)
            success = .False.
        endif

        allocate(tmp_ivec(size(self%neigh_idx)))
        if(self%me == root) tmp_ivec = self%neigh_idx
        call MPI_Bcast(tmp_ivec(:), size(self%neigh_idx), MPI_INTEGER, &
                       root, comm, ierr(6))
        if(any(tmp_ivec /= self%neigh_idx)) then
            write (*,*) self%me, "neigh_idx", self%neigh_idx
            write (*,*) self%me, "tmp_ivec", tmp_ivec
            call error_msg("neigh_idx doesn't match", abort=.True.)
            success = .False.
        endif
        deallocate(tmp_ivec)

        if(self%me == root) tmp_i = size(self%neigh_conn)
        call MPI_Bcast(tmp_i, 1, MPI_INTEGER, root, comm, ierr(7))
        if(tmp_i /= size(self%neigh_conn)) then
            call error_msg("size(neigh_conn) doesn't match", abort=.True.)
            success = .False.
        endif
        s1 = size(self%neigh_conn, dim=1)
        s2 = size(self%neigh_conn, dim=2)
        allocate(tmp_rmtx(s1, s2))
        if(self%me == root) tmp_rmtx = self%neigh_conn
        call MPI_Bcast(tmp_rmtx(1:s1,1:s2), size(tmp_rmtx), MPI_REAL8, &
                                      root, comm, ierr(8))
        if(mtx_norm(tmp_rmtx - self%neigh_conn) >  1d-11) then
            call error_msg("neigh_conn doesn't match", abort=.True.)
            success = .False.
        endif
        deallocate(tmp_rmtx)

        ! compare conn types
        if(self%me == root) tmp_i = size(self%conn_type)
        call MPI_Bcast(tmp_i, 1, MPI_INTEGER, root, comm, ierr(9))
        if(tmp_i /= size(self%conn_type)) then
            call error_msg("size(conn_type) doesn't match", abort=.True.)
            success = .False.
        endif

        allocate(tmp_i4vec(size(self%conn_type)))
        if(self%me == root) tmp_i4vec = self%conn_type
        call MPI_Bcast(tmp_i4vec(:), size(tmp_i4vec), MPI_INTEGER4, &
                       root, comm, ierr(10))
        if(any(tmp_i4vec /= self%conn_type)) then
            call error_msg("conn_type doesn't match", abort=.True.)
            success = .False.
        endif
        deallocate(tmp_i4vec)

        if(.not. success) then
            call error_msg("Atom test failed", abort=.True.)
        endif

        call check_ierr(ierr, self%me, "compare atoms")
    end function compare_to_root
end module
