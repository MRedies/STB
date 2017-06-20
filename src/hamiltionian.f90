module Class_hamiltionian
    use m_config
    use output
    use Class_unit_cell
    use m_npy
    use mpi
    implicit none
    
    type hamil
        real(8)         :: E_s !> onsite eigenenergy
        real(8)         :: t_nn !> nearest neighbour hopping
        real(8)         :: I !> stoner parameter
        complex(8), allocatable    :: del_H(:,:)
        integer(4)      :: nProcs
        integer(4)      :: me
        type(unit_cell) :: UC !> unit cell
        type(units)     :: units
    contains

        procedure :: Bcast_hamil            => Bcast_hamil
        procedure :: setup_H                => setup_H
        procedure :: calc_eigenvalues       => calc_eigenvalues
        procedure :: calc_single_eigenvalue => calc_single_eigenvalue
        procedure :: set_EigenE             => set_EigenE
        procedure :: set_hopping            => set_hopping
        procedure :: set_Stoner             => set_Stoner
        procedure :: setup_Stoner_mtx       => setup_Stoner_mtx
        procedure :: set_derivative_k       => set_derivative_k
        procedure :: set_deriv_FD           =>  set_deriv_FD
        procedure :: calc_deriv_elem        => calc_deriv_elem
        procedure :: calc_berry_tensor_elem => calc_berry_tensor_elem
        procedure :: calc_berry_z           => calc_berry_z
        procedure :: compare_derivative     => compare_derivative
    end type hamil

contains
    subroutine compare_derivative(self, k)
        implicit none 
        class(hamil)                :: self
        real(8), intent(in)         :: k(3)
        complex(8), allocatable     :: fd_H(:,:)
        integer(4)                  :: N, k_idx 
        
        N = 2 * self%UC%num_atoms
        allocate(fd_H(N,N))
        !allocate(self%del_H(N,N))
       
        do k_idx =  1,2
            call self%set_deriv_FD(k, k_idx, fd_H)
            call self%set_derivative_k(k, k_idx)
            
            if(cnorm2(reshape(fd_H - self%del_H, [N*N])) >= 1d-8) then
                write (*,*) "mist"
                stop
            else
                write (*,*) "Cool"
            endif

        enddo
        deallocate(fd_H)
        !deallocate(self%del_H)
    end subroutine compare_derivative

    subroutine setup_H(self,k,H)
        implicit none
        class(hamil)        :: self
        real(8), intent(in) :: k(3)
        complex(8), intent(inout) :: H(:,:)

        if(k(3) /= 0d0) then
            write (*,*) "K_z is non-zero. Abort.", k
            stop
        endif
       
        H =  0d0

        call self%set_EigenE(H)
        call self%set_hopping(k,H)
        call self%set_Stoner(H)
    end subroutine setup_H

    function init_hamil(cfg) result(self)
        implicit none
        type(CFG_t)    :: cfg
        type(hamil)    :: self
        real(8)        :: tmp
        integer(4)     :: ierr
        
        call MPI_Comm_size(MPI_COMM_WORLD, self%nProcs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)
        
        self%units = init_units(cfg, self%me)
        self%UC    = init_unit(cfg)

        if(self%me ==  0) then 
            call CFG_get(cfg, "hamil%E_s", tmp)
            self%E_s =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%t_nn", tmp)
            self%t_nn =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%I", tmp)
            self%I =  tmp * self%units%energy
        endif
        call self%Bcast_hamil()
    end function init_hamil

    subroutine Bcast_hamil(self)
        implicit none
        class(hamil)          :: self
        integer(4), parameter :: num_cast =  3
        integer(4)            :: ierr(num_cast)

        call MPI_Bcast(self%E_s,  1, MPI_REAL8, root, MPI_COMM_WORLD, ierr(1))
        call MPI_Bcast(self%t_nn, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr(2))
        call MPI_Bcast(self%I,    1, MPI_REAL8, root, MPI_COMM_WORLD, ierr(3))
        call check_ierr(ierr, self%me, "Hamiltionian check err")
    end subroutine

    subroutine set_Stoner(self,H)
        implicit none
        class(hamil), intent(in)   :: self
        complex(8), intent(inout)  :: H(:,:)
        complex(8)                 :: S(2,2) !> Stonermatrix
        integer(4)                 :: i, i_up, i_dw

        do i =  1,self%UC%num_atoms
            i_up =  i
            i_dw =  i +  self%UC%num_atoms 

            call self%setup_Stoner_mtx(i,S)

            H(i_up,i_up) = H(i_up, i_up) +  S(1,1)
            H(i_up,i_dw) = H(i_up, i_dw) +  S(1,2)
            H(i_dw,i_up) = H(i_dw, i_up) +  S(2,1)
            H(i_dw,i_dw) = H(i_dw, i_dw) +  S(2,2)
        enddo

    end subroutine set_Stoner

    subroutine setup_Stoner_mtx(self,i,S)
        implicit none
        class(hamil), intent(in) :: self
        integer(4), intent(in)   :: i
        complex(8), intent(inout):: S(2,2)
        real(8)                  :: m(3), fac
        
        m = self%UC%atoms(i)%get_m_cart()
        fac =  - 0.5d0 *  self%I
        
        S = fac * ( m(1) * sigma_x &
                  + m(2) * sigma_y &
                  + m(3) * sigma_z)
    end subroutine setup_Stoner_mtx

    subroutine set_EigenE(self,H)
        implicit none
        class(hamil), intent(in) :: self
        complex(8), intent(inout):: H(:,:)
        integer(4) :: i

        do i =  1,size(H,dim=1)
            H(i,i) = H(i,i) + self%E_s
        enddo
    end subroutine set_EigenE 

    subroutine set_hopping(self,k, H)
        implicit none
        class(hamil), intent(in)          :: self 
        real(8), intent(in)               :: k(3)
        complex(8), intent(inout)         :: H(:,:)
        integer(4)                        :: i, i_d, j,&
                                             j_d, conn
        real(8)                           :: k_dot_r
        complex(8)                        :: new

        ! Spin up
        do i = 1,self%UC%num_atoms
            do conn = 1,self%UC%atoms(i)%n_neigh
                j =  self%UC%atoms(i)%neigh_idx(conn)
                k_dot_r =  dot_product(k, self%UC%atoms(i)%neigh_conn(conn,:))
                
                new = exp(i_unit * k_dot_r) * self%t_nn
                H(i,j) =  H(i,j) + new
                H(j,i) =  H(j,i) + conjg(new)
            enddo
        enddo

        ! Spin down
        do i = 1,self%UC%num_atoms
            i_d =  i + self%UC%num_atoms
            do conn = 1,self%UC%atoms(i)%n_neigh

                j      = self%UC%atoms(i)%neigh_idx(conn)
                j_d = j + self%UC%num_atoms
                
                k_dot_r =  dot_product(k, self%UC%atoms(i)%neigh_conn(conn,:))
                new =  exp(i_unit *  k_dot_r) * self%t_nn
                H(i_d, j_d) = H(i_d, j_d) + new
                H(j_d, i_d) = H(j_d, i_d) + conjg(new)
            enddo
        enddo

    end subroutine set_hopping


    subroutine set_deriv_FD(self, k, k_idx, del_H)
        implicit none
        class(hamil), intent(in) :: self
        real(8), intent(in)      :: k(3)
        integer(4), intent(in)   :: k_idx
        complex(8), allocatable     :: H_forw(:,:), H_back(:,:), del_H(:,:)
        real(8) :: k_forw(3), k_back(3)
        real(8), parameter :: delta_k =  1d-6
        integer(4)         :: N

        N = 2 * self%UC%num_atoms
        allocate(H_back(N,N))
        allocate(H_forw(N,N))
        
        if(k(3) /= 0) then
            write (*,*) "K_z not zero in set_derivative_k"
            stop
        endif
    
        del_H = 0d0
        k_forw = k
        k_back = k
        k_forw(k_idx) = k_forw(k_idx) + 0.5d0 * delta_k
        k_back(k_idx) = k_back(k_idx) - 0.5d0 * delta_k

        call self%setup_H(k_forw, H_forw)
        call self%setup_H(k_back, H_back)
        del_H =  (H_forw - H_back) / delta_k

        deallocate(H_back)
        deallocate(H_forw)
    end subroutine set_deriv_FD

    subroutine set_derivative_k(self, k, k_idx)
        implicit none
        class(hamil)              :: self
        integer(4), intent(in)    :: k_idx
        real(8), intent(in)       :: k(3)
        real(8)                   :: r(3), k_dot_r
        complex(8)                :: forw, back
        integer(4)                :: i, j, conn, i_d, j_d

        if(k(3) /= 0) then
            write (*,*) "K_z not zero in set_derivative_k"
            stop
        endif
    
        self%del_H = 0d0
        do i = 1,self%UC%num_atoms
            i_d =  i + self%UC%num_atoms
            do conn = 1,self%UC%atoms(i)%n_neigh

                j   = self%UC%atoms(i)%neigh_idx(conn)
                j_d = j + self%UC%num_atoms
                r   = self%UC%atoms(i)%neigh_conn(conn,:)

                k_dot_r = dot_product(k, r)
                forw    = i_unit * r(k_idx) * self%t_nn &
                          * exp(i_unit * k_dot_r)
                
                r       = - r
                k_dot_r = - k_dot_r
                back    = i_unit * r(k_idx) * self%t_nn &
                          * exp(i_unit * k_dot_r)

                !Spin up
                self%del_H(i,j)     = self%del_H(i,j) + forw 
                self%del_H(j,i)     = self%del_H(j,i) + back
                !Spin down
                self%del_H(i_d,j_d) = self%del_H(i_d,j_d) + forw
                self%del_H(j_d,i_d) = self%del_H(j_d,i_d) + back
            enddo
        enddo

    end subroutine set_derivative_k

    function calc_deriv_elem(self, psi_nk, psi_mk, k, k_idx) result(elem)
        implicit none
        class(hamil)               :: self
        complex(8), intent(in)     :: psi_nk(:), psi_mk(:)
        real(8), intent(in)        :: k(3)
        integer(4), intent(in)     :: k_idx
        complex(8)                 :: elem
        !complex(8), allocatable    :: tmp_vec(:)
        integer(4)                 :: n

        n =  size(psi_nk)
        
        call self%set_derivative_k(k, k_idx)
        
        !verkackte zgemv does not work don't even think about it
        !tmp_vec =  matmul(self%del_H, psi_mk)

        ! citing intel ref here:
        ! If vector_a is of type complex, the result 
        ! value is SUM (CONJG ( vector_a)* vector_b).
        elem =  dot_product(psi_nk, matmul(self%del_H, psi_mk))

        
    end function

    subroutine calc_berry_tensor_elem(self, k_i, k_j, k, omega)
        implicit none 
        class(hamil)                       :: self
        integer(4), intent(in)             :: k_i, k_j
        real(8), intent(in)                :: k(3)
        real(8), allocatable               :: omega(:) !> \f$ \Omega_{ij}^n\f$
        real(8), allocatable               :: eig_val(:), rwork(:)
        complex(8), allocatable            :: H(:,:), work(:)
        complex(8)                         :: summe, term
        complex(8)  :: fac
        integer(4), allocatable :: iwork(:) 
        integer(4)  :: n_dim, n, m, info, lwork, lrwork, liwork

        n_dim = 2 * self%UC%num_atoms
        allocate(H(n_dim,n_dim))
        !allocate(self%del_H(n_dim,n_dim))
        allocate(eig_val(n_dim))
        
        if(.not. allocated(omega))then
            allocate(omega(n_dim))
        endif

        H = 0d0
        call self%setup_H(k, H)
        lwork  = 4 * n_dim**2 +  2*n_dim 
        lrwork = 8 * n_dim**2 +  5*n_dim + 1
        liwork = 10* n_dim +  3
        allocate(work(lwork))
        allocate(rwork(lrwork))
        allocate(iwork(liwork))
        call zheevd('V', 'L', n_dim, H, n_dim, eig_val, &
                    work, lwork, rwork, lrwork, iwork, liwork, info)
        if(info /= 0) then
            write (*,*) "ZHEEVD in berry calculation failed"
        endif
        do n = 1,n_dim
            summe = 0d0
            do m = 1,n_dim
                if(n /= m) then
                    fac = 1d0 / ((eig_val(n) - eig_val(m))**2 + i_unit * small_imag)

                    term =   self%calc_deriv_elem(H(:,n), H(:,m), k, k_i) &
                           * self%calc_deriv_elem(H(:,m), H(:,n), k, k_j)

                    summe =  summe +  fac * term
                endif
            enddo
            omega(n) = -2d0 * aimag(summe)
        enddo

        deallocate(eig_val)
        deallocate(H)
        !deallocate(self%del_H)
    end subroutine calc_berry_tensor_elem

    subroutine calc_berry_z(self,k, z_comp)
        implicit none
        class(hamil), intent(in)            :: self
        real(8), intent(in)                 :: k(3)
        real(8), allocatable                :: z_comp(:) !> \f$ \Omega^n_z \f$
        real(8), allocatable                :: tmp(:)

        call self%calc_berry_tensor_elem(1,2,k, tmp)
        
        if(.not. allocated(z_comp)) then
            allocate(z_comp(size(tmp)))
        endif

        ! 0.5 *  sigma_xy
        z_comp = 0.5d0 *  tmp

        ! -0.5 *  sigma_yx
        call self%calc_berry_tensor_elem(2,1,k,tmp)
        z_comp = z_comp - 0.5d0 * tmp

        deallocate(tmp)
    end subroutine calc_berry_z

    Subroutine  calc_eigenvalues(self, k_list, eig_val)
        Implicit None
        class(hamil)                      :: self
        real(8), intent(in)               :: k_list(:,:)
        real(8), allocatable,intent(out)  :: eig_val(:,:)
        real(8)                           :: k(3)
        complex(8), allocatable           :: H(:,:)
        integer(4) :: i, N, LWMAX, info
        real(8), allocatable              :: RWORK(:)
        complex(8), allocatable           :: WORK(:)
        integer(4), allocatable           :: IWORK(:)
        N =  2 * self%UC%num_atoms
        LWMAX =  10*N
        allocate(eig_val(N, size(k_list, 2)))
        allocate(H(N,N))
        allocate(RWORK(LWMAX))
        allocate(IWORK(LWMAX))
        allocate(WORK(LWMAX))
        
        do i = 1,size(k_list,2)
            k =  k_list(:,i)
            call self%setup_H(k, H)
            
            call zheevd('N', 'U', N, H, N, eig_val(:,i), WORK, LWMAX, &
                                RWORK, LWMAX, IWORK, LWMAX, info)
            if( info /= 0) then
                write (*,*) "ZHEEV failed: ", info
                stop
            endif
        enddo

        deallocate(H)
        deallocate(RWORK)
        deallocate(WORK)
    End Subroutine calc_eigenvalues

    subroutine calc_single_eigenvalue(self, k, eig_val)
        implicit none
        class(hamil)                      :: self
        real(8), intent(in)               :: k(3)
        real(8), allocatable, intent(out) :: eig_val(:)
        complex(8), allocatable           :: H(:,:), work(:)
        real(8), allocatable              :: rwork(:)
        integer(4), allocatable           :: iwork(:)
        integer(4)                        :: N, LWMAX, info

        N = 2 * self%UC%num_atoms
        LWMAX = 10 * N
        allocate(eig_val(N))
        allocate(H(N,N))
        allocate(work(LWMAX))
        allocate(iwork(LWMAX))
        allocate(rwork(LWMAX))

        call self%setup_H(k, H)
        call zheevd('N', 'U', N, H, N, eig_val, work, LWMAX, &
                     RWORK, LWMAX, IWORK, LWMAX, info)
        
        if(info /=  0) write (*,*) "ZHEEVD failed:", info

        deallocate(H)
        deallocate(work)
        deallocate(iwork)
        deallocate(rwork)
    end subroutine calc_single_eigenvalue

end module
