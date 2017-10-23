module Class_hamiltionian
    use m_config
    use output
    use Class_unit_cell
    use m_npy
    use mpi
    implicit none

    type hamil
        real(8)         :: E_s !> onsite eigenenergy
        real(8)         :: E_A, E_B !> onsite energies for A and B sites in honeycomb
        real(8)         :: Vss_sig !> nearest neighbour hopping
        real(8)         :: Vsp_sig, Vpp_pi, Vpp_sig
        real(8)         :: E_soc
        real(8)         :: t_2 !> amplitude for 2nd nearest neighbour hopping
        real(8)         :: phi_2 !> polar angle of 2nd nearest neighbour hopping in rad
        real(8)         :: t_so !> Rashba spin orb
        real(8)         :: lambda !> local exchange
        real(8)         :: lambda_nl !> non-local exchange (not implemented yet)
        complex(8), allocatable    :: del_H(:,:)
        integer(4)      :: nProcs
        integer(4)      :: me
        integer(4)      :: num_orb, num_up
        type(unit_cell) :: UC !> unit cell
        type(units)     :: units
    contains

        procedure :: Bcast_hamil                => Bcast_hamil
        procedure :: setup_H                    => setup_H
        procedure :: calc_eigenvalues           => calc_eigenvalues
        procedure :: calc_single_eigenvalue     => calc_single_eigenvalue
        procedure :: set_EigenE                 => set_EigenE
        procedure :: set_hopping                => set_hopping
        procedure :: set_snd_hopping            => set_snd_hopping
        procedure :: set_loc_exch               => set_loc_exch
        procedure :: setup_Stoner_mtx           => setup_Stoner_mtx
        procedure :: set_derivative_k           => set_derivative_k
        procedure :: set_rashba_SO              => set_rashba_SO
        procedure :: set_deriv_FD               => set_deriv_FD
        procedure :: calc_berry_z               => calc_berry_z
        procedure :: calc_velo_mtx              => calc_velo_mtx
        procedure :: calc_eig_and_velo          => calc_eig_and_velo
        procedure :: compare_derivative         => compare_derivative
        procedure :: set_derivative_hopping     => set_derivative_hopping 
        procedure :: set_derivative_snd_hopping => set_derivative_snd_hopping
        procedure :: set_derivative_rashba_so   => set_derivative_rashba_so
        procedure :: set_hopp_mtx               => set_hopp_mtx
        procedure :: free_ham                   => free_ham
    end type hamil

contains
    subroutine free_ham(self)
        implicit none
    class(hamil)     :: self

        if(allocated(self%del_H)) deallocate(self%del_H)
        call self%UC%free_uc()
    end subroutine free_ham

    subroutine compare_derivative(self, k)
        implicit none 
    class(hamil)                :: self
        real(8), intent(in)         :: k(3)
        complex(8), allocatable     :: fd_H(:,:)
        integer(4)                  :: N, k_idx, info

        N = 2 * self%num_up
        allocate(fd_H(N,N))
        if(.not. allocated(self%del_H)) then
            allocate(self%del_H(N,N))
        endif

        do k_idx =  1,2
            call self%set_deriv_FD(k, k_idx, fd_H)
            call self%set_derivative_k(k, k_idx)
            
            call test_herm(self%del_H, tag="del_H")
            call test_herm(fd_H, tag="FD H")

            if(cnorm2(reshape(fd_H - self%del_H, [N*N])) >= 1d-8) then
                write (*,*) "FD comp failed"
                write (*,*) "FD"
                call save_npy("output/dbg/fd_H.npy",fd_H) 
                !call print_mtx(fd_H)
                write (*,*) "analytic"
                !call print_mtx(self%del_H)
                call save_npy("output/dbg/del_H.npy", self%del_H) 
                !write (*,*) "diff"
                !call print_mtx(self%del_H -  fd_H)
                call MPI_Abort(MPI_COMM_WORLD, 0, info)
            else
                write (*,*) "great"
            endif
            
        enddo
        deallocate(fd_H)
    end subroutine compare_derivative

    subroutine setup_H(self,k,H)
        implicit none
        class(hamil)              :: self
        real(8), intent(in)       :: k(3)
        complex(8), intent(inout) :: H(:,:)
        logical                   :: has_E, has_hopp

        if(k(3) /= 0d0) then
            write (*,*) "K_z is non-zero. Abort.", k
            stop
        endif

        H =  0d0

        has_E = (self%E_s /= 0) .or. (self%E_A /= 0) .or. (self%E_B /= 0)
        if(has_E) call self%set_EigenE(H)

        has_hopp =   (self%Vss_sig /= 0d0) .or. (self%Vsp_sig /= 0d0) &
                .or. (self%Vpp_sig /= 0d0) .or. (self%Vpp_pi /= 0d0)
        if(has_hopp) call self%set_hopping(k,H)
        
        if(self%t_2       /= 0d0) call self%set_snd_hopping(k,H)
        if(self%t_so      /= 0d0) call self%set_rashba_SO(k,h)
        if(self%lambda    /= 0d0) call self%set_loc_exch(H)
        
    end subroutine setup_H

    subroutine test_herm(H, tag)
        implicit none
        complex(8), intent(in)      :: H(:,:)
        integer(4)                  :: n, info
        character(len=*), optional  :: tag

        n = size(H, dim=1)

        if((my_norm2(reshape( real(H - transpose(conjg(H))),[n*n]))/(n**2) > 1d-10) .or.&
            (my_norm2(reshape(aimag(H - transpose(conjg(H))),[n*n]))/(n**2) > 1d-10)) then
            if(present(tag)) then
                write (*,*) "nope - > ", tag
            else
                write (*,*) "nope"
            endif
            call MPI_Abort(MPI_COMM_WORLD, 0, info)
        endif
    end subroutine test_herm

    function init_hamil(cfg) result(self)
        implicit none
        type(CFG_t)    :: cfg
        type(hamil)    :: self
        real(8)        :: tmp
        integer(4)     :: ierr, n

        call MPI_Comm_size(MPI_COMM_WORLD, self%nProcs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)

        self%units = init_units(cfg, self%me)
        self%UC    = init_unit(cfg)

        if(self%me ==  0) then 
            call CFG_get(cfg, "hamil%Vss_sig", tmp)
            self%Vss_sig =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%Vsp_sig", tmp)
            self%Vsp_sig =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%Vpp_pi", tmp)
            self%Vpp_pi =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%Vpp_sig", tmp)
            self%Vpp_sig =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%t_2", tmp)
            self%t_2 =  tmp * self%units%energy

            call CFG_get(cfg, "hamil%phi_2", self%phi_2)

            call CFG_get(cfg, "hamil%t_so", tmp)
            self%t_so =  tmp * self%units%energy

            call CFG_get(cfg, "hamil%lambda", tmp)
            self%lambda =  tmp * self%units%energy

            call CFG_get(cfg, "hamil%lambda_nl", tmp)
            self%lambda_nl =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%E_s", tmp)
            self%E_s =  tmp * self%units%energy
            
            call CFG_get(cfg, "hamil%E_A", tmp)
            self%E_A =  tmp * self%units%energy
            call CFG_get(cfg, "hamil%E_B", tmp)
            self%E_B =  tmp * self%units%energy

            call CFG_get(cfg, "hamil%n", n)
            self%num_orb =  2 * n +  1
            self%num_up  =  self%num_orb *  self%UC%num_atoms
        endif
        call self%Bcast_hamil()
    end function init_hamil

    subroutine Bcast_hamil(self)
        implicit none
    class(hamil)          :: self
        integer(4), parameter :: num_cast =  11
        integer(4)            :: ierr(num_cast)

        call MPI_Bcast(self%E_s,       1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(1))
        call MPI_Bcast(self%E_A,       1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(2))
        call MPI_Bcast(self%E_B,       1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(3))
        call MPI_Bcast(self%Vss_sig,   1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(4))
        call MPI_Bcast(self%t_2,       1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(5))
        call MPI_Bcast(self%phi_2,     1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(6))
        call MPI_Bcast(self%t_so,      1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(7))
        call MPI_Bcast(self%lambda,    1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(8))
        call MPI_Bcast(self%lambda_nl, 1, MPI_REAL8,    root, MPI_COMM_WORLD, ierr(9))
        call MPI_Bcast(self%num_orb,   1, MPI_INTEGER4, root, MPI_COMM_WORLD, ierr(10))
        call MPI_Bcast(self%num_up,    1, MPI_INTEGER4, root, MPI_COMM_WORLD, ierr(11))

        call check_ierr(ierr, self%me, "Hamiltionian check err")
    end subroutine

    subroutine set_loc_exch(self,H)
        implicit none
    class(hamil), intent(in)   :: self
        complex(8), intent(inout)  :: H(:,:)
        complex(8)                 :: S(2,2) !> Stonermatrix
        integer(4)                 :: i, i_up, i_dw, atm, m, j, j_up, j_dw

        m = self%num_orb - 1
        atm = 1

        do i =  1,self%num_up,self%num_orb
            i_up =  i
            i_dw =  i +  self%num_up 

            call self%setup_Stoner_mtx(atm,S)

            do j = 0,self%num_orb-1
                j_up = j + i_up
                j_dw = j + i_dw
                H(j_up,j_up) = H(j_up, j_up) + S(1,1)
                H(j_up,j_dw) = H(j_up, j_dw) + S(1,2)
                H(j_dw,j_up) = H(j_dw, j_up) + S(2,1)
                H(j_dw,j_dw) = H(j_dw, j_dw) + S(2,2)
            enddo
            atm = atm + 1
        enddo
    end subroutine set_loc_exch

    subroutine setup_Stoner_mtx(self,i,S)
        implicit none
    class(hamil), intent(in) :: self
        integer(4), intent(in)   :: i
        complex(8), intent(inout):: S(2,2)
        real(8)                  :: m(3), fac

        m = self%UC%atoms(i)%get_m_cart()
        fac =  - 0.5d0 *  self%lambda

        S = fac * ( m(1) * sigma_x &
            + m(2) * sigma_y &
            + m(3) * sigma_z)
    end subroutine setup_Stoner_mtx

    subroutine set_EigenE(self,H)
        implicit none
    class(hamil), intent(in) :: self
        complex(8), intent(inout):: H(:,:)
        integer(4) :: i,id, ierr

        do i =  1,self%num_up
            id = i +  self%num_up
            select case(self%UC%atoms(i)%site_type)
            case(A_site)
                H(i,  i)  = H(i,   i) + self%E_A
                H(id, id) = H(id, id) + self%E_A
            case(B_site)
                H(i,  i)  = H(i,   i) + self%E_B
                H(id, id) = H(id, id) + self%E_B
            case(no_site)
                H(i,  i)  = H(i,   i) + self%E_s 
                H(id, id) = H(id, id) + self%E_s
            case default
                write (*,*) i, "unknown site type: ", self%UC%atoms(i)%site_type
                call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
            end select

        enddo
    end subroutine set_EigenE 

    subroutine set_hopping(self,k, H)
        implicit none
    class(hamil), intent(in)          :: self 
        real(8), intent(in)               :: k(3)
        complex(8), intent(inout)         :: H(:,:)
        integer(4)                        :: i, i_d, j,&
            j_d, conn, m, cnt
        real(8)                           :: k_dot_r, hopp_mtx(self%num_orb, self%num_orb)
        complex(8)                        :: new(self%num_orb, self%num_orb)

        call self%set_hopp_mtx(hopp_mtx)
        m = self%num_orb - 1

        ! Spin up
        cnt = 1
        do i = 1,self%num_up,self%num_orb
            do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
                if(self%UC%atoms(cnt)%conn_type(conn) == nn_conn) then
                    j =  self%UC%atoms(cnt)%neigh_idx(conn)
                    k_dot_r =  dot_product(k, self%UC%atoms(cnt)%neigh_conn(conn,:))

                    new = exp(i_unit * k_dot_r) * hopp_mtx
                    H(i:i+m, j:j+m) =  H(i:i+m,j:j+m) + new
                    H(j:j+m, i:i+m) =  H(j:j+m,i:i+m) + conjg(new)
                endif
            enddo
            cnt = cnt + 1
        enddo

        ! Spin down
        cnt = 1
        do i = 1,self%num_up,self%num_orb
            i_d =  i + self%num_up
            do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
                if(self%UC%atoms(cnt)%conn_type(conn) == nn_conn) then
                    j      = self%UC%atoms(cnt)%neigh_idx(conn)
                    j_d = j + self%num_up

                    k_dot_r =  dot_product(k, self%UC%atoms(cnt)%neigh_conn(conn,:))

                    new =  exp(i_unit *  k_dot_r) * hopp_mtx
                    H(i_d:i_d+m, j_d:j_d+m) = H(i_d:i_d+m, j_d:j_d+m) + new
                    H(j_d:j_d+m, i_d:i_d+m) = H(j_d:j_d+m, i_d:i_d+m) + conjg(new)
                endif
            enddo
            cnt = cnt + 1
        enddo
    end subroutine set_hopping

    subroutine set_SOC(self, H)
        implicit none
    class(hamil), intent(in)              :: self
        complex(8), intent(inout)         :: H(:,:)
        integer(4)                        :: ierr, i, i_d, m

        if(self%num_orb /= 4) then
            if(self%me == root) write (*,*) "SOC only implemented for p-orbitals"
            call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
        endif
        
        m = self%num_orb - 1

        if(self%num_orb == 4) then
            do i=1,self%num_up,self%num_orb
                i_d = i + self%num_up

                H(i:  i+m,   i:  i+m)   = H(i:  i+m,   i:  i+m)   + self%E_soc * Lz
                H(i_d:i_d+m, i_d:i_d+m) = H(i_d:i_d+m, i_d:i_d+m) - self%E_soc * Lz
                H(i:  i+m,   i_d:i_d+m) = H(i:  i+m,   i_d:i_d+m) + self%E_soc * LxpILy 
                H(i_d:i_d+m, i:  i+m)   = H(i_d:i_d+m, i:  i+m)   + self%E_soc * LxmILy
            enddo
        endif
    end subroutine set_SOC

    subroutine set_snd_hopping(self,k, H)
        implicit none
    class(hamil), intent(in)              :: self 
        real(8), intent(in)               :: k(3)
        complex(8), intent(inout)         :: H(:,:)
        integer(4)                        :: i, i_d, j, j_d, conn
        real(8)                           :: k_dot_r
        complex(8)                        :: forw, back, t_full

        t_full =  self%t_2 * exp(i_unit * self%phi_2)
        
        ! Spin up
        do i = 1,self%num_up
            i_d =  i + self%num_up
            do conn = 1,size(self%UC%atoms(i)%neigh_idx)
                if(self%UC%atoms(i)%conn_type(conn) == snd_nn_conn) then
                    j =  self%UC%atoms(i)%neigh_idx(conn)
                    j_d = j + self%num_up
                    
                    k_dot_r =  dot_product(k, self%UC%atoms(i)%neigh_conn(conn,:))
                    forw =  exp(i_unit * k_dot_r) * t_full
                    back =  conjg(forw)

                    H(i,j) =  H(i,j) + forw 
                    H(j,i) =  H(j,i) + back

                    H(i_d, j_d) = H(i_d, j_d) + forw 
                    H(j_d, i_d) = H(j_d, i_d) + back
                endif
            enddo
        enddo

    end subroutine set_snd_hopping

    subroutine set_rashba_SO(self, k, H)
        implicit none
    class(hamil), intent(in)    :: self
        real(8), intent(in)         :: k(3)
        complex(8), intent(inout)   :: H(:,:)
        integer(4)                  :: i, conn, j, i_d, j_d
        real(8)                     :: k_dot_r, d_ij(3)
        complex(8)                  :: new(2,2)

        do i =  1, self%num_up
            i_d =  i + self%num_up
            do conn =  1,size(self%UC%atoms(i)%neigh_idx)
                if(self%UC%atoms(i)%conn_type(conn) == nn_conn) then
                    j =  self%UC%atoms(i)%neigh_idx(conn)
                    j_d = j + self%num_up

                    k_dot_r =  dot_product(k, self%UC%atoms(i)%neigh_conn(conn,:))

                    ! set and normalize d_ij pointing from j to i
                    d_ij =  - self%UC%atoms(i)%neigh_conn(conn,:)
                    d_ij =  d_ij / my_norm2(d_ij)

                    ! hopping from i to j
                    new = exp(i_unit *  k_dot_r) * i_unit * self%t_so &
                        * (sigma_x * d_ij(2) -  sigma_y * d_ij(1))


                    H(i, j) = H(i, j) + new(1,       1)
                    H(j, i) = H(j, i) + conjg(new(1, 1))

                    H(i,   j_d) = H(i,   j_d) + new(1,       2)
                    H(j_d, i)   = H(j_d, i)   + conjg(new(1, 2))

                    H(i_d, j)   = H(i_d, j) + new(2,         1)
                    H(j,   i_d) = H(j, i_d) + conjg(new(2, 1))

                    H(i_d, j_d) = H(i_d, j_d) + new(2,       2)
                    H(j_d, i_d) = H(j_d, i_d) + conjg(new(2, 2))
                endif
            enddo
        enddo
    end subroutine set_rashba_SO

    subroutine set_hopp_mtx(self, hopp_mtx)
        implicit none
    class(hamil), intent(in)     :: self
        real(8)                  :: hopp_mtx(self%num_orb, self%num_orb)
        integer(4)               :: ierr
        
        if(self%num_orb ==  1) then
            hopp_mtx(1,1) =  self%Vss_sig
        elseif(self%num_orb == 4) then
            if(self%me == root) write (*,*) &
                "P orbitals not implemented yet. See comments for necessary adaptations"
            call MPI_Abort(MPI_COMM_WORLD, 0, ierr)

            ! 1. rotations away from R =  R * e_z (See
            ! 2. figure out why this matrix isn't hermitian (why is Vsp_sig complex)
            !hopp_mtx =  transpose(reshape(&
                ![self%Vss_sig,  0d0,         0d0,         self%Vsp_sig, &
                !0d0,            self%Vpp_pi, 0d0,         0d0,          &
                !0d0,            0d0,         self%Vpp_pi, 0d0,          &
                !- self%Vsp_sig, 0d0,         0d0,         self%Vpp_sig], [4,4]))
        endif
    end subroutine

    subroutine set_deriv_FD(self, k, k_idx, del_H)
        implicit none
    class(hamil), intent(in) :: self
        real(8), intent(in)      :: k(3)
        integer(4), intent(in)   :: k_idx
        complex(8), allocatable     :: H_forw(:,:), H_back(:,:), del_H(:,:)
        real(8) :: k_forw(3), k_back(3)
        real(8), parameter :: delta_k =  1d-6
        integer(4)         :: N

        N = 2 * self%num_up
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
    class(hamil)                  :: self
        integer(4), intent(in)    :: k_idx
        real(8), intent(in)       :: k(3)
        logical                   :: has_hopp

        if(k(3) /= 0) then
            write (*,*) "K_z not zero in set_derivative_k"
            stop
        endif

        self%del_H = 0d0
        has_hopp =   (self%Vss_sig /= 0d0) .or. (self%Vsp_sig /= 0d0) &
                .or. (self%Vpp_sig /= 0d0) .or. (self%Vpp_pi /= 0d0)
        if(has_hopp) call self%set_derivative_hopping(k, k_idx)

        if(self%t_so /= 0d0) call self%set_derivative_rashba_so(k, k_idx)
        if(self%t_2  /= 0d0) call self%set_derivative_snd_hopping(k, k_idx)

    end subroutine set_derivative_k

    subroutine set_derivative_hopping(self, k, k_idx)
        implicit none
    class(hamil)             :: self
        real(8), intent(in)      :: k(3)
        integer(4), intent(in)   :: k_idx
        real(8)                   :: r(3), k_dot_r, hopp_mtx(self%num_orb, self%num_orb)
        complex(8)                :: forw(self%num_orb, self%num_orb), back(self%num_orb, self%num_orb)
        integer(4)                :: i, j, conn, i_d, j_d, m, cnt

        call self%set_hopp_mtx(hopp_mtx)
        m =  self%num_orb - 1

        !$omp parallel do default(shared) &
        !$omp& private(i_d, conn, j, j_d, r, k_dot_r, forw, back, i)&
        !$omp& schedule(static)
        do cnt = 1,self%UC%num_atoms
            i   =  1 + (cnt-1)*self%num_orb
            i_d =  i + self%num_up
            do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
                if(self%UC%atoms(cnt)%conn_type(conn) == nn_conn) then
                    j   = self%UC%atoms(cnt)%neigh_idx(conn)
                    j_d = j + self%num_up
                    r   = self%UC%atoms(cnt)%neigh_conn(conn,:)

                    k_dot_r = dot_product(k, r)
                    forw    = i_unit * r(k_idx) * hopp_mtx * exp(i_unit * k_dot_r)
                    back    =  transpose(conjg(forw))             

                    !Spin up
                    self%del_H(i:i+m,j:j+m)     = self%del_H(i:i+m,j:j+m) + forw 
                    self%del_H(j:j+m,i:i+m)     = self%del_H(j:j+m,i:i+m) + back

                    !Spin down
                    self%del_H(i_d:i_d+m,j_d:j_d+m) = self%del_H(i_d:i_d+m,j_d:j_d+m) + forw
                    self%del_H(j_d:j_d+m,i_d:i_d+m) = self%del_H(j_d:j_d+m,i_d:i_d+m) + back
                endif
            enddo
        enddo
    end subroutine set_derivative_hopping
    
    subroutine set_derivative_snd_hopping(self, k, k_idx)
        implicit none
        class(hamil)              :: self
        real(8), intent(in)       :: k(3)
        integer(4), intent(in)    :: k_idx
        real(8)                   :: r(3), k_dot_r
        complex(8)                :: forw, back, t_full
        integer(4)                :: i, j, conn, i_d, j_d
        
        t_full =  self%t_2 * exp(i_unit * self%phi_2)

        !$omp parallel do default(shared) &
        !$omp& private(i_d, conn, j, j_d, r, k_dot_r, forw, back)&
        !$omp& schedule(static)
        do i = 1,self%num_up
            i_d =  i + self%num_up
            do conn = 1,size(self%UC%atoms(i)%neigh_idx)
                if(self%UC%atoms(i)%conn_type(conn) == snd_nn_conn) then
                    j   = self%UC%atoms(i)%neigh_idx(conn)
                    j_d = j + self%num_up
                    r   = self%UC%atoms(i)%neigh_conn(conn,:)

                    k_dot_r = dot_product(k, r)
                    forw    = i_unit * r(k_idx) * t_full * exp(i_unit * k_dot_r)
                    back =  conjg(forw)                

                    !Spin up
                    self%del_H(i,j)     = self%del_H(i,j) + forw 
                    self%del_H(j,i)     = self%del_H(j,i) + back

                    !Spin down
                    self%del_H(i_d,j_d) = self%del_H(i_d,j_d) + forw
                    self%del_H(j_d,i_d) = self%del_H(j_d,i_d) + back
                endif
            enddo
        enddo
    end subroutine set_derivative_snd_hopping

    subroutine set_derivative_rashba_so(self, k, k_idx)
        implicit none
    class(hamil)            :: self
        real(8), intent(in)     :: k(3)
        integer(4), intent(in)  :: k_idx
        real(8)                 :: r(3), d_ij(3), k_dot_r
        integer(4)              :: i, j, i_d, j_d, conn
        complex(8)              :: e_z_sigma(2,2), forw(2,2), back(2,2)

        !$omp  parallel do default(shared) & 
        !$omp& private(i_d, conn, j, j_d, r, d_ij, k_dot_r, e_z_sigma, forw, back)&
        !$omp& schedule(static)
        do i = 1,self%num_up
            i_d = i +  self%num_up
            do conn = 1,size(self%UC%atoms(i)%neigh_idx)
                if(self%UC%atoms(i)%conn_type(conn) == nn_conn) then
                    j   =  self%UC%atoms(i)%neigh_idx(conn)
                    j_d = j + self%num_up 

                    r    = self%UC%atoms(i)%neigh_conn(conn,:)
                    d_ij = - r / my_norm2(r)

                    k_dot_r =  dot_product(k, r)

                    e_z_sigma = sigma_x * d_ij(2) - sigma_y * d_ij(1)
                    forw =  - self%t_so * e_z_sigma * r(k_idx) * exp(i_unit * k_dot_r)
                    back =  conjg(forw)

                    self%del_H(i,j)      = self%del_H(i,j)     + forw(1,1)
                    self%del_H(j,i)      = self%del_H(j,i)     + back(1,1)

                    self%del_H(i,j_d)    = self%del_H(i,j_d)   + forw(1,2)
                    self%del_H(j_d,i)    = self%del_H(j_d,i)   + back(1,2)

                    self%del_H(i_d,j)    = self%del_H(i_d,j)   + forw(2,1)
                    self%del_H(j,i_d)    = self%del_H(j,i_d)   + back(2,1)

                    self%del_H(i_d, j_d) = self%del_H(i_d,j_d) + forw(2,2)
                    self%del_H(j_d, i_d) = self%del_H(j_d,i_d) + back(2,2)
                endif
            enddo
        enddo
    end subroutine set_derivative_rashba_so

    subroutine calc_velo_mtx(self, k, derive_idx, eig_vec_mtx, ret)
        implicit none
    class(hamil)                    :: self
        real(8), intent(in)             :: k(3)
        integer(4), intent(in)          :: derive_idx
        complex(8), intent(in)          :: eig_vec_mtx(:,:)
        complex(8), allocatable         :: ret(:,:), tmp(:,:)
        integer(4)   :: n_dim

        n_dim = 2 * self%num_up
        allocate(tmp(n_dim, n_dim))

        call self%set_derivative_k(k, derive_idx)
        if(allocated(ret))then
            if(size(ret,1) /= n_dim .or. size(ret,2) /= n_dim) then
                deallocate(ret)
            endif
        endif
        if(.not. allocated(ret)) allocate(ret(n_dim, n_dim))

        call zgemm('N', 'N', n_dim, n_dim, n_dim, &
            c_1, self%del_H, n_dim,&
            eig_vec_mtx, n_dim,&
            c_0, tmp, n_dim)
        call zgemm('C', 'N', n_dim, n_dim, n_dim, &
            c_1, eig_vec_mtx, n_dim, &
            tmp, n_dim, &
            c_0, ret, n_dim)
        deallocate(tmp)
    end subroutine calc_velo_mtx

    subroutine calc_eig_and_velo(self, k, eig_val, del_kx, del_ky)
        implicit none
    class(hamil)             :: self
        real(8), intent(in)      :: k(3)
        real(8)                  :: eig_val(:)
        complex(8), allocatable  :: eig_vec(:,:), del_kx(:,:), del_ky(:,:), work(:)
        real(8), allocatable     :: rwork(:)
        integer(4), allocatable  :: iwork(:)
        integer(4)   :: n_dim, lwork, lrwork, liwork, info

        n_dim = 2 * self%num_up
        if(.not. allocated(eig_vec)) allocate(eig_vec(n_dim,n_dim))

        eig_vec = 0d0
        call self%setup_H(k, eig_vec)

        call calc_zheevd_size('V', eig_vec, eig_val, lwork, lrwork, liwork)
        allocate(work(lwork), stat=info)
        allocate(rwork(lrwork))
        allocate(iwork(liwork))


        call zheevd('V', 'L', n_dim, eig_vec, n_dim, eig_val, &
            work, lwork, rwork, lrwork, iwork, liwork, info)
        if(info /= 0) then
            write (*,*) "ZHEEVD in berry calculation failed"
        endif
        
        call self%calc_velo_mtx(k, 1, eig_vec, del_kx)
        call self%calc_velo_mtx(k, 2, eig_vec, del_ky)
        deallocate(work)
        deallocate(rwork)
        deallocate(iwork)
    end subroutine calc_eig_and_velo

    subroutine calc_berry_z(self, z_comp, eig_val, x_mtx, y_mtx)
        implicit none
    class(hamil)             :: self
        real(8)                  :: eig_val(:), dE
        real(8)                  :: z_comp(:) !> \f$ \Omega^n_z \f$
        complex(8)               :: x_mtx(:,:), y_mtx(:,:)
        complex(8) :: fac
        integer(4)   :: n_dim, n, m

        n_dim = 2 * self%num_up
        z_comp =  0d0
        do n = 1,n_dim
            do m = 1,n_dim
                if(n /= m) then
                    dE =  eig_val(n) - eig_val(m)
                    ! fac =  Re[1/(dE + ieta)]^2
                    fac =  dE**2/(dE**2 + eta_sq)**2
                    z_comp(n) = z_comp(n) - 2d0 &
                        * aimag(fac * x_mtx(n,m) * y_mtx(m,n))
                endif
            enddo
        enddo

    end subroutine calc_berry_z

    Subroutine  calc_eigenvalues(self, k_list, eig_val)
        Implicit None
    class(hamil)                      :: self
        real(8), intent(in)               :: k_list(:,:)
        real(8), allocatable,intent(out)  :: eig_val(:,:)
        real(8)                           :: k(3)
        complex(8), allocatable           :: H(:,:)
        integer(4) :: i, N, lwork, lrwork, liwork, info, percentage
        real(8), allocatable              :: RWORK(:)
        complex(8), allocatable           :: WORK(:)
        integer(4), allocatable           :: IWORK(:)

        N =  2 * self%num_up
        allocate(eig_val(N, size(k_list, 2)))
        allocate(H(N,N))

        call calc_zheevd_size('N', H, eig_val(:,1), lwork, lrwork, liwork)
        allocate(RWORK(lrwork))
        allocate(IWORK(liwork))
        allocate(WORK(lwork))

        percentage = 0
        do i = 1,size(k_list,2)
            k =  k_list(:,i)
            call self%setup_H(k, H)

            call zheevd('N', 'U', N, H, N, eig_val(:,i), WORK, lwork, &
                RWORK, lrwork, IWORK, liwork, info)
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
        real(8)             , intent(out) :: eig_val(:)
        complex(8), allocatable           :: H(:,:), work(:)
        real(8), allocatable              :: rwork(:)
        integer(4), allocatable           :: iwork(:)
        integer(4)                        :: N, lwork, lrwork, liwork, info

        N = 2 * self%num_up

        !allocate(eig_val(N))
        allocate(H(N,N))
        call calc_zheevd_size('N', H, eig_val, lwork, lrwork, liwork)
        allocate(work(lwork))
        allocate(iwork(liwork))
        allocate(rwork(lrwork))

        call self%setup_H(k, H)
        call zheevd('N', 'U', N, H, N, eig_val, work, lwork, &
            RWORK, lrwork, IWORK, liwork, info)

        if(info /=  0) write (*,*) "ZHEEVD failed:", info

        deallocate(H)
        deallocate(work)
        deallocate(iwork)
        deallocate(rwork)
    end subroutine calc_single_eigenvalue

end module
