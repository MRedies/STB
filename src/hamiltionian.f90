module Class_hamiltionian
   use m_config
   use output
   use Class_unit_cell
   use m_npy
   use mpi
   use MYPI
   use Constants
   implicit none

   type hamil
      real(8), allocatable :: E_fermi(:) !> Fermi lvl
      real(8)              :: temp !> temperature used in fermi-dirac
      real(8)         :: E_s !> onsite eigenenergy
      real(8)         :: E_A, E_B !> onsite energies for A and B sites in honeycomb
      real(8)         :: E_p(3)
      real(8)         :: Vss_sig !> nearest neighbour hopping for s-orbital
      real(8)         :: Vpp_pi, Vpp_sig !> nearest neigh hopping for p-orbitals
      real(8)         :: V2pp_pi, V2pp_sig !> 2nd n-neigh hopping for p-orbitals
      real(8)         :: eta_soc
      real(8)         :: t_2 !> amplitude for 2nd nearest neighbour hopping
      real(8)         :: phi_2 !> polar angle of 2nd nearest neighbour hopping in rad
      real(8)         :: t_so !> Rashba spin orb
      real(8)         :: lambda !> local exchange
      real(8)         :: HB1, HB2, HB_eta !> parameters for hongbins model
      real(8)         :: lambda_KM !> parameter for Kane Mele term
      real(8)         :: gamma !> broadening, Greens function, sigma_xx
      real(8), allocatable       :: drop_Vx_layers(:), drop_Vy_layers(:)
      complex(8), allocatable    :: del_H(:,:)
      character(len=300)   :: prefix
      integer         :: nProcs
      integer         :: me
      integer         :: num_orb, num_up
      logical      :: test_run !> should unit tests be performed
      type(unit_cell) :: UC !> unit cell
      type(units)     :: units
   contains
      procedure :: fermi_distr                    => fermi_distr
      procedure :: set_fermi                      => set_fermi
      procedure :: write_fermi                    => write_fermi
      procedure :: Bcast_hamil                    => Bcast_hamil
      procedure :: setup_H                        => setup_H
      procedure :: calc_eigenvalues               => calc_eigenvalues
      procedure :: calc_single_eigenvalue         => calc_single_eigenvalue
      procedure :: set_EigenE                     => set_EigenE
      procedure :: set_p_energy                   => set_p_energy
      procedure :: set_hopping                    => set_hopping
      procedure :: set_haldane_hopping            => set_haldane_hopping
      procedure :: set_loc_exch                   => set_loc_exch
      procedure :: setup_Stoner_mtx               => setup_Stoner_mtx
      procedure :: set_derivative_k               => set_derivative_k
      procedure :: set_rashba_SO                  => set_rashba_SO
      procedure :: set_deriv_FD                   => set_deriv_FD
      procedure :: calc_berry_z                   => calc_berry_z
      procedure :: calc_berry_diag                   => calc_berry_diag
      procedure :: calc_berry_diag_sea            => calc_berry_diag_sea
      procedure :: calc_berry_diag_surf           => calc_berry_diag_surf
      procedure :: calc_fac_sea                   => calc_fac_sea
      procedure :: calc_fac_surf                  => calc_fac_surf
      procedure :: calc_fac_diag                  => calc_fac_diag
      procedure :: calc_velo_mtx                  => calc_velo_mtx
      procedure :: calc_right_pert_velo_mtx       => calc_right_pert_velo_mtx
      procedure :: calc_left_pert_velo_mtx        => calc_left_pert_velo_mtx
      procedure :: calc_eig_and_velo              => calc_eig_and_velo
      procedure :: calc_exch_firstord             => calc_exch_firstord
      procedure :: compare_derivative             => compare_derivative
      procedure :: set_derivative_hopping         => set_derivative_hopping
      procedure :: set_derivative_snd_hopping     => set_derivative_snd_hopping
      procedure :: set_derivative_haldane_hopping => set_derivative_haldane_hopping
      procedure :: set_derivative_rashba_so       => set_derivative_rashba_so
      procedure :: set_derivative_KM              => set_derivative_KM
      procedure :: set_hopp_mtx                   => set_hopp_mtx
      procedure :: set_small_SOC                  => set_small_SOC
      procedure :: set_SOC                        => set_SOC
      procedure :: free_ham                       => free_ham
      procedure :: set_snd_hopping                => set_snd_hopping
      procedure :: set_snd_hopp_mtx               => set_snd_hopp_mtx
      procedure :: set_hongbin_hopp               => set_hongbin_hopp
      procedure :: set_hongbin_SOC                => set_hongbin_SOC
      procedure :: set_KaneMele_exch              => set_KaneMele_exch
      procedure :: z_layer_states                 => z_layer_states
      procedure :: drop_layer_derivative          => drop_layer_derivative
   end type hamil

contains
   subroutine set_fermi(self, cfg)
      use mpi
      implicit none
      class(hamil)         :: self
      class(CFG_t)           :: cfg
      real(8)                :: tmp(3)
      integer                :: ierr
      integer                :: n_steps
   
      if(root == self%me) then
         call CFG_get(cfg, "berry%E_fermi", tmp)
      endif
      call MPI_Bcast(tmp, 3, MPI_REAL8, root, MPI_COMM_WORLD, ierr)
      n_steps = nint(tmp(3))
      tmp =  tmp *  self%units%energy
   
      call linspace(tmp(1), tmp(2), n_steps, self%E_fermi)
      call self%write_fermi()
   end subroutine set_fermi

   subroutine write_fermi(self)
      implicit none
      class(hamil)         :: self
      if(self%me ==  root) then
         call save_npy(trim(self%prefix) //  "fermi.npy", self%E_fermi)
      endif
   end subroutine write_fermi

   function fermi_distr(self, E, n_ferm) result(ferm)
      implicit none
      class(hamil), intent(in)    :: self
      real(8), intent(in)           :: E
      integer   , intent(in)        :: n_ferm
      real(8)                       :: ferm, exp_term
   
      exp_term =  (E - self%E_fermi(n_ferm)) /&
                 (boltzmann_const * self%temp)
      !cutting off at exp(x) =  10^16
      ! which is at the machine eps
      if(exp_term > 36d0) then
         ferm = 0d0
      elseif(exp_term < -36d0) then
         ferm = 1d0
      else
         ferm = 1d0 / (exp(exp_term) + 1d0)
      endif
   end function fermi_distr

   function z_layer_states(self) result(z)
      implicit none
      class(hamil), intent(in)      :: self
      real(8)                       :: z(2 * self%num_up)
      integer                       :: n_up, n_down, i_atm

      i_atm = 1
      do n_up = 1, self%num_up, self%num_orb
         n_down = n_up + self%num_up
         z(n_up:   n_up   + self%num_orb -1) = self%UC%atoms(i_atm)%pos(3)
         z(n_down: n_down + self%num_orb -1) = self%UC%atoms(i_atm)%pos(3)

         i_atm = i_atm + 1
      enddo
   end function z_layer_states

   subroutine drop_layer_derivative(self, k_idx)
      implicit None
      class(hamil)                  :: self
      integer, intent(in)           :: k_idx
      real(8), allocatable          :: layers(:)
      logical                       :: mask(2 * self%num_up)
      real(8)                       :: z(2 * self%num_up), t_start, t_stop
      real(8), parameter            :: eps = 1e-6
      integer                       :: i, m, ierr

      t_start = MPI_Wtime()
      if(k_idx == 1) then
         layers = self%drop_Vx_layers
      elseif(k_idx == 2) then
         layers = self%drop_Vy_layers
      else
         write (*,*) "Dropping kidx needs to be 1 or 2"
         call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
      endif

      if(size(layers) > 0) then
         write (*,*) "droping some layers"
         z = self%z_layer_states()
         do i = 1, size(layers)
            mask = abs(layers(i) - z) < eps

            do m = 1, 2*self%num_up
               where(mask) self%del_H(m,:) = 0d0
               where(mask) self%del_H(:,m) = 0d0
            enddo
         enddo
      endif

      t_stop = MPI_Wtime()
   end subroutine drop_layer_derivative

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
      integer                     :: N, k_idx

      N = 2 * self%num_up
      write (*,*) N
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
            call error_msg("FD comp failed")
            write (*,*) "FD"
            call save_npy("output/dbg/fd_H.npy",fd_H)
            write (*,*) "analytic"
            call save_npy("output/dbg/del_H.npy", self%del_H)
            !write (*,*) "diff"
            call error_msg("Not hermitian", abort=.True.)
         else
            call error_msg("All hermitian", c_green)
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
      ! on-site eigenenergy
      has_E = (self%E_s /= 0) .or. (self%E_A /= 0) .or. (self%E_B /= 0)
      if(has_E) call self%set_EigenE(H)

      if(any(self%E_p /= 0)) call self%set_p_energy(H)

      !hopping
      has_hopp =   (self%Vss_sig /= 0d0) &
                 .or. (self%Vpp_sig /= 0d0) &
                 .or. (self%Vpp_pi /= 0d0)
      if(has_hopp) call self%set_hopping(k,H)

      !2nd hopping
      has_hopp =   (self%V2pp_sig /= 0d0) &
                 .or. (self%V2pp_pi  /= 0d0)
      if(has_hopp) call self%set_snd_hopping(k,H)

      if(self%t_2       /= 0d0) call self%set_haldane_hopping(k,H)
      if(self%t_so      /= 0d0) call self%set_rashba_SO(k,H)
      if(self%eta_soc   /= 0d0) call self%set_SOC(H)
      if(self%lambda    /= 0d0) call self%set_loc_exch(H)
      if(self%lambda_KM /= 0d0) call self%set_KaneMele_exch(k,H)

      if(self%HB1 /= 0d0 .or. self%HB2 /= 0d0) then
         call self%set_hongbin_hopp(k, H)
      endif
      if(self%HB_eta /= 0d0) call self%set_hongbin_SOC(H)
   end subroutine setup_H

   subroutine test_herm(H, tag, verbose)
      implicit none
      complex(8), intent(in)      :: H(:,:)
      integer                     :: n
      character(len=*), optional  :: tag
      logical, optional           :: verbose
      logical                     :: p_verbose

      if(present(verbose)) then
         p_verbose = verbose
      else
         p_verbose = .False.
      endif

      n = size(H, dim=1)

      if((my_norm2(reshape( real(H - transpose(conjg(H))),[n*n]))/(n**2) > 1d-10) .or.&
         (my_norm2(reshape(aimag(H - transpose(conjg(H))),[n*n]))/(n**2) > 1d-10)) then
         if(present(tag)) then
            write (*,*) "Tag = ", tag
         endif

         call error_msg("not hermitian", abort=.True.)
      else
         if(p_verbose) write (*,*) tag, " is hermitian"
      endif
   end subroutine test_herm

   function init_hamil(cfg) result(self)
      implicit none
      type(CFG_t)    :: cfg
      type(hamil)    :: self
      real(8)        :: tmp
      integer        :: ierr
      integer        :: n, n_arr

      call MPI_Comm_size(MPI_COMM_WORLD, self%nProcs, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)

      self%units = init_units(cfg, self%me)
      self%UC    = init_unit(cfg)

      if(self%me ==  0) then
         call CFG_get(cfg, "berry%temperature", tmp)
         self%temp = tmp * self%units%temperature

         call CFG_get(cfg, "hamil%Vss_sig", tmp)
         self%Vss_sig =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%Vpp_pi", tmp)
         self%Vpp_pi =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%Vpp_sig", tmp)
         self%Vpp_sig =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%V2pp_pi", tmp)
         self%V2pp_pi =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%V2pp_sig", tmp)
         self%V2pp_sig =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%t_2", tmp)
         self%t_2 =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%phi_2", self%phi_2)

         call CFG_get(cfg, "hamil%t_so", tmp)
         self%t_so =  tmp * self%units%energy
         !write(*,*) "Rasbha spin orbit:", self%t_so,tmp

         call CFG_get(cfg, "hamil%eta_soc", tmp)
         self%eta_soc =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%lambda", tmp)
         self%lambda =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%E_s", tmp)
         self%E_s =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%E_A", tmp)
         self%E_A =  tmp * self%units%energy
         call CFG_get(cfg, "hamil%E_B", tmp)
         self%E_B =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%E_p", self%E_p)
         self%E_p =  self%E_p* self%units%energy

         call CFG_get(cfg, "hamil%n", n)
         self%num_orb =  2 * n +  1
         self%num_up  =  self%num_orb *  self%UC%num_atoms

         call CFG_get(cfg, "hamil%HB1", tmp)
         self%HB1 =  tmp * self%units%energy
         call CFG_get(cfg, "hamil%HB2", tmp)
         self%HB2 =  tmp * self%units%energy
         call CFG_get(cfg, "hamil%HB_eta", tmp)
         self%HB_eta =  tmp * self%units%energy

         call CFG_get(cfg, "hamil%lambda_KM", tmp)
         self%lambda_KM =  tmp * self%units%energy
         
         call CFG_get(cfg, "hamil%gamma", tmp)
         self%gamma =  tmp * self%units%energy

         call CFG_get(cfg, "general%test_run", self%test_run)

         call CFG_get_size(cfg, "layer_dropout%Vx", n_arr)
         allocate(self%drop_Vx_layers(n_arr))
         call CFG_get(cfg, "layer_dropout%Vx", self%drop_Vx_layers)

         call CFG_get_size(cfg, "layer_dropout%Vy", n_arr)
         allocate(self%drop_Vy_layers(n_arr))
         call CFG_get(cfg, "layer_dropout%Vy", self%drop_Vy_layers)
         call CFG_get(cfg, "output%band_prefix", self%prefix)
         if(self%prefix(len_trim(self%prefix):len_trim(self%prefix)) /=  "/") then
            self%prefix =  trim(self%prefix) // "/"
         endif
      endif
      call self%Bcast_hamil()
   end function init_hamil
   
   subroutine Bcast_hamil(self)
      implicit none
      class(hamil)          :: self
      integer   , parameter :: num_cast = 27
      integer               :: ierr(num_cast), Vx_len, Vy_len

      call MPI_Bcast(self%E_s,      1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(1))
      call MPI_Bcast(self%E_A,      1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(2))
      call MPI_Bcast(self%E_B,      1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(3))
      call MPI_Bcast(self%E_p,      3,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(4))
      call MPI_Bcast(self%Vss_sig,  1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(5))
      call MPI_Bcast(self%Vpp_pi,   1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(6))
      call MPI_Bcast(self%Vpp_sig,  1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(7))
      call MPI_Bcast(self%V2pp_pi,  1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(8))
      call MPI_Bcast(self%V2pp_sig, 1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(9))
      call MPI_Bcast(self%t_2,      1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(10))
      call MPI_Bcast(self%phi_2,    1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(11))
      call MPI_Bcast(self%t_so,     1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(12))
      call MPI_Bcast(self%lambda,   1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(13))
      call MPI_Bcast(self%eta_soc,  1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(14))
      call MPI_Bcast(self%num_orb,  1,              MYPI_INT,    &
                     root,          MPI_COMM_WORLD, ierr(15))
      call MPI_Bcast(self%num_up,   1,              MYPI_INT,    &
                     root,          MPI_COMM_WORLD, ierr(16))
      call MPI_Bcast(self%test_run, 1,              MPI_LOGICAL, &
                     root,         MPI_COMM_WORLD, ierr(17))

      call MPI_Bcast(self%HB1,    1,              MPI_REAL8, &
                     root,        MPI_COMM_WORLD, ierr(18))
      call MPI_Bcast(self%HB2,    1,              MPI_REAL8, &
                     root,        MPI_COMM_WORLD, ierr(19))
      call MPI_Bcast(self%HB_eta, 1,              MPI_REAL8, &
                     root,        MPI_COMM_WORLD, ierr(20))
      ! allocate and share Vx_dropout
      if(self%me == root) Vx_len = size(self%drop_Vx_layers)
      call MPI_Bcast(Vx_len, 1, MPI_INTEGER, &
                     root, MPI_COMM_WORLD, ierr(21))
      call MPI_Bcast(self%prefix,  300,          MPI_CHARACTER, &
                     root,        MPI_COMM_WORLD, ierr(22))
      ! allocate and share Vy_dropout
      if(self%me ==root) Vy_len = size(self%drop_Vy_layers)
      call MPI_Bcast(Vy_len, 1, MPI_INTEGER, &
                     root, MPI_COMM_WORLD, ierr(23))

      if(self%me /= root) allocate(self%drop_Vy_layers(Vy_len))
      call MPI_Bcast(self%drop_Vy_layers, Vy_len, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(24))
      call MPI_Bcast(self%lambda_KM, 1,              MPI_REAL8, &
                     root,        MPI_COMM_WORLD, ierr(25))
      if(self%me /= root) allocate(self%drop_Vx_layers(Vx_len))
      call MPI_Bcast(self%drop_Vx_layers, Vx_len, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(26))
      call MPI_Bcast(self%gamma,      1,              MPI_REAL8,   &
                     root,          MPI_COMM_WORLD, ierr(27))
      call check_ierr(ierr, self%me, "Hamiltionian check err")
   end subroutine

    subroutine set_loc_exch(self,H)
      implicit none
      class(hamil), intent(in)   :: self
      complex(8), intent(inout)  :: H(:,:)
      complex(8)                 :: S(2,2) !> Stonermatrix
      integer                    :: i, i_up, i_dw, atm, m, j, j_up, j_dw

      m = self%num_orb - 1
      atm = 1

      do i =  1,self%num_up,self%num_orb
         i_up =  i
         i_dw =  i +  self%num_up

         call self%setup_Stoner_mtx(atm,S)

         do j = 0,self%num_orb-1
            j_up = i_up + j
            j_dw = i_dw + j
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
      integer   , intent(in)   :: i
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
      integer    :: i, id, m

      m = self%num_orb - 1
      do i =  1,self%num_up, self%num_orb
         id = i +  self%num_up
         select case(self%UC%atoms(i)%site_type)
         case(A_site)
            H(i:i+m,   i:i+m)   = H(i:i+m,   i:i+m)   + self%E_A
            H(id:id+m, id:id+m) = H(id:id+m, id:id+m) + self%E_A
         case(B_site)
            H(i:i+m,   i:i+m)   = H(i:i+m,   i:i+m)   + self%E_B
            H(id:id+m, id:id+m) = H(id:id+m, id:id+m) + self%E_B
         case(no_site)
            H(i:i+m,   i:i+m)   = H(i:i+m,   i:i+m)   + self%E_s
            H(id:id+m, id:id+m) = H(id:id+m, id:id+m) + self%E_s
         case default
            write (*,*) i, "unknown site type: ", self%UC%atoms(i)%site_type
            call error_msg("Aborting due to unknown type site", abort=.True.)
         end select

      enddo
   end subroutine set_EigenE

   subroutine set_p_energy(self, H)
      implicit none
      class(hamil), intent(in)      :: self
      complex(8), intent(inout) :: H(:,:)
      integer                   :: i, N, p_idx

      if(self%num_orb == 3) then
         N = 2 * self%num_up

         do i=1,N
            p_idx = mod((i-1),3) + 1
            H(i,i) =  H(i,i) + self%E_p(p_idx)
         enddo
      else
         call error_msg("p_energy only for p-orbitals (duh)!", abort=.True.)
      endif
   end subroutine set_p_energy

   subroutine set_hopping(self,k, H)
      implicit none
      class(hamil), intent(in)          :: self
      real(8), intent(in)               :: k(3)
      complex(8), intent(inout)         :: H(:,:)
      integer                           :: i, i_d, j,&
                                           j_d, conn, m, cnt, n_idx
      real(8)                           :: k_dot_r, hopp_mtx(self%num_orb, self%num_orb), R(3)
      complex(8)                        :: new(self%num_orb, self%num_orb)

      m = self%num_orb - 1

      ! Spin up
      cnt = 1
      do i = 1,self%num_up,self%num_orb
         do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
            if(self%UC%atoms(cnt)%conn_type(conn) == nn_conn) then
               n_idx =  self%UC%atoms(cnt)%neigh_idx(conn)
               j     =  self%num_orb * (n_idx - 1) + 1

               R =  self%UC%atoms(cnt)%neigh_conn(conn,:)
               call self%set_hopp_mtx(R, hopp_mtx)
               k_dot_r =  dot_product(k, R)

               new = exp(i_unit * k_dot_r) * hopp_mtx
               H(i:i+m, j:j+m) =  H(i:i+m,j:j+m) + new
               H(j:j+m, i:i+m) =  H(j:j+m,i:i+m) + conjg(new)
            endif
         enddo
         cnt = cnt + 1
      enddo

      ! Spin down
      cnt = 1
      do i_d = self%num_up+1,2*self%num_up,self%num_orb
         do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
            if(self%UC%atoms(cnt)%conn_type(conn) == nn_conn) then
               n_idx =  self%UC%atoms(cnt)%neigh_idx(conn)
               j_d   =  self%num_orb * (n_idx - 1) + 1 + self%num_up

               R =  self%UC%atoms(cnt)%neigh_conn(conn,:)
               call self%set_hopp_mtx(R, hopp_mtx)
               k_dot_r =  dot_product(k, R)

               new =  exp(i_unit *  k_dot_r) * hopp_mtx
               H(i_d:i_d+m, j_d:j_d+m) = H(i_d:i_d+m, j_d:j_d+m) + new
               H(j_d:j_d+m, i_d:i_d+m) = H(j_d:j_d+m, i_d:i_d+m) + conjg(new)
            endif
         enddo
         cnt = cnt + 1
      enddo
   end subroutine set_hopping

   subroutine set_hongbin_hopp(self, k, H)
      implicit none
      class(hamil), intent(in)             :: self
      real(8), intent(in)              :: k(3)
      complex(8), intent(inout)        :: H(6,6)
      real(8), parameter               :: A = 1d0
      real(8)                          :: kx, ky, kz

      kx = k(1)
      ky = k(2)
      kz = k(3)
      if(kz /= 0d0) call error_msg("We only do 2d here", abort=.True.)

      H(1,1) = H(1,1) - 2d0 * self%HB1 * (cos(ky) + A)
      H(2,2) = H(2,2) - 2d0 * self%HB1 * (cos(kx) + A)
      H(3,3) = H(3,3) - 2d0 * self%HB1 * (cos(kx) + cos(ky))

      H(2,1) = H(2,1) + 4d0 * self%HB2 * sin(kx) * sin(ky)
      H(1,2) = conjg(H(2,1))

      H(4:6,4:6) = H(1:3,1:3)
   end subroutine set_hongbin_hopp

   subroutine set_hongbin_SOC(self, H)
      implicit none
      class(hamil), intent(in)          :: self
      complex(8), intent(inout)     :: H(6,6)
      integer, parameter            :: x_or_z = 1
      complex(8)                    :: soc_mtx(3,3)

      soc_mtx      =   c_0
      if(x_or_z ==  1) then
         soc_mtx(2,1) = - c_i
         soc_mtx(1,2) =   c_i
      elseif(x_or_z == 3) then
         soc_mtx(3,2) = - c_i
         soc_mtx(2,3) =   c_i
      else
         call error_msg("Hongbin hamiltonian only in x or z direction", &
                        abort=.True.)
      endif

      H(1:3,1:3) = H(1:3,1:3) + self%HB_eta * soc_mtx
      H(4:6,4:6) = H(4:6,4:6) + self%HB_eta * soc_mtx
   end subroutine set_hongbin_SOC

   subroutine set_snd_hopping(self,k, H)
      implicit none
      class(hamil), intent(in)              :: self
      real(8), intent(in)               :: k(3)
      complex(8), intent(inout)         :: H(:,:)
      integer                           :: i, i_d, j,&
                                           j_d, conn, m, cnt, n_idx
      real(8)                           :: k_dot_r, hopp_mtx(self%num_orb, self%num_orb), R(3)
      complex(8)                        :: new(self%num_orb, self%num_orb)

      m = self%num_orb - 1

      ! Spin up
      cnt = 1
      do i = 1,self%num_up,self%num_orb
         do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
            if(self%UC%atoms(cnt)%conn_type(conn) == snd_nn_conn) then
               n_idx =  self%UC%atoms(cnt)%neigh_idx(conn)
               j     =  self%num_orb * (n_idx - 1) + 1

               R =  self%UC%atoms(cnt)%neigh_conn(conn,:)
               call self%set_snd_hopp_mtx(R, hopp_mtx)
               k_dot_r =  dot_product(k, R)

               new = exp(i_unit * k_dot_r) * hopp_mtx
               H(i:i+m, j:j+m) =  H(i:i+m,j:j+m) + new
               H(j:j+m, i:i+m) =  H(j:j+m,i:i+m) + conjg(new)
            endif
         enddo
         cnt = cnt + 1
      enddo

      ! Spin down
      cnt = 1
      do i_d = self%num_up+1,2*self%num_up,self%num_orb
         do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
            if(self%UC%atoms(cnt)%conn_type(conn) == snd_nn_conn) then
               n_idx =  self%UC%atoms(cnt)%neigh_idx(conn)
               j_d   =  self%num_orb * (n_idx - 1) + 1 + self%num_up

               R =  self%UC%atoms(cnt)%neigh_conn(conn,:)
               call self%set_snd_hopp_mtx(R, hopp_mtx)
               k_dot_r =  dot_product(k, R)

               new =  exp(i_unit *  k_dot_r) * hopp_mtx
               H(i_d:i_d+m, j_d:j_d+m) = H(i_d:i_d+m, j_d:j_d+m) + new
               H(j_d:j_d+m, i_d:i_d+m) = H(j_d:j_d+m, i_d:i_d+m) + conjg(new)
            endif
         enddo
         cnt = cnt + 1
      enddo
   end subroutine set_snd_hopping

   subroutine set_p_hopp_mtx(R, Vpp_sig, Vpp_pi, hopp_mtx)
      implicit none
      real(8), intent(in)      :: R(3), Vpp_sig, Vpp_pi !> real-space connection between atoms
      real(8), intent(out)     :: hopp_mtx(3,3) !> hopping matrix
      real(8)                  :: l, m, n !> directional cosines

      l = R(1)/my_norm2(R)
      m = R(2)/my_norm2(R)
      n = R(3)/my_norm2(R)

      !cyclic permutations of H_x^x = l^2 Vpp_sig + (1-l^2) Vpp_pi
      hopp_mtx(1,1) = l * l * Vpp_sig + (1d0 - l * l) * Vpp_pi
      hopp_mtx(2,2) = m * m * Vpp_sig + (1d0 - m * m) * Vpp_pi
      hopp_mtx(3,3) = n * n * Vpp_sig + (1d0 - n * n) * Vpp_pi

      !cyclic permutations of H_x^y = l*m Vpp_sig - l*m Vpp_pi
      hopp_mtx(1,2) = l * m * (Vpp_sig - Vpp_pi)
      hopp_mtx(2,3) = m * n * (Vpp_sig - Vpp_pi)
      hopp_mtx(3,1) = n * l * (Vpp_sig - Vpp_pi)

      !cyclic permutations of H_x^z =  l*n (Vpp_sig -  Vpp_pi)
      hopp_mtx(1,3) = l * n * (Vpp_sig - Vpp_pi)
      hopp_mtx(2,1) = m * l * (Vpp_sig - Vpp_pi)
      hopp_mtx(3,2) = n * m * (Vpp_sig - Vpp_pi)
   end subroutine set_p_hopp_mtx

   subroutine set_SOC(self, H)
      implicit none
      class(hamil), intent(in)              :: self
      complex(8), intent(inout)         :: H(:,:)
      complex(8)                        :: loc_H(2,2)
      integer                           :: i_u, i_d, mu, nu, ms, ns, i_atm

      if(self%num_orb /= 3) then
         call error_msg("SOC only implemented for p-orbitals", abort=.True.)
      else
         i_atm = 1
         do i_u = 1,self%num_up,self%num_orb
            i_d =  i_u + self%num_up
            do mu = 1,self%num_orb
               do nu = 1,self%num_orb
                  call self%set_small_SOC(mu,nu, loc_H)

                  ms = mu - 1 ! mu-shift
                  ns = nu - 1 ! nu-shift

                  H(i_u+ms, i_u+ns) = H(i_u+ms, i_u+ns) + loc_H(1, 1)
                  H(i_d+ms, i_u+ns) = H(i_d+ms, i_u+ns) + loc_H(2, 1)
                  H(i_u+ms, i_d+ns) = H(i_u+ms, i_d+ns) + loc_H(1, 2)
                  H(i_d+ms, i_d+ns) = H(i_d+ms, i_d+ns) + loc_H(2, 2)
               enddo
            enddo
            i_atm = i_atm + 1
         enddo
      endif
   end subroutine set_SOC

   subroutine set_rot_SO(atm, U)
      implicit none
      type(atom), intent(in)   :: atm
      complex(8), intent(out)  :: U(2,2)
      complex(8)                  :: t_half, p_half

      t_half =  0.5d0 * atm%m_theta
      p_half =  0.5d0 * atm%m_phi

      U(1,1) =  exp(-i_unit * p_half) * cos(t_half)
      U(2,1) =  exp( i_unit * p_half) * sin(t_half)
      U(1,2) = -exp(-i_unit * p_half) * sin(t_half)
      U(2,2) =  exp( i_unit * p_half) * cos(t_half)
   end subroutine set_rot_SO

   subroutine set_small_SOC(self, mu, nu, H)
      implicit none
      class(hamil), intent(in)       :: self
      integer   , intent(in)     :: mu, nu
      complex(8), intent(out)    :: H(2,2)

      H(1,1) =   self%eta_soc *  Lz(mu,nu)
      H(2,1) =   self%eta_soc * (Lx(mu,nu) + i_unit * Ly(mu,nu))
      H(1,2) =   self%eta_soc * (Lx(mu,nu) - i_unit * Ly(mu,nu))
      H(2,2) = - self%eta_soc *  Lz(mu,nu)
   end subroutine set_small_SOC

   subroutine set_haldane_hopping(self,k, H)
      implicit none
      class(hamil), intent(in)              :: self
      real(8), intent(in)               :: k(3)
      complex(8), intent(inout)         :: H(:,:)
      integer                           :: i, i_d, j, j_d, conn
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

   end subroutine set_haldane_hopping

   subroutine set_KaneMele_exch(self, k, H)
    implicit none
    class(hamil), intent(in)    :: self
    real(8),intent(in)          ::k(3)
    complex(8), intent(inout)   :: H(:,:)
    real(8)                     ::a(3,2)
    real(8)                     ::k_n(3)
    integer                     ::conn,i,i_d,j,j_d
    real(8)                     ::KM,lambda_KM,x,y,f

    if(self%num_orb /= 1) call error_msg("Kane-Mele only for s-oritals", abort=.True.)

    do i = 1,self%num_up
         i_d =  i + self%num_up
         do conn = 1,size(self%UC%atoms(i)%neigh_idx)
            if(self%UC%atoms(i)%conn_type(conn) == snd_nn_conn) then
               a(:,size(self%UC%atoms(i)%neigh_idx) + 1 - conn) = self%UC%atoms(i)%neigh_conn(conn,:) ! connection vectors 
            endif
         enddo
      k_n = k/my_norm2(k)
      if (my_norm2(a(:,2))>pos_eps) then
         a(:,2) = a(:,2)/my_norm2(a(:,2))
      endif
      if (my_norm2(a(:,1))>pos_eps) then
         a(:,1) = a(:,1)/my_norm2(a(:,1))
      endif
      y = 0.5*dot_product(k_n,a(:,2) + a(:,1))
      x = 0.5*dot_product(k_n,a(:,1) - a(:,2))
      f = 4*(2*sin(y)*cos(x) - sin(2*x))
      KM = lambda_KM*f
      H(i,i) = H(i,i) + KM
      H(i_d,i_d) = H(i_d,i_d) - KM
      do conn = 1,size(self%UC%atoms(i)%neigh_idx)
         if(self%UC%atoms(i)%conn_type(conn) == nn_conn) then
            j =  self%UC%atoms(i)%neigh_idx(conn)
            j_d = j + self%num_up
            H(j,j) = H(j,j) - KM
            H(j_d,j_d) = H(j_d,j_d) + KM
         endif
      enddo
    enddo
    end subroutine set_KaneMele_exch

   subroutine set_rashba_SO(self, k, H)
      implicit none
      class(hamil), intent(in)    :: self
      real(8), intent(in)         :: k(3)
      complex(8), intent(inout)   :: H(:,:)
      integer                     :: i, conn, j, i_d, j_d
      real(8)                     :: k_dot_r, d_ij(3)
      complex(8)                  :: new(2,2)

      if(self%num_orb /= 1) call error_msg("Rashba only for s-oritals", abort=.True.)

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

   subroutine set_hopp_mtx(self, R, hopp_mtx)
      implicit none
      class(hamil), intent(in)     :: self
      real(8), intent(in)      :: R(3)
      real(8)                  :: hopp_mtx(self%num_orb, self%num_orb)

      if(self%num_orb ==  1) then
         hopp_mtx(1,1) =  self%Vss_sig
      elseif(self%num_orb == 3) then
         call set_p_hopp_mtx(R, self%Vpp_sig, self%Vpp_pi, hopp_mtx)
      endif
   end subroutine set_hopp_mtx

   subroutine set_snd_hopp_mtx(self, R, hopp_mtx)
      implicit none
      class(hamil), intent(in)     :: self
      real(8), intent(in)      :: R(3)
      real(8)                  :: hopp_mtx(self%num_orb, self%num_orb)

      if(self%num_orb ==  1) then
         call error_msg("snd hopping not implemented for s", abort=.True.)
      elseif(self%num_orb == 3) then
         call set_p_hopp_mtx(R, self%V2pp_sig, self%V2pp_pi, hopp_mtx)
      endif
   end subroutine set_snd_hopp_mtx

   subroutine set_deriv_FD(self, k, k_idx, del_H)
      implicit none
      class(hamil), intent(in) :: self
      real(8), intent(in)      :: k(3)
      integer   , intent(in)   :: k_idx
      complex(8), allocatable     :: H_forw(:,:), H_back(:,:), del_H(:,:)
      real(8) :: k_forw(3), k_back(3)
      real(8), parameter :: delta_k =  1d-6
      integer            :: N

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
      integer   , intent(in)    :: k_idx
      real(8), intent(in)       :: k(3)
      logical                   :: has_hopp, has_hong

      if(k(3) /= 0) then
         write (*,*) "K_z not zero in set_derivative_k"
         stop
      endif

      self%del_H = 0d0
      has_hopp =   (self%Vss_sig /= 0d0)  &
                 .or. (self%Vpp_sig /= 0d0)  &
                 .or. (self%Vpp_pi /= 0d0)
      if(has_hopp) call self%set_derivative_hopping(k, k_idx)

      if(self%t_so /= 0d0) call self%set_derivative_rashba_so(k, k_idx)
      if(self%t_2  /= 0d0) call self%set_derivative_haldane_hopping(k, k_idx)
      if(self%lambda_KM /= 0d0) call self%set_derivative_KM(k)
      has_hopp =  (self%V2pp_sig /= 0d0)  &
                 .or. (self%V2pp_pi /= 0d0)
      if(has_hopp) call self%set_derivative_snd_hopping(k, k_idx)

      has_hong =   (self%HB1 /= 0d0) &
                 .or. (self%HB2 /= 0d0) &
                 .or. (self%HB_eta /= 0d0)
      if(has_hong) call self%set_deriv_FD(k, k_idx, self%del_H)

      ! drop layers if necessary
      call self%drop_layer_derivative(k_idx)
   end subroutine set_derivative_k


   subroutine set_derivative_KM(self, k )
    implicit none
    class(hamil)                :: self
    real(8),intent(in)          ::k(3)
    real(8)                     ::a(3,2)
    real(8)                     ::k_n(3)
    integer                     ::conn,i,i_d,j,j_d
    real(8)                     ::abs_k,KM_deriv,lambda_KM,x,y,f_deriv

    lambda_KM = self%lambda_KM

    if(self%num_orb /= 1) call error_msg("Kane-Mele only for s-oritals", abort=.True.)
    lambda_KM = self%lambda_KM
    do i = 1,self%num_up
        i_d =  i + self%num_up
        self%del_H(i,i) = self%del_H(i,i) + KM_deriv
        self%del_H(i_d,i_d) = self%del_H(i_d,i_d) - KM_deriv
        do conn = 1,size(self%UC%atoms(i)%neigh_idx)
            if(self%UC%atoms(i)%conn_type(conn) == snd_nn_conn) then
               a(:,size(self%UC%atoms(i)%neigh_idx) + 1 - conn) = self%UC%atoms(i)%neigh_conn(conn,:) ! connection vectors
               k_n = k/my_norm2(k)
               if (my_norm2(a(:,2))>pos_eps) then
                  a(:,2) = a(:,2)/my_norm2(a(:,2))
               endif
               if (my_norm2(a(:,1))>pos_eps) then
                  a(:,1) = a(:,1)/my_norm2(a(:,1))
               endif
               !write (*,*) "y: " ,y
               write (*,*) "a2" ,a(:,2)
               !write (*,*) "x: " ,x
               write (*,*) "a1: " ,a(:,1)
               y = 0.5*dot_product(k_n,a(:,2) + a(:,1))
               x = 0.5*dot_product(k_n,a(:,1) - a(:,2))
               abs_k = my_norm2(k)
               f_deriv = 8*(2*cos(y)*cos(x)*y/abs_k - 2*sin(y)*sin(x)*x/abs_k - 2*x/abs_k*cos(2*x))
               KM_deriv = lambda_KM*f_deriv
               j =  self%UC%atoms(i)%neigh_idx(conn)
               j_d = j + self%num_up
               self%del_H(j,j) = self%del_H(j,j) - KM_deriv
               self%del_H(j_d,j_d) = self%del_H(j_d,j_d) + KM_deriv
            endif
        enddo
    enddo

    end subroutine set_derivative_KM
   subroutine set_derivative_hopping(self, k, k_idx)
      implicit none
      class(hamil)             :: self
      real(8), intent(in)      :: k(3)
      integer   , intent(in)   :: k_idx
      real(8)                   :: r(3), k_dot_r, hopp_mtx(self%num_orb, self%num_orb)
      complex(8)                :: forw(self%num_orb, self%num_orb), back(self%num_orb, self%num_orb)
      integer                   :: i, j, conn, i_d, j_d, m, cnt, n_idx

      m =  self%num_orb - 1

!$OMP         parallel do default(shared) &
!$OMP        & private(i_d, conn, j, j_d, r, k_dot_r, forw, back, i)&
!$OMP        & schedule(static)
      do cnt = 1,self%UC%num_atoms
         i   =  1 + (cnt-1)*self%num_orb
         i_d =  i + self%num_up
         do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
            if(self%UC%atoms(cnt)%conn_type(conn) == nn_conn) then
               n_idx = self%UC%atoms(cnt)%neigh_idx(conn)
               j     = self%num_orb * (n_idx - 1) + 1
               j_d = j + self%num_up

               r   = self%UC%atoms(cnt)%neigh_conn(conn,:)
               call self%set_hopp_mtx(r, hopp_mtx)
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
      class(hamil)             :: self
      real(8), intent(in)      :: k(3)
      integer   , intent(in)   :: k_idx
      real(8)                   :: r(3), k_dot_r, hopp_mtx(self%num_orb, self%num_orb)
      complex(8)                :: forw(self%num_orb, self%num_orb), back(self%num_orb, self%num_orb)
      integer                   :: i, j, conn, i_d, j_d, m, cnt, n_idx

      m =  self%num_orb - 1

!$OMP         parallel do default(shared) &
!$OMP        & private(i_d, conn, j, j_d, r, k_dot_r, forw, back, i)&
!$OMP        & schedule(static)
      do cnt = 1,self%UC%num_atoms
         i   =  1 + (cnt-1)*self%num_orb
         i_d =  i + self%num_up
         do conn = 1,size(self%UC%atoms(cnt)%neigh_idx)
            if(self%UC%atoms(cnt)%conn_type(conn) == snd_nn_conn) then
               n_idx = self%UC%atoms(cnt)%neigh_idx(conn)
               j     = self%num_orb * (n_idx - 1) + 1
               j_d = j + self%num_up

               r   = self%UC%atoms(cnt)%neigh_conn(conn,:)
               call self%set_snd_hopp_mtx(r, hopp_mtx)
               k_dot_r = dot_product(k, r)

               forw    = i_unit * r(k_idx) * hopp_mtx * exp(i_unit * k_dot_r)
               back    = transpose(conjg(forw))

               !Spin up
               self%del_H(i:i+m,j:j+m)     = self%del_H(i:i+m,j:j+m) + forw
               self%del_H(j:j+m,i:i+m)     = self%del_H(j:j+m,i:i+m) + back

               !Spin down
               self%del_H(i_d:i_d+m,j_d:j_d+m) = self%del_H(i_d:i_d+m,j_d:j_d+m) + forw
               self%del_H(j_d:j_d+m,i_d:i_d+m) = self%del_H(j_d:j_d+m,i_d:i_d+m) + back
            endif
         enddo
      enddo
   end subroutine set_derivative_snd_hopping

   subroutine set_derivative_haldane_hopping(self, k, k_idx)
      implicit none
      class(hamil)              :: self
      real(8), intent(in)       :: k(3)
      integer   , intent(in)    :: k_idx
      real(8)                   :: r(3), k_dot_r
      complex(8)                :: forw, back, t_full
      integer                   :: i, j, conn, i_d, j_d

      t_full =  self%t_2 * exp(i_unit * self%phi_2)

!$OMP         parallel do default(shared) &
!$OMP        & private(i_d, conn, j, j_d, r, k_dot_r, forw, back)&
!$OMP        & schedule(static)
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
   end subroutine set_derivative_haldane_hopping

   subroutine set_derivative_rashba_so(self, k, k_idx)
      implicit none
      class(hamil)            :: self
      real(8), intent(in)     :: k(3)
      integer   , intent(in)  :: k_idx
      real(8)                 :: r(3), d_ij(3), k_dot_r
      integer                 :: i, j, i_d, j_d, conn
      complex(8)              :: e_z_sigma(2,2), forw(2,2), back(2,2)

!$OMP          parallel do default(shared) &
!$OMP        & private(i_d, conn, j, j_d, r, d_ij, k_dot_r, e_z_sigma, forw, back)&
!$OMP        & schedule(static)
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

   subroutine calc_exch_firstord(self,eig_vec_mtx,eig_val,H_xc_1)
      implicit none
      class(hamil)                    :: self
      complex(8), intent(in)          :: eig_vec_mtx(:,:)
      real(8),intent(in)              :: eig_val(:)
      complex(8), allocatable         :: H_xc_1(:,:)
      complex(8), allocatable         :: temp(:,:),ret(:,:),H_temp(:,:)
      complex(8)                      :: theta(2),phi(2),theta_nc,theta_col,phi_nc,phi_col,dE,fac,Efac
      integer                         :: i,i_d,j,j_u,j_d,n_dim
      logical                         :: full
      n_dim = 2 * self%num_up
      if(allocated(H_temp)) deallocate(H_temp)
      if(allocated(temp)) deallocate(temp)
      if(allocated(ret)) deallocate(ret)
      allocate(H_temp(n_dim, n_dim))
      allocate(temp(n_dim, n_dim))
      allocate(ret(n_dim, n_dim))
      H_temp=(0d0,0d0)
      temp=(0d0,0d0)
      ret = (0d0,0d0)
      theta = self%UC%anticol_theta
      phi = self%UC%anticol_phi
      theta_nc = theta(2)
      theta_col = theta(1)
      phi_nc = phi(2)
      phi_col = phi(1)
      fac = 0.5d0 * self%lambda * PI/180d0
      i = 1
      i_d = i + self%num_up
      j = i + self%num_orb
      j_d = i_d + self%num_orb

      full = .False.
      do i =  1, self%num_up,self%num_orb
        i_d =  i + self%num_up
        do j = 0,self%num_orb-1
            j_u =  i + j
            j_d = i_d + j
            if(full) then
                temp(i,i)     =  fac*(cos(theta_col + theta_nc/2d0) - cos(theta_col))
                temp(i_d,i_d) = -fac*(cos(theta_col + theta_nc/2d0) - cos(theta_col))
                temp(j_u,j_u)     =  fac*(cos(theta_col - theta_nc/2d0) - cos(theta_col))
                temp(j_d,j_d) = -fac*(cos(theta_col - theta_nc/2d0) - cos(theta_col))
                temp(i,i_d)   =  fac*(sin(theta_col + theta_nc/2d0) - sin(theta_col))*exp(-i_unit*(phi_col+phi_nc/2d0))
                temp(i_d,i)   =  fac*(sin(theta_col + theta_nc/2d0) - sin(theta_col))*exp( i_unit*(phi_col+phi_nc/2d0))
                temp(j_u,j_d)   =  fac*(sin(theta_col - theta_nc/2d0) - sin(theta_col))*exp(-i_unit*(phi_col-phi_nc/2d0))
                temp(j_d,j_u)   =  fac*(sin(theta_col - theta_nc/2d0) - sin(theta_col))*exp( i_unit*(phi_col-phi_nc/2d0))
            else
                temp(i,i)     = -fac*sin(theta_col)*theta_nc/2d0
                temp(i_d,i_d) =  fac*sin(theta_col)*theta_nc/2d0
                temp(j_u,j_u)     =  fac*sin(theta_col)*theta_nc/2d0
                temp(j_d,j_d) = -fac*sin(theta_col)*theta_nc/2d0
                temp(i,i_d)   =  fac*cos(theta_col)*theta_nc/2d0*exp(-i_unit*(phi_col+phi_nc/2d0))
                temp(i_d,i)   =  fac*cos(theta_col)*theta_nc/2d0*exp( i_unit*(phi_col+phi_nc/2d0))
                temp(j_u,j_d)   = -fac*cos(theta_col)*theta_nc/2d0*exp(-i_unit*(phi_col-phi_nc/2d0))
                temp(j_d,j_u)   = -fac*cos(theta_col)*theta_nc/2d0*exp( i_unit*(phi_col-phi_nc/2d0))
            endif
        enddo
      enddo
      !rotate hxc into the eigenbasis of the hamiltonian with E_dagger H E
      call zgemm('N', 'N', n_dim, n_dim, n_dim, &
                  c_1, temp, n_dim,&
                  eig_vec_mtx, n_dim,&
                  c_0, ret, n_dim)
      deallocate(temp)
      call zgemm('C', 'N', n_dim, n_dim, n_dim, &
                  c_1, eig_vec_mtx, n_dim,&
                  ret, n_dim,&
                  c_0, H_temp, n_dim)
      deallocate(ret)
      do i=1,n_dim
         do j=1,n_dim
            if(i /= j) then
               !the sign here is tested, it is correct this way, 
               !also the energy factor would scale the outcome by one order of magn.
               dE = (eig_val(i)-eig_val(j))! + (H_temp(i,i)-H_temp(j,j))
               Efac =  dE/(dE + eta_sq)**2
               H_temp(i,j) = H_temp(i,j)*Efac
            else if(i==j) then
               H_temp(i,j) = (0d0,0d0)
            endif
         enddo
      enddo
      if(allocated(H_xc_1))then
         if(size(H_xc_1,1) /= n_dim .or. size(H_xc_1,2) /= n_dim) then
            deallocate(H_xc_1)
         endif
      endif
      if(.not. allocated(H_xc_1)) allocate(H_xc_1(n_dim, n_dim))
      H_xc_1 = (0d0,0d0)
      H_xc_1 = H_temp
      deallocate(H_temp)
   end subroutine calc_exch_firstord

   subroutine calc_left_pert_velo_mtx(self, k, derive_idx, eig_vec_mtx,eig_val, ret)
      implicit none
      class(hamil)                    :: self
      real(8), intent(in)             :: k(3),eig_val(:)
      integer   , intent(in)          :: derive_idx
      complex(8), intent(in)          :: eig_vec_mtx(:,:)
      complex(8), allocatable         :: ret(:,:), tmp(:,:),H_xc_1(:,:)
      integer                         :: n_dim
      n_dim = 2 * self%num_up
      allocate(tmp(n_dim, n_dim))
      tmp=(0d0,0d0)
      call self%calc_exch_firstord(eig_vec_mtx,eig_val,H_xc_1)
      call self%calc_velo_mtx(k,derive_idx,eig_vec_mtx,tmp)
      if(allocated(ret))then
         if(size(ret,1) /= n_dim .or. size(ret,2) /= n_dim) then
            deallocate(ret)
         endif
      endif
      if(.not. allocated(ret)) allocate(ret(n_dim, n_dim))
      ret=(0d0,0d0)
      call zgemm('C', 'N', n_dim, n_dim, n_dim, &
                  c_1, H_xc_1, n_dim, &
                  tmp, n_dim, &
                  c_0, ret, n_dim)
      deallocate(tmp)
      deallocate(H_xc_1)
   end subroutine calc_left_pert_velo_mtx

   subroutine calc_right_pert_velo_mtx(self, k, derive_idx, eig_vec_mtx,eig_val, ret)
      implicit none
      class(hamil)                    :: self
      real(8), intent(in)             :: k(3),eig_val(:)
      integer   , intent(in)          :: derive_idx
      complex(8), intent(in)          :: eig_vec_mtx(:,:)
      complex(8), allocatable         :: ret(:,:), tmp(:,:),H_xc_1(:,:)
      integer                         :: n_dim
      n_dim = 2 * self%num_up
      allocate(tmp(n_dim, n_dim))
      tmp=(0d0,0d0)
      call self%calc_exch_firstord(eig_vec_mtx,eig_val,H_xc_1)
      call self%calc_velo_mtx(k,derive_idx,eig_vec_mtx,tmp)
      if(allocated(ret))then
         if(size(ret,1) /= n_dim .or. size(ret,2) /= n_dim) then
            deallocate(ret)
         endif
      endif
      if(.not. allocated(ret)) allocate(ret(n_dim, n_dim))
      ret = (0d0,0d0)
      call zgemm('N', 'N', n_dim, n_dim, n_dim, &
                    c_1, tmp, n_dim,&
                    H_xc_1, n_dim,&
                    c_0, ret, n_dim)
      deallocate(tmp)
      deallocate(H_xc_1)
   end subroutine calc_right_pert_velo_mtx

   subroutine calc_velo_mtx(self, k, derive_idx, eig_vec_mtx, ret)
      implicit none
      class(hamil)                        :: self
      real(8), intent(in)             :: k(3)
      integer   , intent(in)          :: derive_idx
      complex(8), intent(in)          :: eig_vec_mtx(:,:)
      complex(8), allocatable         :: ret(:,:), tmp(:,:)
      integer                         :: n_dim, ierr(3) = 0
      n_dim = 2 * self%num_up
      allocate(tmp(n_dim, n_dim), stat=ierr(1))
      tmp=(0d0,0d0)
      if(.not. allocated(self%del_H)) allocate(self%del_H(n_dim, n_dim),stat=ierr(2))
      call self%set_derivative_k(k, derive_idx)
      call zgemm('N', 'N', n_dim, n_dim, n_dim, &
                 c_1, self%del_H, n_dim,&
                 eig_vec_mtx, n_dim,&
                 c_0, tmp, n_dim)
      deallocate(self%del_H)
      if(allocated(ret))then
         if(size(ret,1) /= n_dim .or. size(ret,2) /= n_dim) then
            deallocate(ret)
         endif
      endif
      if(.not. allocated(ret)) allocate(ret(n_dim, n_dim),stat=ierr(3))
      call check_ierr(ierr, me_in=self%me, msg=["failed alloc in calc_velo_mtx"])
      ret=(0d0,0d0)
      call zgemm('C', 'N', n_dim, n_dim, n_dim, &
                 c_1, eig_vec_mtx, n_dim, &
                 tmp, n_dim, &
                 c_0, ret, n_dim)

      deallocate(tmp)
   end subroutine calc_velo_mtx

   subroutine calc_eig_and_velo(self, k, eig_val, del_kx, del_ky,pert_log)
      implicit none
      class(hamil)             :: self
      real(8), intent(in)      :: k(3)
      real(8)                  :: eig_val(:)
      complex(8), allocatable  :: eig_vec(:,:), del_kx(:,:), del_ky(:,:), work(:), temp(:,:)
      real(8), allocatable     :: rwork(:)
      integer   , allocatable  :: iwork(:)
      integer, intent(in)      :: pert_log
      character (300)          :: elem_file!, folder = "/p/project/cjiff40/kipp1/output/spinspiral/100QSpirals/100_SW_PiHalf_5/"
      integer      :: n_dim, lwork, lrwork, liwork, info
      integer      :: ierr(3)
      n_dim = 2 * self%num_up
      if(.not. allocated(eig_vec)) allocate(eig_vec(n_dim,n_dim))
      if(.not. allocated(temp)) allocate(temp(n_dim,n_dim))
      eig_vec = (0d0,0d0)
      call self%setup_H(k, eig_vec)
      call calc_zheevd_size('V', eig_vec, eig_val, lwork, lrwork, liwork)
      temp = eig_vec
      allocate(work(lwork), stat=ierr(1))
      allocate(rwork(lrwork), stat=ierr(2))
      allocate(iwork(liwork), stat=ierr(3))
      call check_ierr(ierr, me_in=self%me, msg=[" tried to allocate in zheevd"])
      call zheevd('V', 'L', n_dim, eig_vec, n_dim, eig_val, &
                  work, lwork, rwork, lrwork, iwork, liwork, info)
      if(info /= 0) then
         write (*,*) "ZHEEVD in berry calculation failed, trying ZHEEV", self%me
         write (elem_file, "(A,I0.5,A)") "ham", self%me ,".npy"
         call save_npy(trim(self%prefix) //elem_file,temp)
         write (elem_file, "(A,I0.5,A)") "kpoint", self%me ,".npy"
         call save_npy(trim(self%prefix) //elem_file,k)
         eig_vec = temp
         call zheev('V', 'L', n_dim, eig_vec, n_dim, eig_val, &
                     work, lwork, rwork, info)
         if(info /= 0) then
            write (*,*) "ZHEEVD in berry calculation failed, trying ZHEEV", self%me
         endif
      endif
      deallocate(work, rwork, iwork)
      if     (pert_log==0) then
         call self%calc_velo_mtx(k, 1, eig_vec, del_kx)
         call self%calc_velo_mtx(k, 2, eig_vec, del_ky)
      else if(pert_log==1) then
         call self%calc_left_pert_velo_mtx(k, 1, eig_vec, eig_val, del_kx)
         call self%calc_velo_mtx(k, 2, eig_vec, del_ky)
      else if(pert_log==2) then
         call self%calc_right_pert_velo_mtx(k, 1, eig_vec, eig_val, del_kx)
         call self%calc_velo_mtx(k, 2, eig_vec, del_ky)
      else if(pert_log==3) then
         call self%calc_velo_mtx(k, 1, eig_vec, del_kx)
         call self%calc_left_pert_velo_mtx(k, 2, eig_vec, eig_val, del_ky)
      else if(pert_log==4) then
         call self%calc_velo_mtx(k, 1, eig_vec, del_kx)
         call self%calc_right_pert_velo_mtx(k, 2, eig_vec, eig_val, del_ky)
      endif
      deallocate(eig_vec)
   end subroutine calc_eig_and_velo

   subroutine calc_berry_z(self, z_comp, eig_val, x_mtx, y_mtx)
      implicit none
      class(hamil)             :: self
      real(8)                  :: eig_val(:), dE
      real(8)                  :: z_comp(:) !> \f$ \Omega^n_z \f$
      complex(8)               :: x_mtx(:,:), y_mtx(:,:)
      complex(8) :: fac
      integer    :: n_dim, n, m

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

   subroutine calc_berry_diag(self, z_comp, eig_val, x_mtx)
      implicit none
      class(hamil)             :: self
      real(8)                  :: z_comp(:), eig_val(:), fac, ferm !> \f$ \Omega^n_z \f$
      complex(8)               :: x_mtx(:,:)
      integer    :: n_dim, n, m, n_fermi
      n_dim = 2 * self%num_up
      z_comp =  0d0
      do n_fermi = 1,size(self%E_fermi)
         do n = 1,n_dim
            do m = 1,n_dim
               !if(n /= m) then
                  ferm  =  self%fermi_distr(eig_val(n), n_fermi)
                  call self%calc_fac_diag(eig_val(n), eig_val(m), self%E_fermi(n_fermi),fac)
                  z_comp(n_fermi) = z_comp(n_fermi) + 1d0/(Pi) *&
                              ferm * fac * real(x_mtx(n,m) * x_mtx(m,n))
               !endif
            enddo
         enddo
      enddo

   end subroutine calc_berry_diag

   subroutine calc_berry_diag_surf(self, z_comp, eig_val, x_mtx, y_mtx)
      implicit none
      class(hamil)             :: self
      real(8)                  :: z_comp(:), eig_val(:), fac, ferm !> \f$ \Omega^n_z \f$
      complex(8)               :: x_mtx(:,:), y_mtx(:,:)
      integer    :: n_dim, n, m, n_fermi
      n_dim = 2 * self%num_up
      z_comp =  0d0
      do n_fermi = 1,size(self%E_fermi)
         do n = 1,n_dim
            do m = 1,n_dim
               if(n /= m) then
                  ferm  =  self%fermi_distr(eig_val(n), n_fermi)
                  call self%calc_fac_surf(eig_val(n), eig_val(m), self%E_fermi(n_fermi),fac)
                  z_comp(n_fermi) = z_comp(n_fermi) + 1d0/(2d0*Pi) *&
                              ferm * fac * aimag(x_mtx(n,m) * y_mtx(m,n))
               endif
            enddo
         enddo
      enddo

   end subroutine calc_berry_diag_surf

   subroutine calc_berry_diag_sea(self, z_comp, eig_val, x_mtx, y_mtx)
      implicit none
      class(hamil)             :: self
      real(8)                  :: z_comp(:), eig_val(:), fac, ferm !> \f$ \Omega^n_z \f$
      complex(8)               :: x_mtx(:,:), y_mtx(:,:)
      integer    :: n_dim, n, m, n_fermi
      n_dim = 2 * self%num_up
      z_comp =  0d0
      do n_fermi = 1,size(self%E_fermi)
         do n = 1,n_dim
            do m = 1,n_dim
               if(n /= m) then
                  ferm  =  self%fermi_distr(eig_val(n), n_fermi)
                  call self%calc_fac_sea(eig_val(n), eig_val(m), self%E_fermi(n_fermi),fac)
                  z_comp(n_fermi) = z_comp(n_fermi) + 1d0/Pi *&
                              ferm * fac * aimag(x_mtx(n,m) * y_mtx(m,n))
               endif
            enddo
         enddo
      enddo
   
   end subroutine calc_berry_diag_sea

   subroutine calc_fac_diag(self, e_n, e_m, E_f, fac)
      implicit none
      class(hamil)             :: self
      real(8)                  :: e_n, e_m, dE, E_f, gamma
      real(8)                  :: fac!> \f$ \Omega^n_z \f$
   
      gamma = self%gamma
      fac =  0d0
      dE =  e_m - e_n
      fac =  gamma**2/(((E_f-e_n)**2+gamma**2)*((E_f-e_m)**2+gamma**2))
      fac =  0.5d0 * fac
   end subroutine calc_fac_diag

   subroutine calc_fac_surf(self, e_n, e_m, E_f, fac)
      implicit none
      class(hamil)             :: self
      real(8)                  :: e_n, e_m, dE, E_f, gamma
      real(8)                  :: fac!> \f$ \Omega^n_z \f$
      integer    :: n_dim
   
      gamma = self%gamma
      n_dim = 2 * self%num_up
      fac =  0d0
      dE =  e_m - e_n
      fac =  dE*gamma/(((E_f-e_n)**2+gamma**2)*((E_f-e_m)**2+gamma**2))
      fac = - 0.5d0 * fac
      !fac = - 1d0/(2d0*PI) * fac
   end subroutine calc_fac_surf
   
   subroutine calc_fac_sea(self, e_n, e_m, E_f, fac)
      implicit none
      class(hamil)             :: self
      real(8)                  :: e_n,e_m, dE, E_f, gamma, fac
      complex(8)               :: denom, numer
      integer    :: n_dim
   
      gamma = self%gamma
      n_dim = 2 * self%num_up
      fac =  0d0
      dE =  e_m - e_n
      denom = i_unit*gamma + E_f - e_n
      numer = i_unit*gamma + E_f - e_m
      fac =  gamma/(dE*((E_f-e_m)**2+gamma**2))&
             - dE**2/(dE**2 + eta_sq)**2*aimag(log(numer&
                                  /denom))
      !fac = 1d0/PI * fac
   end subroutine calc_fac_sea
   

   Subroutine  calc_eigenvalues(self, k_list, eig_val)
      Implicit None
      class(hamil)                      :: self
      real(8), intent(in)               :: k_list(:,:)
      real(8), allocatable,intent(out)  :: eig_val(:,:)
      real(8)                           :: k(3)
      complex(8), allocatable           :: H(:,:)
      integer    :: i, N, lwork, lrwork, liwork, info
      real(8), allocatable              :: RWORK(:)
      complex(8), allocatable           :: WORK(:)
      integer   , allocatable           :: IWORK(:)

      N =  2 * self%num_up
      allocate(eig_val(N, size(k_list, 2)))
      allocate(H(N,N))
      call calc_zheevd_size('N', H, eig_val(:,1), lwork, lrwork, liwork)
      allocate(RWORK(lrwork))
      allocate(IWORK(liwork))
      allocate(WORK(lwork))

      call MPI_Barrier(MPI_COMM_WORLD, info)
      do i = 1,size(k_list,2)
         k =  k_list(:,i)
         call self%setup_H(k, H)
         
         call zheevd('N', 'U', N, H, N, eig_val(:,i), WORK, lwork, &
                     RWORK, lrwork, IWORK, liwork, info)
         if( info /= 0) then
            write (*,*) "ZHEEVD failed: ", info
            call error_msg("Aborting now", abort=.True.)
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
      integer   , allocatable           :: iwork(:)
      integer                           :: N, lwork, lrwork, liwork, info

      N = 2 * self%num_up

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
