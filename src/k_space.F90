module Class_k_space
   use m_config
   use m_npy
   use mpi
   use Class_hamiltionian
   use Class_helper
   use MYPI
   use ieee_arithmetic
   implicit none

   type k_space
      real(8), allocatable :: new_k_pts(:,:), all_k_pts(:,:)
      real(8), allocatable :: int_DOS(:) !> integrated Density of states
      real(8), allocatable :: E_DOS(:)
      real(8), allocatable :: E_fermi(:) !> Fermi lvl
      real(8)              :: DOS_gamma !> broadening \f$ \Gamma \f$ used in
      !> DOS calculations
      real(8)              :: DOS_lower !> lower energy bound for DOS calc
      real(8)              :: DOS_upper !> upper energy bound for DOS calc
      real(8)              :: temp !> temperature used in fermi-dirac
      !> used in DOS calculations
      integer              :: DOS_num_k_pts !> number of kpts per dim
      integer              :: berry_num_k_pts !> number of ks per dim in berry calc
      integer              :: ACA_num_k_pts
      integer              :: num_DOS_pts!> number of points on E grid
      integer              :: num_k_pts !> number of k_pts per segment
      integer              :: num_plot_pts !> points per dimension
      integer              :: nProcs !> number of MPI Processes
      integer              :: me !> MPI rank
      integer              :: berry_iter !> number of grid refinements
      integer              :: kpts_per_step !> new kpts per step and Proc
      real(8)              :: k_shift(3) !> shift of brillouine-zone
      real(8)              :: berry_conv_crit !> convergance criterion for berry integration
      real(8), allocatable :: weights(:) !> weights for integration
      integer, allocatable :: elem_nodes(:,:) !> elements in triangulation
      real(8), allocatable :: refine_weights(:)
      real(8), allocatable :: k1_param(:) !> 1st k_space param
      real(8), allocatable :: k2_param(:) !> 2nd k_space param
      character(len=300)   :: filling, prefix, chosen_weights
      character(len=6)     :: ada_mode
      logical              :: perform_pad !> should the k-grid be padded, to match cores
      logical              :: calc_hall !> should hall conductivity be calculated
      logical              :: calc_orbmag !> should orbital magnetism be calculated
      logical              :: test_run !> should unit tests be performed
      !logical              :: pert_log !>should berry be calculated in first order perturbation theory
      type(hamil)          :: ham
      type(units)          :: units
   contains

      procedure :: vol_k_space_para       => vol_k_space_para
      procedure :: calc_and_print_band    => calc_and_print_band
      procedure :: setup_k_path_rel       => setup_k_path_rel
      procedure :: setup_k_path_abs       => setup_k_path_abs
      procedure :: setup_k_grid           => setup_k_grid
      procedure :: lorentzian             => lorentzian
      procedure :: calc_pdos              => calc_pdos
      procedure :: find_fermi             => find_fermi
      procedure :: set_fermi              => set_fermi
      procedure :: fermi_distr            => fermi_distr
      procedure :: write_fermi            => write_fermi
      procedure :: calc_and_print_dos     => calc_and_print_dos
      procedure :: calc_berry_quantities  => calc_berry_quantities
      procedure :: setup_inte_grid_para   => setup_inte_grid_para
      procedure :: setup_inte_grid_hex    => setup_inte_grid_hex
      procedure :: setup_berry_inte_grid  => setup_berry_inte_grid
      procedure :: set_weights_ksp        => set_weights_ksp
      procedure :: test_integration       => test_integration
      procedure :: Bcast_k_space          => Bcast_k_space
      procedure :: hex_border_x           => hex_border_x
      procedure :: vol_k_hex              => vol_k_hex
      procedure :: free_ksp               => free_ksp
      procedure :: find_E_max             => find_E_max
      procedure :: area_of_elem           => area_of_elem
      procedure :: centeroid_of_triang    => centeroid_of_triang
      procedure :: pad_k_points_init      => pad_k_points_init
      procedure :: new_pt                 => new_pt
      procedure :: on_hex_border          => on_hex_border
      procedure :: in_points              => in_points
      procedure :: append_kpts            => append_kpts
      procedure :: add_kpts_iter          => add_kpts_iter
      procedure :: random_pt_hex          => random_pt_hex
      procedure :: set_hall_weights       => set_hall_weights
      procedure :: set_orbmag_weights     => set_orbmag_weights
      procedure :: integrate_hall         => integrate_hall
      procedure :: integrate_orbmag       => integrate_orbmag
      procedure :: finalize_hall          => finalize_hall
      procedure :: finalize_hall_surf     => finalize_hall_surf
      procedure :: finalize_orbmag        => finalize_orbmag
      procedure :: process_hall           => process_hall
      procedure :: process_hall_surf      => process_hall_surf
      procedure :: process_orbmag         => process_orbmag
      procedure :: calc_new_berry_points  => calc_new_berry_points
      procedure :: calc_new_kidx          => calc_new_kidx
      procedure :: calc_orbmag_z_singleK  => calc_orbmag_z_singleK
      !procedure :: setup_A_mtx            => setup_A_mtx
      procedure :: save_grid              => save_grid
      procedure :: calc_ACA               => calc_ACA
      procedure :: calc_ACA_singleK       => calc_ACA_singleK
      procedure :: calc_S_singleK         => calc_S_singleK
      procedure :: calc_local_l           => calc_local_l
      procedure :: plot_omega_square      => plot_omega_square
   end type k_space

   interface
      subroutine run_triang(k_pts, ret_elem)
         real(8), intent(in)              :: k_pts(:,:)
         integer, allocatable             :: ret_elem(:,:)
      end subroutine run_triang
   end interface

contains
   subroutine free_ksp(self)
      implicit none
      class(k_space)              :: self

      if(allocated(self%new_k_pts)) deallocate(self%new_k_pts)
      if(allocated(self%int_DOS)) deallocate(self%int_DOS)
      if(allocated(self%E_DOS)) deallocate(self%E_DOS)
      if(allocated(self%E_fermi)) deallocate(self%E_fermi)
      if(allocated(self%k1_param)) deallocate(self%k1_param)
      if(allocated(self%k2_param)) deallocate(self%k2_param)
      call self%ham%free_ham()
   end subroutine free_ksp

   Subroutine  calc_and_print_band(self)
      use mpi
      Implicit None
      class(k_space)                :: self
      integer                       :: first, last, N
      integer                       :: send_count, ierr
      integer   , allocatable       :: num_elems(:), offsets(:)
      real(8), allocatable          :: eig_val(:,:), sec_eig_val(:,:), k_pts_sec(:,:)

      if(trim(self%filling) ==  "path_rel") then
         call self%setup_k_path_rel()
      else if(trim(self%filling) == "path_abs") then
         call self%setup_k_path_abs()
      else if(trim(self%filling) == "grid") then
         call self%setup_k_grid()
      else
         call error_msg("Filling not known", abort=.True.)
      endif

      call my_section(self%me, self%nProcs, size(self%new_k_pts, 2), first, last)
      allocate(k_pts_sec(3, last - first + 1))
      k_pts_sec = self%new_k_pts(:,first:last)

      call self%ham%calc_eigenvalues(k_pts_sec, sec_eig_val)

      N = 2 *  self%ham%num_up
      allocate(eig_val(N, size(self%new_k_pts,2)))
      allocate(num_elems(self%nProcs))
      allocate(offsets(self%nProcs))
      call sections(self%nProcs, size(self%new_k_pts, 2), num_elems, offsets)
      num_elems =  num_elems * N
      offsets   =  offsets   * N

      send_count =  N *  size(k_pts_sec, 2)

      call MPI_Gatherv(sec_eig_val, send_count, MPI_REAL8, &
                       eig_val,     num_elems,  offsets,   MPI_REAL8,&
                       root,        MPI_COMM_WORLD, ierr)

      if(self%me == root) then
         call save_npy(trim(self%prefix) //  "band_k.npy", self%new_k_pts / self%units%inv_length)
         call save_npy(trim(self%prefix) //  "band_E.npy", eig_val / self%units%energy)
         call save_npy(trim(self%prefix) //  "lattice.npy", &
                       self%ham%UC%lattice / self%units%length)
         call save_npy(trim(self%prefix) //  "rez_lattice.npy", &
                       self%ham%UC%rez_lattice / self%units%inv_length)
      endif

      deallocate(eig_val)
      deallocate(num_elems)
      deallocate(offsets)
      deallocate(sec_eig_val)
      deallocate(self%new_k_pts)
   End Subroutine calc_and_print_band

   subroutine calc_pdos(self, E, PDOS)
      use mpi
      implicit none
      class(k_space)          :: self
      real(8), intent(in)     :: E(:)
      real(8), intent(out)    :: PDOS(:,:)
      real(8), allocatable    :: RWORK(:), eig_val(:), loc_PDOS(:,:)
      real(8)                 :: lor
      complex(8), allocatable :: H(:,:), WORK(:)
      integer   , allocatable :: IWORK(:)
      integer     :: first, last, ierr
      integer     :: k_idx, E_idx, j, m, N, info
      integer     :: lwork, liwork, lrwork, percentage

      N =  2 * self%ham%num_up
      allocate(H(N,N))
      allocate(eig_val(N))

      call calc_zheevd_size('V', H, eig_val, lwork, lrwork, liwork)
      if(self%me ==  0) then
         write (*,*) "shape(H) =  ", shape(H)
         write (*,*) "lwork =  ", lwork
         write (*,*) "lrwork = ", lrwork
         write (*,*) "liwork = ", liwork
      endif

      allocate(WORK(lwork))
      allocate(RWORK(lrwork))
      allocate(IWORK(liwork))

      allocate(loc_PDOS(N, self%num_DOS_pts))

      loc_PDOS =  0d0
      PDOS = 0d0

      call my_section(self%me, self%nProcs, size(self%new_k_pts,2), first, last)
      if(self%me == root) write (*,*) "DOS kpts", size(self%new_k_pts,2)
      percentage = 0
      do k_idx=first, last
         if(self%me == root) then
            if(percentage /=  nint(10d0 * k_idx / (1d0 * last))) then
               percentage =  nint(10d0 * k_idx / (1d0 * last))
               write (*,"(I5,A)") 10*percentage, "%"
            endif
         endif
         call self%ham%setup_H(self%new_k_pts(:,k_idx), H)
         call zheevd('V', 'U', N, H, N, eig_val, WORK, lwork, &
                     RWORK, lrwork, IWORK, liwork, info)
         if( info /= 0) then
            call error_msg("ZHEEVD (with vectors) failed: ", abort=.True.)
         endif

!$OMP             parallel do private(j) default(shared)
         do m = 1,N
            do j = 1,N
               ! calc absolute of eigen_vectors
               H(j,m) =  H(j,m) *  conjg(H(j,m))
            enddo
         enddo

!$OMP             parallel do private(m,j,lor) default(shared)
         do E_idx =  1,self%num_DOS_pts
            ! eigenvectors are stored column-wise
            ! m-th eigenvalue
            ! j-th component of
            do m =  1,N
               lor = self%lorentzian(E(E_idx) - eig_val(m))
               do j = 1,N
                  loc_PDOS(j, E_idx) = loc_PDOS(j, E_idx) &
                                       + lor * H(j,m)
               enddo
            enddo
         enddo
      enddo
      loc_PDOS =  loc_PDOS / real(size(self%new_k_pts, 2))
      call MPI_Reduce(loc_PDOS,  PDOS,    size(loc_PDOS), &
                      MPI_REAL8, MPI_SUM, root, &
                      MPI_COMM_WORLD, ierr)

      deallocate(loc_PDOS)
      deallocate(WORK)
      deallocate(IWORK)
      deallocate(RWORK)
   end subroutine

   subroutine calc_and_print_dos(self)
      implicit none
      class(k_space)       :: self
      real(8), allocatable :: DOS(:), PDOS(:,:), up(:), down(:)
      real(8)              :: dE
      integer              :: i, num_up

      if(trim(self%ham%UC%uc_type) == "square_2d") then
         call self%setup_inte_grid_para(self%DOS_num_k_pts)
      elseif(trim(self%ham%UC%uc_type) == "file_square") then
         call self%setup_inte_grid_para(self%DOS_num_k_pts)
      elseif(trim(self%ham%UC%uc_type) == "honey_2d") then
         call self%setup_inte_grid_hex(self%DOS_num_k_pts)
      elseif(trim(self%ham%UC%uc_type) == "honey_line") then
         call self%setup_inte_grid_para(self%DOS_num_k_pts)
      else
         call error_msg("DOS k-grid not known", abort=.True.)
      endif

      num_up = self%ham%num_up
      allocate(PDOS(2*num_up, self%num_DOS_pts))

      call linspace(self%DOS_lower, self%DOS_upper, self%num_DOS_pts, self%E_DOS)

      if(self%me ==  root) then
         allocate(DOS(self%num_DOS_pts))
         allocate(up(self%num_DOS_pts))
         allocate(down(self%num_DOS_pts))
      endif

      call self%calc_pdos(self%E_DOS, PDOS)

      if(self%me == root) then
         DOS  = sum(PDOS,1)
         up   = sum(PDOS(1:num_up, :),1)
         down = sum(PDOS(num_up+1:2*num_up, :),1)

         call save_npy(trim(self%prefix) //  "DOS_E.npy", self%E_DOS / self%units%energy)
         ! unit of DOS is per energy
         call save_npy(trim(self%prefix) //  "DOS.npy", DOS * self%units%energy)
         call save_npy(trim(self%prefix) //  "DOS_partial.npy", PDOS * self%units%energy)
         call save_npy(trim(self%prefix) //  "DOS_up.npy", up *  self%units%energy)
         call save_npy(trim(self%prefix) //  "DOS_down.npy", down *  self%units%energy)

         allocate(self%int_DOS(self%num_DOS_pts))

         if(size(self%E_DOS) >=  2) then
            dE =  self%E_DOS(2) - self%E_DOS(1)
         else
            call error_msg("Can't perform integration. Only one point", &
                           abort=.True.)
         endif
         self%int_DOS(1) =  0d0
         do i =  2,size(self%E_DOS)
            self%int_DOS(i) =  self%int_DOS(i-1) &
                              + 0.5d0 * dE * (DOS(i-1) +  DOS(i))
         enddo
         ! integrated DOS ist unitless
         call save_npy(trim(self%prefix) // "DOS_integrated.npy", self%int_DOS)
      endif

      deallocate(self%new_k_pts)
      deallocate(PDOS)
      if(self%me == root) then
         deallocate(DOS)
         deallocate(up)
         deallocate(down)
      endif
   end subroutine calc_and_print_dos

   function init_k_space(cfg) result(self)
      use mpi
      implicit none
      type(k_space)         :: self
      type(CFG_t)           :: cfg
      real(8)               :: tmp
      !logical               :: logtmp
      integer               :: sz
      integer               :: ierr

      call MPI_Comm_size(MPI_COMM_WORLD, self%nProcs, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)

      self%units = init_units(cfg, self%me)
      self%ham   = init_hamil(cfg)

      if(self%me ==  0) then
         call CFG_get(cfg, "grid%k_shift", self%k_shift)

         call CFG_get(cfg, "output%band_prefix", self%prefix)
         if(self%prefix(len_trim(self%prefix):len_trim(self%prefix)) /=  "/") then
            self%prefix =  trim(self%prefix) // "/"
         endif
         call create_dir(self%prefix)

         call CFG_get(cfg, "band%filling", self%filling)

         call CFG_get(cfg, "dos%delta_broadening", tmp)
         self%DOS_gamma =  tmp * self%units%energy

         call CFG_get(cfg, "dos%num_points", self%num_DOS_pts)

         call self%ham%UC%save_unit_cell(trim(self%prefix))

         call CFG_get_size(cfg, "band%k_x", sz)
         allocate(self%k1_param(sz))
         call CFG_get_size(cfg, "band%k_y", sz)
         allocate(self%k2_param(sz))

         call CFG_get(cfg, "band%k_x", self%k1_param)
         call CFG_get(cfg, "band%k_y", self%k2_param)
         call CFG_get(cfg, "band%num_points", self%num_k_pts)

         call CFG_get(cfg, "dos%k_pts_per_dim", self%DOS_num_k_pts)

         call CFG_get(cfg, "dos%lower_E_bound", tmp)
         self%DOS_lower =  tmp * self%units%energy

         call CFG_get(cfg, "dos%upper_E_bound", tmp)
         self%DOS_upper =  tmp * self%units%energy
         
         !call CFG_get(cfg,"berry%pert_log",self%pert_log)

         call CFG_get(cfg, "berry%k_pts_per_dim", self%berry_num_k_pts)
         call CFG_get(cfg, "berry%temperature", tmp)
         self%temp = tmp * self%units%temperature

         call CFG_get(cfg, "berry%refinement_iter", self%berry_iter)
         call CFG_get(cfg, "berry%kpts_per_step", self%kpts_per_step)
         call CFG_get(cfg, "berry%conv_criterion", self%berry_conv_crit)
         call CFG_get(cfg, "berry%perform_pad", self%perform_pad)
         call CFG_get(cfg, "berry%calc_hall", self%calc_hall)
         call CFG_get(cfg, "berry%calc_orbmag", self%calc_orbmag)
         call CFG_get(cfg, "berry%weights", self%chosen_weights)
         call CFG_get(cfg, "berry%adaptive_mode", self%ada_mode)

         call CFG_get(cfg, "ACA%num_kpts", self%ACA_num_k_pts)
         call CFG_get(cfg, "plot%num_plot_points", self%num_plot_pts)

         call CFG_get(cfg, "general%test_run", self%test_run)
      endif
      call self%Bcast_k_space()
   end function init_k_space

   subroutine Bcast_k_space(self)
      use mpi
      class(k_space)             :: self
      integer, parameter     :: num_cast =  26
      integer                :: ierr(num_cast)
      integer                :: sz(2)
      ierr =  0

      if(self%me ==  root) then
         sz(1) = size(self%k1_param)
         sz(2) = size(self%k2_param)
      endif

      call MPI_Bcast(self%prefix,  300,          MPI_CHARACTER, &
                     root,                    MPI_COMM_WORLD, ierr(1))
      call MPI_Bcast(self%filling, 300,          MPI_CHARACTER, &
                     root,                    MPI_COMM_WORLD, ierr(2))

      call MPI_Bcast(self%DOS_gamma,   1,            MPI_REAL8,   &
                     root,                        MPI_COMM_WORLD, ierr(3))
      call MPI_Bcast(self%num_DOS_pts, 1,            MYPI_INT, &
                     root,                        MPI_COMM_WORLD, ierr(4))

      call MPI_Bcast(sz, 2,            MYPI_INT, &
                     root,          MPI_COMM_WORLD, ierr(5))

      if(self%me /= root) then
         allocate(self%k1_param(sz(1)))
         allocate(self%k2_param(sz(2)))
      endif
      call MPI_Bcast(self%k1_param,      sz(1),          MPI_REAL8,    &
                     root,               MPI_COMM_WORLD, ierr(6))
      call MPI_Bcast(self%k2_param,      sz(2),          MPI_REAL8,    &
                     root,               MPI_COMM_WORLD, ierr(7))
      call MPI_Bcast(self%num_k_pts,     1,              MYPI_INT, &
                     root,               MPI_COMM_WORLD, ierr(8))
      call MPI_Bcast(self%DOS_num_k_pts, 1,              MYPI_INT, &
                     root,               MPI_COMM_WORLD, ierr(9))
      call MPI_Bcast(self%DOS_lower,     1,              MPI_REAL8, &
                     root,               MPI_COMM_WORLD, ierr(10))
      call MPI_Bcast(self%DOS_upper,     1,              MPI_REAL8, &
                     root,               MPI_COMM_WORLD, ierr(11))

      ! Berry parameter               
      call MPI_Bcast(self%berry_num_k_pts, 1,            MYPI_INT,   &
                     root,                            MPI_COMM_WORLD, ierr(12))
      call MPI_Bcast(self%temp,            1,            MPI_REAL8,     &
                     root,                            MPI_COMM_WORLD, ierr(13))
      call MPI_Bcast(self%berry_iter,      1,            MYPI_INT,   &
                     root,                            MPI_COMM_WORLD, ierr(14))
      call MPI_Bcast(self%kpts_per_step,   1,            MYPI_INT,   &
                     root,                            MPI_COMM_WORLD, ierr(15))
      call MPI_Bcast(self%k_shift,         3,            MPI_REAL8,     &
                     root,                            MPI_COMM_WORLD, ierr(16))
      call MPI_Bcast(self%berry_conv_crit, 1,            MPI_REAL8,     &
                     root,                            MPI_COMM_WORLD, ierr(17))
      call MPI_Bcast(self%perform_pad,     1,            MPI_LOGICAL,   &
                     root,                            MPI_COMM_WORLD, ierr(18))
      call MPI_Bcast(self%calc_hall,       1,            MPI_LOGICAL,   &
                     root,                            MPI_COMM_WORLD, ierr(19))
      call MPI_Bcast(self%calc_orbmag,     1,            MPI_LOGICAL,   &
                     root,                            MPI_COMM_WORLD, ierr(20))
      call MPI_Bcast(self%chosen_weights,  300,          MPI_CHARACTER, &
                     root,                            MPI_COMM_WORLD, ierr(21))
      call MPI_Bcast(self%ada_mode,    6,          MPI_CHARACTER, &
                     root,                            MPI_COMM_WORLD, ierr(22))

      call MPI_Bcast(self%test_run,      1,              MPI_LOGICAL, &
                     root,              MPI_COMM_WORLD, ierr(23))
      call MPI_Bcast(self%ACA_num_k_pts, 1,              MYPI_INT,    &
                     root,              MPI_COMM_WORLD, ierr(24))
      call MPI_Bcast(self%num_plot_pts,  1,              MYPI_INT,    &
                     root,              MPI_COMM_WORLD, ierr(25))
      !call MPI_Bcast(self%pert_log, 1,            MPI_LOGICAL,   &
      !               root,                            MPI_COMM_WORLD, ierr(26))

      call check_ierr(ierr, self%me, "Ksp Bcast")
   end subroutine Bcast_k_space

   subroutine setup_k_grid(self)
      implicit none
      class(k_space)       :: self
      real(8), allocatable :: kx_points(:), ky_points(:)
      real(8), allocatable :: kx_grid(:,:), ky_grid(:,:)
      integer              :: sz_x, sz_y, i, j

      sz_x =  NINT(self%k1_param(3))
      sz_y =  NINT(self%k2_param(3))

      self%k1_param(1:2) =  self%k1_param(1:2) * self%units%inv_length
      self%k2_param(1:2) =  self%k2_param(1:2) * self%units%inv_length

      allocate(kx_grid(sz_x, sz_y))
      allocate(ky_grid(sz_x, sz_y))
      allocate(kx_points(sz_x))
      allocate(ky_points(sz_y))

      call linspace(self%k1_param(1), self%k1_param(2), sz_x , kx_points)
      call linspace(self%k2_param(1), self%k2_param(2), sz_y, ky_points)

      do j =  0,sz_y-1
         kx_grid(:,j+1) =  kx_points
      enddo

      do i =  0,sz_x-1
         ky_grid(i+1,:) =  ky_points
      enddo

      allocate(self%new_k_pts(3, sz_x* sz_y))
      self%new_k_pts(1,:) =  reshape(kx_grid, (/sz_x * sz_y /))
      self%new_k_pts(2,:) =  reshape(ky_grid, (/sz_x * sz_y /))
      self%new_k_pts(3,:) =  0.0d0

      deallocate(kx_grid)
      deallocate(ky_grid)
      deallocate(kx_points)
      deallocate(ky_points)

   end subroutine setup_k_grid

   subroutine setup_inte_grid_para(self, n_k, padding)
      implicit none
      class(k_space)        :: self
      integer, intent(in):: n_k
      real(8), allocatable  :: ls(:)
      real(8)               :: k1(3), k2(3)
      integer               :: i, j, cnt
      logical, optional     :: padding

      if(allocated(self%new_k_pts)) deallocate(self%new_k_pts)
      allocate(self%new_k_pts(3,n_k**2))

      k1 =  0d0
      k2 =  0d0
      k1(1:2) =  self%ham%UC%rez_lattice(:,1)
      k2(1:2) =  self%ham%UC%rez_lattice(:,2)

      call linspace(0d0, 1d0, n_k, ls)

      self%new_k_pts = 0d0

      cnt        = 1
      do i = 1,n_k
         do j =  1,n_k
            self%new_k_pts(:,cnt) = ls(i) *  k1 +  ls(j) *  k2
            cnt =  cnt + 1
         enddo
      enddo
      if(.not. present(padding) .and. self%perform_pad) then
         call run_triang(self%new_k_pts, self%elem_nodes)
         call self%pad_k_points_init()
      elseif(present(padding)) then
         if(padding) then
            call self%pad_k_points_init()
            call run_triang(self%new_k_pts, self%elem_nodes)
         endif
      endif

      forall(i = 1:size(self%new_k_pts,2)) self%new_k_pts(:,i) = &
         self%new_k_pts(:,i) + my_norm2(k1) * self%k_shift
   end subroutine setup_inte_grid_para

   subroutine setup_inte_grid_hex(self, n_k)
      implicit none
      class(k_space)         :: self
      integer, intent(in) :: n_k
      real(8), allocatable   :: x(:), y(:)
      real(8)                :: den, l, a
      integer                :: cnt_k, start, halt, my_n, i

      l = my_norm2(self%ham%UC%rez_lattice(:,1))
      a = l / (2.0 * cos(deg_30))
      den =  (1.0*n_k)/a

      call linspace(-0.5*(l), 0.5*l, nint(den * l),y)
      cnt_k =  0
      do i =  1,size(y)
         cnt_k =  cnt_k + nint(2 * self%hex_border_x(y(i)) * den)
      enddo

      if(allocated(self%new_k_pts)) deallocate(self%new_k_pts)
      allocate(self%new_k_pts(3, cnt_k))

      self%new_k_pts = 0d0

      start = 1
      do i =  1,size(y)
         my_n = nint(2*self%hex_border_x(y(i)) * den)
         call linspace(-self%hex_border_x(y(i)), self%hex_border_x(y(i)), my_n, x)

         halt = start + size(x) - 1

         self%new_k_pts(1,start:halt) = x
         self%new_k_pts(2,start:halt) = y(i)

         start = halt + 1
      enddo

      call run_triang(self%new_k_pts, self%elem_nodes)
      if(self%perform_pad) call self%pad_k_points_init()

      forall(i = 1:size(self%new_k_pts,2)) self%new_k_pts(:,i) = &
         self%new_k_pts(:,i) + l * self%k_shift
   end subroutine setup_inte_grid_hex

   subroutine test_integration(self, iter)
      implicit none
      class(k_space)     :: self
      integer            :: i, iter
      real(8)            :: integral, kpt(3), f
      character(len=300) :: filename

      call run_triang(self%all_k_pts, self%elem_nodes)
      call self%set_weights_ksp()
      integral =  0d0
      do i =  1, size(self%all_k_pts,2)
         kpt =  self%all_k_pts(:,i)
         f = test_func(kpt(1:2))
         integral =  integral + self%weights(i) * f
      enddo
      if(self%me ==  root) write (*,*) iter, "integration =  ", integral

      if(self%me == root) then
         write (filename, "(A, I0.5, A)") "k_points_", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), self%all_k_pts)

         write (filename, "(A, I0.5, A)") "elems_", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), self%elem_nodes)

         write (filename, "(A, I0.5, A)") "integral_", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), [integral])
      endif
   end subroutine test_integration

   subroutine set_weights_ksp(self)
      implicit none
      class(k_space)   :: self
      real(8)          :: A_proj, vec1(3), vec2(3)
      integer          :: i, j, k_idx

      if(allocated(self%weights)) then
         deallocate(self%weights)
      endif
      allocate(self%weights(size(self%all_k_pts,2)))

      self%weights = 0d0
      do i = 1,size(self%elem_nodes,1)
         vec1 = self%all_k_pts(:,self%elem_nodes(i,1)) - self%all_k_pts(:,self%elem_nodes(i,2))
         vec2 = self%all_k_pts(:,self%elem_nodes(i,1)) - self%all_k_pts(:,self%elem_nodes(i,3))
         A_proj =  0.1666666666666d0 * my_norm2(cross_prod(vec1, vec2))
         do j =  1,3
            k_idx = self%elem_nodes(i,j)
            self%weights(k_idx) = self%weights(k_idx) +  A_proj
         enddo
      enddo
   end subroutine set_weights_ksp

   function hex_border_x(self, y) result(x)
      implicit none
      class(k_space)      :: self
      real(8), intent(in) :: y
      real(8)             :: m, b, l, a, x

      l = my_norm2(self%ham%UC%rez_lattice(:,1))
      a = l / (2.0 * cos(deg_30))
      m =  - 2.0 * sin(deg_60)
      b =  - m * a

      x = (abs(y)-b)/m
   end function hex_border_x

   subroutine setup_k_path_abs(self)
      implicit none
      class(k_space)        :: self
      integer               :: n_pts, n_sec, start, halt, i
      real(8), allocatable  :: tmp(:)

      self%k1_param =  self%k1_param * self%units%inv_length
      self%k2_param =  self%k2_param * self%units%inv_length

      n_pts =  self%num_k_pts
      n_sec =  size(self%k1_param) - 1

      allocate(self%new_k_pts(3, n_sec * (self%num_k_pts-1) + 1))
      self%new_k_pts(3,:) =  0d0

      start = 1
      do i =  1,n_sec
         halt =  start +  n_pts - 1

         call linspace(self%k1_param(i), self%k1_param(i+1), n_pts, tmp)
         self%new_k_pts(1,start:halt) = tmp
         call linspace(self%k2_param(i), self%k2_param(i+1), n_pts, tmp)
         self%new_k_pts(2,start:halt) = tmp
         start =  halt
      enddo

   end subroutine setup_k_path_abs

   subroutine setup_k_path_rel(self)
      implicit none
      class(k_space)        :: self
      real(8), allocatable :: c1_sec(:), c2_sec(:)
      real(8), dimension(3)              :: k1, k2
      integer    ::  n_pts, n_sec,i,j, start, halt, cnt

      n_sec =  size(self%k1_param) - 1
      n_pts =  self%num_k_pts

      allocate(self%new_k_pts(3,n_sec * (n_pts - 1) + 1))
      allocate(c1_sec(n_sec))
      allocate(c2_sec(n_sec))

      k1(1:2) =  self%ham%UC%rez_lattice(:,1)
      k1(3)   =  0d0
      k2(1:2) =  self%ham%UC%rez_lattice(:,2)
      k2(3)   =  0d0

      write (*,*) "k1 = ", k1
      write (*,*) "k2 = ", k2

      start =  1
      do i =  1, n_sec
         halt =  start +  n_pts - 1
         ! linear combination of k-vec in current
         ! section.

         call linspace(self%k1_param(i),self% k1_param(i+1), n_pts, c1_sec)
         call linspace(self%k2_param(i), self%k2_param(i+1), n_pts, c2_sec)

         cnt =  1
         do j =  start,halt
            self%new_k_pts(:,j) =  c1_sec(cnt) *  k1 +  c2_sec(cnt) *  k2
            cnt =  cnt + 1
         enddo
         start =  halt
      enddo

   end subroutine setup_k_path_rel

   function vol_k_space_para(self) result(vol)
      implicit none
      class(k_space)         :: self
      real(8)                :: vol, k1(3), k2(3)

      k1      = 0d0
      k2      = 0d0
      k1(1:2) = self%ham%UC%rez_lattice(:,1)
      k2(1:2) = self%ham%UC%rez_lattice(:,2)

      vol = my_norm2(cross_prod(k1,k2))
   end function vol_k_space_para

   function vol_k_hex(self) result(vol)
      implicit none
      class(k_space)    :: self
      real(8)           :: a, l, vol

      l = my_norm2(self%ham%UC%rez_lattice(:,1))
      a = l / (2.0 * cos(deg_30))

      vol =  3d0 * a * a * sin(deg_60)
   end function vol_k_hex

   subroutine calc_berry_quantities(self,pert_log)
      use mpi
      implicit none
      class(k_space)          :: self
      real(8), allocatable    :: eig_val_all(:,:), eig_val_new(:,:),&
                                 hall(:), hall_old(:), omega_z_all(:,:), omega_z_new(:,:),&
                                 hall_surf(:), hall_surf_old(:), omega_surf_all(:,:), omega_surf_new(:,:),&
                                 hall_sea(:), hall_sea_old(:), omega_sea_all(:,:), omega_sea_new(:,:),&
                                 orbmag(:), orbmag_old(:), Q_L_all(:,:), Q_IC_all(:,:), &
                                 Q_L_new(:,:), Q_IC_new(:,:), orbmag_L(:), orbmag_IC(:)
      real(8)                  :: start, factor
      integer   , allocatable  :: kidx_all(:), kidx_new(:)
      integer     :: N_k, num_up, iter, n_ferm,nProcs
      integer     :: all_err(13), info
      character(len=300)       :: msg, surf_name = "hall_cond_surf", sea_name = "hall_cond_sea"
      logical                  :: done_hall = .True., done_orbmag = .True.,  done_hall_surf = .True.,  done_hall_sea = .True.
      logical, intent(in)      :: pert_log
      call self%setup_berry_inte_grid()
      N_k = size(self%new_k_pts, 2)
      num_up =  self%ham%num_up
      n_ferm = size(self%E_fermi)
      nProcs = self%nProcs
      all_err = 0

      allocate(self%ham%del_H(2*num_up, 2*num_up), stat=all_err(1))
      allocate(hall_old(size(self%E_fermi)),       stat=all_err(2))
      allocate(hall(size(self%E_fermi)),           stat=all_err(3))
      allocate(hall_surf_old(size(self%E_fermi)),       stat=all_err(2))
      allocate(hall_surf(size(self%E_fermi)),           stat=all_err(3))
      allocate(hall_sea_old(size(self%E_fermi)),       stat=all_err(2))
      allocate(hall_sea(size(self%E_fermi)),           stat=all_err(3))
      allocate(orbmag_old(size(self%E_fermi)),     stat=all_err(4))
      allocate(orbmag(size(self%E_fermi)),         stat=all_err(5))
      allocate(orbmag_L(size(self%E_fermi)),       stat=all_err(6))
      allocate(orbmag_IC(size(self%E_fermi)),      stat=all_err(7))
      allocate(eig_val_all(2*num_up, 0),           stat=all_err(8))
      allocate(omega_z_all(2*num_up, 0),           stat=all_err(9))
      allocate(omega_surf_all(2*num_up, 0),           stat=all_err(9))
      allocate(omega_sea_all(2*num_up, 0),           stat=all_err(9))
      allocate(Q_L_all(n_ferm, 0),                 stat=all_err(10))
      allocate(Q_IC_all(n_ferm, 0),                stat=all_err(11))
      allocate(kidx_all(0),                        stat=all_err(12))
      allocate(self%all_k_pts(3,0),                stat=all_err(13))

      hall    =  1e35
      hall_surf    =  1e35
      hall_sea    =  1e35
      orbmag =  1e35
      call check_ierr(all_err, self%me, "allocate hall vars")

      iter =  1
      start = MPI_Wtime()
      do iter =1,self%berry_iter
         if(self%me == root) write (*,*) "Time: ", MPI_Wtime() -  start
         call self%calc_new_berry_points(eig_val_new, omega_z_new, omega_surf_new, omega_sea_new, Q_L_new, Q_IC_new,pert_log)
         call self%calc_new_kidx(kidx_new)

         ! concat to old ones
         call append_kidx(kidx_all, kidx_new)
         call self%append_kpts()
         call append_eigval(eig_val_all, eig_val_new)
         if(self%me ==  root) write (*,*) self%me, "post appending"

         if(self%calc_hall)   call append_quantitiy(omega_z_all, omega_z_new)
         if(self%calc_hall)   call append_quantitiy(omega_surf_all, omega_surf_new)
         if(self%calc_hall)   call append_quantitiy(omega_sea_all, omega_sea_new)
         if(self%calc_orbmag) then
            call append_quantitiy(Q_L_all, Q_L_new)
            call append_quantitiy(Q_IC_all, Q_IC_new)
         endif
         if(self%calc_hall) then
            hall_old = hall
            hall_surf_old = hall_surf
            hall_sea_old = hall_sea
            call self%integrate_hall(kidx_all, omega_z_all, eig_val_all, hall)
            call self%integrate_hall(kidx_all, omega_surf_all, eig_val_all, hall_surf)
            call self%integrate_hall(kidx_all, omega_sea_all, eig_val_all, hall_sea)

            ! save current iteration and check if converged
            done_hall =  self%process_hall(hall, hall_old, iter, omega_z_all)
            done_hall_surf =  self%process_hall_surf(hall_surf, hall_surf_old, iter, omega_surf_all,surf_name)
            done_hall_sea =  self%process_hall_surf(hall_sea, hall_sea_old, iter, omega_sea_all,sea_name)
         endif
         if(done_hall .and. trim(self%chosen_weights) == "hall") then
            call error_msg("Switched to orbmag-weights", &
                           p_color=c_green, abort=.False.)
            self%chosen_weights = "orbmag"
         endif
         if(self%calc_orbmag) then
            orbmag_old = orbmag
            call self%integrate_orbmag(kidx_all, Q_L_all, Q_IC_all, orbmag, orbmag_L, orbmag_IC)

            ! save current iteration and check if converged
            factor = self%ham%UC%calc_area() / self%units%mag_dipol
            done_orbmag = self%process_orbmag(orbmag*factor, orbmag_old*factor, &
                                              orbmag_L*factor, orbmag_IC*factor,&
                                              iter)
         endif

         if(done_orbmag .and. trim(self%chosen_weights) == "orbmag") then
            call error_msg("Switched to hall-weights", &
                           p_color=c_green, abort=.False.)
            self%chosen_weights = "hall"
         endif

         ! Stop if both converged
         if(done_hall .and. done_orbmag) exit
         if(trim(self%chosen_weights) == "hall")then
            if(.not.self%calc_hall) then
               call error_msg("Must calculate hall to use it as weights", abort=.True.)
            endif

            call self%set_hall_weights(omega_z_all, kidx_all)
         elseif(trim(self%chosen_weights) == "orbmag")then
            if(.not.self%calc_orbmag) then
               call error_msg("Must calculate orbmag to use it as weights", abort=.True.)
            endif

            call self%set_orbmag_weights(Q_L_all +  Q_IC_all, kidx_all)
         else
            call error_msg("weights unknown", abort=.True.)
         endif
         call save_grid(self,iter)
         call self%add_kpts_iter(self%kpts_per_step*self%nProcs, self%new_k_pts)
      enddo
      if(self%calc_hall) then
         call self%finalize_hall(hall,omega_z_all)
         call self%finalize_hall_surf(hall_surf,omega_surf_all,surf_name)
         call self%finalize_hall_surf(hall_sea,omega_sea_all,sea_name)
      endif
      if(self%calc_orbmag) call self%finalize_orbmag(orbmag, orbmag_L, orbmag_IC)

      if(allocated(self%new_k_pts)) deallocate(self%new_k_pts)
      deallocate(self%ham%del_H, hall_old, eig_val_all, omega_z_all, &
                 hall_surf, hall_surf_old, omega_surf_all, kidx_all, self%all_k_pts, &
                 hall_sea, hall_sea_old, omega_sea_all,hall, stat=info, errmsg=msg)

   end subroutine calc_berry_quantities

   subroutine calc_new_kidx(self, kidx_new)
      implicit none
      class(k_space)               :: self
      integer   , allocatable      :: kidx_new(:)
      integer                      :: cnt, N_k, k_idx
      integer                      :: first, last
      integer                      :: err(1)

      N_k = size(self%new_k_pts, 2)
      call my_section(self%me, self%nProcs, N_k, first, last)
      allocate(kidx_new(last-first+1),      stat=err(1))
      call check_ierr(err, self%me, "new kidx")

      cnt =  1
      do k_idx = first, last
         kidx_new(cnt) =  k_idx + size(self%all_k_pts,2)
         cnt = cnt + 1
      enddo
   end subroutine calc_new_kidx

   subroutine calc_new_berry_points(self, eig_val_new, omega_z_new, omega_surf_new, omega_sea_new, Q_L_new, Q_IC_new, pert_log)
      use mpi
      implicit none
      class(k_space)            :: self
      integer                   :: N_k, cnt, k_idx, num_up, n_ferm,pert_idx
      integer                   :: first, last, err(3), me, ierr
      real(8)                   :: tmp
      real(8)                   :: k(3)
      real(8), allocatable      :: eig_val_new(:,:), omega_z_new(:,:), omega_surf_new(:,:),&
                                   omega_sea_new(:,:),omega_z_pert_new(:), Q_L_new(:,:),&
                                   Q_IC_new(:,:)
      complex(8), allocatable   :: del_kx(:,:), del_ky(:,:)
      logical, intent(in)       :: pert_log
      tmp = 0d0
      N_k = size(self%new_k_pts, 2)
      n_ferm =  size(self%E_fermi)
      num_up =  self%ham%num_up

      call my_section(self%me, self%nProcs, N_k, first, last)

      err =  0
      allocate(eig_val_new(2*num_up, last-first+1), stat=err(1))
      if(self%calc_hall)   allocate(omega_z_new(2*num_up, last-first+1), stat=err(2))
      if(self%calc_hall)   allocate(omega_surf_new(2*num_up, last-first+1), stat=err(2))
      if(self%calc_hall)   allocate(omega_sea_new(2*num_up, last-first+1), stat=err(2))
      if(self%calc_orbmag) allocate(Q_L_new(n_ferm,        last-first+1), stat=err(3))
      if(self%calc_orbmag) allocate(Q_IC_new(n_ferm,        last-first+1), stat=err(3))

      call check_ierr(err, self%me, " new chunk alloc")
      call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
      ! calculate
      cnt =  1
      if(pert_log) then
         do k_idx = first, last
            k = self%new_k_pts(:,k_idx)
            call self%ham%calc_eig_and_velo(k, eig_val_new(:,cnt), del_kx, del_ky,0)
            if(self%calc_hall) then
               call self%ham%calc_berry_z(omega_z_new(:,cnt),&
                                       eig_val_new(:,cnt), del_kx, del_ky)
               if(allocated(omega_z_pert_new)) deallocate(omega_z_pert_new)
               allocate(omega_z_pert_new(2*num_up), stat=err(2))
               omega_z_pert_new=0d0
               !allocation for omega_z_pert_new can be done here, since ham%calc_berry_z sets to zero
               do pert_idx=1,4
                  !write (*,*) "calc_new_berrypoints", me, pert_log, pert_idx
                  if(allocated(del_kx)) deallocate(del_kx)
                  if(allocated(del_ky)) deallocate(del_ky)
                  call self%ham%calc_eig_and_velo(k, eig_val_new(:,cnt), del_kx, del_ky,pert_idx)
                  call self%ham%calc_berry_z(omega_z_pert_new,&
                                             eig_val_new(:,cnt), del_kx, del_ky)
                  omega_z_new(:,cnt) = omega_z_new(:,cnt) + omega_z_pert_new
               enddo
            endif
            if(self%calc_orbmag) then
               call self%calc_orbmag_z_singleK(Q_L_new(:,cnt), Q_IC_new(:,cnt), &
                                            eig_val_new(:,cnt), del_kx, del_ky)
            endif
            cnt = cnt + 1
         enddo
      else if(.not. pert_log) then
         do k_idx = first, last
            k = self%new_k_pts(:,k_idx)
            call self%ham%calc_eig_and_velo(k, eig_val_new(:,cnt), del_kx, del_ky,0)
         
            if(self%calc_hall) then
               call self%ham%calc_berry_z(omega_z_new(:,cnt),&
                                       eig_val_new(:,cnt), del_kx, del_ky)
               call self%ham%calc_berry_diag_surf(omega_surf_new(:,cnt),&
                                       eig_val_new(:,cnt), del_kx)
               call self%ham%calc_berry_diag_sea(omega_sea_new(:,cnt),&
                                       eig_val_new(:,cnt), del_kx)
               !call self%ham%calc_berry_z(omega_xx_new(:,cnt),&
               !                        eig_val_new(:,cnt), del_kx, del_kx)
            endif
         
            if(self%calc_orbmag) then
               call self%calc_orbmag_z_singleK(Q_L_new(:,cnt), Q_IC_new(:,cnt), &
                                            eig_val_new(:,cnt), del_kx, del_ky)
            endif
            cnt = cnt + 1
         enddo
      endif
      if(allocated(del_kx)) deallocate(del_kx)
      if(allocated(del_ky)) deallocate(del_ky)
      if(allocated(omega_z_pert_new)) deallocate(omega_z_pert_new)
   end subroutine calc_new_berry_points

   function process_hall(self, var, var_old, iter, varall) result(cancel)
      implicit none
      class(k_space)                 :: self
      real(8), intent(in)            :: var(:), var_old(:), varall(:,:)
      integer   , intent(in)         :: iter
      character(len=*), parameter    :: var_name = "hall_cond"
      character(len=300)             :: filename
      logical                        :: cancel
      real(8)                        :: rel_error

      cancel = .False.

      if(self%me == root) then
         if(any(ieee_is_nan(var))) then
            write (*,*) "hall is nan"
         endif
         if(any(ieee_is_nan(var_old))) then
            write (*,*) "hall_old is nan"
         endif
      endif

      ! save current iteration data
      if(self%me == root) then
         write (filename, "(A,I0.5,A)") trim(var_name) // "_iter=", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), var)

         call save_npy(trim(self%prefix) // trim(var_name) //  "_E.npy", &
                       self%E_fermi / self%units%energy)
         if (iter == self%berry_iter) then
            call save_npy(trim(self%prefix) // "unitcell_"// trim(filename), varall)
         endif
      endif

      ! check for convergence
      rel_error = my_norm2(var - var_old) &
                  / (self%kpts_per_step * self%nProcs * my_norm2(var))!/ (1d0*size(var))

      if(self%me == root) then
         write (*,"(I5,A,A,A,I7,A,ES10.3)") iter, " var: ", var_name, &
            " nkpts ", size(self%all_k_pts,2),&
            " err ", rel_error
      endif

      if(rel_error < self%berry_conv_crit) then
         if(self%me == root) write (*,*) "Converged " // trim(var_name) //   " interation"
         cancel = .True.
      endif
   end function process_hall
   function process_hall_surf(self, var, var_old, iter, varall, var_name) result(cancel)
      implicit none
      class(k_space)                 :: self
      real(8), intent(in)            :: var(:), var_old(:), varall(:,:)
      integer   , intent(in)         :: iter
      character(len=*), intent(in)   :: var_name
      character(len=300)             :: filename
      logical                        :: cancel
      real(8)                        :: rel_error

      cancel = .False.

      if(self%me == root) then
         if(any(ieee_is_nan(var))) then
            write (*,*) "hall is nan"
         endif
         if(any(ieee_is_nan(var_old))) then
            write (*,*) "hall_old is nan"
         endif
      endif

      ! save current iteration data
      if(self%me == root) then
         write (filename, "(A,I0.5,A)") trim(var_name) // "_iter=", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), var)

         call save_npy(trim(self%prefix) // trim(var_name) //  "_E.npy", &
                       self%E_fermi / self%units%energy)
         if (self%ham%UC%num_atoms==2) then
            call save_npy(trim(self%prefix) // "unitcell_"// trim(filename), varall)
         endif
      endif

      ! check for convergence
      rel_error = my_norm2(var - var_old) &
                  / (self%kpts_per_step * self%nProcs * my_norm2(var))!/ (1d0*size(var))

      if(self%me == root) then
         write (*,"(I5,A,A,A,I7,A,ES10.3)") iter, " var: ", var_name, &
            " nkpts ", size(self%all_k_pts,2),&
            " err ", rel_error
      endif

      if(rel_error < self%berry_conv_crit) then
         if(self%me == root) write (*,*) "Converged " // trim(var_name) //   " interation"
         cancel = .True.
      endif
   end function process_hall_surf

   function process_orbmag(self, orbmag, orbmag_old, &
                           orbmag_L, orbmag_IC, iter) result(cancel)
      implicit none
      class(k_space)                 :: self
      real(8), intent(in)            :: orbmag(:), orbmag_old(:),&
                                        orbmag_L(:), orbmag_IC(:)
      integer   , intent(in)         :: iter
      character(len=*), parameter    :: var_name = "orbmag   "
      character(len=300)             :: filename
      logical                        :: cancel
      real(8)                        :: rel_error

      cancel = .False.

      ! save current iteration data
      if(self%me == root) then
         write (filename, "(A,I0.5,A)") "orbmag_iter=", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), orbmag)

         write (filename, "(A,I0.5,A)") "orbmag_L_iter=", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), orbmag_L)

         write (filename, "(A,I0.5,A)") "orbmag_IC_iter=", iter, ".npy"
         call save_npy(trim(self%prefix) // trim(filename), orbmag_IC)

         call save_npy(trim(self%prefix) // "orbmag_E.npy", &
                       self%E_fermi / self%units%energy)
      endif

      ! check for convergence
      rel_error = my_norm2(orbmag - orbmag_old) &
                  / (self%kpts_per_step * self%nProcs * my_norm2(orbmag))!/ (1d0*size(orbmag))

      if(self%me == root) then
         write (*,"(I5,A,A,A,I7,A,ES10.3)") iter, " var: ", var_name, &
            " nkpts ", size(self%all_k_pts,2),&
            " err ", rel_error
      endif

      if(rel_error < self%berry_conv_crit) then
         if(self%me == root) write (*,*) "Converged " // trim(var_name) //   " interation"
         cancel = .True.
      endif
   end function process_orbmag

   subroutine finalize_hall_surf(self, var,varall,var_name)
      implicit none
      class(k_space)              :: self
      real(8), intent(in)         :: var(:)
      real(8), intent(in)         :: varall(:,:)
      character(len=*), intent(in) :: var_name
      character(len=300) :: elem_file

      if(self%me == root) then
         write (*,*) size(self%all_k_pts,2), &
            "saving hall_cond with questionable unit"
         write(*,*) var_name
         write (elem_file, "(A,A)") var_name ,".npy"
         call save_npy(trim(self%prefix) // elem_file, var)
         write (elem_file, "(A,A)") var_name ,"_uc.npy"
         call save_npy(trim(self%prefix) // elem_file, varall)
         write (elem_file, "(A,A)") var_name ,"_E.npy"
         call save_npy(trim(self%prefix) // elem_file, &
                       self%E_fermi / self%units%energy)
      endif
   end subroutine finalize_hall_surf

   subroutine finalize_hall(self, var,varall)
      implicit none
      class(k_space)              :: self
      real(8), intent(in)         :: var(:)
      real(8), intent(in)         :: varall(:,:)

      if(self%me == root) then
         write (*,*) size(self%all_k_pts,2), &
            "saving hall_cond with questionable unit"

         call save_npy(trim(self%prefix) // "hall_cond.npy", var)
         call save_npy(trim(self%prefix) // "hall_cond_uc.npy", varall)
         call save_npy(trim(self%prefix) // "hall_cond_E.npy", &
                       self%E_fermi / self%units%energy)
      endif
   end subroutine finalize_hall

   subroutine finalize_orbmag(self, orbmag, orbmag_L, orbmag_IC)
      implicit none
      class(k_space)          :: self
      real(8), intent(in)     :: orbmag(:), orbmag_L(:), orbmag_IC(:)
      real(8)                 :: area

      area =  self%ham%UC%calc_area()

      if(self%me == root) then
         write (*,*) "saving orbmag"

         call save_npy(trim(self%prefix) // "orbmag.npy", &
                       orbmag * area / self%units%mag_dipol)
         call save_npy(trim(self%prefix) // "orbmag_L.npy", &
                       orbmag_L * area / self%units%mag_dipol)
         call save_npy(trim(self%prefix) // "orbmag_IC.npy", &
                       orbmag_IC * area / self%units%mag_dipol)

         call save_npy(trim(self%prefix) // "orbmag_E.npy", &
                       self%E_fermi / self%units%energy)
      endif
   end subroutine finalize_orbmag

   subroutine integrate_hall(self, kidx_all, omega_z_all, eig_val_all, hall)
      use mpi
      implicit none
      class(k_space)          :: self
      integer   , intent(in)  :: kidx_all(:)
      real(8), intent(in)     :: eig_val_all(:,:), omega_z_all(:,:)
      real(8), allocatable    :: hall(:)
      real(8)                 :: ferm
      integer                 :: loc_idx, n, n_hall, k_idx
      integer                 :: ierr(2), all_err(1)

      all_err = 0
      if(allocated(hall)) deallocate(hall)
      allocate(hall(size(self%E_fermi)), stat=all_err(1))
      call check_ierr(all_err, self%me, "integrate hall allocation")

      !run triangulation
      call run_triang(self%all_k_pts, self%elem_nodes)
      call self%set_weights_ksp()

      !perform integration with all points
      hall =  0d0

      do loc_idx = 1,size(kidx_all)
         k_idx =  kidx_all(loc_idx)
         do n_hall =  1,size(hall)
            n_loop: do n = 1,size(omega_z_all,1)
               ferm  =  self%fermi_distr(eig_val_all(n, loc_idx), n_hall)
               if(ferm /=  0d0) then
                  hall(n_hall) = hall(n_hall) + &
                                 self%weights(k_idx) * omega_z_all(n, loc_idx) * ferm
               else
                  exit n_loop
               endif
            enddo n_loop
         enddo
      enddo

      hall = hall / (2d0*PI)

      ! Allreduce is not suitable for convergence criteria
      ierr = 0
      if(self%me == root) then
         call MPI_Reduce(MPI_IN_PLACE, hall, size(hall), MPI_REAL8, &
                         MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      else
         call MPI_Reduce(hall, hall, size(hall), MPI_REAL8, &
                         MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      endif
      call MPI_Bcast(hall, size(hall), MPI_REAL8, root, &
                     MPI_COMM_WORLD, ierr(2))
      call check_ierr(ierr, self%me, "Hall conductance")
   end subroutine integrate_hall

   subroutine integrate_orbmag(self, Q_kidx_all, Q_L_all, Q_IC_all, orb_mag, orbmag_L, orbmag_IC)
      use mpi
      implicit none
      class(k_space)          :: self
      integer   , intent(in)  :: Q_kidx_all(:)
      real(8), intent(in)     :: Q_L_all(:, :), Q_IC_all(:,:)
      real(8), allocatable    :: orb_mag(:), orbmag_L(:), orbmag_IC(:)
      integer                 :: n_ferm, loc_idx, k_idx
      integer                 :: ierr(4)

      if(allocated(orb_mag))then
         if(size(orb_mag) /= size(self%E_fermi)) deallocate(orb_mag)
      endif
      if(.not. allocated(orb_mag)) allocate(orb_mag(size(self%E_fermi)))
      if(allocated(orbmag_L))then
         if(size(orbmag_L) /= size(self%E_fermi)) deallocate(orbmag_L)
      endif
      if(.not. allocated(orbmag_L)) allocate(orbmag_L(size(self%E_fermi)))
      if(allocated(orbmag_IC))then
         if(size(orbmag_IC) /= size(self%E_fermi)) deallocate(orbmag_IC)
      endif
      if(.not. allocated(orbmag_IC)) allocate(orbmag_IC(size(self%E_fermi)))

      !run triangulation
      call run_triang(self%all_k_pts, self%elem_nodes)
      call self%set_weights_ksp()

      orb_mag   = 0d0
      orbmag_L  = 0d0
      orbmag_IC = 0d0

      do n_ferm = 1, size(orb_mag)
         do loc_idx =  1,size(Q_kidx_all)
            k_idx =  Q_kidx_all(loc_idx)
            orbmag_L(n_ferm) = orbmag_L(n_ferm) &
                               + self%weights(k_idx) * Q_L_all(n_ferm, loc_idx)
            orbmag_IC(n_ferm) = orbmag_IC(n_ferm) &
                                + self%weights(k_idx) * Q_IC_all(n_ferm, loc_idx)
         enddo
      enddo

      orbmag_L  = orbmag_L  / (2d0 * speed_of_light * (2d0*PI)**2)
      orbmag_IC = orbmag_IC / (2d0 * speed_of_light * (2d0*PI)**2)

      ! reduce & bcast local term
      if(self%me == root) then
         call MPI_Reduce(MPI_IN_PLACE, orbmag_L, size(orbmag_L),&
                         MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      else
         call MPI_Reduce(orbmag_L, orbmag_L, size(orbmag_L), MPI_REAL8,&
                         MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      endif
      call MPI_Bcast(orbmag_L, size(orbmag_L), MPI_REAL8, root, &
                     MPI_COMM_WORLD, ierr(2))

      ! reduce & bcast itinerant term
      if(self%me == root) then
         call MPI_Reduce(MPI_IN_PLACE, orbmag_IC, size(orbmag_IC), &
                         MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(3))
      else
         call MPI_Reduce(orbmag_IC, orbmag_IC, size(orbmag_IC), &
                         MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(3))
      endif
      call MPI_Bcast(orbmag_IC, size(orbmag_IC), MPI_REAL8, root, &
                     MPI_COMM_WORLD, ierr(4))

      orb_mag    = orbmag_L + orbmag_IC
      call check_ierr(ierr, me_in=self%me, info="orbmag")
   end subroutine integrate_orbmag

   subroutine set_hall_weights(self, omega_z_all, kidx_all)
      use mpi
      implicit none
      class(k_space)         :: self
      integer   , intent(in) :: kidx_all(:)
      real(8)                :: omega_z_all(:,:)
      integer                :: i, node, k_idx, loc_idx
      integer                :: n_elem
      integer                :: ierr(2), error(2) = [0,0]
      character(len=300)     :: msg = ""

      n_elem = size(self%elem_nodes,1)
      error = 0
      if(allocated(self%refine_weights)) then
         if(size(self%refine_weights) /= n_elem) then
            deallocate(self%refine_weights, stat=error(1), errmsg=msg)
         endif
      endif
      if(.not. allocated(self%refine_weights)) then
         allocate(self%refine_weights(n_elem), stat=error(2))
      endif
      call check_ierr(error, self%me, "hall: set_refine_weights errors")

      self%refine_weights = 0d0

      do i=1,n_elem
         do node=1,3
            k_idx = self%elem_nodes(i,node)
            loc_idx =  find_list_idx(kidx_all, k_idx)

            if(loc_idx > 0) then
               self%refine_weights(i) = self%refine_weights(i) &
                                        + self%weights(k_idx) &
                                        * sum(abs(omega_z_all(:,loc_idx)))
            endif
         enddo
      enddo

      ierr = 0
      if(self%me == root) then
         call MPI_Reduce(MPI_IN_PLACE, self%refine_weights, n_elem,&
                         MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      else
         call MPI_Reduce(self%refine_weights,self%refine_weights, n_elem, &
                         MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      endif
      call MPI_Bcast(self%refine_weights, n_elem, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(2))
      call check_ierr(ierr, self%me, " hall: set_refine_weights")

   end subroutine set_hall_weights

   subroutine set_orbmag_weights(self, Q_all, Q_kidx_all)
      use mpi
      implicit none
      class(k_space)            :: self
      real(8), intent(in)       :: Q_all(:,:)
      integer   , intent(in)    :: Q_kidx_all(:)
      integer                   :: i, node, k_idx, loc_idx
      integer                   :: ierr(2), error(2), n_elem
      character(len=300)        :: msg

      n_elem = size(self%elem_nodes,1)
      error = 0
      if(allocated(self%refine_weights)) then
         if(size(self%refine_weights) /= n_elem) then
            deallocate(self%refine_weights, stat=error(1), errmsg=msg)
         endif
      endif
      if(.not. allocated(self%refine_weights)) then
         allocate(self%refine_weights(n_elem), stat=error(2))
      endif
      call check_ierr(error, self%me, " orb mag: set_refine_weights errors")

      self%refine_weights = 0d0
      do i=1,n_elem
         do node=1,3
            k_idx = self%elem_nodes(i,node)
            loc_idx =  find_list_idx(Q_kidx_all, k_idx)

            if(loc_idx > 0) then
               self%refine_weights(i) = self%refine_weights(i) &
                                        + self%weights(k_idx) * sum(abs(Q_all(:,loc_idx)))
            endif
         enddo
      enddo

      ierr =  0
      if(self%me == root) then
         call MPI_Reduce(MPI_IN_PLACE, self%refine_weights, n_elem,&
                         MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      else
         call MPI_Reduce(self%refine_weights,self%refine_weights, n_elem, &
                         MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
      endif
      call MPI_Bcast(self%refine_weights, n_elem, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(2))
      call check_ierr(ierr, self%me, " orb_mag: set_refine_weights")

   end subroutine set_orbmag_weights

   function setup_A_mtx(Vx_mtx, Vy_mtx) result(A_mtx)
      implicit none
      !class(k_space), intent(in) :: self
      complex(8), intent(in)     :: Vx_mtx(:,:), Vy_mtx(:,:)
      real(8), allocatable       :: A_mtx(:,:)
      integer                    :: m, n
      real(8)                    :: t_start, t_stop

      t_start = MPI_Wtime()

      allocate(A_mtx(size(Vx_mtx,1), size(Vx_mtx,2)))

!$OMP         parallel do private(n) collapse(2) default(shared)
      do m = 1,size(Vx_mtx,1)
         do n = 1,size(Vx_mtx,2)
            A_mtx(n, m) = aimag(Vx_mtx(m, n) *  Vy_mtx(n, m))
         enddo
      enddo

      t_stop = MPI_Wtime()
   end function setup_A_mtx

   subroutine calc_orbmag_z_singleK(self, Q_L, Q_IC, eig_val, Vx_mtx, Vy_mtx)
      implicit none
      class(k_space)           :: self
      real(8)                  :: f_nk, dE, Ef, Q_L(:), Q_IC(:), eig_val(:)
      complex(8)               :: Vx_mtx(:,:), Vy_mtx(:,:)
      real(8), allocatable     :: A_mtx(:,:)
      integer                  :: m, n, n_ferm

      A_mtx = setup_A_mtx(Vx_mtx, Vy_mtx)

      Q_L  = 0
      Q_IC = 0

!$OMP         parallel do private(n, m, Ef, f_nk, dE) collapse(2) default(shared)
      do n_ferm = 1, size(Q_L)
         n_loop: do n = 1, size(A_mtx,1)
            Ef =  self%E_fermi(n_ferm)
            f_nk =  self%fermi_distr(eig_val(n), n_ferm)

            if(f_nk /= 0d0) then
               do m = 1, size(A_mtx,1)
                  dE =  eig_val(m) -  eig_val(n)

                  if(abs(dE) >=  eta) then
                     Q_L(n_ferm)  = Q_L(n_ferm)  &
                                    +       f_nk * A_mtx(m,n)/dE
                     Q_IC(n_ferm) = Q_IC(n_ferm) &
                                    - 2d0 * f_nk * (Ef - eig_val(n)) * A_mtx(m,n)/(dE**2)
                  endif ! m != n
               enddo !m
            endif
         enddo n_loop !n
      enddo ! n_ferm

      deallocate(A_mtx)

   end subroutine calc_orbmag_z_singleK

   subroutine append_eigval(eig_val_all, eig_val_new)
      implicit none
      real(8), allocatable    :: eig_val_all(:,:), eig_val_new(:,:), tmp(:,:)
      integer    :: vec_sz, num_k_all, num_k_new, i,j
      integer    ::ierr(6)

      ierr = 0

      vec_sz =  size(eig_val_new, 1)
      if(.not. allocated(eig_val_all)) allocate(eig_val_all(vec_sz,0), stat=ierr(1))
      num_k_all =  size(eig_val_all, 2)
      num_k_new =  size(eig_val_new, 2)

      allocate(tmp(vec_sz, num_k_all), stat=ierr(2))
      forall(i=1:vec_sz, j=1:num_k_all) tmp(i,j) =  eig_val_all(i,j)
      deallocate(eig_val_all, stat=ierr(3))

      allocate(eig_val_all(vec_sz, num_k_all + num_k_new), stat=ierr(4))
      forall(i=1:vec_sz, j=1:num_k_all) eig_val_all(i,j) = tmp(i,j)
      deallocate(tmp, stat=ierr(5))

      forall(i=1:vec_sz, j=1:num_k_new) eig_val_all(i,j+num_k_all) &
         = eig_val_new(i,j)
      deallocate(eig_val_new, stat=ierr(6))
      call check_ierr(ierr, info="append_eigval")
   end subroutine append_eigval

   subroutine append_kidx(kidx_all, kidx_new)
      implicit none
      integer, allocatable   :: kidx_all(:), kidx_new(:), tmp(:)
      integer                :: i

      if(allocated(kidx_all)) then
         !copy to tmp
         allocate(tmp(size(kidx_all)))
         forall(i=1:size(kidx_all)) tmp(i) = kidx_all(i)

         !reallocate kidx_all
         deallocate(kidx_all)
         allocate(kidx_all(size(tmp) + size(kidx_new)))

         !copy and deallocate
         kidx_all(1:size(tmp)) = tmp
         kidx_all(size(tmp)+1:size(kidx_all)) =  kidx_new
         deallocate(tmp, kidx_new)
      else
         allocate(kidx_all(size(kidx_new)))
         kidx_all =  kidx_new
         deallocate(kidx_new)
      endif
   end subroutine append_kidx

   subroutine append_quantitiy(quantity_all, quantity_new)
      implicit none
      real(8), allocatable      :: quantity_new(:,:), quantity_all(:,:), tmp_z(:,:)
      integer                 :: vec_sz, num_k_old, num_k_new
      integer                 :: ierr(2)

      num_k_old = size(quantity_all, 2)
      num_k_new = size(quantity_new, 2)
      vec_sz    = size(quantity_new, 1)
      ierr      = 0

      allocate(tmp_z(vec_sz, num_k_old),stat=ierr(1))
      tmp_z = quantity_all
      deallocate(quantity_all)
      allocate(quantity_all(vec_sz, num_k_old + num_k_new), stat=ierr(2))
      call check_ierr(ierr, info=" append_quantitiy")

      quantity_all(:,1:num_k_old) = tmp_z
      quantity_all(:,num_k_old+1:num_k_old+num_k_new) = quantity_new
      deallocate(tmp_z)
      deallocate(quantity_new)
   end subroutine append_quantitiy

   subroutine setup_berry_inte_grid(self)
      implicit none
      class(k_space)           :: self

      if(allocated(self%new_k_pts) )then
         deallocate(self%new_k_pts)
      endif

      if(trim(self%ham%UC%uc_type) == "square_2d") then
         call self%setup_inte_grid_para(self%berry_num_k_pts)
      elseif(trim(self%ham%UC%uc_type) == "file_square") then
         call self%setup_inte_grid_para(self%berry_num_k_pts)
      elseif(trim(self%ham%UC%uc_type) == "honey_2d") then
         call self%setup_inte_grid_hex(self%berry_num_k_pts)
      elseif(trim(self%ham%UC%uc_type) == "honey_line") then
         call self%setup_inte_grid_para(self%berry_num_k_pts)
      else
         call error_msg("berry k-grid not known", abort=.True.)
      endif

   end subroutine setup_berry_inte_grid

   subroutine set_fermi(self, cfg)
      use mpi
      implicit none
      class(k_space)         :: self
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

   subroutine find_fermi(self, cfg)
      use mpi
      implicit none
      class(k_space)         :: self
      class(CFG_t)           :: cfg
      real(8)            :: target, delta_old, delta_new
      integer            :: i
      integer            :: ierr

      if(self%me ==  root)then
         call CFG_get(cfg, "dos%fermi_fill", target)
      endif
      call MPI_Bcast(self%E_fermi, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr)

      target = target * self%int_DOS(size(self%int_DOS))
      i =  1

      delta_old = -1d0
      delta_new = -1d0

      do while((delta_old * delta_new) > 0)
         delta_old = target - self%int_DOS(i)
         delta_new = target - self%int_DOS(i+1)

         if(i+1 <= size(self%int_DOS)) then
            i = i + 1
         else
            write (*,*) "Required filling not in DOS range"
            stop
         endif
      enddo

      self%E_fermi = self%E_DOS(i)
      call self%write_fermi()
   end subroutine find_fermi

   subroutine write_fermi(self)
      implicit none
      class(k_space)         :: self
      if(self%me ==  root) then
         call save_npy(trim(self%prefix) //  "fermi.npy", self%E_fermi)
      endif
   end subroutine write_fermi

   function lorentzian(self, x) result(lor)
      implicit none
      class(k_space), intent(in)    :: self
      real(8), intent(in)           :: x
      real(8)                       :: lor

      lor =  self%DOS_gamma &
            / (PI * (x**2 +  self%DOS_gamma**2))
   end function lorentzian

   function fermi_distr(self, E, n_ferm) result(ferm)
      implicit none
      class(k_space), intent(in)    :: self
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

   function find_E_max(self) result(c)
      implicit none
      class(k_space), intent(in)   :: self
      real(8)                      :: l, u, c
      real(8), parameter           :: tol = 1d-6, tar = 1d-12
      integer                      :: Emax, cnt

      l = - 25d0 * self%ham%Vss_sig
      u = 25d0 * self%ham%Vss_sig
      Emax  = size(self%E_fermi)

      if((self%fermi_distr(l, Emax) - tar) &
         * (self%fermi_distr(u, Emax) - tar) > 0d0) then
         write (*,*) "Emax bisection failed. So crossing in window"
         stop
      endif

      c =  0.5 * (u + l)
      cnt =  0
      do while(abs((self%fermi_distr(c, Emax) - tar)/tar) > tol)
         if(sign(1d0,self%fermi_distr(c, Emax) - tar) ==&
            sign(1d0,self%fermi_distr(l,Emax)- tar)) then
            l = c
         else
            u = c
         endif
         c =  0.5 * (u + l)
         cnt =  cnt + 1
      enddo
   end function find_E_max

   pure function area_of_elem(self, kpts, idx) result(area)
      implicit none
      class(k_space), intent(in)     :: self
      integer   , intent(in)     :: idx
      real(8), intent(in)        :: kpts(:,:)
      real(8)                    :: area
      real(8)    :: vec1(3), vec2(3)

      vec1 = kpts(:,self%elem_nodes(idx,1)) - kpts(:,self%elem_nodes(idx,2))
      vec2 = kpts(:,self%elem_nodes(idx,1)) - kpts(:,self%elem_nodes(idx,3))
      area =  0.5d0 * my_norm2(cross_prod(vec1, vec2))
   end function area_of_elem

   function new_pt(self, idx, kpts) result(pt)
      implicit none
      class(k_space), intent(in)   ::self
      integer   , intent(in)   :: idx
      real(8), intent(in)      :: kpts(:,:)
      real(8)                      :: pt(2), start(2), vec(3), len, len_max

      vec = kpts(:,self%elem_nodes(idx,1)) - kpts(:,self%elem_nodes(idx,2))
      len =  my_norm2(vec)
      len_max =  len
      start =  kpts(1:2,self%elem_nodes(idx,2))
      pt = start + 0.5 * vec(1:2)

      vec = kpts(:,self%elem_nodes(idx,1)) - kpts(:,self%elem_nodes(idx,3))
      len =  my_norm2(vec)
      if(len >  len_max) then
         len_max =  len
         start =  kpts(1:2,self%elem_nodes(idx,3))
         pt = start + 0.5 * vec(1:2)
      endif

      vec = kpts(:,self%elem_nodes(idx,2)) - kpts(:,self%elem_nodes(idx,3))
      len =  my_norm2(vec)
      if(len >  len_max) then
         len_max =  len
         start =  kpts(1:2,self%elem_nodes(idx,3))
         pt = start + 0.5 * vec(1:2)
      endif

   end function new_pt

   function centeroid_of_triang(self, idx, k_pts) result(centeroid)
      implicit none
      class(k_space), intent(in)  :: self
      integer   , intent(in)      :: idx
      integer                     :: i
      real(8), intent(in)         :: k_pts(:,:)
      real(8)                     :: centeroid(2)

      centeroid =  0d0
      do i = 1,3
         centeroid(1) =  centeroid(1) + k_pts(1,self%elem_nodes(idx,i))
         centeroid(2) =  centeroid(2) + k_pts(2,self%elem_nodes(idx,i))
      enddo
      centeroid = centeroid * 0.33333333333d0
   end function centeroid_of_triang

   subroutine pad_k_points_init(self)
      implicit none
      class(k_space)                :: self
      integer     :: rest, i,j, cnt, n_kpts, n_elem, per_proc
      integer   , allocatable :: sort(:)
      real(8), allocatable    :: areas(:), new_ks(:,:), tmp(:,:)
      real(8)                 :: cand(2)

      n_kpts =  size(self%new_k_pts,2)
      if(self%me == root) write (*,*) "nkpts =  ", n_kpts

      per_proc = CEILING((1d0*n_kpts)/(1d0*self%nProcs))
      rest =  self%nProcs * per_proc - n_kpts

      if(self%me == root) write (*,*) "per_proc", per_proc
      if(self%me == root) write (*,*) "rest", rest

      if(rest /= 0) then
         !find biggest triangles
         n_elem =  size(self%elem_nodes,1)
         if(rest >  n_elem) then
            call error_msg("Not enough elements for inital padding", abort=.True.)
         endif

         allocate(areas(n_elem))
         allocate(new_ks(3,rest))

         forall(i = 1:n_elem) areas(i) = self%area_of_elem(self%new_k_pts, i)
         call qargsort(areas, sort)

         !calculate new elements
         new_ks =  0d0
         i = n_elem
         do cnt = 1,rest
            cand = self%new_pt(sort(i), self%new_k_pts)
            do while(self%in_points(cand, new_ks(1:2,1:cnt-1)))
               i = i - 1
               if(i < 1) then
                  call error_msg("Not enough elemes for refinement", abort=.True.)
               endif
               cand = self%new_pt(sort(i), self%new_k_pts)
            enddo
            new_ks(1:2,cnt) = cand
         enddo

         !extent list
         allocate(tmp(3,n_kpts))
         forall(i=1:3, j=1:n_kpts) tmp(i,j) =  self%new_k_pts(i,j)
         deallocate(self%new_k_pts)
         allocate(self%new_k_pts(3, n_kpts +  rest))
         forall(i=1:3, j=1:n_kpts) self%new_k_pts(i,j) = tmp(i,j)
         deallocate(tmp)

         !add new ones
         forall(i=1:3, j=1:rest) self%new_k_pts(i,n_kpts+j) = new_ks(i,j)
         call run_triang(self%new_k_pts, self%elem_nodes)
         if(self%me == root)write (*,*) "elem_post pad", shape(self%elem_nodes)
      endif
   end subroutine pad_k_points_init

   function in_points(self, pt, list) result(inside)
      implicit none
      class(k_space)       :: self
      real(8), intent(in)  :: pt(2), list(:,:)
      logical              :: inside
      real(8)              :: l
      integer              :: i

      inside =  .False.
      l = my_norm2(self%ham%UC%rez_lattice(:,1))

      do i = 1,size(list,2)
         if(my_norm2(pt - list(:,i)) < l * pos_eps) then
            inside = .True.
         endif
      enddo
   end function in_points

   function on_hex_border(self, idx) result(on_border)
      implicit none
      class(k_space), intent(in)   :: self
      integer                      :: idx
      real(8)                      :: pt(2)
      logical                      :: on_border
      real(8)                      :: l

      pt =  self%new_k_pts(1:2,idx)
      l = my_norm2(self%ham%UC%rez_lattice(:,1))

      if(abs(abs(pt(2)) - 0.5d0 * l) < pos_eps * l) then
         on_border = .True.
      else
         on_border = abs(abs(pt(1)) - self%hex_border_x(pt(2))) < pos_eps * l
      endif
   end function on_hex_border

   subroutine add_kpts_iter(self, n_new, new_ks)
      implicit none
      class(k_space)              :: self
      integer, intent(in)     :: n_new
      real(8), allocatable    :: new_ks(:,:), areas(:)
      integer                 :: n_elem, i, cnt, area_cnt, weight_cnt, num_kpts
      integer   , allocatable :: sort_weight(:), sort_area(:)

      if(allocated(new_ks)) then
         if(size(new_ks,1) /= 3 .or. size(new_ks,2) /= n_new) then
            deallocate(new_ks)
         endif
      endif
      if(.not. allocated(new_ks)) allocate(new_ks(3, n_new))
      n_elem = size(self%elem_nodes,1)

      allocate(areas(n_elem))
      forall(i = 1:n_elem) areas(i) = self%area_of_elem(self%all_k_pts, i)

      call qargsort(areas, sort_area)
      call qargsort(self%refine_weights, sort_weight)

      new_ks =  0d0
      i = n_elem
      area_cnt   = n_elem
      weight_cnt = n_elem
      num_kpts   = size(self%all_k_pts, 2)

      do cnt = 1,n_new
         if(trim(self%ada_mode) == "area") then
            new_ks(1:2, cnt) = self%centeroid_of_triang(sort_area(i), self%all_k_pts)
         elseif(trim(self%ada_mode) == "weight") then
            new_ks(1:2, cnt) = self%centeroid_of_triang(sort_weight(i), self%all_k_pts)
         elseif(trim(self%ada_mode) == "mixed") then
            if(mod(cnt + num_kpts,2) == 0) then
               new_ks(1:2, cnt) = self%centeroid_of_triang(sort_weight(weight_cnt), self%all_k_pts)
               weight_cnt = weight_cnt - 1
            else
               new_ks(1:2, cnt) = self%centeroid_of_triang(sort_area(area_cnt), self%all_k_pts)
               area_cnt = area_cnt - 1
            endif
         else
            call error_msg("adaptation not known", abort=.True.)
         endif
         i = i - 1
         if(i == 0) then
            call error_msg("Not enough elements", abort=.True.)
         endif
      enddo
      deallocate(areas)
   end subroutine add_kpts_iter

   subroutine append_kpts(self)
      implicit none
      class(k_space)         :: self
      real(8), allocatable   :: tmp(:,:)
      integer                :: old_sz, new_sz, i,j
      integer                :: ierr(2)

      if(.not. allocated(self%all_k_pts)) allocate(self%all_k_pts(3,0))
      old_sz = size(self%all_k_pts, 2)
      new_sz = size(self%new_k_pts, 2)
      ierr   = 0

      allocate(tmp(3, old_sz), stat=ierr(1))

      forall(i=1:3, j=1:old_sz) tmp(i,j) = self%all_k_pts(i,j)

      deallocate(self%all_k_pts)
      allocate(self%all_k_pts(3, old_sz + new_sz), stat=ierr(2))
      call check_ierr(ierr, self%me, " append_kpts")

      forall(i=1:3, j=1:old_sz) self%all_k_pts(i,j)        = tmp(i,j)
      forall(i=1:3, j=1:new_sz) self%all_k_pts(i,old_sz+j) = self%new_k_pts(i,j)

      deallocate(tmp)
      deallocate(self%new_k_pts)
   end subroutine append_kpts

   subroutine calc_ACA(self)
      use mpi
      implicit none
      class(k_space)              :: self
      real(8), allocatable    :: m(:), S(:), l_space(:), eig_val(:), RWORK(:)
      real(8)                 :: area, t_start, t_stop
      complex(8), allocatable :: H(:,:), WORK(:)
      integer                 :: N_k, N, &
                                 first, last, k_idx, info, ierr
      integer                 :: lwork, lrwork, liwork
      integer, allocatable    :: IWORK(:)

      if(self%ham%num_orb /= 3) then
         call error_msg("ACA only for p-orbitals", abort=.True.)
      endif

      N = 2 * self%ham%num_up
      allocate(H(N,N))
      allocate(m(size(self%E_fermi)))
      allocate(S(size(self%E_fermi)))
      allocate(eig_val(N))
      call calc_zheevd_size('V', H, eig_val, lwork, lrwork, liwork)
      allocate(WORK(lwork))
      allocate(RWORK(lrwork))
      allocate(IWORK(liwork))

      m = 0d0
      S = 0d0

      call linspace(0d0, 1d0, self%ACA_num_k_pts, l_space)
      t_start = MPI_Wtime()

      if( trim(self%ham%UC%uc_type) == "square_2d" &
         .or.trim(self%ham%UC%uc_type) == "file_square" ) then
         call self%setup_inte_grid_para(self%ACA_num_k_pts, padding=.False.)
         N_k = size(self%new_k_pts, 2)
         call my_section(self%me, self%nProcs, N_k, first, last)

         do k_idx = first, last
            if(self%me == root) write (*,*) "Go kpt", k_idx, "of", last, date_time()
            call self%ham%setup_H(self%new_k_pts(:,k_idx), H)
            if(self%me == root) write (*,*) "Done setup   ", date_time()

            call zheevd('V', 'U', N, H, N, eig_val, WORK, lwork, &
                        RWORK, lrwork, IWORK, liwork, info)
            if(self%me == root) write (*,*) "Done zheevd    ", date_time()
            if(info /= 0) call error_msg("zheevd fail", abort=.True.)

            ! write (*,*) "k = ", self%new_k_pts(:,k_idx)
            m = m + self%calc_ACA_singleK(H, eig_val)
            if(self%me == root) write (*,*) "Done ACA     ", date_time()

            S = S + self%calc_S_singleK(H, eig_val)
            if(self%me == root) write (*,*) "Done magnetization   ", date_time()
            ! write(*,*) "current m = ", m/k_idx
         enddo

         area =  self%ham%UC%calc_area()
         m =  m / (N_k * area)
         S = S / N_k

         if(self%me == root) then
            call MPI_Reduce(MPI_IN_PLACE, m, size(m), MPI_REAL8, &
                            MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(MPI_IN_PLACE, S, size(S), MPI_REAL8, &
                            MPI_SUM, root, MPI_COMM_WORLD, ierr)
         else
            call MPI_Reduce(m, m, size(m), MPI_REAL8, &
                            MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(S, S, size(S), MPI_REAL8, &
                            MPI_SUM, root, MPI_COMM_WORLD, ierr)
         endif

         if(self%me ==  root) then
            call error_msg("Wrote ACA orbmag with questionable unit", &
                           p_color=c_green)
            call save_npy(trim(self%prefix) // "orbmag_ACA.npy", m)
            call save_npy(trim(self%prefix) // "orbmag_E.npy", &
                          self%E_fermi / self%units%energy)

            call save_npy(trim(self%prefix) // "spinmag.npy", S)
            call save_npy(trim(self%prefix) // "spinmag_E.npy", &
                          self%E_fermi / self%units%energy)
         endif
      else
         call error_msg("ACA implemented for square only.", abort=.True.)
      endif
      t_stop = MPI_Wtime()
      if(self%me == root) write (*,*) "ACA time = ", t_stop - t_start
      deallocate(H, WORK, RWORK, IWORK, m)
   end subroutine calc_ACA

   function calc_ACA_singleK(self, eig_vec, eig_val) result(m)
      implicit none
      class(k_space), intent(in)     :: self
      complex(8), intent(in)     :: eig_vec(:,:)
      real(8), intent(in)        :: eig_val(:)
      real(8)                    :: m(size(self%E_fermi))
      integer                    :: n_ferm

      m        = 0d0
!$OMP         parallel do default(shared)
      do n_ferm = 1,size(self%E_fermi)
         m(n_ferm) = m(n_ferm) + sum(self%calc_local_l(eig_vec, eig_val, n_ferm))
      enddo
!$OMP         end parallel do

      ! we drop the - 0.5 so we are directly in mu_b
      ! m = - 0.5d0 * m
      m = - m
   end function calc_ACA_singleK

   function calc_S_singleK(self, eig_vec, eig_val) result(S)
      implicit none
      class(k_space), intent(in)     :: self
      complex(8), intent(in)     :: eig_vec(:,:)
      real(8), intent(in)        :: eig_val(:)
      real(8)                    :: S(size(self%E_fermi)), f
      integer                    :: n_up, n_ferm, i

      n_up = self%ham%num_up
      S    = 0d0

!$OMP         parallel do default(shared) private(i, f)
      do n_ferm = 1,size(self%E_fermi)
         do i = 1,size(eig_val)
            f = self%fermi_distr(eig_val(i), n_ferm)
            S(n_ferm) = S(n_ferm) &
                        + f * (  sum(abs(eig_vec(     1:  n_up, i))**2) &
                               - sum(abs(eig_vec(n_up+1:2*n_up, i))**2) )
         enddo
      enddo
!$OMP         end parallel do
   end function calc_S_singleK

   function calc_local_l(self, eig_vec, eig_val, n_ferm) result(loc_l)
      implicit none
      class(k_space), intent(in)     :: self
      complex(8), intent(in)         :: eig_vec(:,:)
      real(8), intent(in)            :: eig_val(:)
      integer, intent(in)            :: n_ferm
      real(8)                        :: loc_l(self%ham%UC%num_atoms), f
      integer                        :: i, s, n_stat, v_u, v_d

      n_stat = 2 * self%ham%num_up
      loc_l  = 0d0
      do i = 1,self%ham%UC%num_atoms
         ! vector index up
         v_u = 3 * i - 2
         ! vector index down
         v_d = v_u + self%ham%num_up

         do s = 1,n_stat
            ! get fermi_factor
            f = self%fermi_distr(eig_val(s), n_ferm)
            loc_l(i) = loc_l(i) + f * calc_l(eig_vec(v_u:v_u+2, s))
            loc_l(i) = loc_l(i) + f * calc_l(eig_vec(v_d:v_d+2, s))
         enddo
      enddo
   end function calc_local_l

   function calc_l(p) result(l)
      implicit none
      complex(8), intent(in)  :: p(3)
      complex(8)              :: l_prime
      real(8)                 :: l
      complex(8), parameter   :: m_l(3,3) &
                                 = transpose(reshape( &
                                             [c_0, c_i, c_0, &
                                              -c_i, c_0, c_0, &
                                              c_0, c_0, c_0], &
                                             shape(m_l)))
      complex(8)              :: rhs(3)

      rhs     = matmul(m_l, p)
      ! remember complex dot_product = sum(conjg(a) * b)
      l_prime = dot_product(p, rhs)

      if(abs(aimag(l_prime)) > 1d-11) then
         call error_msg("l is imaginary", abort=.True.)
      endif
      l = real(l_prime)
   end function calc_l

   function test_func(kpt) result(ret)
      implicit none
      real(8)                :: kpt(2), ret, d

      d =  my_norm2(kpt)
      ret = exp(-10d0*d)
   end function test_func

   function random_pt_hex(self) result(ret)
      implicit none
      class(k_space)            :: self
      real(8)                   :: ret(2), l

      l = my_norm2(self%ham%UC%rez_lattice(:,1))

      call random_number(ret)
      ret =  2d0 * l * (ret - 0.5d0)

      do while( (abs(ret(2)) > 0.5 * l) .or. &
               (abs(ret(1)) > self%hex_border_x(ret(2))))
         call random_number(ret)
         ret =  2d0 * l * (ret - 0.5d0)
      enddo
   end function random_pt_hex

   subroutine save_grid(self, iter)
      implicit none
      class(k_space), intent(in)  :: self
      integer, intent(in)     :: iter
      character(len=300)      :: k_file, elem_file

      if(self%me == root) then
        write (elem_file, "(A,I0.5,A)") trim(self%prefix) // "elem_iter=", iter, ".npy"
        call save_npy(trim(elem_file), self%elem_nodes) 
        !if (self%ham%UC%num_atoms==2) then 
        if (iter == self%berry_iter) then
            write (k_file,    "(A,I0.5,A)") trim(self%prefix) // "kpts_iter=", iter, ".npy"
            call save_npy(trim(k_file),    self%all_k_pts)
        endif
      endif

   end subroutine save_grid

   subroutine plot_omega_square(self)
      implicit none
      class(k_space)               :: self
      real(8), allocatable     :: eig_val(:,:), omega_z(:,:), omega_surf(:,:), omega_sea(:,:), Q_L(:,:), Q_IC(:,:)
      logical                  :: tmp_ch, tmp_co

      if(self%nProcs /= 1) call error_msg("Plot only for 1 process", abort=.True.)

      tmp_ch           = self%calc_hall
      tmp_co           = self%calc_orbmag
      self%calc_hall   = .True.
      self%calc_orbmag = .False.

      call self%setup_inte_grid_para(self%num_plot_pts, padding=.False.)
      call self%calc_new_berry_points(eig_val, omega_z, omega_surf, omega_sea, Q_L, Q_IC,.False.)

      call save_npy(trim(self%prefix) // "k_grid.npy", self%new_k_pts)
      call save_npy(trim(self%prefix) // "omega_z.npy", omega_z)

      deallocate(self%new_k_pts)
      self%calc_hall   = tmp_ch
      self%calc_orbmag = tmp_co
   end subroutine plot_omega_square

   subroutine basis_trafo(in_EV, trafo_mtx, out_EV)
      implicit none
      complex(8), intent(in)   :: in_EV(:)
      complex(8), intent(in)   :: trafo_mtx(3,3)
      complex(8)               :: out_EV(:)
      integer                  :: i

      out_EV = 0d0
      do i = 1,size(in_EV),3
         call zgemv('N',      3,        3,    c_1, trafo_mtx, 3, &
                    in_EV(i: i+2),     1,    &
                    c_0,     out_EV(i: i+2), 1)
      enddo
   end subroutine basis_trafo

   subroutine p_real_to_cmplx(real_EV, cmplx_EV)
      implicit none
      complex(8), intent(in)   :: real_EV(:)
      complex(8)               :: cmplx_EV(:)

      call basis_trafo(real_EV, BT_real_to_cmplx, cmplx_EV)
   end subroutine p_real_to_cmplx

   subroutine p_cmplx_to_real(cmplx_EV, real_EV)
      implicit none
      complex(8), intent(in)   :: cmplx_EV(:)
      complex(8)               :: real_EV(:)

      call basis_trafo(cmplx_EV, BT_cmplx_to_real, real_EV)
   end subroutine p_cmplx_to_real
end module
