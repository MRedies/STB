module Class_k_space
    use m_config
    use m_npy
    use mpi
    use Class_hamiltionian
    use Class_helper
    implicit none

    type k_space
        real(8), allocatable :: new_k_pts(:,:), all_k_pts(:,:)
        real(8), allocatable :: int_DOS(:) !> integrated Density of states  
        real(8), allocatable :: E_DOS(:)
        real(8), allocatable :: E_fermi(:) !> Fermi lvl
        real(8) :: DOS_gamma !> broadening \f$ \Gamma \f$ used in
        !> DOS calculations
        real(8) :: DOS_lower !> lower energy bound for DOS calc
        real(8) :: DOS_upper !> upper energy bound for DOS calc
        real(8) :: temp !> temperature used in fermi-dirac
        integer(4) :: DOS_num_k_pts !> number of kpts per dim 
        !> used in DOS calculations
        integer(4) :: berry_num_k_pts !> number of ks per dim in berry calc
        integer(4) :: num_DOS_pts!> number of points on E grid
        integer(4) :: num_k_pts !> number of k_pts per segment
        integer(4) :: nProcs !> number of MPI Processes
        integer(4) :: me !> MPI rank
        integer(4) :: laplace_iter !> number of laplace iterations
        integer(4) :: berry_iter !> number of grid refinements
        integer(4) :: kpts_per_step !> new kpts per step and Proc
        real(8)    :: berry_k_shift(3) !> shift of brillouine-zone
        real(8)    :: berry_conv_crit !> convergance criterion for berry integration
        real(8), allocatable :: weights(:) !> weights for integration
        integer(4), allocatable :: elem_nodes(:,:) !> elements in triangulation
        real(8), allocatable :: hall_weights(:)
        real(8), allocatable :: k1_param(:) !> 1st k_space param
        real(8), allocatable :: k2_param(:) !> 2nd k_space param
        character(len=300)    :: filling, prefix
        logical      :: perform_dos_integration !> param to toggle dos integr.
        logical      :: perform_pad !> should the k-grid be padded, to match cores
        type(hamil)  :: ham
        type(units)  :: units
    contains

        procedure :: vol_k_space_parallelo  => vol_k_space_parallelo
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
        procedure :: calc_hall_conductance  => calc_hall_conductance
        procedure :: setup_inte_grid_square => setup_inte_grid_square
        procedure :: setup_inte_grid_hex    => setup_inte_grid_hex
        procedure :: setup_berry_inte_grid  => setup_berry_inte_grid
        procedure :: set_weights_ksp        => set_weights_ksp
        procedure :: test_integration       => test_integration
        procedure :: Bcast_k_space          => Bcast_k_space
        procedure :: hex_border_x           => hex_border_x
        procedure :: vol_k_hex              => vol_k_hex
        procedure :: plot_omega             => plot_omega
        procedure :: free_ksp               => free_ksp
        procedure :: find_E_max             => find_E_max
        procedure :: area_of_elem           => area_of_elem
        procedure :: centeroid_of_elem      => centeroid_of_elem
        procedure :: pad_k_points_init      => pad_k_points_init
        procedure :: new_pt                 => new_pt
        procedure :: on_hex_border          => on_hex_border
        procedure :: hex_laplace_smoother   => hex_laplace_smoother
        procedure :: in_points              => in_points
        procedure :: append_kpts            => append_kpts
        procedure :: add_kpts_iter          => add_kpts_iter
        procedure :: random_pt_hex          => random_pt_hex
        procedure :: set_hall_weights       => set_hall_weights
        procedure :: integrate_hall         => integrate_hall
    end type k_space

    type :: r8arr
        real(8), allocatable :: arr(:)
    end type r8arr

    interface
        subroutine run_triang(k_pts, ret_elem)
            real(8), intent(in)              :: k_pts(:,:)
            integer(4), allocatable          :: ret_elem(:,:)
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
        Implicit None
    class(k_space)                :: self 
        integer(4)                    :: first, last, N, send_count, ierr
        integer(4), allocatable       :: num_elems(:), offsets(:)
        real(8), allocatable          :: eig_val(:,:), sec_eig_val(:,:), k_pts_sec(:,:)


        if(trim(self%filling) ==  "path_rel") then
            call self%setup_k_path_rel()
        else if(trim(self%filling) == "path_abs") then
            call self%setup_k_path_abs()
        else if(trim(self%filling) == "grid") then
            call self%setup_k_grid()
        else
            write (*,*) "Filling not known"
            stop
        endif

        call my_section(self%me, self%nProcs, size(self%new_k_pts, 2), first, last)
        allocate(k_pts_sec(3, last - first + 1))
        k_pts_sec = self%new_k_pts(:,first:last)

        call self%ham%calc_eigenvalues(k_pts_sec, sec_eig_val)


        N = 2 *  self%ham%UC%num_atoms 
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
        implicit none
    class(k_space)          :: self
        real(8), intent(in)     :: E(:)
        real(8), intent(out)    :: PDOS(:,:)
        real(8), allocatable    :: RWORK(:), eig_val(:), loc_PDOS(:,:)
        real(8)                 :: lor
        complex(8), allocatable :: H(:,:), WORK(:)
        integer(4), allocatable :: IWORK(:)
        integer(4)  :: k_idx, E_idx, j, m, N, info, first, last, ierr, num_atoms
        integer(4)  :: lwork, liwork, lrwork

        N =  2 * self%ham%UC%num_atoms
        allocate(H(N,N))
        allocate(eig_val(N))

        call calc_zheevd_size('V', H, eig_val, lwork, lrwork, liwork)
        allocate(WORK(lwork))
        allocate(RWORK(lrwork))
        allocate(IWORK(liwork))

        num_atoms =  self%ham%UC%num_atoms
        allocate(loc_PDOS(2*num_atoms, self%num_DOS_pts))

        loc_PDOS =  0d0
        PDOS = 0d0

        call my_section(self%me, self%nProcs, size(self%new_k_pts,2), first, last)
        do k_idx=first, last
            call self%ham%setup_H(self%new_k_pts(:,k_idx), H)
            call zheevd('V', 'U', N, H, N, eig_val, WORK, lwork, &
                RWORK, lrwork, IWORK, liwork, info)
            if( info /= 0) then
                write (*,*) self%me, ": ZHEEVD (with vectors) failed: ", info
                stop
            endif

            !$omp parallel do private(j) default(shared)
            do m = 1,N
                do j = 1,N
                    ! calc absolute of eigen_vectors
                    H(j,m) =  H(j,m) *  conjg(H(j,m)) 
                enddo
            enddo


            !$omp parallel do private(m,j,lor) default(shared)
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
        write (*,*) "ksz", size(self%new_k_pts, 2)
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
        integer(4)           :: i, num_atoms, info


        if(trim(self%ham%UC%uc_type) == "square_2d") then
            call self%setup_inte_grid_square(self%DOS_num_k_pts)
        elseif(trim(self%ham%UC%uc_type) == "honey_2d") then
            call self%setup_inte_grid_hex(self%DOS_num_k_pts)
        else
            if(self%me ==  root) write (*,*) "DOS k-grid not known"
            call MPI_Abort(MPI_COMM_WORLD, 0, info)
        endif

        num_atoms =  self%ham%UC%num_atoms
        allocate(PDOS(2*num_atoms, self%num_DOS_pts))

        call linspace(self%DOS_lower, self%DOS_upper, self%num_DOS_pts, self%E_DOS)

        if(self%me ==  root) then
            allocate(DOS(self%num_DOS_pts))
            allocate(up(self%num_DOS_pts))
            allocate(down(self%num_DOS_pts))
        endif

        call self%calc_pdos(self%E_DOS, PDOS)

        if(self%me == root) then
            DOS  = sum(PDOS,1)
            up   = sum(PDOS(1:num_atoms, :),1)
            down = sum(PDOS(num_atoms+1:2*num_atoms, :),1)

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
                write (*,*) "Can't perform integration. Only one point"
                stop
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
        implicit none
        type(k_space)         :: self
        type(CFG_t)           :: cfg
        real(8)               :: tmp
        integer(4)            :: sz, ierr

        call MPI_Comm_size(MPI_COMM_WORLD, self%nProcs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)

        self%units = init_units(cfg, self%me)
        self%ham   = init_hamil(cfg)

        if(self%me ==  0) then 
            call CFG_get(cfg, "output%band_prefix", self%prefix)
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

            call CFG_get(cfg, "berry%k_pts_per_dim", self%berry_num_k_pts)
            call CFG_get(cfg, "berry%temperature", tmp)
            self%temp = tmp * self%units%temperature

            call CFG_get(cfg, "berry%laplace_iter", self%laplace_iter)
            call CFG_get(cfg, "berry%refinement_iter", self%berry_iter)
            call CFG_get(cfg, "berry%kpts_per_step", self%kpts_per_step)
            call CFG_get(cfg, "berry%k_shift", self%berry_k_shift)
            call CFG_get(cfg, "berry%conv_criterion", self%berry_conv_crit)
            call CFG_get(cfg, "berry%perform_pad", self%perform_pad)
        endif
        call self%Bcast_k_space()
    end function init_k_space

    subroutine Bcast_k_space(self)
    class(k_space)         :: self
        integer(4), parameter  :: num_cast =  18
        integer(4)             :: ierr(num_cast), sz(2)

        if(self%me ==  root) then
            sz(1) = size(self%k1_param)
            sz(2) = size(self%k2_param)
        endif

        call MPI_Bcast(self%prefix,  300,            MPI_CHARACTER, &
            root,          MPI_COMM_WORLD, ierr(1))
        call MPI_Bcast(self%filling, 300,            MPI_CHARACTER, &
            root,          MPI_COMM_WORLD, ierr(2))

        call MPI_Bcast(self%DOS_gamma,   1,              MPI_REAL8,    &
            root,            MPI_COMM_WORLD, ierr(3))
        call MPI_Bcast(self%num_DOS_pts, 1,              MPI_INTEGER4, &
            root,            MPI_COMM_WORLD, ierr(4))

        call MPI_Bcast(sz,   2,              MPI_INTEGER4, &
            root, MPI_COMM_WORLD, ierr(5))

        if(self%me /= root) then
            allocate(self%k1_param(sz(1)))
            allocate(self%k2_param(sz(2)))
        endif
        call MPI_Bcast(self%k1_param,      sz(1),          MPI_REAL8,    &
            root,               MPI_COMM_WORLD, ierr(5))
        call MPI_Bcast(self%k2_param,      sz(2),          MPI_REAL8,    &
            root,               MPI_COMM_WORLD, ierr(6))
        call MPI_Bcast(self%num_k_pts,     1,              MPI_INTEGER4, &
            root,               MPI_COMM_WORLD, ierr(7))
        call MPI_Bcast(self%DOS_num_k_pts, 1,              MPI_INTEGER4, &
            root,               MPI_COMM_WORLD, ierr(8))
        call MPI_Bcast(self%DOS_lower,     1,              MPI_REAL8, &
            root,               MPI_COMM_WORLD, ierr(9))
        call MPI_Bcast(self%DOS_upper,     1,              MPI_REAL8, &
            root,               MPI_COMM_WORLD, ierr(10))

        ! Berry parameter
        call MPI_Bcast(self%berry_num_k_pts, 1,              MPI_INTEGER4, &
            root,                            MPI_COMM_WORLD, ierr(11))
        call MPI_Bcast(self%temp,            1,              MPI_REAL8,    &
            root,                            MPI_COMM_WORLD, ierr(12))
        call MPI_Bcast(self%laplace_iter,    1,              MPI_INTEGER4, &
            root,                            MPI_COMM_WORLD, ierr(13))
        call MPI_Bcast(self%berry_iter,      1,              MPI_INTEGER4, &
            root,                            MPI_COMM_WORLD, ierr(14))
        call MPI_Bcast(self%kpts_per_step,   1,              MPI_INTEGER4, &
            root,                            MPI_COMM_WORLD, ierr(15))
        call MPI_Bcast(self%berry_k_shift,   3,              MPI_REAL8,    &
            root,                            MPI_COMM_WORLD, ierr(16))
        call MPI_Bcast(self%berry_conv_crit, 1,              MPI_REAL8,    &
            root,                            MPI_COMM_WORLD, ierr(17))
        call MPI_Bcast(self%perform_pad,     1,              MPI_LOGICAL,  &
            root,                            MPI_COMM_WORLD, ierr(18))

        call check_ierr(ierr, self%me, "Ksp Bcast")
    end subroutine Bcast_k_space

    subroutine setup_k_grid(self)
        implicit none
    class(k_space)       :: self
        real(8), allocatable :: kx_points(:), ky_points(:)
        real(8), allocatable :: kx_grid(:,:), ky_grid(:,:)
        integer(4)           :: sz_x, sz_y, i, j


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

    subroutine setup_inte_grid_square(self, n_k)
        implicit none
    class(k_space)        :: self
        integer(4), intent(in):: n_k
        real(8), allocatable  :: ls(:), ls_help(:)
        real(8)               :: k1(3), k2(3)
        integer(4)            :: i, j, cnt

        if(allocated(self%new_k_pts)) deallocate(self%new_k_pts)
        allocate(self%new_k_pts(3,n_k**2))

        k1 =  0d0
        k2 =  0d0 
        k1(1:2) =  self%ham%UC%rez_lattice(:,1)
        k2(1:2) =  self%ham%UC%rez_lattice(:,2)

        call linspace(0d0, 1d0, n_k+1, ls_help)
        allocate(ls(n_k))

        ls         = ls_help(1:n_k)
        self%new_k_pts = 0d0
        cnt        = 1

        do i = 1,n_k
            do j =  1,n_k
                self%new_k_pts(:,cnt) = ls(i) *  k1 +  ls(j) *  k2
                cnt =  cnt + 1
            enddo
        enddo
    end subroutine setup_inte_grid_square

    subroutine setup_inte_grid_hex(self, n_k)
        implicit none
    class(k_space)         :: self
        integer(4), intent(in) :: n_k
        real(8), allocatable   :: x(:), y(:)
        real(8)                :: den, l, a
        real(8), parameter     :: deg_30 = 30.0/180.0*PI, deg_60 = 60.0/180.0*PI
        integer(4)             :: cnt_k, start, halt, my_n, i

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
                self%new_k_pts(:,i) + l * self%berry_k_shift
    end subroutine setup_inte_grid_hex

    subroutine test_integration(self, iter)
        implicit none
    class(k_space)     :: self
        integer(4)         :: i, iter
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
        integer(4)       :: i, j, k_idx

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
        real(8), parameter  :: deg_30 = 30.0/180.0*PI, deg_60 = 60.0/180.0*PI

        l = my_norm2(self%ham%UC%rez_lattice(:,1))
        a = l / (2.0 * cos(deg_30))
        m =  - 2.0 * sin(deg_60)
        b =  - m * a

        x = (abs(y)-b)/m
    end function hex_border_x

    subroutine setup_k_path_abs(self)
        implicit none
    class(k_space)        :: self
        integer(4)            :: n_pts, n_sec, start, halt, i
        real(8), allocatable  :: tmp(:)



        self%k1_param =  self%k1_param * self%units%inv_length
        self%k2_param =  self%k2_param * self%units%inv_length

        n_pts =  self%num_k_pts
        n_sec =  size(self%k1_param)-1


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
        integer(4) ::  n_pts, n_sec,i,j, start, halt, cnt

        n_sec =  size(self%k1_param) - 1
        n_pts =  self%num_k_pts 

        allocate(self%new_k_pts(3,n_sec * (n_pts - 1) + 1))
        allocate(c1_sec(n_sec))
        allocate(c2_sec(n_sec))

        k1(1:2) =  self%ham%UC%rez_lattice(:,1)
        k1(3)   =  0d0
        k2(1:2) =  self%ham%UC%rez_lattice(:,2)
        k2(3)   =  0d0

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

    function vol_k_space_parallelo(self) result(vol)
        implicit none
    class(k_space)         :: self
        real(8)                :: vol, k1(3), k2(3)

        k1      = 0d0
        k2      = 0d0
        k1(1:2) = self%ham%UC%rez_lattice(:,1)
        k2(1:2) = self%ham%UC%rez_lattice(:,2)

        vol = my_norm2(cross_prod(k1,k2))
    end function vol_k_space_parallelo

    function vol_k_hex(self) result(vol)
        implicit none
    class(k_space)    :: self
        real(8), parameter:: deg_30 = 30.0/180.0*PI, deg_60 = 60.0/180.0*PI
        real(8)           :: a, l, vol

        l = my_norm2(self%ham%UC%rez_lattice(:,1))
        a = l / (2.0 * cos(deg_30))

        vol =  3d0 * a * a * sin(deg_60)
    end function vol_k_hex

    subroutine calc_hall_conductance(self, hall)
        implicit none
    class(k_space)          :: self
        real(8)                 :: k(3), Emax
        real(8), allocatable    :: eig_val_all(:,:), eig_val_new(:,:),&
            hall(:), hall_old(:)
        integer(4), allocatable :: omega_kidx_all(:), omega_kidx_new(:)
        type(r8arr), allocatable:: omega_z_all(:), omega_z_new(:)
        integer(4)  :: N_k, k_idx, first, last, n_atm,&
            iter, cnt
        character(len=300)      :: filename

        call self%setup_berry_inte_grid()
        N_k = size(self%new_k_pts, 2)

        N_k = size(self%new_k_pts, 2)
        n_atm =  self%ham%UC%num_atoms
        allocate(self%ham%del_H(2*n_atm,2*n_atm))
        allocate(hall_old(size(self%E_fermi)))
        allocate(eig_val_all(2*n_atm,0))
        allocate(omega_z_all(0))
        allocate(omega_kidx_all(0))
        allocate(self%all_k_pts(3,0))

        iter =  1
        do iter =1,self%berry_iter
            N_k = size(self%new_k_pts, 2)

            call my_section(self%me, self%nProcs, N_k, first, last)
            allocate(eig_val_new(2*n_atm,last-first+1))
            allocate(omega_z_new(last-first+1))
            allocate(omega_kidx_new(last-first+1))

            ! calculate 
            cnt =  1
            do k_idx = first, last
                !if(self%me == root) write (*,*) k_idx, "of", last
                k = self%new_k_pts(:,k_idx)
                omega_kidx_new(cnt) =  k_idx + size(self%all_k_pts,2)
                call self%ham%calc_berry_z(k, omega_z_new(cnt)%arr,&
                                                   eig_val_new(:,cnt))
                cnt = cnt + 1
            enddo
            
            !write (filename, "(A,I0.5,A,I0.5,A)") "pre_omega_kidx_new=", self%me, "_iter=", iter, ".npy"
            !call save_npy(trim(self%prefix) // trim(filename), omega_kidx_new)

            ! concat to old ones
            call add_new_kpts_omega(omega_kidx_all, omega_kidx_new,&
                omega_z_all, omega_z_new)
            call self%append_kpts()
            call append_eigval(eig_val_all, eig_val_new)
            ! integrate hall conductance

            if(allocated(hall))then
                hall_old = hall
            else 
                hall_old =  1d35
            endif

            call self%integrate_hall(omega_kidx_all, omega_z_all, eig_val_all, hall)

            if(my_norm2(hall - hall_old)/(1d0*size(hall)) &
                                             < self%berry_conv_crit) then
                if(self%me == root) write (*,*) "Converged berry interation"
                exit
            else
                if(self%me == root) write (*,*) iter, "nkpts", size(self%all_k_pts,2), "err", &
                         my_norm2(hall - hall_old)/(1d0*size(hall))
            endif

            if(self%me == root) then
                write (filename, "(A,I0.5,A,I0.5,A)") "kpoint_proc=", self%me, "_iter=", iter, ".npy"
                call save_npy(trim(self%prefix) // trim(filename), self%all_k_pts)
                
                write (filename, "(A,I0.5,A,I0.5,A)") "elem_proc=", self%me, "_iter=", iter, ".npy"
                call save_npy(trim(self%prefix) // trim(filename), self%elem_nodes)
                
                write (filename, "(A,I0.5,A,I0.5,A)") "kweigh_proc=", self%me, "_iter=", iter, ".npy"
                call save_npy(trim(self%prefix) // trim(filename), self%weights)
                
                write (filename, "(A,I0.5,A,I0.5,A)") "hallweigh_proc=", self%me, "_iter=", iter, ".npy"
                call save_npy(trim(self%prefix) // trim(filename), self%hall_weights)
                
                write (filename, "(A,I0.5,A,I0.5,A)") "hall_proc=", self%me, "_iter=", iter, ".npy"
                call save_npy(trim(self%prefix) // trim(filename), hall)
                
                write (filename, "(A,I0.5,A,I0.5,A)") "nkpts_proc=", self%me, "_iter=", iter, ".npy"
                call save_npy(trim(self%prefix) // trim(filename), [size(self%all_k_pts,2)])

                write (filename, "(A,I0.5,A,I0.5,A)") "hallold_proc=", self%me, "_iter=", iter, ".npy"
                call save_npy(trim(self%prefix) // trim(filename), hall_old)
            endif
            
            write (filename, "(A,I0.5,A,I0.5,A)") "omega_kidx_all=", self%me, "_iter=", iter, ".npy"
            call save_npy(trim(self%prefix) // trim(filename), omega_kidx_all)

            write (filename, "(A,I0.5,A,I0.5,A)") "eigval_proc=", self%me, "_iter=", iter, ".npy"
            call save_npy(trim(self%prefix) // trim(filename), eig_val_all)
            
            write (filename, "(A,I0.5,A,I0.5,A)") "omKidx_proc=", self%me, "_iter=", iter, ".npy"
            call save_npy(trim(self%prefix) // trim(filename), omega_kidx_all)
            
            call self%set_hall_weights(omega_z_all, omega_kidx_all)
            call self%add_kpts_iter(self%kpts_per_step*self%nProcs, self%new_k_pts)
        enddo
        if(self%me == root) then
            write (*,*) iter, size(self%all_k_pts,2), &
                "saving hall cond with questionable unit"

            call save_npy(trim(self%prefix) // "hall_cond.npy", hall)
            call save_npy(trim(self%prefix) // "hall_E.npy", &
                self%E_fermi / self%units%energy)
        endif
        deallocate(hall)
        deallocate(self%ham%del_H)
        deallocate(self%new_k_pts)
    end subroutine calc_hall_conductance

    subroutine integrate_hall(self, omega_kidx_all, omega_z_all, eig_val_all, hall)
        implicit none
    class(k_space)          :: self
        integer(4), intent(in)  :: omega_kidx_all(:)
        type(r8arr), intent(in) :: omega_z_all(:)
        real(8), intent(in)     :: eig_val_all(:,:)
        real(8), allocatable    :: hall(:), omega_z(:)
        integer(4)              :: loc_idx, n, n_hall, k_idx, ierr(2)

        if(allocated(hall)) deallocate(hall)
        allocate(hall(size(self%E_fermi)))

        !run triangulation
        call run_triang(self%all_k_pts, self%elem_nodes)
        if(self%me == root)write (*,*) "elem integr", shape(self%elem_nodes)
        call self%set_weights_ksp()
        
        !perform integration with all points
        hall =  0d0

        do loc_idx = 1,size(omega_kidx_all)
            k_idx =  omega_kidx_all(loc_idx)

            omega_z =  omega_z_all(loc_idx)%arr
            do n = 1,size(omega_z)
                do n_hall =  1,size(hall)
                    hall(n_hall) = hall(n_hall) + &
                        self%weights(k_idx) * omega_z(n) *&
                        self%fermi_distr(eig_val_all(n, loc_idx), n_hall)
                enddo
            enddo
            deallocate(omega_z)
        enddo

        hall = hall / (2d0*PI)
        
        ! Allreduce is not suitable for convergence criteria
        if(self%me == root) then
            call MPI_Reduce(MPI_IN_PLACE, hall, size(hall), MPI_REAL8, MPI_SUM, &
                            root, MPI_COMM_WORLD, ierr(1))
        else
            call MPI_Reduce(hall, hall, size(hall), MPI_REAL8, MPI_SUM, &
                            root, MPI_COMM_WORLD, ierr(1))
        endif
        call MPI_Bcast(hall, size(hall), MPI_REAL8, root, MPI_COMM_WORLD, ierr(2))
        call check_ierr(ierr, self%me, "Hall conductance")
    end subroutine integrate_hall

    subroutine set_hall_weights(self, omega_z_all, omega_kidx_all)
        implicit none
    class(k_space)         :: self
        integer(4), intent(in) :: omega_kidx_all(:)
        type(r8arr), intent(in):: omega_z_all(:)
        integer(4)             :: i, n_elem, node, k_idx, loc_idx, ierr(2)
        real(8)                :: kpt(3), om_max
        character(len=300)     :: filename

        n_elem = size(self%elem_nodes,1)
        if(allocated(self%hall_weights)) then
            if(size(self%hall_weights) /= n_elem) then
                deallocate(self%hall_weights)
            endif
        endif
        if(.not. allocated(self%hall_weights)) then
            allocate(self%hall_weights(n_elem))
        endif

        self%hall_weights = 0d0
        om_max =  0d0

        do i=1,n_elem
            do node=1,3
                k_idx = self%elem_nodes(i,node)
                loc_idx =  find_list_idx(omega_kidx_all, k_idx)

                if(loc_idx > 0) then
                    kpt = self%all_k_pts(:,k_idx)
                    self%hall_weights(i) = self%hall_weights(i) &
                        + self%weights(k_idx) &!* abs(test_func(kpt(1:2)))
                        * sum(abs(omega_z_all(loc_idx)%arr))
                    !if(om_max < sum(abs(omega_z_all(loc_idx)%arr))) &
                       !om_max =  sum(abs(omega_z_all(loc_idx)%arr))
                endif
            enddo
        enddo

        !do i = 1,self%nProcs
            !if(self%me == i-1) write (*,*) self%me, "->", om_max
            !call MPI_Barrier(MPI_COMM_WORLD, ierr(1))
        !enddo


        if(self%me == root) then
            call MPI_Reduce(MPI_IN_PLACE, self%hall_weights, n_elem,&
                            MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
        else
            call MPI_Reduce(self%hall_weights,self%hall_weights, n_elem, &
                            MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr(1))
        endif
        call MPI_Bcast(self%hall_weights, n_elem, MPI_REAL8, &
                        root, MPI_COMM_WORLD, ierr(2))
        call check_ierr(ierr, self%me, "set_hall_weights")
        
        !write (filename, "(I0.5,A)") self%me,"_loc_ks.npy"
        !call save_npy(trim(self%prefix) // trim(filename), omega_kidx_all)
        
        !write (filename, "(I0.5,A)") self%me,"_hall_w.npy"
        !call save_npy(trim(self%prefix) // trim(filename), self%hall_weights)
        
        !write (filename, "(I0.5,A)") self%me,"_k_w.npy"
        !call save_npy(trim(self%prefix) // trim(filename), self%weights)
    end subroutine set_hall_weights
    

    subroutine append_eigval(eig_val_all, eig_val_new)
        implicit none
        real(8), allocatable    :: eig_val_all(:,:), eig_val_new(:,:), tmp(:,:)
        integer(4) :: vec_sz, num_k_all, num_k_new, i,j

        vec_sz =  size(eig_val_new, 1)
        if(.not. allocated(eig_val_all)) allocate(eig_val_all(vec_sz,0))
        num_k_all =  size(eig_val_all, 2)
        num_k_new =  size(eig_val_new, 2)

        allocate(tmp(vec_sz, num_k_all))
        forall(i=1:vec_sz, j=1:num_k_all) tmp(i,j) =  eig_val_all(i,j)
        deallocate(eig_val_all)

        allocate(eig_val_all(vec_sz, num_k_all + num_k_new))
        forall(i=1:vec_sz, j=1:num_k_all) eig_val_all(i,j) = tmp(i,j)
        deallocate(tmp)

        forall(i=1:vec_sz, j=1:num_k_new) eig_val_all(i,j+num_k_all) &
                = eig_val_new(i,j)
        deallocate(eig_val_new)
    end subroutine append_eigval

    subroutine add_new_kpts_omega(omega_kidx_all, omega_kidx_new, &
            omega_z_all, omega_z_new)
        implicit none
        integer(4), allocatable   :: omega_kidx_all(:),omega_kidx_new(:) 
        type(r8arr), allocatable  :: omega_z_new(:), omega_z_all(:), tmp_z(:)
        integer(4)                :: n_old, i

        if(allocated(omega_kidx_all)) then
            omega_kidx_all =  [omega_kidx_all, omega_kidx_new]
        else
            omega_kidx_all =  omega_kidx_new
        endif
        deallocate(omega_kidx_new)

        if(.not. allocated(omega_z_all)) then
            allocate(omega_z_all(0))
        endif

        n_old =  size(omega_z_all)
        allocate(tmp_z(n_old))

        do i = 1,n_old
            tmp_z(i) = omega_z_all(i)
            deallocate(omega_z_all(i)%arr)
        enddo

        deallocate(omega_z_all)
        allocate(omega_z_all(n_old + size(omega_z_new)))

        do i = 1,n_old
            omega_z_all(i) =  tmp_z(i)
            deallocate(tmp_z(i)%arr)
        enddo
        deallocate(tmp_z)

        do i = 1,size(omega_z_new)
            omega_z_all(i+n_old) = omega_z_new(i)
            deallocate(omega_z_new(i)%arr)
        enddo
        deallocate(omega_z_new)
    end subroutine add_new_kpts_omega

    subroutine setup_berry_inte_grid(self)
        implicit none
    class(k_space)           :: self
        integer(4)           :: ierr 

        if(allocated(self%new_k_pts) )then
            deallocate(self%new_k_pts)
        endif

        if(trim(self%ham%UC%uc_type) == "square_2d") then
            call self%setup_inte_grid_square(self%berry_num_k_pts)
        elseif(trim(self%ham%UC%uc_type) == "honey_2d") then
            call self%setup_inte_grid_hex(self%berry_num_k_pts)
        else
            if(self%me ==  root) write (*,*) "berry k-grid not known"
            call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
        endif

    end subroutine setup_berry_inte_grid

    subroutine plot_omega(self)
        implicit none
    class(k_space)         :: self
        real(8), allocatable   :: omega_z(:,:), tmp_vec(:), sec_omega_z(:,:),&
                                  eig_val(:)
        real(8)                :: k(3)
        integer(4)  :: N, k_idx, send_count, first, last, ierr, cnt, n_atm
        integer(4), allocatable:: num_elems(:), offsets(:)
        integer(4)  :: dim_sz

        dim_sz =  self%berry_num_k_pts
        n_atm =  self%ham%UC%num_atoms
        allocate(self%ham%del_H(2*n_atm,2*n_atm))
        allocate(num_elems(self%nProcs))
        allocate(offsets(self%nProcs))
        N = 2* self%ham%UC%num_atoms
        call self%setup_inte_grid_hex(dim_sz)

        if(self%me ==  root) write (*,*) "nkpts =  ", size(self%new_k_pts, 2)

        allocate(eig_val(N))
        allocate(omega_z(N, size(self%new_k_pts, 2)))
        allocate(tmp_vec(N))

        call sections(self%nProcs, size(self%new_k_pts, 2), num_elems, offsets)
        call my_section(self%me, self%nProcs, size(self%new_k_pts, 2), first, last)
        num_elems =  num_elems * N
        offsets   =  offsets   * N
        send_count =  N *  (last - first + 1)
        allocate(sec_omega_z(N, send_count))

        cnt =  1
        do k_idx = first,last
            write (*,*) "k_ind" , k_idx
            k = self%new_k_pts(:,k_idx)

            call self%ham%calc_berry_z(k, tmp_vec, eig_val)

            sec_omega_z(:,cnt) =  tmp_vec 
            cnt = cnt + 1
        enddo

        call MPI_Gatherv(sec_omega_z, send_count, MPI_REAL8, &
            omega_z,    num_elems, offsets, MPI_REAL8, &
            root, MPI_COMM_WORLD, ierr)

        if(self%me ==  root) then
            call save_npy(trim(self%prefix) //  "omega_xy_z.npy", omega_z)
            call save_npy(trim(self%prefix) //  "omega_xy_k.npy", self%new_k_pts)
            write (*,*) "Berry curvature saved unitless"
        endif
        deallocate(self%ham%del_H)
        deallocate(eig_val)
    end subroutine plot_omega 

    subroutine set_fermi(self, cfg)
        implicit none
    class(k_space)         :: self
    class(CFG_t)           :: cfg
        real(8)                :: tmp(3)
        integer(4)             :: ierr, n_steps

        if(root == self%me) then
            call CFG_get(cfg, "dos%E_fermi", tmp)
        endif
        call MPI_Bcast(tmp, 3, MPI_REAL8, root, MPI_COMM_WORLD, ierr)
        n_steps = nint(tmp(3))
        tmp =  tmp *  self%units%energy

        call linspace(tmp(1), tmp(2), n_steps, self%E_fermi)
        call self%write_fermi()
    end subroutine set_fermi

    subroutine find_fermi(self, cfg)
        implicit none
    class(k_space)         :: self
    class(CFG_t)           :: cfg
        real(8)                :: target, delta_old, delta_new
        integer(4)             :: i, ierr

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
        real(8)                :: fermi(1)
        if(self%me ==  root) then 
            fermi =  self%E_fermi

            call save_npy(trim(self%prefix) //  "fermi.npy", fermi)
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
        integer(4), intent(in)        :: n_ferm
        real(8)                       :: ferm, exp_term


        exp_term =  (E - self%E_fermi(n_ferm)) /&
                     (boltzmann_const * self%temp)
        if(exp_term > 700d0) then
            ferm =  0d0
        else 
            ferm = 1d0 / (exp(exp_term) + 1d0)
        endif
    end function fermi_distr

    function find_E_max(self) result(c)
        implicit none
    class(k_space), intent(in)   :: self
        real(8)                      :: l, u, c
        real(8), parameter           :: tol = 1d-6, tar = 1d-12
        integer(4)                   :: Emax, cnt

        l = - 25d0 * self%ham%t_nn
        u = 25d0 * self%ham%t_nn
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
        integer(4), intent(in)     :: idx
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
        integer(4), intent(in)   :: idx
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

    function centeroid_of_elem(self, idx, k_pts) result(centeroid)
        implicit none
    class(k_space), intent(in)  :: self
        integer(4), intent(in)      :: idx
        integer(4)                  :: i
        real(8), intent(in)         :: k_pts(:,:)
        real(8)                     :: centeroid(2)

        centeroid =  0d0 
        do i = 1,3
            centeroid(1) =  centeroid(1) + k_pts(1,self%elem_nodes(idx,i))
            centeroid(2) =  centeroid(2) + k_pts(2,self%elem_nodes(idx,i))
        enddo
        centeroid = centeroid * 0.33333333333d0
    end function centeroid_of_elem

    subroutine pad_k_points_init(self)
        implicit none
    class(k_space)                :: self
        integer(4)  :: rest, i,j, cnt, n_kpts, n_elem, ierr, per_proc
        integer(4), allocatable :: sort(:)
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
                write (*,*) "Not enough elements for inital padding: ", n_elem
                stop
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
                        write (*,*) "not enough elems for refinement"
                        call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
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
        integer(4)           :: i

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
        integer(4)                   :: idx
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

    subroutine hex_laplace_smoother(self)
        implicit none
    class(k_space)              :: self
        real(8), allocatable        :: shift(:,:), n_neigh(:)
        integer(4)                  :: i, j, point, neigh
        integer(4), parameter       :: N = 4

        allocate(shift(2, size(self%new_k_pts,2)))
        allocate(n_neigh(size(self%new_k_pts,2)))


        do i =  1,self%laplace_iter
            shift = 0d0
            neigh = 0d0
            do j = 1,size(self%elem_nodes,1)
                do point = 1,3
                    do neigh = 1,3
                        if(.not. self%on_hex_border(self%elem_nodes(j,point)))then
                            shift(:,self%elem_nodes(j,point)) = shift(:,self%elem_nodes(j,point)) + &
                                (self%new_k_pts(1:2,self%elem_nodes(j,neigh)) - self%new_k_pts(1:2,self%elem_nodes(j,point)))
                            n_neigh(self%elem_nodes(j,point)) =  n_neigh(self%elem_nodes(j,point)) + 1d0
                        endif
                    enddo
                enddo
            enddo
            do j = 1,size(shift,2)
                if(n_neigh(j) >  1d-6) then
                    self%new_k_pts(1:2,j) =  self%new_k_pts(1:2,j) + shift(:,j) / n_neigh(j)
                endif

            enddo
        enddo

    end subroutine hex_laplace_smoother

    subroutine add_kpts_iter(self, n_new, new_ks)
        implicit none
    class(k_space)          :: self
        integer(4), intent(in)  :: n_new
        real(8), allocatable    :: new_ks(:,:)
        integer(4)              :: n_elem, i, cnt, ierr
        integer(4), allocatable :: sort(:)
        character(len=300)      :: filename

        if(allocated(new_ks)) then
            if(size(new_ks,1) /= 3 .or. size(new_ks,2) /= n_new) then
                deallocate(new_ks)
            endif
        endif
        if(.not. allocated(new_ks)) allocate(new_ks(3, n_new))
        n_elem = size(self%elem_nodes,1)

        call qargsort(self%hall_weights, sort)
        if(self%me == root) write (*,*) "size sort", shape(sort)

        !write (filename, "(I0.4,A,I0.9,A)") self%me,&
                                !"_sorting_", size(sort), ".npy"
        !call save_npy(trim(self%prefix) // filename, sort)
        
        !write (filename, "(I0.4,A,I0.9,A)") self%me,&
                                !"_hall_weights_", size(sort), ".npy"
        !call save_npy(trim(self%prefix) // filename, self%hall_weights)

        !write (filename, "(I0.4,A,I0.9,A)") self%me,&
                                !"_k_weights_", size(sort), ".npy"
        !call save_npy(trim(self%prefix) // filename, self%weights)


        new_ks =  0d0
        i = n_elem

        do cnt = 1,n_new
            new_ks(1:2, cnt) = self%centeroid_of_elem(sort(i), self%all_k_pts)
            i = i - 1
            if(i == 0) then
                write (*,*) "Not enough elements"
                call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
            endif
        enddo
    end subroutine add_kpts_iter

    subroutine append_kpts(self)
        implicit none
    class(k_space)         :: self
        real(8), allocatable   :: tmp(:,:)
        integer(4)             :: old_sz, new_sz, i,j

        if(.not. allocated(self%all_k_pts)) allocate(self%all_k_pts(3,0))
        old_sz = size(self%all_k_pts, 2)
        new_sz = size(self%new_k_pts, 2)

        allocate(tmp(3, old_sz))

        forall(i=1:3, j=1:old_sz) tmp(i,j) = self%all_k_pts(i,j)

        deallocate(self%all_k_pts)
        allocate(self%all_k_pts(3, old_sz + new_sz))

        forall(i=1:3, j=1:old_sz) self%all_k_pts(i,j)        = tmp(i,j)
        forall(i=1:3, j=1:new_sz) self%all_k_pts(i,old_sz+j) = self%new_k_pts(i,j)

        deallocate(tmp)
        deallocate(self%new_k_pts)
    end subroutine append_kpts

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

end module

