module Class_k_space
    use m_config
    use m_npy
    !use Class_unit_cell 
    use Class_hamiltionian
    use Class_helper
    implicit none

    type k_space
        real(8), allocatable :: k_pts(:,:)
        real(8), allocatable :: int_DOS(:) !> integrated Density of states  
        real(8), allocatable :: E_DOS(:)
        real(8) :: DOS_gamma !> broadening \f$ \Gamma \f$ used in
        !> DOS calculations
        real(8) :: DOS_lower !> lower energy bound for DOS calc
        real(8) :: DOS_upper !> upper energy bound for DOS calc
        real(8) :: E_fermi !> Fermi lvl
        real(8) :: temp !> temperature used in fermi-dirac
        integer(4) :: DOS_num_k_pts !> number of kpts per dim 
        !> used in DOS calculations
        integer(4) :: berry_num_k_pts !> number of ks per dim in berry calc
        integer(4) :: num_DOS_pts!> number of points on E grid
        integer(4) :: num_k_pts !> number of k_pts per segment
        real(8), allocatable :: k1_param(:) !> 1st k_space param
        real(8), allocatable :: k2_param(:) !> 2nd k_space param
        character(len=300)    :: filling, prefix
        logical :: perform_dos_integration !> param to toggle dos integr.
        type(hamil)           :: ham
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
    end type k_space 

contains
    Subroutine  calc_and_print_band(self, cfg)
        Implicit None
        class(k_space)                :: self 
        class(CFG_t)                  :: cfg
        character(len=300)            :: npz_file
        real(8), dimension(:,:), allocatable    :: eig_val

        npz_file = trim(self%prefix) // ".npz"

        if(trim(self%filling) ==  "path_rel") then
            call self%setup_k_path_rel()
        else if(trim(self%filling) == "path_abs") then
            call self%setup_k_path_abs(cfg)
        else if(trim(self%filling) == "grid") then
            call self%setup_k_grid(cfg)
        else
            write (*,*) "Filling not known"
            stop
        endif

        call self%ham%calc_eigenvalues(self%k_pts, eig_val)
        call add_npz(npz_file, "band_k", self%k_pts)
        call add_npz(npz_file, "band_E", eig_val)
        call add_npz(npz_file, "lattice", self%ham%UC%lattice)
        call add_npz(npz_file, "rez_lattice", self%ham%UC%rez_lattice)
        call add_npz(npz_file, "band_num_kpts", (/ self%num_k_pts /))

        deallocate(self%k_pts)
        deallocate(eig_val)
    End Subroutine calc_and_print_band

    subroutine calc_pdos(self, E, PDOS)
        implicit none
        class(k_space)          :: self
        real(8), intent(in)     :: E(:)
        real(8), intent(out)    :: PDOS(:,:)
        real(8), allocatable    :: RWORK(:), eig_val(:)
        complex(8), allocatable :: H(:,:), WORK(:)
        integer(4), allocatable :: IWORK(:)
        integer(4)              :: k_idx, E_idx, j, m, N, LWMAX, info 

        N =  2 * self%ham%UC%num_atoms
        if(N > 4) then 
            LWMAX =  4 * N*N
        else
            LWMAX =  200
        endif

        allocate(H(N,N))
        allocate(eig_val(N))
        allocate(WORK(LWMAX))
        allocate(RWORK(LWMAX))
        allocate(IWORK(LWMAX))

        PDOS =  0d0
        do k_idx=1,size(self%k_pts, 2)
            call self%ham%setup_H(self%k_pts(:,k_idx), H)
            call zheevd('V', 'U', N, H, N, eig_val, WORK, LWMAX, &
                RWORK, LWMAX, IWORK, LWMAX, info)
            if( info /= 0) then
                write (*,*) "ZHEEVD (with vectors) failed: ", info
                stop
            endif

            do E_idx =  1,self%num_DOS_pts 
                ! eigenvectors are stored column-wise
                ! m-th eigenvalue
                ! j-th component of 
                do m =  1,N
                    do j = 1,N
                        PDOS(j, E_idx) = PDOS(j, E_idx) &
                                       + self%lorentzian(E(E_idx) - eig_val(m)) &
                                       * H(j,m) * conjg(H(j,m))
                    enddo
                enddo
            enddo
        enddo
        PDOS =  PDOS / real(size(self%k_pts,2))

        deallocate(WORK)
        deallocate(IWORK)
        deallocate(RWORK)
    end subroutine


    subroutine calc_and_print_dos(self)
        implicit none
        class(k_space)       :: self
        character(len=300)   :: npz_file
        real(8), allocatable :: DOS(:), PDOS(:,:), up(:), down(:)
        real(8)              :: dE
        integer(4)           :: i, num_atoms


        npz_file = trim(self%prefix) // ".npz"
        call self%setup_inte_grid_square(self%DOS_num_k_pts)
        num_atoms =  self%ham%UC%num_atoms
        allocate(PDOS(2*num_atoms, self%num_DOS_pts))

        self%E_DOS =  linspace(self%DOS_lower, self%DOS_upper, self%num_DOS_pts)
        allocate(DOS(self%num_DOS_pts))
        allocate(up(self%num_DOS_pts))
        allocate(down(self%num_DOS_pts))

        call self%calc_pdos(self%E_DOS, PDOS)

        DOS  = sum(PDOS,1)
        up   = sum(PDOS(1:num_atoms, :),1)
        down = sum(PDOS(num_atoms+1:2*num_atoms, :),1)

        call add_npz(npz_file, "DOS_E",       self%E_DOS)
        call add_npz(npz_file, "DOS",         DOS)
        call add_npz(npz_file, "DOS_partial", PDOS)
        call add_npz(npz_file, "DOS_up",      up)
        call add_npz(npz_file, "DOS_down",    down)

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
        call add_npz(npz_file, "DOS_integrated", self%int_DOS)
        
        deallocate(self%k_pts)
        deallocate(DOS)
        deallocate(PDOS)
        deallocate(up)
        deallocate(down)
    end subroutine calc_and_print_dos

    function init_k_space(cfg) result(self)
        implicit none
        type(k_space)         :: self
        type(CFG_t)           :: cfg
        real(8)               :: tmp
        integer(4)            :: sz
        character(len=300)    :: npz_file

        self%ham =  init_hamil(cfg)

        call CFG_get(cfg, "output%band_prefix", self%prefix)
        call CFG_get(cfg, "band%filling", self%filling)

        call CFG_get(cfg, "dos%delta_broadening", tmp)
        self%DOS_gamma =  tmp *  get_unit_conv("energy", cfg)

        call CFG_get(cfg, "dos%num_points", self%num_DOS_pts)

        npz_file = trim(self%prefix) // ".npz"
        call self%ham%UC%save_unit_cell(npz_file)

        call CFG_get_size(cfg, "band%k_x", sz)
        allocate(self%k1_param(sz))
        call CFG_get_size(cfg, "band%k_y", sz)
        allocate(self%k2_param(sz))

        call CFG_get(cfg, "band%k_x", self%k1_param)
        call CFG_get(cfg, "band%k_y", self%k2_param)
        call CFG_get(cfg, "band%num_points", self%num_k_pts)

        call CFG_get(cfg, "dos%k_pts_per_dim", self%DOS_num_k_pts)
        call CFG_get(cfg, "dos%lower_E_bound", tmp)
        self%DOS_lower =  tmp * get_unit_conv("energy", cfg)
        call CFG_get(cfg, "dos%upper_E_bound", tmp)
        self%DOS_upper =  tmp * get_unit_conv("energy", cfg)

        call CFG_get(cfg, "berry%k_pts_per_dim", self%berry_num_k_pts)
        call CFG_get(cfg, "berry%temperature", tmp)
        self%temp = tmp * get_unit_conv("temperature", cfg)
    end function init_k_space

    subroutine setup_k_grid(self, cfg)
        implicit none
        class(k_space)           :: self
        type(CFG_t)          :: cfg
        real(8), allocatable :: kx_points(:), ky_points(:)
        real(8), allocatable :: kx_grid(:,:), ky_grid(:,:)
        integer(4)           :: sz_x, sz_y, i, j


        sz_x =  NINT(self%k1_param(3))
        sz_y =  NINT(self%k2_param(3))

        self%k1_param(1:2) =  self%k1_param(1:2) &
            * get_unit_conv("inv_length",cfg)
        self%k2_param(1:2) =  self%k2_param(1:2) &
            * get_unit_conv("inv_length",cfg)

        allocate(kx_grid(sz_x, sz_y))
        allocate(ky_grid(sz_x, sz_y))
        allocate(kx_points(sz_x))
        allocate(ky_points(sz_y))

        kx_points =  linspace(self%k1_param(1), self%k1_param(2), sz_x)
        ky_points =  linspace(self%k2_param(1), self%k2_param(2), sz_y)

        do j =  0,sz_y-1
            kx_grid(:,j+1) =  kx_points
        enddo

        do i =  0,sz_x-1
            ky_grid(i+1,:) =  ky_points
        enddo

        allocate(self%k_pts(3, sz_x* sz_y))
        self%k_pts(1,:) =  reshape(kx_grid, (/sz_x * sz_y /))
        self%k_pts(2,:) =  reshape(ky_grid, (/sz_x * sz_y /))
        self%k_pts(3,:) =  0.0d0

        deallocate(kx_grid)
        deallocate(ky_grid)
        deallocate(kx_points)
        deallocate(ky_points)
        write (*,*) "K_grid sz: ", sz_x, sz_y

    end subroutine setup_k_grid

    subroutine setup_inte_grid_square(self, n_k)
        implicit none
        class(k_space)        :: self
        integer(4), intent(in):: n_k
        real(8), allocatable  :: ls(:), ls_help(:)
        real(8)               :: k1(3), k2(3)
        integer(4)            :: i, j, cnt

        allocate(self%k_pts(3,n_k**2))

        k1 =  0d0
        k2 =  0d0 
        k1(1:2) =  self%ham%UC%rez_lattice(:,1)
        k2(1:2) =  self%ham%UC%rez_lattice(:,2)

        ls_help    = linspace(0d0, 1d0, n_k+1)
        ls         = ls_help(1:n_k)
        self%k_pts = 0d0
        cnt        = 1

        do i = 1,n_k
            do j =  1,n_k
                self%k_pts(:,cnt) = ls(i) *  k1 +  ls(j) *  k2
                cnt =  cnt + 1
            enddo
        enddo
    end subroutine setup_inte_grid_square

    subroutine setup_k_path_abs(self, cfg)
        implicit none
        class(k_space)        :: self
        class(CFG_t)          :: cfg
        integer(4)            :: n_pts, n_sec, start, halt, i


        self%k1_param =  self%k1_param * get_unit_conv("inv_length",cfg)
        self%k2_param =  self%k2_param * get_unit_conv("inv_length",cfg)
        n_pts =  self%num_k_pts
        write (*,*) "k1_param: ", self%k1_param
        write (*,*) "k2_param: ", self%k2_param

        n_sec =  size(self%k1_param)-1
        write (*,*) "N_sec: ", n_sec


        allocate(self%k_pts(3, n_sec * (self%num_k_pts-1) + 1))
        self%k_pts(3,:) =  0d0

        start = 1
        do i =  1,n_sec
            halt =  start +  n_pts - 1

            self%k_pts(1,start:halt) &
                = linspace(self%k1_param(i), self%k1_param(i+1), n_pts)
            self%k_pts(2,start:halt) &
                = linspace(self%k2_param(i), self%k2_param(i+1), n_pts)
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

        allocate(self%k_pts(3,n_sec * (n_pts - 1) + 1))
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

            c1_sec = linspace(self%k1_param(i),self% k1_param(i+1), n_pts)
            c2_sec = linspace(self%k2_param(i), self%k2_param(i+1), n_pts)

            cnt =  1
            do j =  start,halt
                self%k_pts(:,j) =  c1_sec(cnt) *  k1 +  c2_sec(cnt) *  k2
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

        vol = norm2(cross_prod(k1,k2))
    end function vol_k_space_parallelo

    function calc_hall_conductance(self) result(hall)
        implicit none
        class(k_space)       :: self
        real(8)              :: hall, V_k, k(3)
        real(8), allocatable :: eig_val(:), omega_z(:), omega_plot(:,:)
        integer(4)           :: N_k, n, k_idx
        character(len=300)   :: npz_file

        if(allocated(self%k_pts) )then
            deallocate(self%k_pts)
        endif
        call self%setup_inte_grid_square(self%berry_num_k_pts)
        V_k = self%vol_k_space_parallelo()
        N_k = size(self%k_pts, 2)

        allocate(omega_plot(2*self%ham%UC%num_atoms, N_k))

        hall = 0d0
        do k_idx = 1,N_k
            k = self%k_pts(:,k_idx)
            call self%ham%calc_berry_z(k, omega_z, eig_val)
            omega_plot(:,k_idx) =  omega_z
            do n = 1,2*self%ham%UC%num_atoms
                hall = hall + omega_z(n) * self%fermi_distr(eig_val(n))
            enddo
        enddo
        npz_file = trim(self%prefix) // ".npz"
        write (*,*) "Wrote to: ", npz_file
        call add_npz(npz_file, "berry_plot", omega_plot)
        call add_npz(npz_file, "berry_k", self%k_pts)
        hall = hall * V_k/real(N_k)
        hall = hall / (2d0*PI)
    end function calc_hall_conductance

    subroutine set_fermi(self, cfg)
        implicit none
        class(k_space)         :: self
        class(CFG_t)           :: cfg
        real(8)                :: tmp

        call CFG_get(cfg, "dos%E_fermi", tmp)
        self%E_fermi =  tmp * get_unit_conv("energy", cfg)

        call self%write_fermi()
    end subroutine set_fermi

    subroutine find_fermi(self, cfg)
        implicit none
        class(k_space)         :: self
        class(CFG_t)           :: cfg
        real(8)                :: target, delta_old, delta_new
        integer(4)             :: i

        call CFG_get(cfg, "dos%fermi_fill", target)
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
        character(len=300)     :: npz_file
        real(8)                :: fermi(1)

        fermi =  self%E_fermi

        npz_file = trim(self%prefix) // ".npz"
        call add_npz(npz_file, "E_fermi", fermi)
    end subroutine write_fermi

    function lorentzian(self, x) result(lor)
        implicit none
        class(k_space), intent(in)    :: self
        real(8), intent(in)           :: x
        real(8)                       :: lor

        lor =  self%DOS_gamma &
            / (PI * (x**2 +  self%DOS_gamma**2))
    end function lorentzian

    function fermi_distr(self, E) result(ferm)
        implicit none
        class(k_space), intent(in)    :: self
        real(8), intent(in)           :: E
        real(8)                       :: ferm

        ferm = 1d0 / (exp((E-self%E_fermi)&
                      /(boltzmann_const * self%temp)) + 1d0)
    end function fermi_distr

end module

