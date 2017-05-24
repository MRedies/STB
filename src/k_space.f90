module Class_k_space
    use m_config
    use m_npy
    !use Class_unit_cell 
    use Class_hamiltionian
    use Class_math_helper
    implicit none

    type k_space
        real(8), dimension(:,:), allocatable          :: k_pts
        real(8) :: DOS_gamma !> broadening \f$ \Gamma \f$ used in
        !> DOS calculations
        real(8) :: DOS_lower !> lower energy bound for DOS calc
        real(8) :: DOS_upper !> upper energy bound for DOS calc
        integer(4) :: DOS_num_k_pts !> number of kpts per dim 
        !> used in DOS calculations
        integer(4) :: num_DOS_pts!> number of points on E grid
        integer(4) :: num_k_pts !> number of k_pts per segment
        real(8), allocatable :: k1_param(:) !> 1st k_space param
        real(8), allocatable :: k2_param(:) !> 2nd k_space param
        character(len=300)    :: filling, prefix
        logical :: perform_dos_integration !> param to toggle dos integr.
        type(hamil)           :: ham
    contains
        procedure :: calc_and_print_band => calc_and_print_band
        procedure :: setup_k_path_rel    => setup_k_path_rel
        procedure :: setup_k_path_abs    => setup_k_path_abs
        procedure :: setup_k_grid        => setup_k_grid
        procedure :: lorentzian          => lorentzian
        procedure :: calc_dos            => calc_dos
        procedure :: calc_pdos           => calc_pdos
        procedure :: calc_and_print_dos  => calc_and_print_dos
        procedure :: setup_DOS_grid_square => &
            setup_DOS_grid_square
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
        real(8), intent(in)     :: E
        real(8), intent(out)    :: PDOS(:)
        character(len=300)      :: npz_file
        real(8), allocatable    :: RWORK(:), eig_val(:)
        complex(8), allocatable :: H(:,:), WORK(:)
        real(8)                 :: lower, upper, k(3)
        integer(4), allocatable :: IWORK(:)
        integer(4)              :: k_idx, j, m, N, LWMAX, info 

        npz_file = trim(self%prefix) // ".npz"
        N =  2 * self%ham%UC%num_atoms
        PDOS =  0d0
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

        do k_idx=1,size(self%k_pts, 2)
            call self%ham%setup_H(self%k_pts(:,k_idx), H)
            call zheevd('V', 'U', N, H, N, eig_val, WORK, LWMAX, &
                RWORK, LWMAX, IWORK, LWMAX, info)
            if( info /= 0) then
                write (*,*) "ZHEEVD (with vectors) failed: ", info
                stop
            endif

            ! eigenvectors are stored colum-wise
            do m =  1,N
                do j = 1,N
                    PDOS(j) = PDOS(j) + self%lorentzian(E - eig_val(m)) &
                                      * H(j,m) * conjg(H(j,m))
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
        real(8), allocatable :: E(:), DOS(:), int_DOS(:), PDOS(:,:), up(:), down(:)
        real(8)              :: dE
        integer(4)           :: i, num_atoms


        npz_file = trim(self%prefix) // ".npz"
        call self%setup_DOS_grid_square()
        num_atoms =  self%ham%UC%num_atoms
        allocate(PDOS(2*num_atoms, self%num_DOS_pts))

        E =  linspace(self%DOS_lower, self%DOS_upper, self%num_DOS_pts)
        allocate(DOS(self%num_DOS_pts))
        allocate(up(self%num_DOS_pts))
        allocate(down(self%num_DOS_pts))
        !allocate(DOS,  mold=E)
        !allocate(up,   mold=E)
        !allocate(down, mold=E)

        do i =  1, self%num_DOS_pts 
            write (*,*) i, "/", self%num_DOS_pts
            call self%calc_pdos(E(i), PDOS(:,i))
        enddo

        DOS  = sum(PDOS,1)
        up   = sum(PDOS(1:num_atoms, :),1)
        down = sum(PDOS(num_atoms+1:2*num_atoms, :),1)

        call add_npz(npz_file, "DOS_E",       E)
        call add_npz(npz_file, "DOS",         DOS)
        call add_npz(npz_file, "DOS_partial", PDOS)
        call add_npz(npz_file, "DOS_up",      up)
        call add_npz(npz_file, "DOS_down",    down)

        if(self%perform_dos_integration) then
            allocate(int_DOS(self%num_DOS_pts))
            !allocate(int_DOS, mold=DOS)

            if(size(E) >=  2) then 
                dE =  E(2) - E(1)
            else 
                write (*,*) "Can't perform integration. Only one point"
                stop
            endif
            int_DOS(1) =  0d0
            do i =  2,size(E)
                int_DOS(i) =  int_DOS(i-1) &
                    + 0.5d0 * dE * (DOS(i-1) +  DOS(i))
            enddo
            call add_npz(npz_file, "DOS_integrated", int_DOS)

        endif

    end subroutine calc_and_print_dos

    subroutine calc_dos(self, E, eig_val, DOS)
        implicit none
        class(k_space)         :: self
        real(8), intent(in)    :: E(:), eig_val(:,:)
        real(8), allocatable, intent(out)   :: DOS(:)
        integer(4)             :: i, j, k, cnt 
        real(8)                :: N

        !> the DOS ist calculated using the formular:
        !> \f$ \frac{1}{N} \sum_i \delta(\epsilon - \epsilon_i )\f$

        N = size(eig_val, dim=1)

        allocate(DOS(size(E)))
        DOS =  0d0 

        do i =  1,size(E)
            cnt =  0d0 
            do j = 1,size(eig_val,dim=1)
                do k = 1,size(eig_val,dim=2)
                    cnt =  cnt +  1
                    DOS(i) = DOS(i) &
                        + self%lorentzian(E(i) - eig_val(j,k))
                enddo
            enddo
        enddo
        DOS = DOS / N
    end subroutine calc_dos 

    function init_k_space(cfg) result(k)
        implicit none
        type(k_space)         :: k
        type(CFG_t)           :: cfg
        real(8)               :: tmp
        integer(4)            :: sz
        character(len=300)    :: npz_file

        k%ham =  init_hamil(cfg)

        call CFG_get(cfg, "output%band_prefix", k%prefix)
        call CFG_get(cfg, "band%filling", k%filling)

        call CFG_get(cfg, "dos%delta_broadening", tmp)
        k%DOS_gamma =  tmp *  get_unit_conv("energy", cfg)

        call CFG_get(cfg, "dos%num_points", k%num_DOS_pts)

        npz_file = trim(k%prefix) // ".npz"
        call k%ham%UC%save_unit_cell(npz_file)

        call CFG_get_size(cfg, "band%k_x", sz)
        allocate(k%k1_param(sz))
        call CFG_get_size(cfg, "band%k_y", sz)
        allocate(k%k2_param(sz))

        call CFG_get(cfg, "band%k_x", k%k1_param)
        call CFG_get(cfg, "band%k_y", k%k2_param)
        call CFG_get(cfg, "band%num_points", k%num_k_pts)

        call CFG_get(cfg, "dos%k_pts_per_dim", k%DOS_num_k_pts)
        call CFG_get(cfg, "dos%perform_integration", &
            k%perform_dos_integration)
        call CFG_get(cfg, "dos%lower_E_bound", tmp)
        k%DOS_lower =  tmp * get_unit_conv("energy", cfg)
        call CFG_get(cfg, "dos%upper_E_bound", tmp)
        k%DOS_upper =  tmp * get_unit_conv("energy", cfg)

    end function init_k_space

    subroutine setup_k_grid(self, cfg)
        implicit none
        class(k_space)           :: self
        type(CFG_t)          :: cfg
        real(8)              :: k_mtx(2,2)
        real(8), allocatable :: kx_points(:), ky_points(:)
        real(8), allocatable :: kx_grid(:,:), ky_grid(:,:), RHS(:,:)
        integer(4)           :: sz_x, sz_y, i,j, cnt


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

    subroutine setup_DOS_grid_square(self)
        implicit none
        class(k_space)        :: self
        real(8), allocatable  :: ls(:), ls_help(:)
        real(8)               :: k1(3), k2(3)
        integer(4)            :: n_k, i, j, cnt

        n_k =  self%DOS_num_k_pts
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
    end subroutine setup_DOS_grid_square

    subroutine setup_k_path_abs(self, cfg)
        implicit none
        class(k_space)        :: self
        class(CFG_t)          :: cfg
        integer(4)            :: n, n_pts, n_sec, start, halt, i


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
        integer(4) :: n, n_pts, n_sec,i,j, start, halt, cnt

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


    Function  linspace(start, halt, n) result(x)
        Implicit None
        real(8), intent(in)    :: start, halt
        integer(4), intent(in) :: n
        real(8), dimension(n)  :: x
        real(8)                :: step, curr 
        integer(4)             :: i

        step =  (halt - start) /  (n-1)
        curr =  start

        do i = 1,n
            x(i) =  curr
            curr =  curr +  step
        enddo


    End Function 

    function lorentzian(self, x) result(lor)
        implicit none
        class(k_space), intent(in)    :: self
        real(8), intent(in)           :: x
        real(8)                       :: lor


        lor =  self%DOS_gamma &
            / (PI * (x**2 +  self%DOS_gamma**2))
    end function lorentzian

end module

