module Class_k_space
    use m_config
    use m_npy
    !use Class_unit_cell 
    use Class_hamiltionian
    implicit none

    type k_space
        real(8), dimension(:,:), allocatable          :: k_pts
        character(len=100), dimension(:), allocatable :: label 
        character(len=300)                            :: prefix
        type(hamil)                                   :: ham
    contains
        procedure :: calc_and_print_band => calc_and_print_band
        procedure :: setup_k_path_rel    => setup_k_path_rel
        procedure :: setup_k_path_abs    => setup_k_path_abs
        procedure :: setup_k_grid        => setup_k_grid
    end type k_space 

contains
    Subroutine  calc_and_print_band(this)
        Implicit None
        class(k_space)                :: this 
        character(len=300)            :: npz_file
        real(8), dimension(:,:), allocatable    :: eig_val

        npz_file = trim(this%prefix) // ".npz"

        call this%ham%calc_eigenvalues(this%k_pts, eig_val)
        call add_npz(npz_file, "k", this%k_pts)
        call add_npz(npz_file, "E", eig_val)
        call add_npz(npz_file, "lattice", this%ham%UC%lattice)
        call add_npz(npz_file, "rez_lattice", this%ham%UC%rez_lattice)
        
        call this%ham%UC%save_unit_cell(npz_file)

        deallocate(eig_val)
    End Subroutine calc_and_print_band

    function init_k_space(cfg) result(k)
        implicit none
        type(k_space)         :: k
        type(CFG_t)           :: cfg
        character(len=300)                 :: filling

        k%ham =  init_hamil(cfg)

        call CFG_get(cfg, "output%band_prefix", k%prefix)
        call CFG_get(cfg, "kspace%filling", filling)

        if(trim(filling) ==  "path_rel") then
            call k%setup_k_path_rel(cfg)
        else if(trim(filling) == "path_abs") then
            call k%setup_k_path_abs(cfg)
        else if(trim(filling) == "grid") then
            call k%setup_k_grid(cfg)
        endif
    end function init_k_space

    subroutine setup_k_grid(this, cfg)
        implicit none
        class(k_space)          :: this
        type(CFG_t)             :: cfg
        real(8), dimension(3)   :: kx_para, ky_para 
        real(8), dimension(2,2) :: k_mtx
        real(8), dimension(:), allocatable   :: kx_points, ky_points
        real(8), dimension(:,:), allocatable :: kx_grid, ky_grid, RHS
        integer(4)              :: sz_x, sz_y, i,j, cnt

        call CFG_get(cfg, "kspace%k_x", kx_para)
        call CFG_get(cfg, "kspace%k_y", ky_para)

        sz_x =  NINT(kx_para(3))
        sz_y =  NINT(ky_para(3))

        kx_para(1:2) =  kx_para(1:2) * get_unit_conv("inv_length",cfg)
        ky_para(1:2) =  ky_para(1:2) * get_unit_conv("inv_length",cfg)

        allocate(kx_grid(sz_x, sz_y))
        allocate(ky_grid(sz_x, sz_y))
        allocate(kx_points(sz_x))
        allocate(ky_points(sz_y))

        kx_points =  linspace(kx_para(1), kx_para(2), sz_x)
        ky_points =  linspace(ky_para(1), ky_para(2), sz_y)

        do j =  0,sz_y-1
            kx_grid(:,j+1) =  kx_points
        enddo

        do i =  0,sz_x-1
            ky_grid(i+1,:) =  ky_points
        enddo

        allocate(this%k_pts(3, sz_x* sz_y))
        this%k_pts(1,:) =  reshape(kx_grid, (/sz_x * sz_y /))
        this%k_pts(2,:) =  reshape(ky_grid, (/sz_x * sz_y /))
        this%k_pts(3,:) =  0.0d0

        deallocate(kx_grid)
        deallocate(ky_grid)
        deallocate(kx_points)
        deallocate(ky_points)
        write (*,*) "K_grid sz: ", sz_x, sz_y

    end subroutine setup_k_grid

    subroutine setup_k_path_abs(this, cfg)
        implicit none
        class(k_space)        :: this
        type(CFG_t)           :: cfg
        real(8), allocatable  :: k1(:), k2(:)
        integer(4)            :: n, n_pts, n_sec, start, halt, i

        call CFG_get_size(cfg, "kspace%k_x", n)
        n_sec =  n-1
        allocate(k1(n))
        allocate(k2(n))
        
        call CFG_get(cfg, "kspace%k_x", k1)
        call CFG_get(cfg, "kspace%k_y", k2)
        call CFG_get(cfg, "kspace%num_points", n_pts)

        allocate(this%k_pts(3, n_sec * (n_pts-1) + 1))
        this%k_pts(3,:) =  0d0

        start = 1
        do i =  1,n_sec
            halt =  start +  n_pts - 1

            this%k_pts(1,start:halt) = linspace(k1(i), k1(i+1), n_pts)
            this%k_pts(2,start:halt) = linspace(k2(i), k2(i+1), n_pts)
            start =  halt
        enddo
        
    end subroutine setup_k_path_abs 

    subroutine setup_k_path_rel(this, cfg)
        implicit none
        class(k_space)        :: this
        type(CFG_t)           :: cfg
        real(8), dimension(:), allocatable :: c1, c2, c1_sec, c2_sec
        real(8), dimension(3)              :: k1, k2
        integer(4) :: n, n_pts, n_sec,i,j, start, halt, cnt


        call CFG_get_size(cfg, "kspace%k_label", n)
        allocate(this%label(n))
        call CFG_get(cfg, "kspace%k_label", this%label)

        call CFG_get_size(cfg, "kspace%k_x", n)
        allocate(c1(n))
        allocate(c2(n))
        n_sec =  n - 1
        call CFG_get(cfg, "kspace%k_x", c1(:))
        call CFG_get(cfg, "kspace%k_y", c2(:))
        call CFG_get(cfg, "kspace%num_points", n_pts)

        allocate(this%k_pts(3,n_sec * (n_pts - 1) + 1))
        allocate(c1_sec(n_sec))
        allocate(c2_sec(n_sec))

        k1(1:2) =  this%ham%UC%rez_lattice(:,1)
        k1(3)   =  0d0
        k2(1:2) =  this%ham%UC%rez_lattice(:,2)
        k2(3)   =  0d0

        start =  1
        do i =  1, n_sec
            halt =  start +  n_pts - 1
            ! linear combination of k-vec in current 
            ! section.
            c1_sec =   linspace(c1(i), c1(i+1), n_pts) 
            c2_sec =   linspace(c2(i), c2(i+1), n_pts)
            
            cnt =  1
            do j =  start,halt
                this%k_pts(:,j) =  c1_sec(cnt) *  k1 +  c2_sec(cnt) *  k2
                cnt =  cnt + 1
            enddo
            start =  halt
        enddo

        deallocate(c1)
        deallocate(c1_sec)
        deallocate(c2)
        deallocate(c2_sec)
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

end module

