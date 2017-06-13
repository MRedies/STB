module Class_k_space
    use m_config
    use m_npy
    use mpi
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
        integer(4) :: nProcs !> number of MPI Processes
        integer(4) :: me !> MPI ranki
        real(8), allocatable :: k1_param(:) !> 1st k_space param
        real(8), allocatable :: k2_param(:) !> 2nd k_space param
        character(len=300)    :: filling, prefix
        logical      :: perform_dos_integration !> param to toggle dos integr.
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
        procedure :: Bcast_k_space          => Bcast_k_space
        procedure :: plot_omega             => plot_omega
    end type k_space 

contains
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

        call my_section(self%me, self%nProcs, size(self%k_pts, 2), first, last)
        allocate(k_pts_sec(3, last - first + 1))
        k_pts_sec = self%k_pts(:,first:last)

        call self%ham%calc_eigenvalues(k_pts_sec, sec_eig_val)
        
        
        N = 2 *  self%ham%UC%num_atoms 
        allocate(eig_val(N, size(self%k_pts,2)))
        allocate(num_elems(self%nProcs))
        allocate(offsets(self%nProcs))
        call sections(self%nProcs, size(self%k_pts, 2), num_elems, offsets)
        num_elems =  num_elems * N
        offsets   =  offsets   * N
    
        send_count =  N *  size(k_pts_sec, 2)

        call MPI_Gatherv(sec_eig_val, send_count, MPI_REAL8, &
                         eig_val,     num_elems,  offsets,   MPI_REAL8,&
                         root,        MPI_COMM_WORLD, ierr)
        
        if(self%me == root) then 
            call save_npy(trim(self%prefix) //  "band_k.npy", self%k_pts)
            call save_npy(trim(self%prefix) //  "band_E.npy", eig_val)
            call save_npy(trim(self%prefix) //  "lattice.npy", self%ham%UC%lattice)
            call save_npy(trim(self%prefix) //  "rez_lattice.npy", self%ham%UC%rez_lattice)
        endif
       
        deallocate(eig_val)
        deallocate(num_elems)
        deallocate(offsets)
        deallocate(sec_eig_val)
        deallocate(self%k_pts)
    End Subroutine calc_and_print_band

    subroutine calc_pdos(self, E, PDOS)
        implicit none
        class(k_space)          :: self
        real(8), intent(in)     :: E(:)
        real(8), intent(out)    :: PDOS(:,:)
        real(8), allocatable    :: RWORK(:), eig_val(:), loc_PDOS(:,:)
        complex(8), allocatable :: H(:,:), WORK(:)
        integer(4), allocatable :: IWORK(:)
        integer(4)  :: k_idx, E_idx, j, m, N, LWMAX, info, first, last, ierr, num_atoms 

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
        
        num_atoms =  self%ham%UC%num_atoms
        allocate(loc_PDOS(2*num_atoms, self%num_DOS_pts))

        loc_PDOS =  0d0
        PDOS = 0d0

        call my_section(self%me, self%nProcs, size(self%k_pts,2), first, last)
        do k_idx=first, last
            call self%ham%setup_H(self%k_pts(:,k_idx), H)
            call zheevd('V', 'U', N, H, N, eig_val, WORK, LWMAX, &
                RWORK, LWMAX, IWORK, LWMAX, info)
            if( info /= 0) then
                write (*,*) self%me, ": ZHEEVD (with vectors) failed: ", info
                stop
            endif

            do E_idx =  1,self%num_DOS_pts 
                ! eigenvectors are stored column-wise
                ! m-th eigenvalue
                ! j-th component of 
                do m =  1,N
                    do j = 1,N
                        loc_PDOS(j, E_idx) = loc_PDOS(j, E_idx) &
                                           + self%lorentzian(E(E_idx) - eig_val(m)) &
                                           * H(j,m) * conjg(H(j,m))
                    enddo
                enddo
            enddo
        enddo
        loc_PDOS =  loc_PDOS / real(self%DOS_num_k_pts**2)
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
        integer(4)           :: i, num_atoms


        call self%setup_inte_grid_square(self%DOS_num_k_pts)
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

            call save_npy(trim(self%prefix) //  "DOS_E.npy", self%E_DOS)
            call save_npy(trim(self%prefix) //  "DOS.npy", DOS)
            call save_npy(trim(self%prefix) //  "DOS_partial.npy", PDOS)
            call save_npy(trim(self%prefix) //  "DOS_up.npy", up)
            call save_npy(trim(self%prefix) //  "DOS_down.npy", down)

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
            call save_npy(trim(self%prefix) // "DOS_integrated.npy", self%int_DOS)
        endif 
        
        deallocate(self%k_pts)
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
        endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call self%Bcast_k_space()
    end function init_k_space

    subroutine Bcast_k_space(self)
        class(k_space)         :: self
        integer(4), parameter  :: num_cast =  12
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
                       root,                 MPI_COMM_WORLD, ierr(11))
        call MPI_Bcast(self%temp,            1,              MPI_REAL8,    &
                       root,                 MPI_COMM_WORLD, ierr(12))
       
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
        
        call linspace(0d0, 1d0, n_k+1, ls_help)
        allocate(ls(n_k))

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

    subroutine setup_k_path_abs(self)
        implicit none
        class(k_space)        :: self
        integer(4)            :: n_pts, n_sec, start, halt, i
        real(8), allocatable  :: tmp(:)



        self%k1_param =  self%k1_param * self%units%inv_length
        self%k2_param =  self%k2_param * self%units%inv_length

        n_pts =  self%num_k_pts
        n_sec =  size(self%k1_param)-1


        allocate(self%k_pts(3, n_sec * (self%num_k_pts-1) + 1))
        self%k_pts(3,:) =  0d0

        start = 1
        do i =  1,n_sec
            halt =  start +  n_pts - 1

            call linspace(self%k1_param(i), self%k1_param(i+1), n_pts, tmp)
            self%k_pts(1,start:halt) = tmp
            call linspace(self%k2_param(i), self%k2_param(i+1), n_pts, tmp)
            self%k_pts(2,start:halt) = tmp
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

            call linspace(self%k1_param(i),self% k1_param(i+1), n_pts, c1_sec)
            call linspace(self%k2_param(i), self%k2_param(i+1), n_pts, c2_sec)

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

        vol = my_norm2(cross_prod(k1,k2))
    end function vol_k_space_parallelo

    subroutine calc_hall_conductance(self, ret)
        implicit none
        class(k_space)       :: self
        real(8)              :: hall, V_k, k(3), ret
        real(8), allocatable :: eig_val(:), omega_z(:)
        integer(4)           :: N_k, n, k_idx, first, last, ierr, n_atm

        if(allocated(self%k_pts) )then
            deallocate(self%k_pts)
        endif

        call self%setup_inte_grid_square(self%berry_num_k_pts)
        V_k = self%vol_k_space_parallelo()
        N_k = size(self%k_pts, 2)

        call my_section(self%me, self%nProcs, N_k, first, last)
        
        hall = 0d0

        n_atm =  self%ham%UC%num_atoms
        allocate(self%ham%del_H(2*n_atm,2*n_atm))
        do k_idx = first, last
            k = self%k_pts(:,k_idx)
            call self%ham%calc_berry_z(k, omega_z)
            call self%ham%calc_single_eigenvalue(k, eig_val)
            
            do n = 1,2*self%ham%UC%num_atoms
                hall = hall + omega_z(n) * self%fermi_distr(eig_val(n))
            enddo
        enddo
        deallocate(self%ham%del_H)

        hall = hall * V_k/real(N_k)
        hall = hall / (2d0*PI)
        call MPI_Reduce(hall, ret, 1, &
                        MPI_REAL8, MPI_Sum, &
                        root, MPI_COMM_WORLD, ierr)

        if(self%me == root) then
            call save_npy(trim(self%prefix) // "hall_cond.npy", (/ret /))
        endif
    end subroutine calc_hall_conductance

    subroutine plot_omega(self)
        implicit none
        class(k_space)         :: self
        real(8), allocatable   :: omega_z(:,:), tmp_vec(:), sec_omega_z(:,:)
        real(8)                :: k(3)
        integer(4)             :: N, k_idx, send_count, first, last, ierr, cnt
        integer(4), allocatable:: num_elems(:), offsets(:)
        integer(4), parameter  :: dim_sz = 100
        
        allocate(num_elems(self%nProcs))
        allocate(offsets(self%nProcs))
        N = 2* self%ham%UC%num_atoms
        call self%setup_inte_grid_square(dim_sz)
        allocate(omega_z(N, size(self%k_pts, 2)))
        allocate(tmp_vec(N))
        
        call sections(self%nProcs, size(self%k_pts, 2), num_elems, offsets)
        call my_section(self%me, self%nProcs, size(self%k_pts, 2), first, last)
        num_elems =  num_elems * N
        offsets   =  offsets   * N
        send_count =  N *  (last - first + 1)
        allocate(sec_omega_z(N, send_count))

        !do i =  0,self%nProcs-1
            !if(self%me ==  i) &
                !write (*,*) self%me, first, last, send_count/N, offsets(i+1)/N
            !call MPI_Barrier(MPI_COMM_WORLD, ierr)
        !enddo
        
        cnt =  1
        !write (*,*) self%me, "Pre calc"
        do k_idx = first,last
            k = self%k_pts(:,k_idx)
            
            !omega_z
            !call self%ham%calc_berry_z(k, tmp_vec)

            !omega_xy
            call self%ham%calc_berry_tensor_elem(1,2, k, tmp_vec)
            
            sec_omega_z(:,cnt) =  tmp_vec 
            cnt = cnt + 1
        enddo
        !write (*,*) self%me, "Post calc"

        call MPI_Gatherv(sec_omega_z, send_count, MPI_REAL8, &
                         omega_z,    num_elems, offsets, MPI_REAL8, &
                         root, MPI_COMM_WORLD, ierr)
            
        if(self%me ==  root) then
            call save_npy(trim(self%prefix) //  "omega_xy_z.npy", omega_z)
            call save_npy(trim(self%prefix) //  "omega_xy_k.npy", self%k_pts)
        endif
        cnt =  1
        do k_idx = first,last
            k = self%k_pts(:,k_idx)
            !omega_xy
            call self%ham%calc_berry_tensor_elem(2,1, k, tmp_vec)
            
            sec_omega_z(:,cnt) =  tmp_vec 
            cnt = cnt + 1
        enddo
        !write (*,*) self%me, "Post calc"

        call MPI_Gatherv(sec_omega_z, send_count, MPI_REAL8, &
                         omega_z,    num_elems, offsets, MPI_REAL8, &
                         root, MPI_COMM_WORLD, ierr)
            
        if(self%me ==  root) then
            call save_npy(trim(self%prefix) //  "omega_yx_z.npy", omega_z)
            call save_npy(trim(self%prefix) //  "omega_yx_k.npy", self%k_pts)
        endif
    end subroutine plot_omega 

    subroutine set_fermi(self, cfg)
        implicit none
        class(k_space)         :: self
        class(CFG_t)           :: cfg
        real(8)                :: tmp
        integer(4)             :: ierr

        if(root == self%me) then
            call CFG_get(cfg, "dos%E_fermi", tmp)
        endif
        call MPI_Bcast(tmp, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr)

        self%E_fermi =  tmp * self%units%energy
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

    function fermi_distr(self, E) result(ferm)
        implicit none
        class(k_space), intent(in)    :: self
        real(8), intent(in)           :: E
        real(8)                       :: ferm

        ferm = 1d0 / (exp((E-self%E_fermi)&
            /(boltzmann_const * self%temp)) + 1d0)
    end function fermi_distr

end module

