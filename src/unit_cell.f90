module Class_unit_cell
    use Class_atom
    use Class_helper
    use m_config
    use output
    use m_npy
    use mpi
    use Constants
    use class_Units

    implicit none

    type unit_cell
        real(8), public :: lattice(2,2) !> translation vectors
        !> of the real-space lattice. First index: element of vector
        !> Second index: 1 or second vector
        real(8), public :: rez_lattice(2,2) !> translation vectors
        !> of the reciprocal lattice. Indexs same as lattice
        ! number of non-redundant atoms pre unit cell
        integer(4) :: num_atoms  !> number of non-redundant atmos in a unit cell
        integer(4) :: atom_per_dim !> atoms along the radius of the unit_cell
        integer(4) :: nProcs
        integer(4) :: me
        real(8) :: lattice_constant !> lattice constant in atomic units
        real(8) :: t_nn !> hopping paramater passed for connection 
        real(8) :: eps !> threshold for positional accuracy
        real(8) :: ferro_phi, ferro_theta
        type(atom), dimension(:), allocatable :: atoms !> array containing all atoms
        type(units)       :: units
        character(len=25) :: uc_type !> indicates shape of unitcell
        character(len=25) :: mag_type !> indicates type of magnetization
    contains

        procedure :: get_num_atoms               => get_num_atoms
        procedure :: setup_square                => setup_square
        procedure :: setup_single_hex            => setup_single_hex
        procedure :: in_cell                     => in_cell
        procedure :: setup_gen_conn              => setup_gen_conn
        procedure :: get_atoms                   => get_atoms
        procedure :: gen_find_neigh              => gen_find_neigh
        procedure :: save_unit_cell              => save_unit_cell
        procedure :: set_mag_ferro               => set_mag_ferro
        procedure :: set_mag_random              => set_mag_random
        procedure :: set_mag_x_spiral_square     => set_mag_x_spiral_square
        procedure :: set_mag_linrot_skrym_square => set_mag_linrot_skrym_square
        procedure :: Bcast_UC                    => Bcast_UC
        procedure :: setup_honey                 => setup_honey
    end type unit_cell
contains
    function angle(a ,b) result(ang)
        implicit none
        real(8), dimension(2), intent(in)   :: a,b 
        real(8)                             :: ang
        ang                                       = dot_product(a,b)/(my_norm2(a)* my_norm2(b))
        ang                                       = 180.0d0 / PI *  acos(ang)
    end function angle
    
    function init_unit(cfg) result(self)
        implicit none
        type(CFG_t)       :: cfg !> config file as read by m_config
        type(unit_cell)   :: self
        integer(4), parameter           :: lwork =  20
        real(8)                         :: work(lwork), tmp 
        integer(4), dimension(2)        :: ipiv
        integer(4)                      :: info, ierr
        
        call MPI_Comm_size(MPI_COMM_WORLD, self%nProcs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)
        
        self%units = init_units(cfg, self%me)
        
        
        if(self%me ==  0) then 
            call CFG_get(cfg, "grid%epsilon", tmp)
            self%eps =  tmp * self%units%length
            
            call CFG_get(cfg, "hamil%t_nn", tmp)
            self%t_nn =  tmp * self%units%energy

            call CFG_get(cfg, "grid%mag_type", self%mag_type)
            call CFG_get(cfg, "grid%unit_cell_type", self%uc_type)

            call CFG_get(cfg, "grid%lattice_constant", tmp)
            self%lattice_constant = tmp * self%units%length

            call CFG_get(cfg, "grid%atoms_per_dim", self%atom_per_dim)

            call CFG_get(cfg, "grid%ferro_phi",   self%ferro_phi)
            call CFG_get(cfg, "grid%ferro_theta", self%ferro_theta)
        endif
        call self%Bcast_UC()



        if(trim(self%uc_type) == "square_2d") then
            call init_unit_square(self)
        else if(trim(self%uc_type) == "honey_2d") then
            call init_unit_honey(self)
        else
            write (*,*) self%me, ": Cell type unknown"
            stop
        endif


        ! calculate reciprocal grid
        self%rez_lattice =  transpose(self%lattice)
        call dgetrf(2,2, self%rez_lattice, 2, ipiv, info)
        if(info /= 0) then
            write (*,*) self%me, ": LU-decomp of lattice vectors failed"
            stop
        endif

        call dgetri(2, self%rez_lattice, 2, ipiv, work, lwork, info)
        if(info /= 0) then
            write (*,*) self%me, ": Inversion of lattice vectors failed"
            stop
        endif
        self%rez_lattice =  2 *  PI * self%rez_lattice


    end function

    subroutine Bcast_UC(self)
        implicit none
        class(unit_cell)              :: self
        integer(4), parameter         :: num_cast = 8
        integer(4)                    :: ierr(num_cast)
        
        call MPI_Bcast(self%eps,              1,              MPI_REAL8,     &
                       root,                  MPI_COMM_WORLD, ierr(1))
        call MPI_Bcast(self%t_nn,             1,              MPI_REAL8,     &
                       root,                  MPI_COMM_WORLD, ierr(2))
        call MPI_Bcast(self%mag_type,         25,             MPI_CHARACTER, &
                       root,                  MPI_COMM_WORLD, ierr(3))
        call MPI_Bcast(self%uc_type,          25,             MPI_CHARACTER, &
                       root,                  MPI_COMM_WORLD, ierr(4))
        call MPI_Bcast(self%lattice_constant, 1,              MPI_REAL8, &
                       root,                  MPI_COMM_WORLD, ierr(5))
        call MPI_Bcast(self%atom_per_dim,     1,              MPI_INTEGER4, &
                       root,                  MPI_COMM_WORLD, ierr(6))
        
        call MPI_Bcast(self%ferro_phi,   1,              MPI_REAL8, &
                       root,             MPI_COMM_WORLD, ierr(7))
        call MPI_Bcast(self%ferro_theta, 1,              MPI_REAL8, &
                       root,             MPI_COMM_WORLD, ierr(8))
        
        call check_ierr(ierr, self%me, "Unit cell check err")
    end subroutine Bcast_UC

    subroutine init_unit_square(ret)
        implicit none
        type(unit_cell), intent(inout) :: ret
        real(8)                        :: conn_mtx(2,3), transl_mtx(2,3)
        
        ret%num_atoms = ret%atom_per_dim * ret%atom_per_dim
        allocate(ret%atoms(ret%num_atoms))
    
        call ret%setup_square()

        conn_mtx(1,:) =  (/ ret%lattice_constant, 0d0, 0d0 /)
        conn_mtx(2,:) =  (/ 0d0, ret%lattice_constant, 0d0 /)
        
        transl_mtx = 0d0
        transl_mtx(1,1:2) = ret%lattice(:,1)
        transl_mtx(2,1:2) = ret%lattice(:,2)

        call ret%setup_gen_conn(conn_mtx, transl_mtx)    

        if(trim(ret%mag_type) ==  "x_spiral") then
            call ret%set_mag_x_spiral_square()
        else if(trim(ret%mag_type) == "ferro") then
            call ret%set_mag_ferro()
        else if(trim(ret%mag_type) == "lin_skyrm") then
            call ret%set_mag_linrot_skrym_square()
        else if(trim(ret%mag_type) == "random") then
            call ret%set_mag_random()
        else
            write (*,*) "Mag_type not known"
            stop
        endif
    end subroutine init_unit_square

    subroutine init_unit_honey(ret)
        implicit none
        type(unit_cell), intent(inout)   :: ret
        real(8)  :: transl_mtx(3,3), l, base_len_uc, pos(3), conn_mtx(3,3)
        real(8), allocatable             :: grid(:,:), hexagon(:,:)
        real(8), parameter               :: deg_30 =  30.0 * PI / 180.0
        real(8), parameter               :: deg_60 =  60.0 * PI / 180.0
        integer(4)                       :: apd, cnt, i
        

        apd         = ret%atom_per_dim
        base_len_uc = ret%lattice_constant * apd
        l           = 2 *  cos(deg_30) * base_len_uc
    
        transl_mtx(1, :) =  l *  [1d0,   0d0,           0d0]
        transl_mtx(2, :) =  l *  [0.5d0, sin(deg_60),   0d0]
        transl_mtx(3, :) =  l *  [0.5d0, - sin(deg_60), 0d0]

        ret%lattice(:,1) =  transl_mtx(1,1:2)
        ret%lattice(:,2) =  transl_mtx(2,1:2)

        ret%num_atoms = calc_num_atoms_non_red_honey(apd)
        allocate(ret%atoms(ret%num_atoms))
        allocate(hexagon(ret%num_atoms, 3))
        hexagon = 0d0
        call gen_honey_grid(ret%lattice_constant, apd, grid)

        cnt =  1
        do i =  1,size(grid,1)
            pos =  grid(i,:)
            if(in_hexagon(pos, base_len_uc)) then
                if(.not. already_in_red(pos, hexagon, cnt-1, transl_mtx)) then
                    hexagon(cnt,:) =  grid(i,:)
                    cnt =  cnt + 1
                endif
            endif
        enddo
        call ret%setup_honey(hexagon)
        
        if(trim(ret%mag_type) == "ferro") then
            call ret%set_mag_ferro()
        else if(trim(ret%mag_type) == "random") then
            call ret%set_mag_random()
        else
            write (*,*) "Mag_type not known"
            stop
        endif

        ! only one kind of atom from honey-comb unit cell needed
        ! the other comes through complex conjugate
        conn_mtx(1, :) =  ret%lattice_constant * [0d0,          1d0,           0d0]
        conn_mtx(2, :) =  ret%lattice_constant * [cos(deg_30),  - sin(deg_30), 0d0]
        conn_mtx(3, :) =  ret%lattice_constant * [-cos(deg_30), - sin(deg_30), 0d0]
        
        
        call ret%setup_gen_conn(conn_mtx, transl_mtx)  
    end subroutine init_unit_honey

    subroutine set_mag_ferro(self)
        implicit none
        class(unit_cell)    :: self
        integer(4)          :: i

        do i = 1,self%num_atoms
            call self%atoms(i)%set_sphere(self%ferro_phi, self%ferro_theta)
        enddo
    end subroutine set_mag_ferro


    subroutine set_mag_x_spiral_square(self)
        implicit none
        class(unit_cell)                 :: self 
        real(8)        :: alpha, rel_xpos
        integer(4)     :: i

        do i =  1,self%num_atoms
            rel_xpos =  self%atoms(i)%pos(1) / self%lattice_constant
            alpha =  rel_xpos *  2*PI / self%atom_per_dim
            if(alpha <= PI) then
                self%atoms(i)%m_theta = alpha
                self%atoms(i)%m_phi   = 0d0
            else
                self%atoms(i)%m_theta = 2*PI - alpha 
                self%atoms(i)%m_phi   = PI
            endif
        enddo

    end subroutine set_mag_x_spiral_square

    subroutine set_mag_random(self)
        implicit none
        class(unit_cell)       :: self
        integer(4)             :: i
        real(8)                :: phi, theta, r(2)

        do i =  1,self%num_atoms
            call random_number(r)

            phi   = r(1) * 2d0 * PI
            theta = r(2) * PI
            call self%atoms(i)%set_sphere(phi, theta)
        enddo
    end subroutine set_mag_random

    subroutine set_mag_linrot_skrym_square(self)
        implicit none
        class(unit_cell)     :: self
        real(8), parameter   :: e_z(3) = (/0,0,1/)
        real(8) :: R(3,3), center(3), conn(3), n(3),m(3), radius, alpha 
        integer(4)           :: i

        ! Nagaosa style unit cell. Ferromagnetic border only to the left
        radius = 0.5d0 * self%lattice_constant * self%atom_per_dim
        center = (/radius, radius, 0d0/)
        alpha =  0d0 
        do i =  1,self%num_atoms
            conn  = center - self%atoms(i)%pos
            if(my_norm2(conn) > 1d-6 * self%lattice_constant &
                    .and. my_norm2(conn) <= radius + self%eps) then 
                n     = cross_prod(conn, e_z)
                alpha =  PI * (1d0 -  my_norm2(conn) / radius)
                R     = R_mtx(alpha, n)
                ! center of skyrmion point down
                m     = matmul(R,  e_z)
            else if(my_norm2(conn) <= 1d-6 * self%lattice_constant) then 
                m =  - e_z
            else
                m =  e_z
            endif


            call self%atoms(i)%set_m_cart(m(1), m(2), m(3))
        enddo

    end subroutine set_mag_linrot_skrym_square 

    subroutine save_unit_cell(self, folder)
        implicit none
        class(unit_cell)        :: self
        character(len=*)        :: folder
        real(8), allocatable    :: x(:), y(:), z(:), phi(:), theta(:)
        integer(4)              :: i
        
        allocate(x(self%num_atoms))
        allocate(y(self%num_atoms))
        allocate(z(self%num_atoms))
        allocate(phi(self%num_atoms))
        allocate(theta(self%num_atoms))

        do i =  1,self%num_atoms
            x(i)     = self%atoms(i)%pos(1)
            y(i)     = self%atoms(i)%pos(2)
            z(i)     = self%atoms(i)%pos(3)
            phi(i)   = self%atoms(i)%m_phi
            theta(i) = self%atoms(i)%m_theta

            !if(self%atoms(i)%n_neigh > 0) then
                !write(conn_name, '(a,i4.4,a,a)') 'conn_', i, '.npy'
                !call save_npy(folder // trim(conn_name), self%atoms(i)%neigh_conn)
                
                !write(conn_name, '(a,i4.4,a,a)') 'connidx_', i, '.npy'
                !call save_npy(folder // trim(conn_name), self%atoms(i)%neigh_idx)
            !endif

        enddo

        call save_npy(folder // "pos_x.npy", x / self%units%length)
        call save_npy(folder // "pos_y.npy", y / self%units%length)
        call save_npy(folder // "pos_z.npy", z / self%units%length)
        call save_npy(folder // "m_phi.npy", phi)
        call save_npy(folder // "m_theta.npy", theta)

    end subroutine save_unit_cell

    subroutine setup_single_hex(self)
        implicit none
        class(unit_cell), intent(inout)   :: self
        real(8)                           :: base_len
        real(8), dimension(3,3)           :: base_vecs

        self%atoms(1) =  init_ferro_z((/0d0, 0d0, 0d0/))
        allocate(self%atoms(1)%neigh_idx(3))
        allocate(self%atoms(1)%neigh_conn(3,3))

        self%atoms(1)%n_neigh   = 3
        self%atoms(1)%neigh_idx = (/ 1,1,1 /)

        base_len = self%lattice_constant
        base_vecs(1, :) = (/ 1d0,   0d0,                  0d0 /)
        base_vecs(2, :) = (/ 0.5d0, sin(60d0/180d0 * PI), 0d0 /)
        base_vecs(3, :) = (/-0.5d0, sin(60d0/180d0 * PI), 0d0 /)
        base_vecs =  base_vecs *  base_len
        self%atoms(1)%neigh_conn =  base_vecs 
    
        self%lattice(:,1) = base_vecs(1,1:2)
        self%lattice(:,2) = base_vecs(2,1:2)
    end subroutine

    subroutine setup_square(self)
        implicit none
        class(unit_cell), intent(inout)  :: self
        integer(4)                       :: i, j, cnt
        real(8) :: pos(3)

        cnt =  1
        do i = 0, self%atom_per_dim-1
            do j = 0, self%atom_per_dim-1
                pos             = (/i,j,0/) * self%lattice_constant 
                self%atoms(cnt) = init_ferro_z(pos)
                cnt             = cnt + 1
            enddo
        enddo

        self%lattice(:,1) =  (/ 1d0, 0d0 /) * self%atom_per_dim &
                                            *  self%lattice_constant
        self%lattice(:,2) =  (/ 0d0, 1d0 /) * self%atom_per_dim &
                                            *  self%lattice_constant
    end subroutine setup_square

    subroutine setup_honey(self, hexagon)
        implicit none
        class(unit_cell), intent(inout)  :: self
        real(8), intent(in)              :: hexagon(:,:)
        real(8)                          :: pos(3)
        integer(4)                       :: i

        do i =  1, size(hexagon, dim=1)
            pos           =  hexagon(i,:)
            self%atoms(i) =  init_ferro_z(pos)
        enddo
    end subroutine setup_honey

    subroutine setup_gen_conn(self, conn_mtx, transl_mtx)
        implicit none
        class(unit_cell)    :: self
        real(8), intent(in) :: conn_mtx(:,:) !> Matrix containing
        !> real-space connections. The first index inidcates
        !> the connection vector, the second the vector element
        real(8), intent(in) :: transl_mtx(:,:) !> Matrix containing
        !> real-space translation vectors. Notation as in conn_mtx
        integer(4)                        :: i, j, cnt, candidate, n_conn, n_found
        integer(4), allocatable :: neigh_cand(:)
        real(8)  :: start_pos(3), conn(3)
        logical, allocatable :: found_conn(:)

        n_conn =  size(conn_mtx, 1)
        allocate(found_conn(n_conn))
        allocate(neigh_cand(n_conn))

        do i =  1, self%num_atoms
            start_pos             =  self%atoms(i)%pos

            n_found    = 0
            found_conn = .False.
            cnt        = 1
            neigh_cand =  - 1

            do j =  1,n_conn
                conn =  conn_mtx(j,:)
                candidate = self%gen_find_neigh(start_pos, conn, transl_mtx)

                if(candidate /= - 1) then
                    !self%atoms(i)%neigh_idx(cnt) = candidate
                    found_conn(j) =  .True.
                    neigh_cand(j) =  candidate
                    n_found = n_found + 1
                    cnt     = cnt + 1 
                endif
            enddo



            allocate(self%atoms(i)%neigh_idx(n_found))
            allocate(self%atoms(i)%neigh_conn(n_found, 3))
            self%atoms(i)%n_neigh =  n_found
            
            cnt =  1
            do j = 1,n_conn
                if(found_conn(j)) then
                    self%atoms(i)%neigh_conn(cnt,:) =  conn_mtx(j,:)
                    self%atoms(i)%neigh_idx(cnt)    = neigh_cand(j)
                    cnt =  cnt + 1
                endif
            enddo
        enddo
        deallocate(found_conn)
        deallocate(neigh_cand)
    end subroutine setup_gen_conn

    function gen_find_neigh(self, start, conn, transl_mtx) result(neigh)
        implicit none
        class(unit_cell), intent(in)         :: self
        real(8), intent(in) :: start(3) !> starting position in RS
        real(8), intent(in) :: conn(3) !> RS connection
        real(8), intent(in) :: transl_mtx(:,:) !> RS translation vectors to next unit cell
        !> The vectors are save as columns in the matrix:
        !> The first index indicates the vector
        !> The second index indicates the element of the vector
        integer(4)  :: neigh, idx, n_transl, i
        real(8) :: new(3)
        
        n_transl =  size(transl_mtx, dim=1)

        neigh =  - 1
        idx =  self%in_cell(start, conn)
        if(idx /= - 1) then
            neigh =  idx
            return
        else 
            do i = 1, n_transl
                new =  conn + transl_mtx(i,:)
                idx =  self%in_cell(start, new)
                if (idx /=  - 1) then
                    neigh =  idx
                    return
                endif
                
                new =  conn -  transl_mtx(i,:)
                idx =  self%in_cell(start, new)
                if (idx /=  - 1) then
                    neigh =  idx
                    return
                endif
            enddo
        endif
    end function gen_find_neigh

    function already_in_red(pos, hex, till, transl_mtx) result(inside)
        implicit none
        real(8), intent(in)    :: pos(3), hex(:,:), transl_mtx(:,:)
        integer(4), intent(in) :: till
        logical                :: inside
        real(8)                :: new(3), delta
        integer(4)             :: n_transl, i, trl

        n_transl = size(transl_mtx, dim = 1)
        inside   = .False.
      

        outer: do i =  1, till
            do trl =  1, n_transl
                new   =  pos + transl_mtx(trl,:)
                delta =  my_norm2(hex(i,:) -  new)
                
                if(delta <= pos_eps) then
                    inside = .True.
                    exit outer
                endif
                new   =  pos - transl_mtx(trl,:)
                delta =  my_norm2(hex(i,:) -  new)
                
                if(delta <= pos_eps) then
                    inside = .True.
                    exit outer
                endif
            enddo
        enddo outer
    end function already_in_red 

    function in_cell(self, start, conn) result(idx)
        ! if position is in hexagon the corresponding index is
        ! returned, else - 1
        implicit none
        class(unit_cell), intent(in)          :: self
        real(8), intent(in) :: start(3) !> RS start position
        real(8), intent(in) :: conn(3) !> RZ connection
        real(8) :: new(3), delta
        integer(4) :: idx 
        integer(4) :: i

        new =  start +  conn

        idx =  -1
        do i =  1, self%num_atoms

            delta =  my_norm2(new -  self%atoms(i)%pos)

            if(delta < self%eps) then
                idx =  i
                exit
            endif
        enddo
    end function in_cell
    
    subroutine calc_num_atoms_full_honey(n, n_atm, side)
        implicit none
        integer(4), intent(in)  :: n
        integer(4), intent(out) :: n_atm, side
        integer(4)              :: i

        side =  2
        n_atm =  0

        do i = 1,n
            if(mod(i,3) == 0) then
                n_atm =  n_atm + side
                side  =  side + 2
            else
                n_atm =  n_atm + side - 1
            endif
        enddo
        n_atm =  n_atm * 6
    end subroutine calc_num_atoms_full_honey

    function calc_num_atoms_non_red_honey(n) result(n_atm)
        implicit none
        integer(4), intent(in)   :: n
        integer(4)               :: inner, next_side, n_atm

        call calc_num_atoms_full_honey(n-1, inner, next_side)

        if(mod(n,3) == 0) then
            n_atm =  inner +  3 * next_side
        else
            n_atm =  inner +  3 * next_side - 4
        endif
    end function calc_num_atoms_non_red_honey

    function in_hexagon(pos, a) result(inside)
        implicit none
        real(8), intent(in)         :: pos(3), a
        real(8)                     :: m, b
        real(8), parameter          :: deg_30 =  30d0/180d0 * PI
        logical                     :: inside

        m =  - tan(deg_30)
        inside =  .True.

        if(abs(pos(1)) > a * (cos(deg_30) + pos_eps)) then
            inside = .False.
        endif

        b =  a * (1d0 +  pos_eps)
        if(.not. (abs(pos(2)) <= m * abs(pos(1)) + b )) then
            inside = .False.
        endif
    end function in_hexagon

    subroutine gen_hexa_grid(l, origin, max_ind, grid)
        implicit none
        real(8), intent(in)     :: l, origin(3)
        integer(4), intent(in)  :: max_ind
        real(8), allocatable    :: grid(:,:)
        real(8)                 :: v1(3), v2(3)
        real(8), parameter      :: deg_30 =  30d0/180d0 * PI
        integer(4)              :: cnt, i, j

        if(.not. allocated(grid)) then
            allocate(grid((2*max_ind + 1)**2,3))
        endif

        v1 =  [l,       0d0,             0d0]
        v2 =  [0.5d0*l, cos(deg_30) * l, 0d0]

        cnt = 1
        do i = -max_ind,max_ind
            do j= - max_ind,max_ind
                grid(cnt,:) = origin + i * v1 + j * v2
                cnt =  cnt + 1
            enddo
        enddo
    end subroutine gen_hexa_grid

    subroutine gen_honey_grid(a, max_ind, grid)
        implicit none
        real(8), intent(in)        :: a
        integer(4), intent(in)     :: max_ind
        real(8), allocatable       :: grid(:,:), tmp(:,:)
        real(8)                    :: l, origin(3)
        real(8), parameter         :: deg_30 =  30d0/180d0 * PI
        integer(4)                 :: n

        n = 2 * max_ind + 1
        l =  2d0 * cos(deg_30) * a

        if(.not. allocated(grid)) then
            allocate(grid(2 * n**2, 3))
        endif

        origin =  [0d0, a, 0d0]
        call gen_hexa_grid(l, origin,  max_ind, tmp)
        grid(1:n**2,:) = tmp
        
        origin =   [0d0,-a, 0d0]
        call gen_hexa_grid(l, origin, max_ind, tmp)
        grid(n**2 + 1:2 * n**2,:) = tmp
    end subroutine gen_honey_grid

    function get_num_atoms(self) result(num)
        implicit none
        class(unit_cell), intent(in) :: self
        integer(4) :: num
        num = self%num_atoms 
    end function get_num_atoms

    function get_atoms(self) result(ret)
        implicit none
        class(unit_cell), intent(in)            :: self
        type(atom), dimension(:), allocatable   :: ret

        ret =  self%atoms        
    end function get_atoms 

    function rot_z_deg(deg) result(rot)
        implicit none 
        real(8), intent(in)        :: deg
        real(8), dimension(3,3)    :: rot
        real(8)                    :: bog

        bog =  deg * PI /  180.0d0

        rot      = 0.0d0
        rot(1,1) = cos(bog)
        rot(1,2) = sin(bog)
        rot(2,1) = - sin(bog)
        rot(2,2) = cos(bog)
        rot(3,3) = 1.0d0
    end function rot_z_deg

    function R_mtx(theta, vec) result(R)
        implicit none
        real(8), intent(in)    :: theta!> rotation angle
        real(8), intent(in)    :: vec(3) !> vector to rotate AROUND
        real(8), parameter     :: Iden(3,3) &
                               = reshape((/1,0,0,& !this works only for
                                           0,1,0,& !symm matrcies
                                           0,0,1/), (/3,3/)) !fort-order...
        real(8)  :: R(3,3), u(3,1), u_x_u(3,3), u_x(3,3)
        
        u(:,1) =  vec / my_norm2(vec)
        
        u_x_u =  matmul(u, transpose(u))

        u_x =  0d0
        u_x(2,1) =  u(3,1)
        u_x(3,1) = -u(2,1)
        u_x(1,2) = -u(3,1)
        u_x(3,2) =  u(1,1)
        u_x(1,3) =  u(2,1)
        u_x(2,3) = -u(1,1)

        R =  cos(theta)    * Iden &
          +  sin(theta)    * u_x  &
          + (1-cos(theta)) * u_x_u

    end function R_mtx
end module

