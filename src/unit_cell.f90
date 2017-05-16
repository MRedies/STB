module Class_unit_cell
    use Class_atom
    use m_config
    use output
    use m_npy
    use Constants

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
        real(8) :: lattice_constant !> lattice constant in atomic units
        real(8) :: t_nn !> hopping paramater passed for connection 
        real(8) :: eps !> threshold for positional accuracy
        type(atom), dimension(:), allocatable :: atoms !> array containing all atoms
        character(len=25) :: uc_type !> ind
    contains
        procedure :: get_num_atoms       => get_num_atoms
!        procedure :: setup_hexagon       => setup_hexagon
        procedure :: setup_square        => setup_square
        procedure :: setup_single_hex    => setup_single_hex
        procedure :: in_cell             => in_cell
        procedure :: setup_gen_conn      => setup_gen_conn
        procedure :: get_atoms           => get_atoms
        procedure :: gen_find_neigh      => gen_find_neigh
    end type unit_cell
contains
    function angle(a ,b) result(ang)
        implicit none
        real(8), dimension(2), intent(in)   :: a,b 
        real(8)                             :: ang
        ang =  dot_product(a,b)/(norm2(a)* norm2(b))
        ang =  180.0d0 / PI *  acos(ang)
    end function angle
    
    function init_unit(cfg) result(ret)
        implicit none
        type(CFG_t)       :: cfg !> config file as read by m_config
        type(unit_cell)   :: ret
        integer(4), parameter           :: lwork =  20
        real(8)                         :: work(lwork), tmp 
        integer(4), dimension(2)        :: ipiv
        integer(4)                      :: info, i      
        
        call CFG_get(cfg, "grid%epsilon", tmp)
        ret%eps =  tmp * get_unit_conv("length", cfg)
        
        call CFG_get(cfg, "hamil%t_nn", tmp)
        ret%t_nn =  tmp * get_unit_conv("energy", cfg)

        call CFG_get(cfg, "grid%unit_cell_type", ret%uc_type)
        if(trim(ret%uc_type) == "square_2d") then
            call init_unit_square(cfg, ret)
        else
            write (*,*) "Cell type unknown"
            stop
        endif

        ! calculate reciprocal grid
        ret%rez_lattice =  transpose(ret%lattice)
        call dgetrf(2,2, ret%rez_lattice, 2, ipiv, info)
        if(info /= 0) then
            write (*,*) "LU-decomp of lattice vectors failed"
            stop
        endif

        call dgetri(2, ret%rez_lattice, 2, ipiv, work, lwork, info)
        if(info /= 0) then
            write (*,*) "Inversion of lattice vectors failed"
            stop
        endif
        ret%rez_lattice =  2 *  PI * ret%rez_lattice
    end function

    !function init_unit_hex(cfg) result(ret)
        !implicit none
        !type(unit_cell)              :: ret
        !type(CFG_t),  intent(inout)  :: cfg
        !integer(4)                   :: hex_sz,i 
        !real(8)                      :: tmp

        !call CFG_get(cfg, "grid%unit_cell_dim", tmp)
        !ret%unit_cell_dim = tmp * get_unit_conv("length", cfg)

        !call CFG_get(cfg, "hamil%E_s", tmp)
        !ret%E_s =  tmp * get_unit_conv("energy", cfg)

        !call CFG_get(cfg, "hamil%t_nn", tmp)
        !ret%t_nn =  tmp * get_unit_conv("energy", cfg)

        !call CFG_get(cfg, "grid%hexagon_size", ret%atom_per_dim)
        !if(ret%atom_per_dim >=  1) then
            !ret%num_atoms =  calc_num_atoms(ret%atom_per_dim)
            !allocate(ret%atoms(ret%num_atoms))

            !call ret%setup_hexagon()

            !call ret%setup_conn_1D_layer()
            !call ret%setup_lattice_vec()
        !else if(ret%atom_per_dim ==  0) then
            !ret%num_atoms =  1
            !allocate(ret%atoms(1))
            
            !call ret%setup_single_hex()
        !endif
    !end function init_unit_hex

    subroutine init_unit_square(cfg, ret)
        implicit none
        type(unit_cell), intent(inout) :: ret
        type(CFG_t),  intent(inout)    :: cfg
        integer(4)                     :: i
        real(8)                        :: tmp, conn_mtx(2,3), transl_mtx(2,3)
        
        call CFG_get(cfg, "grid%lattice_constant", tmp)
        ret%lattice_constant = tmp * get_unit_conv("length", cfg)

        call CFG_get(cfg, "grid%atoms_per_dim", ret%atom_per_dim)
        ret%num_atoms = ret%atom_per_dim * ret%atom_per_dim
        
        allocate(ret%atoms(ret%num_atoms))
    
        call ret%setup_square()

        conn_mtx(1,:) =  (/ ret%lattice_constant, 0d0, 0d0 /)
        conn_mtx(2,:) =  (/ 0d0, ret%lattice_constant, 0d0 /)
        
        transl_mtx = 0d0
        transl_mtx(1,1:2) = ret%lattice(:,1)
        transl_mtx(2,1:2) = ret%lattice(:,2)

        call ret%setup_gen_conn(conn_mtx, transl_mtx)
    
    end subroutine init_unit_square 

    subroutine setup_single_hex(self)
        implicit none
        class(unit_cell), intent(inout)   :: self
        real(8)                           :: base_len
        real(8), dimension(3,3)           :: base_vecs
        integer(4), parameter             :: lwork =  20
        real(8), dimension(lwork)         :: work 
        integer(4), dimension(2)          :: ipiv
        integer(4)                        :: info, i, j

        self%atoms(1) =  init_ferro((/0d0, 0d0, 0d0/))
        allocate(self%atoms(1)%neigh_idx(3))
        allocate(self%atoms(1)%hopping(3))
        allocate(self%atoms(1)%neigh_conn(3,3))

        self%atoms(1)%hopping   = self%t_nn
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
                self%atoms(cnt) = init_ferro(pos)
                cnt             = cnt + 1
            enddo
        enddo

        self%lattice(:,1) =  (/ 1d0, 0d0 /) * self%atom_per_dim &
                                            *  self%lattice_constant
        self%lattice(:,2) =  (/ 0d0, 1d0 /) * self%atom_per_dim &
                                            *  self%lattice_constant
    end subroutine setup_square

    !Subroutine  setup_hexagon(self)
        !Implicit None
        !class(unit_cell), intent(inout)   :: self
        !real(8), dimension(3)          :: start, pos, halt, dir
        !integer(4)                        :: cnt, row
        !type(atom)                        :: test

        !cnt =  1
        !! sweep from top to (including) middle
        !start =  (/0, self%atom_per_dim /)
        !halt  =  (/self%atom_per_dim , self%atom_per_dim/)
        !dir   =  (/1,0 /)

        !do row = self%atom_per_dim,0,-1
            !pos =  start 
            !do while(any(pos /= halt))
                !self%atoms(cnt) = init_ferro(pos)
                !cnt =  cnt + 1
                !pos =  pos + dir
            !end do

            !start = start +  (/-1, -1/)
            !halt  = halt +  (/0, -1/)
        !end do


        !! sweep after middle downwards
        !start =  (/-self%atom_per_dim, -1/)
        !halt  =  (/self%atom_per_dim -1, -1/)

        !do row =  -1,-(self%atom_per_dim-1), - 1
            !pos =  start
            !do while(any(pos /= halt))
                !self%atoms(cnt) = init_ferro(pos)
                !cnt =  cnt +  1
                !pos =  pos +  dir
            !enddo

            !start =  start +  (/0, -1 /)
            !halt  =  halt  +  (/-1, -1 /)
        !enddo

    !End Subroutine setup_hexagon
    
    subroutine setup_gen_conn(self, conn_mtx, transl_mtx)
        implicit none
        class(unit_cell)    :: self
        real(8), intent(in) :: conn_mtx(:,:) !> Matrix containing
        !> real-space connections. The first index inidcates
        !> the connection vector, the second the vector element
        real(8), intent(in) :: transl_mtx(:,:) !> Matrix containing
        !> real-space translation vectors. Notation as in conn_mtx
        integer(4)                        :: i, j, n_conn
        real(8)  :: start_pos(3)

        n_conn =  size(conn_mtx, dim=1)

        do i =  1, self%num_atoms
            allocate(self%atoms(i)%neigh_idx(n_conn))
            allocate(self%atoms(i)%hopping(n_conn))
            allocate(self%atoms(i)%neigh_conn, mold=conn_mtx)

            self%atoms(i)%hopping =  self%t_nn
            self%atoms(i)%n_neigh =  n_conn
            start_pos             =  self%atoms(i)%pos

            do j =  1,n_conn
                self%atoms(i)%neigh_idx(j) = &
                   self%gen_find_neigh(start_pos, conn_mtx(j,:), transl_mtx)
            enddo
            self%atoms(i)%neigh_conn =  conn_mtx
        enddo
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

        idx =  self%in_cell(start, conn)
        if(idx /= - 1) then
            neigh =  idx
            return
        else 
            do i = 1, n_transl
                idx =  self%in_cell(start, conn + transl_mtx(i,:))
                if (idx /=  - 1) then
                    neigh =  idx
                    return
                endif
                
                idx =  self%in_cell(start, conn - transl_mtx(i,:))
                if (idx /=  - 1) then
                    neigh =  idx
                    return
                endif
            enddo
        endif

        write (*,*) "Couldn't find a generalized neigbour"
        stop
        
    end function gen_find_neigh

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

            delta =  norm2(new -  self%atoms(i)%pos)

            if(delta < self%eps) then
                idx =  i
                exit
            endif
        enddo
    end function in_cell

    function calc_num_atoms(atom_per_dim) result(num_atoms)
        implicit none
        integer(4), intent(in)       :: atom_per_dim
        integer(4)                   :: num_atoms

        if (atom_per_dim > 0) then
            num_atoms =  3 * atom_per_dim *  atom_per_dim  
        else if(atom_per_dim ==  0) then
            num_atoms =  1
        else
            write(*,*) "Invalid atom_per_dim"
            stop 1 
        end if
    end function calc_num_atoms

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

    function get_unit_conv(field_name, cfg) result(factor)
        implicit none
        character(len=*), intent(in)        :: field_name
        type(CFG_t)                         :: cfg
        real(8)                             :: factor
        character(len=300)                  :: unit_name

        call CFG_get(cfg, "units%" // trim(field_name), unit_name)

        select case(trim(unit_name))
        case ("a0")
            factor =  1.0d0
        case ("eV")
            factor =  0.03674932d0
        case ("a0^-1") 
            factor = 1.0d0
        case default
            write (*,*) "Unit unknown"
            stop
        end select
    end function get_unit_conv

end module

