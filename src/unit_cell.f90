module Class_unit_cell
    use Class_atom
    use m_config
    use output
    implicit none
    real(8), parameter     :: PI     = 3.14159265359d0
    complex(8), parameter :: i_unit = cmplx(0d0, 1d0)

    type unit_cell
        real(8), public, dimension(2,2) :: lattice, rez_lattice
        ! number of non-redundant atoms pre unit cell
        integer(4), private      :: num_atoms, hex_size
        real(8), private         :: unit_cell_dim, E_s, in_plane_hopping
        type(atom), private, dimension(:), allocatable   :: atoms 
    contains
        procedure :: get_num_atoms       => get_num_atoms
        procedure :: setup_hexagon       => setup_hexagon
        procedure :: find_neigh          => find_neigh
        procedure :: in_hexagon          => in_hexagon
        procedure :: setup_conn_1D_layer => setup_conn_1D_layer
        procedure :: get_atoms           => get_atoms
            procedure :: get_ham             => get_ham
        procedure :: calc_eigenvalues    =>  calc_eigenvalues
        procedure :: setup_lattice_vec   =>  setup_lattice_vec
    end type unit_cell
contains
    function angle(a ,b) result(ang)
        implicit none
        real(8), dimension(2), intent(in)   :: a,b 
        real(8)                             :: ang
        ang =  dot_product(a,b)/(norm2(a)* norm2(b))
        ang =  180.0d0 / PI *  acos(ang)
    end function angle
    
    subroutine setup_lattice_vec(this)
        implicit none
        class(unit_cell), intent(inout) :: this 
        real(8)                         :: ucd
        integer(4), parameter           :: lwork =  20
        real(8), dimension(lwork)       :: work 
        integer(4), dimension(2)        :: ipiv
        integer(4)                      :: info, i, j
        
        ucd =  this%unit_cell_dim 
        this%lattice(:,1) =  ucd * (/1.5d0,  sin(60d0 * PI/180d0) /)
        this%lattice(:,2) =  ucd * (/1.5d0, -sin(60d0 * PI/180d0) /)

        ! preform inversion 
        this%rez_lattice =  transpose(this%lattice)
        call dgetrf(2,2, this%rez_lattice, 2, ipiv, info)
        if(info /= 0) then
            write (*,*) "LU-decomp of lattice vectors failed"
            stop
        endif

        call dgetri(2, this%rez_lattice, 2, ipiv, work, lwork, info)
        if(info /= 0) then
            write (*,*) "Inversion of lattice vectors failed"
            stop
        endif
        this%rez_lattice =  2 *  PI * this%rez_lattice
        

        do i =  1,2
            write (*,*) "Real space", i
            call print_mtx(this%lattice(:,i))
        enddo
        write (*,*)  "#############"
        do i =  1,2
            write (*,*) "Rez: ", i
            call print_mtx(this%rez_lattice(:,i))
        enddo

        
    end subroutine setup_lattice_vec
    
    function calc_eigenvalues(this, k_list) result(eig_val)
        implicit none
        class(unit_cell)                :: this
        real(8), dimension(:,:), intent(in)         :: k_list
        real(8), dimension(:,:), allocatable        :: eig_val
        real(8), dimension(3)                       :: k
        complex(8), dimension(:,:), allocatable     :: H
        integer(4) :: i, N, LWMAX, info
        real(8), dimension(:), allocatable          :: RWORK
        complex(8), dimension(:), allocatable       :: WORK 

        N =  2 * this%num_atoms
        LWMAX =  10*N
        allocate(eig_val(size(k_list, 2), N))
        allocate(H(N,N))
        allocate(RWORK(LWMAX))
        allocate(WORK(LWMAX))
        
        do i = 1,size(k_list,2)
            k =  k_list(:,i)
            call this%get_ham(k, H)
            
            call zheev('V', 'U', N, H, N, eig_val(i,:), WORK, LWMAX, RWORK, info)
            if( info /= 0) then
                write (*,*) "ZHEEV failed: ", info
                stop
            endif

        enddo

        deallocate(H)
        deallocate(RWORK)
        deallocate(WORK)

    end function calc_eigenvalues



    subroutine get_ham(this,k, ham)
        implicit none
        class(unit_cell), intent(in)              :: this 
        real(8), dimension(3), intent(in)         :: k
        complex(8), dimension(:,:), intent(out)   :: ham
        integer(4)  :: i, i_up, i_d, j, j_up, j_d, conn
        real(8)                                   :: k_dot_r

        if(k(3) /= 0d0) then
            write (*,*) "K_z is non-zero. Abort."
            stop
        endif

        ham =  0d0
        
        do i =  1,2*this%num_atoms
            ham(i,i) =  this%E_s 
        enddo

        ! Spin up
        do i = 1,this%num_atoms
            do conn = 1,this%atoms(i)%n_neigh
                j =  this%atoms(i)%neigh(conn)
                k_dot_r =  dot_product(k, this%atoms(i)%neigh_conn(conn,:))
                ham(i,j) = ham(i,j) + exp(i_unit * k_dot_r) &
                                    * this%atoms(i)%hopping(conn)
                ham(j,i) = conjg(ham(i,j))
            enddo
        enddo

        ! Spin down
        do i = 1,this%num_atoms
            i_d =  i + this%num_atoms
            do conn = 1,this%atoms(i)%n_neigh

                j      = this%atoms(i)%neigh(conn)
                j_d = j + this%num_atoms
                
                k_dot_r =  dot_product(k, this%atoms(i)%neigh_conn(conn,:))
                ham(i_d,j_d) = ham(i_d,j_d) + exp(i_unit * k_dot_r) &
                                            * this%atoms(i)%hopping(conn)  
                ham(j_d,i_d) = conjg(ham(i_d,j_d))
            enddo
        enddo
    end subroutine get_ham

    function init_unit(cfg) result(ret)
        implicit none
        type(unit_cell)              :: ret
        type(CFG_t),  intent(inout)  :: cfg
        integer(4)                   :: hex_sz
        real(8)                      :: tmp

        call CFG_get(cfg, "grid%hexagon_size", ret%hex_size)
        ret%num_atoms =  calc_num_atoms(ret%hex_size)
        allocate(ret%atoms(ret%num_atoms))

        call CFG_get(cfg, "grid%unit_cell_dim", tmp)
        ret%unit_cell_dim = tmp * get_unit_conv("length", cfg)

        call CFG_get(cfg, "hamil%E_s", tmp)
        ret%E_s =  tmp * get_unit_conv("energy", cfg)

        call CFG_get(cfg, "hamil%in_plane_hopping", tmp)
        ret%in_plane_hopping =  tmp * get_unit_conv("energy", cfg)

        call ret%setup_hexagon()
        call ret%setup_conn_1D_layer()
        call ret%setup_lattice_vec()

    end function init_unit

    Subroutine  setup_hexagon(this)
        Implicit None
        class(unit_cell), intent(inout)   :: this
        integer(4), dimension(2)          :: start, pos, halt, dir
        integer(4)                        :: cnt, row
        type(atom)                        :: test

        cnt =  1
        ! sweep from top to (including) middle
        start =  (/0, this%hex_size /)
        halt  =  (/this%hex_size , this%hex_size/)
        dir   =  (/1,0 /)

        do row = this%hex_size,0,-1
            pos =  start 
            do while(any(pos /= halt))
                this%atoms(cnt) = init_ferro(pos)
                cnt =  cnt + 1
                pos =  pos + dir
            end do
            
            start = start +  (/-1, -1/)
            halt  = halt +  (/0, -1/)
        end do


        ! sweep after middle downwards
        start =  (/-this%hex_size, -1/)
        halt  =  (/this%hex_size -1, -1/)

        do row =  -1,-(this%hex_size-1), - 1
            pos =  start
            do while(any(pos /= halt))
                this%atoms(cnt) = init_ferro(pos)
                cnt =  cnt +  1
                pos =  pos +  dir
            enddo

            start =  start +  (/0, -1 /)
            halt  =  halt  +  (/-1, -1 /)
        enddo

    End Subroutine setup_hexagon

    Subroutine  setup_conn_1D_layer(this)
        Implicit None
        class(unit_cell), intent(inout)       :: this 
        integer(4), dimension(2), parameter   :: conn1 =  (/0, 1 /)
        integer(4), dimension(2), parameter   :: conn2 =  (/1, 0 /)
        integer(4), dimension(2), parameter   :: conn3 =  (/1, 1 /)
        integer(4), dimension(2)              :: start_pos 
        integer(4)                            :: i,j
        real(8)                               :: base_len
        real(8), dimension(3,3)               :: base_vecs

        base_len =  this%unit_cell_dim / this%hex_size
        base_vecs(1, :) = (/ 1d0,   0d0,                  0d0 /)
        base_vecs(2, :) = (/ 0.5d0, sin(60d0/180d0 * PI), 0d0 /)
        base_vecs(3, :) = (/-0.5d0, sin(60d0/180d0 * PI), 0d0 /)
        base_vecs =  base_vecs *  base_len

        write (*,*) "base_vecs(1,:)"
        call print_mtx(base_vecs(1,:))

        do i =  1,this%num_atoms
            allocate(this%atoms(i)%neigh(3))
            allocate(this%atoms(i)%hopping(3))
            allocate(this%atoms(i)%neigh_conn(3,3))
            
            this%atoms(i)%hopping = this%in_plane_hopping
            this%atoms(i)%n_neigh =  3
            start_pos             = this%atoms(i)%pos
            
            this%atoms(i)%neigh(1) = this%find_neigh(start_pos, conn1)
            this%atoms(i)%neigh(2) = this%find_neigh(start_pos, conn2)
            this%atoms(i)%neigh(3) = this%find_neigh(start_pos, conn3)

            this%atoms(i)%neigh_conn =  base_vecs 
        enddo

    End Subroutine setup_conn_1D_layer

    function find_neigh(this, start, conn) result(neigh)
        implicit none
        class(unit_cell), intent(in)           :: this
        integer(4), dimension(2), intent(in)   :: start, conn 
        integer(4)                             :: neigh
        integer(4), dimension(2)               :: trans1, trans2, new
        integer(4)                             :: i,j

        trans1 =  (/this%hex_size, - this%hex_size /)
        trans2 =  (/2*this%hex_size, this%hex_size /)

        if(this%in_hexagon(start+conn) /= -1) then
            neigh = this%in_hexagon(start +  conn)
            return
        else
            do i = -1,1
                do j =  -1,1
                    new =  start  +  conn &
                        +  i *  trans1 +  j *  trans2
                    if(this%in_hexagon(new) /= - 1) then
                        neigh =  this%in_hexagon(new)
                        return
                    endif
                enddo
            enddo
        endif

        write (*,*) "Couldn't find a neighbour"
        stop 2

    end function find_neigh

    function in_hexagon(this, pos) result(idx)
        ! if position is in hexagon the corresponding index is
        ! returned, else - 1
        implicit none
        class(unit_cell), intent(in)          :: this
        integer(4), dimension(2), intent(in)  :: pos
        integer(4)                            :: idx 
        integer(4)                            :: i

        idx =  -1
        do i =  1, this%num_atoms
            if(all(pos == this%atoms(i)%pos)) then
                idx =  i
                exit
            endif
        enddo
    end function in_hexagon

    function calc_num_atoms(hex_size) result(num_atoms)
        implicit none
        integer(4), intent(in)       :: hex_size
        integer(4)                   :: num_atoms

        if (hex_size > 0) then
            num_atoms =  3 * hex_size *  hex_size  
        else if(hex_size ==  0) then
            num_atoms =  1
        else
            write(*,*) "Invalid hex_size"
            stop 1 
        end if
    end function calc_num_atoms

    function get_num_atoms(this) result(num)
        implicit none
        Class(unit_cell), intent(in) :: this
        integer(4) :: num
        num = this%num_atoms 
    end function get_num_atoms

    function get_atoms(this) result(ret)
        implicit none
        Class(unit_cell), intent(in)            :: this
        type(atom), dimension(:), allocatable   :: ret

        ret =  this%atoms        
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
            case ("a0^-1") =  1.0d0
                factor = 1.0d0
            case default
                write (*,*) "Unit unknown"
                stop
        end select
    end function get_unit_conv

end module

