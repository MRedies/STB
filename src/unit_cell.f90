module Class_unit_cell
    use Class_atom
    implicit none

    type unit_cell
        ! number of non-redundant atoms pre unit cell
        private
        integer(4)                              :: num_atoms, hex_size
        type(atom), dimension(:), allocatable   :: atoms 
    contains
        procedure :: get_num_atoms => get_num_atoms
        procedure :: init          => init
        procedure :: setup_hexagon => setup_hexagon
        procedure :: find_neigh    => find_neigh
        procedure :: in_hexagon    => in_hexagon
        procedure :: setup_conn    => setup_conn
        procedure :: get_atoms     =>  get_atoms
    end type unit_cell
contains
    Subroutine  init(this, hex_sz)
        ! Stupid Fortran style constructor
        implicit none
        class(unit_cell), intent(out)   :: this
        integer(4), intent(in)          :: hex_sz

        this%num_atoms =  calc_num_atoms(hex_sz)
        allocate(this%atoms(this%num_atoms))

        this%hex_size =  hex_sz
        
        call this%setup_hexagon()
        call this%setup_conn()

    End Subroutine init

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

    Subroutine  setup_conn(this)
        Implicit None
        class(unit_cell), intent(inout)       :: this 
        integer(4), dimension(2), parameter   :: conn1 =  (/0, 1 /)
        integer(4), dimension(2), parameter   :: conn2 =  (/1, 0 /)
        integer(4), dimension(2), parameter   :: conn3 =  (/1, 1 /)
        integer(4), dimension(2)              :: start_pos
        integer(4)                            :: i

        do i =  1,this%num_atoms
            start_pos =  this%atoms(i)%pos
            this%atoms(i)%neigh(1) = this%find_neigh(start_pos, conn1)
            this%atoms(i)%neigh(2) = this%find_neigh(start_pos, conn2)
            this%atoms(i)%neigh(3) = this%find_neigh(start_pos, conn3)
        enddo


    End Subroutine setup_conn

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
            num_atoms =  0
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
end module

