module Class_unit_cell
    implicit none
    private
    public: : unit_cell, get_num_atom

    type unit_cell
        ! number of non-redundant atoms pre unit cell
        integer(4)                  :: num_atoms, hex_size
        real(8)                     :: phi, theta
        integer(4), dimension(3)    :: neighbours
        real(8),    dimension(3)    :: pos
    contains
        procedure :: get_num_atoms =>  get_num_atoms
    end type unit_cell
contains
        function get_num_atoms(this) result(num)
            Class(unit_cell), intent(in) :: this
            integer(4) :: num
            num = this%num_atoms 
        end function get_num_atoms

        function calc_num_atoms(hex_size) result(num_atoms)
            integer(4), intent(in)       :: hex_size
            integer(4), intent(out)      :: num_atoms

            if (hex_size /= 1) then
                num_atoms =  3 *  hex_size**2 -  5 *  hex_size +  2
            else
                num_atoms =  1
            end if
        end function calc_num_atoms

end module

