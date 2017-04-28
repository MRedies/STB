program STB
    use Class_unit_cell
    implicit none
    integer(4)      :: num, i
    type(unit_cell) :: UC
    type(atom), dimension(:), allocatable   :: atoms 

    call UC%init(200)
    atoms =  UC%get_atoms()
    write (*,*) "Index    ", "connected"
    do i =  1, UC%get_num_atoms()
        write (*,*) i, atoms(i)%neigh
    enddo
end program STB

