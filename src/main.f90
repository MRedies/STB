program STB
    use Class_unit_cell
    use m_config
    !use mkl_service
    implicit none
    integer(4)                              :: num, i
    type(unit_cell)                         :: UC
    type(atom), dimension(:), allocatable   :: atoms 
    type(CFG_t)                             :: input_cfg
    character(len=300)                      :: arg
   
    call CFG_update_from_arguments(input_cfg)
    
    !call UC%init(200)
    !atoms =  UC%get_atoms()
    !write (*,*) "Index    ", "connected"
    !do i =  1, UC%get_num_atoms()
        !write (*,*) i, atoms(i)%neigh
    !enddo
    !call test_herm()
contains
end program STB

