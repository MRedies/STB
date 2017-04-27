program STB
    use Class_unit_cell
    implicit none
    integer(4)      :: num  
    type(unit_cell) :: bla
    bla =  unit_cell(3,2.3,4.5)
    num =  bla%get_num_atom()
    
    write(*,*) "Num: ", num 
    
end program STB

