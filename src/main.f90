program STB
    use Class_unit_cell
    use m_config
    use output
    
    implicit none
    integer(4)                              :: num, i
    type(unit_cell)                         :: UC
    type(atom), dimension(:), allocatable   :: atoms 
    type(CFG_t)                             :: input_cfg
    character(len=300)                      :: arg
    real(8), dimension(3,1)                :: k
    real(8), dimension(:,:), allocatable    :: eig_val
    real(8), dimension(4)                   :: A
    complex(8), dimension(:,:), allocatable:: H

    do i =  1,size(k,2)
        k(:,i) = (/0.5d0 +  0.0001d0 *  i, 0d0, 0d0/)
    enddo
       
    call CFG_update_from_arguments(input_cfg)
    call add_full_cfg(input_cfg)
    
    UC =  init_unit(input_cfg)
    
    allocate(H(2*UC%get_num_atoms(), 2*UC%get_num_atoms()))
    H =  0d0

    call UC%get_ham(k(:,1), H)

    open(unit=2, file="bla.txt")
    call print_mtx(2,H)
    close(2)
    eig_val =  UC%calc_eigenvalues(k)

    do i =  1,2*UC%get_num_atoms()
        write (*,*) i, eig_val(1,i)
    enddo
    
contains
    Subroutine  add_full_cfg(cfg)
        Implicit None
        type(CFG_t)            :: cfg 

        call CFG_add(cfg, "units%length", "none", "")
        call CFG_add(cfg, "units%energy", "none", "")
        call CFG_add(cfg, "units%inv_energy", "none", "")

        call CFG_add(cfg, "hamil%in_plane_hopping", 0.0d0, "")
        call CFG_add(cfg, "hamil%E_s", 0.0d0, "")

        call CFG_add(cfg, "grid%hexagon_size", 0, "")
        call CFG_add(cfg, "grid%unit_cell_dim", 0d0, "")

        call CFG_add(cfg, "output%outfile", "foo/bar", "")
    End Subroutine add_full_cfg
end program STB

