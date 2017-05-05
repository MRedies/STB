program STB
    use Class_k_space
    use m_config
    use output
    
    implicit none
    integer(4)                              :: num, i
    type(k_space)                           :: Ksp
    type(atom), dimension(:), allocatable   :: atoms 
    type(CFG_t)                             :: input_cfg
    character(len=300)                      :: arg
    real(8), dimension(:,:), allocatable    :: eig_val
    real(8), dimension(4)                   :: A
    complex(8), dimension(:,:), allocatable:: H

    call CFG_update_from_arguments(input_cfg)
    call add_full_cfg(input_cfg)
    call CFG_write(input_cfg, "stdout")
    
    Ksp =  init_k_space(input_cfg)
    call Ksp%calc_and_print_band() 
    !do i =  1,2*Ksp%UC%get_num_atoms()
        !write (*,*) i, eig_val(1,i)
    !enddo
    
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

        call CFG_add(cfg, "kspace%k_label", (/ ""/), "", &
                    dynamic_size=.true.)
        call CFG_add(cfg, "kspace%k_x", (/1.0d0/), "", dynamic_size=.true.)
        call CFG_add(cfg, "kspace%k_y", (/1.0d0/), "", dynamic_size=.true.)
        call CFG_add(cfg, "kspace%num_points", 0, "")
        call CFG_add(cfg, "kspace%filling", "","")

        call CFG_add(cfg, "output%outfile", "foo/bar", "")
        call CFG_add(cfg, "output%band_prefix", "bar/foo","")
    End Subroutine add_full_cfg
end program STB

