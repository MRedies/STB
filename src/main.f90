program STB
    use Class_k_space
    use m_config
    use output
    use Constants
    
    implicit none
    type(k_space)                           :: Ksp
    type(CFG_t)                             :: cfg
    character(len=25)                       :: fermi_type
    logical :: perform_band, perform_dos

    call CFG_update_from_arguments(cfg)
    call add_full_cfg(cfg)
    !call CFG_write(cfg, "stdout")
    
    call CFG_get(cfg, "band%perform_band", perform_band)
    call CFG_get(cfg, "dos%perform_dos",   perform_dos)
    call CFG_get(cfg, "dos%fermi_type", fermi_type) 

    Ksp =  init_k_space(cfg)
    if(perform_band) then
        call Ksp%calc_and_print_band(cfg) 
    endif
    if(perform_dos) then
        call Ksp%calc_and_print_dos()

        ! Only set Fermi energy if DOS was performed
        if(trim(fermi_type) == "fixed") then
            call Ksp%set_fermi(cfg)
        else if(trim(fermi_type) == "filling") then
            call Ksp%find_fermi(cfg)
        endif
        
        call Ksp%write_fermi()
    endif
contains
    Subroutine  add_full_cfg(cfg)
        Implicit None
        type(CFG_t)            :: cfg 

        call CFG_add(cfg, "units%length",     "none", "")
        call CFG_add(cfg, "units%energy",     "none", "")
        call CFG_add(cfg, "units%inv_energy", "none", "")
        call CFG_add(cfg, "units%inv_length", "none", "")

        call CFG_add(cfg, "hamil%t_nn", 0.0d0, "")
        call CFG_add(cfg, "hamil%E_s",  0.0d0, "")
        call CFG_add(cfg, "hamil%I",    0d0,   "")

        call CFG_add(cfg, "grid%atoms_per_dim", -1, "")
        call CFG_add(cfg, "grid%unit_cell_type","","")
        call CFG_add(cfg, "grid%lattice_constant", 0d0, "")
        call CFG_add(cfg, "grid%epsilon", 1d-6, "")
        call CFG_add(cfg, "grid%mag_type", "", "")

        call CFG_add(cfg, "band%perform_band",  .False., "")
        call CFG_add(cfg, "band%k_label",    (/ ""/),   "",&
                          dynamic_size=.true.)
        call CFG_add(cfg, "band%k_x",        (/1.0d0/), "",&
                          dynamic_size=.true.)
        call CFG_add(cfg, "band%k_y",        (/1.0d0/), "",&
                          dynamic_size=.true.)
        call CFG_add(cfg, "band%num_points", 0,         "")
        call CFG_add(cfg, "band%filling",    "",        "")

        call CFG_add(cfg, "dos%perform_dos",      .False., "")
        call CFG_add(cfg, "dos%k_pts_per_dim",    0,       "")
        call CFG_add(cfg, "dos%delta_broadening", 0d0,     "")
        call CFG_add(cfg, "dos%num_points",       300,     "")
        call CFG_add(cfg, "dos%lower_E_bound",    0d0,     "")
        call CFG_add(cfg, "dos%upper_E_bound",    0d0,     "")
        call CFG_add(cfg, "dos%fermi_type",       "",      "")
        call CFG_add(cfg, "dos%E_fermi",          0d0,     "")
        call CFG_add(cfg, "dos%fermi_fill",       0.5d0,   "")

        call CFG_add(cfg, "output%band_prefix", "bar/foo","")
    End Subroutine add_full_cfg
end program STB

