program STB
    use Class_k_space
    use m_config
    use output
    use mpi
    use Constants
    use Class_unit_cell
    
    implicit none
    type(k_space)      :: Ksp
    type(CFG_t)        :: cfg
    character(len=25)  :: fermi_type
    character(len=*), parameter :: time_fmt =  "(A,F10.3,A)"
    integer(4)         :: ierr, me
    real(8)            :: start, halt
    logical :: perform_band, perform_dos, calc_hall
    real(8), allocatable :: hall_cond(:)
    complex(8)           :: A(100,100), b(100), x(100)
    real(8)              :: tmp_M(100,100), tmp_vec(100)

  
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)

    start =  MPI_Wtime()

    if(me ==  root)then
        call CFG_update_from_arguments(cfg)
        call add_full_cfg(cfg)
        
        call CFG_get(cfg, "band%perform_band", perform_band)
        call CFG_get(cfg, "dos%perform_dos",   perform_dos)
        call CFG_get(cfg, "dos%fermi_type", fermi_type) 
        call CFG_get(cfg, "berry%calc_hall", calc_hall)
    endif
    
    call MPI_Bcast(perform_band, 1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(perform_dos,  1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(fermi_type,   25, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(calc_hall,    1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
    
    Ksp =  init_k_space(cfg)
    if(me == root) write (*,*) "num atm", Ksp%ham%UC%num_atoms

    halt =  MPI_Wtime()
    if(root ==  me) then
        write (*,time_fmt) "Init: ", halt-start, "s"
    endif
    
    if(perform_band) then
        call Ksp%calc_and_print_band() 
    endif

    if(trim(fermi_type) == "fixed") then
        call Ksp%set_fermi(cfg)
    endif


    if(perform_dos) then
        call Ksp%calc_and_print_dos()

        ! Only set Fermi energy relative if DOS was performed
        if(trim(fermi_type) == "filling") then
            call Ksp%find_fermi(cfg)
        endif
        
    endif
    
    if(calc_hall) then
        call Ksp%calc_hall_conductance(hall_cond)
    endif
    halt = MPI_Wtime()
    if(root ==  me) then
        write (*,time_fmt) "Total: ", halt-start, "s"
    endif
    call MPI_Finalize(ierr)
contains
    Subroutine  add_full_cfg(cfg)
        Implicit None
        type(CFG_t)            :: cfg 

        call CFG_add(cfg, "units%length",     "none", "")
        call CFG_add(cfg, "units%energy",     "none", "")
        call CFG_add(cfg, "units%inv_energy", "none", "")
        call CFG_add(cfg, "units%inv_length", "none", "")
        call CFG_add(cfg, "units%temperature", "none", "")

        call CFG_add(cfg, "hamil%t_nn",      0.0d0, "")
        call CFG_add(cfg, "hamil%t_so",      0d0,   "")
        call CFG_add(cfg, "hamil%E_s",       0.0d0, "")
        call CFG_add(cfg, "hamil%lambda",    0d0,   "")
        call CFG_add(cfg, "hamil%lambda_nl", 0d0,   "")

        call CFG_add(cfg, "grid%atoms_per_dim",    -1,   "")
        call CFG_add(cfg, "grid%unit_cell_type",   "",   "")
        call CFG_add(cfg, "grid%lattice_constant", 0d0,  "")
        call CFG_add(cfg, "grid%epsilon",          1d-6, "")
        call CFG_add(cfg, "grid%mag_type",         "",   "")
        call CFG_add(cfg, "grid%ferro_phi",        0d0,  "")
        call CFG_add(cfg, "grid%ferro_theta",      0d0,  "")
        call CFG_add(cfg, "grid%atan_fac",         0d0,  "")

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
        call CFG_add(cfg, "dos%E_fermi",          [0d0, 0d0, 0d0],   "")
        call CFG_add(cfg, "dos%fermi_fill",       0.5d0,   "")

        call CFG_add(cfg, "berry%calc_hall", .False., "")
        call CFG_add(cfg, "berry%k_pts_per_dim", 25, "")
        call CFG_add(cfg, "berry%temperature", 1d-5, "")

        call CFG_add(cfg, "output%band_prefix", "bar/foo","")
    End Subroutine add_full_cfg
end program STB

