program STB
    use Class_k_space
    use m_config
    use m_npy
    use output
    use mpi
    use Constants
    use Class_unit_cell

    implicit none
    type(k_space)                   :: Ksp
    type(CFG_t)                     :: cfg
    character(len=25)               :: fermi_type
    character(len=*), parameter     :: time_fmt =  "(A,F10.3,A)"
    integer(4)                      :: ierr, me, n_inp, n_files, seed_sz, start_idx, end_idx, cnt
    integer(4), allocatable         :: seed(:)
    real(8)                         :: start, halt
    logical                         :: perform_band, perform_dos, calc_hall, calc_orbmag
    character(len=300), allocatable :: inp_files(:)

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
    
    call random_seed(size = seed_sz)
    allocate(seed(seed_sz))
    seed =  7
    call random_seed(put=seed)
    
    call get_inp_files(n_files, inp_files)

    call MPI_Bcast(n_files, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    do n_inp = 1, n_files
        start =  MPI_Wtime()

        if(me ==  root)then
            write (*,*) "running: ", trim(inp_files(n_inp))
            call CFG_read_file(cfg, trim(inp_files(n_inp)))
            call add_full_cfg(cfg)

            call CFG_get(cfg, "band%perform_band", perform_band)
            call CFG_get(cfg, "dos%perform_dos",   perform_dos)
            call CFG_get(cfg, "berry%fermi_type", fermi_type) 
            call CFG_get(cfg, "berry%calc_hall", calc_hall)
            call CFG_get(cfg, "berry%calc_orbmag", calc_orbmag)
        endif

        call MPI_Bcast(perform_band, 1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(perform_dos,  1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(fermi_type,   25, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(calc_hall,    1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(calc_orbmag,  1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)

        Ksp =  init_k_space(cfg)
        if(me == root) call save_cfg(cfg)
        
        !call Ksp%plot_omega()
        if(me == root) write (*,*) "num atm", Ksp%ham%UC%num_atoms

        halt =  MPI_Wtime()
        if(root ==  me) then
            write (*,time_fmt) "Init: ", halt-start, "s"
        endif

        if(perform_band) then
            if(root == me) write (*,*) "started Band"
            call Ksp%calc_and_print_band() 
        endif

        if(trim(fermi_type) == "fixed") then
            call Ksp%set_fermi(cfg)
        endif

        if(perform_dos) then
            if(root == me) write (*,*) "started DOS"
            call Ksp%calc_and_print_dos()

            ! Only set Fermi energy relative if DOS was performed
            if(trim(fermi_type) == "filling") then
                call Ksp%find_fermi(cfg)
            endif

        endif

        if(calc_hall .or. calc_orbmag) then
            if(root == me) write (*,*) "started Berry"
            call Ksp%calc_berry_quantities()
        endif

        halt = MPI_Wtime()
        if(root ==  me) then
            write (*,time_fmt) "Total: ", halt-start, "s"
        endif
        call Ksp%free_ksp()
        if(me == root) call CFG_clear(cfg)
    enddo

    call MPI_Finalize(ierr)
contains
    subroutine get_inp_files(n_files, inp_files)
        implicit none
        integer(4), intent(out)  :: n_files 
        character(len=300), allocatable :: inp_files(:)
        integer(4)               :: me, ierr
        character(len=300)       :: base_str, tmp_str, start_str, end_str, n_files_str
        
        call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
    
        if(me == root) then
            if(command_argument_count() == 1) then
                n_files = 1
                allocate(inp_files(1))
                call get_command_argument(1, inp_files(1))
            else if(command_argument_count() == 2) then
                call get_command_argument(2, n_files_str)
                read(n_files_str,*) n_files
                allocate(inp_files(n_files))
                call get_command_argument(1, base_str)

                do n_inp = 0,n_files-1
                    write(tmp_str, "(I6)") n_inp
                    inp_files(n_inp+1) = trim(base_str) // trim(adjustl(tmp_str)) // ".cfg"
                enddo
            else if(command_argument_count() == 3) then
                call get_command_argument(2, start_str)
                call get_command_argument(3, end_str)
                read(start_str, *) start_idx
                read(end_str, * ) end_idx
                n_files =  end_idx -  start_idx + 1
                allocate(inp_files(n_files))
                call get_command_argument(1, base_str)

                cnt = 1
                do n_inp = start_idx, end_idx
                    write(tmp_str, "(I6)") n_inp
                    inp_files(cnt) = trim(base_str) // trim(adjustl(tmp_str)) // ".cfg"
                    cnt = cnt + 1
                enddo
            endif
        endif

    end subroutine get_inp_files
        
    Subroutine  add_full_cfg(cfg)
        Implicit None
        type(CFG_t)            :: cfg 

        call CFG_add(cfg, "units%length",     "none", "")
        call CFG_add(cfg, "units%energy",     "none", "")
        call CFG_add(cfg, "units%inv_energy", "none", "")
        call CFG_add(cfg, "units%inv_length", "none", "")
        call CFG_add(cfg, "units%temperature", "none", "")
        call CFG_add(cfg, "units%mag_dipol", "none", "")

        call CFG_add(cfg, "hamil%Vss_sig",   0d0, "")
        call CFG_add(cfg, "hamil%Vsp_sig",   0d0, "")
        call CFG_add(cfg, "hamil%Vpp_sig",   0d0, "")
        call CFG_add(cfg, "hamil%Vpp_pi",    0d0, "")

        call CFG_add(cfg, "hamil%t_2",       0d0, "")
        call CFG_add(cfg, "hamil%phi_2",     0d0, "")
        call CFG_add(cfg, "hamil%t_so",      0d0, "")
        call CFG_add(cfg, "hamil%E_s",       0d0, "")
        call CFG_add(cfg, "hamil%E_A",       0d0, "")
        call CFG_add(cfg, "hamil%E_B",       0d0, "")
        call CFG_add(cfg, "hamil%lambda",    0d0, "")
        call CFG_add(cfg, "hamil%lambda_nl", 0d0, "")
        call CFG_add(cfg, "hamil%n_max",     0,   "")

        call CFG_add(cfg, "grid%atoms_per_dim",    -1,    "")
        call CFG_add(cfg, "grid%unit_cell_type",   "",    "")
        call CFG_add(cfg, "grid%lattice_constant", 0d0,   "")
        call CFG_add(cfg, "grid%epsilon",          1d-6,  "")
        call CFG_add(cfg, "grid%mag_type",         "",    "")
        call CFG_add(cfg, "grid%ferro_phi",        0d0,   "")
        call CFG_add(cfg, "grid%ferro_theta",      0d0,   "")
        call CFG_add(cfg, "grid%atan_fac",         0d0,   "")
        call CFG_add(cfg, "grid%skyrm_middle",     0.5d0, "")
        call CFG_add(cfg, "grid%dblatan_width",    0d0,   "")
        call CFG_add(cfg, "grid%mag_randomness",   0d0,   "")
        call CFG_add(cfg, "grid%number_of_layers", 1,     "")
        call CFG_add(cfg, "grid%stacking_type",    "",    "")
        call CFG_add(cfg, "grid%layer_height",     1d0,   "")
        call CFG_add(cfg, "grid%mag_file",         "",    "")

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
        call CFG_add(cfg, "berry%fermi_type",       "",      "")
        call CFG_add(cfg, "berry%E_fermi",          [0d0, 0d0, 0d0],   "")
        call CFG_add(cfg, "dos%fermi_fill",       0.5d0,   "")

        call CFG_add(cfg, "berry%calc_hall", .False., "")
        call CFG_add(cfg, "berry%calc_orbmag", .False., "")
        call CFG_add(cfg, "berry%k_pts_per_dim", 25, "")
        call CFG_add(cfg, "berry%temperature", 1d-5, "")
        call CFG_add(cfg, "berry%laplace_iter", 0, "")
        call CFG_add(cfg, "berry%refinement_iter", 1,"")
        call CFG_add(cfg, "berry%kpts_per_step", 1, "")
        call CFG_add(cfg, "grid%k_shift", [0d0, 0d0, 0d0], "")
        call CFG_add(cfg, "berry%conv_criterion", 0d0, "")
        call CFG_add(cfg, "berry%perform_pad", .True., "")
        call CFG_add(cfg, "berry%weights", "", "")

        call CFG_add(cfg, "output%band_prefix", "bar/foo","")
    End Subroutine add_full_cfg

    subroutine save_cfg(cfg)
        Implicit None
        type(CFG_t)            :: cfg
        character(len=300)     :: prefix

        call CFG_add(cfg, "calculation%start_time", date_time(), "")
        call CFG_get(cfg, "output%band_prefix", prefix)
        call CFG_write(cfg, trim(prefix) // "setup.cfg")

        end subroutine save_cfg
end program STB

