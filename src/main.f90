program STB
   use Class_k_space
   use m_config
   use m_npy
   use output
   use mpi
   use Constants
   use Class_unit_cell
   use mypi

   implicit none
   type(k_space)                   :: Ksp
   character(len=*), parameter     :: time_fmt =  "(A,F10.3,A)"
   integer                         :: n_inp, n_files, seed_sz, start_idx, end_idx, cnt
   integer   , allocatable         :: seed(:)
   integer                         :: ierr, me
   character(len=300), allocatable :: inp_files(:)
 
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)

   call random_seed(size = seed_sz)
   allocate(seed(seed_sz))
   seed =  7
   call random_seed(put=seed)

   call get_inp_files(n_files, inp_files)

   call MPI_Bcast(n_files, 1, MYPI_INT, root, MPI_COMM_WORLD, ierr)
   
   do n_inp = 1, n_files
      if(me == root) write (*,*) "started at ", date_time()
      call process_file(inp_files(n_inp))
   enddo

   call MPI_Finalize(ierr)
contains
   subroutine process_file(inp_file)
      implicit none
      character(len=300), intent(in) :: inp_file
      real(8)                        :: start, halt
      integer                        :: me, ierr
      logical                        :: perform_band, perform_dos, calc_hall,&
                                        calc_orbmag, perform_ACA,plot_omega,pert_log
      type(CFG_t)                     :: cfg
      character(len=25)               :: fermi_type

      call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
      start =  MPI_Wtime()

      if(me ==  root)then
         write (*,*) "running: ", trim(inp_file)
         call CFG_read_file(cfg, trim(inp_file))

         call add_full_cfg(cfg)

         call CFG_get(cfg, "band%perform_band", perform_band)
         call CFG_get(cfg, "dos%perform_dos",   perform_dos)
         call CFG_get(cfg, "berry%fermi_type",  fermi_type)
         call CFG_get(cfg, "berry%pert_log",  pert_log)
         call CFG_get(cfg, "berry%calc_hall",   calc_hall)
         call CFG_get(cfg, "berry%calc_orbmag", calc_orbmag)
         call CFG_get(cfg, "ACA%perform_ACA",   perform_ACA)
         call CFG_get(cfg, "plot%plot_omega",   plot_omega)
      endif

      call MPI_Bcast(perform_band, 1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(perform_dos,  1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(fermi_type,   25, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(calc_hall,    1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(calc_orbmag,  1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(perform_ACA,  1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(plot_omega,   1,  MPI_LOGICAL,   root, MPI_COMM_WORLD, ierr)

      Ksp =  init_k_space(cfg)
      if(me == root) call save_cfg(cfg)

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
         call Ksp%calc_berry_quantities(pert_log)
      endif

      if(perform_ACA) then
         if(root == me) write (*,*) "Started ACA"
         call Ksp%calc_ACA()
      endif

      if(plot_omega) then
         call error_msg("plotting Berry curvature...", p_color=c_green)
         call Ksp%plot_omega_square()
      endif

      halt = MPI_Wtime()
      if(root ==  me) then
         write (*,time_fmt) "Total: ", halt-start, "s"
      endif
      call Ksp%free_ksp()
      if(me == root) call CFG_clear(cfg)

   end subroutine process_file

   subroutine get_inp_files(n_files, inp_files)
      implicit none
      integer, intent(out)     :: n_files
      character(len=300), allocatable :: inp_files(:)
      integer                  :: me, ierr, i
      character(len=300)       :: base_str, tmp_str, start_str, end_str, n_files_str

      call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)

      if(me == root) then
         if(command_argument_count() == 0) then
            call error_msg("need some input", abort=.True.)
         elseif(command_argument_count() == 1) then
            n_files = 1
            allocate(inp_files(1))
            call get_command_argument(1, inp_files(1))
         elseif(command_argument_count() == 2) then
            call get_command_argument(2, n_files_str)
            read(n_files_str,*) n_files
            allocate(inp_files(n_files))
            call get_command_argument(1, base_str)

            do n_inp = 0,n_files-1
               write(tmp_str, "(I6)") n_inp
               inp_files(n_inp+1) = trim(base_str) // trim(adjustl(tmp_str)) // ".cfg"
            enddo
         elseif(command_argument_count() == 3) then
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
      call MPI_Bcast(n_files, 1, MYPI_INT, root, MPI_COMM_WORLD, ierr)
      if(me /= root) allocate(inp_files(n_files))

      do i=1,n_files
         call MPI_Bcast(inp_files(i), 300, MPI_CHARACTER,&
                        root, MPI_COMM_WORLD, ierr)
      enddo

   end subroutine get_inp_files

   Subroutine  add_full_cfg(cfg)
      Implicit None
      type(CFG_t)            :: cfg
      real(8), allocatable   :: empty_array(:)
      allocate(empty_array(0))

      call CFG_add(cfg, "units%length",     "none", "unit of length")
      call CFG_add(cfg, "units%energy",     "none", "unit of Ener")
      call CFG_add(cfg, "units%inv_energy", "none", "unit of inverse Ener")
      call CFG_add(cfg, "units%inv_length", "none", "unit of inverse length")
      call CFG_add(cfg, "units%temperature", "none", "unit of temperature")
      call CFG_add(cfg, "units%mag_dipol", "none", "unit of magnetic dipol")

      call CFG_add(cfg, "hamil%Vss_sig",  0d0, "s-state hopping")
      call CFG_add(cfg, "hamil%Vpp_sig",  0d0, "p-state del_l=0 hopping")
      call CFG_add(cfg, "hamil%Vpp_pi",   0d0, "p-state del_l=+/-1 hopping")
      call CFG_add(cfg, "hamil%V2pp_sig", 0d0, "2nd nearest neigh")
      call CFG_add(cfg, "hamil%V2pp_pi",  0d0, "2nd nearest neigh")

      call CFG_add(cfg, "hamil%t_2",      0d0,     "2nd nearest neighbour hopping")
      call CFG_add(cfg, "hamil%phi_2",    0d0,     "2nd nearest neighbourh hopping phase")
      call CFG_add(cfg, "hamil%t_so",     0d0,     "rashba spin-orbit")
      call CFG_add(cfg, "hamil%eta_soc",  0d0,     "p-state SOC")
      call CFG_add(cfg, "hamil%E_s",      0d0,     "EigenE")
      call CFG_add(cfg, "hamil%E_A",      0d0,     "A-site eigenE")
      call CFG_add(cfg, "hamil%E_B",      0d0,     "B-site eigenE")
      call CFG_add(cfg, "hamil%E_p",      [0d0,    0d0, 0d0], "Energies of px, py, pz")
      call CFG_add(cfg, "hamil%lambda",   0d0,     "xc-splitting")
      call CFG_add(cfg, "hamil%n",        0,       "n=0 -> s-states; n=1 -> p-states")
      call CFG_add(cfg, "hamil%molecule", .False., "should this be calculated as an isolated molecule")

      call CFG_add(cfg, "hamil%HB1",    0d0, "Hongbin nearest neigh hopping")
      call CFG_add(cfg, "hamil%HB2",    0d0, "Hongbin snd nearest neigh")
      call CFG_add(cfg, "hamil%HB_eta", 0d0, "Hongbin SOC")
      call CFG_add(cfg, "hamil%lambda_KM", 0d0, "Kane Mele parameter")

      call CFG_add(cfg, "grid%atoms_per_dim",    -1,    "")
      call CFG_add(cfg, "grid%unit_cell_type",   "",    "")
      call CFG_add(cfg, "grid%lattice_constant", 0d0,   "")
      call CFG_add(cfg, "grid%epsilon",          1d-6,  "positional accurary for lattice")
      call CFG_add(cfg, "grid%mag_type",         "",    "")
      call CFG_add(cfg, "grid%winding_number",   1,     "winding number/topological charge")
      call CFG_add(cfg, "grid%ferro_phi",        0d0,   "ferro mag polar angle")
      call CFG_add(cfg, "grid%ferro_theta",      0d0,   "ferro mag azimut angle")
      call CFG_add(cfg, "grid%atan_fac",         0d0,   "wall steepness")
      call CFG_add(cfg, "grid%skyrm_middle",     0.5d0, "skyrm position")
      call CFG_add(cfg, "grid%dblatan_width",    0d0,   "plateau width")
      call CFG_add(cfg, "grid%mag_file",         "",    "mag input file")
      call CFG_add(cfg, "grid%anticol_phi",empty_array, "anticollinear polar angle", dynamic_size=.True.)
      call CFG_add(cfg, "grid%anticol_theta",empty_array, "anticollinear azimutal angle", dynamic_size=.True.)

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
      call CFG_add(cfg, "dos%k_pts_per_dim",    0,       "density of k-grid")
      call CFG_add(cfg, "dos%delta_broadening", 0d0,     "DOS broadening")
      call CFG_add(cfg, "dos%num_points",       300,     "")
      call CFG_add(cfg, "dos%lower_E_bound",    0d0,     "")
      call CFG_add(cfg, "dos%upper_E_bound",    0d0,     "")
      call CFG_add(cfg, "berry%fermi_type", "fixed",      "")
      call CFG_add(cfg, "berry%E_fermi",          [-10d0, 10d0, 300d0],   "")
      call CFG_add(cfg, "dos%fermi_fill",       0.5d0,   "")
      call CFG_add(cfg, "berry%pert_log", .False., "")

      call CFG_add(cfg, "berry%calc_hall", .False., "")
      call CFG_add(cfg, "berry%pert_log", .False., "")
      call CFG_add(cfg, "berry%calc_orbmag", .False., "")
      call CFG_add(cfg, "berry%k_pts_per_dim", 25, "inital density of k-grid")
      call CFG_add(cfg, "berry%temperature", 1d-5, "")
      call CFG_add(cfg, "berry%refinement_iter", 1,"number of refinement iter")
      call CFG_add(cfg, "berry%kpts_per_step", 1, "")
      call CFG_add(cfg, "grid%k_shift", [0d0, 0d0, 0d0], "shift brellouin zone")
      call CFG_add(cfg, "berry%conv_criterion", 0d0, "")
      call CFG_add(cfg, "berry%perform_pad", .True., "padding to use all procs")
      call CFG_add(cfg, "berry%weights", "hall", "use this quantity for refinement")
      call CFG_add(cfg, "berry%adaptive_mode", "area", "how to choose new points")

      call CFG_add(cfg, "output%band_prefix", "bar/foo","folder")

      call CFG_add(cfg, "general%test_run", .False., "performed tests")

      call CFG_add(cfg, "ACA%perform_ACA", .False., "perform ACA calculation")
      call CFG_add(cfg, "ACA%num_kpts", 0, "number of ACA k-points")

      call CFG_add(cfg, "layer_dropout%Vx", empty_array, "Vx dropout", dynamic_size=.True.)
      call CFG_add(cfg, "layer_dropout%Vy", empty_array, "Vy dropout", dynamic_size=.True.)

      call CFG_add(cfg, "plot%num_plot_points", 0, "number of points used for plot")
      call CFG_add(cfg, "plot%plot_omega", .False., "perform berry curvature plot")
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
