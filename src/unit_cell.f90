module Class_unit_cell
   use Class_atom
   use Class_helper
   use m_config
   use output
   use m_npy
   use mpi
   use mypi
   use Constants
   use class_Units

   implicit none

   enum, bind(c)
      enumerator :: nn_conn, snd_nn_conn, oop_conn
   end enum

   type unit_cell
      real(8), public :: lattice(2, 2) !> translation vectors
      !> of the real-space lattice. First index: element of vector
      !> Second index: 1 or second vector
      real(8), public :: rez_lattice(2, 2) !> translation vectors
      !> of the reciprocal lattice. Indexs same as lattice
      ! number of non-redundant atoms pre unit cell
      integer    :: num_atoms  !> number of non-redundant atoms in a unit cell
      integer    :: atom_per_dim !> atoms along the radius of the unit_cell
      integer    :: nProcs
      integer    :: me
      integer    :: n_wind !> winding number for lin_rot
      integer, allocatable    :: wavevector(:)
      real(8) :: lattice_constant !> lattice constant in atomic units
      real(8) :: eps !> threshold for positional accuracy
      real(8) :: ferro_phi, ferro_theta
      real(8) :: cone_angle
      real(8), allocatable:: anticol_phi(:), anticol_theta(:), m0_A(:), m0_B(:) !> the angles for anticollinear setups, one
      real(8), allocatable:: axis(:) !> the angles for anticollinear setups, one
      real(8) :: atan_factor !> how fast do we change the border wall
      real(8) :: atan_pref !> prefactor for atan
      real(8) :: dblatan_dist !> width of the atan plateau
      real(8) :: dblatan_pref !> prefactor for angle turning
      real(8) :: skyrm_middle !> position of inplane
      type(atom), dimension(:), allocatable :: atoms !> array containing all atoms
      type(units)       :: units
      character(len=25) :: uc_type !> indicates shape of unitcell
      character(len=25) :: mag_type !> indicates type of magnetization
      character(len=25) :: spiral_type !> indicates type of magnetization
      character(len=300):: mag_file
      logical :: molecule !> should we have a k-space or not?
      logical      :: test_run !> should unit tests be performed
      logical         :: pert_log ! berry in first order pert.
   contains

      procedure :: init_unit_honey_hexa => init_unit_honey_hexa
      procedure :: init_unit_square => init_unit_square
      procedure :: init_unit_honey_line => init_unit_honey_line
      procedure :: get_num_atoms => get_num_atoms
      procedure :: setup_square => setup_square
      procedure :: setup_single_hex => setup_single_hex
      procedure :: in_cell => in_cell
      procedure :: setup_gen_conn => setup_gen_conn
      procedure :: get_atoms => get_atoms
      procedure :: gen_find_neigh => gen_find_neigh
      procedure :: save_unit_cell => save_unit_cell
      procedure :: set_mag_ferro => set_mag_ferro
      procedure :: set_mag_anticol => set_mag_anticol
      procedure :: set_mag_random => set_mag_random
      procedure :: set_mag_x_spiral_square => set_mag_x_spiral_square
      procedure :: set_mag_linrot_1D_spiral => set_mag_linrot_1D_spiral
      procedure :: set_mag_linrot_1D_spiral_honey => set_mag_linrot_1D_spiral_honey
      procedure :: set_mag_linrot_1D_spiral_m0_anticol => set_mag_linrot_1D_spiral_m0_anticol
      procedure :: set_mag_linrot_1D_spiral_m0_cone => set_mag_linrot_1D_spiral_m0_cone
      procedure :: set_mag_linrot_skyrm => set_mag_linrot_skyrm
      procedure :: set_mag_atan_skyrm => set_mag_atan_skyrm
      procedure :: set_mag_atan_skyrm_honey => set_mag_atan_skyrm_honey
      procedure :: set_mag_dblatan_skyrm => set_mag_dblatan_skyrm
      procedure :: set_mag_dblatan_skyrm_honey => set_mag_dblatan_skyrm_honey
      procedure :: set_mag_linrot_skrym_square => set_mag_linrot_skrym_square
      procedure :: set_mag_linrot_skrym_honey => set_mag_linrot_skrym_honey
      procedure :: set_honey_snd_nearest => set_honey_snd_nearest
      procedure :: set_honey_snd_nearest_line => set_honey_snd_nearest_line
      procedure :: find_lattice_vectors => find_lattice_vectors
      procedure :: find_conn_vectors => find_conn_vectors
      procedure :: set_mag_site => set_mag_site
      procedure :: Bcast_UC => Bcast_UC
      procedure :: setup_honey => setup_honey
      procedure :: make_hexagon => make_hexagon
      procedure :: make_honeycomb_line => make_honeycomb_line
      procedure :: free_uc => free_uc
      procedure :: init_file_square => init_file_square
      procedure :: run_tests => run_tests
      procedure :: calc_area => calc_area
   end type unit_cell
contains
   subroutine free_uc(self)
      implicit none
      class(unit_cell)  :: self
      integer           :: i

      if (allocated(self%atoms)) then
         do i = 1, self%num_atoms
            call self%atoms(i)%free_atm()
         enddo
         deallocate (self%atoms)
      endif
   end subroutine free_uc

   function angle(a, b) result(ang)
      implicit none
      real(8), dimension(2), intent(in)   :: a, b
      real(8)                             :: ang
      ang = dot_product(a, b)/(my_norm2(a)*my_norm2(b))
      ang = 180.0d0/PI*acos(ang)
   end function angle

   function init_unit(cfg) result(self)
      implicit none
      type(CFG_t)       :: cfg !> config file as read by m_config
      type(unit_cell)   :: self
      integer, parameter           :: lwork = 20
      real(8)                         :: work(lwork), tmp
      integer, dimension(2)        :: ipiv
      integer                         :: info
      integer                         :: ierr
      integer                         :: anticol_size, wavevector_size
      logical                         :: tmp_log

      call MPI_Comm_size(MPI_COMM_WORLD, self%nProcs, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)

      self%units = init_units(cfg, self%me)

      if (self%me == 0) then
         call CFG_get(cfg, "berry%pert_log", tmp_log)
         self%pert_log = tmp_log
         call CFG_get(cfg, "grid%epsilon", tmp)
         self%eps = tmp*self%units%length

         call CFG_get(cfg, "hamil%molecule", self%molecule)
         if (self%molecule) then
            call error_msg("Using Molecule mode!", p_color=c_green)
         endif

         call CFG_get(cfg, "grid%mag_type", self%mag_type)
         call CFG_get_size(cfg, "grid%anticol_phi", anticol_size)
         allocate (self%anticol_phi(anticol_size)) !allocate phi
         call CFG_get(cfg, "grid%anticol_phi", self%anticol_phi)
         call CFG_get_size(cfg, "grid%anticol_theta", anticol_size)
         allocate (self%anticol_theta(anticol_size)) !allocate theta
         call CFG_get(cfg, "grid%anticol_theta", self%anticol_theta)

         call CFG_get(cfg, "grid%winding_number", self%n_wind)
         call CFG_get(cfg, "grid%unit_cell_type", self%uc_type)
         call CFG_get_size(cfg, "grid%wavevector", wavevector_size)
         allocate (self%wavevector(wavevector_size))
         call CFG_get(cfg, "grid%wavevector", self%wavevector)
         call CFG_get_size(cfg, "grid%axis", wavevector_size)
         allocate (self%axis(wavevector_size))
         call CFG_get(cfg, "grid%axis", self%axis)
         call CFG_get(cfg, "grid%cone_angle", self%cone_angle)
         call CFG_get(cfg, "grid%spiral_type", self%spiral_type)
         call CFG_get(cfg, "grid%lattice_constant", tmp)
         self%lattice_constant = tmp*self%units%length

         call CFG_get(cfg, "grid%atoms_per_dim", self%atom_per_dim)

         call CFG_get(cfg, "grid%ferro_phi", self%ferro_phi)
         call CFG_get(cfg, "grid%ferro_theta", self%ferro_theta)
         call CFG_get(cfg, "grid%atan_fac", self%atan_factor)
         call CFG_get(cfg, "grid%atan_pref", self%atan_pref)
         call CFG_get(cfg, "grid%skyrm_middle", self%skyrm_middle)
         call CFG_get(cfg, "grid%dblatan_width", self%dblatan_dist)
         call CFG_get(cfg, "grid%dblatan_pref", self%dblatan_pref)

         call CFG_get(cfg, "grid%mag_file", self%mag_file)
         call CFG_get(cfg, "general%test_run", self%test_run)


      endif

      call self%Bcast_UC()

      if (trim(self%uc_type) == "square_2d") then
         call self%init_unit_square()
      else if (trim(self%uc_type) == "honey_2d") then
         call self%init_unit_honey_hexa()
      else if (trim(self%uc_type) == "honey_line") then
         call self%init_unit_honey_line()
      else if (trim(self%uc_type) == "file_square") then
         call self%init_file_square()
      else
         write (*, *) self%me, ": Cell type unknown"
         stop
      endif

      ! calculate reciprocal grid
      self%rez_lattice = transpose(self%lattice)
      call dgetrf(2, 2, self%rez_lattice, 2, ipiv, info)
      if (info /= 0) then
         write (*, *) self%me, ": LU-decomp of lattice vectors failed", info
         stop
      endif

      call dgetri(2, self%rez_lattice, 2, ipiv, work, lwork, info)
      if (info /= 0) then
         write (*, *) self%me, ": Inversion of lattice vectors failed"
         stop
      endif
      self%rez_lattice = 2*PI*self%rez_lattice

      if (self%test_run) call self%run_tests()
   end function

   subroutine Bcast_UC(self)
      use mpi
      implicit none
      class(unit_cell)           :: self
      integer, parameter         :: num_cast = 26
      integer                    :: ierr(num_cast)
      integer                    :: anticol_size_phi, wsize, asize
      integer                    :: anticol_size_theta

      if (self%me == root) then
         anticol_size_phi = size(self%anticol_phi)
         anticol_size_theta = size(self%anticol_theta)
         wsize = size(self%wavevector)
         asize = size(self%axis)
      endif
      call MPI_Bcast(self%eps, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(1))
      call MPI_Bcast(self%mag_type, 25, MPI_CHARACTER, &
                     root, MPI_COMM_WORLD, ierr(2))
      call MPI_Bcast(self%uc_type, 25, MPI_CHARACTER, &
                     root, MPI_COMM_WORLD, ierr(3))
      call MPI_Bcast(self%lattice_constant, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(4))
      call MPI_Bcast(self%atom_per_dim, 1, MYPI_INT, &
                     root, MPI_COMM_WORLD, ierr(5))

      call MPI_Bcast(self%ferro_phi, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(6))
      call MPI_Bcast(self%ferro_theta, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(7))
      call MPI_Bcast(self%atan_factor, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(8))
      call MPI_Bcast(self%dblatan_dist, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(9))
      call MPI_Bcast(self%skyrm_middle, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(10))

      call MPI_Bcast(self%n_wind, 1, MYPI_INT, &
                     root, MPI_COMM_WORLD, ierr(11))
      call MPI_Bcast(self%molecule, 1, MPI_LOGICAL, &
                     root, MPI_COMM_WORLD, ierr(12))
      call MPI_Bcast(self%test_run, 1, MPI_LOGICAL, &
                     root, MPI_COMM_WORLD, ierr(13))

      call MPI_Bcast(anticol_size_phi, 1, MYPI_INT, &
                     root, MPI_COMM_WORLD, ierr(14))
      call MPI_Bcast(anticol_size_theta, 1, MYPI_INT, &
                     root, MPI_COMM_WORLD, ierr(15))
      if (self%me /= root) then
         allocate (self%anticol_phi(anticol_size_phi))
         allocate (self%anticol_theta(anticol_size_theta))
      endif
      call MPI_Bcast(self%anticol_theta, anticol_size_theta, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(16))
      call MPI_Bcast(self%anticol_phi, anticol_size_phi, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(17))

      call MPI_Bcast(self%pert_log, 1, MPI_LOGICAL, &
                     root, MPI_COMM_WORLD, ierr(18))

      call MPI_Bcast(wsize, 1, MYPI_INT, root, MPI_COMM_WORLD, ierr(19))
      call MPI_Bcast(asize, 1, MYPI_INT, root, MPI_COMM_WORLD, ierr(20))
      if (self%me /= root) then
         allocate (self%wavevector(wsize))
         allocate (self%axis(asize))
      endif
      call MPI_Bcast(self%wavevector, wsize, MYPI_INT, root, MPI_COMM_WORLD, ierr(21))
      call MPI_Bcast(self%axis, asize, MPI_REAL8, root, MPI_COMM_WORLD, ierr(22))
      call MPI_Bcast(self%cone_angle, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr(23))
      call MPI_Bcast(self%spiral_type, 25, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr(24))

      call MPI_Bcast(self%dblatan_pref, 1, MPI_REAL8, &
                     root, MPI_COMM_WORLD, ierr(25))
      call MPI_Bcast(self%atan_pref, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr(26))
      call check_ierr(ierr, self%me, "Unit cell check err")
   end subroutine Bcast_UC

   subroutine init_unit_square(self)
      implicit none
      class(unit_cell), intent(inout) :: self
      real(8)                         :: conn_mtx(5, 3), transl_mtx(4, 3)
      integer(4)                      :: conn_types(size(conn_mtx, dim=1))

      self%num_atoms = self%atom_per_dim**2
      allocate (self%atoms(self%num_atoms))

      call self%setup_square()

      conn_mtx(1, :) = self%lattice_constant*[1d0, 0d0, 0d0]
      conn_mtx(2, :) = self%lattice_constant*[0d0, 1d0, 0d0]
      conn_mtx(3, :) = self%lattice_constant*[0d0, 0d0, 1d0]
      conn_mtx(4, :) = self%lattice_constant*[1d0, 1d0, 0d0]
      conn_mtx(5, :) = self%lattice_constant*[1d0, -1d0, 0d0]
      conn_types = [nn_conn, nn_conn, nn_conn, snd_nn_conn, snd_nn_conn]

      transl_mtx = 0d0
      transl_mtx(1, 1:2) = self%lattice(:, 1)
      transl_mtx(2, 1:2) = self%lattice(:, 2)
      transl_mtx(3, 1:2) = self%lattice(:, 1) + self%lattice(:, 2)
      transl_mtx(4, 1:2) = self%lattice(:, 1) - self%lattice(:, 2)

      call self%setup_gen_conn(conn_mtx, conn_types, transl_mtx)

      if (trim(self%mag_type) == "x_spiral") then
         call self%set_mag_x_spiral_square()
      else if (trim(self%mag_type) == "ferro") then
         call self%set_mag_ferro()
      else if (trim(self%mag_type) == "lin_skyrm") then
         call self%set_mag_linrot_skrym_square()
      else if (trim(self%mag_type) == "random") then
         call self%set_mag_random()
      else if (trim(self%mag_type) == "anticol") then
         call self%set_mag_anticol()
      else
         write (*, *) "Mag_type = ", trim(self%mag_type)
         call error_msg("mag_type not known", abort=.True.)
      endif
   end subroutine init_unit_square

   subroutine init_file_square(self)
      use mpi
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8)                           :: conn_mtx(3, 3)
      real(8), allocatable              :: transl_mtx(:, :), m(:, :), pos(:, :)
      integer                           :: n(3), i, n_transl
      integer                           :: info
      character(len=300)                :: garb

      if (self%me == root) then
         write (*, *) "file =  ", trim(self%mag_file)
         open (unit=21, file=trim(self%mag_file))
         read (21, *) garb, n(1), n(2), n(3)
         write (*, *) n
      endif
      call MPI_Bcast(n, 3, MYPI_INT, root, MPI_COMM_WORLD, info)
      self%num_atoms = n(1)*n(2)*n(3)

      allocate (self%atoms(self%num_atoms))
      allocate (m(3, self%num_atoms))
      allocate (pos(3, self%num_atoms))

      do i = 1, self%num_atoms
         if (self%me == root) then
            read (21, *) pos(1, i), pos(2, i), pos(3, i), m(1, i), m(2, i), m(3, i)
         endif
      enddo

      call MPI_Bcast(pos, int(3*self%num_atoms, 4), MPI_REAL8, &
                     root, MPI_COMM_WORLD, info)
      call MPI_Bcast(m, int(3*self%num_atoms, 4), MPI_REAL8, &
                     root, MPI_COMM_WORLD, info)

      pos = pos*self%lattice_constant

      do i = 1, self%num_atoms
         self%atoms(i) = init_ferro_z(pos(:, i))
         call self%atoms(i)%set_m_cart(m(1, i), m(2, i), m(3, i))
      enddo

      if (self%me == root) read (21, *) garb, n_transl
      call MPI_Bcast(n_transl, 1, MYPI_INT, root, MPI_COMM_WORLD, info)
      allocate (transl_mtx(n_transl, 3))

      if (self%me == root) then
         do i = 1, n_transl
            read (21, *) transl_mtx(i, 1), transl_mtx(i, 2), transl_mtx(i, 3)
         enddo
         close (21)
      endif

      !if we want a molecule, ensure that no wrap-around is found
      if (self%molecule) transl_mtx = transl_mtx*10d0

      call MPI_Bcast(transl_mtx, int(3*n_transl, 4), MPI_REAL8, root, MPI_COMM_WORLD, info)

      conn_mtx(1, :) = (/self%lattice_constant, 0d0, 0d0/)
      conn_mtx(2, :) = (/0d0, self%lattice_constant, 0d0/)
      conn_mtx(3, :) = (/0d0, 0d0, self%lattice_constant/)

      self%lattice(:, 1) = self%lattice_constant*transl_mtx(:, 1)
      self%lattice(:, 2) = self%lattice_constant*transl_mtx(:, 2)

      call self%setup_gen_conn(conn_mtx, [nn_conn, nn_conn, nn_conn], transl_mtx)
      deallocate (m, pos)
   end subroutine init_file_square

   subroutine make_hexagon(self, hexagon, site_type)
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8)                          :: transl_mtx(3, 3), base_len_uc, l, pos(3)
      real(8), allocatable             :: hexagon(:, :), grid(:, :)
      integer                          :: num_atoms, cnt, apd, i
      integer, allocatable             :: site_type(:)

      apd = self%atom_per_dim
      base_len_uc = self%lattice_constant*apd
      l = 2*cos(deg_30)*base_len_uc

      transl_mtx(1, :) = l*[1d0, 0d0, 0d0]
      transl_mtx(2, :) = l*[0.5d0, sin(deg_60), 0d0]
      transl_mtx(3, :) = l*[0.5d0, -sin(deg_60), 0d0]
      num_atoms = calc_num_atoms_non_red_honey(apd)

      allocate (hexagon(num_atoms, 3))
      allocate (site_type(num_atoms))
      hexagon = 0d0
      call gen_honey_grid(self%lattice_constant, apd, grid)

      cnt = 1
      do i = 1, size(grid, 1)
         pos = grid(i, :)
         if (in_hexagon(pos, base_len_uc)) then
            if (.not. already_in_red(pos, hexagon, cnt - 1, transl_mtx)) then
               hexagon(cnt, :) = grid(i, :)
               if (i <= (2*apd + 1)**2) then
                  site_type(cnt) = A_site
               else
                  site_type(cnt) = B_site
               endif
               cnt = cnt + 1
            endif
         endif
      enddo
      deallocate (grid)
   end subroutine make_hexagon

   subroutine make_honeycomb_line(self, line, site_type)
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8), allocatable              :: line(:, :), conn_vecs(:, :)
      integer, allocatable              :: site_type(:)
      real(8)                           :: shift_mtx(3, 3), conn_mtx(3, 3), transf_mtx(3, 3), base_len_uc, posA(3), &
                                           posB(3), posC(3), posD(3),pos(3), conn_vec_1(3), conn_vec_2(3), l
      integer                           :: i, ii, ierr

      if (mod(self%num_atoms, 2) /= 0) then
         write (*, *) "number of atoms in honey_comb line has to be even"
         call MPI_Abort(MPI_COMM_WORLD, 23, ierr)
      endif
      base_len_uc = self%lattice_constant
      l = 2*cos(deg_30)*base_len_uc
      !conn to nearest neighbors
      conn_mtx(1, :) = self%lattice_constant*[0d0, 1d0, 0d0]!1
      conn_mtx(2, :) = self%lattice_constant*[cos(deg_30), -sin(deg_30), 0d0]!2
      conn_mtx(3, :) = self%lattice_constant*[-cos(deg_30), -sin(deg_30), 0d0]!3
      !conn to next honey unit cell
      shift_mtx(1, :) = l*[1d0, 0d0, 0d0]!1
      shift_mtx(2, :) = l*[0.5d0, sin(deg_60), 0d0]!2
      shift_mtx(3, :) = l*[0.5d0, -sin(deg_60), 0d0]!3
      !change of coordinates, matmul(transpose(shift_mtx),wavevector) = matmul(transpose(conn_mtx),matmul(transf_mtx,wavevector)
      transf_mtx(1, :) = [1d0, 2d0, -1d0]
      transf_mtx(2, :) = [2d0, 1d0, 1d0]
      transf_mtx(3, :) = [0d0, 0d0, 0d0]
      !so far only spirals along connection vectors
      if (trim(self%mag_type) == "1Dspiral") then
         call self%find_conn_vectors(conn_vecs)
         conn_vec_1 = conn_vecs(1, :)
         conn_vec_2 = conn_vecs(2, :)
      else
         conn_vec_1 = shift_mtx(1, :)
         conn_vec_2 = shift_mtx(2, :)
      endif
      allocate (line(self%num_atoms, 3))
      allocate (site_type(self%num_atoms))
      posA = -0.5*self%lattice_constant*[0d0,1d0,0d0]!conn_mtx(3, :)
      posB = 0.5*self%lattice_constant*[0d0,1d0,0d0]!-conn_mtx(2, :)
      posC = posA + conn_mtx(3, :)
      posD = posB - conn_mtx(3, :)
      pos = -1d0*(self%atom_per_dim - 1)/2d0*conn_vec_1
      do i = 1, self%atom_per_dim
         ii = 4*(i - 1)
         line(ii + 1, :) = posA + pos
         site_type(ii + 1) = A_site
         line(ii + 2, :) = posB + pos
         site_type(ii + 2) = B_site
         line(ii + 3, :) = posC + pos
         site_type(ii + 3) = B_site
         line(ii + 4, :) = posD + pos
         site_type(ii + 4) = A_site
         pos = pos + conn_vec_1
      enddo
   end subroutine make_honeycomb_line

   subroutine find_lattice_vectors(self,lattice)
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8)                           :: temp(3), shift_mtx(3, 3), wave_proj(3), proj, check, l
      real(8), allocatable              :: lattice(:, :)
      integer                           :: i

      l = 2*cos(deg_30)*self%lattice_constant
      shift_mtx(1, :) = l*[1d0, 0d0, 0d0]!1
      shift_mtx(2, :) = l*[0.5d0, sin(deg_60), 0d0]!2
      shift_mtx(3, :) = l*[0.5d0, -sin(deg_60), 0d0]!3
      allocate(lattice(2, 3))
      temp = matmul(transpose(shift_mtx), self%wavevector)
      lattice(1, :) = self%atom_per_dim*temp
      !construct second perpendicular lattice vector of correct length
      temp = cross_prod(temp,[0d0,0d0,1d0])
      !check if temp is parallel to any of the shift_mtx
      if (norm2(cross_prod(temp,shift_mtx(1, :))) < pos_eps) then
         lattice(2, :) = shift_mtx(1, :)
      elseif (norm2(cross_prod(temp,shift_mtx(2, :))) < pos_eps) then
         lattice(2, :) = shift_mtx(2, :)
      elseif (norm2(cross_prod(temp,shift_mtx(3, :))) < pos_eps) then
         lattice(2, :) = shift_mtx(3, :)
      else
         wave_proj = matmul(shift_mtx,temp)
         check = 100d0
         do i=1, 3
            proj = abs(wave_proj(i))
            if (proj < check .AND. proj > pos_eps) then
               check = proj
            endif
         enddo
         wave_proj = wave_proj/check
         do i=1, 3
            if (wave_proj(i)-int(wave_proj(i)) > pos_eps) then
               write(*,*) "Coefficients of 2nd lattice vector are not integer!", wave_proj
            endif
         enddo
         lattice(2, :) = matmul(transpose(shift_mtx),wave_proj)
      endif
   end subroutine find_lattice_vectors

   subroutine find_conn_vectors(self,conn_vecs)
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8)                           :: conn_mtx(3, 3), shift_mtx(3, 3), conn_proj(3) ,conn_vec_1(3), conn_vec_2(3) &
                                           , l
      real(8), allocatable              :: conn_vecs(:,:)
      integer                           :: ierr
      
      l = 2*cos(deg_30)*self%lattice_constant
      !conn to next honey neighbor
      conn_mtx(1, :) = self%lattice_constant*[0d0, 1d0, 0d0]!1
      conn_mtx(2, :) = self%lattice_constant*[cos(deg_30), -sin(deg_30), 0d0]!2
      conn_mtx(3, :) = self%lattice_constant*[-cos(deg_30), -sin(deg_30), 0d0]!3
      !conn to next honey unit cell
      shift_mtx(1, :) = l*[1d0, 0d0, 0d0]!1
      shift_mtx(2, :) = l*[0.5d0, sin(deg_60), 0d0]!2
      shift_mtx(3, :) = l*[0.5d0, -sin(deg_60), 0d0]!3
      conn_vec_1 = matmul(transpose(shift_mtx), self%wavevector)
      if (norm2(conn_vec_1) < pos_eps) then
         write(*,*) "q-vector is zero!", conn_vec_1
         call MPI_Abort(MPI_COMM_WORLD, 23, ierr)
      endif
      conn_proj = matmul(conn_mtx, conn_vec_1)
      if (abs(conn_proj(1) - conn_proj(2)) < pos_eps) then
         conn_vec_2 = conn_mtx(3, :)
      elseif (abs(conn_proj(3) - conn_proj(2)) < pos_eps) then
         conn_vec_2 = conn_mtx(1, :)
      elseif (abs(conn_proj(1) - conn_proj(3)) < pos_eps) then
         conn_vec_2 = conn_mtx(2, :)
      elseif (conn_proj(1) > conn_proj(2) .AND. conn_proj(1) > conn_proj(3)) then
         conn_vec_2 = conn_mtx(1, :)
      elseif (conn_proj(2) > conn_proj(1) .AND. conn_proj(2) > conn_proj(3)) then
         conn_vec_2 = conn_mtx(2, :)
      elseif (conn_proj(3) > conn_proj(2) .AND. conn_proj(3) > conn_proj(1)) then
         conn_vec_2 = conn_mtx(3, :)
      endif
      allocate(conn_vecs(2, 3))
      conn_vecs(1, :) = conn_vec_1
      conn_vecs(2, :) = conn_vec_2
   end subroutine find_conn_vectors

   subroutine init_unit_honey_line(self)
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8)                           :: transl_mtx(3, 3), conn_mtx(3, 3), shift_mtx(3, 3)
      real(8)                           :: base_len_uc, l
      real(8), allocatable              :: lattice(:, :), line(:, :)
      integer, allocatable              :: site_type(:)
      integer                           :: apd, check_idx
      apd = self%atom_per_dim
      self%num_atoms = calc_num_atoms_line_honey(apd)
      base_len_uc = self%lattice_constant
      l = 2*cos(deg_30)*base_len_uc
      !conn to next honey unit cell
      shift_mtx(1, :) = l*[1d0, 0d0, 0d0]!1
      shift_mtx(2, :) = l*[0.5d0, sin(deg_60), 0d0]!2
      shift_mtx(3, :) = l*[0.5d0, -sin(deg_60), 0d0]!3
      !spiral uc lat vecs
      call self%find_lattice_vectors(lattice)
      self%lattice(:, :) = lattice(:, 1:2)
      if (dot_product(self%lattice(1, :),self%lattice(2, :)) > pos_eps) then
         write(*,*) "Lattice vectors are not orthogonal!", dot_product(self%lattice(1, :),self%lattice(2, :))
      endif
      !if we want a molecule, ensure that no wrap-around is found
      if (self%molecule) transl_mtx = transl_mtx*10d0
      allocate (self%atoms(self%num_atoms))
      ! only one kind of atoms of the honey-comb unit cell needs
      ! the other comes through complex conjugate
      !translates to neighboring atoms
      conn_mtx(1, :) = self%lattice_constant*[0d0, 1d0, 0d0]
      conn_mtx(2, :) = self%lattice_constant*[cos(deg_30), -sin(deg_30), 0d0]
      conn_mtx(3, :) = self%lattice_constant*[-cos(deg_30), -sin(deg_30), 0d0]
      !translates to neighboring unit cells
      check_idx = 0
      transl_mtx(1,:) = matmul(transpose(shift_mtx),self%wavevector)
      transl_mtx(2,:) = lattice(2, :)
      transl_mtx(3,:) = lattice(1, :) - lattice(2, :)
      call self%make_honeycomb_line(line, site_type)
      call self%setup_honey(line, site_type)
      call self%setup_gen_conn(conn_mtx, [nn_conn, nn_conn, nn_conn], transl_mtx)
      call self%set_honey_snd_nearest_line(transl_mtx)

      if (trim(self%mag_type) == "ferro_uiaeuiaeuia") then
         call self%set_mag_ferro()
      else if (trim(self%mag_type) == "1Dspiral") then
         call self%set_mag_linrot_1D_spiral_honey()
      else
         write (*, *) "Mag_type = ", trim(self%mag_type)
         call error_msg("mag_type not known", abort=.True.)
      endif

   end subroutine init_unit_honey_line

   subroutine init_unit_honey_hexa(self)
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8)  :: transl_mtx(3, 3), l, base_len_uc, conn_mtx(3, 3)
      real(8), allocatable             :: hexagon(:, :)
      integer, allocatable             :: site_type(:)
      integer                          :: apd

      apd = self%atom_per_dim
      base_len_uc = self%lattice_constant*apd
      l = 2*cos(deg_30)*base_len_uc

      transl_mtx(1, :) = l*[1d0, 0d0, 0d0]
      transl_mtx(2, :) = l*[0.5d0, sin(deg_60), 0d0]
      transl_mtx(3, :) = l*[0.5d0, -sin(deg_60), 0d0]

      !if we want a molecule, ensure that no wrap-around is found
      if (self%molecule) transl_mtx = transl_mtx*10d0

      self%lattice(:, 1) = transl_mtx(1, 1:2)
      self%lattice(:, 2) = transl_mtx(2, 1:2)

      self%num_atoms = calc_num_atoms_non_red_honey(apd)
      allocate (self%atoms(self%num_atoms))

      call self%make_hexagon(hexagon, site_type)
      call self%setup_honey(hexagon, site_type)

      ! only one kind of atom from honey-comb unit cell needed
      ! the other comes through complex conjugate
      conn_mtx(1, :) = self%lattice_constant*[0d0, 1d0, 0d0]
      conn_mtx(2, :) = self%lattice_constant*[cos(deg_30), -sin(deg_30), 0d0]
      conn_mtx(3, :) = self%lattice_constant*[-cos(deg_30), -sin(deg_30), 0d0]

      call self%setup_gen_conn(conn_mtx, [nn_conn, nn_conn, nn_conn], transl_mtx)
      call self%set_honey_snd_nearest()

      if (trim(self%mag_type) == "ferro") then
         call self%set_mag_ferro()
      else if (trim(self%mag_type) == "lin_skyrm") then
         call self%set_mag_linrot_skrym_honey()
      else if (trim(self%mag_type) == "atan_skyrm") then
         call self%set_mag_atan_skyrm_honey()
      else if (trim(self%mag_type) == "dblatan_skyrm") then
         call self%set_mag_dblatan_skyrm_honey()
      else if (trim(self%mag_type) == "random") then
         call self%set_mag_random()
      else if (trim(self%mag_type) == "anticol") then
         call self%set_mag_anticol()
      else if (trim(self%mag_type) == "1Dspiral") then
         call self%set_mag_linrot_1D_spiral_honey()
      else
         write (*, *) "Mag_type = ", trim(self%mag_type)
         call error_msg("mag_type not known", abort=.True.)
      endif

      !setup_gen_conn block was here before!

      deallocate (hexagon, site_type)
   end subroutine init_unit_honey_hexa

   subroutine set_honey_snd_nearest_line(self,transl_mtx)
      implicit none
      class(unit_cell)        :: self
      integer                 :: i, j, cand
      real(8)                 :: l, conn_mtx_A(3, 3), conn_mtx_B(3, 3), start_pos(3), &
                                 conn(3), conn_storage(3, 3)
      real(8), allocatable    :: tmp(:, :)
      real(8), intent(in)     :: transl_mtx(3, 3)
      integer                 :: idx(3), curr_size
      l = 2d0*cos(deg_30)*self%lattice_constant
      !transl_mtx(1, :) = apd*l*[1d0, 0d0, 0d0]
      !transl_mtx(2, :) = apd*l*[0.5d0, sin(deg_60), 0d0]
      !transl_mtx(3, :) = apd*l*[0.5d0, -sin(deg_60), 0d0]

      !only clockwise connections
      conn_mtx_A(1, :) = l*[-1d0, 0d0, 0d0]
      conn_mtx_A(2, :) = l*[0.5d0, sin(deg_60), 0d0]
      conn_mtx_A(3, :) = l*[0.5d0, -sin(deg_60), 0d0]

      conn_mtx_B(1, :) = -l*[-1d0, 0d0, 0d0]
      conn_mtx_B(2, :) = -l*[0.5d0, sin(deg_60), 0d0]
      conn_mtx_B(3, :) = -l*[0.5d0, -sin(deg_60), 0d0]

      do i = 1, self%num_atoms
         start_pos = self%atoms(i)%pos

         do j = 1, 3
            if (self%atoms(i)%site_type == A_site) then
               conn = conn_mtx_A(j, :)
            elseif (self%atoms(i)%site_type == B_site) then
               conn = conn_mtx_B(j, :)
            else
               call error_msg("2nd nearest only in honey", abort=.True.)
            endif

            cand = self%gen_find_neigh(start_pos, conn, transl_mtx)
            if (cand /= -1) then
               idx(j) = cand
               conn_storage(j, :) = conn
            else
               if (self%me == root) write (*, *) "couldn't make a match"
            endif
         enddo
         !append 1D-arrays
         self%atoms(i)%neigh_idx = [self%atoms(i)%neigh_idx, idx]
         self%atoms(i)%conn_type = [self%atoms(i)%conn_type, [snd_nn_conn, snd_nn_conn, snd_nn_conn]]
         curr_size = size(self%atoms(i)%neigh_conn, 1)
         allocate (tmp(curr_size + 3, 3))
         tmp(1:curr_size, :) = self%atoms(i)%neigh_conn
         tmp(curr_size + 1:curr_size + 3, :) = conn_storage
         call move_alloc(tmp, self%atoms(i)%neigh_conn)
      enddo

   end subroutine set_honey_snd_nearest_line

   subroutine set_honey_snd_nearest(self)
      implicit none
      class(unit_cell)        :: self
      integer                 :: i, j, cand, apd
      real(8)                 :: l, conn_mtx_A(3, 3), conn_mtx_B(3, 3), start_pos(3), &
                                 conn(3), transl_mtx(3, 3), conn_storage(3, 3)
      real(8), allocatable    :: tmp(:, :)
      integer                 :: idx(3), curr_size
      apd = self%atom_per_dim
      l = 2d0*cos(deg_30)*self%lattice_constant
      transl_mtx(1, :) = apd*l*[1d0, 0d0, 0d0]
      transl_mtx(2, :) = apd*l*[0.5d0, sin(deg_60), 0d0]
      transl_mtx(3, :) = apd*l*[0.5d0, -sin(deg_60), 0d0]

      !only clockwise connections
      conn_mtx_A(1, :) = l*[-1d0, 0d0, 0d0]
      conn_mtx_A(2, :) = l*[0.5d0, sin(deg_60), 0d0]
      conn_mtx_A(3, :) = l*[0.5d0, -sin(deg_60), 0d0]

      conn_mtx_B(1, :) = -l*[-1d0, 0d0, 0d0]
      conn_mtx_B(2, :) = -l*[0.5d0, sin(deg_60), 0d0]
      conn_mtx_B(3, :) = -l*[0.5d0, -sin(deg_60), 0d0]

      do i = 1, self%num_atoms
         start_pos = self%atoms(i)%pos

         do j = 1, 3
            if (self%atoms(i)%site_type == A_site) then
               conn = conn_mtx_A(j, :)
            elseif (self%atoms(i)%site_type == B_site) then
               conn = conn_mtx_B(j, :)
            else
               call error_msg("2nd nearest only in honey", abort=.True.)
            endif

            cand = self%gen_find_neigh(start_pos, conn, transl_mtx)
            if (cand /= -1) then
               idx(j) = cand
               conn_storage(j, :) = conn
            else
               if (self%me == root) write (*, *) "couldn't make a match"
            endif
         enddo
         !append 1D-arrays
         self%atoms(i)%neigh_idx = [self%atoms(i)%neigh_idx, idx]
         self%atoms(i)%conn_type = [self%atoms(i)%conn_type, [snd_nn_conn, snd_nn_conn, snd_nn_conn]]
         curr_size = size(self%atoms(i)%neigh_conn, 1)
         allocate (tmp(curr_size + 3, 3))
         tmp(1:curr_size, :) = self%atoms(i)%neigh_conn
         tmp(curr_size + 1:curr_size + 3, :) = conn_storage
         call move_alloc(tmp, self%atoms(i)%neigh_conn)
      enddo

   end subroutine set_honey_snd_nearest

   function clockwise(v1, v2, v3) result(clock)
      implicit none
      real(8), intent(in)  :: v1(3), v2(3), v3(3)
      real(8)              :: tmp
      logical              :: clock

      tmp = (v2(1) - v1(1))*(v2(2) + v1(2))
      tmp = tmp + (v3(1) - v2(1))*(v3(2) + v2(2))
      tmp = tmp + (v1(1) - v3(1))*(v1(2) + v3(2))

      clock = tmp > 0d0
   end function clockwise

   subroutine set_mag_ferro(self)
      implicit none
      class(unit_cell)    :: self
      integer             :: i

      do i = 1, self%num_atoms
         call self%atoms(i)%set_sphere(self%ferro_phi, self%ferro_theta)
      enddo
   end subroutine set_mag_ferro

   subroutine set_mag_anticol(self)
      implicit none
      class(unit_cell)        :: self
      real(8)                 :: phi, theta, phi_nc, phi_col, theta_nc, theta_col
      integer                 :: i
      if (self%me == root) write (*, *) "This is the collinear perturbation!"
      if (size(self%anticol_phi) /= self%num_atoms &
          .or. size(self%anticol_theta) /= self%num_atoms) then
         call error_msg("sizes of anticol_phi and anticol_theta not consistent with num_atoms", abort=.True.)
      else
         if (size(self%anticol_theta) == 2 &
             .and. size(self%anticol_phi) == 2) then
            phi_col = self%anticol_phi(1)
            phi_nc = self%anticol_phi(2)
            theta_col = self%anticol_theta(1)
            theta_nc = self%anticol_theta(2)
            if (self%pert_log) then
               do i = 1, self%num_atoms
                  call self%atoms(i)%set_sphere(phi_col, theta_col)
               enddo
            else
               do i = 1, self%num_atoms
                  phi = phi_col - (-1)**i*phi_nc/2d0
                  theta = theta_col - (-1)**i*theta_nc/2d0
                  call self%atoms(i)%set_sphere(phi, theta)
               enddo
            endif
         else
            do i = 1, self%num_atoms
               call self%atoms(i)%set_sphere(self%anticol_phi(i), self%anticol_theta(i))
            enddo
         endif
      endif
   end subroutine set_mag_anticol

   subroutine set_mag_linrot_1D_spiral_m0_anticol(self)
      implicit none
      class(unit_cell)        :: self
      real(8)                 :: phi_nc, phi_col, theta_nc, theta_col, thetaA, thetaB, phiA, phiB
      if (mod(self%num_atoms, size(self%anticol_phi)) == 0 &
          .and. mod(self%num_atoms, size(self%anticol_theta)) == 0 &
          .and. size(self%anticol_theta) == 2 &
          .and. size(self%anticol_phi) == 2) then
         allocate (self%m0_A(3))
         allocate (self%m0_B(3))
         phi_col = self%anticol_phi(1)
         phi_nc = self%anticol_phi(2)
         theta_col = self%anticol_theta(1)
         theta_nc = self%anticol_theta(2)
         phiA = phi_col + phi_nc/2d0
         phiB = phi_col - phi_nc/2d0
         thetaA = theta_col + theta_nc/2d0
         thetaB = theta_col - theta_nc/2d0
         self%m0_A(1) = sin(thetaA)*cos(phiA)
         self%m0_A(2) = sin(thetaA)*sin(phiA)
         self%m0_A(3) = cos(thetaA)
         self%m0_B(1) = sin(thetaB)*cos(phiB)
         self%m0_B(2) = sin(thetaB)*sin(phiB)
         self%m0_B(3) = cos(thetaB)
         !call self%atoms(i)%set_sphere(phi,theta)
      else
         call error_msg("sizes of anticol_phi and anticol_theta not consistent with num_atoms", abort=.True.)
      endif
   end subroutine set_mag_linrot_1D_spiral_m0_anticol

   subroutine set_mag_linrot_1D_spiral_m0_cone(self)
      implicit none
      class(unit_cell)        :: self
      real(8)                 :: G(3, 3), axis(3), perp_axis(3), m0(3)
      axis = self%axis
      perp_axis = cross_prod(axis, [0d0, 0d0, 1d0])
      G = R_mtx(self%cone_angle, perp_axis)
      m0 = matmul(G, axis)
      self%m0_A(1) = m0(1)
      self%m0_A(2) = m0(2)
      self%m0_A(3) = m0(3)
      self%m0_B(1) = m0(1)
      self%m0_B(2) = m0(2)
      self%m0_B(3) = m0(3)
   end subroutine set_mag_linrot_1D_spiral_m0_cone
   subroutine set_mag_x_spiral_square(self)
      implicit none
      class(unit_cell)                 :: self
      real(8)        :: alpha, rel_xpos
      integer        :: i

      do i = 1, self%num_atoms
         rel_xpos = self%atoms(i)%pos(1)/self%lattice_constant
         alpha = rel_xpos*2*PI/self%atom_per_dim
         if (alpha <= PI) then
            self%atoms(i)%m_theta = alpha
            self%atoms(i)%m_phi = 0d0
         else
            self%atoms(i)%m_theta = 2*PI - alpha
            self%atoms(i)%m_phi = PI
         endif
      enddo

   end subroutine set_mag_x_spiral_square

   subroutine set_mag_random(self)
      implicit none
      class(unit_cell)       :: self
      integer                :: i
      real(8)                :: phi, theta, r(2)

      do i = 1, self%num_atoms
         call random_number(r)

         phi = r(1)*2d0*PI
         theta = r(2)*PI
         call self%atoms(i)%set_sphere(phi, theta)
      enddo
   end subroutine set_mag_random

   subroutine set_mag_linrot_skrym_square(self)
      implicit none
      class(unit_cell)     :: self
      real(8)              :: radius, center(3)

      ! Nagaosa style unit cell. Ferromagnetic border only to the left
      radius = 0.5d0*self%lattice_constant*self%atom_per_dim
      center = (/radius, radius, 0d0/)

      call self%set_mag_linrot_skyrm(center, radius)
   end subroutine set_mag_linrot_skrym_square

   subroutine set_mag_linrot_skrym_honey(self)
      implicit none
      class(unit_cell)      :: self
      real(8), parameter    :: center(3) = [0d0, 0d0, 0d0]
      real(8)               :: radius

      radius = 0.5d0*my_norm2(self%lattice(:, 1))
      call self%set_mag_linrot_skyrm(center, radius)

   end subroutine set_mag_linrot_skrym_honey

   subroutine set_mag_atan_skyrm_honey(self)
      implicit none
      class(unit_cell)    :: self
      real(8), parameter    :: center(3) = [0d0, 0d0, 0d0]
      real(8)               :: radius

      radius = 0.5d0*my_norm2(self%lattice(:, 1))
      call self%set_mag_atan_skyrm(center, radius)
   end subroutine

   subroutine set_mag_dblatan_skyrm_honey(self)
      implicit none
      class(unit_cell)    :: self
      real(8), parameter    :: center(3) = [0d0, 0d0, 0d0]
      real(8)               :: radius

      radius = 0.5d0*my_norm2(self%lattice(:, 1))
      call self%set_mag_dblatan_skyrm(center, radius)
   end subroutine set_mag_dblatan_skyrm_honey

   subroutine set_mag_linrot_skyrm(self, center, radius)
      implicit none

      class(unit_cell)                      :: self
      real(8), intent(in)           :: center(3), radius
      real(8), parameter            :: e_z(3) = [0, 0, 1]
      real(8)                               :: R(3, 3), conn(3), n(3), m(3), k(3), alpha
      integer                               :: i

      alpha = 0d0
      do i = 1, self%num_atoms
         conn = center - self%atoms(i)%pos
         if (my_norm2(conn) > pos_eps*self%lattice_constant &
             .and. my_norm2(conn) <= radius + pos_eps) then
            k = n_times_phi(conn, self%n_wind)
            n = cross_prod(k, e_z)

            alpha = PI*(1d0 - my_norm2(conn)/radius)
            R = R_mtx(alpha, n)
            ! center of skyrmion point down
            m = matmul(R, e_z)
         else if (my_norm2(conn) <= pos_eps*self%lattice_constant) then
            m = -e_z
         else
            m = e_z
         endif

         call self%atoms(i)%set_m_cart(m(1), m(2), m(3))
      enddo
   end subroutine set_mag_linrot_skyrm

   subroutine set_mag_atan_skyrm(self, center, radius)
      implicit none
      class(unit_cell)      :: self
      real(8), intent(in)   :: center(3), radius
      real(8), parameter    :: e_z(3) = [0, 0, 1]
      real(8)  :: R(3, 3), conn(3), n(3), m(3), alpha, y_min, y_max, x0, x, a, scaling
      integer               :: i

      a = self%atan_factor
      x0 = self%skyrm_middle*radius

      y_min = -atan(-a*x0) + 0.5d0*PI
      y_max = -atan(a*(radius - x0)) + 0.5d0*PI
      scaling = PI/(y_max - y_min)

      alpha = 0d0
      do i = 1, self%num_atoms
         conn = center - self%atoms(i)%pos
         if (my_norm2(conn) > pos_eps*self%lattice_constant &
             .and. my_norm2(conn) <= radius + pos_eps) then
            n = cross_prod(conn, e_z)

            x = my_norm2(conn)
            alpha = PI - scaling*(-atan(a*(x - x0)) + 0.5d0*PI - y_min)

            alpha = self%atan_pref * alpha

            R = R_mtx(alpha, n)
            ! center of skyrmion point down
            m = matmul(R, e_z)
         else if (my_norm2(conn) <= pos_eps*self%lattice_constant) then
            m = -e_z
         else
            m = e_z
         endif
         call self%atoms(i)%set_m_cart(m(1), m(2), m(3))
      enddo
   end subroutine set_mag_atan_skyrm

   subroutine set_mag_dblatan_skyrm(self, center, radius)
      implicit none
      class(unit_cell)      :: self
      real(8), intent(in)   :: center(3), radius
      real(8), parameter    :: e_z(3) = [0, 0, 1]
      real(8)  :: R(3, 3), conn(3), n(3), m(3), alpha, alp_min, alp_max, a, d, x
      integer               :: i

      a = self%atan_factor
      d = self%dblatan_dist
      !x0      = 0.5d0 * radius

      alp_min = (atan(a*(0.5d0 - d*0.5d0)) &
                 + atan(a*(0.5d0 + d*0.5d0)))
      alp_max = (atan(a*(-0.5d0 - d*0.5d0)) &
                 + atan(a*(-0.5d0 + d*0.5d0)))

      alpha = 0d0
      do i = 1, self%num_atoms
         conn = center - self%atoms(i)%pos
         if (my_norm2(conn) > pos_eps*self%lattice_constant &
             .and. my_norm2(conn) <= radius + pos_eps) then
            n = cross_prod(conn, e_z)

            x = my_norm2(conn)/radius
            alpha = atan(a*(x - 0.5d0 - 0.5d0*d)) &
                    + atan(a*(x - 0.5d0 + 0.5d0*d))

            alpha = alpha - alp_min
            alpha = alpha/(alp_max - alp_min)
            alpha = self%dblatan_pref*PI*alpha

            R = R_mtx(alpha, n)
            ! center of skyrmion point down
            m = matmul(R, e_z)
         else if (my_norm2(conn) <= pos_eps*self%lattice_constant) then
            m = -e_z
         else
            m = e_z
         endif
         call self%atoms(i)%set_m_cart(m(1), m(2), m(3))
      enddo

   end subroutine set_mag_dblatan_skyrm

   subroutine set_mag_linrot_1D_spiral_honey(self)
      implicit none
      class(unit_cell)      :: self
      real(8), parameter    :: center(3) = [0d0, 0d0, 0d0]
      real(8)               :: UC_l

      UC_l = my_norm2(self%lattice(:, 1))
      if (self%spiral_type == "anticol") then
         call self%set_mag_linrot_1D_spiral_m0_anticol()
         call self%set_mag_linrot_1D_spiral(center, UC_l)
      elseif (self%spiral_type == "cone") then
         call self%set_mag_linrot_1D_spiral_m0_cone()
         call self%set_mag_linrot_1D_spiral(center, UC_l)
      endif
   end subroutine set_mag_linrot_1D_spiral_honey
   subroutine set_mag_site(self, ii, j, center, UC_l)
      implicit none
      class(unit_cell)    :: self
      integer, intent(in) :: ii, j
      real(8), intent(in) :: center(3), UC_l
      integer             :: site_type, i
      real(8)             :: conn(3), phase_fac, x, l, R(3,3), shift_mtx(3,3), m(3), wavevector(3) &
                             , wavevector_len, wavelength, psi
      
      l = 2*cos(deg_30)*self%lattice_constant
      shift_mtx(1, :) = l*[1d0, 0d0, 0d0]!1
      shift_mtx(2, :) = l*[0.5d0, sin(deg_60), 0d0]!2
      shift_mtx(3, :) = l*[0.5d0, -sin(deg_60), 0d0]!3
      wavevector = matmul(transpose(shift_mtx), self%wavevector)
      wavevector_len = my_norm2(wavevector)
      wavevector = wavevector/wavevector_len
      wavelength = UC_l/(1d0*self%n_wind)
      psi = 2d0*PI/wavelength
      i = ii + j
      site_type = self%atoms(i)%site_type
      conn = self%atoms(i)%pos! - self%atoms(j)%pos
      phase_fac = 0d0!
      if (self%atoms(i)%site_type /= self%atoms(j)%site_type) then
         write(*,*) "Site types do not agree!"
      endif
      if (site_type == 0) then
         x = dot_product(wavevector,conn - phase_fac)
         R = R_mtx(psi*x, self%axis)
         m = matmul(R, self%m0_A)
      elseif (site_type == 1) then
         x = dot_product(wavevector,conn - phase_fac)
         R = R_mtx(psi*x, self%axis)
         m = matmul(R, self%m0_B)
      endif
      call self%atoms(i)%set_m_cart(m(1), m(2), m(3))
   end subroutine set_mag_site

   subroutine set_mag_linrot_1D_spiral(self, center, UC_l)
      implicit none
      class(unit_cell)    :: self
      real(8), intent(in) :: center(3), UC_l
      integer             :: i, ii, j
      do i = 1, self%atom_per_dim
         ii = 4*(i-1)
         do j = 1, 4
            call self%set_mag_site(ii, j, center, UC_l)
         enddo
      enddo
   end subroutine set_mag_linrot_1D_spiral

   subroutine save_unit_cell(self, folder)
      implicit none
      class(unit_cell)        :: self
      character(len=*)        :: folder
      real(8), allocatable    :: x(:), y(:), z(:), phi(:), theta(:)
      integer                 :: i, n_neigh
      integer, allocatable :: neigh(:, :)
      integer, allocatable    :: site_type(:), conn_type(:, :)

      allocate (x(self%num_atoms))
      allocate (y(self%num_atoms))
      allocate (z(self%num_atoms))
      allocate (phi(self%num_atoms))
      allocate (theta(self%num_atoms))
      allocate (site_type(self%num_atoms))
      allocate (neigh(self%num_atoms, 20))
      allocate (conn_type(self%num_atoms, 20))

      neigh = -1
      conn_type = -1
      do i = 1, self%num_atoms

         x(i) = self%atoms(i)%pos(1)
         y(i) = self%atoms(i)%pos(2)
         z(i) = self%atoms(i)%pos(3)

         if (trim(self%mag_type) == "anticol") then
            phi(i) = self%atoms(i)%m_phi
            theta(i) = self%atoms(i)%m_theta
            !phi(i) = self%anticol_phi(i)
            !theta(i) = self%anticol_theta(i)
         else
            phi(i) = self%atoms(i)%m_phi
            theta(i) = self%atoms(i)%m_theta
         endif
         site_type(i) = self%atoms(i)%site_type

         n_neigh = size(self%atoms(i)%neigh_idx)
         neigh(i, 1:n_neigh) = self%atoms(i)%neigh_idx
         conn_type(i, 1:n_neigh) = self%atoms(i)%conn_type
      enddo

      call save_npy(folder//"pos_x.npy", x/self%units%length)
      call save_npy(folder//"pos_y.npy", y/self%units%length)
      call save_npy(folder//"pos_z.npy", z/self%units%length)
      call save_npy(folder//"m_phi.npy", phi)
      call save_npy(folder//"m_theta.npy", theta)
      call save_npy(folder//"site_type.npy", site_type)
      call save_npy(folder//"neigh.npy", neigh)
      call save_npy(folder//"conn_type.npy", conn_type)
      if (trim(self%mag_type) == "1Dspiral") then
         call save_npy(folder//"1Dspiralaxis.npy", self%axis)
         call save_npy(folder//"1Dspiralwavevector.npy", self%wavevector)
         !call save_npy(folder//"1Dspiralconeangle.npy", [self%cone_angle])
      endif
      call save_npy(folder//"lattice.npy", &
                    self%lattice/self%units%length)
      call save_npy(folder//"rez_lattice.npy", &
                    self%rez_lattice/self%units%inv_length)
   end subroutine save_unit_cell

   subroutine setup_single_hex(self)
      implicit none
      class(unit_cell), intent(inout)   :: self
      real(8)                           :: base_len
      real(8), dimension(3, 3)           :: base_vecs

      self%atoms(1) = init_ferro_z((/0d0, 0d0, 0d0/))
      allocate (self%atoms(1)%neigh_idx(3))
      allocate (self%atoms(1)%neigh_conn(3, 3))

      self%atoms(1)%neigh_idx = (/1, 1, 1/)

      base_len = self%lattice_constant
      base_vecs(1, :) = (/1d0, 0d0, 0d0/)
      base_vecs(2, :) = (/0.5d0, sin(60d0/180d0*PI), 0d0/)
      base_vecs(3, :) = (/-0.5d0, sin(60d0/180d0*PI), 0d0/)
      base_vecs = base_vecs*base_len
      self%atoms(1)%neigh_conn = base_vecs

      self%lattice(:, 1) = base_vecs(1, 1:2)
      self%lattice(:, 2) = base_vecs(2, 1:2)
   end subroutine

   subroutine setup_square(self)
      implicit none
      class(unit_cell), intent(inout)  :: self
      integer                          :: i, j, cnt
      real(8) :: pos(3)

      cnt = 1
      do i = 0, self%atom_per_dim - 1
         do j = 0, self%atom_per_dim - 1
            pos = [i, j, 0]*self%lattice_constant
            self%atoms(cnt) = init_ferro_z(pos)
            cnt = cnt + 1
         enddo
      enddo

      self%lattice(:, 1) = (/1d0, 0d0/)*self%atom_per_dim &
                           *self%lattice_constant
      self%lattice(:, 2) = (/0d0, 1d0/)*self%atom_per_dim &
                           *self%lattice_constant

      !if we want a molecule, ensure that no wrap-around is found
      if (self%molecule) self%lattice = self%lattice*10d0
   end subroutine setup_square

   subroutine setup_honey(self, hexagon, site_type)
      implicit none
      class(unit_cell), intent(inout)  :: self
      real(8), intent(in)              :: hexagon(:, :)
      integer, intent(in)              :: site_type(:)
      real(8)                          :: pos(3)
      integer                          :: i

      do i = 1, size(hexagon, dim=1)
         pos = hexagon(i, :)
         self%atoms(i) = init_ferro_z(pos, site=site_type(i))
      enddo
   end subroutine setup_honey

   subroutine setup_gen_conn(self, conn_mtx, conn_type, transl_mtx)
      implicit none
      class(unit_cell)    :: self
      real(8), intent(in) :: conn_mtx(:, :) !> Matrix containing
      !> real-space connections. The first index inidcates
      !> the connection vector, the second the vector element
      integer(4), intent(in) :: conn_type(:)
      real(8), intent(in) :: transl_mtx(:, :) !> Matrix containing
      !> real-space translation vectors. Notation as in conn_mtx
      integer                 :: i, j, cnt, candidate, n_conn, n_found
      integer, allocatable :: neigh(:)
      real(8)  :: start_pos(3), conn(3)
      logical, allocatable :: found_conn(:)

      n_conn = size(conn_mtx, 1)
      if (n_conn /= size(conn_type)) then
         call error_msg("number of connections have to agree", abort=.True.)
      endif

      allocate (found_conn(n_conn))
      allocate (neigh(n_conn))
      !$omp parallel do default(shared) schedule(dynamic)&
      !$omp& private(start_pos, n_found, found_conn, cnt, neigh, j, conn, &
      !$omp& candidate)
      do i = 1, self%num_atoms
         start_pos = self%atoms(i)%pos

         n_found = 0
         found_conn = .False.
         cnt = 1
         neigh = -1

         do j = 1, n_conn
            conn = conn_mtx(j, :)
            candidate = self%gen_find_neigh(start_pos, conn, transl_mtx)

            if (candidate /= -1) then
               found_conn(j) = .True.
               neigh(j) = candidate
               n_found = n_found + 1
               cnt = cnt + 1
            endif
         enddo
         allocate (self%atoms(i)%neigh_idx(n_found))
         allocate (self%atoms(i)%neigh_conn(n_found, 3))
         allocate (self%atoms(i)%conn_type(n_found))
         cnt = 1
         do j = 1, n_conn
            if (found_conn(j)) then
               self%atoms(i)%neigh_conn(cnt, :) = conn_mtx(j, :)
               self%atoms(i)%conn_type(cnt) = conn_type(j)
               self%atoms(i)%neigh_idx(cnt) = neigh(j)
               cnt = cnt + 1
            endif
         enddo
      enddo
      deallocate (found_conn)
      deallocate (neigh)
   end subroutine setup_gen_conn

   function gen_find_neigh(self, start, conn, transl_mtx) result(neigh)
      implicit none
      class(unit_cell), intent(in) :: self
      real(8), intent(in)          :: start(3) !> starting position in RS
      real(8), intent(in)          :: conn(3) !> RS connection
      real(8), intent(in)          :: transl_mtx(:, :) !> RS translation vectors to next unit cell
      !> The vectors are save as columns in the matrix:
      !> The first index indicates the vector
      !> The second index indicates the element of the vector
      integer     :: neigh, idx, n_transl, i
      real(8)     :: new(3)

      n_transl = size(transl_mtx, dim=1)

      neigh = -1
      idx = self%in_cell(start, conn)
      if (idx /= -1) then
         neigh = idx
         return
      else
         do i = 1, n_transl
            new = conn + transl_mtx(i, :)
            idx = self%in_cell(start, new)
            if (idx /= -1) then
               neigh = idx
               return
            endif
            new = conn - transl_mtx(i, :)
            idx = self%in_cell(start, new)

            if (idx /= -1) then
               neigh = idx
               return
            endif
         enddo
      endif
   end function gen_find_neigh

   function already_in_red(pos, hex, till, transl_mtx) result(inside)
      implicit none
      real(8), intent(in)    :: pos(3), hex(:, :), transl_mtx(:, :)
      integer, intent(in) :: till
      logical                :: inside
      real(8)                :: new(3), delta_vec(3), delta
      integer                :: n_transl, i, trl

      n_transl = size(transl_mtx, dim=1)
      inside = .False.

      outer: do i = 1, till
         do trl = 1, n_transl

            new = pos + transl_mtx(trl, :)
            delta_vec = hex(i, :) - new
            delta = my_norm2(delta_vec)

            if (delta <= pos_eps) then
               inside = .True.
               exit outer
            endif

            new = pos - transl_mtx(trl, :)
            delta_vec = hex(i, :) - new
            delta = sqrt(dot_product(delta_vec, delta_vec))
            !delta     = my_norm2(delta_vec)

            if (delta <= pos_eps) then
               inside = .True.
               exit outer
            endif
         enddo
      enddo outer
   end function already_in_red

   function in_cell(self, start, conn) result(idx)
      ! if position is in hexagon the corresponding index is
      ! returned, else - 1
      implicit none
      class(unit_cell), intent(in)          :: self
      real(8), intent(in) :: start(3) !> RS start position
      real(8), intent(in) :: conn(3) !> RZ connection
      real(8) :: new(3), delta_vec(3), delta
      real(8), parameter :: repl_eps = 1d-8
      integer    :: idx
      integer    :: i

      new = start + conn

      idx = -1
      do i = 1, self%num_atoms
         delta_vec = new - self%atoms(i)%pos
         ! inlining could help
         delta = sqrt(dot_product(delta_vec, delta_vec))
         !delta     = my_norm2(delta_vec)

         !if(delta < self%eps) then
         if (delta < repl_eps) then
            idx = i
            exit
         endif
      enddo
   end function in_cell

   subroutine calc_num_atoms_full_honey(n, n_atm, side)
      implicit none
      integer, intent(in)  :: n
      integer, intent(out) :: n_atm, side
      integer                 :: i

      side = 2
      n_atm = 0

      do i = 1, n
         if (mod(i, 3) == 0) then
            n_atm = n_atm + side
            side = side + 2
         else
            n_atm = n_atm + side - 1
         endif
      enddo
      n_atm = n_atm*6
   end subroutine calc_num_atoms_full_honey

   function calc_num_atoms_non_red_honey(n) result(n_atm)
      implicit none
      integer, intent(in)   :: n
      integer                  :: inner, next_side, n_atm

      call calc_num_atoms_full_honey(n - 1, inner, next_side)

      if (mod(n, 3) == 0) then
         n_atm = inner + 3*next_side
      else
         n_atm = inner + 3*next_side - 4
      endif
   end function calc_num_atoms_non_red_honey

   function calc_num_atoms_line_honey(n) result(n_atm)
      implicit none
      integer, intent(in)   :: n
      integer                  :: n_atm

      n_atm = 4*n

   end function calc_num_atoms_line_honey

   function in_hexagon(pos, a) result(inside)
      implicit none
      real(8), intent(in)         :: pos(3), a
      real(8)                     :: m, b
      logical                     :: inside

      m = -tan(deg_30)
      inside = .True.

      if (abs(pos(1)) > a*(cos(deg_30) + pos_eps)) then
         inside = .False.
      endif

      b = a*(1d0 + pos_eps)
      if (.not. (abs(pos(2)) <= m*abs(pos(1)) + b)) then
         inside = .False.
      endif
   end function in_hexagon

   subroutine gen_hexa_grid(l, origin, max_ind, grid)
      implicit none
      real(8), intent(in)     :: l, origin(3)
      integer, intent(in)  :: max_ind
      real(8), allocatable    :: grid(:, :)
      real(8)                 :: v1(3), v2(3)
      integer                 :: cnt, i, j

      if (.not. allocated(grid)) then
         allocate (grid((2*max_ind + 1)**2, 3))
      endif

      v1 = [l, 0d0, 0d0]
      v2 = [0.5d0*l, cos(deg_30)*l, 0d0]

      cnt = 1
      do i = -max_ind, max_ind
         do j = -max_ind, max_ind
            grid(cnt, :) = origin + i*v1 + j*v2
            cnt = cnt + 1
         enddo
      enddo
   end subroutine gen_hexa_grid

   subroutine gen_honey_grid(a, max_ind, grid)
      implicit none
      real(8), intent(in)        :: a
      integer, intent(in)     :: max_ind
      real(8), allocatable       :: grid(:, :), tmp(:, :)
      real(8)                    :: l, origin(3)
      integer                    :: n

      n = 2*max_ind + 1
      l = 2d0*cos(deg_30)*a

      if (.not. allocated(grid)) then
         allocate (grid(2*n**2, 3))
      endif

      origin = [0d0, a, 0d0]
      call gen_hexa_grid(l, origin, max_ind, tmp)
      grid(1:n**2, :) = tmp

      origin = [0d0, -a, 0d0]
      call gen_hexa_grid(l, origin, max_ind, tmp)
      grid(n**2 + 1:2*n**2, :) = tmp
   end subroutine gen_honey_grid

   function get_num_atoms(self) result(num)
      implicit none
      class(unit_cell), intent(in) :: self
      integer    :: num
      num = self%num_atoms
   end function get_num_atoms

   function get_atoms(self) result(ret)
      implicit none
      class(unit_cell), intent(in)            :: self
      type(atom), dimension(:), allocatable   :: ret

      ret = self%atoms
   end function get_atoms

   function calc_area(self) result(area)
      implicit none
      class(unit_cell), intent(in)            :: self
      real(8)                                 :: area, base_len_uc

      base_len_uc = self%lattice_constant*self%atom_per_dim
      if (trim(self%uc_type) == "honey_2d") then
         area = 1.5d0*sqrt(3d0)*base_len_uc**2
      else
         area = base_len_uc**2
      endif
   end function calc_area

   function rot_z_deg(deg) result(rot)
      implicit none
      real(8), intent(in)        :: deg
      real(8), dimension(3, 3)    :: rot
      real(8)                    :: bog

      bog = deg*PI/180.0d0

      rot = 0.0d0
      rot(1, 1) = cos(bog)
      rot(1, 2) = sin(bog)
      rot(2, 1) = -sin(bog)
      rot(2, 2) = cos(bog)
      rot(3, 3) = 1.0d0
   end function rot_z_deg

   function R_mtx(theta, vec) result(R)
      implicit none
      real(8), intent(in)    :: theta!> rotation angle
      real(8), intent(in)    :: vec(3) !> vector to rotate AROUND
      real(8), parameter     :: Iden(3, 3) &
                                = reshape((/1, 0, 0, & !this works only for
                                            0, 1, 0, & !symm matrcies
                                            0, 0, 1/), (/3, 3/)) !fort-order...
      real(8)  :: R(3, 3), u(3, 1), u_x_u(3, 3), u_x(3, 3)

      u(:, 1) = vec/my_norm2(vec)

      u_x_u = matmul(u, transpose(u))

      u_x = 0d0
      u_x(2, 1) = u(3, 1)
      u_x(3, 1) = -u(2, 1)
      u_x(1, 2) = -u(3, 1)
      u_x(3, 2) = u(1, 1)
      u_x(1, 3) = u(2, 1)
      u_x(2, 3) = -u(1, 1)

      R = cos(theta)*Iden &
          + sin(theta)*u_x &
          + (1 - cos(theta))*u_x_u

   end function R_mtx

   function n_times_phi(x, n) result(y)
      implicit none
      real(8), intent(in)   :: x(3)
      integer, intent(in) :: n
      real(8)               :: y(3), r, phi, theta

      r = my_norm2(x)
      theta = acos(x(3)/r)
      phi = atan2(x(2), x(1))

      phi = n*phi

      y(1) = r*sin(theta)*cos(phi)
      y(2) = r*sin(theta)*sin(phi)
      y(3) = r*cos(theta)
   end function n_times_phi

   subroutine run_tests(self)
      use mpi
      implicit none
      class(unit_cell), intent(in)   :: self
      integer                        :: i, ierr
      logical                        :: passed, tmp

      passed = .True.
      !compare perturbation logical
      if (self%me == root) tmp = self%pert_log
      call MPI_Bcast(tmp, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
      if (tmp .NEQV. self%pert_log) then
         call error_msg("pert_log doesn't match", abort=.True.)
         passed = .False.
      endif
      do i = 1, self%num_atoms
         passed = passed .and. self%atoms(i)%compare_to_root()
      enddo

      if (passed) then
         call error_msg("Atom test passed", p_color=c_green, abort=.False.)
      endif
   end subroutine run_tests

end module

