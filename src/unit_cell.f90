module Class_unit_cell
    implicit none
    private

    type, public        :: unit_cell
        ! number of non-redundant atoms pre unit cell
        integer*4       :: num_atoms
        real*8          :: phi, theta
        integer*4(3)    :: neighbours
        real*8(2)       :: pos
    contains
    end type unit_cell
contains

