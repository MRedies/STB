module Class_append_funcs
    use mpi
    use m_config
    use m_npy
    use Class_units
    use stdlib_io_npy, only: npy_load
    implicit none
    type collect_quantities
        real(8), allocatable ::  int_DOS_collect(:,:)
        real(8), allocatable ::  DOS_collect(:,:)
        real(8), allocatable ::  up_collect(:,:)
        real(8), allocatable ::  down_collect(:,:)
        real(8), allocatable ::  spins_collect(:,:)
        integer, allocatable ::  sample_idx(:)
        integer              ::  me,me_sample,color
        character(len=300)   :: prefix
        type(units)          :: units
    contains
        procedure :: add_DOS_collect => add_DOS_collect
        procedure :: save_DOS_collect => save_DOS_collect
        procedure :: add_spins_collect => add_spins_collect
        procedure :: add_sample_idx => add_sample_idx
        procedure :: save_sample_idx => save_sample_idx
        procedure :: save_spins_collect => save_spins_collect
        procedure :: add_to_arr1D_int => add_to_arr1D_int
        procedure :: add_to_arr1D_real => add_to_arr1D_real
        procedure :: add_to_arr2D_real => add_to_arr2D_real
        procedure :: add_2D_to_arr2D_real => add_2D_to_arr2D_real
    end type collect_quantities
    
    contains
        function init_collect_quantities(cfg,prefix,sample_comm,color) result(self)
            use mpi
            implicit none
            type(collect_quantities) :: self
            type(CFG_t)              :: cfg
            integer                  :: ierr,sample_comm,color
            character(len=300)       :: prefix
    
            call MPI_Comm_rank(MPI_COMM_WORLD, self%me, ierr)
            call MPI_Comm_rank(sample_comm, self%me_sample, ierr)
            self%units = init_units(cfg, self%me)
            self%prefix = trim(prefix)
            self%color = color
        end function init_collect_quantities

        subroutine add_DOS_collect(self, DOS, up, down, int_DOS)
            use mpi
            implicit none
            class(collect_quantities)           :: self
            real(8), intent(in)                 :: DOS(:), up(:), down(:), int_DOS(:)
    
            if(self%me_sample==root) then
                if(.NOT. allocated(self%DOS_collect)) then
                    allocate(self%DOS_collect(1,size(DOS)))
                    self%DOS_collect(1,:) = DOS
                else
                    call self%add_to_arr2D_real(self%DOS_collect,DOS)
                endif
                if(.NOT. allocated(self%up_collect)) then
                    allocate(self%up_collect(1,size(up)))
                    self%up_collect(1,:) = up
                else
                    call self%add_to_arr2D_real(self%up_collect,up)
                endif
                if(.NOT. allocated(self%down_collect)) then
                    allocate(self%down_collect(1,size(down)))
                    self%down_collect(1,:) = down
                else
                    call self%add_to_arr2D_real(self%down_collect,down)
                endif
                if(.NOT. allocated(self%int_DOS_collect)) then
                    allocate(self%int_DOS_collect(1,size(int_DOS)))
                    self%int_DOS_collect(1,:) = int_DOS
                else
                    call self%add_to_arr2D_real(self%int_DOS_collect,int_DOS)
                endif
            endif
        end subroutine
        
        subroutine add_spins_collect(self, spins)
            use mpi
            implicit none
            class(collect_quantities)           :: self
            integer                             :: i
            integer,allocatable                 :: isize(:)
            real(8), intent(in)                 :: spins(:,:)
            allocate(isize(2))
            if(self%me_sample==root) then
                if(.NOT. allocated(self%spins_collect)) then
                    isize = shape(spins)
                    allocate(self%spins_collect(isize(1),isize(2)))
                    do i=1,isize(1)
                        self%spins_collect(i,:) = spins(i,:)
                    enddo
                else
                    call self%add_2D_to_arr2D_real(self%spins_collect,spins)
                endif
            endif
        end subroutine

        subroutine add_sample_idx(self, idx)
            use mpi
            implicit none
            class(collect_quantities)           :: self
            integer, intent(in)                 :: idx

            if(self%me_sample==root) then
                if(.NOT. allocated(self%sample_idx)) then
                    allocate(self%sample_idx(1))
                    self%sample_idx = idx
                else
                    call self%add_to_arr1D_int(self%sample_idx,idx)
                endif
            endif
        end subroutine
        
        subroutine save_DOS_collect(self)
            use mpi
            implicit none
            class(collect_quantities)           :: self
            character(len=300)                  :: filename
            
    
            if(self%me_sample ==  root) then
                write (filename,  "(A,I0.6,A)") "DOS_collect=", self%color,".npy"
                call save_npy(trim(self%prefix) //  trim(filename), self%DOS_collect * self%units%energy)
                write (filename,  "(A,I0.6,A)") "int_DOS_collect=", self%color,".npy"
                call save_npy(trim(self%prefix) //  trim(filename), self%int_DOS_collect * self%units%energy)
                write (filename,  "(A,I0.6,A)") "up_collect=", self%color,".npy"
                call save_npy(trim(self%prefix) //  trim(filename), self%up_collect * self%units%energy)
                write (filename,  "(A,I0.6,A)") "down_collect=", self%color,".npy"
                call save_npy(trim(self%prefix) //  trim(filename), self%down_collect * self%units%energy)
            endif
        end subroutine

        subroutine save_spins_collect(self)
            use mpi
            implicit none
            class(collect_quantities)           :: self
            character(len=300)                  :: filename
            
    
            if(self%me_sample ==  root) then
                write (filename,  "(A,I0.6,A)") "spins_collect=", self%color,".npy"
                call save_npy(trim(self%prefix) //  trim(filename), self%spins_collect)
            endif
        end subroutine

        subroutine save_sample_idx(self)
            use mpi
            implicit none
            class(collect_quantities)           :: self
            character(len=300)                  :: filename
            
    
            if(self%me_sample ==  root) then
                write (filename,  "(A,I0.6,A)") "sample_idx_collect=", self%color,".npy"
                call save_npy(trim(self%prefix) //  trim(filename), self%sample_idx)
            endif
        end subroutine

        subroutine add_to_arr1D_int(self,list, element)
            implicit none
            class(collect_quantities)           :: self
            integer                             :: i,isize
            integer, intent(in)                 :: element
            integer, allocatable, intent(inout) :: list(:)
            integer, allocatable                :: clist(:)
    
            if(allocated(list)) then
                isize = size(list)
                allocate(clist(isize+1))
                do i=1,isize          
                clist(i) = list(i)
                end do
                clist(isize+1) = element
    
                deallocate(list)
                call move_alloc(clist, list)
    
            else
                allocate(list(1))
                list(1) = element
            end if
    
        end subroutine add_to_arr1D_int

        subroutine add_to_arr1D_real(self,list, element)
            implicit none
            class(collect_quantities)           :: self
            integer                             :: i,isize
            real(8), intent(in)                 :: element
            real(8), allocatable, intent(inout) :: list(:)
            real(8), allocatable                :: clist(:)
    
            if(allocated(list)) then
                isize = size(list)
                allocate(clist(isize+1))
                do i=1,isize          
                clist(i) = list(i)
                end do
                clist(isize+1) = element
    
                deallocate(list)
                call move_alloc(clist, list)
    
            else
                allocate(list(1))
                list(1) = element
            end if
    
        end subroutine add_to_arr1D_real
    
        subroutine add_to_arr2D_real(self,list, element)
            implicit none
            class(collect_quantities)           :: self
            integer                             :: i
            integer,allocatable                 :: isize(:)
            real(8)             , intent(in)    :: element(:)
            real(8), allocatable, intent(inout) :: list(:,:)
            real(8), allocatable                :: clist(:,:)
        
        allocate(isize(2))

        if(allocated(list)) then
            isize = shape(list)
            allocate(clist(isize(1)+1,isize(2)))
            do i=1,isize(1)
            clist(i,:) = list(i,:)
            end do
            if (size(element)==isize(2)) then
            clist(isize(1)+1,:) = element
            else
                write(*,*) "Append shape does not agree!"
            endif
            deallocate(list)
            call move_alloc(clist, list)
    
        else
            allocate(list(1,1))
            list(1,:) = element
        end if
    end subroutine add_to_arr2D_real

    subroutine add_2D_to_arr2D_real(self,list, element)
        implicit none
        class(collect_quantities)           :: self
        integer                             :: i
        integer,allocatable                 :: isize(:),esize(:)
        real(8)             , intent(in)    :: element(:,:)
        real(8), allocatable, intent(inout) :: list(:,:)
        real(8), allocatable                :: clist(:,:)
    
    allocate(isize(2))
    allocate(esize(2))

    if(allocated(list)) then
        isize = shape(list)
        esize = shape(element)
        allocate(clist(isize(1)+esize(1),isize(2)))
        do i=1,isize(1)
        clist(i,:) = list(i,:)
        end do
        if (esize(2)==isize(2)) then
            do i=1,esize(1)
                clist(isize(1)+i,:) = element(i,:)
            enddo
        else
            write(*,*) "Append shape does not agree!"
        endif
        deallocate(list)
        call move_alloc(clist, list)

    else
        esize = shape(element)
        allocate(list(esize(1),esize(2)))
        do i=1,esize(1)
            list(i,:) = element(i,:)
        enddo
    end if
    
    end subroutine add_2D_to_arr2D_real
end module