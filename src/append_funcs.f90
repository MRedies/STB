module Class_append_funcs
    use mpi
    implicit none
    
    contains
    subroutine add_to_arr1D(list, element)
        implicit none
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
  
    end subroutine add_to_arr1D
  
    subroutine add_to_arr2D(list, element)
        implicit none
        integer                             :: i
        integer,allocatable                 :: isize(:)
        real(8), allocatable, intent(in)    :: element
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
            write(*,*) "Append shape does agree!"
        endif
        deallocate(list)
        call move_alloc(clist, list)
  
    else
         allocate(list(1,1))
         list(1,1) = element
    end if
  
  end subroutine add_to_arr2D
end module