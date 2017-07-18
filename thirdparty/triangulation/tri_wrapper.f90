subroutine run_triang(k_pts, ret_elem)
    use output
    use iso_c_binding, only : C_INT, C_DOUBLE
    implicit none
    real(8), intent(in)              :: k_pts(:,:)
    integer(4), allocatable          :: ret_elem(:,:)
    integer(kind=C_INT)              :: n_nodes
    real(kind=C_DOUBLE), allocatable :: x(:), y(:)
    integer(kind=C_INT), allocatable :: elem(:)
    integer(kind=C_INT)              :: n_elem, ierr
    integer, allocatable ::  resh(:,:)

    interface
        subroutine my_tri (x, y, n_nodes, elem, n_elem, ierr ) bind ( C, name = "my_tri" )
            use iso_c_binding, only : C_INT, C_DOUBLE 
            real(kind = C_DOUBLE) :: x(*), y(*)
            integer(kind = C_INT) :: n_nodes, n_elem, ierr, elem(*)

        end subroutine my_tri
    end interface

    n_nodes =  size(k_pts, 2)
    allocate(x(n_nodes))
    allocate(y(n_nodes))
    n_elem =  2 * n_nodes 
    allocate(elem(3*n_elem))
    allocate(resh(n_elem, 3))

    x =  k_pts(1,:)
    y =  k_pts(2,:)

    call my_tri(x,y, n_nodes, elem, n_elem, ierr)
    if(ierr /= 0) then
        write (*,*) "Error: Triangulation failed", ierr
        stop 0
    endif


    resh =  reshape(elem, [2*n_nodes,3])
   
    if(allocated(ret_elem))then
        if(size(ret_elem,1) /= n_elem) deallocate(ret_elem)
    endif
    if(.not. allocated(ret_elem)) allocate(ret_elem(n_elem,3))

    ret_elem = resh(1:n_elem, :)

    deallocate(x)
    deallocate(y)
    deallocate(elem)
    deallocate(resh)
end subroutine run_triang
