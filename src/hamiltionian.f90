module Class_hamiltionian
    use m_config
    use Class_unit_cell
    implicit none
    
    type hamil
        real(8)         :: E_s !> onsite eigenenergy
        real(8)         :: t_nn !> nearest neighbour hopping
        real(8)         :: I !> stoner parameter
        type(unit_cell) :: UC !> unit cell
    contains
        procedure :: setup_H             => setup_H
        procedure :: calc_eigenvalues    => calc_eigenvalues
        procedure :: set_EigenE          => set_EigenE
        procedure :: set_hopping         => set_hopping
    end type hamil

contains
    subroutine setup_H(self,k,H)
        implicit none
        class(hamil)        :: self
        real(8), intent(in) :: k(3)
        complex(8), intent(inout) :: H(:,:)
        
        if(k(3) /= 0d0) then
            write (*,*) "K_z is non-zero. Abort."
            stop
        endif
       
        call reset_H(H)
        call self%set_EigenE(H)
        call self%set_hopping(k,H)
    end subroutine setup_H

    function init_hamil(cfg) result(ret)
        implicit none
        type(CFG_t)    :: cfg
        type(hamil)    :: ret
        real(8)        :: tmp

        ret%UC =  init_unit(cfg)
        call CFG_get(cfg, "hamil%E_s", tmp)
        ret%E_s =  tmp * get_unit_conv("energy", cfg)
        
        call CFG_get(cfg, "hamil%t_nn", tmp)
        ret%t_nn =  tmp * get_unit_conv("energy", cfg)
        
        call CFG_get(cfg, "hamil%I", tmp)
        ret%I =  tmp * get_unit_conv("energy", cfg)
    end function init_hamil

    subroutine set_EigenE(self,H)
        implicit none
        class(hamil), intent(in) :: self
        complex(8), intent(inout):: H(:,:)
        integer(4) :: i

        do i =  1,size(H,dim=1)
            H(i,i) = H(i,i) + self%E_s
        enddo
    end subroutine set_EigenE 

    subroutine reset_H(H)
        implicit none
        complex(8), intent(inout):: H(:,:)

        H =  0d0
    end subroutine reset_H
    
    subroutine set_hopping(self,k, H)
        implicit none
        class(hamil), intent(in)          :: self 
        real(8), intent(in)               :: k(3)
        complex(8), intent(inout)         :: H(:,:)
        integer(4)                        :: i, i_up, i_d, j,&
                                             j_up, j_d, conn
        real(8)                           :: k_dot_r
        complex(8)                        :: new

        ! Spin up
        do i = 1,self%UC%num_atoms
            do conn = 1,self%UC%atoms(i)%n_neigh
                j =  self%UC%atoms(i)%neigh_idx(conn)
                k_dot_r =  dot_product(k, self%UC%atoms(i)%neigh_conn(conn,:))
                
                new = exp(i_unit * k_dot_r) * self%UC%atoms(i)%hopping(conn)
                H(i,j) =  H(i,j) + new
                H(j,i) =  H(j,i) + conjg(new)
            enddo
        enddo

        ! Spin down
        do i = 1,self%UC%num_atoms
            i_d =  i + self%UC%num_atoms
            do conn = 1,self%UC%atoms(i)%n_neigh

                j      = self%UC%atoms(i)%neigh_idx(conn)
                j_d = j + self%UC%num_atoms
                
                k_dot_r =  dot_product(k, self%UC%atoms(i)%neigh_conn(conn,:))
                new =  exp(i_unit *  k_dot_r) * self%UC%atoms(i)%hopping(conn)
                H(i_d, j_d) = H(i_d, j_d) + new
                H(j_d, i_d) = H(j_d, i_d) + conjg(new)
            enddo
        enddo

    end subroutine set_hopping
    
    Subroutine  calc_eigenvalues(self, k_list, eig_val)
        Implicit None
        class(hamil)                      :: self
        real(8), intent(in)               :: k_list(:,:)
        real(8), allocatable,intent(out)  :: eig_val(:,:)
        real(8)                           :: k(3)
        complex(8), allocatable           :: H(:,:)
        integer(4) :: i, N, LWMAX, info
        real(8), allocatable              :: RWORK(:), tmp_out(:)
        complex(8), allocatable           :: WORK(:) 
        character(len=20) :: filename
        N =  2 * self%UC%num_atoms
        LWMAX =  10*N
        allocate(eig_val(size(k_list, 2), N))
        allocate(H(N,N))
        allocate(RWORK(LWMAX))
        allocate(WORK(LWMAX))
        allocate(tmp_out(N))
        
        do i = 1,size(k_list,2)
            k =  k_list(:,i)
            call self%setup_H(k, H)
            
            call zheev('V', 'U', N, H, N, tmp_out, WORK, LWMAX, RWORK, info)
            eig_val(i,:) =  tmp_out 
            if( info /= 0) then
                write (*,*) "ZHEEV failed: ", info
                stop
            endif
        enddo

        deallocate(H)
        deallocate(RWORK)
        deallocate(WORK)
        deallocate(tmp_out)
    End Subroutine calc_eigenvalues

end module
