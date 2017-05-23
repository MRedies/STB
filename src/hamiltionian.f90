module Class_hamiltionian
    use m_config
    use output
    use Class_unit_cell
    use m_npy
    USE, INTRINSIC :: IEEE_ARITHMETIC
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
        procedure :: set_Stoner          => set_Stoner
        procedure :: setup_Stoner_mtx    => setup_Stoner_mtx
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
       
        H =  0d0

        call self%set_EigenE(H)
        call self%set_hopping(k,H)
        call self%set_Stoner(H)
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

    subroutine set_Stoner(self,H)
        implicit none
        class(hamil), intent(in)   :: self
        complex(8), intent(inout)  :: H(:,:)
        complex(8)                 :: S(2,2) !> Stonermatrix
        integer(4)                 :: i, i_up, i_dw

        do i =  1,self%UC%num_atoms
            i_up =  i
            i_dw =  i +  self%UC%num_atoms 

            call self%setup_Stoner_mtx(i,S)

            H(i_up,i_up) = H(i_up, i_up) +  S(1,1)
            H(i_up,i_dw) = H(i_up, i_dw) +  S(1,2)
            H(i_dw,i_up) = H(i_dw, i_up) +  S(2,1)
            H(i_dw,i_dw) = H(i_dw, i_dw) +  S(2,2)
        enddo

    end subroutine set_Stoner

    subroutine setup_Stoner_mtx(self,i,S)
        implicit none
        class(hamil), intent(in) :: self
        integer(4), intent(in)   :: i
        complex(8), intent(inout):: S(2,2)
        real(8)                  :: m(3), fac
        
        m = self%UC%atoms(i)%get_m_cart()
        fac =  - 0.5d0 *  self%I
        
        S = fac * ( m(1) * sigma_x &
                  + m(2) * sigma_y &
                  + m(3) * sigma_z)
    end subroutine setup_Stoner_mtx

    subroutine set_EigenE(self,H)
        implicit none
        class(hamil), intent(in) :: self
        complex(8), intent(inout):: H(:,:)
        integer(4) :: i

        do i =  1,size(H,dim=1)
            H(i,i) = H(i,i) + self%E_s
        enddo
    end subroutine set_EigenE 

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
        integer(4), allocatable           :: IWORK(:)
        character(len=20) :: filename
        N =  2 * self%UC%num_atoms
        LWMAX =  10*N
        allocate(eig_val(size(k_list, 2), N))
        allocate(H(N,N))
        allocate(RWORK(LWMAX))
        allocate(IWORK(LWMAX))
        allocate(WORK(LWMAX))
        allocate(tmp_out(N))
        
        do i = 1,size(k_list,2)
            k =  k_list(:,i)
            call self%setup_H(k, H)
            
            !call zheev('V', 'U', N, H, N, tmp_out, WORK, LWMAX, RWORK, info)
            call zheevd('N', 'U', N, H, N, tmp_out, WORK, LWMAX, &
                                RWORK, LWMAX, IWORK, LWMAX, info)
            if( info == 0) then
                eig_val(i,:) =  tmp_out 
            else
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
