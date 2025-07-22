  subroutine dyn_measure(nt,gt0up,g0tup,gttup,g00up)
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF  
    implicit none
    
    !g0t_(i,j) =-<c^{dagger}(tau)_j c_i>
    !gt0_(i,j) = <c_i(tau) c^{dagger}_j>

    !	arguments.
    type(gfunc), intent(in) :: g00up, g0tup, gt0up, gttup
    integer, intent(in) :: nt
    
    !       local
    integer :: i, j, no_i, no_j, nu_i, nu_j, imj, i_0, j_0, i_n, j_n
    complex(dp) :: zi, zj

    !g0t_(i,j) = -<c^{dagger}(tau)_j c_i>
    !gt0_(i,j) =  <c_i(tau) c^{dagger}_j>
    ! Attention to this definition!!!
    ! As in projector version code, I define g0t_(i,j) = <c^{dagger}(tau)_j c_i>, has NO sgn in front
    !
    ! How to use it?
    ! for example, we calculate
    ! <S^+_i(tau) S^-_j(0)> = < cup^+_i(tau) cdn_i(tau) cdn^+_j(0) cup_j(0) >
    !                       = < cup^+_i(tau) cup_j(0) > < cdn_i(tau) cdn^+_j(0) >
    !                       = -g0tup(j,i) gt0dn(i,j)

    if( dble(phase) < 0.d0 ) then
        sgn = -1.d0
    else
        sgn = 1.d0
    end if
    zphi = dcmplx(sgn/dble(phase), 0.d0) * phase

    gtau0_tau = czero
    zspsm_tau = czero
    znn_tau = czero
    do j = 1, latt%nsites
        nu_j = latt%list(j,1)
        no_j = latt%list(j,2)
        do i = 1, latt%nsites
            nu_i = latt%list(i,1)
            no_i = latt%list(i,2)
            imj  = latt%imj(nu_i,nu_j)
!            gtau0_tau(imj,no_i,no_j,nt) = gtau0_tau(imj,no_i,no_j,nt) + gt0up%blk1(i,j)
!            zspsm_tau(imj,no_i,no_j,nt) = zspsm_tau(imj,no_i,no_j,nt) - g0tup%blk1(j,i)*gt0up%blk1(i,j)*dcmplx(1.d0-1.d0/dble(nflr*nflr), 0.d0)
            znn_tau(imj,no_i,no_j,nt) = znn_tau(imj,no_i,no_j,nt) + dcmplx( (0.5d0-dble(gttup%blk1(i,i)))*(0.5d0-dble(g00up%blk1(j,j))) ) &
                               - dcmplx(1.d0/dble(nflr),0.d0)*g0tup%blk1(j,i)*gt0up%blk1(i,j)
        end do
    end do
!    gtau0_tau_bin = gtau0_tau_bin + gtau0_tau*zphi
!    zspsm_tau_bin = zspsm_tau_bin + zspsm_tau*zphi
    znn_tau_bin(:,:,:,nt) = znn_tau_bin(:,:,:,nt) + znn_tau(:,:,:,nt)*zphi

!    zbb_tau = czero
!    zb_tau = czero
!    do nu_j = 1, lq
!        j_0 = nu_j  ! for 1d chain
!        j_n = latt%nnlist(j_0,1) ! for 1d, nnlist should give site j+1
!        zj = expar(j_0,1,xmag,flux_x,flux_y,dimer)
!        zb_tau(nu_j,nt) = zb_tau(nu_j,nt) + (zj*(-gttup%blk1(j_n,j_0)) + dconjg(zj)*(-gttup%blk1(j_0,j_n)))
!        do nu_i = 1, lq
!            i_0 = nu_i  ! for 1d chain
!            i_n = latt%nnlist(i_0,1)
!            zi = expar(i_0,1,xmag,flux_x,flux_y,dimer)
!            imj = latt%imj(nu_i,nu_j)
!            zbb_tau(imj,nt) = zbb_tau(imj,nt) + (zj*(-g00up%blk1(j_n,j_0)) + dconjg(zj)*(-g00up%blk1(j_0,j_n)))*(zi*(-gttup%blk1(i_n,i_0)) + dconjg(zi)*(-gttup%blk1(i_0,i_n))) &
!                                          + ( -zi*zj*g0tup%blk1(j_n,i_0)*gt0up%blk1(i_n,j_0) - dconjg(zi)*dconjg(zj)*g0tup%blk1(j_0,i_n)*gt0up%blk1(i_0,j_n) &
!                                          -   zi*dconjg(zj)*g0tup%blk1(j_0,i_0)*gt0up%blk1(i_n,j_n) - dconjg(zi)*zj*g0tup%blk1(j_n,i_n)*gt0up%blk1(i_0,j_0) )*dcmplx(1.d0/dble(nflr),0.d0)
!        end do
!    end do
    zbb_tau_bin = zbb_tau_bin + zbb_tau*zphi
    zb_tau_bin = zb_tau_bin + zb_tau*zphi
  end subroutine dyn_measure
