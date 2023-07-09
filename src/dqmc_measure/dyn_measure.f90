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
    
    !!!allocate( g00do(ndim,ndim), g0tdo(ndim,ndim), gt0do(ndim,ndim), gttdo(ndim,ndim) )
    !!!do j = 1,ndim
    !!!   xj = 1.d0
    !!!   if (list(j,2) == 1 ) xj= -1.d0
    !!!   do i = 1,ndim
    !!!      xi = 1.d0
    !!!      if (list(i,2) == 1 ) xi = -1.d0
    !!!      gt0do(i,j)  = dcmplx(xi*xj,0.d0)*dconjg( g0tup(j,i) )
    !!!      g0tdo(i,j)  = dcmplx(xi*xj,0.d0)*dconjg( gt0up(j,i) )
    !!!      gttdo(i,j)  = dcmplx(xi*xj,0.d0)* ( Imat(i,j) - dconjg( gttup(j,i) ) )
    !!!      g00do(i,j)  = dcmplx(xi*xj,0.d0)* ( Imat(i,j) - dconjg( g00up(j,i) ) )
    !!!   enddo
    !!!enddo

    if( dble(phase) < 0.d0 ) then
        sgn = -1.d0
    else
        sgn = 1.d0
    end if
    zphi = dcmplx(sgn/dble(phase), 0.d0) * phase

    gtau0_orb1_tau = czero
    zspsm_orb1_tau = czero
    znn_orb1_tau = czero
    do j = 1, latt%nsites
        nu_j = latt%list(j,1)
        no_j = latt%list(j,2)
        do i = 1, latt%nsites
            nu_i = latt%list(i,1)
            no_i = latt%list(i,2)
            imj  = latt%imj(nu_i,nu_j)
!            gtau0_orb1_tau(imj,no_i,no_j,nt) = gtau0_orb1_tau(imj,no_i,no_j,nt) + gt0up%orb1(i,j)
!            zspsm_orb1_tau(imj,no_i,no_j,nt) = zspsm_orb1_tau(imj,no_i,no_j,nt) - g0tup%orb1(j,i)*gt0up%orb1(i,j)*dcmplx(1.d0-1.d0/dble(nflr*nflr), 0.d0)
            znn_orb1_tau(imj,no_i,no_j,nt) = znn_orb1_tau(imj,no_i,no_j,nt) + dcmplx( (0.5d0-dble(gttup%orb1(i,i)))*(0.5d0-dble(g00up%orb1(j,j))) ) &
                               - dcmplx(1.d0/dble(nflr),0.d0)*g0tup%orb1(j,i)*gt0up%orb1(i,j)
        end do
    end do
!    gtau0_orb1_tau_bin = gtau0_orb1_tau_bin + gtau0_orb1_tau*zphi
!    zspsm_orb1_tau_bin = zspsm_orb1_tau_bin + zspsm_orb1_tau*zphi
    znn_orb1_tau_bin = znn_orb1_tau_bin + znn_orb1_tau*zphi

!    zbb_orb1_tau = czero
!    zb_orb1_tau = czero
!    do nu_j = 1, lq
!#IFDEF HONEYCOMB
!        j_0 = 2*nu_j - 1 ! A site
!        j_n = latt%nnlist(j_0,1) ! B site
!#ENDIF
!#IFDEF SQUARE
!        j_0 = nu_j  ! A site
!        j_n = latt%nnlist(j_0,1) ! B site
!#ENDIF
!        zj = expar(j_0,1,xmag,flux_x,flux_y,dimer)
!        zb_orb1_tau(nu_j,nt) = zb_orb1_tau(nu_j,nt) + (zj*(-gttup%orb1(j_n,j_0)) + dconjg(zj)*(-gttup%orb1(j_0,j_n)))
!        do nu_i = 1, lq
!#IFDEF HONEYCOMB
!            i_0 = 2*nu_i - 1 ! A site
!            i_n = latt%nnlist(i_0,1) ! B site
!#ENDIF
!#IFDEF SQUARE
!            i_0 = nu_i ! A site
!            i_n = latt%nnlist(i_0,1) ! B site
!#ENDIF
!            zi = expar(i_0,1,xmag,flux_x,flux_y,dimer)
!            imj = latt%imj(nu_i,nu_j)
!            zbb_orb1_tau(imj,nt) = zbb_orb1_tau(imj,nt) + (zj*(-g00up%orb1(j_n,j_0)) + dconjg(zj)*(-g00up%orb1(j_0,j_n)))*(zi*(-gttup%orb1(i_n,i_0)) + dconjg(zi)*(-gttup%orb1(i_0,i_n))) &
!                                          + ( -zi*zj*g0tup%orb1(j_n,i_0)*gt0up%orb1(i_n,j_0) - dconjg(zi)*dconjg(zj)*g0tup%orb1(j_0,i_n)*gt0up%orb1(i_0,j_n) &
!                                          -   zi*dconjg(zj)*g0tup%orb1(j_0,i_0)*gt0up%orb1(i_n,j_n) - dconjg(zi)*zj*g0tup%orb1(j_n,i_n)*gt0up%orb1(i_0,j_0) )*dcmplx(1.d0/dble(nflr),0.d0)
!        end do
!    end do
    zbb_orb1_tau_bin = zbb_orb1_tau_bin + zbb_orb1_tau*zphi
    zb_orb1_tau_bin = zb_orb1_tau_bin + zb_orb1_tau*zphi
  end subroutine dyn_measure
