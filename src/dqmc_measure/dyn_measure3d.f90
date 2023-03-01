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
#IFDEF TIMING
    real(dp) :: starttime21, endtime21
#ENDIF
    
    !       local
    integer :: i, j, no_i, no_j, nu_i, nu_j, imj

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
#IFDEF TIMING
    starttime21 = omp_get_wtime()
#ENDIF

    if( dble(phase) < 0.d0 ) then
        sgn = -1.d0
    else
        sgn = 1.d0
    end if
    zphi = dcmplx(sgn/dble(phase), 0.d0) * phase

    gtau0 = czero
    zspsm_tau = czero
    znn_tau = czero
    do j = 1, ndim
        nu_j = latt%list(j,1)
        no_j = latt%list(j,2)
        do i = 1, ndim
            nu_i = latt%list(i,1)
            no_i = latt%list(i,2)
            imj  = latt%imj(nu_i,nu_j)
            gtau0(imj,nt) = gtau0(imj,nt) + gt0up%orb1(i,j)
            zspsm_tau(imj,nt) = zspsm_tau(imj,nt) - g0tup%orb1(j,i)*gt0up%orb1(i,j)*dcmplx(1.d0-1.d0/dble(nflr*nflr), 0.d0)
            znn_tau(imj,nt) = znn_tau(imj,nt) + dcmplx( (0.5d0-dble(gttup%orb1(i,i)))*(0.5d0-dble(g00up%orb1(j,j))) ) &
                               - dcmplx(1.d0/dble(nflr),0.d0)*g0tup%orb1(j,i)*gt0up%orb1(i,j)
        end do
    end do
    gtau0_bin = gtau0_bin + gtau0*zphi
    zspsm_tau_bin = zspsm_tau_bin + zspsm_tau*zphi
    znn_tau_bin = znn_tau_bin + znn_tau*zphi
#IFDEF TIMING
    endtime21 = omp_get_wtime()
    timecalculation(5)=timecalculation(5)+endtime21-starttime21
#ENDIF

  end subroutine dyn_measure
