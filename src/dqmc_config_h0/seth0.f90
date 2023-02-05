  subroutine seth0(h0mat, rt, mu)
    use model_para, only: latt, ndim, xmag, flux_x, flux_y, dimer, dtau
    implicit none
    
    real(dp), intent(in) :: rt, mu
    complex(dp), dimension(ndim,ndim), intent(out) :: h0mat
    
    ! local
    integer :: i, i_0, i_n, nf
    complex(dp) :: z1

    h0mat = czero

    ! nearest hopping
    do i = 1, latt%nn_lf
      i_0 = latt%nnlf_list(i) ! A site
      do nf = 1, latt%nn_nf
        i_n = latt%nnlist(i_0,nf)
        z1 = dcmplx(-rt,0.d0)*expar(i_0,nf,xmag,flux_x,flux_y,dimer)
        h0mat(i_0,i_n)  = h0mat(i_0,i_n) + z1
        h0mat(i_n,i_0)  = h0mat(i_n,i_0) + dconjg(z1)
      end do
    end do
    ! add chemical potential
    do i = 1, ndim
        h0mat(i,i) = h0mat(i,i) - dcmplx(mu, 0.d0)
    end do
  end subroutine seth0
