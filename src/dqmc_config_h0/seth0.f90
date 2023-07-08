  subroutine seth0(h0mat, rt1, mu1, xmag1, flux_x1, flux_y1, dimer1)
    use model_para, only: latt, ndim, dtau
    implicit none
    
    real(dp), intent(in) :: rt1, mu1, xmag1, flux_x1, flux_y1, dimer1
    complex(dp), dimension(latt%nsites,latt%nsites), intent(out) :: h0mat
    
    ! local
    integer :: i, i_0, i_n, nf
    complex(dp) :: z1

    h0mat = czero

    ! nearest hopping
    do i = 1, latt%nn_lf
      i_0 = latt%nnlf_list(i) ! A site
      do nf = 1, latt%nn_nf
        i_n = latt%nnlist(i_0,nf)
        z1 = dcmplx(-rt1,0.d0)*expar(i_0,nf,xmag1,flux_x1,flux_y1,dimer1)
        h0mat(i_0,i_n)  = h0mat(i_0,i_n) + z1
        h0mat(i_n,i_0)  = h0mat(i_n,i_0) + dconjg(z1)
      end do
    end do
    ! add chemical potential
    do i = 1, latt%nsites
        h0mat(i,i) = h0mat(i,i) - dcmplx(mu1, 0.d0)
    end do
  end subroutine seth0
