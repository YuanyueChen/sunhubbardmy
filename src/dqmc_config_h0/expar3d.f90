  ! do not forgett to declare expar as complex in calling programm.
  function expar(i,nf,xmag,flux_x,flux_y,dimer)
    ! For honeycomb lattice, nn hopping phase
    use model_para, only: latt
    implicit none
    complex(dp) :: expar
    integer, intent(in) :: i, nf
    real(dp), intent(in) :: xmag, flux_x, flux_y,dimer
  
    !  local 
    integer :: nax, nay
    real(dp) :: x, x1, xmag1, avec(2), rvec(2)

    expar = dcmplx(1.d0,0.d0)
    
  end function expar
