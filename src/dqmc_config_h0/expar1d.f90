  ! do not forgett to declare expar as complex in calling programm.
  function expar(i,nf,xmag,flux_x,flux_y,dimer)
    ! For chain lattice, no hopping phase. flux_y is meaningless. 
    use model_para, only: latt
    implicit none
    complex(dp) :: expar
    integer, intent(in) :: i, nf
    real(dp), intent(in) :: xmag, flux_x, flux_y, dimer
  
    !  local 
    integer :: nax, nay
    real(dp) :: x, x1, xmag1, avec(1), rvec(1)

    nax = latt%nnveci(1,nf,i)
  
    !! xmag is magnetic xmag per plaquette.
    ! flux  is the twisting of boundary condition in  x-direction.
    ! both xmag and flux are in units of flux quantum.

    !     xmg.
    x = 0.d0 ! bulk hopping phase
    x1 = 0.d0 ! boundary hopping phase
    
    expar = exp( dcmplx(0.d0, x+x1) )
    
    !     flux.
    ! shift k point with vector flux_x*a1_p/l1
    avec = flux_x*latt%a1_p/dble(latt%l1)
    rvec = dble(nax)*latt%a1_p
    expar = expar*exp( dcmplx(0.d0, 2.d0*pi*( avec(1)*rvec(1))) )

    !     dimerization.
    if( nf == 1 ) then
        expar = expar*dcmplx(1.d0+dimer,0.d0)
    end if
    
  end function expar
