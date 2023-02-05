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

    nax = latt%nnveci(1,nf,i)
    nay = latt%nnveci(2,nf,i)
  
    ! xmag is magnetic xmag per plaquette.
    ! flux  is the twisting of boundary condition in  x-direction.
    ! both xmag and flux are in units of flux quantum.
    
    !write(6,*) 'lq in thop_mag: ', lq
  
  
    !     uses landau gauge to compute the matix element 
    !     c^{dagger}_i c_j exp(2 pi i / phi_0 \int_i^j a dl),  j = i + (nax,nay)
    !     a(x) = -b(x_2,0,0) with bondary conditions.
  
    xmag1  = xmag/dble(latt%ncell)*(latt%a1_p(1)*latt%a2_p(2) - latt%a1_p(2)*latt%a2_p(1))

    x = -xmag1*latt%nnvecr(1,nf,i)*(2.d0*latt%rcord(i,2) + latt%nnvecr(2,nf,i))*pi
  
    x1 = 0.d0
    if ( ( ( latt%cord(i,2) == latt%l2  .and. nay > 0 ) .or. ( latt%cord(i,2) == 1  .and. nay < 0 ) ) .and. nax == 0)  then
       x1 = 2.0*pi * xmag1 * dble(sign(1,nay)*latt%l2) * latt%a2_p(2) * &
                             ( (latt%rcord(i,2)+latt%nnvecr(2,nf,i))*latt%a2_p(1) + (latt%rcord(i,1)+latt%nnvecr(1,nf,i))*latt%a1_p(1) )
    else if ( ( ( latt%cord(i,1) == latt%l1  .and. nax > 0 ) .or. ( latt%cord(i,1) == 1  .and. nax < 0 ) ) .and. nay == 0)  then
       x1 = 2.0*pi * xmag1 * dble(sign(1,nax)*latt%l1) * latt%a1_p(2) * &
                             ( (latt%rcord(i,2)+latt%nnvecr(2,nf,i))*latt%a2_p(1) + (latt%rcord(i,1)+latt%nnvecr(1,nf,i))*latt%a1_p(1) )
    else if ( ( ( latt%cord(i,1) == latt%l1  .and. nax > 0 ) .or. ( latt%cord(i,1) == 1  .and. nax < 0 ) ) .and. &
              ( ( latt%cord(i,2) == latt%l2  .and. nay > 0 ) .or. ( latt%cord(i,2) == 1  .and. nay < 0 ) ) ) then
       x1 = 2.0*pi * xmag1 * ( dble(sign(1,nax)*latt%l1) * latt%a1_p(2) + dble(sign(1,nay)*latt%l2) * latt%a2_p(2) ) * &
                             ( (latt%rcord(i,2)+latt%nnvecr(2,nf,i))*latt%a2_p(1) + (latt%rcord(i,1)+latt%nnvecr(1,nf,i))*latt%a1_p(1) )
    end if
    
    expar = exp( dcmplx(0.d0, x+x1) )
    
    !     flux.
    ! shift k point with vector flux_x*a1_p/l1 + flux_y*a2_p/l2
    avec = flux_x*latt%a1_p/dble(latt%l1) + flux_y*latt%a2_p/dble(latt%l2)
    rvec = dble(nax)*latt%a1_p + dble(nay)*latt%a2_p
    expar = expar*exp( dcmplx(0.d0, 2.d0*pi*( avec(1)*rvec(1) + avec(2)*rvec(2) )) )

    !     dimerization.
    if( nf == 1 ) then
        expar = expar*dcmplx(1.d0+dimer,0.d0)
    end if
    
  end function expar
