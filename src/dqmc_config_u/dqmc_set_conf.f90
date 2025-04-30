  subroutine dqmc_set_conf( this, lq, ltrot, nu, lproju, alpha, theta, u, dtau )
    implicit none
    class(uconf) :: this
    integer :: lq, ltrot
    logical :: lproju
    real(dp) :: nu, alpha, theta, u, dtau

    ! local
    integer :: i,j,is,m,iflip, isp
    complex(dp) :: z0, z1, z2, z3, z1h
    real(dp) :: a, d  ! Parameters for exact HS transformation
    real(dp) :: xp, xm ! Arguments for acos/acosh in exact HS

    this%lq = latt%ncell
    this%nsites = latt%nsites
    this%ltrot = ltrot

    ! Determine the number of components
    if (lproju ) then
        if( dble(int(nu))-nu == 0.d0 ) then ! if nu is integer
            this%lcomp = nflr/2+1+int( 0.5d0*abs(nu) )
        else ! if nu is not integer
            stop "Error!!! In the projection case, please set nu to be integer !!! "
        end if
    else
#ifdef EXACTHS
        ! Determine the number of components based on nflr
        select case (nflr)
        case (2)
            this%lcomp = 2
        case (4, 6)
            this%lcomp = 4
        case default
            ! use general 4-component HS transformation
            write(fout,'(a)') " warning: the exact HS transformations for SU(N>6) Hubbard models are not implemented, we turn to general 4-component HS transformation. " 
            this%lcomp = 4
        end select
#else
        this%lcomp = 4
#endif
    end if
    
    ! set etal, alpha and gaml which together define the Hubbard-Stratonovich (HS) transformation
    allocate(this%etal(this%lcomp), this%gaml(this%lcomp))
#ifdef EXACTHS
    ! For exact HS, set parameters based on nflr and sign of u
    select case (nflr)
    case (2)  ! SU(2) case
        ! set eta: here -1,1 is mapped to indices 1,2
        this%etal = (/-1.d0,1.d0/)

        ! Set alpha
        if (u .ge. 0.d0) then
            this%alpha = dacos(dexp(-dtau*dabs(u)*0.5d0))  ! Repulsive case
        else
            this%alpha = dacosh(dexp(dtau*dabs(u)*0.5d0))  ! Attractive case
        end if
        
        ! Set gamma - SU(2) exact HS do not use gamma so both are 1
        this%gaml(1:2) = 1.d0

    case (4, 6)  ! SU(4) or SU(6) case
        ! Set eta: here -2,-1,1,2 is mapped to indices 1,2,3,4
        a = dexp(-0.5d0 * dtau * u)
        d = dsqrt(8.d0 + a**2 * (3.d0 + a**2)**2)
        xp = (a+2.d0*a**3+a**5+(a**2-1.d0)*d)/4.d0
        xm = (a+2.d0*a**3+a**5-(a**2-1.d0)*d)/4.d0
        if (u .ge. 0.d0) then
            this%etal = (/-dacos(xm), -dacos(xp), dacos(xp), dacos(xm)/)
        else
            this%etal = (/-dacosh(xm), -dacosh(xp), dacosh(xp), dacosh(xm)/)
        end if
        
        ! Set alpha - SU(4/6) exact HS do not use alpha so all are 1
        this%alpha = 1.d0
        
        ! Set gamma
        this%gaml = (/  (a*(3.d0+a**2)+d)/d, (-a*(3.d0+a**2)+d)/d, &
                       (-a*(3.d0+a**2)+d)/d,  (a*(3.d0+a**2)+d)/d  /)
        
    case default
        ! use general 4-component HS transformation
        ! Set eta: here -2,-1,1,2 is mapped to indices 1,2,3,4
        this%etal(1:4) = (/  - dsqrt(2.d0 * ( 3.d0 + dsqrt(6.d0) ) ), &
                             - dsqrt(2.d0 * ( 3.d0 - dsqrt(6.d0) ) ), &
                               dsqrt(2.d0 * ( 3.d0 - dsqrt(6.d0) ) ), &
                               dsqrt(2.d0 * ( 3.d0 + dsqrt(6.d0) ) ) /)

        ! set alpha
        this%alpha = dsqrt( dtau*dabs(u)*0.5d0 )

        ! Set gamma
        this%gaml(1:4) =  (/ 1.d0 - dsqrt(6.d0)/3.d0, &
                             1.d0 + dsqrt(6.d0)/3.d0, &
                             1.d0 + dsqrt(6.d0)/3.d0, &
                             1.d0 - dsqrt(6.d0)/3.d0 /)
    end select
#else
    ! general 4-component HS transformation
    ! Set eta: here -2,-1,1,2 is mapped to indices 1,2,3,4
    this%etal(1:4) = (/  - dsqrt(2.d0 * ( 3.d0 + dsqrt(6.d0) ) ), &
                         - dsqrt(2.d0 * ( 3.d0 - dsqrt(6.d0) ) ), &
                           dsqrt(2.d0 * ( 3.d0 - dsqrt(6.d0) ) ), &
                           dsqrt(2.d0 * ( 3.d0 + dsqrt(6.d0) ) ) /)
    
    ! set alpha
    this%alpha = dsqrt( dtau*dabs(u)*0.5d0 )
    
    ! Set gamma
    this%gaml(1:4) =  (/ 1.d0 - dsqrt(6.d0)/3.d0, &
                         1.d0 + dsqrt(6.d0)/3.d0, &
                         1.d0 + dsqrt(6.d0)/3.d0, &
                         1.d0 - dsqrt(6.d0)/3.d0 /)
#endif

    call allocate_zvfunc(this%bmat, this%lcomp)
    call allocate_zvfunc(this%bmat_inv, this%lcomp)
    call allocate_gfunc( this%delta_bmat, this%lcomp-1, this%lcomp)
    allocate( this%phase(this%lcomp) )
    allocate( this%phase_ratio(this%lcomp-1, this%lcomp) )

    do is = 1, this%lcomp
      if (lproju ) then
          this%phase(is) = exp( dcmplx( 0.d0, pi*dble(-2*is)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
          this%bmat%blk1(is) = exp( dcmplx(0.d0, pi*dble(2*is)/dble(this%lcomp)) )
          this%bmat_inv%blk1(is) = exp( dcmplx(0.d0,-pi*dble(2*is)/dble(this%lcomp)) )
      else
          if( u .ge. 0.d0 ) then
              this%phase(is) = exp( dcmplx(0.d0, this%alpha*(-this%etal(is))*0.5d0*(dble(nflr)+nu)) )
              this%bmat%blk1(is) =      exp( dcmplx(0.d0, this%alpha*this%etal(is)) )
              this%bmat_inv%blk1(is) =  exp( dcmplx(0.d0,-this%alpha*this%etal(is)) )
          else
              this%phase(is) = exp( dcmplx(this%alpha*(-this%etal(is))*0.5d0*(dble(nflr)+nu), 0.d0) ) ! this is not a "phase" actually, use it (or this%phase_ratio) to calculate ratio but do not use it to set phase
              this%bmat%blk1(is) =      exp( dcmplx( this%alpha*this%etal(is), 0.d0) )
              this%bmat_inv%blk1(is) =  exp( dcmplx(-this%alpha*this%etal(is), 0.d0) )
          end if
      end if
#IFDEF PLEVEL2
      if( irank == 0 ) then
          write(fout,'(a)') " in uconf%set_conf "
          write(fout,'(a,i3,a)') "bmat%blk1(",is," ) = "
          write(fout,'(2f9.5)') this%bmat%blk1(is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_inv%blk1(",is," ) = "
          write(fout,'(2f9.5)') this%bmat_inv%blk1(is)
          write(fout,'(a)') " "
      end if
#ENDIF
    end do
    do is = 1, this%lcomp
    do iflip = 1, this%lcomp-1
      isp = mod(is+iflip-1,this%lcomp) + 1
      if(lproju) then
          this%phase_ratio(iflip,is) = exp( dcmplx( 0.d0, pi*dble(2*is-2*isp)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
          this%delta_bmat%blk1(iflip,is) = exp( dcmplx(0.d0, pi*dble(2*isp-2*is)/dble(this%lcomp)) ) - cone
      else
          if( u .ge. 0.d0 ) then
              this%phase_ratio(iflip,is) = exp( dcmplx(0.d0, this%alpha*(this%etal(is)-this%etal(isp))*0.5d0*(dble(nflr)+nu)) )
              this%delta_bmat%blk1(iflip,is) = exp( dcmplx(0.d0, this%alpha*(this%etal(isp)-this%etal(is))) ) - cone
          else
              this%phase_ratio(iflip,is) = exp( dcmplx(this%alpha*(this%etal(is)-this%etal(isp))*0.5d0*(dble(nflr)+nu),0.d0) )
              this%delta_bmat%blk1(iflip,is) = exp( dcmplx(this%alpha*(this%etal(isp)-this%etal(is)),0.d0) ) - cone
          end if
      end if
#IFDEF PLEVEL2
      if( irank == 0 ) then
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a)') "delta_bmat%blk1(",iflip," ,",is," ) = "
          write(fout,'(2f9.5)') this%delta_bmat%blk1(iflip,is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a,2f16.8)') "phase_ratio(",iflip," ,",is," ) = ", this%phase_ratio(iflip,is)
      end if
#ENDIF
    end do
    end do
  end subroutine dqmc_set_conf
