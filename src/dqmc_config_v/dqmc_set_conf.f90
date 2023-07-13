  subroutine dqmc_set_conf( this, lq, ltrot, nu, lprojv, alpha, theta, v, dtau )
    implicit none
    class(vconf) :: this
    integer :: lq, ltrot
    logical :: lprojv
    real(dp) :: nu, alpha, theta, v, dtau

    ! local
    integer :: i,j,is,m,iflip, isp
    complex(dp) :: z0, z1, z2, z3, z1h
    this%lq = latt%ncell
    this%nsites = latt%nsites
    this%nn_nf = latt%nn_nf
    this%ltrot = ltrot
    this%lcomp = 2
    this%lambda = 0.5d0*exp(-dtau*v*0.25d0)
    this%alpha = dcmplx( acosh(exp(dtau*v*0.5d0)), 0.d0)
    call allocate_zvfunc( this%bmat_p, this%lcomp )
    call allocate_zvfunc( this%bmat_m, this%lcomp )
    call allocate_zvfunc( this%bmat_p_inv, this%lcomp )
    call allocate_zvfunc( this%bmat_m_inv, this%lcomp )
    call allocate_gfunc( this%delta_bmat_p, this%lcomp-1, this%lcomp )
    call allocate_gfunc( this%delta_bmat_m, this%lcomp-1, this%lcomp )
    allocate( this%phase_ratio(this%lcomp-1,this%lcomp) )
    allocate( this%phase(this%lcomp) )

    do is = 1, this%lcomp
      if (lprojv ) then
          this%bmat_p%orb1(is)     = exp( dcmplx(0.d0, pi*dble(2*is)/dble(this%lcomp)) )
          this%bmat_p_inv%orb1(is) = exp( dcmplx(0.d0,-pi*dble(2*is)/dble(this%lcomp)) )
          this%bmat_m%orb1(is)     = exp( dcmplx(0.d0, pi*dble(2*is)/dble(this%lcomp)) )
          this%bmat_m_inv%orb1(is) = exp( dcmplx(0.d0,-pi*dble(2*is)/dble(this%lcomp)) )
      else
          this%bmat_p%orb1(is)     =  exp(  this%alpha*etal(is) )
          this%bmat_p_inv%orb1(is) =  exp( -this%alpha*etal(is) )
          this%bmat_m%orb1(is)     =  exp( -this%alpha*etal(is) ) 
          this%bmat_m_inv%orb1(is) =  exp(  this%alpha*etal(is) )
      end if
      if( irank == 0 ) then
          write(fout,'(a)') " in vconf%set_conf "
          write(fout,'(a,i3,a)') "bmat_p%orb1(",is," ) = "
          write(fout,'(2f9.5)') this%bmat_p%orb1(is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_p_inv%orb1(",is," ) = "
          write(fout,'(2f9.5)') this%bmat_p_inv%orb1(is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_m%orb1(",is," ) = "
          write(fout,'(2f9.5)') this%bmat_m%orb1(is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_m_inv1%orb1(",is," ) = "
          write(fout,'(2f9.5)') this%bmat_m_inv%orb1(is)
          write(fout,'(a)') " "

      end if
    end do
    do is = 1, this%lcomp
    do iflip = 1, this%lcomp-1
      isp = mod(is+iflip-1,this%lcomp) + 1
      if(lprojv) then
          this%phase_ratio(iflip,is) = exp( dcmplx( 0.d0, pi*dble(2*is-2*isp)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
          this%delta_bmat_p%orb1(iflip,is) = exp( dcmplx(0.d0, pi*dble(2*isp-2*is)/dble(this%lcomp)) ) - cone
          this%delta_bmat_m%orb1(iflip,is) = exp( this%alpha*(etal(isp)-etal(is)) ) - cone
      else
          this%phase_ratio(iflip,is) = exp( dcmplx(0.d0, this%alpha*(etal(is)-etal(isp))*(dble(nflr)+nu)) )
          this%delta_bmat_p%orb1(iflip,is) = exp(  this%alpha*(etal(isp)-etal(is)) ) - cone
          this%delta_bmat_m%orb1(iflip,is) = exp( -this%alpha*(etal(isp)-etal(is)) ) - cone
      end if
      if( irank == 0 ) then
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a)') "delta_bmat_p%orb1(",iflip," ,",is," ) = "
          write(fout,'(2f9.5)') this%delta_bmat_p%orb1(iflip,is)
          write(fout,'(a,i3,a,i3,a)') "delta_bmat_m%orb1(",iflip," ,",is," ) = "
          write(fout,'(2f9.5)') this%delta_bmat_m%orb1(iflip,is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a,2f16.8)') "phase_ratio(",iflip," ,",is," ) = ", this%phase_ratio(iflip,is)
      end if
    end do
    end do
  end subroutine dqmc_set_conf
