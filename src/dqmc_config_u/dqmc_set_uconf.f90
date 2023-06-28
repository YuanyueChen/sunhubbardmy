  subroutine dqmc_set_uconf( this, lq, ltrot, nu, lproju, alpha, theta, u, dtau )
    implicit none
    class(uconf) :: this
    integer :: lq, ltrot
    logical :: lproju
    real(dp) :: nu, alpha, theta, u, dtau

    ! local
    integer :: i,j,is,m,iflip, isp
    complex(dp) :: z0, z1, z2, z3, z1h
    this%lq = latt%ncell
    this%nsites = latt%nsites
    this%ltrot = ltrot
    if (lproju ) then
        if( dble(int(nu))-nu == 0.d0 ) then ! if nu is integer
            this%lcomp = nflr/2+1+int( 0.5d0*abs(nu) )
        else ! if nu is not integer
            stop "Error!!! In the projection case, please set nu to be integer !!! "
        end if
    else
        this%lcomp = 4
    end if
    this%alpha_u = dsqrt( dtau*u*0.5d0 )
    allocate( this%bmat_u_orb1(this%lcomp) )
    allocate( this%bmat_u_orb1_inv(this%lcomp) )
    allocate( this%delta_bmat_u_orb1(this%lcomp-1,this%lcomp) )
    allocate( this%phase(this%lcomp) )
    allocate( this%phase_ratio(this%lcomp-1, this%lcomp) )

    do is = 1, this%lcomp
      if (lproju ) then
          this%phase(is) = exp( dcmplx( 0.d0, pi*dble(-2*is)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
          this%bmat_u_orb1(is) = exp( dcmplx(0.d0, pi*dble(2*is)/dble(this%lcomp)) )
          this%bmat_u_orb1_inv(is) = exp( dcmplx(0.d0,-pi*dble(2*is)/dble(this%lcomp)) )
      else
          this%phase(is) = exp( dcmplx(0.d0, this%alpha_u*(-etal(is))*0.5d0*(dble(nflr)+nu)) )
          this%bmat_u_orb1(is) =      exp( dcmplx(0.d0, this%alpha_u*etal(is)) )
          this%bmat_u_orb1_inv(is) =  exp( dcmplx(0.d0,-this%alpha_u*etal(is)) )
      end if
#IFDEF PLEVEL2
      if( irank == 0 ) then
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_u_orb1(",is," ) = "
          write(fout,'(2f9.5)') this%bmat_u_orb1(is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_u_orb1_inv(",is," ) = "
          write(fout,'(2f9.5)') this%bmat_u_orb1_inv(is)
          write(fout,'(a)') " "
      end if
#ENDIF
    end do
    do is = 1, this%lcomp
    do iflip = 1, this%lcomp-1
      isp = mod(is+iflip-1,this%lcomp) + 1
      if(lproju) then
          this%phase_ratio(iflip,is) = exp( dcmplx( 0.d0, pi*dble(2*is-2*isp)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
          this%delta_bmat_u_orb1(iflip,is) = exp( dcmplx(0.d0, pi*dble(2*isp-2*is)/dble(this%lcomp)) ) - cone
      else
          this%phase_ratio(iflip,is) = exp( dcmplx(0.d0, this%alpha_u*(etal(is)-etal(isp))*0.5d0*(dble(nflr)+nu)) )
          this%delta_bmat_u_orb1(iflip,is) = exp( dcmplx(0.d0, this%alpha_u*(etal(isp)-etal(is))) ) - cone
      end if
#IFDEF PLEVEL2
      if( irank == 0 ) then
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a)') "delta_bmat_u_orb1(",iflip," ,",is," ) = "
          write(fout,'(2f9.5)') this%delta_bmat_u_orb1(iflip,is)
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a,2f16.8)') "phase_ratio(",iflip," ,",is," ) = ", this%phase_ratio(iflip,is)
      end if
#ENDIF
    end do
    end do
  end subroutine dqmc_set_uconf
