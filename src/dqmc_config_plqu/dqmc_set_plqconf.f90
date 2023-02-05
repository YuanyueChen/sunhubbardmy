  subroutine dqmc_set_plqconf( this, lq, ltrot, nu, lprojplq, alpha, theta, plqu, dtau )
    use model_para, only: latt
    implicit none
    class(plqconf) :: this
    integer :: lq, ltrot
    logical :: lprojplq
    real(dp) :: nu, alpha, theta, plqu, dtau

    ! local
    complex(dp), allocatable, dimension(:,:) :: vstmat_tmp, umat_tmp
    real(dp), allocatable, dimension(:) :: eig_tmp
    integer :: i,j,is,m,iflip, isp
    complex(dp) :: z0, z1, z2, z3, z1h
    this%lq = lq
    this%ltrot = ltrot
    if (lprojplq ) then
#IFDEF HONEYCOMB
        if( dble(int(nu))-nu == 0.d0 ) then ! if nu is integer
            this%lcomp = latt%nn_nf*int( nflr+abs(nu) ) + 1
#ELIF SQUARE
        if( dble(int(dble(latt%nn_nf)*nu)) - dble(latt%nn_nf)*nu == 0.d0 ) then !
            this%lcomp = int( dble(latt%nn_nf)*( 0.5d0*dble(nflr)+abs(nu) ) ) + 1
#ELSE
        if( dble(int(dble(latt%nn_nf)*nu)) - dble(latt%nn_nf)*nu == 0.d0 ) then !
            this%lcomp = int( dble(latt%nn_nf)*( 0.5d0*dble(nflr)+abs(nu) ) ) + 1
#ENDIF
        else ! if nu is not integer
            stop "Error!!! In the projection case, please set nu to be integer !!! "
        end if
    else
        this%lcomp = 4
    end if
    this%alpha_plqu = dsqrt( dtau*plqu )
    allocate( this%bmat_plqu_orb1(latt%z_plq,latt%z_plq,this%lcomp) )
    allocate( this%bmat_plqu_orb1_inv(latt%z_plq,latt%z_plq,this%lcomp) )
    allocate( this%delta_bmat_plqu_orb1(latt%z_plq,latt%z_plq,this%lcomp-1,this%lcomp) )
    allocate( this%phase(this%lcomp) )
    allocate( this%phase_ratio(this%lcomp-1, this%lcomp) )

    allocate( vstmat_tmp(latt%z_plq,latt%z_plq) )
    allocate( umat_tmp(latt%z_plq,latt%z_plq) )
    allocate( eig_tmp(latt%z_plq) )

    z0 = dcmplx(1.d0/dble(latt%nn_nf), 0.d0)
    z1  = dcmplx( alpha, 0.d0)*exp( dcmplx(0.d0, -theta) )
    z1h = dcmplx( alpha, 0.d0)*exp( dcmplx(0.d0,  theta) )
    if( latt%z_plq == 6 ) then ! honeycomb lattice
        vstmat_tmp(1:latt%z_plq,1:latt%z_plq) = reshape( &
        (/ z0,    z1,    czero, czero, czero,-z1  &
          ,z1h,   z0,   -z1h,   czero, czero, czero &
          ,czero,-z1,    z0,    z1,    czero, czero &
          ,czero, czero, z1h,   z0,   -z1h,   czero &
          ,czero, czero, czero,-z1,    z0,    z1 &
          ,-z1h,  czero, czero, czero, z1h,   z0   /) &
          ,(/latt%z_plq,latt%z_plq/))
    else if( latt%z_plq == 4 ) then ! square lattice
        vstmat_tmp(1:latt%z_plq,1:latt%z_plq) = reshape( &
        (/ z0,    z1,    czero,-z1  &
          ,z1h,   z0,   -z1h,   czero &
          ,czero,-z1,    z0,    z1 &
          ,-z1h,  czero, z1h,   z0   /) &
          ,(/latt%z_plq,latt%z_plq/))
    end if
    call s_eig_he(latt%z_plq,latt%z_plq,vstmat_tmp,eig_tmp,umat_tmp)
    if( irank == 0 ) then
        write(fout,'(a)') " "
        write(fout,'(a)') "vstmat_tmp(:,:) = "
        do i = 1, latt%z_plq
            write(fout,'(6(2f9.5,2x))') vstmat_tmp(i,:)
        end do
        write(fout,'(a)') " "
        write(fout,'(a)') "umat_tmp(:,:) = "
        do i = 1, latt%z_plq
            write(fout,'(6(2f9.5,2x))') umat_tmp(i,:)
        end do
        write(fout,'(a)') " "
        write(fout,'(a)') "eig_tmp(:) = "
        do i = 1, latt%z_plq
            write(fout,'(f9.5)') eig_tmp(i)
        end do
    end if

    do is = 1, this%lcomp
      do i = 1, latt%z_plq
         do j = 1, latt%z_plq
            z0 = czero
            z1 = czero
            z2 = czero
            z3 = czero
            do m = 1, latt%z_plq
               if (lprojplq ) then
                 z0 = z0 +         umat_tmp(i,m)  * exp( dcmplx(0.d0, eig_tmp(m)*pi*dble(2*latt%nn_nf*is)/dble(this%lcomp)) ) * dconjg(umat_tmp(j,m))
                 z1 = z1 +  dconjg(umat_tmp(i,m)) * exp( dcmplx(0.d0, eig_tmp(m)*pi*dble(2*latt%nn_nf*is)/dble(this%lcomp)) ) *        umat_tmp(j,m)
                 z2 = z2 +         umat_tmp(i,m)  * exp( dcmplx(0.d0,-eig_tmp(m)*pi*dble(2*latt%nn_nf*is)/dble(this%lcomp)) ) * dconjg(umat_tmp(j,m))
                 z3 = z3 +  dconjg(umat_tmp(i,m)) * exp( dcmplx(0.d0,-eig_tmp(m)*pi*dble(2*latt%nn_nf*is)/dble(this%lcomp)) ) *        umat_tmp(j,m)
               else
                 z0 = z0 +         umat_tmp(i,m)  * exp( dcmplx(0.d0, eig_tmp(m)*this%alpha_plqu*etal(is)) ) * dconjg(umat_tmp(j,m))
                 z1 = z1 +  dconjg(umat_tmp(i,m)) * exp( dcmplx(0.d0, eig_tmp(m)*this%alpha_plqu*etal(is)) ) *        umat_tmp(j,m)
                 z2 = z2 +         umat_tmp(i,m)  * exp( dcmplx(0.d0,-eig_tmp(m)*this%alpha_plqu*etal(is)) ) * dconjg(umat_tmp(j,m))
                 z3 = z3 +  dconjg(umat_tmp(i,m)) * exp( dcmplx(0.d0,-eig_tmp(m)*this%alpha_plqu*etal(is)) ) *        umat_tmp(j,m)
               end if
            enddo
            this%bmat_plqu_orb1(i,j,is) = z0
            this%bmat_plqu_orb1_inv(i,j,is) = z2
         enddo
      enddo
      if(lprojplq) then
#IFDEF HONEYCOMB
          this%phase(is) = exp( dcmplx( 0.d0, pi*dble(-2*latt%nn_nf*is)*(dble(nflr)+nu)/dble(this%lcomp) ) )
#ELIF SQUARE
          this%phase(is) = exp( dcmplx( 0.d0, pi*dble(-2*latt%nn_nf*is)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
#ELSE
          this%phase(is) = exp( dcmplx( 0.d0, pi*dble(-2*latt%nn_nf*is)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
#ENDIF
      else
#IFDEF HONEYCOMB
          this%phase(is) = exp( dcmplx(0.d0, this%alpha_plqu*(-etal(is))*(dble(nflr)+nu)) )
#ELIF SQUARE
          this%phase(is) = exp( dcmplx(0.d0, this%alpha_plqu*(-etal(is))*0.5d0*(dble(nflr)+nu)) )
#ELSE
          this%phase(is) = exp( dcmplx(0.d0, this%alpha_plqu*(-etal(is))*0.5d0*(dble(nflr)+nu)) )
#ENDIF
      end if
      if( irank == 0 ) then
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_plqu_orb1(:,:,",is," ) = "
          do i = 1, latt%z_plq
                write(fout,'(6(2f9.5,2x))') this%bmat_plqu_orb1(i,:,is)
          end do
          write(fout,'(a)') " "
          write(fout,'(a,i3,a)') "bmat_plqu_orb1_inv(:,:,",is," ) = "
          do i = 1, latt%z_plq
                write(fout,'(6(2f9.5,2x))') this%bmat_plqu_orb1_inv(i,:,is)
          end do
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,2f16.8)') "phase(",is," ) = ", this%phase(is)
      end if
    end do
    do is = 1, this%lcomp
    do iflip = 1, this%lcomp-1
      isp = mod(is+iflip-1,this%lcomp) + 1
      do i = 1, latt%z_plq
         do j = 1, latt%z_plq
            z0 = czero
            z1 = czero
            do m = 1, latt%z_plq
               if(lprojplq) then
                 z0 = z0 +         umat_tmp(i,m)  * exp( dcmplx(0.d0, eig_tmp(m)*pi*dble(2*latt%nn_nf*isp-2*latt%nn_nf*is)/dble(this%lcomp)) ) * dconjg(umat_tmp(j,m))
                 z1 = z1 +  dconjg(umat_tmp(i,m)) * exp( dcmplx(0.d0, eig_tmp(m)*pi*dble(2*latt%nn_nf*isp-2*latt%nn_nf*is)/dble(this%lcomp)) ) *        umat_tmp(j,m)
               else
                 z0 = z0 +         umat_tmp(i,m)  * exp( dcmplx(0.d0, eig_tmp(m)*this%alpha_plqu*(etal(isp)-etal(is))) ) * dconjg(umat_tmp(j,m))
                 z1 = z1 +  dconjg(umat_tmp(i,m)) * exp( dcmplx(0.d0, eig_tmp(m)*this%alpha_plqu*(etal(isp)-etal(is))) ) *        umat_tmp(j,m)
               end if
            enddo
            if( i .ne. j ) then
                this%delta_bmat_plqu_orb1(i,j,iflip,is) = z0
            else
                this%delta_bmat_plqu_orb1(i,j,iflip,is) = z0 - cone
            end if
         enddo
      enddo
      if(lprojplq) then
#IFDEF HONEYCOMB
          this%phase_ratio(iflip,is) = exp( dcmplx( 0.d0, pi*dble(2*latt%nn_nf*is-2*latt%nn_nf*isp)*(dble(nflr)+nu)/dble(this%lcomp) ) )
#ELIF SQUARE
          this%phase_ratio(iflip,is) = exp( dcmplx( 0.d0, pi*dble(2*latt%nn_nf*is-2*latt%nn_nf*isp)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
#ELSE
          this%phase_ratio(iflip,is) = exp( dcmplx( 0.d0, pi*dble(2*latt%nn_nf*is-2*latt%nn_nf*isp)*0.5d0*(dble(nflr)+nu)/dble(this%lcomp) ) )
#ENDIF
      else
#IFDEF HONEYCOMB
          this%phase_ratio(iflip,is) = exp( dcmplx(0.d0, this%alpha_plqu*(etal(is)-etal(isp))*(dble(nflr)+nu)) )
#ELIF SQUARE
          this%phase_ratio(iflip,is) = exp( dcmplx(0.d0, this%alpha_plqu*(etal(is)-etal(isp))*0.5d0*(dble(nflr)+nu)) )
#ELSE
          this%phase_ratio(iflip,is) = exp( dcmplx(0.d0, this%alpha_plqu*(etal(is)-etal(isp))*0.5d0*(dble(nflr)+nu)) )
#ENDIF
      end if
      if( irank == 0 ) then
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a)') "delta_bmat_plqu_orb1(:,:,",iflip," ,",is," ) = "
          do i = 1, latt%z_plq
                write(fout,'(6(2f9.5,2x))') this%delta_bmat_plqu_orb1(i,:,iflip,is)
          end do
          write(fout,'(a)') " "
          write(fout,'(a,i3,a,i3,a,2f16.8)') "phase_ratio(",iflip," ,",is," ) = ", this%phase_ratio(iflip,is)
      end if
    end do
    end do

    deallocate( vstmat_tmp )
    deallocate( umat_tmp )
    deallocate( eig_tmp )
  end subroutine dqmc_set_plqconf
