  subroutine dqmc_set_h0conf(this, lq, ltrot, rt, mu )
    use model_para, only: latt, xmag, flux_x, flux_y, dimer, dtau
    implicit none
    class(h0conf) :: this
    integer :: lq, ltrot
    real(dp) :: rt, mu

    ! local
    integer :: i,j,m
    integer, parameter :: nch = 2
    complex(dp) ::  z, z0, z1

    complex(dp), allocatable ::  h0mat(:,:), umat_tmp(:,:)
    real    (dp), allocatable ::   eig_tmp(:)
    integer info
    integer :: nf, iu, ist, i0, i1

    this%lq = lq
    this%ltrot = ltrot
    this%rt = rt
    this%mu = mu
#IFDEF BREAKUP_T
    allocate( this%urt(latt%nn_nf*latt%nn_lf,2,2) )
    allocate( this%urtm1(latt%nn_nf*latt%nn_lf,2,2) )

    allocate ( h0mat(nch,nch), umat_tmp(nch,nch), eig_tmp(nch)  )
    do nf = 1, latt%nn_nf
        do iu = 1,latt%nn_lf
           ist = iu + (nf - 1)*latt%nn_lf
           i0 = latt%nnlf_list(iu) ! A site
           i1 = latt%nnlist(i0,nf) ! B site
           
           h0mat = czero
           
           z = dcmplx(-rt,0.d0)*expar(i0,nf,xmag,flux_x,flux_y,dimer)
           h0mat(1,2) = z
           h0mat(2,1) = dconjg(z)
           h0mat(1,1) = dcmplx(-mu/dble(latt%nn_nf),0.d0) ! each site has been counted latt%nn_nf times
           h0mat(2,2) = dcmplx(-mu/dble(latt%nn_nf),0.d0)
           
           call s_eig_he(nch,nch,h0mat,eig_tmp,umat_tmp)
           
           do i = 1,nch
              do j = 1,nch
                 z0 = czero
                 z1 = czero
                 do m = 1,nch
                    z0 = z0 +  umat_tmp(i,m) * exp(-0.5d0*dtau *eig_tmp(m)) * dconjg(umat_tmp(j,m))
                    z1 = z1 +  umat_tmp(i,m) * exp( 0.5d0*dtau *eig_tmp(m)) * dconjg(umat_tmp(j,m))
                 enddo
                 this%urt  (ist,i,j) = z0
                 this%urtm1(ist,i,j) = z1
              enddo
           enddo
        enddo
    enddo

    deallocate ( h0mat, umat_tmp, eig_tmp )
#ELSE
    allocate( this%urt(ndim,ndim) )
    allocate( this%urtm1(ndim,ndim) )

    allocate ( h0mat(ndim,ndim), umat_tmp(ndim,ndim), eig_tmp(ndim)  )
    
    call seth0(h0mat,rt, mu)

    call s_eig_he(ndim,ndim,h0mat,eig_tmp,umat_tmp)
          
    do i = 1,ndim
       do j = 1,ndim
          z0 = czero
          z1 = czero
          do m = 1,ndim
             z0 = z0 +  umat_tmp(i,m) * exp(-0.5d0*dtau *eig_tmp(m)) * dconjg(umat_tmp(j,m))
             z1 = z1 +  umat_tmp(i,m) * exp( 0.5d0*dtau *eig_tmp(m)) * dconjg(umat_tmp(j,m))
          enddo
          this%urt  (i,j) = z0        !!! the exponential matrix of hopping term
          this%urtm1(i,j) = z1        !!! the exponential matrix of hopping term ( conjugate) 
       enddo
    enddo
#IFDEF TEST
    write(fout,*)
    write(fout,*) ' h0mat(:,:) = '
    do i = 1, ndim
        write(fout,'(18(2e14.6,2x))') h0mat(i,:)
    end do
    write(fout,*)
    write(fout,*) ' this%urt(:,:) = '
    do i = 1, ndim
        write(fout,'(18(2e14.6,2x))') this%urt(i,:)
    end do
    write(fout,*)
    write(fout,*) ' eig_tmp(:) = '
    do i = 1, ndim
        write(fout,'(f16.8)') eig_tmp(i)
    end do
#ENDIF
    deallocate ( h0mat, umat_tmp, eig_tmp  )
#ENDIF
  
  end subroutine dqmc_set_h0conf
