subroutine dqmc_set_h0conf(this, lq, ltrot, rt, mu )
    use model_para, only: latt, xmag, flux_x, flux_y, dimer, dtau
    implicit none
    class(h0conf) :: this
    integer :: lq, ltrot
    real(dp) :: rt, mu

    ! local
    integer :: i,j,m
    integer, parameter :: nch = 2
    real(dp) :: en_free
    complex(dp) ::  z, z0, z1

    integer info
    integer :: nf, iu, ist, i0, i1
    complex(dp), allocatable ::  hmat_tmp(:,:), umat_tmp(:,:)
    real(dp), allocatable :: eig_tmp(:)

! case 1, use trotter
#IFDEF BREAKUP_T
    this%lq = lq
    this%ltrot = ltrot
    this%rt = rt
    this%mu = mu
    allocate( this%urt(latt%nn_nf*latt%nn_lf,2,2) )
    allocate( this%urtm1(latt%nn_nf*latt%nn_lf,2,2) )

    allocate( hmat_tmp(nch,nch), umat_tmp(nch,nch), eig_tmp(nch)  )
    allocate( this%h0mat(latt%nsites,latt%nsites) ) 
    this%h0mat = czero
    do nf = 1, latt%nn_nf
        do iu = 1,latt%nn_lf
           ist = iu + (nf - 1)*latt%nn_lf
           i0 = latt%nnlf_list(iu) ! A site
           i1 = latt%nnlist(i0,nf) ! B site
           
           hmat_tmp = czero
           
           z = dcmplx(-rt,0.d0)*expar(i0,nf,xmag,flux_x,flux_y,dimer)
           hmat_tmp(1,2) = z
           hmat_tmp(2,1) = dconjg(z)
           hmat_tmp(1,1) = dcmplx(-mu/dble(latt%nn_nf),0.d0) ! each site has been counted latt%nn_nf times
           hmat_tmp(2,2) = dcmplx(-mu/dble(latt%nn_nf),0.d0)

           this%h0mat(i0,i1) = this%h0mat(i0,i1) + z
           this%h0mat(i1,i0) = this%h0mat(i1,i0) + dconjg(z)
           this%h0mat(i0,i0) = this%h0mat(i0,i0) + dcmplx(-mu/dble(latt%nn_nf),0.d0) ! each site has been counted latt%nn_nf times
           this%h0mat(i1,i1) = this%h0mat(i1,i1) + dcmplx(-mu/dble(latt%nn_nf),0.d0) ! each site has been counted latt%nn_nf times
           
           call s_eig_he(nch,nch,hmat_tmp,eig_tmp,umat_tmp)
           
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

    deallocate ( hmat_tmp, umat_tmp, eig_tmp )
! case 2, use fft
#ELIF FFT
    complex(dp), allocatable :: htmp(:,:), htmpC(:,:), h0vec(:)
    type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
    type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim2
#IFDEF CUBIC
    complex(dp), allocatable :: h0k(:)
    integer :: Status, Stride(4), dimm(3)
#ELIF HONEYCOMB
    complex(dp), allocatable :: h0k(:,:,:)
    integer :: Status, Stride(3), dimm(2)
    real(dp) :: eig2(2)
#ELIF SQUARE
    complex(dp), allocatable :: h0k(:)
    integer :: Status, Stride(3), dimm(2)
#ELIF CHAIN
    complex(dp), allocatable :: h0k(:)
    integer :: Status, Stride(2), dimm(1)
#ENDIF
    this%lq = lq
    this%ltrot = ltrot
    this%rt = rt
    this%mu = mu
    allocate( this%h0mat(latt%nsites,latt%nsites) )
    allocate( eig_tmp(latt%nsites) )
    allocate( h0vec(latt%nsites*latt%nsites) )
    call seth0(this%h0mat, rt, mu, xmag, flux_x, flux_y, dimer)

  ! case 2.0, fft for chain
  ! H = U D U^+ => D = U^+ H U
#IFDEF CHAIN
  dimm = (/latt%l1/)
  allocate( this%exph0k(lq), this%exph0kinv(lq) )
  allocate( h0k(lq) )
#IFDEF TEST
  write(fout, '(a)') 'before fft, h0mat = '
  do i = 1, lq
      write(fout, '(16f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  h0vec = reshape( this%h0mat, (/lq*lq/) )
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 1, dimm )
  Status = DftiCreateDescriptor(Desc_Handle_Dim2, DFTI_DOUBLE, DFTI_COMPLEX, 1, dimm )

  ! perform lq two-dimensional transforms along 1st dimension of h0mat
  ! U^+ H
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, h0vec )
  ! perform lq two-dimensional transforms along 2nd dimension of h0mat
  ! (U^+ H) U
  ! Stride is used to locate the data
  ! e.g. suppose our data has structure data(i,j1), where j1 are the index for 2nd dimension
  !      then the address is i + s1*i + s2*j1
  Stride(1) = 0; Stride(2) = lq;
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_NUMBER_OF_TRANSFORMS, lq )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_DISTANCE, 1 )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_DISTANCE, 1 )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_STRIDES, Stride )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_STRIDES, Stride )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim2 )
  Status = DftiComputeBackward( Desc_Handle_Dim2, h0vec )
  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )
  Status = DftiFreeDescriptor( Desc_Handle_Dim2 )

  this%h0mat = reshape( h0vec, (/lq,lq/) )
  do i = 1, lq
      h0k(i) = this%h0mat(i,i)
      this%exph0k(i)    = exp( -dcmplx(0.5d0*dtau,0.d0)*this%h0mat(i,i) )
      this%exph0kinv(i) = exp(  dcmplx(0.5d0*dtau,0.d0)*this%h0mat(i,i) )
  end do
#IFDEF TEST
  write(fout, '(a)') 'after fft, h0mat = '
  do i = 1, lq
      write(fout, '(16f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  if( irank.eq.0 ) then
      write(fout,*)
      write(fout,'(a)') 'fft get eigenvalues = '
      do i = 1, lq
          write(fout,'(e16.8)') dble(h0k(i))
      end do

      en_free = 0.d0
      do i = 1, lq
          if (dble(h0k(i))<=0.d0) then
              en_free = en_free + dble(h0k(i))
          end if
      end do
      write(fout,*)
      write(fout,'(a,e16.8)') ' fft get en_free = ', en_free
  end if

  ! case 2.1, fft for square
  ! H = U D U^+ => D = U^+ H U
#ELIF SQUARE
  dimm = (/latt%l1, latt%l2/)
  allocate( this%exph0k(lq), this%exph0kinv(lq) )
  allocate( h0k(lq) )
#IFDEF TEST
  write(fout, '(a)') 'before fft, h0mat = '
  do i = 1, lq
      write(fout, '(16f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  h0vec = reshape( this%h0mat, (/lq*lq/) )
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm )
  Status = DftiCreateDescriptor(Desc_Handle_Dim2, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm )

  ! perform lq two-dimensional transforms along 1st dimension of h0mat
  ! U^+ H
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, h0vec )
  ! perform lq two-dimensional transforms along 2nd dimension of h0mat
  ! (U^+ H) U
  ! Stride is used to locate the data
  ! e.g. suppose our data has structure data(i,j1,j2), where j1,j2 are the index for 2nd dimension
  !      then the address is i + s1*i + s2*j1 + s3*j2
  Stride(1) = 0; Stride(2) = lq; Stride(3) = latt%l1*lq;
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_NUMBER_OF_TRANSFORMS, lq )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_DISTANCE, 1 )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_DISTANCE, 1 )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_STRIDES, Stride )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_STRIDES, Stride )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim2 )
  Status = DftiComputeBackward( Desc_Handle_Dim2, h0vec )
  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )
  Status = DftiFreeDescriptor( Desc_Handle_Dim2 )

  this%h0mat = reshape( h0vec, (/lq,lq/) )
  do i = 1, lq
      h0k(i) = this%h0mat(i,i)
      this%exph0k(i)    = exp( -dcmplx(0.5d0*dtau,0.d0)*this%h0mat(i,i) )
      this%exph0kinv(i) = exp(  dcmplx(0.5d0*dtau,0.d0)*this%h0mat(i,i) )
  end do
#IFDEF TEST
  write(fout, '(a)') 'after fft, h0mat = '
  do i = 1, lq
      write(fout, '(16f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  if( irank.eq.0 ) then
      write(fout,*)
      write(fout,'(a)') 'fft get eigenvalues = '
      do i = 1, lq
          write(fout,'(e16.8)') dble(h0k(i))
      end do

      en_free = 0.d0
      do i = 1, lq
          if (dble(h0k(i))<=0.d0) then
              en_free = en_free + dble(h0k(i))
          end if
      end do
      write(fout,*)
      write(fout,'(a,e16.8)') ' fft get en_free = ', en_free
  end if

  ! case 2.2, fft for cubic
#ELIF CUBIC
  dimm(:) = (/latt%l1,latt%l2,latt%l3/)
  allocate( this%exph0k(lq) )
  allocate( this%exph0kinv(lq) )
  allocate( h0k(lq) )
#IFDEF TEST
  write(fout, '(a)') 'before fft, h0mat = '
  do i = 1, lq
      write(fout, '(16f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  h0vec = reshape( this%h0mat, (/lq*lq/) )
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 3, dimm )
  Status = DftiCreateDescriptor(Desc_Handle_Dim2, DFTI_DOUBLE, DFTI_COMPLEX, 3, dimm )

  ! perform lq two-dimensional transforms along 1st dimension of h0mat
  ! U^+ H
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, h0vec )
  ! perform lq three-dimensional transforms along 2nd dimension of h0mat
  ! (U^+ H) U
  ! Stride is used to locate the data
  ! e.g. suppose our data have structure data(i,j1,j2,j3), where j1,j2,j3 are the index for 2nd dimension
  !      then the address is i + s1*i + s2*j1 + s3*j2 + s4*j3
  Stride(1) = 0; Stride(2) = lq; Stride(3) = latt%l1*lq; Stride(4) = latt%l1*latt%l2*lq;
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_NUMBER_OF_TRANSFORMS, lq )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_DISTANCE, 1 )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_DISTANCE, 1 )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_STRIDES, Stride )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_STRIDES, Stride )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim2 )
  Status = DftiComputeBackward( Desc_Handle_Dim2, h0vec )
  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )
  Status = DftiFreeDescriptor( Desc_Handle_Dim2 )

  this%h0mat = reshape( h0vec, (/lq,lq/) )
  do i = 1, lq
      h0k(i) = this%h0mat(i,i)
      this%exph0k(i)    = exp( -dcmplx(0.5d0*dtau,0.d0)*this%h0mat(i,i) )
      this%exph0kinv(i) = exp(  dcmplx(0.5d0*dtau,0.d0)*this%h0mat(i,i) )
  end do
#IFDEF TEST
  write(fout, '(a)') 'after fft, h0mat = '
  do i = 1, lq
      write(fout, '(16f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  if( irank.eq.0 ) then
      write(fout,*)
      write(fout,'(a)') 'fft get eigenvalues = '
      do i = 1, lq
          write(fout,'(e16.8)') dble(h0k(i))
      end do

      en_free = 0.d0
      do i = 1, lq
          if (dble(h0k(i))<=0.d0) then
              en_free = en_free + dble(h0k(i))
          end if
      end do
      write(fout,*)
      write(fout,'(a,e16.8)') ' fft get en_free = ', en_free
  end if

  ! case 2.3, fft for honeycomb
#ELIF HONEYCOMB
  ! fft for honeycomb case is more tedious
  dimm(:) = (/latt%l1, latt%l2/)
  allocate( this%exph0k(2,2,lq) )
  allocate( this%exph0kinv(2,2,lq) )
  allocate( h0k(2,2,lq) )
  allocate( htmp(lq*2,lq*2) )
  allocate( htmpC(lq*2,lq*2) )
#IFDEF TEST
  write(fout, '(a)') 'before fft, h0mat = '
  do i = 1, lq*2
      write(fout, '(18f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  ! we first reshape h0mat, such that both the 1st and 2nd dimension have basis:
  ! all sites of sublattice A, all sites of sublattice B
  do j = 1, 2*lq, 2
      do i = 1, 2*lq, 2
          htmp((i+1)/2,    (j+1)/2)    = this%h0mat(i,   j)    ! AA
          htmp((i+1)/2+lq, (j+1)/2+lq) = this%h0mat(i+1, j+1)  ! BB
          htmp((i+1)/2,    (j+1)/2+lq) = this%h0mat(i,   j+1)  ! AB
          htmp((i+1)/2+lq, (j+1)/2)    = this%h0mat(i+1, j)    ! BA
      end do
  end do

  h0vec = reshape( htmp, (/lq*2*lq*2/) )
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm )
  Status = DftiCreateDescriptor(Desc_Handle_Dim2, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm )

  ! perform lq*4 two-dimensional transforms along 1st dimension of htmp
  ! U^+ H
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, lq*4 )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, h0vec )

  ! due to the special storage for honeycomb lattice, to perform fft along 2nd dimension of htmp,
  ! we first perform hermitian conjugate on htmp, then perform fft along 1st dimension of htmp^+
  htmp = reshape( h0vec, (/lq*2,lq*2/) )
  htmpC = dconjg( transpose( htmp ) )
  h0vec = reshape( htmpC, (/lq*2*lq*2/) )
  ! perform lq*4 two-dimensional transforms along 1st dimension of htmp^+
  ! U^+ (U^+ H)^+ = ( U^+ H U )^+
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_NUMBER_OF_TRANSFORMS, lq*4 )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim2, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim2 )
  Status = DftiComputeForward( Desc_Handle_Dim2, h0vec )
  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )
  Status = DftiFreeDescriptor( Desc_Handle_Dim2 )

  htmpC = reshape( h0vec, (/lq*2,lq*2/) )
  htmp = dconjg( transpose( htmpC ) )
  ! reshape back to h0mat
  do i = 1, 2*lq, 2
      this%h0mat(i,  :) = htmp((i+1)/2,    :)
      this%h0mat(i+1,:) = htmp((i+1)/2+lq, :)
  end do
  do j = 1, 2*lq, 2
      do i = 1, 2*lq, 2
          this%h0mat(i,   j)   =  htmp((i+1)/2,    (j+1)/2)    ! AA
          this%h0mat(i+1, j+1) =  htmp((i+1)/2+lq, (j+1)/2+lq) ! BB
          this%h0mat(i,   j+1) =  htmp((i+1)/2,    (j+1)/2+lq) ! AB
          this%h0mat(i+1, j)   =  htmp((i+1)/2+lq, (j+1)/2)    ! BA
      end do
  end do

  allocate( hmat_tmp(2,2), umat_tmp(2,2) )
  do i0 = 1, lq
      h0k(1:2,1:2,i0) = this%h0mat(2*i0-1:2*i0,2*i0-1:2*i0)
      hmat_tmp(1:2,1:2) = this%h0mat(2*i0-1:2*i0,2*i0-1:2*i0)
      call s_eig_he(2,2,hmat_tmp,eig2,umat_tmp)
      eig_tmp(2*i0-1:2*i0) = eig2(1:2)
      do i = 1, 2
         do j = 1, 2
            z0 = czero
            z1 = czero
            do m = 1, 2
               z0 = z0 +  umat_tmp(i,m) * exp(-0.5d0*dtau *eig2(m)) * dconjg(umat_tmp(j,m))
               z1 = z1 +  umat_tmp(i,m) * exp( 0.5d0*dtau *eig2(m)) * dconjg(umat_tmp(j,m))
            enddo
            this%exph0k(i,j,i0)    = z0
            this%exph0kinv(i,j,i0) = z1
         enddo
      enddo
  end do
  deallocate( hmat_tmp, umat_tmp )
#IFDEF TEST
  write(fout, '(a)') 'after fft, h0mat = '
  do i = 1, lq*2
      write(fout, '(18f8.3)') dble(this%h0mat(i,:))
  end do
#ENDIF
  if( irank.eq.0 ) then
      write(fout,*)
      write(fout,'(a)') 'fft get eigenvalues = '
      do i = 1, lq*2
          write(fout,'(e16.8)') eig_tmp(i)
      end do

      en_free = 0.d0
      do i = 1, lq*2
          if (eig_tmp(i)<=0.d0) then
              en_free = en_free + eig_tmp(i)
          end if
      end do
      write(fout,*)
      write(fout,'(a,e16.8)') ' fft get en_free = ', en_free
  end if
  deallocate( htmp, htmpC )

#ENDIF
  deallocate( h0k )
  deallocate( h0vec )
  deallocate( eig_tmp )

! case 3, use zgemm
#ELSE
    this%lq = lq
    this%ltrot = ltrot
    this%rt = rt
    this%mu = mu
    allocate( this%urt(latt%nsites,latt%nsites) )
    allocate( this%urtm1(latt%nsites,latt%nsites) )

    allocate ( this%h0mat(latt%nsites,latt%nsites), umat_tmp(latt%nsites,latt%nsites), eig_tmp(latt%nsites)  )
    
    call seth0(this%h0mat, rt, mu, xmag, flux_x, flux_y, dimer)

    call s_eig_he(latt%nsites,latt%nsites,this%h0mat,eig_tmp,umat_tmp)
          
    do i = 1,latt%nsites
       do j = 1,latt%nsites
          z0 = czero
          z1 = czero
          do m = 1,latt%nsites
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
    do i = 1, latt%nsites
        write(fout,'(18(2e14.6,2x))') this%h0mat(i,:)
    end do
    write(fout,*)
    write(fout,*) ' this%urt(:,:) = '
    do i = 1, latt%nsites
        write(fout,'(18(2e14.6,2x))') this%urt(i,:)
    end do
#ENDIF
    if( irank.eq.0 ) then
        write(fout,*)
        write(fout,'(a)') ' s_eig_he get eigenvalues = '
        do i = 1, latt%nsites
            write(fout,'(e16.8)') eig_tmp(i)
        end do

        en_free = 0.d0
        do i = 1, latt%nsites
            if (eig_tmp(i)<=0.d0) then
                en_free = en_free + eig_tmp(i)
            end if
        end do
        write(fout,*)
        write(fout,'(a,e16.8)') ' s_eig_he get en_free = ', en_free
    end if

    deallocate ( umat_tmp, eig_tmp  )
#ENDIF
  
  end subroutine dqmc_set_h0conf
