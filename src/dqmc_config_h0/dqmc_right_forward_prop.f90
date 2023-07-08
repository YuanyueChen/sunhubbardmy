subroutine dqmc_right_forward_prop(this, nf, gmat)
! gmat*B(t+dt,t)

#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  use dqmc_basic_data
  implicit none
  class(h0conf) :: this
  integer, intent(in) :: nf
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF 
  type(gfunc) :: gmat
  integer :: n1
! case 1, use trotter
#IFDEF BREAKUP_T
  complex(dp), allocatable, dimension(:) :: v1, v2
  integer :: i, ist, i1, i2, j
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  n1 = size(gmat%orb1,1)
  allocate( v1(n1), v2(n1) )
  if (rt.gt.zero) then
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, j, v1, v2 )
!$OMP DO
      do i = 1,latt%nn_lf
         ist = i + (nf - 1)*latt%nn_lf
         i1 = latt%nnlf_list(i) ! A site
         i2 = latt%nnlist(i1,nf)
         do j = 1, n1
            v1(j) = gmat%orb1(j,i1)*this%urt(ist,1,1) + gmat%orb1(j,i2)*this%urt(ist,2,1)
            v2(j) = gmat%orb1(j,i1)*this%urt(ist,1,2) + gmat%orb1(j,i2)*this%urt(ist,2,2) 
         enddo
         do j = 1, n1
            gmat%orb1(j,i1) = v1(j)
            gmat%orb1(j,i2) = v2(j)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
  endif
  deallocate(v1, v2)

! case 2, use fft
#ELIF DEFINED(FFT)
  ! G B = ( (G B)^+ )^+
  ! (G B)^+ = B^+ G^+ = U D^+ U^+ G^+

  integer :: Status, i, j
  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
  complex(dp), dimension(:), allocatable :: fftmp
  complex(dp), dimension(:,:), allocatable :: gtmpC

#IFDEF SQUARE
  integer :: dimm(2)
#IFDEF TIMING 
  call cpu_time_now(starttime)
#ENDIF
  n1 = size(gmat%orb1,1)
  dimm(:) = (/latt%l1,latt%l2/)
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm)
  allocate( fftmp(lq*n1) )

  ! to perform fft transforms along 2nd dimension of gmat, we can first hermitian conjugate it
  ! and perform fft transforms along 1st dimension of gmat^+
  ! namely to calcualte G B = G U D U^+, we calculate ( G B )^+ = B^+ G^+ = U D^+ U^+ G^+
  allocate( gtmpC(lq,n1) )
  gtmpC = dconjg( transpose( gmat%orb1 ) )

  ! perform n1 two-dimensional transforms along 1st dimension of gmat^+
  ! U^+ G^+
  fftmp = reshape( gtmpC, (/lq*n1/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, n1 )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, fftmp )

  gtmpC = reshape( fftmp, (/lq,n1/) )
  do j = 1, n1
      do i = 1, lq
          ! D^+ (U^+ G^+)
          gtmpC(i,j) = gtmpC(i,j) * dconjg( this%exph0k(i) )
      end do
  end do

  ! perform n1 two-dimensional transforms along 1st dimension of gmat^+
  ! U (D^+ U^+ G^+)
  fftmp = reshape( gtmpC, (/lq*n1/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeBackward( Desc_Handle_Dim1, fftmp )

  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )

  gtmpC = reshape( fftmp, (/lq,n1/) )
  gmat%orb1 = dconjg( transpose( gtmpC ) )

  deallocate( fftmp )
  deallocate( gtmpC )
#ELIF CUBIC
  integer :: dimm(3)
#IFDEF TIMING 
  call cpu_time_now(starttime)
#ENDIF
  n1 = size(gmat%orb1,1)
  dimm(:) = (/latt%l1,latt%l2,latt%l3/)
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 3, dimm)
  allocate( fftmp(lq*n1) )

  ! to perform fft transforms along 2nd dimension of gmat, we can first hermitian conjugate it
  ! and perform fft transforms along 1st dimension of gmat^+
  allocate( gtmpC(lq,n1) )
  gtmpC = dconjg( transpose( gmat%orb1 ) )

  ! perform n1 three-dimensional transforms along 1st dimension of gmat^+
  fftmp = reshape( gtmpC, (/lq*n1/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, n1 )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, fftmp )

  gtmpC = reshape( fftmp, (/lq,n1/) )
  do j = 1, n1
      do i = 1, lq
          gtmpC(i,j) = gtmpC(i,j) * dconjg( this%exph0k(i) )
      end do
  end do

  ! perform n1 three-dimensional transforms along 1st dimension of gmat^+
  fftmp = reshape( gtmpC, (/lq*n1/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeBackward( Desc_Handle_Dim1, fftmp )

  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )

  gtmpC = reshape( fftmp, (/lq,n1/) )
  gmat%orb1 = dconjg( transpose( gtmpC ) )

  deallocate( fftmp )
  deallocate( gtmpC )
#ELIF HONEYCOMB
  integer :: i1, i2, dimm(2)
  complex(dp), dimension(:,:), allocatable :: gtmp
  complex(dp), dimension(:), allocatable :: v1, v2
  n1 = size(gmat%orb1,1)
#IFDEF TIMING 
  call cpu_time_now(starttime)
#ENDIF
  dimm(:) = (/latt%l1,latt%l2/)
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm)
  allocate( fftmp(lq*2*n1) )
  allocate( v1(n1), v2(n1) )
  ! we first reshape gmat, such that the 2nd dimension has basis:
  ! all sites of sublattice A, all sites of sublattice B
  allocate( gtmp(n1,lq*2) )
  do i = 1, 2*lq, 2
      gtmp(:, (i+1)/2   ) = gmat%orb1(:, i  )
      gtmp(:, (i+1)/2+lq) = gmat%orb1(:, i+1)
  end do

  ! to perform fft transforms along 2nd dimension of gmat, we can first hermitian conjugate it
  ! and perform fft transforms along 1st dimension of gmat^+
  allocate( gtmpC(lq*2,n1) )
  gtmpC = dconjg( transpose( gtmp ) )

  ! perform 2*n1 two-dimensional transforms along 1st dimension of gmat^+, factor 2 comes from sublattices
  fftmp = reshape( gtmpC, (/lq*2*n1/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, 2*n1 )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, fftmp )

  gtmpC = reshape( fftmp, (/lq*2,n1/) )
  do i = 1, lq
      i1 = i
      i2 = i + lq
      ! D^+ (U^+ G^+)
      v1(:) = dconjg( this%exph0k(1,1,i) )*gtmpC(i1,:) + dconjg( this%exph0k(2,1,i) )*gtmpC(i2,:)
      v2(:) = dconjg( this%exph0k(1,2,i) )*gtmpC(i1,:) + dconjg( this%exph0k(2,2,i) )*gtmpC(i2,:)
      gtmpC(i1,:) = v1(:)
      gtmpC(i2,:) = v2(:)
  end do

  ! perform 2*n1 two-dimensional transforms along 1st dimension of gmat^+, factor 2 comes from sublattices
  fftmp = reshape( gtmpC, (/lq*2*n1/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeBackward( Desc_Handle_Dim1, fftmp )

  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )

  gtmpC = reshape( fftmp, (/lq*2,n1/) )
  gtmp = dconjg( transpose( gtmpC ) )

  ! reshape back to gmat
  do i = 1, 2*lq, 2
      gmat%orb1(:, i  ) = gtmp(:, (i+1)/2   )
      gmat%orb1(:, i+1) = gtmp(:, (i+1)/2+lq)
  end do

  deallocate( fftmp )
  deallocate( gtmp )
  deallocate( gtmpC )
  deallocate( v1, v2 )
#ENDIF

! case 3, use zgemm
#ELSE
  type(gfunc) :: Atmp
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  n1 = size(gmat%orb1,1)
  call allocate_gfunc( Atmp, n1, ndim )
  if (rt.gt.zero) then
      call zgemm('n','n',n1,ndim,ndim,cone,gmat%orb1,n1,this%urt,ndim,czero,Atmp%orb1,n1)
      gmat = Atmp
  endif
#ENDIF
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(9)=timecalculation(9)+endtime-starttime
#ENDIF

end subroutine dqmc_right_forward_prop
