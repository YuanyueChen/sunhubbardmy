subroutine dqmc_left_forward_prop(this, nf, gmat)
! B(t+dt,t)*gmat

#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  use dqmc_basic_data
  implicit none
  class(h0conf) :: this
  integer, intent(in) :: nf
  type(gfunc) :: gmat
  integer :: n2
#IFDEF TIMING
  real(dp) :: starttime3, endtime3
#ENDIF
! case 1, use trotter
#IFDEF BREAKUP_T
  complex(dp), allocatable, dimension(:) :: v1, v2
  integer :: i, ist, i1, i2, j
#IFDEF TIMING
  starttime3 = omp_get_wtime()
#ENDIF
  n2 = size(gmat%orb1,2)
  allocate( v1(n2), v2(n2) )
  if (rt.gt.zero) then
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, j, v1, v2 )
!$OMP DO
      do i = 1,latt%nn_lf
         ist = i + (nf - 1)*latt%nn_lf
         i1 = latt%nnlf_list(i) ! A site
         i2 = latt%nnlist(i1,nf)
         do j = 1, n2
            v1(j) = this%urt(ist,1,1) * gmat%orb1(i1,j) + this%urt(ist,1,2) * gmat%orb1(i2,j) 
            v2(j) = this%urt(ist,2,1) * gmat%orb1(i1,j) + this%urt(ist,2,2) * gmat%orb1(i2,j) 
         enddo
         do j = 1, n2
            gmat%orb1(i1,j) = v1(j)
            gmat%orb1(i2,j) = v2(j)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
  endif
  deallocate(v1, v2)

! case 2, use fft
#ELIF DEFINED(FFT)
  ! B G = U D U^+ G

  integer :: Status, i, j
  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
  complex(dp), dimension(:), allocatable :: fftmp

#IFDEF SQUARE
  integer :: dimm(2)
  n2 = size(gmat%orb1,2)
  dimm(:) = (/latt%l1,latt%l2/)
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm)
  allocate( fftmp(lq*n2) )

  ! perform n2 two-dimensional transforms along 1st dimension of gmat
  ! U^+ G
  fftmp = reshape( gmat%orb1, (/lq*n2/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, n2 )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, fftmp )

  gmat%orb1 = reshape( fftmp, (/lq,n2/) )
  do j = 1, n2
      do i = 1, lq
          ! D (U^+ G)
          gmat%orb1(i,j) = gmat%orb1(i,j) * this%exph0k(i)
      end do
  end do

  ! perform n2 two-dimensional transforms along 1st dimension of gmat
  ! U (D U^+ G)
  fftmp = reshape( gmat%orb1, (/lq*n2/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeBackward( Desc_Handle_Dim1, fftmp )

  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )

  gmat%orb1 = reshape( fftmp, (/lq,n2/) )

  deallocate( fftmp )
#ELIF CUBIC
  integer :: dimm(3)
  n2 = size(gmat%orb1,2)
  dimm(:) = (/latt%l1,latt%l2,latt%l3/)
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 3, dimm)
  allocate( fftmp(lq*n2) )

  ! perform n2 three-dimensional transforms along 1st dimension of gmat
  fftmp = reshape( gmat%orb1, (/lq*n2/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, n2 )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, fftmp )

  gmat%orb1 = reshape( fftmp, (/lq,n2/) )
  do j = 1, n2
      do i = 1, lq
          gmat%orb1(i,j) = gmat%orb1(i,j) * this%exph0k(i)
      end do
  end do

  ! perform n2 three-dimensional transforms along 1st dimension of gmat
  fftmp = reshape( gmat%orb1, (/lq*n2/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeBackward( Desc_Handle_Dim1, fftmp )

  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )

  gmat%orb1 = reshape( fftmp, (/lq,n2/) )

  deallocate( fftmp )
#ELIF HONEYCOMB
  integer :: i1, i2, dimm(2)
  complex(dp), dimension(:,:), allocatable :: gtmp
  complex(dp), dimension(:), allocatable :: v1, v2
  n2 = size(gmat%orb1,2)
  dimm(:) = (/latt%l1,latt%l2/)
  Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE, DFTI_COMPLEX, 2, dimm)
  allocate( fftmp(lq*2*n2) )
  allocate( gtmp(lq*2,n2) )
  allocate( v1(n2), v2(n2) )
  ! we first reshape gmat, such that the 1st dimension has basis:
  ! all sites of sublattice A, all sites of sublattice B
  do i = 1, 2*lq, 2
      gtmp((i+1)/2,    :) = gmat%orb1(i,   :)
      gtmp((i+1)/2+lq, :) = gmat%orb1(i+1, :)
  end do

  ! perform 2*n2 two-dimensional transforms along 1st dimension of gmat, factor 2 comes from sublattices
  fftmp = reshape( gtmp, (/lq*2*n2/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, 2*n2 )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, lq )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_FORWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeForward( Desc_Handle_Dim1, fftmp )

  gtmp = reshape( fftmp, (/lq*2,n2/) )
  do i = 1, lq
      i1 = i
      i2 = i + lq
      v1(:) = this%exph0k(1,1,i)*gtmp(i1,:) + this%exph0k(1,2,i)*gtmp(i2,:)
      v2(:) = this%exph0k(2,1,i)*gtmp(i1,:) + this%exph0k(2,2,i)*gtmp(i2,:)
      gtmp(i1,:) = v1(:)
      gtmp(i2,:) = v2(:)
  end do

  ! perform 2*n2 two-dimensional transforms along 1st dimension of gmat, factor 2 comes from sublattices
  fftmp = reshape( gtmp, (/lq*2*n2/) )
  Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1.d0/dsqrt(dble(lq)) )
  Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
  Status = DftiComputeBackward( Desc_Handle_Dim1, fftmp )

  Status = DftiFreeDescriptor( Desc_Handle_Dim1 )

  gtmp = reshape( fftmp, (/lq*2,n2/) )

  ! reshape back to gmat
  do i = 1, 2*lq, 2
      gmat%orb1(i,   :) = gtmp((i+1)/2,    :)
      gmat%orb1(i+1, :) = gtmp((i+1)/2+lq, :)
  end do

  deallocate( fftmp )
  deallocate( gtmp )
  deallocate( v1, v2 )
#ENDIF

! case 3, use zgemm
#ELSE
  type(gfunc) :: Atmp
#IFDEF TIMING
  starttime3 = omp_get_wtime()
#ENDIF
  n2 = size(gmat%orb1,2)
  call allocate_gfunc( Atmp, ndim, n2 )
  if (rt.gt.zero) then
      call zgemm('n','n',ndim,n2,ndim,cone,this%urt,ndim,gmat%orb1,ndim,czero,Atmp%orb1,ndim)
      gmat = Atmp
  endif
#ENDIF
#IFDEF TIMING
  endtime3 = omp_get_wtime()
  timecalculation(8)=timecalculation(8)+endtime3-starttime3
#ENDIF
end subroutine dqmc_left_forward_prop
