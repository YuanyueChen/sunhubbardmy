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
  type(gfunc) :: gmat
  integer :: n1
#IFDEF TIMING
  real(dp) :: starttime6, endtime6
#ENDIF
#IFDEF BREAKUP_T
  complex(dp), allocatable, dimension(:) :: v1, v2
  integer :: i, ist, i1, i2, j
#IFDEF TIMING
  starttime6 = omp_get_wtime()
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
#ELSE
  type(gfunc) :: Atmp
#IFDEF TIMING
  starttime6 = omp_get_wtime()
#ENDIF
  n1 = size(gmat%orb1,1)
  call allocate_gfunc( Atmp, n1, ndim )
  if (rt.gt.zero) then
      call zgemm('n','n',n1,ndim,ndim,cone,gmat%orb1,n1,this%urt,ndim,czero,Atmp%orb1,n1)
      gmat = Atmp
  endif
#ENDIF
#IFDEF TIMING
  endtime6 = omp_get_wtime()
  timecalculation(9)=timecalculation(9)+endtime6-starttime6
#ENDIF

end subroutine dqmc_right_forward_prop
