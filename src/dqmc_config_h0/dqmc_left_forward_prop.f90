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
