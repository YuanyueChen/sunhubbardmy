subroutine dqmc_left_forward_prop(this, gmat, ntau, nf, nflag)
  ! B(t+dt,t)*gmat
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  implicit none

  class(vconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau, nf, nflag

  ! local
  integer :: i, j, i1, i2, n2, is, isite
  complex(dp), allocatable, dimension(:) :: v1, v2
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  n2 = size(gmat%blk1,2)
  allocate( v1(n2), v2(n2) )

  ! B = U D U^+
  ! B G = U D U^+ G
  if (nflag.eq.2) then
      ! D U^+ G
!$OMP PARALLEL &
!$OMP PRIVATE ( i, is, i1, i2, j, v1, v2)
!$OMP DO
     do i = 1, latt%nn_lf
        i1 = latt%nnlf_list(i)
        i2 = latt%nnlist(i1,nf)
        do j = 1, n2
           v1(j) = this%ut(1,1) * gmat%blk1(i1,j) + this%ut(1,2) * gmat%blk1(i2,j)
           v2(j) = this%ut(2,1) * gmat%blk1(i1,j) + this%ut(2,2) * gmat%blk1(i2,j) 
        enddo
        is = this%conf(i,nf,ntau)
        do j = 1, n2
           gmat%blk1(i1,j) = this%bmat_p%blk1(is)*v1(j)
           gmat%blk1(i2,j) = this%bmat_m%blk1(is)*v2(j)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif

  if (nflag.eq.1) then
      ! U (D U^+ G)
!$OMP PARALLEL &
!$OMP PRIVATE ( i, i1, i2, j, v1, v2)
!$OMP DO
     do i = 1, latt%nn_lf
        i1 = latt%nnlf_list(i)
        i2 = latt%nnlist(i1,nf)
        do j = 1, n2
           v1(j) = this%u(1,1) * gmat%blk1(i1,j) +  this%u(1,2) * gmat%blk1(i2,j)
           v2(j) = this%u(2,1) * gmat%blk1(i1,j) +  this%u(2,2) * gmat%blk1(i2,j)
        enddo
        do j = 1, n2
           gmat%blk1(i1,j) = v1(j)
           gmat%blk1(i2,j) = v2(j)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif

#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(10)=timecalculation(10)+endtime-starttime
#ENDIF
  deallocate( v2, v1)
end subroutine dqmc_left_forward_prop
