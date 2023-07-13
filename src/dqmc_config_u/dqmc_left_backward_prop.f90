subroutine dqmc_left_backward_prop(this, gmat, ntau )
  ! B(t+dt,t)^-1 *gmat
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  implicit none

  class(uconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau

  ! local
  integer :: i, n2, is, isite
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  n2 = size(gmat%orb1,2)

!$OMP PARALLEL &
!$OMP PRIVATE ( isite, is, i )
!$OMP DO
  do isite = 1, this%nsites
      is = this%conf(isite,ntau)
      do i = 1, n2
          gmat%orb1(isite, i) = this%bmat_inv%orb1(is)*gmat%orb1(isite, i)
      end do
  end do
!$OMP END DO
!$OMP END PARALLEL
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(10)=timecalculation(10)+endtime-starttime
#ENDIF

end subroutine dqmc_left_backward_prop
