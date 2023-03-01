subroutine dqmc_right_forward_prop(this, gmat, ntau)
  ! gmat*B(t+dt,t)
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  implicit none

  class(uconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau

  ! local
  integer :: i, n1, is, jsite
#IFDEF TIMING
  real(dp) :: starttime11, endtime11
#ENDIF
#IFDEF TIMING
  starttime11 = omp_get_wtime()
#ENDIF
  n1 = size(gmat%orb1,1)

!$OMP PARALLEL &
!$OMP PRIVATE ( jsite, is, i )
!$OMP DO
  do jsite = 1, this%nsites
      is = this%conf_u(jsite,ntau)
      do i = 1, n1
          gmat%orb1(i, jsite) = gmat%orb1(i, jsite)*this%bmat_u_orb1(is)
      end do
  end do
!$OMP END DO
!$OMP END PARALLEL
#IFDEF TIMING
      endtime11 = omp_get_wtime()
      timecalculation(11)=timecalculation(11)+endtime11-starttime11
#ENDIF 
end subroutine dqmc_right_forward_prop
