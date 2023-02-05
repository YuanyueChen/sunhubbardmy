subroutine dqmc_left_forward_prop(this, gmat, ntau)
  ! B(t+dt,t)*gmat
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
  n2 = size(gmat%orb1,2)

!$OMP PARALLEL &
!$OMP PRIVATE ( isite, is, i )
!$OMP DO
  do isite = 1, this%nsites
      is = this%conf_u(isite,ntau)
      do i = 1, n2
          gmat%orb1(isite, i) = this%bmat_u_orb1(is)*gmat%orb1(isite, i)
      end do
  end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine dqmc_left_forward_prop
