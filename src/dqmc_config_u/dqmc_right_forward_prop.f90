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

end subroutine dqmc_right_forward_prop
