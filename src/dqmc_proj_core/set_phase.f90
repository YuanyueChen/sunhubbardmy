subroutine set_phase( logweight, w_phase )
  implicit none
  type(zfunc), intent(in) :: logweight
  complex(dp), intent(out) :: w_phase
  
  ! local 
  integer :: i, nt
  w_phase = exp( dcmplx(0.0_dp, dble(nflr)*aimag(logweight%orb1) ) )
  if( lwrapplqu ) then
      do nt = 1, ltrot
          do i = 1, latt%ncell
              w_phase = w_phase*hconf%phase( hconf%conf_plqu(i,nt) )
          end do
      end do
  end if
  if( lwrapu ) then
      do nt = 1, ltrot
          do i = 1, latt%nsites
              w_phase = w_phase*u0conf%phase( u0conf%conf_u(i,nt) )
          end do
      end do
  end if

end subroutine set_phase
