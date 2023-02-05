subroutine dyn_output
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  implicit none

  complex(dp) :: znorm
  character(40) :: ftag
  integer :: i

  include 'mpif.h'
  if( ltau ) then
      znorm = cone / dcmplx( dble(nsweep*isize), 0.d0 )
  end if

  if (irank.eq.0) then
     if(ltau) then
         znorm = cone / dcmplx( dble(nsweep*isize), 0.d0 )
     end if
  endif
end subroutine dyn_output
