subroutine dqmc_update(this, gmat, ntau )

  use spring
  use model_para
  use dqmc_basic_data
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  implicit none
  !arguments
  class(uconf), intent(inout) :: this
  integer,intent(in) :: ntau
  type(gfunc) :: gmat

  !local
  complex(dp) ::  ratio1, ratiotot
  integer :: iflip, is, i0, k, isp, isite, i1, i2
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zvfunc) :: ukmat, vkmat
  type(zfunc) :: sscl
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  call allocate_zvfunc( ukmat, ndim)
  call allocate_zvfunc( vkmat, ndim)

  accm  = 0.d0
  do isite = 1, latt%nsites
      is = this%conf(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      ukmat%orb1 = czero
      vkmat%orb1 = czero
      vkmat%orb1(:) =      - gmat%orb1(isite, :)
      vkmat%orb1(isite) = cone - gmat%orb1(isite, isite)

      ratio1 = vkmat%orb1(isite)*this%delta_bmat%orb1(iflip, is) + cone
      sscl%orb1 = this%delta_bmat%orb1(iflip,is)/ratio1

      ratiotot = ratio1**nflr
      if( lproju ) then
          ratiotot = ratiotot*this%phase_ratio(iflip, is)
      else
          ratiotot = ratiotot*dcmplx( gaml(isp)/gaml(is), 0.d0 ) * this%phase_ratio(iflip, is)
      end if
#IFDEF TEST
      write(fout, '(a,2e16.8)') ' in update_u, ratio1 = ', ratio1
      write(fout, '(a,2e16.8)') ' in update_u, phase_ratio = ', this%phase_ratio(iflip, is)
      write(fout, '(a,2e16.8)') ' in update_u, ratiotot = ', ratiotot
#ENDIF
      
	  ratio_re = dble( ratiotot * phase )/dble( phase )

      ratio_re_abs = ratio_re
      if ( ratio_re .lt. 0.d0 ) ratio_re_abs = - ratio_re

      random = spring_sfmt_stream()

      if ( ratio_re_abs.gt.random ) then

         accm  = accm + 1.d0
         !weight_track = weight_track + log( ratio_re_abs )
	     weight = dsqrt(dble(ratiotot*dconjg(ratiotot)))
	     phase =  phase*ratiotot/dcmplx(weight,0.d0)
#IFDEF TEST
         write(fout, '(a,2e16.8)') ' in update_u, phase = ', phase
#ENDIF

         ! update gmat
         ukmat%orb1(:) = gmat%orb1(:,isite)
         do i2 = 1, ndim
         do i1 = 1, ndim
             gmat%orb1(i1,i2) = gmat%orb1(i1,i2) - sscl%orb1*ukmat%orb1(i1)*vkmat%orb1(i2)
         end do
         end do

         ! flip
         this%conf(isite,ntau) =  isp
      endif
   end do
   main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(12)=timecalculation(12)+endtime-starttime
#ENDIF
end subroutine dqmc_update
