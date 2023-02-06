subroutine dqmc_proj_update_u(this, ntau, ul, ur, ulrinv)

  use spring
  use model_para
  use dqmc_basic_data

  implicit none

  !arguments
  class(uconf), intent(inout) :: this
  integer,intent(in) :: ntau
  type(gfunc) :: ul, ur, ulrinv

  !local
  complex(dp) ::  ratio1, ratiotot
  integer :: iflip, is, i0, k, isp, isite, i1, i2
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zvfunc) :: ukmat, vkmat, vuvkmat, lrukmat, vkmat_tmp
  type(zfunc) :: sscl

  call allocate_zvfunc( ukmat, ne)
  call allocate_zvfunc( vkmat, ne)
  call allocate_zvfunc( vuvkmat, ne)
  call allocate_zvfunc( lrukmat, ne )
  call allocate_zvfunc( vkmat_tmp, ne)

  accm  = 0.d0
  do isite = 1, latt%nsites
      is = this%conf_u(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      vuvkmat%orb1(:) = ur%orb1(isite,:)
      ukmat%orb1(:) = ul%orb1(:,isite)

      vkmat_tmp%orb1(:) = this%delta_bmat_u_orb1(iflip,is)*vuvkmat%orb1(:)

      vkmat%orb1 = czero
      do i1 = 1, ne
          do i2 = 1, ne
              vkmat%orb1(i1) = vkmat%orb1(i1) + vkmat_tmp%orb1(i2)*ulrinv%orb1(i2,i1)
          end do
      end do

      ratio1 = cone
      do i1 = 1, ne
          ratio1 = ratio1 + vkmat%orb1(i1)*ukmat%orb1(i1)
      end do

      sscl%orb1 = cone / ratio1


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

         ! update ur
         ur%orb1(isite,:) = ur%orb1(isite,:) + vkmat_tmp%orb1(:)

         ! update urlinv
         vuvkmat%orb1(:) = sscl%orb1*vkmat%orb1(:)
         lrukmat%orb1 = czero
         do i2 = 1, ne
             do i1 = 1, ne
                 lrukmat%orb1(i1) = lrukmat%orb1(i1) + ulrinv%orb1(i1,i2)*ukmat%orb1(i2)
             end do
         end do

         do i1 = 1, ne
             do i2 = 1, ne
                 ulrinv%orb1(i1,i2) = ulrinv%orb1(i1,i2) - lrukmat%orb1(i1)*vuvkmat%orb1(i2)
             end do
         end do


         ! flip
         this%conf_u(isite,ntau) =  isp
     end if
  end do
  main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )
end subroutine dqmc_proj_update_u
