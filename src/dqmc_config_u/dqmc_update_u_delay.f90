subroutine dqmc_update_u(this, gmat, ntau )

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
  integer :: iflip, is, isp, isite
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zfunc) :: sscl
  complex(dp), allocatable, dimension(:) :: diagg_up
  complex(dp), allocatable, dimension(:) :: diagg_dn
  complex(dp), allocatable, dimension(:,:) :: avec_up, bvec_up
  complex(dp), allocatable, dimension(:,:) :: avec_dn, bvec_dn
  integer :: i, ik, m
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  allocate( diagg_up(ndim) )
  allocate( avec_up(ndim,nublock) )
  allocate( bvec_up(ndim,nublock) )

  accm  = 0.d0
  ik = 0
  ! initial diag G
  do isite = 1, latt%nsites
      diagg_up(isite) = gmat%orb1(isite,isite)
  end do
  ! intial avec, bvec
  avec_up = czero
  bvec_up = czero
  do isite = 1, latt%nsites
      ! delay update: after nublock steps of local update, perform a whole update of Green function
      ! calculate weight ratio, fermion part
      is = this%conf_u(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      ratio1 = (cone - diagg_up(isite))*this%delta_bmat_u_orb1(iflip, is) + cone
      sscl%orb1 = this%delta_bmat_u_orb1(iflip,is)/ratio1

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

         ik = ik + 1
         ! store avec(:,ik) and bvec(:,ik)
         avec_up(:,ik) = gmat%orb1(:,isite)
         bvec_up(:,ik) = gmat%orb1(isite,:)
         do m = 1, ik-1
             avec_up(:,ik) = avec_up(:,ik) + bvec_up(isite,m)*avec_up(:,m)
             bvec_up(:,ik) = bvec_up(:,ik) + avec_up(isite,m)*bvec_up(:,m)
         end do
         avec_up(:,ik) =avec_up(:,ik)*sscl%orb1
         bvec_up(isite,ik)=bvec_up(isite,ik) - cone
         ! update diag G
         do i = 1, ndim
             diagg_up(i) = diagg_up(i) + avec_up(i,ik)*bvec_up(i,ik)
         end do

         ! flip
         this%conf_u(isite,ntau) =  isp
     end if

      if( (ik.eq.nublock) .or. (isite.eq.latt%nsites) ) then
          ik = 0
          ! delay update: update the whole Green function
          call  zgemm('N', 'T', ndim, ndim, nublock, cone, avec_up, ndim, bvec_up, ndim, cone, gmat%orb1, ndim)
          if( isite.lt.latt%nsites) then
              ! initial diag G
              do i = 1, ndim
              diagg_up(i) = gmat%orb1(i,i)
              end do
              ! intial avec, bvec
              avec_up = czero
              bvec_up = czero
          end if
      end if
   end do
   main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(15)=timecalculation(15)+endtime-starttime
#ENDIF
  deallocate( bvec_up )
  deallocate( avec_up )
  deallocate( diagg_up )
end subroutine dqmc_update_u
