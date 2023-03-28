subroutine dqmc_proj_update_u(this, ntau, ul, ur, ulrinv)

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
  type(gfunc) :: ul, ur, ulrinv

  !local
  complex(dp) ::  ratio1, ratiotot
  integer :: iflip, is, isp, isite
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zfunc) :: sscl
  complex(dp), allocatable, dimension(:) :: diagg_up
  complex(dp), allocatable, dimension(:) :: diagg_dn
  complex(dp), allocatable, dimension(:,:) :: avec_up, bvec_up
  complex(dp), allocatable, dimension(:,:) :: avec_dn, bvec_dn
  complex(dp), allocatable, dimension(:,:) :: gmmat, ulrinvul, urrecord
  integer :: i, ik, m
#IFDEF TIMING
  real(dp) :: starttime26, endtime26
#ENDIF
#IFDEF TIMING
  starttime26 = omp_get_wtime()
#ENDIF

  allocate( diagg_up(ndim) )
  allocate( avec_up(ndim,ndim) )
  allocate( bvec_up(ndim,ndim) )
  allocate( gmmat(ndim,ndim) )
  allocate( ulrinvul(ne,ndim) )
  allocate( urrecord(ndim,ne) )
  
  ulrinvul = czero
  gmmat = czero
  urrecord = ur%orb1
  !calculate the G function
  call zgemm('N', 'N', ne, ndim, ne, cone, ulrinv%orb1, ne, ul%orb1, ne, czero, ulrinvul, ne)
  call zgemm('N', 'N', ndim, ndim, ne, cone, ur%orb1, ndim, ulrinvul,ne, czero, gmmat, ndim)
  gmmat(:,:) =-gmmat(:,:)
  do isite = 1, latt%nsites
      gmmat(isite,isite) =cone + gmmat(isite,isite)
  end do
  accm  = 0.d0
  ik = 0
  ! initial diag G
  do isite = 1, latt%nsites
      diagg_up(isite) =gmmat(isite,isite)
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
         avec_up(:,ik) = gmmat(:,isite)
         bvec_up(:,ik) = gmmat(isite,:)
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
         urrecord(isite,:) = urrecord(isite,:) + this%delta_bmat_u_orb1(iflip,is)*urrecord(isite,:)
         ! flip
         this%conf_u(isite,ntau) =  isp
     end if
     
     if(  isite.eq.latt%nsites ) then
          ik = 0
          ! delay update: update the R and the (LR)^-1
          ur%orb1 = urrecord
          call zgemm('N', 'N', ne, ne, ndim, cone, ul%orb1, ne, ur%orb1, ndim, czero, ulrinv%orb1, ne)
          call s_invlu_z(ne, ulrinv%orb1) 
      end if
  end do
  main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )
#IFDEF TIMING
  endtime26 = omp_get_wtime()
  timecalculation(17)=timecalculation(17)+endtime26-starttime26
#ENDIF
  deallocate( bvec_up )
  deallocate( avec_up )
  deallocate( diagg_up )
  deallocate( ulrinvul )
  deallocate( gmmat )
  deallocate( urrecord )
end subroutine dqmc_proj_update_u
