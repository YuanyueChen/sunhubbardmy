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
  complex(dp) ::  ratio1, ratiotot, v, v2, v3
  integer :: iflip, is, isp, isite
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zfunc) :: sscl
  complex(dp), allocatable, dimension(:,:) ::   gammainv_up,  Gcuta_up, gammainvtmp_up, gmattmp1_up, gmattmp2_up
  complex(dp), allocatable, dimension(:,:) ::   gammainv_dn,  Gcuta_dn, gammainvtmp_dn, gmattmp1_dn, gmattmp2_dn
  complex(dp), allocatable, dimension(:,:) ::   Gmat1_up, Gmat2_up
  complex(dp), allocatable, dimension(:,:) ::   Gmat1_dn, Gmat2_dn
  complex(dp), allocatable, dimension(:) ::  Qmat_up, Qmat_dn, v1, v4, cvec_up, cvec_dn, bvec_up, bvec_dn
  complex(dp), allocatable, dimension(:) ::  Gtmp_up, Gtmp_dn
  integer, allocatable, dimension(:) :: pvec_up, pvec_dn
  integer :: i, ik, m, isf
#IFDEF TIMING
  real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  allocate( pvec_up(nublock) )
  allocate( gammainv_up(nublock,nublock) )
  allocate( Gcuta_up(ndim,ndim) )
  allocate(cvec_up(nublock))
  allocate(bvec_up(nublock))
  allocate(v1(nublock))
  allocate(v4(nublock))
  allocate( gmattmp1_up(ndim,nublock) )
  allocate( gmattmp2_up(nublock,ndim) )
  allocate(Gmat1_up(ndim,nublock))


    ! intial pvec,svec,wvec
  accm  = 0.d0
  ik = 0
  pvec_up = 0
  cvec_up = czero
  bvec_up = czero
  v2 = czero
  v3 = czero
  Gcuta_up = gmat%blk1
  gmattmp1_up = czero
  gmattmp2_up = czero
  do isite = 1, latt%nsites
      ! submatrix update: after nublock steps of local update, perform a whole update of Green function
      ! calculate weight ratio, fermion part

      is = this%conf(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1
      ! intial svec, wvec
      cvec_up = czero
      bvec_up = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( isf, i )
!$OMP DO
      do i = 1, ik
        isf = pvec_up(i) 
        cvec_up(i) = gmat%blk1(isite,isf)
        bvec_up(i) = gmat%blk1(isf,isite)
      end do
!$OMP END DO
!$OMP END PARALLEL
      v1 = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( m, i )
!$OMP DO
      do i = 1, ik
        do m = 1, ik
          v1(i) = v1(i) + cvec_up(m)*gammainv_up(m,i)
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      v = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( i )
!$OMP DO
      do i = 1, ik
        v = v + v1(i)*bvec_up(i)
      end do
!$OMP END DO
!$OMP END PARALLEL
      v2 = v + gmat%blk1(isite,isite) - cone 
      ratio1 = -v2 *this%delta_bmat%blk1(iflip, is) + cone
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
         ! store pvec, Qmat, gammainv
         pvec_up(ik) = isite
         v3 = cone/(-v + cone + cone/this%delta_bmat%blk1(iflip, is)-gmat%blk1(isite,isite))
         gammainv_up(ik,ik) = v3
         v4 = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( m, i )
!$OMP DO
          do i = 1, ik-1
            do m = 1, ik-1
              v4(i) = v4(i) + gammainv_up(i,m)*bvec_up(m)
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE ( i )
!$OMP DO
          do i = 1, ik-1
            gammainv_up(ik,i) = v3*v1(i)
            gammainv_up(i,ik) = v4(i)*v3
          end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE ( m, i )
!$OMP DO
          do m = 1, ik-1
            do i = 1, ik-1
              gammainv_up(i,m) = gammainv_up(i,m) + v4(i)*v3*v1(m)
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
          Gcuta_up(:,isite) = Gcuta_up(:,isite) - this%delta_bmat%blk1(iflip, is)/(cone + this%delta_bmat%blk1(iflip, is))*Gcuta_up(:,isite)
          gmattmp1_up(:,ik) = gmat%blk1(:,isite)
         ! flip
         this%conf(isite,ntau) =  isp
     end if

      if( (ik.eq.nublock) .or. (isite.eq.latt%nsites) ) then
          ! submatrix update: update the whole Green function
#IFDEF TIMING
          call cpu_time_now(starttime11)
#ENDIF
!$OMP PARALLEL &
!$OMP PRIVATE ( isf, i )
!$OMP DO
          do i = 1, ik
            isf = pvec_up(i) 
            gmattmp2_up(i,:) = Gcuta_up(isf,:)
          end do
!$OMP END DO
!$OMP END PARALLEL
          call zgemm('N','N', ndim, nublock, ik, cone, gmattmp1_up, ndim, gammainv_up, nublock, czero, Gmat1_up, ndim)
          call zgemm('N','N', ndim, ndim, ik, cone, Gmat1_up, ndim, gmattmp2_up, nublock, cone, Gcuta_up, ndim)
          gmat%blk1 = Gcuta_up
          ik = 0
#IFDEF TIMING
          call cpu_time_now(endtime11)
          timecalculation(18)=timecalculation(18)+endtime11-starttime11
#ENDIF
      end if
   end do
   main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )
  deallocate( gmattmp1_up )
  deallocate( gmattmp2_up )
  deallocate( pvec_up )
  deallocate( Gcuta_up)
  deallocate(cvec_up)
  deallocate(bvec_up)
  deallocate(v1)
  deallocate(v4)
  deallocate( gammainv_up )
  deallocate( Gmat1_up )

#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(16)=timecalculation(16)+endtime-starttime
#ENDIF
end subroutine dqmc_update
