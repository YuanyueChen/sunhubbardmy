subroutine dqmc_proj_update(this, ntau, ul, ur, ulrinv)

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
  complex(dp) ::  ratio1, ratiotot, v, v2, v3
  integer :: iflip, is, isp, isite
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zfunc) :: sscl
  complex(dp), allocatable, dimension(:,:) :: gammainv_up,  Gcuta_up, gammainvtmp_up, gmattmp1_up, gmattmp2_up
  complex(dp), allocatable, dimension(:,:) :: gammainv_dn,  Gcuta_dn, gammainvtmp_dn, gmattmp1_dn, gmattmp2_dn
  complex(dp), allocatable, dimension(:,:) :: Gmat1_up, Gmat2_up
  complex(dp), allocatable, dimension(:,:) :: Gmat1_dn, Gmat2_dn
  complex(dp), allocatable, dimension(:) ::  Qmat_up, Qmat_dn, v1, v4, cvec_up, cvec_dn, bvec_up, bvec_dn
  complex(dp), allocatable, dimension(:) ::  Gtmp_up, Gtmp_dn 
  integer, allocatable, dimension(:) :: pvec_up, pvec_dn
  type(gfunc) :: Fmat_cuta, Fmat, ulrinvul, urrecord,urFcuta_up, urFcuta_dn, uFcutal_up, ulFcutal_dn
  integer :: i, ik, m, isf
#IFDEF TIMING
  real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  allocate( pvec_up(nublock) )
  allocate( gammainv_up(nublock,nublock) )
  allocate(cvec_up(nublock))
  allocate(bvec_up(nublock))
  allocate(v1(nublock))
  allocate(v4(nublock))
  allocate( gmattmp1_up(ndim,nublock) )
  allocate( gmattmp2_up(nublock,ndim) )
  allocate( Gmat1_up(ne,nublock) )
  call allocate_gfunc(urFcuta_up, ndim, ne)
  call allocate_gfunc(uFcutal_up, ne, ndim)
  call allocate_gfunc(Fmat_cuta, ne, ne)
  call allocate_gfunc(Fmat, ndim, ndim)
  call allocate_gfunc(ulrinvul, ne, ndim)
  call allocate_gfunc(urrecord, ndim, ne)

  
  ulrinvul%blk1 = czero
  urrecord%blk1 = ur%blk1
  Fmat_cuta%blk1 = ulrinv%blk1
  accm  = 0.d0
  ik = 0
  pvec_up = 0
  v2 = czero
  v3 = czero
  gammainv_up = czero
  call zgemm('N', 'N', ne, ndim, ne, cone, Fmat_cuta%blk1, ne, ul%blk1, ne, czero, ulrinvul%blk1, ne)
  call zgemm('N','N', ndim, ne, ne, cone, ur%blk1, ndim, Fmat_cuta%blk1, ne, czero, urFcuta_up%blk1, ndim)
  call zgemm('N', 'N', ndim, ndim, ne, cone, ur%blk1, ndim, ulrinvul%blk1, ne, czero, Fmat%blk1, ndim)
  do isite = 1, latt%nsites
      ! delay update: after nublock steps of local update, perform a whole update of Green function
      ! calculate weight ratio, fermion part
      is = this%conf(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1
      !calculate the F function
      cvec_up = czero
      bvec_up = czero
      do i = 1, ik
        isf = pvec_up(i)
        cvec_up(i) = Fmat%blk1(isite,isf)
        bvec_up(i) = Fmat%blk1(isf,isite)
      end do
      v1 = czero
      do i = 1, ik
        do m = 1, ik
          v1(i) = v1(i) + cvec_up(m)*gammainv_up(m,i)
        end do
      end do
      v = czero
      do i = 1, ik
        v = v + v1(i)*bvec_up(i)
      end do
      v2 = Fmat%blk1(isite,isite) - v 
      ratio1 = v2 *this%delta_bmat%blk1(iflip, is) + cone
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
         v3 = cone/(-v  + cone/this%delta_bmat%blk1(iflip, is)+Fmat%blk1(isite,isite))
         gammainv_up(ik,ik) = v3
         v4 = czero
          do i = 1, ik-1
            do m = 1, ik-1
              v4(i) = v4(i) + gammainv_up(i,m)*bvec_up(m)
            end do
          end do
          do i = 1, ik-1
            gammainv_up(ik,i) = -v3*v1(i)
            gammainv_up(i,ik) = -v4(i)*v3
          end do
          do i = 1, ik-1
            do m = 1, ik-1
              gammainv_up(i,m) = gammainv_up(i,m) + v4(i)*v3*v1(m)
            end do
          end do
          urrecord%blk1(isite,:) = urrecord%blk1(isite,:) + this%delta_bmat%blk1(iflip,is)*urrecord%blk1(isite,:)
          gmattmp1_up(:,ik) = ulrinvul%blk1(:,isite)
          gmattmp2_up(ik,:) = urFcuta_up%blk1(isite,:)
         ! flip
         this%conf(isite,ntau) =  isp
     end if
      
     if( (ik.eq.nublock) .and. (isite.lt.latt%nsites) ) then
#IFDEF TIMING
          call cpu_time_now(starttime11)
#ENDIF 

          !call zgemm('N','N', ne, ndim, ne, cone, Fmat_cuta%blk1, ne, ul%blk1, ne, czero, uFcutal_up%blk1, ne)
          call zgemm('N','N', ne, nublock, ik, cone, gmattmp1_up, ne, gammainv_up, nublock, czero, Gmat1_up, ne)
          call zgemm('N','N', ne, ne, ik, -cone, Gmat1_up, ne, gmattmp2_up, nublock, cone, Fmat_cuta%blk1, ne)
          call zgemm('N', 'N', ne, ndim, ne, cone, Fmat_cuta%blk1, ne, ul%blk1, ne, czero, ulrinvul%blk1, ne)
          call zgemm('N', 'N', ndim, ndim, ne, cone, ur%blk1, ndim, ulrinvul%blk1, ne, czero, Fmat%blk1, ndim)
#IFDEF TEST
         write(fout, '(a,2e16.8)') ' after update the G, Fmat = ' 
         do i = 1, ndim
           write(fout, '(18(2e12.4))')  Fmat%blk1(i,:)
         end do
         write(fout, '(a,2e16.8)') ' after update the G, gmattmp1 = '
          do i = 1, ndim
            write(fout, '(18(2e12.4))')  gmattmp1_up(i,:)
          end do
          write(fout, '(a,2e16.8)') ' after update the G, gmattmp2 = '
          do i = 1, nublock
            write(fout, '(18(2e12.4))')  gmattmp2_up(i,:)
          end do
         write(fout, '(a,2e16.8)') ' after update the G, pvec = '
          write(fout, '(18(2i4))')  pvec_up
#ENDIF
          ik = 0
          gammainv_up = czero
#IFDEF TIMING
          call cpu_time_now(endtime11)
          timecalculation(20)=timecalculation(20)+endtime11-starttime11
#ENDIF
      end if

     if(  isite.eq.latt%nsites ) then
#IFDEF TIMING
          call cpu_time_now(starttime11)
#ENDIF 
#IFDEF TEST
         write(fout, '(a,2e16.8)') 'update the R ' 
#ENDIF
          ik = 0
          ! delay update: update the R and the (LR)^-1
          ! update the R
          ur = urrecord
          ulrinv = Fmat_cuta
#IFDEF TIMING
          call cpu_time_now(endtime11)
          timecalculation(20)=timecalculation(20)+endtime11-starttime11
#ENDIF
      end if
  end do
  main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )
  call deallocate_gfunc(Fmat_cuta)
  call deallocate_gfunc(Fmat)
  call deallocate_gfunc(ulrinvul)
  call deallocate_gfunc(urrecord)
  call deallocate_gfunc(urFcuta_up)
  call deallocate_gfunc(uFcutal_up)
  deallocate( pvec_up )
  deallocate(cvec_up)
  deallocate(bvec_up)
  deallocate(v1)
  deallocate(v4)
  deallocate( gammainv_up )
  deallocate( gmattmp1_up )
  deallocate( gmattmp2_up )
  deallocate( Gmat1_up )
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(17)=timecalculation(17)+endtime-starttime
#ENDIF
end subroutine dqmc_proj_update
