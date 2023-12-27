subroutine dqmc_update(this, gmat, ntau, nf )

  use spring
  use model_para
  use dqmc_basic_data
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  implicit none
  !arguments
  class(vconf), intent(inout) :: this
  integer,intent(in) :: ntau,nf
  type(gfunc) :: gmat

  !local
  complex(dp) ::  ratio1, ratiotot
  integer :: iflip, is, isp, isite,i1,i2
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zfunc) :: sscl
  complex(dp), allocatable, dimension(:,:) ::   gammainv_up,  Gcuta_up, gammainvtmp_up, gmattmp1_up, gmattmp2_up
  complex(dp), allocatable, dimension(:,:) ::   gammainv_dn,  Gcuta_dn, gammainvtmp_dn, gmattmp1_dn, gmattmp2_dn
  complex(dp), allocatable, dimension(:,:) ::   Gmat1_up, Gmat2_up,v,v2,v3
  complex(dp), allocatable, dimension(:,:) ::   Gmat1_dn, Gmat2_dn,v1, v4,v5, cvec_up, cvec_dn, bvec_up, bvec_dn
  complex(dp), allocatable, dimension(:) ::  Qmat_up, Qmat_dn
  complex(dp), allocatable, dimension(:) ::  Gtmp_up, Gtmp_dn
  integer, allocatable, dimension(:) :: pvec_up, pvec_dn
  integer :: i, ik, m,is1,is2
#IFDEF TIMING
  real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  allocate( pvec_up(nublock) )
  allocate( gammainv_up(nublock,nublock) )
  allocate( Gtmp_up(nublock) )
  allocate( Gcuta_up(ndim,ndim) )
  allocate(cvec_up(2,nublock))
  allocate(bvec_up(nublock,2))
  allocate(v1(2,nublock))
  allocate(v4(nublock,2))
  allocate(v(2,2))
  allocate(v2(2,2))
  allocate(v3(2,2))
  allocate(v5(2,2))
  allocate( gmattmp1_up(ndim,nublock) )
  allocate( gmattmp2_up(nublock,ndim) )
  allocate( Gmat1_up(ndim,nublock) )
    ! intial pvec,svec,wvec
  accm  = 0.d0
  ik = 0
  pvec_up = 0
  cvec_up = czero
  bvec_up = czero
  v2 = czero
  v3 = czero
  v5 = czero
  gammainv_up = czero
  Gcuta_up = gmat%blk1
  gmattmp1_up = czero
  gmattmp2_up = czero
  do isite = 1, latt%nn_lf
      i1 = latt%nnlf_list(isite)
      i2 = latt%nnlist(i1,nf)
      is = this%conf(isite, nf, ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      ! intial svec, wvec
      cvec_up = czero
      bvec_up = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( is1, is2, i )
!$OMP DO
      do i = 0, ik-2,2
        is1 = pvec_up(i+1)
        is2 = pvec_up(i+2)
        cvec_up(1,1+i) = gmat%blk1(i1,is1)
        cvec_up(1,2+i) = gmat%blk1(i1,is2)
        cvec_up(2,1+i) = gmat%blk1(i2,is1)
        cvec_up(2,2+i) = gmat%blk1(i2,is2)
        bvec_up(1+i,1) = gmat%blk1(is1,i1)
        bvec_up(2+i,1) = gmat%blk1(is2,i1)
        bvec_up(1+i,2) = gmat%blk1(is1,i2)
        bvec_up(2+i,2) = gmat%blk1(is2,i2)
      end do
!$OMP END DO
!$OMP END PARALLEL
      v1 = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( m, i )
!$OMP DO
      do i = 0, ik-2,2
        do m = 0, ik-2,2
          v1(1,i+1) = v1(1,i+1) + cvec_up(1,m+1)*gammainv_up(m+1,i+1)+cvec_up(1,m+2)*gammainv_up(m+2,i+1)
          v1(1,i+2) = v1(1,i+2) + cvec_up(1,m+1)*gammainv_up(m+1,i+2)+cvec_up(1,m+2)*gammainv_up(m+2,i+2)
          v1(2,i+1) = v1(2,i+1) + cvec_up(2,m+1)*gammainv_up(m+1,i+1)+cvec_up(2,m+2)*gammainv_up(m+2,i+1)
          v1(2,i+2) = v1(2,i+2) + cvec_up(2,m+1)*gammainv_up(m+1,i+2)+cvec_up(2,m+2)*gammainv_up(m+2,i+2)
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      v = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( i )
!$OMP DO
      do i = 0, ik-2,2
        v(1,1) = v(1,1) + v1(1,i+1)*bvec_up(i+1,1)+v1(1,i+2)*bvec_up(i+2,1)
        v(1,2) = v(1,2) + v1(1,i+1)*bvec_up(i+1,2)+v1(1,i+2)*bvec_up(i+2,2)
        v(2,1) = v(2,1) + v1(2,i+1)*bvec_up(i+1,1)+v1(2,i+2)*bvec_up(i+2,1)
        v(2,2) = v(2,2) + v1(2,i+1)*bvec_up(i+1,2)+v1(2,i+2)*bvec_up(i+2,2)
      end do
!$OMP END DO
!$OMP END PARALLEL
      v2(1,1) = this%delta_bmat_p%blk1(iflip,is)*(-v(1,1) - gmat%blk1(i1,i1) + cone) + cone
      v2(1,2) = this%delta_bmat_m%blk1(iflip,is)*(-v(1,2) - gmat%blk1(i1,i2))
      v2(2,1) = this%delta_bmat_p%blk1(iflip,is)*(-v(2,1) - gmat%blk1(i2,i1))
      v2(2,2) = this%delta_bmat_m%blk1(iflip,is)*(-v(2,2) - gmat%blk1(i2,i2) + cone) + cone

      ratio1 = v2(1,1)*v2(2,2) - v2(1,2)*v2(2,1)
      ratiotot = ratio1
      ratio_re = dble( ratiotot * phase )/dble( phase )

#IFDEF TEST
      write(fout, '(a,2e16.8)') ' in update_v, ratiotot = ', ratiotot
#ENDIF
      
      ratio_re_abs = ratio_re
      if ( ratio_re .lt. 0.d0 ) ratio_re_abs = - ratio_re

      random = spring_sfmt_stream()
      if ( ratio_re_abs.gt.random ) then

         accm  = accm + 1.d0
         !weight_track = weight_track + log( ratio_re_abs )
	     weight = dsqrt(dble(ratiotot*dconjg(ratiotot)))
	     phase =  phase*ratiotot/dcmplx(weight,0.d0)
#IFDEF TEST
         write(fout, '(a,2e16.8)') ' in update_v, phase = ', phase
#ENDIF

         ik = ik + 2
         ! store pvec, Qmat, gammainv
         pvec_up(ik-1) = i1
         pvec_up(ik) = i2
         v5(1,1)=cone + cone/ this%delta_bmat_p%blk1(iflip,is) - gmat%blk1(i1,i1)-v(1,1)
         v5(1,2)=-gmat%blk1(i1,i2)-v(1,2)
         v5(2,1)=-gmat%blk1(i2,i1)-v(2,1)
         v5(2,2)=cone + cone/ this%delta_bmat_m%blk1(iflip,is) - gmat%blk1(i2,i2)-v(2,2)
         v3(1,1)= cone/v5(1,1)+(cone/v5(1,1))*v5(1,2)*(cone/(v5(2,2)-(cone/v5(1,1))*v5(1,2)*v5(2,1)))*(cone/v5(1,1))*v5(2,1)
         v3(1,2)= -(cone/v5(1,1))*v5(1,2)*(cone/(v5(2,2)-(cone/v5(1,1))*v5(1,2)*v5(2,1)))
         v3(2,1)= -(cone/v5(1,1))*v5(2,1)*(cone/(v5(1,1)-(cone/v5(2,2))*v5(2,1)*v5(1,2))) 
         v3(2,2)= cone/(v5(2,2)-(cone/v5(1,1))*v5(1,2)*v5(2,1))

         gammainv_up(ik-1,ik-1) = v3(1,1)
         gammainv_up(ik-1,ik) = v3(1,2)
         gammainv_up(ik,ik-1) = v3(2,1)
         gammainv_up(ik,ik) = v3(2,2)

         v4 = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( m, i )
!$OMP DO
          do i = 0, ik-4,2
            do m = 0, ik-4,2
              v4(i+1,1) = v4(i+1,1) + gammainv_up(i+1,m+1)*bvec_up(m+1,1)+gammainv_up(i+1,m+2)*bvec_up(m+2,1) 
              v4(i+1,2) = v4(i+1,2) + gammainv_up(i+1,m+1)*bvec_up(m+1,2)+gammainv_up(i+1,m+2)*bvec_up(m+2,2)
              v4(i+2,1) = v4(i+2,1) + gammainv_up(i+2,m+1)*bvec_up(m+1,1)+gammainv_up(i+2,m+2)*bvec_up(m+2,1)
              v4(i+2,2) = v4(i+2,2) + gammainv_up(i+2,m+1)*bvec_up(m+1,2)+gammainv_up(i+2,m+2)*bvec_up(m+2,2)
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE ( i )
!$OMP DO
          do i = 0, ik-4,2
            gammainv_up(ik-1,i+1) = gammainv_up(ik-1,i+1) + v3(1,1)*v1(1,i+1)+v3(1,2)*v1(2,i+1)
            gammainv_up(ik-1,i+2) = gammainv_up(ik-1,i+2) + v3(1,1)*v1(1,i+2)+v3(1,2)*v1(2,i+2)
            gammainv_up(ik,i+1) = gammainv_up(ik,i+1) + v3(2,1)*v1(1,i+1)+v3(2,2)*v1(2,i+1)
            gammainv_up(ik,i+2) = gammainv_up(ik,i+2) + v3(2,1)*v1(1,i+2)+v3(2,2)*v1(2,i+2)
            gammainv_up(i+1,ik-1) = gammainv_up(i+1,ik-1) + v4(i+1,1)*v3(1,1)+v4(i+1,2)*v3(2,1)
            gammainv_up(i+2,ik-1) = gammainv_up(i+2,ik-1) +v4(i+2,1)*v3(1,1)+v4(i+2,2)*v3(2,1) 
            gammainv_up(i+1,ik) = gammainv_up(i+1,ik) + v4(i+1,1)*v3(1,2)+v4(i+1,2)*v3(2,2)
            gammainv_up(i+2,ik) = gammainv_up(i+2,ik) + v4(i+2,1)*v3(1,2)+v4(i+2,2)*v3(2,2) 
          end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE ( m, i )
!$OMP DO
          do i = 1, ik-2
            do m = 1, ik-2
              gammainv_up(i,m) = gammainv_up(i,m) + v4(i,1)*gammainv_up(ik-1,m) + v4(i,2)*gammainv_up(ik,m)
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
          Gcuta_up(:,i1) = Gcuta_up(:,i1) - this%delta_bmat_p%blk1(iflip, is)/(cone + this%delta_bmat_p%blk1(iflip, is))*Gcuta_up(:,i1)
          Gcuta_up(:,i2) = Gcuta_up(:,i2) - this%delta_bmat_m%blk1(iflip, is)/(cone + this%delta_bmat_m%blk1(iflip, is))*Gcuta_up(:,i2)
          gmattmp1_up(:,ik-1) = gmat%blk1(:,i1)
          gmattmp1_up(:,ik) = gmat%blk1(:,i2)
#IFDEF TEST
          write(fout, '(a,2e16.8)') ' after accepted, cvec_up = '
          do i = 1, 2
            write(fout, '(18(2e12.4))') cvec_up(i,:)
          end do
          write(fout, '(a,2e16.8)') ' after accepted, bvec_up = '
          do i = 1, ik
            write(fout, '(18(2e12.4))') bvec_up(i,:)
          end do
         write(fout, '(a,2e16.8)') ' after accepted, v3 = '
         do i = 1, 2
           write(fout, '(18(2e12.4))') v3(i,:)
         end do
         write(fout, '(a,2e16.8)') ' after accepted, v4 = '
         do i = 1, ik
           write(fout, '(18(2e12.4))') v4(i,:)
         end do
          write(fout, '(a,2e16.8)') ' after accepted, gammainv_up = '
         do i = 1, ik
           write(fout, '(18(2e12.4))') gammainv_up(i,:)
         end do 
#ENDIF
          this%conf(isite,nf,ntau) =  isp
     end if


      if( (ik.eq.nublock) .or. (isite.eq.latt%nn_lf) ) then
          ! submatrix update: update the whole Green function
#IFDEF TIMING
          call cpu_time_now(starttime11)
#ENDIF 
          do i = 1, ik
              is1 = pvec_up(i) 
              gmattmp2_up(i,:) = Gcuta_up(is1,:)
          end do
          call zgemm('N','N', ndim, nublock, ik, cone, gmattmp1_up, ndim, gammainv_up, nublock, czero, Gmat1_up, ndim)
          call zgemm('N','N', ndim, ndim, ik, cone, Gmat1_up, ndim, gmattmp2_up, nublock, cone, Gcuta_up, ndim)
          gmat%blk1 = Gcuta_up
          ik = 0
          gammainv_up = czero
#IFDEF TIMING
          call cpu_time_now(endtime11)
          timecalculation(18)=timecalculation(18)+endtime11-starttime11
#ENDIF
      end if
   end do
   main_obs(4) = main_obs(4) + dcmplx( accm, latt%nn_lf )
  deallocate( gmattmp1_up )
  deallocate( gmattmp2_up )
  deallocate( Gmat1_up )
  deallocate( pvec_up )
  deallocate( Gcuta_up)
  deallocate(cvec_up)
  deallocate(bvec_up)
  deallocate(v)
  deallocate(v1)
  deallocate(v2)
  deallocate(v3)
  deallocate(v4)
  deallocate(v5)
  deallocate(gammainv_up)
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(16)=timecalculation(16)+endtime-starttime
#ENDIF
end subroutine dqmc_update
