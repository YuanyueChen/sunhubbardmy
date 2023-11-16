subroutine dqmc_proj_update(this, ntau, nf, ul, ur, ulrinv)

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
  type(gfunc) :: ul, ur, ulrinv

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
  type(gfunc) :: Fmat, ulrinvul, urrecord
  integer :: i, ik, m,is1,is2
#IFDEF TIMING
  real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  allocate( pvec_up(nublock) )
  allocate( gammainv_up(nublock,nublock) )
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
  call allocate_gfunc(Fmat, ndim, ndim)
  call allocate_gfunc(ulrinvul, ne, ndim)
  call allocate_gfunc(urrecord, ndim, ne)

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
  gmattmp1_up = czero
  gmattmp2_up = czero
  urrecord%blk1 = ur%blk1
  call zgemm('N', 'N', ne, ndim, ne, cone, ulrinv%blk1, ne, ul%blk1, ne, czero, ulrinvul%blk1, ne)
  call zgemm('N', 'N', ndim, ndim, ne, cone, ur%blk1, ndim, ulrinvul%blk1, ne, czero, Fmat%blk1, ndim)
  do isite = 1, latt%nn_lf
      i1 = latt%nnlf_list(isite)
      i2 = latt%nnlist(i1,nf)
      is = this%conf(isite, nf, ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      ! intial svec, wvec
      cvec_up = czero
      bvec_up = czero
      do i = 0, ik-2,2
        is1 = pvec_up(i+1)
        is2 = pvec_up(i+2)
        cvec_up(1,1+i) = Fmat%blk1(i1,is1)
        cvec_up(1,2+i) = Fmat%blk1(i1,is2)
        cvec_up(2,1+i) = Fmat%blk1(i2,is1)
        cvec_up(2,2+i) = Fmat%blk1(i2,is2)
        bvec_up(1+i,1) = Fmat%blk1(is1,i1)
        bvec_up(2+i,1) = Fmat%blk1(is2,i1)
        bvec_up(1+i,2) = Fmat%blk1(is1,i2)
        bvec_up(2+i,2) = Fmat%blk1(is2,i2)
      end do
      v1 = czero
      do i = 0, ik-2,2
        do m = 0, ik-2,2
          v1(1,i+1) = v1(1,i+1) + cvec_up(1,m+1)*gammainv_up(m+1,i+1)+cvec_up(1,m+2)*gammainv_up(m+2,i+1)
          v1(1,i+2) = v1(1,i+2) + cvec_up(1,m+1)*gammainv_up(m+1,i+2)+cvec_up(1,m+2)*gammainv_up(m+2,i+2)
          v1(2,i+1) = v1(2,i+1) + cvec_up(2,m+1)*gammainv_up(m+1,i+1)+cvec_up(2,m+2)*gammainv_up(m+2,i+1)
          v1(2,i+2) = v1(2,i+2) + cvec_up(2,m+1)*gammainv_up(m+1,i+2)+cvec_up(2,m+2)*gammainv_up(m+2,i+2)
        end do
      end do
      v = czero
      do i = 0, ik-2,2
        v(1,1) = v(1,1) + v1(1,i+1)*bvec_up(i+1,1)+v1(1,i+2)*bvec_up(i+2,1)
        v(1,2) = v(1,2) + v1(1,i+1)*bvec_up(i+1,2)+v1(1,i+2)*bvec_up(i+2,2)
        v(2,1) = v(2,1) + v1(2,i+1)*bvec_up(i+1,1)+v1(2,i+2)*bvec_up(i+2,1)
        v(2,2) = v(2,2) + v1(2,i+1)*bvec_up(i+1,2)+v1(2,i+2)*bvec_up(i+2,2)
      end do
      v2(1,1) = this%delta_bmat_p%blk1(iflip,is)*(-v(1,1) + Fmat%blk1(i1,i1) ) + cone
      v2(1,2) = this%delta_bmat_m%blk1(iflip,is)*(-v(1,2) + Fmat%blk1(i1,i2))
      v2(2,1) = this%delta_bmat_p%blk1(iflip,is)*(-v(2,1) + Fmat%blk1(i2,i1))
      v2(2,2) = this%delta_bmat_m%blk1(iflip,is)*(-v(2,2) + Fmat%blk1(i2,i2) ) + cone

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
         v5(1,1)= cone/this%delta_bmat_p%blk1(iflip,is) + Fmat%blk1(i1,i1)-v(1,1)
         v5(1,2)= Fmat%blk1(i1,i2)-v(1,2)
         v5(2,1)= Fmat%blk1(i2,i1)-v(2,1)
         v5(2,2)= cone/this%delta_bmat_m%blk1(iflip,is) + Fmat%blk1(i2,i2)-v(2,2)
         v3(1,1)= cone/v5(1,1)+(cone/v5(1,1))*v5(1,2)*(cone/(v5(2,2)-(cone/v5(1,1))*v5(1,2)*v5(2,1)))*(cone/v5(1,1))*v5(2,1)
         v3(1,2)= -(cone/v5(1,1))*v5(1,2)*(cone/(v5(2,2)-(cone/v5(1,1))*v5(1,2)*v5(2,1)))
         v3(2,1)= -(cone/v5(1,1))*v5(2,1)*(cone/(v5(1,1)-(cone/v5(2,2))*v5(2,1)*v5(1,2))) 
         v3(2,2)= cone/(v5(2,2)-(cone/v5(1,1))*v5(1,2)*v5(2,1))

         gammainv_up(ik-1,ik-1) = v3(1,1)
         gammainv_up(ik-1,ik) = v3(1,2)
         gammainv_up(ik,ik-1) = v3(2,1)
         gammainv_up(ik,ik) = v3(2,2)

         v4 = czero
          do i = 0, ik-4,2
            do m = 0, ik-4,2
              v4(i+1,1) = v4(i+1,1) + gammainv_up(i+1,m+1)*bvec_up(m+1,1)+gammainv_up(i+1,m+2)*bvec_up(m+2,1) 
              v4(i+1,2) = v4(i+1,2) + gammainv_up(i+1,m+1)*bvec_up(m+1,2)+gammainv_up(i+1,m+2)*bvec_up(m+2,2)
              v4(i+2,1) = v4(i+2,1) + gammainv_up(i+2,m+1)*bvec_up(m+1,1)+gammainv_up(i+2,m+2)*bvec_up(m+2,1)
              v4(i+2,2) = v4(i+2,2) + gammainv_up(i+2,m+1)*bvec_up(m+1,2)+gammainv_up(i+2,m+2)*bvec_up(m+2,2)
            end do
          end do
          do i = 0, ik-4,2
            gammainv_up(ik-1,i+1) = gammainv_up(ik-1,i+1) - v3(1,1)*v1(1,i+1) - v3(1,2)*v1(2,i+1)
            gammainv_up(ik-1,i+2) = gammainv_up(ik-1,i+2) - v3(1,1)*v1(1,i+2) - v3(1,2)*v1(2,i+2)
            gammainv_up(ik,i+1) = gammainv_up(ik,i+1) - v3(2,1)*v1(1,i+1) - v3(2,2)*v1(2,i+1)
            gammainv_up(ik,i+2) = gammainv_up(ik,i+2) - v3(2,1)*v1(1,i+2) - v3(2,2)*v1(2,i+2)
            gammainv_up(i+1,ik-1) = gammainv_up(i+1,ik-1) - v4(i+1,1)*v3(1,1)- v4(i+1,2)*v3(2,1)
            gammainv_up(i+2,ik-1) = gammainv_up(i+2,ik-1) - v4(i+2,1)*v3(1,1)- v4(i+2,2)*v3(2,1) 
            gammainv_up(i+1,ik) = gammainv_up(i+1,ik) - v4(i+1,1)*v3(1,2) - v4(i+1,2)*v3(2,2)
            gammainv_up(i+2,ik) = gammainv_up(i+2,ik) - v4(i+2,1)*v3(1,2) - v4(i+2,2)*v3(2,2) 
          end do
          do i = 1, ik-2
            do m = 1, ik-2
              gammainv_up(i,m) = gammainv_up(i,m) - v4(i,1)*gammainv_up(ik-1,m) - v4(i,2)*gammainv_up(ik,m)
            end do
          end do
          urrecord%blk1(i1,:) = urrecord%blk1(i1,:) + this%delta_bmat_p%blk1(iflip, is)*urrecord%blk1(i1,:)
          urrecord%blk1(i2,:) = urrecord%blk1(i2,:) + this%delta_bmat_m%blk1(iflip, is)*urrecord%blk1(i2,:)
          gmattmp1_up(:,ik-1) = Fmat%blk1(:,i1)
          gmattmp1_up(:,ik) = Fmat%blk1(:,i2)
          gmattmp2_up(ik-1,:) = Fmat%blk1(i1,:)
          gmattmp2_up(ik,:) = Fmat%blk1(i2,:)
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
         write(fout, '(a,2e16.8)') ' after accepted, v1 = '
          do i = 1, 2
            write(fout, '(18(2e12.4))') v1(i,:)
          end do
         write(fout, '(a,2e16.8)') ' after accepted, v = '
         do i = 1, 2
           write(fout, '(18(2e12.4))') v(i,:)
         end do
          write(fout, '(a,2e16.8)') ' after accepted, gammainv_up = '
         do i = 1, ik
           write(fout, '(18(2e12.4))') gammainv_up(i,:)
         end do 
#ENDIF
          this%conf(isite,nf,ntau) =  isp
     end if


      if( (ik.eq.nublock) .and. (isite.lt.latt%nn_lf) ) then
          ! submatrix update: update the whole Green function
#IFDEF TIMING
          call cpu_time_now(starttime11)
#ENDIF 
          call zgemm('N','N', ndim, nublock, ik, cone, gmattmp1_up, ndim, gammainv_up, nublock, czero, Gmat1_up, ndim)
          call zgemm('N','N', ndim, ndim, ik, -cone, Gmat1_up, ndim, gmattmp2_up, nublock, cone, Fmat%blk1, ndim)
          ik = 0
          gammainv_up = czero
#IFDEF TIMING
          call cpu_time_now(endtime11)
          timecalculation(20)=timecalculation(20)+endtime11-starttime11
#ENDIF
      end if
     
      if(  isite.eq.latt%nn_lf ) then
#IFDEF TIMING
          call cpu_time_now(starttime11)
#ENDIF 
          ik = 0
          ! delay update: update the R and the (LR)^-1
          ! update the R
          ur%blk1 = urrecord%blk1
          call zgemm('N', 'N', ne, ne, ndim, cone, ul%blk1, ne, ur%blk1, ndim, czero, ulrinv%blk1, ne)
          call s_invlu_z(ne, ulrinv%blk1) 
#IFDEF TIMING
          call cpu_time_now(endtime11)
          timecalculation(20)=timecalculation(20)+endtime11-starttime11
#ENDIF
      end if
   end do
   main_obs(4) = main_obs(4) + dcmplx( accm, latt%nn_lf )
  call deallocate_gfunc(Fmat)
  call deallocate_gfunc(ulrinvul)
  call deallocate_gfunc(urrecord)
  deallocate( gmattmp1_up )
  deallocate( gmattmp2_up )
  deallocate( Gmat1_up )
  deallocate( pvec_up )
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
  timecalculation(17)=timecalculation(17)+endtime-starttime
#ENDIF
end subroutine dqmc_proj_update
