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
  integer,intent(in) :: ntau, nf
  type(gfunc) :: ul, ur, ulrinv

  !local
  complex(dp) ::  ratio1, ratiotot
  integer :: iflip, is, isp, i, i1, i2
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(gfunc) :: ukmat, vkmat, vuvkmat, vukmat, lrukmat, vkmat_tmp, Dmat
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  call allocate_gfunc( ukmat, ne, 2)
  call allocate_gfunc( vkmat, 2, ne)
  call allocate_gfunc( vuvkmat, 2, ne)
  call allocate_gfunc( vukmat, 2, 2)
  call allocate_gfunc( lrukmat, ne, 2)
  call allocate_gfunc( vkmat_tmp, 2, ne)
  call allocate_gfunc( Dmat, 2, 2)

  accm  = 0.d0
  do i = 1, latt%nn_lf
      i1 = latt%nnlf_list(i)
      i2 = latt%nnlist(i1,nf)
      is = this%conf(i, nf, ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1
      vuvkmat%blk1(1,:) = ur%blk1(i1,:)
      vuvkmat%blk1(2,:) = ur%blk1(i2,:)
      ukmat%blk1(:,1) = ul%blk1(:,i1)
      ukmat%blk1(:,2) = ul%blk1(:,i2)
      Dmat%blk1 = czero
      Dmat%blk1(1,1) = this%delta_bmat_p%blk1(iflip, is)
      Dmat%blk1(2,2) = this%delta_bmat_m%blk1(iflip, is)
      call zgemm('n','n',2,ne,2,cone,Dmat%blk1,2,vuvkmat%blk1,2,czero,vkmat_tmp%blk1,2)  ! vkmat_tmp = D*vuvkmat
      call zgemm('n','n',2,ne,ne,cone,vkmat_tmp%blk1,2,ulrinv%blk1,ne,czero,vkmat%blk1,2)  ! vkmat = vkmat_tmp*ulrinv

      vukmat%blk1 = czero
      vukmat%blk1(1,1) = cone 
      vukmat%blk1(2,2) = cone
      call zgemm('n','n',2,2,ne,cone,vkmat%blk1,2,ukmat%blk1,ne,cone,vukmat%blk1,2)  ! vukmat = I + vkmat*ukmat
      
      call s_inv_det_qr_z(2, vukmat%blk1, ratio1) ! cal det(vukmat) and vukmat^-1

      ratiotot = ratio1

#IFDEF TEST
      write(fout, '(a,2e16.8)') ' in update_v, ratiotot = ', ratiotot
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
         write(fout, '(a,2e16.8)') ' in update_v, phase = ', phase
#ENDIF 

         ! update ur
         ur%blk1(i1, :) = ur%blk1(i1, :) + vkmat_tmp%blk1(1, :)
         ur%blk1(i2, :) = ur%blk1(i2, :) + vkmat_tmp%blk1(2, :)

         ! update urlinv
         call zgemm('n','n',2,ne,2,cone,vukmat%blk1,2,vkmat%blk1,2,czero,vuvkmat%blk1,2)  ! vuvkmat = vukmat*vkmat
         call zgemm('n','n',ne,2,ne,cone,ulrinv%blk1,ne,ukmat%blk1,ne,czero,lrukmat%blk1,ne)  ! lrukmat = ulrinv*ukmat
         call zgemm('n','n',ne,ne,2,-cone,lrukmat%blk1,ne,vuvkmat%blk1,2,cone,ulrinv%blk1,ne)  ! ulrinv = - lrukmat*vuvkmat


         ! flip
         this%conf(i,nf,ntau) =  isp
     end if
  end do
  main_obs(4) = main_obs(4) + dcmplx( accm, dble(latt%nn_lf) )

#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(14)=timecalculation(14)+endtime-starttime
#ENDIF
end subroutine dqmc_proj_update
