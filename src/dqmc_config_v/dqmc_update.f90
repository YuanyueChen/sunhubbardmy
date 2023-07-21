subroutine dqmc_update(this, gmat, ntau, nf)

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
  integer :: iflip, is, isp, i, i1, i2
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(gfunc) :: ukmat, vkmat, svkmat, vukmat, smat, Dmat
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  call allocate_gfunc( ukmat, ndim, 2)
  call allocate_gfunc( vkmat, 2, ndim)
  call allocate_gfunc( svkmat, 2, ndim)
  call allocate_gfunc( vukmat, 2, 2)
  call allocate_gfunc( smat, 2, 2)
  call allocate_gfunc( Dmat, 2, 2)

  accm  = 0.d0
  do i = 1, latt%nn_lf
      i1 = latt%nnlf_list(i)
      i2 = latt%nnlist(i1,nf)
      is = this%conf(i, nf, ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      ukmat%blk1 = czero
      vkmat%blk1 = czero
      vukmat%blk1 = czero
      Dmat%blk1 = czero
      Dmat%blk1(1,1) = this%delta_bmat_p%blk1(iflip, is)
      Dmat%blk1(2,2) = this%delta_bmat_m%blk1(iflip, is)
      ukmat%blk1(i1, 1) = this%delta_bmat_p%blk1(iflip, is)
      ukmat%blk1(i2, 2) = this%delta_bmat_m%blk1(iflip, is)
      vkmat%blk1(1, :) =      - gmat%blk1(i1, :)
      vkmat%blk1(1, i1) = cone - gmat%blk1(i1, i1)
      vkmat%blk1(2, :) =      - gmat%blk1(i2, :)
      vkmat%blk1(2, i2) = cone - gmat%blk1(i2, i2)
      vukmat%blk1(1,1) = cone 
      vukmat%blk1(2,2) = cone
      call zgemm('n','n',2,2,ndim,cone,vkmat%blk1,2,ukmat%blk1,ndim,cone,vukmat%blk1,2)
      call s_inv_det_qr_z(2, vukmat%blk1, ratio1) ! cal det(I+vukmat) and (I+vukmat)^-1
      call zgemm('n','n',2,2,2,cone,Dmat%blk1, 2, vukmat%blk1, 2, czero, smat%blk1, 2)

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

         ! update gmat
         ukmat%blk1(:,1) = gmat%blk1(:,i1)
         ukmat%blk1(:,2) = gmat%blk1(:,i2)

         call zgemm('n','n',2,ndim,2,-cone,smat%blk1,2,vkmat%blk1,2,czero,svkmat%blk1,2)  ! svkmat = -smat*vkmat
         call zgemm('n','n',ndim,ndim,2,cone,ukmat%blk1,ndim,svkmat%blk1,2,cone,gmat%blk1,ndim)  ! gmat + ukmat*svkmat

         ! flip
         this%conf(i,nf,ntau) =  isp
      endif
   end do
   main_obs(4) = main_obs(4) + dcmplx( accm, dble(latt%nn_lf) )
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(12)=timecalculation(12)+endtime-starttime
#ENDIF
end subroutine dqmc_update
