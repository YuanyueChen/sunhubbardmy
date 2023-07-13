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
      vuvkmat%orb1(1,:) = ur%orb1(i1,:)
      vuvkmat%orb1(2,:) = ur%orb1(i2,:)
      ukmat%orb1(:,1) = ul%orb1(:,i1)
      ukmat%orb1(:,2) = ul%orb1(:,i2)
      Dmat%orb1 = czero
      Dmat%orb1(1,1) = this%delta_bmat_p%orb1(iflip, is)
      Dmat%orb1(2,2) = this%delta_bmat_m%orb1(iflip, is)
      call zgemm('n','n',2,ne,2,cone,Dmat%orb1,2,vuvkmat%orb1,2,czero,vkmat_tmp%orb1,2)  ! vkmat_tmp = D*vuvkmat
      call zgemm('n','n',2,ne,ne,cone,vkmat_tmp%orb1,2,ulrinv%orb1,ne,czero,vkmat%orb1,2)  ! vkmat = vkmat_tmp*ulrinv

      vukmat%orb1 = czero
      vukmat%orb1(1,1) = cone 
      vukmat%orb1(2,2) = cone
      call zgemm('n','n',2,2,ne,cone,vkmat%orb1,2,ukmat%orb1,ne,cone,vukmat%orb1,2)  ! vukmat = I + vkmat*ukmat
      
      call s_inv_det_qr_z(2, vukmat%orb1, ratio1) ! cal det(vukmat) and vukmat^-1

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
         ur%orb1(i1, :) = ur%orb1(i1, :) + vkmat_tmp%orb1(1, :)
         ur%orb1(i2, :) = ur%orb1(i2, :) + vkmat_tmp%orb1(2, :)

         ! update urlinv
         call zgemm('n','n',2,ne,2,cone,vukmat%orb1,2,vkmat%orb1,2,czero,vuvkmat%orb1,2)  ! vuvkmat = vukmat*vkmat
         call zgemm('n','n',ne,2,ne,cone,ulrinv%orb1,ne,ukmat%orb1,ne,czero,lrukmat%orb1,ne)  ! lrukmat = ulrinv*ukmat
         call zgemm('n','n',ne,ne,2,-cone,lrukmat%orb1,ne,vuvkmat%orb1,2,cone,ulrinv%orb1,ne)  ! ulrinv = - lrukmat*vuvkmat


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
