subroutine dqmc_proj_update_plqu(this, isite, ntau, ul, ur, ulrinv)

  use spring
  use model_para
  use dqmc_ft_basic_data

  implicit none

  !arguments
  class(plqconf), intent(inout) :: this
  integer,intent(in) :: ntau, isite
  type(gfunc) :: ul, ur, ulrinv

  !local
  complex(dp) ::  ratio1, ratiotot
  integer :: iflip, is, i0, k, isp
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(gfunc) :: ukmat, vkmat, vuvkmat, vukmat, lrukmat, vkmat_tmp

  call allocate_gfunc( ukmat, ne, latt%z_plq)
  call allocate_gfunc( vkmat, latt%z_plq, ne)
  call allocate_gfunc( vuvkmat, latt%z_plq, ne)
  call allocate_gfunc( vukmat, latt%z_plq, latt%z_plq)
  call allocate_gfunc( lrukmat, ne, latt%z_plq)
  call allocate_gfunc( vkmat_tmp, latt%z_plq, ne)

  accm  = 0.d0
      is = this%conf_plqu(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      do i0 = 1, latt%z_plq
          k = latt%plq_cord(i0,isite)
          vuvkmat%orb1(i0,:) = ur%orb1(k,:)
          ukmat%orb1(:,i0) = ul%orb1(:,k)
      end do
      call zgemm('n','n',latt%z_plq,ne,latt%z_plq,cone,this%delta_bmat_plqu_orb1(:,:,iflip,is),latt%z_plq,vuvkmat%orb1,latt%z_plq,czero,vkmat_tmp%orb1,latt%z_plq)  ! vkmat_tmp = D*vuvkmat
      call zgemm('n','n',latt%z_plq,ne,ne,cone,vkmat_tmp%orb1,latt%z_plq,ulrinv%orb1,ne,czero,vkmat%orb1,latt%z_plq)  ! vkmat = vkmat_tmp*ulrinv

      vukmat%orb1 = Imat_plq
      call zgemm('n','n',latt%z_plq,latt%z_plq,ne,cone,vkmat%orb1,latt%z_plq,ukmat%orb1,ne,cone,vukmat%orb1,latt%z_plq)  ! vukmat = I + vkmat*ukmat

      call s_inv_det_qr_z(latt%z_plq, vukmat%orb1, ratio1) ! cal det(vukmat) and vukmat^-1

      ratiotot = ratio1**nflr
      if( lprojplqu ) then
          ratiotot = ratiotot*this%phase_ratio(iflip, is)
      else
          ratiotot = ratiotot*dcmplx( gaml(isp)/gaml(is), 0.d0 ) * this%phase_ratio(iflip, is)
      end if
#IFDEF TEST
      write(fout, '(a,2e16.8)') ' in update_plqu, ratio1 = ', ratio1
      write(fout, '(a,2e16.8)') ' in update_plqu, phase_ratio = ', this%phase_ratio(iflip, is)
      write(fout, '(a,2e16.8)') ' in update_plqu, ratiotot = ', ratiotot
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
         write(fout, '(a,2e16.8)') ' in update_plqu, phase = ', phase
#ENDIF

         ! update ur
         do i0 = 1, latt%z_plq
             k = latt%plq_cord(i0,isite)
             ur%orb1(k, :) = ur%orb1(k, :) + vkmat_tmp%orb1(i0, :)
         end do

         ! update urlinv
         call zgemm('n','n',latt%z_plq,ne,latt%z_plq,cone,vukmat%orb1,latt%z_plq,vkmat%orb1,latt%z_plq,czero,vuvkmat%orb1,latt%z_plq)  ! vuvkmat = vukmat*vkmat
         call zgemm('n','n',ne,latt%z_plq,ne,cone,ulrinv%orb1,ne,ukmat%orb1,ne,czero,lrukmat%orb1,ne)  ! lrukmat = ulrinv*ukmat
         call zgemm('n','n',ne,ne,latt%z_plq,-cone,lrukmat%orb1,ne,vuvkmat%orb1,latt%z_plq,cone,ulrinv%orb1,ne)  ! ulrinv - lrukmat*vuvkmat

         ! flip
         this%conf_plqu(isite,ntau) =  isp
     end if
    main_obs(2) = main_obs(2) + dcmplx( accm, 1.d0 )
  end subroutine dqmc_proj_update_plqu
