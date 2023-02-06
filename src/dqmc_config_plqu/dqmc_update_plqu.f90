subroutine dqmc_update_plqu(this, gmat, isite, ntau )

  use spring
  use model_para
  use dqmc_basic_data

  implicit none

  !arguments
  class(plqconf), intent(inout) :: this
  integer,intent(in) :: ntau, isite
  type(gfunc) :: gmat

  !local
  complex(dp) ::  ratio1, ratiotot
  integer :: iflip, is, i0, k, isp
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(gfunc) :: ukmat, vkmat, svkmat, vukmat, smat

  call allocate_gfunc( ukmat, ndim, latt%z_plq)
  call allocate_gfunc( vkmat, latt%z_plq, ndim)
  call allocate_gfunc( svkmat, latt%z_plq, ndim)
  call allocate_gfunc( vukmat, latt%z_plq, latt%z_plq)
  call allocate_gfunc( smat, latt%z_plq, latt%z_plq)

  accm  = 0.d0
      is = this%conf_plqu(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      ! orb 1
      ukmat%orb1 = czero
      vkmat%orb1 = czero
      do i0 = 1, latt%z_plq
          k = latt%plq_cord(i0,isite)
          ukmat%orb1(k, :) = this%delta_bmat_plqu_orb1(i0, :, iflip, is)
          vkmat%orb1(i0, :) =      - gmat%orb1(k, :)
          vkmat%orb1(i0, k) = cone - gmat%orb1(k, k)
      end do
      vukmat%orb1 = Imat_plq
      call zgemm('n','n',latt%z_plq,latt%z_plq,ndim,cone,vkmat%orb1,latt%z_plq,ukmat%orb1,ndim,cone,vukmat%orb1,latt%z_plq)
      call s_inv_det_qr_z(latt%z_plq, vukmat%orb1, ratio1) ! cal det(I+vukmat) and (I+vukmat)^-1
      call zgemm('n','n',latt%z_plq,latt%z_plq,latt%z_plq,cone,this%delta_bmat_plqu_orb1(:,:,iflip,is), latt%z_plq, vukmat%orb1, latt%z_plq, czero, smat%orb1, latt%z_plq)

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

         ! update gmat
         do i0 = 1, latt%z_plq
             ukmat%orb1(:,i0) = gmat%orb1(:,latt%plq_cord(i0,isite))
         end do
         call zgemm('n','n',latt%z_plq,ndim,latt%z_plq,-cone,smat%orb1,latt%z_plq,vkmat%orb1,latt%z_plq,czero,svkmat%orb1,latt%z_plq)  ! svkmat = -smat*vkmat
         call zgemm('n','n',ndim,ndim,latt%z_plq,cone,ukmat%orb1,ndim,svkmat%orb1,latt%z_plq,cone,gmat%orb1,ndim)  ! gmat + ukmat*svkmat

         ! flip
         this%conf_plqu(isite,ntau) =  isp
      endif
    main_obs(2) = main_obs(2) + dcmplx( accm, 1.d0 )
  end subroutine dqmc_update_plqu
