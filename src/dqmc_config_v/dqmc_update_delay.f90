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
  integer :: iflip, is, isp, i, i1, i2, ik, m
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  complex(dp), allocatable, dimension(:,:) :: smat, gpart, vukmat
  complex(dp), allocatable, dimension(:,:) :: avec
  complex(dp), allocatable, dimension(:,:) :: bvec
  complex(dp), allocatable, dimension(:,:) :: grow

#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  allocate( avec(ndim,nublock) )
  allocate( bvec(ndim,nublock) )
  allocate( grow(ndim,2) )
  allocate( smat(2,2), gpart(2,2), vukmat(2,2) )

  accm  = 0.d0
  ik = 0

  do i = 1, latt%nn_lf
      i1 = latt%nnlf_list(i)
      i2 = latt%nnlist(i1,nf)
      is = this%conf(i, nf, ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      gpart(1,1) = gmat%orb1(i1,i1)
      gpart(1,2) = gmat%orb1(i1,i2)
      gpart(2,1) = gmat%orb1(i2,i1)
      gpart(2,2) = gmat%orb1(i2,i2)
      do m = 0, ik-2, 2
          gpart(1,1) = gpart(1,1) + avec(i1,1+m)*bvec(i1,1+m) + avec(i1,2+m)*bvec(i1,2+m)
          gpart(1,2) = gpart(1,2) + avec(i1,1+m)*bvec(i2,1+m) + avec(i1,2+m)*bvec(i2,2+m)
          gpart(2,1) = gpart(2,1) + avec(i2,1+m)*bvec(i1,1+m) + avec(i2,2+m)*bvec(i1,2+m)
          gpart(2,2) = gpart(2,2) + avec(i2,1+m)*bvec(i2,1+m) + avec(i2,2+m)*bvec(i2,2+m)
      end do
      vukmat(1,1) = (cone - gpart(1,1))*this%delta_bmat_p%orb1(iflip, is) + cone
      vukmat(1,2) =       - gpart(1,2) *this%delta_bmat_m%orb1(iflip, is)
      vukmat(2,1) =       - gpart(2,1) *this%delta_bmat_p%orb1(iflip, is)
      vukmat(2,2) = (cone - gpart(2,2))*this%delta_bmat_m%orb1(iflip, is) + cone
      call s_inv_det_qr_z(2, vukmat, ratio1) ! cal det(I+vukmat) and (I+vukmat)^-1
      smat(1,:) = this%delta_bmat_p%orb1(iflip, is) * vukmat(1,:)
      smat(2,:) = this%delta_bmat_m%orb1(iflip, is) * vukmat(2,:)

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

          grow(:,1) = gmat%orb1(:,i1)
          grow(:,2) = gmat%orb1(:,i2)
          bvec(:,1+ik) = gmat%orb1(i1,:)
          bvec(:,2+ik) = gmat%orb1(i2,:)
          bvec(i1,1+ik) = gmat%orb1(i1,i1) - cone
          bvec(i2,2+ik) = gmat%orb1(i2,i2) - cone
          if ( ratio_re_abs.gt.random ) then
              do m = 0, ik-2, 2
                  grow(:,1) = grow(:,1) + avec(:,1+m)*bvec(i1,1+m) + avec(:,2+m)*bvec(i1,2+m)
                  grow(:,2) = grow(:,2) + avec(:,1+m)*bvec(i2,1+m) + avec(:,2+m)*bvec(i2,2+m)
                  bvec(:,1+ik) = bvec(:,1+ik) + avec(i1,1+m)*bvec(:,1+m) + avec(i1,2+m)*bvec(:,2+m)
                  bvec(:,2+ik) = bvec(:,2+ik) + avec(i2,1+m)*bvec(:,1+m) + avec(i2,2+m)*bvec(:,2+m)
              end do
          end if
          avec(:,1+ik) = grow(:,1)*smat(1,1) + grow(:,2)*smat(2,1)
          avec(:,2+ik) = grow(:,1)*smat(1,2) + grow(:,2)*smat(2,2)
         ! flip
         this%conf(i,nf,ntau) =  isp
         ik = ik + 2
      end if

      if( (ik.eq.nublock) .or. (i.eq.latt%nn_lf) ) then
          ik = 0
          ! delay update: update the whole Green function
          call  zgemm('N', 'T', ndim, ndim, nublock, cone, avec, ndim, bvec, ndim, cone, gmat%orb1, ndim)
          if( i.lt.latt%nn_lf) then
              ! intial avec, bvec
              avec = czero
              bvec = czero
          end if
      end if
   end do
   main_obs(4) = main_obs(4) + dcmplx( accm, latt%nn_lf )
  deallocate( smat, gpart, vukmat )
  deallocate( grow )
  deallocate( bvec )
  deallocate( avec )
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(13)=timecalculation(13)+endtime-starttime
#ENDIF
end subroutine dqmc_update
