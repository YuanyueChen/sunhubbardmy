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
  complex(dp) ::  ratio1, ratiotot, v
  integer :: iflip, is, isp, isite
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  type(zfunc) :: sscl
  complex(dp), allocatable, dimension(:,:) ::   gammainv_up, gammainv_uptmp, Gtmp_up, Gcuta_up
  complex(dp), allocatable, dimension(:,:) ::   gammainv_dn, gammainv_dntmp, Gtmp_dn, Gcuta_dn
  complex(dp), allocatable, dimension(:,:) ::   Amat_up
  complex(dp), allocatable, dimension(:,:) ::   Amat_dn
  complex(dp), allocatable, dimension(:) ::  Qmat_up, Qmat_dn,v1,svec_up,svec_dn,wvec_up,wvec_dn
  integer, allocatable, dimension(:,:) :: pmat_up, pmat_dn
  integer, allocatable, dimension(:) :: pvec_up, pvec_dn
  integer :: i, ik, m
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF
  allocate( pmat_up(nublock,ndim) )
  allocate( pvec_up(nublock) )
  allocate( gammainv_up(nublock,nublock) )
  allocate( Qmat_up(nublock) )
  allocate( Gtmp_up(ndim,ndim) )
  allocate( Gcuta_up(ndim,ndim) )


  accm  = 0.d0
  ik = 0
  ! intial pvec
  pvec_up = 0
  do isite = 1, latt%nsites
      ! submatrix update: after nublock steps of local update, perform a whole update of Green function
      ! calculate weight ratio, fermion part
      is = this%conf_u(isite,ntau)
      iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
      isp = mod(is+iflip-1,this%lcomp) + 1

      if( ik.eq.0 ) then
          ratio1 = (cone - gmat%blk1(isite,isite))*this%delta_bmat%blk1(iflip, is) + cone
      else
        allocate(svec_up(ik))
        allocate(wvec_up(ik))
        ! intial svec, wvec
        svec_up = czero
        wvec_up = czero
        do i = 1, ik
            svec_up(i) = gmat%blk1(isite,pvec_up(i))
            wvec_up(i) = gmat%blk1(pvec_up(i),isite)
        end do
        allocate(gammainv_uptmp(ik,ik))
        gammainv_uptmp = czero
        do i = 1, ik
            do m = 1, ik
                gammainv_uptmp(i,m) = gammainv_up(i,m)
            end do
        end do
        allocate(v1(ik))
        v1 = czero
        do i = 1, ik
            do m = 1, ik
                v1(i) = v1(i) + svec_up(m)*gammainv_uptmp(m,i)
            end do
        end do
        v = czero
        do i = 1, ik
            v = v + v1(i)*wvec_up(i)
        end do
        v = v + gmat%blk1(isite,isite) - cone 
        ratio1 = -v *this%delta_bmat%blk1(iflip, is) + cone
        deallocate(svec_up)
        deallocate(wvec_up)
        deallocate(gammainv_uptmp)
        deallocate(v1)
      end if
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
         pmat_up(ik,isite) = 1
         Qmat_up(ik) = this%delta_bmat%blk1(iflip, is)
         if( ik.eq.1 ) then
             gammainv_up(1,1) = cone/(cone + cone/Qmat_up(1) - gmat%blk1(isite,isite))
         else
             allocate(Amat_up(ik,ik))
             Amat_up = czero
             do i = 1, ik
                do m = 1, ik
                   Amat_up(i,m) = -gmat%blk1(pvec_up(i),pvec_up(m))
                end do
             end do
             do i = 1, ik
                Amat_up(i,i) = Amat_up(i,i)+ cone+ cone/Qmat_up(i)
             end do
             call s_invlu_z(ik, Amat_up)
             do i = 1, ik
                 do m = 1, ik
                   gammainv_up(i,m) = Amat_up(i,m)
                 end do
             end do
              deallocate(Amat_up)
         end if
         
         ! flip
         this%conf_u(isite,ntau) =  isp
     end if

      if( (ik.eq.nublock) .or. (isite.eq.latt%nsites) ) then
          ! submatrix update: update the whole Green function
          allocate(gammainv_uptmp(ik,ik))
          gammainv_uptmp = czero
          do i = 1, ik
              do m = 1, ik
                  gammainv_uptmp(i,m) = gammainv_up(i,m)
              end do
          end do
          Gtmp_up = czero
          do i = 1, ik
              Gtmp_up(pvec_up(i),pvec_up(i)) = -Qmat_up(i)/(cone + Qmat_up(i))
          end do
          do i = 1, ndim
              Gtmp_up(i,i) = Gtmp_up(i,i) + cone
          end do
          Gcuta_up = czero
          call zgemm('N','N', ndim, ndim, ndim, cone, gmat%blk1, ndim, Gtmp_up, ndim, czero, Gcuta_up, ndim)   
          Gtmp_up = czero
          do i = 1, ik
              do m = 1, ik
                  Gtmp_up(pvec_up(i),pvec_up(m)) = gammainv_uptmp(i,m)
              end do
          end do
          call zgemm('N','N', ndim, ndim, ndim, cone, gmat%blk1, ndim, Gtmp_up, ndim, czero, Gtmp_up, ndim)
          do i = 1, ndim
              Gtmp_up(i,i) = Gtmp_up(i,i) + cone
          end do
          call zgemm('N','N', ndim, ndim, ndim, cone, Gtmp_up, ndim, Gcuta_up, ndim, czero, gmat%blk1, ndim)
          deallocate(gammainv_uptmp)
          ik = 0
      end if
   end do
   main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )
  deallocate( pmat_up )
  deallocate( pvec_up )
  deallocate( gammainv_up )
  deallocate( Qmat_up )
  deallocate( Gtmp_up )
  deallocate( Gcuta_up)

#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(16)=timecalculation(16)+endtime-starttime
#ENDIF
end subroutine dqmc_update
