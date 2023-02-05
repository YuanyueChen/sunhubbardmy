    subroutine green_equaltime( n, ndm, ure, dre, vre, vle, dle, ule, gtt, logdetqr, logdetql, logweightf, infoe )
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(in) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      complex(dp), dimension(ndm,ndm), intent(out) :: gtt
      complex(dp), intent(in) :: logdetqr, logdetql
      complex(dp), intent(out) :: logweightf
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin
#IFDEF TEST_sig
      real(dp), allocatable, dimension(:) :: sigval_tmp
#ENDIF
      complex(dp), allocatable, dimension(:,:) :: uutmp, vvtmp, dvvdtmp, Btmp
      complex(dp) :: zlogdet_tmp
      real(dp) :: rlogdet_tmp1, rlogdet_tmp2

      allocate( drmax(ndm), drmin(ndm), dlmax(ndm), dlmin(ndm) )
#IFDEF TEST_sig
      allocate( sigval_tmp(ndm) )
#ENDIF
      allocate( uutmp(ndm,ndm) )
      allocate( vvtmp(ndm,ndm) )
      allocate( dvvdtmp(ndm,ndm) )
      allocate( Btmp(ndm,ndm) )

      ! breakup dre = drmax * drmin
      !         dle = dlmax * dlmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dre,drmax,drmin)
      call s_dvec_min_max(ndm,dle,dlmax,dlmin)

#IFDEF TEST
      write(fout,*)
      write(fout,'(a,i6)') ' in green_equaltime, n = ', n
      write(fout,'(a)') '     dre(i)        drmax(i)        drmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dre(i), drmax(i), drmin(i)
      end do
      write(fout,*)
      write(fout,'(a)') '     dle(i)        dlmax(i)        dlmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dle(i), dlmax(i), dlmin(i)
      end do
#ENDIF

      ! uutmp = ure^-1*ule
      call zgemm('c','n',ndm,ndm,ndm,cone,ure,ndm,ule,ndm,czero,uutmp,ndm)
      ! vvtmp = vre*vle^H
      call zgemm('n','c',ndm,ndm,ndm,cone,vre,ndm,vle,ndm,czero,vvtmp,ndm)

      !! >> g(t,t)
      ! drmax^-1 * ( ule * ure )^-1 dlmax^-1
      ! drmin * ( vre * vle ) * dlmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ! note uutmp^-1 = uutmp ^ +
              dvvdtmp(i,j) = uutmp(i,j) / ( drmax(i)*dlmax(j) ) + vvtmp(i,j) * drmin(i) * dlmin(j)
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL

#IFDEF TEST_sig
      ! for test
      grtmp = dvvdtmp
      call s_svd_zg(ndim, ndim, ndim, grtmp, Atmp, sigval_tmp, Btmp)
      write(fout,'(a)') ' singular value of dvvdtmp^-1 = '
      do i = 1, ndm
          write(fout,'(e24.16)') 1.d0/sigval_tmp(i)
      end do
      deallocate(sigval_tmp)
#ENDIF

      !!!call s_invlu_z(ndm,dvvdtmp)
      call s_inv_logdet_lu_z(ndm,dvvdtmp,zlogdet_tmp)
      call s_v_invd_u( ndm, ule, dlmax, dvvdtmp, Btmp )
      call s_v_invd_uH( ndm, Btmp, drmax, ure, gtt )

      rlogdet_tmp1 = 0.d0
      do i = 1, ndm
          rlogdet_tmp1 = rlogdet_tmp1 + log( dlmax(i) )
      end do
      rlogdet_tmp2 = 0.d0
      do i = 1, ndm
          rlogdet_tmp2 = rlogdet_tmp2 + log( drmax(i) )
      end do
      logweightf = -logdetql + logdetqr + zlogdet_tmp + dcmplx(rlogdet_tmp1+rlogdet_tmp2,0.d0)


      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( uutmp )
      deallocate( vvtmp )
      deallocate( dvvdtmp )
      deallocate( Btmp )
    end subroutine green_equaltime
