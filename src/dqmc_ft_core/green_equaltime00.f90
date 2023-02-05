    subroutine green_equaltime00( n, ndm, vle, dle, ule, gtt, logdetql, logweightf, infoe )
      ! calcultate G(0,0), can save 3 matrix products when compare with green_equaltime
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(in) :: vle, ule
      real(dp), dimension(ndm), intent(in) :: dle
#IFDEF TEST_sig
      real(dp), allocatable, dimension(:) :: sigval_tmp
#ENDIF
      complex(dp), dimension(ndm,ndm), intent(out) :: gtt
      complex(dp), intent(in) :: logdetql
      complex(dp), intent(out) :: logweightf
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      real(dp), allocatable, dimension(:) :: dlmax, dlmin
      complex(dp), allocatable, dimension(:,:) :: dvvdtmp, Btmp
      complex(dp) :: zlogdet_tmp
      real(dp) :: rlogdet_tmp

      allocate( dlmax(ndm), dlmin(ndm) )
#IFDEF TEST_sig
      allocate( sigval_tmp(ndm) )
#ENDIF
      allocate( dvvdtmp(ndm,ndm) )
      allocate( Btmp(ndm,ndm) )

      ! breakup dle = dlmax * dlmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dle,dlmax,dlmin)

#IFDEF TEST
      write(fout,'(a,i6)') ' in green_equaltime00, n = ', n
      write(fout,'(a)') '     dle(i)        dlmax(i)        dlmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dle(i), dlmax(i), dlmin(i)
      end do
#ENDIF

      !! >> g(0,0)
      ! ule^-1 * dlmax^-1
      ! vle * dlmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ! note uutmp^-1 = uutmp ^ +
              dvvdtmp(i,j) = ule(i,j) / dlmax(j) + dconjg( vle(j,i) ) * dlmin(j)
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
      call s_v_invd_u( ndm, ule, dlmax, dvvdtmp, gtt )

      rlogdet_tmp = 0.d0
      do i = 1, ndm
          rlogdet_tmp = rlogdet_tmp + log( dlmax(i) )
      end do
      logweightf = -logdetql + zlogdet_tmp + dcmplx(rlogdet_tmp,0.d0)

      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( dvvdtmp )
      deallocate( Btmp )
    end subroutine green_equaltime00
