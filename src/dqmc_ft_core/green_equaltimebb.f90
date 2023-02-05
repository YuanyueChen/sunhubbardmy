    subroutine green_equaltimebb( n, ndm, ure, dre, vre, gtt, logdetqr, logweightf, infoe )
      ! calcultate G(beta,beta), can save 3 matrix products when compare with green_equaltime
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(in) :: ure, vre
#IFDEF TEST_sig
      real(dp), allocatable, dimension(:) :: sigval_tmp
#ENDIF
      real(dp), dimension(ndm), intent(in) :: dre
      complex(dp), dimension(ndm,ndm), intent(out) :: gtt
      complex(dp), intent(in) :: logdetqr
      complex(dp), intent(out) :: logweightf
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      real(dp), allocatable, dimension(:) :: drmax, drmin
      complex(dp), allocatable, dimension(:,:) :: dvvdtmp, Btmp
      complex(dp) :: zlogdet_tmp
      real(dp) :: rlogdet_tmp

      allocate( drmax(ndm), drmin(ndm) )
#IFDEF TEST_sig
      allocate( sigval_tmp(ndm) )
#ENDIF
      allocate( dvvdtmp(ndm,ndm) )
      allocate( Btmp(ndm,ndm) )

      ! breakup dre = drmax * drmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dre,drmax,drmin)

#IFDEF TEST
      write(fout,*)
      write(fout,'(a,i6)') ' in green_equaltimebb, n = ', n
      write(fout,'(a)') '     dre(i)        drmax(i)        drmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dre(i), drmax(i), drmin(i)
      end do
#ENDIF

      !! >> g(beta,beta)
      ! drmax^-1 * ure^-1
      ! drmin * vre
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ! note uutmp^-1 = uutmp ^ +
              dvvdtmp(i,j) = dconjg(ure(j,i)) / drmax(i) + vre(i,j) * drmin(i)
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

      !call s_invlu_z(ndm,dvvdtmp)
      call s_inv_logdet_lu_z(ndm,dvvdtmp,zlogdet_tmp)
      call s_v_invd_uH( ndm, dvvdtmp, drmax, ure, gtt )

      rlogdet_tmp = 0.d0
      do i = 1, ndm
          rlogdet_tmp = rlogdet_tmp + log( drmax(i) )
      end do
      logweightf =  logdetqr + zlogdet_tmp + dcmplx(rlogdet_tmp,0.d0)

      infoe = 0

      deallocate( drmin )
      deallocate( drmax )
      deallocate( dvvdtmp )
      deallocate( Btmp )
    end subroutine green_equaltimebb
