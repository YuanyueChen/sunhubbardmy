    subroutine green_tau(n, ndm, ure, dre, vre, vle, dle, ule, g00, gt0, g0t, gtt, logdetqr, logdetql, logweightf, infoe )
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(inout) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      complex(dp), dimension(ndm,ndm), intent(out) :: g00, gt0, g0t, gtt
      complex(dp), intent(in) :: logdetqr, logdetql
      complex(dp), intent(out) :: logweightf
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      complex(dp), allocatable, dimension(:,:) :: vlhtr_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin
      complex(dp), allocatable, dimension(:,:) :: uutmp, vvtmp, dvvdtmp, Btmp
      complex(dp) :: zlogdet_tmp
      real(dp) :: rlogdet_tmp1, rlogdet_tmp2

      allocate( vlhtr_tmp(ndm,ndm) )
      allocate( drmax(ndm), drmin(ndm), dlmax(ndm), dlmin(ndm) )
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
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              vlhtr_tmp(i,j) = dconjg(vle(j,i))
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL

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

      !! >> g(t,0)
      call s_v_d_u( ndm, Btmp, drmin, vre, gt0 )


      call s_invlu_z(ndm,vre)       ! vre has been rewritten here.
      call s_invlu_z(ndm,vlhtr_tmp) ! vlhtr_tmp has been rewritten here.
      ! vvtmp = (vre*vle)^-1 = vle^-1 * vre^-1
      call zgemm('n','n',ndm,ndm,ndm,cone,vlhtr_tmp,ndm,vre,ndm,czero,vvtmp,ndm)  ! vvtmp = vle^-1 * vre^-1
      !! >> g(0,0)
      ! dlmax^-1 * ( vre * vle )^-1 drmax^-1
      ! dlmin * ( ule * ure ) * drmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              dvvdtmp(i,j) = vvtmp(i,j) / ( dlmax(i)*drmax(j) ) + dconjg(uutmp(j,i)) * dlmin(i) * drmin(j)
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_invlu_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, vre, drmax, dvvdtmp, Btmp )
      call s_v_invd_u( ndm, Btmp, dlmax, vlhtr_tmp, g00 )

      !! >> g(0,t)
      dlmin(:)=-dlmin(:)  ! note minus sign here
      call s_v_d_uH( ndm, Btmp, dlmin, ule, g0t )

      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( vlhtr_tmp )
      deallocate( uutmp )
      deallocate( vvtmp )
      deallocate( dvvdtmp )
      deallocate( Btmp )

    end subroutine green_tau

    !!!subroutine green_tau(nt, ndm, ure, dre, vre, vle, dle, ule, g00, gt0, g0t, gtt, infoe )
    !!!  implicit none
    !!!  integer, intent(in) :: nt, ndm
    !!!  complex(dp), dimension(ndm,ndm), intent(inout) :: ure, vre, vle, ule
    !!!  real(dp), dimension(ndm), intent(in) :: dre, dle
    !!!  complex(dp), dimension(ndm,ndm), intent(out) :: g00, gt0, g0t, gtt
    !!!  integer, intent(out) :: infoe

    !!!  ! local
    !!!  integer :: i, j
    !!!  complex(dp), allocatable, dimension(:,:) :: vrulmat, v2mat, u2mat, vlurmat, udvmat, vuvmat, uvumat
    !!!  real(dp), allocatable, dimension(:) :: d2vec

    !!!  allocate( vrulmat( 2*ndm, 2*ndm ) )   ! 1
    !!!  allocate(   v2mat( 2*ndm, 2*ndm ) )   ! 2
    !!!  allocate(   u2mat( 2*ndm, 2*ndm ) )   ! 3
    !!!  allocate( vlurmat( 2*ndm, 2*ndm ) )   ! 4
    !!!  allocate(  udvmat( 2*ndm, 2*ndm ) )   ! 5
    !!!  allocate(  vuvmat( 2*ndm, 2*ndm ) )   ! 6
    !!!  allocate(  uvumat( 2*ndm, 2*ndm ) )   ! 7
    !!!  allocate(   d2vec( 2*ndm ) )          ! 8

    !!!  call zgemm('n','n',ndm,ndm,ndm,cone,vre,ndm,vle,ndm,czero,vvtmp,ndm)  ! vvtmp = vre*vle
    !!!  call s_invlu_z(ndm,vvtmp)
    !!!  call zgemm('n','n',ndm,ndm,ndm,cone,ule,ndm,ure,ndm,czero,Atmp,ndm)  ! Atmp = ule*ure
    !!!  call s_invlu_z(ndm,Atmp)
    !!!  udvmat = czero
    !!!  do j = 1, ndm
    !!!      do i = 1, ndm
    !!!          udvmat(i,     j     ) =  vvtmp(i,j)
    !!!          udvmat(i+ndm, j+ndm ) =  Atmp(i,j)
    !!!          if(i.eq.j) then
    !!!              udvmat(i+ndm, j     ) = dcmplx( -dre(i), 0.d0 )
    !!!              udvmat(i    , j+ndm ) = dcmplx(  dle(i), 0.d0 )
    !!!          end if
    !!!      end do
    !!!  end do
    !!!  call s_svd_zg(2*ndm, 2*ndm, 2*ndm, udvmat, u2mat, d2vec, v2mat)

    !!!  ! check d2vec
    !!!  infoe = 0
    !!!  do i = 1, 2*ndm
    !!!      if( d2vec(i) .eq. 0.d0 ) then
    !!!          infoe = -1
    !!!          write(fout,'(a,i5,a)') 'WARNING!!! d2vec(',i, ' ) = 0 in green_tau !!! '
    !!!      end if
    !!!  end do

    !!!  call s_invlu_z(2*ndm, v2mat)
    !!!  call s_invlu_z(2*ndm, u2mat)

    !!!  !! attention here, we are now changing vre, ule, vle, ure
    !!!  call s_invlu_z(ndm,vre)
    !!!  call s_invlu_z(ndm,ule)
    !!!  call s_invlu_z(ndm,vle)
    !!!  call s_invlu_z(ndm,ure)

    !!!  vrulmat = czero
    !!!  vrulmat(1:ndm,1:ndm) = vre(1:ndm,1:ndm)
    !!!  vrulmat(ndm+1:2*ndm, ndm+1:2*ndm) = ule(1:ndm,1:ndm)

    !!!  vlurmat = czero
    !!!  vlurmat(1:ndm,1:ndm) = vle(1:ndm,1:ndm)
    !!!  vlurmat(ndm+1:2*ndm, ndm+1:2*ndm) = ure(1:ndm,1:ndm)

    !!!  call zgemm('n','n',2*ndm,2*ndm,2*ndm,cone,vrulmat,2*ndm,v2mat,  2*ndm,czero,vuvmat,2*ndm)
    !!!  call zgemm('n','n',2*ndm,2*ndm,2*ndm,cone,u2mat,  2*ndm,vlurmat,2*ndm,czero,uvumat,2*ndm)

    !!!  call s_v_invd_u( 2*ndm, vuvmat, d2vec, uvumat, udvmat )

    !!!  do j = 1, ndm
    !!!      do i = 1, ndm
    !!!          g00(i,j) = udvmat(i,j)
    !!!          gt0(i,j) = udvmat(i+ndm,j)
    !!!          g0t(i,j) = udvmat(i,j+ndm)
    !!!          gtt(i,j) = udvmat(i+ndm,j+ndm)
    !!!      end do
    !!!  end do

    !!!  deallocate(   d2vec )   ! 8
    !!!  deallocate(  uvumat )   ! 7
    !!!  deallocate(  vuvmat )   ! 6
    !!!  deallocate(  udvmat )   ! 5
    !!!  deallocate( vlurmat )   ! 4
    !!!  deallocate(   u2mat )   ! 3
    !!!  deallocate(   v2mat )   ! 2
    !!!  deallocate( vrulmat )   ! 1
    !!!end subroutine green_tau
