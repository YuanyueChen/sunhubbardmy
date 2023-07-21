    subroutine dqmc_ft_sweep_start_0b
      implicit none
      integer :: n, i, info
      real(dp) :: tmp

      type(gfunc) :: UR, VR
      type(dfunc) :: DRvec
      type(zfunc) :: logdetQR
      call allocate_gfunc(UR,ndim,ndim)
      call allocate_gfunc(VR,ndim,ndim)
      call allocate_dfunc(DRvec,ndim)

      IF ( nst .gt. 0 ) THEN

      ! at tau = 0
      Ust(0) = Ifmat
      Dst(0) = Ifvec
      Vst(0) = Ifmat
      logdetQst(0) = zfunc(czero)

      do n = 1, nst
          ! at tau = n * tau1
          call dqmc_ft_stablize_0b_qr(n)
#IFDEF TEST
          write(fout, '(a,i4,a)') ' in dqmc_ft_sweep_start_0b, Dst(', n, ')%blk1(:) = '
          write(fout,'(4(e16.8))') Dst(n)%blk1(:)
#ENDIF
      end do
  
      ! at tau = beta
      UR    = Ust(nst)
      DRvec = Dst(nst)
      VR    = Vst(nst)
      logdetQR = logdetQst(nst)
      call green_equaltimebb( nst, ndim, UR%blk1, DRvec%blk1, VR%blk1, gf%blk1, logdetQR%blk1, logweightf%blk1, info )

      call set_phase( logweightf, phase )
!#IFDEF TEST
      if( irank == 0 ) then
          write(fout,'(a,2e16.8)') 'in dqmc_ft_sweep_start_0b, phase = ', phase
      end if
!#ENDIF

      if( info .eq. -1 ) then
          write(fout,'(a)') ' WRONG in sweep_start, exit '
          stop
      end if

#IFDEF TEST_LEVEL3
      write(fout, '(a)') ' After sweep_start_0b, gf%blk1(:,:) = '
      do i = 1, ndim
          write(fout,'(18(f8.4))') dble(gf%blk1(i,:))
      end do
#ENDIF

      ELSE

      gf = Ifmat
      call Bmat_left_forward( ltrot, 1, gf)
      do  i = 1, ndim
          gf%blk1(i,i) = gf%blk1(i,i) + cone
      end do
      call s_invlu_z( ndim, gf%blk1 )

      END IF

      call deallocate_gfunc(UR)
      call deallocate_gfunc(VR)
      call deallocate_dfunc(DRvec)
  
    end subroutine dqmc_ft_sweep_start_0b
