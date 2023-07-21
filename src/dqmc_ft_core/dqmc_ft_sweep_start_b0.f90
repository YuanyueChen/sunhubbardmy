    subroutine dqmc_ft_sweep_start_b0
      implicit none
      integer :: n, i, info
      real(dp) :: tmp

      type(gfunc) :: UL, VL
      type(dfunc) :: DLvec
      type(zfunc) :: logdetQL
      call allocate_gfunc(UL,ndim,ndim)
      call allocate_gfunc(VL,ndim,ndim)
      call allocate_dfunc(DLvec,ndim)

      IF ( nst .gt. 0 ) THEN

      ! at tau = beta
      Ust(nst) = Ifmat
      Dst(nst) = Ifvec
      Vst(nst) = Ifmat
      logdetQst(nst) = zfunc(czero)

      do n = nst, 1, -1
          ! at tau = (n-1) * tau1
          ! calculate B(n*tau1,(n-1)*tau1), and set Vst(:,:,n-1), Dst(:,:,n-1), Ust(:,:,n-1)
          call dqmc_ft_stablize_b0_qr(n)
#IFDEF TEST
          write(fout, '(a,i4,a)') ' in dqmc_ft_sweep_start_b0, Dst(', n, ')%blk1(:) = '
          write(fout,'(4(e16.8))') Dst(n)%blk1(:)
#ENDIF
      end do

      ! at tau = 0
      UL    = Ust(0)
      DLvec = Dst(0)
      VL    = Vst(0)
      logdetQL = logdetQst(0)
      call green_equaltime00( nst, ndim, VL%blk1, DLvec%blk1, UL%blk1, gf%blk1, logdetQL%blk1, logweightf%blk1, info )

      call set_phase( logweightf, phase)
#IFDEF TEST
      if( irank == 0 ) then
          write(fout,'(a,2e16.8)') 'in dqmc_ft_sweep_start_b0, phase = ', phase
      end if
#ENDIF

      if( info .eq. -1 ) then
          write(fout,'(a)') ' WRONG in sweep_start, exit '
          stop
      end if

#IFDEF TEST_LEVEL3
      write(fout, '(a)') ' After sweep_start, gf%blk1(:,:) = '
      do i = 1, ndim
          write(fout,'(18(2e12.4))') gf%blk1(i,:)
      end do
#ENDIF

      ELSE

      gf = Ifmat
      call Bmat_left_forward( ltrot, 1, gf )
      do  i = 1, ndim
          gf%blk1(i,i) = gf%blk1(i,i) + cone
      end do
      call s_invlu_z( ndim, gf%blk1 )

      END IF

      call deallocate_gfunc(UL)
      call deallocate_gfunc(VL)
      call deallocate_dfunc(DLvec)
  
    end subroutine dqmc_ft_sweep_start_b0
