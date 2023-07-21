    subroutine dqmc_proj_stablize_0b_qr(ure,logwtmp)
      implicit none
      type(gfunc), intent(inout) :: ure
      type(rfunc), intent(out) :: logwtmp

      ! local
      type(gfunc) :: Utmp, Vtmp
#IFDEF TEST
      integer :: i
#ENDIF

      call allocate_gfunc(Utmp, ndim, ne)
      call allocate_gfunc(Vtmp, ne, ndim)
      call s_zgeQR_pqmc(ndim, ne, ure%blk1, Utmp%blk1, Vtmp%blk1, logwtmp%blk1 )
      ure = Utmp
#IFDEF TEST_LEVEL3
      write(fout,*)
      write(fout,'(a)') ' in qr 0->beta, ure(:,:) = '
      do i = 1, ndim
          write(fout,'(18(2f7.4))') ure%blk1(i,:)
      end do
      write(fout,*)
      write(fout,'(a)') ' Vtmp(:) = '
      do i = 1, ne
          write(fout,'(18(2f7.4))') Vtmp%blk1(i,:)
      end do
      write(fout,'(a,e24.12)') ' logwtmp%blk1 = ', logwtmp%blk1
#ENDIF
      call deallocate_gfunc( Utmp )
      call deallocate_gfunc( Vtmp )
    end subroutine dqmc_proj_stablize_0b_qr
