    subroutine dqmc_proj_stablize_b0_qr(ule,logwtmp)
      implicit none
      type(gfunc), intent(inout) :: ule
      type(rfunc), intent(out) :: logwtmp

      ! local
      type(gfunc) :: Utmp, Vtmp

      call allocate_gfunc( Utmp, ne, ndim )
      call allocate_gfunc( Vtmp, ne, ne )
      call s_zgeQR_pqmc(ne, ndim, ule%blk1, Utmp%blk1, Vtmp%blk1, logwtmp%blk1 )
      ule = Utmp
      call deallocate_gfunc( Vtmp )
      call deallocate_gfunc( Utmp )
    end subroutine dqmc_proj_stablize_b0_qr
