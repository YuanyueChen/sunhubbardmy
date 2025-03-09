    subroutine dqmc_proj_sweep_start_0b
      implicit none

      integer :: n, nt, nf, nflag, nl, nr, i
      type(zfunc) :: zlogwtmp
      type(rfunc) :: logwtmp

      ! UR at tau = 0
      do nl = 1, ne
         do i = 1,ndim
            UR%blk1(i,nl) = projR%blk1(i,nl)
            Ust(nst)%blk1(i,nl) = dconjg(projL%blk1(i,nl))
         end do
      end do

      logweightf = zfunc(czero)
      logwDV(0) = rfunc(0.d0)
      do n = 1, nst
          call Bmat_left_forward( wrap_step(2,n), wrap_step(1,n), UR )
          call dqmc_proj_stablize_0b_qr(UR,logwtmp)
          !n = nt/nwrap
          Ust(n) = UR
          logwDV(n) = logwDV(n-1) + logwtmp
      end do

      ! now UR at tau = beta

      ! UL at tau = beta
      do nl = 1, ne
          do i = 1, ndim
              UL%blk1(nl,i) = dconjg(projL%blk1(i,nl))
          end do
      end do

      ! calculate ULRINV for calculating green function
      call zgemm('n','n',ne,ne,ndim,cone,UL%blk1,ne,UR%blk1,ndim,czero,ULR%blk1,ne)  ! ULR = UL*UR
      ULRINV = ULR
      !call s_inv_z(ne,ULRINV)
      call s_inv_logdet_lu_z(ne,ULRINV%blk1,zlogwtmp%blk1)
      call set_phase(zlogwtmp, phase )
      logweightf%blk1 = dcmplx(logwDV(nst)%blk1,0.d0) + zlogwtmp%blk1 ! the true weight of fermion determinant
!#IFDEF TEST
      if( irank == 0 ) then
          write(fout,'(a,2e16.8)') 'in dqmc_proj_sweep_start_0b, phase = ', phase
      end if
!#ENDIF

    end subroutine dqmc_proj_sweep_start_0b
