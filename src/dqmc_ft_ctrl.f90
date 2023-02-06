module dqmc_ctrl
  use model_para
  use dqmc_util
  use dqmc_ft_core
  implicit none

  contains

  subroutine dqmc_start
    use spring
  
    ! prepare for the DQMC
    call h0c%set_h0conf(lq, ltrot, rt, mu)                  !h0c:H0
!    call hconf%set_plqconf(lq, ltrot, nu, lprojplqu, alpha, theta, rhub_plq, dtau )                     !plaquette interaction
!    call hconf%input_plqconf
    call u0conf%set_uconf(lq, ltrot, nu, lproju, alpha, theta, rhub, dtau )                             !on-site interaction
    call u0conf%input_uconf

    call allocate_core
    call allocate_obs                                                   !observable
    call dqmc_ft_sweep_start_0b

    if( irank .eq. 0 ) then
        write(fout,'(a)') ' dqmc_ft_sweep_start done '
    end if
  end subroutine dqmc_start

  subroutine dqmc_core_deallocate
    call deallocate_obs
    call deallocate_core
  end subroutine dqmc_core_deallocate

  subroutine dqmc_warmup
    implicit none
    ! warmup
    !lwarmup = .false.
    if( lwarmup ) then
        ! set nwarmup
        nwarmup = nsweep
#IFDEF TEST
        nwarmup = 0
#ENDIF
        if( irank.eq.0 ) then
            write(fout,'(a,i8)') ' nwarmup = ', nwarmup
        end if
        do nsw = 1, nwarmup
            if(llocal) then
                call dqmc_ft_sweep_b0(lupdate=.true., lmeasure_equaltime=.false.)
                call dqmc_ft_sweep_0b(lupdate=.true., lmeasure_equaltime=.false., lmeasure_dyn=.false.)
            end if
        end do
  
        if(irank.eq.0) then
            write(fout, '(a)') 'after warmup : '
        end if
        call wrap_error_print
    end if
  end subroutine dqmc_warmup

  subroutine dqmc_sweep
    implicit none
    include 'mpif.h'
    call obser_init
    do nsw = 1, nsweep
        !! only perform local update, measure equaltime quantities and dyn quantities when turnning off updates
  
        !!call dqmc_ft_sweep_b0(lupdate=.true., lmeasure_equaltime=.true.)
        call dqmc_ft_sweep_b0(lupdate=.true., lmeasure_equaltime=.true.)
  
        ! dynamic measurement
        if( ltau) then
            call push_stage
            call dqmc_ft_sweep_0b(lupdate=.false., lmeasure_equaltime=.false., lmeasure_dyn=ltau)
            call pop_stage
        end if
  
        call dqmc_ft_sweep_0b(lupdate=.true., lmeasure_equaltime=.true., lmeasure_dyn=.false.)
  
        !call dqmc_ft_sweep_0b(lupdate=.true., lmeasure_equaltime=.true., lmeasure_dyn=.false.)
        !!!call dqmc_ft_sweep_0b(lupdate=.true., lmeasure_equaltime=.false., lmeasure_dyn=.false.)
        !!call dqmc_ft_sweep_0b(lupdate=.false., lmeasure_equaltime=.true., lmeasure_dyn=.false.)
  
#IFDEF TEST
        if( irank .eq. 0 ) then
            write(fout,'(a,i4,i4,a)') ' dqmc_sweep ', nbc, nsw,  '  done'
        end if
#ENDIF
    end do
    call equaltime_output  ! reduce
    if(ltau) call dyn_output
  end subroutine dqmc_sweep

end module dqmc_ctrl
