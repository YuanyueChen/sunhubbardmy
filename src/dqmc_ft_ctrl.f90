module dqmc_ctrl
  use model_para
  use dqmc_util
  use dqmc_ft_core
  implicit none

  contains

  subroutine dqmc_start
    use spring
  
    ! prepare for the DQMC

    ! H0 part
    if(lwrapT) then
        call h0c%set_h0conf(lq, ltrot, rt, mu)
    end if

    ! plaquette interaction
    if(lwrapplqu) then
        call hconf%set_plqconf(lq, ltrot, nu, lprojplqu, alpha, theta, rhub_plq, dtau )
        call hconf%input_plqconf
    end if

    ! on-site interaction
    if(lwrapu) then
        call u0conf%set_conf(lq, ltrot, nu, lproju, alpha, theta, rhub, dtau )
        call u0conf%input_conf
    end if

    ! nearest neighbor interaction
    if(lwrapv) then
        call v0conf%set_conf(lq, ltrot, nu, lprojv, alpha, theta, rv, dtau )
        call v0conf%input_conf
    end if

    call allocate_core
    ! allocate observables
    call allocate_obs
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
#IFDEF TIMING
    real(dp) :: starttime, endtime, time1
#ENDIF

    include 'mpif.h'
#IFDEF TIMING
     call cpu_time_now(starttime)
#ENDIF
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
#IFDEF TIMING
     call cpu_time_now(time1)
#ENDIF
    call equaltime_output  ! reduce
    if(ltau) call dyn_output
#IFDEF TIMING
    call cpu_time_now(endtime)
    timecalculation(6)=timecalculation(6)+endtime-time1
    timecalculation(2)=timecalculation(2)+endtime-starttime
#ENDIF
  end subroutine dqmc_sweep

end module dqmc_ctrl
