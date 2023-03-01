module dqmc_util
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  use dqmc_para
  implicit none

  integer :: nbc, nsw
  character (len = 24) :: date_time_string
  real(dp) :: start_time, end_time, time1, time2

  contains

  subroutine wrap_error_print
      implicit none
      include 'mpif.h'
      call mpi_reduce( max_wrap_error, max_wrap_error_tmp, 1, mpi_real8, mpi_max, 0, mpi_comm_world, ierr )
      call mpi_reduce( max_phasediff, max_phasediff_tmp, 1, mpi_real8, mpi_max, 0, mpi_comm_world, ierr )
      max_wrap_error = max_wrap_error_tmp
      max_phasediff = max_phasediff_tmp

      call mpi_reduce( logcount, max_wrap_error_tmp, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr )
      logcount = max_wrap_error_tmp

      call mpi_reduce( avglog10error, max_wrap_error_tmp, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr )
      call mpi_reduce( avgexpphasediff, max_phasediff_tmp, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr )
      avglog10error = max_wrap_error_tmp / logcount
      avgexpphasediff = max_phasediff_tmp / logcount
      if(irank.eq.0) then
          write(fout, '(a,e16.8)') '--------------max_wrap_error  = ', max_wrap_error
          write(fout, '(a,e16.8)') '--------------max_phasediff   = ', max_phasediff
          write(fout, '(a,f16.8)') '--------------avglog10error     = ', avglog10error
          write(fout, '(a,f16.8)') '--------------avgexpphasediff = ', avgexpphasediff
      end if
      max_wrap_error = 0.d0
      max_phasediff = 0.d0
      logcount = 0.d0
      avglog10error = 0.d0
      avgexpphasediff = 0.d0
      if( ltau ) then
          call mpi_reduce( xmax_dyn, xmax_dyn_tmp, 1, mpi_real8, mpi_max, 0, mpi_comm_world, ierr )
          xmax_dyn = xmax_dyn_tmp
          call mpi_reduce( logdyncount, xmax_dyn_tmp, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr )
          logdyncount = xmax_dyn_tmp
          call mpi_reduce( avglog10dynerror, xmax_dyn_tmp, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr )
          avglog10dynerror = xmax_dyn_tmp / logdyncount

          if(irank.eq.0) then
              write(fout,'(a,e16.8)')'--------------xmax_dyn        = ', xmax_dyn
              write(fout,'(a,f16.8)')'--------------avglog10dynerror  = ', avglog10dynerror
          end if
          xmax_dyn = 0.d0 ! in warmup, xmax_dyn is not right, reset it here
          logdyncount = 0.d0
          avglog10dynerror = 0.d0
      end if
  end subroutine wrap_error_print

  subroutine dqmc_initial
  
    integer :: system_time
    integer :: stream_seed
    character (len = 24) :: date_time_string
  
    !================================================
    !%% inital the pseudo random number generator   $
    !------------------------------------------------
    call system_clock(system_time)
    stream_seed = abs( system_time - ( irank * 1981 + 2008 ) * 951049 )
#IFDEF TEST
    stream_seed = abs( 0 - ( irank * 1981 + 2008 ) * 951049 )
    write(fout, '(a,i20)') ' stream_seed = ', stream_seed
#ENDIF
    call spring_sfmt_init(stream_seed)
  
    call fdate(date_time_string)
  
    ! print head
    if(irank.eq.0) then
  
        write(fout,'(a)') ' ===================================================================================='
        write(fout,*)
        write(fout,'(a)') '        The determinant quantum monte carlo (DQMC) package '
        write(fout,*)
        write(fout,'(a)') '            DDD      QQQ       M   M      CCCC                 '
        write(fout,'(a)') '            D  D    Q   Q     M M M M    C                     '
        write(fout,'(a)') '            D   D   Q   Q     M M M M    C                     '
        write(fout,'(a)') '            D  D    Q   Q     M M M M    C                     '
        write(fout,'(a)') '            DDD      QQQ Q   M   M   M    CCCC                 '
        write(fout,'(a)') '                                                               '
        write(fout,*)
        write(fout,*)
        write(fout,'(a)') ' written by Xiao Yan Xu ( wanderxu@gmail.com )                                '
        write(fout,*)
        write(fout,'(a)') ' history: '
        write(fout,*)
        write(fout,'(a)') '     22/01/2016,  version 1.0  '
        write(fout,'(a)') '     05/02/2023,  version 1.x  '
        write(fout,*)
        write(fout,'(a)') ' ------------------------------------------------------------------------------------'
        write(fout,*)
        write(fout,'(a)') ' >>> The simulation start running at '//date_time_string
        if( isize .gt. 1 ) then
            write(fout,'(a,i6,a)') ' >>> Parallelism running with', isize, '  processes'
        else
            write(fout,'(a)') ' >>> Serial running '
        end if
#IFDEF _OPENMP
        !$OMP PARALLEL
        if( OMP_GET_THREAD_NUM() == 0 ) then
            write(fout,'(a,i6)') ' >>> OMP_NUM_THREADS = ', OMP_GET_NUM_THREADS()
        end if
        !$OMP END PARALLEL
#ENDIF
    end if
#IFDEF _OPENMP
    start_time = omp_get_wtime()
#ELSE
    call cpu_time(start_time)
#ENDIF
    main_obs(:) = czero
#IFDEF TIMING
    timecalculation(:)=0.d0
#ENDIF
    max_wrap_error = 0.d0
    max_phasediff = 0.d0
    logcount = 0.d0
    avglog10error = 0.d0
    avgexpphasediff = 0.d0

    if(ltau) then
        xmax_dyn = 0.d0
        logdyncount = 0.d0
        avglog10dynerror = 0.d0
    end if
    call make_tables

    call dqmc_initial_print

  end subroutine dqmc_initial

  subroutine dqmc_timing
    implicit none
    include 'mpif.h'
    !!! --- Timming and outconfc
    if( nbc .eq. 1 )  then
#IFDEF _OPENMP
        time2 = omp_get_wtime()
#ELSE
        call cpu_time(time2)
#ENDIF
        if(irank.eq.0) then
            n_outconf_pace = nint( dble( 3600 * 12 ) / ( time2-time1 ) )
            if( n_outconf_pace .lt. 1 ) n_outconf_pace = 1
            write(fout,'(a,e16.8,a)') ' time for 1 bin: ', time2-time1, ' s'
            write(fout,'(a,i12)') ' n_out_conf_pace = ', n_outconf_pace
        end if
        call mpi_bcast( n_outconf_pace, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr )
    end if

    if( n_outconf_pace .lt. nbin/3 ) then
        if( mod(nbc,n_outconf_pace) .eq. 0 ) then
            if(lwrapplqu) call hconf%output_plqconf
            if(lwrapu)    call u0conf%output_uconf
        end if
    else if( mod( nbc, max(nbin/3,1) ) .eq. 0 ) then
        if(lwrapplqu) call hconf%output_plqconf
        if(lwrapu)    call u0conf%output_uconf
    end if

    if( irank.eq.0 .and. mod(nbc,max(nbin/10,1) ).eq.0 ) then
        write( fout, '(i5,a,i5,a)' ) nbc, '  /', nbin, '   finished '
    end if
    call wrap_error_print
    !!! --- END Timming and outconfc
  end subroutine dqmc_timing

  subroutine dqmc_initial_print
    implicit none
  
    integer :: i
  
    IF(irank.eq.0) THEN
  
    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The input parameters  '
    write(fout,'(a)')' --------------------- '
    write(fout,*)
    write(fout,'(a,f6.2)')    ' t      = ', rt
    write(fout,'(a,f7.3)')    ' mu     = ', mu
    write(fout,'(a,f7.3)')    ' rhub_plq  = ', rhub_plq
    write(fout,*)             ' lprojplqu = ', lprojplqu
    write(fout,'(a,f7.3)')    ' rhub      = ', rhub
    write(fout,*)             ' lproju    = ', lproju
    write(fout,'(a,i4)')      ' la     = ', latt%l1
    write(fout,'(a,i4)')      ' lb     = ', latt%l2
#IFDEF CUBIC
    write(fout,'(a,i4)')      ' lc     = ', latt%l3
#ENDIF
    write(fout,'(a,f6.2)')    ' beta   = ', beta
    write(fout,'(a,f7.3)')    ' dtau   = ', dtau
    write(fout,'(a,i6)')      ' ndim = ', ndim
    write(fout,'(a,i6)')      ' ne   = ', ne
    write(fout,'(a,i6)')      ' nwrap  = ', nwrap
    write(fout,'(a,i6)')      ' nsweep = ', nsweep
    write(fout,'(a,i6)')      ' nbin = ', nbin
    write(fout,'(a,i6)')      ' nublock = ', nublock
    write(fout,'(a,f6.2)')    ' B      = ', xmag
    write(fout,'(a,f8.5)')    ' flux_x = ', flux_x
    write(fout,'(a,f8.5)')    ' flux_y = ', flux_y
    write(fout,'(a,f8.5)')    ' dimer  = ', dimer
    write(fout,'(a,f8.5)')    ' alpha  = ', alpha
    write(fout,'(a,f8.5)')    ' theta  = ', theta
    write(fout,'(a,i4)')      ' nflr   = ', nflr
    write(fout,*)  ' ltau = ', ltau
  
    if( ltau ) then
        write(fout, '(a,i6)') ' ntauin = ', ntauin
        write(fout, '(a,i6)') ' ntdm = ', ntdm
    end if
  
    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' wrapping coordinates  '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')'       wrap_step(1,i)   wrap_step(2,i)   iwrap_nt(nt) '
    do i = 1, nst
        write( fout, '(3i16)') wrap_step(1,i), wrap_step(2,i), iwrap_nt( wrap_step(2,i) )
    end do

    write(fout,'(a)') ' '
#IFDEF HONEYCOMB
    write(fout,'(a)') 'honeycomb_lattice'
#ENDIF
#IFDEF SQUARE
    write(fout,'(a)') 'square_lattice'
#ENDIF
#IFDEF CUBIC
    write(fout,'(a)') 'cubic_lattice'
#ENDIF
    call latt%print_latt(fout)
  
    END IF
  
  end subroutine dqmc_initial_print

  subroutine dqmc_end
    implicit none
    include 'mpif.h'
    if(lwrapplqu) call hconf%output_plqconf
    if(lwrapu)    call u0conf%output_uconf

    call mpi_reduce(main_obs, mpi_main_obs, size(main_obs), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    if(irank.eq.0) then
        if(lwrapplqu)  write(fout,'(a,e16.8)') ' >>> accep_plqu  = ', dble(mpi_main_obs(2))/aimag(mpi_main_obs(2))
        if(lwrapu)     write(fout,'(a,e16.8)') ' >>> accep_u     = ', dble(mpi_main_obs(3))/aimag(mpi_main_obs(3))
    end if

    call deallocate_tables
    if( irank.eq.0 ) then
#IFDEF _OPENMP
        end_time = omp_get_wtime()
#ELSE
        call cpu_time(end_time)
#ENDIF
        call fdate( date_time_string )
        write(fout,*)
        write(fout,'(a,f10.2,a)') ' >>> Total time spent:', end_time-start_time, 's'
        write(fout,'(a)') ' >>> Happy ending at '//date_time_string
        write(fout,*)
#IFDEF TIMING
        write(fout,'(a,f10.2,a)') 'The total time of sweep:', timecalculation(1), 's'  
        write(fout,'(a,f10.2,a)') 'The time of sweep in calculation:', timecalculation(2), 's'
        write(fout,'(a,f10.2,a)') 'The time of sweep in warmup:', timecalculation(3), 's' 
        write(fout,'(a,f10.2,a)') 'The time of measuring equaltime:', timecalculation(4), 's'  
        write(fout,'(a,f10.2,a)') 'The time of measuring unequaltime:', timecalculation(5), 's'
        write(fout,'(a,f10.2,a)') 'The time of outputing each bin:', timecalculation(6), 's'     
        write(fout,'(a,f10.2,a)') 'The total time of numerical stability:', timecalculation(7), 's'
        write(fout,'(a,f10.2,a)') 'The time of H0 left:', timecalculation(8), 's'
        write(fout,'(a,f10.2,a)') 'The time of H0 right:', timecalculation(9), 's'
        write(fout,'(a,f10.2,a)') 'The time of u left:', timecalculation(10), 's'
        write(fout,'(a,f10.2,a)') 'The time of u right:', timecalculation(11), 's'
        write(fout,'(a,f10.2,a)') 'The time of plqu left:', timecalculation(12), 's'
        write(fout,'(a,f10.2,a)') 'The time of plqu right:', timecalculation(13), 's'
        if( timecalculation(14).ne.0.d0 ) then
            write(fout,'(a,f10.2,a)') 'The time of update without delay update:', timecalculation(14), 's'
        end if
        if( timecalculation(15).ne.0.d0 ) then
            write(fout,'(a,f10.2,a)') 'The time of update with delay update:', timecalculation(15), 's'
        end if
#ENDIF      
        write(fout,'(a)') ' The simulation done !!! '
        write(fout,*)
        write(fout,'(a)') '        o         o    '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '        o         o    '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '        o         o    '
    end if
  
    close(fout)
  end subroutine dqmc_end


end module dqmc_util
