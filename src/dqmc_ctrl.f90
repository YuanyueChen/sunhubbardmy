module dqmc_ctrl
  use model_para
  use dqmc_config_h0
  use dqmc_config_plqu
  use dqmc_ft_basic_data
  use dqmc_ft_core
  implicit none


  ! main
  integer :: nbc, nsw
  character (len = 24) :: date_time_string
  real(dp) :: start_time, end_time, time1, time2

  contains

  subroutine make_tables
    implicit none

    include 'mpif.h'

    integer :: i, nwrap_mid
    real(dp) :: rtmp
    complex(dp) :: zw
    logical :: exists

    integer :: la, lb
    real(dp) :: a1_p(2), a2_p(2)

    namelist /model_para/ la, lb, beta, dtau, rt, nu, mu, rhub, rhub_plq, alpha, theta, nflr, lprojplqu, lproju, xmag, flux_x, flux_y, dimer
    namelist /ctrl_para/ lstabilize, nwrap, nsweep, nbin, nublock, ltau, dyntau, obs_eqt_mid_len

    ! default parameters
    la   = 2
    lb   = 2
    beta = 20
    dtau = 0.05d0
    rt   = 1.d0
    nu = 0.d0
    mu   = 0.d0 ! default is half filling
    rhub = 0.0d0
    rhub_plq = 0.0d0
    alpha = 0.2d0
    theta = 0.d0
    nflr = 2
    lprojplqu = .false.
    lproju = .false.
    xmag = 0.d0
    flux_x = 0.d0
    flux_y = 0.d0
    dimer = 0.d0

    lstabilize = .true.
    nwrap = 10
    nsweep = 20
    nbin = 10
    nublock = 16
    ltau = .false.
    dyntau = 0
    obs_eqt_mid_len = 1

    ! read parameters
    if ( irank.eq.0 ) then
        exists = .false.
        inquire (file = 'dqmc.in', exist = exists)
        if ( exists .eqv. .true. ) then
            open(unit=100, file='dqmc.in',status='unknown')
            read(100, model_para)
            read(100, ctrl_para)
            close(100)
        end if
    end if

!#ifdef MPI
    call mpi_bcast( la,           1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( lb,           1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( beta,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( dtau,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rt,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( nu,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( mu,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rhub,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rhub_plq,     1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( alpha,        1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( theta,        1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( nflr,         1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( lprojplqu,    1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( lproju,       1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( xmag,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( flux_x,       1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( flux_y,       1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( dimer,        1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( lstabilize,   1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( nwrap,        1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nsweep,       1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nbin,         1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nublock,      1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( ltau,         1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( dyntau,       1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( obs_eqt_mid_len,  1, mpi_integer,    0, mpi_comm_world, ierr )
    call mpi_barrier(mpi_comm_world,ierr)
!#endif

    ! tune parameters
    theta = theta*pi
    if( dble(int(nu))-nu .ne. 0.d0 ) then ! if nu is not integer
        lprojplqu = .false.
        if( irank == 0 ) then
            write(fout, '(a)') "nu is not integer, lprojplqu is set to be false!"
        end if
    end if
    if( dble(int(2.d0*nu))-2.d0*nu .ne. 0.d0 ) then ! if 2*nu is not integer
        lproju = .false.
        if( irank == 0 ) then
            write(fout, '(a)') "2*nu is not integer, lproju is set to be false!"
        end if
    end if

    if( rt .gt. 0.d0 ) then
        lwrapT = .true.
    else
        lwrapT = .false.
    end if

    if( rhub .ne. 0.d0 ) then
        lwrapu = .true.
    else
        lwrapu = .false.
    end if

    if( rhub_plq .ne. 0.d0 ) then
        lwrapplqu = .true.
    else
        lwrapplqu = .false.
    end if

    if( rhub .ne. 0.d0 ) then
        lupdateu = .true.
    else
        lupdateu = .false.
    end if

    if( rhub_plq .ne. 0.d0 ) then
        lupdateplqu = .true.
    else
        lupdateplqu = .false.
    end if


#IFDEF HONEYCOMB
    ! lattice
    lq = la*lb
    ndim = lq*2
   	a1_p(1) = 0.5d0 ; a1_p(2) = -0.5d0*dsqrt(3.d0)
    a2_p(1) = 0.5d0 ; a2_p(2) =  0.5d0*dsqrt(3.d0)
    call setup_honeycomb(la,lb,a1_p,a2_p,latt)
#ENDIF
#IFDEF SQUARE
    ! lattice
    lq = la*lb
    ndim = lq
   	a1_p(1) = 1.0d0 ; a1_p(2) =  0.0d0
    a2_p(1) = 0.0d0 ; a2_p(2) =  1.0d0
    call setup_square(la,lb,a1_p,a2_p,latt)
#ENDIF
    call print_latt( latt )

    ne = ndim/2+int(nu*dble(ndim)/8.d0)  ! number of particles
    ltrot = nint( beta / dtau )

#IFDEF PIFLUX
    xmag = dble(la*lb)*(a1_p(1)*a2_p(2) - a1_p(2)*a2_p(1))/2.d0
#ENDIF

    ! tune para for delay update
    if( lq/5 .lt. 16) then
        nublock = 4
    else if( lq/5 .lt. 32 ) then
        nublock = 8
    else if( lq/5 .lt. 64 ) then
        nublock = 16
    else if( lq/5 .lt. 256 ) then
        nublock = 32
    else ! equal to or greater than 256
        nublock = 64
    end if

    allocate( iwrap_nt(0:ltrot) )
    iwrap_nt(0:ltrot) = 0
    ! set nst, and wrap_step
    if( ltrot .lt. nwrap ) then
        nwrap = ltrot
        if(irank.eq.0) write(fout,'(a,i3,a)')  " WARNNING, ltrot is less than nwrap ", ltrot, ', then nwrap set to equal to ltrot '
    end if
    if( ltrot .lt. nwrap .or. (.not.lstabilize) ) then
        if(irank.eq.0) write(fout,'(a,i3,a)')  " WARNNING, ltrot is less than nwrap ", ltrot, ', do not need stablization '
        nst = 0
    else
        if( mod(ltrot,nwrap) .eq. 0 ) then
            nst = ltrot/nwrap
            allocate( wrap_step(2,nst) )
            do i = 1, nst
                iwrap_nt(i*nwrap) = i
                wrap_step(1,i) = 1+(i-1)*nwrap
                wrap_step(2,i) = i*nwrap
            end do
        else
            nst = ltrot/nwrap + 1
            allocate( wrap_step(2,nst) )
            do i = 1, nst-1
                iwrap_nt(i*nwrap) = i
                wrap_step(1,i) = 1+(i-1)*nwrap
                wrap_step(2,i) = i*nwrap
            end do
            i = nst
            iwrap_nt(ltrot) = i
            wrap_step(1,i) = (i-1)*nwrap+1
            wrap_step(2,i) = ltrot
        end if
    end if

    if(ltau) then
        ntdm = nint(dyntau/dtau)
        ntauin = (ltrot - ntdm)/2
        xmax_dyn = 0.d0
    end if

    weight_track = 0.d0

    allocate( Imat(ndim,ndim) )
    call s_identity_z(ndim,Imat)
    allocate( Imat_plq(latt%z_plq,latt%z_plq) )
    call s_identity_z(latt%z_plq,Imat_plq)
    allocate( Ivec(ndim) )
    do i = 1, ndim
        Ivec(i) = 1.d0
    end do
    sgnm = 1.d0
    phase = cone

  end subroutine make_tables

  subroutine deallocate_tables
    deallocate( wrap_step )
    deallocate( iwrap_nt )
  end subroutine deallocate_tables

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

  subroutine dqmc_start
    use spring                      !随机数
  
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
        write(fout,'(a)') '        The finite temperature determinant quantum monte carlo (DQMC) package '
        write(fout,*)
        write(fout,'(a)') '            FFFF   TTTTT   DDD      QQQ     M   M      CCCC                    '
        write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
        write(fout,'(a)') '            FFFF     T     D   D   Q   Q   M M M M    C                        '
        write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
        write(fout,'(a)') '            F        T     DDD      QQQ   M   M   M    CCCC                    '
        write(fout,'(a)') '                                      \                                       '
        write(fout,*)
        write(fout,*)
        write(fout,'(a)') ' written by Xiao Yan Xu ( wanderxu@gmail.com )                                '
        write(fout,*)
        write(fout,'(a)') ' history: '
        write(fout,*)
        write(fout,'(a)') '     22/01/2016,  version 1.0  '
        write(fout,*)
        write(fout,'(a)') ' ------------------------------------------------------------------------------------'
        write(fout,*)
        write(fout,'(a)') ' >>> The simulation start running at '//date_time_string
        if( isize .gt. 1 ) then
            write(fout,'(a,i6,a)') ' >>> Parallelism running with', isize, '  processes'
        else
            write(fout,'(a)') ' >>> Serial running '
        end if
    end if
#IFDEF _OPENMP
    start_time = omp_get_wtime()
#ELSE
    call cpu_time(start_time)
#ENDIF
    main_obs(:) = czero                 !检测程序的参数
    max_wrap_error = 0.d0
    max_phasediff = 0.d0
    logcount = 0.d0
    avglog10error = 0.d0
    avgexpphasediff = 0.d0

    if(ltau) then                       !是否做动力学测量
        xmax_dyn = 0.d0
        logdyncount = 0.d0
        avglog10dynerror = 0.d0
    end if
    call make_tables                    !输出表格, 读入dqmc.in

    call dqmc_initial_print             !打印模型参数

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
            call hconf%output_plqconf
            call u0conf%output_uconf
        end if
    else if( mod( nbc, max(nbin/3,1) ) .eq. 0 ) then
        call hconf%output_plqconf
        call u0conf%output_uconf
    end if

    if( irank.eq.0 .and. mod(nbc,max(nbin/10,1) ).eq.0 ) then
        write( fout, '(i5,a,i5,a)' ) nbc, '  /', nbin, '   finished '
    end if
    call wrap_error_print
    !!! --- END Timming and outconfc
  end subroutine dqmc_timing

  subroutine dqmc_initial_print
    implicit none
  
    integer :: i, j, nf, imj
  
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
    write(fout,'(a,f6.2)')    ' beta   = ', beta
    write(fout,'(a,f7.3)')    ' dtau   = ', dtau
    write(fout,'(a,i6)')      ' ndim = ', ndim
    write(fout,'(a,i6)')      ' nwrap  = ', nwrap
    write(fout,'(a,i6)')      ' nsweep = ', nsweep
    write(fout,'(a,i6)')      ' nbin = ', nbin
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
  
    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The lattice sites list '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)') '   i     list(i,:) '
    do i = 1, ndim
        write(fout,'(i6,2i4)') i,  latt%list(i,:)
    end do

    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The lattice sites cord '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)') '   i     cord(i,:) '
    do i = 1, ndim
        write(fout,'(i6,2i4)') i,  latt%cord(i,:)
    end do

    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The lattice nnlist '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)') '   i     nnlist(i,:) '
    do i = 1, ndim
        write(fout,'(i6,10i4)') i,  latt%nnlist(i,:)
    end do

    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The lattice snlist '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)') '   i     snlist(i,:) '
    do i = 1, ndim
        write(fout,'(i6,10i4)') i,  latt%snlist(i,:)
    end do

    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The lattice tnlist '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)') '   i     tnlist(i,:) '
    do i = 1, ndim
        write(fout,'(i6,10i4)') i,  latt%tnlist(i,:)
    end do

    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The lattice plq_cord '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)') '   i     plq_cord(:,i) '
    do i = 1, latt%ncell
        write(fout,'(i6,6i4)') i,  latt%plq_cord(:,i)
    end do

    write(fout,*)
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)')' The momentum listk '
    write(fout,'(a)')' --------------------- '
    write(fout,'(a)') '   i     listk(i,:) '
    do i = 1, latt%ncell
        write(fout,'(i6,2i4)') i,  latt%listk(i,:)
    end do
  
    if( latt%ncell < 50 ) then
    write(fout, *)
    write(fout,'(a)') '----------------------------'
    write(fout,'(a)') ' latt_imj info   '
    write(fout,'(a)') '----------------------------'
    write(fout, '(a)') '   i   j   imj '
    do j = 1, latt%ncell
        do i = 1, latt%ncell
            imj  = latt%imj(i,j)
            write(fout, '(3i4)') i, j, imj
        end do
    end do
    end if
  
  
    END IF
  
  end subroutine dqmc_initial_print

  subroutine dqmc_end
    implicit none
    include 'mpif.h'
    call hconf%output_plqconf
    call u0conf%output_uconf

    call mpi_reduce(main_obs, mpi_main_obs, size(main_obs), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    if(irank.eq.0) then
        if(lwrapplqu)  write(fout,'(a,e16.8)') ' >>> accep_plqu  = ', dble(mpi_main_obs(2))/aimag(mpi_main_obs(2))
        if(lwrapu)     write(fout,'(a,e16.8)') ' >>> accep_u     = ', dble(mpi_main_obs(3))/aimag(mpi_main_obs(3))
    end if

    call deallocate_obs
    call deallocate_core

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

end module dqmc_ctrl
