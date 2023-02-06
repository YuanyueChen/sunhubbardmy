module model_para
  use constants, only: dp, czero, cone, ctwo, cthree, cfour, pi
#IFDEF HONEYCOMB
  use honeycomb_lattice
#ENDIF
#IFDEF SQUARE
  use square_lattice
#ENDIF
#IFDEF CUBIC
  use cubic_lattice
#ENDIF
  use spring
  use dqmc_basic_data
  implicit none

  ! lattice
#IFDEF HONEYCOMB
  type(honeycomb) :: latt
#ENDIF
#IFDEF SQUARE
  type(square) :: latt
#ENDIF
#IFDEF CUBIC
  type(cubic) :: latt
#ENDIF
  integer, save :: ndim
  integer, save :: lq
  integer, save :: norb
  integer, save :: nzz, nnzz

  ! model
  integer, save :: ne
  real(dp), save :: beta
  real(dp), save :: dtau
  real(dp), save :: rt
  real(dp), save :: t2
  real(dp), save :: t3
  real(dp), save :: nu
  real(dp), save :: mu
  real(dp), save :: rhub
  real(dp), save :: rhub_plq
  real(dp), save :: alpha
  real(dp), save :: theta
  logical, save :: lproju
  logical, save :: lprojplqu
  real(dp), save :: xmag
  real(dp), save :: flux_x
  real(dp), save :: flux_y
  real(dp), save :: dimer
  integer, save :: ltrot
  integer, save :: nflr

  ! global para
  integer, parameter :: fout = 50
  integer, save :: irank, isize, ierr
  complex(dp), dimension(10), save :: main_obs, main_obs_recv, mpi_main_obs

  ! dqmc relative
  logical, save :: lwrapT, lwrapu, lwrapplqu, llocal, lstglobal, lworm
  logical, save :: lupdateu, lupdateplqu
  complex(dp), save :: phase, phase_tmp
  real(dp), save :: weight_track
  real(dp), save :: logweightf_old, logweightf_new, logweights_old, logweights_new
  complex(dp), save :: logweightf_up, logweightf_dn
  real(dp), save :: sgnm, nupdate

  ! cal. control
  integer, save :: nbin
  integer, save :: nsweep
  integer, save :: nst
  integer, save :: nwrap
  logical, save :: lstabilize
  logical, save :: lwarmup
  integer, save :: nwarmup
  integer, save :: n_outconf_pace 
  integer, allocatable, dimension(:) :: iwrap_nt
  integer, allocatable, dimension(:,:) :: wrap_step

  ! for dynamical
  logical, save :: ltau
  real(dp), save :: dyntau ! max dynamical tau
  integer, save :: obs_eqt_mid_len
  integer, save :: ntdm
  integer, save :: ntauin


  real(dp), save :: max_wrap_error, max_wrap_error_tmp, max_phasediff, max_phasediff_tmp
  real(dp), save :: xmax_dyn, xmax_dyn_tmp
  real(dp), save :: avglog10error, avgexpphasediff, logcount, logdyncount, avglog10dynerror

  ! delay update
  integer, save :: nublock

  contains

  subroutine make_tables
    implicit none

    include 'mpif.h'

    integer :: i, nwrap_mid
    real(dp) :: rtmp
    complex(dp) :: zw
    logical :: exists

#IFDEF CUBIC
    integer :: la, lb, lc
    real(dp) :: a1_p(3), a2_p(3), a3_p(3)
#ELSE
    integer :: la, lb
    real(dp) :: a1_p(2), a2_p(2)
#ENDIF

#IFDEF CUBIC
    namelist /model_para/ la, lb, lc, beta, dtau, rt, t2, t3, nu, mu, rhub, rhub_plq, alpha, theta, nflr, lprojplqu, lproju, xmag, flux_x, flux_y, dimer
#ELSE
    namelist /model_para/ la, lb, beta, dtau, rt, t2, t3, nu, mu, rhub, rhub_plq, alpha, theta, nflr, lprojplqu, lproju, xmag, flux_x, flux_y, dimer
#ENDIF
    namelist /ctrl_para/ lstabilize, nwrap, nsweep, nbin, nublock, ltau, dyntau, obs_eqt_mid_len

    ! default parameters
    la   = 2
    lb   = 2
#IFDEF CUBIC
    lc   = 2
#ENDIF
    beta = 20
    dtau = 0.05d0
    rt   = 1.d0
    t2   = 0.d0
    t3   = 0.d0
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
#IFDEF CUBIC
    call mpi_bcast( lc,           1, mpi_integer,  0, mpi_comm_world, ierr )
#ENDIF
    call mpi_bcast( beta,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( dtau,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rt,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( t2,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( t3,           1, mpi_real8,    0, mpi_comm_world, ierr )
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
#IFDEF CUBIC
    ! lattice
    lq = la*lb*lc
    ndim = lq
   	a1_p(1) = 1.0d0 ; a1_p(2) =  0.0d0; a1_p(3) = 0.d0;
    a2_p(1) = 0.0d0 ; a2_p(2) =  1.0d0; a2_p(3) = 0.d0;
    a3_p(1) = 0.0d0 ; a3_p(2) =  0.0d0; a3_p(3) = 1.d0;
    call setup_cubic(la,lb,lc,a1_p,a2_p,a3_p,latt)
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

end module model_para
