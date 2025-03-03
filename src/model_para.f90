module model_para
  ! in this module, we define model related parameters and calculation control parameters
  ! The model related parameters include lattice, orbitals, hopping, interaction, ...
  ! The control parameters includes  parameters for mpi, parameters for numerical stablization, parameters for update, ...
  ! It also contains several subroutines:
  ! subroutine make_tables, read parameters from input file, and set up parameters
  ! subroutine deallocate_tables, deallocate parameter arrays
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
  use spring ! random number generator
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
  integer, save :: ndim ! typical dimension of matrix, usually it is the total number of sites of the system
  integer, save :: lq   ! total number of unit cell
  integer, save :: norb ! number of orbital*sublattice (e.g. for single orbital model, 1 for square, 2 for honeycomb)
  integer, save :: nzz, nnzz ! coordination number of the nearest neighbor and next nearest neighbor

  ! model
  integer, save :: ne ! number of particles, ne/ndim is the filling factor
  real(dp), save :: beta ! inverse temperature
  real(dp), save :: dtau ! time step
  real(dp), save :: rt   ! nearest neighbor hopping
  real(dp), save :: t2   ! second nearest neighbor hopping
  real(dp), save :: t3   ! third nearest neighbor hopping
  real(dp), save :: nu   ! reference filling factor (0 for half-filling)
  real(dp), save :: mu   ! chemical potential
  real(dp), save :: rhub  ! Hubbard onsite interaction
  real(dp), save :: rhub_plq ! plaquette interaction
  real(dp), save :: rv       ! nn interaction
  real(dp), save :: alpha    ! parameter in front of asistant hopping term
  real(dp), save :: theta    ! phase factor in asistant hopping term
  logical, save :: lprojv     ! whether using infinite interaction projection for nn interaction
  logical, save :: lproju     ! whether using infinite interaction projection for Hubbard onsite interaction
  logical, save :: lprojplqu  ! whether using infinite interaction projection for plaquette interaction
  real(dp), save :: xmag  ! z-direction magnetic field
  real(dp), save :: flux_x ! twisting flux, magnetic field along x-direction
  real(dp), save :: flux_y ! twisting flux, magnetic field along y-direction
  real(dp), save :: dimer  ! strength of dimerization
  integer, save :: ltrot   ! total number of time slices
  integer, save :: nflr    ! total number of fermion flavors (2 for spin-1/2 fermion systems)
  real(dp), save :: rndness  ! a number control the randomness added to hopping when we generate trial wave function in PQMC

  ! global para
  integer, parameter :: fout = 50  ! tag for main output
  integer, save :: irank ! mpi parameters, which process id
  integer, save :: isize ! mpi parameters, total number process
  integer, save :: ierr  ! mpi parameters, error tag
  complex(dp), dimension(10), save :: main_obs, main_obs_recv, mpi_main_obs ! global observables, monitoring running status
#IFDEF TIMING
  real(dp), dimension(30),save :: timecalculation  ! calculate the time spent in update, measure, sweep and etc.
  real(dp)::starttime1,endtime1 !use these to do the time calculation
#ENDIF
  ! dqmc relative
  logical, save :: lwrapT ! whether wrap hopping (false when rt=0, true when rt \= 0)
  logical, save :: lwrapu ! whether wrap onsite interaction (false when rhub=0, true when rhub\=0)
  logical, save :: lwrapv ! whether wrap nn interaction (false when rv=0, true when rv\=0)
  logical, save :: lwrapplqu ! whether wrap plaquette interaction (false when rhub_plq=0, true when rhub_plq\=0)
  logical, save :: llocal ! whether use local update
  logical, save :: lstglobal ! whether use space-time global update
  logical, save :: lworm     ! whether use worm update
  logical, save :: lupdateu  ! whether update configurations for onsite interaction
  logical, save :: lupdatev  ! whether update configurations for nn interaction
  logical, save :: lupdateplqu ! whether update configurations for plaquette interaction
  complex(dp), save :: phase, phase_tmp ! phase factor of determinant
  real(dp), save :: weight_track ! track the weight of determinant
  real(dp), save :: logweightf_old, logweightf_new ! log weight of fermion part
  real(dp), save :: logweights_old, logweights_new ! log weight of bosonic part

  ! cal. control
  integer, save :: nbin  ! total number of bins
  integer, save :: nsweep ! number of sweeps in each bin
  integer, save :: nst    ! number of stablization from 0 to beta (or from beta to 0)
  integer, save :: nwrap  ! time step for stablization
  logical, save :: lstabilize ! whether use stablization
  logical, save :: lwarmup ! whether use warmup
  integer, save :: nwarmup ! number of sweeps for warmup
  integer, save :: n_outconf_pace  ! pace for output configurations
  integer, allocatable, dimension(:) :: iwrap_nt    ! for the control of stablization, >0 at the time point of stablization
  integer, allocatable, dimension(:,:) :: wrap_step ! for the control of stablization, begin and end time point for each stablization

  ! for dynamical
  logical, save :: ltau ! whether measure unequal-time correlations
  real(dp), save :: dyntau ! max dynamical tau
  integer, save :: obs_eqt_mid_len ! length of equal-time observation time points for projection QMC
  integer, save :: ntdm   ! number of time slices for unequal-time measurements
  integer, save :: ntauin ! start point for unequal-time measurements


  real(dp), save :: max_wrap_error, max_wrap_error_tmp ! record max wrap error
  real(dp), save :: max_phasediff, max_phasediff_tmp   ! record max wrap phase error
  real(dp), save :: xmax_dyn, xmax_dyn_tmp             ! record max dyn error
  real(dp), save :: avglog10error    ! average of log10 of max wrap error
  real(dp), save :: avgexpphasediff  ! average of log10 of max wrap phase error
  real(dp), save :: avglog10dynerror ! average of log10 of max dyn error
  real(dp), save :: logcount      ! counter for calculating average of log10 of max wrap error
  real(dp), save :: logdyncount   ! counter for calculating average of log10 of max dyn error

  ! delay update
  integer, save :: nublock ! size of block when doing delay update, if nublock = 0 in input file, it uses default setting based on system size
  integer, save :: nustock ! size of stock array for submatrix-LR update for projective QMC

  contains

  subroutine make_tables
    implicit none

    include 'mpif.h'

    integer :: i, nwrap_mid
    real(dp) :: rtmp
    complex(dp) :: zw
    logical :: exists

#IFDEF CUBIC
    integer :: la, lb, lc ! lattice size in each direction, for 3d
    real(dp) :: a1_p(3), a2_p(3), a3_p(3) ! primitive vectors, for 3d
#ELSE
    integer :: la, lb ! lattice size in each direction, for 2d
    real(dp) :: a1_p(2), a2_p(2) ! primitive vectors, for 2d
#ENDIF

    ! namelist is used for reading parameters from input file
#IFDEF CUBIC
    namelist /model_para/ la, lb, lc, beta, dtau, rt, t2, t3, nu, mu, rhub, rv, rhub_plq, alpha, theta, nflr, lprojplqu, lprojv, lproju, xmag, flux_x, flux_y, dimer, rndness
#ELSE
    namelist /model_para/ la, lb, beta, dtau, rt, t2, t3, nu, mu, rhub, rv, rhub_plq, alpha, theta, nflr, lprojplqu, lprojv, lproju, xmag, flux_x, flux_y, dimer, rndness
#ENDIF
    namelist /ctrl_para/ lstabilize, nwrap, nsweep, nbin, nublock, nustock, ltau, dyntau, obs_eqt_mid_len
    
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
    rv   = 0.0d0
    rhub_plq = 0.0d0
    alpha = 0.2d0
    theta = 0.d0
    nflr = 2
    lprojplqu = .false.
    lprojv = .false.
    lproju = .false.
    xmag = 0.d0
    flux_x = 0.d0
    flux_y = 0.d0
    dimer = 0.d0
    rndness = 0.d0

    lstabilize = .true.
    nwrap = 10
    nsweep = 20
    nbin = 10
    nublock = 0
    nustock = 0
    ltau = .false.
    dyntau = 0
    obs_eqt_mid_len = 1

    ! read parameters, read only in irank=0 process
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
    ! broadcast parameters to all processes
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
    call mpi_bcast( rv,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rhub_plq,     1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( alpha,        1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( theta,        1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( nflr,         1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( lprojplqu,    1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( lprojv,       1, mpi_logical,  0, mpi_comm_world, ierr )   
    call mpi_bcast( lproju,       1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( xmag,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( flux_x,       1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( flux_y,       1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( dimer,        1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rndness,      1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( lstabilize,   1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( nwrap,        1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nsweep,       1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nbin,         1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nublock,      1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nustock,      1, mpi_integer,  0, mpi_comm_world, ierr )
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
    if( dble(int(nu))-nu .ne. 0.d0 ) then ! if nu is not integer
        lprojv = .false.
        if( irank == 0 ) then
            write(fout, '(a)') "nu is not integer, lprojv is set to be false!"
        end if
    end if
    if( dble(int(2.d0*nu))-2.d0*nu .ne. 0.d0 ) then ! if 2*nu is not integer
        lproju = .false.
        if( irank == 0 ) then
            write(fout, '(a)') "2*nu is not integer, lproju is set to be false!"
        end if
    end if

    if( rt .ne. 0.d0 ) then
        lwrapT = .true.
    else
        lwrapT = .false.
    end if

    if( rhub .ne. 0.d0 ) then
        lwrapu = .true.
    else
        lwrapu = .false.
    end if

    if( rv .ne. 0.d0 ) then
        lwrapv = .true.
    else
        lwrapv = .false.
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

    if( rv .ne. 0.d0 ) then
        lupdatev = .true.
    else
        lupdatev = .false.
    end if

    if( rhub_plq .ne. 0.d0 ) then
        lupdateplqu = .true.
    else
        lupdateplqu = .false.
    end if

    ! set up lattice
#IFDEF HONEYCOMB
    ! lattice
    lq = la*lb
    ndim = lq*2
   	a1_p(1) = 0.5d0 ; a1_p(2) = -0.5d0*dsqrt(3.d0)
    a2_p(1) = 0.5d0 ; a2_p(2) =  0.5d0*dsqrt(3.d0)
    call latt%setup_honeycomb(la,lb,a1_p,a2_p)
#ENDIF
#IFDEF SQUARE
    ! lattice
    lq = la*lb
    ndim = lq
   	a1_p(1) = 1.0d0 ; a1_p(2) =  0.0d0
    a2_p(1) = 0.0d0 ; a2_p(2) =  1.0d0
    call latt%setup_square(la,lb,a1_p,a2_p)
#ENDIF
#IFDEF CUBIC
    ! lattice
    lq = la*lb*lc
    ndim = lq
   	a1_p(1) = 1.0d0 ; a1_p(2) =  0.0d0; a1_p(3) = 0.d0;
    a2_p(1) = 0.0d0 ; a2_p(2) =  1.0d0; a2_p(3) = 0.d0;
    a3_p(1) = 0.0d0 ; a3_p(2) =  0.0d0; a3_p(3) = 1.d0;
    call latt%setup_cubic(la,lb,lc,a1_p,a2_p,a3_p)
#ENDIF

    ne = ndim/2+int(nu*dble(ndim)/8.d0)  ! number of particles
    ltrot = nint( beta / dtau ) ! total number of time slices

#IFDEF PIFLUX
    xmag = dble(la*lb)*(a1_p(1)*a2_p(2) - a1_p(2)*a2_p(1))/2.d0
#ENDIF
#if defined(DELAY) && defined(SUBMATRIX)
#elif defined(DELAY)
    if( nublock .eq. 0 ) then
        ! tune para for delay update, it is based on some experiences
        if( ndim/5 .lt. 16) then
            nublock = 4
        else if( ndim/5 .lt. 32 ) then
            nublock = 8
        else if( ndim/5 .lt. 64 ) then
            nublock = 16
        else if( ndim/5 .lt. 256 ) then
            nublock = 32
        else ! equal to or greater than 256
            nublock = 64
        end if
    end if
#elif defined(SUBMATRIX)
    if( nublock .eq. 0 ) then
        ! tune para for submatrix update, it is based on some experiences
        if( ndim/20 .lt. 32) then
            nublock = 32
        else if( ndim/20 .lt. 64 ) then
            nublock = 64
        else if( ndim/20 .lt. 128 ) then
            nublock = 16
        else if( ndim/20 .lt. 256 ) then
            nublock = 256
        else if( ndim/20 .lt. 512 ) then
            nublock = 512
        else ! equal to or greater than 1024
            nublock = 1024
        end if
    end if
#endif

    if (nustock .eq. 0) then
        nustock = nublock
    end if

    ! set nst, iwrap_nt and wrap_step
    allocate( iwrap_nt(0:ltrot) )
    iwrap_nt(0:ltrot) = 0
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

    ! set ntdm, ntauin
    if(ltau) then
        ntdm = nint(dyntau/dtau)
        ntauin = (ltrot - ntdm)/2
        xmax_dyn = 0.d0
    end if

    weight_track = 0.d0

    allocate( Imat(ndim,ndim) )
    call s_identity_z(ndim,Imat) ! identity matrix
    allocate( Imat_plq(latt%z_plq,latt%z_plq) )
    call s_identity_z(latt%z_plq,Imat_plq) ! identity matrix with dimension latt%z_plq*latt%z_plq
    allocate( Ivec(ndim) ) ! identity vector
    do i = 1, ndim
        Ivec(i) = 1.d0
    end do
    phase = cone ! inital phase factor is one

  end subroutine make_tables

  subroutine deallocate_tables
    deallocate( wrap_step )
    deallocate( iwrap_nt )
    call latt%deallocate_latt
  end subroutine deallocate_tables

end module model_para
