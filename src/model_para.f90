module model_para
  use constants, only: dp, czero, cone, ctwo, cthree, cfour
#IFDEF HONEYCOMB
  use honeycomb_lattice
#ENDIF
#IFDEF SQUARE
  use square_lattice
#ENDIF
#IFDEF CUBIC
  use cubic_lattice
#ENDIF
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

end module model_para
