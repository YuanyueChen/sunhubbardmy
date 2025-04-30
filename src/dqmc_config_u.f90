module dqmc_config_u
  use constants, only: dp, cone, czero, pi
  use model_para, only: fout, ierr, irank, isize, lproju, nflr, latt
  use dqmc_basic_data
  type uconf
  !! gaml(l) * exp( (i) * alpha * etal(l) * n ) * phase(l)
      logical :: lwarmup
      integer :: lq, nsites, ltrot, lcomp
      real(dp) :: alpha
      real(dp), allocatable :: etal(:), gaml(:)
      integer, allocatable :: conf(:,:)
      type(zvfunc) :: bmat
      type(zvfunc) :: bmat_inv
      type(gfunc) :: delta_bmat
      complex(dp), allocatable :: phase(:)
      complex(dp), allocatable :: phase_ratio(:,:)
      contains
        procedure, public :: set_conf           => dqmc_set_conf
        procedure, public :: input_conf         => dqmc_input_conf                !input configuration
        procedure, public :: output_conf        => dqmc_output_conf               !output configuration
        procedure, public :: left_backward_prop    => dqmc_left_backward_prop
        procedure, public :: left_forward_prop     => dqmc_left_forward_prop
        procedure, public :: left_forward_prop_hc  => dqmc_left_forward_prop_hc
        procedure, public :: right_backward_prop   => dqmc_right_backward_prop
        procedure, public :: right_forward_prop    => dqmc_right_forward_prop
        procedure, public :: update           => dqmc_update                !finite temperature
        procedure, public :: proj_update      => dqmc_proj_update           !projector
        final :: deallocate_conf
  end type uconf
  private :: dqmc_set_conf
  private :: dqmc_input_conf
  private :: dqmc_output_conf
  private :: dqmc_left_backward_prop
  private :: dqmc_left_forward_prop
  private :: dqmc_left_forward_prop_hc
  private :: dqmc_right_backward_prop
  private :: dqmc_right_forward_prop
  private :: dqmc_update

  contains

#include 'dqmc_config_u/dqmc_set_conf.f90'
#include 'dqmc_config_u/dqmc_input_conf.f90'
#include 'dqmc_config_u/dqmc_output_conf.f90'
#include 'dqmc_config_u/dqmc_left_backward_prop.f90'
#include 'dqmc_config_u/dqmc_left_forward_prop.f90'
#include 'dqmc_config_u/dqmc_left_forward_prop_hc.f90'
#include 'dqmc_config_u/dqmc_right_backward_prop.f90'
#include 'dqmc_config_u/dqmc_right_forward_prop.f90'
#if defined(DELAY) && defined(SUBMATRIX)
#elif defined(DELAY)
#include 'dqmc_config_u/dqmc_update_delay.f90'
#include 'dqmc_config_u/dqmc_proj_update_delay.f90'
#elif defined(SUBMATRIX)
#include 'dqmc_config_u/dqmc_update_submatrix.f90'
#if defined(DELAYG)
#include 'dqmc_config_u/dqmc_proj_update_submatrix_G.f90'
#else
#include 'dqmc_config_u/dqmc_proj_update_submatrix_LR.f90'
#endif
#else
#include 'dqmc_config_u/dqmc_update.f90'
#include 'dqmc_config_u/dqmc_proj_update.f90'
#endif

  subroutine deallocate_conf( this )
    implicit none
    type(uconf) :: this
    deallocate( this%conf )
    deallocate( this%etal )
    deallocate( this%gaml )
    deallocate( this%phase )
    deallocate( this%phase_ratio )
  end subroutine deallocate_conf


end module dqmc_config_u
