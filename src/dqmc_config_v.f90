module dqmc_config_v
  use constants, only: dp, cone, czero, pi
  use model_para, only: fout, ierr, irank, isize, lprojv, nflr, latt
  use dqmc_basic_data
  type vconf
      logical :: lwarmup
      integer :: lq, nsites, ltrot, lcomp, nn_nf
      real(dp) :: alpha
      real(dp) :: lambda
      integer, allocatable:: conf(:,:,:)
      type(zvfunc) :: bmat_p
      type(zvfunc) :: bmat_p_inv
      type(zvfunc) :: bmat_m
      type(zvfunc) :: bmat_m_inv
      type(gfunc) :: delta_bmat_p
      type(gfunc) :: delta_bmat_m
      complex(dp), allocatable :: phase(:)
      complex(dp), allocatable :: phase_ratio(:,:)
      complex(dp) :: u(2,2) = (/ 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), -1.d0/dsqrt(2.d0)/)
      complex(dp) :: ut(2,2) = (/ 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), -1.d0/dsqrt(2.d0) /)
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
        procedure, public :: proj_update     => dqmc_proj_update           !projector
        final :: deallocate_conf
  end type vconf
  private :: dqmc_set_conf
  private :: dqmc_input_conf
  private :: dqmc_output_conf
  private :: dqmc_left_backward_prop
  private :: dqmc_left_forward_prop
  private :: dqmc_left_forward_prop_hc
  private :: dqmc_right_backward_prop
  private :: dqmc_right_forward_prop
  private :: dqmc_update
  ! for DQMC
  ! here -1,1 is mapped to 1,2
  real(dp), public, parameter :: etal(1:2) = (/-1,1/)
  integer, public, parameter :: nflipl(1:2,1:2)  = reshape( (/ 1, 2, & 
                                                               2, 1  /), (/2,2/) )
  contains
#include 'dqmc_config_v/dqmc_set_conf.f90'
#include 'dqmc_config_v/dqmc_input_conf.f90'
#include 'dqmc_config_v/dqmc_output_conf.f90'
#include 'dqmc_config_v/dqmc_left_backward_prop.f90'
#include 'dqmc_config_v/dqmc_left_forward_prop.f90'
#include 'dqmc_config_v/dqmc_left_forward_prop_hc.f90'
#include 'dqmc_config_v/dqmc_right_backward_prop.f90'
#include 'dqmc_config_v/dqmc_right_forward_prop.f90'
#if defined(DELAY) && defined(SUBMATRIX)
#elif defined(DELAY)
#include 'dqmc_config_v/dqmc_update_delay.f90'
#include 'dqmc_config_v/dqmc_proj_update_delay.f90'
#elif defined(SUBMATRIX1)
#include 'dqmc_config_v/dqmc_update_submatrix.f90'
#include 'dqmc_config_v/dqmc_proj_update.f90'
#else
#include 'dqmc_config_v/dqmc_update.f90'
#include 'dqmc_config_v/dqmc_proj_update.f90'
#endif

  subroutine deallocate_conf( this )
    implicit none
    type(vconf) :: this
    deallocate( this%conf )
    deallocate( this%phase )
    deallocate( this%phase_ratio )
  end subroutine deallocate_conf


end module dqmc_config_v
