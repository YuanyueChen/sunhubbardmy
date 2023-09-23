module dqmc_config_u
  use constants, only: dp, cone, czero, pi
  use model_para, only: fout, ierr, irank, isize, lproju, nflr, latt
  use dqmc_basic_data
  type uconf
      logical :: lwarmup
      integer :: lq, nsites, ltrot, lcomp
      real(dp) :: alpha
      integer, allocatable:: conf(:,:)
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
  ! for DQMC
  ! here -2,-1,1,2 is mapped to 1,2,3,4
  real(dp), public, parameter :: etal(1:4) = (/  - dsqrt(2.d0 * ( 3.d0 + dsqrt(6.d0) ) ), &
                                                 - dsqrt(2.d0 * ( 3.d0 - dsqrt(6.d0) ) ), &
                                                   dsqrt(2.d0 * ( 3.d0 - dsqrt(6.d0) ) ), &
                                                   dsqrt(2.d0 * ( 3.d0 + dsqrt(6.d0) ) ) /)
  real(dp), public, parameter :: gaml(1:4) =  (/ 1.d0 - dsqrt(6.d0)/3.d0, &
                                                 1.d0 + dsqrt(6.d0)/3.d0, &
                                                 1.d0 + dsqrt(6.d0)/3.d0, &
                                                 1.d0 - dsqrt(6.d0)/3.d0 /)
  integer, public, parameter :: nflipl(1:4,1:3)  = reshape( (/ 2, 3, 4, & 
                                                               3, 4, 1, &
                                                               4, 1, 2, &
                                                               1, 2, 3 /), (/4,3/) )
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
  ! 如果同时定义了 DELAY 和 SUBMATRIX，则不包含任何文件
#elif defined(DELAY)
  ! 如果只定义了 DELAY，根据需要包含相关文件
#include 'dqmc_config_u/dqmc_update_delay.f90'
#include 'dqmc_config_u/dqmc_proj_update_delay.f90'
#elif defined(SUBMATRIX)
  ! 如果只定义了 SUBMATRIX，根据需要包含相关文件
#include 'dqmc_config_u/dqmc_update_submatrix.f90'
#include 'dqmc_config_u/dqmc_proj_update.f90'
#else
 !  如果没有定义 DELAY 和 SUBMATRIX，包含默认文件
#include 'dqmc_config_u/dqmc_update.f90'
#include 'dqmc_config_u/dqmc_proj_update.f90'
#endif

  subroutine deallocate_conf( this )
    implicit none
    type(uconf) :: this
    deallocate( this%conf )
    deallocate( this%phase )
    deallocate( this%phase_ratio )
  end subroutine deallocate_conf


end module dqmc_config_u
