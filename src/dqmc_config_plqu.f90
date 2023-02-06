module dqmc_config_plqu
  use constants, only: dp, cone, czero, pi
  use model_para, only: fout, ierr, irank, isize, lprojplqu, nflr
  use dqmc_basic_data
  type plqconf
      logical :: lwarmup
      integer :: lq, ltrot, lcomp
      real(dp) :: alpha_plqu
      integer, allocatable:: conf_plqu(:,:)
      complex(dp), allocatable :: bmat_plqu_orb1(:,:,:)
      complex(dp), allocatable :: bmat_plqu_orb1_inv(:,:,:)
      complex(dp), allocatable :: delta_bmat_plqu_orb1(:,:,:,:)
      complex(dp), allocatable :: phase(:)
      complex(dp), allocatable :: phase_ratio(:,:)
      contains
        procedure, public :: set_plqconf           => dqmc_set_plqconf
        procedure, public :: input_plqconf         => dqmc_input_plqconf
        procedure, public :: output_plqconf        => dqmc_output_plqconf
        procedure, public :: left_backward_prop    => dqmc_left_backward_prop
        procedure, public :: left_forward_prop     => dqmc_left_forward_prop
        procedure, public :: left_forward_prop_hc  => dqmc_left_forward_prop_hc
        procedure, public :: right_backward_prop   => dqmc_right_backward_prop
        procedure, public :: right_forward_prop    => dqmc_right_forward_prop
        procedure, public :: update_plqu           => dqmc_update_plqu
        procedure, public :: proj_update_plqu      => dqmc_proj_update_plqu
        final :: deallocate_plqconf
  end type plqconf
  private :: dqmc_set_plqconf
  private :: dqmc_input_plqconf
  private :: dqmc_output_plqconf
  private :: dqmc_left_backward_prop
  private :: dqmc_left_forward_prop
  private :: dqmc_left_forward_prop_hc
  private :: dqmc_right_backward_prop
  private :: dqmc_right_forward_prop
  private :: dqmc_update_plqu
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
#include 'dqmc_config_plqu/dqmc_set_plqconf.f90'
#include 'dqmc_config_plqu/dqmc_input_plqconf.f90'
#include 'dqmc_config_plqu/dqmc_output_plqconf.f90'
#include 'dqmc_config_plqu/dqmc_left_backward_prop.f90'
#include 'dqmc_config_plqu/dqmc_left_forward_prop.f90'
#include 'dqmc_config_plqu/dqmc_left_forward_prop_hc.f90'
#include 'dqmc_config_plqu/dqmc_right_backward_prop.f90'
#include 'dqmc_config_plqu/dqmc_right_forward_prop.f90'
#include 'dqmc_config_plqu/dqmc_update_plqu.f90'
#include 'dqmc_config_plqu/dqmc_proj_update_plqu.f90'

  subroutine deallocate_plqconf( this )
    implicit none
    type(plqconf) :: this
    deallocate( this%conf_plqu )
    deallocate( this%bmat_plqu_orb1 )
    deallocate( this%bmat_plqu_orb1_inv )
    deallocate( this%delta_bmat_plqu_orb1 )
    deallocate( this%phase )
    deallocate( this%phase_ratio )
  end subroutine deallocate_plqconf


end module dqmc_config_plqu
