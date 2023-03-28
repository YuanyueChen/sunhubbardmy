module dqmc_config_u
  use constants, only: dp, cone, czero, pi
  use model_para, only: fout, ierr, irank, isize, lproju, nflr, latt
  use dqmc_basic_data
  type uconf
      logical :: lwarmup
      integer :: lq, nsites, ltrot, lcomp
      real(dp) :: alpha_u
      integer, allocatable:: conf_u(:,:)
      complex(dp), allocatable :: bmat_u_orb1(:)
      complex(dp), allocatable :: bmat_u_orb1_inv(:)
      complex(dp), allocatable :: delta_bmat_u_orb1(:,:)
      complex(dp), allocatable :: phase(:)
      complex(dp), allocatable :: phase_ratio(:,:)
      contains
        procedure, public :: set_uconf           => dqmc_set_uconf
        procedure, public :: input_uconf         => dqmc_input_uconf                !input configuration
        procedure, public :: output_uconf        => dqmc_output_uconf               !output configuration
        procedure, public :: left_backward_prop    => dqmc_left_backward_prop
        procedure, public :: left_forward_prop     => dqmc_left_forward_prop
        procedure, public :: left_forward_prop_hc  => dqmc_left_forward_prop_hc
        procedure, public :: right_backward_prop   => dqmc_right_backward_prop
        procedure, public :: right_forward_prop    => dqmc_right_forward_prop
        procedure, public :: update_u           => dqmc_update_u                !finite temperature
        procedure, public :: proj_update_u      => dqmc_proj_update_u           !projector
        final :: deallocate_uconf
  end type uconf
  private :: dqmc_set_uconf
  private :: dqmc_input_uconf
  private :: dqmc_output_uconf
  private :: dqmc_left_backward_prop
  private :: dqmc_left_forward_prop
  private :: dqmc_left_forward_prop_hc
  private :: dqmc_right_backward_prop
  private :: dqmc_right_forward_prop
  private :: dqmc_update_u
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
#include 'dqmc_config_u/dqmc_set_uconf.f90'
#include 'dqmc_config_u/dqmc_input_uconf.f90'
#include 'dqmc_config_u/dqmc_output_uconf.f90'
#include 'dqmc_config_u/dqmc_left_backward_prop.f90'
#include 'dqmc_config_u/dqmc_left_forward_prop.f90'
#include 'dqmc_config_u/dqmc_left_forward_prop_hc.f90'
#include 'dqmc_config_u/dqmc_right_backward_prop.f90'
#include 'dqmc_config_u/dqmc_right_forward_prop.f90'
#IFDEF DELAY
#include 'dqmc_config_u/dqmc_update_u_delay.f90'
#ELSE
#include 'dqmc_config_u/dqmc_update_u.f90'
#ENDIF
#IFDEF DELAY
#include 'dqmc_config_u/dqmc_proj_update_u_delay.f90'
#ELSE
#include 'dqmc_config_u/dqmc_proj_update_u.f90'
#ENDIF
  subroutine deallocate_uconf( this )
    implicit none
    type(uconf) :: this
    deallocate( this%conf_u )
    deallocate( this%bmat_u_orb1 )
    deallocate( this%bmat_u_orb1_inv )
    deallocate( this%delta_bmat_u_orb1 )
    deallocate( this%phase )
    deallocate( this%phase_ratio )
  end subroutine deallocate_uconf


end module dqmc_config_u
