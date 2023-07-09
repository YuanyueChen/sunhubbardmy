module dqmc_config_v
  use constants, only: dp, cone, czero, pi
  use model_para, only: fout, ierr, irank, isize, lprojv, nflr, latt
  use dqmc_basic_data
  type vconf
      logical :: lwarmup
      integer :: lq, nsites, ltrot, lcomp, nn_nf
      real(dp) :: alpha_v
      real(dp) :: lambda
      integer, allocatable:: conf_v(:,:,:)
      complex(dp), allocatable :: bmat_v_p_orb1(:)
      complex(dp), allocatable :: bmat_v_p_orb1_inv(:)
      complex(dp), allocatable :: bmat_v_m_orb1(:)
      complex(dp), allocatable :: bmat_v_m_orb1_inv(:)
      complex(dp), allocatable :: delta_bmat_v_p_orb1(:,:)
      complex(dp), allocatable :: delta_bmat_v_m_orb1(:,:)
      complex(dp), allocatable :: phase(:)
      complex(dp), allocatable :: phase_ratio(:,:)
      complex(dp) :: u(2,2) = (/ 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), -1.d0/dsqrt(2.d0)/)
      complex(dp) :: ut(2,2) = (/ 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), 1.d0/dsqrt(2.d0), -1.d0/dsqrt(2.d0) /)
      contains
        procedure, public :: set_vconf           => dqmc_set_vconf
        procedure, public :: input_vconf         => dqmc_input_vconf                !input configuration
        procedure, public :: output_vconf        => dqmc_output_vconf               !output configuration
        procedure, public :: left_backward_prop    => dqmc_left_backward_prop
        procedure, public :: left_forward_prop     => dqmc_left_forward_prop
        procedure, public :: left_forward_prop_hc  => dqmc_left_forward_prop_hc
        procedure, public :: right_backward_prop   => dqmc_right_backward_prop
        procedure, public :: right_forward_prop    => dqmc_right_forward_prop
        procedure, public :: update_v           => dqmc_update_v                !finite temperature
        procedure, public :: proj_update_v     => dqmc_proj_update_v           !projector
        final :: deallocate_vconf
  end type vconf
  private :: dqmc_set_vconf
  private :: dqmc_input_vconf
  private :: dqmc_output_vconf
  private :: dqmc_left_backward_prop
  private :: dqmc_left_forward_prop
  private :: dqmc_left_forward_prop_hc
  private :: dqmc_right_backward_prop
  private :: dqmc_right_forward_prop
  private :: dqmc_update_v
  ! for DQMC
  ! here -1,1 is mapped to 1,2
  real(dp), public, parameter :: etal(1:2) = (/-1,1/)
  integer, public, parameter :: nflipl(1:2,1:2)  = reshape( (/ 1, 2, & 
                                                               2, 1  /), (/2,2/) )
  contains
#include 'dqmc_config_v/dqmc_set_vconf.f90'
#include 'dqmc_config_v/dqmc_input_vconf.f90'
#include 'dqmc_config_v/dqmc_output_vconf.f90'
#include 'dqmc_config_v/dqmc_left_backward_prop.f90'
#include 'dqmc_config_v/dqmc_left_forward_prop.f90'
#include 'dqmc_config_v/dqmc_left_forward_prop_hc.f90'
#include 'dqmc_config_v/dqmc_right_backward_prop.f90'
#include 'dqmc_config_v/dqmc_right_forward_prop.f90'
#IFDEF DELAY
#include 'dqmc_config_v/dqmc_update_v_delay.f90'
#ELSE
#include 'dqmc_config_v/dqmc_update_v.f90'
#ENDIF
#IFDEF DELAY
#include 'dqmc_config_v/dqmc_proj_update_v_delay.f90'
#ELSE
#include 'dqmc_config_v/dqmc_proj_update_v.f90'
#ENDIF
  subroutine deallocate_vconf( this )
    implicit none
    type(vconf) :: this
    deallocate( this%conf_v )
    deallocate( this%bmat_v_p_orb1 )
    deallocate( this%bmat_v_p_orb1_inv )
    deallocate( this%bmat_v_m_orb1 )
    deallocate( this%bmat_v_m_orb1_inv )
    deallocate( this%delta_bmat_v_p_orb1)
    deallocate( this%delta_bmat_v_m_orb1)
    deallocate( this%phase )
    deallocate( this%phase_ratio )
  end subroutine deallocate_vconf


end module dqmc_config_v
