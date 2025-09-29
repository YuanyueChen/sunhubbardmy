module dqmc_config_h0
  use mkl_dfti
  use constants, only: dp, cone, czero, pi, zero
  use model_para, only: fout, ierr, irank, isize, lprojplqu, ndim

  type h0conf
      integer :: lq, ltrot
      real(dp) :: rt, mu, rth
      complex(dp), allocatable :: h0mat(:,:)
#if defined(BREAKUP_T)
      complex(dp), allocatable :: urt(:,:,:), urtm1(:,:,:)
#elif defined(FFT)
#if defined(SQUARE) || defined(CUBIC) || defined(CHAIN)
      complex(dp), allocatable :: exph0k(:)
      complex(dp), allocatable :: exph0kinv(:)
#elif defined(HONEYCOMB)
      complex(dp), allocatable :: exph0k(:,:,:)
      complex(dp), allocatable :: exph0kinv(:,:,:)
#endif
#else
      complex(dp), allocatable :: urt(:,:), urtm1(:,:)
#endif
      contains
        procedure, public :: set_h0conf           => dqmc_set_h0conf
        procedure, public :: left_backward_prop   => dqmc_left_backward_prop
        procedure, public :: left_forward_prop    => dqmc_left_forward_prop
        procedure, public :: left_forward_prop_hc => dqmc_left_forward_prop_hc
        procedure, public :: right_backward_prop  => dqmc_right_backward_prop
        procedure, public :: right_forward_prop   => dqmc_right_forward_prop
        final :: deallocate_h0conf
  end type h0conf
  
  private :: dqmc_set_h0conf
  private :: dqmc_left_backward_prop
  private :: dqmc_left_forward_prop
  private :: dqmc_left_forward_prop_hc
  private :: dqmc_right_backward_prop
  private :: dqmc_right_forward_prop

  contains

#if defined(CUBIC) || defined(R_CUBIC)
#include 'dqmc_config_h0/expar3d.f90'
#elif defined(CHAIN)
#include 'dqmc_config_h0/expar1d.f90'
#else
#include 'dqmc_config_h0/expar.f90'
#endif
#include 'dqmc_config_h0/seth0.f90'
#include 'dqmc_config_h0/set_trial_h0.f90'
#include 'dqmc_config_h0/dqmc_set_h0conf.f90'
#include 'dqmc_config_h0/dqmc_left_backward_prop.f90'
#include 'dqmc_config_h0/dqmc_left_forward_prop.f90'
#include 'dqmc_config_h0/dqmc_left_forward_prop_hc.f90'
#include 'dqmc_config_h0/dqmc_right_backward_prop.f90'
#include 'dqmc_config_h0/dqmc_right_forward_prop.f90'

  subroutine deallocate_h0conf( this )
    implicit none
    type(h0conf) :: this
    deallocate( this%h0mat )
#if defined(FFT)
    deallocate( this%exph0k )
    deallocate( this%exph0kinv )
#else
    deallocate( this%urt )
    deallocate( this%urtm1 )
#endif
  end subroutine deallocate_h0conf

end module dqmc_config_h0
