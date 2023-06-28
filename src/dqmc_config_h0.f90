module dqmc_config_h0
  use mkl_dfti
  use constants, only: dp, cone, czero, pi, zero
  use model_para, only: fout, ierr, irank, isize, lprojplqu, ndim

  type h0conf
      integer :: lq, ltrot
      real(dp) :: rt, mu
      complex(dp), allocatable :: h0mat(:,:)
#IFDEF BREAKUP_T
      complex(dp), allocatable :: urt(:,:,:), urtm1(:,:,:)

#ELIF FFT

#IFDEF SQUARE
      complex(dp), allocatable :: exph0k(:)
      complex(dp), allocatable :: exph0kinv(:)
#ELIF CUBIC
      complex(dp), allocatable :: exph0k(:)
      complex(dp), allocatable :: exph0kinv(:)
#ELIF HONEYCOMB
      complex(dp), allocatable :: exph0k(:,:,:)
      complex(dp), allocatable :: exph0kinv(:,:,:)
#ENDIF

#ELSE
      complex(dp), allocatable :: urt(:,:), urtm1(:,:)
#ENDIF
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

#IFDEF CUBIC
#include 'dqmc_config_h0/expar3d.f90'
#ELSE
#include 'dqmc_config_h0/expar.f90'
#ENDIF
#include 'dqmc_config_h0/seth0.f90'
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
#IFDEF FFT
    deallocate( this%exph0k )
    deallocate( this%exph0kinv )
#ELSE
    deallocate( this%urt )
    deallocate( this%urtm1 )
#ENDIF
  end subroutine deallocate_h0conf

end module dqmc_config_h0
