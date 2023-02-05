module dqmc_measure
  use model_para
  use dqmc_ft_basic_data
  use dqmc_config_h0
  !use dqmc_config_plqu

  integer, save :: nobs, nobst
  real(dp), save :: sgn
  complex(dp), save :: zphi



  complex(dp), allocatable, dimension(:), save :: zspsm_orb1, zspsm_orb1_bin
  complex(dp), allocatable, dimension(:), save :: znn_orb1_bin, znn_orb1
  complex(dp), allocatable, dimension(:), save :: zbb_orb1_bin, zbb_orb1

  complex(dp), save :: energy_bin(10), energy_bin_recv(10)
  complex(dp), save :: zn_orb1_bin, zn_orb1
  complex(dp), save :: zb_orb1_bin, zb_orb1

  complex(dp), allocatable, dimension(:,:), save :: zspsm_tau, zspsm_tau_bin 
  complex(dp), allocatable, dimension(:,:), save :: znn_tau, znn_tau_bin 
  complex(dp), allocatable, dimension(:,:), save :: gtau0, gtau0_bin

  contains
#IFDEF CUBIC
#include 'dqmc_measure/equaltime_measure3d.f90'
#include 'dqmc_measure/equaltime_output3d.f90'
#include 'dqmc_measure/dyn_measure3d.f90'
#include 'dqmc_measure/dyn_output3d.f90'
#ELSE
#include 'dqmc_measure/equaltime_measure.f90'
#include 'dqmc_measure/equaltime_output.f90'
#include 'dqmc_measure/dyn_measure.f90'
#include 'dqmc_measure/dyn_output.f90'
#ENDIF

  subroutine allocate_obs
    implicit none
    allocate( zspsm_orb1_bin(lq) )
    allocate( zspsm_orb1(lq) )
    if(ltau) then
        allocate( gtau0(lq,ltrot))                     !ltrot = time slices
        allocate( gtau0_bin(lq,ltrot))
        allocate( zspsm_tau(lq,ltrot))                 !ltrot = time slices
        allocate( zspsm_tau_bin(lq,ltrot))
        allocate( znn_tau(lq,ltrot))                 !ltrot = time slices
        allocate( znn_tau_bin(lq,ltrot))
    end if

    allocate( znn_orb1_bin(lq) )
    allocate( zbb_orb1_bin(lq) )
    allocate( znn_orb1(lq) )
    allocate( zbb_orb1(lq) )
  end subroutine allocate_obs

  subroutine deallocate_obs
    implicit none
    deallocate( zspsm_orb1_bin )
    deallocate( zspsm_orb1 )
    if(allocated(znn_tau))      deallocate(znn_tau)
    if(allocated(znn_tau_bin))  deallocate(znn_tau_bin)
    if(allocated(zspsm_tau))      deallocate(zspsm_tau)
    if(allocated(zspsm_tau_bin))  deallocate(zspsm_tau_bin)
    if(allocated(gtau0))      deallocate(gtau0)
    if(allocated(gtau0_bin))  deallocate(gtau0_bin)
    deallocate( znn_orb1_bin )
    deallocate( zbb_orb1_bin )
    deallocate( znn_orb1 )
    deallocate( zbb_orb1 )
  end subroutine deallocate_obs

  subroutine obser_init
    implicit none
    nobs = 0
    nobst = 0
    energy_bin = czero
    zspsm_orb1_bin = czero
    znn_orb1_bin = czero
    zn_orb1_bin = czero
    zbb_orb1_bin = czero
    zb_orb1_bin = czero
    if(ltau) then
        gtau0_bin = czero
        zspsm_tau_bin = czero
        znn_tau_bin = czero
    end if
  end subroutine obser_init

end module dqmc_measure
