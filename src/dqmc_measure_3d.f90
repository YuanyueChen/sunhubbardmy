module dqmc_measure
  use model_para
  use dqmc_basic_data
  use dqmc_config_h0
  !use dqmc_config_plqu

  integer, save :: nobs, nobst
  real(dp), save :: sgn
  complex(dp), save :: zphi



  complex(dp), allocatable, dimension(:), save :: zcpcm, zcpcm_bin
  complex(dp), allocatable, dimension(:), save :: zspsm, zspsm_bin
  complex(dp), allocatable, dimension(:), save :: znn_bin, znn
  complex(dp), allocatable, dimension(:), save :: zbb_bin, zbb

  complex(dp), save :: energy_bin(10), energy_bin_recv(10)
  complex(dp), save :: zn_bin, zn
  complex(dp), save :: zb_bin, zb

  complex(dp), allocatable, dimension(:,:), save :: zspsm_tau, zspsm_tau_bin 
  complex(dp), allocatable, dimension(:,:), save :: znn_tau, znn_tau_bin 
  complex(dp), allocatable, dimension(:,:), save :: gtau0, gtau0_bin

  contains
#include 'dqmc_measure/equaltime_measure3d.f90'
#include 'dqmc_measure/equaltime_output3d.f90'
#include 'dqmc_measure/dyn_measure3d.f90'
#include 'dqmc_measure/dyn_output3d.f90'

  subroutine allocate_obs
    implicit none
    allocate( zcpcm_bin(lq) )
    allocate( zcpcm(lq) )
    allocate( zspsm_bin(lq) )
    allocate( zspsm(lq) )
    if(ltau) then
        allocate( gtau0(lq,ntdm+1))
        allocate( gtau0_bin(lq,ntdm+1))
        allocate( zspsm_tau(lq,ntdm+1))
        allocate( zspsm_tau_bin(lq,ntdm+1))
        allocate( znn_tau(lq,ntdm+1))
        allocate( znn_tau_bin(lq,ntdm+1))
    end if

    allocate( znn_bin(lq) )
    allocate( zbb_bin(lq) )
    allocate( znn(lq) )
    allocate( zbb(lq) )
  end subroutine allocate_obs

  subroutine deallocate_obs
    implicit none
    deallocate( zcpcm_bin )
    deallocate( zcpcm )
    deallocate( zspsm_bin )
    deallocate( zspsm )
    if(allocated(znn_tau))      deallocate(znn_tau)
    if(allocated(znn_tau_bin))  deallocate(znn_tau_bin)
    if(allocated(zspsm_tau))      deallocate(zspsm_tau)
    if(allocated(zspsm_tau_bin))  deallocate(zspsm_tau_bin)
    if(allocated(gtau0))      deallocate(gtau0)
    if(allocated(gtau0_bin))  deallocate(gtau0_bin)
    deallocate( znn_bin )
    deallocate( zbb_bin )
    deallocate( znn )
    deallocate( zbb )
  end subroutine deallocate_obs

  subroutine obser_init
    implicit none
    nobs = 0
    nobst = 0
    energy_bin = czero
    zspsm_bin = czero
    zcpcm_bin = czero
    znn_bin = czero
    zn_bin = czero
    zbb_bin = czero
    zb_bin = czero
    if(ltau) then
        gtau0_bin = czero
        zspsm_tau_bin = czero
        znn_tau_bin = czero
    end if
  end subroutine obser_init

end module dqmc_measure
