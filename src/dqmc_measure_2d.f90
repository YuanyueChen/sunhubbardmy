module dqmc_measure
  use model_para
  use dqmc_basic_data
  use dqmc_config_h0
  !use dqmc_config_plqu

  integer, save :: nobs, nobst
  real(dp), save :: sgn
  complex(dp), save :: zphi



  complex(dp), allocatable, dimension(:,:,:), save :: zcpcm, zcpcm_bin
  complex(dp), allocatable, dimension(:,:,:), save :: zspsm, zspsm_bin
  complex(dp), allocatable, dimension(:,:,:), save :: znn_bin, znn
  complex(dp), allocatable, dimension(:), save :: zjj_bin, zjj
  complex(dp), allocatable, dimension(:), save :: zbb_bin, zbb
  complex(dp), allocatable, dimension(:,:,:), save :: pair_onsite, pair_onsite_bin 
  complex(dp), allocatable, dimension(:), save :: pair_nn, pair_nn_bin
  complex(dp), allocatable, dimension(:), save :: pair_sn, pair_sn_bin

  complex(dp), save :: energy_bin(10), energy_bin_recv(10)
  complex(dp), save :: zn_bin(2), zn(2)
  complex(dp), save :: zj_bin, zj
  complex(dp), save :: zb_bin, zb

  complex(dp), allocatable, dimension(:,:,:,:), save :: zspsm_tau, zspsm_tau_bin 
  complex(dp), allocatable, dimension(:,:,:,:), save :: znn_tau, znn_tau_bin 
  complex(dp), allocatable, dimension(:,:,:,:), save :: gtau0_tau, gtau0_tau_bin
  complex(dp), allocatable, dimension(:,:), save :: zbb_tau, zbb_tau_bin
  complex(dp), allocatable, dimension(:,:), save :: zb_tau, zb_tau_bin

  contains
#include 'dqmc_measure/equaltime_measure.f90'
#include 'dqmc_measure/equaltime_output.f90'
#include 'dqmc_measure/dyn_measure.f90'
#include 'dqmc_measure/dyn_output.f90'

  subroutine allocate_obs
    implicit none
    allocate( zcpcm_bin(lq,2,2) )
    allocate( zcpcm(lq,2,2) )
    allocate( zspsm_bin(lq,2,2) )
    allocate( zspsm(lq,2,2) )
    allocate( znn_bin(lq,2,2) )
    allocate( zjj_bin(lq) )
    allocate( zbb_bin(lq) )
    allocate( znn(lq,2,2) )
    allocate( zjj(lq) )
    allocate( zbb(lq) )
    allocate( pair_onsite_bin(lq,2,2) )
    allocate( pair_onsite(lq,2,2) )
    allocate( pair_nn_bin(lq) )
    allocate( pair_sn_bin(lq) )
    allocate( pair_nn(lq) )
    allocate( pair_sn(lq) )
    if(ltau) then
        allocate( zspsm_tau(lq,2,2,ntdm+1) )
        allocate( znn_tau(lq,2,2,ntdm+1) )
        allocate( gtau0_tau(lq,2,2,ntdm+1) )
        allocate( zspsm_tau_bin(lq,2,2,ntdm+1) )
        allocate( znn_tau_bin(lq,2,2,ntdm+1) )
        allocate( gtau0_tau_bin(lq,2,2,ntdm+1) )
        allocate( zbb_tau(lq,ntdm+1) )
        allocate( zbb_tau_bin(lq,ntdm+1) )
        allocate( zb_tau(lq,ntdm+1) )
        allocate( zb_tau_bin(lq,ntdm+1) )
    end if
  end subroutine allocate_obs

  subroutine deallocate_obs
    implicit none
    if(ltau) then
        deallocate( zb_tau_bin )
        deallocate( zb_tau )
        deallocate( zbb_tau_bin )
        deallocate( zbb_tau )
        deallocate( gtau0_tau_bin )
        deallocate( znn_tau_bin )
        deallocate( zspsm_tau_bin )
        deallocate( gtau0_tau )
        deallocate( znn_tau )
        deallocate( zspsm_tau )
    end if
    deallocate( zcpcm_bin )
    deallocate( zcpcm )
    deallocate( zspsm_bin )
    deallocate( zspsm )
    deallocate( znn_bin )
    deallocate( zjj_bin )
    deallocate( zbb_bin )
    deallocate( znn )
    deallocate( zjj )
    deallocate( zbb )
    deallocate( pair_onsite_bin )
    deallocate( pair_onsite )
    deallocate( pair_nn_bin )
    deallocate( pair_sn_bin )
    deallocate( pair_nn )
    deallocate( pair_sn )
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
    zjj_bin = czero
    zj_bin = czero
    zbb_bin = czero
    zb_bin = czero
    pair_onsite_bin = czero
    pair_nn_bin = czero
    pair_sn_bin = czero
    if(ltau) then
        gtau0_tau_bin = czero
        znn_tau_bin = czero
        zspsm_tau_bin = czero
        zbb_tau_bin = czero
        zb_tau_bin = czero
    end if
  end subroutine obser_init

end module dqmc_measure
