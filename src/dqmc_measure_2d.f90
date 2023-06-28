module dqmc_measure
  use model_para
  use dqmc_basic_data
  use dqmc_config_h0
  !use dqmc_config_plqu

  integer, save :: nobs, nobst
  real(dp), save :: sgn
  complex(dp), save :: zphi



  complex(dp), allocatable, dimension(:,:,:), save :: zcpcm_orb1, zcpcm_orb1_bin
  complex(dp), allocatable, dimension(:,:,:), save :: zspsm_orb1, zspsm_orb1_bin
  complex(dp), allocatable, dimension(:,:,:), save :: znn_orb1_bin, znn_orb1
  complex(dp), allocatable, dimension(:), save :: zjj_orb1_bin, zjj_orb1
  complex(dp), allocatable, dimension(:), save :: zbb_orb1_bin, zbb_orb1
  complex(dp), allocatable, dimension(:,:,:), save :: pair_onsite_orb1, pair_onsite_orb1_bin 
  complex(dp), allocatable, dimension(:), save :: pair_nn_orb1, pair_nn_orb1_bin
  complex(dp), allocatable, dimension(:), save :: pair_sn_orb1, pair_sn_orb1_bin

  complex(dp), save :: energy_bin(10), energy_bin_recv(10)
  complex(dp), save :: zn_orb1_bin(2), zn_orb1(2)
  complex(dp), save :: zj_orb1_bin, zj_orb1
  complex(dp), save :: zb_orb1_bin, zb_orb1

  complex(dp), allocatable, dimension(:,:,:,:), save :: zspsm_orb1_tau, zspsm_orb1_tau_bin 
  complex(dp), allocatable, dimension(:,:,:,:), save :: znn_orb1_tau, znn_orb1_tau_bin 
  complex(dp), allocatable, dimension(:,:,:,:), save :: gtau0_orb1_tau, gtau0_orb1_tau_bin
  complex(dp), allocatable, dimension(:,:), save :: zbb_orb1_tau, zbb_orb1_tau_bin
  complex(dp), allocatable, dimension(:,:), save :: zb_orb1_tau, zb_orb1_tau_bin

  contains
#include 'dqmc_measure/equaltime_measure.f90'
#include 'dqmc_measure/equaltime_output.f90'
#include 'dqmc_measure/dyn_measure.f90'
#include 'dqmc_measure/dyn_output.f90'

  subroutine allocate_obs
    implicit none
    allocate( zcpcm_orb1_bin(lq,2,2) )
    allocate( zcpcm_orb1(lq,2,2) )
    allocate( zspsm_orb1_bin(lq,2,2) )
    allocate( zspsm_orb1(lq,2,2) )
    allocate( znn_orb1_bin(lq,2,2) )
    allocate( zjj_orb1_bin(lq) )
    allocate( zbb_orb1_bin(lq) )
    allocate( znn_orb1(lq,2,2) )
    allocate( zjj_orb1(lq) )
    allocate( zbb_orb1(lq) )
    allocate( pair_onsite_orb1_bin(lq,2,2) )
    allocate( pair_onsite_orb1(lq,2,2) )
    allocate( pair_nn_orb1_bin(lq) )
    allocate( pair_sn_orb1_bin(lq) )
    allocate( pair_nn_orb1(lq) )
    allocate( pair_sn_orb1(lq) )
    if(ltau) then
        allocate( zspsm_orb1_tau(lq,2,2,ntdm+1) )
        allocate( znn_orb1_tau(lq,2,2,ntdm+1) )
        allocate( gtau0_orb1_tau(lq,2,2,ntdm+1) )
        allocate( zspsm_orb1_tau_bin(lq,2,2,ntdm+1) )
        allocate( znn_orb1_tau_bin(lq,2,2,ntdm+1) )
        allocate( gtau0_orb1_tau_bin(lq,2,2,ntdm+1) )
        allocate( zbb_orb1_tau(lq,ntdm+1) )
        allocate( zbb_orb1_tau_bin(lq,ntdm+1) )
        allocate( zb_orb1_tau(lq,ntdm+1) )
        allocate( zb_orb1_tau_bin(lq,ntdm+1) )
    end if
  end subroutine allocate_obs

  subroutine deallocate_obs
    implicit none
    if(ltau) then
        deallocate( zb_orb1_tau_bin )
        deallocate( zb_orb1_tau )
        deallocate( zbb_orb1_tau_bin )
        deallocate( zbb_orb1_tau )
        deallocate( gtau0_orb1_tau_bin )
        deallocate( znn_orb1_tau_bin )
        deallocate( zspsm_orb1_tau_bin )
        deallocate( gtau0_orb1_tau )
        deallocate( znn_orb1_tau )
        deallocate( zspsm_orb1_tau )
    end if
    deallocate( zcpcm_orb1_bin )
    deallocate( zcpcm_orb1 )
    deallocate( zspsm_orb1_bin )
    deallocate( zspsm_orb1 )
    deallocate( znn_orb1_bin )
    deallocate( zjj_orb1_bin )
    deallocate( zbb_orb1_bin )
    deallocate( znn_orb1 )
    deallocate( zjj_orb1 )
    deallocate( zbb_orb1 )
    deallocate( pair_onsite_orb1_bin )
    deallocate( pair_onsite_orb1 )
    deallocate( pair_nn_orb1_bin )
    deallocate( pair_sn_orb1_bin )
    deallocate( pair_nn_orb1 )
    deallocate( pair_sn_orb1 )
  end subroutine deallocate_obs

  subroutine obser_init
    implicit none
    nobs = 0
    nobst = 0
    energy_bin = czero
    zspsm_orb1_bin = czero
    zcpcm_orb1_bin = czero
    znn_orb1_bin = czero
    zn_orb1_bin = czero
    zjj_orb1_bin = czero
    zj_orb1_bin = czero
    zbb_orb1_bin = czero
    zb_orb1_bin = czero
    pair_onsite_orb1_bin = czero
    pair_nn_orb1_bin = czero
    pair_sn_orb1_bin = czero
    if(ltau) then
        gtau0_orb1_tau_bin = czero
        znn_orb1_tau_bin = czero
        zspsm_orb1_tau_bin = czero
        zbb_orb1_tau_bin = czero
        zb_orb1_tau_bin = czero
    end if
  end subroutine obser_init

end module dqmc_measure
