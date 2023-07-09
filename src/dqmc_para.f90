module dqmc_para
  use dqmc_config_h0
  use dqmc_config_plqu
  use dqmc_config_u
  use dqmc_config_v
  implicit none

  ! conf.
  type(h0conf), save :: h0c
  type(plqconf), save :: hconf
  type(uconf), save :: u0conf
  type(vconf), save :: v0conf

end module dqmc_para
