#ifdef CUBIC
#include 'dqmc_measure_3d.f90'
#elif CHAIN
#include 'dqmc_measure_1d.f90'
#else
#include 'dqmc_measure_2d.f90'
#endif
