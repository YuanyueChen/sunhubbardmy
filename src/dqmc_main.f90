program dqmc_main
  use constants
  use dqmc_util
  use dqmc_ctrl
  implicit none

  include 'mpif.h'

  call MPI_INIT(ierr)                             
  call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr) 

  if(irank==0) open( unit=fout, file='dqmc.out', status='unknown' )

  call dqmc_initial
  call dqmc_start
  call dqmc_warmup

#IFDEF _OPENMP
  time1 = omp_get_wtime()
#ELSE
  call cpu_time(time1)
#ENDIF

  do nbc =  1, nbin
      call dqmc_sweep
      call dqmc_timing
  end do
  
  call dqmc_core_deallocate
  call dqmc_end

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

end program dqmc_main
