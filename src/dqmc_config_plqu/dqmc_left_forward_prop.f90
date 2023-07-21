subroutine dqmc_left_forward_prop(this, gmat, isite, ntau)
  ! B(t+dt,t)*gmat
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  implicit none

  class(plqconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau, isite

  ! local
  integer :: i, i0, i1, n2, is
  type(gfunc) :: vkmat, vkmat_tmp
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  n2 = size(gmat%blk1,2)
  call allocate_gfunc( vkmat, latt%z_plq, n2 )
  call allocate_gfunc( vkmat_tmp, latt%z_plq, n2 )
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, isite)
    do i = 1, n2
      vkmat%blk1(i0, i) = gmat%blk1(i1, i)
    end do
  end do
  is = this%conf_plqu(isite,ntau)
  call zgemm('n','n',latt%z_plq,n2,latt%z_plq,cone,this%bmat_plqu_blk1(:,:,is),latt%z_plq,vkmat%blk1,latt%z_plq,czero,vkmat_tmp%blk1,latt%z_plq)
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, isite)
    do i = 1, n2
      gmat%blk1(i1, i) = vkmat_tmp%blk1(i0, i)
    end do
  end do
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(12)=timecalculation(12)+endtime-starttime
#ENDIF

end subroutine dqmc_left_forward_prop
