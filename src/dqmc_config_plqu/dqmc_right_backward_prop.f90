subroutine dqmc_right_backward_prop(this,gmat,jsite,ntau)
  ! gmat*B(t+dt,t)^-1
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  implicit none

  class(plqconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau, jsite

  ! local
  integer :: i, i0, i1, n1, is
  type(gfunc) :: ukmat, ukmat_tmp
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  n1 = size(gmat%blk1,1)
  call allocate_gfunc(ukmat, n1, latt%z_plq)
  call allocate_gfunc(ukmat_tmp, n1, latt%z_plq)
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, jsite)
    do i = 1, n1
      ukmat%blk1(i, i0) = gmat%blk1(i, i1)
    end do
  end do
  is = this%conf_plqu(jsite,ntau)
  call zgemm('n','n',n1,latt%z_plq,latt%z_plq,cone,ukmat%blk1,n1,this%bmat_plqu_blk1_inv(:,:,is),latt%z_plq,czero,ukmat_tmp%blk1,n1)
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, jsite)
    do i = 1, n1
      gmat%blk1(i, i1) = ukmat_tmp%blk1(i, i0)
    end do
  end do
#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(13)=timecalculation(13)+endtime-starttime
#ENDIF

end subroutine dqmc_right_backward_prop
