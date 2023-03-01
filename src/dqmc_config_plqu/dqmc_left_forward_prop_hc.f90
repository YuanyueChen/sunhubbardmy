subroutine dqmc_left_forward_prop_hc(this, gmat, isite, ntau)
  ! *gmat
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
  real(dp) :: starttime13, endtime13
#ENDIF
#IFDEF TIMING
  starttime13 = omp_get_wtime()
#ENDIF

  n2 = size(gmat%orb1,2)
  call allocate_gfunc( vkmat, latt%z_plq, n2 )
  call allocate_gfunc( vkmat_tmp, latt%z_plq, n2 )
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, isite)
    do i = 1, n2
      vkmat%orb1(i0, i) = gmat%orb1(i1, i)
    end do
  end do
  is = this%conf_plqu(isite,ntau)
  ! exp(iV)^H = exp(-iV)
  call zgemm('n','n',latt%z_plq,n2,latt%z_plq,cone,this%bmat_plqu_orb1_inv(:,:,is),latt%z_plq,vkmat%orb1,latt%z_plq,czero,vkmat_tmp%orb1,latt%z_plq)
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, isite)
    do i = 1, n2
      gmat%orb1(i1, i) = vkmat_tmp%orb1(i0, i)
    end do
  end do
#IFDEF TIMING
  endtime13 = omp_get_wtime()
  timecalculation(12)=timecalculation(12)+endtime13-starttime13
#ENDIF 

end subroutine dqmc_left_forward_prop_hc
