subroutine dqmc_right_backward_prop(this,gmat,jsite,ntau)
  ! gmat*B(t+dt,t)^-1
  use model_para
  implicit none

  class(plqconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau, jsite

  ! local
  integer :: i, i0, i1, n1, is
  type(gfunc) :: ukmat, ukmat_tmp
  n1 = size(gmat%orb1,1)
  call allocate_gfunc(ukmat, n1, latt%z_plq)
  call allocate_gfunc(ukmat_tmp, n1, latt%z_plq)
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, jsite)
    do i = 1, n1
      ukmat%orb1(i, i0) = gmat%orb1(i, i1)
    end do
  end do
  is = this%conf_plqu(jsite,ntau)
  call zgemm('n','n',n1,latt%z_plq,latt%z_plq,cone,ukmat%orb1,n1,this%bmat_plqu_orb1_inv(:,:,is),latt%z_plq,czero,ukmat_tmp%orb1,n1)
  do i0 = 1, latt%z_plq
    i1 = latt%plq_cord(i0, jsite)
    do i = 1, n1
      gmat%orb1(i, i1) = ukmat_tmp%orb1(i, i0)
    end do
  end do

end subroutine dqmc_right_backward_prop
