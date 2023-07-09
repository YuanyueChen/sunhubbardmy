subroutine dqmc_right_forward_prop(this, gmat, ntau, nf, nflag)
  ! gmat*B(t+dt,t)
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  implicit none

  class(vconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau, nf, nflag

  ! local
  integer :: i, j, i1, i2, is, n1
  complex(dp), allocatable, dimension(:) :: v1, v2
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  n1 = size(gmat%orb1,1)
  allocate( v1(n1), v2(n1) )

  if (nflag.eq.2) then
!$OMP PARALLEL &
!$OMP PRIVATE ( i, i1, i2, j, v1, v2 )
!$OMP DO
     do i = 1,latt%nn_lf
        i1 = latt%nnlf_list(i) 
        i2 = latt%nnlist(i1,nf)
        do j = 1,n1
           v1(j)   =  gmat%orb1(j,i1) * this%u(1,1) + gmat%orb1(j,i2) * this%u(2,1) 
           v2(j)   =  gmat%orb1(j,i1) * this%u(1,2) + gmat%orb1(j,i2) * this%u(2,2) 
        enddo
        do j = 1,n1
           gmat%orb1(j,i1) = v1(j)
           gmat%orb1(j,i2) = v2(j)
        enddo   
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif
   
  if (nflag.eq.1) then
!$OMP PARALLEL &
!$OMP PRIVATE (  i, is, i1, i2, j, v1, v2 )
!$OMP DO
     do i = 1,latt%nn_lf
        i1 = latt%nnlf_list(i)
        i2 = latt%nnlist(i1,nf) 
        is = this%conf_v(i,nf,ntau)
        do j = 1,n1
           gmat%orb1(j,i1) = this%bmat_v_p_orb1(is)*gmat%orb1(j,i1)
           gmat%orb1(j,i2) = this%bmat_v_m_orb1(is)*gmat%orb1(j,i2)
        enddo
        do j = 1,n1
           v1(j) = gmat%orb1(j,i1) * this%ut(1,1) +  gmat%orb1(j,i2) * this%ut(2,1)
           v2(j) = gmat%orb1(j,i1) * this%ut(1,2) +  gmat%orb1(j,i2) * this%ut(2,2)
        enddo
        do j = 1,n1
           gmat%orb1(j,i1) = v1(j)
           gmat%orb1(j,i2) = v2(j)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL     
  endif

#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(11)=timecalculation(11)+endtime-starttime
#ENDIF
  deallocate( v2, v1)
end subroutine dqmc_right_forward_prop
