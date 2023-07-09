subroutine dqmc_left_backward_prop(this, gmat, ntau, nf, nflag)
  ! B(t+dt,t)^-1 *gmat
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use model_para
  implicit none

  class(vconf) :: this
  type(gfunc) :: gmat
  integer, intent(in) :: ntau, nf, nflag

  ! local
  integer :: i, j, i1, i2, n2, is, isite
  complex(dp), allocatable, dimension(:) :: v1, v2
#IFDEF TIMING
  real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
  call cpu_time_now(starttime)
#ENDIF

  n2 = size(gmat%orb1,2)
  allocate( v1(n2), v2(n2) )

  if (nflag.eq.2) then
!$OMP PARALLEL &
!$OMP PRIVATE ( i, i1, i2, j, v1, v2 )
!$OMP DO
       do i = 1,latt%nn_lf
        i1 = latt%nnlf_list(i)
        i2 = latt%nnlist(i1,nf)
        do j = 1,n2
           v1(j)   =  this%ut(1,1) *  gmat%orb1(i1,j) + this%ut(1,2) *  gmat%orb1(i2,j)
           v2(j)   =  this%ut(2,1) *  gmat%orb1(i1,j) + this%ut(2,2) *  gmat%orb1(i2,j) 
        enddo
        do j = 1,n2
            gmat%orb1(i1,j) = v1(j)
            gmat%orb1(i2,j) = v2(j)
        enddo   
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif

  if (nflag.eq.1) then
!$OMP PARALLEL &
!$OMP PRIVATE ( i, is, i1, i2, j, v1, v2)
!$OMP DO
     do i = 1,latt%nn_lf
        i1 = latt%nnlf_list(i)
        i2 = latt%nnlist(i1,nf)
        is = this%conf_v(i,nf,ntau)
        do j = 1,n2
           gmat%orb1(i1,j) =  this%bmat_v_p_orb1_inv(is) * gmat%orb1(i1,j)  
           gmat%orb1(i2,j) =  this%bmat_v_m_orb1_inv(is) * gmat%orb1(i2,j) 
        enddo
        do j = 1,n2
           v1(j)   =  this%u(1,1) * gmat%orb1(i1,j) +  this%u(1,2) * gmat%orb1(i2,j) 
           v2(j)   =  this%u(2,1) * gmat%orb1(i1,j) +  this%u(2,2) * gmat%orb1(i2,j) 
        enddo
        do j = 1,n2
           gmat%orb1(i1,j) = v1(j)
           gmat%orb1(i2,j) = v2(j)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif

#IFDEF TIMING
  call cpu_time_now(endtime)
  timecalculation(10)=timecalculation(10)+endtime-starttime
#ENDIF
  deallocate( v2, v1)
end subroutine dqmc_left_backward_prop
