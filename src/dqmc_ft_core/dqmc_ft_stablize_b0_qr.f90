    subroutine dqmc_ft_stablize_b0_qr(n)
      ! B( beta, (n-1)*tau1 ) = B( beta, n*tau1 ) * B( n*tau1, (n-1)*tau1 )
      implicit none
      integer, intent(in) :: n

      ! local
      integer :: i
      type(dfint) :: jpvt
      type(dfunc) :: Dvec1, Dvec2
      type(gfunc) :: Umat2, Vmat1, Vmat2, Bdtau1, Btmp, Vtmp
#IFDEF TIMING
      real(dp) :: starttime, endtime
#ENDIF

      call allocate_dfint(jpvt,ndim)
      call allocate_dfunc(Dvec1,ndim)
      call allocate_dfunc(Dvec2,ndim)
      call allocate_gfunc(Umat2,ndim,ndim)
      call allocate_gfunc(Vmat1,ndim,ndim)
      call allocate_gfunc(Vmat2,ndim,ndim)
      call allocate_gfunc(Bdtau1,ndim,ndim)
      call allocate_gfunc(Btmp,ndim,ndim)
      call allocate_gfunc(Vtmp,ndim,ndim)
#IFDEF TIMING
      call cpu_time_now(starttime)
#ENDIF
      Bdtau1 = Ust(n)
      ! Bdtau1 = U*B(tau+dtau,tau)
      ! call Bmat_right_forward( wrap_step(2,n), wrap_step(1,n), Bdtau1 )
      call Bmat_left_forward_hc( wrap_step(2,n), wrap_step(1,n), Bdtau1 ) ! we need to get ( U*B(n*tau1,(n-1)*tau1) )^H = B^H * U^H
      if( n.eq.nst ) then
          !call s_zgeQRPT(ndim, ndim, Bdtau1%orb1, Umat2%orb1, Vtmp%orb1, jpvt%orb1) ! we have Vtmp^+ and Umat^+, actually
          call s_zgeQRPT_logdetQ(ndim, ndim, Bdtau1%orb1, Umat2%orb1, Vtmp%orb1, jpvt%orb1, logdetQst(n-1)%orb1) ! we have Vtmp^+ and Umat^+, actually
          do i = 1, ndim
              Dvec2%orb1(i) = Vtmp%orb1(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2%orb1, Vtmp%orb1, Vmat2%orb1 )
          call s_zmp(ndim,ndim,Vmat2%orb1,jpvt%orb1)  ! we have Vmat2^+, actually
      else
          Dvec1 = Dst(n)
          Vmat1 = Vst(n)
          call s_z_x_diag_d(ndim,Bdtau1%orb1,Dvec1%orb1,Btmp%orb1) ! Btmp = Bdtau1_up * Dmat1
          call s_zmcpt(ndim,ndim,Btmp%orb1,jpvt%orb1)
          !call s_zgeQR(ndim,ndim,Btmp%orb1,Umat2%orb1,Vmat2%orb1)
          call s_zgeQR_logdetQ(ndim,ndim,Btmp%orb1,Umat2%orb1,Vmat2%orb1, logdetQst(n-1)%orb1)
          do i = 1, ndim
              Dvec2%orb1(i) = Vmat2%orb1(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2%orb1, Vmat2%orb1, Vtmp%orb1 )
          call s_zpm(ndim,ndim,jpvt%orb1,Vmat1%orb1)
          call zgemm('n','n',ndim,ndim,ndim,cone,Vtmp%orb1,ndim,Vmat1%orb1,ndim,czero,Vmat2%orb1,ndim)
      end if
      Ust(n-1) = Umat2  ! note we have Umat2^+ and Vmat2^+
      Dst(n-1) = Dvec2
      Vst(n-1) = Vmat2
#IFDEF TIMING
      call cpu_time_now(endtime)
      timecalculation(7)=timecalculation(7)+endtime-starttime
#ENDIF

      call deallocate_dfint(jpvt)
      call deallocate_dfunc(Dvec1)
      call deallocate_dfunc(Dvec2)
      call deallocate_gfunc(Umat2)
      call deallocate_gfunc(Vmat1)
      call deallocate_gfunc(Vmat2)
      call deallocate_gfunc(Bdtau1)
      call deallocate_gfunc(Btmp)
      call deallocate_gfunc(Vtmp)
    end subroutine dqmc_ft_stablize_b0_qr
