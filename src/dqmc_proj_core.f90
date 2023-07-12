module dqmc_proj_core
!#DEFINE TEST_sig
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
  use model_para
  use dqmc_para
  use dqmc_measure
  implicit none

  type(gfunc), save :: Ifmat
  type(dfunc), save :: Ifvec
  type(gfunc), save :: gf, gfc ! gfc play as gf_tmp when outside measurement
  type(gfunc), save :: proj, UR, UL, ULR, ULRINV
  type(zfunc), save :: logweightf
  type(gfunc), allocatable, dimension(:), save :: Ust, Ust_tmp
  type(zfunc), allocatable, dimension(:), save :: logdetQst, logdetQst_tmp
  type(rfunc), allocatable, dimension(:), save :: logwDV, logwDV_tmp
  type(gfunc), save :: gt0, g0t, g00

  contains

    subroutine allocate_core
      implicit none
      integer :: i
      call allocate_gfunc(Ifmat,ndim,ndim)
      Ifmat%orb1 = Imat
      call allocate_dfunc(Ifvec,ndim)
      Ifvec%orb1 = Ivec
      call allocate_gfunc(gf,ndim,ndim)
      call allocate_gfunc(gfc,ndim,ndim)
      call allocate_gfunc(proj,ndim,ndim)
      call allocate_gfunc(UR,ndim,ne)
      call allocate_gfunc(UL,ne,ndim)
      call allocate_gfunc(ULR,ne,ne)
      call allocate_gfunc(ULRINV,ne,ne)
      allocate( Ust(0:nst), Ust_tmp(0:nst) )
      allocate( logdetQst(0:nst), logdetQst_tmp(0:nst) )
      allocate( logwDV(0:nst), logwDV_tmp(0:nst) )
      do i = 0, nst
          call allocate_gfunc(Ust(i),ndim,ne)
          call allocate_gfunc(Ust_tmp(i),ndim,ne)
      end do
      call allocate_gfunc( gt0, ndim, ndim )
      call allocate_gfunc( g0t, ndim, ndim )
      call allocate_gfunc( g00, ndim, ndim )

    end subroutine allocate_core

    subroutine deallocate_core
      implicit none
      deallocate( Ust, Ust_tmp)
      deallocate( logdetQst, logdetQst_tmp)
      deallocate( logwDV, logwDV_tmp)
    end subroutine deallocate_core

#include 'dqmc_proj_core/dqmc_proj_stablize_0b_qr.f90'
#include 'dqmc_proj_core/dqmc_proj_stablize_b0_qr.f90'
#include 'dqmc_proj_core/dqmc_proj_sweep_start_0b.f90'  
#include 'dqmc_proj_core/dqmc_proj_sweep_b0.f90'
#include 'dqmc_proj_core/dqmc_proj_sweep_0b.f90'
#include 'dqmc_proj_core/dqmc_proj_dyn.f90'
#include 'dqmc_proj_core/dqmc_set_proj.f90'
#include 'dqmc_proj_core/green_equaltime.f90'
#include 'dqmc_proj_core/set_phase.f90'
  
    subroutine Bmat_left_forward( nt1, nt2, bmat )
      ! B(tau1,tau2) * 
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      type(gfunc), intent(inout) :: bmat

      ! local
      integer :: nt, nf, i, nflag

      do nt = nt2, nt1
          if( lwrapT ) then
              !call mmthr(bmat_up)
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                call h0c%left_forward_prop(nf, bmat)
              end do
#ELSE
              call h0c%left_forward_prop(nf, bmat)
#ENDIF
          end if
          if( lwrapplqu ) then
              !call mmuur( bmat_up, ity, nf, nt, nflag )
              do i = 1, lq
                  call hconf%left_forward_prop(bmat,i,nt)
              end do
          end if
         
          if( lwrapv ) then
                do nf = 1, latt%nn_nf
                    do nflag = 2, 1, -1
                      call v0conf%left_forward_prop(bmat,nt,nf,nflag)
                    end do
                end do
          end if

          if( lwrapu ) then
              call u0conf%left_forward_prop(bmat,nt)
          end if

          if( lwrapT ) then
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                call h0c%left_forward_prop(nf, bmat)
              end do
#ELSE
              call h0c%left_forward_prop(nf, bmat)
#ENDIF
          end if
      end do

    end subroutine Bmat_left_forward

    subroutine Bmat_left_forward_hc( nt1, nt2, bmat )
      ! B(tau1,tau2) * 
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      type(gfunc), intent(inout) :: bmat

      ! local
      integer :: nt, nf, i, nflag

      do nt = nt1, nt2, -1

          if( lwrapT ) then
              !call mmthr_revH(bmat_up)
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                  call h0c%left_forward_prop_hc(nf, bmat)
              end do
#ELSE
              call h0c%left_forward_prop_hc(nf, bmat)
#ENDIF
          end if

          if( lwrapu ) then
              call u0conf%left_forward_prop_hc(bmat,nt)
          end if
          
          if( lwrapv ) then
                do nf = latt%nn_nf, 1, -1
                    do nflag = 2, 1, -1
                      call v0conf%left_forward_prop_hc(bmat,nt,nf,nflag)
                    end do
                end do
          end if

          if( lwrapplqu ) then
              do i = lq, 1, -1
                  !call mmuurH( bmat_up, ity, nf, nt, nflag )
                  call hconf%left_forward_prop_hc(bmat,i,nt)
              end do
          end if

          if( lwrapT ) then
              !call mmthrH(bmat_up)
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                  call h0c%left_forward_prop_hc(nf, bmat)
              end do
#ELSE
              call h0c%left_forward_prop_hc(nf, bmat)
#ENDIF
          end if
      end do
    end subroutine Bmat_left_forward_hc

    subroutine Bmat_right_forward( nt1, nt2, bmat )
      ! * B(tau1,tau2)
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      type(gfunc), intent(inout) :: bmat

      ! local
      integer :: nt, nf, i, nflag

      do nt = nt1, nt2, -1
          ! wrap H0/2
          if( lwrapT ) then
              !call mmthl_rev(bmat_up)
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                  call h0c%right_forward_prop(nf, bmat)
              end do
#ELSE
              call h0c%right_forward_prop(nf, bmat)
#ENDIF
          end if

          if( lwrapu ) then
              call u0conf%right_forward_prop(bmat,nt)
          end if
          
          if( lwrapv ) then
                do nf = latt%nn_nf, 1, -1
                    do nflag = 2, 1, -1
                      call v0conf%right_forward_prop(bmat,nt,nf,nflag)
                    end do
                end do
          end if

          if( lwrapplqu ) then
              !call mmuul( bmat_up, ity, nf, nt, nflag )
              do i = lq, 1, -1
                  call hconf%right_forward_prop(bmat,i,nt)
              end do
          end if
  
          ! wrap H0/2
          if( lwrapT ) then
              !call mmthl(bmat_up)
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                  call h0c%right_forward_prop(nf, bmat)
              end do
#ELSE
              call h0c%right_forward_prop(nf, bmat)
#ENDIF
          end if
      end do
    end subroutine Bmat_right_forward

    subroutine Bmat_right_backward( nt1, nt2, bmat )
      ! *B(tau1,tau2)^-1
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      type(gfunc) :: bmat

      ! local
      integer :: nt, nf, i, nflag

      do nt = nt2, nt1
          ! wrap H0/2
          if( lwrapT ) then
              !call mmthlm1(bmat_up)
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                  call h0c%right_backward_prop(nf, bmat)
              end do
#ELSE
              call h0c%right_backward_prop(nf, bmat)
#ENDIF
          end if

          if( lwrapplqu ) then
              !call mmuulm1( bmat_up, ity, nf, nt, nflag )
              do i = 1, lq
                  call hconf%right_backward_prop(bmat,i,nt)
              end do
          end if
          
          if( lwrapv ) then
                do nf = 1, latt%nn_nf
                    do nflag = 2, 1, -1
                      call v0conf%right_backward_prop(bmat,nt,nf,nflag)
                    end do
                end do
          end if

          if( lwrapu ) then
              call u0conf%right_backward_prop(bmat,nt)
          end if

          ! wrap H0/2
          if( lwrapT ) then
              !call mmthlm1_rev(bmat_up)
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                  call h0c%right_backward_prop(nf, bmat)
              end do
#ELSE
              call h0c%right_backward_prop(nf, bmat)
#ENDIF
          end if
      end do
    end subroutine Bmat_right_backward

    subroutine push_stage
      implicit none
      integer :: i
      do i = 0, nst
          Ust_tmp(i) = Ust(i)
          logdetQst_tmp(i) = logdetQst(i)
          logwDV_tmp(i) = logwDV(i)
      end do
      gfc = gf
    end subroutine push_stage

    subroutine pop_stage
      implicit none
      integer :: i
      do i = 0, nst
          Ust(i) =  Ust_tmp(i)
          logdetQst(i) = logdetQst_tmp(i)
          logwDV(i) = logwDV_tmp(i)
      end do
      gf = gfc
    end subroutine pop_stage

end module dqmc_proj_core
