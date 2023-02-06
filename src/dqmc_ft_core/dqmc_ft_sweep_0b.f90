    subroutine dqmc_ft_sweep_0b(lupdate, lmeasure_equaltime, lmeasure_dyn )
      implicit none
      logical, intent(in) :: lupdate, lmeasure_equaltime, lmeasure_dyn
      ! local variables
      integer :: nt, n, nf, i, j, info, nt1, nt2
      logical :: lterminate 
      real(dp) :: tmp, ratiof, ratiofi
      type(gfunc) :: gf_tmp, g0t_tmp, gt0_tmp, UR, VR, UL, VL
      type(dfunc) :: DRvec, DLvec
      type(zfunc) :: logdetQR, logdetQL, logweightf_tmp
      call allocate_gfunc(gf_tmp,ndim,ndim)
      call allocate_gfunc(g0t_tmp,ndim,ndim)
      call allocate_gfunc(gt0_tmp,ndim,ndim)
      call allocate_gfunc(UR,ndim,ndim)
      call allocate_gfunc(VR,ndim,ndim)
      call allocate_dfunc(DRvec,ndim)
      call allocate_gfunc(UL,ndim,ndim)
      call allocate_gfunc(VL,ndim,ndim)
      call allocate_dfunc(DLvec,ndim)
      ! at tau = 0
      if(nst.gt.0) then
      Ust(0) = Ifmat
      Dst(0) = Ifvec
      Vst(0) = Ifmat
      logdetQst(0) = zfunc(czero)
      end if


!!!!#include "stglobal.f90"

      if( ltau ) then
          g00 = gf
          gt0 = gf
          g0t = gf - Ifmat
      end if
      if( ltau .and. lmeasure_dyn .and. (.not.lupdate) ) then
          call dyn_measure(1,gt0,g0t,gf,g00)
      end if ! if( ltau .and. lmeasure_dyn .and. (.not.lupdate) ) then
  
      do nt = 1, ltrot, 1

#IFDEF TEST
          write(fout,*)
          write(fout,'(a)') " ----------------"
          write(fout, '(a,i8)') ' |=> nt = ',  nt
          write(fout,'(a)') " ----------------"
          write(fout,*)
#ENDIF
          ! wrap H0/2
          if( lwrapT ) then
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                  call h0c%left_forward_prop(nf,gf)
                  call h0c%right_backward_prop(nf,gf)
              end do
#ELSE
              call h0c%left_forward_prop(nf,gf)
              call h0c%right_backward_prop(nf,gf)
#ENDIF
              !call mmthr  (grup)
              !call mmthlm1(grup)
          end if

          ! updateplqu
          if( lwrapplqu ) then
            do i = 1, lq
              !nflag = 5 ! onsite u for f fermion
              !call mmuur  ( grup, ity, nf, nt, nflag )
              !call mmuulm1( grup, ity, nf, nt, nflag )
              call hconf%left_forward_prop(gf,i,nt)
              call hconf%right_backward_prop(gf,i,nt)
              !if(lupdate) call upgradeuf( nt, grup )
              if(lupdate .and. lupdateplqu) call hconf%update_plqu(gf,i,nt)
            end do
          end if

          ! updateu
          if( lwrapu ) then
            call u0conf%left_forward_prop(gf,nt)
            call u0conf%right_backward_prop(gf,nt)
            if(lupdate .and. lupdateu) call u0conf%update_u(gf,nt)
          end if

          ! wrap H0/2
          if( lwrapT ) then
              !call mmthr_rev  (grup)
              !call mmthlm1_rev(grup)
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                  call h0c%left_forward_prop(nf,gf)
                  call h0c%right_backward_prop(nf,gf)
              end do
#ELSE
              call h0c%left_forward_prop(nf,gf)
              call h0c%right_backward_prop(nf,gf)
#ENDIF
          end if

          ! obser
          !!!if( lmeasure_equaltime .and. ( abs(nt-nt_ob) .le. obs_segment_len .or. abs(nt-nt_ob) .ge. (ltrot-obs_segment_len) ) ) then
          if( lmeasure_equaltime ) then
             call equaltime_measure(nt,gf,gfc)
          end if
  
          if ( iwrap_nt(nt) .gt. 0 ) then
              n = iwrap_nt(nt)
              ! at tau = n * tau1
              VL    = Vst(n)
              DLvec = Dst(n)
              UL    = Ust(n)
              logdetQL = logdetQst(n)

              call dqmc_ft_stablize_0b_qr(n)

              UR    = Ust(n)
              DRvec = Dst(n)
              VR    = Vst(n)
              logdetQR = logdetQst(n)

              if( .not. ltau .or. .not. lmeasure_dyn ) then
                  call green_equaltime( n, ndim, UR%orb1, DRvec%orb1, VR%orb1, VL%orb1, DLvec%orb1, UL%orb1, gf_tmp%orb1, logdetQR%orb1, logdetQL%orb1, logweightf_tmp%orb1, info )
              else
              ! only when we need measure dynamical quantities, we will call green_tau
#IFDEF DYNERROR
                  ! B(nt1,nt2) with nt1 >= nt2
                  nt1 = nt
                  nt2 = nt
                  ! G(t',0) = B(t',t) * G(t,0)
                  call Bmat_left_forward( nt1, nt2, gt0 )

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call Bmat_right_backward( nt1, nt2, g0t )
#ENDIF

                  call  green_tau(n, ndim, UR%orb1, DRvec%orb1, VR%orb1, VL%orb1, DLvec%orb1, UL%orb1, g00%orb1, gt0_tmp%orb1,  g0t_tmp%orb1, gf_tmp%orb1, logdetQR%orb1, logdetQL%orb1, logweightf_tmp%orb1, info )
#IFDEF TEST_LEVEL3
                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' progating gt0%orb1(:,:) = '
                  do i = 1, ndim
                      write(fout,'(18(2e12.4))') gt0%orb1(i,:)
                  end do

                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch gt0_tmp%orb1(:,:) = '
                  do i = 1, ndim
                      write(fout,'(18(2e12.4))') gt0_tmp%orb1(i,:)
                  end do

#ENDIF

#IFDEF DYNERROR
                  call s_compare_max_z( ndim, gt0%orb1, gt0_tmp%orb1, xmax_dyn_tmp )
                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
                  logdyncount = logdyncount + 1.d0
                  avglog10dynerror = avglog10dynerror + dlog10( xmax_dyn_tmp )
#ENDIF
                  gt0 = gt0_tmp
#IFDEF TEST
                  write(fout, '(a,e16.8)') 'gt0up, xmax_dyn_tmp = ',  xmax_dyn_tmp
#ENDIF
#IFDEF DYNERROR
                  call s_compare_max_z( ndim, g0t%orb1, g0t_tmp%orb1, xmax_dyn_tmp )
                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
                  logdyncount = logdyncount + 1.d0
                  avglog10dynerror = avglog10dynerror + dlog10( xmax_dyn_tmp )
#ENDIF
                  g0t = g0t_tmp
#IFDEF TEST
                  write(fout, '(a,e16.8)') 'g0tup, xmax_dyn_tmp = ',  xmax_dyn_tmp
#ENDIF
              end if
              call s_compare_max_z( ndim, gf_tmp%orb1, gf%orb1, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp

              call set_phase(logweightf_tmp, phase_tmp)

              max_phasediff_tmp = abs( phase - phase_tmp )
              if( max_phasediff_tmp .gt. max_phasediff ) max_phasediff = max_phasediff_tmp

              logcount = logcount + 1.d0
              avglog10error = avglog10error + dlog10( max_wrap_error_tmp )
              avgexpphasediff = avgexpphasediff + dexp( max_phasediff_tmp )
#IFDEF TEST
              write(fout,'(a,2e16.8)') 'progating phase = ', phase
              write(fout,'(a,2e16.8)') 'scratch phase = ', phase_tmp
#ENDIF
#IFDEF TEST_LEVEL3
              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' progating gf%orb1(:,:) = '
              do i = 1, ndim
                  write(fout,'(18(2e12.4))') gf%orb1(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch gf_tmp%orb1(:,:) = '
              do i = 1, ndim
                  write(fout,'(18(2e12.4))') gf_tmp%orb1(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst(', n, ' )%orb1(:) before wrap = '
              write(fout,'(18(e16.8))') DLvec%orb1(:)

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst(', n, ' )%orb1(:) after wrap = '
              write(fout,'(18(e16.8))') Dst(n)%orb1(:)

#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' gf, max_wrap_error_tmp = ',  max_wrap_error_tmp
              write(fout, '(a,e16.8)') ' gf, max_phasediff = ',  max_phasediff
#ENDIF
              if( info .eq. 0 ) then
                  gf = gf_tmp
                  logweightf = logweightf_tmp
                  phase = phase_tmp
              end if

          end if
  

#include "dqmc_ft_core/dyn.f90"
  
      end do

      call deallocate_gfunc(gf_tmp)
      call deallocate_gfunc(UR)
      call deallocate_gfunc(VR)
      call deallocate_dfunc(DRvec)
      call deallocate_gfunc(UL)
      call deallocate_gfunc(VL)
      call deallocate_dfunc(DLvec)
    end subroutine dqmc_ft_sweep_0b
