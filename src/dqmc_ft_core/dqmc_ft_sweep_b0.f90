    subroutine dqmc_ft_sweep_b0(lupdate, lmeasure_equaltime)
      implicit none
      logical, intent(in) :: lupdate, lmeasure_equaltime
      ! local variables
      integer :: nt, n, nf, i, j, info
      logical :: lterminate 
      real(dp) :: tmp, ratiof, ratiofi
      type(gfunc) :: gf_tmp, UR, VR, UL, VL
      type(dfunc) :: DRvec, DLvec
      type(zfunc) :: logdetQR, logdetQL, logweightf_tmp
#IFDEF TIMING
      real(dp) :: starttime, endtime
#ENDIF

      call allocate_gfunc(gf_tmp,ndim,ndim)
      call allocate_gfunc(UR,ndim,ndim)
      call allocate_gfunc(VR,ndim,ndim)
      call allocate_dfunc(DRvec,ndim)
      call allocate_gfunc(UL,ndim,ndim)
      call allocate_gfunc(VL,ndim,ndim)
      call allocate_dfunc(DLvec,ndim)
#IFDEF TIMING
      call cpu_time_now(starttime)
#ENDIF

      ! at tau = beta
      if(nst.gt.0) then
      Vst(nst) = Ifmat
      Dst(nst) = Ifvec
      Ust(nst) = Ifmat
      logdetQst(nst) = zfunc(czero)
      end if
  
      do nt = ltrot, 1, -1
#IFDEF TEST
          write(fout,*)
          write(fout,'(a)') " ----------------"
          write(fout, '(a,i8)') ' |=> nt = ',  nt
          write(fout,'(a)') " ----------------"
          write(fout,*)
#ENDIF
          !!!! obser
          !!!if( lmeasure_equaltime .and. ( abs(nt-nt_ob) .le. obs_segment_len .or. abs(nt-nt_ob) .ge. (ltrot-obs_segment_len) ) ) then
          if( lmeasure_equaltime ) then
             call equaltime_measure(nt,gf,gfc)
          end if

          ! wrap H0/2
          if( lwrapT ) then
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                  call h0c%right_forward_prop(nf,gf)
                  call h0c%left_backward_prop(nf,gf)
              end do
#ELSE
              call h0c%right_forward_prop(nf,gf)
              call h0c%left_backward_prop(nf,gf)
#ENDIF
              !call mmthl_rev  (gf)
              !call mmthrm1_rev(gf)
          end if
  
          !! update
          ! updateu
          if( lwrapu ) then
            if(lupdate .and. lupdateu) call u0conf%update_u(gf,nt)
            call u0conf%right_forward_prop(gf,nt)
            call u0conf%left_backward_prop(gf,nt)
          end if
          ! updateplqu
          if( lwrapplqu ) then
              !if(lupdate) call upgradeuf( nt, grup)
            do i = lq, 1, -1
              if(lupdate .and. lupdateplqu) call hconf%update_plqu(gf,i,nt)
              call hconf%right_forward_prop(gf,i,nt)
              call hconf%left_backward_prop(gf,i,nt)
              !call mmuul  ( grup, ity, nf, nt, nflag )
              !call mmuurm1( grup, ity, nf, nt, nflag )
            end do
          end if

          ! wrap H0/2
          if( lwrapT ) then
              !call mmthl  (grup)
              !call mmthrm1(grup)
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                  call h0c%right_forward_prop(nf,gf)
                  call h0c%left_backward_prop(nf,gf)
              end do
#ELSE
              call h0c%right_forward_prop(nf,gf)
              call h0c%left_backward_prop(nf,gf)
#ENDIF
          end if

          !if ( (iwrap_nt(nt-1).gt.0 .and. nwrap.eq.1) .or. (iwrap_nt(nt-1) .gt. 0 .and. nt.ne.ltrot) .or. ( nt.eq.1 .and. nst .gt. 0 ) ) then
          !if ( (iwrap_nt(nt-1).gt.0 .and. (nt.ne.ltrot .or. nwrap.eq.1)) .or. (nt.eq.1) ) then
          if ( (iwrap_nt(nt-1).gt.0 .or. nwrap.eq.1) .or. (nt .eq. 1) ) then
              n = iwrap_nt(nt-1)
               ! at tau = n * tau1
              UR    = Ust(n)
              DRvec = Dst(n)
              VR    = Vst(n)
              logdetQR = logdetQst(n)

              call dqmc_ft_stablize_b0_qr(n+1)

              UL    = Ust(n)
              DLvec = Dst(n)
              VL    = Vst(n)
              logdetQL = logdetQst(n)

              !if( .not. ltau ) then
                  call green_equaltime( n, ndim, UR%orb1, DRvec%orb1, VR%orb1, VL%orb1, DLvec%orb1, UL%orb1, gf_tmp%orb1, logdetQR%orb1, logdetQL%orb1, logweightf_tmp%orb1, info )
              !else
              !    call green_tau(n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, g00up, gt0up, g0tup, grtmp, info )
              !end if
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
              write(fout,'(18(e16.8))') DRvec%orb1(:)

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst(', n, ' )%orb1(:) after wrap = '
              write(fout,'(18(e16.8))') Dst(n)%orb1(:)

#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' gf, max_wrap_error_tmp = ',  max_wrap_error_tmp
              write(fout, '(a,e16.8)') ' gf, max_phasediff = ',  max_phasediff
#ENDIF
              ! whether use the scrath gf
              if( info .eq. 0 ) then 
                  gf = gf_tmp
                  logweightf = logweightf_tmp
                  phase = phase_tmp
              else
                  stop " Error in dqmc_ft_sweep_b0 "
              end if

          end if
      end do
#IFDEF TIMING
      call cpu_time_now(endtime)
      timecalculation(1)=timecalculation(1)+endtime-starttime
#ENDIF

      call deallocate_gfunc(gf_tmp)
      call deallocate_gfunc(UR)
      call deallocate_gfunc(VR)
      call deallocate_dfunc(DRvec)
      call deallocate_gfunc(UL)
      call deallocate_gfunc(VL)
      call deallocate_dfunc(DLvec)

    end subroutine dqmc_ft_sweep_b0
