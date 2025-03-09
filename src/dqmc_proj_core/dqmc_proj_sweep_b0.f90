    subroutine dqmc_proj_sweep_b0(lupdate, lmeasure_equaltime)
      implicit none
      logical, intent(in) :: lupdate, lmeasure_equaltime
      ! local variables
      integer :: nt, n, nf, nl, i, nflag
      type(rfunc) :: logwtmp
      type(zfunc) :: zlogwtmp, logweightf_tmp
      type(gfunc) :: ultmp, urtmp, ULRtmp, ULRINVtmp, gftmp, gfctmp
#IFDEF TIMING
      real(dp) :: starttime, endtime, time1, time2
#ENDIF
#IFDEF TIMING
      call cpu_time_now(starttime)
#ENDIF

      call allocate_gfunc(ultmp,ne,ndim)
      call allocate_gfunc(urtmp,ndim,ne)
      call allocate_gfunc(ULRtmp,ne,ne)
      call allocate_gfunc(ULRINVtmp,ne,ne)
      call allocate_gfunc(gftmp,ndim,ndim)
      call allocate_gfunc(gfctmp,ndim,ndim)
      
      ! at this point, Ust(1:nst) contains UR(wrap_step(2,1:nst))
      ! UL at tau = beta
      do nl = 1, ne
          do i = 1, ndim
              UL%blk1(nl,i) = dconjg(projL%blk1(i,nl))
              Ust(0)%blk1(i,nl) = projR%blk1(i,nl)
          end do
      end do
      ! now Ust(0:nst) contains UR([0;wrap_step(2,1:nst)])

      do nt = ltrot, 1, -1
#IFDEF TEST
          write(fout,*)
          write(fout,'(a)') " ----------------"
          write(fout, '(a,i8)') ' |=> nt = ',  nt
          write(fout,'(a)') " ----------------"
          write(fout,*)
#ENDIF
#IFDEF TIMING
          call cpu_time_now(time1)
#ENDIF
          ! obser
          if( lmeasure_equaltime .and. nt > (ltrot/2 - (obs_eqt_mid_len+1)/2) .and. nt <= (ltrot/2 + obs_eqt_mid_len/2) ) then
              call zgemm('n','n',ne,ne,ndim,cone,UL%blk1,ne,UR%blk1,ndim,czero,ULR%blk1,ne)  ! ULR = UL*UR
              ULRINV = ULR
              call s_invlu_z(ne,ULRINV%blk1)
              call green_equaltime(UL%blk1,UR%blk1,ULRINV%blk1,gf%blk1,gfc%blk1)
              call equaltime_measure(nt,gf,gfc)
          end if
#IFDEF TIMING
          call cpu_time_now(time2)
          timecalculation(4)=timecalculation(4)+time2-time1
#ENDIF

          ! wrap H0/2
          if( lwrapT ) then
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                  call h0c%right_forward_prop(nf,UL)
                  call h0c%left_backward_prop(nf,UR)
              end do
#ELSE
              call h0c%right_forward_prop(nf,UL)
              call h0c%left_backward_prop(nf,UR)
#ENDIF
              !call mmthl_rev  (gf)
              !call mmthrm1_rev(gf)
          end if
  
          !! update
          ! updateu
          if( lwrapu ) then
            if(lupdate .and. lupdateu) call u0conf%proj_update(nt,UL,UR,ULRINV)
            call u0conf%right_forward_prop(UL,nt)
            call u0conf%left_backward_prop(UR,nt)
          end if
          ! updatev
          if( lwrapv ) then
            do  nf = latt%nn_nf,1,-1
                nflag = 2
                call v0conf%right_forward_prop(UL,nt,nf,nflag)
                call v0conf%left_backward_prop(UR,nt,nf,nflag)
                if(lupdate .and. lupdatev) call v0conf%proj_update(nt,nf,UL,UR,ULRINV)
                nflag = 1
                call v0conf%right_forward_prop(UL,nt,nf,nflag)
                call v0conf%left_backward_prop(UR,nt,nf,nflag)
            end do
          end if
          ! updateplqu
          if( lwrapplqu ) then
              !if(lupdate) call upgradeuf( nt, grup)
            do i = lq, 1, -1
              if(lupdate .and. lupdateplqu) call hconf%proj_update_plqu(i,nt,UL,UR,ULRINV)
              call hconf%right_forward_prop(UL,i,nt)
              call hconf%left_backward_prop(UR,i,nt)
            end do
          end if

          ! wrap H0/2
          if( lwrapT ) then
              !call mmthl  (grup)
              !call mmthrm1(grup)
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                  call h0c%right_forward_prop(nf,UL)
                  call h0c%left_backward_prop(nf,UR)
              end do
#ELSE
              call h0c%right_forward_prop(nf,UL)
              call h0c%left_backward_prop(nf,UR)
#ENDIF
          end if

#IFDEF TIMING
          call cpu_time_now(time1)
          timecalculation(3)=timecalculation(3)+time1-time2
#ENDIF

          !if ( mod(nt, nwrap) .eq. 0  ) then
          if ( (iwrap_nt(nt-1).gt.0 .or. nwrap.eq.1) .or. (nt .eq. 1) ) then
              !n = nt/nwrap
              n = iwrap_nt(nt-1)

              ! at tau = n * tau1

              ! read UR and store the old for calculate grtmp
#IFDEF WRAPERROR
              urtmp = UR
              ultmp = UL
#ENDIF
              UR = Ust(n)

              ! ortho UL and store the old for calculate grtmp
              call dqmc_proj_stablize_b0_qr(UL,logwtmp)
              if( (n+1) == nst ) then
                  logwDV(n+1) = logwtmp
              else
                  logwDV(n+1) = logwDV(n+2) + logwtmp
              end if

              call zgemm('n','n',ne,ne,ndim,cone,UL%blk1,ne,UR%blk1,ndim,czero,ULR%blk1,ne)  ! ULR = UL*UR
              ULRINV = ULR
              !call s_inv_z(ne,ULRINV)
              call s_inv_logdet_lu_z(ne,ULRINV%blk1,zlogwtmp%blk1)
              if( n == 0 ) then
                  logweightf_tmp%blk1 = dcmplx(logwDV(n+1)%blk1,0.d0) + zlogwtmp%blk1
              else
                  logweightf_tmp%blk1 = dcmplx(logwDV(n+1)%blk1+logwDV(n)%blk1,0.d0) + zlogwtmp%blk1
              end if
              call set_phase(logweightf_tmp, phase_tmp)
              call green_equaltime(UL%blk1,UR%blk1,ULRINV%blk1,gf%blk1,gfc%blk1)
#IFDEF WRAPERROR
              call zgemm('n','n',ne,ne,ndim,cone,ultmp%blk1,ne,urtmp%blk1,ndim,czero,ULRtmp%blk1,ne)  ! ULRtmp = ultmp*urtmp
              ULRINVtmp = ULRtmp
              call s_invlu_z(ne,ULRINVtmp%blk1)
              call green_equaltime(ultmp%blk1,urtmp%blk1,ULRINVtmp%blk1,gftmp%blk1,gfctmp%blk1)
              call s_compare_max_z(ndim,gftmp%blk1,gf%blk1,max_wrap_error_tmp)
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
              max_phasediff_tmp = abs( phase - phase_tmp )
              if( max_phasediff_tmp .gt. max_phasediff ) max_phasediff = max_phasediff_tmp
              logcount = logcount + 1.d0
              avglog10error = avglog10error + dlog10( max_wrap_error_tmp )
              avgexpphasediff = avgexpphasediff + dexp( max_phasediff_tmp )
#ENDIF
#IFDEF TEST
              write(fout,'(a,2e16.8)') 'progating phase = ', phase
              write(fout,'(a,2e16.8)') 'scratch phase = ', phase_tmp
              !write(fout,'(a,2e24.12)') 'progating logweightf = ', logweightf%blk1
              write(fout,'(a,2e24.12)') 'scratch   logweightf = ', logweightf_tmp%blk1
#ENDIF
              phase = phase_tmp
              logweightf = logweightf_tmp
#IFDEF TEST_LEVEL3
              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' progating grup(:,:) = '
              do i = 1, ndim
                  write(fout,'(36(2f7.4))') gftmp%blk1(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch grup(:,:) = '
              do i = 1, ndim
                  write(fout,'(36(2f7.4))') gf%blk1(i,:)
              end do
#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
              ! store UL
              do nl = 1, ne
                 do i = 1, ndim
                    Ust(n)%blk1(i,nl) = UL%blk1(nl,i)
                 end do
              end do
          end if
#IFDEF TIMING
          call cpu_time_now(time2)
          timecalculation(7)=timecalculation(7)+time2-time1
#ENDIF
      end do

      ! at this point, Ust(0:nst-1) contain UL([0;wrap_step(2,1:nst-1)])^T, Ust(nst) reserves outdated UR(wrap_step(2,nst))
      ! Ust = [UL(0)^T, UL(nwrap)^T, UL(2*nwrap)^T, ..., UL((nst-1)*nwrap)^T, outdated UR(ltrot)]

      call deallocate_gfunc(ultmp)
      call deallocate_gfunc(urtmp)
      call deallocate_gfunc(ULRtmp)
      call deallocate_gfunc(ULRINVtmp)
      call deallocate_gfunc(gftmp)
      call deallocate_gfunc(gfctmp)

#IFDEF TIMING
      call cpu_time_now(endtime)
      timecalculation(1)=timecalculation(1)+endtime-starttime
#ENDIF
    end subroutine dqmc_proj_sweep_b0
