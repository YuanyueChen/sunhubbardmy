    subroutine dqmc_proj_sweep_0b(lupdate, lmeasure_equaltime, lmeasure_dyn )
      implicit none
      logical, intent(in) :: lupdate, lmeasure_equaltime, lmeasure_dyn
      ! local variables
      integer :: nt, n, nf, nl, nr, i, nflag
      type(rfunc) :: logwtmp
      type(zfunc) :: zlogwtmp, logweightf_tmp
      type(gfunc) :: ultmp, urtmp, ULRtmp, ULRINVtmp, ul_t, ur_t, gftmp, gfctmp
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
      call allocate_gfunc(ul_t,ne,ndim)
      call allocate_gfunc(ur_t,ndim,ne)
      call allocate_gfunc(gftmp,ndim,ndim)
      call allocate_gfunc(gfctmp,ndim,ndim)

      ! at this point, Ust(0:nst-1) contain UL([0;wrap_step(2,1:nst-1)])^T
      do nl = 1, ne
         do i = 1,ndim
            UR%blk1(i,nl) = projR%blk1(i,nl)
            Ust(nst)%blk1(i,nl) = dconjg(projL%blk1(i,nl))
         end do
      end do
      logwDV(0) = rfunc(0.d0)
      ! now Ust(0:nst) contain UL([0;wrap_step(2,1:nst)])^T
      ! Ust = [UL(0)^T, UL(nwrap)^T, UL(2*nwrap)^T, ..., UL(nst*nwrap)^T]

      do nt = 1, ltrot, 1

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
          ! wrap H0/2
          if( lwrapT ) then
#IFDEF BREAKUP_T
              do nf = 1, latt%nn_nf
                  call h0c%left_forward_prop(nf,UR)
                  call h0c%right_backward_prop(nf,UL)
              end do
#ELSE
              call h0c%left_forward_prop(nf,UR)
              call h0c%right_backward_prop(nf,UL)
#ENDIF
          end if

          ! updateplqu
          if( lwrapplqu ) then
            do i = 1, lq
              call hconf%left_forward_prop(UR,i,nt)
              call hconf%right_backward_prop(UL,i,nt)
              if(lupdate .and. lupdateplqu) call hconf%proj_update_plqu(i,nt,UL,UR,ULRINV )
            end do
          end if
          ! updatev
          if( lwrapv ) then
            do nf = 1, latt%nn_nf
                 nflag = 2
                 call v0conf%left_forward_prop(UR,nt,nf,nflag)
                 call v0conf%right_backward_prop(UL,nt,nf,nflag)
                 if(lupdate .and. lupdatev) call v0conf%proj_update(nt,nf,UL,UR,ULRINV )
                 nflag = 1
                 call v0conf%left_forward_prop(UR,nt,nf,nflag)
                 call v0conf%right_backward_prop(UL,nt,nf,nflag)
            end do
          end if
          ! updateu
          if( lwrapu ) then
            call u0conf%left_forward_prop(UR,nt)
            call u0conf%right_backward_prop(UL,nt)
            if(lupdate .and. lupdateu) call u0conf%proj_update(nt,UL,UR,ULRINV)
          end if

          ! wrap H0/2
          if( lwrapT ) then
              !call mmthr_rev  (grup)
              !call mmthlm1_rev(grup)
#IFDEF BREAKUP_T
              do nf = latt%nn_nf, 1, -1
                  call h0c%left_forward_prop(nf,UR)
                  call h0c%right_backward_prop(nf,UL)
              end do
#ELSE
              call h0c%left_forward_prop(nf,UR)
              call h0c%right_backward_prop(nf,UL)
#ENDIF
          end if
#IFDEF TIMING
          call cpu_time_now(time2)
          timecalculation(3)=timecalculation(3)+time2-time1
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
          call cpu_time_now(time1)
          timecalculation(4)=timecalculation(4)+time1-time2
#ENDIF

          !if ( mod(nt, nwrap) .eq. 0 ) then
          if ( iwrap_nt(nt) .gt. 0 ) then
              !n = nt/nwrap
              n = iwrap_nt(nt)
              ! at tau = n * tau1
              ! read UL
#IFDEF WRAPERROR
              ultmp = UL
              urtmp = UR
#ENDIF
              do nl = 1, ne
                  do nr = 1, ndim
                      UL%blk1(nl,nr) = Ust(n)%blk1(nr,nl)
                  end do
              end do
              call dqmc_proj_stablize_0b_qr(UR,logwtmp)
              logwDV(n) = logwDV(n-1) + logwtmp

              call zgemm('n','n',ne,ne,ndim,cone,UL%blk1,ne,UR%blk1,ndim,czero,ULR%blk1,ne)  ! ULR = UL*UR
              ULRINV = ULR
              !call s_inv_z(ne,ULRINV%blk1)
              call s_inv_logdet_lu_z(ne,ULRINV%blk1,zlogwtmp%blk1)
              if( n== nst ) then
                  logweightf_tmp%blk1 = dcmplx(logwDV(n)%blk1,0.d0)+zlogwtmp%blk1
              else
                  logweightf_tmp%blk1 = dcmplx(logwDV(n)%blk1 + logwDV(n+1)%blk1,0.d0) + zlogwtmp%blk1
              end if
              call set_phase(logweightf_tmp, phase_tmp )
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
              write(fout, '(a,i5,a)') 'nt = ', nt, ' progating gf(:,:) = '
              do i = 1, ndim
                  write(fout,'(36(2f7.4))') gftmp%blk1(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch gf(:,:) = '
              do i = 1, ndim
                  write(fout,'(36(2f7.4))') gf%blk1(i,:)
              end do
#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
              ! store UR
              Ust(n) = UR
          end if
#IFDEF TIMING
          call cpu_time_now(time2)
          timecalculation(7)=timecalculation(7)+time2-time1
#ENDIF

          !write(fout, '(a,i5,a)') 'nt = ', nt, ' ul = '
          !do i = 1, ne
          !    write(fout,'(8(2e12.4))') UL(i,:)
          !end do

          !write(fout, '(a,i5,a)') 'nt = ', nt, ' ur = '
          !do i = 1, ndim
          !    write(fout,'(8(2e12.4))') UR(i,:)
          !end do

          if(ltau .and. nt.eq.ntauin .and. lmeasure_dyn ) then
              ul_t = UL
              ur_t = UR
              call dqmc_proj_dyn(Ust, ul_t, ur_t, xmax_dyn)
              nobst = nobst + 1
          end if
#IFDEF TIMING
          call cpu_time_now(time1)
          timecalculation(5)=timecalculation(5)+time1-time2
#ENDIF
  
      end do

      ! at this point, Ust(1:nst) contain UR(wrap_step(2,1:nst)), Ust(0) reserves outdated UL(0)^T
      ! Ust = [outdated UL(0)^T, UR(nwrap), UR(2*nwrap), ..., UR(ltrot)]

      call deallocate_gfunc(ultmp)
      call deallocate_gfunc(urtmp)
      call deallocate_gfunc(ULRtmp)
      call deallocate_gfunc(ULRINVtmp)
      call deallocate_gfunc(ul_t)
      call deallocate_gfunc(ur_t)
      call deallocate_gfunc(gftmp)
      call deallocate_gfunc(gfctmp)

#IFDEF TIMING
      call cpu_time_now(endtime)
      timecalculation(1)=timecalculation(1)+endtime-starttime
#ENDIF
    end subroutine dqmc_proj_sweep_0b
