subroutine dqmc_proj_dyn(ust, ul, ur, xmax_dyn1)
  implicit none

  ! storage is full with u^{<} (left) propagations.
  type(gfunc), dimension(0:nst), intent(in) :: ust
  type(gfunc), intent(inout) :: ul, ur
  real(dp), intent(inout) :: xmax_dyn1
  ! local.
  type(gfunc) :: gftmpb, gftmp, g00, g0t, gt0, gtt, temp
  type(gfunc):: ulr, ulrinv
  type(rfunc) :: logwtmp
  integer :: nt, ntau, nt_st, i, nt1, ntau1
  real(dp) :: xmax
  call allocate_gfunc( gftmpb, ndim, ndim )
  call allocate_gfunc( gftmp, ndim, ndim )
  call allocate_gfunc( g00, ndim, ndim )
  call allocate_gfunc( g0t, ndim, ndim )
  call allocate_gfunc( gt0, ndim, ndim )
  call allocate_gfunc( gtt, ndim, ndim )
  call allocate_gfunc( temp, ndim, ndim )
  call allocate_gfunc( ulr, ne, ne )
  call allocate_gfunc( ulrinv, ne, ne )
#IFDEF TEST
  write(fout,*) 'starting dyn'
#ENDIF
  do nt = ntauin, ntauin + ntdm -1
#IFDEF TEST
     write(fout,'(a,i6)') 'in dyn, nt = ', nt
#ENDIF
     ! ur is on time slice nt
     ntau = nt - ntauin
     !if ( mod(nt,nwrap).eq.0) then
     if ( iwrap_nt(nt) .gt. 0 ) then
        ! write(fout,'(a,i6)') 'wrapping, nt = ', nt
        ntau = nt - ntauin
        !nt_st = nt/nwrap
        nt_st = iwrap_nt(nt)
        if (ntau.gt.0) then
           ul = ust(nt_st)
        endif
        !write(fout, '(a)') 'ul='
        !do i = 1, ne
        !    write(fout,'(8(2e12.4,2x))') ul(i,:)
        !end do
        !write(fout, '(a)') 'ur='
        !do i = 1, ndim
        !    write(fout,'(8(2e12.4,2x))') ur(i,:)
        !end do
        call zgemm('n','n',ne,ne,ndim,cone,ul%orb1,ne,ur%orb1,ndim,czero,ulr%orb1,ne)  ! ulr = ul*ur
        !write(fout, '(a)') 'ulr='
        !do i = 1, ne
        !    write(fout,'(8(2e12.4,2x))') ulr(i,:)
        !end do
        ulrinv = ulr
        call s_invlu_z(ne,ulrinv%orb1)
        ! compute green functions.
        call green_equaltime(ul%orb1,ur%orb1,ulrinv%orb1,gftmp%orb1,temp%orb1)

        !write(fout,*) 'computed green function at nt = ', nt
        !do i = 1, ndim
        !    write(fout,'(8(2e12.4,2x))') gftmp(i,:)
        !end do
        !you have: gftmp (i,j) = <c_i c^+_j >
        gftmpb = Ifmat - gftmp
        if (ntau.eq.0) then
           g00 = gftmp
           gtt = gftmp
           gt0 = gftmp
           g0t%orb1 = -gftmpb%orb1
           nt1 = 0
           call dyn_measure(nt,gt0,g0t,gtt,g00)
        else
           xmax = 0.d0
           call s_compare_max_z(ndim, gtt%orb1, gftmp%orb1, xmax)
           if (xmax.gt.xmax_dyn1) xmax_dyn1 = xmax
           logdyncount = logdyncount + 1.d0
           avglog10dynerror = avglog10dynerror + dlog10( xmax )
#IFDEF TEST_LEVEL3
           write(fout,'(a,e16.8)') ' xmax_dyn_tmp      = ', xmax
           write(fout,*)
           write(fout, '(a,i5,a)') 'nt = ', nt, ' in dyn, gftmp(:,:) = '
           do i = 1, ndim
               write(fout,'(36(2e12.4,2x))') gftmp%orb1(i,:)
           end do
           write(fout,*)
           write(fout, '(a,i5,a)') 'nt = ', nt, ' in dyn, gtt(:,:) = '
           do i = 1, ndim
               write(fout,'(36(2e12.4,2x))') gtt%orb1(i,:)
           end do
#ENDIF
           gtt = gftmp
           call zgemm('n','n',ndim,ndim,ndim,cone,gftmp%orb1,ndim,gt0%orb1,ndim,czero,temp%orb1,ndim)
           gt0 = temp
           call zgemm('n','n',ndim,ndim,ndim,cone,g0t%orb1,ndim,gftmpb%orb1,ndim,czero,temp%orb1,ndim)
           g0t = temp
        endif
     endif ! ortho.
     ! now propagate to ntau + 1 and call dyn_measure.
     nt1 = nt + 1
     call Bmat_left_forward(nt1,nt1,gt0)
     call Bmat_right_backward(nt1,nt1,g0t)
     call Bmat_right_backward(nt1,nt1,gtt)
     call Bmat_left_forward(nt1,nt1,gtt)
     ntau1 = nt1 - ntauin
     ! write(6,*) 'dyn: calling obsett: ', ntau1
     call dyn_measure(ntau1,gt0,g0t,gtt,g00)

     call Bmat_left_forward(nt1,nt1,ur)

     !if (mod(nt1,nwrap).eq.0 .and. nt1 .ne. (ntauin + ntdm) ) then
     if ( iwrap_nt(nt) .gt. 0 .and. nt1 .ne. (ntauin + ntdm) ) then
        call dqmc_proj_stablize_0b_qr(ur,logwtmp)
        !write(6,*) 'dyn ortho on: ', ntauin, ntauin+ntdm, nt1
     endif
  enddo
  call deallocate_gfunc( gftmpb )
  call deallocate_gfunc( gftmp )
  call deallocate_gfunc( g00 )
  call deallocate_gfunc( g0t )
  call deallocate_gfunc( gt0 )
  call deallocate_gfunc( gtt )
  call deallocate_gfunc( temp )
  call deallocate_gfunc( ulr )
  call deallocate_gfunc( ulrinv )
end subroutine dqmc_proj_dyn
