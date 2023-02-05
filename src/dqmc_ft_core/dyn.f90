          ! dyn
          ! g00up g00dn
          ! gt0up gt0dn
          if( ltau .and. lmeasure_dyn .and. (.not.lupdate) ) then
              if( iwrap_nt(nt) .gt. 0 ) then
                  ! at stablization point, already have gt0,g0t from green_tau
              else
                  ! B(nt1,nt2) with nt1 >= nt2
                  nt1 = nt
                  nt2 = nt
                  ! G(t',0) = B(t',t) * G(t,0)       Gij(t,0) = < C_i(t) C_j(0)^\dagger >  (t>0)
                  call Bmat_left_forward( nt1, nt2, gt0 )

                  ! G(0,t') = G(0,t) * B(t',t)^-1    Gij(0,t) = - < C_j(t)^\dagger C_i(0) >  (t>0)
                  call Bmat_right_backward( nt1, nt2, g0t )
              end if

#IFDEF TEST_LEVEL3
!!!              write(fout,*)
!!!              write(fout, '(a)') ' gt0up(:,:) = '
!!!              do i = 1, ndim
!!!                  write(fout,'(4(2e12.4))') gt0up(i,:)
!!!              end do
!!!
!!!              write(fout,*)
!!!              write(fout, '(a)') ' g0tup(:,:) = '
!!!              do i = 1, ndim
!!!                  write(fout,'(4(2e12.4))') g0tup(i,:)
!!!              end do
!!!              write(fout,*)
!!!              write(fout, '(a)') ' gttup(:,:) = '
!!!              do i = 1, ndim
!!!                  write(fout,'(4(2e12.4))') grup(i,:)
!!!              end do
#ENDIF
              call dyn_measure(nt,gt0,g0t,gf,g00)

          end if ! if( ltau .and. lmeasure_dyn .and. (.not.lupdate) ) then
