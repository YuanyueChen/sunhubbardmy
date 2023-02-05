subroutine dqmc_set_proj
  implicit none

  complex(dp), dimension(:,:), allocatable :: tmp
  real(dp), dimension(:), allocatable ::  wc
  real(dp) :: rt1, en_free, degen
  integer :: i

  allocate(tmp(ndim,ndim), wc(ndim))
  tmp  = czero 
  
  rt1 = rt
  if (rt.le.zero)  rt1 = 1.d0
  call seth0(tmp,rt1,mu)
  !!!do i = 1, ndim
  !!!    ! add mass
  !!!    if( mod(i,2) == 0 ) then
  !!!        tmp(i,i) = tmp(i,i) + dcmplx(0.d0, 0.00001d0)
  !!!    else
  !!!        tmp(i,i) = tmp(i,i) - dcmplx(0.d0, 0.00001d0)
  !!!    end if
  !!!end do

#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' h0mat(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4,2x))') tmp(i,:)
  end do
#ENDIF

  call s_eig_he(ndim,ndim,tmp,wc,proj%orb1)
#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' proj(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4))') proj%orb1(i,:)
  end do
  write(fout,*)
  write(fout,'(a)') ' wc(:) = '
  do i = 1, ndim
      write(fout,'(f10.6)') wc(i)
  end do
#ENDIF

#IFDEF TEST
  tmp = proj%orb1
  call s_invlu_z(ndim,tmp)
  write(fout,*)
  write(fout,'(a)') ' inv(proj)(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4))') tmp(i,:)
  end do
  write(fout,*)
#ENDIF

  en_free = 0.d0
  do i = 1,ne
     en_free = en_free + wc(i)
  enddo
  en_free = en_free

  degen = wc(ne+1) - wc(ne)

  if (irank == 0) then
      write(50,*) 'degen: ', degen
      write(50,*) 'en_free: ', en_free
  end if

  deallocate(tmp)
  deallocate(wc)

end subroutine dqmc_set_proj
