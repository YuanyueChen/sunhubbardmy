subroutine dqmc_set_proj
  implicit none

  type(gfunc) :: proj
  complex(dp), dimension(:,:), allocatable :: tmp
  real(dp), dimension(:), allocatable ::  wc
  real(dp) :: en_free, degen
  real(dp) :: rt1, mu1, xmag1, flux_x1, flux_y1, dimer1, rndness1
  integer :: i

  call allocate_gfunc(proj, ndim, ndim)
  allocate(tmp(ndim,ndim), wc(ndim))
  tmp  = czero 
  
  rt1=rt; mu1=mu; xmag1=xmag; flux_x1=flux_x; flux_y1=flux_y; dimer1=dimer; rndness1=rndness;
  call set_trial_h0(tmp, rt1, mu1, xmag1, flux_x1, flux_y1, dimer1, rndness1)
  !do i = 1, ndim
  !    ! add small mass
  !    if( mod(i,2) == 0 ) then
  !        tmp(i,i) = tmp(i,i) + dcmplx(0.d0, 0.00001d0)
  !    else
  !        tmp(i,i) = tmp(i,i) - dcmplx(0.d0, 0.00001d0)
  !    end if
  !end do

#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' h0mat(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4,2x))') tmp(i,:)
  end do
#ENDIF

  call s_eig_he(ndim,ndim,tmp,wc,proj%blk1)
#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' proj(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4))') proj%blk1(i,:)
  end do
  write(fout,*)
  write(fout,'(a)') ' wc(:) = '
  do i = 1, ndim
      write(fout,'(f10.6)') wc(i)
  end do
#ENDIF

#IFDEF TEST
  tmp = proj%blk1
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

  projL%blk1 = proj%blk1
  projR%blk1 = proj%blk1

  deallocate(tmp)
  deallocate(wc)

end subroutine dqmc_set_proj

subroutine dqmc_set_proj_mixed
  implicit none

  complex(dp), dimension(:,:), allocatable :: tmp
  real(dp), dimension(:), allocatable ::  wc
  real(dp) :: en_free, degen
  real(dp) :: rt1, mu1, xmag1, flux_x1, flux_y1, dimer1, rndness1
  integer :: i

  allocate(tmp(ndim,ndim), wc(ndim))
  tmp  = czero 

  ! ---------------------- set up left trial wave function --------------------- !
  !! trial Hamiltonian
  rt1=rt; mu1=mu; xmag1=xmag; flux_x1=flux_x; flux_y1=flux_y; dimer1=dimer; rndness1=rndness;
  call set_trial_h0(tmp, rt1, mu1, xmag1, flux_x1, flux_y1, dimer1, rndness1)
  !do i = 1, ndim
  !    ! add small mass
  !    if( mod(i,2) == 0 ) then
  !        tmp(i,i) = tmp(i,i) + dcmplx(0.d0, 0.00001d0)
  !    else
  !        tmp(i,i) = tmp(i,i) - dcmplx(0.d0, 0.00001d0)
  !    end if
  !end do

#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' trial_h0mat_L(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4,2x))') tmp(i,:)
  end do
#ENDIF

  !! trial Slater determinant
  call s_eig_he(ndim,ndim,tmp,wc,projL%blk1)
#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' projL(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4))') projL%blk1(i,:)
  end do
  write(fout,*)
  write(fout,'(a)') ' wc(:) = '
  do i = 1, ndim
      write(fout,'(f10.6)') wc(i)
  end do
  tmp = projL%blk1
  call s_invlu_z(ndim,tmp)
  write(fout,*)
  write(fout,'(a)') ' inv(projL)(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4))') tmp(i,:)
  end do
  write(fout,*)
#ENDIF

  !! trial ground-state energy and degeneracy
  en_free = 0.d0
  do i = 1,ne
    en_free = en_free + wc(i)
  enddo
  en_free = en_free

  degen = wc(ne+1) - wc(ne)

  if (irank == 0) then
      write(50,*) 'degen_L: ', degen
      write(50,*) 'en_free_L: ', en_free
  end if

  ! --------------------- set up right trial wave function --------------------- !
  !! trial Hamiltonian
  rt1=rt; mu1=mu; xmag1=xmag; flux_x1=flux_x; flux_y1=flux_y; dimer1=dimer; rndness1=rndness;
  call set_trial_h0(tmp, rt1, mu1, xmag1, flux_x1, flux_y1, dimer1, rndness1)
  !do i = 1, ndim
  !    ! add small mass
  !    if( mod(i,2) == 0 ) then
  !        tmp(i,i) = tmp(i,i) + dcmplx(0.d0, 0.00001d0)
  !    else
  !        tmp(i,i) = tmp(i,i) - dcmplx(0.d0, 0.00001d0)
  !    end if
  !end do

#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' trial_h0mat_R(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4,2x))') tmp(i,:)
  end do
#ENDIF

  !! trial Slater determinant
  call s_eig_he(ndim,ndim,tmp,wc,projR%blk1)
#IFDEF TEST
  write(fout,*)
  write(fout,'(a)') ' projR(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4))') projR%blk1(i,:)
  end do
  write(fout,*)
  write(fout,'(a)') ' wc(:) = '
  do i = 1, ndim
      write(fout,'(f10.6)') wc(i)
  end do
  tmp = projR%blk1
  call s_invlu_z(ndim,tmp)
  write(fout,*)
  write(fout,'(a)') ' inv(projR)(:,:) = '
  do i = 1, ndim
      write(fout,'(36(2f7.4))') tmp(i,:)
  end do
  write(fout,*)
#ENDIF

  !! trial ground-state energy and degeneracy
  en_free = 0.d0
  do i = 1,ne
     en_free = en_free + wc(i)
  enddo
  en_free = en_free

  degen = wc(ne+1) - wc(ne)

  if (irank == 0) then
      write(50,*) 'degen_R: ', degen
      write(50,*) 'en_free_R: ', en_free
  end if



  deallocate(tmp)
  deallocate(wc)

end subroutine dqmc_set_proj_mixed
