subroutine dyn_output
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  implicit none
  include 'mpif.h'

  complex(dp) :: znorm
  character(40) :: ftag
  integer :: i, j, imj, nt

  if( ltau ) then
      znorm = cone / dcmplx( dble(nsweep*isize), 0.d0 )
  end if

  if (irank.eq.0) then
     if(ltau) then
         znorm = cone / dcmplx( dble(nsweep*isize), 0.d0 )
     end if
  endif

!  call mpi_reduce( gtau0_orb1_tau_bin, gtau0_orb1_tau, size(gtau0_orb1_tau), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
!  call mpi_reduce( zspsm_orb1_tau_bin, zspsm_orb1_tau, size(zspsm_orb1_tau), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( znn_orb1_tau_bin, znn_orb1_tau, size(znn_orb1_tau), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
!  call mpi_reduce( zbb_orb1_tau_bin, zbb_orb1_tau, size(zbb_orb1_tau), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
!  call mpi_reduce( zb_orb1_tau_bin, zb_orb1_tau, size(zb_orb1_tau), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )

  if( irank .eq. 0 ) then
!      zspsm_orb1_tau = zspsm_orb1_tau * znorm
!      ftag = "spsm"
!      call fourier_trans_tau_22(zspsm_orb1_tau, ftag)

      znn_orb1_tau = znn_orb1_tau * znorm
      ftag = "nn"
      call fourier_trans_tau_22(znn_orb1_tau, ftag)

!      gtau0_orb1_tau = gtau0_orb1_tau * znorm
!      ftag = "gf"
!      call fourier_trans_tau_22(gtau0_orb1_tau, ftag)

!      zbb_orb1_tau = zbb_orb1_tau * znorm
!      ftag = "bb"
!      call fourier_trans_tau(zbb_orb1_tau, ftag)

!      zb_orb1_tau = zb_orb1_tau * znorm
!      ftag = "b"
!      call fourier_trans_tau(zb_orb1_tau, ftag)

!      do nt = 1, ntdm+1
!          do j = 1, lq
!              do i = 1, lq
!                  imj = latt%imj(i,j)
!                  zbb_orb1_tau(imj,nt) = zbb_orb1_tau(imj,nt) - zb_orb1_tau(i,nt)*zb_orb1_tau(j,1)
!              end do
!          end do
!      end do
!      ftag = "bbdbg"
!      call fourier_trans_tau(zbb_orb1_tau, ftag)
  end if

end subroutine dyn_output

subroutine fourier_trans_tau_22(gr,ftag)
#ifdef _OPENMP
  USE OMP_LIB
#endif
  implicit none
  complex(dp), dimension(:,:,:,:) :: gr
  character (40) :: ftag

  ! local
  integer :: imj, nt, nk, no_i, no_j
  character (80) :: filek, filer
  real(dp) :: xk_p(2), aimj_p(2)
  complex(dp), allocatable , dimension(:,:,:,:) :: gk

  filek = trim(ftag)//trim('_ktau.bin')
  filer = trim(ftag)//trim('_rtau.bin')

  allocate (gk(lq,2,2,ntdm+1))

  gk = dcmplx(0.d0,0.d0)
  open (unit=20,file=filer,status='unknown', action="write", position="append")
  do imj = 1,lq
     do nt = 1,ntdm+1
        do nk = 1,lq
           gk(nk,:,:,nt) = gk(nk,:,:,nt) +  gr(imj,:,:,nt)/latt%zexpiqr(nk, imj)
        enddo
        write(20,'(8e16.8)') gr(imj,1:latt%nsub,1:latt%nsub,nt)
     enddo
  enddo
  close(20)

  gk = gk/dcmplx(dble(lq),0.d0)

  open (unit=20,file=filek,status='unknown', action="write", position="append")
  do nk = 1,lq
     write(20,'(2f16.8)') latt%listk(nk,1)/dble(latt%l1), latt%listk(nk,2)/dble(latt%l2)
     do nt = 1,ntdm+1
         write(20,'(8e16.8)') gk(nk,1:latt%nsub,1:latt%nsub,nt)
     enddo
  enddo
  close(20)
  deallocate (gk)
end subroutine fourier_trans_tau_22

subroutine fourier_trans_tau(gr,ftag)
#ifdef _OPENMP
  USE OMP_LIB
#endif
  implicit none
  complex(dp), dimension(:,:) :: gr
  character (40) :: ftag

  ! local
  integer :: imj, nt, nk, no_i, no_j
  character (80) :: filek, filer
  real(dp) :: xk_p(2), aimj_p(2)
  complex(dp), allocatable , dimension(:,:) :: gk

  filek = trim(ftag)//trim('_ktau.bin')
  filer = trim(ftag)//trim('_rtau.bin')

  allocate (gk(lq,ntdm+1))

  gk = dcmplx(0.d0,0.d0)
  open (unit=20,file=filer,status='unknown', action="write", position="append")
  do imj = 1,lq
     do nt = 1,ntdm+1
        do nk = 1,lq
           gk(nk,nt) = gk(nk,nt) +  gr(imj,nt)/latt%zexpiqr(nk, imj)
        enddo
        write(20,'(8e16.8)') gr(imj,nt)
     enddo
  enddo
  close(20)

  gk = gk/dcmplx(dble(lq),0.d0)

  open (unit=20,file=filek,status='unknown', action="write", position="append")
  do nk = 1,lq
     write(20,'(2f16.8)') latt%listk(nk,1)/dble(latt%l1), latt%listk(nk,2)/dble(latt%l2)
     do nt = 1,ntdm+1
         write(20,'(8e16.8)') gk(nk,nt)
     enddo
  enddo
  close(20)
  deallocate (gk)
end subroutine fourier_trans_tau
