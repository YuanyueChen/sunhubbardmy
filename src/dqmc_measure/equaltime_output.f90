subroutine equaltime_output
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  implicit none

  include 'mpif.h'

  integer :: i, j, imj, nf
  real(dp) :: rnorm
  complex(dp) :: znorm
  character(40) :: ftag

  call mpi_reduce( zcpcm_orb1_bin, zcpcm_orb1, size(zcpcm_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( zspsm_orb1_bin, zspsm_orb1, size(zspsm_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( znn_orb1_bin, znn_orb1, size(znn_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( zn_orb1_bin, zn_orb1, size(zn_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( zjj_orb1_bin, zjj_orb1, size(zjj_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( zj_orb1_bin, zj_orb1, 1, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( zbb_orb1_bin, zbb_orb1, size(zbb_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( zb_orb1_bin, zb_orb1, 1, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( pair_onsite_orb1_bin, pair_onsite_orb1, size(pair_onsite_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( pair_nn_orb1_bin, pair_nn_orb1, size(pair_nn_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( pair_sn_orb1_bin, pair_sn_orb1, size(pair_sn_orb1), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  call mpi_reduce( energy_bin, energy_bin_recv, size(energy_bin_recv), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )

  energy_bin(:) = energy_bin_recv(:)

  if( irank .eq. 0 ) then
      rnorm = dble( isize*nobs )
      znorm = dcmplx( rnorm, 0.d0 )

      zcpcm_orb1 = zcpcm_orb1 / znorm
      ftag = "cpcm_orb1"
      call fourier_trans_eqt_22(zcpcm_orb1, ftag)

      zspsm_orb1 = zspsm_orb1 / znorm
      ftag = "spsm_orb1"
      call fourier_trans_eqt_22(zspsm_orb1, ftag)

      znn_orb1 = znn_orb1 / znorm
      ftag = "nn_orb1"
      call fourier_trans_eqt_22(znn_orb1, ftag)

      zjj_orb1 = zjj_orb1 / znorm
      ftag = "jj_orb1"
      call fourier_trans_eqt(zjj_orb1, ftag)

      zbb_orb1 = zbb_orb1 / znorm
      ftag = "bb_orb1"
      call fourier_trans_eqt(zbb_orb1, ftag)

      pair_onsite_orb1 = pair_onsite_orb1 / znorm
      ftag = "pair_onsite_orb1"
      call fourier_trans_eqt_22(pair_onsite_orb1, ftag)

      pair_nn_orb1 = pair_nn_orb1 / znorm
      ftag = "pair_nn_orb1"
      call fourier_trans_eqt(pair_nn_orb1, ftag)

      pair_sn_orb1 = pair_sn_orb1 / znorm
      ftag = "pair_sn_orb1"
      call fourier_trans_eqt(pair_sn_orb1, ftag)

      ! calculate energy_bin
      energy_bin(:) = energy_bin(:) / rnorm
      open (unit=90,file='energy.bin',status='unknown', action="write", position="append")
      write(90, '(16e16.8)') energy_bin(1:8)
      close(90)

      zn_orb1 = zn_orb1 / znorm
      zj_orb1 = zj_orb1 / znorm
      zb_orb1 = zb_orb1 / znorm
      open (unit=90,file='n.bin',status='unknown', action="write", position="append")
      write(90, '(4e16.8)') zn_orb1
      close(90)
      open (unit=90,file='j.bin',status='unknown', action="write", position="append")
      write(90, '(2e16.8)') zj_orb1
      close(90)
      open (unit=90,file='b.bin',status='unknown', action="write", position="append")
      write(90, '(2e16.8)') zb_orb1
      close(90)

  end if

  call mpi_barrier( mpi_comm_world, ierr )

end subroutine equaltime_output

subroutine fourier_trans_eqt(gr,ftag)
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  implicit none
  complex(dp), dimension(:) :: gr
  character (40) :: ftag

  ! local
  integer :: imj, nk
  real(dp) :: xk_p(2)
  complex(dp), allocatable , dimension(:) :: gk
  character (80) :: filek, filer
  filek = trim(ftag)//trim('_k.bin')
  filer = trim(ftag)//trim('_r.bin')

  allocate (gk(lq))

  open (unit=20,file=filer,status='unknown', action="write", position="append")
  gk = dcmplx(0.d0,0.d0)
  do imj = 1,lq
     do nk = 1,lq
        gk(nk) = gk(nk) +  gr(imj)/latt%zexpiqr(nk,imj)
     enddo
     write(20,'(2e16.8)') gr(imj)
  enddo
  gk = gk/dcmplx(dble(lq),0.d0)

  open (unit=20,file=filek,status='unknown', action="write", position="append")
  do nk = 1,lq
     !xk_p = dble(latt%listk(nk,1))*latt%b1_p + dble(latt%listk(nk,2))*latt%b2_p
     !write(20,'(4e16.8)') xk_p(1), xk_p(2), gk(nk)
     write(20,'(10e16.8)') latt%listk(nk,1)/dble(latt%l1), latt%listk(nk,2)/dble(latt%l2), gk(nk)
  enddo
  close(20)

  deallocate (gk)
end subroutine fourier_trans_eqt

subroutine fourier_trans_eqt_22(gr,ftag)
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  implicit none
  complex(dp), intent(in), dimension(:,:,:) :: gr
  character (40) :: ftag

  ! local
  integer :: imj, nk, nf, nfi, nfj
  real(dp) :: xk_p(2)
  complex(dp), allocatable , dimension(:,:,:) :: gk
  character (80) :: filek, filer
  filek = trim(ftag)//trim('_k.bin')
  filer = trim(ftag)//trim('_r.bin')

  allocate (gk(lq,latt%nsub,latt%nsub))

  open (unit=20,file=filer,status='unknown', action="write", position="append")
  gk = dcmplx(0.d0,0.d0)
  do imj = 1,lq
   do nfj = 1, latt%nsub
    do nfi = 1, latt%nsub
     do nk = 1,lq
       gk(nk,nfi,nfj) = gk(nk,nfi,nfj) +  gr(imj,nfi,nfj)/latt%zexpiqr(nk,imj)
     enddo
    end do
   end do
   write(20,'(10e16.8)') gr(imj,1:latt%nsub,1:latt%nsub)
  end do
  close(20)
  gk = gk/dcmplx(dble(lq),0.d0)

  open (unit=20,file=filek,status='unknown', action="write", position="append")
  do nk = 1,lq
     !xk_p = dble(latt%listk(nk,1))*latt%b1_p + dble(latt%listk(nk,2))*latt%b2_p
     !write(20,'(10e16.8)') xk_p(1), xk_p(2), gk(nk,1:2,1:2)
     write(20,'(10e16.8)') latt%listk(nk,1)/dble(latt%l1), latt%listk(nk,2)/dble(latt%l2), gk(nk,1:latt%nsub,1:latt%nsub)
  enddo
  close(20)

  deallocate (gk)
end subroutine fourier_trans_eqt_22
