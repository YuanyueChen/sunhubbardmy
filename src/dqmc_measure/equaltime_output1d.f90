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
  
    call mpi_reduce( zcpcm_bin, zcpcm, size(zcpcm), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( zspsm_bin, zspsm, size(zspsm), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( znn_bin, znn, size(znn), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( zn_bin, zn, size(zn), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( zjj_bin, zjj, size(zjj), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( zj_bin, zj, 1, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( zbb_bin, zbb, size(zbb), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( zb_bin, zb, 1, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( pair_onsite_bin, pair_onsite, size(pair_onsite), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( pair_nn_bin, pair_nn, size(pair_nn), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( pair_sn_bin, pair_sn, size(pair_sn), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    call mpi_reduce( energy_bin, energy_bin_recv, size(energy_bin_recv), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    
    energy_bin(:) = energy_bin_recv(:)
  
    if( irank .eq. 0 ) then
        rnorm = dble( isize*nobs )
        znorm = dcmplx( rnorm, 0.d0 )
  
        zcpcm = zcpcm / znorm
        ftag = "cpcm"
        call fourier_trans_eqt_22(zcpcm, ftag)
  
        zspsm = zspsm / znorm
        ftag = "spsm"
        call fourier_trans_eqt_22(zspsm, ftag)
  
        znn = znn / znorm
        ftag = "nn"
        call fourier_trans_eqt_22(znn, ftag)
  
        zjj = zjj / znorm
        ftag = "jj"
        call fourier_trans_eqt(zjj, ftag)
  
        zbb = zbb / znorm
        ftag = "bb"
        call fourier_trans_eqt(zbb, ftag)
  
        pair_onsite = pair_onsite / znorm
        ftag = "pair_onsite"
        call fourier_trans_eqt_22(pair_onsite, ftag)
  
        pair_nn = pair_nn / znorm
        ftag = "pair_nn"
        call fourier_trans_eqt(pair_nn, ftag)
  
        pair_sn = pair_sn / znorm
        ftag = "pair_sn"
        call fourier_trans_eqt(pair_sn, ftag)
  
        ! calculate energy_bin
        energy_bin(:) = energy_bin(:) / rnorm
        open (unit=90,file='energy.bin',status='unknown', action="write", position="append")
        write(90, '(20e16.8)') energy_bin(:)
        close(90)
  
        zn = zn / znorm
        zj = zj / znorm
        zb = zb / znorm
        open (unit=90,file='n.bin',status='unknown', action="write", position="append")
        write(90, '(4e16.8)') zn
        close(90)
        open (unit=90,file='j.bin',status='unknown', action="write", position="append")
        write(90, '(2e16.8)') zj
        close(90)
        open (unit=90,file='b.bin',status='unknown', action="write", position="append")
        write(90, '(2e16.8)') zb
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
       write(20,'(10e16.8)') latt%listk(nk,1)/dble(latt%l1), gk(nk)
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
       write(20,'(10e16.8)') latt%listk(nk,1)/dble(latt%l1), gk(nk,1:latt%nsub,1:latt%nsub)
    enddo
    close(20)
  
    deallocate (gk)
  end subroutine fourier_trans_eqt_22
  