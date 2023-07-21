  subroutine equaltime_measure(nt, gf, gfc)
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    implicit none
    integer,intent(in) :: nt
    type(gfunc) :: gf, gfc

    ! local 
    integer :: i, j
    complex(dp) :: zne

    nobs = nobs + 1

    if( dble(phase) < 0.d0 ) then
        sgn = -1.d0
    else
        sgn = 1.d0
    end if
    zphi = dcmplx(sgn/dble(phase), 0.d0) * phase

    !grup (i,j) = < c_i c^+_j >
    !grupc (i,j) = < c^+_i c_j >

    ! get grupc
    do i = 1, ndim
        do j = 1, ndim
            gfc%blk1(j,i) = Imat(i,j) - gf%blk1(i,j)
        end do
    end do

    ! sign
    energy_bin(1) = energy_bin(1) + dcmplx(sgn, 0.d0)

    ! zne
    zne = czero
    do i = 1, ndim
        zne = zne + gfc%blk1(i,i)
    end do
    energy_bin(2) = energy_bin(2) + zne*zphi/dcmplx(dble(lq),0.d0)
    energy_bin(3) = energy_bin(3) + zkint(gf,gfc)*zphi
    energy_bin(4) = energy_bin(4) + zint_u(gf,gfc)*zphi
    call measure_cpcm(gf,gfc)
    zcpcm_bin = zcpcm_bin + zcpcm*zphi
    call measure_spsm(gf,gfc)
    zspsm_bin = zspsm_bin + zspsm*zphi
    call measure_nn(gf,gfc)
    znn_bin = znn_bin + znn*zphi
    zn_bin = zn_bin + zn*zphi
    call measure_bondcorr(gf,gfc)
    zbb_bin = zbb_bin + zbb*zphi
    zb_bin = zb_bin + zb*zphi

  end subroutine equaltime_measure

    complex(dp) function zkint(gf,gfc)
      ! measure kinetic energy
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i, i_0, nf, i_n
      zkint = czero
      ! nn
      do i = 1, latt%nn_lf
          i_0 = latt%nnlf_list(i)
          do nf = 1, latt%nn_nf
              i_n = latt%nnlist(i_0,nf)
              zkint = zkint + dcmplx(-2.d0*dble( gfc%blk1(i_0,i_n)*expar(i_0,nf,xmag,flux_x,flux_y,dimer) ), 0.d0 )
          end do
      end do
    end function zkint

    
    complex(dp) function zint_u(gf,gfc)
      ! onsite u
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i
      zint_u = czero
      do i = 1,  latt%nsites
          zint_u = zint_u + gfc%blk1(i,i)*gfc%blk1(i,i)
      end do
    end function zint_u

    subroutine measure_cpcm(gf,gfc)
      ! <c+ c->
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj
      zcpcm = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              zcpcm(imj) = zcpcm(imj) + gfc%blk1(i,j)
          end do
      end do
    end subroutine measure_cpcm

    subroutine measure_spsm(gf,gfc)
      ! <S+ S->
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj
      zspsm = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              zspsm(imj) = zspsm(imj) + gfc%blk1(i,j)*gf%blk1(i,j)*dcmplx(1.d0-1.d0/dble(nflr*nflr), 0.d0)
          end do
      end do
    end subroutine measure_spsm

    subroutine measure_nn(gf,gfc)
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj
      znn = czero
      zn = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          zn = zn + gfc%blk1(j,j)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              znn(imj) = znn(imj) + gfc%blk1(i,i)*gfc%blk1(j,j) + dcmplx(1.d0/dble(nflr),0.d0)*gfc%blk1(i,j)*gf%blk1(i,j)
          end do
      end do
    end subroutine measure_nn

    subroutine measure_bondcorr(gf,gfc)
      implicit none 
      type(gfunc), intent(in) :: gf, gfc
      integer :: nu_j, j_0, j_n, nu_i, i_0, i_n, imj, nf, fid
      character(40) :: filek
      complex(dp) :: zi, zj, ztmp
      zbb(:) = czero
      zb = czero
      do nu_j = 1, lq
#IFDEF HONEYCOMB
          j_0 = 2*nu_j - 1 ! A site
          j_n = latt%nnlist(j_0,1) ! B site
#ELIF SQUARE
          j_0 = nu_j  ! A site
          j_n = latt%nnlist(j_0,1) ! B site
#ELSE
          j_0 = nu_j  ! A site
          j_n = latt%nnlist(j_0,1) ! B site
#ENDIF
          zj = expar(j_0,1,xmag,flux_x,flux_y,dimer)
          zb = zb + zj*gfc%blk1(j_0,j_n) + dconjg(zj)*gfc%blk1(j_n,j_0)
          do nu_i = 1, lq
#IFDEF HONEYCOMB
              i_0 = 2*nu_i - 1 ! A site
              i_n = latt%nnlist(i_0,1) ! B site
#ELIF SQUARE
              i_0 = nu_i ! A site
              i_n = latt%nnlist(i_0,1) ! B site
#ELSE
              i_0 = nu_i ! A site
              i_n = latt%nnlist(i_0,1) ! B site
#ENDIF
              zi = expar(i_0,1,xmag,flux_x,flux_y,dimer)
              imj = latt%imj(nu_i,nu_j)
              zbb(imj) = zbb(imj) + (zj*gfc%blk1(j_0,j_n) + dconjg(zj)*gfc%blk1(j_n,j_0))*(zi*gfc%blk1(i_0,i_n) + dconjg(zi)*gfc%blk1(i_n,i_0)) &
                                            + ( zi*zj*gfc%blk1(i_0,j_n)*gf%blk1(i_n,j_0) + dconjg(zi)*dconjg(zj)*gfc%blk1(i_n,j_0)*gf%blk1(i_0,j_n) &
                                            +   zi*dconjg(zj)*gfc%blk1(i_0,j_0)*gf%blk1(i_n,j_n) + dconjg(zi)*zj*gfc%blk1(i_n,j_n)*gf%blk1(i_0,j_0) )*dcmplx(1.d0/dble(nflr),0.d0)
          end do
      end do

    end subroutine measure_bondcorr
