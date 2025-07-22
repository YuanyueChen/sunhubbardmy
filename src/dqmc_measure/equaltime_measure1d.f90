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
    energy_bin(2) = energy_bin(2) + zne*zphi/dcmplx(dble(lq),0.d0) ! number density per flavor
    energy_bin(3) = energy_bin(3) + zkint(gf,gfc)*zphi ! kinetic energy
    energy_bin(4) = energy_bin(4) + zint_u(gf,gfc)*zphi ! double occupancy
    energy_bin(5) = energy_bin(5) + zint_v(gf,gfc)*zphi ! nn interaction
    call measure_cpcm(gf,gfc)
    zcpcm_bin = zcpcm_bin + zcpcm*zphi
    call measure_spsm(gf,gfc)
    zspsm_bin = zspsm_bin + zspsm*zphi
    call measure_nn(gf,gfc)
    znn_bin = znn_bin + znn*zphi
    zn_bin = zn_bin + zn*zphi
    call measure_jj(gf,gfc)
    zjj_bin = zjj_bin + zjj*zphi
    zj_bin = zj_bin + zj*zphi
    call measure_bondcorr(gf,gfc)
    zbb_bin = zbb_bin + zbb*zphi
    zb_bin = zb_bin + zb*zphi
    call measure_paircorr(gf,gfc)
    pair_onsite_bin = pair_onsite_bin + pair_onsite*zphi
    pair_nn_bin = pair_nn_bin + pair_nn*zphi
    pair_sn_bin = pair_sn_bin + pair_sn*zphi

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
              ! nflr flavor
              zkint = zkint + dcmplx(-2*nflr*rt*dble( gfc%blk1(i_0,i_n)*expar(i_0,nf,xmag,flux_x,flux_y,dimer) ), 0.d0 )
          end do
      end do
      ! second nearest hopping
      if( t2 .ne. 0.d0 ) then
      do i_0 = 1, latt%nsites
        do nf = 1, 2
          i_n = latt%snlist(i_0,nf)
          ! 0.5 factor to avoid double counting * 2 flavor
          zkint = zkint + dcmplx(-2.d0*t2*dble( gfc%blk1(i_0,i_n) ), 0.d0 )
        end do
      end do
      end if
      ! third nearest hopping
      if( t3 .ne. 0.d0 ) then
      do i_0 = 1, latt%nsites
        do nf = 1, 2
          i_n = latt%tnlist(i_0,nf)
          ! 0.5 factor to avoid double counting * 2 flavor
          zkint = zkint + dcmplx(-2.d0*t3*dble( gfc%blk1(i_0,i_n) ), 0.d0 )
        end do
      end do
      end if
    end function zkint

    
    complex(dp) function zint_u(gf,gfc)
      ! double occupancy
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i
      zint_u = czero
      do i = 1,  latt%nsites
          zint_u = zint_u + gfc%blk1(i,i)*gfc%blk1(i,i)
      end do
    end function zint_u

    complex(dp) function zHsq(gf,gfc)
      ! connected part of <H^2>
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i, j, i_0, j_0, i_n, j_n, nfi, nfj
      zHsq = czero
      do j = 1, latt%nn_lf
          do nfj = 1, latt%nn_nf
              j_0 = latt%nnlf_list(j)
              j_n = latt%nnlist(j_0, nfj)
              do i = 1, latt%nn_lf
                  do nfi = 1, latt%nn_nf
                      i_0 = latt%nnlf_list(i)
                      i_n = latt%nnlist(i_0, nfi)
                      zHsq = zHsq + dcmplx(2.d0*rt*rt,0.d0)*( gfc%blk1(i_0,j_n)*gf%blk1(i_n,j_0) + &
                                                              gfc%blk1(i_0,j_0)*gf%blk1(i_n,j_n) + &
                                                              gfc%blk1(i_n,j_n)*gf%blk1(i_0,j_0) + &
                                                              gfc%blk1(i_n,j_0)*gf%blk1(i_0,j_n) )
                  end do
              end do
          end do
      end do

      do i = 1, latt%nsites
          do j = 1, latt%nn_lf
              do nfj = 1, latt%nn_nf
                  j_0 = latt%nnlf_list(j)
                  j_n = latt%nnlist(j_0, nfj)
                  zHsq = zHsq - dcmplx(2.d0*rt*rhub,0.d0)*gfc%blk1(i,i)*( gfc%blk1(i,j_n)*gf%blk1(i,j_0) + &
                                                                          gfc%blk1(j_0,i)*gf%blk1(j_n,i) + &
                                                                          gfc%blk1(i,j_0)*gf%blk1(i,j_n) + &
                                                                          gfc%blk1(j_n,i)*gf%blk1(j_0,i) )
              end do
          end do
      end do

      do j = 1, latt%nsites
          do i = 1, latt%nsites
              zHsq = zHsq + dcmplx(rhub*rhub,0.d0)*( ( gfc%blk1(i,j)*gf%blk1(i,j) + gfc%blk1(i,i)*gfc%blk1(j,j) )**2 - &
                                                     ( gfc%blk1(i,i)*gfc%blk1(j,j) )**2 )
          end do
      end do
    end function zHsq

    complex(dp) function zint_v(gf,gfc)
      ! measure nn interaction
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i, i_0, nf, i_n
      zint_v = czero
      ! nn
      do i = 1, latt%nn_lf
          i_0 = latt%nnlf_list(i)
          do nf = 1, latt%nn_nf
              i_n = latt%nnlist(i_0,nf)
              zint_v = zint_v + dcmplx(dble( gfc%blk1(i_0,i_n)), 0.d0 )
          end do
      end do
    end function zint_v

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
              zcpcm(imj,no_i,no_j) = zcpcm(imj,no_i,no_j) + gfc%blk1(i,j)
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
              zspsm(imj,no_i,no_j) = zspsm(imj,no_i,no_j) + gfc%blk1(i,j)*gf%blk1(i,j)*dcmplx(1.d0-1.d0/dble(nflr*nflr), 0.d0)
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
          zn(no_j) = zn(no_j) + gfc%blk1(j,j)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              znn(imj,no_i,no_j) = znn(imj,no_i,no_j) + gfc%blk1(i,i)*gfc%blk1(j,j) + dcmplx(1.d0/dble(nflr),0.d0)*gfc%blk1(i,j)*gf%blk1(i,j)
          end do
      end do
    end subroutine measure_nn

    subroutine measure_jj(gf,gfc)
      ! current current correlation
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i, imj, i_0, i_n, j, j_0, j_n
      complex(dp) :: z1, z2, z3, z4
      ! current_i = i * ( c_i^+ c_i+d - c_i+d^+ c_i )
      !   < c_i^+ c_i+d c_j^+ c_j+d >
      ! = < c_i^+ c_i+d > < c_j^+ c_j+d > + < c_i^+ c_j+d > < c_i+d c_j^+ >
      !
      ! < current_i * current_j > = -1 * < ( c_i^+ c_i+d - c_i+d^+ c_i ) * ( c_j^+ c_j+d - c_j+d^+ c_j ) >
      zjj = czero
      do j = 1, lq
          j_0 = j
          j_n = latt%nnlist(j_0,1)
          do i = 1, lq
              imj = latt%imj(i,j)
              i_0 = i
              i_n = latt%nnlist(i_0,1)
              ! disconneted part
              z1 = ( gfc%blk1(i_0,i_n) - gfc%blk1(i_n,i_0) ) * ( gfc%blk1(j_0,j_n) - gfc%blk1(j_n,j_0) )
              ! conneted part
              z3 = gfc%blk1(i_0,j_n)*gf%blk1(i_n,j_0) + gfc%blk1(i_n,j_0)*gf%blk1(i_0,j_n) - gfc%blk1(i_0,j_0)*gf%blk1(i_n,j_n) - gfc%blk1(i_n,j_n)*gf%blk1(i_0,j_0)
              zjj(imj) = zjj(imj) - z1 - dcmplx(1.d0/dble(nflr),0.d0)*z3
          end do
      end do

      ! background
      zj = czero
      do i = 1, lq
          i_0 = i
          i_n = latt%nnlist(i_0,1)
          zj = zj + czi*( gfc%blk1(i_0,i_n) - gfc%blk1(i_n,i_0) )
      end do
    end subroutine measure_jj

    subroutine measure_bondcorr(gf,gfc)
      implicit none 
      type(gfunc), intent(in) :: gf, gfc
      integer :: nu_j, j_0, j_n, nu_i, i_0, i_n, imj, nf, fid
      character(40) :: filek
      complex(dp) :: zi, zj, ztmp
      zbb(:) = czero
      zb = czero
      do nu_j = 1, lq
          j_0 = nu_j  ! A site
          j_n = latt%nnlist(j_0,1) ! B site
          zj = expar(j_0,1,xmag,flux_x,flux_y,dimer)
          zb = zb + zj*gfc%blk1(j_0,j_n) + dconjg(zj)*gfc%blk1(j_n,j_0)
          do nu_i = 1, lq
              i_0 = nu_i ! A site
              i_n = latt%nnlist(i_0,1) ! B site
              zi = expar(i_0,1,xmag,flux_x,flux_y,dimer)
              imj = latt%imj(nu_i,nu_j)
              zbb(imj) = zbb(imj) + (zj*gfc%blk1(j_0,j_n) + dconjg(zj)*gfc%blk1(j_n,j_0))*(zi*gfc%blk1(i_0,i_n) + dconjg(zi)*gfc%blk1(i_n,i_0)) &
                                            + ( zi*zj*gfc%blk1(i_0,j_n)*gf%blk1(i_n,j_0) + dconjg(zi)*dconjg(zj)*gfc%blk1(i_n,j_0)*gf%blk1(i_0,j_n) &
                                            +   zi*dconjg(zj)*gfc%blk1(i_0,j_0)*gf%blk1(i_n,j_n) + dconjg(zi)*zj*gfc%blk1(i_n,j_n)*gf%blk1(i_0,j_0) )*dcmplx(1.d0/dble(nflr),0.d0)
          end do
      end do
    end subroutine measure_bondcorr

    subroutine measure_paircorr(gf,gfc)
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj, j_0, j_n, i_0, i_n
      pair_onsite = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              pair_onsite(imj,no_i,no_j) = pair_onsite(imj,no_i,no_j) + gf%blk1(i,j)*gf%blk1(i,j)
          end do
      end do

      pair_nn = czero
      do nu_j = 1, lq
          j_0 = nu_j  ! A site
          j_n = latt%nnlist(j_0,1) ! B site
          do nu_i = 1, lq
              i_0 = nu_i  ! A site
              i_n = latt%nnlist(i_0,1) ! B site
              imj = latt%imj(nu_i,nu_j)
              pair_nn(imj) = pair_nn(imj) + gf%blk1(i_0,j_0)*gf%blk1(i_n,j_n)
          end do
      end do

      pair_sn = czero
      do nu_j = 1, lq
          j_0 = nu_j  ! A site
          j_n = latt%snlist(j_0,1) ! A site
          do nu_i = 1, lq
              i_0 = nu_i  ! A site
              i_n = latt%snlist(i_0,1) ! A site
              imj = latt%imj(nu_i,nu_j)
              pair_sn(imj) = pair_sn(imj) + gf%blk1(i_0,j_0)*gf%blk1(i_n,j_n)
          end do
      end do
    end subroutine measure_paircorr
