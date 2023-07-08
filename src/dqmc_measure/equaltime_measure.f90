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
#IFDEF TIMING
    real(dp) :: starttime, endtime
#ENDIF
#IFDEF TIMING
    call cpu_time_now(starttime)
#ENDIF

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
            gfc%orb1(j,i) = Imat(i,j) - gf%orb1(i,j)
        end do
    end do

    ! sign
    energy_bin(1) = energy_bin(1) + dcmplx(sgn, 0.d0)

    ! zne
    zne = czero
    do i = 1, ndim
        zne = zne + gfc%orb1(i,i)
    end do
    energy_bin(2) = energy_bin(2) + zne*zphi/dcmplx(dble(lq),0.d0)
    energy_bin(3) = energy_bin(3) + zkint(gf,gfc)*zphi
    energy_bin(4) = energy_bin(4) + zint_u(gf,gfc)*zphi
    energy_bin(5) = energy_bin(5) + zqsq(gf,gfc)*dcmplx(1.d0/9.d0,0.d0)*zphi
    energy_bin(6) = energy_bin(6) + ztsq(gf,gfc)*zphi
    energy_bin(7) = energy_bin(7) + ztq(gf,gfc)*zphi
    call measure_cpcm(gf,gfc)
    zcpcm_orb1_bin = zcpcm_orb1_bin + zcpcm_orb1*zphi
    call measure_spsm(gf,gfc)
    zspsm_orb1_bin = zspsm_orb1_bin + zspsm_orb1*zphi
    call measure_nn(gf,gfc)
    znn_orb1_bin = znn_orb1_bin + znn_orb1*zphi
    zn_orb1_bin = zn_orb1_bin + zn_orb1*zphi
    call measure_jj(gf,gfc)
    zjj_orb1_bin = zjj_orb1_bin + zjj_orb1*zphi
    zj_orb1_bin = zj_orb1_bin + zj_orb1*zphi
    call measure_bondcorr(gf,gfc)
    zbb_orb1_bin = zbb_orb1_bin + zbb_orb1*zphi
    zb_orb1_bin = zb_orb1_bin + zb_orb1*zphi
    call measure_paircorr(gf,gfc)
    pair_onsite_orb1_bin = pair_onsite_orb1_bin + pair_onsite_orb1*zphi
    pair_nn_orb1_bin = pair_nn_orb1_bin + pair_nn_orb1*zphi
    pair_sn_orb1_bin = pair_sn_orb1_bin + pair_sn_orb1*zphi
#IFDEF TIMING
    call cpu_time_now(endtime)
    timecalculation(4)=timecalculation(4)+endtime-starttime
#ENDIF

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
              zkint = zkint + dcmplx(-2.d0*dble( gfc%orb1(i_0,i_n)*expar(i_0,nf,xmag,flux_x,flux_y,dimer) ), 0.d0 )
          end do
      end do
      ! second nearest hopping
      if( t2 .ne. 0.d0 ) then
      do i_0 = 1, latt%nsites
#IFDEF SQUARE
        do nf = 1, 4
#ELIF HONEYCOMB
        do nf = 1, 6
#ENDIF
          i_n = latt%snlist(i_0,nf)
          ! 0.5 factor to avoid double counting
          zkint = zkint + dcmplx(-t2*dble( gfc%orb1(i_0,i_n) ), 0.d0 )
        end do
      end do
      end if
      ! third nearest hopping
      if( t3 .ne. 0.d0 ) then
      do i_0 = 1, latt%nsites
#IFDEF SQUARE
        do nf = 1, 4
#ELIF HONEYCOMB
        do nf = 1, 3
#ENDIF
          i_n = latt%tnlist(i_0,nf)
          ! 0.5 factor to avoid double counting
          zkint = zkint + dcmplx(-t3*dble( gfc%orb1(i_0,i_n) ), 0.d0 )
        end do
      end do
      end if
    end function zkint

    
    complex(dp) function zint_u(gf,gfc)
      ! onsite u
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i
      zint_u = czero
      do i = 1,  latt%nsites
          zint_u = zint_u + gfc%orb1(i,i)*gfc%orb1(i,i)
      end do
    end function zint_u

    complex(dp) function zqsq(gf,gfc)
      ! Q^2 term
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i, i1, i2, isite1, isite2
      complex(dp) :: ztmp
      zqsq = czero
      do i = 1, lq
          ! disconnected part
          ztmp = czero
          do isite1 = 1, latt%z_plq
              i1 = latt%plq_cord(isite1,i)
              ztmp = ztmp + gfc%orb1(i1,i1)
          end do
          zqsq = zqsq + ztmp*ztmp*dcmplx(dble(nflr),0.d0) ! the factor comes from nflr^2/nflr

          ! connected part
          do isite1 = 1,  latt%z_plq
              i1 = latt%plq_cord(isite1,i)
              do isite2 = 1,  latt%z_plq
                  i2 = latt%plq_cord(isite2,i)
                  zqsq = zqsq + gfc%orb1(i1,i2)*gf%orb1(i1,i2)
              end do
          end do
      end do
    end function zqsq

    complex(dp) function ztsq(gf,gfc)
      ! T^2 term
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i, i1, i1a, i2, i2a, isite1, isite2
      complex(dp) :: ztmp, z1, z2, z3, z4
      ztsq = czero
      do i = 1, lq
          ! disconnected part
          ztmp = czero
          do isite1 = 1, latt%z_plq
              i1 = latt%plq_cord(isite1,i)
              i1a = latt%plq_cord(npbc(isite1+1,latt%z_plq),i)
              if( mod(isite1,2) == 0 ) then
                  z1 = -exp( dcmplx(0.d0, theta) )
              else
                  z1 =  exp( dcmplx(0.d0,-theta) )
              end if
              ztmp = ztmp + dcmplx( 2.d0*dble( z1*gfc%orb1(i1,i1a) ), 0.d0)
          end do
          ztsq = ztsq + ztmp*ztmp*dcmplx(dble(nflr),0.d0) ! the factor comes from nflr^2/nflr

          ! connected part
          do isite1 = 1, latt%z_plq
              i1 = latt%plq_cord(isite1,i)
              i1a = latt%plq_cord(npbc(isite1+1,latt%z_plq),i)
              if( mod(isite1,2) == 0 ) then
                  z1 = -exp( dcmplx(0.d0, theta) )
              else
                  z1 =  exp( dcmplx(0.d0,-theta) )
              end if
              z2 = dconjg(z1)
              do isite2 = 1, latt%z_plq
                  i2 = latt%plq_cord(isite2,i)
                  i2a = latt%plq_cord(npbc(isite2+1,latt%z_plq),i)
                  if( mod(isite2,2) == 0 ) then
                      z3 = -exp( dcmplx(0.d0, theta) )
                  else
                      z3 =  exp( dcmplx(0.d0,-theta) )
                  end if
                  z4 = dconjg(z3)
                  ztsq = ztsq + z1*z3*( gfc%orb1(i1,i2a)*gf%orb1(i1a,i2) ) + &
                                z1*z4*( gfc%orb1(i1,i2)*gf%orb1(i1a,i2a) ) + &
                                z2*z3*( gfc%orb1(i1a,i2a)*gf%orb1(i1,i2) ) + &
                                z2*z4*( gfc%orb1(i1a,i2)*gf%orb1(i1,i2a) )
              end do
          end do
      end do
    end function ztsq

    complex(dp) function ztq(gf,gfc)
      ! T^2 term
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: i, i1, i1a, i2, i2a, isite1, isite2
      complex(dp) :: ztmp1, ztmp2, z1, z2
      ztq = czero
      do i = 1, lq
          ! disconnected part
          ztmp1 = czero
          ztmp2 = czero
          do isite1 = 1, latt%z_plq
              i1 = latt%plq_cord(isite1,i)
              i1a = latt%plq_cord(npbc(isite1+1,latt%z_plq),i)
              if( mod(isite1,2) == 0 ) then
                  z1 = -exp( dcmplx(0.d0, theta) )
              else
                  z1 =  exp( dcmplx(0.d0,-theta) )
              end if
              z2 = dconjg(z1)
              ztmp1 = ztmp1 + dcmplx( 2.d0*dble( z1*gfc%orb1(i1,i1a) ), 0.d0)
              ztmp2 = ztmp2 + dcmplx(1.d0/3.d0,0.d0)*gfc%orb1(i1,i1)
          end do
          ztq = ztq + ztmp1*ztmp2*dcmplx(dble(nflr),0.d0) ! the factor comes from nflr^2/nflr

          ! connected part
          do isite1 = 1, latt%z_plq
              i1 = latt%plq_cord(isite1,i)
              i1a = latt%plq_cord(npbc(isite1+1,latt%z_plq),i)
              if( mod(isite1,2) == 0 ) then
                  z1 = -exp( dcmplx(0.d0, theta) )/dcmplx(3.d0,0.d0)
              else
                  z1 =  exp( dcmplx(0.d0,-theta) )/dcmplx(3.d0,0.d0)
              end if
              z2 = dconjg(z1)
              do isite2 = 1, latt%z_plq
                  i2 = latt%plq_cord(isite2,i)
                  i2a = latt%plq_cord(npbc(isite2+1,latt%z_plq),i)
                  ztq = ztq + z1*( gfc%orb1(i1,i2)*gf%orb1(i1a,i2) ) + &
                              z2*( gfc%orb1(i1a,i2)*gf%orb1(i1,i2) )
              end do
          end do
      end do
    end function ztq

    subroutine measure_cpcm(gf,gfc)
      ! <c+ c->
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj
      zcpcm_orb1 = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              zcpcm_orb1(imj,no_i,no_j) = zcpcm_orb1(imj,no_i,no_j) + gfc%orb1(i,j)
          end do
      end do
    end subroutine measure_cpcm

    subroutine measure_spsm(gf,gfc)
      ! <S+ S-> for orb1
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj
      zspsm_orb1 = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              zspsm_orb1(imj,no_i,no_j) = zspsm_orb1(imj,no_i,no_j) + gfc%orb1(i,j)*gf%orb1(i,j)*dcmplx(1.d0-1.d0/dble(nflr*nflr), 0.d0)
          end do
      end do
    end subroutine measure_spsm

    subroutine measure_nn(gf,gfc)
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj
      znn_orb1 = czero
      zn_orb1 = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          zn_orb1(no_j) = zn_orb1(no_j) + gfc%orb1(j,j)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              znn_orb1(imj,no_i,no_j) = znn_orb1(imj,no_i,no_j) + gfc%orb1(i,i)*gfc%orb1(j,j) + dcmplx(1.d0/dble(nflr),0.d0)*gfc%orb1(i,j)*gf%orb1(i,j)
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
      zjj_orb1 = czero
      do j = 1, lq
#IFDEF HONEYCOMB
          j_0 = 2*j - 1
          j_n = latt%snlist(j_0,1)
#ELIF SQUARE
          j_0 = j
          j_n = latt%nnlist(j_0,1)
#ENDIF
          do i = 1, lq
              imj = latt%imj(i,j)
#IFDEF HONEYCOMB
              i_0 = 2*i - 1
              i_n = latt%snlist(i_0,1)
#ELIF SQUARE
              i_0 = i
              i_n = latt%nnlist(i_0,1)
#ENDIF
              ! disconneted part
              z1 = ( gfc%orb1(i_0,i_n) - gfc%orb1(i_n,i_0) ) * ( gfc%orb1(j_0,j_n) - gfc%orb1(j_n,j_0) )
              ! conneted part
              z3 = gfc%orb1(i_0,j_n)*gf%orb1(i_n,j_0) + gfc%orb1(i_n,j_0)*gf%orb1(i_0,j_n) - gfc%orb1(i_0,j_0)*gf%orb1(i_n,j_n) - gfc%orb1(i_n,j_n)*gf%orb1(i_0,j_0)
              zjj_orb1(imj) = zjj_orb1(imj) - z1 - dcmplx(1.d0/dble(nflr),0.d0)*z3
          end do
      end do

      ! background
      zj_orb1 = czero
      do i = 1, lq
#IFDEF HONEYCOMB
          i_0 = 2*i - 1
          i_n = latt%snlist(i_0,1)
#ELIF SQUARE
          i_0 = i
          i_n = latt%nnlist(i_0,1)
#ENDIF
          zj_orb1 = zj_orb1 + czi*( gfc%orb1(i_0,i_n) - gfc%orb1(i_n,i_0) )
      end do
    end subroutine measure_jj

    subroutine measure_bondcorr(gf,gfc)
      implicit none 
      type(gfunc), intent(in) :: gf, gfc
      integer :: nu_j, j_0, j_n, nu_i, i_0, i_n, imj, nf, fid
      character(40) :: filek
      complex(dp) :: zi, zj, ztmp
      zbb_orb1(:) = czero
      zb_orb1 = czero
      do nu_j = 1, lq
#IFDEF HONEYCOMB
          j_0 = 2*nu_j - 1 ! A site
          j_n = latt%nnlist(j_0,1) ! B site
#ELIF SQUARE
          j_0 = nu_j  ! A site
          j_n = latt%nnlist(j_0,1) ! B site
#ENDIF
          zj = expar(j_0,1,xmag,flux_x,flux_y,dimer)
          zb_orb1 = zb_orb1 + zj*gfc%orb1(j_0,j_n) + dconjg(zj)*gfc%orb1(j_n,j_0)
          do nu_i = 1, lq
#IFDEF HONEYCOMB
              i_0 = 2*nu_i - 1 ! A site
              i_n = latt%nnlist(i_0,1) ! B site
#ELIF SQUARE
              i_0 = nu_i ! A site
              i_n = latt%nnlist(i_0,1) ! B site
#ENDIF
              zi = expar(i_0,1,xmag,flux_x,flux_y,dimer)
              imj = latt%imj(nu_i,nu_j)
              zbb_orb1(imj) = zbb_orb1(imj) + (zj*gfc%orb1(j_0,j_n) + dconjg(zj)*gfc%orb1(j_n,j_0))*(zi*gfc%orb1(i_0,i_n) + dconjg(zi)*gfc%orb1(i_n,i_0)) &
                                            + ( zi*zj*gfc%orb1(i_0,j_n)*gf%orb1(i_n,j_0) + dconjg(zi)*dconjg(zj)*gfc%orb1(i_n,j_0)*gf%orb1(i_0,j_n) &
                                            +   zi*dconjg(zj)*gfc%orb1(i_0,j_0)*gf%orb1(i_n,j_n) + dconjg(zi)*zj*gfc%orb1(i_n,j_n)*gf%orb1(i_0,j_0) )*dcmplx(1.d0/dble(nflr),0.d0)
          end do
      end do
#IFDEF HIST
      ztmp = czero
      zi = exp(dcmplx(0.d0,2.d0/3.d0*pi))
      do nu_j = 1, lq
#IFDEF HONEYCOMB
          j_0 = 2*nu_j - 1 ! A site
          !do nf = 1, 3
          !    j_n = latt%nnlist(j_0,nf) ! B site
          !    zj = expar(j_0,nf,xmag,flux_x,flux_y,dimer)
          !    ztmp = ztmp + zi**(nf-1)*(zj*gfc%orb1(j_0,j_n) + dconjg(zj)*gfc%orb1(j_n,j_0))*latt%zexpiqr( (latt%l1-1)/2*latt%l2+latt%l1/3*latt%l2+(latt%l2-1)/2+latt%l2/3+1, nu_j )
          !end do
          j_n = latt%nnlist(j_0,1) ! B site
          zj = expar(j_0,1,xmag,flux_x,flux_y,dimer)
          ztmp = ztmp + (zj*gfc%orb1(j_0,j_n) + dconjg(zj)*gfc%orb1(j_n,j_0))*latt%zexpiqr( (latt%l1-1)/2*latt%l2+latt%l1/3*latt%l2+(latt%l2-1)/2+latt%l2/3+1, nu_j )
          !write(*,'(i6,2e16.8)') (latt%l1-1)/2*latt%l2+latt%l1/3*latt%l2+(latt%l2-1)/2+latt%l2/3+1, latt%zexpiqr( (latt%l1-1)/2*latt%l2+latt%l1/3*latt%l2+(latt%l2-1)/2+latt%l2/3+1, nu_j )
#ELIF SQUARE
          j_0 = nu_j ! A site
          j_n = latt%nnlist(j_0,1) ! B site, x-dir
          zj = expar(j_0,1,xmag,flux_x,flux_y,dimer)
          ztmp = ztmp + dcmplx( dble( (zj*gfc%orb1(j_0,j_n) + dconjg(zj)*gfc%orb1(j_n,j_0))*latt%zexpiqr( latt%l1*latt%l2-latt%l2/2, nu_j ) ), 0.d0 ) ! (pi,0)
          j_n = latt%nnlist(j_0,2) ! B site, y-dir
          zj = expar(j_0,2,xmag,flux_x,flux_y,dimer)
          ztmp = ztmp + dcmplx( 0.d0, dble( (zj*gfc%orb1(j_0,j_n) + dconjg(zj)*gfc%orb1(j_n,j_0))*latt%zexpiqr( latt%l1/2*latt%l2, nu_j ) ) ) ! (0,pi)
          !write(*,'(i6,2e16.8,i6,2e16.8)') latt%l1*latt%l2-latt%l2/2, latt%zexpiqr( latt%l1*latt%l2-latt%l2/2, nu_j ), latt%l1/2*latt%l2, latt%zexpiqr( latt%l1/2*latt%l2, nu_j )
#ENDIF
      end do
      fid = 1001+irank
      write(filek,'(a,i4,a)')"bb_hist",fid,".bin"
      open(unit=fid,file=filek,status='unknown', action="write", position="append")
      write(fid,'(2e16.8)') ztmp/dcmplx( dble(lq), 0.d0 )
#ENDIF
    end subroutine measure_bondcorr

    subroutine measure_paircorr(gf,gfc)
      implicit none
      type(gfunc), intent(in) :: gf, gfc
      integer :: j, nu_j, no_j, i, nu_i, no_i, imj, j_0, j_n, i_0, i_n
      pair_onsite_orb1 = czero
      do j = 1, latt%nsites
          nu_j = latt%list(j,1)
          no_j = latt%list(j,2)
          do i = 1, latt%nsites
              nu_i = latt%list(i,1)
              no_i = latt%list(i,2)
              imj  = latt%imj(nu_i,nu_j)
              pair_onsite_orb1(imj,no_i,no_j) = pair_onsite_orb1(imj,no_i,no_j) + gf%orb1(i,j)*gf%orb1(i,j)
          end do
      end do

      pair_nn_orb1 = czero
      do nu_j = 1, lq
#IFDEF HONEYCOMB
          j_0 = 2*nu_j - 1 ! A site
          j_n = latt%nnlist(j_0,1) ! B site
#ELIF SQUARE
          j_0 = nu_j  ! A site
          j_n = latt%nnlist(j_0,1) ! B site
#ENDIF
          do nu_i = 1, lq
#IFDEF HONEYCOMB
              i_0 = 2*nu_i - 1 ! A site
              i_n = latt%nnlist(i_0,1) ! B site
#ELIF SQUARE
              i_0 = nu_i  ! A site
              i_n = latt%nnlist(i_0,1) ! B site
#ENDIF
              imj = latt%imj(nu_i,nu_j)
              pair_nn_orb1(imj) = pair_nn_orb1(imj) + gf%orb1(i_0,j_0)*gf%orb1(i_n,j_n)
          end do
      end do

      pair_sn_orb1 = czero
      do nu_j = 1, lq
#IFDEF HONEYCOMB
          j_0 = 2*nu_j - 1 ! A site
          j_n = latt%snlist(j_0,1) ! A site
#ELIF SQUARE
          j_0 = nu_j  ! A site
          j_n = latt%snlist(j_0,1) ! A site
#ENDIF
          do nu_i = 1, lq
#IFDEF HONEYCOMB
              i_0 = 2*nu_i - 1 ! A site
              i_n = latt%snlist(i_0,1) ! A site
#ELIF SQUARE
              i_0 = nu_i  ! A site
              i_n = latt%snlist(i_0,1) ! A site
#ENDIF
              imj = latt%imj(nu_i,nu_j)
              pair_sn_orb1(imj) = pair_sn_orb1(imj) + gf%orb1(i_0,j_0)*gf%orb1(i_n,j_n)
          end do
      end do
    end subroutine measure_paircorr
