subroutine dqmc_proj_update(this, ntau, nf, ul, ur, ulrinv)
!! update the configuration using submatrix LR algorithm
!! intemediate vectors are preallocated through an oversize array
!! the matrix-vector multiplications are collected in group of "stocks" and done by calling zgemm
!! update is controlled by nustock
    use spring
    use model_para
    use dqmc_basic_data
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    implicit none

    ! arguments
    class(vconf), intent(inout) :: this
    integer, intent(in) :: ntau, nf
    type(gfunc), intent(inout) :: ul, ur, ulrinv

    ! local
    complex(dp) ::  ratio1, ratiotot, delta
    integer :: iflip, is, isp, i, ibond, ibond2, i1, i2, i3, i4, ik, m, n, x, x1
    real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

    complex(dp), allocatable, dimension(:,:) :: vukmat

    integer, allocatable, dimension(:) :: xvec, svec
    complex(dp), allocatable, dimension(:) :: Dvec
    complex(dp), allocatable, dimension(:,:) :: Gammainv
    complex(dp), allocatable, dimension(:,:) :: Rrow, Rtmp, LDcol ! intermediate matrices when calculating PFDP
    complex(dp), allocatable, dimension(:,:) :: PFP_kk, PFP_ikk, PFP_kik ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: ikktmp ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: kiktmp, Gammainv_kk ! for update Gammainv
    type(gfunc) :: Umat, Utmp, Vmat, Vtmp, VUmat, VUtmp ! for update (LR)^-1
    integer :: nustkup, nstk, istk_s, istk_e, i1_stk, i2_stk
    complex(dp), allocatable, dimension(:,:) :: Rrows, Felms ! stock arrays for preventing zgemv
    type(gfunc) :: ulP, urP ! ul, ur after permutation

#IFDEF TIMING
    real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
    call cpu_time_now(starttime)
#ENDIF

#IFDEF TEST
    write(fout, '(a,2e16.8)') ' in update_v using submatrix LR2-3'
#ENDIF

    ! for delay update (LR)^-1
    call allocate_gfunc( Umat, ne, nustock )
    call allocate_gfunc( Utmp, ne, nustock )
    call allocate_gfunc( Vmat, nustock, ne )
    call allocate_gfunc( Vtmp, nustock, ne )
    call allocate_gfunc( VUmat, nustock, nustock )
    call allocate_gfunc( VUtmp, nustock, nustock )

    ! for calculate update ratio
    allocate( vukmat(2,2) )           ! I + Vcal*D
    allocate( xvec(nustock) )         ! x^{(i)}_{j}, j=(i-1)*k+1,...,i*k, update index in ndim arrays
    allocate( svec(nustock) )         ! s^{(i)}_{j}, j=(i-1)*k+1,...,i*k, update index in nustock arrays
    allocate( Dvec(nustock))          ! Delta_{i}(j), i=1,...,nustock/k, j=1,...,k
    allocate( PFP_kk(2,2) )           ! P^{(i+1)}_{k*N}*F*P^{(i+1)}_{N*k}
    allocate( Gammainv(nustock,nustock) )
    allocate( Gammainv_kk(2,2))       ! (I + Vcal*D)^{-1}

    ! matrices for delay update
    xvec = 0
    svec = 0
    Dvec = czero
    Gammainv = czero
    Vmat%blk1 = czero
    Umat%blk1 = czero

    ! preallocate oversize arrays for intermediate matrices
    allocate( kiktmp(2,nustock) )
    allocate( ikktmp(nustock,2) )
    allocate( PFP_ikk(nustock,2) )
    allocate( PFP_kik(2,nustock) )
    kiktmp = czero
    ikktmp = czero
    PFP_ikk = czero
    PFP_kik = czero

    ! reorder the ul and ur to make the sequential update continuous in memory
    ! urP = P*ur, ulP = ul*P^T, ulrinv is not changed
    call allocate_gfunc( ulP, ne, nustock )
    call allocate_gfunc( urP, nustock, ne )

    ! allocate stock arrays for preventing zgemv
    allocate( Rrows(nustock, ne) )      ! rows of R(LR)^-1
    allocate( Felms(nustock, nustock) ) ! selected rows and columns of F

    accm  = 0.d0 ! acceptance rate
    ik = 0       ! delay step * k
    istk_s = 0   ! starting bond index in the stock arrays
    istk_e = 0   ! ending bond index in the stock arrays
    nstk = 0     ! number of bonds in the stock arrays

    do ibond = 1, latt%nn_lf
        i1 = latt%nnlf_list(ibond)
        i2 = latt%nnlist(i1,nf)
        is = this%conf(ibond, nf, ntau)
        iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
        isp = mod(is+iflip-1,this%lcomp) + 1
        
        !! stock up intermediate matrices
        if (nstk == 0) then
            nustkup = min(nustock, (latt%nn_lf - ibond + 1)*2)
#ifdef TEST
            write(fout, '(a,i4,a,i4)') ' in update_v, stock up intermediate matrices from ibond', ibond, ' to ibond', ibond+nustkup/2-1
#endif
            ! reset ulP and urP
            ulP%blk1 = czero
            urP%blk1 = czero
            i = 0
            do ibond2 = ibond, ibond+nustkup/2-1
                i3 = latt%nnlf_list(ibond2)
                i4 = latt%nnlist(i3,nf)
                ulP%blk1(:,i+1) = ul%blk1(:,i3)
                ulP%blk1(:,i+2) = ul%blk1(:,i4)
                urP%blk1(i+1,:) = ur%blk1(i3,:)
                urP%blk1(i+2,:) = ur%blk1(i4,:)
                i = i + 2
            end do

            ! Rrows = R(LR)^-1
            call zgemm('N', 'N', nustkup, ne, ne, cone, urP%blk1, nustock, ulrinv%blk1, ne, czero, Rrows, nustock)
            ! Felms = R(LR)^-1*L
            call zgemm('N', 'N', nustkup, nustkup, ne, cone, Rrows, nustock, ulP%blk1, ne, czero, Felms, nustock)
            istk_s = ibond
            istk_e = ibond + nustkup/2 - 1
            nstk = nustkup/2
        end if

        ! index in the stock arrays
        i1_stk = (ibond - istk_s)*2 + 1
        i2_stk = i1_stk + 1

        !! obtain new blocks of PFP = P^{(i+1)}_{(i+1)k*N}*F*P^{(i+1)}_{N*(i+1)k}
        !!! obtain PFP_kk
        ! Rrow = P^{(i+1)}_{k*N}*R*(LR)^-1
        ! Lcol = L*P^{(i+1)}_{N*k}
        ! PFP_kk = Rrow * Lcol = P^{(i+1)}_{k*N}*F*P^{(i+1)}_{N*k}
        PFP_kk(:,:) = Felms(i1_stk:i2_stk,i1_stk:i2_stk)
        !
        if ( ik > 0 ) then
            !!! obtain PFP_kik
            ! Umat%blk1(:,1:ik) = L*P^{(i)}_{N*ik}
            ! PFP_kik = Rrow * Umat%blk1(:,1:ik) = P^{(i+1)}_{k*N}*F*P^{(i)}_{N*ik}
            PFP_kik(:,1:ik) = Felms(i1_stk:i2_stk, svec(1:ik))
            !
            !!! obtain PFP_ikk
            ! Vmat%blk1(1:ik,:) = P^{(i)}_{ik*N}*R*ulrinv
            ! PFP_ikk = Vmat%blk1(1:ik,:) * Lcol = P^{(i)}_{ik*N}*F*P^{(i+1)}_{N*k}
            PFP_ikk(1:ik,:) = Felms(svec(1:ik), i1_stk:i2_stk)
        end if

        !! vukmat = I + Vcal*D
        vukmat = czero
        ! vukmat = PFP_kik * Gammainv * PFP_ikk
        if ( ik > 0 ) then
            ! ikktmp = Gammainv(1:ik,1:ik)*PFP_ikk
            call zgemm('N', 'N', ik, 2, ik, cone, Gammainv, nustock, PFP_ikk, nustock, czero, ikktmp, nustock)
            ! vukmat = PFP_kik * ikktmp
            call zgemm('N', 'N', 2, 2, ik, cone, PFP_kik, 2, ikktmp, nustock, cone, vukmat, 2)
        end if
        ! vukmat = Vcal = PFP_kk - PFP_kik * Gammainv * PFP_ikk
        vukmat = PFP_kk - vukmat
        ! vukmat = I + Vcal*D
        vukmat(1,1) = cone + vukmat(1,1)*this%delta_bmat_p%blk1(iflip,is)
        vukmat(2,1) =        vukmat(2,1)*this%delta_bmat_p%blk1(iflip,is)
        vukmat(1,2) =        vukmat(1,2)*this%delta_bmat_m%blk1(iflip,is)
        vukmat(2,2) = cone + vukmat(2,2)*this%delta_bmat_m%blk1(iflip,is)

        ratio1 = vukmat(1,1)*vukmat(2,2) - vukmat(1,2)*vukmat(2,1)
        ratiotot = ratio1
#IFDEF TEST
        write(fout, '(a,2e16.8)') ' in update_v, ratiotot = ', ratiotot
#ENDIF
        ratio_re = dble( ratiotot * phase )/dble( phase )
        ratio_re_abs = ratio_re
        if ( ratio_re .lt. 0.d0 ) ratio_re_abs = - ratio_re

        random = spring_sfmt_stream()
        if ( ratio_re_abs.gt.random ) then

            accm  = accm + 1.d0
            !weight_track = weight_track + log( ratio_re_abs )
            weight = dsqrt(dble(ratiotot*dconjg(ratiotot)))
            phase =  phase*ratiotot/dcmplx(weight,0.d0)
#IFDEF TEST
            write(fout, '(a,2e16.8)') ' in update_v, phase = ', phase
#ENDIF

            ! update Gammainv: from (ik*ik) Gammainv^{(i)} to ((ik+2)*(ik+2)) Gammainv^{(i+1)} 
            !! obtain Gammainv_kk = (D^{(i)}_{k*k})(I + Vcal*D)^-1 = (D^{(i)}_{k*k})*vukmat^-1
            Gammainv_kk(1,1) =  vukmat(2,2)*this%delta_bmat_p%blk1(iflip,is)
            Gammainv_kk(1,2) = -vukmat(1,2)*this%delta_bmat_p%blk1(iflip,is)
            Gammainv_kk(2,1) = -vukmat(2,1)*this%delta_bmat_m%blk1(iflip,is)
            Gammainv_kk(2,2) =  vukmat(1,1)*this%delta_bmat_m%blk1(iflip,is)
            Gammainv_kk = Gammainv_kk/ratio1
            ! update 22 block of Gammainv
            Gammainv(ik+1:ik+2,ik+1:ik+2) = Gammainv_kk(:,:)
            if ( ik > 0 ) then
                !! obtain Gammainv_kik = - Gammainv_kk * PFP_kik * Gammainv^{(i)}
                ! kiktmp = PFP_kik * Gammainv^{(i)}
                call zgemm('N', 'N', 2, ik, ik, cone, PFP_kik, 2, Gammainv, nustock, czero, kiktmp, 2)
                ! update 21 block of Gammainv
                call zgemm('N', 'N', 2, ik, 2, -cone, Gammainv_kk, 2, kiktmp, 2, czero, Gammainv(ik+1,1), nustock)

                !! obtain Gammainv_ikk = - Gammainv^{(i)} * PFP_ikk * Gammainv_kk
                ! we already have ikktmp = Gammainv^{(i)} * PFP_ikk
                ! update 12 block of Gammainv
                call zgemm('N', 'N', ik, 2, 2, -cone, ikktmp, nustock, Gammainv_kk, 2, czero, Gammainv(1,ik+1), nustock)

                !! obtain Gammainv_ik = Gammainv^{(i)} - Gammainv^{(i)} * PFDP_ikk * Gammainv_kik
                ! we already have ikktmp = Gammainv^{(i)} * PFP_ikk
                ! we already have iktmp = Gammainv^{(i)}
                ! Gammainv_kik is sliced from Gammainv_kik
                ! update 11 block of Gammainv
                call zgemm('N', 'N', ik, ik, 2, -cone, ikktmp, nustock, Gammainv(ik+1,1), nustock, cone, Gammainv, nustock)
            end if

            ik = ik + 2
            ! store xvec and Dvec
            xvec(ik-1) = i1
            xvec(ik) = i2
            svec(ik-1) = i1_stk
            svec(ik) = i2_stk
            Dvec(ik-1) = this%delta_bmat_p%blk1(iflip,is)
            Dvec(ik) = this%delta_bmat_m%blk1(iflip,is)
            ! append Rrow, Lcol to Vmat%blk1, Umat%blk1
            Vmat%blk1(ik-1:ik,:) = Rrows(i1_stk:i2_stk,:)
            Umat%blk1(:,ik-1:ik) = ulP%blk1(:,i1_stk:i2_stk)

            ! flip
            this%conf(ibond,nf,ntau) =  isp
        end if

        if( (ibond == istk_e) .or. (ibond == latt%nn_lf) ) then
#IFDEF TIMING
            call cpu_time_now(starttime11)
#ENDIF
#ifdef TEST
            write(fout, '(a,i4)') ' in update_v, update the whole LRinv'
#endif
            ! delay update: update (LR)^-1 and R
            
            !! update (LR)^-1
            ! (LR)^-1 = (LR)^-1 - (LR)^-1 * Umat * (D^-1 + Vmat*Umat)^-1 * Vmat
            ! we already have Umat=L*P^{(i)}_{N*ik} and Vmat = P^{(i)}_{ik*N}*R*(LR)^-1
            ! we already have Gammainv = (D^-1 + Vmat*Umat)^-1
            ! 
            ! obtain Vtmp = Gammainv * Vmat
            call zgemm('N', 'N', ik, ne, ik, cone, Gammainv, nustock, Vmat%blk1, nustock, czero, Vtmp%blk1, nustock)
            ! obtain Utmp = ulrinv * Umat
            call zgemm('N', 'N', ne, ik, ne, cone, ulrinv%blk1, ne, Umat%blk1, ne, czero, Utmp%blk1, ne)
            ! obtain new ulrinv = ulrinv - (ulrinv * Umat) * (Gammainv * Vmat)
            call zgemm('N', 'N', ne, ne, ik, -cone, Utmp%blk1, ne, Vtmp%blk1, nustock, cone, ulrinv%blk1, ne)
            
            !! update R^{(n)} = (I + Delta^{(i)}) R^{(0)}
            do m = 1, ik
                x = xvec(m)
                delta = Dvec(m)
                ur%blk1(x,:) = ur%blk1(x,:) + delta*ur%blk1(x,:)
            end do

            ik = 0
            istk_s = 0
            istk_e = 0
            nstk = 0
            ! initial xvec, svec, Dvec and Gammainv
            xvec = 0
            svec = 0
            Dvec = czero
            Gammainv = czero
            ! initial Vmat, Umat
            Vmat%blk1 = czero
            Umat%blk1 = czero
            ! initial oversize arrays for intermediate matrices
            kiktmp = czero
            ikktmp = czero
            PFP_ikk = czero
            PFP_kik = czero

#IFDEF TIMING
            call cpu_time_now(endtime11)
            timecalculation(20)=timecalculation(20)+endtime11-starttime11
#ENDIF
        end if
    end do
    main_obs(4) = main_obs(4) + dcmplx( accm, latt%nn_lf )

#IFDEF TIMING
    call cpu_time_now(endtime)
    timecalculation(17)=timecalculation(17)+endtime-starttime
#ENDIF
end subroutine dqmc_proj_update
