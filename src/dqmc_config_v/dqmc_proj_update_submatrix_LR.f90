subroutine dqmc_proj_update(this, ntau, nf, ul, ur, ulrinv)
!! update the configuration using submatrix LR1-1
!! intemediate vectors are re-allocated when their size changes
!! the matrix-vector multiplications are done by calling zgemm
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
    integer :: iflip, is, isp, i, i1, i2, ik, ik_mat, m, n, x, x1
    real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

    complex(dp), allocatable, dimension(:,:) :: vukmat

    integer, allocatable, dimension(:) :: xvec
    complex(dp), allocatable, dimension(:) :: Dvec
    complex(dp), allocatable, dimension(:,:) :: Gammainv
    complex(dp), allocatable, dimension(:,:) :: Rrow, Rtmp, LDcol ! intermediate matrices when calculating PFDP
    complex(dp), allocatable, dimension(:,:) :: PFDP_kk, PFDP_ikk, PFDP_kik ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: neiktmp, iknetmp ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: iktmp, ikktmp ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: kiktmp, Gammainv_kk ! for update Gammainv
    type(gfunc) :: Umat, Utmp, Vmat, Vtmp, VUmat, VUtmp ! for update (LR)^-1
    type(gfunc) :: Fmmat, Rmmat

#IFDEF TIMING
    real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
    call cpu_time_now(starttime)
#ENDIF

#IFDEF TEST
    write(fout, '(a,2e16.8)') ' in update_v using submatrix LR'
#ENDIF

    ! for delay update (LR)^-1
    call allocate_gfunc( Umat, ne, nublock )
    call allocate_gfunc( Utmp, ne, nublock )
    call allocate_gfunc( Vmat, nublock, ne )
    call allocate_gfunc( Vtmp, nublock, ne )
    call allocate_gfunc( VUmat, nublock, nublock )
    call allocate_gfunc( VUtmp, nublock, nublock )

    ! for calculate update ratio
    allocate( vukmat(2,2) )           ! I + Vcal*D
    allocate( xvec(nublock) )         ! x^{(i)}_{j}, j=(i-1)*k+1,...,i*k
    allocate( Dvec(nublock))          ! Delta_{i}(j), i=1,...,nublock/k, j=1,...,k
    allocate( PFDP_kk(2,2) )          ! P^{(i+1)}_{k*N}*F*Delta*P^{(i+1)}_{N*k}
    allocate( Rrow(2,ne) )            ! P^{(i)}_{k*N}*R
    allocate( Rtmp(2,ne) )
    allocate( LDcol(ne,2) )           ! L*Delta_{i}*P^{(i)}_{N*k}
    allocate( Gammainv(nublock,nublock) )
    allocate( Gammainv_kk(2,2))       ! (I + Vcal*D)^{-1}

#IFDEF TEST
    write(fout, '(a)') ' in update_v'
    call allocate_gfunc( Fmmat, ndim, ndim )
    call allocate_gfunc( Rmmat, ndim, ne )
    ! Rmmat = R * (LR)^-1
    call zgemm('N', 'N', ndim, ne, ne, cone, ur%blk1, ndim, ulrinv%blk1, ne, czero, Rmmat%blk1, ndim)
    ! Fmmat = Rmmat * L
    call zgemm('N', 'N', ndim, ndim, ne, cone, Rmmat%blk1, ndim, ul%blk1, ne, czero, Fmmat%blk1, ndim)
#ENDIF

    ! matrices for delay update
    xvec = 0
    Dvec = czero
    Gammainv = czero
    Vmat%blk1 = czero
    Umat%blk1 = czero

    accm  = 0.d0 ! acceptance rate
    ik = 0       ! delay step * k
    ik_mat = 0   ! size of iktmp

    do i = 1, latt%nn_lf
        i1 = latt%nnlf_list(i)
        i2 = latt%nnlist(i1,nf)
        is = this%conf(i, nf, ntau)
        iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
        isp = mod(is+iflip-1,this%lcomp) + 1

        !! obtain new blocks of PFDP = P^{(i+1)}_{(i+1)k*N}*F*Delta*P^{(i+1)}_{N*(i+1)k}
        !!! obtain PFDP_kk
        ! Rrow = P^{(i+1)}_{k*N}*R
        Rrow(1,:) = ur%blk1(i1,:)
        Rrow(2,:) = ur%blk1(i2,:)
        ! Rtmp = Rrow * (LR)^-1
        call zgemm('N', 'N', 2, ne, ne, cone, Rrow, 2, ulrinv%blk1, ne, czero, Rtmp, 2)
        ! LDcol = L*Delta*P^{(i+1)}_{N*k}
        LDcol(:,1) = ul%blk1(:,i1)*this%delta_bmat_p%blk1(iflip, is)
        LDcol(:,2) = ul%blk1(:,i2)*this%delta_bmat_m%blk1(iflip, is)
        ! PFDP_kk = Rtmp * LDcol = P^{(i+1)}_{k*N}*F*Delta*P^{(i+1)}_{N*k}
        call zgemm('N', 'N', 2, 2, ne, cone, Rtmp, 2, LDcol, ne, czero, PFDP_kk, 2)
#IFDEF TEST
        ! check PFDP_kk with Fmmat
        do m = 1, 2
            do n = 1, 2
                if (m == 1) then
                    x = i1
                else
                    x = i2
                end if
                if (n == 1) then
                    x1 = i1
                    delta = this%delta_bmat_p%blk1(iflip, is)
                else
                    x1 = i2
                    delta = this%delta_bmat_m%blk1(iflip, is)
                end if
                if ( abs(PFDP_kk(m,n)-Fmmat%blk1(x,x1)*delta) .gt. 1.d-8 ) then
                    write(fout, '(a,6i4)') ' in update_v, PFDP_kk != Fmmat, i, nf, m, n, x, x1 = ', i, nf, m, n, x, x1 
                    write(fout, '(a,2e16.8)') ' in update_v, PFDP_kk(m,n) = ', PFDP_kk(m,n)
                    write(fout, '(a,2e16.8)') ' in update_v, Fmmat(x,x1)*delta = ', Fmmat%blk1(x,x1)*delta
                end if
            end do
        end do
#ENDIF
        !
        if ( ik > 0 ) then
            ! check if we need to reallocate the matrices
            if (ik_mat < ik) then
                ik_mat = ik
                if ( allocated( ikktmp ) ) deallocate( ikktmp )
                allocate( ikktmp(ik,2) )
                if ( allocated( PFDP_ikk ) ) deallocate( PFDP_ikk )
                allocate( PFDP_ikk(ik,2) )
                if ( allocated( PFDP_kik ) ) deallocate( PFDP_kik )
                allocate( PFDP_kik(2,ik) )
            end if
            !!! obtain PFDP_kik
            ! neiktmp = L*Delta^{(i)}*P^{(i)}_{N*ik}
            neiktmp = Umat%blk1(:,1:ik)
            ! PFDP_kik = Rtmp * neiktmp = P^{(i+1)}_{k*N}*F*Delta^{(i)}*P^{(i)}_{N*ik}
            call zgemm('N', 'N', 2, ik, ne, cone, Rtmp, 2, neiktmp, ne, czero, PFDP_kik, 2)
            !
            !!! obtain PFDP_ikk
            ! iknetmp = P^{(i)}_{ik*N}*R*ulrinv
            iknetmp = Vmat%blk1(1:ik,:)
            ! PFDP_ikk = iknetmp * LDcol = P^{(i)}_{ik*N}*F*Delta*P^{(i+1)}_{N*k}
            call zgemm('N', 'N', ik, 2, ne, cone, iknetmp, ik, LDcol, ne, czero, PFDP_ikk, ik)
#IFDEF TEST
            ! check PFDP_kik with Fmmat
            do m = 1, 2
                do n = 1, ik
                    if (m == 1) then
                        x = i1
                    else
                        x = i2
                    end if
                    x1 = xvec(n)
                    delta = Dvec(n)
                    if ( abs(PFDP_kik(m,n)-Fmmat%blk1(x,x1)*delta) .gt. 1.d-8 ) then
                        write(fout, '(a,6i4)') ' in update_v, PFDP_kik != Fmmat, i, nf, m, n, x, x1 = ', i, nf, m, n, x, x1
                        write(fout, '(a,2e16.8)') ' in update_v, PFDP_kik(m,n) = ', PFDP_kik(m,n)
                        write(fout, '(a,2e16.8)') ' in update_v, Fmmat(x,x1)*delta = ', Fmmat%blk1(x,x1)*delta
                    end if
                end do
            end do
            ! check PFDP_ikk with Fmmat
            do m = 1, ik
                do n = 1, 2
                    x = xvec(m)
                    if (n == 1) then
                        x1 = i1
                        delta = this%delta_bmat_p%blk1(iflip, is)
                    else
                        x1 = i2
                        delta = this%delta_bmat_m%blk1(iflip, is)
                    end if
                    if ( abs(PFDP_ikk(m,n)-Fmmat%blk1(x,x1)*delta) .gt. 1.d-8 ) then
                        write(fout, '(a,6i4)') ' in update_v, PFDP_ikk != Fmmat, i, nf, m, n, x, x1 = ', i, nf, m, n, x, x1
                        write(fout, '(a,2e16.8)') ' in update_v, PFDP_ikk(m,n) = ', PFDP_ikk(m,n)
                        write(fout, '(a,2e16.8)') ' in update_v, Fmmat(x,x1)*delta = ', Fmmat%blk1(x,x1)*delta
                    end if
                end do
            end do
#ENDIF
        end if

        !! vukmat = I + Vcal*D
        ! I + PFDP_kk
        vukmat(1,1) = PFDP_kk(1,1) + cone
        vukmat(1,2) = PFDP_kk(1,2)
        vukmat(2,2) = PFDP_kk(2,2) + cone
        vukmat(2,1) = PFDP_kk(2,1)
        if ( ik > 0 ) then
            ! iktmp = Gammainv
            iktmp = Gammainv(1:ik,1:ik)
            ! ikktmp = iktmp*PFDP_ikk
            call zgemm('N', 'N', ik, 2, ik, cone, iktmp, ik, PFDP_ikk, ik, czero, ikktmp, ik)
            ! vukmat = vukmat - PFDP_kik*ikktmp
            call zgemm('N', 'N', 2, 2, ik, -cone, PFDP_kik, 2, ikktmp, ik, cone, vukmat, 2)
        end if

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
            !! obtain Gammainv_kk = (I + Vcal*D)^-1 = vukmat^-1
            Gammainv_kk(1,1) = vukmat(2,2)
            Gammainv_kk(1,2) = -vukmat(1,2)
            Gammainv_kk(2,1) = -vukmat(2,1)
            Gammainv_kk(2,2) = vukmat(1,1)
            Gammainv_kk = Gammainv_kk/ratio1
            ! update 22 block of Gammainv
            Gammainv(ik+1:ik+2,ik+1:ik+2) = Gammainv_kk(:,:)
            if ( ik > 0 ) then
                !! obtain Gammainv_ikk = - Gammainv^{(i)} * PFDP_ikk * Gammainv_kk
                ! we already have ikktmp = Gammainv^{(i)} * PFDP_ikk
                ! PFDP_ikk is useless now
                ! PFDP_ikk = - ikktmp * Gammainv_kk = Gammainv_ikk
                call zgemm('N', 'N', ik, 2, 2, -cone, ikktmp, ik, Gammainv_kk, 2, czero, PFDP_ikk, ik)
                ! update 12 block of Gammainv
                Gammainv(1:ik,ik+1:ik+2) = PFDP_ikk(:,:)

                !! obtain Gammainv_kik = - Gammainv_kk * PFDP_kik * Gammainv^{(i)}
                ! we already have iktmp = Gammainv^{(i)}
                ! kiktmp = PFDP_kik * Gammainv^{(i)}
                kiktmp = PFDP_kik(:,:) ! a new (re-)allocated temporary matrix
                call zgemm('N', 'N', 2, ik, ik, cone, PFDP_kik, 2, iktmp, ik, czero, kiktmp, 2)
                ! PFDP_kik is useless now
                ! PFDP_kik = - Gammainv_kk * kiktmp = Gammainv_kik
                call zgemm('N', 'N', 2, ik, 2, -cone, Gammainv_kk, 2, kiktmp, 2, czero, PFDP_kik, 2)
                ! update 21 block of Gammainv
                Gammainv(ik+1:ik+2,1:ik) = PFDP_kik(:,:)

                !! obtain Gammainv_ik = Gammainv^{(i)} - Gammainv^{(i)} * PFDP_ikk * Gammainv_kik
                ! we already have ikktmp = Gammainv^{(i)} * PFDP_ikk
                ! we already have PFDP_kik = Gammainv_kik
                ! we already have iktmp = Gammainv^{(i)}
                ! iktmp = iktmp - ikktmp * PFDP_kik
                call zgemm('N', 'N', ik, ik, 2, -cone, ikktmp, ik, PFDP_kik, 2, cone, iktmp, ik)
                ! update 11 block of Gammainv
                Gammainv(1:ik,1:ik) = iktmp(:,:)
            end if

            ik = ik + 2
            ! store xvec and Dvec
            xvec(ik-1) = i1
            xvec(ik) = i2
            Dvec(ik-1) = this%delta_bmat_p%blk1(iflip,is)
            Dvec(ik) = this%delta_bmat_m%blk1(iflip,is)
            ! append Rtmp, LDcol to Vmat%blk1, Umat%blk1
            Vmat%blk1(ik-1:ik,:) = Rtmp(:,:)
            Umat%blk1(:,ik-1:ik) = LDcol(:,:)

            ! flip
            this%conf(i,nf,ntau) =  isp
        end if

        if( (ik.eq.nublock) .or. (i.eq.latt%nn_lf) ) then
#IFDEF TIMING
            call cpu_time_now(starttime11)
#ENDIF
            ! delay update: update (LR)^-1 and R
            
            !! update (LR)^-1
            ! (LR)^-1 = (LR)^-1 - (LR)^-1 * Umat * (I + Vmat*Umat)^-1 * Vmat
            ! we already have Umat=L*Delta^{(i)}*P^{(i)}_{N*ik} and Vmat = P^{(i)}_{ik*N}*R*(LR)^-1
            ! we already have Gammainv = (I + Vmat*Umat)^-1
            ! 
            ! obtain Vtmp = Gammainv * Vmat
            call zgemm('N', 'N', nublock, ne, nublock, cone, Gammainv, nublock, Vmat%blk1, nublock, czero, Vtmp%blk1, nublock)
            ! obtain Utmp = ulrinv * Umat
            call zgemm('N', 'N', ne, nublock, ne, cone, ulrinv%blk1, ne, Umat%blk1, ne, czero, Utmp%blk1, ne)
            ! obtain new ulrinv = ulrinv - (ulrinv * Umat) * (Gammainv * Vmat)
            call zgemm('N', 'N', ne, ne, nublock, -cone, Utmp%blk1, ne, Vtmp%blk1, nublock, cone, ulrinv%blk1, ne)
            
            !! update R^{(n)} = (I + Delta^{(i)}) R^{(0)}
            do m = 1, ik
                x = xvec(m)
                delta = Dvec(m)
                ur%blk1(x,:) = ur%blk1(x,:) + delta*ur%blk1(x,:)
            end do
#IFDEF TEST
            ! update Fmmat
            ! Rmmat = R * (LR)^-1
            call zgemm('N', 'N', ndim, ne, ne, cone, ur%blk1, ndim, ulrinv%blk1, ne, czero, Rmmat%blk1, ndim)
            ! Fmmat = Rmmat * L
            call zgemm('N', 'N', ndim, ndim, ne, cone, Rmmat%blk1, ndim, ul%blk1, ne, czero, Fmmat%blk1, ndim)
            write(fout, '(a,2e16.8)') ' in update_v, finish delay update'
#ENDIF

            ik = 0
            ik_mat = 0
            ! initial xvec, Dvec and Gammainv
            xvec = czero
            Dvec = czero
            Gammainv = czero
            ! initial Vmat, Umat
            Vmat%blk1 = czero
            Umat%blk1 = czero

#IFDEF TIMING
            call cpu_time_now(endtime11)
            timecalculation(21)=timecalculation(21)+endtime11-starttime11
#ENDIF
        end if
    end do
    main_obs(4) = main_obs(4) + dcmplx( accm, latt%nn_lf )

#IFDEF TIMING
    call cpu_time_now(endtime)
    timecalculation(15)=timecalculation(15)+endtime-starttime
#ENDIF
end subroutine dqmc_proj_update
