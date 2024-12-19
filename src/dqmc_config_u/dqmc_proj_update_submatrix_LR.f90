subroutine dqmc_proj_update(this, ntau, ul, ur, ulrinv)
!! update the configuration using submatrix LR1-1
!! intemediate vectors are re-allocated when their size changes
!! the matrix-vector multiplications are done by calling zgemv
    use spring
    use model_para
    use dqmc_basic_data
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    implicit none

    ! arguments
    class(uconf), intent(inout) :: this
    integer, intent(in) :: ntau
    type(gfunc), intent(inout) :: ul, ur, ulrinv

    ! local
    complex(dp) ::  ratio1, ratiotot, delta
    integer :: iflip, is, isp, isite, i, ik, ik_mat, m, n, x
    real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

    complex(dp) :: vukmat

    integer, allocatable, dimension(:) :: xvec
    complex(dp), allocatable, dimension(:) :: Dvec
    complex(dp), allocatable, dimension(:,:) :: Gammainv
    complex(dp), allocatable, dimension(:) :: Rrow, Rtmp, LDcol ! intermediate matrices when calculating PFDP
    complex(dp) :: PFDP_kk
    complex(dp), allocatable, dimension(:) :: PFDP_ikk, PFDP_kik ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: neiktmp, iknetmp ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: iktmp ! change size at each step
    complex(dp), allocatable, dimension(:) :: ikktmp ! change size at each step
    complex(dp), allocatable, dimension(:) :: kiktmp ! for update Gammainv
    complex(dp) :: Gammainv_kk ! for update Gammainv
    type(gfunc) :: Umat, Utmp, Vmat, Vtmp, VUmat, VUtmp ! for update (LR)^-1

#IFDEF TIMING
    real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
    call cpu_time_now(starttime)
#ENDIF

#IFDEF TEST
    write(fout, '(a,2e16.8)') ' in update_u using submatrix LR'
#ENDIF

    ! for delay update (LR)^-1
    call allocate_gfunc( Umat, ne, nublock )
    call allocate_gfunc( Utmp, ne, nublock )
    call allocate_gfunc( Vmat, nublock, ne )
    call allocate_gfunc( Vtmp, nublock, ne )
    call allocate_gfunc( VUmat, nublock, nublock )
    call allocate_gfunc( VUtmp, nublock, nublock )

    ! for calculate update ratio
    allocate( xvec(nublock) )         ! x^{(i)}
    allocate( Dvec(nublock) )         ! Delta_{i}, i=1,...,nublock
    allocate( Rrow(ne) )              ! P^{(i)}_{1*N}*R
    allocate( Rtmp(ne) )
    allocate( LDcol(ne) )             ! L*Delta_{i}*P^{(i)}_{N*1}
    allocate( Gammainv(nublock,nublock) )

    ! matrices for delay update
    xvec = 0
    Dvec = czero
    Gammainv = czero
    Vmat%blk1 = czero
    Umat%blk1 = czero

    accm  = 0.d0 ! acceptance rate
    ik = 0       ! delay step * k
    ik_mat = 0   ! size of iktmp

    do isite = 1, latt%nsites
        is = this%conf(isite,ntau)
        iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
        isp = mod(is+iflip-1,this%lcomp) + 1

        !! obtain new blocks of PFDP = P^{(i+1)}_{(i+1)k*N}*F*Delta*P^{(i+1)}_{N*(i+1)k}
        !!! obtain PFDP_kk
        ! Rrow = P^{(i+1)}_{1*N}*R
        Rrow(:) = ur%blk1(isite,:)
        ! Rtmp = Rrow * (LR)^-1
        call zgemv('T', ne, ne, cone, ulrinv%blk1, ne, Rrow, 1, czero, Rtmp, 1)
        ! LDcol = L*Delta*P^{(i+1)}_{N*1}
        LDcol(:) = ul%blk1(:,isite)*this%delta_bmat%blk1(iflip,is)
        ! PFDP_kk = Rtmp * LDcol = P^{(i+1)}_{k*N}*F*Delta*P^{(i+1)}_{N*k}
        PFDP_kk = czero
        do m = 1, ne
            !! Warning: BLAS level 1 operation (inner product)
            PFDP_kk = PFDP_kk + Rtmp(m)*LDcol(m)
        end do
        !
        if ( ik > 0 ) then
            ! check if we need to reallocate the matrices
            if (ik_mat < ik) then
                ik_mat = ik
                if ( allocated( ikktmp ) ) deallocate( ikktmp )
                allocate( ikktmp(ik) )
                if ( allocated( PFDP_ikk ) ) deallocate( PFDP_ikk )
                allocate( PFDP_ikk(ik) )
                if ( allocated( PFDP_kik ) ) deallocate( PFDP_kik )
                allocate( PFDP_kik(ik) )
            end if
            !!! obtain PFDP_kik
            ! neiktmp = L*Delta^{(i)}*P^{(i)}_{N*ik}
            neiktmp = Umat%blk1(:,1:ik)
            ! PFDP_kik = Rtmp * neiktmp = P^{(i+1)}_{k*N}*F*Delta^{(i)}*P^{(i)}_{N*ik}
            call zgemv('T', ne, ik, cone, neiktmp, ne, Rtmp, 1, czero, PFDP_kik, 1)
            !
            !!! obtain PFDP_ikk
            ! iknetmp = P^{(i)}_{ik*N}*R*ulrinv
            iknetmp = Vmat%blk1(1:ik,:)
            ! PFDP_ikk = iknetmp * LDcol = P^{(i)}_{ik*N}*F*Delta*P^{(i+1)}_{N*k}
            call zgemv('N', ik, ne, cone, iknetmp, ik, LDcol, 1, czero, PFDP_ikk, 1)
        end if


        !! vukmat = 1 + Vcal*D
        ! 1 + PFDP_kk
        vukmat = cone + PFDP_kk
        if ( ik > 0 ) then
            ! iktmp = Gammainv
            iktmp = Gammainv(1:ik,1:ik)
            ! ikktmp = iktmp*PFDP_ikk
            call zgemv('N', ik, ik, cone, iktmp, ik, PFDP_ikk, 1, czero, ikktmp, 1)
            ! vukmat = vukmat - PFDP_kik * ikktmp
            do m = 1, ik
                !! warning: BLAS level 1 operation (inner product)
                vukmat = vukmat - PFDP_kik(m)*ikktmp(m)
            end do
        end if
        
        ratio1 = vukmat
        ratiotot = ratio1**nflr
        if( lproju ) then
            ratiotot = ratiotot*this%phase_ratio(iflip, is)
        else
            ratiotot = ratiotot*dcmplx( gaml(isp)/gaml(is), 0.d0 ) * this%phase_ratio(iflip, is)
        end if
#IFDEF TEST
        write(fout, '(a,2e16.8)') ' in update_u, ratio1 = ', ratio1
        write(fout, '(a,2e16.8)') ' in update_u, phase_ratio = ', this%phase_ratio(iflip, is)
        write(fout, '(a,2e16.8)') ' in update_u, ratiotot = ', ratiotot
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
            write(fout, '(a,2e16.8)') ' in update_u, phase = ', phase
#ENDIF

            ! update Gammainv: from (ik*ik) Gammainv^{(i)} to ((ik+1)*(ik+1)) Gammainv^{(i+1)} 
            !! obtain Gammainv_kk = (I + Vcal*D)^-1 = vukmat^-1
            Gammainv_kk = cone/vukmat
            ! update 22 block of Gammainv
            Gammainv(ik+1,ik+1) = Gammainv_kk
            if ( ik > 0 ) then
                !! obtain Gammainv_ikk = - Gammainv^{(i)} * PFDP_ikk * Gammainv_kk
                ! we already have ikktmp = Gammainv^{(i)} * PFDP_ikk
                ! PFDP_ikk is useless now
                ! PFDP_ikk = - ikktmp * Gammainv_kk = Gammainv_ikk
                PFDP_ikk(:) = - ikktmp(:) * Gammainv_kk
                ! update 12 block of Gammainv
                Gammainv(1:ik,ik+1) = PFDP_ikk(:)

                !! obtain Gammainv_kik = - Gammainv_kk * PFDP_kik * Gammainv^{(i)}
                ! we already have iktmp = Gammainv^{(i)}
                ! kiktmp = PFDP_kik * Gammainv^{(i)}
                kiktmp = PFDP_kik(:) ! a new (re-)allocated temporary vector
                call zgemv('T', ik, ik, cone, iktmp, ik, PFDP_kik, 1, czero, kiktmp, 1)
                ! PFDP_kik is useless now
                ! PFDP_kik = - Gammainv_kk * kiktmp = Gammainv_kik
                PFDP_kik(:) = - Gammainv_kk * kiktmp(:)
                ! update 21 block of Gammainv
                Gammainv(ik+1,1:ik) = PFDP_kik(:)

                !! obtain Gammainv_ik = Gammainv^{(i)} - Gammainv^{(i)} * PFDP_ikk * Gammainv_kik
                ! we already have ikktmp = Gammainv^{(i)} * PFDP_ikk
                ! we already have PFDP_kik = Gammainv_kik
                ! Gammainv_ik = Gammainv_ik - ikktmp * PFDP_kik
                ! calculate and update 11 block of Gammainv
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
                do n = 1, ik
                    do m = 1, ik
                        !! warning: BLAS level 1 operation (outer product)
                        Gammainv(m,n) = Gammainv(m,n) - ikktmp(m)*PFDP_kik(n)
                    end do
                end do
!$OMP END DO
!$OMP END PARALLEL
            end if

            ik = ik + 1
            ! store xvec and Dvec
            xvec(ik) = isite
            Dvec(ik) = this%delta_bmat%blk1(iflip,is)
            ! append Rtmp, LDcol to Vmat%blk1, Umat%blk1
            Vmat%blk1(ik,:) = Rtmp(:)
            Umat%blk1(:,ik) = LDcol(:)
            
            ! flip
            this%conf(isite,ntau) =  isp
        end if

        if( (ik.eq.nublock) .or. (isite.eq.latt%nsites) ) then
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
    main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )

#IFDEF TIMING
    call cpu_time_now(endtime)
    timecalculation(15)=timecalculation(15)+endtime-starttime
#ENDIF
end subroutine dqmc_proj_update
