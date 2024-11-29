subroutine dqmc_proj_update(this, ntau, ul, ur, ulrinv)
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
    complex(dp), allocatable, dimension(:,:) :: PFDP
    complex(dp) :: PFDP_kk
    complex(dp), allocatable, dimension(:) :: PFDP_ikk, PFDP_kik
    complex(dp), allocatable, dimension(:) :: Rrow, Rtmp, LDcol
    complex(dp), allocatable, dimension(:,:) :: neiktmp, iknetmp, iknetmp2 ! change size at each step
    complex(dp), allocatable, dimension(:,:) :: iktmp ! change size at each step
    complex(dp), allocatable, dimension(:) :: ikktmp ! change size at each step
    type(gfunc) :: Umat, Utmp, Vmat, Vtmp, VUmat, VUtmp ! for update (LR)^-1

#IFDEF TIMING
    real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
    call cpu_time_now(starttime)
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
    allocate( PFDP(nublock,nublock) ) ! P^{(i)}_{i*N}*F*Delta^{(i)}*P^{(i)}_{N*i}
    allocate( Rrow(ne) )              ! P^{(i)}_{1*N}*R
    allocate( Rtmp(ne) )
    allocate( LDcol(ne) )             ! L*Delta_{i}*P^{(i)}_{N*1}

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
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
        do m = 1, ne
            !! Warning: BLAS level 1 operation
            PFDP_kk = PFDP_kk + Rtmp(m)*LDcol(m)
        end do
!$OMP END DO
!$OMP END PARALLEL
        !
        if ( ik > 0 ) then
            ! check if we need to reallocate the matrices
            if (ik_mat < ik) then
                ik_mat = ik
                if ( allocated( iknetmp2 ) ) deallocate( iknetmp2 )
                allocate( iknetmp2(ik,ne) )
                if ( allocated( PFDP_ikk ) ) deallocate( PFDP_ikk )
                allocate( PFDP_ikk(ik) )
                if ( allocated( iknetmp ) ) deallocate( iknetmp )
                allocate( iknetmp(ik,ne) )
                if ( allocated( neiktmp ) ) deallocate( neiktmp )
                allocate( neiktmp(ne,ik) )
                if ( allocated( PFDP_kik ) ) deallocate( PFDP_kik )
                allocate( PFDP_kik(ik) )
            end if
            !!! obtain PFDP_kik
            ! neiktmp = L*Delta^{(i)}*P^{(i)}_{N*ik}
            ! iknetmp = P^{(i)}_{ik*N}*R
!$OMP PARALLEL &
!$OMP PRIVATE ( n, x, delta )
!$OMP DO
            do n = 1, ik
                x = xvec(n)
                delta = Dvec(n)
                neiktmp(:,n) = ul%blk1(:,x)*delta
                iknetmp(n,:) = ur%blk1(x,:)
            end do
!$OMP END DO
!$OMP END PARALLEL
            ! PFDP_kik = Rtmp * neiktmp = P^{(i+1)}_{k*N}*F*Delta^{(i)}*P^{(i)}_{N*ik}
            call zgemv('T', ne, ik, cone, neiktmp, ne, Rtmp, 1, czero, PFDP_kik, 1)
            !
            !!! obtain PFDP_ikk
            ! iknetmp2 = iknetmp * (LR)^-1
            call zgemm('N', 'N', ik, ne, ne, cone, iknetmp, ik, ulrinv%blk1, ne, czero, iknetmp2, ik)
            ! PFDP_ikk = iknetmp2 * LDcol = P^{(i)}_{ik*N}*F*Delta*P^{(i+1)}_{N*k}
            call zgemv('N', ik, ne, cone, iknetmp2, ik, LDcol, 1, czero, PFDP_ikk, 1)
        end if


        !! vukmat = 1 + Vcal*D
        ! 1 + PFDP_kk
        vukmat = cone + PFDP_kk

        if ( ik > 0 ) then
            ! ikmat = (I + PFDP)^-1
            iktmp = PFDP(1:ik,1:ik)
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
            do m = 1, ik
                iktmp(m,m) = iktmp(m,m) + cone
            end do
!$OMP END DO
!$OMP END PARALLEL
            call s_invlu_z(ik, iktmp)
        
            ! ikktmp = iktmp*PFDP_ikk
            ikktmp = PFDP_ikk(:)
            call zgemv('N', ik, ik, cone, iktmp, ik, PFDP_ikk, 1, czero, ikktmp, 1)
            ! vukmat = vukmat - PFDP_kik * ikktmp
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
            do m = 1, ik
                !! warning: BLAS level 1 operation
                vukmat = vukmat - PFDP_kik(m)*ikktmp(m)
            end do
!$OMP END DO
!$OMP END PARALLEL
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

            ik = ik + 1
            ! store xvec and Dvec
            xvec(ik) = isite
            Dvec(ik) = this%delta_bmat%blk1(iflip,is)
            ! enlarge PFDP
            PFDP(ik,ik) = PFDP_kk
            if ( ik > 1 ) then
                PFDP(1:ik-1,ik) = PFDP_ikk
                PFDP(ik,1:ik-1) = PFDP_kik
            end if
            
            ! flip
            this%conf(isite,ntau) =  isp
        end if

        if( (ik.eq.nublock) .or. (isite.eq.latt%nsites) ) then
#IFDEF TIMING
            call cpu_time_now(starttime11)
#ENDIF     
            ! delay update: update (LR)^-1 and R
            ! note that since (LR)^-1 depends on R, we should R later
            
            !! update (LR)^-1
            ! obtain Umat=L*Delta^{(i)}*P^{(i)}_{N*ik} and Vtmp = P_{ik*N}*R
            Umat%blk1 = czero
            Vtmp%blk1 = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( m, x, delta )
!$OMP DO
            do m = 1, ik
                x = xvec(m)
                delta = Dvec(m)
                Umat%blk1(:,m) = ul%blk1(:,x)*delta
                Vtmp%blk1(m,:) = ur%blk1(x,:)
            end do
!$OMP END DO
!$OMP END PARALLEL
            ! obtain Vmat = Vtmp*ulrinv
            call zgemm('N', 'N', nublock, ne, ne, cone, Vtmp%blk1, nublock, ulrinv%blk1, ne, czero, Vmat%blk1, nublock)
            !
            ! obtain VUmat = (I + Vmat*Umat)^-1
            ! VUmat = Vmat*Umat = PFDP
            VUmat%blk1 = PFDP(:,:)
            ! VUmat = I + Vmat*Umat
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
            do m = 1, nublock
                ! although VUmat may have only ik (< nublock) non-zero elements, we still use nublock here
                ! otherwise the inverse of VUmat may throw error
                ! this would not affect the final result since the two diagonal blocks in VUmat are inversed separately
                VUmat%blk1(m,m) = VUmat%blk1(m,m) + cone
            end do
!$OMP END DO
!$OMP END PARALLEL
            ! VUmat = (I + Vmat*Umat)^-1
            call s_invlu_z(nublock, VUmat%blk1)
            !
            ! obtain Vtmp = VUmat * Vmat
            call zgemm('N', 'N', nublock, ne, nublock, cone, VUmat%blk1, nublock, Vmat%blk1, nublock, czero, Vtmp%blk1, nublock)
            ! obtain Utmp = ulrinv * Umat
            call zgemm('N', 'N', ne, nublock, ne, cone, ulrinv%blk1, ne, Umat%blk1, ne, czero, Utmp%blk1, ne)
            ! obtain new ulrinv = ulrinv - (ulrinv * Umat) * ((I + Vmat*Umat)^-1 * Vmat)
            call zgemm('N', 'N', ne, ne, nublock, -cone, Utmp%blk1, ne, Vtmp%blk1, nublock, cone, ulrinv%blk1, ne)
            
            !! update R^{(n)} = (I + Delta^{(i)}) R^{(0)}
!$OMP PARALLEL &
!$OMP PRIVATE ( m, x, delta )
!$OMP DO
            do m = 1, ik
                x = xvec(m)
                delta = Dvec(m)
                ur%blk1(x,:) = ur%blk1(x,:) + delta*ur%blk1(x,:)
            end do
!$OMP END DO
!$OMP END PARALLEL

            ik = 0
            ik_mat = 0
            ! initial xvec, Dvec and PFvecs
            xvec = czero
            Dvec = czero
            PFDP = czero

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
