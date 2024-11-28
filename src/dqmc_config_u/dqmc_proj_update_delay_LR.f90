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
    integer :: iflip, is, isp, isite, i, ik, m, n, x
    real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

    complex(dp) :: vukmat
    complex(dp), allocatable, dimension(:) :: Rrow, Rtmp

    integer, allocatable, dimension(:) :: xvec
    complex(dp), allocatable, dimension(:) :: Dvec
    complex(dp), allocatable, dimension(:) :: PFvec
    complex(dp), allocatable, dimension(:,:) :: PFvecs
    complex(dp), allocatable, dimension(:,:) :: iktmp
    complex(dp), allocatable, dimension(:) :: ikktmp, kiktmp, kiktmp2
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
    allocate( xvec(nublock) )        ! x^{(i)}
    allocate( Dvec(nublock) )        ! Delta_{i}, i=1,...,nublock
    allocate( PFvecs(nublock,ndim) ) ! P^{(i)}_{i*N}*F
    allocate( PFvec(ndim) )        ! P^{(i)}_{1*N}*F
    allocate( Rrow(ne) )           ! P^{(i)}_{1*N}*R
    allocate( Rtmp(ne) )

    accm  = 0.d0 ! acceptance rate
    ik = 0       ! delay step * k

    do isite = 1, latt%nsites
        is = this%conf(isite,ntau)
        iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
        isp = mod(is+iflip-1,this%lcomp) + 1

        !! obtain PFvec = P^{(i)}_{1*N}*F = P^{(i)}_{1*N}*R*(LR)^-1*L
        ! Rrow = P^{(i)}_{1*N}*R
        Rrow(:) = ur%blk1(isite,:)
        ! Rtmp = Rrow * (LR)^-1
        call zgemv('T', ne, ne, cone, ulrinv%blk1, ne, Rrow, 1, czero, Rtmp, 1)
        ! PFvec = Rtmp * ul%blk1
        call zgemv('T', ne, ndim, cone, ul%blk1, ne, Rtmp, 1, czero, PFvec, 1)

        !! vukmat = 1 + Vcal*D
        ! 1 + PFvec * Delta_{i}P
        vukmat = PFvec(isite)*this%delta_bmat%blk1(iflip,is) + cone

        if ( ik > 0 ) then
            
            ! ikmat = (I + PFvecs * Delta^{(i-1)}P)^-1
            if ( allocated( iktmp ) ) deallocate( iktmp )
            allocate( iktmp(ik,ik) )
!$OMP PARALLEL &
!$OMP PRIVATE ( m, n, x, delta )
!$OMP DO
            do m = 1, ik
                do n = 1, ik
                    x = xvec(n)
                    delta = Dvec(n)
                    iktmp(m,n) = PFvecs(m,x)*delta
                end do
            end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
            do m = 1, ik
                iktmp(m,m) = iktmp(m,m) + cone
            end do
!$OMP END DO
!$OMP END PARALLEL
            call s_invlu_z(ik, iktmp)
        
            ! ikktmp = PFvecs * Delta_{i}P
            if ( allocated( ikktmp) ) deallocate( ikktmp )
            allocate( ikktmp(ik) )
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
            do m = 1, ik
                ikktmp(m) = PFvecs(m,isite)*this%delta_bmat%blk1(iflip,is)
            end do
!$OMP END DO
!$OMP END PARALLEL

            ! kiktmp = PFvec * Delta^{(i-1)}P
            if ( allocated( kiktmp) ) deallocate( kiktmp )
            allocate( kiktmp(ik) )
!$OMP PARALLEL &
!$OMP PRIVATE ( n, x, delta )
!$OMP DO
            do n = 1, ik
                x = xvec(n)
                delta = Dvec(n)
                kiktmp(n) = PFvec(x)*delta
            end do
!$OMP END DO
!$OMP END PARALLEL

            ! vuktmp = vuktmp - kiktmp * iktmp * ikktmp
            if ( allocated( kiktmp2 ) ) deallocate( kiktmp2 )
            allocate( kiktmp2(ik) )
            call zgemv('T', ik, ik, cone, iktmp, ik, kiktmp, 1, czero, kiktmp2, 1)
!$OMP PRIVATE ( m )
!$OMP DO
            do m = 1, ik
                vukmat = vukmat - kiktmp2(m)*ikktmp(m)
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
            ! append PFvec to PFvecs
            PFvecs(ik,:) = PFvec(:)
            
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
            ! obtain Umat=L*Delta^{(i)} and Vtmp = P_{ik*N}*R
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
            ! obtain VUmat = I + Vmat*Umat
            VUmat%blk1 = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
            ! although VUmat may have only ik (< nublock) non-zero elements, we still use nublock here
            ! otherwise the inverse of VUmat may throw error
            ! this would not affect the final result since the two diagonal blocks in VUmat are inversed separately
            do m = 1, nublock
                VUmat%blk1(m,m) = cone
            end do
!$OMP END DO
!$OMP END PARALLEL
            call zgemm('N', 'N', nublock, nublock, ne, cone, Vmat%blk1, nublock, Umat%blk1, ne, cone, VUmat%blk1, nublock)
            ! obtain (VUmat)^-1
            call s_invlu_z(nublock, VUmat%blk1)
            ! obtain Vtmp = VUmat^-1 * Vmat
            call zgemm('N', 'N', nublock, ne, nublock, cone, VUmat%blk1, nublock, Vmat%blk1, nublock, czero, Vtmp%blk1, nublock)
            ! obtain Utmp = ulrinv * Umat
            call zgemm('N', 'N', ne, nublock, ne, cone, ulrinv%blk1, ne, Umat%blk1, ne, czero, Utmp%blk1, ne)
            ! obtain new ulrinv = ulrinv - Utmp * Vtmp
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
            ! initial xvec, Dvec and PFvecs
            xvec = czero
            Dvec = czero
            PFvecs = czero

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
