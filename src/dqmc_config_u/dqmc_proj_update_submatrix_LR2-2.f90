subroutine dqmc_proj_update(this, ntau, ul, ur, ulrinv)
!! update the configuration using submatrix LR2-2
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
    class(uconf), intent(inout) :: this
    integer, intent(in) :: ntau
    type(gfunc), intent(inout) :: ul, ur, ulrinv

    ! local
    complex(dp) ::  ratio1, ratiotot, delta
    integer :: iflip, is, isp, isite, i, ik, m, n, x
    real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

    complex(dp) :: vukmat

    integer, allocatable, dimension(:) :: xvec
    complex(dp), allocatable, dimension(:) :: Dvec
    complex(dp), allocatable, dimension(:,:) :: Gammainv
    complex(dp), allocatable, dimension(:) :: Rrow, Lcol ! intermediate matrices when calculating PFDP
    complex(dp) :: PFP_kk
    complex(dp), allocatable, dimension(:) :: PFP_ikk, PFP_kik ! change size at each step
    complex(dp), allocatable, dimension(:) :: ikktmp ! change size at each step
    complex(dp), allocatable, dimension(:) :: kiktmp ! for update Gammainv
    complex(dp) :: Gammainv_kk ! for update Gammainv
    type(gfunc) :: Umat, Utmp, Vmat, Vtmp, VUmat, VUtmp ! for update (LR)^-1
    integer :: nustkup, nstk, istk_s, istk_e
    complex(dp), allocatable, dimension(:,:) :: Rrows, Lcols, Felms ! stock arrays for preventing zgemv

#IFDEF TIMING
    real(dp) :: starttime, endtime, starttime11, endtime11
#ENDIF
#IFDEF TIMING
    call cpu_time_now(starttime)
#ENDIF

#IFDEF TEST
    write(fout, '(a,2e16.8)') ' in update_u using submatrix LR2'
#ENDIF

    ! for delay update (LR)^-1
    call allocate_gfunc( Umat, ne, nustock )
    call allocate_gfunc( Utmp, ne, nustock )
    call allocate_gfunc( Vmat, nustock, ne )
    call allocate_gfunc( Vtmp, nustock, ne )
    call allocate_gfunc( VUmat, nustock, nustock )
    call allocate_gfunc( VUtmp, nustock, nustock )

    ! for calculate update ratio
    allocate( xvec(nustock) )         ! x^{(i)}
    allocate( Dvec(nustock) )         ! Delta_{i}, i=1,...,nublock
    allocate( Gammainv(nustock,nustock) )

    ! matrices for delay update
    xvec = 0
    Dvec = czero
    Gammainv = czero
    Vmat%blk1 = czero
    Umat%blk1 = czero

    ! preallocate oversize arrays for intermediate matrices
    allocate( kiktmp(nustock) )
    allocate( ikktmp(nustock) )
    allocate( PFP_ikk(nustock) )
    allocate( PFP_kik(nustock) )
    kiktmp = czero
    ikktmp = czero
    PFP_ikk = czero
    PFP_kik = czero

    ! allocate stock arrays for preventing zgemv
    allocate( Rrows(ndim, ne) ) ! rows of R(LR)^-1
    ! allocate( Lcols(ne, ndim) ) ! columns of L
    allocate( Felms(ndim, ndim) ) ! selected rows and columns of F

    accm  = 0.d0 ! acceptance rate
    ik = 0       ! delay step * k
    istk_s = 0   ! starting isite in the stock arrays
    istk_e = 0   ! ending isite in the stock arrays
    nstk = 0     ! number of sites in the stock arrays

    do isite = 1, latt%nsites
        is = this%conf(isite,ntau)
        iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
        isp = mod(is+iflip-1,this%lcomp) + 1

        !! stock up intermediate matrices
        if (nstk == 0) then
            nustkup = min(nustock, latt%nsites - isite + 1)
#ifdef TEST
            write(fout, '(a,i4,a,i4)') ' in update_u, stock up intermediate matrices from ', isite, ' to ', isite+nustkup-1
#endif
            ! Rrows(isite:isite+nustock,:) = R(LR)^-1(isite:isite+nustock,:)
            call zgemm('N', 'N', nustkup, ne, ne, cone, ur%blk1(isite,1), ndim, ulrinv%blk1, ne, czero, Rrows(isite,1), ndim)
            ! Lcols(:,isite:isite+nustkup-1) = L(:,isite:isite+nustkup-1)
            ! Lcols(:,isite:isite+nustkup-1) = ul%blk1(:,isite:isite+nustkup-1)
            ! Felms = R(LR)^-1*L(isite:isite+nustock,isite:isite+nustock)
            ! call zgemm('N', 'N', nustkup, nustkup, ne, cone, Rrows, ndim, Lcols, ne, czero, Felms, nustock)
            call zgemm('N', 'N', nustkup, nustkup, ne, cone, Rrows(isite,1), ndim, ul%blk1(1,isite), ne, czero, Felms(isite,isite), ndim)
            istk_s = isite
            istk_e = isite + nustkup - 1
            nstk = nustkup
        end if

        !! obtain new blocks of PFP = P^{(i+1)}_{(i+1)k*N}*F*P^{(i+1)}_{N*(i+1)k}
        !!! obtain PFDP_kk
        ! Rrow = P^{(i+1)}_{1*N}*R * (LR)^-1
        ! Lcol = L*P^{(i+1)}_{N*1}
        ! PFP_kk = Rrow * Lcol = P^{(i+1)}_{k*N}*F*P^{(i+1)}_{N*k}
        PFP_kk = Felms(isite, isite)
        !
        if ( ik > 0 ) then
            !!! obtain PFP_kik
            ! Umat%blk1(:,1:ik) = L*P^{(i)}_{N*ik}
            ! PFP_kik = Rrow * Umat%blk1(:,1:ik) = P^{(i+1)}_{k*N}*F*P^{(i)}_{N*ik}
            PFP_kik(1:ik) = Felms(isite, xvec(1:ik))
            !
            !!! obtain PFP_ikk
            ! Vmat%blk1(1:ik,:) = P^{(i)}_{ik*N}*R*ulrinv
            ! PFP_ikk = Vmat%blk1(1:ik,:) * Lcol = P^{(i)}_{ik*N}*F*P^{(i+1)}_{N*k}
            PFP_ikk(1:ik) = Felms(xvec(1:ik), isite)
        end if


        !! vukmat = 1 + Vcal*D
        vukmat = czero
        ! vukmat = PFP_kik * Gammainv * PFDP_ikk
        if ( ik > 0 ) then
            ! ikktmp = Gammainv(1:ik,1:ik)*PFP_ikk
            call zgemv('N', ik, ik, cone, Gammainv, nublock, PFP_ikk, 1, czero, ikktmp, 1)
            ! vukmat = PFP_kik * ikktmp
            do m = 1, ik
                !! warning: BLAS level 1 operation (inner product)
                vukmat = vukmat + PFP_kik(m)*ikktmp(m)
            end do
        end if
        ! vukmat = Vcal = PFP_kk - vukmat
        vukmat = PFP_kk - vukmat
        ! vukmat = 1 + Vcal*D
        vukmat = cone + vukmat * this%delta_bmat%blk1(iflip,is)
        
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
            !! obtain Gammainv_kk = (D^{(i)}_{k*k})(I + Vcal*D)^-1 = (D^{(i)}_{k*k})*vukmat^-1
            Gammainv_kk = this%delta_bmat%blk1(iflip,is)/vukmat
            ! update 22 block of Gammainv
            Gammainv(ik+1,ik+1) = Gammainv_kk
            if ( ik > 0 ) then
                !! obtain Gammainv_kik = - Gammainv_kk * PFP_kik * Gammainv^{(i)}
                ! kiktmp = PFP_kik * Gammainv^{(i)}
                call zgemv('T', ik, ik, cone, Gammainv, nublock, PFP_kik, 1, czero, kiktmp, 1)
                ! Gammainv_kik = - Gammainv_kk * kiktmp
                ! update 21 block of Gammainv
                Gammainv(ik+1,1:ik) = - Gammainv_kk * kiktmp(1:ik)

                !! obtain Gammainv_ikk = - Gammainv^{(i)} * PFP_ikk * Gammainv_kk
                ! we already have ikktmp = Gammainv^{(i)} * PFP_ikk
                ! Gammainv_ikk = - ikktmp * Gammainv_kk
                ! update 12 block of Gammainv
                Gammainv(1:ik,ik+1) = - ikktmp(1:ik) * Gammainv_kk

                !! obtain Gammainv_ik = Gammainv^{(i)} - Gammainv^{(i)} * PFP_ikk * Gammainv_kik
                ! we already have ikktmp = Gammainv^{(i)} * PFP_ikk
                ! we already have Gammainv(ik+1,1:ik) = Gammainv_kik
                ! Gammainv_ik = Gammainv_ik - ikktmp * Gammainv(ik+1,1:ik)
                ! calculate and update 11 block of Gammainv
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
                do n = 1, ik
                    do m = 1, ik
                        !! warning: BLAS level 1 operation (outer product)
                        Gammainv(m,n) = Gammainv(m,n) - ikktmp(m)*Gammainv(ik+1,n)
                    end do
                end do
!$OMP END DO
!$OMP END PARALLEL
            end if

            ik = ik + 1
            ! store xvec and Dvec
            xvec(ik) = isite
            Dvec(ik) = this%delta_bmat%blk1(iflip,is)
            ! append Rrow, Lcol to Vmat%blk1, Umat%blk1
            Vmat%blk1(ik,:) = Rrows(isite,:)
            Umat%blk1(:,ik) = ul%blk1(:,isite)
            
            ! flip
            this%conf(isite,ntau) =  isp
        end if

        ! when the stock arrays are running out, update (LR)^-1
        if( (isite == istk_e) .or. (isite == latt%nsites) ) then
#IFDEF TIMING
            call cpu_time_now(starttime11)
#ENDIF     
#ifdef TEST
            write(fout, '(a,i4)') ' in update_u, update the whole LRinv'
#endif
            ! delay update: update (LR)^-1 and R
            
            !! update (LR)^-1
            ! (LR)^-1 = (LR)^-1 - (LR)^-1 * Umat * (D^-1 + Vmat*Umat)^-1 * Vmat
            ! we already have Umat=L*P^{(i)}_{N*ik} and Vmat = P^{(i)}_{ik*N}*R*(LR)^-1
            ! we already have Gammainv = (D^-1 + Vmat*Umat)^-1
            !
            ! obtain Vtmp = Gammainv * Vmat
            call zgemm('N', 'N', ik, ne, ik, cone, Gammainv, nublock, Vmat%blk1, nublock, czero, Vtmp%blk1, nublock)
            ! obtain Utmp = ulrinv * Umat
            call zgemm('N', 'N', ne, ik, ne, cone, ulrinv%blk1, ne, Umat%blk1, ne, czero, Utmp%blk1, ne)
            ! obtain new ulrinv = ulrinv - (ulrinv * Umat) * (Gammainv * Vmat)
            call zgemm('N', 'N', ne, ne, ik, -cone, Utmp%blk1, ne, Vtmp%blk1, nublock, cone, ulrinv%blk1, ne)
            
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
            ! initial xvec, Dvec and Gammainv
            xvec = czero
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
    main_obs(3) = main_obs(3) + dcmplx( accm, latt%nsites )

#IFDEF TIMING
    call cpu_time_now(endtime)
    timecalculation(17)=timecalculation(17)+endtime-starttime
#ENDIF
end subroutine dqmc_proj_update
