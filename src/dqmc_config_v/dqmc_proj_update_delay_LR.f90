subroutine dqmc_proj_update(this, ntau, nf, ul, ur, ulrinv)
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
    integer :: iflip, is, isp, i, i1, i2, ik, ik_mat, m, n, x
    real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

    complex(dp), allocatable, dimension(:,:) :: vukmat
    complex(dp), allocatable, dimension(:,:) :: Rrow, Rtmp

    integer, allocatable, dimension(:) :: xvec
    complex(dp), allocatable, dimension(:) :: Dvec
    complex(dp), allocatable, dimension(:,:) :: PFvec, PFvecs
    complex(dp), allocatable, dimension(:,:) :: iktmp, ikktmp, kiktmp, kiktmp2
    type(gfunc) :: Umat, Utmp, Vmat, Vtmp, VUmat, VUtmp ! for update (LR)^-1
    type(gfunc) :: Fmmat, Rmmat

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
    allocate( vukmat(2,2) )          ! I + Vcal*D
    allocate( xvec(nublock) )        ! x^{(i)}_{j}, j=(i-1)*k+1,...,i*k
    allocate( Dvec(nublock))         ! Delta_{i}(j), i=1,...,nublock/k, j=1,...,k
    allocate( PFvecs(nublock,ndim) ) ! P^{(i)}_{ik*N}*F
    allocate( PFvec(2,ndim) )        ! P^{(i)}_{k*N}*F
    allocate( Rrow(2,ne) )           ! P^{(i)}_{k*N}*R
    allocate( Rtmp(2,ne) )           ! P^{(i)}_{k*N}*R*(LR)^-1

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
    PFvecs = czero
    Vmat%blk1 = czero

    accm  = 0.d0 ! acceptance rate
    ik = 0       ! delay step * k
    ik_mat = 0   ! size of iktmp

    do i = 1, latt%nn_lf
        i1 = latt%nnlf_list(i)
        i2 = latt%nnlist(i1,nf)
        is = this%conf(i, nf, ntau)
        iflip = ceiling( spring_sfmt_stream() * (this%lcomp-1) )
        isp = mod(is+iflip-1,this%lcomp) + 1

        !! obtain PFvec = P^{(i)}_{k*N}*F = P^{(i)}_{k*N}*R*(LR)^-1*L
        ! Rrow = P^{(i)}_{k*N}*R
        Rrow(1,:) = ur%blk1(i1,:)
        Rrow(2,:) = ur%blk1(i2,:)
        ! Rtmp = Rrow * (LR)^-1
        call zgemm('N', 'N', 2, ne, ne, cone, Rrow, 2, ulrinv%blk1, ne, czero, Rtmp, 2)
        ! PFvec = Rtmp * ul%blk1
        call zgemm('N', 'N', 2, ndim, ne, cone, Rtmp, 2, ul%blk1, ne, czero, PFvec, 2)
#IFDEF TEST
        ! check PFvec with Fmmat
        ! PFvec should be two rows of Fmmat
        do m = 1, 2
            do n = 1, ndim
                if (m == 1) then
                    x = i1
                else
                    x = i2
                end if
                if ( abs(PFvec(m,n)-Fmmat%blk1(x,n)) .gt. 1.d-8 ) then
                    write(fout, '(a,5i4)') ' in update_v, PFvec != Fmmat, i, nf, m, x, n = ', i, nf, m, x, n
                    write(fout, '(a,2e16.8)') ' in update_v, PFvec(m,n) = ', PFvec(m,n)
                    write(fout, '(a,2e16.8)') ' in update_v, Fmmat(i_m,n) = ', Fmmat%blk1(x,n)
                end if
            end do
        end do
#ENDIF

        !! vukmat = I + Vcal*D
        ! I + PFvec * Delta_{i}P
        vukmat(1,1) = PFvec(1,i1)*this%delta_bmat_p%blk1(iflip, is) + cone
        vukmat(1,2) = PFvec(1,i2)*this%delta_bmat_m%blk1(iflip, is)
        vukmat(2,1) = PFvec(2,i1)*this%delta_bmat_p%blk1(iflip, is)
        vukmat(2,2) = PFvec(2,i2)*this%delta_bmat_m%blk1(iflip, is) + cone

        if ( ik > 0 ) then

            ! check if we need to reallocate the matrices
            if (ik_mat < ik) then
                ik_mat = ik
                if ( allocated( iktmp ) ) deallocate( iktmp )
                allocate( iktmp(ik,ik) )
                if ( allocated( ikktmp) ) deallocate( ikktmp )
                allocate( ikktmp(ik,2) )
                if ( allocated( kiktmp) ) deallocate( kiktmp )
                allocate( kiktmp(2,ik) )
                if ( allocated( kiktmp2) ) deallocate( kiktmp2 )
                allocate( kiktmp2(2,ik) )
            end if

            ! ikmat = (I + PFvecs * Delta^{(i-1)}P)^-1
!$OMP PARALLEL &
!$OMP PRIVATE ( n, x, delta )
!$OMP DO
            do n = 1, ik
                x = xvec(n)
                delta = Dvec(n)
                iktmp(:,n) = PFvecs(1:ik,x)*delta
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
!$OMP PARALLEL &
!$OMP PRIVATE ( m )
!$OMP DO
            do m = 1, ik
                ikktmp(m,1) = PFvecs(m,i1)*this%delta_bmat_p%blk1(iflip,is)
                ikktmp(m,2) = PFvecs(m,i2)*this%delta_bmat_m%blk1(iflip,is)
            end do
!$OMP END DO
!$OMP END PARALLEL

            ! kiktmp = PFvec * Delta^{(i-1)}P
!$OMP PARALLEL &
!$OMP PRIVATE ( n, x, delta )
!$OMP DO
            do n = 1, ik
                x = xvec(n)
                delta = Dvec(n)
                kiktmp(1,n) = PFvec(1,x)*delta
                kiktmp(2,n) = PFvec(2,x)*delta
            end do
!$OMP END DO
!$OMP END PARALLEL

            ! vuktmp = vuktmp - kiktmp * iktmp * ikktmp
            call zgemm('N', 'N', 2, ik, ik, cone, kiktmp, 2, iktmp, ik, czero, kiktmp2, 2)
            call zgemm('N', 'N', 2, 2, ik, -cone, kiktmp2, 2, ikktmp, ik, cone, vukmat, 2)
        
        end if

        call s_inv_det_qr_z(2, vukmat, ratio1) ! cal det(I+vukmat)
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

            ik = ik + 2
            ! store xvec and Dvec
            xvec(ik-1) = i1
            xvec(ik) = i2
            Dvec(ik-1) = this%delta_bmat_p%blk1(iflip,is)
            Dvec(ik) = this%delta_bmat_m%blk1(iflip,is)
            ! append PFvec to PFvecs
            PFvecs(ik-1,:) = PFvec(1,:)
            PFvecs(ik,:) = PFvec(2,:)
            ! append Rtmp to Vmat
            Vmat%blk1(ik-1:ik,:) = Rtmp(:,:)

            ! flip
            this%conf(i,nf,ntau) =  isp
        end if

        if( (ik.eq.nublock) .or. (i.eq.latt%nn_lf) ) then
#IFDEF TIMING
            call cpu_time_now(starttime11)
#ENDIF
            ! delay update: update (LR)^-1 and R
            ! note that since (LR)^-1 depends on R, we should R later
            
            !! update (LR)^-1
            ! obtain Umat=L*Delta^{(i)}
            Umat%blk1 = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( m, x, delta )
!$OMP DO
            do m = 1, ik
                x = xvec(m)
                delta = Dvec(m)
                Umat%blk1(:,m) = ul%blk1(:,x)*delta
            end do
!$OMP END DO
!$OMP END PARALLEL
            !
            ! obtain VUmat = (I + Vmat*Umat)^-1
            VUmat%blk1 = czero
            ! VUmat = Vmat*Umat = PFvecs * Delta^{(i)}P
!$OMP PARALLEL &
!$OMP PRIVATE ( n, x, delta )
!$OMP DO
            do n = 1, ik
                x = xvec(n)
                delta = Dvec(n)
                VUmat%blk1(:,n) = PFvecs(:,x)*delta
            end do
!$OMP END DO
!$OMP END PARALLEL
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
            ! initial xvec, Dvec and PFvecs
            xvec = czero
            Dvec = czero
            PFvecs = czero
            Vmat%blk1 = czero

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
