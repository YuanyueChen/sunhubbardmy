  subroutine dqmc_input_vconf( this )
    use spring
    implicit none
  
  	include 'mpif.h'
    class(vconf) :: this
    ! local
  	integer  status(mpi_status_size)
    logical :: exists
    integer, dimension(:,:,:),   allocatable ::  itmpv
    integer :: iseed0, i, n, nt, eof, n_re, nf
    real(dp) :: x
  
  	call mpi_comm_size(mpi_comm_world,isize,ierr)
  	call mpi_comm_rank(mpi_comm_world,irank,ierr)
  
    allocate (itmpv(this%nsites,this%nn_nf,this%ltrot) )
  
    allocate( this%conf_v(this%nsites,this%nn_nf,this%ltrot) )
  
    if (irank .eq. 0 ) then
         inquire (file = 'conf_v.in', exist = exists)
         if ( exists .eqv. .true. ) then
#IFDEF BFIELD
             open (unit=30,file='conf_v.in', status='unknown', form='unformatted', access='sequential')
  	       read(30) iseed0
#ELSE
             open (unit=30,file='conf_v.in', status='unknown')
             read(30,*) iseed0
#ENDIF
         else
             iseed0 = 0
         end if
    endif
  
  	if ( irank.eq.0 )  then
  	   if (iseed0.eq.0) then
            ! start from scratch
            this%lwarmup = .true.
            write(fout,'(a)') ' start from scratch, need warnup '
  	      do n = 1,isize - 1
             ! setup node i and send data.
  	         do nt = 1,this%ltrot
                 do i  = 1, this%nsites
                    do nf  = 1, this%nn_nf
                        if( lprojv ) then
                           this%conf_v(i,nf,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                        else
                           this%conf_v(i,nf,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                        end if
                    end do
                 end do
             end do
  	         call mpi_send(this%conf_v, this%nsites*this%nn_nf*this%ltrot,mpi_integer, n, n+1536,mpi_comm_world,ierr)
  	      end do
            !	set node zero.
  	      do nt = 1,this%ltrot
              do i  = 1,this%nsites
                do nf  = 1, this%nn_nf
                    if( lprojv ) then
                        this%conf_v(i,nf,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                    else
                        this%conf_v(i,nf,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                    end if
                end do
  		      end do
          end do
  	   else
            ! read all config from node 0. 
            this%lwarmup = .false.
            write(fout,'(a)') ' start from old conf, do not need warnup '
            !	setup node 0
            do nt = 1,this%ltrot
              do i  = 1,this%nsites
                do nf  = 1, this%nn_nf
                  read(30,*,IOSTAT=eof) this%conf_v(i,nf,nt)
                  if(eof.lt.0) exit 
                end do
                if(eof.lt.0) exit 
              end do
              if(eof.lt.0) exit 
            end do
  	        do n = 1,isize - 1
               do nt = 1,this%ltrot
                 do i  = 1,this%nsites
                    do nf  = 1, this%nn_nf          
                      read(30,*,IOSTAT=eof) itmpv(i,nf,nt)
                      if(eof.lt.0) exit 
                    end do
                    if(eof.lt.0) exit 
                 end do
                 if(eof.lt.0) exit 
               end do
               if(eof.lt.0) exit 
               call mpi_send(itmpv, this%nsites*this%nn_nf*this%ltrot,mpi_integer, n, n+1536,mpi_comm_world,ierr)
            enddo
            ! if we do not have enough configurations, we have to copy configurations from master process
            if( eof .lt. 0 ) then
            write(fout, '(a)') '----------------------------------------------------------------------------------------------------------------!'
            write(fout, '(a)') '|| Note: In inconfc, we do not have enough configurations, we will copy configurations from master process!'
            write(fout, '(a)') '----------------------------------------------------------------------------------------------------------------!'
            do n_re = n, isize-1
                write(fout, '(a,i4,a)') ' In inconfc, irank = ', n_re, ' are copying configurations ... '
                itmpv = this%conf_v
                call mpi_send(itmpv, this%nsites*this%nn_nf*this%ltrot,mpi_integer, n_re, n_re+1536,mpi_comm_world,ierr)
            end do
            end if
            write(fout, '(a)') ' '
  	   endif
  	else
  	   call mpi_recv(this%conf_v, this%nsites*this%nn_nf*this%ltrot, mpi_integer,0,  irank + 1536, mpi_comm_world,status,ierr)
  	endif
  
      if (irank .eq. 0 ) then
         close(30)
      endif
  
      deallocate(itmpv)
  
  end subroutine dqmc_input_vconf
