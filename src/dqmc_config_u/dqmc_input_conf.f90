  subroutine dqmc_input_conf( this )
    use spring
    implicit none
  
  	include 'mpif.h'
    class(uconf) :: this
    ! local
  	integer  status(mpi_status_size)
    logical :: exists
    integer, dimension(:,:),   allocatable ::  itmpu
    integer :: iseed0, i, n, nt, eof, n_re
    real(dp) :: x
  
  	call mpi_comm_size(mpi_comm_world,isize,ierr)
  	call mpi_comm_rank(mpi_comm_world,irank,ierr)
  
    allocate (itmpu(this%nsites, this%ltrot) )
  
    allocate( this%conf(this%nsites, this%ltrot) )
  
    if (irank .eq. 0 ) then
         inquire (file = 'conf_u.in', exist = exists)
         if ( exists .eqv. .true. ) then
#IFDEF BFIELD
             open (unit=30,file='conf_u.in', status='unknown', form='unformatted', access='sequential')
  	       read(30) iseed0
#ELSE
             open (unit=30,file='conf_u.in', status='unknown')
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
                   if( lproju ) then
                       this%conf(i,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                   else
                       this%conf(i,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                   end if
                 end do
             end do
  	         call mpi_send(this%conf, this%nsites*this%ltrot,mpi_integer, n, n+1024,mpi_comm_world,ierr)
  	      end do
            !	set node zero.
  	      do nt = 1,this%ltrot
              do i  = 1,this%nsites
                  if( lproju ) then
                      this%conf(i,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                  else
                      this%conf(i,nt) = ceiling( spring_sfmt_stream() * this%lcomp )
                  end if
  		      end do
          end do
  	   else
            ! read all config from node 0. 
            this%lwarmup = .false.
            write(fout,'(a)') ' start from old conf, do not need warnup '
            !	setup node 0
            do nt = 1,this%ltrot
              do i  = 1,this%nsites
                  read(30,*,IOSTAT=eof) this%conf(i,nt)
                  if(eof.lt.0) exit 
              end do
              if(eof.lt.0) exit 
            end do
  	        do n = 1,isize - 1
               do nt = 1,this%ltrot
                 do i  = 1,this%nsites
                   read(30,*,IOSTAT=eof) itmpu(i,nt)
                   if(eof.lt.0) exit 
                 end do
                 if(eof.lt.0) exit 
               end do
               if(eof.lt.0) exit 
               call mpi_send(itmpu, this%nsites*this%ltrot,mpi_integer, n, n+1024,mpi_comm_world,ierr)
            enddo
            ! if we do not have enough configurations, we have to copy configurations from master process
            if( eof .lt. 0 ) then
            write(fout, '(a)') '----------------------------------------------------------------------------------------------------------------!'
            write(fout, '(a)') '|| Note: In inconfc, we do not have enough configurations, we will copy configurations from master process!'
            write(fout, '(a)') '----------------------------------------------------------------------------------------------------------------!'
            do n_re = n, isize-1
                write(fout, '(a,i4,a)') ' In inconfc, irank = ', n_re, ' are copying configurations ... '
                itmpu = this%conf
                call mpi_send(itmpu, this%nsites*this%ltrot,mpi_integer, n_re, n_re+1024,mpi_comm_world,ierr)
            end do
            end if
            write(fout, '(a)') ' '
  	   endif
  	else
  	   call mpi_recv(this%conf, this%nsites*this%ltrot, mpi_integer,0,  irank + 1024, mpi_comm_world,status,ierr)
  	endif
  
      if (irank .eq. 0 ) then
         close(30)
      endif
  
      deallocate(itmpu)
  
  end subroutine dqmc_input_conf
