  subroutine dqmc_output_uconf(this)
     use model_para
     implicit none
  
  	include 'mpif.h'
  
    class(uconf) :: this

    ! local
  	integer  status(mpi_status_size)
    integer, dimension(:,:),   allocatable ::  itmpu
    integer :: i, n, nt
  
  	call mpi_comm_size(mpi_comm_world,isize,ierr)
  	call mpi_comm_rank(mpi_comm_world,irank,ierr)
  
    allocate (itmpu(this%nsites, this%ltrot) )
  	
    if (irank.eq.0) then
         open (unit=35, file='conf_u.out', status='unknown')
    endif
  
  	if ( irank.ne.0 )  then
  	   call mpi_send(this%conf_u, this%nsites*this%ltrot,mpi_integer, 0, irank+1024,mpi_comm_world,ierr)
  	endif
  
  	if (irank.eq.0)  then
         write(35,*) 1    ! useful for conf_u.in
  	     do nt = 1,this%ltrot
             do i  = 1,this%nsites
                write(35,'(i4)') this%conf_u(i,nt)
             enddo
         enddo
         do n = 1,isize - 1
  	      call mpi_recv(itmpu,this%nsites*this%ltrot, mpi_integer,n, n+1024, mpi_comm_world,status,ierr)
            do nt = 1,this%ltrot
            do i  = 1,this%nsites
                write(35,'(i4)') itmpu(i,nt)
            enddo
            enddo
         enddo
    endif
    if (irank.eq.0) then
       close(35)
    endif
    deallocate(itmpu)
  end subroutine dqmc_output_uconf
