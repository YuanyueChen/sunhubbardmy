  subroutine dqmc_output_conf(this)
     use model_para
     implicit none
  
  	include 'mpif.h'
  
    class(vconf) :: this

    ! local
  	integer  status(mpi_status_size)
    integer, dimension(:,:,:),   allocatable ::  itmpv
    integer :: i, n, nt, nf
  
  	call mpi_comm_size(mpi_comm_world,isize,ierr)
  	call mpi_comm_rank(mpi_comm_world,irank,ierr)
  
    allocate (itmpv(this%nsites, this%nn_nf, this%ltrot) )
  	
    if (irank.eq.0) then
         open (unit=35, file='conf_v.out', status='unknown')
    endif
  
  	if ( irank.ne.0 )  then
  	   call mpi_send(this%conf, this%nsites*this%nn_nf*this%ltrot,mpi_integer, 0, irank+1536,mpi_comm_world,ierr)
  	endif
  
  	if (irank.eq.0)  then
         write(35,*) 1    ! useful for conf_v.in
  	     do nt = 1,this%ltrot
             do i  = 1,this%nsites
               do nf  = 1, this%nn_nf
                  write(35,'(i4)') this%conf(i,nf,nt)
               enddo
             enddo
         enddo
         do n = 1,isize - 1
  	      call mpi_recv(itmpv,this%nsites*this%nn_nf*this%ltrot, mpi_integer,n, n+1536, mpi_comm_world,status,ierr)
            do nt = 1,this%ltrot
            do i  = 1,this%nsites
            do nf  = 1, this%nn_nf
                write(35,'(i4)') itmpv(i,nf,nt)
            enddo
            enddo
            enddo
         enddo
    endif
    if (irank.eq.0) then
       close(35)
    endif
    deallocate(itmpv)
  end subroutine dqmc_output_conf
