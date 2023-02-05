module honeycomb_lattice
  use constants, only: dp
  type honeycomb
     integer :: ncell, nsub, nsites, l1, l2, z_plq, nn_nf, nn_lf
     real(dp) :: a1_p(2), a2_p(2), l1_p(2), l2_p(2), b1_p(2), b2_p(2), bz1_p(2), bz2_p(2)
     integer, allocatable :: ulist(:,:), nnlf_list(:), list(:,:), cord(:,:), invlist(:,:,:), nnlist(:,:), snlist(:,:), tnlist(:,:), plq_cord(:,:), listk(:,:), invlistk(:,:),  imj(:,:), nnveci(:,:,:)
     real(dp), allocatable :: nnvecr(:,:,:), rcord(:,:)
     complex(dp), allocatable :: zexpiqr(:,:)
     contains
       private
       final :: deallocate_honeycomb
  end type honeycomb

  interface iscalar
     module procedure iscalar_ii, iscalar_ir, iscalar_rr
  end interface

  interface xnorm
     module procedure xnorm_i, xnorm_r
  end interface

  contains

  subroutine deallocate_honeycomb( this )
    implicit none
    type(honeycomb) :: this
    if( allocated(this%ulist) ) deallocate( this%ulist )
    if( allocated(this%nnlf_list) ) deallocate( this%nnlf_list )
    if( allocated(this%list) ) deallocate( this%list )
    if( allocated(this%cord) ) deallocate( this%cord )
    if( allocated(this%invlist) ) deallocate( this%invlist )
    if( allocated(this%nnlist) ) deallocate( this%nnlist )
    if( allocated(this%snlist) ) deallocate( this%snlist )
    if( allocated(this%tnlist) ) deallocate( this%tnlist )
    if( allocated(this%plq_cord) ) deallocate( this%plq_cord )
    if( allocated(this%listk) ) deallocate( this%listk )
    if( allocated(this%invlistk) ) deallocate( this%invlistk )
    if( allocated(this%imj) ) deallocate( this%imj )
    if( allocated(this%nnveci) ) deallocate( this%nnveci )
    if( allocated(this%nnvecr) ) deallocate( this%nnvecr )
    if( allocated(this%rcord) ) deallocate( this%rcord )
    if( allocated(this%zexpiqr) ) deallocate( this%zexpiqr )
    write(*,*) " deallocate_honeycomb is performed "
  end subroutine deallocate_honeycomb
  
  subroutine setup_honeycomb(l1, l2, a1_p, a2_p, latt) 

    ! set up honeycomb lattice
    
    implicit real (kind=8) (a-g,o-z)
    implicit integer (h-n)
    
    integer, intent(in) :: l1, l2
    real(dp), dimension(:), intent(in) :: a1_p, a2_p
    type (honeycomb) :: latt

    ! local
    real(dp) :: ai_p(2), xk_p(2), b1_p(2), b2_p(2), bz1_p(2), bz2_p(2), b_p(2)
    real(dp) :: mat(2,2), mat_inv(2,2)

    latt%nsub = 2
    latt%l1 = l1
    latt%l2 = l2
    latt%ncell = l1*l2
    latt%nsites = l1*l2*2
    latt%z_plq = 6
    latt%nn_nf = 3
    latt%nn_lf = l1*l2
    allocate( latt%cord( latt%nsites,2) )
    allocate( latt%rcord(latt%nsites,2) )
    allocate( latt%ulist( latt%ncell,2) )
    allocate( latt%nnlf_list( latt%nn_lf) )
    allocate( latt%list( latt%nsites,2) )
    allocate( latt%invlist(l1, l2, 2) )
    zero = 1.e-6
    latt%l1_p = l1*a1_p
    latt%l2_p = l2*a2_p
    latt%a1_p = a1_p
    latt%a2_p = a2_p

    do i = 1, latt%nn_lf
        latt%nnlf_list(i) = 2*i - 1
    end do

    ind = 0
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            latt%ulist( i2+(i1-1)*l2, 1) = i1
            latt%ulist( i2+(i1-1)*l2, 2) = i2
            do iorb = 1, latt%nsub
                ind = ind + 1
                latt%list(ind, 1) = i2+(i1-1)*l2 ! unit cell
                latt%list(ind, 2) = iorb         ! orb
                latt%cord(ind, 1) = i1
                latt%cord(ind, 2) = i2
                latt%invlist(i1,i2,iorb) = ind
                if( iorb == 1 ) then ! A site
                    latt%rcord(ind, 1) = dble(i1)
                    latt%rcord(ind, 2) = dble(i2)
                else ! B site
                    latt%rcord(ind, 1) = dble(i1) - 1.d0/3.d0
                    latt%rcord(ind, 2) = dble(i2) + 1.d0/3.d0
                end if
            end do
        end do
    end do
!                     A 
!                 *  /    * 
!              B    /        B
!                  a2         
!              *  /   x      *
!                /                                        three kinds of plaq
!              A ----------- A                                                  
!           * / \ *       *     *                               / \ / \
!       B    /   \    B            B                           | B | C |                   |   |   |
!           /     a1                                          / \ / \ /                   a2                     B site
!       *  /       \  *            *                         | C | A |  the lattice:     |   |   |               |
!         /         \                                       / \ / \ /                     a1                     A site
!       A             A            A                       | A | B |                   |   |   |
!           *     *       *     *                           \ / \ /                            
!              B             B

! set nearest neighbor
    allocate( latt%nnvecr(2,3,latt%nsites) )
    allocate( latt%nnveci(2,3,latt%nsites) )
    allocate( latt%nnlist(latt%nsites,3) )
    allocate( latt%snlist(latt%nsites,6) )
    allocate( latt%tnlist(latt%nsites,3) )
    ! A->B
    do i = 1, latt%nsites, 2
        latt%nnvecr(:,1,i) = (/-1.d0/3.d0, 1.d0/3.d0/) ! in the unit of a1_p and a2_p
        latt%nnvecr(:,2,i) = (/ 2.d0/3.d0, 1.d0/3.d0/)
        latt%nnvecr(:,3,i) = (/-1.d0/3.d0,-2.d0/3.d0/)
        latt%nnveci(:,1,i) = (/ 0, 0 /)
        latt%nnveci(:,2,i) = (/ 1, 0 /)
        latt%nnveci(:,3,i) = (/ 0,-1 /)
    end do
    ! B->A
    do i = 2, latt%nsites, 2
        latt%nnvecr(:,1,i) = (/ 1.d0/3.d0,-1.d0/3.d0/) ! in the unit of a1_p and a2_p
        latt%nnvecr(:,2,i) = (/-2.d0/3.d0,-1.d0/3.d0/)
        latt%nnvecr(:,3,i) = (/ 1.d0/3.d0, 2.d0/3.d0/)
        latt%nnveci(:,1,i) = (/ 0, 0 /)
        latt%nnveci(:,2,i) = (/-1, 0 /)
        latt%nnveci(:,3,i) = (/ 0, 1 /)
    end do
    do i  = 1, latt%ncell*latt%nsub
        if( mod(i,2) == 1 ) then ! A site
            latt%nnlist(i,1) = i+1
            latt%nnlist(i,2) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)  , l2 ) )*latt%nsub
            latt%nnlist(i,3) = ( (npbc(latt%cord(i,1)  , l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub
            latt%snlist(i,1) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub - 1
            latt%snlist(i,2) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)  , l2 ) )*latt%nsub - 1
            latt%snlist(i,3) = ( (npbc(latt%cord(i,1)  , l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub - 1
            latt%snlist(i,4) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub - 1
            latt%snlist(i,5) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)  , l2 ) )*latt%nsub - 1
            latt%snlist(i,6) = ( (npbc(latt%cord(i,1)  , l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub - 1
            latt%tnlist(i,1) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub
            latt%tnlist(i,2) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub
            latt%tnlist(i,3) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub
        else                     ! B site
            latt%nnlist(i,1) = i-1
            latt%nnlist(i,2) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)  , l2 ) )*latt%nsub - 1
            latt%nnlist(i,3) = ( (npbc(latt%cord(i,1)  , l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub - 1
            latt%snlist(i,1) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub
            latt%snlist(i,2) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)  , l2 ) )*latt%nsub
            latt%snlist(i,3) = ( (npbc(latt%cord(i,1)  , l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub
            latt%snlist(i,4) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub
            latt%snlist(i,5) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)  , l2 ) )*latt%nsub
            latt%snlist(i,6) = ( (npbc(latt%cord(i,1)  , l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub
            latt%tnlist(i,1) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub - 1
            latt%tnlist(i,2) = ( (npbc(latt%cord(i,1)+1, l1) - 1)*l2 + npbc(latt%cord(i,2)+1, l2 ) )*latt%nsub - 1
            latt%tnlist(i,3) = ( (npbc(latt%cord(i,1)-1, l1) - 1)*l2 + npbc(latt%cord(i,2)-1, l2 ) )*latt%nsub - 1
        end if
    end do

! set plq_cord
    allocate( latt%plq_cord(6,latt%ncell))
    ! count number of A,B,C type plaq first
    id1 = 0; id2 = 0; id3 = 0
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            !i = i2 + (i1-1)*latt%l2
            if( mod(i1+i2,3) == 2 ) then ! A type plaq
                id1 = id1 + 1
            else if( mod(i1+i2,3) == 0 ) then ! C type plaq
                id3 = id3 + 1
            else ! B type plaq
                id2 = id2 + 1
            end if
        end do
    end do
    nd1 = id1; nd2 = id2; nd3 = id3
    id1 = 0; id2 = 0; id3 = 0
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            !i = i2 + (i1-1)*latt%l2
            if( mod(i1+i2,3) == 2 ) then ! A type plaq
                id1 = id1 + 1
                id = id1 
            else if( mod(i1+i2,3) == 0 ) then ! C type plaq
                id3 = id3 + 1
                id = id3 + nd1 + nd2
            else ! B type plaq
                id2 = id2 + 1
                id = id2 + nd1
            end if
            ! with clockwise order, 3-clock-direction is the first one
            latt%plq_cord(1, id) = latt%invlist(i1, i2, 2)
            latt%plq_cord(2, id) = latt%invlist(i1, i2, 1)
            latt%plq_cord(3, id) = latt%invlist(i1, npbc(i2-1,l2), 2)
            latt%plq_cord(4, id) = latt%invlist(npbc(i1-1,l1), npbc(i2-1,l2), 1)
            latt%plq_cord(5, id) = latt%invlist(npbc(i1-1,l1), npbc(i2-1,l2), 2)
            latt%plq_cord(6, id) = latt%invlist(npbc(i1-1,l1), i2, 1)
        end do
    end do
    
    !setting up k-space lattice
    allocate ( latt%listk(latt%ncell,2), latt%invlistk(-l1:l1, -l2:l2) )
    pi   = dacos(-1.d0)
    ! setup the 2x2 matrix to determine  bz1_p, bz2_p
    mat(1,1) = dble(a1_p(1))
    mat(1,2) = dble(a1_p(2))
    mat(2,1) = dble(a2_p(1))
    mat(2,2) = dble(a2_p(2))
    x = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)
    mat_inv(1,1) =  mat(2,2)/x
    mat_inv(2,2) =  mat(1,1)/x
    mat_inv(1,2) = -mat(1,2)/x
    mat_inv(2,1) = -mat(2,1)/x
    bz1_p(1)      = 2.d0*pi*mat_inv(1,1)
    bz1_p(2)      = 2.d0*pi*mat_inv(2,1)
    bz2_p(1)      = 2.d0*pi*mat_inv(1,2)
    bz2_p(2)      = 2.d0*pi*mat_inv(2,2)
    latt%bz1_p = bz1_p
    latt%bz2_p = bz2_p
    

    x =  2.d0*pi / ( iscalar(bz1_p,latt%l1_p) * iscalar(bz2_p,latt%l2_p) -   &
         &           iscalar(bz2_p,latt%l1_p) * iscalar(bz1_p,latt%l2_p)   )
    x = abs(x)

    b1_p = x*( iscalar(bz2_p,latt%l2_p) * bz1_p - iscalar(bz1_p,latt%l2_p) * bz2_p ) 
    b2_p = x*( iscalar(bz1_p,latt%l1_p) * bz2_p - iscalar(bz2_p,latt%l1_p) * bz1_p )


    latt%b1_p  = b1_p
    latt%b2_p  = b2_p

    !write(6,*) 'bz1: ', bz1_p(1), bz1_p(2)
    !write(6,*) 'bz2: ', bz2_p(1), bz2_p(2)
    !write(6,*) 'b1 : ', b1_p(1), b1_p(2)
    !write(6,*) 'b2 : ', b2_p(1), b2_p(2)
    !!!nc = 0
    !!!do m = -l1,l1
    !!!   do n = -l2,l2
    !!!      xk_p = dble(m)* b1_p + dble(n) * b2_p
    !!!      l_f = 1
    !!!      do i = 1,4
    !!!         if (i.eq.1) b_p = bz2_p 
    !!!         if (i.eq.2) b_p = bz1_p  
    !!!         if (i.eq.3) b_p = bz2_p - bz1_p 
    !!!         if (i.eq.4) b_p = bz2_p + bz1_p
    !!!         if  (  iscalar( xk_p, b_p )    .le.  xnorm(b_p)**2/2.d0 + zero   .and.   &
    !!!              & iscalar( xk_p, b_p )    .ge. -xnorm(b_p)**2/2.d0 + zero    ) then
    !!!            l_f = l_f * 1
    !!!         else
    !!!            l_f = 0
    !!!         endif
    !!!      enddo
    !!!      if (l_f .eq. 1) then   
    !!!         !write(11,"(f14.7,2x,f14.7)")  xk_p(1), xk_p(2)
    !!!         nc = nc + 1
    !!!         latt%listk(nc,1) = m
    !!!         latt%listk(nc,2) = n
    !!!         latt%invlistk(m,n) = nc
    !!!      endif
    !!!   enddo
    !!!enddo
    nc = 0
    do m = -(l1+1)/2+1, l1/2 
      do n = -(l2+1)/2+1, l2/2
          nc = nc + 1
          latt%listk(nc,1) = m
          latt%listk(nc,2) = n
          latt%invlistk(m,n) = nc
      end do
    end do
    if (nc.ne.latt%ncell) then 
       write(6,*) 'error ', nc, latt%ncell
       stop
    endif

    !setup imj 
    allocate ( latt%imj(latt%ncell,latt%ncell) )
    ! latt_imj
    do j = 1, latt%ncell
        do i = 1, latt%ncell
            imj_nx = npbc( latt%cord(2*i,1) - latt%cord(2*j,1), latt%l1 )
            imj_ny = npbc( latt%cord(2*i,2) - latt%cord(2*j,2), latt%l2 )
            latt%imj(i,j) = imj_ny + (imj_nx-1)*latt%l2
        end do
    end do

    allocate( latt%zexpiqr(latt%ncell,latt%ncell) )
    do i = 1,latt%ncell
       ai_p = dble(latt%ulist(i,1))*latt%a1_p + dble(latt%ulist(i,2))*latt%a2_p
       do nk = 1,latt%ncell
          xk_p = dble(latt%listk(nk,1))*latt%b1_p + dble(latt%listk(nk,2))*latt%b2_p
          latt%zexpiqr(nk,i) = exp( dcmplx( 0.d0, iscalar(xk_p, ai_p)) )
       enddo
    enddo
    
  end subroutine setup_honeycomb

  integer function npbc(nr,l)
    implicit none
    integer, intent(in) :: nr
    integer, intent(in) :: l
    npbc = nr
    if (nr.gt.l) npbc = nr - l
    if (nr.lt.1) npbc = nr + l
  end function npbc

  integer function inv_k(xk_p,latt) 
    
    implicit none
    real (kind=8)  :: xk_p(2)
    type (honeycomb) :: latt
    
    real (kind=8)     :: xk1_p(2), x, zero
    integer :: nk 

    zero = 1.d-10
    nk = 1
    do 
       xk1_p = latt%listk(nk,1)*latt%b1_p + latt%listk(nk,2)*latt%b2_p
       if (xnorm(xk1_p - xk_p)  < zero ) then
          inv_k = nk
          exit
       elseif (nk < latt%ncell) then 
          nk = nk + 1
       else
          write(6,*) 'error in inv_k lattice_new'
          stop
       endif
    enddo

  end function inv_k
  
  integer function iscalar_ii(i_p, j_p)
    implicit none
    integer, dimension(:) :: i_p, j_p
    integer i
    
    iscalar_ii = 0
    !write(6,*) size(i_p)
    do i = 1,  size(i_p)
      ! write(6,*) i
       iscalar_ii = iscalar_ii + i_p(i)*j_p(i)
    enddo
  end function iscalar_ii

  real (kind=8)  function iscalar_ir(x_p, j_p)
    implicit none
    real (kind=8), dimension(:) ::  x_p
    integer, dimension(:) ::  j_p
    integer i
    
    iscalar_ir = 0.d0
    !write(6,*) size(i_p)
    do i = 1,  size(x_p)
      ! write(6,*) i
       iscalar_ir = iscalar_ir + x_p(i)*dble(j_p(i))
    enddo
  end function iscalar_ir

  real (kind=8)  function iscalar_rr(x_p, y_p)
    implicit none
    real (kind=8), dimension(:) ::  x_p, y_p
    integer i
    
    iscalar_rr = 0.d0
    do i = 1,  size(x_p)
       iscalar_rr = iscalar_rr + x_p(i)*y_p(i)
    enddo
  end function iscalar_rr

  real (kind=8) function xnorm_i(i_p)
    implicit none
    integer, dimension(:) :: i_p
    integer :: i

    xnorm_i = 0.d0
    do i = 1,  size(i_p)
       xnorm_i = xnorm_i + dble(i_p(i)*i_p(i))
    enddo
    xnorm_i = sqrt(xnorm_i)
  end function xnorm_i

  real (kind=8) function xnorm_r(x_p)
    implicit none
    real (kind=8), dimension(:) :: x_p
    integer :: i

    xnorm_r = 0.d0
    do i = 1,  size(x_p)
       xnorm_r = xnorm_r + x_p(i)*x_p(i)
    enddo
    xnorm_r = sqrt(xnorm_r)
  end function xnorm_r

  subroutine print_latt(latt)
    
    implicit real (kind=8) (a-g,o-z)
    implicit integer (h-n)
    
    type (honeycomb) :: latt
    integer          :: i_p(2),nd_p(2)
    real    (kind=8) :: x_p(2)

    open (unit=55,file="latt_info", status = "unknown")
    write(55,*) ' reciprocal vector 1: ', latt%bz1_p(1), latt%bz1_p(2)
    write(55,*) ' reciprocal vector 2: ', latt%bz2_p(1), latt%bz2_p(2)
    write(55,*) ' latt       vector 1: ', latt%a1_p(1), latt%a1_p(2)
    write(55,*) ' latt       vector 2: ', latt%a2_p(1), latt%a2_p(2)
    close(55)
  end subroutine print_latt

end module honeycomb_lattice


module square_lattice
  use constants, only: dp
  type square
     integer :: ncell, nsub, nsites, l1, l2, z_plq, nn_nf, nn_lf
     real(dp) :: a1_p(2), a2_p(2), l1_p(2), l2_p(2), b1_p(2), b2_p(2), bz1_p(2), bz2_p(2)
     integer, allocatable :: ulist(:,:), nnlf_list(:), list(:,:), cord(:,:), invlist(:,:), nnlist(:,:), snlist(:,:), tnlist(:,:), plq_cord(:,:), listk(:,:), invlistk(:,:), imj(:,:), nnveci(:,:,:)
     real(dp), allocatable :: nnvecr(:,:,:), rcord(:,:)
     complex(dp), allocatable :: zexpiqr(:,:)
     contains
       private
       final :: deallocate_square
  end type square

  interface iscalar
     module procedure iscalar_ii, iscalar_ir, iscalar_rr
  end interface

  interface xnorm
     module procedure xnorm_i, xnorm_r
  end interface

  contains

  subroutine deallocate_square( this )
    implicit none
    type(square) :: this
    if( allocated(this%ulist) ) deallocate( this%ulist )
    if( allocated(this%nnlf_list) ) deallocate( this%nnlf_list )
    if( allocated(this%list) ) deallocate( this%list )
    if( allocated(this%cord) ) deallocate( this%cord )
    if( allocated(this%invlist) ) deallocate( this%invlist )
    if( allocated(this%nnlist) ) deallocate( this%nnlist )
    if( allocated(this%snlist) ) deallocate( this%snlist )
    if( allocated(this%tnlist) ) deallocate( this%tnlist )
    if( allocated(this%plq_cord) ) deallocate( this%plq_cord )
    if( allocated(this%listk) ) deallocate( this%listk )
    if( allocated(this%invlistk) ) deallocate( this%invlistk )
    if( allocated(this%imj) ) deallocate( this%imj )
    if( allocated(this%nnveci) ) deallocate( this%nnveci )
    if( allocated(this%nnvecr) ) deallocate( this%nnvecr )
    if( allocated(this%rcord) ) deallocate( this%rcord )
    if( allocated(this%zexpiqr) ) deallocate( this%zexpiqr )
    write(*,*) " deallocate_square is performed "
  end subroutine deallocate_square
  
  subroutine setup_square(l1, l2, a1_p, a2_p, latt) 

    ! set up square lattice
    
    implicit real (kind=8) (a-g,o-z)
    implicit integer (h-n)
    
    integer, intent(in) :: l1, l2
    real(dp), dimension(:), intent(in) :: a1_p, a2_p
    type (square) :: latt

    ! local
    real(dp) :: ai_p(2), xk_p(2), b1_p(2), b2_p(2), bz1_p(2), bz2_p(2), b_p(2)
    real(dp) :: mat(2,2), mat_inv(2,2)

    latt%nsub = 1
    latt%l1 = l1
    latt%l2 = l2
    latt%ncell = l1*l2
    latt%nsites = l1*l2
    latt%z_plq = 4
    latt%nn_nf = 4
    latt%nn_lf = l1*l2/2
    allocate( latt%rcord( latt%nsites,2) )
    allocate( latt%ulist( latt%ncell,2) )
    allocate( latt%nnlf_list( latt%nn_lf) )
    allocate( latt%list( latt%nsites,2) )
    allocate( latt%invlist(l1, l2) )
    zero = 1.e-6
    latt%l1_p = l1*a1_p
    latt%l2_p = l2*a2_p
    latt%a1_p = a1_p
    latt%a2_p = a2_p

    ind = 0
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            latt%ulist( i2+(i1-1)*l2, 1) = i1
            latt%ulist( i2+(i1-1)*l2, 2) = i2
            latt%list( i2+(i1-1)*l2, 1) = i2+(i1-1)*l2
            latt%list( i2+(i1-1)*l2, 2) = 1
            latt%rcord( i2+(i1-1)*l2, 1) = dble(i1)
            latt%rcord( i2+(i1-1)*l2, 2) = dble(i2)
            latt%invlist(i1,i2) = i2+(i1-1)*l2
            if( mod(i1+i2,2) == 0 ) then
                ind = ind + 1
                latt%nnlf_list(ind) = i2+(i1-1)*l2
            end if
        end do
    end do
    latt%cord = latt%ulist
!                                                                               
!              A --------- A --------- A                                                  
!              |           |           |                                         
!              |           |           |                                                                      
!              |           a2          |                                                                                         
!              |           |           |                                                                                    
!              |           |           |                                                                                         
!              A --------- A ----a1--- A                                                                  
!              |           |           |                                         
!              |           |           |                                                                      
!              |           |           |                                                                                         
!              |           |           |                                                                                    
!              |           |           |                                                                                         
!              A --------- A --------- A                                                                  
!                                                                                              

! set nearest neighbor
    allocate( latt%nnvecr(2,4,latt%nsites) )
    allocate( latt%nnveci(2,4,latt%nsites) )
    allocate( latt%nnlist(latt%nsites,4) )
    allocate( latt%snlist(latt%nsites,4) )
    allocate( latt%tnlist(latt%nsites,4) )
    ! A->B
    do i = 1, latt%nsites, 2
        latt%nnvecr(:,1,i) = (/ 1.d0, 0.d0/) ! in the unit of a1_p and a2_p
        latt%nnvecr(:,2,i) = (/ 0.d0, 1.d0/)
        latt%nnvecr(:,3,i) = (/-1.d0, 0.d0/)
        latt%nnvecr(:,4,i) = (/ 0.d0,-1.d0/)

        latt%nnveci(:,1,i) = (/ 1, 0 /)
        latt%nnveci(:,2,i) = (/ 0, 1 /)
        latt%nnveci(:,3,i) = (/-1, 0 /)
        latt%nnveci(:,4,i) = (/ 0,-1 /)
    end do
    ! B->A
    do i = 2, latt%nsites, 2
        latt%nnvecr(:,1,i) = (/ 1.d0, 0.d0/) ! in the unit of a1_p and a2_p
        latt%nnvecr(:,2,i) = (/ 0.d0, 1.d0/)
        latt%nnvecr(:,3,i) = (/-1.d0, 0.d0/)
        latt%nnvecr(:,4,i) = (/ 0.d0,-1.d0/)
        latt%nnveci(:,1,i) = (/ 1, 0 /)
        latt%nnveci(:,2,i) = (/ 0, 1 /)
        latt%nnveci(:,3,i) = (/-1, 0 /)
        latt%nnveci(:,4,i) = (/ 0,-1 /)
    end do
    do i  = 1, latt%ncell
        latt%nnlist(i,1) = ( npbc(latt%ulist(i,1)+1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)  , l2)
        latt%nnlist(i,2) = ( npbc(latt%ulist(i,1)  ,l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2)
        latt%nnlist(i,3) = ( npbc(latt%ulist(i,1)-1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)  , l2)
        latt%nnlist(i,4) = ( npbc(latt%ulist(i,1)  ,l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2)

        latt%snlist(i,1) = ( npbc(latt%ulist(i,1)+1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2)
        latt%snlist(i,2) = ( npbc(latt%ulist(i,1)-1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2)
        latt%snlist(i,3) = ( npbc(latt%ulist(i,1)-1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2)
        latt%snlist(i,4) = ( npbc(latt%ulist(i,1)+1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2)

        latt%tnlist(i,1) = ( npbc(latt%ulist(i,1)+2,l1) - 1 )*l2 + npbc(latt%ulist(i,2)  , l2)
        latt%tnlist(i,2) = ( npbc(latt%ulist(i,1)  ,l1) - 1 )*l2 + npbc(latt%ulist(i,2)+2, l2)
        latt%tnlist(i,3) = ( npbc(latt%ulist(i,1)-2,l1) - 1 )*l2 + npbc(latt%ulist(i,2)  , l2)
        latt%tnlist(i,4) = ( npbc(latt%ulist(i,1)  ,l1) - 1 )*l2 + npbc(latt%ulist(i,2)-2, l2)
    end do

    ! set plq_cord
    allocate( latt%plq_cord(4,latt%ncell))
    ! count number of A,B,C,D type plaq first
    !
    !  ---------------------------------
    !  |       |       |       |       |
    !  |   D   |   C   |   D   |   C   |
    !  |       |       |       |       |
    !  ---------------------------------
    !  |       |       |       |       |
    !  |   A   |   B   |   A   |   B   |
    !  |       |       |       |       |
    !  ---------------------------------
    !  |       |       |       |       |
    !  |   D   |   C   |   D   |   C   |
    !  |       |       |       |       |
    !  ---------------------------------
    !  |       |       |       |       |
    !  |   A   |   B   |   A   |   B   |
    !  |       |       |       |       |
    !  ---------------------------------
    id1 = 0; id2 = 0; id3 = 0; id4 = 0;
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            !i = i2 + (i1-1)*latt%l2
            if( mod(i1,2) == 1 .and. mod(i2,2) == 1 ) then ! A type plaq
                id1 = id1 + 1
            else if( mod(i1,2) == 0 .and. mod(i2,2) == 1 ) then ! B type plaq
                id2 = id2 + 1
            else if( mod(i1,2) == 0 .and. mod(i2,2) == 0 ) then ! C type plaq
                id3 = id3 + 1
            else ! D type plaq
                id4 = id4 + 1
            end if
        end do
    end do
    nd1 = id1; nd2 = id2; nd3 = id3; nd4 = id4;
    id1 = 0; id2 = 0; id3 = 0; id4 = 0;
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            !i = i2 + (i1-1)*latt%l2
            if( mod(i1,2) == 1 .and. mod(i2,2) == 1 ) then ! A type plaq
                id1 = id1 + 1
                id = id1
            else if( mod(i1,2) == 0 .and. mod(i2,2) == 1 ) then ! B type plaq
                id2 = id2 + 1
                id = id2 + nd1
            else if( mod(i1,2) == 0 .and. mod(i2,2) == 0 ) then ! C type plaq
                id3 = id3 + 1
                id = id3 + nd1 + nd2
            else ! D type plaq
                id4 = id4 + 1
                id = id4 + nd1 + nd2 + nd3
            end if

            ! with anti-clockwise order, left-bottom corner is the first one
            latt%plq_cord(1, id) = latt%invlist(i1, i2)
            latt%plq_cord(2, id) = latt%invlist(npbc(i1+1,l1), i2)
            latt%plq_cord(3, id) = latt%invlist(npbc(i1+1,l1), npbc(i2+1,l2))
            latt%plq_cord(4, id) = latt%invlist(i1, npbc(i2+1,l2))
        end do
    end do

    !setting up k-space lattice
    allocate ( latt%listk(latt%ncell,2), latt%invlistk(-l1:l1, -l2:l2) )
    pi   = dacos(-1.d0)
    ! setup the 2x2 matrix to determine  bz1_p, bz2_p
    mat(1,1) = dble(a1_p(1))
    mat(1,2) = dble(a1_p(2))
    mat(2,1) = dble(a2_p(1))
    mat(2,2) = dble(a2_p(2))
    x = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)
    mat_inv(1,1) =  mat(2,2)/x
    mat_inv(2,2) =  mat(1,1)/x
    mat_inv(1,2) = -mat(1,2)/x
    mat_inv(2,1) = -mat(2,1)/x
    bz1_p(1)      = 2.d0*pi*mat_inv(1,1)
    bz1_p(2)      = 2.d0*pi*mat_inv(2,1)
    bz2_p(1)      = 2.d0*pi*mat_inv(1,2)
    bz2_p(2)      = 2.d0*pi*mat_inv(2,2)
    latt%bz1_p = bz1_p
    latt%bz2_p = bz2_p
    

    x =  2.d0*pi / ( iscalar(bz1_p,latt%l1_p) * iscalar(bz2_p,latt%l2_p) -   &
         &           iscalar(bz2_p,latt%l1_p) * iscalar(bz1_p,latt%l2_p)   )
    x = abs(x)

    b1_p = x*( iscalar(bz2_p,latt%l2_p) * bz1_p - iscalar(bz1_p,latt%l2_p) * bz2_p ) 
    b2_p = x*( iscalar(bz1_p,latt%l1_p) * bz2_p - iscalar(bz2_p,latt%l1_p) * bz1_p )


    latt%b1_p  = b1_p
    latt%b2_p  = b2_p

    !write(6,*) 'bz1: ', bz1_p(1), bz1_p(2)
    !write(6,*) 'bz2: ', bz2_p(1), bz2_p(2)
    !write(6,*) 'b1 : ', b1_p(1), b1_p(2)
    !write(6,*) 'b2 : ', b2_p(1), b2_p(2)
    !!!nc = 0
    !!!do m = -l1,l1
    !!!   do n = -l2,l2
    !!!      xk_p = dble(m)* b1_p + dble(n) * b2_p
    !!!      l_f = 1
    !!!      do i = 1,4
    !!!         if (i.eq.1) b_p = bz2_p 
    !!!         if (i.eq.2) b_p = bz1_p  
    !!!         if (i.eq.3) b_p = bz2_p - bz1_p 
    !!!         if (i.eq.4) b_p = bz2_p + bz1_p
    !!!         if  (  iscalar( xk_p, b_p )    .le.  xnorm(b_p)**2/2.d0 + zero   .and.   &
    !!!              & iscalar( xk_p, b_p )    .ge. -xnorm(b_p)**2/2.d0 + zero    ) then
    !!!            l_f = l_f * 1
    !!!         else
    !!!            l_f = 0
    !!!         endif
    !!!      enddo
    !!!      if (l_f .eq. 1) then   
    !!!         !write(11,"(f14.7,2x,f14.7)")  xk_p(1), xk_p(2)
    !!!         nc = nc + 1
    !!!         latt%listk(nc,1) = m
    !!!         latt%listk(nc,2) = n
    !!!         latt%invlistk(m,n) = nc
    !!!      endif
    !!!   enddo
    !!!enddo
    nc = 0
    do m = -(l1+1)/2+1, l1/2 
      do n = -(l2+1)/2+1, l2/2
          nc = nc + 1
          latt%listk(nc,1) = m
          latt%listk(nc,2) = n
          latt%invlistk(m,n) = nc
      end do
    end do
    if (nc.ne.latt%ncell) then 
       write(6,*) 'error ', nc, latt%ncell
       stop
    endif

    !setup imj 
    allocate ( latt%imj(latt%ncell,latt%ncell) )
    ! latt_imj
    do j = 1, latt%ncell
        do i = 1, latt%ncell
            imj_nx = npbc( latt%ulist(i,1) - latt%ulist(j,1), latt%l1 )
            imj_ny = npbc( latt%ulist(i,2) - latt%ulist(j,2), latt%l2 )
            latt%imj(i,j) = imj_ny + (imj_nx-1)*latt%l2
        end do
    end do

    allocate( latt%zexpiqr(latt%ncell,latt%ncell) )
    do i = 1,latt%ncell
       ai_p = dble(latt%ulist(i,1))*latt%a1_p + dble(latt%ulist(i,2))*latt%a2_p
       do nk = 1,latt%ncell
          xk_p = dble(latt%listk(nk,1))*latt%b1_p + dble(latt%listk(nk,2))*latt%b2_p
          latt%zexpiqr(nk,i) = exp( dcmplx( 0.d0, iscalar(xk_p, ai_p)) )
       enddo
    enddo
    
  end subroutine setup_square

  integer function npbc(nr,l)
    implicit none
    integer, intent(in) :: nr
    integer, intent(in) :: l
    npbc = nr
    if (nr.gt.l) npbc = nr - l
    if (nr.lt.1) npbc = nr + l
  end function npbc

  integer function inv_k(xk_p,latt) 
    
    implicit none
    real (kind=8)  :: xk_p(2)
    type (square) :: latt
    
    real (kind=8)     :: xk1_p(2), x, zero
    integer :: nk 

    zero = 1.d-10
    nk = 1
    do 
       xk1_p = latt%listk(nk,1)*latt%b1_p + latt%listk(nk,2)*latt%b2_p
       if (xnorm(xk1_p - xk_p)  < zero ) then
          inv_k = nk
          exit
       elseif (nk < latt%ncell) then 
          nk = nk + 1
       else
          write(6,*) 'error in inv_k lattice_new'
          stop
       endif
    enddo

  end function inv_k
  
  integer function iscalar_ii(i_p, j_p)
    implicit none
    integer, dimension(:) :: i_p, j_p
    integer i
    
    iscalar_ii = 0
    !write(6,*) size(i_p)
    do i = 1,  size(i_p)
      ! write(6,*) i
       iscalar_ii = iscalar_ii + i_p(i)*j_p(i)
    enddo
  end function iscalar_ii

  real (kind=8)  function iscalar_ir(x_p, j_p)
    implicit none
    real (kind=8), dimension(:) ::  x_p
    integer, dimension(:) ::  j_p
    integer i
    
    iscalar_ir = 0.d0
    !write(6,*) size(i_p)
    do i = 1,  size(x_p)
      ! write(6,*) i
       iscalar_ir = iscalar_ir + x_p(i)*dble(j_p(i))
    enddo
  end function iscalar_ir

  real (kind=8)  function iscalar_rr(x_p, y_p)
    implicit none
    real (kind=8), dimension(:) ::  x_p, y_p
    integer i
    
    iscalar_rr = 0.d0
    do i = 1,  size(x_p)
       iscalar_rr = iscalar_rr + x_p(i)*y_p(i)
    enddo
  end function iscalar_rr

  real (kind=8) function xnorm_i(i_p)
    implicit none
    integer, dimension(:) :: i_p
    integer :: i

    xnorm_i = 0.d0
    do i = 1,  size(i_p)
       xnorm_i = xnorm_i + dble(i_p(i)*i_p(i))
    enddo
    xnorm_i = sqrt(xnorm_i)
  end function xnorm_i

  real (kind=8) function xnorm_r(x_p)
    implicit none
    real (kind=8), dimension(:) :: x_p
    integer :: i

    xnorm_r = 0.d0
    do i = 1,  size(x_p)
       xnorm_r = xnorm_r + x_p(i)*x_p(i)
    enddo
    xnorm_r = sqrt(xnorm_r)
  end function xnorm_r

  subroutine print_latt(latt)
    
    implicit real (kind=8) (a-g,o-z)
    implicit integer (h-n)
    
    type (square) :: latt
    integer          :: i_p(2),nd_p(2)
    real    (kind=8) :: x_p(2)

    open (unit=55,file="latt_info", status = "unknown")
    write(55,*) ' reciprocal vector 1: ', latt%bz1_p(1), latt%bz1_p(2)
    write(55,*) ' reciprocal vector 2: ', latt%bz2_p(1), latt%bz2_p(2)
    write(55,*) ' latt       vector 1: ', latt%a1_p(1), latt%a1_p(2)
    write(55,*) ' latt       vector 2: ', latt%a2_p(1), latt%a2_p(2)
    close(55)
  end subroutine print_latt

end module square_lattice


module cubic_lattice
   use constants, only: dp
   type cubic
      integer :: ncell, nsub, nsites, l1, l2, l3, z_plq, nn_nf, nn_lf
      real(dp) :: a1_p(3), a2_p(3), a3_p(3), l1_p(3), l2_p(3), l3_p(3), b1_p(3),  b2_p(3), b3_p(3), bz1_p(3), bz2_p(3), bz3_p(3)
      integer, allocatable :: ulist(:,:), nnlf_list(:), list(:,:), cord(:,:), invlist(:,:,:), nnlist(:,:), snlist(:,:), tnlist(:,:), plq_cord(:,:), listk(:,:), invlistk(:,:,:), imj(:,:), nnveci(:,:,:)
      real(dp), allocatable :: nnvecr(:,:,:), rcord(:,:)
      complex(dp), allocatable :: zexpiqr(:,:)
      contains 
        private
        final :: deallocate_cubic
   end type cubic     

  interface iscalar
     module procedure iscalar_ii, iscalar_ir, iscalar_rr
  end interface

  interface xnorm
     module procedure xnorm_i, xnorm_r
  end interface

  contains

  subroutine deallocate_cubic( this )
    implicit none
    type(cubic) :: this
    if( allocated(this%ulist) ) deallocate( this%ulist )
    if( allocated(this%nnlf_list) ) deallocate( this%nnlf_list )
    if( allocated(this%list) ) deallocate( this%list )
    if( allocated(this%cord) ) deallocate( this%cord )
    if( allocated(this%invlist) ) deallocate( this%invlist )
    if( allocated(this%nnlist) ) deallocate( this%nnlist )
    if( allocated(this%snlist) ) deallocate( this%snlist )
    if( allocated(this%tnlist) ) deallocate( this%tnlist )
    if( allocated(this%plq_cord) ) deallocate( this%plq_cord )
    if( allocated(this%listk) ) deallocate( this%listk )
    if( allocated(this%invlistk) ) deallocate( this%invlistk )
    if( allocated(this%imj) ) deallocate( this%imj )
    if( allocated(this%nnveci) ) deallocate( this%nnveci )
    if( allocated(this%nnvecr) ) deallocate( this%nnvecr )
    if( allocated(this%rcord) ) deallocate( this%rcord )
    if( allocated(this%zexpiqr) ) deallocate( this%zexpiqr )
    write(*,*) " deallocate_cubic is performed "
  end subroutine deallocate_cubic

  subroutine setup_cubic(l1, l2, l3, a1_p, a2_p, a3_p, latt)

    integer, intent(in) :: l1, l2, l3
    real(dp), dimension(:), intent(in) :: a1_p, a2_p, a3_p
    type(cubic) :: latt

    real(dp) :: ai_p(3), xk_p(3), b1_p(3), b2_p(3), b3_p(3), bz1_p(3), bz2_p(3), bz3_p(3), b_p(3) 
    real(dp) :: mat(3,3), mat_inv(3,3)


    latt%nsub = 1               !??
    latt%l1 = l1
    latt%l2 = l2
    latt%l3 = l3
    latt%ncell = l1*l2*l3
    latt%nsites = l1*l2*l3
    latt%z_plq = 8              ! sites contained in one cubic, 可重复计算. 
    latt%nn_nf = 6
    latt%nn_lf = l1*l2*l3/2 
    allocate( latt%rcord( latt%nsites,3) )
    allocate( latt%ulist( latt%ncell,3) )  
    allocate( latt%nnlf_list(latt%nn_lf) )
    allocate( latt%list( latt%nsites, 2) )
    allocate( latt%invlist(l1, l2, l3))
    zero = 1.e-6

    latt%l1_p = l1*a1_p
    latt%l2_p = l2*a2_p
    latt%l3_p = l3*a3_p 

    latt%a1_p = a1_p
    latt%a2_p = a2_p
    latt%a3_p = a3_p

    ind = 0
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            do i3 = 1, latt%l3
                latt%ulist( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 1) = i1
                latt%ulist( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 2) = i2
                latt%ulist( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 3) = i3
                latt%list( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 1) = (i3-1)*(l1*l2)+(i1-1)*l2+i2
                latt%list( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 2) = 1              ! orbital
                latt%rcord( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 1) = dble(i1)
                latt%rcord( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 2) = dble(i2) 
                latt%rcord( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 3) = dble(i3) 
                latt%invlist(i1,i2,i3) = (i3-1)*(l1*l2)+(i1-1)*l2+i2
                if ( mod(i1+i2+i3,2) == 0 ) then                    !??
                    ind = ind + 1
                    latt%nnlf_list(ind) = (i3-1)*(l1*l2)+(i1-1)*l2+i2
                end if 
            end do
        end do
    end do        
    latt%cord = latt%ulist


    !set nearest neighbor
    allocate( latt%nnvecr(3, 6, latt%nsites))
    allocate( latt%nnveci(3, 6, latt%nsites))
    allocate( latt%nnlist(latt%nsites, 6))
    allocate( latt%snlist(latt%nsites, 12))
    allocate( latt%tnlist(latt%nsites, 8))

    ! A->B
    do i = 1, latt%nsites, 2
        latt%nnvecr(:,1,i) = (/ 1.d0, 0.d0, 0.d0/)   ! in the unit of a1_p and a2_p
        latt%nnvecr(:,2,i) = (/ 0.d0, 1.d0, 0.d0/)   
        latt%nnvecr(:,3,i) = (/ 0.d0, 0.d0, 1.d0/)   
        latt%nnvecr(:,4,i) = (/-1.d0, 0.d0, 0.d0/)   
        latt%nnvecr(:,5,i) = (/ 0.d0,-1.d0, 0.d0/)   
        latt%nnvecr(:,6,i) = (/ 0.d0, 0.d0,-1.d0/)   

        latt%nnveci(:, 1, i) = (/ 1, 0, 0/) 
        latt%nnveci(:, 2, i) = (/ 0, 1, 0/)
        latt%nnveci(:, 3, i) = (/ 0, 0, 1/)
        latt%nnveci(:, 4, i) = (/-1, 0, 0/) 
        latt%nnveci(:, 5, i) = (/ 0,-1, 0/) 
        latt%nnveci(:, 6, i) = (/ 0, 0,-1/) 
    end do    

    ! B->A
    do i = 2, latt%nsites, 2
        latt%nnvecr(:,1,i) = (/ 1.d0, 0.d0, 0.d0/)   ! in the unit of a1_p and a2_p
        latt%nnvecr(:,2,i) = (/ 0.d0, 1.d0, 0.d0/)   
        latt%nnvecr(:,3,i) = (/ 0.d0, 0.d0, 1.d0/)   
        latt%nnvecr(:,4,i) = (/-1.d0, 0.d0, 0.d0/)   
        latt%nnvecr(:,5,i) = (/ 0.d0,-1.d0, 0.d0/)   
        latt%nnvecr(:,6,i) = (/ 0.d0, 0.d0,-1.d0/)   

        latt%nnveci(:, 1, i) = (/ 1, 0, 0/) 
        latt%nnveci(:, 2, i) = (/ 0, 1, 0/)
        latt%nnveci(:, 3, i) = (/ 0, 0, 1/)
        latt%nnveci(:, 4, i) = (/-1, 0, 0/) 
        latt%nnveci(:, 5, i) = (/ 0,-1, 0/) 
        latt%nnveci(:, 6, i) = (/ 0, 0,-1/) 
    end do 

    do i = 1, latt%ncell
        ! nearest neighbor
        latt%nnlist(i,1) = ( npbc(latt%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),   l1) - 1 )*l2 + npbc(latt%ulist(i,2),   l2) 
        latt%nnlist(i,2) = ( npbc(latt%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1, l1) - 1 )*l2 + npbc(latt%ulist(i,2),   l2) 
        latt%nnlist(i,3) = ( npbc(latt%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),   l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%nnlist(i,4) = ( npbc(latt%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),   l1) - 1 )*l2 + npbc(latt%ulist(i,2),   l2) 
        latt%nnlist(i,5) = ( npbc(latt%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1, l1) - 1 )*l2 + npbc(latt%ulist(i,2),   l2) 
        latt%nnlist(i,6) = ( npbc(latt%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),   l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 

        ! next nearest neighbor 
        latt%snlist(i,1) = ( npbc(latt%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%snlist(i,2) = ( npbc(latt%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%snlist(i,3) = ( npbc(latt%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 
        latt%snlist(i,4) = ( npbc(latt%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1,l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 

        latt%snlist(i,5) = ( npbc(latt%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%snlist(i,6) = ( npbc(latt%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%snlist(i,7) = ( npbc(latt%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 
        latt%snlist(i,8) = ( npbc(latt%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1),l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 

        latt%snlist(i,9) =  ( npbc(latt%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1,l1) - 1 )*l2 + npbc(latt%ulist(i,2), l2)  
        latt%snlist(i,10) = ( npbc(latt%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1,l1) - 1 )*l2 + npbc(latt%ulist(i,2), l2) 
        latt%snlist(i,11) = ( npbc(latt%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1,l1) - 1 )*l2 + npbc(latt%ulist(i,2), l2) 
        latt%snlist(i,12) = ( npbc(latt%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1,l1) - 1 )*l2 + npbc(latt%ulist(i,2), l2) 

        ! next next nearest neighbor
        latt%tnlist(i,1) = ( npbc(latt%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%tnlist(i,2) = ( npbc(latt%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%tnlist(i,3) = ( npbc(latt%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 
        latt%tnlist(i,4) = ( npbc(latt%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 
        latt%tnlist(i,5) = ( npbc(latt%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%tnlist(i,6) = ( npbc(latt%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)+1, l2) 
        latt%tnlist(i,7) = ( npbc(latt%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)-1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 
        latt%tnlist(i,8) = ( npbc(latt%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(latt%ulist(i,1)+1, l1) - 1 )*l2 + npbc(latt%ulist(i,2)-1, l2) 
    end do


    ! set plq_cord          !??
    allocate( latt%plq_cord(8, latt%ncell))
    id1 = 0; id2 = 0; id3 = 0; id4 = 0; id5 = 0; id6 = 0; id7 = 0; id8 = 0; 
    do i1 = 1, latt%l1
        do i2 = 1, latt%l2
            do i3 = 1, latt%l3
                if( mod(i1,2) == 1 .and. mod(i2,2) == 1 .and. mod(i3,2) == 1 ) then
                    id1 = id1 + 1
                else if ( mod(i1,2) == 1 .and. mod(i2,2) == 1 .and. mod(i3,2) == 0 ) then
                    id2 = id2 + 1
                else if ( mod(i1,2) == 1 .and. mod(i2,2) == 0 .and. mod(i3,2) == 1 ) then
                    id3 = id3 + 1
                else if ( mod(i1,2) == 1 .and. mod(i2,2) == 0 .and. mod(i3,2) == 0 ) then
                    id4 = id4 + 1
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 1 .and. mod(i3,2) == 1) then    
                    id5 = id5 + 1
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 1 .and. mod(i3,2) == 0 ) then
                    id6 = id6 + 1
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 0 .and. mod(i3,2) == 1 ) then
                    id7 = id7 + 1
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 0 .and. mod(i3,2) == 0 ) then 
                    id8 = id8 + 1
                end if
            end do
        end do
    end do



    nd1 = id1; nd2 = id2; nd3 = id3; nd4 = id4; nd5 = id5; nd6 = id6; nd7 = id7; nd8 = id8;
    id1 = 0; id2 = 0; id3 = 0; id4 = 0; id5 = 0; id6 = 0; id7 = 0; id8 = 0;
    do i1 = 1, latt%l1 
        do i2 = 1, latt%l2
            do i3 = 1, latt%l3

                if( mod(i1,2) == 1 .and. mod(i2,2) == 1 .and. mod(i3,2) == 1 ) then
                    id1 = id1 + 1
                    id = id1
                else if ( mod(i1,2) == 1 .and. mod(i2,2) == 1 .and. mod(i3,2) == 0 ) then
                    id2 = id2 + 1
                    id = id2 + nd1
                else if ( mod(i1,2) == 1 .and. mod(i2,2) == 0 .and. mod(i3,2) == 1 ) then
                    id3 = id3 + 1
                    id = id3 + nd1 + nd2
                else if ( mod(i1,2) == 1 .and. mod(i2,2) == 0 .and. mod(i3,2) == 0 ) then
                    id4 = id4 + 1
                    id = id4 + nd1 + nd2 + nd3
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 1 .and. mod(i3,2) == 1) then    
                    id5 = id5 + 1
                    id = id5 + nd1 + nd2 + nd3 + nd4
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 1 .and. mod(i3,2) == 0 ) then
                    id6 = id6 + 1
                    id = id6 + nd1 + nd2 + nd3 + nd4 + nd5
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 0 .and. mod(i3,2) == 1 ) then
                    id7 = id7 + 1
                    id = id7 + nd1 + nd2 + nd3 + nd4 + nd5 + nd6
                else if ( mod(i1,2) == 0 .and. mod(i2,2) == 0 .and. mod(i3,2) == 0 ) then 
                    id8 = id8 + 1
                    id = id8 + nd1 + nd2 + nd3 + nd4 + nd5 + nd6 + nd7
                end if

                latt%plq_cord(1, id) = latt%invlist(i1, i2, i3)         
                latt%plq_cord(2, id) = latt%invlist(npbc(i1+1,l1), i2, i3)
                latt%plq_cord(3, id) = latt%invlist(npbc(i1+1,l1), npbc(i2+1,l2), i3)
                latt%plq_cord(4, id) = latt%invlist(i1, npbc(i2+1,l2), i3)
                latt%plq_cord(5, id) = latt%invlist(i1, i2, npbc(i3+1, l3))         
                latt%plq_cord(6, id) = latt%invlist(npbc(i1+1,l1), i2, npbc(i3+1, l3))
                latt%plq_cord(7, id) = latt%invlist(npbc(i1+1,l1), npbc(i2+1,l2), npbc(i3+1,l3))
                latt%plq_cord(8, id) = latt%invlist(i1, npbc(i2+1,l2), npbc(i3+1,l3))

            end do 
        end do
    end do 


    !setting up k-space lattice
    allocate( latt%listk(latt%ncell,3), latt%invlistk(-l1:l1, -l2:l2, -l3:l3))
    pi = dacos(-1.d0)
    ! setup the 3x3 matrix to determine bz1_p, bz2_p, bz3_p
    mat(1,1) = dble(a1_p(1))
    mat(1,2) = dble(a1_p(2))
    mat(1,3) = dble(a1_p(3))
    mat(2,1) = dble(a2_p(1))
    mat(2,2) = dble(a2_p(2))
    mat(2,3) = dble(a2_p(3))
    mat(3,1) = dble(a3_p(1))
    mat(3,2) = dble(a3_p(2))
    mat(3,3) = dble(a3_p(3))

    x = mat(1,1)*( mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2) ) - mat(1,2)*( mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1) ) + mat(1,3)*( mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1) ) 
    mat_inv(1,1) = 1.d0/x * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2))
    mat_inv(1,2) = 1.d0/x * (mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)) 
    mat_inv(1,3) = 1.d0/x * (mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)) 
    mat_inv(2,1) = 1.d0/x * (mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)) 
    mat_inv(2,2) = 1.d0/x * (mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)) 
    mat_inv(2,3) = 1.d0/x * (mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)) 
    mat_inv(3,1) = 1.d0/x * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)) 
    mat_inv(3,2) = 1.d0/x * (mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)) 
    mat_inv(3,3) = 1.d0/x * (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)) 

    bz1_p(1)     = 2.d0*pi*mat_inv(1,1) 
    bz1_p(2)     = 2.d0*pi*mat_inv(2,1) 
    bz1_p(3)     = 2.d0*pi*mat_inv(3,1)

    bz2_p(1)     = 2.d0*pi*mat_inv(1,2) 
    bz2_p(2)     = 2.d0*pi*mat_inv(2,2) 
    bz2_p(3)     = 2.d0*pi*mat_inv(3,2)

    bz3_p(1)     = 2.d0*pi*mat_inv(1,3) 
    bz3_p(2)     = 2.d0*pi*mat_inv(2,3) 
    bz3_p(3)     = 2.d0*pi*mat_inv(3,3)

    latt%bz1_p   = bz1_p 
    latt%bz2_p   = bz2_p 
    latt%bz3_p   = bz3_p 


    b1_p = bz1_p/(l1*1.d0)
    b2_p = bz2_p/(l1*1.d0)
    b3_p = bz3_p/(l1*1.d0)

    latt%b1_p = b1_p
    latt%b2_p = b2_p
    latt%b3_p = b3_p

    nc = 0
    do l = -(l3+1)/2+1, l3/2
        do m = -(l1+1)/2+1, l1/2
            do n = -(l2+1)/2+1, l2/2
                nc = nc + 1
                latt%listk(nc,1) = m
                latt%listk(nc,2) = n
                latt%listk(nc,3) = l
                latt%invlistk(m,n,l) = nc
            end do
        end do
    end do

    if (nc.ne.latt%ncell) then 
       write(6,*) 'error ', nc, latt%ncell
       stop
    endif


    !setup imj 
    allocate ( latt%imj(latt%ncell,latt%ncell) ) 
    ! latt_imj
    do j = 1, latt%ncell
        do i = 1, latt%ncell
            imj_nx = npbc( latt%ulist(i,1) - latt%ulist(j,1), latt%l1 ) 
            imj_ny = npbc( latt%ulist(i,2) - latt%ulist(j,2), latt%l2 )
            imj_nz = npbc( latt%ulist(i,3) - latt%ulist(j,3), latt%l3 )
            latt%imj(i, j) = (imj_nz - 1)*l1*l2 + imj_ny + (imj_nx-1)*latt%l2  
        end do
    end do



    allocate( latt%zexpiqr(latt%ncell,latt%ncell) )
    do i = 1,latt%ncell
       ai_p = dble(latt%ulist(i,1))*latt%a1_p + dble(latt%ulist(i,2))*latt%a2_p + dble(latt%ulist(i,3))*latt%a3_p
       do nk = 1,latt%ncell
          xk_p = dble(latt%listk(nk,1))*latt%b1_p + dble(latt%listk(nk,2))*latt%b2_p + dble(latt%listk(nk,3))*latt%b3_p
          latt%zexpiqr(nk,i) = exp( dcmplx( 0.d0, iscalar(xk_p, ai_p)) )
       enddo
    enddo

    
  end subroutine setup_cubic


  integer function npbc(nr, l)
    implicit none
    integer, intent(in) :: nr
    integer, intent(in) :: l
    npbc = nr
    if (nr.gt.l) npbc = nr - l
    if (nr.lt.1) npbc = nr + l
  end function npbc



  integer function inv_k(xk_p,latt) 
    
    implicit none
    real (kind=8)  :: xk_p(3)
    type (cubic) :: latt
    
    real (kind=8)     :: xk1_p(3), x, zero
    integer :: nk 

    zero = 1.d-10
    nk = 1
    do 
       xk1_p = latt%listk(nk,1)*latt%b1_p + latt%listk(nk,2)*latt%b2_p + latt%listk(nk,3)*latt%b3_p 
       if (xnorm(xk1_p - xk_p)  < zero ) then
          inv_k = nk
          exit
       elseif (nk < latt%ncell) then 
          nk = nk + 1
       else
          write(6,*) 'error in inv_k lattice_new'
          stop
       endif
    enddo

  end function inv_k


  integer function iscalar_ii(i_p, j_p) 
    implicit none
    integer, dimension(:) :: i_p, j_p
    integer i
    
    iscalar_ii = 0
    do i = 1,  size(i_p)
       iscalar_ii = iscalar_ii + i_p(i)*j_p(i)
    enddo
  end function iscalar_ii



  real (kind=8)  function iscalar_ir(x_p, j_p)
    implicit none
    real (kind=8), dimension(:) ::  x_p
    integer, dimension(:) ::  j_p
    integer i
    
    iscalar_ir = 0.d0
    !write(6,*) size(i_p)
    do i = 1,  size(x_p)
      ! write(6,*) i
       iscalar_ir = iscalar_ir + x_p(i)*dble(j_p(i))
    enddo
  end function iscalar_ir


  real (kind=8)  function iscalar_rr(x_p, y_p)
    implicit none
    real (kind=8), dimension(:) ::  x_p, y_p
    integer i
    
    iscalar_rr = 0.d0
    do i = 1,  size(x_p)
       iscalar_rr = iscalar_rr + x_p(i)*y_p(i)
    enddo
  end function iscalar_rr


  real (kind=8) function xnorm_i(i_p)
    implicit none
    integer, dimension(:) :: i_p
    integer :: i

    xnorm_i = 0.d0
    do i = 1,  size(i_p)
       xnorm_i = xnorm_i + dble(i_p(i)*i_p(i))
    enddo
    xnorm_i = sqrt(xnorm_i)
  end function xnorm_i


  real (kind=8) function xnorm_r(x_p)
    implicit none
    real (kind=8), dimension(:) :: x_p
    integer :: i

    xnorm_r = 0.d0
    do i = 1,  size(x_p)
       xnorm_r = xnorm_r + x_p(i)*x_p(i)
    enddo
    xnorm_r = sqrt(xnorm_r)
  end function xnorm_r

  subroutine print_latt(latt)
    
    implicit real (kind=8) (a-g,o-z)
    implicit integer (h-n)
    
    type (cubic) :: latt
    integer          :: i_p(3),nd_p(3)
    real    (kind=8) :: x_p(3)

    open (unit=55,file="latt_info", status = "unknown")
    write(55,*) ' reciprocal vector 1: ', latt%bz1_p(1), latt%bz1_p(2), latt%bz1_p(3)
    write(55,*) ' reciprocal vector 2: ', latt%bz2_p(1), latt%bz2_p(2), latt%bz2_p(3)
    write(55,*) ' reciprocal vector 3: ', latt%bz3_p(1), latt%bz3_p(2), latt%bz3_p(3)
    write(55,*) ' latt       vector 1: ', latt%a1_p(1), latt%a1_p(2), latt%a1_p(3)
    write(55,*) ' latt       vector 2: ', latt%a2_p(1), latt%a2_p(2), latt%a2_p(3)
    write(55,*) ' latt       vector 3: ', latt%a3_p(1), latt%a3_p(2), latt%a3_p(3)
    close(55)
  end subroutine print_latt

end module cubic_lattice
