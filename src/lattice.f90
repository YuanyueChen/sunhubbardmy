module latt2d_class
  use constants, only: dp, pi
  type, public :: latt2d
    integer :: ncell ! number of unit cells
    integer :: nsub  ! number of sublattice, 1 for square, 2 for honeycomb
    integer :: nsites ! number of sites
    integer :: l1     ! length in direction a1_p
    integer :: l2     ! length in direction a2_p
    integer :: z_plq  ! number of sites for each plaquette, 4 for square, 6 for honeycomb
    integer :: nn_nf  ! decomposition of nearest neighbor bonds, number of families
    integer :: nn_lf  ! decomposition of nearest neighbor bonds, number of bonds in each families
    real(dp) :: a1_p(2), a2_p(2) ! primitive vectors
    real(dp) :: l1_p(2), l2_p(2) ! lattice vectors, length in each direction * corresponding primitive vectors
    real(dp) :: bz1_p(2), bz2_p(2) ! primitive reciprocal vectors
    real(dp) :: b1_p(2), b2_p(2)   ! lattice reciprocal vectors
    integer, allocatable :: ulist(:,:)   ! the unit cell list, dimension(ncell,2), (:,1) is the index in a1_p direction, (:,2) is the index in a2_p direction
    integer, allocatable :: nnlf_list(:) ! decomposition of nearest neighbor bonds, the starting site list, dimension(nnlf)
    integer, allocatable :: list(:,:)    ! site list, dimension(nsite,2), (:,1) is the unit cell index, (:,2) is the sublattice index
    integer, allocatable :: cord(:,:)    ! site unit cell coordination, dimension(nsites,2), (:,1) is the index in a1_p direction (:,2) is the index in a2_p direction
    integer, allocatable :: invlist(:,:,:) ! site invlist, dimension(l1,l2,nsub), given the index in a1_p direction, index in a2_p direction, and index for sublattice, identify the site index
    integer, allocatable :: nnlist(:,:)    ! nearest neighbor list, dimension(nsites,*), * is the number of nearest neighbors for each site
    integer, allocatable :: snlist(:,:)    ! second neighbor list, dimension(nistes,*), * is the number of second neighbors for each site
    integer, allocatable :: tnlist(:,:)    ! third neighbor list, dimension(nsites,*), * is the number of third neighbors for each site
    integer, allocatable :: plq_cord(:,:)  ! plaquette coordination, dimension(*1,*2), *1 is the number of sites in each plaquette, *2 is the number of plaquettes
    integer, allocatable :: listk(:,:)     ! k points list, dimension(ncell,2)
    integer, allocatable :: invlistk(:,:)  ! k points invlist, dimension(-l1:l1,-l2:l2), given the k points coordinate, identify the k point index
    integer, allocatable :: imj(:,:)       ! equivalent unit cell of Ri-Rj, dimension(ncell,ncell), given the index for unit cell Ri and Rj, identify the equivalent unit cell index for Ri-Rj
    integer, allocatable :: nnveci(:,:,:)  ! unit cell difference for nearest neighbor sites, dimension(2,*,nsites), * is the number of nearest neighbors for each site
    real(dp), allocatable :: nnvecr(:,:,:) ! real space difference for nearest neighbor sites, dimension(2,*,nsites), * is the number of nearest neighbors for each site, in the unit of a1_p and a2_p
    real(dp), allocatable :: rcord(:,:)    ! real sapce coordinates for each site, dimension(2,nsites), in the unit of a1_p and a2_p
    complex(dp), allocatable :: zexpiqr(:,:) ! exp(i*qi*Ri), dimension(ncell,ncell)

    contains
      procedure, public :: setup_latt => setup_latt2d_sub
      procedure, public :: print_latt => print_latt2d_sub
      procedure, public :: deallocate_latt => deallocate_latt2d_sub
  end type latt2d

  private :: setup_latt2d_sub
  private :: print_latt2d_sub
  private :: deallocate_latt2d

  contains

  subroutine deallocate_latt2d_sub( this )
    implicit none
    class(latt2d) :: this
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
    !write(*,*) " deallocate_latt2d is performed "
  end subroutine deallocate_latt2d_sub
  
  subroutine setup_latt2d_sub(this, l1, l2, a1_p, a2_p)
    ! set up latt2d
    
    implicit none
    
    class(latt2d) :: this
    integer, intent(in) :: l1, l2
    real(dp), dimension(:), intent(in) :: a1_p, a2_p

    ! local
    integer :: i1, i2, m, n, nc, i, j, nk, imj_nx, imj_ny
    real(dp) :: x
    real(dp) :: ai_p(2), xk_p(2), b1_p(2), b2_p(2), bz1_p(2), bz2_p(2), b_p(2)
    real(dp) :: mat(2,2), mat_inv(2,2)

    this%l1 = l1
    this%l2 = l2
    this%ncell = l1*l2

    this%l1_p = l1*a1_p
    this%l2_p = l2*a2_p
    this%a1_p = a1_p
    this%a2_p = a2_p

    ! setup ulist
    allocate( this%ulist( this%ncell,2) )
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            this%ulist( i2+(i1-1)*l2, 1) = i1
            this%ulist( i2+(i1-1)*l2, 2) = i2
        end do
    end do

    !setting up k-space lattice
    allocate ( this%listk(this%ncell,2), this%invlistk(-l1:l1, -l2:l2) )
    ! get the primitive reciprocal vectors bz1_p, bz2_p
    ! the formula [b1 b2]^T = 2*pi*[a1 a2]^{-1} is used. This formula can be generalized to any dimension.
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
    this%bz1_p = bz1_p
    this%bz2_p = bz2_p

    ! get the lattice reciprocal vectors b1_p, b2_p
    ! the above same formula is used.
    mat(1,1) = dble(this%l1_p(1))
    mat(1,2) = dble(this%l1_p(2))
    mat(2,1) = dble(this%l2_p(1))
    mat(2,2) = dble(this%l2_p(2))
    x = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)
    mat_inv(1,1) =  mat(2,2)/x
    mat_inv(2,2) =  mat(1,1)/x
    mat_inv(1,2) = -mat(1,2)/x
    mat_inv(2,1) = -mat(2,1)/x
    b1_p(1)      = 2.d0*pi*mat_inv(1,1)
    b1_p(2)      = 2.d0*pi*mat_inv(2,1)
    b2_p(1)      = 2.d0*pi*mat_inv(1,2)
    b2_p(2)      = 2.d0*pi*mat_inv(2,2)
    this%b1_p = b1_p
    this%b2_p = b2_p
    
    nc = 0
    do m = -(l1+1)/2+1, l1/2 
      do n = -(l2+1)/2+1, l2/2
        nc = nc + 1
        this%listk(nc,1) = m
        this%listk(nc,2) = n
        this%invlistk(m,n) = nc
      end do
    end do
    if (nc.ne.this%ncell) then 
      write(6,*) 'error ', nc, this%ncell
      stop
    endif

    ! setup imj 
    allocate ( this%imj(this%ncell,this%ncell) )
    do j = 1, this%ncell
      do i = 1, this%ncell
        imj_nx = npbc( this%ulist(i,1) - this%ulist(j,1), this%l1 ) 
        imj_ny = npbc( this%ulist(i,2) - this%ulist(j,2), this%l2 )
        this%imj(i,j) = imj_ny + (imj_nx-1)*this%l2
      end do
    end do

    ! setup zexpiqr
    allocate( this%zexpiqr(this%ncell,this%ncell) )
    do i = 1,this%ncell
      ai_p = dble(this%ulist(i,1))*this%a1_p + dble(this%ulist(i,2))*this%a2_p
      do nk = 1,this%ncell
        xk_p = dble(this%listk(nk,1))*this%b1_p + dble(this%listk(nk,2))*this%b2_p
        this%zexpiqr(nk,i) = exp( dcmplx( 0.d0, xk_p(1)*ai_p(1)+xk_p(2)*ai_p(2) ) )
      enddo
    enddo
    
  end subroutine setup_latt2d_sub

  integer function npbc(nr,l)
    implicit none
    integer, intent(in) :: nr
    integer, intent(in) :: l
    npbc = nr
    if (nr.gt.l) npbc = nr - l
    if (nr.lt.1) npbc = nr + l
  end function npbc

  subroutine print_latt2d_sub(this, fid)
    
    implicit none

    class(latt2d) :: this
    integer,intent(in) :: fid

    integer :: i, j, imj

    write(fid,'(a)') 'lattice information: '
    write(fid,'(a,2f16.12)') ' primitive lattice  vector a1_p(:) =  ', this%a1_p(:)
    write(fid,'(a,2f16.12)') ' primitive lattice  vector a2_p(:) =  ', this%a2_p(:)
    write(fid,'(a,2f16.12)') ' primitive reciprocal vector bz1_p / 2pi = ', this%bz1_p(:)/(2.d0*pi)
    write(fid,'(a,2f16.12)') ' primitive reciprocal vector bz2_p / 2pi = ', this%bz2_p(:)/(2.d0*pi)
    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice sites list '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     list(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,3i4)') i,  this%list(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice sites cord '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     cord(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,3i4)') i,  this%cord(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice nnlist '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     nnlist(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,12i4)') i,  this%nnlist(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice snlist '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     snlist(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,12i4)') i,  this%snlist(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice tnlist '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     tnlist(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,12i4)') i,  this%tnlist(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice plq_cord '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     plq_cord(:,i) '
    do i = 1, this%ncell
        write(fid,'(i6,8i4)') i,  this%plq_cord(:,i)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The momentum listk '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     listk(i,:) '
    do i = 1, this%ncell
        write(fid,'(i6,3i4)') i,  this%listk(i,:)
    end do
  
#IFDEF PLEVEL2
    if( this%ncell < 50 ) then
    write(fid, *)
    write(fid,'(a)') '----------------------------'
    write(fid,'(a)') ' latt_imj info   '
    write(fid,'(a)') '----------------------------'
    write(fid, '(a)') '   i   j   imj '
    do j = 1, this%ncell
        do i = 1, this%ncell
            imj  = this%imj(i,j)
            write(fid, '(3i4)') i, j, imj
        end do
    end do
    end if
#ENDIF
  end subroutine print_latt2d_sub

end module latt2d_class


module honeycomb_lattice
  use constants, only: dp
  use latt2d_class
  type, public, extends(latt2d) :: honeycomb
     contains
       procedure, public :: setup_honeycomb => setup_honeycomb_sub
  end type honeycomb

  private :: setup_honeycomb_sub

  contains

  subroutine setup_honeycomb_sub(this, l1, l2, a1_p, a2_p) 

    ! set up honeycomb lattice
    
    implicit none
    
    class(honeycomb) :: this
    integer, intent(in) :: l1, l2
    real(dp), dimension(:), intent(in) :: a1_p, a2_p

    ! local
    integer :: i1, i2, m, n, nc, i, j, nk, ind, iorb, id, id1, id2, id3, nd1, nd2, nd3
    real(dp) :: x
    real(dp) :: ai_p(2), xk_p(2), b1_p(2), b2_p(2), bz1_p(2), bz2_p(2), b_p(2)
    real(dp) :: mat(2,2), mat_inv(2,2)

    call this%setup_latt(l1, l2, a1_p, a2_p) 

    this%nsub = 2
    this%nsites = l1*l2*2
    this%z_plq = 6
    this%nn_nf = 3
    this%nn_lf = l1*l2

    allocate( this%cord( this%nsites,2) )
    allocate( this%rcord(this%nsites,2) )
    allocate( this%nnlf_list( this%nn_lf) )
    allocate( this%list( this%nsites,2) )
    allocate( this%invlist(l1, l2, this%nsub) )

    do i = 1, this%nn_lf
        this%nnlf_list(i) = 2*i - 1
    end do

    ind = 0
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            do iorb = 1, this%nsub
                ind = ind + 1
                this%list(ind, 1) = i2+(i1-1)*l2 ! unit cell
                this%list(ind, 2) = iorb         ! orb
                this%cord(ind, 1) = i1
                this%cord(ind, 2) = i2
                this%invlist(i1,i2,iorb) = ind
                if( iorb == 1 ) then ! A site
                    this%rcord(ind, 1) = dble(i1)
                    this%rcord(ind, 2) = dble(i2)
                else ! B site
                    this%rcord(ind, 1) = dble(i1) - 1.d0/3.d0
                    this%rcord(ind, 2) = dble(i2) + 1.d0/3.d0
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
    allocate( this%nnvecr(2,3,this%nsites) )
    allocate( this%nnveci(2,3,this%nsites) )
    allocate( this%nnlist(this%nsites,3) )
    allocate( this%snlist(this%nsites,6) )
    allocate( this%tnlist(this%nsites,3) )
    ! A->B
    do i = 1, this%nsites, 2
        this%nnvecr(:,1,i) = (/-1.d0/3.d0, 1.d0/3.d0/) ! in the unit of a1_p and a2_p
        this%nnvecr(:,2,i) = (/ 2.d0/3.d0, 1.d0/3.d0/)
        this%nnvecr(:,3,i) = (/-1.d0/3.d0,-2.d0/3.d0/)
        this%nnveci(:,1,i) = (/ 0, 0 /)
        this%nnveci(:,2,i) = (/ 1, 0 /)
        this%nnveci(:,3,i) = (/ 0,-1 /)
    end do
    ! B->A
    do i = 2, this%nsites, 2
        this%nnvecr(:,1,i) = (/ 1.d0/3.d0,-1.d0/3.d0/) ! in the unit of a1_p and a2_p
        this%nnvecr(:,2,i) = (/-2.d0/3.d0,-1.d0/3.d0/)
        this%nnvecr(:,3,i) = (/ 1.d0/3.d0, 2.d0/3.d0/)
        this%nnveci(:,1,i) = (/ 0, 0 /)
        this%nnveci(:,2,i) = (/-1, 0 /)
        this%nnveci(:,3,i) = (/ 0, 1 /)
    end do
    do i  = 1, this%ncell*this%nsub
        if( mod(i,2) == 1 ) then ! A site
            this%nnlist(i,1) = i+1
            this%nnlist(i,2) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)  , l2 ) )*this%nsub
            this%nnlist(i,3) = ( (npbc(this%cord(i,1)  , l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub
            this%snlist(i,1) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub - 1
            this%snlist(i,2) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)  , l2 ) )*this%nsub - 1
            this%snlist(i,3) = ( (npbc(this%cord(i,1)  , l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub - 1
            this%snlist(i,4) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub - 1
            this%snlist(i,5) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)  , l2 ) )*this%nsub - 1
            this%snlist(i,6) = ( (npbc(this%cord(i,1)  , l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub - 1
            this%tnlist(i,1) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub
            this%tnlist(i,2) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub
            this%tnlist(i,3) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub
        else                     ! B site
            this%nnlist(i,1) = i-1
            this%nnlist(i,2) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)  , l2 ) )*this%nsub - 1
            this%nnlist(i,3) = ( (npbc(this%cord(i,1)  , l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub - 1
            this%snlist(i,1) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub
            this%snlist(i,2) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)  , l2 ) )*this%nsub
            this%snlist(i,3) = ( (npbc(this%cord(i,1)  , l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub
            this%snlist(i,4) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub
            this%snlist(i,5) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)  , l2 ) )*this%nsub
            this%snlist(i,6) = ( (npbc(this%cord(i,1)  , l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub
            this%tnlist(i,1) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub - 1
            this%tnlist(i,2) = ( (npbc(this%cord(i,1)+1, l1) - 1)*l2 + npbc(this%cord(i,2)+1, l2 ) )*this%nsub - 1
            this%tnlist(i,3) = ( (npbc(this%cord(i,1)-1, l1) - 1)*l2 + npbc(this%cord(i,2)-1, l2 ) )*this%nsub - 1
        end if
    end do

! set plq_cord
    allocate( this%plq_cord(6,this%ncell))
    ! count number of A,B,C type plaq first
    id1 = 0; id2 = 0; id3 = 0
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            !i = i2 + (i1-1)*this%l2
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
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            !i = i2 + (i1-1)*this%l2
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
            this%plq_cord(1, id) = this%invlist(i1, i2, 2)
            this%plq_cord(2, id) = this%invlist(i1, i2, 1)
            this%plq_cord(3, id) = this%invlist(i1, npbc(i2-1,l2), 2)
            this%plq_cord(4, id) = this%invlist(npbc(i1-1,l1), npbc(i2-1,l2), 1)
            this%plq_cord(5, id) = this%invlist(npbc(i1-1,l1), npbc(i2-1,l2), 2)
            this%plq_cord(6, id) = this%invlist(npbc(i1-1,l1), i2, 1)
        end do
    end do
    
  end subroutine setup_honeycomb_sub

end module honeycomb_lattice


module square_lattice
  use constants, only: dp
  use latt2d_class
  type, public, extends(latt2d) :: square
     contains
       procedure, public :: setup_square => setup_square_sub
  end type square
  private :: setup_square_sub

  contains

  subroutine setup_square_sub(this, l1, l2, a1_p, a2_p)

    ! set up square lattice
    
    implicit real (kind=8) (a-g,o-z)
    implicit integer (h-n)
    
    integer, intent(in) :: l1, l2
    real(dp), dimension(:), intent(in) :: a1_p, a2_p
    class(square) :: this

    ! local
    real(dp) :: ai_p(2), xk_p(2), b1_p(2), b2_p(2), bz1_p(2), bz2_p(2), b_p(2)
    real(dp) :: mat(2,2), mat_inv(2,2)

    call this%setup_latt(l1, l2, a1_p, a2_p) 

    this%nsub = 1
    this%nsites = l1*l2
    this%z_plq = 4
    this%nn_nf = 4
    this%nn_lf = l1*l2/2

    allocate( this%cord( this%nsites,2) )
    allocate( this%rcord(this%nsites,2) )
    allocate( this%nnlf_list( this%nn_lf) )
    allocate( this%list( this%nsites,2) )
    allocate( this%invlist(l1, l2, this%nsub) )

    ind = 0
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            this%list( i2+(i1-1)*l2, 1) = i2+(i1-1)*l2
            this%list( i2+(i1-1)*l2, 2) = 1
            this%rcord( i2+(i1-1)*l2, 1) = dble(i1)
            this%rcord( i2+(i1-1)*l2, 2) = dble(i2)
            this%invlist(i1,i2,1) = i2+(i1-1)*l2
            if( mod(i1+i2,2) == 0 ) then
                ind = ind + 1
                this%nnlf_list(ind) = i2+(i1-1)*l2
            end if
        end do
    end do
    this%cord = this%ulist
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
    allocate( this%nnvecr(2,4,this%nsites) )
    allocate( this%nnveci(2,4,this%nsites) )
    allocate( this%nnlist(this%nsites,4) )
    allocate( this%snlist(this%nsites,4) )
    allocate( this%tnlist(this%nsites,4) )
    ! A->B
    do i = 1, this%nsites, 2
        this%nnvecr(:,1,i) = (/ 1.d0, 0.d0/) ! in the unit of a1_p and a2_p
        this%nnvecr(:,2,i) = (/ 0.d0, 1.d0/)
        this%nnvecr(:,3,i) = (/-1.d0, 0.d0/)
        this%nnvecr(:,4,i) = (/ 0.d0,-1.d0/)

        this%nnveci(:,1,i) = (/ 1, 0 /)
        this%nnveci(:,2,i) = (/ 0, 1 /)
        this%nnveci(:,3,i) = (/-1, 0 /)
        this%nnveci(:,4,i) = (/ 0,-1 /)
    end do
    ! B->A
    do i = 2, this%nsites, 2
        this%nnvecr(:,1,i) = (/ 1.d0, 0.d0/) ! in the unit of a1_p and a2_p
        this%nnvecr(:,2,i) = (/ 0.d0, 1.d0/)
        this%nnvecr(:,3,i) = (/-1.d0, 0.d0/)
        this%nnvecr(:,4,i) = (/ 0.d0,-1.d0/)
        this%nnveci(:,1,i) = (/ 1, 0 /)
        this%nnveci(:,2,i) = (/ 0, 1 /)
        this%nnveci(:,3,i) = (/-1, 0 /)
        this%nnveci(:,4,i) = (/ 0,-1 /)
    end do
    do i  = 1, this%ncell
        this%nnlist(i,1) = ( npbc(this%ulist(i,1)+1,l1) - 1 )*l2 + npbc(this%ulist(i,2)  , l2)
        this%nnlist(i,2) = ( npbc(this%ulist(i,1)  ,l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2)
        this%nnlist(i,3) = ( npbc(this%ulist(i,1)-1,l1) - 1 )*l2 + npbc(this%ulist(i,2)  , l2)
        this%nnlist(i,4) = ( npbc(this%ulist(i,1)  ,l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2)

        this%snlist(i,1) = ( npbc(this%ulist(i,1)+1,l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2)
        this%snlist(i,2) = ( npbc(this%ulist(i,1)-1,l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2)
        this%snlist(i,3) = ( npbc(this%ulist(i,1)-1,l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2)
        this%snlist(i,4) = ( npbc(this%ulist(i,1)+1,l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2)

        this%tnlist(i,1) = ( npbc(this%ulist(i,1)+2,l1) - 1 )*l2 + npbc(this%ulist(i,2)  , l2)
        this%tnlist(i,2) = ( npbc(this%ulist(i,1)  ,l1) - 1 )*l2 + npbc(this%ulist(i,2)+2, l2)
        this%tnlist(i,3) = ( npbc(this%ulist(i,1)-2,l1) - 1 )*l2 + npbc(this%ulist(i,2)  , l2)
        this%tnlist(i,4) = ( npbc(this%ulist(i,1)  ,l1) - 1 )*l2 + npbc(this%ulist(i,2)-2, l2)
    end do

    ! set plq_cord
    allocate( this%plq_cord(4,this%ncell))
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
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            !i = i2 + (i1-1)*this%l2
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
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            !i = i2 + (i1-1)*this%l2
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
            this%plq_cord(1, id) = this%invlist(i1, i2, 1)
            this%plq_cord(2, id) = this%invlist(npbc(i1+1,l1), i2, 1)
            this%plq_cord(3, id) = this%invlist(npbc(i1+1,l1), npbc(i2+1,l2), 1)
            this%plq_cord(4, id) = this%invlist(i1, npbc(i2+1,l2), 1)
        end do
    end do

  end subroutine setup_square_sub

end module square_lattice

module latt3d_class
  use constants, only: dp, pi
  type, public :: latt3d
    integer :: ncell ! number of unit cells
    integer :: nsub  ! number of sublattice, 1 for cubic
    integer :: nsites ! number of sites
    integer :: l1     ! length in direction a1_p
    integer :: l2     ! length in direction a2_p
    integer :: l3     ! length in direction a3_p
    integer :: z_plq  ! number of sites for each plaquette, 4 for square, 6 for honeycomb
    integer :: nn_nf  ! decomposition of nearest neighbor bonds, number of families
    integer :: nn_lf  ! decomposition of nearest neighbor bonds, number of bonds in each families
    real(dp) :: a1_p(3), a2_p(3), a3_p(3) ! primitive vectors
    real(dp) :: l1_p(3), l2_p(3), l3_p(3) ! lattice vectors, length in each direction * corresponding primitive vectors
    real(dp) :: bz1_p(3), bz2_p(3), bz3_p(3) ! primitive reciprocal vectors
    real(dp) :: b1_p(3), b2_p(3), b3_p(3)   ! lattice reciprocal vectors
    integer, allocatable :: ulist(:,:)   ! the unit cell list, dimension(ncell,3), (:,1) is the index in a1_p direction, (:,2) is the index in a2_p direction, (:,3) is the index in a3_p direction
    integer, allocatable :: nnlf_list(:) ! decomposition of nearest neighbor bonds, the starting site list, dimension(nnlf)
    integer, allocatable :: list(:,:)    ! site list, dimension(nsites,2), (:,1) is the unit cell index, (:,2) is the sublattice index
    integer, allocatable :: cord(:,:)    ! site unit cell coordination, dimension(nsites,3), (:,1) is the index in a1_p direction (:,2) is the index in a2_p direction, (:,3) is the index in a3_p direction
    integer, allocatable :: invlist(:,:,:,:) ! site invlist, dimension(l1,l2,l3,2), given the index in a1_p direction, index in a2_p direction, index in a3_p direction and index for sublattice, identify the site index
    integer, allocatable :: nnlist(:,:)    ! nearest neighbor list, dimension(nsites,*), * is the number of nearest neighbors for each site
    integer, allocatable :: snlist(:,:)    ! second neighbor list, dimension(nistes,*), * is the number of second neighbors for each site
    integer, allocatable :: tnlist(:,:)    ! third neighbor list, dimension(nsites,*), * is the number of third neighbors for each site
    integer, allocatable :: plq_cord(:,:)  ! plaquette coordination, dimension(*1,*2), *1 is the number of sites in each plaquette, *2 is the number of plaquettes
    integer, allocatable :: listk(:,:)     ! k points list, dimension(ncell,3)
    integer, allocatable :: invlistk(:,:,:)  ! k points invlist, dimension(-l1:l1,-l2:l2,-l3:l3), given the k points coordinate, identify the k point index
    integer, allocatable :: imj(:,:)       ! equivalent unit cell of Ri-Rj, dimension(ncell,ncell), given the index for unit cell Ri and Rj, identify the equivalent unit cell index for Ri-Rj
    integer, allocatable :: nnveci(:,:,:)  ! unit cell difference for nearest neighbor sites, dimension(3,*,nsites), * is the number of nearest neighbors for each site
    real(dp), allocatable :: nnvecr(:,:,:) ! real space difference for nearest neighbor sites, dimension(3,*,nsites), * is the number of nearest neighbors for each site, in the unit of a1_p and a2_p
    real(dp), allocatable :: rcord(:,:)    ! real sapce coordinates for each site, dimension(3,nsites), in the unit of a1_p and a2_p
    complex(dp), allocatable :: zexpiqr(:,:) ! exp(i*qi*Ri), dimension(ncell,ncell)

    contains 
      procedure, public :: setup_latt => setup_latt3d_sub
      procedure, public :: print_latt => print_latt3d_sub
      procedure, public :: deallocate_latt => deallocate_latt3d_sub
  end type latt3d

  private :: setup_latt3d_sub
  private :: print_latt3d_sub
  private :: deallocate_latt3d_sub
  private :: deallocate_latt3d

  contains

  subroutine deallocate_latt3d_sub( this )
    implicit none
    class(latt3d) :: this
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
    !write(*,*) " deallocate_latt3d is performed "
  end subroutine deallocate_latt3d_sub

  subroutine setup_latt3d_sub(this, l1, l2, l3, a1_p, a2_p, a3_p)

    class(latt3d) :: this
    integer, intent(in) :: l1, l2, l3
    real(dp), dimension(:), intent(in) :: a1_p, a2_p, a3_p

    real(dp) :: ai_p(3), xk_p(3), b1_p(3), b2_p(3), b3_p(3), bz1_p(3), bz2_p(3), bz3_p(3), b_p(3) 
    real(dp) :: mat(3,3), mat_inv(3,3)


    this%l1 = l1
    this%l2 = l2
    this%l3 = l3
    this%ncell = l1*l2*l3

    this%l1_p = l1*a1_p
    this%l2_p = l2*a2_p
    this%l3_p = l3*a3_p 

    this%a1_p = a1_p
    this%a2_p = a2_p
    this%a3_p = a3_p

    allocate( this%ulist( this%ncell,3) )  
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            do i3 = 1, this%l3
                this%ulist( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 1) = i1
                this%ulist( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 2) = i2
                this%ulist( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 3) = i3
            end do
        end do
    end do        

    !setting up k-space lattice
    allocate( this%listk(this%ncell,3), this%invlistk(-l1:l1, -l2:l2, -l3:l3))
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

    this%bz1_p   = bz1_p 
    this%bz2_p   = bz2_p 
    this%bz3_p   = bz3_p 


    b1_p = bz1_p/(l1*1.d0)
    b2_p = bz2_p/(l2*1.d0)
    b3_p = bz3_p/(l3*1.d0)

    this%b1_p = b1_p
    this%b2_p = b2_p
    this%b3_p = b3_p

    nc = 0
    do l = -(l3+1)/2+1, l3/2
        do m = -(l1+1)/2+1, l1/2
            do n = -(l2+1)/2+1, l2/2
                nc = nc + 1
                this%listk(nc,1) = m
                this%listk(nc,2) = n
                this%listk(nc,3) = l
                this%invlistk(m,n,l) = nc
            end do
        end do
    end do

    if (nc.ne.this%ncell) then 
       write(6,*) 'error ', nc, this%ncell
       stop
    endif


    !setup imj 
    allocate ( this%imj(this%ncell,this%ncell) ) 
    do j = 1, this%ncell
        do i = 1, this%ncell
            imj_nx = npbc( this%ulist(i,1) - this%ulist(j,1), this%l1 ) 
            imj_ny = npbc( this%ulist(i,2) - this%ulist(j,2), this%l2 )
            imj_nz = npbc( this%ulist(i,3) - this%ulist(j,3), this%l3 )
            this%imj(i, j) = (imj_nz - 1)*l1*l2 + imj_ny + (imj_nx-1)*this%l2  
        end do
    end do

    ! setup zexpiqr
    allocate( this%zexpiqr(this%ncell,this%ncell) )
    do i = 1,this%ncell
       ai_p = dble(this%ulist(i,1))*this%a1_p + dble(this%ulist(i,2))*this%a2_p + dble(this%ulist(i,3))*this%a3_p
       do nk = 1,this%ncell
          xk_p = dble(this%listk(nk,1))*this%b1_p + dble(this%listk(nk,2))*this%b2_p + dble(this%listk(nk,3))*this%b3_p
          this%zexpiqr(nk,i) = exp( dcmplx( 0.d0, xk_p(1)*ai_p(1)+xk_p(2)*ai_p(2)+xk_p(3)*ai_p(3) ) )
       enddo
    enddo
    
  end subroutine setup_latt3d_sub

  integer function npbc(nr, l)
    implicit none
    integer, intent(in) :: nr
    integer, intent(in) :: l
    npbc = nr
    if (nr.gt.l) npbc = nr - l
    if (nr.lt.1) npbc = nr + l
  end function npbc


  subroutine print_latt3d_sub(this,fid)
    
    implicit none
    
    class(latt3d) :: this
    integer, intent(in) :: fid

    integer :: i, j, imj

    write(fid,'(a)') 'lattice information: '
    write(fid,'(a,3f16.12)') ' primitive lattice vector a1_p(:) = ', this%a1_p(:)
    write(fid,'(a,3f16.12)') ' primitive lattice vector a2_p(:) = ', this%a2_p(:)
    write(fid,'(a,3f16.12)') ' primitive lattice vector a3_p(:) = ', this%a3_p(:)
    write(fid,'(a,3f16.12)') ' primitive reciprocal vector bz1_p(:) / 2pi = ', this%bz1_p(:)/(2.d0*pi)
    write(fid,'(a,3f16.12)') ' primitive reciprocal vector bz2_p(:) / 2pi = ', this%bz2_p(:)/(2.d0*pi)
    write(fid,'(a,3f16.12)') ' primitive reciprocal vector bz3_p(:) / 2pi = ', this%bz3_p(:)/(2.d0*pi)
    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice sites list '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     list(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,3i4)') i,  this%list(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice sites cord '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     cord(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,3i4)') i,  this%cord(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice nnlist '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     nnlist(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,12i4)') i,  this%nnlist(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice snlist '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     snlist(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,12i4)') i,  this%snlist(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice tnlist '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     tnlist(i,:) '
    do i = 1, this%nsites
        write(fid,'(i6,12i4)') i,  this%tnlist(i,:)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The lattice plq_cord '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     plq_cord(:,i) '
    do i = 1, this%ncell
        write(fid,'(i6,8i4)') i,  this%plq_cord(:,i)
    end do

    write(fid,*)
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)')' The momentum listk '
    write(fid,'(a)')' --------------------- '
    write(fid,'(a)') '   i     listk(i,:) '
    do i = 1, this%ncell
        write(fid,'(i6,3i4)') i,  this%listk(i,:)
    end do
  
#IFDEF PLEVEL2
    if( this%ncell < 50 ) then
    write(fid, *)
    write(fid,'(a)') '----------------------------'
    write(fid,'(a)') ' latt_imj info   '
    write(fid,'(a)') '----------------------------'
    write(fid, '(a)') '   i   j   imj '
    do j = 1, this%ncell
        do i = 1, this%ncell
            imj  = this%imj(i,j)
            write(fid, '(3i4)') i, j, imj
        end do
    end do
    end if
#ENDIF
  end subroutine print_latt3d_sub

end module latt3d_class


module cubic_lattice
   use constants, only: dp
   use latt3d_class
   type, public, extends(latt3d) :: cubic
      contains 
        procedure, public :: setup_cubic => setup_cubic_sub 
   end type cubic     

  contains

  subroutine setup_cubic_sub(this, l1, l2, l3, a1_p, a2_p, a3_p)

    class(cubic) :: this
    integer, intent(in) :: l1, l2, l3
    real(dp), dimension(:), intent(in) :: a1_p, a2_p, a3_p

    real(dp) :: ai_p(3), xk_p(3), b1_p(3), b2_p(3), b3_p(3), bz1_p(3), bz2_p(3), bz3_p(3), b_p(3) 
    real(dp) :: mat(3,3), mat_inv(3,3)

    call this%setup_latt(l1, l2, l3, a1_p, a2_p, a3_p)

    this%nsub = 1
    this%nsites = l1*l2*l3
    this%z_plq = 8              ! sites contained in one cubic, 可重复计算. 
    this%nn_nf = 6
    this%nn_lf = l1*l2*l3/2 
    allocate( this%rcord( this%nsites,3) )
    allocate( this%nnlf_list(this%nn_lf) )
    allocate( this%list( this%nsites, 2) )
    allocate( this%invlist(l1, l2, l3,this%nsub))

    ind = 0
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            do i3 = 1, this%l3
                this%list( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 1) = (i3-1)*(l1*l2)+(i1-1)*l2+i2
                this%list( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 2) = 1              ! orbital
                this%rcord( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 1) = dble(i1)
                this%rcord( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 2) = dble(i2) 
                this%rcord( (i3-1)*(l1*l2)+(i1-1)*l2+i2, 3) = dble(i3) 
                this%invlist(i1,i2,i3,1) = (i3-1)*(l1*l2)+(i1-1)*l2+i2
                if ( mod(i1+i2+i3,2) == 0 ) then                    !??
                    ind = ind + 1
                    this%nnlf_list(ind) = (i3-1)*(l1*l2)+(i1-1)*l2+i2
                end if 
            end do
        end do
    end do        
    this%cord = this%ulist


    !set nearest neighbor
    allocate( this%nnvecr(3, 6, this%nsites))
    allocate( this%nnveci(3, 6, this%nsites))
    allocate( this%nnlist(this%nsites, 6))
    allocate( this%snlist(this%nsites, 12))
    allocate( this%tnlist(this%nsites, 8))

    ! A->B
    do i = 1, this%nsites, 2
        this%nnvecr(:,1,i) = (/ 1.d0, 0.d0, 0.d0/)   ! in the unit of a1_p and a2_p
        this%nnvecr(:,2,i) = (/ 0.d0, 1.d0, 0.d0/)   
        this%nnvecr(:,3,i) = (/ 0.d0, 0.d0, 1.d0/)   
        this%nnvecr(:,4,i) = (/-1.d0, 0.d0, 0.d0/)   
        this%nnvecr(:,5,i) = (/ 0.d0,-1.d0, 0.d0/)   
        this%nnvecr(:,6,i) = (/ 0.d0, 0.d0,-1.d0/)   

        this%nnveci(:, 1, i) = (/ 1, 0, 0/) 
        this%nnveci(:, 2, i) = (/ 0, 1, 0/)
        this%nnveci(:, 3, i) = (/ 0, 0, 1/)
        this%nnveci(:, 4, i) = (/-1, 0, 0/) 
        this%nnveci(:, 5, i) = (/ 0,-1, 0/) 
        this%nnveci(:, 6, i) = (/ 0, 0,-1/) 
    end do    

    ! B->A
    do i = 2, this%nsites, 2
        this%nnvecr(:,1,i) = (/ 1.d0, 0.d0, 0.d0/)   ! in the unit of a1_p and a2_p
        this%nnvecr(:,2,i) = (/ 0.d0, 1.d0, 0.d0/)   
        this%nnvecr(:,3,i) = (/ 0.d0, 0.d0, 1.d0/)   
        this%nnvecr(:,4,i) = (/-1.d0, 0.d0, 0.d0/)   
        this%nnvecr(:,5,i) = (/ 0.d0,-1.d0, 0.d0/)   
        this%nnvecr(:,6,i) = (/ 0.d0, 0.d0,-1.d0/)   

        this%nnveci(:, 1, i) = (/ 1, 0, 0/) 
        this%nnveci(:, 2, i) = (/ 0, 1, 0/)
        this%nnveci(:, 3, i) = (/ 0, 0, 1/)
        this%nnveci(:, 4, i) = (/-1, 0, 0/) 
        this%nnveci(:, 5, i) = (/ 0,-1, 0/) 
        this%nnveci(:, 6, i) = (/ 0, 0,-1/) 
    end do 

    do i = 1, this%ncell
        ! nearest neighbor
        this%nnlist(i,1) = ( npbc(this%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),   l1) - 1 )*l2 + npbc(this%ulist(i,2),   l2) 
        this%nnlist(i,2) = ( npbc(this%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1, l1) - 1 )*l2 + npbc(this%ulist(i,2),   l2) 
        this%nnlist(i,3) = ( npbc(this%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),   l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%nnlist(i,4) = ( npbc(this%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),   l1) - 1 )*l2 + npbc(this%ulist(i,2),   l2) 
        this%nnlist(i,5) = ( npbc(this%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1, l1) - 1 )*l2 + npbc(this%ulist(i,2),   l2) 
        this%nnlist(i,6) = ( npbc(this%ulist(i,3)  ,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),   l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 

        ! next nearest neighbor 
        this%snlist(i,1) = ( npbc(this%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1,l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%snlist(i,2) = ( npbc(this%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1,l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%snlist(i,3) = ( npbc(this%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1,l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 
        this%snlist(i,4) = ( npbc(this%ulist(i,3), l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1,l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 

        this%snlist(i,5) = ( npbc(this%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%snlist(i,6) = ( npbc(this%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%snlist(i,7) = ( npbc(this%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 
        this%snlist(i,8) = ( npbc(this%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1),l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 

        this%snlist(i,9) =  ( npbc(this%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1,l1) - 1 )*l2 + npbc(this%ulist(i,2), l2)  
        this%snlist(i,10) = ( npbc(this%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1,l1) - 1 )*l2 + npbc(this%ulist(i,2), l2) 
        this%snlist(i,11) = ( npbc(this%ulist(i,3)-1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1,l1) - 1 )*l2 + npbc(this%ulist(i,2), l2) 
        this%snlist(i,12) = ( npbc(this%ulist(i,3)+1, l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1,l1) - 1 )*l2 + npbc(this%ulist(i,2), l2) 

        ! next next nearest neighbor
        this%tnlist(i,1) = ( npbc(this%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1, l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%tnlist(i,2) = ( npbc(this%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1, l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%tnlist(i,3) = ( npbc(this%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1, l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 
        this%tnlist(i,4) = ( npbc(this%ulist(i,3)+1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1, l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 
        this%tnlist(i,5) = ( npbc(this%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1, l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%tnlist(i,6) = ( npbc(this%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1, l1) - 1 )*l2 + npbc(this%ulist(i,2)+1, l2) 
        this%tnlist(i,7) = ( npbc(this%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)-1, l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 
        this%tnlist(i,8) = ( npbc(this%ulist(i,3)-1,l3) - 1 )*(l1*l2) + ( npbc(this%ulist(i,1)+1, l1) - 1 )*l2 + npbc(this%ulist(i,2)-1, l2) 
    end do


    ! set plq_cord          !??
    allocate( this%plq_cord(8, this%ncell))
    id1 = 0; id2 = 0; id3 = 0; id4 = 0; id5 = 0; id6 = 0; id7 = 0; id8 = 0; 
    do i1 = 1, this%l1
        do i2 = 1, this%l2
            do i3 = 1, this%l3
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
    do i1 = 1, this%l1 
        do i2 = 1, this%l2
            do i3 = 1, this%l3

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

                this%plq_cord(1, id) = this%invlist(i1, i2, i3,1)         
                this%plq_cord(2, id) = this%invlist(npbc(i1+1,l1), i2, i3,1)
                this%plq_cord(3, id) = this%invlist(npbc(i1+1,l1), npbc(i2+1,l2), i3,1)
                this%plq_cord(4, id) = this%invlist(i1, npbc(i2+1,l2), i3,1)
                this%plq_cord(5, id) = this%invlist(i1, i2, npbc(i3+1, l3),1)         
                this%plq_cord(6, id) = this%invlist(npbc(i1+1,l1), i2, npbc(i3+1, l3),1)
                this%plq_cord(7, id) = this%invlist(npbc(i1+1,l1), npbc(i2+1,l2), npbc(i3+1,l3),1)
                this%plq_cord(8, id) = this%invlist(i1, npbc(i2+1,l2), npbc(i3+1,l3),1)

            end do 
        end do
    end do 

  end subroutine setup_cubic_sub

end module cubic_lattice
