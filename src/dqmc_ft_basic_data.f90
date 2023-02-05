module dqmc_ft_basic_data
  use constants, only: dp, czi
  implicit none

  complex(dp), allocatable, dimension(:,:) :: Imat, Imat_plq
  real(dp), allocatable, dimension(:) :: Ivec

  type gfunc
      complex(dp), allocatable :: orb1(:,:)
      contains
        private
        procedure, pass(this) :: gfcopy
        procedure, pass(this) :: gfadd
        procedure, pass(this) :: gfsubtract
        generic, public :: assignment(=) => gfcopy
        generic, public :: operator(+) => gfadd
        generic, public :: operator(-) => gfsubtract
        final :: deallocate_gfunc
  end type gfunc

  type zvfunc
      complex(dp), allocatable :: orb1(:)
      contains
        private
        procedure, pass(this) :: zvcopy
        procedure, pass(this) :: zvadd
        procedure, pass(this) :: zvsubtract
        generic, public :: assignment(=) => zvcopy
        generic, public :: operator(+) => zvadd
        generic, public :: operator(-) => zvsubtract
        final :: deallocate_zvfunc
  end type zvfunc

  type dfunc
      real(dp), allocatable :: orb1(:)
      contains
        private
        procedure, pass(this) :: dfcopy
        generic, public :: assignment(=) => dfcopy
        final :: deallocate_dfunc
  end type dfunc

  type dfint
      integer, allocatable :: orb1(:)
      contains
        private
        procedure, pass(this) :: dicopy
        generic, public :: assignment(=) => dicopy
        final :: deallocate_dfint
  end type dfint

  type zfunc
      complex(dp) :: orb1
      contains
        private
        procedure, pass(this) :: zcopy
        procedure, pass(this) :: zadd
        generic, public :: assignment(=) => zcopy
        generic, public :: operator(+) => zadd
  end type zfunc

  type rfunc
      real(dp) :: orb1
      contains
        private
        procedure, pass(this) :: rcopy
        procedure, pass(this) :: radd
        generic, public :: assignment(=) => rcopy
        generic, public :: operator(+) => radd
  end type rfunc

  contains

    pure subroutine rcopy(this, from)
      class(rfunc), intent(inout) :: this
      class(rfunc), intent(in) :: from
      this%orb1 = from%orb1
    end subroutine rcopy

    pure function radd(this, that) result(total)
      class(rfunc), intent(in) :: this, that
      type(rfunc) :: total
      total%orb1 = this%orb1 + that%orb1
    end function radd

    pure subroutine zcopy(this, from)
      class(zfunc), intent(inout) :: this
      class(zfunc), intent(in) :: from
      this%orb1 = from%orb1
    end subroutine zcopy

    pure function zadd(this, that) result(total)
      class(zfunc), intent(in) :: this, that
      type(zfunc) :: total
      total%orb1 = this%orb1 + that%orb1
    end function zadd

    pure subroutine gfcopy(this, from)
      class(gfunc), intent(inout) :: this
      class(gfunc), intent(in) :: from
      this%orb1 = from%orb1
    end subroutine gfcopy

    pure function gfadd(this, that) result(total)
      class(gfunc), intent(in) :: this, that
      type(gfunc) :: total
      total%orb1 = this%orb1 + that%orb1
    end function gfadd

    pure function gfsubtract(this, that) result(total)
      class(gfunc), intent(in) :: this, that
      type(gfunc) :: total
      total%orb1 = this%orb1 - that%orb1
    end function gfsubtract

    pure subroutine zvcopy(this, from)
      class(zvfunc), intent(inout) :: this
      class(zvfunc), intent(in) :: from
      this%orb1 = from%orb1
    end subroutine zvcopy

    pure function zvadd(this, that) result(total)
      class(zvfunc), intent(in) :: this, that
      type(zvfunc) :: total
      total%orb1 = this%orb1 + that%orb1
    end function zvadd

    pure function zvsubtract(this, that) result(total)
      class(zvfunc), intent(in) :: this, that
      type(zvfunc) :: total
      total%orb1 = this%orb1 - that%orb1
    end function zvsubtract

    pure subroutine dfcopy(this, from)
      class(dfunc), intent(inout) :: this
      class(dfunc), intent(in) :: from
      this%orb1 = from%orb1
    end subroutine dfcopy

    pure subroutine dicopy(this, from)
      class(dfint), intent(inout) :: this
      class(dfint), intent(in) :: from
      this%orb1 = from%orb1
    end subroutine dicopy

    subroutine allocate_gfunc( gf, n, m )
      implicit none
      integer, intent(in) :: n, m
      type(gfunc) :: gf
      allocate( gf%orb1(n,m) )
    end subroutine allocate_gfunc

    subroutine deallocate_gfunc( this )
      implicit none
      type(gfunc) :: this
      if( allocated(this%orb1) ) deallocate( this%orb1 )
      !write(*,*) " deallocate_gfunc is performed "
    end subroutine deallocate_gfunc

    subroutine allocate_zvfunc( zv, n )
      implicit none
      integer, intent(in) :: n
      type(zvfunc) :: zv
      allocate( zv%orb1(n) )
    end subroutine allocate_zvfunc

    subroutine deallocate_zvfunc( this )
      implicit none
      type(zvfunc) :: this
      if( allocated(this%orb1) ) deallocate( this%orb1 )
      !write(*,*) " deallocate_zvfunc is performed "
    end subroutine deallocate_zvfunc

    subroutine allocate_dfunc( df, n )
      implicit none
      integer, intent(in) :: n
      type(dfunc) :: df
      allocate( df%orb1(n) )
    end subroutine allocate_dfunc

    subroutine deallocate_dfunc( this )
      implicit none
      type(dfunc) :: this
      if( allocated(this%orb1) ) deallocate( this%orb1 )
      !write(*,*) " deallocate_dfunc is performed "
    end subroutine deallocate_dfunc

    subroutine allocate_dfint( df, n )
      implicit none
      integer, intent(in) :: n
      type(dfint) :: df
      allocate( df%orb1(n) )
    end subroutine allocate_dfint

    subroutine deallocate_dfint( this )
      implicit none
      type(dfint) :: this
      if( allocated(this%orb1) ) deallocate( this%orb1 )
      !write(*,*) " deallocate_dfint is performed "
    end subroutine deallocate_dfint

end module dqmc_ft_basic_data
