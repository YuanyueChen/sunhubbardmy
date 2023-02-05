subroutine green_equaltime( ul1, ur1, ulrinv1, gr, grc )
  implicit none
  complex(dp), dimension(ne,ndim), intent(in) :: ul1
  complex(dp), dimension(ndim,ne), intent(in) :: ur1
  complex(dp), dimension(ne,ne), intent(in) :: ulrinv1
  complex(dp), dimension(ndim,ndim), intent(out) :: gr, grc

  ! local
  integer :: i, j
  complex(dp), dimension(:,:), allocatable :: ur1tmp1

  allocate( ur1tmp1( ndim, ne ) )

  call zgemm('n','n',ndim,ne,ne,cone,ur1,ndim,ulrinv1,ne,czero,ur1tmp1,ndim)  ! ur1*ulrinv1, (ndim,ne) * (ne,ne)
  call zgemm('n','n',ndim,ndim,ne,cone,ur1tmp1,ndim,ul1,ne,czero,gr,ndim)  ! ur1tmp1*ul1, (ndim,ne) * (ne,ndim)
  gr  = Imat - gr
  do i = 1, ndim
      do j = 1, ndim
          grc(j,i) = Imat(i,j) - gr(i,j)
      end do
  end do

  deallocate( ur1tmp1 )

end subroutine green_equaltime
