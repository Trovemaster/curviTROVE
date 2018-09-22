function abs_x(x) result(f)

  type(adf_realq), intent(in) :: x
  type(adf_realq) :: f

  f = sqrt(x**2)

end function abs_x



function abs_x1(x) result(f)

  type(adf_realq), intent(in) :: x(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    f(i) = sqrt(x(i)**2)
  enddo

end function abs_x1



function abs_x2(x) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  type(adf_realq) :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = sqrt(x(j,i)**2)
    enddo
  enddo

end function abs_x2
