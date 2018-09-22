function div_xy(x, y) result(f)

  type(adf_realq), intent(in) :: x, y
  type(adf_realq) :: f

  f = x * y**(-1)

end function div_xy



function div_xy_10(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  type(adf_realq), intent(in) :: y
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    f(i) = x(i) * y**(-1)
  enddo

end function div_xy_10



function div_xa(x, y) result(f)

  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: y
  type(adf_realq) :: f

  f = x * y**(-1)

end function div_xa



function div_xa_10(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  real(ark), intent(in) :: y
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    f(i) = x(i) * y**(-1)
  enddo

end function div_xa_10



function div_xa_20(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  real(ark), intent(in) :: y
  type(adf_realq) :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = x(j,i) * y**(-1)
    enddo
  enddo

end function div_xa_20



function div_ax(x, y) result(f)

  real(ark), intent(in) :: x
  type(adf_realq), intent(in) :: y
  type(adf_realq) :: f

  f = x * y**(-1)

end function div_ax
