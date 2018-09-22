function dot_product_xy(x, y) result(f)

  type(adf_realq), intent(in) :: x(:), y(:)
  type(adf_realq) :: f

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'dot_product_xy error: shapes of x and y do not conform'
    stop
  endif

  f = 0.0_ark
  do i=1, size(x)
    f = f + x(i) * y(i)
  enddo

end function dot_product_xy



function dot_product_xa(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  real(ark), intent(in) :: y(:)
  type(adf_realq) :: f

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'dot_product_xa error: shapes of x and y do not conform'
    stop
  endif

  f = 0.0_ark
  do i=1, size(x)
    f = f + x(i) * y(i)
  enddo

end function dot_product_xa



function dot_product_ax(x, y) result(f)

  real(ark), intent(in) :: x(:)
  type(adf_realq), intent(in) :: y(:)
  type(adf_realq) :: f

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'dot_product_ax error: shapes of x and y do not conform'
    stop
  endif

  f = 0.0_ark
  do i=1, size(x)
    f = f + x(i) * y(i)
  enddo

end function dot_product_ax
