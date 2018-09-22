function product_x(x) result(f)

  type(adf_realq), intent(in) :: x(:)
  type(adf_realq) :: f

  integer(ik) :: i

  f = x(1)
  do i=2, size(x)
    f = f * x(i)
  enddo

end function product_x
