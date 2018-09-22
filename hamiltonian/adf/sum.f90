function sum_x1(x) result(f)

  type(adf_realq), intent(in) :: x(:)
  type(adf_realq) :: f

  integer(ik) :: i

  f = 0.0_ark
  do i=1, size(x)
    f = f + x(i)
  enddo

end function sum_x1
