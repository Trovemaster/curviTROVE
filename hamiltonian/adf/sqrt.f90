function sqrt_x(x) result(f)

  type(adf_realq), intent(in) :: x
  type(adf_realq) :: f

  f = x**0.5_ark

end function sqrt_x
