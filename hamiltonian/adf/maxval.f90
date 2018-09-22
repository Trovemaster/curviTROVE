function maxval_x1(x) result(f)

  type(adf_realq), intent(in) :: x(:)
  type(adf_realq) :: f

  integer(ik) :: i, imax
  real(ark) :: vmax

  imax = 1
  vmax = x(1)%v
  do i=1, size(x)
    if (x(i)%v>vmax) then
      vmax = x(i)%v
      imax = i
    endif
  enddo

  f = x(imax)

end function maxval_x1



function maxval_x2(x) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  type(adf_realq) :: f

  integer(ik) :: i, j, imax, jmax
  real(ark) :: vmax

  imax = 1
  jmax = 1
  vmax = x(1,1)%v
  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      if (x(j,i)%v>vmax) then
        vmax = x(j,i)%v
        imax = i
        jmax = j
      endif
    enddo
  enddo

  f = x(jmax,imax)

end function maxval_x2
