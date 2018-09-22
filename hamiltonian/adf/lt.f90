function lt_xy(x, y) result(f)

  type(adf_realq), intent(in) :: x, y
  logical :: f

  if (x%v.lt.y%v) then
    f = .true.
  else
    f = .false.
  endif

end function lt_xy



function lt_xy_10(x, y) result(f)

  type(adf_realq), intent(in) :: x(:), y
  logical :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    if (x(i)%v.lt.y%v) then
      f(i) = .true.
    else
      f(i) = .false.
    endif
  enddo

end function lt_xy_10



function lt_xy_20(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:), y
  logical :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      if (x(j,i)%v.lt.y%v) then
        f(j,i) = .true.
      else
        f(j,i) = .false.
      endif
    enddo
  enddo

end function lt_xy_20



function lt_xy_11(x, y) result(f)

  type(adf_realq), intent(in) :: x(:), y(:)
  logical :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'lt_xy_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    if (x(i)%v.lt.y(i)%v) then
      f(i) = .true.
    else
      f(i) = .false.
    endif
  enddo

end function lt_xy_11



function lt_xy_22(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:), y(:,:)
  logical :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'lt_xy_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      if (x(j,i)%v.lt.y(j,i)%v) then
        f(j,i) = .true.
      else
        f(j,i) = .false.
      endif
    enddo
  enddo

end function lt_xy_22



function lt_xa(x, y) result(f)

  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: y
  logical :: f

  if (x%v.lt.y) then
    f = .true.
  else
    f = .false.
  endif

end function lt_xa



function lt_xa_10(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  real(ark), intent(in) :: y
  logical :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    if (x(i)%v.lt.y) then
      f(i) = .true.
    else
      f(i) = .false.
    endif
  enddo

end function lt_xa_10



function lt_xa_20(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  real(ark), intent(in) :: y
  logical :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      if (x(j,i)%v.lt.y) then
        f(j,i) = .true.
      else
        f(j,i) = .false.
      endif
    enddo
  enddo

end function lt_xa_20



function lt_xa_11(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  real(ark), intent(in) :: y(:)
  logical :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'lt_xa_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    if (x(i)%v.lt.y(i)) then
      f(i) = .true.
    else
      f(i) = .false.
    endif
  enddo

end function lt_xa_11



function lt_xa_22(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  real(ark), intent(in) :: y(:,:)
  logical :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'lt_xa_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      if (x(j,i)%v.lt.y(j,i)) then
        f(j,i) = .true.
      else
        f(j,i) = .false.
      endif
    enddo
  enddo

end function lt_xa_22



function lt_ax(x, y) result(f)

  real(ark), intent(in) :: x
  type(adf_realq), intent(in) :: y
  logical :: f

  if (x.lt.y%v) then
    f = .true.
  else
    f = .false.
  endif

end function lt_ax



function lt_ax_10(x, y) result(f)

  real(ark), intent(in) :: x(:)
  type(adf_realq), intent(in) :: y
  logical :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    if (x(i).lt.y%v) then
      f(i) = .true.
    else
      f(i) = .false.
    endif
  enddo

end function lt_ax_10



function lt_ax_20(x, y) result(f)

  real(ark), intent(in) :: x(:,:)
  type(adf_realq), intent(in) :: y
  logical :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      if (x(j,i).lt.y%v) then
        f(j,i) = .true.
      else
        f(j,i) = .false.
      endif
    enddo
  enddo

end function lt_ax_20



function lt_ax_11(x, y) result(f)

  real(ark), intent(in) :: x(:)
  type(adf_realq), intent(in) :: y(:)
  logical :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'lt_ax_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    if (x(i).lt.y(i)%v) then
      f(i) = .true.
    else
      f(i) = .false.
    endif
  enddo

end function lt_ax_11



function lt_ax_22(x, y) result(f)

  real(ark), intent(in) :: x(:,:)
  type(adf_realq), intent(in) :: y(:,:)
  logical :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'lt_ax_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      if (x(j,i).lt.y(j,i)%v) then
        f(j,i) = .true.
      else
        f(j,i) = .false.
      endif
    enddo
  enddo

end function lt_ax_22
