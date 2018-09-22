function matmul_xy(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:), y(:,:)
  type(adf_realq) :: f(size(x,dim=1),size(y,dim=2))

  integer(ik) :: i, j

  if (size(x,2)/=size(y,1)) then
    write(out, '(/a)') 'matmul_xy error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(y,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = dot_product(x(j,:), y(:,i))
    enddo
  enddo

end function matmul_xy



function matmul_xa(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  real(ark), intent(in) :: y(:,:)
  type(adf_realq) :: f(size(x,dim=1),size(y,dim=2))

  integer(ik) :: i, j

  if (size(x,2)/=size(y,1)) then
    write(out, '(/a)') 'matmul_xa error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(y,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = dot_product(x(j,:), y(:,i))
    enddo
  enddo

end function matmul_xa



function matmul_ax(x, y) result(f)

  real(ark), intent(in) :: x(:,:)
  type(adf_realq), intent(in) :: y(:,:)
  type(adf_realq) :: f(size(x,dim=1),size(y,dim=2))

  integer(ik) :: i, j

  if (size(x,2)/=size(y,1)) then
    write(out, '(/a)') 'matmul_xa error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(y,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = dot_product(x(j,:), y(:,i))
    enddo
  enddo

end function matmul_ax



function matmul_ax_1(x, y) result(f)

  real(ark), intent(in) :: x(:,:)
  type(adf_realq), intent(in) :: y(:)
  type(adf_realq) :: f(size(x,dim=1))

  integer(ik) :: i, j

  if (size(x,2)/=size(y,1)) then
    write(out, '(/a)') 'matmul_xa_1 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=1)
    f(i) = dot_product(x(i,:), y(:))
  enddo

end function matmul_ax_1
