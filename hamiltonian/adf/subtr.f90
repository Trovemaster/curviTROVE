function subtr_xy(x, y) result(f)

  type(adf_realq), intent(in) :: x, y
  type(adf_realq) :: f

  integer(ik) :: iterm

  if (trim(x%except)/='') then
    write(out, '(/a,1x,a)') 'f(x)-g(y) error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'g(y) =', trim(y%sexpr)
    stop
  endif

  if (trim(y%except)/='') then
    write(out, '(/a,1x,a)') 'f(x)-g(y) error:', trim(y%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'g(y) =', trim(y%sexpr)
    stop
  endif

  f%v = x%v - y%v

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1.or.y%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '-', x, y)
#endif

  f%d(1:adf_nterms) = x%d(1:adf_nterms) - y%d(1:adf_nterms)

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '-', x, y)
  enddo
#endif

#ifdef _ADF_OPT1_
  where(abs(f%d(1:adf_nterms))<=small_number)
    f%ifd0(1:adf_nterms) = .true.
  elsewhere
    f%ifd0(1:adf_nterms) = .false.
  endwhere
#endif

#ifdef _ADF_EXPR_
  call expr('-', x, y, f)
#endif

end function subtr_xy



function subtr_xy_11(x, y) result(f)

  type(adf_realq), intent(in) :: x(:), y(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'subtr_xy_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    f(i) = x(i) - y(i)
  enddo

end function subtr_xy_11



function subtr_xy_22(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:), y(:,:)
  type(adf_realq) :: f(size(x,dim=1),size(x,dim=2))

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'subtr_xy_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = x(j,i) - y(j,i)
    enddo
  enddo

end function subtr_xy_22



function subtr_0x(x) result(f)

  type(adf_realq), intent(in) :: x
  type(adf_realq) :: f

  integer(ik) :: iterm

  f%except = x%except

  f%v = -x%v

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

  f%d(1:adf_nterms) = -x%d(1:adf_nterms)

#ifdef _ADF_OPT1_
  where(abs(f%d(1:adf_nterms))<=small_number)
    f%ifd0(1:adf_nterms) = .true.
  elsewhere
    f%ifd0(1:adf_nterms) = .false.
  endwhere
#endif

#ifdef _ADF_EXPR_
  call expr('-', 0, x, f)
#endif

end function subtr_0x



function subtr_0x_1(x) result(f)

  type(adf_realq), intent(in) :: x(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    f(i) = -x(i)
  enddo

end function subtr_0x_1



function subtr_0x_2(x) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  type(adf_realq) :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = -x(j,i)
    enddo
  enddo

end function subtr_0x_2



function subtr_xa(x, y) result(f)

  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: y
  type(adf_realq) :: f

  integer(ik) :: iterm

  if (trim(x%except)/='') then
    write(out, '(/a,1x,a)') 'f(x)-a error:', trim(x%except)
    write(out, '(a,1x,es16.8)') 'a =', y
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    stop
  endif

  f%v = x%v - y

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '-', x, y)
#endif

  !omp parallel do private(iterm) schedule(dynamic)
  do iterm=1, adf_nterms
    if (adf_order(iterm)==0) then
      f%d(iterm) = x%d(iterm) - y
    else
      f%d(iterm) = x%d(iterm)
    endif
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '-', x, y)
  enddo
#endif

#ifdef _ADF_OPT1_
  where(abs(f%d(1:adf_nterms))<=small_number)
    f%ifd0(1:adf_nterms) = .true.
  elsewhere
    f%ifd0(1:adf_nterms) = .false.
  endwhere
#endif

#ifdef _ADF_EXPR_
  call expr('-', x, y, f)
#endif

end function subtr_xa



function subtr_xa_11(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  real(ark), intent(in) :: y(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'subtr_xa_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    f(i) = x(i) - y(i)
  enddo

end function subtr_xa_11



function subtr_xa_22(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  real(ark), intent(in) :: y(:,:)
  type(adf_realq) :: f(size(x,dim=1),size(x,dim=2))

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'subtr_xa_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = x(j,i) - y(j,i)
    enddo
  enddo

end function subtr_xa_22



function subtr_ax(x, y) result(f)

  real(ark), intent(in) :: x
  type(adf_realq), intent(in) :: y
  type(adf_realq) :: f

  integer(ik) :: iterm

  if (trim(y%except)/='') then
    write(out, '(/a,1x,a)') 'a-f(x) error:', trim(y%except)
    write(out, '(a,1x,es16.8)') 'a =', x
    write(out, '(a,1x,a)') 'f(x) =', trim(y%sexpr)
    stop
  endif

  f%v = x - y%v

#ifdef _ADF_OPT2_
  where(y%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '-', x, y)
#endif

  !omp parallel do private(iterm) schedule(dynamic)
  do iterm=1, adf_nterms
    if (adf_order(iterm)==0) then
      f%d(iterm) = x - y%d(iterm)
    else
      f%d(iterm) = -y%d(iterm)
    endif
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '-', x, y)
  enddo
#endif

#ifdef _ADF_OPT1_
  where(abs(f%d(1:adf_nterms))<=small_number)
    f%ifd0(1:adf_nterms) = .true.
  elsewhere
    f%ifd0(1:adf_nterms) = .false.
  endwhere
#endif

#ifdef _ADF_EXPR_
  call expr('-', x, y, f)
#endif

end function subtr_ax



function subtr_ax_11(x, y) result(f)

  real(ark), intent(in) :: x(:)
  type(adf_realq), intent(in) :: y(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'subtr_ax_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    f(i) = x(i) - y(i)
  enddo

end function subtr_ax_11



function subtr_ax_22(x, y) result(f)

  real(ark), intent(in) :: x(:,:)
  type(adf_realq), intent(in) :: y(:,:)
  type(adf_realq) :: f(size(x,dim=1),size(x,dim=2))

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'subtr_ax_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = x(j,i) - y(j,i)
    enddo
  enddo

end function subtr_ax_22

