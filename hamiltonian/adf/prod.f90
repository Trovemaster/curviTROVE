function prod_xy(x, y) result(f)

  type(adf_realq), intent(in) :: x, y
  type(adf_realq) :: f

  integer(ik) :: iterm, ielem, iprod, ind1, ind2, deg1, deg2, l
  real(ark) :: prod

  if (trim(x%except)/='') then
    write(out, '(/a,1x,a)') 'f(x)*g(y) error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'g(y) =', trim(y%sexpr)
    stop
  endif

  if (trim(y%except)/='') then
    write(out, '(/a,1x,a)') 'f(x)*g(y) error:', trim(y%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'g(y) =', trim(y%sexpr)
    stop
  endif

  f%v = x%v * y%v
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1.or.y%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '*', x, y)
#endif

  !omp parallel do private(iterm,ielem,prod,deg1,deg2,ind1,ind2) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
    do ielem=1, adf_prod(iterm)%nelem

      prod = adf_prod(iterm)%chrule(ielem)%coef

      deg1 = adf_prod(iterm)%chrule(ielem)%pow(1)
      deg2 = adf_prod(iterm)%chrule(ielem)%pow(2)

      ind1 = adf_prod(iterm)%chrule(ielem)%ind(1)
      ind2 = adf_prod(iterm)%chrule(ielem)%ind(2)

#ifdef _ADF_OPT1_
      if ( abs(prod)<=small_number .or. x%ifd0(ind1) .or. y%ifd0(ind2) ) cycle
#endif

      if (deg1==1) then
        prod = prod * x%d(ind1)
      elseif (deg1>1) then
        prod = prod * x%d(ind1)**deg1
      else
        write(out, '(/a)') 'prod_xy error: x**0 must not appear in the chain rule expression for derivative of x*y'
        stop
      endif

      if (deg2==1) then
        prod = prod * y%d(ind2)
      elseif (deg2>1) then
        prod = prod * y%d(ind2)**deg2
      else
        write(out, '(/a)') 'prod_xy error: y**0 must not appear in the chain rule expression for derivative of x*y'
        stop
      endif

      f%d(iterm) = f%d(iterm) + prod

    enddo
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '*', x, y)
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
  call expr('*', x, y, f)
#endif

end function prod_xy



function prod_xy_11(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  type(adf_realq), intent(in) :: y(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'prod_xy_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    f(i) = x(i) * y(i)
  enddo

end function prod_xy_11



function prod_xa(x, y) result(f)

  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: y
  type(adf_realq) :: f

  integer(ik) :: iterm

  if (trim(x%except)/='') then
    write(out, '(/a,1x,a)') 'f(x)*a error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    write(out, '(a,1x,es16.8)') 'a =', y
    stop
  endif

  f%v = x%v * y
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '*', x, y)
#endif

  !omp parallel do private(iterm) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
#ifdef _ADF_OPT1_
    if ( abs(y)<=small_number .or. x%ifd0(iterm) ) cycle
#endif
    f%d(iterm) = x%d(iterm) * y
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '*', x, y)
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
  call expr('*', x, y, f)
#endif

end function prod_xa



function prod_xa_10(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  real(ark), intent(in) :: y
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    f(i) = x(i) * y
  enddo

end function prod_xa_10



function prod_xa_20(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  real(ark), intent(in) :: y
  type(adf_realq) :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = x(j,i) * y
    enddo
  enddo

end function prod_xa_20



function prod_xa_11(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  real(ark), intent(in) :: y(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'prod_xa_11 error: shapes of x and a do not conform'
    stop
  endif

  do i=1, size(x)
    f(i) = x(i) * y(i)
  enddo

end function prod_xa_11



function prod_ax(x, y) result(f)

  real(ark), intent(in) :: x
  type(adf_realq), intent(in) :: y
  type(adf_realq) :: f

  integer(ik) :: iterm

  if (trim(y%except)/='') then
    write(out, '(/a,1x,a)') 'a*f(x) error:', trim(y%except)
    write(out, '(a,1x,es16.8)') 'a =', x
    write(out, '(a,1x,a)') 'f(x) =', trim(y%sexpr)
    stop
  endif

  f%v = x * y%v
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(y%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '*', x, y)
#endif

  !omp parallel do private(iterm) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
#ifdef _ADF_OPT1_
    if ( abs(x)<=small_number .or. y%ifd0(iterm) ) cycle
#endif
    f%d(iterm) = y%d(iterm) * x
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '*', x, y)
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
  call expr('*', x, y, f)
#endif

end function prod_ax



function prod_ax_11(x, y) result(f)

  real(ark), intent(in) :: x(:)
  type(adf_realq), intent(in) :: y(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'prod_ax_11 error: shapes of a and y do not conform'
    stop
  endif

  do i=1, size(x)
    f(i) = x(i) * y(i)
  enddo

end function prod_ax_11



function prod_xi(x, y) result(f)

  type(adf_realq), intent(in) :: x
  integer(ik), intent(in) :: y
  type(adf_realq) :: f

  integer(ik) :: iterm

  if (trim(x%except)/='') then
    write(out, '(/a,1x,a)') 'f(x)*i error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    write(out, '(a,1x,i8)') 'i =', y
    stop
  endif

  f%v = x%v * y
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '*', x, y)
#endif

  !omp parallel do private(iterm) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
#ifdef _ADF_OPT1_
    if ( abs(y)==0 .or. x%ifd0(iterm) ) cycle
#endif
    f%d(iterm) = x%d(iterm) * y
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '*', x, y)
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
  call expr('*', x, y, f)
#endif

end function prod_xi



function prod_xi_10(x, y) result(f)

  type(adf_realq), intent(in) :: x(:)
  integer(ik), intent(in) :: y
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    f(i) = x(i) * y
  enddo

end function prod_xi_10



function prod_xi_20(x, y) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  integer(ik), intent(in) :: y
  type(adf_realq) :: f(size(x,1),size(x,2))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(j,i) = x(j,i) * y
    enddo
  enddo

end function prod_xi_20



function prod_ix(x, y) result(f)

  integer(ik), intent(in) :: x
  type(adf_realq), intent(in) :: y
  type(adf_realq) :: f

  integer(ik) :: iterm

  if (trim(y%except)/='') then
    write(out, '(/a,1x,a)') 'i*f(x) error:', trim(y%except)
    write(out, '(a,1x,i8)') 'i =', x
    write(out, '(a,1x,a)') 'f(x) =', trim(y%sexpr)
    stop
  endif

  f%v = x * y%v
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(y%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '*', x, y)
#endif

  !omp parallel do private(iterm) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
#ifdef _ADF_OPT1_
    if ( abs(x)==0 .or. y%ifd0(iterm) ) cycle
#endif
    f%d(iterm) = y%d(iterm) * x
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '*', x, y)
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
  call expr('*', x, y, f)
#endif

end function prod_ix
