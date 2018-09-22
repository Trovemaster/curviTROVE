function ipow_x(x, ipow) result(f)

  type(adf_realq), intent(in) :: x
  integer(ik), intent(in) :: ipow
  type(adf_realq) :: f

  integer(ik) :: iterm, ielem, iprod, ind, order, deg, iorder, i
  real(ark) :: prod, coef, v
  logical :: if0

  if (trim(x%except)/='') then
    write(out, '(/a,i6,1x,a,1x,a)') 'f(x) **', ipow, 'error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    stop
  endif

  if (ipow/=0) then
    f%v = x%v**ipow
  else
    f%v = 1.0_ark
  endif
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '**', x, ipow)
#endif

  !omp parallel do private(iterm,ielem,prod,order,coef,if0,iorder,iprod,deg,ind) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
    do ielem=1, adf_intrf(iterm)%nelem

      coef = adf_intrf(iterm)%chrule(ielem)%coef
      order = adf_intrf(iterm)%chrule(ielem)%order(1)

      if (ipow>0) then

        deg = ipow
        if0 = .false.
        do iorder=1, order
          if ((deg-1)<0) then
            if0 = .true.
            exit
          endif
          coef = coef * deg
          deg = deg - 1
        enddo

        if (if0) then
          prod = 0.0
        elseif (deg==0) then
          prod = coef
        elseif (deg==1) then
          prod = coef * x%v
        else
          prod = coef * x%v**deg
        endif

      elseif (ipow<0) then

        deg = ipow
        do iorder=1, order
          coef = coef * deg
          deg = deg - 1
        enddo
        prod = coef * x%v**deg

      else

        prod = coef

      endif

#ifdef _ADF_OPT1_
      if ( abs(prod)<=small_number .or. any( (/( x%ifd0(adf_intrf(iterm)%chrule(ielem)%ind(iprod)), iprod=2, adf_intrf(iterm)%chrule(ielem)%nprod )/) ) ) cycle
#endif

      do iprod=2, adf_intrf(iterm)%chrule(ielem)%nprod
        deg = adf_intrf(iterm)%chrule(ielem)%pow(iprod)
        ind = adf_intrf(iterm)%chrule(ielem)%ind(iprod)
        if (deg==1) then
          prod = prod * x%d(ind)
        elseif (deg>1) then
          prod = prod * x%d(ind)**deg
        else
          write(out, '(/a,i3)') 'ipow_x error: x**0 must not appear in the chain rule expression for derivative of x**', ipow
          stop
        endif
      enddo

      f%d(iterm) = f%d(iterm) + prod

    enddo
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '**', x, ipow)
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
  call expr('**', x, ipow, f)
#endif

end function ipow_x



function ipow_x_10(x, ipow) result(f)

  type(adf_realq), intent(in) :: x(:)
  integer(ik), intent(in) :: ipow
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  do i=1, size(x)
    f(i) = x(i)**ipow
  enddo

end function ipow_x_10



function ipow_x_11(x, ipow) result(f)

  type(adf_realq), intent(in) :: x(:)
  integer(ik), intent(in) :: ipow(:)
  type(adf_realq) :: f(size(x))

  integer(ik) :: i

  if (size(x)/=size(ipow)) then
    write(out, '(/a)') 'ipow_x_11 error: shapes of x and ipow do not conform'
    stop
  endif

  do i=1, size(x)
    f(i) = x(i)**ipow(i)
  enddo

end function ipow_x_11



function qpow_x(x, pow) result(f)

  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: pow
  type(adf_realq) :: f

  integer(ik) :: iterm, ielem, iprod, ind, order, iorder, i
  real(ark) :: prod, coef, deg

  if (trim(x%except)/='') then
    write(out, '(/a,f8.3,1x,a,1x,a)') 'f(x) **', pow, 'error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    stop
  endif

  if (abs(pow)>small_number_check) then
    f%v = x%v**pow
  else
    f%v = 1.0_ark
  endif
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, '**', x, pow)
#endif

  !omp parallel do private(iterm,ielem,prod,order,coef,iorder,iprod,deg,ind) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
    do ielem=1, adf_intrf(iterm)%nelem

      coef = adf_intrf(iterm)%chrule(ielem)%coef
      order = adf_intrf(iterm)%chrule(ielem)%order(1)

      if (abs(pow)>small_number_check) then

        deg = pow
        do iorder=1, order
          coef = coef * deg
          deg = deg - 1.0_ark
        enddo
        prod = coef * x%v**deg

      else

        prod = coef

      endif

#ifdef _ADF_OPT1_
      if ( abs(prod)<=small_number .or. any( (/( x%ifd0(adf_intrf(iterm)%chrule(ielem)%ind(iprod)), iprod=2, adf_intrf(iterm)%chrule(ielem)%nprod )/) ) ) cycle
#endif

      do iprod=2, adf_intrf(iterm)%chrule(ielem)%nprod
        deg = adf_intrf(iterm)%chrule(ielem)%pow(iprod)
        ind = adf_intrf(iterm)%chrule(ielem)%ind(iprod)
        if (deg==1) then
          prod = prod * x%d(ind)
        elseif (deg>1) then
          prod = prod * x%d(ind)**deg
        else
          write(out, '(/a,f8.3)') 'qpow_x error: x**0 must not appear in the chain rule expression for derivative of x **', pow
          stop
        endif
      enddo

      f%d(iterm) = f%d(iterm) + prod

    enddo
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), '**', x, pow)
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
  call expr('**', x, pow, f)
#endif

end function qpow_x
