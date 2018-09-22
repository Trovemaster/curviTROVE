function asin_x(x) result(f)

  type(adf_realq), intent(in) :: x
  type(adf_realq) :: f

  integer(ik) :: iterm, ielem, iprod, ind, order, deg, l
  real(ark) :: prod, d

  if (trim(x%except)/='') then
    write(out, '(/a,1x,a)') 'asin(f(x)) error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    stop
  endif

  if (abs(x%v)>=1.0_ark-small_number_check.and.abs(x%v)<=1.0_ark+small_number_check) then
#ifdef _ADF_EXPR_
    call expr('asin', x, f)
#endif
    f%v = x%v
    f%d(1:adf_nterms) = x%d(1:adf_nterms)
#ifdef _ADF_OPT2_
    where(x%ivar(1:adf_nvar)==1)
      f%ivar(1:adf_nvar) = 1
    elsewhere
      f%ivar(1:adf_nvar) = 0
    endwhere
#endif
#ifdef _ADF_OPT1_
    where(abs(f%d(1:adf_nterms))<=small_number)
      f%ifd0(1:adf_nterms) = .true.
    elsewhere
      f%ifd0(1:adf_nterms) = .false.
    endwhere
#endif
    f%except = 'derivatives(asin(x=+-1))=Infinity'
    return
  endif

  f%v = asin(x%v)
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, 'asin', x)
#endif

  !omp parallel do private(iterm,ielem,prod,order,d,iprod,deg,ind) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
    do ielem=1, adf_intrf(iterm)%nelem

      prod = adf_intrf(iterm)%chrule(ielem)%coef

      order = adf_intrf(iterm)%chrule(ielem)%order(1)
      d = -deriv_acos(order, x%v)
      prod = prod * d

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
          write(out, '(/a)') 'asin_x error: x**0 must not appear in the chain rule expression for derivative of acos(x)'
          stop
        endif
      enddo

      f%d(iterm) = f%d(iterm) + prod

    enddo
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), 'asin', x)
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
  call expr('asin', x, f)
#endif

end function asin_x
