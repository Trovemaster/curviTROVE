function acos_x(x) result(f)

  type(adf_realq), intent(in) :: x
  type(adf_realq) :: f

  integer(ik) :: iterm, ielem, iprod, ind, order, deg, l
  real(ark) :: prod, d

  if (trim(x%except)/='') then
    write(out, '(/a,1x,a)') 'acos(f(x)) error:', trim(x%except)
    write(out, '(a,1x,a)') 'f(x) =', trim(x%sexpr)
    stop
  endif

  if (abs(x%v)>=1.0_ark-small_number_check.and.abs(x%v)<=1.0_ark+small_number_check) then
#ifdef _ADF_EXPR_
    call expr('acos', x, f)
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
    f%except = 'derivatives(acos(x=+-1))=Infinity'
    return
  endif

  f%v = acos(x%v)
  f%d(1:adf_nterms) = 0.0

#ifdef _ADF_OPT2_
  where(x%ivar(1:adf_nvar)==1)
    f%ivar(1:adf_nvar) = 1
  elsewhere
    f%ivar(1:adf_nvar) = 0
  endwhere
#endif

#ifdef _ADF_NANINF_
  call check_nan_inf(f%v, 'acos', x)
#endif

  !omp parallel do private(iterm,ielem,prod,order,d,iprod,deg,ind) schedule(dynamic)
  do iterm=1, adf_nterms
#ifdef _ADF_OPT2_
    if (any(f%ivar(1:adf_nvar)-adf_terms_ivar(1:adf_nvar,iterm)<0)) cycle
#endif
    do ielem=1, adf_intrf(iterm)%nelem

      prod = adf_intrf(iterm)%chrule(ielem)%coef

      order = adf_intrf(iterm)%chrule(ielem)%order(1)
      d = deriv_acos(order, x%v)
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
          write(out, '(/a)') 'acos_x error: x**0 must not appear in the chain rule expression for derivative of acos(x)'
          stop
        endif
      enddo

      f%d(iterm) = f%d(iterm) + prod

    enddo
  enddo
  !omp end parallel do

#ifdef _ADF_NANINF_
  do iterm=1, adf_nterms
    call check_nan_inf(adf_terms(:,iterm), f%d(iterm), 'acos', x)
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
  call expr('acos', x, f)
#endif

end function acos_x



function deriv_acos(n, x) result(f)
  integer(ik), intent(in) :: n
  real(ark), intent(in) :: x
  real(ark) :: f
  select case(n)
  case(0)
    f = acos(x)
  case(1)
    f = -(1.0_ark/sqrt(1.0_ark - x**2))
  case(2)
    f = -(x/(1.0_ark - x**2)**1.5_ark)
  case(3)
    f = (-1.0_ark - 2*x**2)/(1.0_ark - x**2)**2.5_ark
  case(4)
    f = -((9.0_ark*x + 6.0_ark*x**3)/(1.0_ark - x**2)**3.5_ark)
  case(5)
    f = (-3.0_ark*(3.0_ark + 24.0_ark*x**2 + 8.0_ark*x**4))/(1.0_ark - x**2)**4.5_ark
  case(6)
    f = (-15.0_ark*x*(15.0_ark + 40.0_ark*x**2 + 8.0_ark*x**4))/(1.0_ark - x**2)**5.5_ark
  case(7)
    f = (-45.0_ark*(5.0_ark + 90.0_ark*x**2 + 120.0_ark*x**4 + 16.0_ark*x**6))/(1.0_ark - x**2)**6.5_ark
  case(8)
    f = (-315.0_ark*x*(35.0_ark + 210.0_ark*x**2 + 168.0_ark*x**4 + 16.0_ark*x**6))/(1.0_ark - x**2)**7.5_ark
  case(9)
    f = (-315.0_ark*(35.0_ark + 1120.0_ark*x**2 + 3360.0_ark*x**4 + 1792.0_ark*x**6 + 128.0_ark*x**8))/(1.0_ark - x**2)**8.5_ark
  case(10)
    f = (-2835.0_ark*x*(315.0_ark + 3360.0_ark*x**2 + 6048.0_ark*x**4 + 2304.0_ark*x**6 + 128.0_ark*x**8))/(1.0_ark - x**2)**9.5_ark
  case(11)
    f = (-14175.0_ark*(63.0_ark + 3150.0_ark*x**2 + 16800.0_ark*x**4 + 20160.0_ark*x**6 + 5760.0_ark*x**8 + 256.0_ark*x**10))&
        /(1.0_ark - x**2)**10.5_ark
  case(12)
    f = (-155925.0_ark*x*(693.0_ark + 11550.0_ark*x**2 + 36960.0_ark*x**4 + 31680.0_ark*x**6 + 7040.0_ark*x**8 + 256.0_ark*x**10))&
        /(1.0_ark - x**2)**11.5_ark
  case(13)
    f = (-467775.0_ark*(231.0_ark + 16632.0_ark*x**2 + 138600.0_ark*x**4 + 295680.0_ark*x**6 + 190080.0_ark*x**8 + 33792.0_ark*x**10&
         + 1024.0_ark*x**12))/(1.0_ark - x**2)**12.5_ark
  case(14)
    f = (-6081075.0_ark*x*(3003.0_ark + 72072.0_ark*x**2 + 360360.0_ark*x**4 + 549120.0_ark*x**6 + 274560.0_ark*x**8 + 39936.0_ark*x**10&
         + 1024.0_ark*x**12))/(1.0_ark - x**2)**13.5_ark
  case(15)
    f = (-42567525.0_ark*(429.0_ark + 42042.0_ark*x**2 + 504504.0_ark*x**4 + 1681680.0_ark*x**6 + 1921920.0_ark*x**8 + 768768.0_ark*x**10&
         + 93184.0_ark*x**12 + 2048.0_ark*x**14))/(1.0_ark - x**2)**14.5_ark
  case(16)
    f = (-638512875.0_ark*x*(6435.0_ark + 210210.0_ark*x**2 + 1513512.0_ark*x**4 + 3603600.0_ark*x**6 + 3203200.0_ark*x**8 + 1048320.0_ark*x**10&
         + 107520.0_ark*x**12 + 2048.0_ark*x**14))/(1.0_ark - x**2)**15.5_ark
  case(17)
    f = (-638512875.0_ark*(6435.0_ark + 823680.0_ark*x**2 + 13453440.0_ark*x**4 + 64576512.0_ark*x**6 + 115315200.0_ark*x**8 + 82001920.0_ark*x**10&
         + 22364160.0_ark*x**12 + 1966080.0_ark*x**14 + 32768.0_ark*x**16))/(1.0_ark - x**2)**16.5_ark
  case(18)
    f = (-10854718875.0_ark*x*(109395.0_ark + 4667520.0_ark*x**2 + 45741696.0_ark*x**4 + 156828672.0_ark*x**6 + 217817600.0_ark*x**8&
         + 126730240.0_ark*x**10 + 29245440.0_ark*x**12 + 2228224.0_ark*x**14 + 32768.0_ark*x**16))/(1.0_ark - x**2)**17.5_ark
  case(19)
    f = (-97692469875.0_ark*(12155.0_ark + 1969110.0_ark*x**2 + 42007680.0_ark*x**4 + 274450176.0_ark*x**6 + 705729024.0_ark*x**8&
         + 784143360.0_ark*x**10 + 380190720.0_ark*x**12 + 75202560.0_ark*x**14 + 5013504.0_ark*x**16 + 65536.0_ark*x**18))/(1.0_ark - x**2)**18.5_ark
  case(20)
    f = (-1856156927625.0_ark*x*(230945.0_ark + 12471030.0_ark*x**2 + 159629184.0_ark*x**4 + 744936192.0_ark*x**6 + 1489872384.0_ark*x**8&
         + 1354429440.0_ark*x**10 + 555663360.0_ark*x**12 + 95256576.0_ark*x**14 + 5603328.0_ark*x**16 + 65536.0_ark*x**18))/(1.0_ark - x**2)**19.5_ark
  case default
    write(out, '(/a,1x,i3)') 'deriv_acos error: derivative order for acos(x) exceeds maximum =', 20
    stop
  end select
end function deriv_acos
