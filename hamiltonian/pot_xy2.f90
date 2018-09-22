subroutine poten_xy2_morbid_ADF(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)


  type(adf_realq) :: r12,r32,cosalpha,y1,y3,coro
  type(adf_realq) :: v0,fe1,fe3,fe11,fe33,fe13,fe111,fe333,fe113,fe133,fe1111,fe3333,fe1113,fe1333,fe1133

  integer(ik) :: irank
  real(ark) :: rhoe,re12,aa1,ve,f1,f11,f13,f111,f113,f1111,f1113,f1133,fa1,fa2,fa3,fa4,fa5,fa6,fa7,fa8, &!
              f1a1,f2a1,f3a1,f4a1,f1a11,f2a11,f3a11,f1a13,f2a13,f3a13,f1a111,f2a111,f1a113,f2a113,f1a1111,f1a1113,f1a1133, &!
              f1a3,f2a3,f3a3,f4a3,f1a33,f2a33,f3a33,f1a333,f2a333,f1a133,f2a133,f1a3333,f1a1333,f3,f33,f333,f133,f3333,f1333

  irank = 1

  ! define internal coordinates

  if (any(trim(molec%coord_transform)==(/'XY2_RRHO','XY2_RRHO_SING','XY2_RRHO_ZMAT','XY2_RRHO_TESTSING'/))) then

    r12 = internal(1)
    r32 = internal(2)
    !cosalpha = cos(real(pi,ark)-internal(3))
    cosalpha = -cos(internal(3))

  else

    write(out, '(/a,a,a)') 'poten_xy2_morbid_ADF error: coordinate type = "', trim(molec%coord_transform), '" is not supported'
    stop

  endif

  re12  = func%params(1,irank)
  rhoe  = real(pi,ark) - func%params(2,irank)*real(pi,ark)/180.0_ark
  aa1   = func%params(3,irank)

  ve        = func%params( 4,irank)
  fa1       = func%params( 5,irank)
  fa2       = func%params( 6,irank)
  fa3       = func%params( 7,irank)
  fa4       = func%params( 8,irank)
  fa5       = func%params( 9,irank)
  fa6       = func%params(10,irank)
  fa7       = func%params(11,irank)
  fa8       = func%params(12,irank)
  f1a1      = func%params(13,irank)
  f2a1      = func%params(14,irank)
  f3a1      = func%params(15,irank)
  f4a1      = func%params(16,irank)
  f11       = func%params(17,irank)
  f1a11     = func%params(18,irank)
  f2a11     = func%params(19,irank)
  f3a11     = func%params(20,irank)
  f13       = func%params(21,irank)
  f1a13     = func%params(22,irank)
  f2a13     = func%params(23,irank)
  f3a13     = func%params(24,irank)
  f111      = func%params(25,irank)
  f1a111    = func%params(26,irank)
  f2a111    = func%params(27,irank)
  f113      = func%params(28,irank)
  f1a113    = func%params(29,irank)
  f2a113    = func%params(30,irank)
  f1111     = func%params(31,irank)
  f1a1111   = func%params(32,irank)
  f1113     = func%params(33,irank)
  f1a1113   = func%params(34,irank)
  f1133     = func%params(35,irank)
  f1a1133   = func%params(36,irank)

  f1     = 0.0
  f3     = f1
  f33    = f11
  f333   = f111
  f133   = f113
  f3333  = f1111
  f1333  = f1113

  f1a3    = f1a1
  f2a3    = f2a1
  f3a3    = f3a1
  f4a3    = f4a1
  f1a33   = f1a11
  f2a33   = f2a11
  f3a33   = f3a11
  f1a333  = f1a111
  f2a333  = f2a111
  f1a133  = f1a113
  f2a133  = f2a113
  f1a3333 = f1a1111
  f1a1333 = f1a1113

  if (abs(rhoe)<=epsilon(1.0_rk)) then
    coro=1.0_ark+cosalpha
  else
    coro=cos(rhoe)+cosalpha
  endif

  y1 = 1.0_ark-exp(-aa1*(r12-re12))
  y3 = 1.0_ark-exp(-aa1*(r32-re12))

  v0= ve+fa1*coro+fa2*coro**2+fa3*coro**3+fa4*coro**4+fa5*coro**5+fa6*coro**6+fa7*coro**7+fa8*coro**8
  fe1= f1+f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4
  fe3= f3+f1a3*coro+f2a3*coro**2+f3a3*coro**3+f4a3*coro**4
  fe11= f11+f1a11*coro+f2a11*coro**2+f3a11*coro**3
  fe33= f33+f1a33*coro+f2a33*coro**2+f3a33*coro**3
  fe13= f13+f1a13*coro+f2a13*coro**2+f3a13*coro**3
  fe111= f111+f1a111*coro+f2a111*coro**2
  fe333= f333+f1a333*coro+f2a333*coro**2
  fe113= f113+f1a113*coro+f2a113*coro**2
  fe133= f133+f1a133*coro+f2a133*coro**2
  fe1111= f1111+f1a1111*coro
  fe3333= f3333+f1a3333*coro
  fe1113= f1113+f1a1113*coro
  fe1333= f1333+f1a1333*coro
  fe1133= f1133+f1a1133*coro
  f(irank) = v0+fe1*y1+fe3*y3                            &!
           +fe11*y1**2+fe33*y3**2+fe13*y1*y3           &!
           +fe111*y1**3+fe333*y3**3+fe113*y1**2*y3     &!
           +fe133*y1*y3**2                             &!
           +fe1111*y1**4+fe3333*y3**4+fe1113*y1**3*y3  &!
           +fe1333*y1*y3**3+fe1133*y1**2*y3**2

end subroutine poten_xy2_morbid_ADF
