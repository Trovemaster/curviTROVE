! D2h-symmetry-adapted potential energy function for C2H4-type molecule

subroutine poten_c2h4_88666678_ADF(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: i, irank
  real(ark) :: rad_, re0, re1, alphae, betae, taue, amorse1, amorse2, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7), &!
               gamma1, gamma2, delta1, delta2, delta3, delta4, delta5
  type(adf_realq) :: r0, r1, r2, r3, r4, alpha1, alpha2, alpha3, alpha4, beta1, beta2, delta, db1, db2, &!
                     vshort, vlong, vdamp, r0_2, r0_4, r1_2, r1_4, y(12)

  irank = 1
  rad_ = real(pi,ark)/180.0_ark

  ! equilibrium parameters

  re0     = func%params(1,irank)
  re1     = func%params(2,irank)
  alphae  = func%params(3,irank)*rad_
  betae   = func%params(4,irank)*rad_
  taue    = func%params(5,irank)*rad_
  amorse1 = func%params(6,irank)
  amorse2 = func%params(7,irank)

  ! long-range function parameters

  gamma1 = func%params(8,irank)  ! cc stretch
  gamma2 = func%params(9,irank)  ! ch stretch

  de_str2(1) = func%params(10,irank)  ! cc stretch
  de_str4(1) = func%params(11,irank)  ! cc stretch

  de_str2(2:5) = func%params(12,irank)  ! ch stretch
  de_str4(2:5) = func%params(13,irank)  ! ch stretch

  de_bnd2(1:4) = func%params(14,irank)  ! hcc bend
  de_bnd4(1:4) = func%params(15,irank)  ! hcc bend

  de_bnd2(5:6) = func%params(16,irank)  ! beta bend
  de_bnd4(5:6) = func%params(17,irank)  ! beta bend

  de_bnd2(7) = func%params(18,irank)  ! tau bend
  de_bnd4(7) = func%params(19,irank)  ! tau bend

  ! damping-function parameters

  delta1 = func%params(20,irank)
  delta2 = func%params(21,irank)
  delta3 = func%params(22,irank)
  delta4 = func%params(23,irank)
  delta5 = func%params(24,irank)

  ! define internal coordinates

  if (trim(molec%coord_transform)=='C2H4_2BETA_1TAU') then

    r0 = internal(1)
    r1 = internal(2)
    r2 = internal(3)
    r3 = internal(4)
    r4 = internal(5)
    alpha1 = internal(6)
    alpha2 = internal(7)
    alpha3 = internal(8)
    alpha4 = internal(9)
    beta1 = internal(10)
    beta2 = internal(11)
    delta = internal(12)
    db1 = real(pi,ark) - beta1
    db2 = beta2 - real(pi,ark)

  else

    write(out, '(/a,a,a)') 'poten_c2h4_88666678_ADF error: coordinate type = "', trim(molec%coord_transform), '" is not supported'
    stop

  endif

  ! expansion functions

  y(1) = 1.0_ark-exp(-amorse1*(r0-re0))

  y(2) = 1.0_ark-exp(-amorse2*(r1-re1))
  y(3) = 1.0_ark-exp(-amorse2*(r2-re1))
  y(4) = 1.0_ark-exp(-amorse2*(r3-re1))
  y(5) = 1.0_ark-exp(-amorse2*(r4-re1))

  y(6) = alpha1 - alphae
  y(7) = alpha2 - alphae
  y(8) = alpha3 - alphae
  y(9) = alpha4 - alphae

  y(10) = db1
  y(11) = db2
  y(12) = delta

  ! short-range potential

  vshort = func%params(25,irank) &!
         + c2h4_poten_n1_d8_ADF( y, func%params(26:57,irank) ) &!
         + c2h4_poten_n2_d8_ADF( y, func%params(58:433,irank) ) &!
         + c2h4_poten_n3_d6_ADF( y, func%params(434:1149,irank) ) &!
         + c2h4_poten_n4_d6_ADF( y, func%params(1150:2167,irank) ) &!
         + c2h4_poten_n5_d6_ADF( y, func%params(2168:2755,irank) ) &!
         + c2h4_poten_n6_d6_ADF( y, func%params(2756:2875,irank) ) &!
         + c2h4_poten_n7_d7_ADF( y, func%params(2876:2983,irank) ) &!
         + c2h4_poten_n8_d8_ADF( y, func%params(2984:3050,irank) )

  ! long-range potential

  r0_2 = (r0-re0)**2
  r0_4 = (r0-re0)**4
  r1_2 = (r1-re1)**2 + (r2-re1)**2 + (r3-re1)**2 + (r4-re1)**2
  r1_4 = (r1-re1)**4 + (r2-re1)**4 + (r3-re1)**4 + (r4-re1)**4

  vlong = sum(de_str2(1:5)*y(1:5)**2 + de_str4(1:5)*y(1:5)**4) &!
        + exp( -gamma1*r0_2 -gamma2*r1_2 ) &!
        * sum( de_bnd2(1:7)*y(6:12)**2 + de_bnd4(1:7)*y(6:12)**4 )

  ! damping function

  vdamp = exp( -delta1*r0_2 -2.0_ark*delta1*r0_4 ) &!
        * exp( -delta2*r1_2 -2.0_ark*delta2*r1_4 ) &!
        * exp( -delta3*sum(y(6:9)**2)   -2.0_ark*delta3*sum(y(6:9)**4) ) &!
        * exp( -delta4*sum(y(10:11)**2) -2.0_ark*delta4*sum(y(10:11)**4) ) &!
        * exp( -delta5*sum(y(12:12)**2) -2.0_ark*delta5*sum(y(12:12)**4) )

  ! total potential

  f = vshort*vdamp + vlong

end subroutine poten_c2h4_88666678_ADF


!###############################################################################


subroutine poten_c2h4_886666_ADF(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: i, irank
  real(ark) :: rad_, re0, re1, alphae, betae, taue, amorse1, amorse2, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7), &!
               gamma1, gamma2, delta1, delta2, delta3, delta4, delta5
  type(adf_realq) :: r0, r1, r2, r3, r4, alpha1, alpha2, alpha3, alpha4, beta1, beta2, delta, db1, db2, &!
                     vshort, vlong, vdamp, r0_2, r0_4, r1_2, r1_4, y(12)

  irank = 1
  rad_ = real(pi,ark)/180.0_ark

  ! equilibrium parameters

  re0     = func%params(1,irank)
  re1     = func%params(2,irank)
  alphae  = func%params(3,irank)*rad_
  betae   = func%params(4,irank)*rad_
  taue    = func%params(5,irank)*rad_
  amorse1 = func%params(6,irank)
  amorse2 = func%params(7,irank)

  ! long-range function parameters

  gamma1 = func%params(8,irank)  ! cc stretch
  gamma2 = func%params(9,irank)  ! ch stretch

  de_str2(1) = func%params(10,irank)  ! cc stretch
  de_str4(1) = func%params(11,irank)  ! cc stretch

  de_str2(2:5) = func%params(12,irank)  ! ch stretch
  de_str4(2:5) = func%params(13,irank)  ! ch stretch

  de_bnd2(1:4) = func%params(14,irank)  ! hcc bend
  de_bnd4(1:4) = func%params(15,irank)  ! hcc bend

  de_bnd2(5:6) = func%params(16,irank)  ! beta bend
  de_bnd4(5:6) = func%params(17,irank)  ! beta bend

  de_bnd2(7) = func%params(18,irank)  ! tau bend
  de_bnd4(7) = func%params(19,irank)  ! tau bend

  ! damping-function parameters

  delta1 = func%params(20,irank)
  delta2 = func%params(21,irank)
  delta3 = func%params(22,irank)
  delta4 = func%params(23,irank)
  delta5 = func%params(24,irank)

  ! define internal coordinates

  if (trim(molec%coord_transform)=='C2H4_2BETA_1TAU') then

    r0 = internal(1)
    r1 = internal(2)
    r2 = internal(3)
    r3 = internal(4)
    r4 = internal(5)
    alpha1 = internal(6)
    alpha2 = internal(7)
    alpha3 = internal(8)
    alpha4 = internal(9)
    beta1 = internal(10)
    beta2 = internal(11)
    delta = internal(12)
    db1 = real(pi,ark) - beta1
    db2 = beta2 - real(pi,ark)

  else

    write(out, '(/a,a,a)') 'poten_c2h4_886666_ADF error: coordinate type = "', trim(molec%coord_transform), '" is not supported'
    stop

  endif

  ! expansion functions

  y(1) = 1.0_ark-exp(-amorse1*(r0-re0))

  y(2) = 1.0_ark-exp(-amorse2*(r1-re1))
  y(3) = 1.0_ark-exp(-amorse2*(r2-re1))
  y(4) = 1.0_ark-exp(-amorse2*(r3-re1))
  y(5) = 1.0_ark-exp(-amorse2*(r4-re1))

  y(6) = alpha1 - alphae
  y(7) = alpha2 - alphae
  y(8) = alpha3 - alphae
  y(9) = alpha4 - alphae

  y(10) = db1
  y(11) = db2
  y(12) = delta

  ! short-range potential

  vshort = func%params(25,irank) &!
         + c2h4_poten_n1_d8_ADF( y, func%params(26:57,irank) ) &!
         + c2h4_poten_n2_d8_ADF( y, func%params(58:433,irank) ) &!
         + c2h4_poten_n3_d6_ADF( y, func%params(434:1149,irank) ) &!
         + c2h4_poten_n4_d6_ADF( y, func%params(1150:2167,irank) ) &!
         + c2h4_poten_n5_d6_ADF( y, func%params(2168:2755,irank) ) &!
         + c2h4_poten_n6_d6_ADF( y, func%params(2756:2875,irank) )

  ! long-range potential

  r0_2 = (r0-re0)**2
  r0_4 = (r0-re0)**4
  r1_2 = (r1-re1)**2 + (r2-re1)**2 + (r3-re1)**2 + (r4-re1)**2
  r1_4 = (r1-re1)**4 + (r2-re1)**4 + (r3-re1)**4 + (r4-re1)**4

  vlong = sum(de_str2(1:5)*y(1:5)**2 + de_str4(1:5)*y(1:5)**4) &!
        + exp( -gamma1*r0_2 -gamma2*r1_2 ) &!
        * sum( de_bnd2(1:7)*y(6:12)**2 + de_bnd4(1:7)*y(6:12)**4 )

  ! damping function

  vdamp = exp( -delta1*r0_2 -2.0_ark*delta1*r0_4 ) &!
        * exp( -delta2*r1_2 -2.0_ark*delta2*r1_4 ) &!
        * exp( -delta3*sum(y(6:9)**2)   -2.0_ark*delta3*sum(y(6:9)**4) ) &!
        * exp( -delta4*sum(y(10:11)**2) -2.0_ark*delta4*sum(y(10:11)**4) ) &!
        * exp( -delta5*sum(y(12:12)**2) -2.0_ark*delta5*sum(y(12:12)**4) )

  ! total potential

  f = vshort*vdamp + vlong

end subroutine poten_c2h4_886666_ADF


!###############################################################################


subroutine dipole_c2h4_4m_ADF(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: i, j, natoms, nmodes, irank, iatom, iterm
  real(ark) :: rad_, re0, re1, alphae, betae, taue, amorse1, amorse2, coefs_tmat(3,3,adf_nterms), coefs_tinv(3,3,adf_nterms)
  type(adf_realq) :: r0, r1, r2, r3, r4, alpha1, alpha2, alpha3, alpha4, beta1, beta2, delta, db1, db2, y(12), mu(3), &!
                     cart0(3), xyz(molec%natoms,3), nc1(3), nc2(3), nh1(3), nh2(3), nh3(3), nh4(3), e1(3), e2(3), e3(3), e4(3), &!
                     e_b1u(3), e_b2u(3), e_b3u(3), tmat(3,3), tinv(3,3)

  nmodes = molec%nmodes
  natoms = molec%natoms

  if (.not.present(cart)) then
    write(out, '(/a)') 'dipole_c2h4_4m_ADF error: Cartesian-coordinate derivatives are not present'
    stop
  endif

  ! equilibrium parameters

  irank = 1
  rad_ = real(pi,ark)/180.0_ark

  re0     = func%params(1,irank)
  re1     = func%params(2,irank)
  alphae  = func%params(3,irank)*rad_
  betae   = func%params(4,irank)*rad_
  taue    = func%params(5,irank)*rad_
  amorse1 = func%params(6,irank)
  amorse2 = func%params(7,irank)

  ! define expansion functions

  if (trim(molec%coord_transform)=='C2H4_2BETA_1TAU') then

    ! internal coordinates

    r0 = internal(1)
    r1 = internal(2)
    r2 = internal(3)
    r3 = internal(4)
    r4 = internal(5)
    alpha1 = internal(6)
    alpha2 = internal(7)
    alpha3 = internal(8)
    alpha4 = internal(9)
    beta1 = internal(10)
    beta2 = internal(11)
    delta = internal(12)
    db1 = real(pi,ark) - beta1
    db2 = beta2 - real(pi,ark)

    ! expansion functions

    y(1) = (r0-re0) * exp(-amorse1*(r0-re0)**2)
    y(2) = (r1-re1) * exp(-amorse2*(r1-re1)**2)
    y(3) = (r2-re1) * exp(-amorse2*(r2-re1)**2)
    y(4) = (r3-re1) * exp(-amorse2*(r3-re1)**2)
    y(5) = (r4-re1) * exp(-amorse2*(r4-re1)**2)
    y(6) = alpha1 - alphae
    y(7) = alpha2 - alphae
    y(8) = alpha3 - alphae
    y(9) = alpha4 - alphae
    y(10) = db1
    y(11) = db2
    y(12) = delta

  else

    write(out, '(/a,a,a)') 'dipole_c2h4_4m_ADF error: coordinate type = "', trim(molec%coord_transform), '" is not supported'
    stop

  endif

  ! symmetry-adapted components

  irank = 1
  mu(irank) = c2h4_dipole_b1u_n1_d6_ADF( y, func%params(8:22,irank) ) &!
            + c2h4_dipole_b1u_n2_d6_ADF( y, func%params(23:210,irank) ) &!
            + c2h4_dipole_b1u_n3_d6_ADF( y, func%params(211:898,irank) ) &!
            + c2h4_dipole_b1u_n4_d6_ADF( y, func%params(899:1881,irank) )

  irank = 2
  mu(irank) = c2h4_dipole_b2u_n1_d6_ADF( y, func%params(8:19,irank) ) &!
            + c2h4_dipole_b2u_n2_d6_ADF( y, func%params(20:196,irank) ) &!
            + c2h4_dipole_b2u_n3_d6_ADF( y, func%params(197:876,irank) ) &!
            + c2h4_dipole_b2u_n4_d6_ADF( y, func%params(877:1861,irank) )

  irank = 3
  mu(irank) = c2h4_dipole_b3u_n1_d6_ADF( y, func%params(8:10,irank) ) &!
            + c2h4_dipole_b3u_n2_d6_ADF( y, func%params(11:79,irank) ) &!
            + c2h4_dipole_b3u_n3_d6_ADF( y, func%params(80:512,irank) ) &!
            + c2h4_dipole_b3u_n4_d6_ADF( y, func%params(513:1399,irank) )


  ! *construct transformation matrix from Cartesian to symmetry-adapted projections

  ! choose origin X to be in the middle of CC bond

  cart0 = (cart(2,:)+cart(1,:))*0.5_ark
  do iatom=1, natoms
    xyz(iatom,:) = cart(iatom,:) - cart0(:)
  enddo

  nc1 = xyz(1,:)/sqrt(sum(xyz(1,:)**2))
  nc2 = xyz(2,:)/sqrt(sum(xyz(2,:)**2))
  nh1 = xyz(3,:)/sqrt(sum(xyz(3,:)**2))
  nh2 = xyz(4,:)/sqrt(sum(xyz(4,:)**2))
  nh3 = xyz(5,:)/sqrt(sum(xyz(5,:)**2))
  nh4 = xyz(6,:)/sqrt(sum(xyz(6,:)**2))

  ! vectors perpendicular to XC2H1, XC1H2, XC2H3, and XC2H4 planes

  e1 = -vector_product_ADF(nc1,nh1)
  e2 =  vector_product_ADF(nc1,nh2)
  e3 =  vector_product_ADF(nc2,nh3)
  e4 = -vector_product_ADF(nc2,nh4)

  e_b1u = nc1;                                e_b1u=e_b1u/sqrt(sum(e_b1u**2))
  e_b3u = (e1 + e2) + (e3 + e4);              e_b3u=e_b3u/sqrt(sum(e_b3u**2))
  e_b2u = vector_product_ADF(e_b1u, e_b3u);   e_b2u=e_b2u/sqrt(sum(e_b2u**2))

  tmat(1,:) = e_b1u
  tmat(2,:) = e_b2u
  tmat(3,:) = e_b3u

  ! *find tmat^{-1} using algebraic approach

  do iterm=1, adf_nterms
    do i=1, 3
      do j=1, 3
        coefs_tmat(j,i,iterm) = tmat(j,i)%d(iterm)
      enddo
    enddo
  enddo

  call deriv_invmat(nmodes, adf_nterms, adf_terms, 3, coefs_tmat, coefs_tinv)

  ! switch to ADF approach

  do i=1, 3
    do j=1, 3
      call adf_set_var(tinv(j,i), coefs_tinv(j,i,1:adf_nterms))
    enddo
  enddo

  ! compute Cartesian projections of dipole moment

  do irank=1, 3
    f(irank) = dot_product(tinv(irank,1:3), mu)
  enddo

end subroutine dipole_c2h4_4m_ADF


!###############################################################################


! D2h-symmetry-adapted 1-mode 8-order expansion for C2H4 molecule

function c2h4_poten_n1_d8_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(32)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(32)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 8*r0
f(2) = 8*r0**2
f(3) = 8*r0**3
f(4) = 8*r0**4
f(5) = 8*r0**5
f(6) = 8*r0**6
f(7) = 8*r0**7
f(8) = 8*r0**8
f(9) = 2*r1 + 2*r2 + 2*r3 + 2*r4
f(10) = 2*r1**2 + 2*r2**2 + 2*r3**2 + 2*r4**2
f(11) = 2*r1**3 + 2*r2**3 + 2*r3**3 + 2*r4**3
f(12) = 2*r1**4 + 2*r2**4 + 2*r3**4 + 2*r4**4
f(13) = 2*r1**5 + 2*r2**5 + 2*r3**5 + 2*r4**5
f(14) = 2*r1**6 + 2*r2**6 + 2*r3**6 + 2*r4**6
f(15) = 2*r1**7 + 2*r2**7 + 2*r3**7 + 2*r4**7
f(16) = 2*r1**8 + 2*r2**8 + 2*r3**8 + 2*r4**8
f(17) = 2*a1 + 2*a2 + 2*a3 + 2*a4
f(18) = 2*a1**2 + 2*a2**2 + 2*a3**2 + 2*a4**2
f(19) = 2*a1**3 + 2*a2**3 + 2*a3**3 + 2*a4**3
f(20) = 2*a1**4 + 2*a2**4 + 2*a3**4 + 2*a4**4
f(21) = 2*a1**5 + 2*a2**5 + 2*a3**5 + 2*a4**5
f(22) = 2*a1**6 + 2*a2**6 + 2*a3**6 + 2*a4**6
f(23) = 2*a1**7 + 2*a2**7 + 2*a3**7 + 2*a4**7
f(24) = 2*a1**8 + 2*a2**8 + 2*a3**8 + 2*a4**8
f(25) = 4*b1**2 + 4*b2**2
f(26) = 4*b1**4 + 4*b2**4
f(27) = 4*b1**6 + 4*b2**6
f(28) = 4*b1**8 + 4*b2**8
f(29) = 8*dtau**2
f(30) = 8*dtau**4
f(31) = 8*dtau**6
f(32) = 8*dtau**8
v = sum(f*params)
end function c2h4_poten_n1_d8_ADF


!###############################################################################


! D2h-symmetry-adapted 2-mode 8-order expansion for C2H4 molecule

function c2h4_poten_n2_d8_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(376)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(376)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 2*r0*(r1 + r2 + r3 + r4)
f(2) = 2*r0*(r1**2 + r2**2 + r3**2 + r4**2)
f(3) = 2*r0*(r1**3 + r2**3 + r3**3 + r4**3)
f(4) = 2*r0*(r1**4 + r2**4 + r3**4 + r4**4)
f(5) = 2*r0*(r1**5 + r2**5 + r3**5 + r4**5)
f(6) = 2*r0*(r1**6 + r2**6 + r3**6 + r4**6)
f(7) = 2*r0*(r1**7 + r2**7 + r3**7 + r4**7)
f(8) = 2*r0**2*(r1 + r2 + r3 + r4)
f(9) = 2*r0**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(10) = 2*r0**2*(r1**3 + r2**3 + r3**3 + r4**3)
f(11) = 2*r0**2*(r1**4 + r2**4 + r3**4 + r4**4)
f(12) = 2*r0**2*(r1**5 + r2**5 + r3**5 + r4**5)
f(13) = 2*r0**2*(r1**6 + r2**6 + r3**6 + r4**6)
f(14) = 2*r0**3*(r1 + r2 + r3 + r4)
f(15) = 2*r0**3*(r1**2 + r2**2 + r3**2 + r4**2)
f(16) = 2*r0**3*(r1**3 + r2**3 + r3**3 + r4**3)
f(17) = 2*r0**3*(r1**4 + r2**4 + r3**4 + r4**4)
f(18) = 2*r0**3*(r1**5 + r2**5 + r3**5 + r4**5)
f(19) = 2*r0**4*(r1 + r2 + r3 + r4)
f(20) = 2*r0**4*(r1**2 + r2**2 + r3**2 + r4**2)
f(21) = 2*r0**4*(r1**3 + r2**3 + r3**3 + r4**3)
f(22) = 2*r0**4*(r1**4 + r2**4 + r3**4 + r4**4)
f(23) = 2*r0**5*(r1 + r2 + r3 + r4)
f(24) = 2*r0**5*(r1**2 + r2**2 + r3**2 + r4**2)
f(25) = 2*r0**5*(r1**3 + r2**3 + r3**3 + r4**3)
f(26) = 2*r0**6*(r1 + r2 + r3 + r4)
f(27) = 2*r0**6*(r1**2 + r2**2 + r3**2 + r4**2)
f(28) = 2*r0**7*(r1 + r2 + r3 + r4)
f(29) = 2*r0*(a1 + a2 + a3 + a4)
f(30) = 2*r0*(a1**2 + a2**2 + a3**2 + a4**2)
f(31) = 2*r0*(a1**3 + a2**3 + a3**3 + a4**3)
f(32) = 2*r0*(a1**4 + a2**4 + a3**4 + a4**4)
f(33) = 2*r0*(a1**5 + a2**5 + a3**5 + a4**5)
f(34) = 2*r0*(a1**6 + a2**6 + a3**6 + a4**6)
f(35) = 2*r0*(a1**7 + a2**7 + a3**7 + a4**7)
f(36) = 2*r0**2*(a1 + a2 + a3 + a4)
f(37) = 2*r0**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(38) = 2*r0**2*(a1**3 + a2**3 + a3**3 + a4**3)
f(39) = 2*r0**2*(a1**4 + a2**4 + a3**4 + a4**4)
f(40) = 2*r0**2*(a1**5 + a2**5 + a3**5 + a4**5)
f(41) = 2*r0**2*(a1**6 + a2**6 + a3**6 + a4**6)
f(42) = 2*r0**3*(a1 + a2 + a3 + a4)
f(43) = 2*r0**3*(a1**2 + a2**2 + a3**2 + a4**2)
f(44) = 2*r0**3*(a1**3 + a2**3 + a3**3 + a4**3)
f(45) = 2*r0**3*(a1**4 + a2**4 + a3**4 + a4**4)
f(46) = 2*r0**3*(a1**5 + a2**5 + a3**5 + a4**5)
f(47) = 2*r0**4*(a1 + a2 + a3 + a4)
f(48) = 2*r0**4*(a1**2 + a2**2 + a3**2 + a4**2)
f(49) = 2*r0**4*(a1**3 + a2**3 + a3**3 + a4**3)
f(50) = 2*r0**4*(a1**4 + a2**4 + a3**4 + a4**4)
f(51) = 2*r0**5*(a1 + a2 + a3 + a4)
f(52) = 2*r0**5*(a1**2 + a2**2 + a3**2 + a4**2)
f(53) = 2*r0**5*(a1**3 + a2**3 + a3**3 + a4**3)
f(54) = 2*r0**6*(a1 + a2 + a3 + a4)
f(55) = 2*r0**6*(a1**2 + a2**2 + a3**2 + a4**2)
f(56) = 2*r0**7*(a1 + a2 + a3 + a4)
f(57) = 4*r0*(b1**2 + b2**2)
f(58) = 4*r0*(b1**4 + b2**4)
f(59) = 4*r0*(b1**6 + b2**6)
f(60) = 4*r0**2*(b1**2 + b2**2)
f(61) = 4*r0**2*(b1**4 + b2**4)
f(62) = 4*r0**2*(b1**6 + b2**6)
f(63) = 4*r0**3*(b1**2 + b2**2)
f(64) = 4*r0**3*(b1**4 + b2**4)
f(65) = 4*r0**4*(b1**2 + b2**2)
f(66) = 4*r0**4*(b1**4 + b2**4)
f(67) = 4*r0**5*(b1**2 + b2**2)
f(68) = 4*r0**6*(b1**2 + b2**2)
f(69) = 8*dtau**2*r0
f(70) = 8*dtau**4*r0
f(71) = 8*dtau**6*r0
f(72) = 8*dtau**2*r0**2
f(73) = 8*dtau**4*r0**2
f(74) = 8*dtau**6*r0**2
f(75) = 8*dtau**2*r0**3
f(76) = 8*dtau**4*r0**3
f(77) = 8*dtau**2*r0**4
f(78) = 8*dtau**4*r0**4
f(79) = 8*dtau**2*r0**5
f(80) = 8*dtau**2*r0**6
f(81) = 4*r1*r2 + 4*r3*r4
f(82) = 2*r1**2*r2 + 2*r1*r2**2 + 2*r3**2*r4 + 2*r3*r4**2
f(83) = 2*r1**3*r2 + 2*r1*r2**3 + 2*r3**3*r4 + 2*r3*r4**3
f(84) = 2*r1**4*r2 + 2*r1*r2**4 + 2*r3**4*r4 + 2*r3*r4**4
f(85) = 2*r1**5*r2 + 2*r1*r2**5 + 2*r3**5*r4 + 2*r3*r4**5
f(86) = 2*r1**6*r2 + 2*r1*r2**6 + 2*r3**6*r4 + 2*r3*r4**6
f(87) = 2*r1**7*r2 + 2*r1*r2**7 + 2*r3**7*r4 + 2*r3*r4**7
f(88) = 4*r1**2*r2**2 + 4*r3**2*r4**2
f(89) = 2*r1**3*r2**2 + 2*r1**2*r2**3 + 2*r3**3*r4**2 + 2*r3**2*r4**3
f(90) = 2*r1**4*r2**2 + 2*r1**2*r2**4 + 2*r3**4*r4**2 + 2*r3**2*r4**4
f(91) = 2*r1**5*r2**2 + 2*r1**2*r2**5 + 2*r3**5*r4**2 + 2*r3**2*r4**5
f(92) = 2*r1**6*r2**2 + 2*r1**2*r2**6 + 2*r3**6*r4**2 + 2*r3**2*r4**6
f(93) = 4*r1**3*r2**3 + 4*r3**3*r4**3
f(94) = 2*r1**4*r2**3 + 2*r1**3*r2**4 + 2*r3**4*r4**3 + 2*r3**3*r4**4
f(95) = 2*r1**5*r2**3 + 2*r1**3*r2**5 + 2*r3**5*r4**3 + 2*r3**3*r4**5
f(96) = 4*r1**4*r2**4 + 4*r3**4*r4**4
f(97) = 4*r1*r3 + 4*r2*r4
f(98) = 2*r1**2*r3 + 2*r1*r3**2 + 2*r2**2*r4 + 2*r2*r4**2
f(99) = 2*r1**3*r3 + 2*r1*r3**3 + 2*r2**3*r4 + 2*r2*r4**3
f(100) = 2*r1**4*r3 + 2*r1*r3**4 + 2*r2**4*r4 + 2*r2*r4**4
f(101) = 2*r1**5*r3 + 2*r1*r3**5 + 2*r2**5*r4 + 2*r2*r4**5
f(102) = 2*r1**6*r3 + 2*r1*r3**6 + 2*r2**6*r4 + 2*r2*r4**6
f(103) = 2*r1**7*r3 + 2*r1*r3**7 + 2*r2**7*r4 + 2*r2*r4**7
f(104) = 4*r1**2*r3**2 + 4*r2**2*r4**2
f(105) = 2*r1**3*r3**2 + 2*r1**2*r3**3 + 2*r2**3*r4**2 + 2*r2**2*r4**3
f(106) = 2*r1**4*r3**2 + 2*r1**2*r3**4 + 2*r2**4*r4**2 + 2*r2**2*r4**4
f(107) = 2*r1**5*r3**2 + 2*r1**2*r3**5 + 2*r2**5*r4**2 + 2*r2**2*r4**5
f(108) = 2*r1**6*r3**2 + 2*r1**2*r3**6 + 2*r2**6*r4**2 + 2*r2**2*r4**6
f(109) = 4*r1**3*r3**3 + 4*r2**3*r4**3
f(110) = 2*r1**4*r3**3 + 2*r1**3*r3**4 + 2*r2**4*r4**3 + 2*r2**3*r4**4
f(111) = 2*r1**5*r3**3 + 2*r1**3*r3**5 + 2*r2**5*r4**3 + 2*r2**3*r4**5
f(112) = 4*r1**4*r3**4 + 4*r2**4*r4**4
f(113) = 4*r1*r4 + 4*r2*r3
f(114) = 2*r1**2*r4 + 2*r1*r4**2 + 2*r2**2*r3 + 2*r2*r3**2
f(115) = 2*r1**3*r4 + 2*r1*r4**3 + 2*r2**3*r3 + 2*r2*r3**3
f(116) = 2*r1**4*r4 + 2*r1*r4**4 + 2*r2**4*r3 + 2*r2*r3**4
f(117) = 2*r1**5*r4 + 2*r1*r4**5 + 2*r2**5*r3 + 2*r2*r3**5
f(118) = 2*r1**6*r4 + 2*r1*r4**6 + 2*r2**6*r3 + 2*r2*r3**6
f(119) = 2*r1**7*r4 + 2*r1*r4**7 + 2*r2**7*r3 + 2*r2*r3**7
f(120) = 4*r1**2*r4**2 + 4*r2**2*r3**2
f(121) = 2*r1**3*r4**2 + 2*r1**2*r4**3 + 2*r2**3*r3**2 + 2*r2**2*r3**3
f(122) = 2*r1**4*r4**2 + 2*r1**2*r4**4 + 2*r2**4*r3**2 + 2*r2**2*r3**4
f(123) = 2*r1**5*r4**2 + 2*r1**2*r4**5 + 2*r2**5*r3**2 + 2*r2**2*r3**5
f(124) = 2*r1**6*r4**2 + 2*r1**2*r4**6 + 2*r2**6*r3**2 + 2*r2**2*r3**6
f(125) = 4*r1**3*r4**3 + 4*r2**3*r3**3
f(126) = 2*r1**4*r4**3 + 2*r1**3*r4**4 + 2*r2**4*r3**3 + 2*r2**3*r3**4
f(127) = 2*r1**5*r4**3 + 2*r1**3*r4**5 + 2*r2**5*r3**3 + 2*r2**3*r3**5
f(128) = 4*r1**4*r4**4 + 4*r2**4*r3**4
f(129) = 2*a1*r1 + 2*a2*r2 + 2*a3*r3 + 2*a4*r4
f(130) = 2*a1**2*r1 + 2*a2**2*r2 + 2*a3**2*r3 + 2*a4**2*r4
f(131) = 2*a1**3*r1 + 2*a2**3*r2 + 2*a3**3*r3 + 2*a4**3*r4
f(132) = 2*a1**4*r1 + 2*a2**4*r2 + 2*a3**4*r3 + 2*a4**4*r4
f(133) = 2*a1**5*r1 + 2*a2**5*r2 + 2*a3**5*r3 + 2*a4**5*r4
f(134) = 2*a1**6*r1 + 2*a2**6*r2 + 2*a3**6*r3 + 2*a4**6*r4
f(135) = 2*a1**7*r1 + 2*a2**7*r2 + 2*a3**7*r3 + 2*a4**7*r4
f(136) = 2*a1*r1**2 + 2*a2*r2**2 + 2*a3*r3**2 + 2*a4*r4**2
f(137) = 2*a1**2*r1**2 + 2*a2**2*r2**2 + 2*a3**2*r3**2 + 2*a4**2*r4**2
f(138) = 2*a1**3*r1**2 + 2*a2**3*r2**2 + 2*a3**3*r3**2 + 2*a4**3*r4**2
f(139) = 2*a1**4*r1**2 + 2*a2**4*r2**2 + 2*a3**4*r3**2 + 2*a4**4*r4**2
f(140) = 2*a1**5*r1**2 + 2*a2**5*r2**2 + 2*a3**5*r3**2 + 2*a4**5*r4**2
f(141) = 2*a1**6*r1**2 + 2*a2**6*r2**2 + 2*a3**6*r3**2 + 2*a4**6*r4**2
f(142) = 2*a1*r1**3 + 2*a2*r2**3 + 2*a3*r3**3 + 2*a4*r4**3
f(143) = 2*a1**2*r1**3 + 2*a2**2*r2**3 + 2*a3**2*r3**3 + 2*a4**2*r4**3
f(144) = 2*a1**3*r1**3 + 2*a2**3*r2**3 + 2*a3**3*r3**3 + 2*a4**3*r4**3
f(145) = 2*a1**4*r1**3 + 2*a2**4*r2**3 + 2*a3**4*r3**3 + 2*a4**4*r4**3
f(146) = 2*a1**5*r1**3 + 2*a2**5*r2**3 + 2*a3**5*r3**3 + 2*a4**5*r4**3
f(147) = 2*a1*r1**4 + 2*a2*r2**4 + 2*a3*r3**4 + 2*a4*r4**4
f(148) = 2*a1**2*r1**4 + 2*a2**2*r2**4 + 2*a3**2*r3**4 + 2*a4**2*r4**4
f(149) = 2*a1**3*r1**4 + 2*a2**3*r2**4 + 2*a3**3*r3**4 + 2*a4**3*r4**4
f(150) = 2*a1**4*r1**4 + 2*a2**4*r2**4 + 2*a3**4*r3**4 + 2*a4**4*r4**4
f(151) = 2*a1*r1**5 + 2*a2*r2**5 + 2*a3*r3**5 + 2*a4*r4**5
f(152) = 2*a1**2*r1**5 + 2*a2**2*r2**5 + 2*a3**2*r3**5 + 2*a4**2*r4**5
f(153) = 2*a1**3*r1**5 + 2*a2**3*r2**5 + 2*a3**3*r3**5 + 2*a4**3*r4**5
f(154) = 2*a1*r1**6 + 2*a2*r2**6 + 2*a3*r3**6 + 2*a4*r4**6
f(155) = 2*a1**2*r1**6 + 2*a2**2*r2**6 + 2*a3**2*r3**6 + 2*a4**2*r4**6
f(156) = 2*a1*r1**7 + 2*a2*r2**7 + 2*a3*r3**7 + 2*a4*r4**7
f(157) = 2*a1*r2 + 2*a2*r1 + 2*a3*r4 + 2*a4*r3
f(158) = 2*a1**2*r2 + 2*a2**2*r1 + 2*a3**2*r4 + 2*a4**2*r3
f(159) = 2*a1**3*r2 + 2*a2**3*r1 + 2*a3**3*r4 + 2*a4**3*r3
f(160) = 2*a1**4*r2 + 2*a2**4*r1 + 2*a3**4*r4 + 2*a4**4*r3
f(161) = 2*a1**5*r2 + 2*a2**5*r1 + 2*a3**5*r4 + 2*a4**5*r3
f(162) = 2*a1**6*r2 + 2*a2**6*r1 + 2*a3**6*r4 + 2*a4**6*r3
f(163) = 2*a1**7*r2 + 2*a2**7*r1 + 2*a3**7*r4 + 2*a4**7*r3
f(164) = 2*a1*r2**2 + 2*a2*r1**2 + 2*a3*r4**2 + 2*a4*r3**2
f(165) = 2*a1**2*r2**2 + 2*a2**2*r1**2 + 2*a3**2*r4**2 + 2*a4**2*r3**2
f(166) = 2*a1**3*r2**2 + 2*a2**3*r1**2 + 2*a3**3*r4**2 + 2*a4**3*r3**2
f(167) = 2*a1**4*r2**2 + 2*a2**4*r1**2 + 2*a3**4*r4**2 + 2*a4**4*r3**2
f(168) = 2*a1**5*r2**2 + 2*a2**5*r1**2 + 2*a3**5*r4**2 + 2*a4**5*r3**2
f(169) = 2*a1**6*r2**2 + 2*a2**6*r1**2 + 2*a3**6*r4**2 + 2*a4**6*r3**2
f(170) = 2*a1*r2**3 + 2*a2*r1**3 + 2*a3*r4**3 + 2*a4*r3**3
f(171) = 2*a1**2*r2**3 + 2*a2**2*r1**3 + 2*a3**2*r4**3 + 2*a4**2*r3**3
f(172) = 2*a1**3*r2**3 + 2*a2**3*r1**3 + 2*a3**3*r4**3 + 2*a4**3*r3**3
f(173) = 2*a1**4*r2**3 + 2*a2**4*r1**3 + 2*a3**4*r4**3 + 2*a4**4*r3**3
f(174) = 2*a1**5*r2**3 + 2*a2**5*r1**3 + 2*a3**5*r4**3 + 2*a4**5*r3**3
f(175) = 2*a1*r2**4 + 2*a2*r1**4 + 2*a3*r4**4 + 2*a4*r3**4
f(176) = 2*a1**2*r2**4 + 2*a2**2*r1**4 + 2*a3**2*r4**4 + 2*a4**2*r3**4
f(177) = 2*a1**3*r2**4 + 2*a2**3*r1**4 + 2*a3**3*r4**4 + 2*a4**3*r3**4
f(178) = 2*a1**4*r2**4 + 2*a2**4*r1**4 + 2*a3**4*r4**4 + 2*a4**4*r3**4
f(179) = 2*a1*r2**5 + 2*a2*r1**5 + 2*a3*r4**5 + 2*a4*r3**5
f(180) = 2*a1**2*r2**5 + 2*a2**2*r1**5 + 2*a3**2*r4**5 + 2*a4**2*r3**5
f(181) = 2*a1**3*r2**5 + 2*a2**3*r1**5 + 2*a3**3*r4**5 + 2*a4**3*r3**5
f(182) = 2*a1*r2**6 + 2*a2*r1**6 + 2*a3*r4**6 + 2*a4*r3**6
f(183) = 2*a1**2*r2**6 + 2*a2**2*r1**6 + 2*a3**2*r4**6 + 2*a4**2*r3**6
f(184) = 2*a1*r2**7 + 2*a2*r1**7 + 2*a3*r4**7 + 2*a4*r3**7
f(185) = 2*a1*r3 + 2*a2*r4 + 2*a3*r1 + 2*a4*r2
f(186) = 2*a1**2*r3 + 2*a2**2*r4 + 2*a3**2*r1 + 2*a4**2*r2
f(187) = 2*a1**3*r3 + 2*a2**3*r4 + 2*a3**3*r1 + 2*a4**3*r2
f(188) = 2*a1**4*r3 + 2*a2**4*r4 + 2*a3**4*r1 + 2*a4**4*r2
f(189) = 2*a1**5*r3 + 2*a2**5*r4 + 2*a3**5*r1 + 2*a4**5*r2
f(190) = 2*a1**6*r3 + 2*a2**6*r4 + 2*a3**6*r1 + 2*a4**6*r2
f(191) = 2*a1**7*r3 + 2*a2**7*r4 + 2*a3**7*r1 + 2*a4**7*r2
f(192) = 2*a1*r3**2 + 2*a2*r4**2 + 2*a3*r1**2 + 2*a4*r2**2
f(193) = 2*a1**2*r3**2 + 2*a2**2*r4**2 + 2*a3**2*r1**2 + 2*a4**2*r2**2
f(194) = 2*a1**3*r3**2 + 2*a2**3*r4**2 + 2*a3**3*r1**2 + 2*a4**3*r2**2
f(195) = 2*a1**4*r3**2 + 2*a2**4*r4**2 + 2*a3**4*r1**2 + 2*a4**4*r2**2
f(196) = 2*a1**5*r3**2 + 2*a2**5*r4**2 + 2*a3**5*r1**2 + 2*a4**5*r2**2
f(197) = 2*a1**6*r3**2 + 2*a2**6*r4**2 + 2*a3**6*r1**2 + 2*a4**6*r2**2
f(198) = 2*a1*r3**3 + 2*a2*r4**3 + 2*a3*r1**3 + 2*a4*r2**3
f(199) = 2*a1**2*r3**3 + 2*a2**2*r4**3 + 2*a3**2*r1**3 + 2*a4**2*r2**3
f(200) = 2*a1**3*r3**3 + 2*a2**3*r4**3 + 2*a3**3*r1**3 + 2*a4**3*r2**3
f(201) = 2*a1**4*r3**3 + 2*a2**4*r4**3 + 2*a3**4*r1**3 + 2*a4**4*r2**3
f(202) = 2*a1**5*r3**3 + 2*a2**5*r4**3 + 2*a3**5*r1**3 + 2*a4**5*r2**3
f(203) = 2*a1*r3**4 + 2*a2*r4**4 + 2*a3*r1**4 + 2*a4*r2**4
f(204) = 2*a1**2*r3**4 + 2*a2**2*r4**4 + 2*a3**2*r1**4 + 2*a4**2*r2**4
f(205) = 2*a1**3*r3**4 + 2*a2**3*r4**4 + 2*a3**3*r1**4 + 2*a4**3*r2**4
f(206) = 2*a1**4*r3**4 + 2*a2**4*r4**4 + 2*a3**4*r1**4 + 2*a4**4*r2**4
f(207) = 2*a1*r3**5 + 2*a2*r4**5 + 2*a3*r1**5 + 2*a4*r2**5
f(208) = 2*a1**2*r3**5 + 2*a2**2*r4**5 + 2*a3**2*r1**5 + 2*a4**2*r2**5
f(209) = 2*a1**3*r3**5 + 2*a2**3*r4**5 + 2*a3**3*r1**5 + 2*a4**3*r2**5
f(210) = 2*a1*r3**6 + 2*a2*r4**6 + 2*a3*r1**6 + 2*a4*r2**6
f(211) = 2*a1**2*r3**6 + 2*a2**2*r4**6 + 2*a3**2*r1**6 + 2*a4**2*r2**6
f(212) = 2*a1*r3**7 + 2*a2*r4**7 + 2*a3*r1**7 + 2*a4*r2**7
f(213) = 2*a1*r4 + 2*a2*r3 + 2*a3*r2 + 2*a4*r1
f(214) = 2*a1**2*r4 + 2*a2**2*r3 + 2*a3**2*r2 + 2*a4**2*r1
f(215) = 2*a1**3*r4 + 2*a2**3*r3 + 2*a3**3*r2 + 2*a4**3*r1
f(216) = 2*a1**4*r4 + 2*a2**4*r3 + 2*a3**4*r2 + 2*a4**4*r1
f(217) = 2*a1**5*r4 + 2*a2**5*r3 + 2*a3**5*r2 + 2*a4**5*r1
f(218) = 2*a1**6*r4 + 2*a2**6*r3 + 2*a3**6*r2 + 2*a4**6*r1
f(219) = 2*a1**7*r4 + 2*a2**7*r3 + 2*a3**7*r2 + 2*a4**7*r1
f(220) = 2*a1*r4**2 + 2*a2*r3**2 + 2*a3*r2**2 + 2*a4*r1**2
f(221) = 2*a1**2*r4**2 + 2*a2**2*r3**2 + 2*a3**2*r2**2 + 2*a4**2*r1**2
f(222) = 2*a1**3*r4**2 + 2*a2**3*r3**2 + 2*a3**3*r2**2 + 2*a4**3*r1**2
f(223) = 2*a1**4*r4**2 + 2*a2**4*r3**2 + 2*a3**4*r2**2 + 2*a4**4*r1**2
f(224) = 2*a1**5*r4**2 + 2*a2**5*r3**2 + 2*a3**5*r2**2 + 2*a4**5*r1**2
f(225) = 2*a1**6*r4**2 + 2*a2**6*r3**2 + 2*a3**6*r2**2 + 2*a4**6*r1**2
f(226) = 2*a1*r4**3 + 2*a2*r3**3 + 2*a3*r2**3 + 2*a4*r1**3
f(227) = 2*a1**2*r4**3 + 2*a2**2*r3**3 + 2*a3**2*r2**3 + 2*a4**2*r1**3
f(228) = 2*a1**3*r4**3 + 2*a2**3*r3**3 + 2*a3**3*r2**3 + 2*a4**3*r1**3
f(229) = 2*a1**4*r4**3 + 2*a2**4*r3**3 + 2*a3**4*r2**3 + 2*a4**4*r1**3
f(230) = 2*a1**5*r4**3 + 2*a2**5*r3**3 + 2*a3**5*r2**3 + 2*a4**5*r1**3
f(231) = 2*a1*r4**4 + 2*a2*r3**4 + 2*a3*r2**4 + 2*a4*r1**4
f(232) = 2*a1**2*r4**4 + 2*a2**2*r3**4 + 2*a3**2*r2**4 + 2*a4**2*r1**4
f(233) = 2*a1**3*r4**4 + 2*a2**3*r3**4 + 2*a3**3*r2**4 + 2*a4**3*r1**4
f(234) = 2*a1**4*r4**4 + 2*a2**4*r3**4 + 2*a3**4*r2**4 + 2*a4**4*r1**4
f(235) = 2*a1*r4**5 + 2*a2*r3**5 + 2*a3*r2**5 + 2*a4*r1**5
f(236) = 2*a1**2*r4**5 + 2*a2**2*r3**5 + 2*a3**2*r2**5 + 2*a4**2*r1**5
f(237) = 2*a1**3*r4**5 + 2*a2**3*r3**5 + 2*a3**3*r2**5 + 2*a4**3*r1**5
f(238) = 2*a1*r4**6 + 2*a2*r3**6 + 2*a3*r2**6 + 2*a4*r1**6
f(239) = 2*a1**2*r4**6 + 2*a2**2*r3**6 + 2*a3**2*r2**6 + 2*a4**2*r1**6
f(240) = 2*a1*r4**7 + 2*a2*r3**7 + 2*a3*r2**7 + 2*a4*r1**7
f(241) = 2*b1**2*r1 + 2*b1**2*r2 + 2*b2**2*r3 + 2*b2**2*r4
f(242) = 2*b1**4*r1 + 2*b1**4*r2 + 2*b2**4*r3 + 2*b2**4*r4
f(243) = 2*b1**6*r1 + 2*b1**6*r2 + 2*b2**6*r3 + 2*b2**6*r4
f(244) = 2*b1**2*r1**2 + 2*b1**2*r2**2 + 2*b2**2*r3**2 + 2*b2**2*r4**2
f(245) = 2*b1**4*r1**2 + 2*b1**4*r2**2 + 2*b2**4*r3**2 + 2*b2**4*r4**2
f(246) = 2*b1**6*r1**2 + 2*b1**6*r2**2 + 2*b2**6*r3**2 + 2*b2**6*r4**2
f(247) = 2*b1**2*r1**3 + 2*b1**2*r2**3 + 2*b2**2*r3**3 + 2*b2**2*r4**3
f(248) = 2*b1**4*r1**3 + 2*b1**4*r2**3 + 2*b2**4*r3**3 + 2*b2**4*r4**3
f(249) = 2*b1**2*r1**4 + 2*b1**2*r2**4 + 2*b2**2*r3**4 + 2*b2**2*r4**4
f(250) = 2*b1**4*r1**4 + 2*b1**4*r2**4 + 2*b2**4*r3**4 + 2*b2**4*r4**4
f(251) = 2*b1**2*r1**5 + 2*b1**2*r2**5 + 2*b2**2*r3**5 + 2*b2**2*r4**5
f(252) = 2*b1**2*r1**6 + 2*b1**2*r2**6 + 2*b2**2*r3**6 + 2*b2**2*r4**6
f(253) = 2*b1**2*r3 + 2*b1**2*r4 + 2*b2**2*r1 + 2*b2**2*r2
f(254) = 2*b1**4*r3 + 2*b1**4*r4 + 2*b2**4*r1 + 2*b2**4*r2
f(255) = 2*b1**6*r3 + 2*b1**6*r4 + 2*b2**6*r1 + 2*b2**6*r2
f(256) = 2*b1**2*r3**2 + 2*b1**2*r4**2 + 2*b2**2*r1**2 + 2*b2**2*r2**2
f(257) = 2*b1**4*r3**2 + 2*b1**4*r4**2 + 2*b2**4*r1**2 + 2*b2**4*r2**2
f(258) = 2*b1**6*r3**2 + 2*b1**6*r4**2 + 2*b2**6*r1**2 + 2*b2**6*r2**2
f(259) = 2*b1**2*r3**3 + 2*b1**2*r4**3 + 2*b2**2*r1**3 + 2*b2**2*r2**3
f(260) = 2*b1**4*r3**3 + 2*b1**4*r4**3 + 2*b2**4*r1**3 + 2*b2**4*r2**3
f(261) = 2*b1**2*r3**4 + 2*b1**2*r4**4 + 2*b2**2*r1**4 + 2*b2**2*r2**4
f(262) = 2*b1**4*r3**4 + 2*b1**4*r4**4 + 2*b2**4*r1**4 + 2*b2**4*r2**4
f(263) = 2*b1**2*r3**5 + 2*b1**2*r4**5 + 2*b2**2*r1**5 + 2*b2**2*r2**5
f(264) = 2*b1**2*r3**6 + 2*b1**2*r4**6 + 2*b2**2*r1**6 + 2*b2**2*r2**6
f(265) = 2*dtau**2*(r1 + r2 + r3 + r4)
f(266) = 2*dtau**4*(r1 + r2 + r3 + r4)
f(267) = 2*dtau**6*(r1 + r2 + r3 + r4)
f(268) = 2*dtau**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(269) = 2*dtau**4*(r1**2 + r2**2 + r3**2 + r4**2)
f(270) = 2*dtau**6*(r1**2 + r2**2 + r3**2 + r4**2)
f(271) = 2*dtau**2*(r1**3 + r2**3 + r3**3 + r4**3)
f(272) = 2*dtau**4*(r1**3 + r2**3 + r3**3 + r4**3)
f(273) = 2*dtau**2*(r1**4 + r2**4 + r3**4 + r4**4)
f(274) = 2*dtau**4*(r1**4 + r2**4 + r3**4 + r4**4)
f(275) = 2*dtau**2*(r1**5 + r2**5 + r3**5 + r4**5)
f(276) = 2*dtau**2*(r1**6 + r2**6 + r3**6 + r4**6)
f(277) = 4*a1*a2 + 4*a3*a4
f(278) = 2*a1**2*a2 + 2*a1*a2**2 + 2*a3**2*a4 + 2*a3*a4**2
f(279) = 2*a1**3*a2 + 2*a1*a2**3 + 2*a3**3*a4 + 2*a3*a4**3
f(280) = 2*a1**4*a2 + 2*a1*a2**4 + 2*a3**4*a4 + 2*a3*a4**4
f(281) = 2*a1**5*a2 + 2*a1*a2**5 + 2*a3**5*a4 + 2*a3*a4**5
f(282) = 2*a1**6*a2 + 2*a1*a2**6 + 2*a3**6*a4 + 2*a3*a4**6
f(283) = 2*a1**7*a2 + 2*a1*a2**7 + 2*a3**7*a4 + 2*a3*a4**7
f(284) = 4*a1**2*a2**2 + 4*a3**2*a4**2
f(285) = 2*a1**3*a2**2 + 2*a1**2*a2**3 + 2*a3**3*a4**2 + 2*a3**2*a4**3
f(286) = 2*a1**4*a2**2 + 2*a1**2*a2**4 + 2*a3**4*a4**2 + 2*a3**2*a4**4
f(287) = 2*a1**5*a2**2 + 2*a1**2*a2**5 + 2*a3**5*a4**2 + 2*a3**2*a4**5
f(288) = 2*a1**6*a2**2 + 2*a1**2*a2**6 + 2*a3**6*a4**2 + 2*a3**2*a4**6
f(289) = 4*a1**3*a2**3 + 4*a3**3*a4**3
f(290) = 2*a1**4*a2**3 + 2*a1**3*a2**4 + 2*a3**4*a4**3 + 2*a3**3*a4**4
f(291) = 2*a1**5*a2**3 + 2*a1**3*a2**5 + 2*a3**5*a4**3 + 2*a3**3*a4**5
f(292) = 4*a1**4*a2**4 + 4*a3**4*a4**4
f(293) = 4*a1*a3 + 4*a2*a4
f(294) = 2*a1**2*a3 + 2*a1*a3**2 + 2*a2**2*a4 + 2*a2*a4**2
f(295) = 2*a1**3*a3 + 2*a1*a3**3 + 2*a2**3*a4 + 2*a2*a4**3
f(296) = 2*a1**4*a3 + 2*a1*a3**4 + 2*a2**4*a4 + 2*a2*a4**4
f(297) = 2*a1**5*a3 + 2*a1*a3**5 + 2*a2**5*a4 + 2*a2*a4**5
f(298) = 2*a1**6*a3 + 2*a1*a3**6 + 2*a2**6*a4 + 2*a2*a4**6
f(299) = 2*a1**7*a3 + 2*a1*a3**7 + 2*a2**7*a4 + 2*a2*a4**7
f(300) = 4*a1**2*a3**2 + 4*a2**2*a4**2
f(301) = 2*a1**3*a3**2 + 2*a1**2*a3**3 + 2*a2**3*a4**2 + 2*a2**2*a4**3
f(302) = 2*a1**4*a3**2 + 2*a1**2*a3**4 + 2*a2**4*a4**2 + 2*a2**2*a4**4
f(303) = 2*a1**5*a3**2 + 2*a1**2*a3**5 + 2*a2**5*a4**2 + 2*a2**2*a4**5
f(304) = 2*a1**6*a3**2 + 2*a1**2*a3**6 + 2*a2**6*a4**2 + 2*a2**2*a4**6
f(305) = 4*a1**3*a3**3 + 4*a2**3*a4**3
f(306) = 2*a1**4*a3**3 + 2*a1**3*a3**4 + 2*a2**4*a4**3 + 2*a2**3*a4**4
f(307) = 2*a1**5*a3**3 + 2*a1**3*a3**5 + 2*a2**5*a4**3 + 2*a2**3*a4**5
f(308) = 4*a1**4*a3**4 + 4*a2**4*a4**4
f(309) = 4*a1*a4 + 4*a2*a3
f(310) = 2*a1**2*a4 + 2*a1*a4**2 + 2*a2**2*a3 + 2*a2*a3**2
f(311) = 2*a1**3*a4 + 2*a1*a4**3 + 2*a2**3*a3 + 2*a2*a3**3
f(312) = 2*a1**4*a4 + 2*a1*a4**4 + 2*a2**4*a3 + 2*a2*a3**4
f(313) = 2*a1**5*a4 + 2*a1*a4**5 + 2*a2**5*a3 + 2*a2*a3**5
f(314) = 2*a1**6*a4 + 2*a1*a4**6 + 2*a2**6*a3 + 2*a2*a3**6
f(315) = 2*a1**7*a4 + 2*a1*a4**7 + 2*a2**7*a3 + 2*a2*a3**7
f(316) = 4*a1**2*a4**2 + 4*a2**2*a3**2
f(317) = 2*a1**3*a4**2 + 2*a1**2*a4**3 + 2*a2**3*a3**2 + 2*a2**2*a3**3
f(318) = 2*a1**4*a4**2 + 2*a1**2*a4**4 + 2*a2**4*a3**2 + 2*a2**2*a3**4
f(319) = 2*a1**5*a4**2 + 2*a1**2*a4**5 + 2*a2**5*a3**2 + 2*a2**2*a3**5
f(320) = 2*a1**6*a4**2 + 2*a1**2*a4**6 + 2*a2**6*a3**2 + 2*a2**2*a3**6
f(321) = 4*a1**3*a4**3 + 4*a2**3*a3**3
f(322) = 2*a1**4*a4**3 + 2*a1**3*a4**4 + 2*a2**4*a3**3 + 2*a2**3*a3**4
f(323) = 2*a1**5*a4**3 + 2*a1**3*a4**5 + 2*a2**5*a3**3 + 2*a2**3*a3**5
f(324) = 4*a1**4*a4**4 + 4*a2**4*a3**4
f(325) = 2*a1*b1**2 + 2*a2*b1**2 + 2*a3*b2**2 + 2*a4*b2**2
f(326) = 2*a1*b1**4 + 2*a2*b1**4 + 2*a3*b2**4 + 2*a4*b2**4
f(327) = 2*a1*b1**6 + 2*a2*b1**6 + 2*a3*b2**6 + 2*a4*b2**6
f(328) = 2*a1**2*b1**2 + 2*a2**2*b1**2 + 2*a3**2*b2**2 + 2*a4**2*b2**2
f(329) = 2*a1**2*b1**4 + 2*a2**2*b1**4 + 2*a3**2*b2**4 + 2*a4**2*b2**4
f(330) = 2*a1**2*b1**6 + 2*a2**2*b1**6 + 2*a3**2*b2**6 + 2*a4**2*b2**6
f(331) = 2*a1**3*b1**2 + 2*a2**3*b1**2 + 2*a3**3*b2**2 + 2*a4**3*b2**2
f(332) = 2*a1**3*b1**4 + 2*a2**3*b1**4 + 2*a3**3*b2**4 + 2*a4**3*b2**4
f(333) = 2*a1**4*b1**2 + 2*a2**4*b1**2 + 2*a3**4*b2**2 + 2*a4**4*b2**2
f(334) = 2*a1**4*b1**4 + 2*a2**4*b1**4 + 2*a3**4*b2**4 + 2*a4**4*b2**4
f(335) = 2*a1**5*b1**2 + 2*a2**5*b1**2 + 2*a3**5*b2**2 + 2*a4**5*b2**2
f(336) = 2*a1**6*b1**2 + 2*a2**6*b1**2 + 2*a3**6*b2**2 + 2*a4**6*b2**2
f(337) = 2*a1*b2**2 + 2*a2*b2**2 + 2*a3*b1**2 + 2*a4*b1**2
f(338) = 2*a1*b2**4 + 2*a2*b2**4 + 2*a3*b1**4 + 2*a4*b1**4
f(339) = 2*a1*b2**6 + 2*a2*b2**6 + 2*a3*b1**6 + 2*a4*b1**6
f(340) = 2*a1**2*b2**2 + 2*a2**2*b2**2 + 2*a3**2*b1**2 + 2*a4**2*b1**2
f(341) = 2*a1**2*b2**4 + 2*a2**2*b2**4 + 2*a3**2*b1**4 + 2*a4**2*b1**4
f(342) = 2*a1**2*b2**6 + 2*a2**2*b2**6 + 2*a3**2*b1**6 + 2*a4**2*b1**6
f(343) = 2*a1**3*b2**2 + 2*a2**3*b2**2 + 2*a3**3*b1**2 + 2*a4**3*b1**2
f(344) = 2*a1**3*b2**4 + 2*a2**3*b2**4 + 2*a3**3*b1**4 + 2*a4**3*b1**4
f(345) = 2*a1**4*b2**2 + 2*a2**4*b2**2 + 2*a3**4*b1**2 + 2*a4**4*b1**2
f(346) = 2*a1**4*b2**4 + 2*a2**4*b2**4 + 2*a3**4*b1**4 + 2*a4**4*b1**4
f(347) = 2*a1**5*b2**2 + 2*a2**5*b2**2 + 2*a3**5*b1**2 + 2*a4**5*b1**2
f(348) = 2*a1**6*b2**2 + 2*a2**6*b2**2 + 2*a3**6*b1**2 + 2*a4**6*b1**2
f(349) = 2*dtau**2*(a1 + a2 + a3 + a4)
f(350) = 2*dtau**4*(a1 + a2 + a3 + a4)
f(351) = 2*dtau**6*(a1 + a2 + a3 + a4)
f(352) = 2*dtau**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(353) = 2*dtau**4*(a1**2 + a2**2 + a3**2 + a4**2)
f(354) = 2*dtau**6*(a1**2 + a2**2 + a3**2 + a4**2)
f(355) = 2*dtau**2*(a1**3 + a2**3 + a3**3 + a4**3)
f(356) = 2*dtau**4*(a1**3 + a2**3 + a3**3 + a4**3)
f(357) = 2*dtau**2*(a1**4 + a2**4 + a3**4 + a4**4)
f(358) = 2*dtau**4*(a1**4 + a2**4 + a3**4 + a4**4)
f(359) = 2*dtau**2*(a1**5 + a2**5 + a3**5 + a4**5)
f(360) = 2*dtau**2*(a1**6 + a2**6 + a3**6 + a4**6)
f(361) = 8*b1*b2
f(362) = 4*b1*b2*(b1**2 + b2**2)
f(363) = 4*b1*b2*(b1**4 + b2**4)
f(364) = 4*b1*b2*(b1**6 + b2**6)
f(365) = 8*b1**2*b2**2
f(366) = 4*b1**2*b2**2*(b1**2 + b2**2)
f(367) = 4*b1**2*b2**2*(b1**4 + b2**4)
f(368) = 8*b1**3*b2**3
f(369) = 4*b1**3*b2**3*(b1**2 + b2**2)
f(370) = 8*b1**4*b2**4
f(371) = 4*dtau**2*(b1**2 + b2**2)
f(372) = 4*dtau**4*(b1**2 + b2**2)
f(373) = 4*dtau**6*(b1**2 + b2**2)
f(374) = 4*dtau**2*(b1**4 + b2**4)
f(375) = 4*dtau**4*(b1**4 + b2**4)
f(376) = 4*dtau**2*(b1**6 + b2**6)
v = sum(f*params)
end function c2h4_poten_n2_d8_ADF


!###############################################################################


! D2h-symmetry-adapted 3-mode 6-order expansion for C2H4 molecule

function c2h4_poten_n3_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(716)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(716)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 4*r0*(r1*r2 + r3*r4)
f(2) = 2*r0*(r1**2*r2 + r1*r2**2 + r3**2*r4 + r3*r4**2)
f(3) = 2*r0*(r1**3*r2 + r1*r2**3 + r3**3*r4 + r3*r4**3)
f(4) = 2*r0*(r1**4*r2 + r1*r2**4 + r3**4*r4 + r3*r4**4)
f(5) = 4*r0*(r1**2*r2**2 + r3**2*r4**2)
f(6) = 2*r0*(r1**3*r2**2 + r1**2*r2**3 + r3**3*r4**2 + r3**2*r4**3)
f(7) = 4*r0**2*(r1*r2 + r3*r4)
f(8) = 2*r0**2*(r1**2*r2 + r1*r2**2 + r3**2*r4 + r3*r4**2)
f(9) = 2*r0**2*(r1**3*r2 + r1*r2**3 + r3**3*r4 + r3*r4**3)
f(10) = 4*r0**2*(r1**2*r2**2 + r3**2*r4**2)
f(11) = 4*r0**3*(r1*r2 + r3*r4)
f(12) = 2*r0**3*(r1**2*r2 + r1*r2**2 + r3**2*r4 + r3*r4**2)
f(13) = 4*r0**4*(r1*r2 + r3*r4)
f(14) = 4*r0*(r1*r3 + r2*r4)
f(15) = 2*r0*(r1**2*r3 + r1*r3**2 + r2**2*r4 + r2*r4**2)
f(16) = 2*r0*(r1**3*r3 + r1*r3**3 + r2**3*r4 + r2*r4**3)
f(17) = 2*r0*(r1**4*r3 + r1*r3**4 + r2**4*r4 + r2*r4**4)
f(18) = 4*r0*(r1**2*r3**2 + r2**2*r4**2)
f(19) = 2*r0*(r1**3*r3**2 + r1**2*r3**3 + r2**3*r4**2 + r2**2*r4**3)
f(20) = 4*r0**2*(r1*r3 + r2*r4)
f(21) = 2*r0**2*(r1**2*r3 + r1*r3**2 + r2**2*r4 + r2*r4**2)
f(22) = 2*r0**2*(r1**3*r3 + r1*r3**3 + r2**3*r4 + r2*r4**3)
f(23) = 4*r0**2*(r1**2*r3**2 + r2**2*r4**2)
f(24) = 4*r0**3*(r1*r3 + r2*r4)
f(25) = 2*r0**3*(r1**2*r3 + r1*r3**2 + r2**2*r4 + r2*r4**2)
f(26) = 4*r0**4*(r1*r3 + r2*r4)
f(27) = 4*r0*(r1*r4 + r2*r3)
f(28) = 2*r0*(r1**2*r4 + r1*r4**2 + r2**2*r3 + r2*r3**2)
f(29) = 2*r0*(r1**3*r4 + r1*r4**3 + r2**3*r3 + r2*r3**3)
f(30) = 2*r0*(r1**4*r4 + r1*r4**4 + r2**4*r3 + r2*r3**4)
f(31) = 4*r0*(r1**2*r4**2 + r2**2*r3**2)
f(32) = 2*r0*(r1**3*r4**2 + r1**2*r4**3 + r2**3*r3**2 + r2**2*r3**3)
f(33) = 4*r0**2*(r1*r4 + r2*r3)
f(34) = 2*r0**2*(r1**2*r4 + r1*r4**2 + r2**2*r3 + r2*r3**2)
f(35) = 2*r0**2*(r1**3*r4 + r1*r4**3 + r2**3*r3 + r2*r3**3)
f(36) = 4*r0**2*(r1**2*r4**2 + r2**2*r3**2)
f(37) = 4*r0**3*(r1*r4 + r2*r3)
f(38) = 2*r0**3*(r1**2*r4 + r1*r4**2 + r2**2*r3 + r2*r3**2)
f(39) = 4*r0**4*(r1*r4 + r2*r3)
f(40) = 2*r0*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(41) = 2*r0*(a1**2*r1 + a2**2*r2 + a3**2*r3 + a4**2*r4)
f(42) = 2*r0*(a1**3*r1 + a2**3*r2 + a3**3*r3 + a4**3*r4)
f(43) = 2*r0*(a1**4*r1 + a2**4*r2 + a3**4*r3 + a4**4*r4)
f(44) = 2*r0*(a1*r1**2 + a2*r2**2 + a3*r3**2 + a4*r4**2)
f(45) = 2*r0*(a1**2*r1**2 + a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2)
f(46) = 2*r0*(a1**3*r1**2 + a2**3*r2**2 + a3**3*r3**2 + a4**3*r4**2)
f(47) = 2*r0*(a1*r1**3 + a2*r2**3 + a3*r3**3 + a4*r4**3)
f(48) = 2*r0*(a1**2*r1**3 + a2**2*r2**3 + a3**2*r3**3 + a4**2*r4**3)
f(49) = 2*r0*(a1*r1**4 + a2*r2**4 + a3*r3**4 + a4*r4**4)
f(50) = 2*r0**2*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(51) = 2*r0**2*(a1**2*r1 + a2**2*r2 + a3**2*r3 + a4**2*r4)
f(52) = 2*r0**2*(a1**3*r1 + a2**3*r2 + a3**3*r3 + a4**3*r4)
f(53) = 2*r0**2*(a1*r1**2 + a2*r2**2 + a3*r3**2 + a4*r4**2)
f(54) = 2*r0**2*(a1**2*r1**2 + a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2)
f(55) = 2*r0**2*(a1*r1**3 + a2*r2**3 + a3*r3**3 + a4*r4**3)
f(56) = 2*r0**3*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(57) = 2*r0**3*(a1**2*r1 + a2**2*r2 + a3**2*r3 + a4**2*r4)
f(58) = 2*r0**3*(a1*r1**2 + a2*r2**2 + a3*r3**2 + a4*r4**2)
f(59) = 2*r0**4*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(60) = 2*r0*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(61) = 2*r0*(a1**2*r2 + a2**2*r1 + a3**2*r4 + a4**2*r3)
f(62) = 2*r0*(a1**3*r2 + a2**3*r1 + a3**3*r4 + a4**3*r3)
f(63) = 2*r0*(a1**4*r2 + a2**4*r1 + a3**4*r4 + a4**4*r3)
f(64) = 2*r0*(a1*r2**2 + a2*r1**2 + a3*r4**2 + a4*r3**2)
f(65) = 2*r0*(a1**2*r2**2 + a2**2*r1**2 + a3**2*r4**2 + a4**2*r3**2)
f(66) = 2*r0*(a1**3*r2**2 + a2**3*r1**2 + a3**3*r4**2 + a4**3*r3**2)
f(67) = 2*r0*(a1*r2**3 + a2*r1**3 + a3*r4**3 + a4*r3**3)
f(68) = 2*r0*(a1**2*r2**3 + a2**2*r1**3 + a3**2*r4**3 + a4**2*r3**3)
f(69) = 2*r0*(a1*r2**4 + a2*r1**4 + a3*r4**4 + a4*r3**4)
f(70) = 2*r0**2*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(71) = 2*r0**2*(a1**2*r2 + a2**2*r1 + a3**2*r4 + a4**2*r3)
f(72) = 2*r0**2*(a1**3*r2 + a2**3*r1 + a3**3*r4 + a4**3*r3)
f(73) = 2*r0**2*(a1*r2**2 + a2*r1**2 + a3*r4**2 + a4*r3**2)
f(74) = 2*r0**2*(a1**2*r2**2 + a2**2*r1**2 + a3**2*r4**2 + a4**2*r3**2)
f(75) = 2*r0**2*(a1*r2**3 + a2*r1**3 + a3*r4**3 + a4*r3**3)
f(76) = 2*r0**3*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(77) = 2*r0**3*(a1**2*r2 + a2**2*r1 + a3**2*r4 + a4**2*r3)
f(78) = 2*r0**3*(a1*r2**2 + a2*r1**2 + a3*r4**2 + a4*r3**2)
f(79) = 2*r0**4*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(80) = 2*r0*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(81) = 2*r0*(a1**2*r3 + a2**2*r4 + a3**2*r1 + a4**2*r2)
f(82) = 2*r0*(a1**3*r3 + a2**3*r4 + a3**3*r1 + a4**3*r2)
f(83) = 2*r0*(a1**4*r3 + a2**4*r4 + a3**4*r1 + a4**4*r2)
f(84) = 2*r0*(a1*r3**2 + a2*r4**2 + a3*r1**2 + a4*r2**2)
f(85) = 2*r0*(a1**2*r3**2 + a2**2*r4**2 + a3**2*r1**2 + a4**2*r2**2)
f(86) = 2*r0*(a1**3*r3**2 + a2**3*r4**2 + a3**3*r1**2 + a4**3*r2**2)
f(87) = 2*r0*(a1*r3**3 + a2*r4**3 + a3*r1**3 + a4*r2**3)
f(88) = 2*r0*(a1**2*r3**3 + a2**2*r4**3 + a3**2*r1**3 + a4**2*r2**3)
f(89) = 2*r0*(a1*r3**4 + a2*r4**4 + a3*r1**4 + a4*r2**4)
f(90) = 2*r0**2*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(91) = 2*r0**2*(a1**2*r3 + a2**2*r4 + a3**2*r1 + a4**2*r2)
f(92) = 2*r0**2*(a1**3*r3 + a2**3*r4 + a3**3*r1 + a4**3*r2)
f(93) = 2*r0**2*(a1*r3**2 + a2*r4**2 + a3*r1**2 + a4*r2**2)
f(94) = 2*r0**2*(a1**2*r3**2 + a2**2*r4**2 + a3**2*r1**2 + a4**2*r2**2)
f(95) = 2*r0**2*(a1*r3**3 + a2*r4**3 + a3*r1**3 + a4*r2**3)
f(96) = 2*r0**3*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(97) = 2*r0**3*(a1**2*r3 + a2**2*r4 + a3**2*r1 + a4**2*r2)
f(98) = 2*r0**3*(a1*r3**2 + a2*r4**2 + a3*r1**2 + a4*r2**2)
f(99) = 2*r0**4*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(100) = 2*r0*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(101) = 2*r0*(a1**2*r4 + a2**2*r3 + a3**2*r2 + a4**2*r1)
f(102) = 2*r0*(a1**3*r4 + a2**3*r3 + a3**3*r2 + a4**3*r1)
f(103) = 2*r0*(a1**4*r4 + a2**4*r3 + a3**4*r2 + a4**4*r1)
f(104) = 2*r0*(a1*r4**2 + a2*r3**2 + a3*r2**2 + a4*r1**2)
f(105) = 2*r0*(a1**2*r4**2 + a2**2*r3**2 + a3**2*r2**2 + a4**2*r1**2)
f(106) = 2*r0*(a1**3*r4**2 + a2**3*r3**2 + a3**3*r2**2 + a4**3*r1**2)
f(107) = 2*r0*(a1*r4**3 + a2*r3**3 + a3*r2**3 + a4*r1**3)
f(108) = 2*r0*(a1**2*r4**3 + a2**2*r3**3 + a3**2*r2**3 + a4**2*r1**3)
f(109) = 2*r0*(a1*r4**4 + a2*r3**4 + a3*r2**4 + a4*r1**4)
f(110) = 2*r0**2*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(111) = 2*r0**2*(a1**2*r4 + a2**2*r3 + a3**2*r2 + a4**2*r1)
f(112) = 2*r0**2*(a1**3*r4 + a2**3*r3 + a3**3*r2 + a4**3*r1)
f(113) = 2*r0**2*(a1*r4**2 + a2*r3**2 + a3*r2**2 + a4*r1**2)
f(114) = 2*r0**2*(a1**2*r4**2 + a2**2*r3**2 + a3**2*r2**2 + a4**2*r1**2)
f(115) = 2*r0**2*(a1*r4**3 + a2*r3**3 + a3*r2**3 + a4*r1**3)
f(116) = 2*r0**3*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(117) = 2*r0**3*(a1**2*r4 + a2**2*r3 + a3**2*r2 + a4**2*r1)
f(118) = 2*r0**3*(a1*r4**2 + a2*r3**2 + a3*r2**2 + a4*r1**2)
f(119) = 2*r0**4*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(120) = 2*r0*(b1**2*r1 + b1**2*r2 + b2**2*r3 + b2**2*r4)
f(121) = 2*r0*(b1**4*r1 + b1**4*r2 + b2**4*r3 + b2**4*r4)
f(122) = 2*r0*(b1**2*r1**2 + b1**2*r2**2 + b2**2*r3**2 + b2**2*r4**2)
f(123) = 2*r0*(b1**2*r1**3 + b1**2*r2**3 + b2**2*r3**3 + b2**2*r4**3)
f(124) = 2*r0**2*(b1**2*r1 + b1**2*r2 + b2**2*r3 + b2**2*r4)
f(125) = 2*r0**2*(b1**2*r1**2 + b1**2*r2**2 + b2**2*r3**2 + b2**2*r4**2)
f(126) = 2*r0**3*(b1**2*r1 + b1**2*r2 + b2**2*r3 + b2**2*r4)
f(127) = 2*r0*(b1**2*r3 + b1**2*r4 + b2**2*r1 + b2**2*r2)
f(128) = 2*r0*(b1**4*r3 + b1**4*r4 + b2**4*r1 + b2**4*r2)
f(129) = 2*r0*(b1**2*r3**2 + b1**2*r4**2 + b2**2*r1**2 + b2**2*r2**2)
f(130) = 2*r0*(b1**2*r3**3 + b1**2*r4**3 + b2**2*r1**3 + b2**2*r2**3)
f(131) = 2*r0**2*(b1**2*r3 + b1**2*r4 + b2**2*r1 + b2**2*r2)
f(132) = 2*r0**2*(b1**2*r3**2 + b1**2*r4**2 + b2**2*r1**2 + b2**2*r2**2)
f(133) = 2*r0**3*(b1**2*r3 + b1**2*r4 + b2**2*r1 + b2**2*r2)
f(134) = 2*dtau**2*r0*(r1 + r2 + r3 + r4)
f(135) = 2*dtau**4*r0*(r1 + r2 + r3 + r4)
f(136) = 2*dtau**2*r0*(r1**2 + r2**2 + r3**2 + r4**2)
f(137) = 2*dtau**2*r0*(r1**3 + r2**3 + r3**3 + r4**3)
f(138) = 2*dtau**2*r0**2*(r1 + r2 + r3 + r4)
f(139) = 2*dtau**2*r0**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(140) = 2*dtau**2*r0**3*(r1 + r2 + r3 + r4)
f(141) = 4*r0*(a1*a2 + a3*a4)
f(142) = 2*r0*(a1**2*a2 + a1*a2**2 + a3**2*a4 + a3*a4**2)
f(143) = 2*r0*(a1**3*a2 + a1*a2**3 + a3**3*a4 + a3*a4**3)
f(144) = 2*r0*(a1**4*a2 + a1*a2**4 + a3**4*a4 + a3*a4**4)
f(145) = 4*r0*(a1**2*a2**2 + a3**2*a4**2)
f(146) = 2*r0*(a1**3*a2**2 + a1**2*a2**3 + a3**3*a4**2 + a3**2*a4**3)
f(147) = 4*r0**2*(a1*a2 + a3*a4)
f(148) = 2*r0**2*(a1**2*a2 + a1*a2**2 + a3**2*a4 + a3*a4**2)
f(149) = 2*r0**2*(a1**3*a2 + a1*a2**3 + a3**3*a4 + a3*a4**3)
f(150) = 4*r0**2*(a1**2*a2**2 + a3**2*a4**2)
f(151) = 4*r0**3*(a1*a2 + a3*a4)
f(152) = 2*r0**3*(a1**2*a2 + a1*a2**2 + a3**2*a4 + a3*a4**2)
f(153) = 4*r0**4*(a1*a2 + a3*a4)
f(154) = 4*r0*(a1*a3 + a2*a4)
f(155) = 2*r0*(a1**2*a3 + a1*a3**2 + a2**2*a4 + a2*a4**2)
f(156) = 2*r0*(a1**3*a3 + a1*a3**3 + a2**3*a4 + a2*a4**3)
f(157) = 2*r0*(a1**4*a3 + a1*a3**4 + a2**4*a4 + a2*a4**4)
f(158) = 4*r0*(a1**2*a3**2 + a2**2*a4**2)
f(159) = 2*r0*(a1**3*a3**2 + a1**2*a3**3 + a2**3*a4**2 + a2**2*a4**3)
f(160) = 4*r0**2*(a1*a3 + a2*a4)
f(161) = 2*r0**2*(a1**2*a3 + a1*a3**2 + a2**2*a4 + a2*a4**2)
f(162) = 2*r0**2*(a1**3*a3 + a1*a3**3 + a2**3*a4 + a2*a4**3)
f(163) = 4*r0**2*(a1**2*a3**2 + a2**2*a4**2)
f(164) = 4*r0**3*(a1*a3 + a2*a4)
f(165) = 2*r0**3*(a1**2*a3 + a1*a3**2 + a2**2*a4 + a2*a4**2)
f(166) = 4*r0**4*(a1*a3 + a2*a4)
f(167) = 4*r0*(a1*a4 + a2*a3)
f(168) = 2*r0*(a1**2*a4 + a1*a4**2 + a2**2*a3 + a2*a3**2)
f(169) = 2*r0*(a1**3*a4 + a1*a4**3 + a2**3*a3 + a2*a3**3)
f(170) = 2*r0*(a1**4*a4 + a1*a4**4 + a2**4*a3 + a2*a3**4)
f(171) = 4*r0*(a1**2*a4**2 + a2**2*a3**2)
f(172) = 2*r0*(a1**3*a4**2 + a1**2*a4**3 + a2**3*a3**2 + a2**2*a3**3)
f(173) = 4*r0**2*(a1*a4 + a2*a3)
f(174) = 2*r0**2*(a1**2*a4 + a1*a4**2 + a2**2*a3 + a2*a3**2)
f(175) = 2*r0**2*(a1**3*a4 + a1*a4**3 + a2**3*a3 + a2*a3**3)
f(176) = 4*r0**2*(a1**2*a4**2 + a2**2*a3**2)
f(177) = 4*r0**3*(a1*a4 + a2*a3)
f(178) = 2*r0**3*(a1**2*a4 + a1*a4**2 + a2**2*a3 + a2*a3**2)
f(179) = 4*r0**4*(a1*a4 + a2*a3)
f(180) = 2*r0*(a1*b1**2 + a2*b1**2 + a3*b2**2 + a4*b2**2)
f(181) = 2*r0*(a1*b1**4 + a2*b1**4 + a3*b2**4 + a4*b2**4)
f(182) = 2*r0*(a1**2*b1**2 + a2**2*b1**2 + a3**2*b2**2 + a4**2*b2**2)
f(183) = 2*r0*(a1**3*b1**2 + a2**3*b1**2 + a3**3*b2**2 + a4**3*b2**2)
f(184) = 2*r0**2*(a1*b1**2 + a2*b1**2 + a3*b2**2 + a4*b2**2)
f(185) = 2*r0**2*(a1**2*b1**2 + a2**2*b1**2 + a3**2*b2**2 + a4**2*b2**2)
f(186) = 2*r0**3*(a1*b1**2 + a2*b1**2 + a3*b2**2 + a4*b2**2)
f(187) = 2*r0*(a1*b2**2 + a2*b2**2 + a3*b1**2 + a4*b1**2)
f(188) = 2*r0*(a1*b2**4 + a2*b2**4 + a3*b1**4 + a4*b1**4)
f(189) = 2*r0*(a1**2*b2**2 + a2**2*b2**2 + a3**2*b1**2 + a4**2*b1**2)
f(190) = 2*r0*(a1**3*b2**2 + a2**3*b2**2 + a3**3*b1**2 + a4**3*b1**2)
f(191) = 2*r0**2*(a1*b2**2 + a2*b2**2 + a3*b1**2 + a4*b1**2)
f(192) = 2*r0**2*(a1**2*b2**2 + a2**2*b2**2 + a3**2*b1**2 + a4**2*b1**2)
f(193) = 2*r0**3*(a1*b2**2 + a2*b2**2 + a3*b1**2 + a4*b1**2)
f(194) = 2*dtau**2*r0*(a1 + a2 + a3 + a4)
f(195) = 2*dtau**4*r0*(a1 + a2 + a3 + a4)
f(196) = 2*dtau**2*r0*(a1**2 + a2**2 + a3**2 + a4**2)
f(197) = 2*dtau**2*r0*(a1**3 + a2**3 + a3**3 + a4**3)
f(198) = 2*dtau**2*r0**2*(a1 + a2 + a3 + a4)
f(199) = 2*dtau**2*r0**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(200) = 2*dtau**2*r0**3*(a1 + a2 + a3 + a4)
f(201) = 8*b1*b2*r0
f(202) = 4*b1*b2*r0*(b1**2 + b2**2)
f(203) = 8*b1**2*b2**2*r0
f(204) = 8*b1*b2*r0**2
f(205) = 4*b1*b2*r0**2*(b1**2 + b2**2)
f(206) = 8*b1**2*b2**2*r0**2
f(207) = 8*b1*b2*r0**3
f(208) = 8*b1*b2*r0**4
f(209) = 4*dtau**2*r0*(b1**2 + b2**2)
f(210) = 4*dtau**2*r0**2*(b1**2 + b2**2)
f(211) = 2*r1*r2*r3 + 2*r1*r2*r4 + 2*r1*r3*r4 + 2*r2*r3*r4
f(212) = 2*r1**2*r3*r4 + 2*r1*r2*r3**2 + 2*r1*r2*r4**2 + 2*r2**2*r3*r4
f(213) = 2*r1**3*r3*r4 + 2*r1*r2*r3**3 + 2*r1*r2*r4**3 + 2*r2**3*r3*r4
f(214) = 2*r1**4*r3*r4 + 2*r1*r2*r3**4 + 2*r1*r2*r4**4 + 2*r2**4*r3*r4
f(215) = 2*r1**2*r2*r4 + 2*r1*r2**2*r3 + 2*r1*r3*r4**2 + 2*r2*r3**2*r4
f(216) = 2*r1**2*r2*r4**2 + 2*r1**2*r3*r4**2 + 2*r1*r2**2*r3**2 + 2*r2** &
      2*r3**2*r4
f(217) = 2*r1**3*r3*r4**2 + 2*r1**2*r2*r4**3 + 2*r1*r2**2*r3**3 + 2*r2** &
      3*r3**2*r4
f(218) = 2*r1**3*r2*r4 + 2*r1*r2**3*r3 + 2*r1*r3*r4**3 + 2*r2*r3**3*r4
f(219) = 2*r1**3*r2*r4**2 + 2*r1**2*r3*r4**3 + 2*r1*r2**3*r3**2 + 2*r2** &
      2*r3**3*r4
f(220) = 2*r1**4*r2*r4 + 2*r1*r2**4*r3 + 2*r1*r3*r4**4 + 2*r2*r3**4*r4
f(221) = 2*r1**2*r2*r3 + 2*r1*r2**2*r4 + 2*r1*r3**2*r4 + 2*r2*r3*r4**2
f(222) = 2*r1**2*r2*r3**2 + 2*r1**2*r3**2*r4 + 2*r1*r2**2*r4**2 + 2*r2** &
      2*r3*r4**2
f(223) = 2*r1**3*r3**2*r4 + 2*r1**2*r2*r3**3 + 2*r1*r2**2*r4**3 + 2*r2** &
      3*r3*r4**2
f(224) = 2*r1**2*r2**2*r3 + 2*r1**2*r2**2*r4 + 2*r1*r3**2*r4**2 + 2*r2* &
      r3**2*r4**2
f(225) = 2*r1**2*r2**2*r3**2 + 2*r1**2*r2**2*r4**2 + 2*r1**2*r3**2*r4**2 &
      + 2*r2**2*r3**2*r4**2
f(226) = 2*r1**3*r2**2*r4 + 2*r1**2*r2**3*r3 + 2*r1*r3**2*r4**3 + 2*r2* &
      r3**3*r4**2
f(227) = 2*r1**3*r2*r3 + 2*r1*r2**3*r4 + 2*r1*r3**3*r4 + 2*r2*r3*r4**3
f(228) = 2*r1**3*r2*r3**2 + 2*r1**2*r3**3*r4 + 2*r1*r2**3*r4**2 + 2*r2** &
      2*r3*r4**3
f(229) = 2*r1**3*r2**2*r3 + 2*r1**2*r2**3*r4 + 2*r1*r3**3*r4**2 + 2*r2* &
      r3**2*r4**3
f(230) = 2*r1**4*r2*r3 + 2*r1*r2**4*r4 + 2*r1*r3**4*r4 + 2*r2*r3*r4**4
f(231) = 2*a1*r1*r2 + 2*a2*r1*r2 + 2*a3*r3*r4 + 2*a4*r3*r4
f(232) = 2*a1**2*r1*r2 + 2*a2**2*r1*r2 + 2*a3**2*r3*r4 + 2*a4**2*r3*r4
f(233) = 2*a1**3*r1*r2 + 2*a2**3*r1*r2 + 2*a3**3*r3*r4 + 2*a4**3*r3*r4
f(234) = 2*a1**4*r1*r2 + 2*a2**4*r1*r2 + 2*a3**4*r3*r4 + 2*a4**4*r3*r4
f(235) = 2*a1*r1*r2**2 + 2*a2*r1**2*r2 + 2*a3*r3*r4**2 + 2*a4*r3**2*r4
f(236) = 2*a1**2*r1*r2**2 + 2*a2**2*r1**2*r2 + 2*a3**2*r3*r4**2 + 2*a4** &
      2*r3**2*r4
f(237) = 2*a1**3*r1*r2**2 + 2*a2**3*r1**2*r2 + 2*a3**3*r3*r4**2 + 2*a4** &
      3*r3**2*r4
f(238) = 2*a1*r1*r2**3 + 2*a2*r1**3*r2 + 2*a3*r3*r4**3 + 2*a4*r3**3*r4
f(239) = 2*a1**2*r1*r2**3 + 2*a2**2*r1**3*r2 + 2*a3**2*r3*r4**3 + 2*a4** &
      2*r3**3*r4
f(240) = 2*a1*r1*r2**4 + 2*a2*r1**4*r2 + 2*a3*r3*r4**4 + 2*a4*r3**4*r4
f(241) = 2*a1*r1**2*r2 + 2*a2*r1*r2**2 + 2*a3*r3**2*r4 + 2*a4*r3*r4**2
f(242) = 2*a1**2*r1**2*r2 + 2*a2**2*r1*r2**2 + 2*a3**2*r3**2*r4 + 2*a4** &
      2*r3*r4**2
f(243) = 2*a1**3*r1**2*r2 + 2*a2**3*r1*r2**2 + 2*a3**3*r3**2*r4 + 2*a4** &
      3*r3*r4**2
f(244) = 2*a1*r1**2*r2**2 + 2*a2*r1**2*r2**2 + 2*a3*r3**2*r4**2 + 2*a4* &
      r3**2*r4**2
f(245) = 2*a1**2*r1**2*r2**2 + 2*a2**2*r1**2*r2**2 + 2*a3**2*r3**2*r4**2 &
      + 2*a4**2*r3**2*r4**2
f(246) = 2*a1*r1**2*r2**3 + 2*a2*r1**3*r2**2 + 2*a3*r3**2*r4**3 + 2*a4* &
      r3**3*r4**2
f(247) = 2*a1*r1**3*r2 + 2*a2*r1*r2**3 + 2*a3*r3**3*r4 + 2*a4*r3*r4**3
f(248) = 2*a1**2*r1**3*r2 + 2*a2**2*r1*r2**3 + 2*a3**2*r3**3*r4 + 2*a4** &
      2*r3*r4**3
f(249) = 2*a1*r1**3*r2**2 + 2*a2*r1**2*r2**3 + 2*a3*r3**3*r4**2 + 2*a4* &
      r3**2*r4**3
f(250) = 2*a1*r1**4*r2 + 2*a2*r1*r2**4 + 2*a3*r3**4*r4 + 2*a4*r3*r4**4
f(251) = 2*a1*r3*r4 + 2*a2*r3*r4 + 2*a3*r1*r2 + 2*a4*r1*r2
f(252) = 2*a1**2*r3*r4 + 2*a2**2*r3*r4 + 2*a3**2*r1*r2 + 2*a4**2*r1*r2
f(253) = 2*a1**3*r3*r4 + 2*a2**3*r3*r4 + 2*a3**3*r1*r2 + 2*a4**3*r1*r2
f(254) = 2*a1**4*r3*r4 + 2*a2**4*r3*r4 + 2*a3**4*r1*r2 + 2*a4**4*r1*r2
f(255) = 2*a1*r3*r4**2 + 2*a2*r3**2*r4 + 2*a3*r1*r2**2 + 2*a4*r1**2*r2
f(256) = 2*a1**2*r3*r4**2 + 2*a2**2*r3**2*r4 + 2*a3**2*r1*r2**2 + 2*a4** &
      2*r1**2*r2
f(257) = 2*a1**3*r3*r4**2 + 2*a2**3*r3**2*r4 + 2*a3**3*r1*r2**2 + 2*a4** &
      3*r1**2*r2
f(258) = 2*a1*r3*r4**3 + 2*a2*r3**3*r4 + 2*a3*r1*r2**3 + 2*a4*r1**3*r2
f(259) = 2*a1**2*r3*r4**3 + 2*a2**2*r3**3*r4 + 2*a3**2*r1*r2**3 + 2*a4** &
      2*r1**3*r2
f(260) = 2*a1*r3*r4**4 + 2*a2*r3**4*r4 + 2*a3*r1*r2**4 + 2*a4*r1**4*r2
f(261) = 2*a1*r3**2*r4 + 2*a2*r3*r4**2 + 2*a3*r1**2*r2 + 2*a4*r1*r2**2
f(262) = 2*a1**2*r3**2*r4 + 2*a2**2*r3*r4**2 + 2*a3**2*r1**2*r2 + 2*a4** &
      2*r1*r2**2
f(263) = 2*a1**3*r3**2*r4 + 2*a2**3*r3*r4**2 + 2*a3**3*r1**2*r2 + 2*a4** &
      3*r1*r2**2
f(264) = 2*a1*r3**2*r4**2 + 2*a2*r3**2*r4**2 + 2*a3*r1**2*r2**2 + 2*a4* &
      r1**2*r2**2
f(265) = 2*a1**2*r3**2*r4**2 + 2*a2**2*r3**2*r4**2 + 2*a3**2*r1**2*r2**2 &
      + 2*a4**2*r1**2*r2**2
f(266) = 2*a1*r3**2*r4**3 + 2*a2*r3**3*r4**2 + 2*a3*r1**2*r2**3 + 2*a4* &
      r1**3*r2**2
f(267) = 2*a1*r3**3*r4 + 2*a2*r3*r4**3 + 2*a3*r1**3*r2 + 2*a4*r1*r2**3
f(268) = 2*a1**2*r3**3*r4 + 2*a2**2*r3*r4**3 + 2*a3**2*r1**3*r2 + 2*a4** &
      2*r1*r2**3
f(269) = 2*a1*r3**3*r4**2 + 2*a2*r3**2*r4**3 + 2*a3*r1**3*r2**2 + 2*a4* &
      r1**2*r2**3
f(270) = 2*a1*r3**4*r4 + 2*a2*r3*r4**4 + 2*a3*r1**4*r2 + 2*a4*r1*r2**4
f(271) = 4*b1**2*r1*r2 + 4*b2**2*r3*r4
f(272) = 4*b1**4*r1*r2 + 4*b2**4*r3*r4
f(273) = 2*b1**2*r1**2*r2 + 2*b1**2*r1*r2**2 + 2*b2**2*r3**2*r4 + 2*b2** &
      2*r3*r4**2
f(274) = 2*b1**2*r1**3*r2 + 2*b1**2*r1*r2**3 + 2*b2**2*r3**3*r4 + 2*b2** &
      2*r3*r4**3
f(275) = 4*b1**2*r1**2*r2**2 + 4*b2**2*r3**2*r4**2
f(276) = 4*b1**2*r3*r4 + 4*b2**2*r1*r2
f(277) = 4*b1**4*r3*r4 + 4*b2**4*r1*r2
f(278) = 2*b1**2*r3**2*r4 + 2*b1**2*r3*r4**2 + 2*b2**2*r1**2*r2 + 2*b2** &
      2*r1*r2**2
f(279) = 2*b1**2*r3**3*r4 + 2*b1**2*r3*r4**3 + 2*b2**2*r1**3*r2 + 2*b2** &
      2*r1*r2**3
f(280) = 4*b1**2*r3**2*r4**2 + 4*b2**2*r1**2*r2**2
f(281) = 4*dtau**2*(r1*r2 + r3*r4)
f(282) = 4*dtau**4*(r1*r2 + r3*r4)
f(283) = 2*dtau**2*(r1**2*r2 + r1*r2**2 + r3**2*r4 + r3*r4**2)
f(284) = 2*dtau**2*(r1**3*r2 + r1*r2**3 + r3**3*r4 + r3*r4**3)
f(285) = 4*dtau**2*(r1**2*r2**2 + r3**2*r4**2)
f(286) = 2*a1*r1*r3 + 2*a2*r2*r4 + 2*a3*r1*r3 + 2*a4*r2*r4
f(287) = 2*a1**2*r1*r3 + 2*a2**2*r2*r4 + 2*a3**2*r1*r3 + 2*a4**2*r2*r4
f(288) = 2*a1**3*r1*r3 + 2*a2**3*r2*r4 + 2*a3**3*r1*r3 + 2*a4**3*r2*r4
f(289) = 2*a1**4*r1*r3 + 2*a2**4*r2*r4 + 2*a3**4*r1*r3 + 2*a4**4*r2*r4
f(290) = 2*a1*r1*r3**2 + 2*a2*r2*r4**2 + 2*a3*r1**2*r3 + 2*a4*r2**2*r4
f(291) = 2*a1**2*r1*r3**2 + 2*a2**2*r2*r4**2 + 2*a3**2*r1**2*r3 + 2*a4** &
      2*r2**2*r4
f(292) = 2*a1**3*r1*r3**2 + 2*a2**3*r2*r4**2 + 2*a3**3*r1**2*r3 + 2*a4** &
      3*r2**2*r4
f(293) = 2*a1*r1*r3**3 + 2*a2*r2*r4**3 + 2*a3*r1**3*r3 + 2*a4*r2**3*r4
f(294) = 2*a1**2*r1*r3**3 + 2*a2**2*r2*r4**3 + 2*a3**2*r1**3*r3 + 2*a4** &
      2*r2**3*r4
f(295) = 2*a1*r1*r3**4 + 2*a2*r2*r4**4 + 2*a3*r1**4*r3 + 2*a4*r2**4*r4
f(296) = 2*a1*r1**2*r3 + 2*a2*r2**2*r4 + 2*a3*r1*r3**2 + 2*a4*r2*r4**2
f(297) = 2*a1**2*r1**2*r3 + 2*a2**2*r2**2*r4 + 2*a3**2*r1*r3**2 + 2*a4** &
      2*r2*r4**2
f(298) = 2*a1**3*r1**2*r3 + 2*a2**3*r2**2*r4 + 2*a3**3*r1*r3**2 + 2*a4** &
      3*r2*r4**2
f(299) = 2*a1*r1**2*r3**2 + 2*a2*r2**2*r4**2 + 2*a3*r1**2*r3**2 + 2*a4* &
      r2**2*r4**2
f(300) = 2*a1**2*r1**2*r3**2 + 2*a2**2*r2**2*r4**2 + 2*a3**2*r1**2*r3**2 &
      + 2*a4**2*r2**2*r4**2
f(301) = 2*a1*r1**2*r3**3 + 2*a2*r2**2*r4**3 + 2*a3*r1**3*r3**2 + 2*a4* &
      r2**3*r4**2
f(302) = 2*a1*r1**3*r3 + 2*a2*r2**3*r4 + 2*a3*r1*r3**3 + 2*a4*r2*r4**3
f(303) = 2*a1**2*r1**3*r3 + 2*a2**2*r2**3*r4 + 2*a3**2*r1*r3**3 + 2*a4** &
      2*r2*r4**3
f(304) = 2*a1*r1**3*r3**2 + 2*a2*r2**3*r4**2 + 2*a3*r1**2*r3**3 + 2*a4* &
      r2**2*r4**3
f(305) = 2*a1*r1**4*r3 + 2*a2*r2**4*r4 + 2*a3*r1*r3**4 + 2*a4*r2*r4**4
f(306) = 2*a1*r2*r4 + 2*a2*r1*r3 + 2*a3*r2*r4 + 2*a4*r1*r3
f(307) = 2*a1**2*r2*r4 + 2*a2**2*r1*r3 + 2*a3**2*r2*r4 + 2*a4**2*r1*r3
f(308) = 2*a1**3*r2*r4 + 2*a2**3*r1*r3 + 2*a3**3*r2*r4 + 2*a4**3*r1*r3
f(309) = 2*a1**4*r2*r4 + 2*a2**4*r1*r3 + 2*a3**4*r2*r4 + 2*a4**4*r1*r3
f(310) = 2*a1*r2*r4**2 + 2*a2*r1*r3**2 + 2*a3*r2**2*r4 + 2*a4*r1**2*r3
f(311) = 2*a1**2*r2*r4**2 + 2*a2**2*r1*r3**2 + 2*a3**2*r2**2*r4 + 2*a4** &
      2*r1**2*r3
f(312) = 2*a1**3*r2*r4**2 + 2*a2**3*r1*r3**2 + 2*a3**3*r2**2*r4 + 2*a4** &
      3*r1**2*r3
f(313) = 2*a1*r2*r4**3 + 2*a2*r1*r3**3 + 2*a3*r2**3*r4 + 2*a4*r1**3*r3
f(314) = 2*a1**2*r2*r4**3 + 2*a2**2*r1*r3**3 + 2*a3**2*r2**3*r4 + 2*a4** &
      2*r1**3*r3
f(315) = 2*a1*r2*r4**4 + 2*a2*r1*r3**4 + 2*a3*r2**4*r4 + 2*a4*r1**4*r3
f(316) = 2*a1*r2**2*r4 + 2*a2*r1**2*r3 + 2*a3*r2*r4**2 + 2*a4*r1*r3**2
f(317) = 2*a1**2*r2**2*r4 + 2*a2**2*r1**2*r3 + 2*a3**2*r2*r4**2 + 2*a4** &
      2*r1*r3**2
f(318) = 2*a1**3*r2**2*r4 + 2*a2**3*r1**2*r3 + 2*a3**3*r2*r4**2 + 2*a4** &
      3*r1*r3**2
f(319) = 2*a1*r2**2*r4**2 + 2*a2*r1**2*r3**2 + 2*a3*r2**2*r4**2 + 2*a4* &
      r1**2*r3**2
f(320) = 2*a1**2*r2**2*r4**2 + 2*a2**2*r1**2*r3**2 + 2*a3**2*r2**2*r4**2 &
      + 2*a4**2*r1**2*r3**2
f(321) = 2*a1*r2**2*r4**3 + 2*a2*r1**2*r3**3 + 2*a3*r2**3*r4**2 + 2*a4* &
      r1**3*r3**2
f(322) = 2*a1*r2**3*r4 + 2*a2*r1**3*r3 + 2*a3*r2*r4**3 + 2*a4*r1*r3**3
f(323) = 2*a1**2*r2**3*r4 + 2*a2**2*r1**3*r3 + 2*a3**2*r2*r4**3 + 2*a4** &
      2*r1*r3**3
f(324) = 2*a1*r2**3*r4**2 + 2*a2*r1**3*r3**2 + 2*a3*r2**2*r4**3 + 2*a4* &
      r1**2*r3**3
f(325) = 2*a1*r2**4*r4 + 2*a2*r1**4*r3 + 2*a3*r2*r4**4 + 2*a4*r1*r3**4
f(326) = 2*b1**2*r1*r3 + 2*b1**2*r2*r4 + 2*b2**2*r1*r3 + 2*b2**2*r2*r4
f(327) = 2*b1**4*r1*r3 + 2*b1**4*r2*r4 + 2*b2**4*r1*r3 + 2*b2**4*r2*r4
f(328) = 2*b1**2*r1*r3**2 + 2*b1**2*r2*r4**2 + 2*b2**2*r1**2*r3 + 2*b2** &
      2*r2**2*r4
f(329) = 2*b1**2*r1*r3**3 + 2*b1**2*r2*r4**3 + 2*b2**2*r1**3*r3 + 2*b2** &
      2*r2**3*r4
f(330) = 2*b1**2*r1**2*r3 + 2*b1**2*r2**2*r4 + 2*b2**2*r1*r3**2 + 2*b2** &
      2*r2*r4**2
f(331) = 2*b1**2*r1**2*r3**2 + 2*b1**2*r2**2*r4**2 + 2*b2**2*r1**2*r3**2 &
      + 2*b2**2*r2**2*r4**2
f(332) = 2*b1**2*r1**3*r3 + 2*b1**2*r2**3*r4 + 2*b2**2*r1*r3**3 + 2*b2** &
      2*r2*r4**3
f(333) = 4*dtau**2*(r1*r3 + r2*r4)
f(334) = 4*dtau**4*(r1*r3 + r2*r4)
f(335) = 2*dtau**2*(r1**2*r3 + r1*r3**2 + r2**2*r4 + r2*r4**2)
f(336) = 2*dtau**2*(r1**3*r3 + r1*r3**3 + r2**3*r4 + r2*r4**3)
f(337) = 4*dtau**2*(r1**2*r3**2 + r2**2*r4**2)
f(338) = 2*a1*r1*r4 + 2*a2*r2*r3 + 2*a3*r2*r3 + 2*a4*r1*r4
f(339) = 2*a1**2*r1*r4 + 2*a2**2*r2*r3 + 2*a3**2*r2*r3 + 2*a4**2*r1*r4
f(340) = 2*a1**3*r1*r4 + 2*a2**3*r2*r3 + 2*a3**3*r2*r3 + 2*a4**3*r1*r4
f(341) = 2*a1**4*r1*r4 + 2*a2**4*r2*r3 + 2*a3**4*r2*r3 + 2*a4**4*r1*r4
f(342) = 2*a1*r1*r4**2 + 2*a2*r2*r3**2 + 2*a3*r2**2*r3 + 2*a4*r1**2*r4
f(343) = 2*a1**2*r1*r4**2 + 2*a2**2*r2*r3**2 + 2*a3**2*r2**2*r3 + 2*a4** &
      2*r1**2*r4
f(344) = 2*a1**3*r1*r4**2 + 2*a2**3*r2*r3**2 + 2*a3**3*r2**2*r3 + 2*a4** &
      3*r1**2*r4
f(345) = 2*a1*r1*r4**3 + 2*a2*r2*r3**3 + 2*a3*r2**3*r3 + 2*a4*r1**3*r4
f(346) = 2*a1**2*r1*r4**3 + 2*a2**2*r2*r3**3 + 2*a3**2*r2**3*r3 + 2*a4** &
      2*r1**3*r4
f(347) = 2*a1*r1*r4**4 + 2*a2*r2*r3**4 + 2*a3*r2**4*r3 + 2*a4*r1**4*r4
f(348) = 2*a1*r1**2*r4 + 2*a2*r2**2*r3 + 2*a3*r2*r3**2 + 2*a4*r1*r4**2
f(349) = 2*a1**2*r1**2*r4 + 2*a2**2*r2**2*r3 + 2*a3**2*r2*r3**2 + 2*a4** &
      2*r1*r4**2
f(350) = 2*a1**3*r1**2*r4 + 2*a2**3*r2**2*r3 + 2*a3**3*r2*r3**2 + 2*a4** &
      3*r1*r4**2
f(351) = 2*a1*r1**2*r4**2 + 2*a2*r2**2*r3**2 + 2*a3*r2**2*r3**2 + 2*a4* &
      r1**2*r4**2
f(352) = 2*a1**2*r1**2*r4**2 + 2*a2**2*r2**2*r3**2 + 2*a3**2*r2**2*r3**2 &
      + 2*a4**2*r1**2*r4**2
f(353) = 2*a1*r1**2*r4**3 + 2*a2*r2**2*r3**3 + 2*a3*r2**3*r3**2 + 2*a4* &
      r1**3*r4**2
f(354) = 2*a1*r1**3*r4 + 2*a2*r2**3*r3 + 2*a3*r2*r3**3 + 2*a4*r1*r4**3
f(355) = 2*a1**2*r1**3*r4 + 2*a2**2*r2**3*r3 + 2*a3**2*r2*r3**3 + 2*a4** &
      2*r1*r4**3
f(356) = 2*a1*r1**3*r4**2 + 2*a2*r2**3*r3**2 + 2*a3*r2**2*r3**3 + 2*a4* &
      r1**2*r4**3
f(357) = 2*a1*r1**4*r4 + 2*a2*r2**4*r3 + 2*a3*r2*r3**4 + 2*a4*r1*r4**4
f(358) = 2*a1*r2*r3 + 2*a2*r1*r4 + 2*a3*r1*r4 + 2*a4*r2*r3
f(359) = 2*a1**2*r2*r3 + 2*a2**2*r1*r4 + 2*a3**2*r1*r4 + 2*a4**2*r2*r3
f(360) = 2*a1**3*r2*r3 + 2*a2**3*r1*r4 + 2*a3**3*r1*r4 + 2*a4**3*r2*r3
f(361) = 2*a1**4*r2*r3 + 2*a2**4*r1*r4 + 2*a3**4*r1*r4 + 2*a4**4*r2*r3
f(362) = 2*a1*r2*r3**2 + 2*a2*r1*r4**2 + 2*a3*r1**2*r4 + 2*a4*r2**2*r3
f(363) = 2*a1**2*r2*r3**2 + 2*a2**2*r1*r4**2 + 2*a3**2*r1**2*r4 + 2*a4** &
      2*r2**2*r3
f(364) = 2*a1**3*r2*r3**2 + 2*a2**3*r1*r4**2 + 2*a3**3*r1**2*r4 + 2*a4** &
      3*r2**2*r3
f(365) = 2*a1*r2*r3**3 + 2*a2*r1*r4**3 + 2*a3*r1**3*r4 + 2*a4*r2**3*r3
f(366) = 2*a1**2*r2*r3**3 + 2*a2**2*r1*r4**3 + 2*a3**2*r1**3*r4 + 2*a4** &
      2*r2**3*r3
f(367) = 2*a1*r2*r3**4 + 2*a2*r1*r4**4 + 2*a3*r1**4*r4 + 2*a4*r2**4*r3
f(368) = 2*a1*r2**2*r3 + 2*a2*r1**2*r4 + 2*a3*r1*r4**2 + 2*a4*r2*r3**2
f(369) = 2*a1**2*r2**2*r3 + 2*a2**2*r1**2*r4 + 2*a3**2*r1*r4**2 + 2*a4** &
      2*r2*r3**2
f(370) = 2*a1**3*r2**2*r3 + 2*a2**3*r1**2*r4 + 2*a3**3*r1*r4**2 + 2*a4** &
      3*r2*r3**2
f(371) = 2*a1*r2**2*r3**2 + 2*a2*r1**2*r4**2 + 2*a3*r1**2*r4**2 + 2*a4* &
      r2**2*r3**2
f(372) = 2*a1**2*r2**2*r3**2 + 2*a2**2*r1**2*r4**2 + 2*a3**2*r1**2*r4**2 &
      + 2*a4**2*r2**2*r3**2
f(373) = 2*a1*r2**2*r3**3 + 2*a2*r1**2*r4**3 + 2*a3*r1**3*r4**2 + 2*a4* &
      r2**3*r3**2
f(374) = 2*a1*r2**3*r3 + 2*a2*r1**3*r4 + 2*a3*r1*r4**3 + 2*a4*r2*r3**3
f(375) = 2*a1**2*r2**3*r3 + 2*a2**2*r1**3*r4 + 2*a3**2*r1*r4**3 + 2*a4** &
      2*r2*r3**3
f(376) = 2*a1*r2**3*r3**2 + 2*a2*r1**3*r4**2 + 2*a3*r1**2*r4**3 + 2*a4* &
      r2**2*r3**3
f(377) = 2*a1*r2**4*r3 + 2*a2*r1**4*r4 + 2*a3*r1*r4**4 + 2*a4*r2*r3**4
f(378) = 2*b1**2*r1*r4 + 2*b1**2*r2*r3 + 2*b2**2*r1*r4 + 2*b2**2*r2*r3
f(379) = 2*b1**4*r1*r4 + 2*b1**4*r2*r3 + 2*b2**4*r1*r4 + 2*b2**4*r2*r3
f(380) = 2*b1**2*r1*r4**2 + 2*b1**2*r2*r3**2 + 2*b2**2*r1**2*r4 + 2*b2** &
      2*r2**2*r3
f(381) = 2*b1**2*r1*r4**3 + 2*b1**2*r2*r3**3 + 2*b2**2*r1**3*r4 + 2*b2** &
      2*r2**3*r3
f(382) = 2*b1**2*r1**2*r4 + 2*b1**2*r2**2*r3 + 2*b2**2*r1*r4**2 + 2*b2** &
      2*r2*r3**2
f(383) = 2*b1**2*r1**2*r4**2 + 2*b1**2*r2**2*r3**2 + 2*b2**2*r1**2*r4**2 &
      + 2*b2**2*r2**2*r3**2
f(384) = 2*b1**2*r1**3*r4 + 2*b1**2*r2**3*r3 + 2*b2**2*r1*r4**3 + 2*b2** &
      2*r2*r3**3
f(385) = 4*dtau**2*(r1*r4 + r2*r3)
f(386) = 4*dtau**4*(r1*r4 + r2*r3)
f(387) = 2*dtau**2*(r1**2*r4 + r1*r4**2 + r2**2*r3 + r2*r3**2)
f(388) = 2*dtau**2*(r1**3*r4 + r1*r4**3 + r2**3*r3 + r2*r3**3)
f(389) = 4*dtau**2*(r1**2*r4**2 + r2**2*r3**2)
f(390) = 2*a1*a2*r1 + 2*a1*a2*r2 + 2*a3*a4*r3 + 2*a3*a4*r4
f(391) = 2*a1**2*a2*r2 + 2*a1*a2**2*r1 + 2*a3**2*a4*r4 + 2*a3*a4**2*r3
f(392) = 2*a1**3*a2*r2 + 2*a1*a2**3*r1 + 2*a3**3*a4*r4 + 2*a3*a4**3*r3
f(393) = 2*a1**4*a2*r2 + 2*a1*a2**4*r1 + 2*a3**4*a4*r4 + 2*a3*a4**4*r3
f(394) = 2*a1**2*a2*r1 + 2*a1*a2**2*r2 + 2*a3**2*a4*r3 + 2*a3*a4**2*r4
f(395) = 2*a1**2*a2**2*r1 + 2*a1**2*a2**2*r2 + 2*a3**2*a4**2*r3 + 2*a3** &
      2*a4**2*r4
f(396) = 2*a1**3*a2**2*r2 + 2*a1**2*a2**3*r1 + 2*a3**3*a4**2*r4 + 2*a3** &
      2*a4**3*r3
f(397) = 2*a1**3*a2*r1 + 2*a1*a2**3*r2 + 2*a3**3*a4*r3 + 2*a3*a4**3*r4
f(398) = 2*a1**3*a2**2*r1 + 2*a1**2*a2**3*r2 + 2*a3**3*a4**2*r3 + 2*a3** &
      2*a4**3*r4
f(399) = 2*a1**4*a2*r1 + 2*a1*a2**4*r2 + 2*a3**4*a4*r3 + 2*a3*a4**4*r4
f(400) = 2*a1*a2*r1**2 + 2*a1*a2*r2**2 + 2*a3*a4*r3**2 + 2*a3*a4*r4**2
f(401) = 2*a1**2*a2*r2**2 + 2*a1*a2**2*r1**2 + 2*a3**2*a4*r4**2 + 2*a3* &
      a4**2*r3**2
f(402) = 2*a1**3*a2*r2**2 + 2*a1*a2**3*r1**2 + 2*a3**3*a4*r4**2 + 2*a3* &
      a4**3*r3**2
f(403) = 2*a1**2*a2*r1**2 + 2*a1*a2**2*r2**2 + 2*a3**2*a4*r3**2 + 2*a3* &
      a4**2*r4**2
f(404) = 2*a1**2*a2**2*r1**2 + 2*a1**2*a2**2*r2**2 + 2*a3**2*a4**2*r3**2 &
      + 2*a3**2*a4**2*r4**2
f(405) = 2*a1**3*a2*r1**2 + 2*a1*a2**3*r2**2 + 2*a3**3*a4*r3**2 + 2*a3* &
      a4**3*r4**2
f(406) = 2*a1*a2*r1**3 + 2*a1*a2*r2**3 + 2*a3*a4*r3**3 + 2*a3*a4*r4**3
f(407) = 2*a1**2*a2*r2**3 + 2*a1*a2**2*r1**3 + 2*a3**2*a4*r4**3 + 2*a3* &
      a4**2*r3**3
f(408) = 2*a1**2*a2*r1**3 + 2*a1*a2**2*r2**3 + 2*a3**2*a4*r3**3 + 2*a3* &
      a4**2*r4**3
f(409) = 2*a1*a2*r1**4 + 2*a1*a2*r2**4 + 2*a3*a4*r3**4 + 2*a3*a4*r4**4
f(410) = 2*a1*a3*r1 + 2*a1*a3*r3 + 2*a2*a4*r2 + 2*a2*a4*r4
f(411) = 2*a1**2*a3*r3 + 2*a1*a3**2*r1 + 2*a2**2*a4*r4 + 2*a2*a4**2*r2
f(412) = 2*a1**3*a3*r3 + 2*a1*a3**3*r1 + 2*a2**3*a4*r4 + 2*a2*a4**3*r2
f(413) = 2*a1**4*a3*r3 + 2*a1*a3**4*r1 + 2*a2**4*a4*r4 + 2*a2*a4**4*r2
f(414) = 2*a1**2*a3*r1 + 2*a1*a3**2*r3 + 2*a2**2*a4*r2 + 2*a2*a4**2*r4
f(415) = 2*a1**2*a3**2*r1 + 2*a1**2*a3**2*r3 + 2*a2**2*a4**2*r2 + 2*a2** &
      2*a4**2*r4
f(416) = 2*a1**3*a3**2*r3 + 2*a1**2*a3**3*r1 + 2*a2**3*a4**2*r4 + 2*a2** &
      2*a4**3*r2
f(417) = 2*a1**3*a3*r1 + 2*a1*a3**3*r3 + 2*a2**3*a4*r2 + 2*a2*a4**3*r4
f(418) = 2*a1**3*a3**2*r1 + 2*a1**2*a3**3*r3 + 2*a2**3*a4**2*r2 + 2*a2** &
      2*a4**3*r4
f(419) = 2*a1**4*a3*r1 + 2*a1*a3**4*r3 + 2*a2**4*a4*r2 + 2*a2*a4**4*r4
f(420) = 2*a1*a3*r1**2 + 2*a1*a3*r3**2 + 2*a2*a4*r2**2 + 2*a2*a4*r4**2
f(421) = 2*a1**2*a3*r3**2 + 2*a1*a3**2*r1**2 + 2*a2**2*a4*r4**2 + 2*a2* &
      a4**2*r2**2
f(422) = 2*a1**3*a3*r3**2 + 2*a1*a3**3*r1**2 + 2*a2**3*a4*r4**2 + 2*a2* &
      a4**3*r2**2
f(423) = 2*a1**2*a3*r1**2 + 2*a1*a3**2*r3**2 + 2*a2**2*a4*r2**2 + 2*a2* &
      a4**2*r4**2
f(424) = 2*a1**2*a3**2*r1**2 + 2*a1**2*a3**2*r3**2 + 2*a2**2*a4**2*r2**2 &
      + 2*a2**2*a4**2*r4**2
f(425) = 2*a1**3*a3*r1**2 + 2*a1*a3**3*r3**2 + 2*a2**3*a4*r2**2 + 2*a2* &
      a4**3*r4**2
f(426) = 2*a1*a3*r1**3 + 2*a1*a3*r3**3 + 2*a2*a4*r2**3 + 2*a2*a4*r4**3
f(427) = 2*a1**2*a3*r3**3 + 2*a1*a3**2*r1**3 + 2*a2**2*a4*r4**3 + 2*a2* &
      a4**2*r2**3
f(428) = 2*a1**2*a3*r1**3 + 2*a1*a3**2*r3**3 + 2*a2**2*a4*r2**3 + 2*a2* &
      a4**2*r4**3
f(429) = 2*a1*a3*r1**4 + 2*a1*a3*r3**4 + 2*a2*a4*r2**4 + 2*a2*a4*r4**4
f(430) = 2*a1*a4*r1 + 2*a1*a4*r4 + 2*a2*a3*r2 + 2*a2*a3*r3
f(431) = 2*a1**2*a4*r4 + 2*a1*a4**2*r1 + 2*a2**2*a3*r3 + 2*a2*a3**2*r2
f(432) = 2*a1**3*a4*r4 + 2*a1*a4**3*r1 + 2*a2**3*a3*r3 + 2*a2*a3**3*r2
f(433) = 2*a1**4*a4*r4 + 2*a1*a4**4*r1 + 2*a2**4*a3*r3 + 2*a2*a3**4*r2
f(434) = 2*a1**2*a4*r1 + 2*a1*a4**2*r4 + 2*a2**2*a3*r2 + 2*a2*a3**2*r3
f(435) = 2*a1**2*a4**2*r1 + 2*a1**2*a4**2*r4 + 2*a2**2*a3**2*r2 + 2*a2** &
      2*a3**2*r3
f(436) = 2*a1**3*a4**2*r4 + 2*a1**2*a4**3*r1 + 2*a2**3*a3**2*r3 + 2*a2** &
      2*a3**3*r2
f(437) = 2*a1**3*a4*r1 + 2*a1*a4**3*r4 + 2*a2**3*a3*r2 + 2*a2*a3**3*r3
f(438) = 2*a1**3*a4**2*r1 + 2*a1**2*a4**3*r4 + 2*a2**3*a3**2*r2 + 2*a2** &
      2*a3**3*r3
f(439) = 2*a1**4*a4*r1 + 2*a1*a4**4*r4 + 2*a2**4*a3*r2 + 2*a2*a3**4*r3
f(440) = 2*a1*a4*r1**2 + 2*a1*a4*r4**2 + 2*a2*a3*r2**2 + 2*a2*a3*r3**2
f(441) = 2*a1**2*a4*r4**2 + 2*a1*a4**2*r1**2 + 2*a2**2*a3*r3**2 + 2*a2* &
      a3**2*r2**2
f(442) = 2*a1**3*a4*r4**2 + 2*a1*a4**3*r1**2 + 2*a2**3*a3*r3**2 + 2*a2* &
      a3**3*r2**2
f(443) = 2*a1**2*a4*r1**2 + 2*a1*a4**2*r4**2 + 2*a2**2*a3*r2**2 + 2*a2* &
      a3**2*r3**2
f(444) = 2*a1**2*a4**2*r1**2 + 2*a1**2*a4**2*r4**2 + 2*a2**2*a3**2*r2**2 &
      + 2*a2**2*a3**2*r3**2
f(445) = 2*a1**3*a4*r1**2 + 2*a1*a4**3*r4**2 + 2*a2**3*a3*r2**2 + 2*a2* &
      a3**3*r3**2
f(446) = 2*a1*a4*r1**3 + 2*a1*a4*r4**3 + 2*a2*a3*r2**3 + 2*a2*a3*r3**3
f(447) = 2*a1**2*a4*r4**3 + 2*a1*a4**2*r1**3 + 2*a2**2*a3*r3**3 + 2*a2* &
      a3**2*r2**3
f(448) = 2*a1**2*a4*r1**3 + 2*a1*a4**2*r4**3 + 2*a2**2*a3*r2**3 + 2*a2* &
      a3**2*r3**3
f(449) = 2*a1*a4*r1**4 + 2*a1*a4*r4**4 + 2*a2*a3*r2**4 + 2*a2*a3*r3**4
f(450) = 2*a1*b1**2*r1 + 2*a2*b1**2*r2 + 2*a3*b2**2*r3 + 2*a4*b2**2*r4
f(451) = 2*a1*b1**4*r1 + 2*a2*b1**4*r2 + 2*a3*b2**4*r3 + 2*a4*b2**4*r4
f(452) = 2*a1**2*b1**2*r1 + 2*a2**2*b1**2*r2 + 2*a3**2*b2**2*r3 + 2*a4** &
      2*b2**2*r4
f(453) = 2*a1**3*b1**2*r1 + 2*a2**3*b1**2*r2 + 2*a3**3*b2**2*r3 + 2*a4** &
      3*b2**2*r4
f(454) = 2*a1*b1**2*r1**2 + 2*a2*b1**2*r2**2 + 2*a3*b2**2*r3**2 + 2*a4* &
      b2**2*r4**2
f(455) = 2*a1**2*b1**2*r1**2 + 2*a2**2*b1**2*r2**2 + 2*a3**2*b2**2*r3**2 &
      + 2*a4**2*b2**2*r4**2
f(456) = 2*a1*b1**2*r1**3 + 2*a2*b1**2*r2**3 + 2*a3*b2**2*r3**3 + 2*a4* &
      b2**2*r4**3
f(457) = 2*a1*b2**2*r1 + 2*a2*b2**2*r2 + 2*a3*b1**2*r3 + 2*a4*b1**2*r4
f(458) = 2*a1*b2**4*r1 + 2*a2*b2**4*r2 + 2*a3*b1**4*r3 + 2*a4*b1**4*r4
f(459) = 2*a1**2*b2**2*r1 + 2*a2**2*b2**2*r2 + 2*a3**2*b1**2*r3 + 2*a4** &
      2*b1**2*r4
f(460) = 2*a1**3*b2**2*r1 + 2*a2**3*b2**2*r2 + 2*a3**3*b1**2*r3 + 2*a4** &
      3*b1**2*r4
f(461) = 2*a1*b2**2*r1**2 + 2*a2*b2**2*r2**2 + 2*a3*b1**2*r3**2 + 2*a4* &
      b1**2*r4**2
f(462) = 2*a1**2*b2**2*r1**2 + 2*a2**2*b2**2*r2**2 + 2*a3**2*b1**2*r3**2 &
      + 2*a4**2*b1**2*r4**2
f(463) = 2*a1*b2**2*r1**3 + 2*a2*b2**2*r2**3 + 2*a3*b1**2*r3**3 + 2*a4* &
      b1**2*r4**3
f(464) = 2*dtau**2*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(465) = 2*dtau**4*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(466) = 2*dtau**2*(a1**2*r1 + a2**2*r2 + a3**2*r3 + a4**2*r4)
f(467) = 2*dtau**2*(a1**3*r1 + a2**3*r2 + a3**3*r3 + a4**3*r4)
f(468) = 2*dtau**2*(a1*r1**2 + a2*r2**2 + a3*r3**2 + a4*r4**2)
f(469) = 2*dtau**2*(a1**2*r1**2 + a2**2*r2**2 + a3**2*r3**2 + a4**2*r4** &
      2)
f(470) = 2*dtau**2*(a1*r1**3 + a2*r2**3 + a3*r3**3 + a4*r4**3)
f(471) = 2*a1*a4*r2 + 2*a1*a4*r3 + 2*a2*a3*r1 + 2*a2*a3*r4
f(472) = 2*a1**2*a4*r3 + 2*a1*a4**2*r2 + 2*a2**2*a3*r4 + 2*a2*a3**2*r1
f(473) = 2*a1**3*a4*r3 + 2*a1*a4**3*r2 + 2*a2**3*a3*r4 + 2*a2*a3**3*r1
f(474) = 2*a1**4*a4*r3 + 2*a1*a4**4*r2 + 2*a2**4*a3*r4 + 2*a2*a3**4*r1
f(475) = 2*a1**2*a4*r2 + 2*a1*a4**2*r3 + 2*a2**2*a3*r1 + 2*a2*a3**2*r4
f(476) = 2*a1**2*a4**2*r2 + 2*a1**2*a4**2*r3 + 2*a2**2*a3**2*r1 + 2*a2** &
      2*a3**2*r4
f(477) = 2*a1**3*a4**2*r3 + 2*a1**2*a4**3*r2 + 2*a2**3*a3**2*r4 + 2*a2** &
      2*a3**3*r1
f(478) = 2*a1**3*a4*r2 + 2*a1*a4**3*r3 + 2*a2**3*a3*r1 + 2*a2*a3**3*r4
f(479) = 2*a1**3*a4**2*r2 + 2*a1**2*a4**3*r3 + 2*a2**3*a3**2*r1 + 2*a2** &
      2*a3**3*r4
f(480) = 2*a1**4*a4*r2 + 2*a1*a4**4*r3 + 2*a2**4*a3*r1 + 2*a2*a3**4*r4
f(481) = 2*a1*a4*r2**2 + 2*a1*a4*r3**2 + 2*a2*a3*r1**2 + 2*a2*a3*r4**2
f(482) = 2*a1**2*a4*r3**2 + 2*a1*a4**2*r2**2 + 2*a2**2*a3*r4**2 + 2*a2* &
      a3**2*r1**2
f(483) = 2*a1**3*a4*r3**2 + 2*a1*a4**3*r2**2 + 2*a2**3*a3*r4**2 + 2*a2* &
      a3**3*r1**2
f(484) = 2*a1**2*a4*r2**2 + 2*a1*a4**2*r3**2 + 2*a2**2*a3*r1**2 + 2*a2* &
      a3**2*r4**2
f(485) = 2*a1**2*a4**2*r2**2 + 2*a1**2*a4**2*r3**2 + 2*a2**2*a3**2*r1**2 &
      + 2*a2**2*a3**2*r4**2
f(486) = 2*a1**3*a4*r2**2 + 2*a1*a4**3*r3**2 + 2*a2**3*a3*r1**2 + 2*a2* &
      a3**3*r4**2
f(487) = 2*a1*a4*r2**3 + 2*a1*a4*r3**3 + 2*a2*a3*r1**3 + 2*a2*a3*r4**3
f(488) = 2*a1**2*a4*r3**3 + 2*a1*a4**2*r2**3 + 2*a2**2*a3*r4**3 + 2*a2* &
      a3**2*r1**3
f(489) = 2*a1**2*a4*r2**3 + 2*a1*a4**2*r3**3 + 2*a2**2*a3*r1**3 + 2*a2* &
      a3**2*r4**3
f(490) = 2*a1*a4*r2**4 + 2*a1*a4*r3**4 + 2*a2*a3*r1**4 + 2*a2*a3*r4**4
f(491) = 2*a1*a3*r2 + 2*a1*a3*r4 + 2*a2*a4*r1 + 2*a2*a4*r3
f(492) = 2*a1**2*a3*r4 + 2*a1*a3**2*r2 + 2*a2**2*a4*r3 + 2*a2*a4**2*r1
f(493) = 2*a1**3*a3*r4 + 2*a1*a3**3*r2 + 2*a2**3*a4*r3 + 2*a2*a4**3*r1
f(494) = 2*a1**4*a3*r4 + 2*a1*a3**4*r2 + 2*a2**4*a4*r3 + 2*a2*a4**4*r1
f(495) = 2*a1**2*a3*r2 + 2*a1*a3**2*r4 + 2*a2**2*a4*r1 + 2*a2*a4**2*r3
f(496) = 2*a1**2*a3**2*r2 + 2*a1**2*a3**2*r4 + 2*a2**2*a4**2*r1 + 2*a2** &
      2*a4**2*r3
f(497) = 2*a1**3*a3**2*r4 + 2*a1**2*a3**3*r2 + 2*a2**3*a4**2*r3 + 2*a2** &
      2*a4**3*r1
f(498) = 2*a1**3*a3*r2 + 2*a1*a3**3*r4 + 2*a2**3*a4*r1 + 2*a2*a4**3*r3
f(499) = 2*a1**3*a3**2*r2 + 2*a1**2*a3**3*r4 + 2*a2**3*a4**2*r1 + 2*a2** &
      2*a4**3*r3
f(500) = 2*a1**4*a3*r2 + 2*a1*a3**4*r4 + 2*a2**4*a4*r1 + 2*a2*a4**4*r3
f(501) = 2*a1*a3*r2**2 + 2*a1*a3*r4**2 + 2*a2*a4*r1**2 + 2*a2*a4*r3**2
f(502) = 2*a1**2*a3*r4**2 + 2*a1*a3**2*r2**2 + 2*a2**2*a4*r3**2 + 2*a2* &
      a4**2*r1**2
f(503) = 2*a1**3*a3*r4**2 + 2*a1*a3**3*r2**2 + 2*a2**3*a4*r3**2 + 2*a2* &
      a4**3*r1**2
f(504) = 2*a1**2*a3*r2**2 + 2*a1*a3**2*r4**2 + 2*a2**2*a4*r1**2 + 2*a2* &
      a4**2*r3**2
f(505) = 2*a1**2*a3**2*r2**2 + 2*a1**2*a3**2*r4**2 + 2*a2**2*a4**2*r1**2 &
      + 2*a2**2*a4**2*r3**2
f(506) = 2*a1**3*a3*r2**2 + 2*a1*a3**3*r4**2 + 2*a2**3*a4*r1**2 + 2*a2* &
      a4**3*r3**2
f(507) = 2*a1*a3*r2**3 + 2*a1*a3*r4**3 + 2*a2*a4*r1**3 + 2*a2*a4*r3**3
f(508) = 2*a1**2*a3*r4**3 + 2*a1*a3**2*r2**3 + 2*a2**2*a4*r3**3 + 2*a2* &
      a4**2*r1**3
f(509) = 2*a1**2*a3*r2**3 + 2*a1*a3**2*r4**3 + 2*a2**2*a4*r1**3 + 2*a2* &
      a4**2*r3**3
f(510) = 2*a1*a3*r2**4 + 2*a1*a3*r4**4 + 2*a2*a4*r1**4 + 2*a2*a4*r3**4
f(511) = 2*a1*b1**2*r2 + 2*a2*b1**2*r1 + 2*a3*b2**2*r4 + 2*a4*b2**2*r3
f(512) = 2*a1*b1**4*r2 + 2*a2*b1**4*r1 + 2*a3*b2**4*r4 + 2*a4*b2**4*r3
f(513) = 2*a1**2*b1**2*r2 + 2*a2**2*b1**2*r1 + 2*a3**2*b2**2*r4 + 2*a4** &
      2*b2**2*r3
f(514) = 2*a1**3*b1**2*r2 + 2*a2**3*b1**2*r1 + 2*a3**3*b2**2*r4 + 2*a4** &
      3*b2**2*r3
f(515) = 2*a1*b1**2*r2**2 + 2*a2*b1**2*r1**2 + 2*a3*b2**2*r4**2 + 2*a4* &
      b2**2*r3**2
f(516) = 2*a1**2*b1**2*r2**2 + 2*a2**2*b1**2*r1**2 + 2*a3**2*b2**2*r4**2 &
      + 2*a4**2*b2**2*r3**2
f(517) = 2*a1*b1**2*r2**3 + 2*a2*b1**2*r1**3 + 2*a3*b2**2*r4**3 + 2*a4* &
      b2**2*r3**3
f(518) = 2*a1*b2**2*r2 + 2*a2*b2**2*r1 + 2*a3*b1**2*r4 + 2*a4*b1**2*r3
f(519) = 2*a1*b2**4*r2 + 2*a2*b2**4*r1 + 2*a3*b1**4*r4 + 2*a4*b1**4*r3
f(520) = 2*a1**2*b2**2*r2 + 2*a2**2*b2**2*r1 + 2*a3**2*b1**2*r4 + 2*a4** &
      2*b1**2*r3
f(521) = 2*a1**3*b2**2*r2 + 2*a2**3*b2**2*r1 + 2*a3**3*b1**2*r4 + 2*a4** &
      3*b1**2*r3
f(522) = 2*a1*b2**2*r2**2 + 2*a2*b2**2*r1**2 + 2*a3*b1**2*r4**2 + 2*a4* &
      b1**2*r3**2
f(523) = 2*a1**2*b2**2*r2**2 + 2*a2**2*b2**2*r1**2 + 2*a3**2*b1**2*r4**2 &
      + 2*a4**2*b1**2*r3**2
f(524) = 2*a1*b2**2*r2**3 + 2*a2*b2**2*r1**3 + 2*a3*b1**2*r4**3 + 2*a4* &
      b1**2*r3**3
f(525) = 2*dtau**2*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(526) = 2*dtau**4*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(527) = 2*dtau**2*(a1**2*r2 + a2**2*r1 + a3**2*r4 + a4**2*r3)
f(528) = 2*dtau**2*(a1**3*r2 + a2**3*r1 + a3**3*r4 + a4**3*r3)
f(529) = 2*dtau**2*(a1*r2**2 + a2*r1**2 + a3*r4**2 + a4*r3**2)
f(530) = 2*dtau**2*(a1**2*r2**2 + a2**2*r1**2 + a3**2*r4**2 + a4**2*r3** &
      2)
f(531) = 2*dtau**2*(a1*r2**3 + a2*r1**3 + a3*r4**3 + a4*r3**3)
f(532) = 2*a1*a2*r3 + 2*a1*a2*r4 + 2*a3*a4*r1 + 2*a3*a4*r2
f(533) = 2*a1**2*a2*r4 + 2*a1*a2**2*r3 + 2*a3**2*a4*r2 + 2*a3*a4**2*r1
f(534) = 2*a1**3*a2*r4 + 2*a1*a2**3*r3 + 2*a3**3*a4*r2 + 2*a3*a4**3*r1
f(535) = 2*a1**4*a2*r4 + 2*a1*a2**4*r3 + 2*a3**4*a4*r2 + 2*a3*a4**4*r1
f(536) = 2*a1**2*a2*r3 + 2*a1*a2**2*r4 + 2*a3**2*a4*r1 + 2*a3*a4**2*r2
f(537) = 2*a1**2*a2**2*r3 + 2*a1**2*a2**2*r4 + 2*a3**2*a4**2*r1 + 2*a3** &
      2*a4**2*r2
f(538) = 2*a1**3*a2**2*r4 + 2*a1**2*a2**3*r3 + 2*a3**3*a4**2*r2 + 2*a3** &
      2*a4**3*r1
f(539) = 2*a1**3*a2*r3 + 2*a1*a2**3*r4 + 2*a3**3*a4*r1 + 2*a3*a4**3*r2
f(540) = 2*a1**3*a2**2*r3 + 2*a1**2*a2**3*r4 + 2*a3**3*a4**2*r1 + 2*a3** &
      2*a4**3*r2
f(541) = 2*a1**4*a2*r3 + 2*a1*a2**4*r4 + 2*a3**4*a4*r1 + 2*a3*a4**4*r2
f(542) = 2*a1*a2*r3**2 + 2*a1*a2*r4**2 + 2*a3*a4*r1**2 + 2*a3*a4*r2**2
f(543) = 2*a1**2*a2*r4**2 + 2*a1*a2**2*r3**2 + 2*a3**2*a4*r2**2 + 2*a3* &
      a4**2*r1**2
f(544) = 2*a1**3*a2*r4**2 + 2*a1*a2**3*r3**2 + 2*a3**3*a4*r2**2 + 2*a3* &
      a4**3*r1**2
f(545) = 2*a1**2*a2*r3**2 + 2*a1*a2**2*r4**2 + 2*a3**2*a4*r1**2 + 2*a3* &
      a4**2*r2**2
f(546) = 2*a1**2*a2**2*r3**2 + 2*a1**2*a2**2*r4**2 + 2*a3**2*a4**2*r1**2 &
      + 2*a3**2*a4**2*r2**2
f(547) = 2*a1**3*a2*r3**2 + 2*a1*a2**3*r4**2 + 2*a3**3*a4*r1**2 + 2*a3* &
      a4**3*r2**2
f(548) = 2*a1*a2*r3**3 + 2*a1*a2*r4**3 + 2*a3*a4*r1**3 + 2*a3*a4*r2**3
f(549) = 2*a1**2*a2*r4**3 + 2*a1*a2**2*r3**3 + 2*a3**2*a4*r2**3 + 2*a3* &
      a4**2*r1**3
f(550) = 2*a1**2*a2*r3**3 + 2*a1*a2**2*r4**3 + 2*a3**2*a4*r1**3 + 2*a3* &
      a4**2*r2**3
f(551) = 2*a1*a2*r3**4 + 2*a1*a2*r4**4 + 2*a3*a4*r1**4 + 2*a3*a4*r2**4
f(552) = 2*a1*b2**2*r3 + 2*a2*b2**2*r4 + 2*a3*b1**2*r1 + 2*a4*b1**2*r2
f(553) = 2*a1*b2**4*r3 + 2*a2*b2**4*r4 + 2*a3*b1**4*r1 + 2*a4*b1**4*r2
f(554) = 2*a1**2*b2**2*r3 + 2*a2**2*b2**2*r4 + 2*a3**2*b1**2*r1 + 2*a4** &
      2*b1**2*r2
f(555) = 2*a1**3*b2**2*r3 + 2*a2**3*b2**2*r4 + 2*a3**3*b1**2*r1 + 2*a4** &
      3*b1**2*r2
f(556) = 2*a1*b2**2*r3**2 + 2*a2*b2**2*r4**2 + 2*a3*b1**2*r1**2 + 2*a4* &
      b1**2*r2**2
f(557) = 2*a1**2*b2**2*r3**2 + 2*a2**2*b2**2*r4**2 + 2*a3**2*b1**2*r1**2 &
      + 2*a4**2*b1**2*r2**2
f(558) = 2*a1*b2**2*r3**3 + 2*a2*b2**2*r4**3 + 2*a3*b1**2*r1**3 + 2*a4* &
      b1**2*r2**3
f(559) = 2*a1*b1**2*r3 + 2*a2*b1**2*r4 + 2*a3*b2**2*r1 + 2*a4*b2**2*r2
f(560) = 2*a1*b1**4*r3 + 2*a2*b1**4*r4 + 2*a3*b2**4*r1 + 2*a4*b2**4*r2
f(561) = 2*a1**2*b1**2*r3 + 2*a2**2*b1**2*r4 + 2*a3**2*b2**2*r1 + 2*a4** &
      2*b2**2*r2
f(562) = 2*a1**3*b1**2*r3 + 2*a2**3*b1**2*r4 + 2*a3**3*b2**2*r1 + 2*a4** &
      3*b2**2*r2
f(563) = 2*a1*b1**2*r3**2 + 2*a2*b1**2*r4**2 + 2*a3*b2**2*r1**2 + 2*a4* &
      b2**2*r2**2
f(564) = 2*a1**2*b1**2*r3**2 + 2*a2**2*b1**2*r4**2 + 2*a3**2*b2**2*r1**2 &
      + 2*a4**2*b2**2*r2**2
f(565) = 2*a1*b1**2*r3**3 + 2*a2*b1**2*r4**3 + 2*a3*b2**2*r1**3 + 2*a4* &
      b2**2*r2**3
f(566) = 2*dtau**2*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(567) = 2*dtau**4*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(568) = 2*dtau**2*(a1**2*r3 + a2**2*r4 + a3**2*r1 + a4**2*r2)
f(569) = 2*dtau**2*(a1**3*r3 + a2**3*r4 + a3**3*r1 + a4**3*r2)
f(570) = 2*dtau**2*(a1*r3**2 + a2*r4**2 + a3*r1**2 + a4*r2**2)
f(571) = 2*dtau**2*(a1**2*r3**2 + a2**2*r4**2 + a3**2*r1**2 + a4**2*r2** &
      2)
f(572) = 2*dtau**2*(a1*r3**3 + a2*r4**3 + a3*r1**3 + a4*r2**3)
f(573) = 2*a1*b2**2*r4 + 2*a2*b2**2*r3 + 2*a3*b1**2*r2 + 2*a4*b1**2*r1
f(574) = 2*a1*b2**4*r4 + 2*a2*b2**4*r3 + 2*a3*b1**4*r2 + 2*a4*b1**4*r1
f(575) = 2*a1**2*b2**2*r4 + 2*a2**2*b2**2*r3 + 2*a3**2*b1**2*r2 + 2*a4** &
      2*b1**2*r1
f(576) = 2*a1**3*b2**2*r4 + 2*a2**3*b2**2*r3 + 2*a3**3*b1**2*r2 + 2*a4** &
      3*b1**2*r1
f(577) = 2*a1*b2**2*r4**2 + 2*a2*b2**2*r3**2 + 2*a3*b1**2*r2**2 + 2*a4* &
      b1**2*r1**2
f(578) = 2*a1**2*b2**2*r4**2 + 2*a2**2*b2**2*r3**2 + 2*a3**2*b1**2*r2**2 &
      + 2*a4**2*b1**2*r1**2
f(579) = 2*a1*b2**2*r4**3 + 2*a2*b2**2*r3**3 + 2*a3*b1**2*r2**3 + 2*a4* &
      b1**2*r1**3
f(580) = 2*a1*b1**2*r4 + 2*a2*b1**2*r3 + 2*a3*b2**2*r2 + 2*a4*b2**2*r1
f(581) = 2*a1*b1**4*r4 + 2*a2*b1**4*r3 + 2*a3*b2**4*r2 + 2*a4*b2**4*r1
f(582) = 2*a1**2*b1**2*r4 + 2*a2**2*b1**2*r3 + 2*a3**2*b2**2*r2 + 2*a4** &
      2*b2**2*r1
f(583) = 2*a1**3*b1**2*r4 + 2*a2**3*b1**2*r3 + 2*a3**3*b2**2*r2 + 2*a4** &
      3*b2**2*r1
f(584) = 2*a1*b1**2*r4**2 + 2*a2*b1**2*r3**2 + 2*a3*b2**2*r2**2 + 2*a4* &
      b2**2*r1**2
f(585) = 2*a1**2*b1**2*r4**2 + 2*a2**2*b1**2*r3**2 + 2*a3**2*b2**2*r2**2 &
      + 2*a4**2*b2**2*r1**2
f(586) = 2*a1*b1**2*r4**3 + 2*a2*b1**2*r3**3 + 2*a3*b2**2*r2**3 + 2*a4* &
      b2**2*r1**3
f(587) = 2*dtau**2*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(588) = 2*dtau**4*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(589) = 2*dtau**2*(a1**2*r4 + a2**2*r3 + a3**2*r2 + a4**2*r1)
f(590) = 2*dtau**2*(a1**3*r4 + a2**3*r3 + a3**3*r2 + a4**3*r1)
f(591) = 2*dtau**2*(a1*r4**2 + a2*r3**2 + a3*r2**2 + a4*r1**2)
f(592) = 2*dtau**2*(a1**2*r4**2 + a2**2*r3**2 + a3**2*r2**2 + a4**2*r1** &
      2)
f(593) = 2*dtau**2*(a1*r4**3 + a2*r3**3 + a3*r2**3 + a4*r1**3)
f(594) = 2*b1*b2*(r1 + r2 + r3 + r4)
f(595) = 2*b1*b2*(b1**2*r3 + b1**2*r4 + b2**2*r1 + b2**2*r2)
f(596) = 2*b1**2*b2**2*(r1 + r2 + r3 + r4)
f(597) = 2*b1*b2*(b1**2*r1 + b1**2*r2 + b2**2*r3 + b2**2*r4)
f(598) = 2*b1*b2*(r1**2 + r2**2 + r3**2 + r4**2)
f(599) = 2*b1*b2*(b1**2*r3**2 + b1**2*r4**2 + b2**2*r1**2 + b2**2*r2**2)
f(600) = 2*b1**2*b2**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(601) = 2*b1*b2*(b1**2*r1**2 + b1**2*r2**2 + b2**2*r3**2 + b2**2*r4**2)
f(602) = 2*b1*b2*(r1**3 + r2**3 + r3**3 + r4**3)
f(603) = 2*b1*b2*(r1**4 + r2**4 + r3**4 + r4**4)
f(604) = 2*dtau*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(605) = 2*dtau**3*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(606) = 2*dtau**2*(b1**2*r1 + b1**2*r2 + b2**2*r3 + b2**2*r4)
f(607) = 2*dtau*(b1**3*r1 - b1**3*r2 + b2**3*r3 - b2**3*r4)
f(608) = 2*dtau*(b1*r1**2 - b1*r2**2 + b2*r3**2 - b2*r4**2)
f(609) = 2*dtau**3*(b1*r1**2 - b1*r2**2 + b2*r3**2 - b2*r4**2)
f(610) = 2*dtau**2*(b1**2*r1**2 + b1**2*r2**2 + b2**2*r3**2 + b2**2*r4** &
      2)
f(611) = 2*dtau*(b1**3*r1**2 - b1**3*r2**2 + b2**3*r3**2 - b2**3*r4**2)
f(612) = 2*dtau*(b1*r1**3 - b1*r2**3 + b2*r3**3 - b2*r4**3)
f(613) = 2*dtau*(b1*r1**4 - b1*r2**4 + b2*r3**4 - b2*r4**4)
f(614) = 2*dtau*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(615) = 2*dtau**3*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(616) = 2*dtau**2*(b1**2*r3 + b1**2*r4 + b2**2*r1 + b2**2*r2)
f(617) = 2*dtau*(b1**3*r3 - b1**3*r4 + b2**3*r1 - b2**3*r2)
f(618) = 2*dtau*(b1*r3**2 - b1*r4**2 + b2*r1**2 - b2*r2**2)
f(619) = 2*dtau**3*(b1*r3**2 - b1*r4**2 + b2*r1**2 - b2*r2**2)
f(620) = 2*dtau**2*(b1**2*r3**2 + b1**2*r4**2 + b2**2*r1**2 + b2**2*r2** &
      2)
f(621) = 2*dtau*(b1**3*r3**2 - b1**3*r4**2 + b2**3*r1**2 - b2**3*r2**2)
f(622) = 2*dtau*(b1*r3**3 - b1*r4**3 + b2*r1**3 - b2*r2**3)
f(623) = 2*dtau*(b1*r3**4 - b1*r4**4 + b2*r1**4 - b2*r2**4)
f(624) = 2*a1*a2*a3 + 2*a1*a2*a4 + 2*a1*a3*a4 + 2*a2*a3*a4
f(625) = 2*a1**2*a3*a4 + 2*a1*a2*a3**2 + 2*a1*a2*a4**2 + 2*a2**2*a3*a4
f(626) = 2*a1**3*a3*a4 + 2*a1*a2*a3**3 + 2*a1*a2*a4**3 + 2*a2**3*a3*a4
f(627) = 2*a1**4*a3*a4 + 2*a1*a2*a3**4 + 2*a1*a2*a4**4 + 2*a2**4*a3*a4
f(628) = 2*a1**2*a2*a4 + 2*a1*a2**2*a3 + 2*a1*a3*a4**2 + 2*a2*a3**2*a4
f(629) = 2*a1**2*a2*a4**2 + 2*a1**2*a3*a4**2 + 2*a1*a2**2*a3**2 + 2*a2** &
      2*a3**2*a4
f(630) = 2*a1**3*a3*a4**2 + 2*a1**2*a2*a4**3 + 2*a1*a2**2*a3**3 + 2*a2** &
      3*a3**2*a4
f(631) = 2*a1**3*a2*a4 + 2*a1*a2**3*a3 + 2*a1*a3*a4**3 + 2*a2*a3**3*a4
f(632) = 2*a1**3*a2*a4**2 + 2*a1**2*a3*a4**3 + 2*a1*a2**3*a3**2 + 2*a2** &
      2*a3**3*a4
f(633) = 2*a1**4*a2*a4 + 2*a1*a2**4*a3 + 2*a1*a3*a4**4 + 2*a2*a3**4*a4
f(634) = 2*a1**2*a2*a3 + 2*a1*a2**2*a4 + 2*a1*a3**2*a4 + 2*a2*a3*a4**2
f(635) = 2*a1**2*a2*a3**2 + 2*a1**2*a3**2*a4 + 2*a1*a2**2*a4**2 + 2*a2** &
      2*a3*a4**2
f(636) = 2*a1**3*a3**2*a4 + 2*a1**2*a2*a3**3 + 2*a1*a2**2*a4**3 + 2*a2** &
      3*a3*a4**2
f(637) = 2*a1**2*a2**2*a3 + 2*a1**2*a2**2*a4 + 2*a1*a3**2*a4**2 + 2*a2* &
      a3**2*a4**2
f(638) = 2*a1**2*a2**2*a3**2 + 2*a1**2*a2**2*a4**2 + 2*a1**2*a3**2*a4**2 &
      + 2*a2**2*a3**2*a4**2
f(639) = 2*a1**3*a2**2*a4 + 2*a1**2*a2**3*a3 + 2*a1*a3**2*a4**3 + 2*a2* &
      a3**3*a4**2
f(640) = 2*a1**3*a2*a3 + 2*a1*a2**3*a4 + 2*a1*a3**3*a4 + 2*a2*a3*a4**3
f(641) = 2*a1**3*a2*a3**2 + 2*a1**2*a3**3*a4 + 2*a1*a2**3*a4**2 + 2*a2** &
      2*a3*a4**3
f(642) = 2*a1**3*a2**2*a3 + 2*a1**2*a2**3*a4 + 2*a1*a3**3*a4**2 + 2*a2* &
      a3**2*a4**3
f(643) = 2*a1**4*a2*a3 + 2*a1*a2**4*a4 + 2*a1*a3**4*a4 + 2*a2*a3*a4**4
f(644) = 4*a1*a2*b1**2 + 4*a3*a4*b2**2
f(645) = 4*a1*a2*b1**4 + 4*a3*a4*b2**4
f(646) = 2*a1**2*a2*b1**2 + 2*a1*a2**2*b1**2 + 2*a3**2*a4*b2**2 + 2*a3* &
      a4**2*b2**2
f(647) = 2*a1**3*a2*b1**2 + 2*a1*a2**3*b1**2 + 2*a3**3*a4*b2**2 + 2*a3* &
      a4**3*b2**2
f(648) = 4*a1**2*a2**2*b1**2 + 4*a3**2*a4**2*b2**2
f(649) = 4*a1*a2*b2**2 + 4*a3*a4*b1**2
f(650) = 4*a1*a2*b2**4 + 4*a3*a4*b1**4
f(651) = 2*a1**2*a2*b2**2 + 2*a1*a2**2*b2**2 + 2*a3**2*a4*b1**2 + 2*a3* &
      a4**2*b1**2
f(652) = 2*a1**3*a2*b2**2 + 2*a1*a2**3*b2**2 + 2*a3**3*a4*b1**2 + 2*a3* &
      a4**3*b1**2
f(653) = 4*a1**2*a2**2*b2**2 + 4*a3**2*a4**2*b1**2
f(654) = 4*dtau**2*(a1*a2 + a3*a4)
f(655) = 4*dtau**4*(a1*a2 + a3*a4)
f(656) = 2*dtau**2*(a1**2*a2 + a1*a2**2 + a3**2*a4 + a3*a4**2)
f(657) = 2*dtau**2*(a1**3*a2 + a1*a2**3 + a3**3*a4 + a3*a4**3)
f(658) = 4*dtau**2*(a1**2*a2**2 + a3**2*a4**2)
f(659) = 2*a1*a3*b1**2 + 2*a1*a3*b2**2 + 2*a2*a4*b1**2 + 2*a2*a4*b2**2
f(660) = 2*a1*a3*b1**4 + 2*a1*a3*b2**4 + 2*a2*a4*b1**4 + 2*a2*a4*b2**4
f(661) = 2*a1**2*a3*b2**2 + 2*a1*a3**2*b1**2 + 2*a2**2*a4*b2**2 + 2*a2* &
      a4**2*b1**2
f(662) = 2*a1**3*a3*b2**2 + 2*a1*a3**3*b1**2 + 2*a2**3*a4*b2**2 + 2*a2* &
      a4**3*b1**2
f(663) = 2*a1**2*a3*b1**2 + 2*a1*a3**2*b2**2 + 2*a2**2*a4*b1**2 + 2*a2* &
      a4**2*b2**2
f(664) = 2*a1**2*a3**2*b1**2 + 2*a1**2*a3**2*b2**2 + 2*a2**2*a4**2*b1**2 &
      + 2*a2**2*a4**2*b2**2
f(665) = 2*a1**3*a3*b1**2 + 2*a1*a3**3*b2**2 + 2*a2**3*a4*b1**2 + 2*a2* &
      a4**3*b2**2
f(666) = 4*dtau**2*(a1*a3 + a2*a4)
f(667) = 4*dtau**4*(a1*a3 + a2*a4)
f(668) = 2*dtau**2*(a1**2*a3 + a1*a3**2 + a2**2*a4 + a2*a4**2)
f(669) = 2*dtau**2*(a1**3*a3 + a1*a3**3 + a2**3*a4 + a2*a4**3)
f(670) = 4*dtau**2*(a1**2*a3**2 + a2**2*a4**2)
f(671) = 2*a1*a4*b1**2 + 2*a1*a4*b2**2 + 2*a2*a3*b1**2 + 2*a2*a3*b2**2
f(672) = 2*a1*a4*b1**4 + 2*a1*a4*b2**4 + 2*a2*a3*b1**4 + 2*a2*a3*b2**4
f(673) = 2*a1**2*a4*b2**2 + 2*a1*a4**2*b1**2 + 2*a2**2*a3*b2**2 + 2*a2* &
      a3**2*b1**2
f(674) = 2*a1**3*a4*b2**2 + 2*a1*a4**3*b1**2 + 2*a2**3*a3*b2**2 + 2*a2* &
      a3**3*b1**2
f(675) = 2*a1**2*a4*b1**2 + 2*a1*a4**2*b2**2 + 2*a2**2*a3*b1**2 + 2*a2* &
      a3**2*b2**2
f(676) = 2*a1**2*a4**2*b1**2 + 2*a1**2*a4**2*b2**2 + 2*a2**2*a3**2*b1**2 &
      + 2*a2**2*a3**2*b2**2
f(677) = 2*a1**3*a4*b1**2 + 2*a1*a4**3*b2**2 + 2*a2**3*a3*b1**2 + 2*a2* &
      a3**3*b2**2
f(678) = 4*dtau**2*(a1*a4 + a2*a3)
f(679) = 4*dtau**4*(a1*a4 + a2*a3)
f(680) = 2*dtau**2*(a1**2*a4 + a1*a4**2 + a2**2*a3 + a2*a3**2)
f(681) = 2*dtau**2*(a1**3*a4 + a1*a4**3 + a2**3*a3 + a2*a3**3)
f(682) = 4*dtau**2*(a1**2*a4**2 + a2**2*a3**2)
f(683) = 2*b1*b2*(a1 + a2 + a3 + a4)
f(684) = 2*b1*b2*(a1*b2**2 + a2*b2**2 + a3*b1**2 + a4*b1**2)
f(685) = 2*b1**2*b2**2*(a1 + a2 + a3 + a4)
f(686) = 2*b1*b2*(a1*b1**2 + a2*b1**2 + a3*b2**2 + a4*b2**2)
f(687) = 2*b1*b2*(a1**2 + a2**2 + a3**2 + a4**2)
f(688) = 2*b1*b2*(a1**2*b2**2 + a2**2*b2**2 + a3**2*b1**2 + a4**2*b1**2)
f(689) = 2*b1**2*b2**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(690) = 2*b1*b2*(a1**2*b1**2 + a2**2*b1**2 + a3**2*b2**2 + a4**2*b2**2)
f(691) = 2*b1*b2*(a1**3 + a2**3 + a3**3 + a4**3)
f(692) = 2*b1*b2*(a1**4 + a2**4 + a3**4 + a4**4)
f(693) = 2*dtau*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(694) = 2*dtau**3*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(695) = 2*dtau**2*(a1*b1**2 + a2*b1**2 + a3*b2**2 + a4*b2**2)
f(696) = 2*dtau*(a1*b1**3 - a2*b1**3 + a3*b2**3 - a4*b2**3)
f(697) = 2*dtau*(a1**2*b1 - a2**2*b1 + a3**2*b2 - a4**2*b2)
f(698) = 2*dtau**3*(a1**2*b1 - a2**2*b1 + a3**2*b2 - a4**2*b2)
f(699) = 2*dtau**2*(a1**2*b1**2 + a2**2*b1**2 + a3**2*b2**2 + a4**2*b2** &
      2)
f(700) = 2*dtau*(a1**2*b1**3 - a2**2*b1**3 + a3**2*b2**3 - a4**2*b2**3)
f(701) = 2*dtau*(a1**3*b1 - a2**3*b1 + a3**3*b2 - a4**3*b2)
f(702) = 2*dtau*(a1**4*b1 - a2**4*b1 + a3**4*b2 - a4**4*b2)
f(703) = 2*dtau*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(704) = 2*dtau**3*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(705) = 2*dtau**2*(a1*b2**2 + a2*b2**2 + a3*b1**2 + a4*b1**2)
f(706) = 2*dtau*(a1*b2**3 - a2*b2**3 + a3*b1**3 - a4*b1**3)
f(707) = 2*dtau*(a1**2*b2 - a2**2*b2 + a3**2*b1 - a4**2*b1)
f(708) = 2*dtau**3*(a1**2*b2 - a2**2*b2 + a3**2*b1 - a4**2*b1)
f(709) = 2*dtau**2*(a1**2*b2**2 + a2**2*b2**2 + a3**2*b1**2 + a4**2*b1** &
      2)
f(710) = 2*dtau*(a1**2*b2**3 - a2**2*b2**3 + a3**2*b1**3 - a4**2*b1**3)
f(711) = 2*dtau*(a1**3*b2 - a2**3*b2 + a3**3*b1 - a4**3*b1)
f(712) = 2*dtau*(a1**4*b2 - a2**4*b2 + a3**4*b1 - a4**4*b1)
f(713) = 8*b1*b2*dtau**2
f(714) = 8*b1*b2*dtau**4
f(715) = 4*b1*b2*dtau**2*(b1**2 + b2**2)
f(716) = 8*b1**2*b2**2*dtau**2
v = sum(f*params)
end function c2h4_poten_n3_d6_ADF


!###############################################################################


! D2h-symmetry-adapted 4-mode 6-order expansion for C2H4 molecule

function c2h4_poten_n4_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(1018)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(1018)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 2*r0*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(2) = 2*r0*(r1**2*r3*r4 + r1*r2*r3**2 + r1*r2*r4**2 + r2**2*r3*r4)
f(3) = 2*r0*(r1**3*r3*r4 + r1*r2*r3**3 + r1*r2*r4**3 + r2**3*r3*r4)
f(4) = 2*r0*(r1**2*r2*r4 + r1*r2**2*r3 + r1*r3*r4**2 + r2*r3**2*r4)
f(5) = 2*r0*(r1**2*r2*r4**2 + r1**2*r3*r4**2 + r1*r2**2*r3**2 + r2**2*r3 &
      **2*r4)
f(6) = 2*r0*(r1**3*r2*r4 + r1*r2**3*r3 + r1*r3*r4**3 + r2*r3**3*r4)
f(7) = 2*r0*(r1**2*r2*r3 + r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4**2)
f(8) = 2*r0*(r1**2*r2*r3**2 + r1**2*r3**2*r4 + r1*r2**2*r4**2 + r2**2*r3 &
      *r4**2)
f(9) = 2*r0*(r1**2*r2**2*r3 + r1**2*r2**2*r4 + r1*r3**2*r4**2 + r2*r3**2 &
      *r4**2)
f(10) = 2*r0*(r1**3*r2*r3 + r1*r2**3*r4 + r1*r3**3*r4 + r2*r3*r4**3)
f(11) = 2*r0**2*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(12) = 2*r0**2*(r1**2*r3*r4 + r1*r2*r3**2 + r1*r2*r4**2 + r2**2*r3*r4)
f(13) = 2*r0**2*(r1**2*r2*r4 + r1*r2**2*r3 + r1*r3*r4**2 + r2*r3**2*r4)
f(14) = 2*r0**2*(r1**2*r2*r3 + r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4**2)
f(15) = 2*r0**3*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(16) = 2*r0*(a1*r1*r2 + a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(17) = 2*r0*(a1**2*r1*r2 + a2**2*r1*r2 + a3**2*r3*r4 + a4**2*r3*r4)
f(18) = 2*r0*(a1**3*r1*r2 + a2**3*r1*r2 + a3**3*r3*r4 + a4**3*r3*r4)
f(19) = 2*r0*(a1*r1*r2**2 + a2*r1**2*r2 + a3*r3*r4**2 + a4*r3**2*r4)
f(20) = 2*r0*(a1**2*r1*r2**2 + a2**2*r1**2*r2 + a3**2*r3*r4**2 + a4**2* &
      r3**2*r4)
f(21) = 2*r0*(a1*r1*r2**3 + a2*r1**3*r2 + a3*r3*r4**3 + a4*r3**3*r4)
f(22) = 2*r0*(a1*r1**2*r2 + a2*r1*r2**2 + a3*r3**2*r4 + a4*r3*r4**2)
f(23) = 2*r0*(a1**2*r1**2*r2 + a2**2*r1*r2**2 + a3**2*r3**2*r4 + a4**2* &
      r3*r4**2)
f(24) = 2*r0*(a1*r1**2*r2**2 + a2*r1**2*r2**2 + a3*r3**2*r4**2 + a4*r3** &
      2*r4**2)
f(25) = 2*r0*(a1*r1**3*r2 + a2*r1*r2**3 + a3*r3**3*r4 + a4*r3*r4**3)
f(26) = 2*r0**2*(a1*r1*r2 + a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(27) = 2*r0**2*(a1**2*r1*r2 + a2**2*r1*r2 + a3**2*r3*r4 + a4**2*r3*r4)
f(28) = 2*r0**2*(a1*r1*r2**2 + a2*r1**2*r2 + a3*r3*r4**2 + a4*r3**2*r4)
f(29) = 2*r0**2*(a1*r1**2*r2 + a2*r1*r2**2 + a3*r3**2*r4 + a4*r3*r4**2)
f(30) = 2*r0**3*(a1*r1*r2 + a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(31) = 2*r0*(a1*r3*r4 + a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(32) = 2*r0*(a1**2*r3*r4 + a2**2*r3*r4 + a3**2*r1*r2 + a4**2*r1*r2)
f(33) = 2*r0*(a1**3*r3*r4 + a2**3*r3*r4 + a3**3*r1*r2 + a4**3*r1*r2)
f(34) = 2*r0*(a1*r3*r4**2 + a2*r3**2*r4 + a3*r1*r2**2 + a4*r1**2*r2)
f(35) = 2*r0*(a1**2*r3*r4**2 + a2**2*r3**2*r4 + a3**2*r1*r2**2 + a4**2* &
      r1**2*r2)
f(36) = 2*r0*(a1*r3*r4**3 + a2*r3**3*r4 + a3*r1*r2**3 + a4*r1**3*r2)
f(37) = 2*r0*(a1*r3**2*r4 + a2*r3*r4**2 + a3*r1**2*r2 + a4*r1*r2**2)
f(38) = 2*r0*(a1**2*r3**2*r4 + a2**2*r3*r4**2 + a3**2*r1**2*r2 + a4**2* &
      r1*r2**2)
f(39) = 2*r0*(a1*r3**2*r4**2 + a2*r3**2*r4**2 + a3*r1**2*r2**2 + a4*r1** &
      2*r2**2)
f(40) = 2*r0*(a1*r3**3*r4 + a2*r3*r4**3 + a3*r1**3*r2 + a4*r1*r2**3)
f(41) = 2*r0**2*(a1*r3*r4 + a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(42) = 2*r0**2*(a1**2*r3*r4 + a2**2*r3*r4 + a3**2*r1*r2 + a4**2*r1*r2)
f(43) = 2*r0**2*(a1*r3*r4**2 + a2*r3**2*r4 + a3*r1*r2**2 + a4*r1**2*r2)
f(44) = 2*r0**2*(a1*r3**2*r4 + a2*r3*r4**2 + a3*r1**2*r2 + a4*r1*r2**2)
f(45) = 2*r0**3*(a1*r3*r4 + a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(46) = 4*r0*(b1**2*r1*r2 + b2**2*r3*r4)
f(47) = 2*r0*(b1**2*r1**2*r2 + b1**2*r1*r2**2 + b2**2*r3**2*r4 + b2**2* &
      r3*r4**2)
f(48) = 4*r0**2*(b1**2*r1*r2 + b2**2*r3*r4)
f(49) = 4*r0*(b1**2*r3*r4 + b2**2*r1*r2)
f(50) = 2*r0*(b1**2*r3**2*r4 + b1**2*r3*r4**2 + b2**2*r1**2*r2 + b2**2* &
      r1*r2**2)
f(51) = 4*r0**2*(b1**2*r3*r4 + b2**2*r1*r2)
f(52) = 4*dtau**2*r0*(r1*r2 + r3*r4)
f(53) = 2*dtau**2*r0*(r1**2*r2 + r1*r2**2 + r3**2*r4 + r3*r4**2)
f(54) = 4*dtau**2*r0**2*(r1*r2 + r3*r4)
f(55) = 2*r0*(a1*r1*r3 + a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(56) = 2*r0*(a1**2*r1*r3 + a2**2*r2*r4 + a3**2*r1*r3 + a4**2*r2*r4)
f(57) = 2*r0*(a1**3*r1*r3 + a2**3*r2*r4 + a3**3*r1*r3 + a4**3*r2*r4)
f(58) = 2*r0*(a1*r1*r3**2 + a2*r2*r4**2 + a3*r1**2*r3 + a4*r2**2*r4)
f(59) = 2*r0*(a1**2*r1*r3**2 + a2**2*r2*r4**2 + a3**2*r1**2*r3 + a4**2* &
      r2**2*r4)
f(60) = 2*r0*(a1*r1*r3**3 + a2*r2*r4**3 + a3*r1**3*r3 + a4*r2**3*r4)
f(61) = 2*r0*(a1*r1**2*r3 + a2*r2**2*r4 + a3*r1*r3**2 + a4*r2*r4**2)
f(62) = 2*r0*(a1**2*r1**2*r3 + a2**2*r2**2*r4 + a3**2*r1*r3**2 + a4**2* &
      r2*r4**2)
f(63) = 2*r0*(a1*r1**2*r3**2 + a2*r2**2*r4**2 + a3*r1**2*r3**2 + a4*r2** &
      2*r4**2)
f(64) = 2*r0*(a1*r1**3*r3 + a2*r2**3*r4 + a3*r1*r3**3 + a4*r2*r4**3)
f(65) = 2*r0**2*(a1*r1*r3 + a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(66) = 2*r0**2*(a1**2*r1*r3 + a2**2*r2*r4 + a3**2*r1*r3 + a4**2*r2*r4)
f(67) = 2*r0**2*(a1*r1*r3**2 + a2*r2*r4**2 + a3*r1**2*r3 + a4*r2**2*r4)
f(68) = 2*r0**2*(a1*r1**2*r3 + a2*r2**2*r4 + a3*r1*r3**2 + a4*r2*r4**2)
f(69) = 2*r0**3*(a1*r1*r3 + a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(70) = 2*r0*(a1*r2*r4 + a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(71) = 2*r0*(a1**2*r2*r4 + a2**2*r1*r3 + a3**2*r2*r4 + a4**2*r1*r3)
f(72) = 2*r0*(a1**3*r2*r4 + a2**3*r1*r3 + a3**3*r2*r4 + a4**3*r1*r3)
f(73) = 2*r0*(a1*r2*r4**2 + a2*r1*r3**2 + a3*r2**2*r4 + a4*r1**2*r3)
f(74) = 2*r0*(a1**2*r2*r4**2 + a2**2*r1*r3**2 + a3**2*r2**2*r4 + a4**2* &
      r1**2*r3)
f(75) = 2*r0*(a1*r2*r4**3 + a2*r1*r3**3 + a3*r2**3*r4 + a4*r1**3*r3)
f(76) = 2*r0*(a1*r2**2*r4 + a2*r1**2*r3 + a3*r2*r4**2 + a4*r1*r3**2)
f(77) = 2*r0*(a1**2*r2**2*r4 + a2**2*r1**2*r3 + a3**2*r2*r4**2 + a4**2* &
      r1*r3**2)
f(78) = 2*r0*(a1*r2**2*r4**2 + a2*r1**2*r3**2 + a3*r2**2*r4**2 + a4*r1** &
      2*r3**2)
f(79) = 2*r0*(a1*r2**3*r4 + a2*r1**3*r3 + a3*r2*r4**3 + a4*r1*r3**3)
f(80) = 2*r0**2*(a1*r2*r4 + a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(81) = 2*r0**2*(a1**2*r2*r4 + a2**2*r1*r3 + a3**2*r2*r4 + a4**2*r1*r3)
f(82) = 2*r0**2*(a1*r2*r4**2 + a2*r1*r3**2 + a3*r2**2*r4 + a4*r1**2*r3)
f(83) = 2*r0**2*(a1*r2**2*r4 + a2*r1**2*r3 + a3*r2*r4**2 + a4*r1*r3**2)
f(84) = 2*r0**3*(a1*r2*r4 + a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(85) = 2*r0*(b1**2*r1*r3 + b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2*r4)
f(86) = 2*r0*(b1**2*r1*r3**2 + b1**2*r2*r4**2 + b2**2*r1**2*r3 + b2**2* &
      r2**2*r4)
f(87) = 2*r0*(b1**2*r1**2*r3 + b1**2*r2**2*r4 + b2**2*r1*r3**2 + b2**2* &
      r2*r4**2)
f(88) = 2*r0**2*(b1**2*r1*r3 + b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2*r4)
f(89) = 4*dtau**2*r0*(r1*r3 + r2*r4)
f(90) = 2*dtau**2*r0*(r1**2*r3 + r1*r3**2 + r2**2*r4 + r2*r4**2)
f(91) = 4*dtau**2*r0**2*(r1*r3 + r2*r4)
f(92) = 2*r0*(a1*r1*r4 + a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(93) = 2*r0*(a1**2*r1*r4 + a2**2*r2*r3 + a3**2*r2*r3 + a4**2*r1*r4)
f(94) = 2*r0*(a1**3*r1*r4 + a2**3*r2*r3 + a3**3*r2*r3 + a4**3*r1*r4)
f(95) = 2*r0*(a1*r1*r4**2 + a2*r2*r3**2 + a3*r2**2*r3 + a4*r1**2*r4)
f(96) = 2*r0*(a1**2*r1*r4**2 + a2**2*r2*r3**2 + a3**2*r2**2*r3 + a4**2* &
      r1**2*r4)
f(97) = 2*r0*(a1*r1*r4**3 + a2*r2*r3**3 + a3*r2**3*r3 + a4*r1**3*r4)
f(98) = 2*r0*(a1*r1**2*r4 + a2*r2**2*r3 + a3*r2*r3**2 + a4*r1*r4**2)
f(99) = 2*r0*(a1**2*r1**2*r4 + a2**2*r2**2*r3 + a3**2*r2*r3**2 + a4**2* &
      r1*r4**2)
f(100) = 2*r0*(a1*r1**2*r4**2 + a2*r2**2*r3**2 + a3*r2**2*r3**2 + a4*r1 &
      **2*r4**2)
f(101) = 2*r0*(a1*r1**3*r4 + a2*r2**3*r3 + a3*r2*r3**3 + a4*r1*r4**3)
f(102) = 2*r0**2*(a1*r1*r4 + a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(103) = 2*r0**2*(a1**2*r1*r4 + a2**2*r2*r3 + a3**2*r2*r3 + a4**2*r1*r4)
f(104) = 2*r0**2*(a1*r1*r4**2 + a2*r2*r3**2 + a3*r2**2*r3 + a4*r1**2*r4)
f(105) = 2*r0**2*(a1*r1**2*r4 + a2*r2**2*r3 + a3*r2*r3**2 + a4*r1*r4**2)
f(106) = 2*r0**3*(a1*r1*r4 + a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(107) = 2*r0*(a1*r2*r3 + a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(108) = 2*r0*(a1**2*r2*r3 + a2**2*r1*r4 + a3**2*r1*r4 + a4**2*r2*r3)
f(109) = 2*r0*(a1**3*r2*r3 + a2**3*r1*r4 + a3**3*r1*r4 + a4**3*r2*r3)
f(110) = 2*r0*(a1*r2*r3**2 + a2*r1*r4**2 + a3*r1**2*r4 + a4*r2**2*r3)
f(111) = 2*r0*(a1**2*r2*r3**2 + a2**2*r1*r4**2 + a3**2*r1**2*r4 + a4**2* &
      r2**2*r3)
f(112) = 2*r0*(a1*r2*r3**3 + a2*r1*r4**3 + a3*r1**3*r4 + a4*r2**3*r3)
f(113) = 2*r0*(a1*r2**2*r3 + a2*r1**2*r4 + a3*r1*r4**2 + a4*r2*r3**2)
f(114) = 2*r0*(a1**2*r2**2*r3 + a2**2*r1**2*r4 + a3**2*r1*r4**2 + a4**2* &
      r2*r3**2)
f(115) = 2*r0*(a1*r2**2*r3**2 + a2*r1**2*r4**2 + a3*r1**2*r4**2 + a4*r2 &
      **2*r3**2)
f(116) = 2*r0*(a1*r2**3*r3 + a2*r1**3*r4 + a3*r1*r4**3 + a4*r2*r3**3)
f(117) = 2*r0**2*(a1*r2*r3 + a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(118) = 2*r0**2*(a1**2*r2*r3 + a2**2*r1*r4 + a3**2*r1*r4 + a4**2*r2*r3)
f(119) = 2*r0**2*(a1*r2*r3**2 + a2*r1*r4**2 + a3*r1**2*r4 + a4*r2**2*r3)
f(120) = 2*r0**2*(a1*r2**2*r3 + a2*r1**2*r4 + a3*r1*r4**2 + a4*r2*r3**2)
f(121) = 2*r0**3*(a1*r2*r3 + a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(122) = 2*r0*(b1**2*r1*r4 + b1**2*r2*r3 + b2**2*r1*r4 + b2**2*r2*r3)
f(123) = 2*r0*(b1**2*r1*r4**2 + b1**2*r2*r3**2 + b2**2*r1**2*r4 + b2**2* &
      r2**2*r3)
f(124) = 2*r0*(b1**2*r1**2*r4 + b1**2*r2**2*r3 + b2**2*r1*r4**2 + b2**2* &
      r2*r3**2)
f(125) = 2*r0**2*(b1**2*r1*r4 + b1**2*r2*r3 + b2**2*r1*r4 + b2**2*r2*r3)
f(126) = 4*dtau**2*r0*(r1*r4 + r2*r3)
f(127) = 2*dtau**2*r0*(r1**2*r4 + r1*r4**2 + r2**2*r3 + r2*r3**2)
f(128) = 4*dtau**2*r0**2*(r1*r4 + r2*r3)
f(129) = 2*r0*(a1*a2*r1 + a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(130) = 2*r0*(a1**2*a2*r2 + a1*a2**2*r1 + a3**2*a4*r4 + a3*a4**2*r3)
f(131) = 2*r0*(a1**3*a2*r2 + a1*a2**3*r1 + a3**3*a4*r4 + a3*a4**3*r3)
f(132) = 2*r0*(a1**2*a2*r1 + a1*a2**2*r2 + a3**2*a4*r3 + a3*a4**2*r4)
f(133) = 2*r0*(a1**2*a2**2*r1 + a1**2*a2**2*r2 + a3**2*a4**2*r3 + a3**2* &
      a4**2*r4)
f(134) = 2*r0*(a1**3*a2*r1 + a1*a2**3*r2 + a3**3*a4*r3 + a3*a4**3*r4)
f(135) = 2*r0*(a1*a2*r1**2 + a1*a2*r2**2 + a3*a4*r3**2 + a3*a4*r4**2)
f(136) = 2*r0*(a1**2*a2*r2**2 + a1*a2**2*r1**2 + a3**2*a4*r4**2 + a3*a4 &
      **2*r3**2)
f(137) = 2*r0*(a1**2*a2*r1**2 + a1*a2**2*r2**2 + a3**2*a4*r3**2 + a3*a4 &
      **2*r4**2)
f(138) = 2*r0*(a1*a2*r1**3 + a1*a2*r2**3 + a3*a4*r3**3 + a3*a4*r4**3)
f(139) = 2*r0**2*(a1*a2*r1 + a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(140) = 2*r0**2*(a1**2*a2*r2 + a1*a2**2*r1 + a3**2*a4*r4 + a3*a4**2*r3)
f(141) = 2*r0**2*(a1**2*a2*r1 + a1*a2**2*r2 + a3**2*a4*r3 + a3*a4**2*r4)
f(142) = 2*r0**2*(a1*a2*r1**2 + a1*a2*r2**2 + a3*a4*r3**2 + a3*a4*r4**2)
f(143) = 2*r0**3*(a1*a2*r1 + a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(144) = 2*r0*(a1*a3*r1 + a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(145) = 2*r0*(a1**2*a3*r3 + a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2*r2)
f(146) = 2*r0*(a1**3*a3*r3 + a1*a3**3*r1 + a2**3*a4*r4 + a2*a4**3*r2)
f(147) = 2*r0*(a1**2*a3*r1 + a1*a3**2*r3 + a2**2*a4*r2 + a2*a4**2*r4)
f(148) = 2*r0*(a1**2*a3**2*r1 + a1**2*a3**2*r3 + a2**2*a4**2*r2 + a2**2* &
      a4**2*r4)
f(149) = 2*r0*(a1**3*a3*r1 + a1*a3**3*r3 + a2**3*a4*r2 + a2*a4**3*r4)
f(150) = 2*r0*(a1*a3*r1**2 + a1*a3*r3**2 + a2*a4*r2**2 + a2*a4*r4**2)
f(151) = 2*r0*(a1**2*a3*r3**2 + a1*a3**2*r1**2 + a2**2*a4*r4**2 + a2*a4 &
      **2*r2**2)
f(152) = 2*r0*(a1**2*a3*r1**2 + a1*a3**2*r3**2 + a2**2*a4*r2**2 + a2*a4 &
      **2*r4**2)
f(153) = 2*r0*(a1*a3*r1**3 + a1*a3*r3**3 + a2*a4*r2**3 + a2*a4*r4**3)
f(154) = 2*r0**2*(a1*a3*r1 + a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(155) = 2*r0**2*(a1**2*a3*r3 + a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2*r2)
f(156) = 2*r0**2*(a1**2*a3*r1 + a1*a3**2*r3 + a2**2*a4*r2 + a2*a4**2*r4)
f(157) = 2*r0**2*(a1*a3*r1**2 + a1*a3*r3**2 + a2*a4*r2**2 + a2*a4*r4**2)
f(158) = 2*r0**3*(a1*a3*r1 + a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(159) = 2*r0*(a1*a4*r1 + a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(160) = 2*r0*(a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 + a2*a3**2*r2)
f(161) = 2*r0*(a1**3*a4*r4 + a1*a4**3*r1 + a2**3*a3*r3 + a2*a3**3*r2)
f(162) = 2*r0*(a1**2*a4*r1 + a1*a4**2*r4 + a2**2*a3*r2 + a2*a3**2*r3)
f(163) = 2*r0*(a1**2*a4**2*r1 + a1**2*a4**2*r4 + a2**2*a3**2*r2 + a2**2* &
      a3**2*r3)
f(164) = 2*r0*(a1**3*a4*r1 + a1*a4**3*r4 + a2**3*a3*r2 + a2*a3**3*r3)
f(165) = 2*r0*(a1*a4*r1**2 + a1*a4*r4**2 + a2*a3*r2**2 + a2*a3*r3**2)
f(166) = 2*r0*(a1**2*a4*r4**2 + a1*a4**2*r1**2 + a2**2*a3*r3**2 + a2*a3 &
      **2*r2**2)
f(167) = 2*r0*(a1**2*a4*r1**2 + a1*a4**2*r4**2 + a2**2*a3*r2**2 + a2*a3 &
      **2*r3**2)
f(168) = 2*r0*(a1*a4*r1**3 + a1*a4*r4**3 + a2*a3*r2**3 + a2*a3*r3**3)
f(169) = 2*r0**2*(a1*a4*r1 + a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(170) = 2*r0**2*(a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 + a2*a3**2*r2)
f(171) = 2*r0**2*(a1**2*a4*r1 + a1*a4**2*r4 + a2**2*a3*r2 + a2*a3**2*r3)
f(172) = 2*r0**2*(a1*a4*r1**2 + a1*a4*r4**2 + a2*a3*r2**2 + a2*a3*r3**2)
f(173) = 2*r0**3*(a1*a4*r1 + a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(174) = 2*r0*(a1*b1**2*r1 + a2*b1**2*r2 + a3*b2**2*r3 + a4*b2**2*r4)
f(175) = 2*r0*(a1**2*b1**2*r1 + a2**2*b1**2*r2 + a3**2*b2**2*r3 + a4**2* &
      b2**2*r4)
f(176) = 2*r0*(a1*b1**2*r1**2 + a2*b1**2*r2**2 + a3*b2**2*r3**2 + a4*b2 &
      **2*r4**2)
f(177) = 2*r0**2*(a1*b1**2*r1 + a2*b1**2*r2 + a3*b2**2*r3 + a4*b2**2*r4)
f(178) = 2*r0*(a1*b2**2*r1 + a2*b2**2*r2 + a3*b1**2*r3 + a4*b1**2*r4)
f(179) = 2*r0*(a1**2*b2**2*r1 + a2**2*b2**2*r2 + a3**2*b1**2*r3 + a4**2* &
      b1**2*r4)
f(180) = 2*r0*(a1*b2**2*r1**2 + a2*b2**2*r2**2 + a3*b1**2*r3**2 + a4*b1 &
      **2*r4**2)
f(181) = 2*r0**2*(a1*b2**2*r1 + a2*b2**2*r2 + a3*b1**2*r3 + a4*b1**2*r4)
f(182) = 2*dtau**2*r0*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(183) = 2*dtau**2*r0*(a1**2*r1 + a2**2*r2 + a3**2*r3 + a4**2*r4)
f(184) = 2*dtau**2*r0*(a1*r1**2 + a2*r2**2 + a3*r3**2 + a4*r4**2)
f(185) = 2*dtau**2*r0**2*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(186) = 2*r0*(a1*a4*r2 + a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(187) = 2*r0*(a1**2*a4*r3 + a1*a4**2*r2 + a2**2*a3*r4 + a2*a3**2*r1)
f(188) = 2*r0*(a1**3*a4*r3 + a1*a4**3*r2 + a2**3*a3*r4 + a2*a3**3*r1)
f(189) = 2*r0*(a1**2*a4*r2 + a1*a4**2*r3 + a2**2*a3*r1 + a2*a3**2*r4)
f(190) = 2*r0*(a1**2*a4**2*r2 + a1**2*a4**2*r3 + a2**2*a3**2*r1 + a2**2* &
      a3**2*r4)
f(191) = 2*r0*(a1**3*a4*r2 + a1*a4**3*r3 + a2**3*a3*r1 + a2*a3**3*r4)
f(192) = 2*r0*(a1*a4*r2**2 + a1*a4*r3**2 + a2*a3*r1**2 + a2*a3*r4**2)
f(193) = 2*r0*(a1**2*a4*r3**2 + a1*a4**2*r2**2 + a2**2*a3*r4**2 + a2*a3 &
      **2*r1**2)
f(194) = 2*r0*(a1**2*a4*r2**2 + a1*a4**2*r3**2 + a2**2*a3*r1**2 + a2*a3 &
      **2*r4**2)
f(195) = 2*r0*(a1*a4*r2**3 + a1*a4*r3**3 + a2*a3*r1**3 + a2*a3*r4**3)
f(196) = 2*r0**2*(a1*a4*r2 + a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(197) = 2*r0**2*(a1**2*a4*r3 + a1*a4**2*r2 + a2**2*a3*r4 + a2*a3**2*r1)
f(198) = 2*r0**2*(a1**2*a4*r2 + a1*a4**2*r3 + a2**2*a3*r1 + a2*a3**2*r4)
f(199) = 2*r0**2*(a1*a4*r2**2 + a1*a4*r3**2 + a2*a3*r1**2 + a2*a3*r4**2)
f(200) = 2*r0**3*(a1*a4*r2 + a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(201) = 2*r0*(a1*a3*r2 + a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(202) = 2*r0*(a1**2*a3*r4 + a1*a3**2*r2 + a2**2*a4*r3 + a2*a4**2*r1)
f(203) = 2*r0*(a1**3*a3*r4 + a1*a3**3*r2 + a2**3*a4*r3 + a2*a4**3*r1)
f(204) = 2*r0*(a1**2*a3*r2 + a1*a3**2*r4 + a2**2*a4*r1 + a2*a4**2*r3)
f(205) = 2*r0*(a1**2*a3**2*r2 + a1**2*a3**2*r4 + a2**2*a4**2*r1 + a2**2* &
      a4**2*r3)
f(206) = 2*r0*(a1**3*a3*r2 + a1*a3**3*r4 + a2**3*a4*r1 + a2*a4**3*r3)
f(207) = 2*r0*(a1*a3*r2**2 + a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3**2)
f(208) = 2*r0*(a1**2*a3*r4**2 + a1*a3**2*r2**2 + a2**2*a4*r3**2 + a2*a4 &
      **2*r1**2)
f(209) = 2*r0*(a1**2*a3*r2**2 + a1*a3**2*r4**2 + a2**2*a4*r1**2 + a2*a4 &
      **2*r3**2)
f(210) = 2*r0*(a1*a3*r2**3 + a1*a3*r4**3 + a2*a4*r1**3 + a2*a4*r3**3)
f(211) = 2*r0**2*(a1*a3*r2 + a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(212) = 2*r0**2*(a1**2*a3*r4 + a1*a3**2*r2 + a2**2*a4*r3 + a2*a4**2*r1)
f(213) = 2*r0**2*(a1**2*a3*r2 + a1*a3**2*r4 + a2**2*a4*r1 + a2*a4**2*r3)
f(214) = 2*r0**2*(a1*a3*r2**2 + a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3**2)
f(215) = 2*r0**3*(a1*a3*r2 + a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(216) = 2*r0*(a1*b1**2*r2 + a2*b1**2*r1 + a3*b2**2*r4 + a4*b2**2*r3)
f(217) = 2*r0*(a1**2*b1**2*r2 + a2**2*b1**2*r1 + a3**2*b2**2*r4 + a4**2* &
      b2**2*r3)
f(218) = 2*r0*(a1*b1**2*r2**2 + a2*b1**2*r1**2 + a3*b2**2*r4**2 + a4*b2 &
      **2*r3**2)
f(219) = 2*r0**2*(a1*b1**2*r2 + a2*b1**2*r1 + a3*b2**2*r4 + a4*b2**2*r3)
f(220) = 2*r0*(a1*b2**2*r2 + a2*b2**2*r1 + a3*b1**2*r4 + a4*b1**2*r3)
f(221) = 2*r0*(a1**2*b2**2*r2 + a2**2*b2**2*r1 + a3**2*b1**2*r4 + a4**2* &
      b1**2*r3)
f(222) = 2*r0*(a1*b2**2*r2**2 + a2*b2**2*r1**2 + a3*b1**2*r4**2 + a4*b1 &
      **2*r3**2)
f(223) = 2*r0**2*(a1*b2**2*r2 + a2*b2**2*r1 + a3*b1**2*r4 + a4*b1**2*r3)
f(224) = 2*dtau**2*r0*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(225) = 2*dtau**2*r0*(a1**2*r2 + a2**2*r1 + a3**2*r4 + a4**2*r3)
f(226) = 2*dtau**2*r0*(a1*r2**2 + a2*r1**2 + a3*r4**2 + a4*r3**2)
f(227) = 2*dtau**2*r0**2*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(228) = 2*r0*(a1*a2*r3 + a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(229) = 2*r0*(a1**2*a2*r4 + a1*a2**2*r3 + a3**2*a4*r2 + a3*a4**2*r1)
f(230) = 2*r0*(a1**3*a2*r4 + a1*a2**3*r3 + a3**3*a4*r2 + a3*a4**3*r1)
f(231) = 2*r0*(a1**2*a2*r3 + a1*a2**2*r4 + a3**2*a4*r1 + a3*a4**2*r2)
f(232) = 2*r0*(a1**2*a2**2*r3 + a1**2*a2**2*r4 + a3**2*a4**2*r1 + a3**2* &
      a4**2*r2)
f(233) = 2*r0*(a1**3*a2*r3 + a1*a2**3*r4 + a3**3*a4*r1 + a3*a4**3*r2)
f(234) = 2*r0*(a1*a2*r3**2 + a1*a2*r4**2 + a3*a4*r1**2 + a3*a4*r2**2)
f(235) = 2*r0*(a1**2*a2*r4**2 + a1*a2**2*r3**2 + a3**2*a4*r2**2 + a3*a4 &
      **2*r1**2)
f(236) = 2*r0*(a1**2*a2*r3**2 + a1*a2**2*r4**2 + a3**2*a4*r1**2 + a3*a4 &
      **2*r2**2)
f(237) = 2*r0*(a1*a2*r3**3 + a1*a2*r4**3 + a3*a4*r1**3 + a3*a4*r2**3)
f(238) = 2*r0**2*(a1*a2*r3 + a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(239) = 2*r0**2*(a1**2*a2*r4 + a1*a2**2*r3 + a3**2*a4*r2 + a3*a4**2*r1)
f(240) = 2*r0**2*(a1**2*a2*r3 + a1*a2**2*r4 + a3**2*a4*r1 + a3*a4**2*r2)
f(241) = 2*r0**2*(a1*a2*r3**2 + a1*a2*r4**2 + a3*a4*r1**2 + a3*a4*r2**2)
f(242) = 2*r0**3*(a1*a2*r3 + a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(243) = 2*r0*(a1*b2**2*r3 + a2*b2**2*r4 + a3*b1**2*r1 + a4*b1**2*r2)
f(244) = 2*r0*(a1**2*b2**2*r3 + a2**2*b2**2*r4 + a3**2*b1**2*r1 + a4**2* &
      b1**2*r2)
f(245) = 2*r0*(a1*b2**2*r3**2 + a2*b2**2*r4**2 + a3*b1**2*r1**2 + a4*b1 &
      **2*r2**2)
f(246) = 2*r0**2*(a1*b2**2*r3 + a2*b2**2*r4 + a3*b1**2*r1 + a4*b1**2*r2)
f(247) = 2*r0*(a1*b1**2*r3 + a2*b1**2*r4 + a3*b2**2*r1 + a4*b2**2*r2)
f(248) = 2*r0*(a1**2*b1**2*r3 + a2**2*b1**2*r4 + a3**2*b2**2*r1 + a4**2* &
      b2**2*r2)
f(249) = 2*r0*(a1*b1**2*r3**2 + a2*b1**2*r4**2 + a3*b2**2*r1**2 + a4*b2 &
      **2*r2**2)
f(250) = 2*r0**2*(a1*b1**2*r3 + a2*b1**2*r4 + a3*b2**2*r1 + a4*b2**2*r2)
f(251) = 2*dtau**2*r0*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(252) = 2*dtau**2*r0*(a1**2*r3 + a2**2*r4 + a3**2*r1 + a4**2*r2)
f(253) = 2*dtau**2*r0*(a1*r3**2 + a2*r4**2 + a3*r1**2 + a4*r2**2)
f(254) = 2*dtau**2*r0**2*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(255) = 2*r0*(a1*b2**2*r4 + a2*b2**2*r3 + a3*b1**2*r2 + a4*b1**2*r1)
f(256) = 2*r0*(a1**2*b2**2*r4 + a2**2*b2**2*r3 + a3**2*b1**2*r2 + a4**2* &
      b1**2*r1)
f(257) = 2*r0*(a1*b2**2*r4**2 + a2*b2**2*r3**2 + a3*b1**2*r2**2 + a4*b1 &
      **2*r1**2)
f(258) = 2*r0**2*(a1*b2**2*r4 + a2*b2**2*r3 + a3*b1**2*r2 + a4*b1**2*r1)
f(259) = 2*r0*(a1*b1**2*r4 + a2*b1**2*r3 + a3*b2**2*r2 + a4*b2**2*r1)
f(260) = 2*r0*(a1**2*b1**2*r4 + a2**2*b1**2*r3 + a3**2*b2**2*r2 + a4**2* &
      b2**2*r1)
f(261) = 2*r0*(a1*b1**2*r4**2 + a2*b1**2*r3**2 + a3*b2**2*r2**2 + a4*b2 &
      **2*r1**2)
f(262) = 2*r0**2*(a1*b1**2*r4 + a2*b1**2*r3 + a3*b2**2*r2 + a4*b2**2*r1)
f(263) = 2*dtau**2*r0*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(264) = 2*dtau**2*r0*(a1**2*r4 + a2**2*r3 + a3**2*r2 + a4**2*r1)
f(265) = 2*dtau**2*r0*(a1*r4**2 + a2*r3**2 + a3*r2**2 + a4*r1**2)
f(266) = 2*dtau**2*r0**2*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(267) = 2*b1*b2*r0*(r1 + r2 + r3 + r4)
f(268) = 2*b1*b2*r0*(b1**2*r3 + b1**2*r4 + b2**2*r1 + b2**2*r2)
f(269) = 2*b1**2*b2**2*r0*(r1 + r2 + r3 + r4)
f(270) = 2*b1*b2*r0*(b1**2*r1 + b1**2*r2 + b2**2*r3 + b2**2*r4)
f(271) = 2*b1*b2*r0*(r1**2 + r2**2 + r3**2 + r4**2)
f(272) = 2*b1*b2*r0*(r1**3 + r2**3 + r3**3 + r4**3)
f(273) = 2*b1*b2*r0**2*(r1 + r2 + r3 + r4)
f(274) = 2*b1*b2*r0**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(275) = 2*b1*b2*r0**3*(r1 + r2 + r3 + r4)
f(276) = 2*dtau*r0*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(277) = 2*dtau**3*r0*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(278) = 2*dtau**2*r0*(b1**2*r1 + b1**2*r2 + b2**2*r3 + b2**2*r4)
f(279) = 2*dtau*r0*(b1**3*r1 - b1**3*r2 + b2**3*r3 - b2**3*r4)
f(280) = 2*dtau*r0*(b1*r1**2 - b1*r2**2 + b2*r3**2 - b2*r4**2)
f(281) = 2*dtau*r0*(b1*r1**3 - b1*r2**3 + b2*r3**3 - b2*r4**3)
f(282) = 2*dtau*r0**2*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(283) = 2*dtau*r0**2*(b1*r1**2 - b1*r2**2 + b2*r3**2 - b2*r4**2)
f(284) = 2*dtau*r0**3*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(285) = 2*dtau*r0*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(286) = 2*dtau**3*r0*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(287) = 2*dtau**2*r0*(b1**2*r3 + b1**2*r4 + b2**2*r1 + b2**2*r2)
f(288) = 2*dtau*r0*(b1**3*r3 - b1**3*r4 + b2**3*r1 - b2**3*r2)
f(289) = 2*dtau*r0*(b1*r3**2 - b1*r4**2 + b2*r1**2 - b2*r2**2)
f(290) = 2*dtau*r0*(b1*r3**3 - b1*r4**3 + b2*r1**3 - b2*r2**3)
f(291) = 2*dtau*r0**2*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(292) = 2*dtau*r0**2*(b1*r3**2 - b1*r4**2 + b2*r1**2 - b2*r2**2)
f(293) = 2*dtau*r0**3*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(294) = 2*r0*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(295) = 2*r0*(a1**2*a3*a4 + a1*a2*a3**2 + a1*a2*a4**2 + a2**2*a3*a4)
f(296) = 2*r0*(a1**3*a3*a4 + a1*a2*a3**3 + a1*a2*a4**3 + a2**3*a3*a4)
f(297) = 2*r0*(a1**2*a2*a4 + a1*a2**2*a3 + a1*a3*a4**2 + a2*a3**2*a4)
f(298) = 2*r0*(a1**2*a2*a4**2 + a1**2*a3*a4**2 + a1*a2**2*a3**2 + a2**2* &
      a3**2*a4)
f(299) = 2*r0*(a1**3*a2*a4 + a1*a2**3*a3 + a1*a3*a4**3 + a2*a3**3*a4)
f(300) = 2*r0*(a1**2*a2*a3 + a1*a2**2*a4 + a1*a3**2*a4 + a2*a3*a4**2)
f(301) = 2*r0*(a1**2*a2*a3**2 + a1**2*a3**2*a4 + a1*a2**2*a4**2 + a2**2* &
      a3*a4**2)
f(302) = 2*r0*(a1**2*a2**2*a3 + a1**2*a2**2*a4 + a1*a3**2*a4**2 + a2*a3 &
      **2*a4**2)
f(303) = 2*r0*(a1**3*a2*a3 + a1*a2**3*a4 + a1*a3**3*a4 + a2*a3*a4**3)
f(304) = 2*r0**2*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(305) = 2*r0**2*(a1**2*a3*a4 + a1*a2*a3**2 + a1*a2*a4**2 + a2**2*a3*a4)
f(306) = 2*r0**2*(a1**2*a2*a4 + a1*a2**2*a3 + a1*a3*a4**2 + a2*a3**2*a4)
f(307) = 2*r0**2*(a1**2*a2*a3 + a1*a2**2*a4 + a1*a3**2*a4 + a2*a3*a4**2)
f(308) = 2*r0**3*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(309) = 4*r0*(a1*a2*b1**2 + a3*a4*b2**2)
f(310) = 2*r0*(a1**2*a2*b1**2 + a1*a2**2*b1**2 + a3**2*a4*b2**2 + a3*a4 &
      **2*b2**2)
f(311) = 4*r0**2*(a1*a2*b1**2 + a3*a4*b2**2)
f(312) = 4*r0*(a1*a2*b2**2 + a3*a4*b1**2)
f(313) = 2*r0*(a1**2*a2*b2**2 + a1*a2**2*b2**2 + a3**2*a4*b1**2 + a3*a4 &
      **2*b1**2)
f(314) = 4*r0**2*(a1*a2*b2**2 + a3*a4*b1**2)
f(315) = 4*dtau**2*r0*(a1*a2 + a3*a4)
f(316) = 2*dtau**2*r0*(a1**2*a2 + a1*a2**2 + a3**2*a4 + a3*a4**2)
f(317) = 4*dtau**2*r0**2*(a1*a2 + a3*a4)
f(318) = 2*r0*(a1*a3*b1**2 + a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2**2)
f(319) = 2*r0*(a1**2*a3*b2**2 + a1*a3**2*b1**2 + a2**2*a4*b2**2 + a2*a4 &
      **2*b1**2)
f(320) = 2*r0*(a1**2*a3*b1**2 + a1*a3**2*b2**2 + a2**2*a4*b1**2 + a2*a4 &
      **2*b2**2)
f(321) = 2*r0**2*(a1*a3*b1**2 + a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2**2)
f(322) = 4*dtau**2*r0*(a1*a3 + a2*a4)
f(323) = 2*dtau**2*r0*(a1**2*a3 + a1*a3**2 + a2**2*a4 + a2*a4**2)
f(324) = 4*dtau**2*r0**2*(a1*a3 + a2*a4)
f(325) = 2*r0*(a1*a4*b1**2 + a1*a4*b2**2 + a2*a3*b1**2 + a2*a3*b2**2)
f(326) = 2*r0*(a1**2*a4*b2**2 + a1*a4**2*b1**2 + a2**2*a3*b2**2 + a2*a3 &
      **2*b1**2)
f(327) = 2*r0*(a1**2*a4*b1**2 + a1*a4**2*b2**2 + a2**2*a3*b1**2 + a2*a3 &
      **2*b2**2)
f(328) = 2*r0**2*(a1*a4*b1**2 + a1*a4*b2**2 + a2*a3*b1**2 + a2*a3*b2**2)
f(329) = 4*dtau**2*r0*(a1*a4 + a2*a3)
f(330) = 2*dtau**2*r0*(a1**2*a4 + a1*a4**2 + a2**2*a3 + a2*a3**2)
f(331) = 4*dtau**2*r0**2*(a1*a4 + a2*a3)
f(332) = 2*b1*b2*r0*(a1 + a2 + a3 + a4)
f(333) = 2*b1*b2*r0*(a1*b2**2 + a2*b2**2 + a3*b1**2 + a4*b1**2)
f(334) = 2*b1**2*b2**2*r0*(a1 + a2 + a3 + a4)
f(335) = 2*b1*b2*r0*(a1*b1**2 + a2*b1**2 + a3*b2**2 + a4*b2**2)
f(336) = 2*b1*b2*r0*(a1**2 + a2**2 + a3**2 + a4**2)
f(337) = 2*b1*b2*r0*(a1**3 + a2**3 + a3**3 + a4**3)
f(338) = 2*b1*b2*r0**2*(a1 + a2 + a3 + a4)
f(339) = 2*b1*b2*r0**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(340) = 2*b1*b2*r0**3*(a1 + a2 + a3 + a4)
f(341) = 2*dtau*r0*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(342) = 2*dtau**3*r0*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(343) = 2*dtau**2*r0*(a1*b1**2 + a2*b1**2 + a3*b2**2 + a4*b2**2)
f(344) = 2*dtau*r0*(a1*b1**3 - a2*b1**3 + a3*b2**3 - a4*b2**3)
f(345) = 2*dtau*r0*(a1**2*b1 - a2**2*b1 + a3**2*b2 - a4**2*b2)
f(346) = 2*dtau*r0*(a1**3*b1 - a2**3*b1 + a3**3*b2 - a4**3*b2)
f(347) = 2*dtau*r0**2*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(348) = 2*dtau*r0**2*(a1**2*b1 - a2**2*b1 + a3**2*b2 - a4**2*b2)
f(349) = 2*dtau*r0**3*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(350) = 2*dtau*r0*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(351) = 2*dtau**3*r0*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(352) = 2*dtau**2*r0*(a1*b2**2 + a2*b2**2 + a3*b1**2 + a4*b1**2)
f(353) = 2*dtau*r0*(a1*b2**3 - a2*b2**3 + a3*b1**3 - a4*b1**3)
f(354) = 2*dtau*r0*(a1**2*b2 - a2**2*b2 + a3**2*b1 - a4**2*b1)
f(355) = 2*dtau*r0*(a1**3*b2 - a2**3*b2 + a3**3*b1 - a4**3*b1)
f(356) = 2*dtau*r0**2*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(357) = 2*dtau*r0**2*(a1**2*b2 - a2**2*b2 + a3**2*b1 - a4**2*b1)
f(358) = 2*dtau*r0**3*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(359) = 8*b1*b2*dtau**2*r0
f(360) = 8*b1*b2*dtau**2*r0**2
f(361) = 8*r1*r2*r3*r4
f(362) = 2*r1*r2*r3*r4*(r1 + r2 + r3 + r4)
f(363) = 2*r1*r2*r3*r4*(r1**2 + r2**2 + r3**2 + r4**2)
f(364) = 4*r1*r2*r3*r4*(r1*r2 + r3*r4)
f(365) = 4*r1*r2*r3*r4*(r1*r3 + r2*r4)
f(366) = 4*r1*r2*r3*r4*(r1*r4 + r2*r3)
f(367) = 2*a1*r1*r2*r3 + 2*a2*r1*r2*r4 + 2*a3*r1*r3*r4 + 2*a4*r2*r3*r4
f(368) = 2*a1**2*r1*r2*r3 + 2*a2**2*r1*r2*r4 + 2*a3**2*r1*r3*r4 + 2*a4** &
      2*r2*r3*r4
f(369) = 2*a1**3*r1*r2*r3 + 2*a2**3*r1*r2*r4 + 2*a3**3*r1*r3*r4 + 2*a4** &
      3*r2*r3*r4
f(370) = 2*a1*r1*r2*r3**2 + 2*a2*r1*r2*r4**2 + 2*a3*r1**2*r3*r4 + 2*a4* &
      r2**2*r3*r4
f(371) = 2*a1**2*r1*r2*r3**2 + 2*a2**2*r1*r2*r4**2 + 2*a3**2*r1**2*r3*r4 &
      + 2*a4**2*r2**2*r3*r4
f(372) = 2*a1*r1*r2*r3**3 + 2*a2*r1*r2*r4**3 + 2*a3*r1**3*r3*r4 + 2*a4* &
      r2**3*r3*r4
f(373) = 2*a1*r1*r2**2*r3 + 2*a2*r1**2*r2*r4 + 2*a3*r1*r3*r4**2 + 2*a4* &
      r2*r3**2*r4
f(374) = 2*a1**2*r1*r2**2*r3 + 2*a2**2*r1**2*r2*r4 + 2*a3**2*r1*r3*r4**2 &
      + 2*a4**2*r2*r3**2*r4
f(375) = 2*a1*r1*r2**2*r3**2 + 2*a2*r1**2*r2*r4**2 + 2*a3*r1**2*r3*r4**2 &
      + 2*a4*r2**2*r3**2*r4
f(376) = 2*a1*r1*r2**3*r3 + 2*a2*r1**3*r2*r4 + 2*a3*r1*r3*r4**3 + 2*a4* &
      r2*r3**3*r4
f(377) = 2*a1*r1**2*r2*r3 + 2*a2*r1*r2**2*r4 + 2*a3*r1*r3**2*r4 + 2*a4* &
      r2*r3*r4**2
f(378) = 2*a1**2*r1**2*r2*r3 + 2*a2**2*r1*r2**2*r4 + 2*a3**2*r1*r3**2*r4 &
      + 2*a4**2*r2*r3*r4**2
f(379) = 2*a1*r1**2*r2*r3**2 + 2*a2*r1*r2**2*r4**2 + 2*a3*r1**2*r3**2*r4 &
      + 2*a4*r2**2*r3*r4**2
f(380) = 2*a1*r1**2*r2**2*r3 + 2*a2*r1**2*r2**2*r4 + 2*a3*r1*r3**2*r4**2 &
      + 2*a4*r2*r3**2*r4**2
f(381) = 2*a1*r1**3*r2*r3 + 2*a2*r1*r2**3*r4 + 2*a3*r1*r3**3*r4 + 2*a4* &
      r2*r3*r4**3
f(382) = 2*a1*r1*r2*r4 + 2*a2*r1*r2*r3 + 2*a3*r2*r3*r4 + 2*a4*r1*r3*r4
f(383) = 2*a1**2*r1*r2*r4 + 2*a2**2*r1*r2*r3 + 2*a3**2*r2*r3*r4 + 2*a4** &
      2*r1*r3*r4
f(384) = 2*a1**3*r1*r2*r4 + 2*a2**3*r1*r2*r3 + 2*a3**3*r2*r3*r4 + 2*a4** &
      3*r1*r3*r4
f(385) = 2*a1*r1*r2*r4**2 + 2*a2*r1*r2*r3**2 + 2*a3*r2**2*r3*r4 + 2*a4* &
      r1**2*r3*r4
f(386) = 2*a1**2*r1*r2*r4**2 + 2*a2**2*r1*r2*r3**2 + 2*a3**2*r2**2*r3*r4 &
      + 2*a4**2*r1**2*r3*r4
f(387) = 2*a1*r1*r2*r4**3 + 2*a2*r1*r2*r3**3 + 2*a3*r2**3*r3*r4 + 2*a4* &
      r1**3*r3*r4
f(388) = 2*a1*r1**2*r2*r4 + 2*a2*r1*r2**2*r3 + 2*a3*r2*r3**2*r4 + 2*a4* &
      r1*r3*r4**2
f(389) = 2*a1**2*r1**2*r2*r4 + 2*a2**2*r1*r2**2*r3 + 2*a3**2*r2*r3**2*r4 &
      + 2*a4**2*r1*r3*r4**2
f(390) = 2*a1*r1**2*r2*r4**2 + 2*a2*r1*r2**2*r3**2 + 2*a3*r2**2*r3**2*r4 &
      + 2*a4*r1**2*r3*r4**2
f(391) = 2*a1*r1**3*r2*r4 + 2*a2*r1*r2**3*r3 + 2*a3*r2*r3**3*r4 + 2*a4* &
      r1*r3*r4**3
f(392) = 2*a1*r1*r2**2*r4 + 2*a2*r1**2*r2*r3 + 2*a3*r2*r3*r4**2 + 2*a4* &
      r1*r3**2*r4
f(393) = 2*a1**2*r1*r2**2*r4 + 2*a2**2*r1**2*r2*r3 + 2*a3**2*r2*r3*r4**2 &
      + 2*a4**2*r1*r3**2*r4
f(394) = 2*a1*r1*r2**2*r4**2 + 2*a2*r1**2*r2*r3**2 + 2*a3*r2**2*r3*r4**2 &
      + 2*a4*r1**2*r3**2*r4
f(395) = 2*a1*r1**2*r2**2*r4 + 2*a2*r1**2*r2**2*r3 + 2*a3*r2*r3**2*r4**2 &
      + 2*a4*r1*r3**2*r4**2
f(396) = 2*a1*r1*r2**3*r4 + 2*a2*r1**3*r2*r3 + 2*a3*r2*r3*r4**3 + 2*a4* &
      r1*r3**3*r4
f(397) = 2*a1*r1*r3*r4 + 2*a2*r2*r3*r4 + 2*a3*r1*r2*r3 + 2*a4*r1*r2*r4
f(398) = 2*a1**2*r1*r3*r4 + 2*a2**2*r2*r3*r4 + 2*a3**2*r1*r2*r3 + 2*a4** &
      2*r1*r2*r4
f(399) = 2*a1**3*r1*r3*r4 + 2*a2**3*r2*r3*r4 + 2*a3**3*r1*r2*r3 + 2*a4** &
      3*r1*r2*r4
f(400) = 2*a1*r1**2*r3*r4 + 2*a2*r2**2*r3*r4 + 2*a3*r1*r2*r3**2 + 2*a4* &
      r1*r2*r4**2
f(401) = 2*a1**2*r1**2*r3*r4 + 2*a2**2*r2**2*r3*r4 + 2*a3**2*r1*r2*r3**2 &
      + 2*a4**2*r1*r2*r4**2
f(402) = 2*a1*r1**3*r3*r4 + 2*a2*r2**3*r3*r4 + 2*a3*r1*r2*r3**3 + 2*a4* &
      r1*r2*r4**3
f(403) = 2*a1*r1*r3*r4**2 + 2*a2*r2*r3**2*r4 + 2*a3*r1*r2**2*r3 + 2*a4* &
      r1**2*r2*r4
f(404) = 2*a1**2*r1*r3*r4**2 + 2*a2**2*r2*r3**2*r4 + 2*a3**2*r1*r2**2*r3 &
      + 2*a4**2*r1**2*r2*r4
f(405) = 2*a1*r1**2*r3*r4**2 + 2*a2*r2**2*r3**2*r4 + 2*a3*r1*r2**2*r3**2 &
      + 2*a4*r1**2*r2*r4**2
f(406) = 2*a1*r1*r3*r4**3 + 2*a2*r2*r3**3*r4 + 2*a3*r1*r2**3*r3 + 2*a4* &
      r1**3*r2*r4
f(407) = 2*a1*r1*r3**2*r4 + 2*a2*r2*r3*r4**2 + 2*a3*r1**2*r2*r3 + 2*a4* &
      r1*r2**2*r4
f(408) = 2*a1**2*r1*r3**2*r4 + 2*a2**2*r2*r3*r4**2 + 2*a3**2*r1**2*r2*r3 &
      + 2*a4**2*r1*r2**2*r4
f(409) = 2*a1*r1**2*r3**2*r4 + 2*a2*r2**2*r3*r4**2 + 2*a3*r1**2*r2*r3**2 &
      + 2*a4*r1*r2**2*r4**2
f(410) = 2*a1*r1*r3**2*r4**2 + 2*a2*r2*r3**2*r4**2 + 2*a3*r1**2*r2**2*r3 &
      + 2*a4*r1**2*r2**2*r4
f(411) = 2*a1*r1*r3**3*r4 + 2*a2*r2*r3*r4**3 + 2*a3*r1**3*r2*r3 + 2*a4* &
      r1*r2**3*r4
f(412) = 2*a1*r2*r3*r4 + 2*a2*r1*r3*r4 + 2*a3*r1*r2*r4 + 2*a4*r1*r2*r3
f(413) = 2*a1**2*r2*r3*r4 + 2*a2**2*r1*r3*r4 + 2*a3**2*r1*r2*r4 + 2*a4** &
      2*r1*r2*r3
f(414) = 2*a1**3*r2*r3*r4 + 2*a2**3*r1*r3*r4 + 2*a3**3*r1*r2*r4 + 2*a4** &
      3*r1*r2*r3
f(415) = 2*a1*r2**2*r3*r4 + 2*a2*r1**2*r3*r4 + 2*a3*r1*r2*r4**2 + 2*a4* &
      r1*r2*r3**2
f(416) = 2*a1**2*r2**2*r3*r4 + 2*a2**2*r1**2*r3*r4 + 2*a3**2*r1*r2*r4**2 &
      + 2*a4**2*r1*r2*r3**2
f(417) = 2*a1*r2**3*r3*r4 + 2*a2*r1**3*r3*r4 + 2*a3*r1*r2*r4**3 + 2*a4* &
      r1*r2*r3**3
f(418) = 2*a1*r2*r3**2*r4 + 2*a2*r1*r3*r4**2 + 2*a3*r1**2*r2*r4 + 2*a4* &
      r1*r2**2*r3
f(419) = 2*a1**2*r2*r3**2*r4 + 2*a2**2*r1*r3*r4**2 + 2*a3**2*r1**2*r2*r4 &
      + 2*a4**2*r1*r2**2*r3
f(420) = 2*a1*r2**2*r3**2*r4 + 2*a2*r1**2*r3*r4**2 + 2*a3*r1**2*r2*r4**2 &
      + 2*a4*r1*r2**2*r3**2
f(421) = 2*a1*r2*r3**3*r4 + 2*a2*r1*r3*r4**3 + 2*a3*r1**3*r2*r4 + 2*a4* &
      r1*r2**3*r3
f(422) = 2*a1*r2*r3*r4**2 + 2*a2*r1*r3**2*r4 + 2*a3*r1*r2**2*r4 + 2*a4* &
      r1**2*r2*r3
f(423) = 2*a1**2*r2*r3*r4**2 + 2*a2**2*r1*r3**2*r4 + 2*a3**2*r1*r2**2*r4 &
      + 2*a4**2*r1**2*r2*r3
f(424) = 2*a1*r2**2*r3*r4**2 + 2*a2*r1**2*r3**2*r4 + 2*a3*r1*r2**2*r4**2 &
      + 2*a4*r1**2*r2*r3**2
f(425) = 2*a1*r2*r3**2*r4**2 + 2*a2*r1*r3**2*r4**2 + 2*a3*r1**2*r2**2*r4 &
      + 2*a4*r1**2*r2**2*r3
f(426) = 2*a1*r2*r3*r4**3 + 2*a2*r1*r3**3*r4 + 2*a3*r1*r2**3*r4 + 2*a4* &
      r1**3*r2*r3
f(427) = 2*b1**2*r1*r2*r3 + 2*b1**2*r1*r2*r4 + 2*b2**2*r1*r3*r4 + 2*b2** &
      2*r2*r3*r4
f(428) = 2*b1**2*r1*r2*r3**2 + 2*b1**2*r1*r2*r4**2 + 2*b2**2*r1**2*r3*r4 &
      + 2*b2**2*r2**2*r3*r4
f(429) = 2*b1**2*r1**2*r2*r4 + 2*b1**2*r1*r2**2*r3 + 2*b2**2*r1*r3*r4**2 &
      + 2*b2**2*r2*r3**2*r4
f(430) = 2*b1**2*r1**2*r2*r3 + 2*b1**2*r1*r2**2*r4 + 2*b2**2*r1*r3**2*r4 &
      + 2*b2**2*r2*r3*r4**2
f(431) = 2*b1**2*r1*r3*r4 + 2*b1**2*r2*r3*r4 + 2*b2**2*r1*r2*r3 + 2*b2** &
      2*r1*r2*r4
f(432) = 2*b1**2*r1**2*r3*r4 + 2*b1**2*r2**2*r3*r4 + 2*b2**2*r1*r2*r3**2 &
      + 2*b2**2*r1*r2*r4**2
f(433) = 2*b1**2*r1*r3*r4**2 + 2*b1**2*r2*r3**2*r4 + 2*b2**2*r1**2*r2*r4 &
      + 2*b2**2*r1*r2**2*r3
f(434) = 2*b1**2*r1*r3**2*r4 + 2*b1**2*r2*r3*r4**2 + 2*b2**2*r1**2*r2*r3 &
      + 2*b2**2*r1*r2**2*r4
f(435) = 2*dtau**2*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(436) = 2*dtau**2*(r1**2*r3*r4 + r1*r2*r3**2 + r1*r2*r4**2 + r2**2*r3* &
      r4)
f(437) = 2*dtau**2*(r1**2*r2*r4 + r1*r2**2*r3 + r1*r3*r4**2 + r2*r3**2* &
      r4)
f(438) = 2*dtau**2*(r1**2*r2*r3 + r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4** &
      2)
f(439) = 4*a1*a2*r1*r2 + 4*a3*a4*r3*r4
f(440) = 2*a1**2*a2*r1*r2 + 2*a1*a2**2*r1*r2 + 2*a3**2*a4*r3*r4 + 2*a3* &
      a4**2*r3*r4
f(441) = 2*a1**3*a2*r1*r2 + 2*a1*a2**3*r1*r2 + 2*a3**3*a4*r3*r4 + 2*a3* &
      a4**3*r3*r4
f(442) = 4*a1**2*a2**2*r1*r2 + 4*a3**2*a4**2*r3*r4
f(443) = 2*a1*a2*r1**2*r2 + 2*a1*a2*r1*r2**2 + 2*a3*a4*r3**2*r4 + 2*a3* &
      a4*r3*r4**2
f(444) = 2*a1**2*a2*r1**2*r2 + 2*a1*a2**2*r1*r2**2 + 2*a3**2*a4*r3**2*r4 &
      + 2*a3*a4**2*r3*r4**2
f(445) = 2*a1**2*a2*r1*r2**2 + 2*a1*a2**2*r1**2*r2 + 2*a3**2*a4*r3*r4**2 &
      + 2*a3*a4**2*r3**2*r4
f(446) = 2*a1*a2*r1**3*r2 + 2*a1*a2*r1*r2**3 + 2*a3*a4*r3**3*r4 + 2*a3* &
      a4*r3*r4**3
f(447) = 4*a1*a2*r1**2*r2**2 + 4*a3*a4*r3**2*r4**2
f(448) = 2*a1*a3*r1*r2 + 2*a1*a3*r3*r4 + 2*a2*a4*r1*r2 + 2*a2*a4*r3*r4
f(449) = 2*a1**2*a3*r3*r4 + 2*a1*a3**2*r1*r2 + 2*a2**2*a4*r3*r4 + 2*a2* &
      a4**2*r1*r2
f(450) = 2*a1**3*a3*r3*r4 + 2*a1*a3**3*r1*r2 + 2*a2**3*a4*r3*r4 + 2*a2* &
      a4**3*r1*r2
f(451) = 2*a1**2*a3*r1*r2 + 2*a1*a3**2*r3*r4 + 2*a2**2*a4*r1*r2 + 2*a2* &
      a4**2*r3*r4
f(452) = 2*a1**2*a3**2*r1*r2 + 2*a1**2*a3**2*r3*r4 + 2*a2**2*a4**2*r1*r2 &
      + 2*a2**2*a4**2*r3*r4
f(453) = 2*a1**3*a3*r1*r2 + 2*a1*a3**3*r3*r4 + 2*a2**3*a4*r1*r2 + 2*a2* &
      a4**3*r3*r4
f(454) = 2*a1*a3*r1*r2**2 + 2*a1*a3*r3*r4**2 + 2*a2*a4*r1**2*r2 + 2*a2* &
      a4*r3**2*r4
f(455) = 2*a1**2*a3*r3*r4**2 + 2*a1*a3**2*r1*r2**2 + 2*a2**2*a4*r3**2*r4 &
      + 2*a2*a4**2*r1**2*r2
f(456) = 2*a1**2*a3*r1*r2**2 + 2*a1*a3**2*r3*r4**2 + 2*a2**2*a4*r1**2*r2 &
      + 2*a2*a4**2*r3**2*r4
f(457) = 2*a1*a3*r1*r2**3 + 2*a1*a3*r3*r4**3 + 2*a2*a4*r1**3*r2 + 2*a2* &
      a4*r3**3*r4
f(458) = 2*a1*a3*r1**2*r2 + 2*a1*a3*r3**2*r4 + 2*a2*a4*r1*r2**2 + 2*a2* &
      a4*r3*r4**2
f(459) = 2*a1**2*a3*r3**2*r4 + 2*a1*a3**2*r1**2*r2 + 2*a2**2*a4*r3*r4**2 &
      + 2*a2*a4**2*r1*r2**2
f(460) = 2*a1**2*a3*r1**2*r2 + 2*a1*a3**2*r3**2*r4 + 2*a2**2*a4*r1*r2**2 &
      + 2*a2*a4**2*r3*r4**2
f(461) = 2*a1*a3*r1**2*r2**2 + 2*a1*a3*r3**2*r4**2 + 2*a2*a4*r1**2*r2**2 &
      + 2*a2*a4*r3**2*r4**2
f(462) = 2*a1*a3*r1**3*r2 + 2*a1*a3*r3**3*r4 + 2*a2*a4*r1*r2**3 + 2*a2* &
      a4*r3*r4**3
f(463) = 2*a1*a4*r1*r2 + 2*a1*a4*r3*r4 + 2*a2*a3*r1*r2 + 2*a2*a3*r3*r4
f(464) = 2*a1**2*a4*r3*r4 + 2*a1*a4**2*r1*r2 + 2*a2**2*a3*r3*r4 + 2*a2* &
      a3**2*r1*r2
f(465) = 2*a1**3*a4*r3*r4 + 2*a1*a4**3*r1*r2 + 2*a2**3*a3*r3*r4 + 2*a2* &
      a3**3*r1*r2
f(466) = 2*a1**2*a4*r1*r2 + 2*a1*a4**2*r3*r4 + 2*a2**2*a3*r1*r2 + 2*a2* &
      a3**2*r3*r4
f(467) = 2*a1**2*a4**2*r1*r2 + 2*a1**2*a4**2*r3*r4 + 2*a2**2*a3**2*r1*r2 &
      + 2*a2**2*a3**2*r3*r4
f(468) = 2*a1**3*a4*r1*r2 + 2*a1*a4**3*r3*r4 + 2*a2**3*a3*r1*r2 + 2*a2* &
      a3**3*r3*r4
f(469) = 2*a1*a4*r1*r2**2 + 2*a1*a4*r3**2*r4 + 2*a2*a3*r1**2*r2 + 2*a2* &
      a3*r3*r4**2
f(470) = 2*a1**2*a4*r3**2*r4 + 2*a1*a4**2*r1*r2**2 + 2*a2**2*a3*r3*r4**2 &
      + 2*a2*a3**2*r1**2*r2
f(471) = 2*a1**2*a4*r1*r2**2 + 2*a1*a4**2*r3**2*r4 + 2*a2**2*a3*r1**2*r2 &
      + 2*a2*a3**2*r3*r4**2
f(472) = 2*a1*a4*r1*r2**3 + 2*a1*a4*r3**3*r4 + 2*a2*a3*r1**3*r2 + 2*a2* &
      a3*r3*r4**3
f(473) = 2*a1*a4*r1**2*r2 + 2*a1*a4*r3*r4**2 + 2*a2*a3*r1*r2**2 + 2*a2* &
      a3*r3**2*r4
f(474) = 2*a1**2*a4*r3*r4**2 + 2*a1*a4**2*r1**2*r2 + 2*a2**2*a3*r3**2*r4 &
      + 2*a2*a3**2*r1*r2**2
f(475) = 2*a1**2*a4*r1**2*r2 + 2*a1*a4**2*r3*r4**2 + 2*a2**2*a3*r1*r2**2 &
      + 2*a2*a3**2*r3**2*r4
f(476) = 2*a1*a4*r1**2*r2**2 + 2*a1*a4*r3**2*r4**2 + 2*a2*a3*r1**2*r2**2 &
      + 2*a2*a3*r3**2*r4**2
f(477) = 2*a1*a4*r1**3*r2 + 2*a1*a4*r3*r4**3 + 2*a2*a3*r1*r2**3 + 2*a2* &
      a3*r3**3*r4
f(478) = 2*a1*b1**2*r1*r2 + 2*a2*b1**2*r1*r2 + 2*a3*b2**2*r3*r4 + 2*a4* &
      b2**2*r3*r4
f(479) = 2*a1**2*b1**2*r1*r2 + 2*a2**2*b1**2*r1*r2 + 2*a3**2*b2**2*r3*r4 &
      + 2*a4**2*b2**2*r3*r4
f(480) = 2*a1*b1**2*r1*r2**2 + 2*a2*b1**2*r1**2*r2 + 2*a3*b2**2*r3*r4**2 &
      + 2*a4*b2**2*r3**2*r4
f(481) = 2*a1*b1**2*r1**2*r2 + 2*a2*b1**2*r1*r2**2 + 2*a3*b2**2*r3**2*r4 &
      + 2*a4*b2**2*r3*r4**2
f(482) = 2*a1*b2**2*r1*r2 + 2*a2*b2**2*r1*r2 + 2*a3*b1**2*r3*r4 + 2*a4* &
      b1**2*r3*r4
f(483) = 2*a1**2*b2**2*r1*r2 + 2*a2**2*b2**2*r1*r2 + 2*a3**2*b1**2*r3*r4 &
      + 2*a4**2*b1**2*r3*r4
f(484) = 2*a1*b2**2*r1*r2**2 + 2*a2*b2**2*r1**2*r2 + 2*a3*b1**2*r3*r4**2 &
      + 2*a4*b1**2*r3**2*r4
f(485) = 2*a1*b2**2*r1**2*r2 + 2*a2*b2**2*r1*r2**2 + 2*a3*b1**2*r3**2*r4 &
      + 2*a4*b1**2*r3*r4**2
f(486) = 2*dtau**2*(a1*r1*r2 + a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(487) = 2*dtau**2*(a1**2*r1*r2 + a2**2*r1*r2 + a3**2*r3*r4 + a4**2*r3* &
      r4)
f(488) = 2*dtau**2*(a1*r1*r2**2 + a2*r1**2*r2 + a3*r3*r4**2 + a4*r3**2* &
      r4)
f(489) = 2*dtau**2*(a1*r1**2*r2 + a2*r1*r2**2 + a3*r3**2*r4 + a4*r3*r4** &
      2)
f(490) = 4*a1*a2*r3*r4 + 4*a3*a4*r1*r2
f(491) = 2*a1**2*a2*r3*r4 + 2*a1*a2**2*r3*r4 + 2*a3**2*a4*r1*r2 + 2*a3* &
      a4**2*r1*r2
f(492) = 2*a1**3*a2*r3*r4 + 2*a1*a2**3*r3*r4 + 2*a3**3*a4*r1*r2 + 2*a3* &
      a4**3*r1*r2
f(493) = 4*a1**2*a2**2*r3*r4 + 4*a3**2*a4**2*r1*r2
f(494) = 2*a1*a2*r3**2*r4 + 2*a1*a2*r3*r4**2 + 2*a3*a4*r1**2*r2 + 2*a3* &
      a4*r1*r2**2
f(495) = 2*a1**2*a2*r3**2*r4 + 2*a1*a2**2*r3*r4**2 + 2*a3**2*a4*r1**2*r2 &
      + 2*a3*a4**2*r1*r2**2
f(496) = 2*a1**2*a2*r3*r4**2 + 2*a1*a2**2*r3**2*r4 + 2*a3**2*a4*r1*r2**2 &
      + 2*a3*a4**2*r1**2*r2
f(497) = 2*a1*a2*r3**3*r4 + 2*a1*a2*r3*r4**3 + 2*a3*a4*r1**3*r2 + 2*a3* &
      a4*r1*r2**3
f(498) = 4*a1*a2*r3**2*r4**2 + 4*a3*a4*r1**2*r2**2
f(499) = 2*a1*b2**2*r3*r4 + 2*a2*b2**2*r3*r4 + 2*a3*b1**2*r1*r2 + 2*a4* &
      b1**2*r1*r2
f(500) = 2*a1**2*b2**2*r3*r4 + 2*a2**2*b2**2*r3*r4 + 2*a3**2*b1**2*r1*r2 &
      + 2*a4**2*b1**2*r1*r2
f(501) = 2*a1*b2**2*r3*r4**2 + 2*a2*b2**2*r3**2*r4 + 2*a3*b1**2*r1*r2**2 &
      + 2*a4*b1**2*r1**2*r2
f(502) = 2*a1*b2**2*r3**2*r4 + 2*a2*b2**2*r3*r4**2 + 2*a3*b1**2*r1**2*r2 &
      + 2*a4*b1**2*r1*r2**2
f(503) = 2*a1*b1**2*r3*r4 + 2*a2*b1**2*r3*r4 + 2*a3*b2**2*r1*r2 + 2*a4* &
      b2**2*r1*r2
f(504) = 2*a1**2*b1**2*r3*r4 + 2*a2**2*b1**2*r3*r4 + 2*a3**2*b2**2*r1*r2 &
      + 2*a4**2*b2**2*r1*r2
f(505) = 2*a1*b1**2*r3*r4**2 + 2*a2*b1**2*r3**2*r4 + 2*a3*b2**2*r1*r2**2 &
      + 2*a4*b2**2*r1**2*r2
f(506) = 2*a1*b1**2*r3**2*r4 + 2*a2*b1**2*r3*r4**2 + 2*a3*b2**2*r1**2*r2 &
      + 2*a4*b2**2*r1*r2**2
f(507) = 2*dtau**2*(a1*r3*r4 + a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(508) = 2*dtau**2*(a1**2*r3*r4 + a2**2*r3*r4 + a3**2*r1*r2 + a4**2*r1* &
      r2)
f(509) = 2*dtau**2*(a1*r3*r4**2 + a2*r3**2*r4 + a3*r1*r2**2 + a4*r1**2* &
      r2)
f(510) = 2*dtau**2*(a1*r3**2*r4 + a2*r3*r4**2 + a3*r1**2*r2 + a4*r1*r2** &
      2)
f(511) = 4*b1*b2*(r1*r2 + r3*r4)
f(512) = 4*b1*b2*(b1**2*r3*r4 + b2**2*r1*r2)
f(513) = 4*b1**2*b2**2*(r1*r2 + r3*r4)
f(514) = 4*b1*b2*(b1**2*r1*r2 + b2**2*r3*r4)
f(515) = 2*b1*b2*(r1**2*r2 + r1*r2**2 + r3**2*r4 + r3*r4**2)
f(516) = 2*b1*b2*(r1**3*r2 + r1*r2**3 + r3**3*r4 + r3*r4**3)
f(517) = 4*b1*b2*(r1**2*r2**2 + r3**2*r4**2)
f(518) = 4*dtau**2*(b1**2*r1*r2 + b2**2*r3*r4)
f(519) = 2*dtau*(-b1*r1**2*r2 + b1*r1*r2**2 - b2*r3**2*r4 + b2*r3*r4**2)
f(520) = 2*dtau*(-b1*r1**3*r2 + b1*r1*r2**3 - b2*r3**3*r4 + b2*r3*r4**3)
f(521) = 2*dtau*(b1*r1**2*r2 - b1*r1*r2**2 + b2*r3**2*r4 - b2*r3*r4**2)
f(522) = 2*dtau*(b1*r1**3*r2 - b1*r1*r2**3 + b2*r3**3*r4 - b2*r3*r4**3)
f(523) = 4*dtau**2*(b1**2*r3*r4 + b2**2*r1*r2)
f(524) = 2*dtau*(-b1*r3**2*r4 + b1*r3*r4**2 - b2*r1**2*r2 + b2*r1*r2**2)
f(525) = 2*dtau*(-b1*r3**3*r4 + b1*r3*r4**3 - b2*r1**3*r2 + b2*r1*r2**3)
f(526) = 2*dtau*(b1*r3**2*r4 - b1*r3*r4**2 + b2*r1**2*r2 - b2*r1*r2**2)
f(527) = 2*dtau*(b1*r3**3*r4 - b1*r3*r4**3 + b2*r1**3*r2 - b2*r1*r2**3)
f(528) = 2*a1*a2*r1*r3 + 2*a1*a2*r2*r4 + 2*a3*a4*r1*r3 + 2*a3*a4*r2*r4
f(529) = 2*a1**2*a2*r2*r4 + 2*a1*a2**2*r1*r3 + 2*a3**2*a4*r2*r4 + 2*a3* &
      a4**2*r1*r3
f(530) = 2*a1**3*a2*r2*r4 + 2*a1*a2**3*r1*r3 + 2*a3**3*a4*r2*r4 + 2*a3* &
      a4**3*r1*r3
f(531) = 2*a1**2*a2*r1*r3 + 2*a1*a2**2*r2*r4 + 2*a3**2*a4*r1*r3 + 2*a3* &
      a4**2*r2*r4
f(532) = 2*a1**2*a2**2*r1*r3 + 2*a1**2*a2**2*r2*r4 + 2*a3**2*a4**2*r1*r3 &
      + 2*a3**2*a4**2*r2*r4
f(533) = 2*a1**3*a2*r1*r3 + 2*a1*a2**3*r2*r4 + 2*a3**3*a4*r1*r3 + 2*a3* &
      a4**3*r2*r4
f(534) = 2*a1*a2*r1*r3**2 + 2*a1*a2*r2*r4**2 + 2*a3*a4*r1**2*r3 + 2*a3* &
      a4*r2**2*r4
f(535) = 2*a1**2*a2*r2*r4**2 + 2*a1*a2**2*r1*r3**2 + 2*a3**2*a4*r2**2*r4 &
      + 2*a3*a4**2*r1**2*r3
f(536) = 2*a1**2*a2*r1*r3**2 + 2*a1*a2**2*r2*r4**2 + 2*a3**2*a4*r1**2*r3 &
      + 2*a3*a4**2*r2**2*r4
f(537) = 2*a1*a2*r1*r3**3 + 2*a1*a2*r2*r4**3 + 2*a3*a4*r1**3*r3 + 2*a3* &
      a4*r2**3*r4
f(538) = 2*a1*a2*r1**2*r3 + 2*a1*a2*r2**2*r4 + 2*a3*a4*r1*r3**2 + 2*a3* &
      a4*r2*r4**2
f(539) = 2*a1**2*a2*r2**2*r4 + 2*a1*a2**2*r1**2*r3 + 2*a3**2*a4*r2*r4**2 &
      + 2*a3*a4**2*r1*r3**2
f(540) = 2*a1**2*a2*r1**2*r3 + 2*a1*a2**2*r2**2*r4 + 2*a3**2*a4*r1*r3**2 &
      + 2*a3*a4**2*r2*r4**2
f(541) = 2*a1*a2*r1**2*r3**2 + 2*a1*a2*r2**2*r4**2 + 2*a3*a4*r1**2*r3**2 &
      + 2*a3*a4*r2**2*r4**2
f(542) = 2*a1*a2*r1**3*r3 + 2*a1*a2*r2**3*r4 + 2*a3*a4*r1*r3**3 + 2*a3* &
      a4*r2*r4**3
f(543) = 4*a1*a3*r1*r3 + 4*a2*a4*r2*r4
f(544) = 2*a1**2*a3*r1*r3 + 2*a1*a3**2*r1*r3 + 2*a2**2*a4*r2*r4 + 2*a2* &
      a4**2*r2*r4
f(545) = 2*a1**3*a3*r1*r3 + 2*a1*a3**3*r1*r3 + 2*a2**3*a4*r2*r4 + 2*a2* &
      a4**3*r2*r4
f(546) = 4*a1**2*a3**2*r1*r3 + 4*a2**2*a4**2*r2*r4
f(547) = 2*a1*a3*r1**2*r3 + 2*a1*a3*r1*r3**2 + 2*a2*a4*r2**2*r4 + 2*a2* &
      a4*r2*r4**2
f(548) = 2*a1**2*a3*r1**2*r3 + 2*a1*a3**2*r1*r3**2 + 2*a2**2*a4*r2**2*r4 &
      + 2*a2*a4**2*r2*r4**2
f(549) = 2*a1**2*a3*r1*r3**2 + 2*a1*a3**2*r1**2*r3 + 2*a2**2*a4*r2*r4**2 &
      + 2*a2*a4**2*r2**2*r4
f(550) = 2*a1*a3*r1**3*r3 + 2*a1*a3*r1*r3**3 + 2*a2*a4*r2**3*r4 + 2*a2* &
      a4*r2*r4**3
f(551) = 4*a1*a3*r1**2*r3**2 + 4*a2*a4*r2**2*r4**2
f(552) = 2*a1*a4*r1*r3 + 2*a1*a4*r2*r4 + 2*a2*a3*r1*r3 + 2*a2*a3*r2*r4
f(553) = 2*a1**2*a4*r2*r4 + 2*a1*a4**2*r1*r3 + 2*a2**2*a3*r1*r3 + 2*a2* &
      a3**2*r2*r4
f(554) = 2*a1**3*a4*r2*r4 + 2*a1*a4**3*r1*r3 + 2*a2**3*a3*r1*r3 + 2*a2* &
      a3**3*r2*r4
f(555) = 2*a1**2*a4*r1*r3 + 2*a1*a4**2*r2*r4 + 2*a2**2*a3*r2*r4 + 2*a2* &
      a3**2*r1*r3
f(556) = 2*a1**2*a4**2*r1*r3 + 2*a1**2*a4**2*r2*r4 + 2*a2**2*a3**2*r1*r3 &
      + 2*a2**2*a3**2*r2*r4
f(557) = 2*a1**3*a4*r1*r3 + 2*a1*a4**3*r2*r4 + 2*a2**3*a3*r2*r4 + 2*a2* &
      a3**3*r1*r3
f(558) = 2*a1*a4*r1*r3**2 + 2*a1*a4*r2**2*r4 + 2*a2*a3*r1**2*r3 + 2*a2* &
      a3*r2*r4**2
f(559) = 2*a1**2*a4*r2**2*r4 + 2*a1*a4**2*r1*r3**2 + 2*a2**2*a3*r1**2*r3 &
      + 2*a2*a3**2*r2*r4**2
f(560) = 2*a1**2*a4*r1*r3**2 + 2*a1*a4**2*r2**2*r4 + 2*a2**2*a3*r2*r4**2 &
      + 2*a2*a3**2*r1**2*r3
f(561) = 2*a1*a4*r1*r3**3 + 2*a1*a4*r2**3*r4 + 2*a2*a3*r1**3*r3 + 2*a2* &
      a3*r2*r4**3
f(562) = 2*a1*a4*r1**2*r3 + 2*a1*a4*r2*r4**2 + 2*a2*a3*r1*r3**2 + 2*a2* &
      a3*r2**2*r4
f(563) = 2*a1**2*a4*r2*r4**2 + 2*a1*a4**2*r1**2*r3 + 2*a2**2*a3*r1*r3**2 &
      + 2*a2*a3**2*r2**2*r4
f(564) = 2*a1**2*a4*r1**2*r3 + 2*a1*a4**2*r2*r4**2 + 2*a2**2*a3*r2**2*r4 &
      + 2*a2*a3**2*r1*r3**2
f(565) = 2*a1*a4*r1**2*r3**2 + 2*a1*a4*r2**2*r4**2 + 2*a2*a3*r1**2*r3**2 &
      + 2*a2*a3*r2**2*r4**2
f(566) = 2*a1*a4*r1**3*r3 + 2*a1*a4*r2*r4**3 + 2*a2*a3*r1*r3**3 + 2*a2* &
      a3*r2**3*r4
f(567) = 2*a1*b1**2*r1*r3 + 2*a2*b1**2*r2*r4 + 2*a3*b2**2*r1*r3 + 2*a4* &
      b2**2*r2*r4
f(568) = 2*a1**2*b1**2*r1*r3 + 2*a2**2*b1**2*r2*r4 + 2*a3**2*b2**2*r1*r3 &
      + 2*a4**2*b2**2*r2*r4
f(569) = 2*a1*b1**2*r1*r3**2 + 2*a2*b1**2*r2*r4**2 + 2*a3*b2**2*r1**2*r3 &
      + 2*a4*b2**2*r2**2*r4
f(570) = 2*a1*b1**2*r1**2*r3 + 2*a2*b1**2*r2**2*r4 + 2*a3*b2**2*r1*r3**2 &
      + 2*a4*b2**2*r2*r4**2
f(571) = 2*a1*b2**2*r1*r3 + 2*a2*b2**2*r2*r4 + 2*a3*b1**2*r1*r3 + 2*a4* &
      b1**2*r2*r4
f(572) = 2*a1**2*b2**2*r1*r3 + 2*a2**2*b2**2*r2*r4 + 2*a3**2*b1**2*r1*r3 &
      + 2*a4**2*b1**2*r2*r4
f(573) = 2*a1*b2**2*r1*r3**2 + 2*a2*b2**2*r2*r4**2 + 2*a3*b1**2*r1**2*r3 &
      + 2*a4*b1**2*r2**2*r4
f(574) = 2*a1*b2**2*r1**2*r3 + 2*a2*b2**2*r2**2*r4 + 2*a3*b1**2*r1*r3**2 &
      + 2*a4*b1**2*r2*r4**2
f(575) = 2*dtau**2*(a1*r1*r3 + a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(576) = 2*dtau**2*(a1**2*r1*r3 + a2**2*r2*r4 + a3**2*r1*r3 + a4**2*r2* &
      r4)
f(577) = 2*dtau**2*(a1*r1*r3**2 + a2*r2*r4**2 + a3*r1**2*r3 + a4*r2**2* &
      r4)
f(578) = 2*dtau**2*(a1*r1**2*r3 + a2*r2**2*r4 + a3*r1*r3**2 + a4*r2*r4** &
      2)
f(579) = 4*a1*a3*r2*r4 + 4*a2*a4*r1*r3
f(580) = 2*a1**2*a3*r2*r4 + 2*a1*a3**2*r2*r4 + 2*a2**2*a4*r1*r3 + 2*a2* &
      a4**2*r1*r3
f(581) = 2*a1**3*a3*r2*r4 + 2*a1*a3**3*r2*r4 + 2*a2**3*a4*r1*r3 + 2*a2* &
      a4**3*r1*r3
f(582) = 4*a1**2*a3**2*r2*r4 + 4*a2**2*a4**2*r1*r3
f(583) = 2*a1*a3*r2**2*r4 + 2*a1*a3*r2*r4**2 + 2*a2*a4*r1**2*r3 + 2*a2* &
      a4*r1*r3**2
f(584) = 2*a1**2*a3*r2**2*r4 + 2*a1*a3**2*r2*r4**2 + 2*a2**2*a4*r1**2*r3 &
      + 2*a2*a4**2*r1*r3**2
f(585) = 2*a1**2*a3*r2*r4**2 + 2*a1*a3**2*r2**2*r4 + 2*a2**2*a4*r1*r3**2 &
      + 2*a2*a4**2*r1**2*r3
f(586) = 2*a1*a3*r2**3*r4 + 2*a1*a3*r2*r4**3 + 2*a2*a4*r1**3*r3 + 2*a2* &
      a4*r1*r3**3
f(587) = 4*a1*a3*r2**2*r4**2 + 4*a2*a4*r1**2*r3**2
f(588) = 2*a1*b1**2*r2*r4 + 2*a2*b1**2*r1*r3 + 2*a3*b2**2*r2*r4 + 2*a4* &
      b2**2*r1*r3
f(589) = 2*a1**2*b1**2*r2*r4 + 2*a2**2*b1**2*r1*r3 + 2*a3**2*b2**2*r2*r4 &
      + 2*a4**2*b2**2*r1*r3
f(590) = 2*a1*b1**2*r2*r4**2 + 2*a2*b1**2*r1*r3**2 + 2*a3*b2**2*r2**2*r4 &
      + 2*a4*b2**2*r1**2*r3
f(591) = 2*a1*b1**2*r2**2*r4 + 2*a2*b1**2*r1**2*r3 + 2*a3*b2**2*r2*r4**2 &
      + 2*a4*b2**2*r1*r3**2
f(592) = 2*a1*b2**2*r2*r4 + 2*a2*b2**2*r1*r3 + 2*a3*b1**2*r2*r4 + 2*a4* &
      b1**2*r1*r3
f(593) = 2*a1**2*b2**2*r2*r4 + 2*a2**2*b2**2*r1*r3 + 2*a3**2*b1**2*r2*r4 &
      + 2*a4**2*b1**2*r1*r3
f(594) = 2*a1*b2**2*r2*r4**2 + 2*a2*b2**2*r1*r3**2 + 2*a3*b1**2*r2**2*r4 &
      + 2*a4*b1**2*r1**2*r3
f(595) = 2*a1*b2**2*r2**2*r4 + 2*a2*b2**2*r1**2*r3 + 2*a3*b1**2*r2*r4**2 &
      + 2*a4*b1**2*r1*r3**2
f(596) = 2*dtau**2*(a1*r2*r4 + a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(597) = 2*dtau**2*(a1**2*r2*r4 + a2**2*r1*r3 + a3**2*r2*r4 + a4**2*r1* &
      r3)
f(598) = 2*dtau**2*(a1*r2*r4**2 + a2*r1*r3**2 + a3*r2**2*r4 + a4*r1**2* &
      r3)
f(599) = 2*dtau**2*(a1*r2**2*r4 + a2*r1**2*r3 + a3*r2*r4**2 + a4*r1*r3** &
      2)
f(600) = 4*b1*b2*(r1*r3 + r2*r4)
f(601) = 2*b1*b2*(b1**2*r1*r3 + b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2*r4)
f(602) = 4*b1**2*b2**2*(r1*r3 + r2*r4)
f(603) = 2*b1*b2*(r1**2*r3 + r1*r3**2 + r2**2*r4 + r2*r4**2)
f(604) = 2*b1*b2*(r1**3*r3 + r1*r3**3 + r2**3*r4 + r2*r4**3)
f(605) = 4*b1*b2*(r1**2*r3**2 + r2**2*r4**2)
f(606) = 2*dtau*(b1*r1*r3 - b1*r2*r4 + b2*r1*r3 - b2*r2*r4)
f(607) = 2*dtau**3*(b1*r1*r3 - b1*r2*r4 + b2*r1*r3 - b2*r2*r4)
f(608) = 2*dtau**2*(b1**2*r1*r3 + b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2* &
      r4)
f(609) = 2*dtau*(b1**3*r1*r3 - b1**3*r2*r4 + b2**3*r1*r3 - b2**3*r2*r4)
f(610) = 2*dtau*(b1*r1*r3**2 - b1*r2*r4**2 + b2*r1**2*r3 - b2*r2**2*r4)
f(611) = 2*dtau*(b1*r1*r3**3 - b1*r2*r4**3 + b2*r1**3*r3 - b2*r2**3*r4)
f(612) = 2*dtau*(b1*r1**2*r3 - b1*r2**2*r4 + b2*r1*r3**2 - b2*r2*r4**2)
f(613) = 2*dtau*(b1*r1**2*r3**2 - b1*r2**2*r4**2 + b2*r1**2*r3**2 - b2* &
      r2**2*r4**2)
f(614) = 2*dtau*(b1*r1**3*r3 - b1*r2**3*r4 + b2*r1*r3**3 - b2*r2*r4**3)
f(615) = 2*a1*a2*r1*r4 + 2*a1*a2*r2*r3 + 2*a3*a4*r1*r4 + 2*a3*a4*r2*r3
f(616) = 2*a1**2*a2*r2*r3 + 2*a1*a2**2*r1*r4 + 2*a3**2*a4*r1*r4 + 2*a3* &
      a4**2*r2*r3
f(617) = 2*a1**3*a2*r2*r3 + 2*a1*a2**3*r1*r4 + 2*a3**3*a4*r1*r4 + 2*a3* &
      a4**3*r2*r3
f(618) = 2*a1**2*a2*r1*r4 + 2*a1*a2**2*r2*r3 + 2*a3**2*a4*r2*r3 + 2*a3* &
      a4**2*r1*r4
f(619) = 2*a1**2*a2**2*r1*r4 + 2*a1**2*a2**2*r2*r3 + 2*a3**2*a4**2*r1*r4 &
      + 2*a3**2*a4**2*r2*r3
f(620) = 2*a1**3*a2*r1*r4 + 2*a1*a2**3*r2*r3 + 2*a3**3*a4*r2*r3 + 2*a3* &
      a4**3*r1*r4
f(621) = 2*a1*a2*r1*r4**2 + 2*a1*a2*r2*r3**2 + 2*a3*a4*r1**2*r4 + 2*a3* &
      a4*r2**2*r3
f(622) = 2*a1**2*a2*r2*r3**2 + 2*a1*a2**2*r1*r4**2 + 2*a3**2*a4*r1**2*r4 &
      + 2*a3*a4**2*r2**2*r3
f(623) = 2*a1**2*a2*r1*r4**2 + 2*a1*a2**2*r2*r3**2 + 2*a3**2*a4*r2**2*r3 &
      + 2*a3*a4**2*r1**2*r4
f(624) = 2*a1*a2*r1*r4**3 + 2*a1*a2*r2*r3**3 + 2*a3*a4*r1**3*r4 + 2*a3* &
      a4*r2**3*r3
f(625) = 2*a1*a2*r1**2*r4 + 2*a1*a2*r2**2*r3 + 2*a3*a4*r1*r4**2 + 2*a3* &
      a4*r2*r3**2
f(626) = 2*a1**2*a2*r2**2*r3 + 2*a1*a2**2*r1**2*r4 + 2*a3**2*a4*r1*r4**2 &
      + 2*a3*a4**2*r2*r3**2
f(627) = 2*a1**2*a2*r1**2*r4 + 2*a1*a2**2*r2**2*r3 + 2*a3**2*a4*r2*r3**2 &
      + 2*a3*a4**2*r1*r4**2
f(628) = 2*a1*a2*r1**2*r4**2 + 2*a1*a2*r2**2*r3**2 + 2*a3*a4*r1**2*r4**2 &
      + 2*a3*a4*r2**2*r3**2
f(629) = 2*a1*a2*r1**3*r4 + 2*a1*a2*r2**3*r3 + 2*a3*a4*r1*r4**3 + 2*a3* &
      a4*r2*r3**3
f(630) = 2*a1*a3*r1*r4 + 2*a1*a3*r2*r3 + 2*a2*a4*r1*r4 + 2*a2*a4*r2*r3
f(631) = 2*a1**2*a3*r2*r3 + 2*a1*a3**2*r1*r4 + 2*a2**2*a4*r1*r4 + 2*a2* &
      a4**2*r2*r3
f(632) = 2*a1**3*a3*r2*r3 + 2*a1*a3**3*r1*r4 + 2*a2**3*a4*r1*r4 + 2*a2* &
      a4**3*r2*r3
f(633) = 2*a1**2*a3*r1*r4 + 2*a1*a3**2*r2*r3 + 2*a2**2*a4*r2*r3 + 2*a2* &
      a4**2*r1*r4
f(634) = 2*a1**2*a3**2*r1*r4 + 2*a1**2*a3**2*r2*r3 + 2*a2**2*a4**2*r1*r4 &
      + 2*a2**2*a4**2*r2*r3
f(635) = 2*a1**3*a3*r1*r4 + 2*a1*a3**3*r2*r3 + 2*a2**3*a4*r2*r3 + 2*a2* &
      a4**3*r1*r4
f(636) = 2*a1*a3*r1*r4**2 + 2*a1*a3*r2**2*r3 + 2*a2*a4*r1**2*r4 + 2*a2* &
      a4*r2*r3**2
f(637) = 2*a1**2*a3*r2**2*r3 + 2*a1*a3**2*r1*r4**2 + 2*a2**2*a4*r1**2*r4 &
      + 2*a2*a4**2*r2*r3**2
f(638) = 2*a1**2*a3*r1*r4**2 + 2*a1*a3**2*r2**2*r3 + 2*a2**2*a4*r2*r3**2 &
      + 2*a2*a4**2*r1**2*r4
f(639) = 2*a1*a3*r1*r4**3 + 2*a1*a3*r2**3*r3 + 2*a2*a4*r1**3*r4 + 2*a2* &
      a4*r2*r3**3
f(640) = 2*a1*a3*r1**2*r4 + 2*a1*a3*r2*r3**2 + 2*a2*a4*r1*r4**2 + 2*a2* &
      a4*r2**2*r3
f(641) = 2*a1**2*a3*r2*r3**2 + 2*a1*a3**2*r1**2*r4 + 2*a2**2*a4*r1*r4**2 &
      + 2*a2*a4**2*r2**2*r3
f(642) = 2*a1**2*a3*r1**2*r4 + 2*a1*a3**2*r2*r3**2 + 2*a2**2*a4*r2**2*r3 &
      + 2*a2*a4**2*r1*r4**2
f(643) = 2*a1*a3*r1**2*r4**2 + 2*a1*a3*r2**2*r3**2 + 2*a2*a4*r1**2*r4**2 &
      + 2*a2*a4*r2**2*r3**2
f(644) = 2*a1*a3*r1**3*r4 + 2*a1*a3*r2*r3**3 + 2*a2*a4*r1*r4**3 + 2*a2* &
      a4*r2**3*r3
f(645) = 4*a1*a4*r1*r4 + 4*a2*a3*r2*r3
f(646) = 2*a1**2*a4*r1*r4 + 2*a1*a4**2*r1*r4 + 2*a2**2*a3*r2*r3 + 2*a2* &
      a3**2*r2*r3
f(647) = 2*a1**3*a4*r1*r4 + 2*a1*a4**3*r1*r4 + 2*a2**3*a3*r2*r3 + 2*a2* &
      a3**3*r2*r3
f(648) = 4*a1**2*a4**2*r1*r4 + 4*a2**2*a3**2*r2*r3
f(649) = 2*a1*a4*r1**2*r4 + 2*a1*a4*r1*r4**2 + 2*a2*a3*r2**2*r3 + 2*a2* &
      a3*r2*r3**2
f(650) = 2*a1**2*a4*r1**2*r4 + 2*a1*a4**2*r1*r4**2 + 2*a2**2*a3*r2**2*r3 &
      + 2*a2*a3**2*r2*r3**2
f(651) = 2*a1**2*a4*r1*r4**2 + 2*a1*a4**2*r1**2*r4 + 2*a2**2*a3*r2*r3**2 &
      + 2*a2*a3**2*r2**2*r3
f(652) = 2*a1*a4*r1**3*r4 + 2*a1*a4*r1*r4**3 + 2*a2*a3*r2**3*r3 + 2*a2* &
      a3*r2*r3**3
f(653) = 4*a1*a4*r1**2*r4**2 + 4*a2*a3*r2**2*r3**2
f(654) = 2*a1*b1**2*r1*r4 + 2*a2*b1**2*r2*r3 + 2*a3*b2**2*r2*r3 + 2*a4* &
      b2**2*r1*r4
f(655) = 2*a1**2*b1**2*r1*r4 + 2*a2**2*b1**2*r2*r3 + 2*a3**2*b2**2*r2*r3 &
      + 2*a4**2*b2**2*r1*r4
f(656) = 2*a1*b1**2*r1*r4**2 + 2*a2*b1**2*r2*r3**2 + 2*a3*b2**2*r2**2*r3 &
      + 2*a4*b2**2*r1**2*r4
f(657) = 2*a1*b1**2*r1**2*r4 + 2*a2*b1**2*r2**2*r3 + 2*a3*b2**2*r2*r3**2 &
      + 2*a4*b2**2*r1*r4**2
f(658) = 2*a1*b2**2*r1*r4 + 2*a2*b2**2*r2*r3 + 2*a3*b1**2*r2*r3 + 2*a4* &
      b1**2*r1*r4
f(659) = 2*a1**2*b2**2*r1*r4 + 2*a2**2*b2**2*r2*r3 + 2*a3**2*b1**2*r2*r3 &
      + 2*a4**2*b1**2*r1*r4
f(660) = 2*a1*b2**2*r1*r4**2 + 2*a2*b2**2*r2*r3**2 + 2*a3*b1**2*r2**2*r3 &
      + 2*a4*b1**2*r1**2*r4
f(661) = 2*a1*b2**2*r1**2*r4 + 2*a2*b2**2*r2**2*r3 + 2*a3*b1**2*r2*r3**2 &
      + 2*a4*b1**2*r1*r4**2
f(662) = 2*dtau**2*(a1*r1*r4 + a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(663) = 2*dtau**2*(a1**2*r1*r4 + a2**2*r2*r3 + a3**2*r2*r3 + a4**2*r1* &
      r4)
f(664) = 2*dtau**2*(a1*r1*r4**2 + a2*r2*r3**2 + a3*r2**2*r3 + a4*r1**2* &
      r4)
f(665) = 2*dtau**2*(a1*r1**2*r4 + a2*r2**2*r3 + a3*r2*r3**2 + a4*r1*r4** &
      2)
f(666) = 4*a1*a4*r2*r3 + 4*a2*a3*r1*r4
f(667) = 2*a1**2*a4*r2*r3 + 2*a1*a4**2*r2*r3 + 2*a2**2*a3*r1*r4 + 2*a2* &
      a3**2*r1*r4
f(668) = 2*a1**3*a4*r2*r3 + 2*a1*a4**3*r2*r3 + 2*a2**3*a3*r1*r4 + 2*a2* &
      a3**3*r1*r4
f(669) = 4*a1**2*a4**2*r2*r3 + 4*a2**2*a3**2*r1*r4
f(670) = 2*a1*a4*r2**2*r3 + 2*a1*a4*r2*r3**2 + 2*a2*a3*r1**2*r4 + 2*a2* &
      a3*r1*r4**2
f(671) = 2*a1**2*a4*r2**2*r3 + 2*a1*a4**2*r2*r3**2 + 2*a2**2*a3*r1**2*r4 &
      + 2*a2*a3**2*r1*r4**2
f(672) = 2*a1**2*a4*r2*r3**2 + 2*a1*a4**2*r2**2*r3 + 2*a2**2*a3*r1*r4**2 &
      + 2*a2*a3**2*r1**2*r4
f(673) = 2*a1*a4*r2**3*r3 + 2*a1*a4*r2*r3**3 + 2*a2*a3*r1**3*r4 + 2*a2* &
      a3*r1*r4**3
f(674) = 4*a1*a4*r2**2*r3**2 + 4*a2*a3*r1**2*r4**2
f(675) = 2*a1*b1**2*r2*r3 + 2*a2*b1**2*r1*r4 + 2*a3*b2**2*r1*r4 + 2*a4* &
      b2**2*r2*r3
f(676) = 2*a1**2*b1**2*r2*r3 + 2*a2**2*b1**2*r1*r4 + 2*a3**2*b2**2*r1*r4 &
      + 2*a4**2*b2**2*r2*r3
f(677) = 2*a1*b1**2*r2*r3**2 + 2*a2*b1**2*r1*r4**2 + 2*a3*b2**2*r1**2*r4 &
      + 2*a4*b2**2*r2**2*r3
f(678) = 2*a1*b1**2*r2**2*r3 + 2*a2*b1**2*r1**2*r4 + 2*a3*b2**2*r1*r4**2 &
      + 2*a4*b2**2*r2*r3**2
f(679) = 2*a1*b2**2*r2*r3 + 2*a2*b2**2*r1*r4 + 2*a3*b1**2*r1*r4 + 2*a4* &
      b1**2*r2*r3
f(680) = 2*a1**2*b2**2*r2*r3 + 2*a2**2*b2**2*r1*r4 + 2*a3**2*b1**2*r1*r4 &
      + 2*a4**2*b1**2*r2*r3
f(681) = 2*a1*b2**2*r2*r3**2 + 2*a2*b2**2*r1*r4**2 + 2*a3*b1**2*r1**2*r4 &
      + 2*a4*b1**2*r2**2*r3
f(682) = 2*a1*b2**2*r2**2*r3 + 2*a2*b2**2*r1**2*r4 + 2*a3*b1**2*r1*r4**2 &
      + 2*a4*b1**2*r2*r3**2
f(683) = 2*dtau**2*(a1*r2*r3 + a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(684) = 2*dtau**2*(a1**2*r2*r3 + a2**2*r1*r4 + a3**2*r1*r4 + a4**2*r2* &
      r3)
f(685) = 2*dtau**2*(a1*r2*r3**2 + a2*r1*r4**2 + a3*r1**2*r4 + a4*r2**2* &
      r3)
f(686) = 2*dtau**2*(a1*r2**2*r3 + a2*r1**2*r4 + a3*r1*r4**2 + a4*r2*r3** &
      2)
f(687) = 4*b1*b2*(r1*r4 + r2*r3)
f(688) = 2*b1*b2*(b1**2*r1*r4 + b1**2*r2*r3 + b2**2*r1*r4 + b2**2*r2*r3)
f(689) = 4*b1**2*b2**2*(r1*r4 + r2*r3)
f(690) = 2*b1*b2*(r1**2*r4 + r1*r4**2 + r2**2*r3 + r2*r3**2)
f(691) = 2*b1*b2*(r1**3*r4 + r1*r4**3 + r2**3*r3 + r2*r3**3)
f(692) = 4*b1*b2*(r1**2*r4**2 + r2**2*r3**2)
f(693) = 2*dtau*(b1*r1*r4 - b1*r2*r3 - b2*r1*r4 + b2*r2*r3)
f(694) = 2*dtau**3*(b1*r1*r4 - b1*r2*r3 - b2*r1*r4 + b2*r2*r3)
f(695) = 2*dtau**2*(b1**2*r1*r4 + b1**2*r2*r3 + b2**2*r1*r4 + b2**2*r2* &
      r3)
f(696) = 2*dtau*(b1**3*r1*r4 - b1**3*r2*r3 - b2**3*r1*r4 + b2**3*r2*r3)
f(697) = 2*dtau*(b1*r1*r4**2 - b1*r2*r3**2 - b2*r1**2*r4 + b2*r2**2*r3)
f(698) = 2*dtau*(b1*r1*r4**3 - b1*r2*r3**3 - b2*r1**3*r4 + b2*r2**3*r3)
f(699) = 2*dtau*(b1*r1**2*r4 - b1*r2**2*r3 - b2*r1*r4**2 + b2*r2*r3**2)
f(700) = 2*dtau*(b1*r1**2*r4**2 - b1*r2**2*r3**2 - b2*r1**2*r4**2 + b2* &
      r2**2*r3**2)
f(701) = 2*dtau*(b1*r1**3*r4 - b1*r2**3*r3 - b2*r1*r4**3 + b2*r2*r3**3)
f(702) = 2*a1*a2*a3*r1 + 2*a1*a2*a4*r2 + 2*a1*a3*a4*r3 + 2*a2*a3*a4*r4
f(703) = 2*a1**2*a3*a4*r3 + 2*a1*a2*a3**2*r1 + 2*a1*a2*a4**2*r2 + 2*a2** &
      2*a3*a4*r4
f(704) = 2*a1**3*a3*a4*r3 + 2*a1*a2*a3**3*r1 + 2*a1*a2*a4**3*r2 + 2*a2** &
      3*a3*a4*r4
f(705) = 2*a1**2*a2*a4*r2 + 2*a1*a2**2*a3*r1 + 2*a1*a3*a4**2*r3 + 2*a2* &
      a3**2*a4*r4
f(706) = 2*a1**2*a2*a4**2*r2 + 2*a1**2*a3*a4**2*r3 + 2*a1*a2**2*a3**2*r1 &
      + 2*a2**2*a3**2*a4*r4
f(707) = 2*a1**3*a2*a4*r2 + 2*a1*a2**3*a3*r1 + 2*a1*a3*a4**3*r3 + 2*a2* &
      a3**3*a4*r4
f(708) = 2*a1**2*a2*a3*r1 + 2*a1*a2**2*a4*r2 + 2*a1*a3**2*a4*r3 + 2*a2* &
      a3*a4**2*r4
f(709) = 2*a1**2*a2*a3**2*r1 + 2*a1**2*a3**2*a4*r3 + 2*a1*a2**2*a4**2*r2 &
      + 2*a2**2*a3*a4**2*r4
f(710) = 2*a1**2*a2**2*a3*r1 + 2*a1**2*a2**2*a4*r2 + 2*a1*a3**2*a4**2*r3 &
      + 2*a2*a3**2*a4**2*r4
f(711) = 2*a1**3*a2*a3*r1 + 2*a1*a2**3*a4*r2 + 2*a1*a3**3*a4*r3 + 2*a2* &
      a3*a4**3*r4
f(712) = 2*a1*a2*a3*r1**2 + 2*a1*a2*a4*r2**2 + 2*a1*a3*a4*r3**2 + 2*a2* &
      a3*a4*r4**2
f(713) = 2*a1**2*a3*a4*r3**2 + 2*a1*a2*a3**2*r1**2 + 2*a1*a2*a4**2*r2**2 &
      + 2*a2**2*a3*a4*r4**2
f(714) = 2*a1**2*a2*a4*r2**2 + 2*a1*a2**2*a3*r1**2 + 2*a1*a3*a4**2*r3**2 &
      + 2*a2*a3**2*a4*r4**2
f(715) = 2*a1**2*a2*a3*r1**2 + 2*a1*a2**2*a4*r2**2 + 2*a1*a3**2*a4*r3**2 &
      + 2*a2*a3*a4**2*r4**2
f(716) = 2*a1*a2*a3*r1**3 + 2*a1*a2*a4*r2**3 + 2*a1*a3*a4*r3**3 + 2*a2* &
      a3*a4*r4**3
f(717) = 2*a1*a2*a3*r2 + 2*a1*a2*a4*r1 + 2*a1*a3*a4*r4 + 2*a2*a3*a4*r3
f(718) = 2*a1**2*a3*a4*r4 + 2*a1*a2*a3**2*r2 + 2*a1*a2*a4**2*r1 + 2*a2** &
      2*a3*a4*r3
f(719) = 2*a1**3*a3*a4*r4 + 2*a1*a2*a3**3*r2 + 2*a1*a2*a4**3*r1 + 2*a2** &
      3*a3*a4*r3
f(720) = 2*a1**2*a2*a3*r2 + 2*a1*a2**2*a4*r1 + 2*a1*a3**2*a4*r4 + 2*a2* &
      a3*a4**2*r3
f(721) = 2*a1**2*a2*a3**2*r2 + 2*a1**2*a3**2*a4*r4 + 2*a1*a2**2*a4**2*r1 &
      + 2*a2**2*a3*a4**2*r3
f(722) = 2*a1**3*a2*a3*r2 + 2*a1*a2**3*a4*r1 + 2*a1*a3**3*a4*r4 + 2*a2* &
      a3*a4**3*r3
f(723) = 2*a1**2*a2*a4*r1 + 2*a1*a2**2*a3*r2 + 2*a1*a3*a4**2*r4 + 2*a2* &
      a3**2*a4*r3
f(724) = 2*a1**2*a2*a4**2*r1 + 2*a1**2*a3*a4**2*r4 + 2*a1*a2**2*a3**2*r2 &
      + 2*a2**2*a3**2*a4*r3
f(725) = 2*a1**2*a2**2*a3*r2 + 2*a1**2*a2**2*a4*r1 + 2*a1*a3**2*a4**2*r4 &
      + 2*a2*a3**2*a4**2*r3
f(726) = 2*a1**3*a2*a4*r1 + 2*a1*a2**3*a3*r2 + 2*a1*a3*a4**3*r4 + 2*a2* &
      a3**3*a4*r3
f(727) = 2*a1*a2*a3*r2**2 + 2*a1*a2*a4*r1**2 + 2*a1*a3*a4*r4**2 + 2*a2* &
      a3*a4*r3**2
f(728) = 2*a1**2*a3*a4*r4**2 + 2*a1*a2*a3**2*r2**2 + 2*a1*a2*a4**2*r1**2 &
      + 2*a2**2*a3*a4*r3**2
f(729) = 2*a1**2*a2*a3*r2**2 + 2*a1*a2**2*a4*r1**2 + 2*a1*a3**2*a4*r4**2 &
      + 2*a2*a3*a4**2*r3**2
f(730) = 2*a1**2*a2*a4*r1**2 + 2*a1*a2**2*a3*r2**2 + 2*a1*a3*a4**2*r4**2 &
      + 2*a2*a3**2*a4*r3**2
f(731) = 2*a1*a2*a3*r2**3 + 2*a1*a2*a4*r1**3 + 2*a1*a3*a4*r4**3 + 2*a2* &
      a3*a4*r3**3
f(732) = 2*a1*a2*b1**2*r1 + 2*a1*a2*b1**2*r2 + 2*a3*a4*b2**2*r3 + 2*a3* &
      a4*b2**2*r4
f(733) = 2*a1**2*a2*b1**2*r2 + 2*a1*a2**2*b1**2*r1 + 2*a3**2*a4*b2**2*r4 &
      + 2*a3*a4**2*b2**2*r3
f(734) = 2*a1**2*a2*b1**2*r1 + 2*a1*a2**2*b1**2*r2 + 2*a3**2*a4*b2**2*r3 &
      + 2*a3*a4**2*b2**2*r4
f(735) = 2*a1*a2*b1**2*r1**2 + 2*a1*a2*b1**2*r2**2 + 2*a3*a4*b2**2*r3**2 &
      + 2*a3*a4*b2**2*r4**2
f(736) = 2*a1*a2*b2**2*r1 + 2*a1*a2*b2**2*r2 + 2*a3*a4*b1**2*r3 + 2*a3* &
      a4*b1**2*r4
f(737) = 2*a1**2*a2*b2**2*r2 + 2*a1*a2**2*b2**2*r1 + 2*a3**2*a4*b1**2*r4 &
      + 2*a3*a4**2*b1**2*r3
f(738) = 2*a1**2*a2*b2**2*r1 + 2*a1*a2**2*b2**2*r2 + 2*a3**2*a4*b1**2*r3 &
      + 2*a3*a4**2*b1**2*r4
f(739) = 2*a1*a2*b2**2*r1**2 + 2*a1*a2*b2**2*r2**2 + 2*a3*a4*b1**2*r3**2 &
      + 2*a3*a4*b1**2*r4**2
f(740) = 2*dtau**2*(a1*a2*r1 + a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(741) = 2*dtau**2*(a1**2*a2*r2 + a1*a2**2*r1 + a3**2*a4*r4 + a3*a4**2* &
      r3)
f(742) = 2*dtau**2*(a1**2*a2*r1 + a1*a2**2*r2 + a3**2*a4*r3 + a3*a4**2* &
      r4)
f(743) = 2*dtau**2*(a1*a2*r1**2 + a1*a2*r2**2 + a3*a4*r3**2 + a3*a4*r4** &
      2)
f(744) = 2*a1*a2*a3*r3 + 2*a1*a2*a4*r4 + 2*a1*a3*a4*r1 + 2*a2*a3*a4*r2
f(745) = 2*a1**2*a2*a4*r4 + 2*a1*a2**2*a3*r3 + 2*a1*a3*a4**2*r1 + 2*a2* &
      a3**2*a4*r2
f(746) = 2*a1**3*a2*a4*r4 + 2*a1*a2**3*a3*r3 + 2*a1*a3*a4**3*r1 + 2*a2* &
      a3**3*a4*r2
f(747) = 2*a1**2*a2*a3*r3 + 2*a1*a2**2*a4*r4 + 2*a1*a3**2*a4*r1 + 2*a2* &
      a3*a4**2*r2
f(748) = 2*a1**2*a2**2*a3*r3 + 2*a1**2*a2**2*a4*r4 + 2*a1*a3**2*a4**2*r1 &
      + 2*a2*a3**2*a4**2*r2
f(749) = 2*a1**3*a2*a3*r3 + 2*a1*a2**3*a4*r4 + 2*a1*a3**3*a4*r1 + 2*a2* &
      a3*a4**3*r2
f(750) = 2*a1**2*a3*a4*r1 + 2*a1*a2*a3**2*r3 + 2*a1*a2*a4**2*r4 + 2*a2** &
      2*a3*a4*r2
f(751) = 2*a1**2*a2*a4**2*r4 + 2*a1**2*a3*a4**2*r1 + 2*a1*a2**2*a3**2*r3 &
      + 2*a2**2*a3**2*a4*r2
f(752) = 2*a1**2*a2*a3**2*r3 + 2*a1**2*a3**2*a4*r1 + 2*a1*a2**2*a4**2*r4 &
      + 2*a2**2*a3*a4**2*r2
f(753) = 2*a1**3*a3*a4*r1 + 2*a1*a2*a3**3*r3 + 2*a1*a2*a4**3*r4 + 2*a2** &
      3*a3*a4*r2
f(754) = 2*a1*a2*a3*r3**2 + 2*a1*a2*a4*r4**2 + 2*a1*a3*a4*r1**2 + 2*a2* &
      a3*a4*r2**2
f(755) = 2*a1**2*a2*a4*r4**2 + 2*a1*a2**2*a3*r3**2 + 2*a1*a3*a4**2*r1**2 &
      + 2*a2*a3**2*a4*r2**2
f(756) = 2*a1**2*a2*a3*r3**2 + 2*a1*a2**2*a4*r4**2 + 2*a1*a3**2*a4*r1**2 &
      + 2*a2*a3*a4**2*r2**2
f(757) = 2*a1**2*a3*a4*r1**2 + 2*a1*a2*a3**2*r3**2 + 2*a1*a2*a4**2*r4**2 &
      + 2*a2**2*a3*a4*r2**2
f(758) = 2*a1*a2*a3*r3**3 + 2*a1*a2*a4*r4**3 + 2*a1*a3*a4*r1**3 + 2*a2* &
      a3*a4*r2**3
f(759) = 2*a1*a3*b1**2*r1 + 2*a1*a3*b2**2*r3 + 2*a2*a4*b1**2*r2 + 2*a2* &
      a4*b2**2*r4
f(760) = 2*a1**2*a3*b2**2*r3 + 2*a1*a3**2*b1**2*r1 + 2*a2**2*a4*b2**2*r4 &
      + 2*a2*a4**2*b1**2*r2
f(761) = 2*a1**2*a3*b1**2*r1 + 2*a1*a3**2*b2**2*r3 + 2*a2**2*a4*b1**2*r2 &
      + 2*a2*a4**2*b2**2*r4
f(762) = 2*a1*a3*b1**2*r1**2 + 2*a1*a3*b2**2*r3**2 + 2*a2*a4*b1**2*r2**2 &
      + 2*a2*a4*b2**2*r4**2
f(763) = 2*a1*a3*b1**2*r3 + 2*a1*a3*b2**2*r1 + 2*a2*a4*b1**2*r4 + 2*a2* &
      a4*b2**2*r2
f(764) = 2*a1**2*a3*b1**2*r3 + 2*a1*a3**2*b2**2*r1 + 2*a2**2*a4*b1**2*r4 &
      + 2*a2*a4**2*b2**2*r2
f(765) = 2*a1**2*a3*b2**2*r1 + 2*a1*a3**2*b1**2*r3 + 2*a2**2*a4*b2**2*r2 &
      + 2*a2*a4**2*b1**2*r4
f(766) = 2*a1*a3*b1**2*r3**2 + 2*a1*a3*b2**2*r1**2 + 2*a2*a4*b1**2*r4**2 &
      + 2*a2*a4*b2**2*r2**2
f(767) = 2*dtau**2*(a1*a3*r1 + a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(768) = 2*dtau**2*(a1**2*a3*r3 + a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2* &
      r2)
f(769) = 2*dtau**2*(a1**2*a3*r1 + a1*a3**2*r3 + a2**2*a4*r2 + a2*a4**2* &
      r4)
f(770) = 2*dtau**2*(a1*a3*r1**2 + a1*a3*r3**2 + a2*a4*r2**2 + a2*a4*r4** &
      2)
f(771) = 2*a1*a4*b1**2*r1 + 2*a1*a4*b2**2*r4 + 2*a2*a3*b1**2*r2 + 2*a2* &
      a3*b2**2*r3
f(772) = 2*a1**2*a4*b2**2*r4 + 2*a1*a4**2*b1**2*r1 + 2*a2**2*a3*b2**2*r3 &
      + 2*a2*a3**2*b1**2*r2
f(773) = 2*a1**2*a4*b1**2*r1 + 2*a1*a4**2*b2**2*r4 + 2*a2**2*a3*b1**2*r2 &
      + 2*a2*a3**2*b2**2*r3
f(774) = 2*a1*a4*b1**2*r1**2 + 2*a1*a4*b2**2*r4**2 + 2*a2*a3*b1**2*r2**2 &
      + 2*a2*a3*b2**2*r3**2
f(775) = 2*a1*a4*b1**2*r4 + 2*a1*a4*b2**2*r1 + 2*a2*a3*b1**2*r3 + 2*a2* &
      a3*b2**2*r2
f(776) = 2*a1**2*a4*b1**2*r4 + 2*a1*a4**2*b2**2*r1 + 2*a2**2*a3*b1**2*r3 &
      + 2*a2*a3**2*b2**2*r2
f(777) = 2*a1**2*a4*b2**2*r1 + 2*a1*a4**2*b1**2*r4 + 2*a2**2*a3*b2**2*r2 &
      + 2*a2*a3**2*b1**2*r3
f(778) = 2*a1*a4*b1**2*r4**2 + 2*a1*a4*b2**2*r1**2 + 2*a2*a3*b1**2*r3**2 &
      + 2*a2*a3*b2**2*r2**2
f(779) = 2*dtau**2*(a1*a4*r1 + a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(780) = 2*dtau**2*(a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 + a2*a3**2* &
      r2)
f(781) = 2*dtau**2*(a1**2*a4*r1 + a1*a4**2*r4 + a2**2*a3*r2 + a2*a3**2* &
      r3)
f(782) = 2*dtau**2*(a1*a4*r1**2 + a1*a4*r4**2 + a2*a3*r2**2 + a2*a3*r3** &
      2)
f(783) = 2*b1*b2*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(784) = 2*b1*b2*(a1*b2**2*r1 + a2*b2**2*r2 + a3*b1**2*r3 + a4*b1**2*r4)
f(785) = 2*b1**2*b2**2*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(786) = 2*b1*b2*(a1*b1**2*r1 + a2*b1**2*r2 + a3*b2**2*r3 + a4*b2**2*r4)
f(787) = 2*b1*b2*(a1**2*r1 + a2**2*r2 + a3**2*r3 + a4**2*r4)
f(788) = 2*b1*b2*(a1**3*r1 + a2**3*r2 + a3**3*r3 + a4**3*r4)
f(789) = 2*b1*b2*(a1*r1**2 + a2*r2**2 + a3*r3**2 + a4*r4**2)
f(790) = 2*b1*b2*(a1**2*r1**2 + a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2)
f(791) = 2*b1*b2*(a1*r1**3 + a2*r2**3 + a3*r3**3 + a4*r4**3)
f(792) = 2*dtau*(a1*b1*r1 - a2*b1*r2 + a3*b2*r3 - a4*b2*r4)
f(793) = 2*dtau**3*(a1*b1*r1 - a2*b1*r2 + a3*b2*r3 - a4*b2*r4)
f(794) = 2*dtau**2*(a1*b1**2*r1 + a2*b1**2*r2 + a3*b2**2*r3 + a4*b2**2* &
      r4)
f(795) = 2*dtau*(a1*b1**3*r1 - a2*b1**3*r2 + a3*b2**3*r3 - a4*b2**3*r4)
f(796) = 2*dtau*(a1**2*b1*r1 - a2**2*b1*r2 + a3**2*b2*r3 - a4**2*b2*r4)
f(797) = 2*dtau*(a1**3*b1*r1 - a2**3*b1*r2 + a3**3*b2*r3 - a4**3*b2*r4)
f(798) = 2*dtau*(a1*b1*r1**2 - a2*b1*r2**2 + a3*b2*r3**2 - a4*b2*r4**2)
f(799) = 2*dtau*(a1**2*b1*r1**2 - a2**2*b1*r2**2 + a3**2*b2*r3**2 - a4** &
      2*b2*r4**2)
f(800) = 2*dtau*(a1*b1*r1**3 - a2*b1*r2**3 + a3*b2*r3**3 - a4*b2*r4**3)
f(801) = 2*dtau*(a1*b2*r1 - a2*b2*r2 + a3*b1*r3 - a4*b1*r4)
f(802) = 2*dtau**3*(a1*b2*r1 - a2*b2*r2 + a3*b1*r3 - a4*b1*r4)
f(803) = 2*dtau**2*(a1*b2**2*r1 + a2*b2**2*r2 + a3*b1**2*r3 + a4*b1**2* &
      r4)
f(804) = 2*dtau*(a1*b2**3*r1 - a2*b2**3*r2 + a3*b1**3*r3 - a4*b1**3*r4)
f(805) = 2*dtau*(a1**2*b2*r1 - a2**2*b2*r2 + a3**2*b1*r3 - a4**2*b1*r4)
f(806) = 2*dtau*(a1**3*b2*r1 - a2**3*b2*r2 + a3**3*b1*r3 - a4**3*b1*r4)
f(807) = 2*dtau*(a1*b2*r1**2 - a2*b2*r2**2 + a3*b1*r3**2 - a4*b1*r4**2)
f(808) = 2*dtau*(a1**2*b2*r1**2 - a2**2*b2*r2**2 + a3**2*b1*r3**2 - a4** &
      2*b1*r4**2)
f(809) = 2*dtau*(a1*b2*r1**3 - a2*b2*r2**3 + a3*b1*r3**3 - a4*b1*r4**3)
f(810) = 2*a1*a2*a3*r4 + 2*a1*a2*a4*r3 + 2*a1*a3*a4*r2 + 2*a2*a3*a4*r1
f(811) = 2*a1**2*a2*a3*r4 + 2*a1*a2**2*a4*r3 + 2*a1*a3**2*a4*r2 + 2*a2* &
      a3*a4**2*r1
f(812) = 2*a1**3*a2*a3*r4 + 2*a1*a2**3*a4*r3 + 2*a1*a3**3*a4*r2 + 2*a2* &
      a3*a4**3*r1
f(813) = 2*a1**2*a2*a4*r3 + 2*a1*a2**2*a3*r4 + 2*a1*a3*a4**2*r2 + 2*a2* &
      a3**2*a4*r1
f(814) = 2*a1**2*a2**2*a3*r4 + 2*a1**2*a2**2*a4*r3 + 2*a1*a3**2*a4**2*r2 &
      + 2*a2*a3**2*a4**2*r1
f(815) = 2*a1**3*a2*a4*r3 + 2*a1*a2**3*a3*r4 + 2*a1*a3*a4**3*r2 + 2*a2* &
      a3**3*a4*r1
f(816) = 2*a1**2*a3*a4*r2 + 2*a1*a2*a3**2*r4 + 2*a1*a2*a4**2*r3 + 2*a2** &
      2*a3*a4*r1
f(817) = 2*a1**2*a2*a3**2*r4 + 2*a1**2*a3**2*a4*r2 + 2*a1*a2**2*a4**2*r3 &
      + 2*a2**2*a3*a4**2*r1
f(818) = 2*a1**2*a2*a4**2*r3 + 2*a1**2*a3*a4**2*r2 + 2*a1*a2**2*a3**2*r4 &
      + 2*a2**2*a3**2*a4*r1
f(819) = 2*a1**3*a3*a4*r2 + 2*a1*a2*a3**3*r4 + 2*a1*a2*a4**3*r3 + 2*a2** &
      3*a3*a4*r1
f(820) = 2*a1*a2*a3*r4**2 + 2*a1*a2*a4*r3**2 + 2*a1*a3*a4*r2**2 + 2*a2* &
      a3*a4*r1**2
f(821) = 2*a1**2*a2*a3*r4**2 + 2*a1*a2**2*a4*r3**2 + 2*a1*a3**2*a4*r2**2 &
      + 2*a2*a3*a4**2*r1**2
f(822) = 2*a1**2*a2*a4*r3**2 + 2*a1*a2**2*a3*r4**2 + 2*a1*a3*a4**2*r2**2 &
      + 2*a2*a3**2*a4*r1**2
f(823) = 2*a1**2*a3*a4*r2**2 + 2*a1*a2*a3**2*r4**2 + 2*a1*a2*a4**2*r3**2 &
      + 2*a2**2*a3*a4*r1**2
f(824) = 2*a1*a2*a3*r4**3 + 2*a1*a2*a4*r3**3 + 2*a1*a3*a4*r2**3 + 2*a2* &
      a3*a4*r1**3
f(825) = 2*a1*a4*b1**2*r2 + 2*a1*a4*b2**2*r3 + 2*a2*a3*b1**2*r1 + 2*a2* &
      a3*b2**2*r4
f(826) = 2*a1**2*a4*b2**2*r3 + 2*a1*a4**2*b1**2*r2 + 2*a2**2*a3*b2**2*r4 &
      + 2*a2*a3**2*b1**2*r1
f(827) = 2*a1**2*a4*b1**2*r2 + 2*a1*a4**2*b2**2*r3 + 2*a2**2*a3*b1**2*r1 &
      + 2*a2*a3**2*b2**2*r4
f(828) = 2*a1*a4*b1**2*r2**2 + 2*a1*a4*b2**2*r3**2 + 2*a2*a3*b1**2*r1**2 &
      + 2*a2*a3*b2**2*r4**2
f(829) = 2*a1*a4*b1**2*r3 + 2*a1*a4*b2**2*r2 + 2*a2*a3*b1**2*r4 + 2*a2* &
      a3*b2**2*r1
f(830) = 2*a1**2*a4*b1**2*r3 + 2*a1*a4**2*b2**2*r2 + 2*a2**2*a3*b1**2*r4 &
      + 2*a2*a3**2*b2**2*r1
f(831) = 2*a1**2*a4*b2**2*r2 + 2*a1*a4**2*b1**2*r3 + 2*a2**2*a3*b2**2*r1 &
      + 2*a2*a3**2*b1**2*r4
f(832) = 2*a1*a4*b1**2*r3**2 + 2*a1*a4*b2**2*r2**2 + 2*a2*a3*b1**2*r4**2 &
      + 2*a2*a3*b2**2*r1**2
f(833) = 2*dtau**2*(a1*a4*r2 + a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(834) = 2*dtau**2*(a1**2*a4*r3 + a1*a4**2*r2 + a2**2*a3*r4 + a2*a3**2* &
      r1)
f(835) = 2*dtau**2*(a1**2*a4*r2 + a1*a4**2*r3 + a2**2*a3*r1 + a2*a3**2* &
      r4)
f(836) = 2*dtau**2*(a1*a4*r2**2 + a1*a4*r3**2 + a2*a3*r1**2 + a2*a3*r4** &
      2)
f(837) = 2*a1*a3*b1**2*r2 + 2*a1*a3*b2**2*r4 + 2*a2*a4*b1**2*r1 + 2*a2* &
      a4*b2**2*r3
f(838) = 2*a1**2*a3*b2**2*r4 + 2*a1*a3**2*b1**2*r2 + 2*a2**2*a4*b2**2*r3 &
      + 2*a2*a4**2*b1**2*r1
f(839) = 2*a1**2*a3*b1**2*r2 + 2*a1*a3**2*b2**2*r4 + 2*a2**2*a4*b1**2*r1 &
      + 2*a2*a4**2*b2**2*r3
f(840) = 2*a1*a3*b1**2*r2**2 + 2*a1*a3*b2**2*r4**2 + 2*a2*a4*b1**2*r1**2 &
      + 2*a2*a4*b2**2*r3**2
f(841) = 2*a1*a3*b1**2*r4 + 2*a1*a3*b2**2*r2 + 2*a2*a4*b1**2*r3 + 2*a2* &
      a4*b2**2*r1
f(842) = 2*a1**2*a3*b1**2*r4 + 2*a1*a3**2*b2**2*r2 + 2*a2**2*a4*b1**2*r3 &
      + 2*a2*a4**2*b2**2*r1
f(843) = 2*a1**2*a3*b2**2*r2 + 2*a1*a3**2*b1**2*r4 + 2*a2**2*a4*b2**2*r1 &
      + 2*a2*a4**2*b1**2*r3
f(844) = 2*a1*a3*b1**2*r4**2 + 2*a1*a3*b2**2*r2**2 + 2*a2*a4*b1**2*r3**2 &
      + 2*a2*a4*b2**2*r1**2
f(845) = 2*dtau**2*(a1*a3*r2 + a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(846) = 2*dtau**2*(a1**2*a3*r4 + a1*a3**2*r2 + a2**2*a4*r3 + a2*a4**2* &
      r1)
f(847) = 2*dtau**2*(a1**2*a3*r2 + a1*a3**2*r4 + a2**2*a4*r1 + a2*a4**2* &
      r3)
f(848) = 2*dtau**2*(a1*a3*r2**2 + a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3** &
      2)
f(849) = 2*b1*b2*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(850) = 2*b1*b2*(a1*b2**2*r2 + a2*b2**2*r1 + a3*b1**2*r4 + a4*b1**2*r3)
f(851) = 2*b1**2*b2**2*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(852) = 2*b1*b2*(a1*b1**2*r2 + a2*b1**2*r1 + a3*b2**2*r4 + a4*b2**2*r3)
f(853) = 2*b1*b2*(a1**2*r2 + a2**2*r1 + a3**2*r4 + a4**2*r3)
f(854) = 2*b1*b2*(a1**3*r2 + a2**3*r1 + a3**3*r4 + a4**3*r3)
f(855) = 2*b1*b2*(a1*r2**2 + a2*r1**2 + a3*r4**2 + a4*r3**2)
f(856) = 2*b1*b2*(a1**2*r2**2 + a2**2*r1**2 + a3**2*r4**2 + a4**2*r3**2)
f(857) = 2*b1*b2*(a1*r2**3 + a2*r1**3 + a3*r4**3 + a4*r3**3)
f(858) = 2*dtau*(-a1*b1*r2 + a2*b1*r1 - a3*b2*r4 + a4*b2*r3)
f(859) = 2*dtau**3*(-a1*b1*r2 + a2*b1*r1 - a3*b2*r4 + a4*b2*r3)
f(860) = 2*dtau**2*(a1*b1**2*r2 + a2*b1**2*r1 + a3*b2**2*r4 + a4*b2**2* &
      r3)
f(861) = 2*dtau*(-a1*b1**3*r2 + a2*b1**3*r1 - a3*b2**3*r4 + a4*b2**3*r3)
f(862) = 2*dtau*(-a1**2*b1*r2 + a2**2*b1*r1 - a3**2*b2*r4 + a4**2*b2*r3)
f(863) = 2*dtau*(-a1**3*b1*r2 + a2**3*b1*r1 - a3**3*b2*r4 + a4**3*b2*r3)
f(864) = 2*dtau*(-a1*b1*r2**2 + a2*b1*r1**2 - a3*b2*r4**2 + a4*b2*r3**2)
f(865) = 2*dtau*(-a1**2*b1*r2**2 + a2**2*b1*r1**2 - a3**2*b2*r4**2 + a4 &
      **2*b2*r3**2)
f(866) = 2*dtau*(-a1*b1*r2**3 + a2*b1*r1**3 - a3*b2*r4**3 + a4*b2*r3**3)
f(867) = 2*dtau*(-a1*b2*r2 + a2*b2*r1 - a3*b1*r4 + a4*b1*r3)
f(868) = 2*dtau**3*(-a1*b2*r2 + a2*b2*r1 - a3*b1*r4 + a4*b1*r3)
f(869) = 2*dtau**2*(a1*b2**2*r2 + a2*b2**2*r1 + a3*b1**2*r4 + a4*b1**2* &
      r3)
f(870) = 2*dtau*(-a1*b2**3*r2 + a2*b2**3*r1 - a3*b1**3*r4 + a4*b1**3*r3)
f(871) = 2*dtau*(-a1**2*b2*r2 + a2**2*b2*r1 - a3**2*b1*r4 + a4**2*b1*r3)
f(872) = 2*dtau*(-a1**3*b2*r2 + a2**3*b2*r1 - a3**3*b1*r4 + a4**3*b1*r3)
f(873) = 2*dtau*(-a1*b2*r2**2 + a2*b2*r1**2 - a3*b1*r4**2 + a4*b1*r3**2)
f(874) = 2*dtau*(-a1**2*b2*r2**2 + a2**2*b2*r1**2 - a3**2*b1*r4**2 + a4 &
      **2*b1*r3**2)
f(875) = 2*dtau*(-a1*b2*r2**3 + a2*b2*r1**3 - a3*b1*r4**3 + a4*b1*r3**3)
f(876) = 2*a1*a2*b2**2*r3 + 2*a1*a2*b2**2*r4 + 2*a3*a4*b1**2*r1 + 2*a3* &
      a4*b1**2*r2
f(877) = 2*a1**2*a2*b2**2*r4 + 2*a1*a2**2*b2**2*r3 + 2*a3**2*a4*b1**2*r2 &
      + 2*a3*a4**2*b1**2*r1
f(878) = 2*a1**2*a2*b2**2*r3 + 2*a1*a2**2*b2**2*r4 + 2*a3**2*a4*b1**2*r1 &
      + 2*a3*a4**2*b1**2*r2
f(879) = 2*a1*a2*b2**2*r3**2 + 2*a1*a2*b2**2*r4**2 + 2*a3*a4*b1**2*r1**2 &
      + 2*a3*a4*b1**2*r2**2
f(880) = 2*a1*a2*b1**2*r3 + 2*a1*a2*b1**2*r4 + 2*a3*a4*b2**2*r1 + 2*a3* &
      a4*b2**2*r2
f(881) = 2*a1**2*a2*b1**2*r4 + 2*a1*a2**2*b1**2*r3 + 2*a3**2*a4*b2**2*r2 &
      + 2*a3*a4**2*b2**2*r1
f(882) = 2*a1**2*a2*b1**2*r3 + 2*a1*a2**2*b1**2*r4 + 2*a3**2*a4*b2**2*r1 &
      + 2*a3*a4**2*b2**2*r2
f(883) = 2*a1*a2*b1**2*r3**2 + 2*a1*a2*b1**2*r4**2 + 2*a3*a4*b2**2*r1**2 &
      + 2*a3*a4*b2**2*r2**2
f(884) = 2*dtau**2*(a1*a2*r3 + a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(885) = 2*dtau**2*(a1**2*a2*r4 + a1*a2**2*r3 + a3**2*a4*r2 + a3*a4**2* &
      r1)
f(886) = 2*dtau**2*(a1**2*a2*r3 + a1*a2**2*r4 + a3**2*a4*r1 + a3*a4**2* &
      r2)
f(887) = 2*dtau**2*(a1*a2*r3**2 + a1*a2*r4**2 + a3*a4*r1**2 + a3*a4*r2** &
      2)
f(888) = 2*b1*b2*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(889) = 2*b1*b2*(a1*b1**2*r3 + a2*b1**2*r4 + a3*b2**2*r1 + a4*b2**2*r2)
f(890) = 2*b1**2*b2**2*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(891) = 2*b1*b2*(a1*b2**2*r3 + a2*b2**2*r4 + a3*b1**2*r1 + a4*b1**2*r2)
f(892) = 2*b1*b2*(a1**2*r3 + a2**2*r4 + a3**2*r1 + a4**2*r2)
f(893) = 2*b1*b2*(a1**3*r3 + a2**3*r4 + a3**3*r1 + a4**3*r2)
f(894) = 2*b1*b2*(a1*r3**2 + a2*r4**2 + a3*r1**2 + a4*r2**2)
f(895) = 2*b1*b2*(a1**2*r3**2 + a2**2*r4**2 + a3**2*r1**2 + a4**2*r2**2)
f(896) = 2*b1*b2*(a1*r3**3 + a2*r4**3 + a3*r1**3 + a4*r2**3)
f(897) = 2*dtau*(a1*b2*r3 - a2*b2*r4 + a3*b1*r1 - a4*b1*r2)
f(898) = 2*dtau**3*(a1*b2*r3 - a2*b2*r4 + a3*b1*r1 - a4*b1*r2)
f(899) = 2*dtau**2*(a1*b2**2*r3 + a2*b2**2*r4 + a3*b1**2*r1 + a4*b1**2* &
      r2)
f(900) = 2*dtau*(a1*b2**3*r3 - a2*b2**3*r4 + a3*b1**3*r1 - a4*b1**3*r2)
f(901) = 2*dtau*(a1**2*b2*r3 - a2**2*b2*r4 + a3**2*b1*r1 - a4**2*b1*r2)
f(902) = 2*dtau*(a1**3*b2*r3 - a2**3*b2*r4 + a3**3*b1*r1 - a4**3*b1*r2)
f(903) = 2*dtau*(a1*b2*r3**2 - a2*b2*r4**2 + a3*b1*r1**2 - a4*b1*r2**2)
f(904) = 2*dtau*(a1**2*b2*r3**2 - a2**2*b2*r4**2 + a3**2*b1*r1**2 - a4** &
      2*b1*r2**2)
f(905) = 2*dtau*(a1*b2*r3**3 - a2*b2*r4**3 + a3*b1*r1**3 - a4*b1*r2**3)
f(906) = 2*dtau*(a1*b1*r3 - a2*b1*r4 + a3*b2*r1 - a4*b2*r2)
f(907) = 2*dtau**3*(a1*b1*r3 - a2*b1*r4 + a3*b2*r1 - a4*b2*r2)
f(908) = 2*dtau**2*(a1*b1**2*r3 + a2*b1**2*r4 + a3*b2**2*r1 + a4*b2**2* &
      r2)
f(909) = 2*dtau*(a1*b1**3*r3 - a2*b1**3*r4 + a3*b2**3*r1 - a4*b2**3*r2)
f(910) = 2*dtau*(a1**2*b1*r3 - a2**2*b1*r4 + a3**2*b2*r1 - a4**2*b2*r2)
f(911) = 2*dtau*(a1**3*b1*r3 - a2**3*b1*r4 + a3**3*b2*r1 - a4**3*b2*r2)
f(912) = 2*dtau*(a1*b1*r3**2 - a2*b1*r4**2 + a3*b2*r1**2 - a4*b2*r2**2)
f(913) = 2*dtau*(a1**2*b1*r3**2 - a2**2*b1*r4**2 + a3**2*b2*r1**2 - a4** &
      2*b2*r2**2)
f(914) = 2*dtau*(a1*b1*r3**3 - a2*b1*r4**3 + a3*b2*r1**3 - a4*b2*r2**3)
f(915) = 2*b1*b2*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(916) = 2*b1*b2*(a1*b1**2*r4 + a2*b1**2*r3 + a3*b2**2*r2 + a4*b2**2*r1)
f(917) = 2*b1**2*b2**2*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(918) = 2*b1*b2*(a1*b2**2*r4 + a2*b2**2*r3 + a3*b1**2*r2 + a4*b1**2*r1)
f(919) = 2*b1*b2*(a1**2*r4 + a2**2*r3 + a3**2*r2 + a4**2*r1)
f(920) = 2*b1*b2*(a1**3*r4 + a2**3*r3 + a3**3*r2 + a4**3*r1)
f(921) = 2*b1*b2*(a1*r4**2 + a2*r3**2 + a3*r2**2 + a4*r1**2)
f(922) = 2*b1*b2*(a1**2*r4**2 + a2**2*r3**2 + a3**2*r2**2 + a4**2*r1**2)
f(923) = 2*b1*b2*(a1*r4**3 + a2*r3**3 + a3*r2**3 + a4*r1**3)
f(924) = 2*dtau*(-a1*b2*r4 + a2*b2*r3 - a3*b1*r2 + a4*b1*r1)
f(925) = 2*dtau**3*(-a1*b2*r4 + a2*b2*r3 - a3*b1*r2 + a4*b1*r1)
f(926) = 2*dtau**2*(a1*b2**2*r4 + a2*b2**2*r3 + a3*b1**2*r2 + a4*b1**2* &
      r1)
f(927) = 2*dtau*(-a1*b2**3*r4 + a2*b2**3*r3 - a3*b1**3*r2 + a4*b1**3*r1)
f(928) = 2*dtau*(-a1**2*b2*r4 + a2**2*b2*r3 - a3**2*b1*r2 + a4**2*b1*r1)
f(929) = 2*dtau*(-a1**3*b2*r4 + a2**3*b2*r3 - a3**3*b1*r2 + a4**3*b1*r1)
f(930) = 2*dtau*(-a1*b2*r4**2 + a2*b2*r3**2 - a3*b1*r2**2 + a4*b1*r1**2)
f(931) = 2*dtau*(-a1**2*b2*r4**2 + a2**2*b2*r3**2 - a3**2*b1*r2**2 + a4 &
      **2*b1*r1**2)
f(932) = 2*dtau*(-a1*b2*r4**3 + a2*b2*r3**3 - a3*b1*r2**3 + a4*b1*r1**3)
f(933) = 2*dtau*(-a1*b1*r4 + a2*b1*r3 - a3*b2*r2 + a4*b2*r1)
f(934) = 2*dtau**3*(-a1*b1*r4 + a2*b1*r3 - a3*b2*r2 + a4*b2*r1)
f(935) = 2*dtau**2*(a1*b1**2*r4 + a2*b1**2*r3 + a3*b2**2*r2 + a4*b2**2* &
      r1)
f(936) = 2*dtau*(-a1*b1**3*r4 + a2*b1**3*r3 - a3*b2**3*r2 + a4*b2**3*r1)
f(937) = 2*dtau*(-a1**2*b1*r4 + a2**2*b1*r3 - a3**2*b2*r2 + a4**2*b2*r1)
f(938) = 2*dtau*(-a1**3*b1*r4 + a2**3*b1*r3 - a3**3*b2*r2 + a4**3*b2*r1)
f(939) = 2*dtau*(-a1*b1*r4**2 + a2*b1*r3**2 - a3*b2*r2**2 + a4*b2*r1**2)
f(940) = 2*dtau*(-a1**2*b1*r4**2 + a2**2*b1*r3**2 - a3**2*b2*r2**2 + a4 &
      **2*b2*r1**2)
f(941) = 2*dtau*(-a1*b1*r4**3 + a2*b1*r3**3 - a3*b2*r2**3 + a4*b2*r1**3)
f(942) = 2*b1*b2*dtau**2*(r1 + r2 + r3 + r4)
f(943) = 2*b1*b2*dtau*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(944) = 2*b1*b2*dtau*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(945) = 2*b1*b2*dtau**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(946) = 2*b1*b2*dtau*(b1*r3**2 - b1*r4**2 + b2*r1**2 - b2*r2**2)
f(947) = 2*b1*b2*dtau*(b1*r1**2 - b1*r2**2 + b2*r3**2 - b2*r4**2)
f(948) = 8*a1*a2*a3*a4
f(949) = 2*a1*a2*a3*a4*(a1 + a2 + a3 + a4)
f(950) = 2*a1*a2*a3*a4*(a1**2 + a2**2 + a3**2 + a4**2)
f(951) = 4*a1*a2*a3*a4*(a1*a2 + a3*a4)
f(952) = 4*a1*a2*a3*a4*(a1*a3 + a2*a4)
f(953) = 4*a1*a2*a3*a4*(a1*a4 + a2*a3)
f(954) = 2*a1*a2*a3*b1**2 + 2*a1*a2*a4*b1**2 + 2*a1*a3*a4*b2**2 + 2*a2* &
      a3*a4*b2**2
f(955) = 2*a1**2*a3*a4*b2**2 + 2*a1*a2*a3**2*b1**2 + 2*a1*a2*a4**2*b1**2 &
      + 2*a2**2*a3*a4*b2**2
f(956) = 2*a1**2*a2*a4*b1**2 + 2*a1*a2**2*a3*b1**2 + 2*a1*a3*a4**2*b2**2 &
      + 2*a2*a3**2*a4*b2**2
f(957) = 2*a1**2*a2*a3*b1**2 + 2*a1*a2**2*a4*b1**2 + 2*a1*a3**2*a4*b2**2 &
      + 2*a2*a3*a4**2*b2**2
f(958) = 2*a1*a2*a3*b2**2 + 2*a1*a2*a4*b2**2 + 2*a1*a3*a4*b1**2 + 2*a2* &
      a3*a4*b1**2
f(959) = 2*a1**2*a3*a4*b1**2 + 2*a1*a2*a3**2*b2**2 + 2*a1*a2*a4**2*b2**2 &
      + 2*a2**2*a3*a4*b1**2
f(960) = 2*a1**2*a2*a4*b2**2 + 2*a1*a2**2*a3*b2**2 + 2*a1*a3*a4**2*b1**2 &
      + 2*a2*a3**2*a4*b1**2
f(961) = 2*a1**2*a2*a3*b2**2 + 2*a1*a2**2*a4*b2**2 + 2*a1*a3**2*a4*b1**2 &
      + 2*a2*a3*a4**2*b1**2
f(962) = 2*dtau**2*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(963) = 2*dtau**2*(a1**2*a3*a4 + a1*a2*a3**2 + a1*a2*a4**2 + a2**2*a3* &
      a4)
f(964) = 2*dtau**2*(a1**2*a2*a4 + a1*a2**2*a3 + a1*a3*a4**2 + a2*a3**2* &
      a4)
f(965) = 2*dtau**2*(a1**2*a2*a3 + a1*a2**2*a4 + a1*a3**2*a4 + a2*a3*a4** &
      2)
f(966) = 4*b1*b2*(a1*a2 + a3*a4)
f(967) = 4*b1*b2*(a1*a2*b2**2 + a3*a4*b1**2)
f(968) = 4*b1**2*b2**2*(a1*a2 + a3*a4)
f(969) = 4*b1*b2*(a1*a2*b1**2 + a3*a4*b2**2)
f(970) = 2*b1*b2*(a1**2*a2 + a1*a2**2 + a3**2*a4 + a3*a4**2)
f(971) = 2*b1*b2*(a1**3*a2 + a1*a2**3 + a3**3*a4 + a3*a4**3)
f(972) = 4*b1*b2*(a1**2*a2**2 + a3**2*a4**2)
f(973) = 4*dtau**2*(a1*a2*b1**2 + a3*a4*b2**2)
f(974) = 2*dtau*(-a1**2*a2*b1 + a1*a2**2*b1 - a3**2*a4*b2 + a3*a4**2*b2)
f(975) = 2*dtau*(-a1**3*a2*b1 + a1*a2**3*b1 - a3**3*a4*b2 + a3*a4**3*b2)
f(976) = 2*dtau*(a1**2*a2*b1 - a1*a2**2*b1 + a3**2*a4*b2 - a3*a4**2*b2)
f(977) = 2*dtau*(a1**3*a2*b1 - a1*a2**3*b1 + a3**3*a4*b2 - a3*a4**3*b2)
f(978) = 4*dtau**2*(a1*a2*b2**2 + a3*a4*b1**2)
f(979) = 2*dtau*(-a1**2*a2*b2 + a1*a2**2*b2 - a3**2*a4*b1 + a3*a4**2*b1)
f(980) = 2*dtau*(-a1**3*a2*b2 + a1*a2**3*b2 - a3**3*a4*b1 + a3*a4**3*b1)
f(981) = 2*dtau*(a1**2*a2*b2 - a1*a2**2*b2 + a3**2*a4*b1 - a3*a4**2*b1)
f(982) = 2*dtau*(a1**3*a2*b2 - a1*a2**3*b2 + a3**3*a4*b1 - a3*a4**3*b1)
f(983) = 4*b1*b2*(a1*a3 + a2*a4)
f(984) = 2*b1*b2*(a1*a3*b1**2 + a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2**2)
f(985) = 4*b1**2*b2**2*(a1*a3 + a2*a4)
f(986) = 2*b1*b2*(a1**2*a3 + a1*a3**2 + a2**2*a4 + a2*a4**2)
f(987) = 2*b1*b2*(a1**3*a3 + a1*a3**3 + a2**3*a4 + a2*a4**3)
f(988) = 4*b1*b2*(a1**2*a3**2 + a2**2*a4**2)
f(989) = 2*dtau*(a1*a3*b1 + a1*a3*b2 - a2*a4*b1 - a2*a4*b2)
f(990) = 2*dtau**3*(a1*a3*b1 + a1*a3*b2 - a2*a4*b1 - a2*a4*b2)
f(991) = 2*dtau**2*(a1*a3*b1**2 + a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2** &
      2)
f(992) = 2*dtau*(a1*a3*b1**3 + a1*a3*b2**3 - a2*a4*b1**3 - a2*a4*b2**3)
f(993) = 2*dtau*(a1**2*a3*b2 + a1*a3**2*b1 - a2**2*a4*b2 - a2*a4**2*b1)
f(994) = 2*dtau*(a1**3*a3*b2 + a1*a3**3*b1 - a2**3*a4*b2 - a2*a4**3*b1)
f(995) = 2*dtau*(a1**2*a3*b1 + a1*a3**2*b2 - a2**2*a4*b1 - a2*a4**2*b2)
f(996) = 2*dtau*(a1**2*a3**2*b1 + a1**2*a3**2*b2 - a2**2*a4**2*b1 - a2** &
      2*a4**2*b2)
f(997) = 2*dtau*(a1**3*a3*b1 + a1*a3**3*b2 - a2**3*a4*b1 - a2*a4**3*b2)
f(998) = 4*b1*b2*(a1*a4 + a2*a3)
f(999) = 2*b1*b2*(a1*a4*b1**2 + a1*a4*b2**2 + a2*a3*b1**2 + a2*a3*b2**2)
f(1000) = 4*b1**2*b2**2*(a1*a4 + a2*a3)
f(1001) = 2*b1*b2*(a1**2*a4 + a1*a4**2 + a2**2*a3 + a2*a3**2)
f(1002) = 2*b1*b2*(a1**3*a4 + a1*a4**3 + a2**3*a3 + a2*a3**3)
f(1003) = 4*b1*b2*(a1**2*a4**2 + a2**2*a3**2)
f(1004) = 2*dtau*(a1*a4*b1 - a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
f(1005) = 2*dtau**3*(a1*a4*b1 - a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
f(1006) = 2*dtau**2*(a1*a4*b1**2 + a1*a4*b2**2 + a2*a3*b1**2 + a2*a3*b2 &
      **2)
f(1007) = 2*dtau*(a1*a4*b1**3 - a1*a4*b2**3 - a2*a3*b1**3 + a2*a3*b2**3)
f(1008) = 2*dtau*(-a1**2*a4*b2 + a1*a4**2*b1 + a2**2*a3*b2 - a2*a3**2*b1 &
      )
f(1009) = 2*dtau*(-a1**3*a4*b2 + a1*a4**3*b1 + a2**3*a3*b2 - a2*a3**3*b1 &
      )
f(1010) = 2*dtau*(a1**2*a4*b1 - a1*a4**2*b2 - a2**2*a3*b1 + a2*a3**2*b2)
f(1011) = 2*dtau*(a1**2*a4**2*b1 - a1**2*a4**2*b2 - a2**2*a3**2*b1 + a2 &
      **2*a3**2*b2)
f(1012) = 2*dtau*(a1**3*a4*b1 - a1*a4**3*b2 - a2**3*a3*b1 + a2*a3**3*b2)
f(1013) = 2*b1*b2*dtau**2*(a1 + a2 + a3 + a4)
f(1014) = 2*b1*b2*dtau*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(1015) = 2*b1*b2*dtau*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(1016) = 2*b1*b2*dtau**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(1017) = 2*b1*b2*dtau*(a1**2*b2 - a2**2*b2 + a3**2*b1 - a4**2*b1)
f(1018) = 2*b1*b2*dtau*(a1**2*b1 - a2**2*b1 + a3**2*b2 - a4**2*b2)
v = sum(f*params)
end function c2h4_poten_n4_d6_ADF


!###############################################################################


! D2h-symmetry-adapted 5-mode 6-order expansion for C2H4 molecule

function c2h4_poten_n5_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(588)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(588)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 8*r0*r1*r2*r3*r4
f(2) = 2*r0*r1*r2*r3*r4*(r1 + r2 + r3 + r4)
f(3) = 8*r0**2*r1*r2*r3*r4
f(4) = 2*r0*(a1*r1*r2*r3 + a2*r1*r2*r4 + a3*r1*r3*r4 + a4*r2*r3*r4)
f(5) = 2*r0*(a1**2*r1*r2*r3 + a2**2*r1*r2*r4 + a3**2*r1*r3*r4 + a4**2*r2 &
      *r3*r4)
f(6) = 2*r0*(a1*r1*r2*r3**2 + a2*r1*r2*r4**2 + a3*r1**2*r3*r4 + a4*r2**2 &
      *r3*r4)
f(7) = 2*r0*(a1*r1*r2**2*r3 + a2*r1**2*r2*r4 + a3*r1*r3*r4**2 + a4*r2*r3 &
      **2*r4)
f(8) = 2*r0*(a1*r1**2*r2*r3 + a2*r1*r2**2*r4 + a3*r1*r3**2*r4 + a4*r2*r3 &
      *r4**2)
f(9) = 2*r0**2*(a1*r1*r2*r3 + a2*r1*r2*r4 + a3*r1*r3*r4 + a4*r2*r3*r4)
f(10) = 2*r0*(a1*r1*r2*r4 + a2*r1*r2*r3 + a3*r2*r3*r4 + a4*r1*r3*r4)
f(11) = 2*r0*(a1**2*r1*r2*r4 + a2**2*r1*r2*r3 + a3**2*r2*r3*r4 + a4**2* &
      r1*r3*r4)
f(12) = 2*r0*(a1*r1*r2*r4**2 + a2*r1*r2*r3**2 + a3*r2**2*r3*r4 + a4*r1** &
      2*r3*r4)
f(13) = 2*r0*(a1*r1**2*r2*r4 + a2*r1*r2**2*r3 + a3*r2*r3**2*r4 + a4*r1* &
      r3*r4**2)
f(14) = 2*r0*(a1*r1*r2**2*r4 + a2*r1**2*r2*r3 + a3*r2*r3*r4**2 + a4*r1* &
      r3**2*r4)
f(15) = 2*r0**2*(a1*r1*r2*r4 + a2*r1*r2*r3 + a3*r2*r3*r4 + a4*r1*r3*r4)
f(16) = 2*r0*(a1*r1*r3*r4 + a2*r2*r3*r4 + a3*r1*r2*r3 + a4*r1*r2*r4)
f(17) = 2*r0*(a1**2*r1*r3*r4 + a2**2*r2*r3*r4 + a3**2*r1*r2*r3 + a4**2* &
      r1*r2*r4)
f(18) = 2*r0*(a1*r1**2*r3*r4 + a2*r2**2*r3*r4 + a3*r1*r2*r3**2 + a4*r1* &
      r2*r4**2)
f(19) = 2*r0*(a1*r1*r3*r4**2 + a2*r2*r3**2*r4 + a3*r1*r2**2*r3 + a4*r1** &
      2*r2*r4)
f(20) = 2*r0*(a1*r1*r3**2*r4 + a2*r2*r3*r4**2 + a3*r1**2*r2*r3 + a4*r1* &
      r2**2*r4)
f(21) = 2*r0**2*(a1*r1*r3*r4 + a2*r2*r3*r4 + a3*r1*r2*r3 + a4*r1*r2*r4)
f(22) = 2*r0*(a1*r2*r3*r4 + a2*r1*r3*r4 + a3*r1*r2*r4 + a4*r1*r2*r3)
f(23) = 2*r0*(a1**2*r2*r3*r4 + a2**2*r1*r3*r4 + a3**2*r1*r2*r4 + a4**2* &
      r1*r2*r3)
f(24) = 2*r0*(a1*r2**2*r3*r4 + a2*r1**2*r3*r4 + a3*r1*r2*r4**2 + a4*r1* &
      r2*r3**2)
f(25) = 2*r0*(a1*r2*r3**2*r4 + a2*r1*r3*r4**2 + a3*r1**2*r2*r4 + a4*r1* &
      r2**2*r3)
f(26) = 2*r0*(a1*r2*r3*r4**2 + a2*r1*r3**2*r4 + a3*r1*r2**2*r4 + a4*r1** &
      2*r2*r3)
f(27) = 2*r0**2*(a1*r2*r3*r4 + a2*r1*r3*r4 + a3*r1*r2*r4 + a4*r1*r2*r3)
f(28) = 2*r0*(b1**2*r1*r2*r3 + b1**2*r1*r2*r4 + b2**2*r1*r3*r4 + b2**2* &
      r2*r3*r4)
f(29) = 2*r0*(b1**2*r1*r3*r4 + b1**2*r2*r3*r4 + b2**2*r1*r2*r3 + b2**2* &
      r1*r2*r4)
f(30) = 2*dtau**2*r0*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(31) = 4*r0*(a1*a2*r1*r2 + a3*a4*r3*r4)
f(32) = 2*r0*(a1**2*a2*r1*r2 + a1*a2**2*r1*r2 + a3**2*a4*r3*r4 + a3*a4** &
      2*r3*r4)
f(33) = 2*r0*(a1*a2*r1**2*r2 + a1*a2*r1*r2**2 + a3*a4*r3**2*r4 + a3*a4* &
      r3*r4**2)
f(34) = 4*r0**2*(a1*a2*r1*r2 + a3*a4*r3*r4)
f(35) = 2*r0*(a1*a3*r1*r2 + a1*a3*r3*r4 + a2*a4*r1*r2 + a2*a4*r3*r4)
f(36) = 2*r0*(a1**2*a3*r3*r4 + a1*a3**2*r1*r2 + a2**2*a4*r3*r4 + a2*a4** &
      2*r1*r2)
f(37) = 2*r0*(a1**2*a3*r1*r2 + a1*a3**2*r3*r4 + a2**2*a4*r1*r2 + a2*a4** &
      2*r3*r4)
f(38) = 2*r0*(a1*a3*r1*r2**2 + a1*a3*r3*r4**2 + a2*a4*r1**2*r2 + a2*a4* &
      r3**2*r4)
f(39) = 2*r0*(a1*a3*r1**2*r2 + a1*a3*r3**2*r4 + a2*a4*r1*r2**2 + a2*a4* &
      r3*r4**2)
f(40) = 2*r0**2*(a1*a3*r1*r2 + a1*a3*r3*r4 + a2*a4*r1*r2 + a2*a4*r3*r4)
f(41) = 2*r0*(a1*a4*r1*r2 + a1*a4*r3*r4 + a2*a3*r1*r2 + a2*a3*r3*r4)
f(42) = 2*r0*(a1**2*a4*r3*r4 + a1*a4**2*r1*r2 + a2**2*a3*r3*r4 + a2*a3** &
      2*r1*r2)
f(43) = 2*r0*(a1**2*a4*r1*r2 + a1*a4**2*r3*r4 + a2**2*a3*r1*r2 + a2*a3** &
      2*r3*r4)
f(44) = 2*r0*(a1*a4*r1*r2**2 + a1*a4*r3**2*r4 + a2*a3*r1**2*r2 + a2*a3* &
      r3*r4**2)
f(45) = 2*r0*(a1*a4*r1**2*r2 + a1*a4*r3*r4**2 + a2*a3*r1*r2**2 + a2*a3* &
      r3**2*r4)
f(46) = 2*r0**2*(a1*a4*r1*r2 + a1*a4*r3*r4 + a2*a3*r1*r2 + a2*a3*r3*r4)
f(47) = 2*r0*(a1*b1**2*r1*r2 + a2*b1**2*r1*r2 + a3*b2**2*r3*r4 + a4*b2** &
      2*r3*r4)
f(48) = 2*r0*(a1*b2**2*r1*r2 + a2*b2**2*r1*r2 + a3*b1**2*r3*r4 + a4*b1** &
      2*r3*r4)
f(49) = 2*dtau**2*r0*(a1*r1*r2 + a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(50) = 4*r0*(a1*a2*r3*r4 + a3*a4*r1*r2)
f(51) = 2*r0*(a1**2*a2*r3*r4 + a1*a2**2*r3*r4 + a3**2*a4*r1*r2 + a3*a4** &
      2*r1*r2)
f(52) = 2*r0*(a1*a2*r3**2*r4 + a1*a2*r3*r4**2 + a3*a4*r1**2*r2 + a3*a4* &
      r1*r2**2)
f(53) = 4*r0**2*(a1*a2*r3*r4 + a3*a4*r1*r2)
f(54) = 2*r0*(a1*b2**2*r3*r4 + a2*b2**2*r3*r4 + a3*b1**2*r1*r2 + a4*b1** &
      2*r1*r2)
f(55) = 2*r0*(a1*b1**2*r3*r4 + a2*b1**2*r3*r4 + a3*b2**2*r1*r2 + a4*b2** &
      2*r1*r2)
f(56) = 2*dtau**2*r0*(a1*r3*r4 + a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(57) = 4*b1*b2*r0*(r1*r2 + r3*r4)
f(58) = 2*b1*b2*r0*(r1**2*r2 + r1*r2**2 + r3**2*r4 + r3*r4**2)
f(59) = 4*b1*b2*r0**2*(r1*r2 + r3*r4)
f(60) = 2*dtau*r0*(-b1*r1**2*r2 + b1*r1*r2**2 - b2*r3**2*r4 + b2*r3*r4** &
      2)
f(61) = 2*dtau*r0*(b1*r1**2*r2 - b1*r1*r2**2 + b2*r3**2*r4 - b2*r3*r4**2 &
      )
f(62) = 2*dtau*r0*(-b1*r3**2*r4 + b1*r3*r4**2 - b2*r1**2*r2 + b2*r1*r2** &
      2)
f(63) = 2*dtau*r0*(b1*r3**2*r4 - b1*r3*r4**2 + b2*r1**2*r2 - b2*r1*r2**2 &
      )
f(64) = 2*r0*(a1*a2*r1*r3 + a1*a2*r2*r4 + a3*a4*r1*r3 + a3*a4*r2*r4)
f(65) = 2*r0*(a1**2*a2*r2*r4 + a1*a2**2*r1*r3 + a3**2*a4*r2*r4 + a3*a4** &
      2*r1*r3)
f(66) = 2*r0*(a1**2*a2*r1*r3 + a1*a2**2*r2*r4 + a3**2*a4*r1*r3 + a3*a4** &
      2*r2*r4)
f(67) = 2*r0*(a1*a2*r1*r3**2 + a1*a2*r2*r4**2 + a3*a4*r1**2*r3 + a3*a4* &
      r2**2*r4)
f(68) = 2*r0*(a1*a2*r1**2*r3 + a1*a2*r2**2*r4 + a3*a4*r1*r3**2 + a3*a4* &
      r2*r4**2)
f(69) = 2*r0**2*(a1*a2*r1*r3 + a1*a2*r2*r4 + a3*a4*r1*r3 + a3*a4*r2*r4)
f(70) = 4*r0*(a1*a3*r1*r3 + a2*a4*r2*r4)
f(71) = 2*r0*(a1**2*a3*r1*r3 + a1*a3**2*r1*r3 + a2**2*a4*r2*r4 + a2*a4** &
      2*r2*r4)
f(72) = 2*r0*(a1*a3*r1**2*r3 + a1*a3*r1*r3**2 + a2*a4*r2**2*r4 + a2*a4* &
      r2*r4**2)
f(73) = 4*r0**2*(a1*a3*r1*r3 + a2*a4*r2*r4)
f(74) = 2*r0*(a1*a4*r1*r3 + a1*a4*r2*r4 + a2*a3*r1*r3 + a2*a3*r2*r4)
f(75) = 2*r0*(a1**2*a4*r2*r4 + a1*a4**2*r1*r3 + a2**2*a3*r1*r3 + a2*a3** &
      2*r2*r4)
f(76) = 2*r0*(a1**2*a4*r1*r3 + a1*a4**2*r2*r4 + a2**2*a3*r2*r4 + a2*a3** &
      2*r1*r3)
f(77) = 2*r0*(a1*a4*r1*r3**2 + a1*a4*r2**2*r4 + a2*a3*r1**2*r3 + a2*a3* &
      r2*r4**2)
f(78) = 2*r0*(a1*a4*r1**2*r3 + a1*a4*r2*r4**2 + a2*a3*r1*r3**2 + a2*a3* &
      r2**2*r4)
f(79) = 2*r0**2*(a1*a4*r1*r3 + a1*a4*r2*r4 + a2*a3*r1*r3 + a2*a3*r2*r4)
f(80) = 2*r0*(a1*b1**2*r1*r3 + a2*b1**2*r2*r4 + a3*b2**2*r1*r3 + a4*b2** &
      2*r2*r4)
f(81) = 2*r0*(a1*b2**2*r1*r3 + a2*b2**2*r2*r4 + a3*b1**2*r1*r3 + a4*b1** &
      2*r2*r4)
f(82) = 2*dtau**2*r0*(a1*r1*r3 + a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(83) = 4*r0*(a1*a3*r2*r4 + a2*a4*r1*r3)
f(84) = 2*r0*(a1**2*a3*r2*r4 + a1*a3**2*r2*r4 + a2**2*a4*r1*r3 + a2*a4** &
      2*r1*r3)
f(85) = 2*r0*(a1*a3*r2**2*r4 + a1*a3*r2*r4**2 + a2*a4*r1**2*r3 + a2*a4* &
      r1*r3**2)
f(86) = 4*r0**2*(a1*a3*r2*r4 + a2*a4*r1*r3)
f(87) = 2*r0*(a1*b1**2*r2*r4 + a2*b1**2*r1*r3 + a3*b2**2*r2*r4 + a4*b2** &
      2*r1*r3)
f(88) = 2*r0*(a1*b2**2*r2*r4 + a2*b2**2*r1*r3 + a3*b1**2*r2*r4 + a4*b1** &
      2*r1*r3)
f(89) = 2*dtau**2*r0*(a1*r2*r4 + a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(90) = 4*b1*b2*r0*(r1*r3 + r2*r4)
f(91) = 2*b1*b2*r0*(r1**2*r3 + r1*r3**2 + r2**2*r4 + r2*r4**2)
f(92) = 4*b1*b2*r0**2*(r1*r3 + r2*r4)
f(93) = 2*dtau*r0*(b1*r1*r3 - b1*r2*r4 + b2*r1*r3 - b2*r2*r4)
f(94) = 2*dtau*r0*(b1*r1*r3**2 - b1*r2*r4**2 + b2*r1**2*r3 - b2*r2**2*r4 &
      )
f(95) = 2*dtau*r0*(b1*r1**2*r3 - b1*r2**2*r4 + b2*r1*r3**2 - b2*r2*r4**2 &
      )
f(96) = 2*dtau*r0**2*(b1*r1*r3 - b1*r2*r4 + b2*r1*r3 - b2*r2*r4)
f(97) = 2*r0*(a1*a2*r1*r4 + a1*a2*r2*r3 + a3*a4*r1*r4 + a3*a4*r2*r3)
f(98) = 2*r0*(a1**2*a2*r2*r3 + a1*a2**2*r1*r4 + a3**2*a4*r1*r4 + a3*a4** &
      2*r2*r3)
f(99) = 2*r0*(a1**2*a2*r1*r4 + a1*a2**2*r2*r3 + a3**2*a4*r2*r3 + a3*a4** &
      2*r1*r4)
f(100) = 2*r0*(a1*a2*r1*r4**2 + a1*a2*r2*r3**2 + a3*a4*r1**2*r4 + a3*a4* &
      r2**2*r3)
f(101) = 2*r0*(a1*a2*r1**2*r4 + a1*a2*r2**2*r3 + a3*a4*r1*r4**2 + a3*a4* &
      r2*r3**2)
f(102) = 2*r0**2*(a1*a2*r1*r4 + a1*a2*r2*r3 + a3*a4*r1*r4 + a3*a4*r2*r3)
f(103) = 2*r0*(a1*a3*r1*r4 + a1*a3*r2*r3 + a2*a4*r1*r4 + a2*a4*r2*r3)
f(104) = 2*r0*(a1**2*a3*r2*r3 + a1*a3**2*r1*r4 + a2**2*a4*r1*r4 + a2*a4 &
      **2*r2*r3)
f(105) = 2*r0*(a1**2*a3*r1*r4 + a1*a3**2*r2*r3 + a2**2*a4*r2*r3 + a2*a4 &
      **2*r1*r4)
f(106) = 2*r0*(a1*a3*r1*r4**2 + a1*a3*r2**2*r3 + a2*a4*r1**2*r4 + a2*a4* &
      r2*r3**2)
f(107) = 2*r0*(a1*a3*r1**2*r4 + a1*a3*r2*r3**2 + a2*a4*r1*r4**2 + a2*a4* &
      r2**2*r3)
f(108) = 2*r0**2*(a1*a3*r1*r4 + a1*a3*r2*r3 + a2*a4*r1*r4 + a2*a4*r2*r3)
f(109) = 4*r0*(a1*a4*r1*r4 + a2*a3*r2*r3)
f(110) = 2*r0*(a1**2*a4*r1*r4 + a1*a4**2*r1*r4 + a2**2*a3*r2*r3 + a2*a3 &
      **2*r2*r3)
f(111) = 2*r0*(a1*a4*r1**2*r4 + a1*a4*r1*r4**2 + a2*a3*r2**2*r3 + a2*a3* &
      r2*r3**2)
f(112) = 4*r0**2*(a1*a4*r1*r4 + a2*a3*r2*r3)
f(113) = 2*r0*(a1*b1**2*r1*r4 + a2*b1**2*r2*r3 + a3*b2**2*r2*r3 + a4*b2 &
      **2*r1*r4)
f(114) = 2*r0*(a1*b2**2*r1*r4 + a2*b2**2*r2*r3 + a3*b1**2*r2*r3 + a4*b1 &
      **2*r1*r4)
f(115) = 2*dtau**2*r0*(a1*r1*r4 + a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(116) = 4*r0*(a1*a4*r2*r3 + a2*a3*r1*r4)
f(117) = 2*r0*(a1**2*a4*r2*r3 + a1*a4**2*r2*r3 + a2**2*a3*r1*r4 + a2*a3 &
      **2*r1*r4)
f(118) = 2*r0*(a1*a4*r2**2*r3 + a1*a4*r2*r3**2 + a2*a3*r1**2*r4 + a2*a3* &
      r1*r4**2)
f(119) = 4*r0**2*(a1*a4*r2*r3 + a2*a3*r1*r4)
f(120) = 2*r0*(a1*b1**2*r2*r3 + a2*b1**2*r1*r4 + a3*b2**2*r1*r4 + a4*b2 &
      **2*r2*r3)
f(121) = 2*r0*(a1*b2**2*r2*r3 + a2*b2**2*r1*r4 + a3*b1**2*r1*r4 + a4*b1 &
      **2*r2*r3)
f(122) = 2*dtau**2*r0*(a1*r2*r3 + a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(123) = 4*b1*b2*r0*(r1*r4 + r2*r3)
f(124) = 2*b1*b2*r0*(r1**2*r4 + r1*r4**2 + r2**2*r3 + r2*r3**2)
f(125) = 4*b1*b2*r0**2*(r1*r4 + r2*r3)
f(126) = 2*dtau*r0*(b1*r1*r4 - b1*r2*r3 - b2*r1*r4 + b2*r2*r3)
f(127) = 2*dtau*r0*(b1*r1*r4**2 - b1*r2*r3**2 - b2*r1**2*r4 + b2*r2**2* &
      r3)
f(128) = 2*dtau*r0*(b1*r1**2*r4 - b1*r2**2*r3 - b2*r1*r4**2 + b2*r2*r3** &
      2)
f(129) = 2*dtau*r0**2*(b1*r1*r4 - b1*r2*r3 - b2*r1*r4 + b2*r2*r3)
f(130) = 2*r0*(a1*a2*a3*r1 + a1*a2*a4*r2 + a1*a3*a4*r3 + a2*a3*a4*r4)
f(131) = 2*r0*(a1**2*a3*a4*r3 + a1*a2*a3**2*r1 + a1*a2*a4**2*r2 + a2**2* &
      a3*a4*r4)
f(132) = 2*r0*(a1**2*a2*a4*r2 + a1*a2**2*a3*r1 + a1*a3*a4**2*r3 + a2*a3 &
      **2*a4*r4)
f(133) = 2*r0*(a1**2*a2*a3*r1 + a1*a2**2*a4*r2 + a1*a3**2*a4*r3 + a2*a3* &
      a4**2*r4)
f(134) = 2*r0*(a1*a2*a3*r1**2 + a1*a2*a4*r2**2 + a1*a3*a4*r3**2 + a2*a3* &
      a4*r4**2)
f(135) = 2*r0**2*(a1*a2*a3*r1 + a1*a2*a4*r2 + a1*a3*a4*r3 + a2*a3*a4*r4)
f(136) = 2*r0*(a1*a2*a3*r2 + a1*a2*a4*r1 + a1*a3*a4*r4 + a2*a3*a4*r3)
f(137) = 2*r0*(a1**2*a3*a4*r4 + a1*a2*a3**2*r2 + a1*a2*a4**2*r1 + a2**2* &
      a3*a4*r3)
f(138) = 2*r0*(a1**2*a2*a3*r2 + a1*a2**2*a4*r1 + a1*a3**2*a4*r4 + a2*a3* &
      a4**2*r3)
f(139) = 2*r0*(a1**2*a2*a4*r1 + a1*a2**2*a3*r2 + a1*a3*a4**2*r4 + a2*a3 &
      **2*a4*r3)
f(140) = 2*r0*(a1*a2*a3*r2**2 + a1*a2*a4*r1**2 + a1*a3*a4*r4**2 + a2*a3* &
      a4*r3**2)
f(141) = 2*r0**2*(a1*a2*a3*r2 + a1*a2*a4*r1 + a1*a3*a4*r4 + a2*a3*a4*r3)
f(142) = 2*r0*(a1*a2*b1**2*r1 + a1*a2*b1**2*r2 + a3*a4*b2**2*r3 + a3*a4* &
      b2**2*r4)
f(143) = 2*r0*(a1*a2*b2**2*r1 + a1*a2*b2**2*r2 + a3*a4*b1**2*r3 + a3*a4* &
      b1**2*r4)
f(144) = 2*dtau**2*r0*(a1*a2*r1 + a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(145) = 2*r0*(a1*a2*a3*r3 + a1*a2*a4*r4 + a1*a3*a4*r1 + a2*a3*a4*r2)
f(146) = 2*r0*(a1**2*a2*a4*r4 + a1*a2**2*a3*r3 + a1*a3*a4**2*r1 + a2*a3 &
      **2*a4*r2)
f(147) = 2*r0*(a1**2*a2*a3*r3 + a1*a2**2*a4*r4 + a1*a3**2*a4*r1 + a2*a3* &
      a4**2*r2)
f(148) = 2*r0*(a1**2*a3*a4*r1 + a1*a2*a3**2*r3 + a1*a2*a4**2*r4 + a2**2* &
      a3*a4*r2)
f(149) = 2*r0*(a1*a2*a3*r3**2 + a1*a2*a4*r4**2 + a1*a3*a4*r1**2 + a2*a3* &
      a4*r2**2)
f(150) = 2*r0**2*(a1*a2*a3*r3 + a1*a2*a4*r4 + a1*a3*a4*r1 + a2*a3*a4*r2)
f(151) = 2*r0*(a1*a3*b1**2*r1 + a1*a3*b2**2*r3 + a2*a4*b1**2*r2 + a2*a4* &
      b2**2*r4)
f(152) = 2*r0*(a1*a3*b1**2*r3 + a1*a3*b2**2*r1 + a2*a4*b1**2*r4 + a2*a4* &
      b2**2*r2)
f(153) = 2*dtau**2*r0*(a1*a3*r1 + a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(154) = 2*r0*(a1*a4*b1**2*r1 + a1*a4*b2**2*r4 + a2*a3*b1**2*r2 + a2*a3* &
      b2**2*r3)
f(155) = 2*r0*(a1*a4*b1**2*r4 + a1*a4*b2**2*r1 + a2*a3*b1**2*r3 + a2*a3* &
      b2**2*r2)
f(156) = 2*dtau**2*r0*(a1*a4*r1 + a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(157) = 2*b1*b2*r0*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(158) = 2*b1*b2*r0*(a1**2*r1 + a2**2*r2 + a3**2*r3 + a4**2*r4)
f(159) = 2*b1*b2*r0*(a1*r1**2 + a2*r2**2 + a3*r3**2 + a4*r4**2)
f(160) = 2*b1*b2*r0**2*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(161) = 2*dtau*r0*(a1*b1*r1 - a2*b1*r2 + a3*b2*r3 - a4*b2*r4)
f(162) = 2*dtau*r0*(a1**2*b1*r1 - a2**2*b1*r2 + a3**2*b2*r3 - a4**2*b2* &
      r4)
f(163) = 2*dtau*r0*(a1*b1*r1**2 - a2*b1*r2**2 + a3*b2*r3**2 - a4*b2*r4** &
      2)
f(164) = 2*dtau*r0**2*(a1*b1*r1 - a2*b1*r2 + a3*b2*r3 - a4*b2*r4)
f(165) = 2*dtau*r0*(a1*b2*r1 - a2*b2*r2 + a3*b1*r3 - a4*b1*r4)
f(166) = 2*dtau*r0*(a1**2*b2*r1 - a2**2*b2*r2 + a3**2*b1*r3 - a4**2*b1* &
      r4)
f(167) = 2*dtau*r0*(a1*b2*r1**2 - a2*b2*r2**2 + a3*b1*r3**2 - a4*b1*r4** &
      2)
f(168) = 2*dtau*r0**2*(a1*b2*r1 - a2*b2*r2 + a3*b1*r3 - a4*b1*r4)
f(169) = 2*r0*(a1*a2*a3*r4 + a1*a2*a4*r3 + a1*a3*a4*r2 + a2*a3*a4*r1)
f(170) = 2*r0*(a1**2*a2*a3*r4 + a1*a2**2*a4*r3 + a1*a3**2*a4*r2 + a2*a3* &
      a4**2*r1)
f(171) = 2*r0*(a1**2*a2*a4*r3 + a1*a2**2*a3*r4 + a1*a3*a4**2*r2 + a2*a3 &
      **2*a4*r1)
f(172) = 2*r0*(a1**2*a3*a4*r2 + a1*a2*a3**2*r4 + a1*a2*a4**2*r3 + a2**2* &
      a3*a4*r1)
f(173) = 2*r0*(a1*a2*a3*r4**2 + a1*a2*a4*r3**2 + a1*a3*a4*r2**2 + a2*a3* &
      a4*r1**2)
f(174) = 2*r0**2*(a1*a2*a3*r4 + a1*a2*a4*r3 + a1*a3*a4*r2 + a2*a3*a4*r1)
f(175) = 2*r0*(a1*a4*b1**2*r2 + a1*a4*b2**2*r3 + a2*a3*b1**2*r1 + a2*a3* &
      b2**2*r4)
f(176) = 2*r0*(a1*a4*b1**2*r3 + a1*a4*b2**2*r2 + a2*a3*b1**2*r4 + a2*a3* &
      b2**2*r1)
f(177) = 2*dtau**2*r0*(a1*a4*r2 + a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(178) = 2*r0*(a1*a3*b1**2*r2 + a1*a3*b2**2*r4 + a2*a4*b1**2*r1 + a2*a4* &
      b2**2*r3)
f(179) = 2*r0*(a1*a3*b1**2*r4 + a1*a3*b2**2*r2 + a2*a4*b1**2*r3 + a2*a4* &
      b2**2*r1)
f(180) = 2*dtau**2*r0*(a1*a3*r2 + a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(181) = 2*b1*b2*r0*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(182) = 2*b1*b2*r0*(a1**2*r2 + a2**2*r1 + a3**2*r4 + a4**2*r3)
f(183) = 2*b1*b2*r0*(a1*r2**2 + a2*r1**2 + a3*r4**2 + a4*r3**2)
f(184) = 2*b1*b2*r0**2*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(185) = 2*dtau*r0*(-a1*b1*r2 + a2*b1*r1 - a3*b2*r4 + a4*b2*r3)
f(186) = 2*dtau*r0*(-a1**2*b1*r2 + a2**2*b1*r1 - a3**2*b2*r4 + a4**2*b2* &
      r3)
f(187) = 2*dtau*r0*(-a1*b1*r2**2 + a2*b1*r1**2 - a3*b2*r4**2 + a4*b2*r3 &
      **2)
f(188) = 2*dtau*r0**2*(-a1*b1*r2 + a2*b1*r1 - a3*b2*r4 + a4*b2*r3)
f(189) = 2*dtau*r0*(-a1*b2*r2 + a2*b2*r1 - a3*b1*r4 + a4*b1*r3)
f(190) = 2*dtau*r0*(-a1**2*b2*r2 + a2**2*b2*r1 - a3**2*b1*r4 + a4**2*b1* &
      r3)
f(191) = 2*dtau*r0*(-a1*b2*r2**2 + a2*b2*r1**2 - a3*b1*r4**2 + a4*b1*r3 &
      **2)
f(192) = 2*dtau*r0**2*(-a1*b2*r2 + a2*b2*r1 - a3*b1*r4 + a4*b1*r3)
f(193) = 2*r0*(a1*a2*b2**2*r3 + a1*a2*b2**2*r4 + a3*a4*b1**2*r1 + a3*a4* &
      b1**2*r2)
f(194) = 2*r0*(a1*a2*b1**2*r3 + a1*a2*b1**2*r4 + a3*a4*b2**2*r1 + a3*a4* &
      b2**2*r2)
f(195) = 2*dtau**2*r0*(a1*a2*r3 + a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(196) = 2*b1*b2*r0*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(197) = 2*b1*b2*r0*(a1**2*r3 + a2**2*r4 + a3**2*r1 + a4**2*r2)
f(198) = 2*b1*b2*r0*(a1*r3**2 + a2*r4**2 + a3*r1**2 + a4*r2**2)
f(199) = 2*b1*b2*r0**2*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(200) = 2*dtau*r0*(a1*b2*r3 - a2*b2*r4 + a3*b1*r1 - a4*b1*r2)
f(201) = 2*dtau*r0*(a1**2*b2*r3 - a2**2*b2*r4 + a3**2*b1*r1 - a4**2*b1* &
      r2)
f(202) = 2*dtau*r0*(a1*b2*r3**2 - a2*b2*r4**2 + a3*b1*r1**2 - a4*b1*r2** &
      2)
f(203) = 2*dtau*r0**2*(a1*b2*r3 - a2*b2*r4 + a3*b1*r1 - a4*b1*r2)
f(204) = 2*dtau*r0*(a1*b1*r3 - a2*b1*r4 + a3*b2*r1 - a4*b2*r2)
f(205) = 2*dtau*r0*(a1**2*b1*r3 - a2**2*b1*r4 + a3**2*b2*r1 - a4**2*b2* &
      r2)
f(206) = 2*dtau*r0*(a1*b1*r3**2 - a2*b1*r4**2 + a3*b2*r1**2 - a4*b2*r2** &
      2)
f(207) = 2*dtau*r0**2*(a1*b1*r3 - a2*b1*r4 + a3*b2*r1 - a4*b2*r2)
f(208) = 2*b1*b2*r0*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(209) = 2*b1*b2*r0*(a1**2*r4 + a2**2*r3 + a3**2*r2 + a4**2*r1)
f(210) = 2*b1*b2*r0*(a1*r4**2 + a2*r3**2 + a3*r2**2 + a4*r1**2)
f(211) = 2*b1*b2*r0**2*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(212) = 2*dtau*r0*(-a1*b2*r4 + a2*b2*r3 - a3*b1*r2 + a4*b1*r1)
f(213) = 2*dtau*r0*(-a1**2*b2*r4 + a2**2*b2*r3 - a3**2*b1*r2 + a4**2*b1* &
      r1)
f(214) = 2*dtau*r0*(-a1*b2*r4**2 + a2*b2*r3**2 - a3*b1*r2**2 + a4*b1*r1 &
      **2)
f(215) = 2*dtau*r0**2*(-a1*b2*r4 + a2*b2*r3 - a3*b1*r2 + a4*b1*r1)
f(216) = 2*dtau*r0*(-a1*b1*r4 + a2*b1*r3 - a3*b2*r2 + a4*b2*r1)
f(217) = 2*dtau*r0*(-a1**2*b1*r4 + a2**2*b1*r3 - a3**2*b2*r2 + a4**2*b2* &
      r1)
f(218) = 2*dtau*r0*(-a1*b1*r4**2 + a2*b1*r3**2 - a3*b2*r2**2 + a4*b2*r1 &
      **2)
f(219) = 2*dtau*r0**2*(-a1*b1*r4 + a2*b1*r3 - a3*b2*r2 + a4*b2*r1)
f(220) = 2*b1*b2*dtau**2*r0*(r1 + r2 + r3 + r4)
f(221) = 2*b1*b2*dtau*r0*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(222) = 2*b1*b2*dtau*r0*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(223) = 8*a1*a2*a3*a4*r0
f(224) = 2*a1*a2*a3*a4*r0*(a1 + a2 + a3 + a4)
f(225) = 8*a1*a2*a3*a4*r0**2
f(226) = 2*r0*(a1*a2*a3*b1**2 + a1*a2*a4*b1**2 + a1*a3*a4*b2**2 + a2*a3* &
      a4*b2**2)
f(227) = 2*r0*(a1*a2*a3*b2**2 + a1*a2*a4*b2**2 + a1*a3*a4*b1**2 + a2*a3* &
      a4*b1**2)
f(228) = 2*dtau**2*r0*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(229) = 4*b1*b2*r0*(a1*a2 + a3*a4)
f(230) = 2*b1*b2*r0*(a1**2*a2 + a1*a2**2 + a3**2*a4 + a3*a4**2)
f(231) = 4*b1*b2*r0**2*(a1*a2 + a3*a4)
f(232) = 2*dtau*r0*(-a1**2*a2*b1 + a1*a2**2*b1 - a3**2*a4*b2 + a3*a4**2* &
      b2)
f(233) = 2*dtau*r0*(a1**2*a2*b1 - a1*a2**2*b1 + a3**2*a4*b2 - a3*a4**2* &
      b2)
f(234) = 2*dtau*r0*(-a1**2*a2*b2 + a1*a2**2*b2 - a3**2*a4*b1 + a3*a4**2* &
      b1)
f(235) = 2*dtau*r0*(a1**2*a2*b2 - a1*a2**2*b2 + a3**2*a4*b1 - a3*a4**2* &
      b1)
f(236) = 4*b1*b2*r0*(a1*a3 + a2*a4)
f(237) = 2*b1*b2*r0*(a1**2*a3 + a1*a3**2 + a2**2*a4 + a2*a4**2)
f(238) = 4*b1*b2*r0**2*(a1*a3 + a2*a4)
f(239) = 2*dtau*r0*(a1*a3*b1 + a1*a3*b2 - a2*a4*b1 - a2*a4*b2)
f(240) = 2*dtau*r0*(a1**2*a3*b2 + a1*a3**2*b1 - a2**2*a4*b2 - a2*a4**2* &
      b1)
f(241) = 2*dtau*r0*(a1**2*a3*b1 + a1*a3**2*b2 - a2**2*a4*b1 - a2*a4**2* &
      b2)
f(242) = 2*dtau*r0**2*(a1*a3*b1 + a1*a3*b2 - a2*a4*b1 - a2*a4*b2)
f(243) = 4*b1*b2*r0*(a1*a4 + a2*a3)
f(244) = 2*b1*b2*r0*(a1**2*a4 + a1*a4**2 + a2**2*a3 + a2*a3**2)
f(245) = 4*b1*b2*r0**2*(a1*a4 + a2*a3)
f(246) = 2*dtau*r0*(a1*a4*b1 - a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
f(247) = 2*dtau*r0*(-a1**2*a4*b2 + a1*a4**2*b1 + a2**2*a3*b2 - a2*a3**2* &
      b1)
f(248) = 2*dtau*r0*(a1**2*a4*b1 - a1*a4**2*b2 - a2**2*a3*b1 + a2*a3**2* &
      b2)
f(249) = 2*dtau*r0**2*(a1*a4*b1 - a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
f(250) = 2*b1*b2*dtau**2*r0*(a1 + a2 + a3 + a4)
f(251) = 2*b1*b2*dtau*r0*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(252) = 2*b1*b2*dtau*r0*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(253) = 2*r1*r2*r3*r4*(a1 + a2 + a3 + a4)
f(254) = 2*r1*r2*r3*r4*(a1**2 + a2**2 + a3**2 + a4**2)
f(255) = 2*r1*r2*r3*r4*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(256) = 2*r1*r2*r3*r4*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(257) = 2*r1*r2*r3*r4*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(258) = 2*r1*r2*r3*r4*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(259) = 4*r1*r2*r3*r4*(b1**2 + b2**2)
f(260) = 8*dtau**2*r1*r2*r3*r4
f(261) = 2*a1*a2*r1*r2*r3 + 2*a1*a2*r1*r2*r4 + 2*a3*a4*r1*r3*r4 + 2*a3* &
      a4*r2*r3*r4
f(262) = 2*a1**2*a2*r1*r2*r4 + 2*a1*a2**2*r1*r2*r3 + 2*a3**2*a4*r2*r3*r4 &
      + 2*a3*a4**2*r1*r3*r4
f(263) = 2*a1**2*a2*r1*r2*r3 + 2*a1*a2**2*r1*r2*r4 + 2*a3**2*a4*r1*r3*r4 &
      + 2*a3*a4**2*r2*r3*r4
f(264) = 2*a1*a2*r1*r2*r3**2 + 2*a1*a2*r1*r2*r4**2 + 2*a3*a4*r1**2*r3*r4 &
      + 2*a3*a4*r2**2*r3*r4
f(265) = 2*a1*a2*r1**2*r2*r4 + 2*a1*a2*r1*r2**2*r3 + 2*a3*a4*r1*r3*r4**2 &
      + 2*a3*a4*r2*r3**2*r4
f(266) = 2*a1*a2*r1**2*r2*r3 + 2*a1*a2*r1*r2**2*r4 + 2*a3*a4*r1*r3**2*r4 &
      + 2*a3*a4*r2*r3*r4**2
f(267) = 2*a1*a3*r1*r2*r3 + 2*a1*a3*r1*r3*r4 + 2*a2*a4*r1*r2*r4 + 2*a2* &
      a4*r2*r3*r4
f(268) = 2*a1**2*a3*r1*r3*r4 + 2*a1*a3**2*r1*r2*r3 + 2*a2**2*a4*r2*r3*r4 &
      + 2*a2*a4**2*r1*r2*r4
f(269) = 2*a1**2*a3*r1*r2*r3 + 2*a1*a3**2*r1*r3*r4 + 2*a2**2*a4*r1*r2*r4 &
      + 2*a2*a4**2*r2*r3*r4
f(270) = 2*a1*a3*r1**2*r3*r4 + 2*a1*a3*r1*r2*r3**2 + 2*a2*a4*r1*r2*r4**2 &
      + 2*a2*a4*r2**2*r3*r4
f(271) = 2*a1*a3*r1*r2**2*r3 + 2*a1*a3*r1*r3*r4**2 + 2*a2*a4*r1**2*r2*r4 &
      + 2*a2*a4*r2*r3**2*r4
f(272) = 2*a1*a3*r1**2*r2*r3 + 2*a1*a3*r1*r3**2*r4 + 2*a2*a4*r1*r2**2*r4 &
      + 2*a2*a4*r2*r3*r4**2
f(273) = 2*a1*a4*r1*r2*r3 + 2*a1*a4*r2*r3*r4 + 2*a2*a3*r1*r2*r4 + 2*a2* &
      a3*r1*r3*r4
f(274) = 2*a1**2*a4*r2*r3*r4 + 2*a1*a4**2*r1*r2*r3 + 2*a2**2*a3*r1*r3*r4 &
      + 2*a2*a3**2*r1*r2*r4
f(275) = 2*a1**2*a4*r1*r2*r3 + 2*a1*a4**2*r2*r3*r4 + 2*a2**2*a3*r1*r2*r4 &
      + 2*a2*a3**2*r1*r3*r4
f(276) = 2*a1*a4*r1*r2*r3**2 + 2*a1*a4*r2**2*r3*r4 + 2*a2*a3*r1**2*r3*r4 &
      + 2*a2*a3*r1*r2*r4**2
f(277) = 2*a1*a4*r1*r2**2*r3 + 2*a1*a4*r2*r3**2*r4 + 2*a2*a3*r1**2*r2*r4 &
      + 2*a2*a3*r1*r3*r4**2
f(278) = 2*a1*a4*r1**2*r2*r3 + 2*a1*a4*r2*r3*r4**2 + 2*a2*a3*r1*r2**2*r4 &
      + 2*a2*a3*r1*r3**2*r4
f(279) = 2*a1*b1**2*r1*r2*r3 + 2*a2*b1**2*r1*r2*r4 + 2*a3*b2**2*r1*r3*r4 &
      + 2*a4*b2**2*r2*r3*r4
f(280) = 2*a1*b2**2*r1*r2*r3 + 2*a2*b2**2*r1*r2*r4 + 2*a3*b1**2*r1*r3*r4 &
      + 2*a4*b1**2*r2*r3*r4
f(281) = 2*dtau**2*(a1*r1*r2*r3 + a2*r1*r2*r4 + a3*r1*r3*r4 + a4*r2*r3* &
      r4)
f(282) = 2*a1*a4*r1*r2*r4 + 2*a1*a4*r1*r3*r4 + 2*a2*a3*r1*r2*r3 + 2*a2* &
      a3*r2*r3*r4
f(283) = 2*a1**2*a4*r1*r3*r4 + 2*a1*a4**2*r1*r2*r4 + 2*a2**2*a3*r2*r3*r4 &
      + 2*a2*a3**2*r1*r2*r3
f(284) = 2*a1**2*a4*r1*r2*r4 + 2*a1*a4**2*r1*r3*r4 + 2*a2**2*a3*r1*r2*r3 &
      + 2*a2*a3**2*r2*r3*r4
f(285) = 2*a1*a4*r1**2*r3*r4 + 2*a1*a4*r1*r2*r4**2 + 2*a2*a3*r1*r2*r3**2 &
      + 2*a2*a3*r2**2*r3*r4
f(286) = 2*a1*a4*r1**2*r2*r4 + 2*a1*a4*r1*r3*r4**2 + 2*a2*a3*r1*r2**2*r3 &
      + 2*a2*a3*r2*r3**2*r4
f(287) = 2*a1*a4*r1*r2**2*r4 + 2*a1*a4*r1*r3**2*r4 + 2*a2*a3*r1**2*r2*r3 &
      + 2*a2*a3*r2*r3*r4**2
f(288) = 2*a1*a3*r1*r2*r4 + 2*a1*a3*r2*r3*r4 + 2*a2*a4*r1*r2*r3 + 2*a2* &
      a4*r1*r3*r4
f(289) = 2*a1**2*a3*r2*r3*r4 + 2*a1*a3**2*r1*r2*r4 + 2*a2**2*a4*r1*r3*r4 &
      + 2*a2*a4**2*r1*r2*r3
f(290) = 2*a1**2*a3*r1*r2*r4 + 2*a1*a3**2*r2*r3*r4 + 2*a2**2*a4*r1*r2*r3 &
      + 2*a2*a4**2*r1*r3*r4
f(291) = 2*a1*a3*r1*r2*r4**2 + 2*a1*a3*r2**2*r3*r4 + 2*a2*a4*r1**2*r3*r4 &
      + 2*a2*a4*r1*r2*r3**2
f(292) = 2*a1*a3*r1**2*r2*r4 + 2*a1*a3*r2*r3**2*r4 + 2*a2*a4*r1*r2**2*r3 &
      + 2*a2*a4*r1*r3*r4**2
f(293) = 2*a1*a3*r1*r2**2*r4 + 2*a1*a3*r2*r3*r4**2 + 2*a2*a4*r1**2*r2*r3 &
      + 2*a2*a4*r1*r3**2*r4
f(294) = 2*a1*b1**2*r1*r2*r4 + 2*a2*b1**2*r1*r2*r3 + 2*a3*b2**2*r2*r3*r4 &
      + 2*a4*b2**2*r1*r3*r4
f(295) = 2*a1*b2**2*r1*r2*r4 + 2*a2*b2**2*r1*r2*r3 + 2*a3*b1**2*r2*r3*r4 &
      + 2*a4*b1**2*r1*r3*r4
f(296) = 2*dtau**2*(a1*r1*r2*r4 + a2*r1*r2*r3 + a3*r2*r3*r4 + a4*r1*r3* &
      r4)
f(297) = 2*a1*a2*r1*r3*r4 + 2*a1*a2*r2*r3*r4 + 2*a3*a4*r1*r2*r3 + 2*a3* &
      a4*r1*r2*r4
f(298) = 2*a1**2*a2*r2*r3*r4 + 2*a1*a2**2*r1*r3*r4 + 2*a3**2*a4*r1*r2*r4 &
      + 2*a3*a4**2*r1*r2*r3
f(299) = 2*a1**2*a2*r1*r3*r4 + 2*a1*a2**2*r2*r3*r4 + 2*a3**2*a4*r1*r2*r3 &
      + 2*a3*a4**2*r1*r2*r4
f(300) = 2*a1*a2*r1**2*r3*r4 + 2*a1*a2*r2**2*r3*r4 + 2*a3*a4*r1*r2*r3**2 &
      + 2*a3*a4*r1*r2*r4**2
f(301) = 2*a1*a2*r1*r3*r4**2 + 2*a1*a2*r2*r3**2*r4 + 2*a3*a4*r1**2*r2*r4 &
      + 2*a3*a4*r1*r2**2*r3
f(302) = 2*a1*a2*r1*r3**2*r4 + 2*a1*a2*r2*r3*r4**2 + 2*a3*a4*r1**2*r2*r3 &
      + 2*a3*a4*r1*r2**2*r4
f(303) = 2*a1*b2**2*r1*r3*r4 + 2*a2*b2**2*r2*r3*r4 + 2*a3*b1**2*r1*r2*r3 &
      + 2*a4*b1**2*r1*r2*r4
f(304) = 2*a1*b1**2*r1*r3*r4 + 2*a2*b1**2*r2*r3*r4 + 2*a3*b2**2*r1*r2*r3 &
      + 2*a4*b2**2*r1*r2*r4
f(305) = 2*dtau**2*(a1*r1*r3*r4 + a2*r2*r3*r4 + a3*r1*r2*r3 + a4*r1*r2* &
      r4)
f(306) = 2*a1*b2**2*r2*r3*r4 + 2*a2*b2**2*r1*r3*r4 + 2*a3*b1**2*r1*r2*r4 &
      + 2*a4*b1**2*r1*r2*r3
f(307) = 2*a1*b1**2*r2*r3*r4 + 2*a2*b1**2*r1*r3*r4 + 2*a3*b2**2*r1*r2*r4 &
      + 2*a4*b2**2*r1*r2*r3
f(308) = 2*dtau**2*(a1*r2*r3*r4 + a2*r1*r3*r4 + a3*r1*r2*r4 + a4*r1*r2* &
      r3)
f(309) = 2*b1*b2*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(310) = 2*b1*b2*(r1**2*r3*r4 + r1*r2*r3**2 + r1*r2*r4**2 + r2**2*r3*r4)
f(311) = 2*b1*b2*(r1**2*r2*r4 + r1*r2**2*r3 + r1*r3*r4**2 + r2*r3**2*r4)
f(312) = 2*b1*b2*(r1**2*r2*r3 + r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4**2)
f(313) = 2*dtau*(b1*r1*r2*r3 - b1*r1*r2*r4 + b2*r1*r3*r4 - b2*r2*r3*r4)
f(314) = 2*dtau*(b1*r1*r2*r3**2 - b1*r1*r2*r4**2 + b2*r1**2*r3*r4 - b2* &
      r2**2*r3*r4)
f(315) = 2*dtau*(-b1*r1**2*r2*r4 + b1*r1*r2**2*r3 + b2*r1*r3*r4**2 - b2* &
      r2*r3**2*r4)
f(316) = 2*dtau*(b1*r1**2*r2*r3 - b1*r1*r2**2*r4 + b2*r1*r3**2*r4 - b2* &
      r2*r3*r4**2)
f(317) = 2*dtau*(b1*r1*r3*r4 - b1*r2*r3*r4 + b2*r1*r2*r3 - b2*r1*r2*r4)
f(318) = 2*dtau*(b1*r1**2*r3*r4 - b1*r2**2*r3*r4 + b2*r1*r2*r3**2 - b2* &
      r1*r2*r4**2)
f(319) = 2*dtau*(b1*r1*r3*r4**2 - b1*r2*r3**2*r4 - b2*r1**2*r2*r4 + b2* &
      r1*r2**2*r3)
f(320) = 2*dtau*(b1*r1*r3**2*r4 - b1*r2*r3*r4**2 + b2*r1**2*r2*r3 - b2* &
      r1*r2**2*r4)
f(321) = 2*a1*a2*a3*r1*r2 + 2*a1*a2*a4*r1*r2 + 2*a1*a3*a4*r3*r4 + 2*a2* &
      a3*a4*r3*r4
f(322) = 2*a1**2*a3*a4*r3*r4 + 2*a1*a2*a3**2*r1*r2 + 2*a1*a2*a4**2*r1*r2 &
      + 2*a2**2*a3*a4*r3*r4
f(323) = 2*a1**2*a2*a4*r1*r2 + 2*a1*a2**2*a3*r1*r2 + 2*a1*a3*a4**2*r3*r4 &
      + 2*a2*a3**2*a4*r3*r4
f(324) = 2*a1**2*a2*a3*r1*r2 + 2*a1*a2**2*a4*r1*r2 + 2*a1*a3**2*a4*r3*r4 &
      + 2*a2*a3*a4**2*r3*r4
f(325) = 2*a1*a2*a3*r1*r2**2 + 2*a1*a2*a4*r1**2*r2 + 2*a1*a3*a4*r3*r4**2 &
      + 2*a2*a3*a4*r3**2*r4
f(326) = 2*a1*a2*a3*r1**2*r2 + 2*a1*a2*a4*r1*r2**2 + 2*a1*a3*a4*r3**2*r4 &
      + 2*a2*a3*a4*r3*r4**2
f(327) = 4*a1*a2*b1**2*r1*r2 + 4*a3*a4*b2**2*r3*r4
f(328) = 4*a1*a2*b2**2*r1*r2 + 4*a3*a4*b1**2*r3*r4
f(329) = 4*dtau**2*(a1*a2*r1*r2 + a3*a4*r3*r4)
f(330) = 2*a1*a2*a3*r3*r4 + 2*a1*a2*a4*r3*r4 + 2*a1*a3*a4*r1*r2 + 2*a2* &
      a3*a4*r1*r2
f(331) = 2*a1**2*a2*a4*r3*r4 + 2*a1*a2**2*a3*r3*r4 + 2*a1*a3*a4**2*r1*r2 &
      + 2*a2*a3**2*a4*r1*r2
f(332) = 2*a1**2*a2*a3*r3*r4 + 2*a1*a2**2*a4*r3*r4 + 2*a1*a3**2*a4*r1*r2 &
      + 2*a2*a3*a4**2*r1*r2
f(333) = 2*a1**2*a3*a4*r1*r2 + 2*a1*a2*a3**2*r3*r4 + 2*a1*a2*a4**2*r3*r4 &
      + 2*a2**2*a3*a4*r1*r2
f(334) = 2*a1*a2*a3*r3*r4**2 + 2*a1*a2*a4*r3**2*r4 + 2*a1*a3*a4*r1*r2**2 &
      + 2*a2*a3*a4*r1**2*r2
f(335) = 2*a1*a2*a3*r3**2*r4 + 2*a1*a2*a4*r3*r4**2 + 2*a1*a3*a4*r1**2*r2 &
      + 2*a2*a3*a4*r1*r2**2
f(336) = 2*a1*a3*b1**2*r1*r2 + 2*a1*a3*b2**2*r3*r4 + 2*a2*a4*b1**2*r1*r2 &
      + 2*a2*a4*b2**2*r3*r4
f(337) = 2*a1*a3*b1**2*r3*r4 + 2*a1*a3*b2**2*r1*r2 + 2*a2*a4*b1**2*r3*r4 &
      + 2*a2*a4*b2**2*r1*r2
f(338) = 2*dtau**2*(a1*a3*r1*r2 + a1*a3*r3*r4 + a2*a4*r1*r2 + a2*a4*r3* &
      r4)
f(339) = 2*a1*a4*b1**2*r1*r2 + 2*a1*a4*b2**2*r3*r4 + 2*a2*a3*b1**2*r1*r2 &
      + 2*a2*a3*b2**2*r3*r4
f(340) = 2*a1*a4*b1**2*r3*r4 + 2*a1*a4*b2**2*r1*r2 + 2*a2*a3*b1**2*r3*r4 &
      + 2*a2*a3*b2**2*r1*r2
f(341) = 2*dtau**2*(a1*a4*r1*r2 + a1*a4*r3*r4 + a2*a3*r1*r2 + a2*a3*r3* &
      r4)
f(342) = 2*b1*b2*(a1*r1*r2 + a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(343) = 2*b1*b2*(a1**2*r1*r2 + a2**2*r1*r2 + a3**2*r3*r4 + a4**2*r3*r4)
f(344) = 2*b1*b2*(a1*r1*r2**2 + a2*r1**2*r2 + a3*r3*r4**2 + a4*r3**2*r4)
f(345) = 2*b1*b2*(a1*r1**2*r2 + a2*r1*r2**2 + a3*r3**2*r4 + a4*r3*r4**2)
f(346) = 2*dtau*(a1*b1*r1*r2 - a2*b1*r1*r2 + a3*b2*r3*r4 - a4*b2*r3*r4)
f(347) = 2*dtau*(a1**2*b1*r1*r2 - a2**2*b1*r1*r2 + a3**2*b2*r3*r4 - a4** &
      2*b2*r3*r4)
f(348) = 2*dtau*(a1*b1*r1*r2**2 - a2*b1*r1**2*r2 + a3*b2*r3*r4**2 - a4* &
      b2*r3**2*r4)
f(349) = 2*dtau*(a1*b1*r1**2*r2 - a2*b1*r1*r2**2 + a3*b2*r3**2*r4 - a4* &
      b2*r3*r4**2)
f(350) = 2*dtau*(a1*b2*r1*r2 - a2*b2*r1*r2 + a3*b1*r3*r4 - a4*b1*r3*r4)
f(351) = 2*dtau*(a1**2*b2*r1*r2 - a2**2*b2*r1*r2 + a3**2*b1*r3*r4 - a4** &
      2*b1*r3*r4)
f(352) = 2*dtau*(a1*b2*r1*r2**2 - a2*b2*r1**2*r2 + a3*b1*r3*r4**2 - a4* &
      b1*r3**2*r4)
f(353) = 2*dtau*(a1*b2*r1**2*r2 - a2*b2*r1*r2**2 + a3*b1*r3**2*r4 - a4* &
      b1*r3*r4**2)
f(354) = 4*a1*a2*b2**2*r3*r4 + 4*a3*a4*b1**2*r1*r2
f(355) = 4*a1*a2*b1**2*r3*r4 + 4*a3*a4*b2**2*r1*r2
f(356) = 4*dtau**2*(a1*a2*r3*r4 + a3*a4*r1*r2)
f(357) = 2*b1*b2*(a1*r3*r4 + a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(358) = 2*b1*b2*(a1**2*r3*r4 + a2**2*r3*r4 + a3**2*r1*r2 + a4**2*r1*r2)
f(359) = 2*b1*b2*(a1*r3*r4**2 + a2*r3**2*r4 + a3*r1*r2**2 + a4*r1**2*r2)
f(360) = 2*b1*b2*(a1*r3**2*r4 + a2*r3*r4**2 + a3*r1**2*r2 + a4*r1*r2**2)
f(361) = 2*dtau*(a1*b2*r3*r4 - a2*b2*r3*r4 + a3*b1*r1*r2 - a4*b1*r1*r2)
f(362) = 2*dtau*(a1**2*b2*r3*r4 - a2**2*b2*r3*r4 + a3**2*b1*r1*r2 - a4** &
      2*b1*r1*r2)
f(363) = 2*dtau*(a1*b2*r3*r4**2 - a2*b2*r3**2*r4 + a3*b1*r1*r2**2 - a4* &
      b1*r1**2*r2)
f(364) = 2*dtau*(a1*b2*r3**2*r4 - a2*b2*r3*r4**2 + a3*b1*r1**2*r2 - a4* &
      b1*r1*r2**2)
f(365) = 2*dtau*(a1*b1*r3*r4 - a2*b1*r3*r4 + a3*b2*r1*r2 - a4*b2*r1*r2)
f(366) = 2*dtau*(a1**2*b1*r3*r4 - a2**2*b1*r3*r4 + a3**2*b2*r1*r2 - a4** &
      2*b2*r1*r2)
f(367) = 2*dtau*(a1*b1*r3*r4**2 - a2*b1*r3**2*r4 + a3*b2*r1*r2**2 - a4* &
      b2*r1**2*r2)
f(368) = 2*dtau*(a1*b1*r3**2*r4 - a2*b1*r3*r4**2 + a3*b2*r1**2*r2 - a4* &
      b2*r1*r2**2)
f(369) = 4*b1*b2*dtau**2*(r1*r2 + r3*r4)
f(370) = 2*a1*a2*a3*r1*r3 + 2*a1*a2*a4*r2*r4 + 2*a1*a3*a4*r1*r3 + 2*a2* &
      a3*a4*r2*r4
f(371) = 2*a1**2*a3*a4*r1*r3 + 2*a1*a2*a3**2*r1*r3 + 2*a1*a2*a4**2*r2*r4 &
      + 2*a2**2*a3*a4*r2*r4
f(372) = 2*a1**2*a2*a4*r2*r4 + 2*a1*a2**2*a3*r1*r3 + 2*a1*a3*a4**2*r1*r3 &
      + 2*a2*a3**2*a4*r2*r4
f(373) = 2*a1**2*a2*a3*r1*r3 + 2*a1*a2**2*a4*r2*r4 + 2*a1*a3**2*a4*r1*r3 &
      + 2*a2*a3*a4**2*r2*r4
f(374) = 2*a1*a2*a3*r1*r3**2 + 2*a1*a2*a4*r2*r4**2 + 2*a1*a3*a4*r1**2*r3 &
      + 2*a2*a3*a4*r2**2*r4
f(375) = 2*a1*a2*a3*r1**2*r3 + 2*a1*a2*a4*r2**2*r4 + 2*a1*a3*a4*r1*r3**2 &
      + 2*a2*a3*a4*r2*r4**2
f(376) = 2*a1*a2*a3*r2*r4 + 2*a1*a2*a4*r1*r3 + 2*a1*a3*a4*r2*r4 + 2*a2* &
      a3*a4*r1*r3
f(377) = 2*a1**2*a3*a4*r2*r4 + 2*a1*a2*a3**2*r2*r4 + 2*a1*a2*a4**2*r1*r3 &
      + 2*a2**2*a3*a4*r1*r3
f(378) = 2*a1**2*a2*a3*r2*r4 + 2*a1*a2**2*a4*r1*r3 + 2*a1*a3**2*a4*r2*r4 &
      + 2*a2*a3*a4**2*r1*r3
f(379) = 2*a1**2*a2*a4*r1*r3 + 2*a1*a2**2*a3*r2*r4 + 2*a1*a3*a4**2*r2*r4 &
      + 2*a2*a3**2*a4*r1*r3
f(380) = 2*a1*a2*a3*r2*r4**2 + 2*a1*a2*a4*r1*r3**2 + 2*a1*a3*a4*r2**2*r4 &
      + 2*a2*a3*a4*r1**2*r3
f(381) = 2*a1*a2*a3*r2**2*r4 + 2*a1*a2*a4*r1**2*r3 + 2*a1*a3*a4*r2*r4**2 &
      + 2*a2*a3*a4*r1*r3**2
f(382) = 2*a1*a2*b1**2*r1*r3 + 2*a1*a2*b1**2*r2*r4 + 2*a3*a4*b2**2*r1*r3 &
      + 2*a3*a4*b2**2*r2*r4
f(383) = 2*a1*a2*b2**2*r1*r3 + 2*a1*a2*b2**2*r2*r4 + 2*a3*a4*b1**2*r1*r3 &
      + 2*a3*a4*b1**2*r2*r4
f(384) = 2*dtau**2*(a1*a2*r1*r3 + a1*a2*r2*r4 + a3*a4*r1*r3 + a3*a4*r2* &
      r4)
f(385) = 2*a1*a3*b1**2*r1*r3 + 2*a1*a3*b2**2*r1*r3 + 2*a2*a4*b1**2*r2*r4 &
      + 2*a2*a4*b2**2*r2*r4
f(386) = 4*dtau**2*(a1*a3*r1*r3 + a2*a4*r2*r4)
f(387) = 2*a1*a4*b1**2*r1*r3 + 2*a1*a4*b2**2*r2*r4 + 2*a2*a3*b1**2*r2*r4 &
      + 2*a2*a3*b2**2*r1*r3
f(388) = 2*a1*a4*b1**2*r2*r4 + 2*a1*a4*b2**2*r1*r3 + 2*a2*a3*b1**2*r1*r3 &
      + 2*a2*a3*b2**2*r2*r4
f(389) = 2*dtau**2*(a1*a4*r1*r3 + a1*a4*r2*r4 + a2*a3*r1*r3 + a2*a3*r2* &
      r4)
f(390) = 2*b1*b2*(a1*r1*r3 + a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(391) = 2*b1*b2*(a1**2*r1*r3 + a2**2*r2*r4 + a3**2*r1*r3 + a4**2*r2*r4)
f(392) = 2*b1*b2*(a1*r1*r3**2 + a2*r2*r4**2 + a3*r1**2*r3 + a4*r2**2*r4)
f(393) = 2*b1*b2*(a1*r1**2*r3 + a2*r2**2*r4 + a3*r1*r3**2 + a4*r2*r4**2)
f(394) = 2*dtau*(a1*b1*r1*r3 - a2*b1*r2*r4 + a3*b2*r1*r3 - a4*b2*r2*r4)
f(395) = 2*dtau*(a1**2*b1*r1*r3 - a2**2*b1*r2*r4 + a3**2*b2*r1*r3 - a4** &
      2*b2*r2*r4)
f(396) = 2*dtau*(a1*b1*r1*r3**2 - a2*b1*r2*r4**2 + a3*b2*r1**2*r3 - a4* &
      b2*r2**2*r4)
f(397) = 2*dtau*(a1*b1*r1**2*r3 - a2*b1*r2**2*r4 + a3*b2*r1*r3**2 - a4* &
      b2*r2*r4**2)
f(398) = 2*dtau*(a1*b2*r1*r3 - a2*b2*r2*r4 + a3*b1*r1*r3 - a4*b1*r2*r4)
f(399) = 2*dtau*(a1**2*b2*r1*r3 - a2**2*b2*r2*r4 + a3**2*b1*r1*r3 - a4** &
      2*b1*r2*r4)
f(400) = 2*dtau*(a1*b2*r1*r3**2 - a2*b2*r2*r4**2 + a3*b1*r1**2*r3 - a4* &
      b1*r2**2*r4)
f(401) = 2*dtau*(a1*b2*r1**2*r3 - a2*b2*r2**2*r4 + a3*b1*r1*r3**2 - a4* &
      b1*r2*r4**2)
f(402) = 2*a1*a3*b1**2*r2*r4 + 2*a1*a3*b2**2*r2*r4 + 2*a2*a4*b1**2*r1*r3 &
      + 2*a2*a4*b2**2*r1*r3
f(403) = 4*dtau**2*(a1*a3*r2*r4 + a2*a4*r1*r3)
f(404) = 2*b1*b2*(a1*r2*r4 + a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(405) = 2*b1*b2*(a1**2*r2*r4 + a2**2*r1*r3 + a3**2*r2*r4 + a4**2*r1*r3)
f(406) = 2*b1*b2*(a1*r2*r4**2 + a2*r1*r3**2 + a3*r2**2*r4 + a4*r1**2*r3)
f(407) = 2*b1*b2*(a1*r2**2*r4 + a2*r1**2*r3 + a3*r2*r4**2 + a4*r1*r3**2)
f(408) = 2*dtau*(-a1*b1*r2*r4 + a2*b1*r1*r3 - a3*b2*r2*r4 + a4*b2*r1*r3)
f(409) = 2*dtau*(-a1**2*b1*r2*r4 + a2**2*b1*r1*r3 - a3**2*b2*r2*r4 + a4 &
      **2*b2*r1*r3)
f(410) = 2*dtau*(-a1*b1*r2*r4**2 + a2*b1*r1*r3**2 - a3*b2*r2**2*r4 + a4* &
      b2*r1**2*r3)
f(411) = 2*dtau*(-a1*b1*r2**2*r4 + a2*b1*r1**2*r3 - a3*b2*r2*r4**2 + a4* &
      b2*r1*r3**2)
f(412) = 2*dtau*(-a1*b2*r2*r4 + a2*b2*r1*r3 - a3*b1*r2*r4 + a4*b1*r1*r3)
f(413) = 2*dtau*(-a1**2*b2*r2*r4 + a2**2*b2*r1*r3 - a3**2*b1*r2*r4 + a4 &
      **2*b1*r1*r3)
f(414) = 2*dtau*(-a1*b2*r2*r4**2 + a2*b2*r1*r3**2 - a3*b1*r2**2*r4 + a4* &
      b1*r1**2*r3)
f(415) = 2*dtau*(-a1*b2*r2**2*r4 + a2*b2*r1**2*r3 - a3*b1*r2*r4**2 + a4* &
      b1*r1*r3**2)
f(416) = 4*b1*b2*dtau**2*(r1*r3 + r2*r4)
f(417) = 2*b1*b2*dtau*(b1*r1*r3 - b1*r2*r4 + b2*r1*r3 - b2*r2*r4)
f(418) = 2*a1*a2*a3*r1*r4 + 2*a1*a2*a4*r2*r3 + 2*a1*a3*a4*r2*r3 + 2*a2* &
      a3*a4*r1*r4
f(419) = 2*a1**2*a3*a4*r2*r3 + 2*a1*a2*a3**2*r1*r4 + 2*a1*a2*a4**2*r2*r3 &
      + 2*a2**2*a3*a4*r1*r4
f(420) = 2*a1**2*a2*a4*r2*r3 + 2*a1*a2**2*a3*r1*r4 + 2*a1*a3*a4**2*r2*r3 &
      + 2*a2*a3**2*a4*r1*r4
f(421) = 2*a1**2*a2*a3*r1*r4 + 2*a1*a2**2*a4*r2*r3 + 2*a1*a3**2*a4*r2*r3 &
      + 2*a2*a3*a4**2*r1*r4
f(422) = 2*a1*a2*a3*r1*r4**2 + 2*a1*a2*a4*r2*r3**2 + 2*a1*a3*a4*r2**2*r3 &
      + 2*a2*a3*a4*r1**2*r4
f(423) = 2*a1*a2*a3*r1**2*r4 + 2*a1*a2*a4*r2**2*r3 + 2*a1*a3*a4*r2*r3**2 &
      + 2*a2*a3*a4*r1*r4**2
f(424) = 2*a1*a2*a3*r2*r3 + 2*a1*a2*a4*r1*r4 + 2*a1*a3*a4*r1*r4 + 2*a2* &
      a3*a4*r2*r3
f(425) = 2*a1**2*a3*a4*r1*r4 + 2*a1*a2*a3**2*r2*r3 + 2*a1*a2*a4**2*r1*r4 &
      + 2*a2**2*a3*a4*r2*r3
f(426) = 2*a1**2*a2*a3*r2*r3 + 2*a1*a2**2*a4*r1*r4 + 2*a1*a3**2*a4*r1*r4 &
      + 2*a2*a3*a4**2*r2*r3
f(427) = 2*a1**2*a2*a4*r1*r4 + 2*a1*a2**2*a3*r2*r3 + 2*a1*a3*a4**2*r1*r4 &
      + 2*a2*a3**2*a4*r2*r3
f(428) = 2*a1*a2*a3*r2*r3**2 + 2*a1*a2*a4*r1*r4**2 + 2*a1*a3*a4*r1**2*r4 &
      + 2*a2*a3*a4*r2**2*r3
f(429) = 2*a1*a2*a3*r2**2*r3 + 2*a1*a2*a4*r1**2*r4 + 2*a1*a3*a4*r1*r4**2 &
      + 2*a2*a3*a4*r2*r3**2
f(430) = 2*a1*a2*b1**2*r1*r4 + 2*a1*a2*b1**2*r2*r3 + 2*a3*a4*b2**2*r1*r4 &
      + 2*a3*a4*b2**2*r2*r3
f(431) = 2*a1*a2*b2**2*r1*r4 + 2*a1*a2*b2**2*r2*r3 + 2*a3*a4*b1**2*r1*r4 &
      + 2*a3*a4*b1**2*r2*r3
f(432) = 2*dtau**2*(a1*a2*r1*r4 + a1*a2*r2*r3 + a3*a4*r1*r4 + a3*a4*r2* &
      r3)
f(433) = 2*a1*a3*b1**2*r1*r4 + 2*a1*a3*b2**2*r2*r3 + 2*a2*a4*b1**2*r2*r3 &
      + 2*a2*a4*b2**2*r1*r4
f(434) = 2*a1*a3*b1**2*r2*r3 + 2*a1*a3*b2**2*r1*r4 + 2*a2*a4*b1**2*r1*r4 &
      + 2*a2*a4*b2**2*r2*r3
f(435) = 2*dtau**2*(a1*a3*r1*r4 + a1*a3*r2*r3 + a2*a4*r1*r4 + a2*a4*r2* &
      r3)
f(436) = 2*a1*a4*b1**2*r1*r4 + 2*a1*a4*b2**2*r1*r4 + 2*a2*a3*b1**2*r2*r3 &
      + 2*a2*a3*b2**2*r2*r3
f(437) = 4*dtau**2*(a1*a4*r1*r4 + a2*a3*r2*r3)
f(438) = 2*b1*b2*(a1*r1*r4 + a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(439) = 2*b1*b2*(a1**2*r1*r4 + a2**2*r2*r3 + a3**2*r2*r3 + a4**2*r1*r4)
f(440) = 2*b1*b2*(a1*r1*r4**2 + a2*r2*r3**2 + a3*r2**2*r3 + a4*r1**2*r4)
f(441) = 2*b1*b2*(a1*r1**2*r4 + a2*r2**2*r3 + a3*r2*r3**2 + a4*r1*r4**2)
f(442) = 2*dtau*(a1*b1*r1*r4 - a2*b1*r2*r3 + a3*b2*r2*r3 - a4*b2*r1*r4)
f(443) = 2*dtau*(a1**2*b1*r1*r4 - a2**2*b1*r2*r3 + a3**2*b2*r2*r3 - a4** &
      2*b2*r1*r4)
f(444) = 2*dtau*(a1*b1*r1*r4**2 - a2*b1*r2*r3**2 + a3*b2*r2**2*r3 - a4* &
      b2*r1**2*r4)
f(445) = 2*dtau*(a1*b1*r1**2*r4 - a2*b1*r2**2*r3 + a3*b2*r2*r3**2 - a4* &
      b2*r1*r4**2)
f(446) = 2*dtau*(a1*b2*r1*r4 - a2*b2*r2*r3 + a3*b1*r2*r3 - a4*b1*r1*r4)
f(447) = 2*dtau*(a1**2*b2*r1*r4 - a2**2*b2*r2*r3 + a3**2*b1*r2*r3 - a4** &
      2*b1*r1*r4)
f(448) = 2*dtau*(a1*b2*r1*r4**2 - a2*b2*r2*r3**2 + a3*b1*r2**2*r3 - a4* &
      b1*r1**2*r4)
f(449) = 2*dtau*(a1*b2*r1**2*r4 - a2*b2*r2**2*r3 + a3*b1*r2*r3**2 - a4* &
      b1*r1*r4**2)
f(450) = 2*a1*a4*b1**2*r2*r3 + 2*a1*a4*b2**2*r2*r3 + 2*a2*a3*b1**2*r1*r4 &
      + 2*a2*a3*b2**2*r1*r4
f(451) = 4*dtau**2*(a1*a4*r2*r3 + a2*a3*r1*r4)
f(452) = 2*b1*b2*(a1*r2*r3 + a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(453) = 2*b1*b2*(a1**2*r2*r3 + a2**2*r1*r4 + a3**2*r1*r4 + a4**2*r2*r3)
f(454) = 2*b1*b2*(a1*r2*r3**2 + a2*r1*r4**2 + a3*r1**2*r4 + a4*r2**2*r3)
f(455) = 2*b1*b2*(a1*r2**2*r3 + a2*r1**2*r4 + a3*r1*r4**2 + a4*r2*r3**2)
f(456) = 2*dtau*(-a1*b1*r2*r3 + a2*b1*r1*r4 - a3*b2*r1*r4 + a4*b2*r2*r3)
f(457) = 2*dtau*(-a1**2*b1*r2*r3 + a2**2*b1*r1*r4 - a3**2*b2*r1*r4 + a4 &
      **2*b2*r2*r3)
f(458) = 2*dtau*(-a1*b1*r2*r3**2 + a2*b1*r1*r4**2 - a3*b2*r1**2*r4 + a4* &
      b2*r2**2*r3)
f(459) = 2*dtau*(-a1*b1*r2**2*r3 + a2*b1*r1**2*r4 - a3*b2*r1*r4**2 + a4* &
      b2*r2*r3**2)
f(460) = 2*dtau*(-a1*b2*r2*r3 + a2*b2*r1*r4 - a3*b1*r1*r4 + a4*b1*r2*r3)
f(461) = 2*dtau*(-a1**2*b2*r2*r3 + a2**2*b2*r1*r4 - a3**2*b1*r1*r4 + a4 &
      **2*b1*r2*r3)
f(462) = 2*dtau*(-a1*b2*r2*r3**2 + a2*b2*r1*r4**2 - a3*b1*r1**2*r4 + a4* &
      b1*r2**2*r3)
f(463) = 2*dtau*(-a1*b2*r2**2*r3 + a2*b2*r1**2*r4 - a3*b1*r1*r4**2 + a4* &
      b1*r2*r3**2)
f(464) = 4*b1*b2*dtau**2*(r1*r4 + r2*r3)
f(465) = 2*b1*b2*dtau*(-b1*r1*r4 + b1*r2*r3 + b2*r1*r4 - b2*r2*r3)
f(466) = 2*b1*b2*dtau*(b1*r1*r4 - b1*r2*r3 - b2*r1*r4 + b2*r2*r3)
f(467) = 2*a1*a2*a3*a4*(r1 + r2 + r3 + r4)
f(468) = 2*a1*a2*a3*a4*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(469) = 2*a1*a2*a3*a4*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(470) = 2*a1*a2*a3*a4*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(471) = 2*a1*a2*a3*a4*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(472) = 2*a1*a2*a3*a4*(r1**2 + r2**2 + r3**2 + r4**2)
f(473) = 2*a1*a2*a3*b1**2*r1 + 2*a1*a2*a4*b1**2*r2 + 2*a1*a3*a4*b2**2*r3 &
      + 2*a2*a3*a4*b2**2*r4
f(474) = 2*a1*a2*a3*b2**2*r1 + 2*a1*a2*a4*b2**2*r2 + 2*a1*a3*a4*b1**2*r3 &
      + 2*a2*a3*a4*b1**2*r4
f(475) = 2*dtau**2*(a1*a2*a3*r1 + a1*a2*a4*r2 + a1*a3*a4*r3 + a2*a3*a4* &
      r4)
f(476) = 2*a1*a2*a3*b1**2*r2 + 2*a1*a2*a4*b1**2*r1 + 2*a1*a3*a4*b2**2*r4 &
      + 2*a2*a3*a4*b2**2*r3
f(477) = 2*a1*a2*a3*b2**2*r2 + 2*a1*a2*a4*b2**2*r1 + 2*a1*a3*a4*b1**2*r4 &
      + 2*a2*a3*a4*b1**2*r3
f(478) = 2*dtau**2*(a1*a2*a3*r2 + a1*a2*a4*r1 + a1*a3*a4*r4 + a2*a3*a4* &
      r3)
f(479) = 2*b1*b2*(a1*a2*r1 + a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(480) = 2*b1*b2*(a1**2*a2*r2 + a1*a2**2*r1 + a3**2*a4*r4 + a3*a4**2*r3)
f(481) = 2*b1*b2*(a1**2*a2*r1 + a1*a2**2*r2 + a3**2*a4*r3 + a3*a4**2*r4)
f(482) = 2*b1*b2*(a1*a2*r1**2 + a1*a2*r2**2 + a3*a4*r3**2 + a3*a4*r4**2)
f(483) = 2*dtau*(a1*a2*b1*r1 - a1*a2*b1*r2 + a3*a4*b2*r3 - a3*a4*b2*r4)
f(484) = 2*dtau*(-a1**2*a2*b1*r2 + a1*a2**2*b1*r1 - a3**2*a4*b2*r4 + a3* &
      a4**2*b2*r3)
f(485) = 2*dtau*(a1**2*a2*b1*r1 - a1*a2**2*b1*r2 + a3**2*a4*b2*r3 - a3* &
      a4**2*b2*r4)
f(486) = 2*dtau*(a1*a2*b1*r1**2 - a1*a2*b1*r2**2 + a3*a4*b2*r3**2 - a3* &
      a4*b2*r4**2)
f(487) = 2*dtau*(a1*a2*b2*r1 - a1*a2*b2*r2 + a3*a4*b1*r3 - a3*a4*b1*r4)
f(488) = 2*dtau*(-a1**2*a2*b2*r2 + a1*a2**2*b2*r1 - a3**2*a4*b1*r4 + a3* &
      a4**2*b1*r3)
f(489) = 2*dtau*(a1**2*a2*b2*r1 - a1*a2**2*b2*r2 + a3**2*a4*b1*r3 - a3* &
      a4**2*b1*r4)
f(490) = 2*dtau*(a1*a2*b2*r1**2 - a1*a2*b2*r2**2 + a3*a4*b1*r3**2 - a3* &
      a4*b1*r4**2)
f(491) = 2*a1*a2*a3*b2**2*r3 + 2*a1*a2*a4*b2**2*r4 + 2*a1*a3*a4*b1**2*r1 &
      + 2*a2*a3*a4*b1**2*r2
f(492) = 2*a1*a2*a3*b1**2*r3 + 2*a1*a2*a4*b1**2*r4 + 2*a1*a3*a4*b2**2*r1 &
      + 2*a2*a3*a4*b2**2*r2
f(493) = 2*dtau**2*(a1*a2*a3*r3 + a1*a2*a4*r4 + a1*a3*a4*r1 + a2*a3*a4* &
      r2)
f(494) = 2*b1*b2*(a1*a3*r1 + a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(495) = 2*b1*b2*(a1**2*a3*r3 + a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2*r2)
f(496) = 2*b1*b2*(a1**2*a3*r1 + a1*a3**2*r3 + a2**2*a4*r2 + a2*a4**2*r4)
f(497) = 2*b1*b2*(a1*a3*r1**2 + a1*a3*r3**2 + a2*a4*r2**2 + a2*a4*r4**2)
f(498) = 2*dtau*(a1*a3*b1*r1 + a1*a3*b2*r3 - a2*a4*b1*r2 - a2*a4*b2*r4)
f(499) = 2*dtau*(a1**2*a3*b2*r3 + a1*a3**2*b1*r1 - a2**2*a4*b2*r4 - a2* &
      a4**2*b1*r2)
f(500) = 2*dtau*(a1**2*a3*b1*r1 + a1*a3**2*b2*r3 - a2**2*a4*b1*r2 - a2* &
      a4**2*b2*r4)
f(501) = 2*dtau*(a1*a3*b1*r1**2 + a1*a3*b2*r3**2 - a2*a4*b1*r2**2 - a2* &
      a4*b2*r4**2)
f(502) = 2*dtau*(a1*a3*b1*r3 + a1*a3*b2*r1 - a2*a4*b1*r4 - a2*a4*b2*r2)
f(503) = 2*dtau*(a1**2*a3*b1*r3 + a1*a3**2*b2*r1 - a2**2*a4*b1*r4 - a2* &
      a4**2*b2*r2)
f(504) = 2*dtau*(a1**2*a3*b2*r1 + a1*a3**2*b1*r3 - a2**2*a4*b2*r2 - a2* &
      a4**2*b1*r4)
f(505) = 2*dtau*(a1*a3*b1*r3**2 + a1*a3*b2*r1**2 - a2*a4*b1*r4**2 - a2* &
      a4*b2*r2**2)
f(506) = 2*b1*b2*(a1*a4*r1 + a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(507) = 2*b1*b2*(a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 + a2*a3**2*r2)
f(508) = 2*b1*b2*(a1**2*a4*r1 + a1*a4**2*r4 + a2**2*a3*r2 + a2*a3**2*r3)
f(509) = 2*b1*b2*(a1*a4*r1**2 + a1*a4*r4**2 + a2*a3*r2**2 + a2*a3*r3**2)
f(510) = 2*dtau*(a1*a4*b1*r1 - a1*a4*b2*r4 - a2*a3*b1*r2 + a2*a3*b2*r3)
f(511) = 2*dtau*(-a1**2*a4*b2*r4 + a1*a4**2*b1*r1 + a2**2*a3*b2*r3 - a2* &
      a3**2*b1*r2)
f(512) = 2*dtau*(a1**2*a4*b1*r1 - a1*a4**2*b2*r4 - a2**2*a3*b1*r2 + a2* &
      a3**2*b2*r3)
f(513) = 2*dtau*(a1*a4*b1*r1**2 - a1*a4*b2*r4**2 - a2*a3*b1*r2**2 + a2* &
      a3*b2*r3**2)
f(514) = 2*dtau*(-a1*a4*b1*r4 + a1*a4*b2*r1 + a2*a3*b1*r3 - a2*a3*b2*r2)
f(515) = 2*dtau*(-a1**2*a4*b1*r4 + a1*a4**2*b2*r1 + a2**2*a3*b1*r3 - a2* &
      a3**2*b2*r2)
f(516) = 2*dtau*(a1**2*a4*b2*r1 - a1*a4**2*b1*r4 - a2**2*a3*b2*r2 + a2* &
      a3**2*b1*r3)
f(517) = 2*dtau*(-a1*a4*b1*r4**2 + a1*a4*b2*r1**2 + a2*a3*b1*r3**2 - a2* &
      a3*b2*r2**2)
f(518) = 2*b1*b2*dtau**2*(a1*r1 + a2*r2 + a3*r3 + a4*r4)
f(519) = 2*b1*b2*dtau*(a1*b2*r1 - a2*b2*r2 + a3*b1*r3 - a4*b1*r4)
f(520) = 2*b1*b2*dtau*(a1*b1*r1 - a2*b1*r2 + a3*b2*r3 - a4*b2*r4)
f(521) = 2*a1*a2*a3*b2**2*r4 + 2*a1*a2*a4*b2**2*r3 + 2*a1*a3*a4*b1**2*r2 &
      + 2*a2*a3*a4*b1**2*r1
f(522) = 2*a1*a2*a3*b1**2*r4 + 2*a1*a2*a4*b1**2*r3 + 2*a1*a3*a4*b2**2*r2 &
      + 2*a2*a3*a4*b2**2*r1
f(523) = 2*dtau**2*(a1*a2*a3*r4 + a1*a2*a4*r3 + a1*a3*a4*r2 + a2*a3*a4* &
      r1)
f(524) = 2*b1*b2*(a1*a4*r2 + a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(525) = 2*b1*b2*(a1**2*a4*r3 + a1*a4**2*r2 + a2**2*a3*r4 + a2*a3**2*r1)
f(526) = 2*b1*b2*(a1**2*a4*r2 + a1*a4**2*r3 + a2**2*a3*r1 + a2*a3**2*r4)
f(527) = 2*b1*b2*(a1*a4*r2**2 + a1*a4*r3**2 + a2*a3*r1**2 + a2*a3*r4**2)
f(528) = 2*dtau*(-a1*a4*b1*r2 + a1*a4*b2*r3 + a2*a3*b1*r1 - a2*a3*b2*r4)
f(529) = 2*dtau*(a1**2*a4*b2*r3 - a1*a4**2*b1*r2 - a2**2*a3*b2*r4 + a2* &
      a3**2*b1*r1)
f(530) = 2*dtau*(-a1**2*a4*b1*r2 + a1*a4**2*b2*r3 + a2**2*a3*b1*r1 - a2* &
      a3**2*b2*r4)
f(531) = 2*dtau*(-a1*a4*b1*r2**2 + a1*a4*b2*r3**2 + a2*a3*b1*r1**2 - a2* &
      a3*b2*r4**2)
f(532) = 2*dtau*(a1*a4*b1*r3 - a1*a4*b2*r2 - a2*a3*b1*r4 + a2*a3*b2*r1)
f(533) = 2*dtau*(a1**2*a4*b1*r3 - a1*a4**2*b2*r2 - a2**2*a3*b1*r4 + a2* &
      a3**2*b2*r1)
f(534) = 2*dtau*(-a1**2*a4*b2*r2 + a1*a4**2*b1*r3 + a2**2*a3*b2*r1 - a2* &
      a3**2*b1*r4)
f(535) = 2*dtau*(a1*a4*b1*r3**2 - a1*a4*b2*r2**2 - a2*a3*b1*r4**2 + a2* &
      a3*b2*r1**2)
f(536) = 2*b1*b2*(a1*a3*r2 + a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(537) = 2*b1*b2*(a1**2*a3*r4 + a1*a3**2*r2 + a2**2*a4*r3 + a2*a4**2*r1)
f(538) = 2*b1*b2*(a1**2*a3*r2 + a1*a3**2*r4 + a2**2*a4*r1 + a2*a4**2*r3)
f(539) = 2*b1*b2*(a1*a3*r2**2 + a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3**2)
f(540) = 2*dtau*(-a1*a3*b1*r2 - a1*a3*b2*r4 + a2*a4*b1*r1 + a2*a4*b2*r3)
f(541) = 2*dtau*(-a1**2*a3*b2*r4 - a1*a3**2*b1*r2 + a2**2*a4*b2*r3 + a2* &
      a4**2*b1*r1)
f(542) = 2*dtau*(-a1**2*a3*b1*r2 - a1*a3**2*b2*r4 + a2**2*a4*b1*r1 + a2* &
      a4**2*b2*r3)
f(543) = 2*dtau*(-a1*a3*b1*r2**2 - a1*a3*b2*r4**2 + a2*a4*b1*r1**2 + a2* &
      a4*b2*r3**2)
f(544) = 2*dtau*(-a1*a3*b1*r4 - a1*a3*b2*r2 + a2*a4*b1*r3 + a2*a4*b2*r1)
f(545) = 2*dtau*(-a1**2*a3*b1*r4 - a1*a3**2*b2*r2 + a2**2*a4*b1*r3 + a2* &
      a4**2*b2*r1)
f(546) = 2*dtau*(-a1**2*a3*b2*r2 - a1*a3**2*b1*r4 + a2**2*a4*b2*r1 + a2* &
      a4**2*b1*r3)
f(547) = 2*dtau*(-a1*a3*b1*r4**2 - a1*a3*b2*r2**2 + a2*a4*b1*r3**2 + a2* &
      a4*b2*r1**2)
f(548) = 2*b1*b2*dtau**2*(a1*r2 + a2*r1 + a3*r4 + a4*r3)
f(549) = 2*b1*b2*dtau*(-a1*b2*r2 + a2*b2*r1 - a3*b1*r4 + a4*b1*r3)
f(550) = 2*b1*b2*dtau*(-a1*b1*r2 + a2*b1*r1 - a3*b2*r4 + a4*b2*r3)
f(551) = 2*b1*b2*(a1*a2*r3 + a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(552) = 2*b1*b2*(a1**2*a2*r4 + a1*a2**2*r3 + a3**2*a4*r2 + a3*a4**2*r1)
f(553) = 2*b1*b2*(a1**2*a2*r3 + a1*a2**2*r4 + a3**2*a4*r1 + a3*a4**2*r2)
f(554) = 2*b1*b2*(a1*a2*r3**2 + a1*a2*r4**2 + a3*a4*r1**2 + a3*a4*r2**2)
f(555) = 2*dtau*(a1*a2*b2*r3 - a1*a2*b2*r4 + a3*a4*b1*r1 - a3*a4*b1*r2)
f(556) = 2*dtau*(-a1**2*a2*b2*r4 + a1*a2**2*b2*r3 - a3**2*a4*b1*r2 + a3* &
      a4**2*b1*r1)
f(557) = 2*dtau*(a1**2*a2*b2*r3 - a1*a2**2*b2*r4 + a3**2*a4*b1*r1 - a3* &
      a4**2*b1*r2)
f(558) = 2*dtau*(a1*a2*b2*r3**2 - a1*a2*b2*r4**2 + a3*a4*b1*r1**2 - a3* &
      a4*b1*r2**2)
f(559) = 2*dtau*(a1*a2*b1*r3 - a1*a2*b1*r4 + a3*a4*b2*r1 - a3*a4*b2*r2)
f(560) = 2*dtau*(-a1**2*a2*b1*r4 + a1*a2**2*b1*r3 - a3**2*a4*b2*r2 + a3* &
      a4**2*b2*r1)
f(561) = 2*dtau*(a1**2*a2*b1*r3 - a1*a2**2*b1*r4 + a3**2*a4*b2*r1 - a3* &
      a4**2*b2*r2)
f(562) = 2*dtau*(a1*a2*b1*r3**2 - a1*a2*b1*r4**2 + a3*a4*b2*r1**2 - a3* &
      a4*b2*r2**2)
f(563) = 2*b1*b2*dtau**2*(a1*r3 + a2*r4 + a3*r1 + a4*r2)
f(564) = 2*b1*b2*dtau*(a1*b1*r3 - a2*b1*r4 + a3*b2*r1 - a4*b2*r2)
f(565) = 2*b1*b2*dtau*(a1*b2*r3 - a2*b2*r4 + a3*b1*r1 - a4*b1*r2)
f(566) = 2*b1*b2*dtau**2*(a1*r4 + a2*r3 + a3*r2 + a4*r1)
f(567) = 2*b1*b2*dtau*(-a1*b1*r4 + a2*b1*r3 - a3*b2*r2 + a4*b2*r1)
f(568) = 2*b1*b2*dtau*(-a1*b2*r4 + a2*b2*r3 - a3*b1*r2 + a4*b1*r1)
f(569) = 4*a1*a2*a3*a4*(b1**2 + b2**2)
f(570) = 8*a1*a2*a3*a4*dtau**2
f(571) = 2*b1*b2*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(572) = 2*b1*b2*(a1**2*a3*a4 + a1*a2*a3**2 + a1*a2*a4**2 + a2**2*a3*a4)
f(573) = 2*b1*b2*(a1**2*a2*a4 + a1*a2**2*a3 + a1*a3*a4**2 + a2*a3**2*a4)
f(574) = 2*b1*b2*(a1**2*a2*a3 + a1*a2**2*a4 + a1*a3**2*a4 + a2*a3*a4**2)
f(575) = 2*dtau*(a1*a2*a3*b1 - a1*a2*a4*b1 + a1*a3*a4*b2 - a2*a3*a4*b2)
f(576) = 2*dtau*(a1**2*a3*a4*b2 + a1*a2*a3**2*b1 - a1*a2*a4**2*b1 - a2** &
      2*a3*a4*b2)
f(577) = 2*dtau*(-a1**2*a2*a4*b1 + a1*a2**2*a3*b1 + a1*a3*a4**2*b2 - a2* &
      a3**2*a4*b2)
f(578) = 2*dtau*(a1**2*a2*a3*b1 - a1*a2**2*a4*b1 + a1*a3**2*a4*b2 - a2* &
      a3*a4**2*b2)
f(579) = 2*dtau*(a1*a2*a3*b2 - a1*a2*a4*b2 + a1*a3*a4*b1 - a2*a3*a4*b1)
f(580) = 2*dtau*(a1**2*a3*a4*b1 + a1*a2*a3**2*b2 - a1*a2*a4**2*b2 - a2** &
      2*a3*a4*b1)
f(581) = 2*dtau*(-a1**2*a2*a4*b2 + a1*a2**2*a3*b2 + a1*a3*a4**2*b1 - a2* &
      a3**2*a4*b1)
f(582) = 2*dtau*(a1**2*a2*a3*b2 - a1*a2**2*a4*b2 + a1*a3**2*a4*b1 - a2* &
      a3*a4**2*b1)
f(583) = 4*b1*b2*dtau**2*(a1*a2 + a3*a4)
f(584) = 4*b1*b2*dtau**2*(a1*a3 + a2*a4)
f(585) = 2*b1*b2*dtau*(a1*a3*b1 + a1*a3*b2 - a2*a4*b1 - a2*a4*b2)
f(586) = 4*b1*b2*dtau**2*(a1*a4 + a2*a3)
f(587) = 2*b1*b2*dtau*(-a1*a4*b1 + a1*a4*b2 + a2*a3*b1 - a2*a3*b2)
f(588) = 2*b1*b2*dtau*(a1*a4*b1 - a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
v = sum(f*params)
end function c2h4_poten_n5_d6_ADF


!###############################################################################


! D2h-symmetry-adapted 6-mode 6-order expansion for C2H4 molecule

function c2h4_poten_n6_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(120)
type(adf_realq):: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(120)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 2*r0*r1*r2*r3*r4*(a1 + a2 + a3 + a4)
f(2) = 2*r0*(a1*a2*r1*r2*r3 + a1*a2*r1*r2*r4 + a3*a4*r1*r3*r4 + a3*a4*r2 &
      *r3*r4)
f(3) = 2*r0*(a1*a3*r1*r2*r3 + a1*a3*r1*r3*r4 + a2*a4*r1*r2*r4 + a2*a4*r2 &
      *r3*r4)
f(4) = 2*r0*(a1*a4*r1*r2*r3 + a1*a4*r2*r3*r4 + a2*a3*r1*r2*r4 + a2*a3*r1 &
      *r3*r4)
f(5) = 2*r0*(a1*a4*r1*r2*r4 + a1*a4*r1*r3*r4 + a2*a3*r1*r2*r3 + a2*a3*r2 &
      *r3*r4)
f(6) = 2*r0*(a1*a3*r1*r2*r4 + a1*a3*r2*r3*r4 + a2*a4*r1*r2*r3 + a2*a4*r1 &
      *r3*r4)
f(7) = 2*r0*(a1*a2*r1*r3*r4 + a1*a2*r2*r3*r4 + a3*a4*r1*r2*r3 + a3*a4*r1 &
      *r2*r4)
f(8) = 2*b1*b2*r0*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(9) = 2*dtau*r0*(b1*r1*r2*r3 - b1*r1*r2*r4 + b2*r1*r3*r4 - b2*r2*r3*r4)
f(10) = 2*dtau*r0*(b1*r1*r3*r4 - b1*r2*r3*r4 + b2*r1*r2*r3 - b2*r1*r2*r4 &
      )
f(11) = 2*r0*(a1*a2*a3*r1*r2 + a1*a2*a4*r1*r2 + a1*a3*a4*r3*r4 + a2*a3* &
      a4*r3*r4)
f(12) = 2*r0*(a1*a2*a3*r3*r4 + a1*a2*a4*r3*r4 + a1*a3*a4*r1*r2 + a2*a3* &
      a4*r1*r2)
f(13) = 2*b1*b2*r0*(a1*r1*r2 + a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(14) = 2*dtau*r0*(a1*b1*r1*r2 - a2*b1*r1*r2 + a3*b2*r3*r4 - a4*b2*r3*r4 &
      )
f(15) = 2*dtau*r0*(a1*b2*r1*r2 - a2*b2*r1*r2 + a3*b1*r3*r4 - a4*b1*r3*r4 &
      )
f(16) = 2*b1*b2*r0*(a1*r3*r4 + a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(17) = 2*dtau*r0*(a1*b2*r3*r4 - a2*b2*r3*r4 + a3*b1*r1*r2 - a4*b1*r1*r2 &
      )
f(18) = 2*dtau*r0*(a1*b1*r3*r4 - a2*b1*r3*r4 + a3*b2*r1*r2 - a4*b2*r1*r2 &
      )
f(19) = 2*r0*(a1*a2*a3*r1*r3 + a1*a2*a4*r2*r4 + a1*a3*a4*r1*r3 + a2*a3* &
      a4*r2*r4)
f(20) = 2*r0*(a1*a2*a3*r2*r4 + a1*a2*a4*r1*r3 + a1*a3*a4*r2*r4 + a2*a3* &
      a4*r1*r3)
f(21) = 2*b1*b2*r0*(a1*r1*r3 + a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(22) = 2*dtau*r0*(a1*b1*r1*r3 - a2*b1*r2*r4 + a3*b2*r1*r3 - a4*b2*r2*r4 &
      )
f(23) = 2*dtau*r0*(a1*b2*r1*r3 - a2*b2*r2*r4 + a3*b1*r1*r3 - a4*b1*r2*r4 &
      )
f(24) = 2*b1*b2*r0*(a1*r2*r4 + a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(25) = 2*dtau*r0*(-a1*b1*r2*r4 + a2*b1*r1*r3 - a3*b2*r2*r4 + a4*b2*r1* &
      r3)
f(26) = 2*dtau*r0*(-a1*b2*r2*r4 + a2*b2*r1*r3 - a3*b1*r2*r4 + a4*b1*r1* &
      r3)
f(27) = 2*r0*(a1*a2*a3*r1*r4 + a1*a2*a4*r2*r3 + a1*a3*a4*r2*r3 + a2*a3* &
      a4*r1*r4)
f(28) = 2*r0*(a1*a2*a3*r2*r3 + a1*a2*a4*r1*r4 + a1*a3*a4*r1*r4 + a2*a3* &
      a4*r2*r3)
f(29) = 2*b1*b2*r0*(a1*r1*r4 + a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(30) = 2*dtau*r0*(a1*b1*r1*r4 - a2*b1*r2*r3 + a3*b2*r2*r3 - a4*b2*r1*r4 &
      )
f(31) = 2*dtau*r0*(a1*b2*r1*r4 - a2*b2*r2*r3 + a3*b1*r2*r3 - a4*b1*r1*r4 &
      )
f(32) = 2*b1*b2*r0*(a1*r2*r3 + a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(33) = 2*dtau*r0*(-a1*b1*r2*r3 + a2*b1*r1*r4 - a3*b2*r1*r4 + a4*b2*r2* &
      r3)
f(34) = 2*dtau*r0*(-a1*b2*r2*r3 + a2*b2*r1*r4 - a3*b1*r1*r4 + a4*b1*r2* &
      r3)
f(35) = 2*a1*a2*a3*a4*r0*(r1 + r2 + r3 + r4)
f(36) = 2*b1*b2*r0*(a1*a2*r1 + a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(37) = 2*dtau*r0*(a1*a2*b1*r1 - a1*a2*b1*r2 + a3*a4*b2*r3 - a3*a4*b2*r4 &
      )
f(38) = 2*dtau*r0*(a1*a2*b2*r1 - a1*a2*b2*r2 + a3*a4*b1*r3 - a3*a4*b1*r4 &
      )
f(39) = 2*b1*b2*r0*(a1*a3*r1 + a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(40) = 2*dtau*r0*(a1*a3*b1*r1 + a1*a3*b2*r3 - a2*a4*b1*r2 - a2*a4*b2*r4 &
      )
f(41) = 2*dtau*r0*(a1*a3*b1*r3 + a1*a3*b2*r1 - a2*a4*b1*r4 - a2*a4*b2*r2 &
      )
f(42) = 2*b1*b2*r0*(a1*a4*r1 + a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(43) = 2*dtau*r0*(a1*a4*b1*r1 - a1*a4*b2*r4 - a2*a3*b1*r2 + a2*a3*b2*r3 &
      )
f(44) = 2*dtau*r0*(-a1*a4*b1*r4 + a1*a4*b2*r1 + a2*a3*b1*r3 - a2*a3*b2* &
      r2)
f(45) = 2*b1*b2*r0*(a1*a4*r2 + a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(46) = 2*dtau*r0*(-a1*a4*b1*r2 + a1*a4*b2*r3 + a2*a3*b1*r1 - a2*a3*b2* &
      r4)
f(47) = 2*dtau*r0*(a1*a4*b1*r3 - a1*a4*b2*r2 - a2*a3*b1*r4 + a2*a3*b2*r1 &
      )
f(48) = 2*b1*b2*r0*(a1*a3*r2 + a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(49) = 2*dtau*r0*(-a1*a3*b1*r2 - a1*a3*b2*r4 + a2*a4*b1*r1 + a2*a4*b2* &
      r3)
f(50) = 2*dtau*r0*(-a1*a3*b1*r4 - a1*a3*b2*r2 + a2*a4*b1*r3 + a2*a4*b2* &
      r1)
f(51) = 2*b1*b2*r0*(a1*a2*r3 + a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(52) = 2*dtau*r0*(a1*a2*b2*r3 - a1*a2*b2*r4 + a3*a4*b1*r1 - a3*a4*b1*r2 &
      )
f(53) = 2*dtau*r0*(a1*a2*b1*r3 - a1*a2*b1*r4 + a3*a4*b2*r1 - a3*a4*b2*r2 &
      )
f(54) = 2*b1*b2*r0*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(55) = 2*dtau*r0*(a1*a2*a3*b1 - a1*a2*a4*b1 + a1*a3*a4*b2 - a2*a3*a4*b2 &
      )
f(56) = 2*dtau*r0*(a1*a2*a3*b2 - a1*a2*a4*b2 + a1*a3*a4*b1 - a2*a3*a4*b1 &
      )
f(57) = 4*r1*r2*r3*r4*(a1*a2 + a3*a4)
f(58) = 4*r1*r2*r3*r4*(a1*a3 + a2*a4)
f(59) = 4*r1*r2*r3*r4*(a1*a4 + a2*a3)
f(60) = 8*b1*b2*r1*r2*r3*r4
f(61) = 2*a1*a2*a3*r1*r2*r3 + 2*a1*a2*a4*r1*r2*r4 + 2*a1*a3*a4*r1*r3*r4 &
      + 2*a2*a3*a4*r2*r3*r4
f(62) = 2*a1*a2*a3*r1*r2*r4 + 2*a1*a2*a4*r1*r2*r3 + 2*a1*a3*a4*r2*r3*r4 &
      + 2*a2*a3*a4*r1*r3*r4
f(63) = 2*a1*a2*a3*r1*r3*r4 + 2*a1*a2*a4*r2*r3*r4 + 2*a1*a3*a4*r1*r2*r3 &
      + 2*a2*a3*a4*r1*r2*r4
f(64) = 2*b1*b2*(a1*r1*r2*r3 + a2*r1*r2*r4 + a3*r1*r3*r4 + a4*r2*r3*r4)
f(65) = 2*dtau*(a1*b1*r1*r2*r3 - a2*b1*r1*r2*r4 + a3*b2*r1*r3*r4 - a4*b2 &
      *r2*r3*r4)
f(66) = 2*dtau*(a1*b2*r1*r2*r3 - a2*b2*r1*r2*r4 + a3*b1*r1*r3*r4 - a4*b1 &
      *r2*r3*r4)
f(67) = 2*a1*a2*a3*r2*r3*r4 + 2*a1*a2*a4*r1*r3*r4 + 2*a1*a3*a4*r1*r2*r4 &
      + 2*a2*a3*a4*r1*r2*r3
f(68) = 2*b1*b2*(a1*r1*r2*r4 + a2*r1*r2*r3 + a3*r2*r3*r4 + a4*r1*r3*r4)
f(69) = 2*dtau*(-a1*b1*r1*r2*r4 + a2*b1*r1*r2*r3 - a3*b2*r2*r3*r4 + a4* &
      b2*r1*r3*r4)
f(70) = 2*dtau*(-a1*b2*r1*r2*r4 + a2*b2*r1*r2*r3 - a3*b1*r2*r3*r4 + a4* &
      b1*r1*r3*r4)
f(71) = 2*b1*b2*(a1*r1*r3*r4 + a2*r2*r3*r4 + a3*r1*r2*r3 + a4*r1*r2*r4)
f(72) = 2*dtau*(a1*b2*r1*r3*r4 - a2*b2*r2*r3*r4 + a3*b1*r1*r2*r3 - a4*b1 &
      *r1*r2*r4)
f(73) = 2*dtau*(a1*b1*r1*r3*r4 - a2*b1*r2*r3*r4 + a3*b2*r1*r2*r3 - a4*b2 &
      *r1*r2*r4)
f(74) = 2*b1*b2*(a1*r2*r3*r4 + a2*r1*r3*r4 + a3*r1*r2*r4 + a4*r1*r2*r3)
f(75) = 2*dtau*(-a1*b2*r2*r3*r4 + a2*b2*r1*r3*r4 - a3*b1*r1*r2*r4 + a4* &
      b1*r1*r2*r3)
f(76) = 2*dtau*(-a1*b1*r2*r3*r4 + a2*b1*r1*r3*r4 - a3*b2*r1*r2*r4 + a4* &
      b2*r1*r2*r3)
f(77) = 4*a1*a2*a3*a4*(r1*r2 + r3*r4)
f(78) = 4*b1*b2*(a1*a2*r1*r2 + a3*a4*r3*r4)
f(79) = 2*b1*b2*(a1*a3*r1*r2 + a1*a3*r3*r4 + a2*a4*r1*r2 + a2*a4*r3*r4)
f(80) = 2*dtau*(a1*a3*b1*r1*r2 + a1*a3*b2*r3*r4 - a2*a4*b1*r1*r2 - a2*a4 &
      *b2*r3*r4)
f(81) = 2*dtau*(a1*a3*b1*r3*r4 + a1*a3*b2*r1*r2 - a2*a4*b1*r3*r4 - a2*a4 &
      *b2*r1*r2)
f(82) = 2*b1*b2*(a1*a4*r1*r2 + a1*a4*r3*r4 + a2*a3*r1*r2 + a2*a3*r3*r4)
f(83) = 2*dtau*(a1*a4*b1*r1*r2 - a1*a4*b2*r3*r4 - a2*a3*b1*r1*r2 + a2*a3 &
      *b2*r3*r4)
f(84) = 2*dtau*(-a1*a4*b1*r3*r4 + a1*a4*b2*r1*r2 + a2*a3*b1*r3*r4 - a2* &
      a3*b2*r1*r2)
f(85) = 4*b1*b2*(a1*a2*r3*r4 + a3*a4*r1*r2)
f(86) = 4*a1*a2*a3*a4*(r1*r3 + r2*r4)
f(87) = 2*b1*b2*(a1*a2*r1*r3 + a1*a2*r2*r4 + a3*a4*r1*r3 + a3*a4*r2*r4)
f(88) = 2*dtau*(a1*a2*b1*r1*r3 - a1*a2*b1*r2*r4 + a3*a4*b2*r1*r3 - a3*a4 &
      *b2*r2*r4)
f(89) = 2*dtau*(a1*a2*b2*r1*r3 - a1*a2*b2*r2*r4 + a3*a4*b1*r1*r3 - a3*a4 &
      *b1*r2*r4)
f(90) = 4*b1*b2*(a1*a3*r1*r3 + a2*a4*r2*r4)
f(91) = 2*dtau*(a1*a3*b1*r1*r3 + a1*a3*b2*r1*r3 - a2*a4*b1*r2*r4 - a2*a4 &
      *b2*r2*r4)
f(92) = 2*b1*b2*(a1*a4*r1*r3 + a1*a4*r2*r4 + a2*a3*r1*r3 + a2*a3*r2*r4)
f(93) = 2*dtau*(a1*a4*b1*r1*r3 - a1*a4*b2*r2*r4 - a2*a3*b1*r2*r4 + a2*a3 &
      *b2*r1*r3)
f(94) = 2*dtau*(-a1*a4*b1*r2*r4 + a1*a4*b2*r1*r3 + a2*a3*b1*r1*r3 - a2* &
      a3*b2*r2*r4)
f(95) = 4*b1*b2*(a1*a3*r2*r4 + a2*a4*r1*r3)
f(96) = 2*dtau*(-a1*a3*b1*r2*r4 - a1*a3*b2*r2*r4 + a2*a4*b1*r1*r3 + a2* &
      a4*b2*r1*r3)
f(97) = 4*a1*a2*a3*a4*(r1*r4 + r2*r3)
f(98) = 2*b1*b2*(a1*a2*r1*r4 + a1*a2*r2*r3 + a3*a4*r1*r4 + a3*a4*r2*r3)
f(99) = 2*dtau*(a1*a2*b1*r1*r4 - a1*a2*b1*r2*r3 - a3*a4*b2*r1*r4 + a3*a4 &
      *b2*r2*r3)
f(100) = 2*dtau*(a1*a2*b2*r1*r4 - a1*a2*b2*r2*r3 - a3*a4*b1*r1*r4 + a3* &
      a4*b1*r2*r3)
f(101) = 2*b1*b2*(a1*a3*r1*r4 + a1*a3*r2*r3 + a2*a4*r1*r4 + a2*a4*r2*r3)
f(102) = 2*dtau*(a1*a3*b1*r1*r4 + a1*a3*b2*r2*r3 - a2*a4*b1*r2*r3 - a2* &
      a4*b2*r1*r4)
f(103) = 2*dtau*(a1*a3*b1*r2*r3 + a1*a3*b2*r1*r4 - a2*a4*b1*r1*r4 - a2* &
      a4*b2*r2*r3)
f(104) = 4*b1*b2*(a1*a4*r1*r4 + a2*a3*r2*r3)
f(105) = 2*dtau*(a1*a4*b1*r1*r4 - a1*a4*b2*r1*r4 - a2*a3*b1*r2*r3 + a2* &
      a3*b2*r2*r3)
f(106) = 4*b1*b2*(a1*a4*r2*r3 + a2*a3*r1*r4)
f(107) = 2*dtau*(-a1*a4*b1*r2*r3 + a1*a4*b2*r2*r3 + a2*a3*b1*r1*r4 - a2* &
      a3*b2*r1*r4)
f(108) = 2*b1*b2*(a1*a2*a3*r1 + a1*a2*a4*r2 + a1*a3*a4*r3 + a2*a3*a4*r4)
f(109) = 2*dtau*(a1*a2*a3*b1*r1 - a1*a2*a4*b1*r2 + a1*a3*a4*b2*r3 - a2* &
      a3*a4*b2*r4)
f(110) = 2*dtau*(a1*a2*a3*b2*r1 - a1*a2*a4*b2*r2 + a1*a3*a4*b1*r3 - a2* &
      a3*a4*b1*r4)
f(111) = 2*b1*b2*(a1*a2*a3*r2 + a1*a2*a4*r1 + a1*a3*a4*r4 + a2*a3*a4*r3)
f(112) = 2*dtau*(-a1*a2*a3*b1*r2 + a1*a2*a4*b1*r1 - a1*a3*a4*b2*r4 + a2* &
      a3*a4*b2*r3)
f(113) = 2*dtau*(-a1*a2*a3*b2*r2 + a1*a2*a4*b2*r1 - a1*a3*a4*b1*r4 + a2* &
      a3*a4*b1*r3)
f(114) = 2*b1*b2*(a1*a2*a3*r3 + a1*a2*a4*r4 + a1*a3*a4*r1 + a2*a3*a4*r2)
f(115) = 2*dtau*(a1*a2*a3*b2*r3 - a1*a2*a4*b2*r4 + a1*a3*a4*b1*r1 - a2* &
      a3*a4*b1*r2)
f(116) = 2*dtau*(a1*a2*a3*b1*r3 - a1*a2*a4*b1*r4 + a1*a3*a4*b2*r1 - a2* &
      a3*a4*b2*r2)
f(117) = 2*b1*b2*(a1*a2*a3*r4 + a1*a2*a4*r3 + a1*a3*a4*r2 + a2*a3*a4*r1)
f(118) = 2*dtau*(-a1*a2*a3*b2*r4 + a1*a2*a4*b2*r3 - a1*a3*a4*b1*r2 + a2* &
      a3*a4*b1*r1)
f(119) = 2*dtau*(-a1*a2*a3*b1*r4 + a1*a2*a4*b1*r3 - a1*a3*a4*b2*r2 + a2* &
      a3*a4*b2*r1)
f(120) = 8*a1*a2*a3*a4*b1*b2
v = sum(f*params)
end function c2h4_poten_n6_d6_ADF


!###############################################################################


! D2h-symmetry-adapted 7-mode 7-order expansion for C2H4 molecule

function c2h4_poten_n7_d7_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(108)
type(adf_realq):: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(108)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 4*r0*r1*r2*r3*r4*(a1*a2 + a3*a4)
f(2) = 4*r0*r1*r2*r3*r4*(a1*a3 + a2*a4)
f(3) = 4*r0*r1*r2*r3*r4*(a1*a4 + a2*a3)
f(4) = 8*b1*b2*r0*r1*r2*r3*r4
f(5) = 2*r0*(a1*a2*a3*r1*r2*r3 + a1*a2*a4*r1*r2*r4 + a1*a3*a4*r1*r3*r4 + &
      a2*a3*a4*r2*r3*r4)
f(6) = 2*r0*(a1*a2*a3*r1*r2*r4 + a1*a2*a4*r1*r2*r3 + a1*a3*a4*r2*r3*r4 + &
      a2*a3*a4*r1*r3*r4)
f(7) = 2*r0*(a1*a2*a3*r1*r3*r4 + a1*a2*a4*r2*r3*r4 + a1*a3*a4*r1*r2*r3 + &
      a2*a3*a4*r1*r2*r4)
f(8) = 2*b1*b2*r0*(a1*r1*r2*r3 + a2*r1*r2*r4 + a3*r1*r3*r4 + a4*r2*r3*r4 &
      )
f(9) = 2*dtau*r0*(a1*b1*r1*r2*r3 - a2*b1*r1*r2*r4 + a3*b2*r1*r3*r4 - a4* &
      b2*r2*r3*r4)
f(10) = 2*dtau*r0*(a1*b2*r1*r2*r3 - a2*b2*r1*r2*r4 + a3*b1*r1*r3*r4 - a4 &
      *b1*r2*r3*r4)
f(11) = 2*r0*(a1*a2*a3*r2*r3*r4 + a1*a2*a4*r1*r3*r4 + a1*a3*a4*r1*r2*r4 &
      + a2*a3*a4*r1*r2*r3)
f(12) = 2*b1*b2*r0*(a1*r1*r2*r4 + a2*r1*r2*r3 + a3*r2*r3*r4 + a4*r1*r3* &
      r4)
f(13) = 2*dtau*r0*(-a1*b1*r1*r2*r4 + a2*b1*r1*r2*r3 - a3*b2*r2*r3*r4 + &
      a4*b2*r1*r3*r4)
f(14) = 2*dtau*r0*(-a1*b2*r1*r2*r4 + a2*b2*r1*r2*r3 - a3*b1*r2*r3*r4 + &
      a4*b1*r1*r3*r4)
f(15) = 2*b1*b2*r0*(a1*r1*r3*r4 + a2*r2*r3*r4 + a3*r1*r2*r3 + a4*r1*r2* &
      r4)
f(16) = 2*dtau*r0*(a1*b2*r1*r3*r4 - a2*b2*r2*r3*r4 + a3*b1*r1*r2*r3 - a4 &
      *b1*r1*r2*r4)
f(17) = 2*dtau*r0*(a1*b1*r1*r3*r4 - a2*b1*r2*r3*r4 + a3*b2*r1*r2*r3 - a4 &
      *b2*r1*r2*r4)
f(18) = 2*b1*b2*r0*(a1*r2*r3*r4 + a2*r1*r3*r4 + a3*r1*r2*r4 + a4*r1*r2* &
      r3)
f(19) = 2*dtau*r0*(-a1*b2*r2*r3*r4 + a2*b2*r1*r3*r4 - a3*b1*r1*r2*r4 + &
      a4*b1*r1*r2*r3)
f(20) = 2*dtau*r0*(-a1*b1*r2*r3*r4 + a2*b1*r1*r3*r4 - a3*b2*r1*r2*r4 + &
      a4*b2*r1*r2*r3)
f(21) = 4*a1*a2*a3*a4*r0*(r1*r2 + r3*r4)
f(22) = 4*b1*b2*r0*(a1*a2*r1*r2 + a3*a4*r3*r4)
f(23) = 2*b1*b2*r0*(a1*a3*r1*r2 + a1*a3*r3*r4 + a2*a4*r1*r2 + a2*a4*r3* &
      r4)
f(24) = 2*dtau*r0*(a1*a3*b1*r1*r2 + a1*a3*b2*r3*r4 - a2*a4*b1*r1*r2 - a2 &
      *a4*b2*r3*r4)
f(25) = 2*dtau*r0*(a1*a3*b1*r3*r4 + a1*a3*b2*r1*r2 - a2*a4*b1*r3*r4 - a2 &
      *a4*b2*r1*r2)
f(26) = 2*b1*b2*r0*(a1*a4*r1*r2 + a1*a4*r3*r4 + a2*a3*r1*r2 + a2*a3*r3* &
      r4)
f(27) = 2*dtau*r0*(a1*a4*b1*r1*r2 - a1*a4*b2*r3*r4 - a2*a3*b1*r1*r2 + a2 &
      *a3*b2*r3*r4)
f(28) = 2*dtau*r0*(-a1*a4*b1*r3*r4 + a1*a4*b2*r1*r2 + a2*a3*b1*r3*r4 - &
      a2*a3*b2*r1*r2)
f(29) = 4*b1*b2*r0*(a1*a2*r3*r4 + a3*a4*r1*r2)
f(30) = 4*a1*a2*a3*a4*r0*(r1*r3 + r2*r4)
f(31) = 2*b1*b2*r0*(a1*a2*r1*r3 + a1*a2*r2*r4 + a3*a4*r1*r3 + a3*a4*r2* &
      r4)
f(32) = 2*dtau*r0*(a1*a2*b1*r1*r3 - a1*a2*b1*r2*r4 + a3*a4*b2*r1*r3 - a3 &
      *a4*b2*r2*r4)
f(33) = 2*dtau*r0*(a1*a2*b2*r1*r3 - a1*a2*b2*r2*r4 + a3*a4*b1*r1*r3 - a3 &
      *a4*b1*r2*r4)
f(34) = 4*b1*b2*r0*(a1*a3*r1*r3 + a2*a4*r2*r4)
f(35) = 2*dtau*r0*(a1*a3*b1*r1*r3 + a1*a3*b2*r1*r3 - a2*a4*b1*r2*r4 - a2 &
      *a4*b2*r2*r4)
f(36) = 2*b1*b2*r0*(a1*a4*r1*r3 + a1*a4*r2*r4 + a2*a3*r1*r3 + a2*a3*r2* &
      r4)
f(37) = 2*dtau*r0*(a1*a4*b1*r1*r3 - a1*a4*b2*r2*r4 - a2*a3*b1*r2*r4 + a2 &
      *a3*b2*r1*r3)
f(38) = 2*dtau*r0*(-a1*a4*b1*r2*r4 + a1*a4*b2*r1*r3 + a2*a3*b1*r1*r3 - &
      a2*a3*b2*r2*r4)
f(39) = 4*b1*b2*r0*(a1*a3*r2*r4 + a2*a4*r1*r3)
f(40) = 2*dtau*r0*(-a1*a3*b1*r2*r4 - a1*a3*b2*r2*r4 + a2*a4*b1*r1*r3 + &
      a2*a4*b2*r1*r3)
f(41) = 4*a1*a2*a3*a4*r0*(r1*r4 + r2*r3)
f(42) = 2*b1*b2*r0*(a1*a2*r1*r4 + a1*a2*r2*r3 + a3*a4*r1*r4 + a3*a4*r2* &
      r3)
f(43) = 2*dtau*r0*(a1*a2*b1*r1*r4 - a1*a2*b1*r2*r3 - a3*a4*b2*r1*r4 + a3 &
      *a4*b2*r2*r3)
f(44) = 2*dtau*r0*(a1*a2*b2*r1*r4 - a1*a2*b2*r2*r3 - a3*a4*b1*r1*r4 + a3 &
      *a4*b1*r2*r3)
f(45) = 2*b1*b2*r0*(a1*a3*r1*r4 + a1*a3*r2*r3 + a2*a4*r1*r4 + a2*a4*r2* &
      r3)
f(46) = 2*dtau*r0*(a1*a3*b1*r1*r4 + a1*a3*b2*r2*r3 - a2*a4*b1*r2*r3 - a2 &
      *a4*b2*r1*r4)
f(47) = 2*dtau*r0*(a1*a3*b1*r2*r3 + a1*a3*b2*r1*r4 - a2*a4*b1*r1*r4 - a2 &
      *a4*b2*r2*r3)
f(48) = 4*b1*b2*r0*(a1*a4*r1*r4 + a2*a3*r2*r3)
f(49) = 2*dtau*r0*(a1*a4*b1*r1*r4 - a1*a4*b2*r1*r4 - a2*a3*b1*r2*r3 + a2 &
      *a3*b2*r2*r3)
f(50) = 4*b1*b2*r0*(a1*a4*r2*r3 + a2*a3*r1*r4)
f(51) = 2*dtau*r0*(-a1*a4*b1*r2*r3 + a1*a4*b2*r2*r3 + a2*a3*b1*r1*r4 - &
      a2*a3*b2*r1*r4)
f(52) = 2*b1*b2*r0*(a1*a2*a3*r1 + a1*a2*a4*r2 + a1*a3*a4*r3 + a2*a3*a4* &
      r4)
f(53) = 2*dtau*r0*(a1*a2*a3*b1*r1 - a1*a2*a4*b1*r2 + a1*a3*a4*b2*r3 - a2 &
      *a3*a4*b2*r4)
f(54) = 2*dtau*r0*(a1*a2*a3*b2*r1 - a1*a2*a4*b2*r2 + a1*a3*a4*b1*r3 - a2 &
      *a3*a4*b1*r4)
f(55) = 2*b1*b2*r0*(a1*a2*a3*r2 + a1*a2*a4*r1 + a1*a3*a4*r4 + a2*a3*a4* &
      r3)
f(56) = 2*dtau*r0*(-a1*a2*a3*b1*r2 + a1*a2*a4*b1*r1 - a1*a3*a4*b2*r4 + &
      a2*a3*a4*b2*r3)
f(57) = 2*dtau*r0*(-a1*a2*a3*b2*r2 + a1*a2*a4*b2*r1 - a1*a3*a4*b1*r4 + &
      a2*a3*a4*b1*r3)
f(58) = 2*b1*b2*r0*(a1*a2*a3*r3 + a1*a2*a4*r4 + a1*a3*a4*r1 + a2*a3*a4* &
      r2)
f(59) = 2*dtau*r0*(a1*a2*a3*b2*r3 - a1*a2*a4*b2*r4 + a1*a3*a4*b1*r1 - a2 &
      *a3*a4*b1*r2)
f(60) = 2*dtau*r0*(a1*a2*a3*b1*r3 - a1*a2*a4*b1*r4 + a1*a3*a4*b2*r1 - a2 &
      *a3*a4*b2*r2)
f(61) = 2*b1*b2*r0*(a1*a2*a3*r4 + a1*a2*a4*r3 + a1*a3*a4*r2 + a2*a3*a4* &
      r1)
f(62) = 2*dtau*r0*(-a1*a2*a3*b2*r4 + a1*a2*a4*b2*r3 - a1*a3*a4*b1*r2 + &
      a2*a3*a4*b1*r1)
f(63) = 2*dtau*r0*(-a1*a2*a3*b1*r4 + a1*a2*a4*b1*r3 - a1*a3*a4*b2*r2 + &
      a2*a3*a4*b2*r1)
f(64) = 8*a1*a2*a3*a4*b1*b2*r0
f(65) = 2*r1*r2*r3*r4*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(66) = 2*b1*b2*r1*r2*r3*r4*(a1 + a2 + a3 + a4)
f(67) = 2*dtau*r1*r2*r3*r4*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(68) = 2*dtau*r1*r2*r3*r4*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(69) = 2*a1*a2*a3*a4*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(70) = 2*b1*b2*(a1*a2*r1*r2*r3 + a1*a2*r1*r2*r4 + a3*a4*r1*r3*r4 + a3* &
      a4*r2*r3*r4)
f(71) = 2*dtau*(a1*a2*b1*r1*r2*r3 - a1*a2*b1*r1*r2*r4 + a3*a4*b2*r1*r3* &
      r4 - a3*a4*b2*r2*r3*r4)
f(72) = 2*dtau*(a1*a2*b2*r1*r2*r3 - a1*a2*b2*r1*r2*r4 + a3*a4*b1*r1*r3* &
      r4 - a3*a4*b1*r2*r3*r4)
f(73) = 2*b1*b2*(a1*a3*r1*r2*r3 + a1*a3*r1*r3*r4 + a2*a4*r1*r2*r4 + a2* &
      a4*r2*r3*r4)
f(74) = 2*dtau*(a1*a3*b1*r1*r2*r3 + a1*a3*b2*r1*r3*r4 - a2*a4*b1*r1*r2* &
      r4 - a2*a4*b2*r2*r3*r4)
f(75) = 2*dtau*(a1*a3*b1*r1*r3*r4 + a1*a3*b2*r1*r2*r3 - a2*a4*b1*r2*r3* &
      r4 - a2*a4*b2*r1*r2*r4)
f(76) = 2*b1*b2*(a1*a4*r1*r2*r3 + a1*a4*r2*r3*r4 + a2*a3*r1*r2*r4 + a2* &
      a3*r1*r3*r4)
f(77) = 2*dtau*(a1*a4*b1*r1*r2*r3 - a1*a4*b2*r2*r3*r4 - a2*a3*b1*r1*r2* &
      r4 + a2*a3*b2*r1*r3*r4)
f(78) = 2*dtau*(-a1*a4*b1*r2*r3*r4 + a1*a4*b2*r1*r2*r3 + a2*a3*b1*r1*r3* &
      r4 - a2*a3*b2*r1*r2*r4)
f(79) = 2*b1*b2*(a1*a4*r1*r2*r4 + a1*a4*r1*r3*r4 + a2*a3*r1*r2*r3 + a2* &
      a3*r2*r3*r4)
f(80) = 2*dtau*(-a1*a4*b1*r1*r2*r4 + a1*a4*b2*r1*r3*r4 + a2*a3*b1*r1*r2* &
      r3 - a2*a3*b2*r2*r3*r4)
f(81) = 2*dtau*(a1*a4*b1*r1*r3*r4 - a1*a4*b2*r1*r2*r4 - a2*a3*b1*r2*r3* &
      r4 + a2*a3*b2*r1*r2*r3)
f(82) = 2*b1*b2*(a1*a3*r1*r2*r4 + a1*a3*r2*r3*r4 + a2*a4*r1*r2*r3 + a2* &
      a4*r1*r3*r4)
f(83) = 2*dtau*(-a1*a3*b1*r1*r2*r4 - a1*a3*b2*r2*r3*r4 + a2*a4*b1*r1*r2* &
      r3 + a2*a4*b2*r1*r3*r4)
f(84) = 2*dtau*(-a1*a3*b1*r2*r3*r4 - a1*a3*b2*r1*r2*r4 + a2*a4*b1*r1*r3* &
      r4 + a2*a4*b2*r1*r2*r3)
f(85) = 2*b1*b2*(a1*a2*r1*r3*r4 + a1*a2*r2*r3*r4 + a3*a4*r1*r2*r3 + a3* &
      a4*r1*r2*r4)
f(86) = 2*dtau*(a1*a2*b2*r1*r3*r4 - a1*a2*b2*r2*r3*r4 + a3*a4*b1*r1*r2* &
      r3 - a3*a4*b1*r1*r2*r4)
f(87) = 2*dtau*(a1*a2*b1*r1*r3*r4 - a1*a2*b1*r2*r3*r4 + a3*a4*b2*r1*r2* &
      r3 - a3*a4*b2*r1*r2*r4)
f(88) = 2*b1*b2*(a1*a2*a3*r1*r2 + a1*a2*a4*r1*r2 + a1*a3*a4*r3*r4 + a2* &
      a3*a4*r3*r4)
f(89) = 2*dtau*(a1*a2*a3*b1*r1*r2 - a1*a2*a4*b1*r1*r2 + a1*a3*a4*b2*r3* &
      r4 - a2*a3*a4*b2*r3*r4)
f(90) = 2*dtau*(a1*a2*a3*b2*r1*r2 - a1*a2*a4*b2*r1*r2 + a1*a3*a4*b1*r3* &
      r4 - a2*a3*a4*b1*r3*r4)
f(91) = 2*b1*b2*(a1*a2*a3*r3*r4 + a1*a2*a4*r3*r4 + a1*a3*a4*r1*r2 + a2* &
      a3*a4*r1*r2)
f(92) = 2*dtau*(a1*a2*a3*b2*r3*r4 - a1*a2*a4*b2*r3*r4 + a1*a3*a4*b1*r1* &
      r2 - a2*a3*a4*b1*r1*r2)
f(93) = 2*dtau*(a1*a2*a3*b1*r3*r4 - a1*a2*a4*b1*r3*r4 + a1*a3*a4*b2*r1* &
      r2 - a2*a3*a4*b2*r1*r2)
f(94) = 2*b1*b2*(a1*a2*a3*r1*r3 + a1*a2*a4*r2*r4 + a1*a3*a4*r1*r3 + a2* &
      a3*a4*r2*r4)
f(95) = 2*dtau*(a1*a2*a3*b1*r1*r3 - a1*a2*a4*b1*r2*r4 + a1*a3*a4*b2*r1* &
      r3 - a2*a3*a4*b2*r2*r4)
f(96) = 2*dtau*(a1*a2*a3*b2*r1*r3 - a1*a2*a4*b2*r2*r4 + a1*a3*a4*b1*r1* &
      r3 - a2*a3*a4*b1*r2*r4)
f(97) = 2*b1*b2*(a1*a2*a3*r2*r4 + a1*a2*a4*r1*r3 + a1*a3*a4*r2*r4 + a2* &
      a3*a4*r1*r3)
f(98) = 2*dtau*(-a1*a2*a3*b1*r2*r4 + a1*a2*a4*b1*r1*r3 - a1*a3*a4*b2*r2* &
      r4 + a2*a3*a4*b2*r1*r3)
f(99) = 2*dtau*(-a1*a2*a3*b2*r2*r4 + a1*a2*a4*b2*r1*r3 - a1*a3*a4*b1*r2* &
      r4 + a2*a3*a4*b1*r1*r3)
f(100) = 2*b1*b2*(a1*a2*a3*r1*r4 + a1*a2*a4*r2*r3 + a1*a3*a4*r2*r3 + a2* &
      a3*a4*r1*r4)
f(101) = 2*dtau*(a1*a2*a3*b1*r1*r4 - a1*a2*a4*b1*r2*r3 + a1*a3*a4*b2*r2* &
      r3 - a2*a3*a4*b2*r1*r4)
f(102) = 2*dtau*(a1*a2*a3*b2*r1*r4 - a1*a2*a4*b2*r2*r3 + a1*a3*a4*b1*r2* &
      r3 - a2*a3*a4*b1*r1*r4)
f(103) = 2*b1*b2*(a1*a2*a3*r2*r3 + a1*a2*a4*r1*r4 + a1*a3*a4*r1*r4 + a2* &
      a3*a4*r2*r3)
f(104) = 2*dtau*(-a1*a2*a3*b1*r2*r3 + a1*a2*a4*b1*r1*r4 - a1*a3*a4*b2*r1 &
      *r4 + a2*a3*a4*b2*r2*r3)
f(105) = 2*dtau*(-a1*a2*a3*b2*r2*r3 + a1*a2*a4*b2*r1*r4 - a1*a3*a4*b1*r1 &
      *r4 + a2*a3*a4*b1*r2*r3)
f(106) = 2*a1*a2*a3*a4*b1*b2*(r1 + r2 + r3 + r4)
f(107) = 2*a1*a2*a3*a4*dtau*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(108) = 2*a1*a2*a3*a4*dtau*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
v = sum(f*params)
end function c2h4_poten_n7_d7_ADF


!###############################################################################


! D2h-symmetry-adapted 8-mode 8-order expansion for C2H4 molecule

function c2h4_poten_n8_d8_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(67)
type(adf_realq):: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(67)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 2*r0*r1*r2*r3*r4*(a1*a2*a3 + a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(2) = 2*b1*b2*r0*r1*r2*r3*r4*(a1 + a2 + a3 + a4)
f(3) = 2*dtau*r0*r1*r2*r3*r4*(a1*b1 - a2*b1 + a3*b2 - a4*b2)
f(4) = 2*dtau*r0*r1*r2*r3*r4*(a1*b2 - a2*b2 + a3*b1 - a4*b1)
f(5) = 2*a1*a2*a3*a4*r0*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(6) = 2*b1*b2*r0*(a1*a2*r1*r2*r3 + a1*a2*r1*r2*r4 + a3*a4*r1*r3*r4 + a3 &
      *a4*r2*r3*r4)
f(7) = 2*dtau*r0*(a1*a2*b1*r1*r2*r3 - a1*a2*b1*r1*r2*r4 + a3*a4*b2*r1*r3 &
      *r4 - a3*a4*b2*r2*r3*r4)
f(8) = 2*dtau*r0*(a1*a2*b2*r1*r2*r3 - a1*a2*b2*r1*r2*r4 + a3*a4*b1*r1*r3 &
      *r4 - a3*a4*b1*r2*r3*r4)
f(9) = 2*b1*b2*r0*(a1*a3*r1*r2*r3 + a1*a3*r1*r3*r4 + a2*a4*r1*r2*r4 + a2 &
      *a4*r2*r3*r4)
f(10) = 2*dtau*r0*(a1*a3*b1*r1*r2*r3 + a1*a3*b2*r1*r3*r4 - a2*a4*b1*r1* &
      r2*r4 - a2*a4*b2*r2*r3*r4)
f(11) = 2*dtau*r0*(a1*a3*b1*r1*r3*r4 + a1*a3*b2*r1*r2*r3 - a2*a4*b1*r2* &
      r3*r4 - a2*a4*b2*r1*r2*r4)
f(12) = 2*b1*b2*r0*(a1*a4*r1*r2*r3 + a1*a4*r2*r3*r4 + a2*a3*r1*r2*r4 + &
      a2*a3*r1*r3*r4)
f(13) = 2*dtau*r0*(a1*a4*b1*r1*r2*r3 - a1*a4*b2*r2*r3*r4 - a2*a3*b1*r1* &
      r2*r4 + a2*a3*b2*r1*r3*r4)
f(14) = 2*dtau*r0*(-a1*a4*b1*r2*r3*r4 + a1*a4*b2*r1*r2*r3 + a2*a3*b1*r1* &
      r3*r4 - a2*a3*b2*r1*r2*r4)
f(15) = 2*b1*b2*r0*(a1*a4*r1*r2*r4 + a1*a4*r1*r3*r4 + a2*a3*r1*r2*r3 + &
      a2*a3*r2*r3*r4)
f(16) = 2*dtau*r0*(-a1*a4*b1*r1*r2*r4 + a1*a4*b2*r1*r3*r4 + a2*a3*b1*r1* &
      r2*r3 - a2*a3*b2*r2*r3*r4)
f(17) = 2*dtau*r0*(a1*a4*b1*r1*r3*r4 - a1*a4*b2*r1*r2*r4 - a2*a3*b1*r2* &
      r3*r4 + a2*a3*b2*r1*r2*r3)
f(18) = 2*b1*b2*r0*(a1*a3*r1*r2*r4 + a1*a3*r2*r3*r4 + a2*a4*r1*r2*r3 + &
      a2*a4*r1*r3*r4)
f(19) = 2*dtau*r0*(-a1*a3*b1*r1*r2*r4 - a1*a3*b2*r2*r3*r4 + a2*a4*b1*r1* &
      r2*r3 + a2*a4*b2*r1*r3*r4)
f(20) = 2*dtau*r0*(-a1*a3*b1*r2*r3*r4 - a1*a3*b2*r1*r2*r4 + a2*a4*b1*r1* &
      r3*r4 + a2*a4*b2*r1*r2*r3)
f(21) = 2*b1*b2*r0*(a1*a2*r1*r3*r4 + a1*a2*r2*r3*r4 + a3*a4*r1*r2*r3 + &
      a3*a4*r1*r2*r4)
f(22) = 2*dtau*r0*(a1*a2*b2*r1*r3*r4 - a1*a2*b2*r2*r3*r4 + a3*a4*b1*r1* &
      r2*r3 - a3*a4*b1*r1*r2*r4)
f(23) = 2*dtau*r0*(a1*a2*b1*r1*r3*r4 - a1*a2*b1*r2*r3*r4 + a3*a4*b2*r1* &
      r2*r3 - a3*a4*b2*r1*r2*r4)
f(24) = 2*b1*b2*r0*(a1*a2*a3*r1*r2 + a1*a2*a4*r1*r2 + a1*a3*a4*r3*r4 + &
      a2*a3*a4*r3*r4)
f(25) = 2*dtau*r0*(a1*a2*a3*b1*r1*r2 - a1*a2*a4*b1*r1*r2 + a1*a3*a4*b2* &
      r3*r4 - a2*a3*a4*b2*r3*r4)
f(26) = 2*dtau*r0*(a1*a2*a3*b2*r1*r2 - a1*a2*a4*b2*r1*r2 + a1*a3*a4*b1* &
      r3*r4 - a2*a3*a4*b1*r3*r4)
f(27) = 2*b1*b2*r0*(a1*a2*a3*r3*r4 + a1*a2*a4*r3*r4 + a1*a3*a4*r1*r2 + &
      a2*a3*a4*r1*r2)
f(28) = 2*dtau*r0*(a1*a2*a3*b2*r3*r4 - a1*a2*a4*b2*r3*r4 + a1*a3*a4*b1* &
      r1*r2 - a2*a3*a4*b1*r1*r2)
f(29) = 2*dtau*r0*(a1*a2*a3*b1*r3*r4 - a1*a2*a4*b1*r3*r4 + a1*a3*a4*b2* &
      r1*r2 - a2*a3*a4*b2*r1*r2)
f(30) = 2*b1*b2*r0*(a1*a2*a3*r1*r3 + a1*a2*a4*r2*r4 + a1*a3*a4*r1*r3 + &
      a2*a3*a4*r2*r4)
f(31) = 2*dtau*r0*(a1*a2*a3*b1*r1*r3 - a1*a2*a4*b1*r2*r4 + a1*a3*a4*b2* &
      r1*r3 - a2*a3*a4*b2*r2*r4)
f(32) = 2*dtau*r0*(a1*a2*a3*b2*r1*r3 - a1*a2*a4*b2*r2*r4 + a1*a3*a4*b1* &
      r1*r3 - a2*a3*a4*b1*r2*r4)
f(33) = 2*b1*b2*r0*(a1*a2*a3*r2*r4 + a1*a2*a4*r1*r3 + a1*a3*a4*r2*r4 + &
      a2*a3*a4*r1*r3)
f(34) = 2*dtau*r0*(-a1*a2*a3*b1*r2*r4 + a1*a2*a4*b1*r1*r3 - a1*a3*a4*b2* &
      r2*r4 + a2*a3*a4*b2*r1*r3)
f(35) = 2*dtau*r0*(-a1*a2*a3*b2*r2*r4 + a1*a2*a4*b2*r1*r3 - a1*a3*a4*b1* &
      r2*r4 + a2*a3*a4*b1*r1*r3)
f(36) = 2*b1*b2*r0*(a1*a2*a3*r1*r4 + a1*a2*a4*r2*r3 + a1*a3*a4*r2*r3 + &
      a2*a3*a4*r1*r4)
f(37) = 2*dtau*r0*(a1*a2*a3*b1*r1*r4 - a1*a2*a4*b1*r2*r3 + a1*a3*a4*b2* &
      r2*r3 - a2*a3*a4*b2*r1*r4)
f(38) = 2*dtau*r0*(a1*a2*a3*b2*r1*r4 - a1*a2*a4*b2*r2*r3 + a1*a3*a4*b1* &
      r2*r3 - a2*a3*a4*b1*r1*r4)
f(39) = 2*b1*b2*r0*(a1*a2*a3*r2*r3 + a1*a2*a4*r1*r4 + a1*a3*a4*r1*r4 + &
      a2*a3*a4*r2*r3)
f(40) = 2*dtau*r0*(-a1*a2*a3*b1*r2*r3 + a1*a2*a4*b1*r1*r4 - a1*a3*a4*b2* &
      r1*r4 + a2*a3*a4*b2*r2*r3)
f(41) = 2*dtau*r0*(-a1*a2*a3*b2*r2*r3 + a1*a2*a4*b2*r1*r4 - a1*a3*a4*b1* &
      r1*r4 + a2*a3*a4*b1*r2*r3)
f(42) = 2*a1*a2*a3*a4*b1*b2*r0*(r1 + r2 + r3 + r4)
f(43) = 2*a1*a2*a3*a4*dtau*r0*(b1*r1 - b1*r2 + b2*r3 - b2*r4)
f(44) = 2*a1*a2*a3*a4*dtau*r0*(b1*r3 - b1*r4 + b2*r1 - b2*r2)
f(45) = 8*a1*a2*a3*a4*r1*r2*r3*r4
f(46) = 4*b1*b2*r1*r2*r3*r4*(a1*a2 + a3*a4)
f(47) = 4*b1*b2*r1*r2*r3*r4*(a1*a3 + a2*a4)
f(48) = 2*dtau*r1*r2*r3*r4*(a1*a3*b1 + a1*a3*b2 - a2*a4*b1 - a2*a4*b2)
f(49) = 4*b1*b2*r1*r2*r3*r4*(a1*a4 + a2*a3)
f(50) = 2*dtau*r1*r2*r3*r4*(a1*a4*b1 - a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
f(51) = 2*b1*b2*(a1*a2*a3*r1*r2*r3 + a1*a2*a4*r1*r2*r4 + a1*a3*a4*r1*r3* &
      r4 + a2*a3*a4*r2*r3*r4)
f(52) = 2*dtau*(a1*a2*a3*b1*r1*r2*r3 - a1*a2*a4*b1*r1*r2*r4 + a1*a3*a4* &
      b2*r1*r3*r4 - a2*a3*a4*b2*r2*r3*r4)
f(53) = 2*dtau*(a1*a2*a3*b2*r1*r2*r3 - a1*a2*a4*b2*r1*r2*r4 + a1*a3*a4* &
      b1*r1*r3*r4 - a2*a3*a4*b1*r2*r3*r4)
f(54) = 2*b1*b2*(a1*a2*a3*r1*r2*r4 + a1*a2*a4*r1*r2*r3 + a1*a3*a4*r2*r3* &
      r4 + a2*a3*a4*r1*r3*r4)
f(55) = 2*dtau*(-a1*a2*a3*b1*r1*r2*r4 + a1*a2*a4*b1*r1*r2*r3 - a1*a3*a4* &
      b2*r2*r3*r4 + a2*a3*a4*b2*r1*r3*r4)
f(56) = 2*dtau*(-a1*a2*a3*b2*r1*r2*r4 + a1*a2*a4*b2*r1*r2*r3 - a1*a3*a4* &
      b1*r2*r3*r4 + a2*a3*a4*b1*r1*r3*r4)
f(57) = 2*b1*b2*(a1*a2*a3*r1*r3*r4 + a1*a2*a4*r2*r3*r4 + a1*a3*a4*r1*r2* &
      r3 + a2*a3*a4*r1*r2*r4)
f(58) = 2*dtau*(a1*a2*a3*b2*r1*r3*r4 - a1*a2*a4*b2*r2*r3*r4 + a1*a3*a4* &
      b1*r1*r2*r3 - a2*a3*a4*b1*r1*r2*r4)
f(59) = 2*dtau*(a1*a2*a3*b1*r1*r3*r4 - a1*a2*a4*b1*r2*r3*r4 + a1*a3*a4* &
      b2*r1*r2*r3 - a2*a3*a4*b2*r1*r2*r4)
f(60) = 2*b1*b2*(a1*a2*a3*r2*r3*r4 + a1*a2*a4*r1*r3*r4 + a1*a3*a4*r1*r2* &
      r4 + a2*a3*a4*r1*r2*r3)
f(61) = 2*dtau*(-a1*a2*a3*b2*r2*r3*r4 + a1*a2*a4*b2*r1*r3*r4 - a1*a3*a4* &
      b1*r1*r2*r4 + a2*a3*a4*b1*r1*r2*r3)
f(62) = 2*dtau*(-a1*a2*a3*b1*r2*r3*r4 + a1*a2*a4*b1*r1*r3*r4 - a1*a3*a4* &
      b2*r1*r2*r4 + a2*a3*a4*b2*r1*r2*r3)
f(63) = 4*a1*a2*a3*a4*b1*b2*(r1*r2 + r3*r4)
f(64) = 4*a1*a2*a3*a4*b1*b2*(r1*r3 + r2*r4)
f(65) = 2*a1*a2*a3*a4*dtau*(b1*r1*r3 - b1*r2*r4 + b2*r1*r3 - b2*r2*r4)
f(66) = 4*a1*a2*a3*a4*b1*b2*(r1*r4 + r2*r3)
f(67) = 2*a1*a2*a3*a4*dtau*(b1*r1*r4 - b1*r2*r3 - b2*r1*r4 + b2*r2*r3)
v = sum(f*params)
end function c2h4_poten_n8_d8_ADF


!###############################################################################


function c2h4_dipole_b1u_n1_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(15)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(15)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = -r1 - r2 + r3 + r4
f(2) = r1**2 + r2**2 - r3**2 - r4**2
f(3) = r1**3 + r2**3 - r3**3 - r4**3
f(4) = -r1**4 - r2**4 + r3**4 + r4**4
f(5) = -r1**5 - r2**5 + r3**5 + r4**5
f(6) = -r1**6 - r2**6 + r3**6 + r4**6
f(7) = -a1 - a2 + a3 + a4
f(8) = -a1**2 - a2**2 + a3**2 + a4**2
f(9) = a1**3 + a2**3 - a3**3 - a4**3
f(10) = -a1**4 - a2**4 + a3**4 + a4**4
f(11) = a1**5 + a2**5 - a3**5 - a4**5
f(12) = a1**6 + a2**6 - a3**6 - a4**6
f(13) = b1**2 - b2**2
f(14) = b1**4 - b2**4
f(15) = b1**6 - b2**6
v = sum(f*params)
end function c2h4_dipole_b1u_n1_d6_ADF


!###############################################################################


function c2h4_dipole_b1u_n2_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(188)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(188)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(-r1 - r2 + r3 + r4)
f(2) = r0*(-r1**2 - r2**2 + r3**2 + r4**2)
f(3) = r0*(-r1**3 - r2**3 + r3**3 + r4**3)
f(4) = r0*(r1**4 + r2**4 - r3**4 - r4**4)
f(5) = r0*(r1**5 + r2**5 - r3**5 - r4**5)
f(6) = r0**2*(r1 + r2 - r3 - r4)
f(7) = r0**2*(r1**2 + r2**2 - r3**2 - r4**2)
f(8) = r0**2*(r1**3 + r2**3 - r3**3 - r4**3)
f(9) = r0**2*(-r1**4 - r2**4 + r3**4 + r4**4)
f(10) = r0**3*(-r1 - r2 + r3 + r4)
f(11) = r0**3*(-r1**2 - r2**2 + r3**2 + r4**2)
f(12) = r0**3*(-r1**3 - r2**3 + r3**3 + r4**3)
f(13) = r0**4*(r1 + r2 - r3 - r4)
f(14) = r0**4*(-r1**2 - r2**2 + r3**2 + r4**2)
f(15) = r0**5*(r1 + r2 - r3 - r4)
f(16) = r0*(-a1 - a2 + a3 + a4)
f(17) = r0*(-a1**2 - a2**2 + a3**2 + a4**2)
f(18) = r0*(a1**3 + a2**3 - a3**3 - a4**3)
f(19) = r0*(-a1**4 - a2**4 + a3**4 + a4**4)
f(20) = r0*(a1**5 + a2**5 - a3**5 - a4**5)
f(21) = r0**2*(-a1 - a2 + a3 + a4)
f(22) = r0**2*(-a1**2 - a2**2 + a3**2 + a4**2)
f(23) = r0**2*(a1**3 + a2**3 - a3**3 - a4**3)
f(24) = r0**2*(-a1**4 - a2**4 + a3**4 + a4**4)
f(25) = r0**3*(a1 + a2 - a3 - a4)
f(26) = r0**3*(a1**2 + a2**2 - a3**2 - a4**2)
f(27) = r0**3*(-a1**3 - a2**3 + a3**3 + a4**3)
f(28) = r0**4*(a1 + a2 - a3 - a4)
f(29) = r0**4*(a1**2 + a2**2 - a3**2 - a4**2)
f(30) = r0**5*(a1 + a2 - a3 - a4)
f(31) = r0*(b1**2 - b2**2)
f(32) = r0*(b1**4 - b2**4)
f(33) = r0**2*(-b1**2 + b2**2)
f(34) = r0**2*(b1**4 - b2**4)
f(35) = r0**3*(-b1**2 + b2**2)
f(36) = r0**4*(-b1**2 + b2**2)
f(37) = -r1*r2 + r3*r4
f(38) = -r1**2*r2 - r1*r2**2 + r3**2*r4 + r3*r4**2
f(39) = -r1**3*r2 - r1*r2**3 + r3**3*r4 + r3*r4**3
f(40) = -r1**4*r2 - r1*r2**4 + r3**4*r4 + r3*r4**4
f(41) = -r1**5*r2 - r1*r2**5 + r3**5*r4 + r3*r4**5
f(42) = r1**2*r2**2 - r3**2*r4**2
f(43) = r1**3*r2**2 + r1**2*r2**3 - r3**3*r4**2 - r3**2*r4**3
f(44) = r1**4*r2**2 + r1**2*r2**4 - r3**4*r4**2 - r3**2*r4**4
f(45) = -r1**3*r2**2 - r1**2*r2**3 + r3**3*r4**2 + r3**2*r4**3
f(46) = r1**3*r2**3 - r3**3*r4**3
f(47) = -r1**4*r2**2 - r1**2*r2**4 + r3**4*r4**2 + r3**2*r4**4
f(48) = -r1**2*r3 + r1*r3**2 - r2**2*r4 + r2*r4**2
f(49) = r1**3*r3 - r1*r3**3 + r2**3*r4 - r2*r4**3
f(50) = r1**4*r3 - r1*r3**4 + r2**4*r4 - r2*r4**4
f(51) = r1**5*r3 - r1*r3**5 + r2**5*r4 - r2*r4**5
f(52) = r1**3*r3**2 - r1**2*r3**3 + r2**3*r4**2 - r2**2*r4**3
f(53) = -r1**4*r3**2 + r1**2*r3**4 - r2**4*r4**2 + r2**2*r4**4
f(54) = -r1**3*r3**2 + r1**2*r3**3 - r2**3*r4**2 + r2**2*r4**3
f(55) = r1**4*r3**2 - r1**2*r3**4 + r2**4*r4**2 - r2**2*r4**4
f(56) = r1**2*r4 - r1*r4**2 + r2**2*r3 - r2*r3**2
f(57) = -r1**3*r4 + r1*r4**3 - r2**3*r3 + r2*r3**3
f(58) = -r1**4*r4 + r1*r4**4 - r2**4*r3 + r2*r3**4
f(59) = r1**5*r4 - r1*r4**5 + r2**5*r3 - r2*r3**5
f(60) = -r1**3*r4**2 + r1**2*r4**3 - r2**3*r3**2 + r2**2*r3**3
f(61) = -r1**4*r4**2 + r1**2*r4**4 - r2**4*r3**2 + r2**2*r3**4
f(62) = r1**4*r4 - r1*r4**4 + r2**4*r3 - r2*r3**4
f(63) = a1*r1 + a2*r2 - a3*r3 - a4*r4
f(64) = -a1**2*r1 - a2**2*r2 + a3**2*r3 + a4**2*r4
f(65) = a1**3*r1 + a2**3*r2 - a3**3*r3 - a4**3*r4
f(66) = a1**4*r1 + a2**4*r2 - a3**4*r3 - a4**4*r4
f(67) = a1**5*r1 + a2**5*r2 - a3**5*r3 - a4**5*r4
f(68) = -a1*r1**2 - a2*r2**2 + a3*r3**2 + a4*r4**2
f(69) = -a1**2*r1**2 - a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2
f(70) = -a1**3*r1**2 - a2**3*r2**2 + a3**3*r3**2 + a4**3*r4**2
f(71) = -a1**4*r1**2 - a2**4*r2**2 + a3**4*r3**2 + a4**4*r4**2
f(72) = -a1*r1**3 - a2*r2**3 + a3*r3**3 + a4*r4**3
f(73) = a1**2*r1**3 + a2**2*r2**3 - a3**2*r3**3 - a4**2*r4**3
f(74) = a1**3*r1**3 + a2**3*r2**3 - a3**3*r3**3 - a4**3*r4**3
f(75) = -a1*r1**4 - a2*r2**4 + a3*r3**4 + a4*r4**4
f(76) = -a1**2*r1**4 - a2**2*r2**4 + a3**2*r3**4 + a4**2*r4**4
f(77) = -a1*r1**5 - a2*r2**5 + a3*r3**5 + a4*r4**5
f(78) = a1*r2 + a2*r1 - a3*r4 - a4*r3
f(79) = -a1**2*r2 - a2**2*r1 + a3**2*r4 + a4**2*r3
f(80) = a1**3*r2 + a2**3*r1 - a3**3*r4 - a4**3*r3
f(81) = -a1**4*r2 - a2**4*r1 + a3**4*r4 + a4**4*r3
f(82) = -a1**5*r2 - a2**5*r1 + a3**5*r4 + a4**5*r3
f(83) = -a1*r2**2 - a2*r1**2 + a3*r4**2 + a4*r3**2
f(84) = -a1**2*r2**2 - a2**2*r1**2 + a3**2*r4**2 + a4**2*r3**2
f(85) = a1**3*r2**2 + a2**3*r1**2 - a3**3*r4**2 - a4**3*r3**2
f(86) = -a1**4*r2**2 - a2**4*r1**2 + a3**4*r4**2 + a4**4*r3**2
f(87) = a1*r2**3 + a2*r1**3 - a3*r4**3 - a4*r3**3
f(88) = -a1**2*r2**3 - a2**2*r1**3 + a3**2*r4**3 + a4**2*r3**3
f(89) = a1**3*r2**3 + a2**3*r1**3 - a3**3*r4**3 - a4**3*r3**3
f(90) = a1*r2**4 + a2*r1**4 - a3*r4**4 - a4*r3**4
f(91) = -a1**2*r2**4 - a2**2*r1**4 + a3**2*r4**4 + a4**2*r3**4
f(92) = -a1*r2**5 - a2*r1**5 + a3*r4**5 + a4*r3**5
f(93) = -a1*r3 - a2*r4 + a3*r1 + a4*r2
f(94) = a1**2*r3 + a2**2*r4 - a3**2*r1 - a4**2*r2
f(95) = -a1**3*r3 - a2**3*r4 + a3**3*r1 + a4**3*r2
f(96) = a1**4*r3 + a2**4*r4 - a3**4*r1 - a4**4*r2
f(97) = a1**5*r3 + a2**5*r4 - a3**5*r1 - a4**5*r2
f(98) = -a1*r3**2 - a2*r4**2 + a3*r1**2 + a4*r2**2
f(99) = a1**2*r3**2 + a2**2*r4**2 - a3**2*r1**2 - a4**2*r2**2
f(100) = a1**3*r3**2 + a2**3*r4**2 - a3**3*r1**2 - a4**3*r2**2
f(101) = -a1**4*r3**2 - a2**4*r4**2 + a3**4*r1**2 + a4**4*r2**2
f(102) = a1*r3**3 + a2*r4**3 - a3*r1**3 - a4*r2**3
f(103) = -a1**2*r3**3 - a2**2*r4**3 + a3**2*r1**3 + a4**2*r2**3
f(104) = -a1**3*r3**3 - a2**3*r4**3 + a3**3*r1**3 + a4**3*r2**3
f(105) = a1*r3**4 + a2*r4**4 - a3*r1**4 - a4*r2**4
f(106) = -a1**2*r3**4 - a2**2*r4**4 + a3**2*r1**4 + a4**2*r2**4
f(107) = -a1*r3**5 - a2*r4**5 + a3*r1**5 + a4*r2**5
f(108) = a1*r4 + a2*r3 - a3*r2 - a4*r1
f(109) = -a1**2*r4 - a2**2*r3 + a3**2*r2 + a4**2*r1
f(110) = -a1**3*r4 - a2**3*r3 + a3**3*r2 + a4**3*r1
f(111) = -a1**4*r4 - a2**4*r3 + a3**4*r2 + a4**4*r1
f(112) = -a1**5*r4 - a2**5*r3 + a3**5*r2 + a4**5*r1
f(113) = a1*r4**2 + a2*r3**2 - a3*r2**2 - a4*r1**2
f(114) = -a1**2*r4**2 - a2**2*r3**2 + a3**2*r2**2 + a4**2*r1**2
f(115) = -a1**3*r4**2 - a2**3*r3**2 + a3**3*r2**2 + a4**3*r1**2
f(116) = -a1**4*r4**2 - a2**4*r3**2 + a3**4*r2**2 + a4**4*r1**2
f(117) = -a1*r4**3 - a2*r3**3 + a3*r2**3 + a4*r1**3
f(118) = -a1**2*r4**3 - a2**2*r3**3 + a3**2*r2**3 + a4**2*r1**3
f(119) = a1**3*r4**3 + a2**3*r3**3 - a3**3*r2**3 - a4**3*r1**3
f(120) = -a1*r4**4 - a2*r3**4 + a3*r2**4 + a4*r1**4
f(121) = -a1**2*r4**4 - a2**2*r3**4 + a3**2*r2**4 + a4**2*r1**4
f(122) = a1*r4**5 + a2*r3**5 - a3*r2**5 - a4*r1**5
f(123) = -b1**2*r1 - b1**2*r2 + b2**2*r3 + b2**2*r4
f(124) = -b1**4*r1 - b1**4*r2 + b2**4*r3 + b2**4*r4
f(125) = b1**2*r1**2 + b1**2*r2**2 - b2**2*r3**2 - b2**2*r4**2
f(126) = -b1**4*r1**2 - b1**4*r2**2 + b2**4*r3**2 + b2**4*r4**2
f(127) = b1**2*r1**3 + b1**2*r2**3 - b2**2*r3**3 - b2**2*r4**3
f(128) = b1**2*r1**4 + b1**2*r2**4 - b2**2*r3**4 - b2**2*r4**4
f(129) = b1**2*r3 + b1**2*r4 - b2**2*r1 - b2**2*r2
f(130) = b1**4*r3 + b1**4*r4 - b2**4*r1 - b2**4*r2
f(131) = b1**2*r3**2 + b1**2*r4**2 - b2**2*r1**2 - b2**2*r2**2
f(132) = b1**4*r3**2 + b1**4*r4**2 - b2**4*r1**2 - b2**4*r2**2
f(133) = b1**2*r3**3 + b1**2*r4**3 - b2**2*r1**3 - b2**2*r2**3
f(134) = -b1**2*r3**4 - b1**2*r4**4 + b2**2*r1**4 + b2**2*r2**4
f(135) = dtau**2*(r1 + r2 - r3 - r4)
f(136) = dtau**4*(r1 + r2 - r3 - r4)
f(137) = dtau**2*(-r1**2 - r2**2 + r3**2 + r4**2)
f(138) = dtau**4*(r1**2 + r2**2 - r3**2 - r4**2)
f(139) = dtau**2*(-r1**3 - r2**3 + r3**3 + r4**3)
f(140) = dtau**2*(r1**4 + r2**4 - r3**4 - r4**4)
f(141) = a1*a2 - a3*a4
f(142) = -a1**2*a2 - a1*a2**2 + a3**2*a4 + a3*a4**2
f(143) = -a1**3*a2 - a1*a2**3 + a3**3*a4 + a3*a4**3
f(144) = -a1**4*a2 - a1*a2**4 + a3**4*a4 + a3*a4**4
f(145) = a1**5*a2 + a1*a2**5 - a3**5*a4 - a3*a4**5
f(146) = -a1**2*a2**2 + a3**2*a4**2
f(147) = a1**3*a2**2 + a1**2*a2**3 - a3**3*a4**2 - a3**2*a4**3
f(148) = a1**4*a2**2 + a1**2*a2**4 - a3**4*a4**2 - a3**2*a4**4
f(149) = a1**3*a2**3 - a3**3*a4**3
f(150) = -a1**5*a2 - a1*a2**5 + a3**5*a4 + a3*a4**5
f(151) = -a1**2*a3 + a1*a3**2 - a2**2*a4 + a2*a4**2
f(152) = a1**3*a3 - a1*a3**3 + a2**3*a4 - a2*a4**3
f(153) = a1**4*a3 - a1*a3**4 + a2**4*a4 - a2*a4**4
f(154) = a1**5*a3 - a1*a3**5 + a2**5*a4 - a2*a4**5
f(155) = -a1**3*a3**2 + a1**2*a3**3 - a2**3*a4**2 + a2**2*a4**3
f(156) = -a1**4*a3**2 + a1**2*a3**4 - a2**4*a4**2 + a2**2*a4**4
f(157) = a1**3*a3**2 - a1**2*a3**3 + a2**3*a4**2 - a2**2*a4**3
f(158) = a1**4*a3**2 - a1**2*a3**4 + a2**4*a4**2 - a2**2*a4**4
f(159) = -a1**2*a4 + a1*a4**2 - a2**2*a3 + a2*a3**2
f(160) = a1**3*a4 - a1*a4**3 + a2**3*a3 - a2*a3**3
f(161) = a1**4*a4 - a1*a4**4 + a2**4*a3 - a2*a3**4
f(162) = -a1**5*a4 + a1*a4**5 - a2**5*a3 + a2*a3**5
f(163) = -a1**3*a4**2 + a1**2*a4**3 - a2**3*a3**2 + a2**2*a3**3
f(164) = a1**4*a4**2 - a1**2*a4**4 + a2**4*a3**2 - a2**2*a3**4
f(165) = -a1*b1**2 - a2*b1**2 + a3*b2**2 + a4*b2**2
f(166) = a1*b1**4 + a2*b1**4 - a3*b2**4 - a4*b2**4
f(167) = -a1**2*b1**2 - a2**2*b1**2 + a3**2*b2**2 + a4**2*b2**2
f(168) = -a1**2*b1**4 - a2**2*b1**4 + a3**2*b2**4 + a4**2*b2**4
f(169) = a1**3*b1**2 + a2**3*b1**2 - a3**3*b2**2 - a4**3*b2**2
f(170) = -a1**4*b1**2 - a2**4*b1**2 + a3**4*b2**2 + a4**4*b2**2
f(171) = -a1*b2**2 - a2*b2**2 + a3*b1**2 + a4*b1**2
f(172) = a1*b2**4 + a2*b2**4 - a3*b1**4 - a4*b1**4
f(173) = -a1**2*b2**2 - a2**2*b2**2 + a3**2*b1**2 + a4**2*b1**2
f(174) = -a1**2*b2**4 - a2**2*b2**4 + a3**2*b1**4 + a4**2*b1**4
f(175) = a1**3*b2**2 + a2**3*b2**2 - a3**3*b1**2 - a4**3*b1**2
f(176) = -a1**4*b2**2 - a2**4*b2**2 + a3**4*b1**2 + a4**4*b1**2
f(177) = dtau**2*(a1 + a2 - a3 - a4)
f(178) = dtau**4*(-a1 - a2 + a3 + a4)
f(179) = dtau**2*(a1**2 + a2**2 - a3**2 - a4**2)
f(180) = dtau**4*(-a1**2 - a2**2 + a3**2 + a4**2)
f(181) = dtau**2*(-a1**3 - a2**3 + a3**3 + a4**3)
f(182) = dtau**2*(a1**4 + a2**4 - a3**4 - a4**4)
f(183) = b1*b2*(b1**2 - b2**2)
f(184) = b1*b2*(-b1**4 + b2**4)
f(185) = b1**2*b2**2*(-b1**2 + b2**2)
f(186) = dtau**2*(-b1**2 + b2**2)
f(187) = dtau**4*(-b1**2 + b2**2)
f(188) = dtau**2*(-b1**4 + b2**4)
v = sum(f*params)
end function c2h4_dipole_b1u_n2_d6_ADF


!###############################################################################


function c2h4_dipole_b1u_n3_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(688)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(688)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(r1*r2 - r3*r4)
f(2) = r0*(r1**2*r2 + r1*r2**2 - r3**2*r4 - r3*r4**2)
f(3) = r0*(-r1**3*r2 - r1*r2**3 + r3**3*r4 + r3*r4**3)
f(4) = r0*(-r1**4*r2 - r1*r2**4 + r3**4*r4 + r3*r4**4)
f(5) = r0*(-r1**2*r2**2 + r3**2*r4**2)
f(6) = r0*(r1**3*r2**2 + r1**2*r2**3 - r3**3*r4**2 - r3**2*r4**3)
f(7) = r0*(r1**3*r2 + r1*r2**3 - r3**3*r4 - r3*r4**3)
f(8) = r0*(-r1**3*r2**2 - r1**2*r2**3 + r3**3*r4**2 + r3**2*r4**3)
f(9) = r0*(r1**4*r2 + r1*r2**4 - r3**4*r4 - r3*r4**4)
f(10) = r0**2*(r1*r2 - r3*r4)
f(11) = r0**2*(r1**2*r2 + r1*r2**2 - r3**2*r4 - r3*r4**2)
f(12) = r0**2*(-r1**3*r2 - r1*r2**3 + r3**3*r4 + r3*r4**3)
f(13) = r0**2*(-r1**2*r2 - r1*r2**2 + r3**2*r4 + r3*r4**2)
f(14) = r0**2*(r1**2*r2**2 - r3**2*r4**2)
f(15) = r0**3*(-r1*r2 + r3*r4)
f(16) = r0**3*(r1**2*r2 + r1*r2**2 - r3**2*r4 - r3*r4**2)
f(17) = r0**4*(-r1*r2 + r3*r4)
f(18) = r0*(r1**2*r3 - r1*r3**2 + r2**2*r4 - r2*r4**2)
f(19) = r0*(-r1**3*r3 + r1*r3**3 - r2**3*r4 + r2*r4**3)
f(20) = r0*(-r1**4*r3 + r1*r3**4 - r2**4*r4 + r2*r4**4)
f(21) = r0*(-r1**3*r3**2 + r1**2*r3**3 - r2**3*r4**2 + r2**2*r4**3)
f(22) = r0*(r1**3*r3 - r1*r3**3 + r2**3*r4 - r2*r4**3)
f(23) = r0*(r1**3*r3**2 - r1**2*r3**3 + r2**3*r4**2 - r2**2*r4**3)
f(24) = r0**2*(-r1**2*r3 + r1*r3**2 - r2**2*r4 + r2*r4**2)
f(25) = r0**2*(-r1**3*r3 + r1*r3**3 - r2**3*r4 + r2*r4**3)
f(26) = r0**2*(r1**3*r3 - r1*r3**3 + r2**3*r4 - r2*r4**3)
f(27) = r0**3*(r1**2*r3 - r1*r3**2 + r2**2*r4 - r2*r4**2)
f(28) = r0*(r1**2*r4 - r1*r4**2 + r2**2*r3 - r2*r3**2)
f(29) = r0*(-r1**3*r4 + r1*r4**3 - r2**3*r3 + r2*r3**3)
f(30) = r0*(-r1**4*r4 + r1*r4**4 - r2**4*r3 + r2*r3**4)
f(31) = r0*(-r1**3*r4**2 + r1**2*r4**3 - r2**3*r3**2 + r2**2*r3**3)
f(32) = r0**2*(-r1**2*r4 + r1*r4**2 - r2**2*r3 + r2*r3**2)
f(33) = r0**2*(-r1**3*r4 + r1*r4**3 - r2**3*r3 + r2*r3**3)
f(34) = r0**3*(r1**2*r4 - r1*r4**2 + r2**2*r3 - r2*r3**2)
f(35) = r0*(a1*r1 + a2*r2 - a3*r3 - a4*r4)
f(36) = r0*(-a1**2*r1 - a2**2*r2 + a3**2*r3 + a4**2*r4)
f(37) = r0*(a1**3*r1 + a2**3*r2 - a3**3*r3 - a4**3*r4)
f(38) = r0*(-a1**4*r1 - a2**4*r2 + a3**4*r3 + a4**4*r4)
f(39) = r0*(-a1*r1**2 - a2*r2**2 + a3*r3**2 + a4*r4**2)
f(40) = r0*(-a1**2*r1**2 - a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2)
f(41) = r0*(-a1**3*r1**2 - a2**3*r2**2 + a3**3*r3**2 + a4**3*r4**2)
f(42) = r0*(-a1*r1**3 - a2*r2**3 + a3*r3**3 + a4*r4**3)
f(43) = r0*(-a1**2*r1**3 - a2**2*r2**3 + a3**2*r3**3 + a4**2*r4**3)
f(44) = r0*(-a1*r1**4 - a2*r2**4 + a3*r3**4 + a4*r4**4)
f(45) = r0**2*(a1*r1 + a2*r2 - a3*r3 - a4*r4)
f(46) = r0**2*(a1**2*r1 + a2**2*r2 - a3**2*r3 - a4**2*r4)
f(47) = r0**2*(-a1**3*r1 - a2**3*r2 + a3**3*r3 + a4**3*r4)
f(48) = r0**2*(-a1*r1**2 - a2*r2**2 + a3*r3**2 + a4*r4**2)
f(49) = r0**2*(-a1**2*r1**2 - a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2)
f(50) = r0**2*(a1*r1**3 + a2*r2**3 - a3*r3**3 - a4*r4**3)
f(51) = r0**3*(-a1*r1 - a2*r2 + a3*r3 + a4*r4)
f(52) = r0**3*(a1**2*r1 + a2**2*r2 - a3**2*r3 - a4**2*r4)
f(53) = r0**3*(a1*r1**2 + a2*r2**2 - a3*r3**2 - a4*r4**2)
f(54) = r0**4*(-a1*r1 - a2*r2 + a3*r3 + a4*r4)
f(55) = r0*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(56) = r0*(-a1**2*r2 - a2**2*r1 + a3**2*r4 + a4**2*r3)
f(57) = r0*(-a1**3*r2 - a2**3*r1 + a3**3*r4 + a4**3*r3)
f(58) = r0*(a1**4*r2 + a2**4*r1 - a3**4*r4 - a4**4*r3)
f(59) = r0*(a1*r2**2 + a2*r1**2 - a3*r4**2 - a4*r3**2)
f(60) = r0*(-a1**2*r2**2 - a2**2*r1**2 + a3**2*r4**2 + a4**2*r3**2)
f(61) = r0*(-a1**3*r2**2 - a2**3*r1**2 + a3**3*r4**2 + a4**3*r3**2)
f(62) = r0*(a1*r2**3 + a2*r1**3 - a3*r4**3 - a4*r3**3)
f(63) = r0*(a1**2*r2**3 + a2**2*r1**3 - a3**2*r4**3 - a4**2*r3**3)
f(64) = r0*(-a1*r2**4 - a2*r1**4 + a3*r4**4 + a4*r3**4)
f(65) = r0**2*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(66) = r0**2*(-a1**2*r2 - a2**2*r1 + a3**2*r4 + a4**2*r3)
f(67) = r0**2*(-a1**3*r2 - a2**3*r1 + a3**3*r4 + a4**3*r3)
f(68) = r0**2*(-a1*r2**2 - a2*r1**2 + a3*r4**2 + a4*r3**2)
f(69) = r0**2*(a1**2*r2**2 + a2**2*r1**2 - a3**2*r4**2 - a4**2*r3**2)
f(70) = r0**2*(-a1*r2**3 - a2*r1**3 + a3*r4**3 + a4*r3**3)
f(71) = r0**3*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(72) = r0**3*(-a1**2*r2 - a2**2*r1 + a3**2*r4 + a4**2*r3)
f(73) = r0**3*(a1*r2**2 + a2*r1**2 - a3*r4**2 - a4*r3**2)
f(74) = r0**4*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(75) = r0*(-a1*r3 - a2*r4 + a3*r1 + a4*r2)
f(76) = r0*(-a1**2*r3 - a2**2*r4 + a3**2*r1 + a4**2*r2)
f(77) = r0*(-a1**3*r3 - a2**3*r4 + a3**3*r1 + a4**3*r2)
f(78) = r0*(a1**4*r3 + a2**4*r4 - a3**4*r1 - a4**4*r2)
f(79) = r0*(a1*r3**2 + a2*r4**2 - a3*r1**2 - a4*r2**2)
f(80) = r0*(a1**2*r3**2 + a2**2*r4**2 - a3**2*r1**2 - a4**2*r2**2)
f(81) = r0*(a1**3*r3**2 + a2**3*r4**2 - a3**3*r1**2 - a4**3*r2**2)
f(82) = r0*(a1*r3**3 + a2*r4**3 - a3*r1**3 - a4*r2**3)
f(83) = r0*(-a1**2*r3**3 - a2**2*r4**3 + a3**2*r1**3 + a4**2*r2**3)
f(84) = r0*(-a1*r3**4 - a2*r4**4 + a3*r1**4 + a4*r2**4)
f(85) = r0**2*(a1*r3 + a2*r4 - a3*r1 - a4*r2)
f(86) = r0**2*(a1**2*r3 + a2**2*r4 - a3**2*r1 - a4**2*r2)
f(87) = r0**2*(a1**3*r3 + a2**3*r4 - a3**3*r1 - a4**3*r2)
f(88) = r0**2*(-a1*r3**2 - a2*r4**2 + a3*r1**2 + a4*r2**2)
f(89) = r0**2*(a1**2*r3**2 + a2**2*r4**2 - a3**2*r1**2 - a4**2*r2**2)
f(90) = r0**2*(a1*r3**3 + a2*r4**3 - a3*r1**3 - a4*r2**3)
f(91) = r0**3*(-a1*r3 - a2*r4 + a3*r1 + a4*r2)
f(92) = r0**3*(-a1**2*r3 - a2**2*r4 + a3**2*r1 + a4**2*r2)
f(93) = r0**3*(a1*r3**2 + a2*r4**2 - a3*r1**2 - a4*r2**2)
f(94) = r0**4*(a1*r3 + a2*r4 - a3*r1 - a4*r2)
f(95) = r0*(-a1*r4 - a2*r3 + a3*r2 + a4*r1)
f(96) = r0*(a1**2*r4 + a2**2*r3 - a3**2*r2 - a4**2*r1)
f(97) = r0*(-a1**3*r4 - a2**3*r3 + a3**3*r2 + a4**3*r1)
f(98) = r0*(-a1**4*r4 - a2**4*r3 + a3**4*r2 + a4**4*r1)
f(99) = r0*(-a1*r4**2 - a2*r3**2 + a3*r2**2 + a4*r1**2)
f(100) = r0*(a1**2*r4**2 + a2**2*r3**2 - a3**2*r2**2 - a4**2*r1**2)
f(101) = r0*(-a1**3*r4**2 - a2**3*r3**2 + a3**3*r2**2 + a4**3*r1**2)
f(102) = r0*(-a1*r4**3 - a2*r3**3 + a3*r2**3 + a4*r1**3)
f(103) = r0*(a1**2*r4**3 + a2**2*r3**3 - a3**2*r2**3 - a4**2*r1**3)
f(104) = r0*(a1*r4**4 + a2*r3**4 - a3*r2**4 - a4*r1**4)
f(105) = r0**2*(-a1*r4 - a2*r3 + a3*r2 + a4*r1)
f(106) = r0**2*(a1**2*r4 + a2**2*r3 - a3**2*r2 - a4**2*r1)
f(107) = r0**2*(a1**3*r4 + a2**3*r3 - a3**3*r2 - a4**3*r1)
f(108) = r0**2*(a1*r4**2 + a2*r3**2 - a3*r2**2 - a4*r1**2)
f(109) = r0**2*(a1**2*r4**2 + a2**2*r3**2 - a3**2*r2**2 - a4**2*r1**2)
f(110) = r0**2*(-a1*r4**3 - a2*r3**3 + a3*r2**3 + a4*r1**3)
f(111) = r0**3*(a1*r4 + a2*r3 - a3*r2 - a4*r1)
f(112) = r0**3*(a1**2*r4 + a2**2*r3 - a3**2*r2 - a4**2*r1)
f(113) = r0**3*(-a1*r4**2 - a2*r3**2 + a3*r2**2 + a4*r1**2)
f(114) = r0**4*(-a1*r4 - a2*r3 + a3*r2 + a4*r1)
f(115) = r0*(-b1**2*r1 - b1**2*r2 + b2**2*r3 + b2**2*r4)
f(116) = r0*(b1**4*r1 + b1**4*r2 - b2**4*r3 - b2**4*r4)
f(117) = r0*(-b1**2*r1**2 - b1**2*r2**2 + b2**2*r3**2 + b2**2*r4**2)
f(118) = r0*(-b1**2*r1**3 - b1**2*r2**3 + b2**2*r3**3 + b2**2*r4**3)
f(119) = r0**2*(b1**2*r1 + b1**2*r2 - b2**2*r3 - b2**2*r4)
f(120) = r0**2*(b1**2*r1**2 + b1**2*r2**2 - b2**2*r3**2 - b2**2*r4**2)
f(121) = r0**3*(-b1**2*r1 - b1**2*r2 + b2**2*r3 + b2**2*r4)
f(122) = r0*(b1**2*r3 + b1**2*r4 - b2**2*r1 - b2**2*r2)
f(123) = r0*(-b1**4*r3 - b1**4*r4 + b2**4*r1 + b2**4*r2)
f(124) = r0*(-b1**2*r3**2 - b1**2*r4**2 + b2**2*r1**2 + b2**2*r2**2)
f(125) = r0*(-b1**2*r3**3 - b1**2*r4**3 + b2**2*r1**3 + b2**2*r2**3)
f(126) = r0**2*(b1**2*r3 + b1**2*r4 - b2**2*r1 - b2**2*r2)
f(127) = r0**2*(b1**2*r3**2 + b1**2*r4**2 - b2**2*r1**2 - b2**2*r2**2)
f(128) = r0**3*(-b1**2*r3 - b1**2*r4 + b2**2*r1 + b2**2*r2)
f(129) = dtau**2*r0*(r1 + r2 - r3 - r4)
f(130) = dtau**4*r0*(-r1 - r2 + r3 + r4)
f(131) = dtau**2*r0*(-r1**2 - r2**2 + r3**2 + r4**2)
f(132) = dtau**2*r0*(-r1**3 - r2**3 + r3**3 + r4**3)
f(133) = dtau**2*r0**2*(-r1 - r2 + r3 + r4)
f(134) = dtau**2*r0**2*(-r1**2 - r2**2 + r3**2 + r4**2)
f(135) = dtau**2*r0**3*(r1 + r2 - r3 - r4)
f(136) = r0*(a1*a2 - a3*a4)
f(137) = r0*(-a1**2*a2 - a1*a2**2 + a3**2*a4 + a3*a4**2)
f(138) = r0*(a1**3*a2 + a1*a2**3 - a3**3*a4 - a3*a4**3)
f(139) = r0*(-a1**4*a2 - a1*a2**4 + a3**4*a4 + a3*a4**4)
f(140) = r0*(a1**2*a2 + a1*a2**2 - a3**2*a4 - a3*a4**2)
f(141) = r0*(a1**2*a2**2 - a3**2*a4**2)
f(142) = r0*(-a1**3*a2**2 - a1**2*a2**3 + a3**3*a4**2 + a3**2*a4**3)
f(143) = r0*(a1**3*a2**2 + a1**2*a2**3 - a3**3*a4**2 - a3**2*a4**3)
f(144) = r0*(a1**4*a2 + a1*a2**4 - a3**4*a4 - a3*a4**4)
f(145) = r0**2*(a1*a2 - a3*a4)
f(146) = r0**2*(a1**2*a2 + a1*a2**2 - a3**2*a4 - a3*a4**2)
f(147) = r0**2*(-a1**3*a2 - a1*a2**3 + a3**3*a4 + a3*a4**3)
f(148) = r0**2*(a1**2*a2**2 - a3**2*a4**2)
f(149) = r0**3*(-a1*a2 + a3*a4)
f(150) = r0**3*(-a1**2*a2 - a1*a2**2 + a3**2*a4 + a3*a4**2)
f(151) = r0**4*(-a1*a2 + a3*a4)
f(152) = r0*(-a1**2*a3 + a1*a3**2 - a2**2*a4 + a2*a4**2)
f(153) = r0*(a1**3*a3 - a1*a3**3 + a2**3*a4 - a2*a4**3)
f(154) = r0*(a1**4*a3 - a1*a3**4 + a2**4*a4 - a2*a4**4)
f(155) = r0*(-a1**3*a3**2 + a1**2*a3**3 - a2**3*a4**2 + a2**2*a4**3)
f(156) = r0*(-a1**3*a3 + a1*a3**3 - a2**3*a4 + a2*a4**3)
f(157) = r0**2*(a1**2*a3 - a1*a3**2 + a2**2*a4 - a2*a4**2)
f(158) = r0**2*(-a1**3*a3 + a1*a3**3 - a2**3*a4 + a2*a4**3)
f(159) = r0**3*(a1**2*a3 - a1*a3**2 + a2**2*a4 - a2*a4**2)
f(160) = r0*(-a1**2*a4 + a1*a4**2 - a2**2*a3 + a2*a3**2)
f(161) = r0*(a1**3*a4 - a1*a4**3 + a2**3*a3 - a2*a3**3)
f(162) = r0*(a1**4*a4 - a1*a4**4 + a2**4*a3 - a2*a3**4)
f(163) = r0*(a1**2*a4 - a1*a4**2 + a2**2*a3 - a2*a3**2)
f(164) = r0*(-a1**3*a4**2 + a1**2*a4**3 - a2**3*a3**2 + a2**2*a3**3)
f(165) = r0*(-a1**4*a4 + a1*a4**4 - a2**4*a3 + a2*a3**4)
f(166) = r0**2*(-a1**2*a4 + a1*a4**2 - a2**2*a3 + a2*a3**2)
f(167) = r0**2*(-a1**3*a4 + a1*a4**3 - a2**3*a3 + a2*a3**3)
f(168) = r0**3*(a1**2*a4 - a1*a4**2 + a2**2*a3 - a2*a3**2)
f(169) = r0**3*(-a1**2*a4 + a1*a4**2 - a2**2*a3 + a2*a3**2)
f(170) = r0*(a1*b1**2 + a2*b1**2 - a3*b2**2 - a4*b2**2)
f(171) = r0*(-a1*b1**4 - a2*b1**4 + a3*b2**4 + a4*b2**4)
f(172) = r0*(a1**2*b1**2 + a2**2*b1**2 - a3**2*b2**2 - a4**2*b2**2)
f(173) = r0*(-a1**3*b1**2 - a2**3*b1**2 + a3**3*b2**2 + a4**3*b2**2)
f(174) = r0**2*(-a1*b1**2 - a2*b1**2 + a3*b2**2 + a4*b2**2)
f(175) = r0**2*(-a1**2*b1**2 - a2**2*b1**2 + a3**2*b2**2 + a4**2*b2**2)
f(176) = r0**3*(a1*b1**2 + a2*b1**2 - a3*b2**2 - a4*b2**2)
f(177) = r0*(a1*b2**2 + a2*b2**2 - a3*b1**2 - a4*b1**2)
f(178) = r0*(a1*b2**4 + a2*b2**4 - a3*b1**4 - a4*b1**4)
f(179) = r0*(a1**2*b2**2 + a2**2*b2**2 - a3**2*b1**2 - a4**2*b1**2)
f(180) = r0*(a1**3*b2**2 + a2**3*b2**2 - a3**3*b1**2 - a4**3*b1**2)
f(181) = r0**2*(a1*b2**2 + a2*b2**2 - a3*b1**2 - a4*b1**2)
f(182) = r0**2*(a1**2*b2**2 + a2**2*b2**2 - a3**2*b1**2 - a4**2*b1**2)
f(183) = r0**3*(-a1*b2**2 - a2*b2**2 + a3*b1**2 + a4*b1**2)
f(184) = dtau**2*r0*(a1 + a2 - a3 - a4)
f(185) = dtau**4*r0*(a1 + a2 - a3 - a4)
f(186) = dtau**2*r0*(a1**2 + a2**2 - a3**2 - a4**2)
f(187) = dtau**2*r0*(a1**3 + a2**3 - a3**3 - a4**3)
f(188) = dtau**2*r0**2*(a1 + a2 - a3 - a4)
f(189) = dtau**2*r0**2*(a1**2 + a2**2 - a3**2 - a4**2)
f(190) = dtau**2*r0**3*(a1 + a2 - a3 - a4)
f(191) = b1*b2*r0*(b1**2 - b2**2)
f(192) = b1*b2*r0**2*(-b1**2 + b2**2)
f(193) = dtau**2*r0*(b1**2 - b2**2)
f(194) = dtau**2*r0**2*(-b1**2 + b2**2)
f(195) = r1*r2*r3 + r1*r2*r4 - r1*r3*r4 - r2*r3*r4
f(196) = -r1**2*r3*r4 + r1*r2*r3**2 + r1*r2*r4**2 - r2**2*r3*r4
f(197) = -r1**3*r3*r4 + r1*r2*r3**3 + r1*r2*r4**3 - r2**3*r3*r4
f(198) = -r1**4*r3*r4 + r1*r2*r3**4 + r1*r2*r4**4 - r2**4*r3*r4
f(199) = r1**2*r2*r4 + r1*r2**2*r3 - r1*r3*r4**2 - r2*r3**2*r4
f(200) = r1**2*r2*r4**2 - r1**2*r3*r4**2 + r1*r2**2*r3**2 - r2**2*r3**2* &
      r4
f(201) = r1**3*r3*r4**2 - r1**2*r2*r4**3 - r1*r2**2*r3**3 + r2**3*r3**2* &
      r4
f(202) = -r1**3*r2*r4 - r1*r2**3*r3 + r1*r3*r4**3 + r2*r3**3*r4
f(203) = -r1**3*r2*r4**2 + r1**2*r3*r4**3 - r1*r2**3*r3**2 + r2**2*r3**3 &
      *r4
f(204) = -r1**4*r2*r4 - r1*r2**4*r3 + r1*r3*r4**4 + r2*r3**4*r4
f(205) = -r1**2*r2*r3 - r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4**2
f(206) = -r1**2*r2*r3**2 + r1**2*r3**2*r4 - r1*r2**2*r4**2 + r2**2*r3*r4 &
      **2
f(207) = -r1**3*r3**2*r4 + r1**2*r2*r3**3 + r1*r2**2*r4**3 - r2**3*r3*r4 &
      **2
f(208) = -r1**2*r2**2*r3 - r1**2*r2**2*r4 + r1*r3**2*r4**2 + r2*r3**2*r4 &
      **2
f(209) = r1**2*r2**2*r3**2 + r1**2*r2**2*r4**2 - r1**2*r3**2*r4**2 - r2 &
      **2*r3**2*r4**2
f(210) = r1**3*r2**2*r4 + r1**2*r2**3*r3 - r1*r3**2*r4**3 - r2*r3**3*r4 &
      **2
f(211) = r1**3*r2*r3 + r1*r2**3*r4 - r1*r3**3*r4 - r2*r3*r4**3
f(212) = r1**3*r2*r3**2 - r1**2*r3**3*r4 + r1*r2**3*r4**2 - r2**2*r3*r4 &
      **3
f(213) = -r1**3*r2**2*r3 - r1**2*r2**3*r4 + r1*r3**3*r4**2 + r2*r3**2*r4 &
      **3
f(214) = -r1**4*r2*r3 - r1*r2**4*r4 + r1*r3**4*r4 + r2*r3*r4**4
f(215) = a1*r1*r2 + a2*r1*r2 - a3*r3*r4 - a4*r3*r4
f(216) = -a1**2*r1*r2 - a2**2*r1*r2 + a3**2*r3*r4 + a4**2*r3*r4
f(217) = a1**3*r1*r2 + a2**3*r1*r2 - a3**3*r3*r4 - a4**3*r3*r4
f(218) = a1**4*r1*r2 + a2**4*r1*r2 - a3**4*r3*r4 - a4**4*r3*r4
f(219) = a1*r1*r2**2 + a2*r1**2*r2 - a3*r3*r4**2 - a4*r3**2*r4
f(220) = -a1**2*r1*r2**2 - a2**2*r1**2*r2 + a3**2*r3*r4**2 + a4**2*r3**2 &
      *r4
f(221) = -a1**3*r1*r2**2 - a2**3*r1**2*r2 + a3**3*r3*r4**2 + a4**3*r3**2 &
      *r4
f(222) = a1*r1*r2**3 + a2*r1**3*r2 - a3*r3*r4**3 - a4*r3**3*r4
f(223) = -a1**2*r1*r2**3 - a2**2*r1**3*r2 + a3**2*r3*r4**3 + a4**2*r3**3 &
      *r4
f(224) = -a1*r1*r2**4 - a2*r1**4*r2 + a3*r3*r4**4 + a4*r3**4*r4
f(225) = -a1*r1**2*r2 - a2*r1*r2**2 + a3*r3**2*r4 + a4*r3*r4**2
f(226) = a1**2*r1**2*r2 + a2**2*r1*r2**2 - a3**2*r3**2*r4 - a4**2*r3*r4 &
      **2
f(227) = a1**3*r1**2*r2 + a2**3*r1*r2**2 - a3**3*r3**2*r4 - a4**3*r3*r4 &
      **2
f(228) = -a1*r1**2*r2**2 - a2*r1**2*r2**2 + a3*r3**2*r4**2 + a4*r3**2*r4 &
      **2
f(229) = -a1**2*r1**2*r2**2 - a2**2*r1**2*r2**2 + a3**2*r3**2*r4**2 + a4 &
      **2*r3**2*r4**2
f(230) = a1*r1**2*r2**3 + a2*r1**3*r2**2 - a3*r3**2*r4**3 - a4*r3**3*r4 &
      **2
f(231) = -a1*r1**3*r2 - a2*r1*r2**3 + a3*r3**3*r4 + a4*r3*r4**3
f(232) = a1**2*r1**3*r2 + a2**2*r1*r2**3 - a3**2*r3**3*r4 - a4**2*r3*r4 &
      **3
f(233) = -a1*r1**3*r2**2 - a2*r1**2*r2**3 + a3*r3**3*r4**2 + a4*r3**2*r4 &
      **3
f(234) = -a1*r1**4*r2 - a2*r1*r2**4 + a3*r3**4*r4 + a4*r3*r4**4
f(235) = -a1*r3*r4 - a2*r3*r4 + a3*r1*r2 + a4*r1*r2
f(236) = a1**2*r3*r4 + a2**2*r3*r4 - a3**2*r1*r2 - a4**2*r1*r2
f(237) = -a1**3*r3*r4 - a2**3*r3*r4 + a3**3*r1*r2 + a4**3*r1*r2
f(238) = -a1**4*r3*r4 - a2**4*r3*r4 + a3**4*r1*r2 + a4**4*r1*r2
f(239) = a1*r3*r4**2 + a2*r3**2*r4 - a3*r1*r2**2 - a4*r1**2*r2
f(240) = -a1**2*r3*r4**2 - a2**2*r3**2*r4 + a3**2*r1*r2**2 + a4**2*r1**2 &
      *r2
f(241) = -a1**3*r3*r4**2 - a2**3*r3**2*r4 + a3**3*r1*r2**2 + a4**3*r1**2 &
      *r2
f(242) = -a1*r3*r4**3 - a2*r3**3*r4 + a3*r1*r2**3 + a4*r1**3*r2
f(243) = a1**2*r3*r4**3 + a2**2*r3**3*r4 - a3**2*r1*r2**3 - a4**2*r1**3* &
      r2
f(244) = -a1*r3*r4**4 - a2*r3**4*r4 + a3*r1*r2**4 + a4*r1**4*r2
f(245) = -a1*r3**2*r4 - a2*r3*r4**2 + a3*r1**2*r2 + a4*r1*r2**2
f(246) = a1**2*r3**2*r4 + a2**2*r3*r4**2 - a3**2*r1**2*r2 - a4**2*r1*r2 &
      **2
f(247) = a1**3*r3**2*r4 + a2**3*r3*r4**2 - a3**3*r1**2*r2 - a4**3*r1*r2 &
      **2
f(248) = a1*r3**2*r4**2 + a2*r3**2*r4**2 - a3*r1**2*r2**2 - a4*r1**2*r2 &
      **2
f(249) = a1**2*r3**2*r4**2 + a2**2*r3**2*r4**2 - a3**2*r1**2*r2**2 - a4 &
      **2*r1**2*r2**2
f(250) = a1*r3**2*r4**3 + a2*r3**3*r4**2 - a3*r1**2*r2**3 - a4*r1**3*r2 &
      **2
f(251) = a1*r3**3*r4 + a2*r3*r4**3 - a3*r1**3*r2 - a4*r1*r2**3
f(252) = -a1**2*r3**3*r4 - a2**2*r3*r4**3 + a3**2*r1**3*r2 + a4**2*r1*r2 &
      **3
f(253) = -a1*r3**3*r4**2 - a2*r3**2*r4**3 + a3*r1**3*r2**2 + a4*r1**2*r2 &
      **3
f(254) = -a1*r3**4*r4 - a2*r3*r4**4 + a3*r1**4*r2 + a4*r1*r2**4
f(255) = b1**2*r1*r2 - b2**2*r3*r4
f(256) = b1**4*r1*r2 - b2**4*r3*r4
f(257) = b1**2*r1**2*r2 + b1**2*r1*r2**2 - b2**2*r3**2*r4 - b2**2*r3*r4 &
      **2
f(258) = -b1**2*r1**3*r2 - b1**2*r1*r2**3 + b2**2*r3**3*r4 + b2**2*r3*r4 &
      **3
f(259) = -b1**2*r1**2*r2**2 + b2**2*r3**2*r4**2
f(260) = -b1**2*r3*r4 + b2**2*r1*r2
f(261) = b1**4*r3*r4 - b2**4*r1*r2
f(262) = b1**2*r3**2*r4 + b1**2*r3*r4**2 - b2**2*r1**2*r2 - b2**2*r1*r2 &
      **2
f(263) = -b1**2*r3**3*r4 - b1**2*r3*r4**3 + b2**2*r1**3*r2 + b2**2*r1*r2 &
      **3
f(264) = -b1**2*r3**2*r4**2 + b2**2*r1**2*r2**2
f(265) = dtau**2*(-r1*r2 + r3*r4)
f(266) = dtau**4*(r1*r2 - r3*r4)
f(267) = dtau**2*(-r1**2*r2 - r1*r2**2 + r3**2*r4 + r3*r4**2)
f(268) = dtau**2*(r1**3*r2 + r1*r2**3 - r3**3*r4 - r3*r4**3)
f(269) = dtau**2*(-r1**2*r2**2 + r3**2*r4**2)
f(270) = a1*r1*r3 + a2*r2*r4 - a3*r1*r3 - a4*r2*r4
f(271) = a1**2*r1*r3 + a2**2*r2*r4 - a3**2*r1*r3 - a4**2*r2*r4
f(272) = a1**3*r1*r3 + a2**3*r2*r4 - a3**3*r1*r3 - a4**3*r2*r4
f(273) = a1**4*r1*r3 + a2**4*r2*r4 - a3**4*r1*r3 - a4**4*r2*r4
f(274) = -a1*r1*r3**2 - a2*r2*r4**2 + a3*r1**2*r3 + a4*r2**2*r4
f(275) = -a1**2*r1*r3**2 - a2**2*r2*r4**2 + a3**2*r1**2*r3 + a4**2*r2**2 &
      *r4
f(276) = a1**3*r1*r3**2 + a2**3*r2*r4**2 - a3**3*r1**2*r3 - a4**3*r2**2* &
      r4
f(277) = -a1*r1*r3**3 - a2*r2*r4**3 + a3*r1**3*r3 + a4*r2**3*r4
f(278) = -a1**2*r1*r3**3 - a2**2*r2*r4**3 + a3**2*r1**3*r3 + a4**2*r2**3 &
      *r4
f(279) = a1*r1*r3**4 + a2*r2*r4**4 - a3*r1**4*r3 - a4*r2**4*r4
f(280) = a1*r1**2*r3 + a2*r2**2*r4 - a3*r1*r3**2 - a4*r2*r4**2
f(281) = a1**2*r1**2*r3 + a2**2*r2**2*r4 - a3**2*r1*r3**2 - a4**2*r2*r4 &
      **2
f(282) = -a1**3*r1**2*r3 - a2**3*r2**2*r4 + a3**3*r1*r3**2 + a4**3*r2*r4 &
      **2
f(283) = a1*r1**2*r3**2 + a2*r2**2*r4**2 - a3*r1**2*r3**2 - a4*r2**2*r4 &
      **2
f(284) = a1**2*r1**2*r3**2 + a2**2*r2**2*r4**2 - a3**2*r1**2*r3**2 - a4 &
      **2*r2**2*r4**2
f(285) = a1*r1**2*r3**3 + a2*r2**2*r4**3 - a3*r1**3*r3**2 - a4*r2**3*r4 &
      **2
f(286) = -a1*r1**3*r3 - a2*r2**3*r4 + a3*r1*r3**3 + a4*r2*r4**3
f(287) = -a1**2*r1**3*r3 - a2**2*r2**3*r4 + a3**2*r1*r3**3 + a4**2*r2*r4 &
      **3
f(288) = a1*r1**3*r3**2 + a2*r2**3*r4**2 - a3*r1**2*r3**3 - a4*r2**2*r4 &
      **3
f(289) = -a1*r1**4*r3 - a2*r2**4*r4 + a3*r1*r3**4 + a4*r2*r4**4
f(290) = a1*r2*r4 + a2*r1*r3 - a3*r2*r4 - a4*r1*r3
f(291) = -a1**2*r2*r4 - a2**2*r1*r3 + a3**2*r2*r4 + a4**2*r1*r3
f(292) = -a1**3*r2*r4 - a2**3*r1*r3 + a3**3*r2*r4 + a4**3*r1*r3
f(293) = -a1**4*r2*r4 - a2**4*r1*r3 + a3**4*r2*r4 + a4**4*r1*r3
f(294) = -a1*r2*r4**2 - a2*r1*r3**2 + a3*r2**2*r4 + a4*r1**2*r3
f(295) = -a1**2*r2*r4**2 - a2**2*r1*r3**2 + a3**2*r2**2*r4 + a4**2*r1**2 &
      *r3
f(296) = a1**3*r2*r4**2 + a2**3*r1*r3**2 - a3**3*r2**2*r4 - a4**3*r1**2* &
      r3
f(297) = -a1*r2*r4**3 - a2*r1*r3**3 + a3*r2**3*r4 + a4*r1**3*r3
f(298) = -a1**2*r2*r4**3 - a2**2*r1*r3**3 + a3**2*r2**3*r4 + a4**2*r1**3 &
      *r3
f(299) = -a1*r2*r4**4 - a2*r1*r3**4 + a3*r2**4*r4 + a4*r1**4*r3
f(300) = a1*r2**2*r4 + a2*r1**2*r3 - a3*r2*r4**2 - a4*r1*r3**2
f(301) = a1**2*r2**2*r4 + a2**2*r1**2*r3 - a3**2*r2*r4**2 - a4**2*r1*r3 &
      **2
f(302) = -a1**3*r2**2*r4 - a2**3*r1**2*r3 + a3**3*r2*r4**2 + a4**3*r1*r3 &
      **2
f(303) = a1*r2**2*r4**2 + a2*r1**2*r3**2 - a3*r2**2*r4**2 - a4*r1**2*r3 &
      **2
f(304) = a1**2*r2**2*r4**2 + a2**2*r1**2*r3**2 - a3**2*r2**2*r4**2 - a4 &
      **2*r1**2*r3**2
f(305) = a1*r2**2*r4**3 + a2*r1**2*r3**3 - a3*r2**3*r4**2 - a4*r1**3*r3 &
      **2
f(306) = -a1*r2**3*r4 - a2*r1**3*r3 + a3*r2*r4**3 + a4*r1*r3**3
f(307) = -a1**2*r2**3*r4 - a2**2*r1**3*r3 + a3**2*r2*r4**3 + a4**2*r1*r3 &
      **3
f(308) = a1*r2**3*r4**2 + a2*r1**3*r3**2 - a3*r2**2*r4**3 - a4*r1**2*r3 &
      **3
f(309) = a1*r2**4*r4 + a2*r1**4*r3 - a3*r2*r4**4 - a4*r1*r3**4
f(310) = -b1**2*r1*r3 - b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2*r4
f(311) = b1**4*r1*r3 + b1**4*r2*r4 - b2**4*r1*r3 - b2**4*r2*r4
f(312) = -b1**2*r1*r3**2 - b1**2*r2*r4**2 + b2**2*r1**2*r3 + b2**2*r2**2 &
      *r4
f(313) = b1**2*r1*r3**3 + b1**2*r2*r4**3 - b2**2*r1**3*r3 - b2**2*r2**3* &
      r4
f(314) = b1**2*r1**2*r3 + b1**2*r2**2*r4 - b2**2*r1*r3**2 - b2**2*r2*r4 &
      **2
f(315) = -b1**2*r1**2*r3**2 - b1**2*r2**2*r4**2 + b2**2*r1**2*r3**2 + b2 &
      **2*r2**2*r4**2
f(316) = -b1**2*r1**3*r3 - b1**2*r2**3*r4 + b2**2*r1*r3**3 + b2**2*r2*r4 &
      **3
f(317) = dtau**2*(r1**2*r3 - r1*r3**2 + r2**2*r4 - r2*r4**2)
f(318) = dtau**2*(-r1**3*r3 + r1*r3**3 - r2**3*r4 + r2*r4**3)
f(319) = a1*r1*r4 + a2*r2*r3 - a3*r2*r3 - a4*r1*r4
f(320) = -a1**2*r1*r4 - a2**2*r2*r3 + a3**2*r2*r3 + a4**2*r1*r4
f(321) = -a1**3*r1*r4 - a2**3*r2*r3 + a3**3*r2*r3 + a4**3*r1*r4
f(322) = a1**4*r1*r4 + a2**4*r2*r3 - a3**4*r2*r3 - a4**4*r1*r4
f(323) = -a1*r1*r4**2 - a2*r2*r3**2 + a3*r2**2*r3 + a4*r1**2*r4
f(324) = -a1**2*r1*r4**2 - a2**2*r2*r3**2 + a3**2*r2**2*r3 + a4**2*r1**2 &
      *r4
f(325) = a1**3*r1*r4**2 + a2**3*r2*r3**2 - a3**3*r2**2*r3 - a4**3*r1**2* &
      r4
f(326) = -a1*r1*r4**3 - a2*r2*r3**3 + a3*r2**3*r3 + a4*r1**3*r4
f(327) = -a1**2*r1*r4**3 - a2**2*r2*r3**3 + a3**2*r2**3*r3 + a4**2*r1**3 &
      *r4
f(328) = a1*r1*r4**4 + a2*r2*r3**4 - a3*r2**4*r3 - a4*r1**4*r4
f(329) = -a1*r1**2*r4 - a2*r2**2*r3 + a3*r2*r3**2 + a4*r1*r4**2
f(330) = -a1**2*r1**2*r4 - a2**2*r2**2*r3 + a3**2*r2*r3**2 + a4**2*r1*r4 &
      **2
f(331) = a1**3*r1**2*r4 + a2**3*r2**2*r3 - a3**3*r2*r3**2 - a4**3*r1*r4 &
      **2
f(332) = a1*r1**2*r4**2 + a2*r2**2*r3**2 - a3*r2**2*r3**2 - a4*r1**2*r4 &
      **2
f(333) = a1**2*r1**2*r4**2 + a2**2*r2**2*r3**2 - a3**2*r2**2*r3**2 - a4 &
      **2*r1**2*r4**2
f(334) = a1*r1**2*r4**3 + a2*r2**2*r3**3 - a3*r2**3*r3**2 - a4*r1**3*r4 &
      **2
f(335) = a1*r1**3*r4 + a2*r2**3*r3 - a3*r2*r3**3 - a4*r1*r4**3
f(336) = a1**2*r1**3*r4 + a2**2*r2**3*r3 - a3**2*r2*r3**3 - a4**2*r1*r4 &
      **3
f(337) = a1*r1**3*r4**2 + a2*r2**3*r3**2 - a3*r2**2*r3**3 - a4*r1**2*r4 &
      **3
f(338) = -a1*r1**4*r4 - a2*r2**4*r3 + a3*r2*r3**4 + a4*r1*r4**4
f(339) = a1*r2*r3 + a2*r1*r4 - a3*r1*r4 - a4*r2*r3
f(340) = a1**2*r2*r3 + a2**2*r1*r4 - a3**2*r1*r4 - a4**2*r2*r3
f(341) = a1**3*r2*r3 + a2**3*r1*r4 - a3**3*r1*r4 - a4**3*r2*r3
f(342) = -a1**4*r2*r3 - a2**4*r1*r4 + a3**4*r1*r4 + a4**4*r2*r3
f(343) = -a1*r2*r3**2 - a2*r1*r4**2 + a3*r1**2*r4 + a4*r2**2*r3
f(344) = -a1**2*r2*r3**2 - a2**2*r1*r4**2 + a3**2*r1**2*r4 + a4**2*r2**2 &
      *r3
f(345) = a1**3*r2*r3**2 + a2**3*r1*r4**2 - a3**3*r1**2*r4 - a4**3*r2**2* &
      r3
f(346) = -a1*r2*r3**3 - a2*r1*r4**3 + a3*r1**3*r4 + a4*r2**3*r3
f(347) = -a1**2*r2*r3**3 - a2**2*r1*r4**3 + a3**2*r1**3*r4 + a4**2*r2**3 &
      *r3
f(348) = -a1*r2*r3**4 - a2*r1*r4**4 + a3*r1**4*r4 + a4*r2**4*r3
f(349) = -a1*r2**2*r3 - a2*r1**2*r4 + a3*r1*r4**2 + a4*r2*r3**2
f(350) = -a1**2*r2**2*r3 - a2**2*r1**2*r4 + a3**2*r1*r4**2 + a4**2*r2*r3 &
      **2
f(351) = a1**3*r2**2*r3 + a2**3*r1**2*r4 - a3**3*r1*r4**2 - a4**3*r2*r3 &
      **2
f(352) = a1*r2**2*r3**2 + a2*r1**2*r4**2 - a3*r1**2*r4**2 - a4*r2**2*r3 &
      **2
f(353) = a1**2*r2**2*r3**2 + a2**2*r1**2*r4**2 - a3**2*r1**2*r4**2 - a4 &
      **2*r2**2*r3**2
f(354) = a1*r2**2*r3**3 + a2*r1**2*r4**3 - a3*r1**3*r4**2 - a4*r2**3*r3 &
      **2
f(355) = a1*r2**3*r3 + a2*r1**3*r4 - a3*r1*r4**3 - a4*r2*r3**3
f(356) = a1**2*r2**3*r3 + a2**2*r1**3*r4 - a3**2*r1*r4**3 - a4**2*r2*r3 &
      **3
f(357) = a1*r2**3*r3**2 + a2*r1**3*r4**2 - a3*r1**2*r4**3 - a4*r2**2*r3 &
      **3
f(358) = a1*r2**4*r3 + a2*r1**4*r4 - a3*r1*r4**4 - a4*r2*r3**4
f(359) = b1**2*r1*r4 + b1**2*r2*r3 - b2**2*r1*r4 - b2**2*r2*r3
f(360) = -b1**4*r1*r4 - b1**4*r2*r3 + b2**4*r1*r4 + b2**4*r2*r3
f(361) = -b1**2*r1*r4**2 - b1**2*r2*r3**2 + b2**2*r1**2*r4 + b2**2*r2**2 &
      *r3
f(362) = -b1**2*r1*r4**3 - b1**2*r2*r3**3 + b2**2*r1**3*r4 + b2**2*r2**3 &
      *r3
f(363) = b1**2*r1**2*r4 + b1**2*r2**2*r3 - b2**2*r1*r4**2 - b2**2*r2*r3 &
      **2
f(364) = b1**2*r1**2*r4**2 + b1**2*r2**2*r3**2 - b2**2*r1**2*r4**2 - b2 &
      **2*r2**2*r3**2
f(365) = b1**2*r1**3*r4 + b1**2*r2**3*r3 - b2**2*r1*r4**3 - b2**2*r2*r3 &
      **3
f(366) = dtau**2*(-r1**2*r4 + r1*r4**2 - r2**2*r3 + r2*r3**2)
f(367) = dtau**2*(r1**3*r4 - r1*r4**3 + r2**3*r3 - r2*r3**3)
f(368) = -a1*a2*r1 - a1*a2*r2 + a3*a4*r3 + a3*a4*r4
f(369) = -a1**2*a2*r2 - a1*a2**2*r1 + a3**2*a4*r4 + a3*a4**2*r3
f(370) = a1**3*a2*r2 + a1*a2**3*r1 - a3**3*a4*r4 - a3*a4**3*r3
f(371) = -a1**4*a2*r2 - a1*a2**4*r1 + a3**4*a4*r4 + a3*a4**4*r3
f(372) = a1**2*a2*r1 + a1*a2**2*r2 - a3**2*a4*r3 - a3*a4**2*r4
f(373) = a1**2*a2**2*r1 + a1**2*a2**2*r2 - a3**2*a4**2*r3 - a3**2*a4**2* &
      r4
f(374) = -a1**3*a2**2*r2 - a1**2*a2**3*r1 + a3**3*a4**2*r4 + a3**2*a4**3 &
      *r3
f(375) = -a1**3*a2*r1 - a1*a2**3*r2 + a3**3*a4*r3 + a3*a4**3*r4
f(376) = a1**3*a2**2*r1 + a1**2*a2**3*r2 - a3**3*a4**2*r3 - a3**2*a4**3* &
      r4
f(377) = -a1**4*a2*r1 - a1*a2**4*r2 + a3**4*a4*r3 + a3*a4**4*r4
f(378) = a1*a2*r1**2 + a1*a2*r2**2 - a3*a4*r3**2 - a3*a4*r4**2
f(379) = -a1**2*a2*r2**2 - a1*a2**2*r1**2 + a3**2*a4*r4**2 + a3*a4**2*r3 &
      **2
f(380) = -a1**3*a2*r2**2 - a1*a2**3*r1**2 + a3**3*a4*r4**2 + a3*a4**3*r3 &
      **2
f(381) = -a1**2*a2*r1**2 - a1*a2**2*r2**2 + a3**2*a4*r3**2 + a3*a4**2*r4 &
      **2
f(382) = a1**2*a2**2*r1**2 + a1**2*a2**2*r2**2 - a3**2*a4**2*r3**2 - a3 &
      **2*a4**2*r4**2
f(383) = -a1**3*a2*r1**2 - a1*a2**3*r2**2 + a3**3*a4*r3**2 + a3*a4**3*r4 &
      **2
f(384) = a1*a2*r1**3 + a1*a2*r2**3 - a3*a4*r3**3 - a3*a4*r4**3
f(385) = -a1**2*a2*r2**3 - a1*a2**2*r1**3 + a3**2*a4*r4**3 + a3*a4**2*r3 &
      **3
f(386) = -a1**2*a2*r1**3 - a1*a2**2*r2**3 + a3**2*a4*r3**3 + a3*a4**2*r4 &
      **3
f(387) = a1*a2*r1**4 + a1*a2*r2**4 - a3*a4*r3**4 - a3*a4*r4**4
f(388) = -a1*a3*r1 + a1*a3*r3 - a2*a4*r2 + a2*a4*r4
f(389) = -a1**2*a3*r3 + a1*a3**2*r1 - a2**2*a4*r4 + a2*a4**2*r2
f(390) = a1**3*a3*r3 - a1*a3**3*r1 + a2**3*a4*r4 - a2*a4**3*r2
f(391) = -a1**4*a3*r3 + a1*a3**4*r1 - a2**4*a4*r4 + a2*a4**4*r2
f(392) = -a1**2*a3*r1 + a1*a3**2*r3 - a2**2*a4*r2 + a2*a4**2*r4
f(393) = a1**2*a3**2*r1 - a1**2*a3**2*r3 + a2**2*a4**2*r2 - a2**2*a4**2* &
      r4
f(394) = -a1**3*a3**2*r3 + a1**2*a3**3*r1 - a2**3*a4**2*r4 + a2**2*a4**3 &
      *r2
f(395) = -a1**3*a3*r1 + a1*a3**3*r3 - a2**3*a4*r2 + a2*a4**3*r4
f(396) = a1**3*a3**2*r1 - a1**2*a3**3*r3 + a2**3*a4**2*r2 - a2**2*a4**3* &
      r4
f(397) = -a1**4*a3*r1 + a1*a3**4*r3 - a2**4*a4*r2 + a2*a4**4*r4
f(398) = -a1*a3*r1**2 + a1*a3*r3**2 - a2*a4*r2**2 + a2*a4*r4**2
f(399) = -a1**2*a3*r3**2 + a1*a3**2*r1**2 - a2**2*a4*r4**2 + a2*a4**2*r2 &
      **2
f(400) = -a1**3*a3*r3**2 + a1*a3**3*r1**2 - a2**3*a4*r4**2 + a2*a4**3*r2 &
      **2
f(401) = -a1**2*a3*r1**2 + a1*a3**2*r3**2 - a2**2*a4*r2**2 + a2*a4**2*r4 &
      **2
f(402) = a1**2*a3**2*r1**2 - a1**2*a3**2*r3**2 + a2**2*a4**2*r2**2 - a2 &
      **2*a4**2*r4**2
f(403) = a1**3*a3*r1**2 - a1*a3**3*r3**2 + a2**3*a4*r2**2 - a2*a4**3*r4 &
      **2
f(404) = a1*a3*r1**3 - a1*a3*r3**3 + a2*a4*r2**3 - a2*a4*r4**3
f(405) = a1**2*a3*r3**3 - a1*a3**2*r1**3 + a2**2*a4*r4**3 - a2*a4**2*r2 &
      **3
f(406) = a1**2*a3*r1**3 - a1*a3**2*r3**3 + a2**2*a4*r2**3 - a2*a4**2*r4 &
      **3
f(407) = a1*a3*r1**4 - a1*a3*r3**4 + a2*a4*r2**4 - a2*a4*r4**4
f(408) = a1*a4*r1 - a1*a4*r4 + a2*a3*r2 - a2*a3*r3
f(409) = -a1**2*a4*r4 + a1*a4**2*r1 - a2**2*a3*r3 + a2*a3**2*r2
f(410) = a1**3*a4*r4 - a1*a4**3*r1 + a2**3*a3*r3 - a2*a3**3*r2
f(411) = -a1**4*a4*r4 + a1*a4**4*r1 - a2**4*a3*r3 + a2*a3**4*r2
f(412) = a1**2*a4*r1 - a1*a4**2*r4 + a2**2*a3*r2 - a2*a3**2*r3
f(413) = -a1**2*a4**2*r1 + a1**2*a4**2*r4 - a2**2*a3**2*r2 + a2**2*a3**2 &
      *r3
f(414) = a1**3*a4**2*r4 - a1**2*a4**3*r1 + a2**3*a3**2*r3 - a2**2*a3**3* &
      r2
f(415) = a1**3*a4*r1 - a1*a4**3*r4 + a2**3*a3*r2 - a2*a3**3*r3
f(416) = -a1**3*a4**2*r1 + a1**2*a4**3*r4 - a2**3*a3**2*r2 + a2**2*a3**3 &
      *r3
f(417) = -a1**4*a4*r1 + a1*a4**4*r4 - a2**4*a3*r2 + a2*a3**4*r3
f(418) = a1*a4*r1**2 - a1*a4*r4**2 + a2*a3*r2**2 - a2*a3*r3**2
f(419) = a1**2*a4*r4**2 - a1*a4**2*r1**2 + a2**2*a3*r3**2 - a2*a3**2*r2 &
      **2
f(420) = a1**3*a4*r4**2 - a1*a4**3*r1**2 + a2**3*a3*r3**2 - a2*a3**3*r2 &
      **2
f(421) = a1**2*a4*r1**2 - a1*a4**2*r4**2 + a2**2*a3*r2**2 - a2*a3**2*r3 &
      **2
f(422) = -a1**2*a4**2*r1**2 + a1**2*a4**2*r4**2 - a2**2*a3**2*r2**2 + a2 &
      **2*a3**2*r3**2
f(423) = -a1**3*a4*r1**2 + a1*a4**3*r4**2 - a2**3*a3*r2**2 + a2*a3**3*r3 &
      **2
f(424) = -a1*a4*r1**3 + a1*a4*r4**3 - a2*a3*r2**3 + a2*a3*r3**3
f(425) = -a1**2*a4*r4**3 + a1*a4**2*r1**3 - a2**2*a3*r3**3 + a2*a3**2*r2 &
      **3
f(426) = -a1**2*a4*r1**3 + a1*a4**2*r4**3 - a2**2*a3*r2**3 + a2*a3**2*r3 &
      **3
f(427) = a1*a4*r1**4 - a1*a4*r4**4 + a2*a3*r2**4 - a2*a3*r3**4
f(428) = a1*b1**2*r1 + a2*b1**2*r2 - a3*b2**2*r3 - a4*b2**2*r4
f(429) = -a1*b1**4*r1 - a2*b1**4*r2 + a3*b2**4*r3 + a4*b2**4*r4
f(430) = -a1**2*b1**2*r1 - a2**2*b1**2*r2 + a3**2*b2**2*r3 + a4**2*b2**2 &
      *r4
f(431) = -a1**3*b1**2*r1 - a2**3*b1**2*r2 + a3**3*b2**2*r3 + a4**3*b2**2 &
      *r4
f(432) = -a1*b1**2*r1**2 - a2*b1**2*r2**2 + a3*b2**2*r3**2 + a4*b2**2*r4 &
      **2
f(433) = a1**2*b1**2*r1**2 + a2**2*b1**2*r2**2 - a3**2*b2**2*r3**2 - a4 &
      **2*b2**2*r4**2
f(434) = a1*b1**2*r1**3 + a2*b1**2*r2**3 - a3*b2**2*r3**3 - a4*b2**2*r4 &
      **3
f(435) = -a1*b2**2*r1 - a2*b2**2*r2 + a3*b1**2*r3 + a4*b1**2*r4
f(436) = -a1*b2**4*r1 - a2*b2**4*r2 + a3*b1**4*r3 + a4*b1**4*r4
f(437) = -a1**2*b2**2*r1 - a2**2*b2**2*r2 + a3**2*b1**2*r3 + a4**2*b1**2 &
      *r4
f(438) = a1**3*b2**2*r1 + a2**3*b2**2*r2 - a3**3*b1**2*r3 - a4**3*b1**2* &
      r4
f(439) = a1*b2**2*r1**2 + a2*b2**2*r2**2 - a3*b1**2*r3**2 - a4*b1**2*r4 &
      **2
f(440) = -a1**2*b2**2*r1**2 - a2**2*b2**2*r2**2 + a3**2*b1**2*r3**2 + a4 &
      **2*b1**2*r4**2
f(441) = -a1*b2**2*r1**3 - a2*b2**2*r2**3 + a3*b1**2*r3**3 + a4*b1**2*r4 &
      **3
f(442) = dtau**2*(-a1*r1 - a2*r2 + a3*r3 + a4*r4)
f(443) = dtau**4*(-a1*r1 - a2*r2 + a3*r3 + a4*r4)
f(444) = dtau**2*(-a1**2*r1 - a2**2*r2 + a3**2*r3 + a4**2*r4)
f(445) = dtau**2*(-a1**3*r1 - a2**3*r2 + a3**3*r3 + a4**3*r4)
f(446) = dtau**2*(-a1*r1**2 - a2*r2**2 + a3*r3**2 + a4*r4**2)
f(447) = dtau**2*(-a1**2*r1**2 - a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2 &
      )
f(448) = dtau**2*(a1*r1**3 + a2*r2**3 - a3*r3**3 - a4*r4**3)
f(449) = -a1*a4*r2 + a1*a4*r3 - a2*a3*r1 + a2*a3*r4
f(450) = -a1**2*a4*r3 + a1*a4**2*r2 - a2**2*a3*r4 + a2*a3**2*r1
f(451) = a1**3*a4*r3 - a1*a4**3*r2 + a2**3*a3*r4 - a2*a3**3*r1
f(452) = a1**4*a4*r3 - a1*a4**4*r2 + a2**4*a3*r4 - a2*a3**4*r1
f(453) = a1**2*a4*r2 - a1*a4**2*r3 + a2**2*a3*r1 - a2*a3**2*r4
f(454) = -a1**2*a4**2*r2 + a1**2*a4**2*r3 - a2**2*a3**2*r1 + a2**2*a3**2 &
      *r4
f(455) = a1**3*a4**2*r3 - a1**2*a4**3*r2 + a2**3*a3**2*r4 - a2**2*a3**3* &
      r1
f(456) = a1**3*a4*r2 - a1*a4**3*r3 + a2**3*a3*r1 - a2*a3**3*r4
f(457) = -a1**3*a4**2*r2 + a1**2*a4**3*r3 - a2**3*a3**2*r1 + a2**2*a3**3 &
      *r4
f(458) = a1**4*a4*r2 - a1*a4**4*r3 + a2**4*a3*r1 - a2*a3**4*r4
f(459) = -a1*a4*r2**2 + a1*a4*r3**2 - a2*a3*r1**2 + a2*a3*r4**2
f(460) = -a1**2*a4*r3**2 + a1*a4**2*r2**2 - a2**2*a3*r4**2 + a2*a3**2*r1 &
      **2
f(461) = -a1**3*a4*r3**2 + a1*a4**3*r2**2 - a2**3*a3*r4**2 + a2*a3**3*r1 &
      **2
f(462) = -a1**2*a4*r2**2 + a1*a4**2*r3**2 - a2**2*a3*r1**2 + a2*a3**2*r4 &
      **2
f(463) = a1**2*a4**2*r2**2 - a1**2*a4**2*r3**2 + a2**2*a3**2*r1**2 - a2 &
      **2*a3**2*r4**2
f(464) = a1**3*a4*r2**2 - a1*a4**3*r3**2 + a2**3*a3*r1**2 - a2*a3**3*r4 &
      **2
f(465) = a1*a4*r2**3 - a1*a4*r3**3 + a2*a3*r1**3 - a2*a3*r4**3
f(466) = a1**2*a4*r3**3 - a1*a4**2*r2**3 + a2**2*a3*r4**3 - a2*a3**2*r1 &
      **3
f(467) = a1**2*a4*r2**3 - a1*a4**2*r3**3 + a2**2*a3*r1**3 - a2*a3**2*r4 &
      **3
f(468) = -a1*a4*r2**4 + a1*a4*r3**4 - a2*a3*r1**4 + a2*a3*r4**4
f(469) = a1*a3*r2 - a1*a3*r4 + a2*a4*r1 - a2*a4*r3
f(470) = -a1**2*a3*r4 + a1*a3**2*r2 - a2**2*a4*r3 + a2*a4**2*r1
f(471) = a1**3*a3*r4 - a1*a3**3*r2 + a2**3*a4*r3 - a2*a4**3*r1
f(472) = a1**4*a3*r4 - a1*a3**4*r2 + a2**4*a4*r3 - a2*a4**4*r1
f(473) = -a1**2*a3*r2 + a1*a3**2*r4 - a2**2*a4*r1 + a2*a4**2*r3
f(474) = a1**2*a3**2*r2 - a1**2*a3**2*r4 + a2**2*a4**2*r1 - a2**2*a4**2* &
      r3
f(475) = -a1**3*a3**2*r4 + a1**2*a3**3*r2 - a2**3*a4**2*r3 + a2**2*a4**3 &
      *r1
f(476) = -a1**3*a3*r2 + a1*a3**3*r4 - a2**3*a4*r1 + a2*a4**3*r3
f(477) = a1**3*a3**2*r2 - a1**2*a3**3*r4 + a2**3*a4**2*r1 - a2**2*a4**3* &
      r3
f(478) = a1**4*a3*r2 - a1*a3**4*r4 + a2**4*a4*r1 - a2*a4**4*r3
f(479) = a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 - a2*a4*r3**2
f(480) = a1**2*a3*r4**2 - a1*a3**2*r2**2 + a2**2*a4*r3**2 - a2*a4**2*r1 &
      **2
f(481) = a1**3*a3*r4**2 - a1*a3**3*r2**2 + a2**3*a4*r3**2 - a2*a4**3*r1 &
      **2
f(482) = a1**2*a3*r2**2 - a1*a3**2*r4**2 + a2**2*a4*r1**2 - a2*a4**2*r3 &
      **2
f(483) = -a1**2*a3**2*r2**2 + a1**2*a3**2*r4**2 - a2**2*a4**2*r1**2 + a2 &
      **2*a4**2*r3**2
f(484) = -a1**3*a3*r2**2 + a1*a3**3*r4**2 - a2**3*a4*r1**2 + a2*a4**3*r3 &
      **2
f(485) = -a1*a3*r2**3 + a1*a3*r4**3 - a2*a4*r1**3 + a2*a4*r3**3
f(486) = -a1**2*a3*r4**3 + a1*a3**2*r2**3 - a2**2*a4*r3**3 + a2*a4**2*r1 &
      **3
f(487) = -a1**2*a3*r2**3 + a1*a3**2*r4**3 - a2**2*a4*r1**3 + a2*a4**2*r3 &
      **3
f(488) = -a1*a3*r2**4 + a1*a3*r4**4 - a2*a4*r1**4 + a2*a4*r3**4
f(489) = -a1*b1**2*r2 - a2*b1**2*r1 + a3*b2**2*r4 + a4*b2**2*r3
f(490) = a1*b1**4*r2 + a2*b1**4*r1 - a3*b2**4*r4 - a4*b2**4*r3
f(491) = a1**2*b1**2*r2 + a2**2*b1**2*r1 - a3**2*b2**2*r4 - a4**2*b2**2* &
      r3
f(492) = a1**3*b1**2*r2 + a2**3*b1**2*r1 - a3**3*b2**2*r4 - a4**3*b2**2* &
      r3
f(493) = a1*b1**2*r2**2 + a2*b1**2*r1**2 - a3*b2**2*r4**2 - a4*b2**2*r3 &
      **2
f(494) = -a1**2*b1**2*r2**2 - a2**2*b1**2*r1**2 + a3**2*b2**2*r4**2 + a4 &
      **2*b2**2*r3**2
f(495) = -a1*b1**2*r2**3 - a2*b1**2*r1**3 + a3*b2**2*r4**3 + a4*b2**2*r3 &
      **3
f(496) = a1*b2**2*r2 + a2*b2**2*r1 - a3*b1**2*r4 - a4*b1**2*r3
f(497) = -a1*b2**4*r2 - a2*b2**4*r1 + a3*b1**4*r4 + a4*b1**4*r3
f(498) = -a1**2*b2**2*r2 - a2**2*b2**2*r1 + a3**2*b1**2*r4 + a4**2*b1**2 &
      *r3
f(499) = a1**3*b2**2*r2 + a2**3*b2**2*r1 - a3**3*b1**2*r4 - a4**3*b1**2* &
      r3
f(500) = -a1*b2**2*r2**2 - a2*b2**2*r1**2 + a3*b1**2*r4**2 + a4*b1**2*r3 &
      **2
f(501) = a1**2*b2**2*r2**2 + a2**2*b2**2*r1**2 - a3**2*b1**2*r4**2 - a4 &
      **2*b1**2*r3**2
f(502) = a1*b2**2*r2**3 + a2*b2**2*r1**3 - a3*b1**2*r4**3 - a4*b1**2*r3 &
      **3
f(503) = dtau**2*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(504) = dtau**4*(-a1*r2 - a2*r1 + a3*r4 + a4*r3)
f(505) = dtau**2*(a1**2*r2 + a2**2*r1 - a3**2*r4 - a4**2*r3)
f(506) = dtau**2*(a1**3*r2 + a2**3*r1 - a3**3*r4 - a4**3*r3)
f(507) = dtau**2*(a1*r2**2 + a2*r1**2 - a3*r4**2 - a4*r3**2)
f(508) = dtau**2*(-a1**2*r2**2 - a2**2*r1**2 + a3**2*r4**2 + a4**2*r3**2 &
      )
f(509) = dtau**2*(-a1*r2**3 - a2*r1**3 + a3*r4**3 + a4*r3**3)
f(510) = -a1*a2*r3 - a1*a2*r4 + a3*a4*r1 + a3*a4*r2
f(511) = a1**2*a2*r4 + a1*a2**2*r3 - a3**2*a4*r2 - a3*a4**2*r1
f(512) = a1**3*a2*r4 + a1*a2**3*r3 - a3**3*a4*r2 - a3*a4**3*r1
f(513) = a1**4*a2*r4 + a1*a2**4*r3 - a3**4*a4*r2 - a3*a4**4*r1
f(514) = -a1**2*a2*r3 - a1*a2**2*r4 + a3**2*a4*r1 + a3*a4**2*r2
f(515) = -a1**2*a2**2*r3 - a1**2*a2**2*r4 + a3**2*a4**2*r1 + a3**2*a4**2 &
      *r2
f(516) = a1**3*a2**2*r4 + a1**2*a2**3*r3 - a3**3*a4**2*r2 - a3**2*a4**3* &
      r1
f(517) = -a1**3*a2*r3 - a1*a2**3*r4 + a3**3*a4*r1 + a3*a4**3*r2
f(518) = -a1**3*a2**2*r3 - a1**2*a2**3*r4 + a3**3*a4**2*r1 + a3**2*a4**3 &
      *r2
f(519) = a1**4*a2*r3 + a1*a2**4*r4 - a3**4*a4*r1 - a3*a4**4*r2
f(520) = a1*a2*r3**2 + a1*a2*r4**2 - a3*a4*r1**2 - a3*a4*r2**2
f(521) = a1**2*a2*r4**2 + a1*a2**2*r3**2 - a3**2*a4*r2**2 - a3*a4**2*r1 &
      **2
f(522) = -a1**3*a2*r4**2 - a1*a2**3*r3**2 + a3**3*a4*r2**2 + a3*a4**3*r1 &
      **2
f(523) = a1**2*a2*r3**2 + a1*a2**2*r4**2 - a3**2*a4*r1**2 - a3*a4**2*r2 &
      **2
f(524) = a1**2*a2**2*r3**2 + a1**2*a2**2*r4**2 - a3**2*a4**2*r1**2 - a3 &
      **2*a4**2*r2**2
f(525) = -a1**3*a2*r3**2 - a1*a2**3*r4**2 + a3**3*a4*r1**2 + a3*a4**3*r2 &
      **2
f(526) = a1*a2*r3**3 + a1*a2*r4**3 - a3*a4*r1**3 - a3*a4*r2**3
f(527) = a1**2*a2*r4**3 + a1*a2**2*r3**3 - a3**2*a4*r2**3 - a3*a4**2*r1 &
      **3
f(528) = a1**2*a2*r3**3 + a1*a2**2*r4**3 - a3**2*a4*r1**3 - a3*a4**2*r2 &
      **3
f(529) = a1*a2*r3**4 + a1*a2*r4**4 - a3*a4*r1**4 - a3*a4*r2**4
f(530) = a1*b2**2*r3 + a2*b2**2*r4 - a3*b1**2*r1 - a4*b1**2*r2
f(531) = -a1*b2**4*r3 - a2*b2**4*r4 + a3*b1**4*r1 + a4*b1**4*r2
f(532) = -a1**2*b2**2*r3 - a2**2*b2**2*r4 + a3**2*b1**2*r1 + a4**2*b1**2 &
      *r2
f(533) = a1**3*b2**2*r3 + a2**3*b2**2*r4 - a3**3*b1**2*r1 - a4**3*b1**2* &
      r2
f(534) = -a1*b2**2*r3**2 - a2*b2**2*r4**2 + a3*b1**2*r1**2 + a4*b1**2*r2 &
      **2
f(535) = a1**2*b2**2*r3**2 + a2**2*b2**2*r4**2 - a3**2*b1**2*r1**2 - a4 &
      **2*b1**2*r2**2
f(536) = a1*b2**2*r3**3 + a2*b2**2*r4**3 - a3*b1**2*r1**3 - a4*b1**2*r2 &
      **3
f(537) = -a1*b1**2*r3 - a2*b1**2*r4 + a3*b2**2*r1 + a4*b2**2*r2
f(538) = a1*b1**4*r3 + a2*b1**4*r4 - a3*b2**4*r1 - a4*b2**4*r2
f(539) = a1**2*b1**2*r3 + a2**2*b1**2*r4 - a3**2*b2**2*r1 - a4**2*b2**2* &
      r2
f(540) = -a1**3*b1**2*r3 - a2**3*b1**2*r4 + a3**3*b2**2*r1 + a4**3*b2**2 &
      *r2
f(541) = a1*b1**2*r3**2 + a2*b1**2*r4**2 - a3*b2**2*r1**2 - a4*b2**2*r2 &
      **2
f(542) = -a1**2*b1**2*r3**2 - a2**2*b1**2*r4**2 + a3**2*b2**2*r1**2 + a4 &
      **2*b2**2*r2**2
f(543) = -a1*b1**2*r3**3 - a2*b1**2*r4**3 + a3*b2**2*r1**3 + a4*b2**2*r2 &
      **3
f(544) = dtau**2*(a1*r3 + a2*r4 - a3*r1 - a4*r2)
f(545) = dtau**4*(a1*r3 + a2*r4 - a3*r1 - a4*r2)
f(546) = dtau**2*(a1**2*r3 + a2**2*r4 - a3**2*r1 - a4**2*r2)
f(547) = dtau**2*(a1**3*r3 + a2**3*r4 - a3**3*r1 - a4**3*r2)
f(548) = dtau**2*(-a1*r3**2 - a2*r4**2 + a3*r1**2 + a4*r2**2)
f(549) = dtau**2*(-a1**2*r3**2 - a2**2*r4**2 + a3**2*r1**2 + a4**2*r2**2 &
      )
f(550) = dtau**2*(a1*r3**3 + a2*r4**3 - a3*r1**3 - a4*r2**3)
f(551) = -a1*b2**2*r4 - a2*b2**2*r3 + a3*b1**2*r2 + a4*b1**2*r1
f(552) = a1*b2**4*r4 + a2*b2**4*r3 - a3*b1**4*r2 - a4*b1**4*r1
f(553) = a1**2*b2**2*r4 + a2**2*b2**2*r3 - a3**2*b1**2*r2 - a4**2*b1**2* &
      r1
f(554) = -a1**3*b2**2*r4 - a2**3*b2**2*r3 + a3**3*b1**2*r2 + a4**3*b1**2 &
      *r1
f(555) = a1*b2**2*r4**2 + a2*b2**2*r3**2 - a3*b1**2*r2**2 - a4*b1**2*r1 &
      **2
f(556) = -a1**2*b2**2*r4**2 - a2**2*b2**2*r3**2 + a3**2*b1**2*r2**2 + a4 &
      **2*b1**2*r1**2
f(557) = -a1*b2**2*r4**3 - a2*b2**2*r3**3 + a3*b1**2*r2**3 + a4*b1**2*r1 &
      **3
f(558) = a1*b1**2*r4 + a2*b1**2*r3 - a3*b2**2*r2 - a4*b2**2*r1
f(559) = -a1*b1**4*r4 - a2*b1**4*r3 + a3*b2**4*r2 + a4*b2**4*r1
f(560) = -a1**2*b1**2*r4 - a2**2*b1**2*r3 + a3**2*b2**2*r2 + a4**2*b2**2 &
      *r1
f(561) = a1**3*b1**2*r4 + a2**3*b1**2*r3 - a3**3*b2**2*r2 - a4**3*b2**2* &
      r1
f(562) = -a1*b1**2*r4**2 - a2*b1**2*r3**2 + a3*b2**2*r2**2 + a4*b2**2*r1 &
      **2
f(563) = a1**2*b1**2*r4**2 + a2**2*b1**2*r3**2 - a3**2*b2**2*r2**2 - a4 &
      **2*b2**2*r1**2
f(564) = a1*b1**2*r4**3 + a2*b1**2*r3**3 - a3*b2**2*r2**3 - a4*b2**2*r1 &
      **3
f(565) = dtau**2*(-a1*r4 - a2*r3 + a3*r2 + a4*r1)
f(566) = dtau**4*(a1*r4 + a2*r3 - a3*r2 - a4*r1)
f(567) = dtau**2*(a1**2*r4 + a2**2*r3 - a3**2*r2 - a4**2*r1)
f(568) = dtau**2*(a1**3*r4 + a2**3*r3 - a3**3*r2 - a4**3*r1)
f(569) = dtau**2*(a1*r4**2 + a2*r3**2 - a3*r2**2 - a4*r1**2)
f(570) = dtau**2*(a1**2*r4**2 + a2**2*r3**2 - a3**2*r2**2 - a4**2*r1**2)
f(571) = dtau**2*(-a1*r4**3 - a2*r3**3 + a3*r2**3 + a4*r1**3)
f(572) = b1*b2*(-r1 - r2 + r3 + r4)
f(573) = b1*b2*(-b1**2*r3 - b1**2*r4 + b2**2*r1 + b2**2*r2)
f(574) = b1**2*b2**2*(-r1 - r2 + r3 + r4)
f(575) = b1*b2*(b1**2*r1 + b1**2*r2 - b2**2*r3 - b2**2*r4)
f(576) = b1*b2*(r1**2 + r2**2 - r3**2 - r4**2)
f(577) = b1*b2*(b1**2*r3**2 + b1**2*r4**2 - b2**2*r1**2 - b2**2*r2**2)
f(578) = b1**2*b2**2*(r1**2 + r2**2 - r3**2 - r4**2)
f(579) = b1*b2*(-b1**2*r1**2 - b1**2*r2**2 + b2**2*r3**2 + b2**2*r4**2)
f(580) = b1*b2*(r1**3 + r2**3 - r3**3 - r4**3)
f(581) = b1*b2*(-r1**4 - r2**4 + r3**4 + r4**4)
f(582) = dtau*(b1*r1 - b1*r2 - b2*r3 + b2*r4)
f(583) = dtau**3*(b1*r1 - b1*r2 - b2*r3 + b2*r4)
f(584) = dtau**2*(b1**2*r1 + b1**2*r2 - b2**2*r3 - b2**2*r4)
f(585) = dtau*(-b1**3*r1 + b1**3*r2 + b2**3*r3 - b2**3*r4)
f(586) = dtau*(-b1*r1**2 + b1*r2**2 + b2*r3**2 - b2*r4**2)
f(587) = dtau**3*(-b1*r1**2 + b1*r2**2 + b2*r3**2 - b2*r4**2)
f(588) = dtau**2*(-b1**2*r1**2 - b1**2*r2**2 + b2**2*r3**2 + b2**2*r4**2 &
      )
f(589) = dtau*(-b1**3*r1**2 + b1**3*r2**2 + b2**3*r3**2 - b2**3*r4**2)
f(590) = dtau*(b1*r1**3 - b1*r2**3 - b2*r3**3 + b2*r4**3)
f(591) = dtau*(-b1*r1**4 + b1*r2**4 + b2*r3**4 - b2*r4**4)
f(592) = dtau*(-b1*r3 + b1*r4 + b2*r1 - b2*r2)
f(593) = dtau**3*(b1*r3 - b1*r4 - b2*r1 + b2*r2)
f(594) = dtau**2*(-b1**2*r3 - b1**2*r4 + b2**2*r1 + b2**2*r2)
f(595) = dtau*(b1**3*r3 - b1**3*r4 - b2**3*r1 + b2**3*r2)
f(596) = dtau*(-b1*r3**2 + b1*r4**2 + b2*r1**2 - b2*r2**2)
f(597) = dtau**3*(b1*r3**2 - b1*r4**2 - b2*r1**2 + b2*r2**2)
f(598) = dtau**2*(-b1**2*r3**2 - b1**2*r4**2 + b2**2*r1**2 + b2**2*r2**2 &
      )
f(599) = dtau*(b1**3*r3**2 - b1**3*r4**2 - b2**3*r1**2 + b2**3*r2**2)
f(600) = dtau*(b1*r3**3 - b1*r4**3 - b2*r1**3 + b2*r2**3)
f(601) = dtau*(-b1*r3**4 + b1*r4**4 + b2*r1**4 - b2*r2**4)
f(602) = -a1*a2*a3 - a1*a2*a4 + a1*a3*a4 + a2*a3*a4
f(603) = a1**2*a3*a4 - a1*a2*a3**2 - a1*a2*a4**2 + a2**2*a3*a4
f(604) = -a1**3*a3*a4 + a1*a2*a3**3 + a1*a2*a4**3 - a2**3*a3*a4
f(605) = -a1**4*a3*a4 + a1*a2*a3**4 + a1*a2*a4**4 - a2**4*a3*a4
f(606) = -a1**2*a2*a4 - a1*a2**2*a3 + a1*a3*a4**2 + a2*a3**2*a4
f(607) = a1**2*a2*a4**2 - a1**2*a3*a4**2 + a1*a2**2*a3**2 - a2**2*a3**2* &
      a4
f(608) = -a1**3*a3*a4**2 + a1**2*a2*a4**3 + a1*a2**2*a3**3 - a2**3*a3**2 &
      *a4
f(609) = -a1**3*a2*a4 - a1*a2**3*a3 + a1*a3*a4**3 + a2*a3**3*a4
f(610) = a1**3*a2*a4**2 - a1**2*a3*a4**3 + a1*a2**3*a3**2 - a2**2*a3**3* &
      a4
f(611) = -a1**4*a2*a4 - a1*a2**4*a3 + a1*a3*a4**4 + a2*a3**4*a4
f(612) = a1**2*a2*a3 + a1*a2**2*a4 - a1*a3**2*a4 - a2*a3*a4**2
f(613) = -a1**2*a2*a3**2 + a1**2*a3**2*a4 - a1*a2**2*a4**2 + a2**2*a3*a4 &
      **2
f(614) = -a1**3*a3**2*a4 + a1**2*a2*a3**3 + a1*a2**2*a4**3 - a2**3*a3*a4 &
      **2
f(615) = a1**2*a2**2*a3 + a1**2*a2**2*a4 - a1*a3**2*a4**2 - a2*a3**2*a4 &
      **2
f(616) = -a1**2*a2**2*a3**2 - a1**2*a2**2*a4**2 + a1**2*a3**2*a4**2 + a2 &
      **2*a3**2*a4**2
f(617) = -a1**3*a2**2*a4 - a1**2*a2**3*a3 + a1*a3**2*a4**3 + a2*a3**3*a4 &
      **2
f(618) = a1**3*a2*a3 + a1*a2**3*a4 - a1*a3**3*a4 - a2*a3*a4**3
f(619) = a1**3*a2*a3**2 - a1**2*a3**3*a4 + a1*a2**3*a4**2 - a2**2*a3*a4 &
      **3
f(620) = a1**3*a2**2*a3 + a1**2*a2**3*a4 - a1*a3**3*a4**2 - a2*a3**2*a4 &
      **3
f(621) = -a1**4*a2*a3 - a1*a2**4*a4 + a1*a3**4*a4 + a2*a3*a4**4
f(622) = a1*a2*b1**2 - a3*a4*b2**2
f(623) = a1*a2*b1**4 - a3*a4*b2**4
f(624) = a1**2*a2*b1**2 + a1*a2**2*b1**2 - a3**2*a4*b2**2 - a3*a4**2*b2 &
      **2
f(625) = -a1**3*a2*b1**2 - a1*a2**3*b1**2 + a3**3*a4*b2**2 + a3*a4**3*b2 &
      **2
f(626) = a1**2*a2**2*b1**2 - a3**2*a4**2*b2**2
f(627) = a1*a2*b2**2 - a3*a4*b1**2
f(628) = -a1*a2*b2**4 + a3*a4*b1**4
f(629) = -a1**2*a2*b2**2 - a1*a2**2*b2**2 + a3**2*a4*b1**2 + a3*a4**2*b1 &
      **2
f(630) = -a1**3*a2*b2**2 - a1*a2**3*b2**2 + a3**3*a4*b1**2 + a3*a4**3*b1 &
      **2
f(631) = a1**2*a2**2*b2**2 - a3**2*a4**2*b1**2
f(632) = a1**3*a2*b2**2 + a1*a2**3*b2**2 - a3**3*a4*b1**2 - a3*a4**3*b1 &
      **2
f(633) = dtau**2*(-a1*a2 + a3*a4)
f(634) = dtau**4*(a1*a2 - a3*a4)
f(635) = dtau**2*(-a1**2*a2 - a1*a2**2 + a3**2*a4 + a3*a4**2)
f(636) = dtau**2*(a1**3*a2 + a1*a2**3 - a3**3*a4 - a3*a4**3)
f(637) = dtau**2*(a1**2*a2 + a1*a2**2 - a3**2*a4 - a3*a4**2)
f(638) = dtau**2*(-a1**2*a2**2 + a3**2*a4**2)
f(639) = -a1*a3*b1**2 + a1*a3*b2**2 - a2*a4*b1**2 + a2*a4*b2**2
f(640) = a1*a3*b1**4 - a1*a3*b2**4 + a2*a4*b1**4 - a2*a4*b2**4
f(641) = -a1**2*a3*b2**2 + a1*a3**2*b1**2 - a2**2*a4*b2**2 + a2*a4**2*b1 &
      **2
f(642) = a1**3*a3*b2**2 - a1*a3**3*b1**2 + a2**3*a4*b2**2 - a2*a4**3*b1 &
      **2
f(643) = a1**2*a3*b1**2 - a1*a3**2*b2**2 + a2**2*a4*b1**2 - a2*a4**2*b2 &
      **2
f(644) = -a1**2*a3**2*b1**2 + a1**2*a3**2*b2**2 - a2**2*a4**2*b1**2 + a2 &
      **2*a4**2*b2**2
f(645) = -a1**3*a3*b1**2 + a1*a3**3*b2**2 - a2**3*a4*b1**2 + a2*a4**3*b2 &
      **2
f(646) = dtau**2*(-a1**2*a3 + a1*a3**2 - a2**2*a4 + a2*a4**2)
f(647) = dtau**2*(a1**3*a3 - a1*a3**3 + a2**3*a4 - a2*a4**3)
f(648) = a1*a4*b1**2 - a1*a4*b2**2 + a2*a3*b1**2 - a2*a3*b2**2
f(649) = -a1*a4*b1**4 + a1*a4*b2**4 - a2*a3*b1**4 + a2*a3*b2**4
f(650) = a1**2*a4*b2**2 - a1*a4**2*b1**2 + a2**2*a3*b2**2 - a2*a3**2*b1 &
      **2
f(651) = -a1**3*a4*b2**2 + a1*a4**3*b1**2 - a2**3*a3*b2**2 + a2*a3**3*b1 &
      **2
f(652) = -a1**2*a4*b1**2 + a1*a4**2*b2**2 - a2**2*a3*b1**2 + a2*a3**2*b2 &
      **2
f(653) = a1**2*a4**2*b1**2 - a1**2*a4**2*b2**2 + a2**2*a3**2*b1**2 - a2 &
      **2*a3**2*b2**2
f(654) = a1**3*a4*b1**2 - a1*a4**3*b2**2 + a2**3*a3*b1**2 - a2*a3**3*b2 &
      **2
f(655) = dtau**2*(-a1**2*a4 + a1*a4**2 - a2**2*a3 + a2*a3**2)
f(656) = dtau**2*(-a1**3*a4 + a1*a4**3 - a2**3*a3 + a2*a3**3)
f(657) = b1*b2*(-a1 - a2 + a3 + a4)
f(658) = b1*b2*(-a1*b2**2 - a2*b2**2 + a3*b1**2 + a4*b1**2)
f(659) = b1**2*b2**2*(a1 + a2 - a3 - a4)
f(660) = b1*b2*(a1*b1**2 + a2*b1**2 - a3*b2**2 - a4*b2**2)
f(661) = b1*b2*(-a1**2 - a2**2 + a3**2 + a4**2)
f(662) = b1*b2*(-a1**2*b2**2 - a2**2*b2**2 + a3**2*b1**2 + a4**2*b1**2)
f(663) = b1**2*b2**2*(a1**2 + a2**2 - a3**2 - a4**2)
f(664) = b1*b2*(a1**2*b1**2 + a2**2*b1**2 - a3**2*b2**2 - a4**2*b2**2)
f(665) = b1*b2*(-a1**3 - a2**3 + a3**3 + a4**3)
f(666) = b1*b2*(a1**4 + a2**4 - a3**4 - a4**4)
f(667) = dtau*(-a1*b1 + a2*b1 + a3*b2 - a4*b2)
f(668) = dtau**3*(-a1*b1 + a2*b1 + a3*b2 - a4*b2)
f(669) = dtau**2*(a1*b1**2 + a2*b1**2 - a3*b2**2 - a4*b2**2)
f(670) = dtau*(a1*b1**3 - a2*b1**3 - a3*b2**3 + a4*b2**3)
f(671) = dtau*(a1**2*b1 - a2**2*b1 - a3**2*b2 + a4**2*b2)
f(672) = dtau**3*(a1**2*b1 - a2**2*b1 - a3**2*b2 + a4**2*b2)
f(673) = dtau**2*(a1**2*b1**2 + a2**2*b1**2 - a3**2*b2**2 - a4**2*b2**2)
f(674) = dtau*(-a1**2*b1**3 + a2**2*b1**3 + a3**2*b2**3 - a4**2*b2**3)
f(675) = dtau*(a1**3*b1 - a2**3*b1 - a3**3*b2 + a4**3*b2)
f(676) = dtau*(-a1**4*b1 + a2**4*b1 + a3**4*b2 - a4**4*b2)
f(677) = dtau*(-a1*b2 + a2*b2 + a3*b1 - a4*b1)
f(678) = dtau**3*(a1*b2 - a2*b2 - a3*b1 + a4*b1)
f(679) = dtau**2*(-a1*b2**2 - a2*b2**2 + a3*b1**2 + a4*b1**2)
f(680) = dtau*(a1*b2**3 - a2*b2**3 - a3*b1**3 + a4*b1**3)
f(681) = dtau*(-a1**2*b2 + a2**2*b2 + a3**2*b1 - a4**2*b1)
f(682) = dtau**3*(-a1**2*b2 + a2**2*b2 + a3**2*b1 - a4**2*b1)
f(683) = dtau**2*(-a1**2*b2**2 - a2**2*b2**2 + a3**2*b1**2 + a4**2*b1**2 &
      )
f(684) = dtau*(-a1**2*b2**3 + a2**2*b2**3 + a3**2*b1**3 - a4**2*b1**3)
f(685) = dtau*(-a1**3*b2 + a2**3*b2 + a3**3*b1 - a4**3*b1)
f(686) = dtau*(a1**4*b2 - a2**4*b2 - a3**4*b1 + a4**4*b1)
f(687) = b1*b2*dtau**2*(-b1**2 + b2**2)
f(688) = b1*b2*dtau**2*(b1**2 - b2**2)
v = sum(f*params)
end function c2h4_dipole_b1u_n3_d6_ADF


!###############################################################################


function c2h4_dipole_b1u_n4_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(983)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(983)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(r1*r2*r3 + r1*r2*r4 - r1*r3*r4 - r2*r3*r4)
f(2) = r0*(r1**2*r3*r4 - r1*r2*r3**2 - r1*r2*r4**2 + r2**2*r3*r4)
f(3) = r0*(-r1**3*r3*r4 + r1*r2*r3**3 + r1*r2*r4**3 - r2**3*r3*r4)
f(4) = r0*(-r1**2*r2*r4 - r1*r2**2*r3 + r1*r3*r4**2 + r2*r3**2*r4)
f(5) = r0*(-r1**2*r2*r4**2 + r1**2*r3*r4**2 - r1*r2**2*r3**2 + r2**2*r3 &
      **2*r4)
f(6) = r0*(-r1**3*r2*r4 - r1*r2**3*r3 + r1*r3*r4**3 + r2*r3**3*r4)
f(7) = r0*(-r1**2*r2*r3 - r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4**2)
f(8) = r0*(-r1**2*r2*r3**2 + r1**2*r3**2*r4 - r1*r2**2*r4**2 + r2**2*r3* &
      r4**2)
f(9) = r0*(-r1**2*r2**2*r3 - r1**2*r2**2*r4 + r1*r3**2*r4**2 + r2*r3**2* &
      r4**2)
f(10) = r0*(-r1**3*r2*r3 - r1*r2**3*r4 + r1*r3**3*r4 + r2*r3*r4**3)
f(11) = r0**2*(-r1*r2*r3 - r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(12) = r0**2*(r1**2*r3*r4 - r1*r2*r3**2 - r1*r2*r4**2 + r2**2*r3*r4)
f(13) = r0**2*(-r1**2*r2*r4 - r1*r2**2*r3 + r1*r3*r4**2 + r2*r3**2*r4)
f(14) = r0**2*(-r1**2*r2*r3 - r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4**2)
f(15) = r0**3*(r1*r2*r3 + r1*r2*r4 - r1*r3*r4 - r2*r3*r4)
f(16) = r0*(a1*r1*r2 + a2*r1*r2 - a3*r3*r4 - a4*r3*r4)
f(17) = r0*(a1**2*r1*r2 + a2**2*r1*r2 - a3**2*r3*r4 - a4**2*r3*r4)
f(18) = r0*(a1**3*r1*r2 + a2**3*r1*r2 - a3**3*r3*r4 - a4**3*r3*r4)
f(19) = r0*(-a1*r1*r2**2 - a2*r1**2*r2 + a3*r3*r4**2 + a4*r3**2*r4)
f(20) = r0*(-a1**2*r1*r2**2 - a2**2*r1**2*r2 + a3**2*r3*r4**2 + a4**2*r3 &
      **2*r4)
f(21) = r0*(-a1*r1*r2**3 - a2*r1**3*r2 + a3*r3*r4**3 + a4*r3**3*r4)
f(22) = r0*(-a1*r1**2*r2 - a2*r1*r2**2 + a3*r3**2*r4 + a4*r3*r4**2)
f(23) = r0*(-a1**2*r1**2*r2 - a2**2*r1*r2**2 + a3**2*r3**2*r4 + a4**2*r3 &
      *r4**2)
f(24) = r0*(a1*r1**2*r2**2 + a2*r1**2*r2**2 - a3*r3**2*r4**2 - a4*r3**2* &
      r4**2)
f(25) = r0*(-a1*r1**3*r2 - a2*r1*r2**3 + a3*r3**3*r4 + a4*r3*r4**3)
f(26) = r0**2*(-a1*r1*r2 - a2*r1*r2 + a3*r3*r4 + a4*r3*r4)
f(27) = r0**2*(-a1**2*r1*r2 - a2**2*r1*r2 + a3**2*r3*r4 + a4**2*r3*r4)
f(28) = r0**2*(a1*r1*r2**2 + a2*r1**2*r2 - a3*r3*r4**2 - a4*r3**2*r4)
f(29) = r0**2*(a1*r1**2*r2 + a2*r1*r2**2 - a3*r3**2*r4 - a4*r3*r4**2)
f(30) = r0**3*(a1*r1*r2 + a2*r1*r2 - a3*r3*r4 - a4*r3*r4)
f(31) = r0*(a1*r3*r4 + a2*r3*r4 - a3*r1*r2 - a4*r1*r2)
f(32) = r0*(-a1**2*r3*r4 - a2**2*r3*r4 + a3**2*r1*r2 + a4**2*r1*r2)
f(33) = r0*(-a1**3*r3*r4 - a2**3*r3*r4 + a3**3*r1*r2 + a4**3*r1*r2)
f(34) = r0*(-a1*r3*r4**2 - a2*r3**2*r4 + a3*r1*r2**2 + a4*r1**2*r2)
f(35) = r0*(a1**2*r3*r4**2 + a2**2*r3**2*r4 - a3**2*r1*r2**2 - a4**2*r1 &
      **2*r2)
f(36) = r0*(-a1*r3*r4**3 - a2*r3**3*r4 + a3*r1*r2**3 + a4*r1**3*r2)
f(37) = r0*(-a1*r3**2*r4 - a2*r3*r4**2 + a3*r1**2*r2 + a4*r1*r2**2)
f(38) = r0*(a1**2*r3**2*r4 + a2**2*r3*r4**2 - a3**2*r1**2*r2 - a4**2*r1* &
      r2**2)
f(39) = r0*(a1*r3**2*r4**2 + a2*r3**2*r4**2 - a3*r1**2*r2**2 - a4*r1**2* &
      r2**2)
f(40) = r0*(-a1*r3**3*r4 - a2*r3*r4**3 + a3*r1**3*r2 + a4*r1*r2**3)
f(41) = r0**2*(-a1*r3*r4 - a2*r3*r4 + a3*r1*r2 + a4*r1*r2)
f(42) = r0**2*(-a1**2*r3*r4 - a2**2*r3*r4 + a3**2*r1*r2 + a4**2*r1*r2)
f(43) = r0**2*(a1*r3*r4**2 + a2*r3**2*r4 - a3*r1*r2**2 - a4*r1**2*r2)
f(44) = r0**2*(a1*r3**2*r4 + a2*r3*r4**2 - a3*r1**2*r2 - a4*r1*r2**2)
f(45) = r0**3*(a1*r3*r4 + a2*r3*r4 - a3*r1*r2 - a4*r1*r2)
f(46) = r0*(-b1**2*r1*r2 + b2**2*r3*r4)
f(47) = r0*(b1**2*r1**2*r2 + b1**2*r1*r2**2 - b2**2*r3**2*r4 - b2**2*r3* &
      r4**2)
f(48) = r0**2*(-b1**2*r1*r2 + b2**2*r3*r4)
f(49) = r0*(-b1**2*r3*r4 + b2**2*r1*r2)
f(50) = r0*(b1**2*r3**2*r4 + b1**2*r3*r4**2 - b2**2*r1**2*r2 - b2**2*r1* &
      r2**2)
f(51) = r0**2*(-b1**2*r3*r4 + b2**2*r1*r2)
f(52) = dtau**2*r0*(r1*r2 - r3*r4)
f(53) = dtau**2*r0*(r1**2*r2 + r1*r2**2 - r3**2*r4 - r3*r4**2)
f(54) = dtau**2*r0*(-r1**2*r2 - r1*r2**2 + r3**2*r4 + r3*r4**2)
f(55) = dtau**2*r0**2*(r1*r2 - r3*r4)
f(56) = r0*(a1*r1*r3 + a2*r2*r4 - a3*r1*r3 - a4*r2*r4)
f(57) = r0*(-a1**2*r1*r3 - a2**2*r2*r4 + a3**2*r1*r3 + a4**2*r2*r4)
f(58) = r0*(a1**3*r1*r3 + a2**3*r2*r4 - a3**3*r1*r3 - a4**3*r2*r4)
f(59) = r0*(-a1*r1*r3**2 - a2*r2*r4**2 + a3*r1**2*r3 + a4*r2**2*r4)
f(60) = r0*(-a1**2*r1*r3**2 - a2**2*r2*r4**2 + a3**2*r1**2*r3 + a4**2*r2 &
      **2*r4)
f(61) = r0*(a1*r1*r3**3 + a2*r2*r4**3 - a3*r1**3*r3 - a4*r2**3*r4)
f(62) = r0*(-a1*r1**2*r3 - a2*r2**2*r4 + a3*r1*r3**2 + a4*r2*r4**2)
f(63) = r0*(-a1**2*r1**2*r3 - a2**2*r2**2*r4 + a3**2*r1*r3**2 + a4**2*r2 &
      *r4**2)
f(64) = r0*(a1*r1**2*r3**2 + a2*r2**2*r4**2 - a3*r1**2*r3**2 - a4*r2**2* &
      r4**2)
f(65) = r0*(a1*r1**3*r3 + a2*r2**3*r4 - a3*r1*r3**3 - a4*r2*r4**3)
f(66) = r0**2*(-a1*r1*r3 - a2*r2*r4 + a3*r1*r3 + a4*r2*r4)
f(67) = r0**2*(-a1**2*r1*r3 - a2**2*r2*r4 + a3**2*r1*r3 + a4**2*r2*r4)
f(68) = r0**2*(a1*r1*r3**2 + a2*r2*r4**2 - a3*r1**2*r3 - a4*r2**2*r4)
f(69) = r0**2*(a1*r1**2*r3 + a2*r2**2*r4 - a3*r1*r3**2 - a4*r2*r4**2)
f(70) = r0**3*(a1*r1*r3 + a2*r2*r4 - a3*r1*r3 - a4*r2*r4)
f(71) = r0*(a1*r2*r4 + a2*r1*r3 - a3*r2*r4 - a4*r1*r3)
f(72) = r0*(a1**2*r2*r4 + a2**2*r1*r3 - a3**2*r2*r4 - a4**2*r1*r3)
f(73) = r0*(-a1**3*r2*r4 - a2**3*r1*r3 + a3**3*r2*r4 + a4**3*r1*r3)
f(74) = r0*(-a1*r2*r4**2 - a2*r1*r3**2 + a3*r2**2*r4 + a4*r1**2*r3)
f(75) = r0*(a1**2*r2*r4**2 + a2**2*r1*r3**2 - a3**2*r2**2*r4 - a4**2*r1 &
      **2*r3)
f(76) = r0*(a1*r2*r4**3 + a2*r1*r3**3 - a3*r2**3*r4 - a4*r1**3*r3)
f(77) = r0*(-a1*r2**2*r4 - a2*r1**2*r3 + a3*r2*r4**2 + a4*r1*r3**2)
f(78) = r0*(a1**2*r2**2*r4 + a2**2*r1**2*r3 - a3**2*r2*r4**2 - a4**2*r1* &
      r3**2)
f(79) = r0*(a1*r2**2*r4**2 + a2*r1**2*r3**2 - a3*r2**2*r4**2 - a4*r1**2* &
      r3**2)
f(80) = r0*(a1*r2**3*r4 + a2*r1**3*r3 - a3*r2*r4**3 - a4*r1*r3**3)
f(81) = r0**2*(-a1*r2*r4 - a2*r1*r3 + a3*r2*r4 + a4*r1*r3)
f(82) = r0**2*(a1**2*r2*r4 + a2**2*r1*r3 - a3**2*r2*r4 - a4**2*r1*r3)
f(83) = r0**2*(a1*r2*r4**2 + a2*r1*r3**2 - a3*r2**2*r4 - a4*r1**2*r3)
f(84) = r0**2*(a1*r2**2*r4 + a2*r1**2*r3 - a3*r2*r4**2 - a4*r1*r3**2)
f(85) = r0**3*(a1*r2*r4 + a2*r1*r3 - a3*r2*r4 - a4*r1*r3)
f(86) = r0*(-b1**2*r1*r3 - b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2*r4)
f(87) = r0*(b1**2*r1*r3**2 + b1**2*r2*r4**2 - b2**2*r1**2*r3 - b2**2*r2 &
      **2*r4)
f(88) = r0*(b1**2*r1**2*r3 + b1**2*r2**2*r4 - b2**2*r1*r3**2 - b2**2*r2* &
      r4**2)
f(89) = r0**2*(b1**2*r1*r3 + b1**2*r2*r4 - b2**2*r1*r3 - b2**2*r2*r4)
f(90) = dtau**2*r0*(r1**2*r3 - r1*r3**2 + r2**2*r4 - r2*r4**2)
f(91) = dtau**2*r0*(-r1**2*r3 + r1*r3**2 - r2**2*r4 + r2*r4**2)
f(92) = r0*(-a1*r1*r4 - a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(93) = r0*(-a1**2*r1*r4 - a2**2*r2*r3 + a3**2*r2*r3 + a4**2*r1*r4)
f(94) = r0*(a1**3*r1*r4 + a2**3*r2*r3 - a3**3*r2*r3 - a4**3*r1*r4)
f(95) = r0*(a1*r1*r4**2 + a2*r2*r3**2 - a3*r2**2*r3 - a4*r1**2*r4)
f(96) = r0*(a1**2*r1*r4**2 + a2**2*r2*r3**2 - a3**2*r2**2*r3 - a4**2*r1 &
      **2*r4)
f(97) = r0*(-a1*r1*r4**3 - a2*r2*r3**3 + a3*r2**3*r3 + a4*r1**3*r4)
f(98) = r0*(a1*r1**2*r4 + a2*r2**2*r3 - a3*r2*r3**2 - a4*r1*r4**2)
f(99) = r0*(a1**2*r1**2*r4 + a2**2*r2**2*r3 - a3**2*r2*r3**2 - a4**2*r1* &
      r4**2)
f(100) = r0*(a1*r1**2*r4**2 + a2*r2**2*r3**2 - a3*r2**2*r3**2 - a4*r1**2 &
      *r4**2)
f(101) = r0*(-a1*r1**3*r4 - a2*r2**3*r3 + a3*r2*r3**3 + a4*r1*r4**3)
f(102) = r0**2*(a1*r1*r4 + a2*r2*r3 - a3*r2*r3 - a4*r1*r4)
f(103) = r0**2*(a1**2*r1*r4 + a2**2*r2*r3 - a3**2*r2*r3 - a4**2*r1*r4)
f(104) = r0**2*(-a1*r1*r4**2 - a2*r2*r3**2 + a3*r2**2*r3 + a4*r1**2*r4)
f(105) = r0**2*(-a1*r1**2*r4 - a2*r2**2*r3 + a3*r2*r3**2 + a4*r1*r4**2)
f(106) = r0**3*(-a1*r1*r4 - a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(107) = r0*(-a1*r2*r3 - a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(108) = r0*(a1**2*r2*r3 + a2**2*r1*r4 - a3**2*r1*r4 - a4**2*r2*r3)
f(109) = r0*(-a1**3*r2*r3 - a2**3*r1*r4 + a3**3*r1*r4 + a4**3*r2*r3)
f(110) = r0*(a1*r2*r3**2 + a2*r1*r4**2 - a3*r1**2*r4 - a4*r2**2*r3)
f(111) = r0*(-a1**2*r2*r3**2 - a2**2*r1*r4**2 + a3**2*r1**2*r4 + a4**2* &
      r2**2*r3)
f(112) = r0*(-a1*r2*r3**3 - a2*r1*r4**3 + a3*r1**3*r4 + a4*r2**3*r3)
f(113) = r0*(a1*r2**2*r3 + a2*r1**2*r4 - a3*r1*r4**2 - a4*r2*r3**2)
f(114) = r0*(-a1**2*r2**2*r3 - a2**2*r1**2*r4 + a3**2*r1*r4**2 + a4**2* &
      r2*r3**2)
f(115) = r0*(a1*r2**2*r3**2 + a2*r1**2*r4**2 - a3*r1**2*r4**2 - a4*r2**2 &
      *r3**2)
f(116) = r0*(-a1*r2**3*r3 - a2*r1**3*r4 + a3*r1*r4**3 + a4*r2*r3**3)
f(117) = r0**2*(a1*r2*r3 + a2*r1*r4 - a3*r1*r4 - a4*r2*r3)
f(118) = r0**2*(-a1**2*r2*r3 - a2**2*r1*r4 + a3**2*r1*r4 + a4**2*r2*r3)
f(119) = r0**2*(-a1*r2*r3**2 - a2*r1*r4**2 + a3*r1**2*r4 + a4*r2**2*r3)
f(120) = r0**2*(-a1*r2**2*r3 - a2*r1**2*r4 + a3*r1*r4**2 + a4*r2*r3**2)
f(121) = r0**3*(-a1*r2*r3 - a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(122) = r0*(b1**2*r1*r4 + b1**2*r2*r3 - b2**2*r1*r4 - b2**2*r2*r3)
f(123) = r0*(-b1**2*r1*r4**2 - b1**2*r2*r3**2 + b2**2*r1**2*r4 + b2**2* &
      r2**2*r3)
f(124) = r0*(-b1**2*r1**2*r4 - b1**2*r2**2*r3 + b2**2*r1*r4**2 + b2**2* &
      r2*r3**2)
f(125) = r0**2*(-b1**2*r1*r4 - b1**2*r2*r3 + b2**2*r1*r4 + b2**2*r2*r3)
f(126) = dtau**2*r0*(r1**2*r4 - r1*r4**2 + r2**2*r3 - r2*r3**2)
f(127) = dtau**2*r0*(-r1**2*r4 + r1*r4**2 - r2**2*r3 + r2*r3**2)
f(128) = r0*(-a1*a2*r1 - a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(129) = r0*(-a1**2*a2*r2 - a1*a2**2*r1 + a3**2*a4*r4 + a3*a4**2*r3)
f(130) = r0*(-a1**3*a2*r2 - a1*a2**3*r1 + a3**3*a4*r4 + a3*a4**3*r3)
f(131) = r0*(-a1**2*a2*r1 - a1*a2**2*r2 + a3**2*a4*r3 + a3*a4**2*r4)
f(132) = r0*(a1**2*a2**2*r1 + a1**2*a2**2*r2 - a3**2*a4**2*r3 - a3**2*a4 &
      **2*r4)
f(133) = r0*(-a1**3*a2*r1 - a1*a2**3*r2 + a3**3*a4*r3 + a3*a4**3*r4)
f(134) = r0*(-a1*a2*r1**2 - a1*a2*r2**2 + a3*a4*r3**2 + a3*a4*r4**2)
f(135) = r0*(a1**2*a2*r2**2 + a1*a2**2*r1**2 - a3**2*a4*r4**2 - a3*a4**2 &
      *r3**2)
f(136) = r0*(-a1**2*a2*r1**2 - a1*a2**2*r2**2 + a3**2*a4*r3**2 + a3*a4** &
      2*r4**2)
f(137) = r0*(a1*a2*r1**3 + a1*a2*r2**3 - a3*a4*r3**3 - a3*a4*r4**3)
f(138) = r0**2*(-a1*a2*r1 - a1*a2*r2 + a3*a4*r3 + a3*a4*r4)
f(139) = r0**2*(a1**2*a2*r2 + a1*a2**2*r1 - a3**2*a4*r4 - a3*a4**2*r3)
f(140) = r0**2*(-a1**2*a2*r1 - a1*a2**2*r2 + a3**2*a4*r3 + a3*a4**2*r4)
f(141) = r0**2*(a1*a2*r1**2 + a1*a2*r2**2 - a3*a4*r3**2 - a3*a4*r4**2)
f(142) = r0**3*(a1*a2*r1 + a1*a2*r2 - a3*a4*r3 - a3*a4*r4)
f(143) = r0*(-a1*a3*r1 + a1*a3*r3 - a2*a4*r2 + a2*a4*r4)
f(144) = r0*(a1**2*a3*r3 - a1*a3**2*r1 + a2**2*a4*r4 - a2*a4**2*r2)
f(145) = r0*(-a1**3*a3*r3 + a1*a3**3*r1 - a2**3*a4*r4 + a2*a4**3*r2)
f(146) = r0*(a1**2*a3*r1 - a1*a3**2*r3 + a2**2*a4*r2 - a2*a4**2*r4)
f(147) = r0*(a1**2*a3**2*r1 - a1**2*a3**2*r3 + a2**2*a4**2*r2 - a2**2*a4 &
      **2*r4)
f(148) = r0*(a1**3*a3*r1 - a1*a3**3*r3 + a2**3*a4*r2 - a2*a4**3*r4)
f(149) = r0*(a1*a3*r1**2 - a1*a3*r3**2 + a2*a4*r2**2 - a2*a4*r4**2)
f(150) = r0*(-a1**2*a3*r3**2 + a1*a3**2*r1**2 - a2**2*a4*r4**2 + a2*a4** &
      2*r2**2)
f(151) = r0*(a1**2*a3*r1**2 - a1*a3**2*r3**2 + a2**2*a4*r2**2 - a2*a4**2 &
      *r4**2)
f(152) = r0*(a1*a3*r1**3 - a1*a3*r3**3 + a2*a4*r2**3 - a2*a4*r4**3)
f(153) = r0**2*(a1*a3*r1 - a1*a3*r3 + a2*a4*r2 - a2*a4*r4)
f(154) = r0**2*(-a1**2*a3*r3 + a1*a3**2*r1 - a2**2*a4*r4 + a2*a4**2*r2)
f(155) = r0**2*(a1**2*a3*r1 - a1*a3**2*r3 + a2**2*a4*r2 - a2*a4**2*r4)
f(156) = r0**2*(a1*a3*r1**2 - a1*a3*r3**2 + a2*a4*r2**2 - a2*a4*r4**2)
f(157) = r0**3*(a1*a3*r1 - a1*a3*r3 + a2*a4*r2 - a2*a4*r4)
f(158) = r0*(-a1*a4*r1 + a1*a4*r4 - a2*a3*r2 + a2*a3*r3)
f(159) = r0*(-a1**2*a4*r4 + a1*a4**2*r1 - a2**2*a3*r3 + a2*a3**2*r2)
f(160) = r0*(a1**3*a4*r4 - a1*a4**3*r1 + a2**3*a3*r3 - a2*a3**3*r2)
f(161) = r0*(a1**2*a4*r1 - a1*a4**2*r4 + a2**2*a3*r2 - a2*a3**2*r3)
f(162) = r0*(a1**2*a4**2*r1 - a1**2*a4**2*r4 + a2**2*a3**2*r2 - a2**2*a3 &
      **2*r3)
f(163) = r0*(a1**3*a4*r1 - a1*a4**3*r4 + a2**3*a3*r2 - a2*a3**3*r3)
f(164) = r0*(a1*a4*r1**2 - a1*a4*r4**2 + a2*a3*r2**2 - a2*a3*r3**2)
f(165) = r0*(a1**2*a4*r4**2 - a1*a4**2*r1**2 + a2**2*a3*r3**2 - a2*a3**2 &
      *r2**2)
f(166) = r0*(a1**2*a4*r1**2 - a1*a4**2*r4**2 + a2**2*a3*r2**2 - a2*a3**2 &
      *r3**2)
f(167) = r0*(a1*a4*r1**3 - a1*a4*r4**3 + a2*a3*r2**3 - a2*a3*r3**3)
f(168) = r0**2*(a1*a4*r1 - a1*a4*r4 + a2*a3*r2 - a2*a3*r3)
f(169) = r0**2*(a1**2*a4*r4 - a1*a4**2*r1 + a2**2*a3*r3 - a2*a3**2*r2)
f(170) = r0**2*(a1**2*a4*r1 - a1*a4**2*r4 + a2**2*a3*r2 - a2*a3**2*r3)
f(171) = r0**2*(a1*a4*r1**2 - a1*a4*r4**2 + a2*a3*r2**2 - a2*a3*r3**2)
f(172) = r0**3*(-a1*a4*r1 + a1*a4*r4 - a2*a3*r2 + a2*a3*r3)
f(173) = r0*(a1*b1**2*r1 + a2*b1**2*r2 - a3*b2**2*r3 - a4*b2**2*r4)
f(174) = r0*(-a1**2*b1**2*r1 - a2**2*b1**2*r2 + a3**2*b2**2*r3 + a4**2* &
      b2**2*r4)
f(175) = r0*(a1*b1**2*r1**2 + a2*b1**2*r2**2 - a3*b2**2*r3**2 - a4*b2**2 &
      *r4**2)
f(176) = r0**2*(a1*b1**2*r1 + a2*b1**2*r2 - a3*b2**2*r3 - a4*b2**2*r4)
f(177) = r0*(-a1*b2**2*r1 - a2*b2**2*r2 + a3*b1**2*r3 + a4*b1**2*r4)
f(178) = r0*(a1**2*b2**2*r1 + a2**2*b2**2*r2 - a3**2*b1**2*r3 - a4**2*b1 &
      **2*r4)
f(179) = r0*(a1*b2**2*r1**2 + a2*b2**2*r2**2 - a3*b1**2*r3**2 - a4*b1**2 &
      *r4**2)
f(180) = r0**2*(a1*b2**2*r1 + a2*b2**2*r2 - a3*b1**2*r3 - a4*b1**2*r4)
f(181) = dtau**2*r0*(-a1*r1 - a2*r2 + a3*r3 + a4*r4)
f(182) = dtau**2*r0*(a1**2*r1 + a2**2*r2 - a3**2*r3 - a4**2*r4)
f(183) = dtau**2*r0*(a1*r1**2 + a2*r2**2 - a3*r3**2 - a4*r4**2)
f(184) = dtau**2*r0**2*(-a1*r1 - a2*r2 + a3*r3 + a4*r4)
f(185) = r0*(-a1*a4*r2 + a1*a4*r3 - a2*a3*r1 + a2*a3*r4)
f(186) = r0*(a1**2*a4*r3 - a1*a4**2*r2 + a2**2*a3*r4 - a2*a3**2*r1)
f(187) = r0*(-a1**3*a4*r3 + a1*a4**3*r2 - a2**3*a3*r4 + a2*a3**3*r1)
f(188) = r0*(a1**2*a4*r2 - a1*a4**2*r3 + a2**2*a3*r1 - a2*a3**2*r4)
f(189) = r0*(a1**2*a4**2*r2 - a1**2*a4**2*r3 + a2**2*a3**2*r1 - a2**2*a3 &
      **2*r4)
f(190) = r0*(a1**3*a4*r2 - a1*a4**3*r3 + a2**3*a3*r1 - a2*a3**3*r4)
f(191) = r0*(a1*a4*r2**2 - a1*a4*r3**2 + a2*a3*r1**2 - a2*a3*r4**2)
f(192) = r0*(-a1**2*a4*r3**2 + a1*a4**2*r2**2 - a2**2*a3*r4**2 + a2*a3** &
      2*r1**2)
f(193) = r0*(-a1**2*a4*r2**2 + a1*a4**2*r3**2 - a2**2*a3*r1**2 + a2*a3** &
      2*r4**2)
f(194) = r0*(-a1*a4*r2**3 + a1*a4*r3**3 - a2*a3*r1**3 + a2*a3*r4**3)
f(195) = r0**2*(a1*a4*r2 - a1*a4*r3 + a2*a3*r1 - a2*a3*r4)
f(196) = r0**2*(-a1**2*a4*r3 + a1*a4**2*r2 - a2**2*a3*r4 + a2*a3**2*r1)
f(197) = r0**2*(-a1**2*a4*r2 + a1*a4**2*r3 - a2**2*a3*r1 + a2*a3**2*r4)
f(198) = r0**2*(a1*a4*r2**2 - a1*a4*r3**2 + a2*a3*r1**2 - a2*a3*r4**2)
f(199) = r0**3*(a1*a4*r2 - a1*a4*r3 + a2*a3*r1 - a2*a3*r4)
f(200) = r0*(-a1*a3*r2 + a1*a3*r4 - a2*a4*r1 + a2*a4*r3)
f(201) = r0*(-a1**2*a3*r4 + a1*a3**2*r2 - a2**2*a4*r3 + a2*a4**2*r1)
f(202) = r0*(a1**3*a3*r4 - a1*a3**3*r2 + a2**3*a4*r3 - a2*a4**3*r1)
f(203) = r0*(a1**2*a3*r2 - a1*a3**2*r4 + a2**2*a4*r1 - a2*a4**2*r3)
f(204) = r0*(a1**2*a3**2*r2 - a1**2*a3**2*r4 + a2**2*a4**2*r1 - a2**2*a4 &
      **2*r3)
f(205) = r0*(a1**3*a3*r2 - a1*a3**3*r4 + a2**3*a4*r1 - a2*a4**3*r3)
f(206) = r0*(a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 - a2*a4*r3**2)
f(207) = r0*(a1**2*a3*r4**2 - a1*a3**2*r2**2 + a2**2*a4*r3**2 - a2*a4**2 &
      *r1**2)
f(208) = r0*(-a1**2*a3*r2**2 + a1*a3**2*r4**2 - a2**2*a4*r1**2 + a2*a4** &
      2*r3**2)
f(209) = r0*(-a1*a3*r2**3 + a1*a3*r4**3 - a2*a4*r1**3 + a2*a4*r3**3)
f(210) = r0**2*(a1*a3*r2 - a1*a3*r4 + a2*a4*r1 - a2*a4*r3)
f(211) = r0**2*(a1**2*a3*r4 - a1*a3**2*r2 + a2**2*a4*r3 - a2*a4**2*r1)
f(212) = r0**2*(-a1**2*a3*r2 + a1*a3**2*r4 - a2**2*a4*r1 + a2*a4**2*r3)
f(213) = r0**2*(a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 - a2*a4*r3**2)
f(214) = r0**3*(-a1*a3*r2 + a1*a3*r4 - a2*a4*r1 + a2*a4*r3)
f(215) = r0*(a1*b1**2*r2 + a2*b1**2*r1 - a3*b2**2*r4 - a4*b2**2*r3)
f(216) = r0*(a1**2*b1**2*r2 + a2**2*b1**2*r1 - a3**2*b2**2*r4 - a4**2*b2 &
      **2*r3)
f(217) = r0*(a1*b1**2*r2**2 + a2*b1**2*r1**2 - a3*b2**2*r4**2 - a4*b2**2 &
      *r3**2)
f(218) = r0**2*(a1*b1**2*r2 + a2*b1**2*r1 - a3*b2**2*r4 - a4*b2**2*r3)
f(219) = r0*(a1*b2**2*r2 + a2*b2**2*r1 - a3*b1**2*r4 - a4*b1**2*r3)
f(220) = r0*(-a1**2*b2**2*r2 - a2**2*b2**2*r1 + a3**2*b1**2*r4 + a4**2* &
      b1**2*r3)
f(221) = r0*(a1*b2**2*r2**2 + a2*b2**2*r1**2 - a3*b1**2*r4**2 - a4*b1**2 &
      *r3**2)
f(222) = r0**2*(a1*b2**2*r2 + a2*b2**2*r1 - a3*b1**2*r4 - a4*b1**2*r3)
f(223) = dtau**2*r0*(-a1*r2 - a2*r1 + a3*r4 + a4*r3)
f(224) = dtau**2*r0*(-a1**2*r2 - a2**2*r1 + a3**2*r4 + a4**2*r3)
f(225) = dtau**2*r0*(a1*r2**2 + a2*r1**2 - a3*r4**2 - a4*r3**2)
f(226) = dtau**2*r0**2*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(227) = r0*(-a1*a2*r3 - a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(228) = r0*(a1**2*a2*r4 + a1*a2**2*r3 - a3**2*a4*r2 - a3*a4**2*r1)
f(229) = r0*(-a1**3*a2*r4 - a1*a2**3*r3 + a3**3*a4*r2 + a3*a4**3*r1)
f(230) = r0*(-a1**2*a2*r3 - a1*a2**2*r4 + a3**2*a4*r1 + a3*a4**2*r2)
f(231) = r0*(-a1**2*a2**2*r3 - a1**2*a2**2*r4 + a3**2*a4**2*r1 + a3**2* &
      a4**2*r2)
f(232) = r0*(a1**3*a2*r3 + a1*a2**3*r4 - a3**3*a4*r1 - a3*a4**3*r2)
f(233) = r0*(-a1*a2*r3**2 - a1*a2*r4**2 + a3*a4*r1**2 + a3*a4*r2**2)
f(234) = r0*(-a1**2*a2*r4**2 - a1*a2**2*r3**2 + a3**2*a4*r2**2 + a3*a4** &
      2*r1**2)
f(235) = r0*(a1**2*a2*r3**2 + a1*a2**2*r4**2 - a3**2*a4*r1**2 - a3*a4**2 &
      *r2**2)
f(236) = r0*(a1*a2*r3**3 + a1*a2*r4**3 - a3*a4*r1**3 - a3*a4*r2**3)
f(237) = r0**2*(-a1*a2*r3 - a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(238) = r0**2*(-a1**2*a2*r4 - a1*a2**2*r3 + a3**2*a4*r2 + a3*a4**2*r1)
f(239) = r0**2*(a1**2*a2*r3 + a1*a2**2*r4 - a3**2*a4*r1 - a3*a4**2*r2)
f(240) = r0**2*(-a1*a2*r3**2 - a1*a2*r4**2 + a3*a4*r1**2 + a3*a4*r2**2)
f(241) = r0**3*(a1*a2*r3 + a1*a2*r4 - a3*a4*r1 - a3*a4*r2)
f(242) = r0*(-a1*b2**2*r3 - a2*b2**2*r4 + a3*b1**2*r1 + a4*b1**2*r2)
f(243) = r0*(-a1**2*b2**2*r3 - a2**2*b2**2*r4 + a3**2*b1**2*r1 + a4**2* &
      b1**2*r2)
f(244) = r0*(-a1*b2**2*r3**2 - a2*b2**2*r4**2 + a3*b1**2*r1**2 + a4*b1** &
      2*r2**2)
f(245) = r0**2*(-a1*b2**2*r3 - a2*b2**2*r4 + a3*b1**2*r1 + a4*b1**2*r2)
f(246) = r0*(-a1*b1**2*r3 - a2*b1**2*r4 + a3*b2**2*r1 + a4*b2**2*r2)
f(247) = r0*(a1**2*b1**2*r3 + a2**2*b1**2*r4 - a3**2*b2**2*r1 - a4**2*b2 &
      **2*r2)
f(248) = r0*(a1*b1**2*r3**2 + a2*b1**2*r4**2 - a3*b2**2*r1**2 - a4*b2**2 &
      *r2**2)
f(249) = r0**2*(a1*b1**2*r3 + a2*b1**2*r4 - a3*b2**2*r1 - a4*b2**2*r2)
f(250) = dtau**2*r0*(a1*r3 + a2*r4 - a3*r1 - a4*r2)
f(251) = dtau**2*r0*(a1**2*r3 + a2**2*r4 - a3**2*r1 - a4**2*r2)
f(252) = dtau**2*r0*(-a1*r3**2 - a2*r4**2 + a3*r1**2 + a4*r2**2)
f(253) = dtau**2*r0**2*(-a1*r3 - a2*r4 + a3*r1 + a4*r2)
f(254) = r0*(a1*b2**2*r4 + a2*b2**2*r3 - a3*b1**2*r2 - a4*b1**2*r1)
f(255) = r0*(a1**2*b2**2*r4 + a2**2*b2**2*r3 - a3**2*b1**2*r2 - a4**2*b1 &
      **2*r1)
f(256) = r0*(a1*b2**2*r4**2 + a2*b2**2*r3**2 - a3*b1**2*r2**2 - a4*b1**2 &
      *r1**2)
f(257) = r0**2*(a1*b2**2*r4 + a2*b2**2*r3 - a3*b1**2*r2 - a4*b1**2*r1)
f(258) = r0*(-a1*b1**2*r4 - a2*b1**2*r3 + a3*b2**2*r2 + a4*b2**2*r1)
f(259) = r0*(-a1**2*b1**2*r4 - a2**2*b1**2*r3 + a3**2*b2**2*r2 + a4**2* &
      b2**2*r1)
f(260) = r0*(a1*b1**2*r4**2 + a2*b1**2*r3**2 - a3*b2**2*r2**2 - a4*b2**2 &
      *r1**2)
f(261) = r0**2*(a1*b1**2*r4 + a2*b1**2*r3 - a3*b2**2*r2 - a4*b2**2*r1)
f(262) = dtau**2*r0*(a1*r4 + a2*r3 - a3*r2 - a4*r1)
f(263) = dtau**2*r0*(-a1**2*r4 - a2**2*r3 + a3**2*r2 + a4**2*r1)
f(264) = dtau**2*r0*(-a1*r4**2 - a2*r3**2 + a3*r2**2 + a4*r1**2)
f(265) = dtau**2*r0**2*(-a1*r4 - a2*r3 + a3*r2 + a4*r1)
f(266) = b1*b2*r0*(-r1 - r2 + r3 + r4)
f(267) = b1*b2*r0*(-b1**2*r3 - b1**2*r4 + b2**2*r1 + b2**2*r2)
f(268) = b1**2*b2**2*r0*(-r1 - r2 + r3 + r4)
f(269) = b1*b2*r0*(b1**2*r1 + b1**2*r2 - b2**2*r3 - b2**2*r4)
f(270) = b1*b2*r0*(r1**2 + r2**2 - r3**2 - r4**2)
f(271) = b1*b2*r0*(-r1**3 - r2**3 + r3**3 + r4**3)
f(272) = b1*b2*r0**2*(r1 + r2 - r3 - r4)
f(273) = b1*b2*r0**2*(r1**2 + r2**2 - r3**2 - r4**2)
f(274) = b1*b2*r0**3*(r1 + r2 - r3 - r4)
f(275) = dtau*r0*(-b1*r1 + b1*r2 + b2*r3 - b2*r4)
f(276) = dtau**3*r0*(-b1*r1 + b1*r2 + b2*r3 - b2*r4)
f(277) = dtau**2*r0*(b1**2*r1 + b1**2*r2 - b2**2*r3 - b2**2*r4)
f(278) = dtau*r0*(-b1**3*r1 + b1**3*r2 + b2**3*r3 - b2**3*r4)
f(279) = dtau*r0*(b1*r1**2 - b1*r2**2 - b2*r3**2 + b2*r4**2)
f(280) = dtau*r0*(-b1*r1**3 + b1*r2**3 + b2*r3**3 - b2*r4**3)
f(281) = dtau*r0**2*(b1*r1 - b1*r2 - b2*r3 + b2*r4)
f(282) = dtau*r0**2*(-b1*r1**2 + b1*r2**2 + b2*r3**2 - b2*r4**2)
f(283) = dtau*r0**3*(-b1*r1 + b1*r2 + b2*r3 - b2*r4)
f(284) = dtau*r0*(b1*r3 - b1*r4 - b2*r1 + b2*r2)
f(285) = dtau**3*r0*(b1*r3 - b1*r4 - b2*r1 + b2*r2)
f(286) = dtau**2*r0*(-b1**2*r3 - b1**2*r4 + b2**2*r1 + b2**2*r2)
f(287) = dtau*r0*(-b1**3*r3 + b1**3*r4 + b2**3*r1 - b2**3*r2)
f(288) = dtau*r0*(-b1*r3**2 + b1*r4**2 + b2*r1**2 - b2*r2**2)
f(289) = dtau*r0*(-b1*r3**3 + b1*r4**3 + b2*r1**3 - b2*r2**3)
f(290) = dtau*r0**2*(-b1*r3 + b1*r4 + b2*r1 - b2*r2)
f(291) = dtau*r0**2*(b1*r3**2 - b1*r4**2 - b2*r1**2 + b2*r2**2)
f(292) = dtau*r0**3*(b1*r3 - b1*r4 - b2*r1 + b2*r2)
f(293) = r0*(a1*a2*a3 + a1*a2*a4 - a1*a3*a4 - a2*a3*a4)
f(294) = r0*(-a1**2*a3*a4 + a1*a2*a3**2 + a1*a2*a4**2 - a2**2*a3*a4)
f(295) = r0*(-a1**3*a3*a4 + a1*a2*a3**3 + a1*a2*a4**3 - a2**3*a3*a4)
f(296) = r0*(a1**2*a2*a4 + a1*a2**2*a3 - a1*a3*a4**2 - a2*a3**2*a4)
f(297) = r0*(-a1**2*a2*a4**2 + a1**2*a3*a4**2 - a1*a2**2*a3**2 + a2**2* &
      a3**2*a4)
f(298) = r0*(-a1**3*a2*a4 - a1*a2**3*a3 + a1*a3*a4**3 + a2*a3**3*a4)
f(299) = r0*(a1**2*a2*a3 + a1*a2**2*a4 - a1*a3**2*a4 - a2*a3*a4**2)
f(300) = r0*(-a1**2*a2*a3**2 + a1**2*a3**2*a4 - a1*a2**2*a4**2 + a2**2* &
      a3*a4**2)
f(301) = r0*(-a1**2*a2**2*a3 - a1**2*a2**2*a4 + a1*a3**2*a4**2 + a2*a3** &
      2*a4**2)
f(302) = r0*(-a1**3*a2*a3 - a1*a2**3*a4 + a1*a3**3*a4 + a2*a3*a4**3)
f(303) = r0**2*(a1*a2*a3 + a1*a2*a4 - a1*a3*a4 - a2*a3*a4)
f(304) = r0**2*(a1**2*a3*a4 - a1*a2*a3**2 - a1*a2*a4**2 + a2**2*a3*a4)
f(305) = r0**2*(-a1**2*a2*a4 - a1*a2**2*a3 + a1*a3*a4**2 + a2*a3**2*a4)
f(306) = r0**2*(-a1**2*a2*a3 - a1*a2**2*a4 + a1*a3**2*a4 + a2*a3*a4**2)
f(307) = r0**3*(-a1*a2*a3 - a1*a2*a4 + a1*a3*a4 + a2*a3*a4)
f(308) = r0*(-a1*a2*b1**2 + a3*a4*b2**2)
f(309) = r0*(-a1**2*a2*b1**2 - a1*a2**2*b1**2 + a3**2*a4*b2**2 + a3*a4** &
      2*b2**2)
f(310) = r0**2*(-a1*a2*b1**2 + a3*a4*b2**2)
f(311) = r0*(a1*a2*b2**2 - a3*a4*b1**2)
f(312) = r0*(-a1**2*a2*b2**2 - a1*a2**2*b2**2 + a3**2*a4*b1**2 + a3*a4** &
      2*b1**2)
f(313) = r0**2*(-a1*a2*b2**2 + a3*a4*b1**2)
f(314) = dtau**2*r0*(a1*a2 - a3*a4)
f(315) = dtau**2*r0*(-a1**2*a2 - a1*a2**2 + a3**2*a4 + a3*a4**2)
f(316) = dtau**2*r0**2*(a1*a2 - a3*a4)
f(317) = r0*(a1*a3*b1**2 - a1*a3*b2**2 + a2*a4*b1**2 - a2*a4*b2**2)
f(318) = r0*(-a1**2*a3*b2**2 + a1*a3**2*b1**2 - a2**2*a4*b2**2 + a2*a4** &
      2*b1**2)
f(319) = r0*(-a1**2*a3*b1**2 + a1*a3**2*b2**2 - a2**2*a4*b1**2 + a2*a4** &
      2*b2**2)
f(320) = r0**2*(-a1*a3*b1**2 + a1*a3*b2**2 - a2*a4*b1**2 + a2*a4*b2**2)
f(321) = dtau**2*r0*(a1**2*a3 - a1*a3**2 + a2**2*a4 - a2*a4**2)
f(322) = r0*(a1*a4*b1**2 - a1*a4*b2**2 + a2*a3*b1**2 - a2*a3*b2**2)
f(323) = r0*(-a1**2*a4*b2**2 + a1*a4**2*b1**2 - a2**2*a3*b2**2 + a2*a3** &
      2*b1**2)
f(324) = r0*(a1**2*a4*b1**2 - a1*a4**2*b2**2 + a2**2*a3*b1**2 - a2*a3**2 &
      *b2**2)
f(325) = r0**2*(a1*a4*b1**2 - a1*a4*b2**2 + a2*a3*b1**2 - a2*a3*b2**2)
f(326) = dtau**2*r0*(a1**2*a4 - a1*a4**2 + a2**2*a3 - a2*a3**2)
f(327) = b1*b2*r0*(-a1 - a2 + a3 + a4)
f(328) = b1*b2*r0*(-a1*b2**2 - a2*b2**2 + a3*b1**2 + a4*b1**2)
f(329) = b1**2*b2**2*r0*(a1 + a2 - a3 - a4)
f(330) = b1*b2*r0*(-a1*b1**2 - a2*b1**2 + a3*b2**2 + a4*b2**2)
f(331) = b1*b2*r0*(-a1**2 - a2**2 + a3**2 + a4**2)
f(332) = b1*b2*r0*(-a1**3 - a2**3 + a3**3 + a4**3)
f(333) = b1*b2*r0**2*(-a1 - a2 + a3 + a4)
f(334) = b1*b2*r0**2*(-a1**2 - a2**2 + a3**2 + a4**2)
f(335) = b1*b2*r0**3*(a1 + a2 - a3 - a4)
f(336) = dtau*r0*(-a1*b1 + a2*b1 + a3*b2 - a4*b2)
f(337) = dtau**3*r0*(a1*b1 - a2*b1 - a3*b2 + a4*b2)
f(338) = dtau**2*r0*(-a1*b1**2 - a2*b1**2 + a3*b2**2 + a4*b2**2)
f(339) = dtau*r0*(a1*b1**3 - a2*b1**3 - a3*b2**3 + a4*b2**3)
f(340) = dtau*r0*(a1**2*b1 - a2**2*b1 - a3**2*b2 + a4**2*b2)
f(341) = dtau*r0*(-a1**3*b1 + a2**3*b1 + a3**3*b2 - a4**3*b2)
f(342) = dtau*r0**2*(-a1*b1 + a2*b1 + a3*b2 - a4*b2)
f(343) = dtau*r0**2*(-a1**2*b1 + a2**2*b1 + a3**2*b2 - a4**2*b2)
f(344) = dtau*r0**3*(a1*b1 - a2*b1 - a3*b2 + a4*b2)
f(345) = dtau*r0*(-a1*b2 + a2*b2 + a3*b1 - a4*b1)
f(346) = dtau**3*r0*(a1*b2 - a2*b2 - a3*b1 + a4*b1)
f(347) = dtau**2*r0*(-a1*b2**2 - a2*b2**2 + a3*b1**2 + a4*b1**2)
f(348) = dtau*r0*(a1*b2**3 - a2*b2**3 - a3*b1**3 + a4*b1**3)
f(349) = dtau*r0*(a1**2*b2 - a2**2*b2 - a3**2*b1 + a4**2*b1)
f(350) = dtau*r0*(-a1**3*b2 + a2**3*b2 + a3**3*b1 - a4**3*b1)
f(351) = dtau*r0**2*(a1*b2 - a2*b2 - a3*b1 + a4*b1)
f(352) = dtau*r0**2*(a1**2*b2 - a2**2*b2 - a3**2*b1 + a4**2*b1)
f(353) = dtau*r0**3*(-a1*b2 + a2*b2 + a3*b1 - a4*b1)
f(354) = r1*r2*r3*r4*(r1 + r2 - r3 - r4)
f(355) = r1*r2*r3*r4*(r1**2 + r2**2 - r3**2 - r4**2)
f(356) = r1*r2*r3*r4*(r1*r2 - r3*r4)
f(357) = a1*r1*r2*r3 + a2*r1*r2*r4 - a3*r1*r3*r4 - a4*r2*r3*r4
f(358) = -a1**2*r1*r2*r3 - a2**2*r1*r2*r4 + a3**2*r1*r3*r4 + a4**2*r2*r3 &
      *r4
f(359) = a1**3*r1*r2*r3 + a2**3*r1*r2*r4 - a3**3*r1*r3*r4 - a4**3*r2*r3* &
      r4
f(360) = a1*r1*r2*r3**2 + a2*r1*r2*r4**2 - a3*r1**2*r3*r4 - a4*r2**2*r3* &
      r4
f(361) = -a1**2*r1*r2*r3**2 - a2**2*r1*r2*r4**2 + a3**2*r1**2*r3*r4 + a4 &
      **2*r2**2*r3*r4
f(362) = -a1*r1*r2*r3**3 - a2*r1*r2*r4**3 + a3*r1**3*r3*r4 + a4*r2**3*r3 &
      *r4
f(363) = a1*r1*r2**2*r3 + a2*r1**2*r2*r4 - a3*r1*r3*r4**2 - a4*r2*r3**2* &
      r4
f(364) = a1**2*r1*r2**2*r3 + a2**2*r1**2*r2*r4 - a3**2*r1*r3*r4**2 - a4 &
      **2*r2*r3**2*r4
f(365) = -a1*r1*r2**2*r3**2 - a2*r1**2*r2*r4**2 + a3*r1**2*r3*r4**2 + a4 &
      *r2**2*r3**2*r4
f(366) = a1*r1*r2**3*r3 + a2*r1**3*r2*r4 - a3*r1*r3*r4**3 - a4*r2*r3**3* &
      r4
f(367) = a1*r1**2*r2*r3 + a2*r1*r2**2*r4 - a3*r1*r3**2*r4 - a4*r2*r3*r4 &
      **2
f(368) = a1**2*r1**2*r2*r3 + a2**2*r1*r2**2*r4 - a3**2*r1*r3**2*r4 - a4 &
      **2*r2*r3*r4**2
f(369) = a1*r1**2*r2*r3**2 + a2*r1*r2**2*r4**2 - a3*r1**2*r3**2*r4 - a4* &
      r2**2*r3*r4**2
f(370) = a1*r1**2*r2**2*r3 + a2*r1**2*r2**2*r4 - a3*r1*r3**2*r4**2 - a4* &
      r2*r3**2*r4**2
f(371) = a1*r1**3*r2*r3 + a2*r1*r2**3*r4 - a3*r1*r3**3*r4 - a4*r2*r3*r4 &
      **3
f(372) = -a1*r1*r2*r4 - a2*r1*r2*r3 + a3*r2*r3*r4 + a4*r1*r3*r4
f(373) = a1**2*r1*r2*r4 + a2**2*r1*r2*r3 - a3**2*r2*r3*r4 - a4**2*r1*r3* &
      r4
f(374) = -a1**3*r1*r2*r4 - a2**3*r1*r2*r3 + a3**3*r2*r3*r4 + a4**3*r1*r3 &
      *r4
f(375) = a1*r1*r2*r4**2 + a2*r1*r2*r3**2 - a3*r2**2*r3*r4 - a4*r1**2*r3* &
      r4
f(376) = a1**2*r1*r2*r4**2 + a2**2*r1*r2*r3**2 - a3**2*r2**2*r3*r4 - a4 &
      **2*r1**2*r3*r4
f(377) = -a1*r1*r2*r4**3 - a2*r1*r2*r3**3 + a3*r2**3*r3*r4 + a4*r1**3*r3 &
      *r4
f(378) = a1*r1**2*r2*r4 + a2*r1*r2**2*r3 - a3*r2*r3**2*r4 - a4*r1*r3*r4 &
      **2
f(379) = -a1**2*r1**2*r2*r4 - a2**2*r1*r2**2*r3 + a3**2*r2*r3**2*r4 + a4 &
      **2*r1*r3*r4**2
f(380) = -a1*r1**2*r2*r4**2 - a2*r1*r2**2*r3**2 + a3*r2**2*r3**2*r4 + a4 &
      *r1**2*r3*r4**2
f(381) = a1*r1**3*r2*r4 + a2*r1*r2**3*r3 - a3*r2*r3**3*r4 - a4*r1*r3*r4 &
      **3
f(382) = a1*r1*r2**2*r4 + a2*r1**2*r2*r3 - a3*r2*r3*r4**2 - a4*r1*r3**2* &
      r4
f(383) = -a1**2*r1*r2**2*r4 - a2**2*r1**2*r2*r3 + a3**2*r2*r3*r4**2 + a4 &
      **2*r1*r3**2*r4
f(384) = a1*r1*r2**2*r4**2 + a2*r1**2*r2*r3**2 - a3*r2**2*r3*r4**2 - a4* &
      r1**2*r3**2*r4
f(385) = -a1*r1**2*r2**2*r4 - a2*r1**2*r2**2*r3 + a3*r2*r3**2*r4**2 + a4 &
      *r1*r3**2*r4**2
f(386) = a1*r1*r2**3*r4 + a2*r1**3*r2*r3 - a3*r2*r3*r4**3 - a4*r1*r3**3* &
      r4
f(387) = a1*r1*r3*r4 + a2*r2*r3*r4 - a3*r1*r2*r3 - a4*r1*r2*r4
f(388) = a1**2*r1*r3*r4 + a2**2*r2*r3*r4 - a3**2*r1*r2*r3 - a4**2*r1*r2* &
      r4
f(389) = -a1**3*r1*r3*r4 - a2**3*r2*r3*r4 + a3**3*r1*r2*r3 + a4**3*r1*r2 &
      *r4
f(390) = a1*r1**2*r3*r4 + a2*r2**2*r3*r4 - a3*r1*r2*r3**2 - a4*r1*r2*r4 &
      **2
f(391) = -a1**2*r1**2*r3*r4 - a2**2*r2**2*r3*r4 + a3**2*r1*r2*r3**2 + a4 &
      **2*r1*r2*r4**2
f(392) = -a1*r1**3*r3*r4 - a2*r2**3*r3*r4 + a3*r1*r2*r3**3 + a4*r1*r2*r4 &
      **3
f(393) = a1*r1*r3*r4**2 + a2*r2*r3**2*r4 - a3*r1*r2**2*r3 - a4*r1**2*r2* &
      r4
f(394) = -a1**2*r1*r3*r4**2 - a2**2*r2*r3**2*r4 + a3**2*r1*r2**2*r3 + a4 &
      **2*r1**2*r2*r4
f(395) = a1*r1**2*r3*r4**2 + a2*r2**2*r3**2*r4 - a3*r1*r2**2*r3**2 - a4* &
      r1**2*r2*r4**2
f(396) = -a1*r1*r3*r4**3 - a2*r2*r3**3*r4 + a3*r1*r2**3*r3 + a4*r1**3*r2 &
      *r4
f(397) = a1*r1*r3**2*r4 + a2*r2*r3*r4**2 - a3*r1**2*r2*r3 - a4*r1*r2**2* &
      r4
f(398) = a1**2*r1*r3**2*r4 + a2**2*r2*r3*r4**2 - a3**2*r1**2*r2*r3 - a4 &
      **2*r1*r2**2*r4
f(399) = -a1*r1**2*r3**2*r4 - a2*r2**2*r3*r4**2 + a3*r1**2*r2*r3**2 + a4 &
      *r1*r2**2*r4**2
f(400) = -a1*r1*r3**2*r4**2 - a2*r2*r3**2*r4**2 + a3*r1**2*r2**2*r3 + a4 &
      *r1**2*r2**2*r4
f(401) = a1*r1*r3**3*r4 + a2*r2*r3*r4**3 - a3*r1**3*r2*r3 - a4*r1*r2**3* &
      r4
f(402) = a1*r2*r3*r4 + a2*r1*r3*r4 - a3*r1*r2*r4 - a4*r1*r2*r3
f(403) = -a1**2*r2*r3*r4 - a2**2*r1*r3*r4 + a3**2*r1*r2*r4 + a4**2*r1*r2 &
      *r3
f(404) = a1**3*r2*r3*r4 + a2**3*r1*r3*r4 - a3**3*r1*r2*r4 - a4**3*r1*r2* &
      r3
f(405) = a1*r2**2*r3*r4 + a2*r1**2*r3*r4 - a3*r1*r2*r4**2 - a4*r1*r2*r3 &
      **2
f(406) = a1**2*r2**2*r3*r4 + a2**2*r1**2*r3*r4 - a3**2*r1*r2*r4**2 - a4 &
      **2*r1*r2*r3**2
f(407) = a1*r2**3*r3*r4 + a2*r1**3*r3*r4 - a3*r1*r2*r4**3 - a4*r1*r2*r3 &
      **3
f(408) = -a1*r2*r3**2*r4 - a2*r1*r3*r4**2 + a3*r1**2*r2*r4 + a4*r1*r2**2 &
      *r3
f(409) = -a1**2*r2*r3**2*r4 - a2**2*r1*r3*r4**2 + a3**2*r1**2*r2*r4 + a4 &
      **2*r1*r2**2*r3
f(410) = a1*r2**2*r3**2*r4 + a2*r1**2*r3*r4**2 - a3*r1**2*r2*r4**2 - a4* &
      r1*r2**2*r3**2
f(411) = a1*r2*r3**3*r4 + a2*r1*r3*r4**3 - a3*r1**3*r2*r4 - a4*r1*r2**3* &
      r3
f(412) = a1*r2*r3*r4**2 + a2*r1*r3**2*r4 - a3*r1*r2**2*r4 - a4*r1**2*r2* &
      r3
f(413) = a1**2*r2*r3*r4**2 + a2**2*r1*r3**2*r4 - a3**2*r1*r2**2*r4 - a4 &
      **2*r1**2*r2*r3
f(414) = -a1*r2**2*r3*r4**2 - a2*r1**2*r3**2*r4 + a3*r1*r2**2*r4**2 + a4 &
      *r1**2*r2*r3**2
f(415) = -a1*r2*r3**2*r4**2 - a2*r1*r3**2*r4**2 + a3*r1**2*r2**2*r4 + a4 &
      *r1**2*r2**2*r3
f(416) = a1*r2*r3*r4**3 + a2*r1*r3**3*r4 - a3*r1*r2**3*r4 - a4*r1**3*r2* &
      r3
f(417) = -b1**2*r1*r2*r3 - b1**2*r1*r2*r4 + b2**2*r1*r3*r4 + b2**2*r2*r3 &
      *r4
f(418) = b1**2*r1*r2*r3**2 + b1**2*r1*r2*r4**2 - b2**2*r1**2*r3*r4 - b2 &
      **2*r2**2*r3*r4
f(419) = b1**2*r1**2*r2*r4 + b1**2*r1*r2**2*r3 - b2**2*r1*r3*r4**2 - b2 &
      **2*r2*r3**2*r4
f(420) = -b1**2*r1**2*r2*r3 - b1**2*r1*r2**2*r4 + b2**2*r1*r3**2*r4 + b2 &
      **2*r2*r3*r4**2
f(421) = b1**2*r1*r3*r4 + b1**2*r2*r3*r4 - b2**2*r1*r2*r3 - b2**2*r1*r2* &
      r4
f(422) = b1**2*r1**2*r3*r4 + b1**2*r2**2*r3*r4 - b2**2*r1*r2*r3**2 - b2 &
      **2*r1*r2*r4**2
f(423) = -b1**2*r1*r3*r4**2 - b1**2*r2*r3**2*r4 + b2**2*r1**2*r2*r4 + b2 &
      **2*r1*r2**2*r3
f(424) = -b1**2*r1*r3**2*r4 - b1**2*r2*r3*r4**2 + b2**2*r1**2*r2*r3 + b2 &
      **2*r1*r2**2*r4
f(425) = dtau**2*(-r1*r2*r3 - r1*r2*r4 + r1*r3*r4 + r2*r3*r4)
f(426) = dtau**2*(r1**2*r3*r4 - r1*r2*r3**2 - r1*r2*r4**2 + r2**2*r3*r4)
f(427) = dtau**2*(-r1**2*r2*r4 - r1*r2**2*r3 + r1*r3*r4**2 + r2*r3**2*r4 &
      )
f(428) = dtau**2*(-r1**2*r2*r3 - r1*r2**2*r4 + r1*r3**2*r4 + r2*r3*r4**2 &
      )
f(429) = a1*a2*r1*r2 - a3*a4*r3*r4
f(430) = a1**2*a2*r1*r2 + a1*a2**2*r1*r2 - a3**2*a4*r3*r4 - a3*a4**2*r3* &
      r4
f(431) = -a1**3*a2*r1*r2 - a1*a2**3*r1*r2 + a3**3*a4*r3*r4 + a3*a4**3*r3 &
      *r4
f(432) = -a1**2*a2**2*r1*r2 + a3**2*a4**2*r3*r4
f(433) = -a1*a2*r1**2*r2 - a1*a2*r1*r2**2 + a3*a4*r3**2*r4 + a3*a4*r3*r4 &
      **2
f(434) = a1**2*a2*r1**2*r2 + a1*a2**2*r1*r2**2 - a3**2*a4*r3**2*r4 - a3* &
      a4**2*r3*r4**2
f(435) = a1**2*a2*r1*r2**2 + a1*a2**2*r1**2*r2 - a3**2*a4*r3*r4**2 - a3* &
      a4**2*r3**2*r4
f(436) = a1*a2*r1**3*r2 + a1*a2*r1*r2**3 - a3*a4*r3**3*r4 - a3*a4*r3*r4 &
      **3
f(437) = a1*a2*r1**2*r2**2 - a3*a4*r3**2*r4**2
f(438) = -a1*a3*r1*r2 + a1*a3*r3*r4 - a2*a4*r1*r2 + a2*a4*r3*r4
f(439) = a1**2*a3*r3*r4 - a1*a3**2*r1*r2 + a2**2*a4*r3*r4 - a2*a4**2*r1* &
      r2
f(440) = a1**3*a3*r3*r4 - a1*a3**3*r1*r2 + a2**3*a4*r3*r4 - a2*a4**3*r1* &
      r2
f(441) = -a1**2*a3*r1*r2 + a1*a3**2*r3*r4 - a2**2*a4*r1*r2 + a2*a4**2*r3 &
      *r4
f(442) = a1**2*a3**2*r1*r2 - a1**2*a3**2*r3*r4 + a2**2*a4**2*r1*r2 - a2 &
      **2*a4**2*r3*r4
f(443) = -a1**3*a3*r1*r2 + a1*a3**3*r3*r4 - a2**3*a4*r1*r2 + a2*a4**3*r3 &
      *r4
f(444) = -a1*a3*r1*r2**2 + a1*a3*r3*r4**2 - a2*a4*r1**2*r2 + a2*a4*r3**2 &
      *r4
f(445) = -a1**2*a3*r3*r4**2 + a1*a3**2*r1*r2**2 - a2**2*a4*r3**2*r4 + a2 &
      *a4**2*r1**2*r2
f(446) = a1**2*a3*r1*r2**2 - a1*a3**2*r3*r4**2 + a2**2*a4*r1**2*r2 - a2* &
      a4**2*r3**2*r4
f(447) = -a1*a3*r1*r2**3 + a1*a3*r3*r4**3 - a2*a4*r1**3*r2 + a2*a4*r3**3 &
      *r4
f(448) = -a1*a3*r1**2*r2 + a1*a3*r3**2*r4 - a2*a4*r1*r2**2 + a2*a4*r3*r4 &
      **2
f(449) = -a1**2*a3*r3**2*r4 + a1*a3**2*r1**2*r2 - a2**2*a4*r3*r4**2 + a2 &
      *a4**2*r1*r2**2
f(450) = a1**2*a3*r1**2*r2 - a1*a3**2*r3**2*r4 + a2**2*a4*r1*r2**2 - a2* &
      a4**2*r3*r4**2
f(451) = a1*a3*r1**2*r2**2 - a1*a3*r3**2*r4**2 + a2*a4*r1**2*r2**2 - a2* &
      a4*r3**2*r4**2
f(452) = -a1*a3*r1**3*r2 + a1*a3*r3**3*r4 - a2*a4*r1*r2**3 + a2*a4*r3*r4 &
      **3
f(453) = -a1*a4*r1*r2 + a1*a4*r3*r4 - a2*a3*r1*r2 + a2*a3*r3*r4
f(454) = a1**2*a4*r3*r4 - a1*a4**2*r1*r2 + a2**2*a3*r3*r4 - a2*a3**2*r1* &
      r2
f(455) = a1**3*a4*r3*r4 - a1*a4**3*r1*r2 + a2**3*a3*r3*r4 - a2*a3**3*r1* &
      r2
f(456) = -a1**2*a4*r1*r2 + a1*a4**2*r3*r4 - a2**2*a3*r1*r2 + a2*a3**2*r3 &
      *r4
f(457) = -a1**2*a4**2*r1*r2 + a1**2*a4**2*r3*r4 - a2**2*a3**2*r1*r2 + a2 &
      **2*a3**2*r3*r4
f(458) = -a1**3*a4*r1*r2 + a1*a4**3*r3*r4 - a2**3*a3*r1*r2 + a2*a3**3*r3 &
      *r4
f(459) = -a1*a4*r1*r2**2 + a1*a4*r3**2*r4 - a2*a3*r1**2*r2 + a2*a3*r3*r4 &
      **2
f(460) = -a1**2*a4*r3**2*r4 + a1*a4**2*r1*r2**2 - a2**2*a3*r3*r4**2 + a2 &
      *a3**2*r1**2*r2
f(461) = -a1**2*a4*r1*r2**2 + a1*a4**2*r3**2*r4 - a2**2*a3*r1**2*r2 + a2 &
      *a3**2*r3*r4**2
f(462) = a1*a4*r1*r2**3 - a1*a4*r3**3*r4 + a2*a3*r1**3*r2 - a2*a3*r3*r4 &
      **3
f(463) = -a1*a4*r1**2*r2 + a1*a4*r3*r4**2 - a2*a3*r1*r2**2 + a2*a3*r3**2 &
      *r4
f(464) = -a1**2*a4*r3*r4**2 + a1*a4**2*r1**2*r2 - a2**2*a3*r3**2*r4 + a2 &
      *a3**2*r1*r2**2
f(465) = -a1**2*a4*r1**2*r2 + a1*a4**2*r3*r4**2 - a2**2*a3*r1*r2**2 + a2 &
      *a3**2*r3**2*r4
f(466) = a1*a4*r1**2*r2**2 - a1*a4*r3**2*r4**2 + a2*a3*r1**2*r2**2 - a2* &
      a3*r3**2*r4**2
f(467) = -a1*a4*r1**3*r2 + a1*a4*r3*r4**3 - a2*a3*r1*r2**3 + a2*a3*r3**3 &
      *r4
f(468) = -a1*b1**2*r1*r2 - a2*b1**2*r1*r2 + a3*b2**2*r3*r4 + a4*b2**2*r3 &
      *r4
f(469) = -a1**2*b1**2*r1*r2 - a2**2*b1**2*r1*r2 + a3**2*b2**2*r3*r4 + a4 &
      **2*b2**2*r3*r4
f(470) = a1*b1**2*r1*r2**2 + a2*b1**2*r1**2*r2 - a3*b2**2*r3*r4**2 - a4* &
      b2**2*r3**2*r4
f(471) = a1*b1**2*r1**2*r2 + a2*b1**2*r1*r2**2 - a3*b2**2*r3**2*r4 - a4* &
      b2**2*r3*r4**2
f(472) = a1*b2**2*r1*r2 + a2*b2**2*r1*r2 - a3*b1**2*r3*r4 - a4*b1**2*r3* &
      r4
f(473) = -a1**2*b2**2*r1*r2 - a2**2*b2**2*r1*r2 + a3**2*b1**2*r3*r4 + a4 &
      **2*b1**2*r3*r4
f(474) = a1*b2**2*r1*r2**2 + a2*b2**2*r1**2*r2 - a3*b1**2*r3*r4**2 - a4* &
      b1**2*r3**2*r4
f(475) = a1*b2**2*r1**2*r2 + a2*b2**2*r1*r2**2 - a3*b1**2*r3**2*r4 - a4* &
      b1**2*r3*r4**2
f(476) = dtau**2*(a1*r1*r2 + a2*r1*r2 - a3*r3*r4 - a4*r3*r4)
f(477) = dtau**2*(-a1**2*r1*r2 - a2**2*r1*r2 + a3**2*r3*r4 + a4**2*r3*r4 &
      )
f(478) = dtau**2*(-a1*r1*r2**2 - a2*r1**2*r2 + a3*r3*r4**2 + a4*r3**2*r4 &
      )
f(479) = dtau**2*(-a1*r1**2*r2 - a2*r1*r2**2 + a3*r3**2*r4 + a4*r3*r4**2 &
      )
f(480) = a1*a2*r3*r4 - a3*a4*r1*r2
f(481) = a1**2*a2*r3*r4 + a1*a2**2*r3*r4 - a3**2*a4*r1*r2 - a3*a4**2*r1* &
      r2
f(482) = a1**3*a2*r3*r4 + a1*a2**3*r3*r4 - a3**3*a4*r1*r2 - a3*a4**3*r1* &
      r2
f(483) = -a1**2*a2**2*r3*r4 + a3**2*a4**2*r1*r2
f(484) = a1*a2*r3**2*r4 + a1*a2*r3*r4**2 - a3*a4*r1**2*r2 - a3*a4*r1*r2 &
      **2
f(485) = -a1**2*a2*r3**2*r4 - a1*a2**2*r3*r4**2 + a3**2*a4*r1**2*r2 + a3 &
      *a4**2*r1*r2**2
f(486) = -a1**2*a2*r3*r4**2 - a1*a2**2*r3**2*r4 + a3**2*a4*r1*r2**2 + a3 &
      *a4**2*r1**2*r2
f(487) = a1*a2*r3**3*r4 + a1*a2*r3*r4**3 - a3*a4*r1**3*r2 - a3*a4*r1*r2 &
      **3
f(488) = a1*a2*r3**2*r4**2 - a3*a4*r1**2*r2**2
f(489) = -a1*a2*r3**3*r4 - a1*a2*r3*r4**3 + a3*a4*r1**3*r2 + a3*a4*r1*r2 &
      **3
f(490) = -a1*b2**2*r3*r4 - a2*b2**2*r3*r4 + a3*b1**2*r1*r2 + a4*b1**2*r1 &
      *r2
f(491) = -a1**2*b2**2*r3*r4 - a2**2*b2**2*r3*r4 + a3**2*b1**2*r1*r2 + a4 &
      **2*b1**2*r1*r2
f(492) = a1*b2**2*r3*r4**2 + a2*b2**2*r3**2*r4 - a3*b1**2*r1*r2**2 - a4* &
      b1**2*r1**2*r2
f(493) = a1*b2**2*r3**2*r4 + a2*b2**2*r3*r4**2 - a3*b1**2*r1**2*r2 - a4* &
      b1**2*r1*r2**2
f(494) = a1*b1**2*r3*r4 + a2*b1**2*r3*r4 - a3*b2**2*r1*r2 - a4*b2**2*r1* &
      r2
f(495) = -a1**2*b1**2*r3*r4 - a2**2*b1**2*r3*r4 + a3**2*b2**2*r1*r2 + a4 &
      **2*b2**2*r1*r2
f(496) = a1*b1**2*r3*r4**2 + a2*b1**2*r3**2*r4 - a3*b2**2*r1*r2**2 - a4* &
      b2**2*r1**2*r2
f(497) = a1*b1**2*r3**2*r4 + a2*b1**2*r3*r4**2 - a3*b2**2*r1**2*r2 - a4* &
      b2**2*r1*r2**2
f(498) = dtau**2*(a1*r3*r4 + a2*r3*r4 - a3*r1*r2 - a4*r1*r2)
f(499) = dtau**2*(-a1**2*r3*r4 - a2**2*r3*r4 + a3**2*r1*r2 + a4**2*r1*r2 &
      )
f(500) = dtau**2*(-a1*r3*r4**2 - a2*r3**2*r4 + a3*r1*r2**2 + a4*r1**2*r2 &
      )
f(501) = dtau**2*(-a1*r3**2*r4 - a2*r3*r4**2 + a3*r1**2*r2 + a4*r1*r2**2 &
      )
f(502) = b1*b2*(-r1*r2 + r3*r4)
f(503) = b1*b2*(b1**2*r3*r4 - b2**2*r1*r2)
f(504) = b1**2*b2**2*(r1*r2 - r3*r4)
f(505) = b1*b2*(-b1**2*r1*r2 + b2**2*r3*r4)
f(506) = b1*b2*(r1**2*r2 + r1*r2**2 - r3**2*r4 - r3*r4**2)
f(507) = b1*b2*(r1**3*r2 + r1*r2**3 - r3**3*r4 - r3*r4**3)
f(508) = b1*b2*(r1**2*r2**2 - r3**2*r4**2)
f(509) = dtau**2*(b1**2*r1*r2 - b2**2*r3*r4)
f(510) = dtau*(-b1*r1**2*r2 + b1*r1*r2**2 + b2*r3**2*r4 - b2*r3*r4**2)
f(511) = dtau*(-b1*r1**3*r2 + b1*r1*r2**3 + b2*r3**3*r4 - b2*r3*r4**3)
f(512) = dtau*(b1*r1**2*r2 - b1*r1*r2**2 - b2*r3**2*r4 + b2*r3*r4**2)
f(513) = dtau*(b1*r1**3*r2 - b1*r1*r2**3 - b2*r3**3*r4 + b2*r3*r4**3)
f(514) = dtau**2*(b1**2*r3*r4 - b2**2*r1*r2)
f(515) = dtau*(b1*r3**2*r4 - b1*r3*r4**2 - b2*r1**2*r2 + b2*r1*r2**2)
f(516) = dtau*(-b1*r3**3*r4 + b1*r3*r4**3 + b2*r1**3*r2 - b2*r1*r2**3)
f(517) = dtau*(-b1*r3**2*r4 + b1*r3*r4**2 + b2*r1**2*r2 - b2*r1*r2**2)
f(518) = dtau*(b1*r3**3*r4 - b1*r3*r4**3 - b2*r1**3*r2 + b2*r1*r2**3)
f(519) = a1*a2*r1*r3 + a1*a2*r2*r4 - a3*a4*r1*r3 - a3*a4*r2*r4
f(520) = -a1**2*a2*r2*r4 - a1*a2**2*r1*r3 + a3**2*a4*r2*r4 + a3*a4**2*r1 &
      *r3
f(521) = a1**3*a2*r2*r4 + a1*a2**3*r1*r3 - a3**3*a4*r2*r4 - a3*a4**3*r1* &
      r3
f(522) = a1**2*a2*r1*r3 + a1*a2**2*r2*r4 - a3**2*a4*r1*r3 - a3*a4**2*r2* &
      r4
f(523) = -a1**2*a2**2*r1*r3 - a1**2*a2**2*r2*r4 + a3**2*a4**2*r1*r3 + a3 &
      **2*a4**2*r2*r4
f(524) = -a1**3*a2*r1*r3 - a1*a2**3*r2*r4 + a3**3*a4*r1*r3 + a3*a4**3*r2 &
      *r4
f(525) = -a1*a2*r1*r3**2 - a1*a2*r2*r4**2 + a3*a4*r1**2*r3 + a3*a4*r2**2 &
      *r4
f(526) = -a1**2*a2*r2*r4**2 - a1*a2**2*r1*r3**2 + a3**2*a4*r2**2*r4 + a3 &
      *a4**2*r1**2*r3
f(527) = a1**2*a2*r1*r3**2 + a1*a2**2*r2*r4**2 - a3**2*a4*r1**2*r3 - a3* &
      a4**2*r2**2*r4
f(528) = a1*a2*r1*r3**3 + a1*a2*r2*r4**3 - a3*a4*r1**3*r3 - a3*a4*r2**3* &
      r4
f(529) = a1*a2*r1**2*r3 + a1*a2*r2**2*r4 - a3*a4*r1*r3**2 - a3*a4*r2*r4 &
      **2
f(530) = -a1**2*a2*r2**2*r4 - a1*a2**2*r1**2*r3 + a3**2*a4*r2*r4**2 + a3 &
      *a4**2*r1*r3**2
f(531) = a1**2*a2*r1**2*r3 + a1*a2**2*r2**2*r4 - a3**2*a4*r1*r3**2 - a3* &
      a4**2*r2*r4**2
f(532) = -a1*a2*r1**2*r3**2 - a1*a2*r2**2*r4**2 + a3*a4*r1**2*r3**2 + a3 &
      *a4*r2**2*r4**2
f(533) = a1*a2*r1**3*r3 + a1*a2*r2**3*r4 - a3*a4*r1*r3**3 - a3*a4*r2*r4 &
      **3
f(534) = -a1**2*a3*r1*r3 + a1*a3**2*r1*r3 - a2**2*a4*r2*r4 + a2*a4**2*r2 &
      *r4
f(535) = -a1**3*a3*r1*r3 + a1*a3**3*r1*r3 - a2**3*a4*r2*r4 + a2*a4**3*r2 &
      *r4
f(536) = -a1*a3*r1**2*r3 + a1*a3*r1*r3**2 - a2*a4*r2**2*r4 + a2*a4*r2*r4 &
      **2
f(537) = a1**2*a3*r1**2*r3 - a1*a3**2*r1*r3**2 + a2**2*a4*r2**2*r4 - a2* &
      a4**2*r2*r4**2
f(538) = a1**2*a3*r1*r3**2 - a1*a3**2*r1**2*r3 + a2**2*a4*r2*r4**2 - a2* &
      a4**2*r2**2*r4
f(539) = a1*a3*r1**3*r3 - a1*a3*r1*r3**3 + a2*a4*r2**3*r4 - a2*a4*r2*r4 &
      **3
f(540) = a1*a3*r1**2*r3 - a1*a3*r1*r3**2 + a2*a4*r2**2*r4 - a2*a4*r2*r4 &
      **2
f(541) = -a1*a4*r1*r3 + a1*a4*r2*r4 + a2*a3*r1*r3 - a2*a3*r2*r4
f(542) = -a1**2*a4*r2*r4 + a1*a4**2*r1*r3 - a2**2*a3*r1*r3 + a2*a3**2*r2 &
      *r4
f(543) = -a1**3*a4*r2*r4 + a1*a4**3*r1*r3 - a2**3*a3*r1*r3 + a2*a3**3*r2 &
      *r4
f(544) = a1**2*a4*r1*r3 - a1*a4**2*r2*r4 + a2**2*a3*r2*r4 - a2*a3**2*r1* &
      r3
f(545) = -a1**2*a4**2*r1*r3 + a1**2*a4**2*r2*r4 + a2**2*a3**2*r1*r3 - a2 &
      **2*a3**2*r2*r4
f(546) = -a1**3*a4*r1*r3 + a1*a4**3*r2*r4 - a2**3*a3*r2*r4 + a2*a3**3*r1 &
      *r3
f(547) = a1*a4*r1*r3**2 - a1*a4*r2**2*r4 - a2*a3*r1**2*r3 + a2*a3*r2*r4 &
      **2
f(548) = -a1**2*a4*r2**2*r4 + a1*a4**2*r1*r3**2 - a2**2*a3*r1**2*r3 + a2 &
      *a3**2*r2*r4**2
f(549) = a1**2*a4*r1*r3**2 - a1*a4**2*r2**2*r4 + a2**2*a3*r2*r4**2 - a2* &
      a3**2*r1**2*r3
f(550) = -a1*a4*r1*r3**3 + a1*a4*r2**3*r4 + a2*a3*r1**3*r3 - a2*a3*r2*r4 &
      **3
f(551) = a1*a4*r1**2*r3 - a1*a4*r2*r4**2 - a2*a3*r1*r3**2 + a2*a3*r2**2* &
      r4
f(552) = -a1**2*a4*r2*r4**2 + a1*a4**2*r1**2*r3 - a2**2*a3*r1*r3**2 + a2 &
      *a3**2*r2**2*r4
f(553) = a1**2*a4*r1**2*r3 - a1*a4**2*r2*r4**2 + a2**2*a3*r2**2*r4 - a2* &
      a3**2*r1*r3**2
f(554) = -a1*a4*r1**2*r3**2 + a1*a4*r2**2*r4**2 + a2*a3*r1**2*r3**2 - a2 &
      *a3*r2**2*r4**2
f(555) = -a1*a4*r1**3*r3 + a1*a4*r2*r4**3 + a2*a3*r1*r3**3 - a2*a3*r2**3 &
      *r4
f(556) = a1*b1**2*r1*r3 + a2*b1**2*r2*r4 - a3*b2**2*r1*r3 - a4*b2**2*r2* &
      r4
f(557) = a1**2*b1**2*r1*r3 + a2**2*b1**2*r2*r4 - a3**2*b2**2*r1*r3 - a4 &
      **2*b2**2*r2*r4
f(558) = a1*b1**2*r1*r3**2 + a2*b1**2*r2*r4**2 - a3*b2**2*r1**2*r3 - a4* &
      b2**2*r2**2*r4
f(559) = a1*b1**2*r1**2*r3 + a2*b1**2*r2**2*r4 - a3*b2**2*r1*r3**2 - a4* &
      b2**2*r2*r4**2
f(560) = a1*b2**2*r1*r3 + a2*b2**2*r2*r4 - a3*b1**2*r1*r3 - a4*b1**2*r2* &
      r4
f(561) = -a1**2*b2**2*r1*r3 - a2**2*b2**2*r2*r4 + a3**2*b1**2*r1*r3 + a4 &
      **2*b1**2*r2*r4
f(562) = a1*b2**2*r1*r3**2 + a2*b2**2*r2*r4**2 - a3*b1**2*r1**2*r3 - a4* &
      b1**2*r2**2*r4
f(563) = a1*b2**2*r1**2*r3 + a2*b2**2*r2**2*r4 - a3*b1**2*r1*r3**2 - a4* &
      b1**2*r2*r4**2
f(564) = dtau**2*(a1*r1*r3 + a2*r2*r4 - a3*r1*r3 - a4*r2*r4)
f(565) = dtau**2*(a1**2*r1*r3 + a2**2*r2*r4 - a3**2*r1*r3 - a4**2*r2*r4)
f(566) = dtau**2*(-a1*r1*r3**2 - a2*r2*r4**2 + a3*r1**2*r3 + a4*r2**2*r4 &
      )
f(567) = dtau**2*(-a1*r1**2*r3 - a2*r2**2*r4 + a3*r1*r3**2 + a4*r2*r4**2 &
      )
f(568) = -a1**2*a3*r2*r4 + a1*a3**2*r2*r4 - a2**2*a4*r1*r3 + a2*a4**2*r1 &
      *r3
f(569) = a1**3*a3*r2*r4 - a1*a3**3*r2*r4 + a2**3*a4*r1*r3 - a2*a4**3*r1* &
      r3
f(570) = a1*a3*r2**2*r4 - a1*a3*r2*r4**2 + a2*a4*r1**2*r3 - a2*a4*r1*r3 &
      **2
f(571) = -a1**2*a3*r2**2*r4 + a1*a3**2*r2*r4**2 - a2**2*a4*r1**2*r3 + a2 &
      *a4**2*r1*r3**2
f(572) = -a1**2*a3*r2*r4**2 + a1*a3**2*r2**2*r4 - a2**2*a4*r1*r3**2 + a2 &
      *a4**2*r1**2*r3
f(573) = -a1*a3*r2**3*r4 + a1*a3*r2*r4**3 - a2*a4*r1**3*r3 + a2*a4*r1*r3 &
      **3
f(574) = -a1*a3*r2**2*r4 + a1*a3*r2*r4**2 - a2*a4*r1**2*r3 + a2*a4*r1*r3 &
      **2
f(575) = -a1*b1**2*r2*r4 - a2*b1**2*r1*r3 + a3*b2**2*r2*r4 + a4*b2**2*r1 &
      *r3
f(576) = -a1**2*b1**2*r2*r4 - a2**2*b1**2*r1*r3 + a3**2*b2**2*r2*r4 + a4 &
      **2*b2**2*r1*r3
f(577) = a1*b1**2*r2*r4**2 + a2*b1**2*r1*r3**2 - a3*b2**2*r2**2*r4 - a4* &
      b2**2*r1**2*r3
f(578) = a1*b1**2*r2**2*r4 + a2*b1**2*r1**2*r3 - a3*b2**2*r2*r4**2 - a4* &
      b2**2*r1*r3**2
f(579) = -a1*b2**2*r2*r4 - a2*b2**2*r1*r3 + a3*b1**2*r2*r4 + a4*b1**2*r1 &
      *r3
f(580) = a1**2*b2**2*r2*r4 + a2**2*b2**2*r1*r3 - a3**2*b1**2*r2*r4 - a4 &
      **2*b1**2*r1*r3
f(581) = -a1*b2**2*r2*r4**2 - a2*b2**2*r1*r3**2 + a3*b1**2*r2**2*r4 + a4 &
      *b1**2*r1**2*r3
f(582) = a1*b2**2*r2**2*r4 + a2*b2**2*r1**2*r3 - a3*b1**2*r2*r4**2 - a4* &
      b1**2*r1*r3**2
f(583) = dtau**2*(a1*r2*r4 + a2*r1*r3 - a3*r2*r4 - a4*r1*r3)
f(584) = dtau**2*(a1**2*r2*r4 + a2**2*r1*r3 - a3**2*r2*r4 - a4**2*r1*r3)
f(585) = dtau**2*(-a1*r2*r4**2 - a2*r1*r3**2 + a3*r2**2*r4 + a4*r1**2*r3 &
      )
f(586) = dtau**2*(-a1*r2**2*r4 - a2*r1**2*r3 + a3*r2*r4**2 + a4*r1*r3**2 &
      )
f(587) = b1*b2*(-b1**2*r1*r3 - b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2*r4)
f(588) = b1*b2*(r1**2*r3 - r1*r3**2 + r2**2*r4 - r2*r4**2)
f(589) = b1*b2*(r1**3*r3 - r1*r3**3 + r2**3*r4 - r2*r4**3)
f(590) = dtau*(-b1*r1*r3 + b1*r2*r4 + b2*r1*r3 - b2*r2*r4)
f(591) = dtau**3*(b1*r1*r3 - b1*r2*r4 - b2*r1*r3 + b2*r2*r4)
f(592) = dtau**2*(-b1**2*r1*r3 - b1**2*r2*r4 + b2**2*r1*r3 + b2**2*r2*r4 &
      )
f(593) = dtau*(b1**3*r1*r3 - b1**3*r2*r4 - b2**3*r1*r3 + b2**3*r2*r4)
f(594) = dtau*(b1*r1*r3**2 - b1*r2*r4**2 - b2*r1**2*r3 + b2*r2**2*r4)
f(595) = dtau*(-b1*r1*r3**3 + b1*r2*r4**3 + b2*r1**3*r3 - b2*r2**3*r4)
f(596) = dtau*(b1*r1**2*r3 - b1*r2**2*r4 - b2*r1*r3**2 + b2*r2*r4**2)
f(597) = dtau*(-b1*r1**2*r3**2 + b1*r2**2*r4**2 + b2*r1**2*r3**2 - b2*r2 &
      **2*r4**2)
f(598) = dtau*(b1*r1**3*r3 - b1*r2**3*r4 - b2*r1*r3**3 + b2*r2*r4**3)
f(599) = -a1*a2*r1*r4 - a1*a2*r2*r3 + a3*a4*r1*r4 + a3*a4*r2*r3
f(600) = -a1**2*a2*r2*r3 - a1*a2**2*r1*r4 + a3**2*a4*r1*r4 + a3*a4**2*r2 &
      *r3
f(601) = a1**3*a2*r2*r3 + a1*a2**3*r1*r4 - a3**3*a4*r1*r4 - a3*a4**3*r2* &
      r3
f(602) = a1**2*a2*r1*r4 + a1*a2**2*r2*r3 - a3**2*a4*r2*r3 - a3*a4**2*r1* &
      r4
f(603) = a1**2*a2**2*r1*r4 + a1**2*a2**2*r2*r3 - a3**2*a4**2*r1*r4 - a3 &
      **2*a4**2*r2*r3
f(604) = -a1**3*a2*r1*r4 - a1*a2**3*r2*r3 + a3**3*a4*r2*r3 + a3*a4**3*r1 &
      *r4
f(605) = a1*a2*r1*r4**2 + a1*a2*r2*r3**2 - a3*a4*r1**2*r4 - a3*a4*r2**2* &
      r3
f(606) = a1**2*a2*r2*r3**2 + a1*a2**2*r1*r4**2 - a3**2*a4*r1**2*r4 - a3* &
      a4**2*r2**2*r3
f(607) = a1**2*a2*r1*r4**2 + a1*a2**2*r2*r3**2 - a3**2*a4*r2**2*r3 - a3* &
      a4**2*r1**2*r4
f(608) = a1*a2*r1*r4**3 + a1*a2*r2*r3**3 - a3*a4*r1**3*r4 - a3*a4*r2**3* &
      r3
f(609) = a1*a2*r1**2*r4 + a1*a2*r2**2*r3 - a3*a4*r1*r4**2 - a3*a4*r2*r3 &
      **2
f(610) = a1**2*a2*r2**2*r3 + a1*a2**2*r1**2*r4 - a3**2*a4*r1*r4**2 - a3* &
      a4**2*r2*r3**2
f(611) = a1**2*a2*r1**2*r4 + a1*a2**2*r2**2*r3 - a3**2*a4*r2*r3**2 - a3* &
      a4**2*r1*r4**2
f(612) = a1*a2*r1**2*r4**2 + a1*a2*r2**2*r3**2 - a3*a4*r1**2*r4**2 - a3* &
      a4*r2**2*r3**2
f(613) = -a1*a2*r1**3*r4 - a1*a2*r2**3*r3 + a3*a4*r1*r4**3 + a3*a4*r2*r3 &
      **3
f(614) = -a1*a3*r1*r4 + a1*a3*r2*r3 + a2*a4*r1*r4 - a2*a4*r2*r3
f(615) = a1**2*a3*r2*r3 - a1*a3**2*r1*r4 + a2**2*a4*r1*r4 - a2*a4**2*r2* &
      r3
f(616) = -a1**3*a3*r2*r3 + a1*a3**3*r1*r4 - a2**3*a4*r1*r4 + a2*a4**3*r2 &
      *r3
f(617) = a1**2*a3*r1*r4 - a1*a3**2*r2*r3 + a2**2*a4*r2*r3 - a2*a4**2*r1* &
      r4
f(618) = -a1**2*a3**2*r1*r4 + a1**2*a3**2*r2*r3 + a2**2*a4**2*r1*r4 - a2 &
      **2*a4**2*r2*r3
f(619) = -a1**3*a3*r1*r4 + a1*a3**3*r2*r3 - a2**3*a4*r2*r3 + a2*a4**3*r1 &
      *r4
f(620) = a1*a3*r1*r4**2 - a1*a3*r2**2*r3 - a2*a4*r1**2*r4 + a2*a4*r2*r3 &
      **2
f(621) = a1**2*a3*r2**2*r3 - a1*a3**2*r1*r4**2 + a2**2*a4*r1**2*r4 - a2* &
      a4**2*r2*r3**2
f(622) = -a1**2*a3*r1*r4**2 + a1*a3**2*r2**2*r3 - a2**2*a4*r2*r3**2 + a2 &
      *a4**2*r1**2*r4
f(623) = a1*a3*r1*r4**3 - a1*a3*r2**3*r3 - a2*a4*r1**3*r4 + a2*a4*r2*r3 &
      **3
f(624) = a1*a3*r1**2*r4 - a1*a3*r2*r3**2 - a2*a4*r1*r4**2 + a2*a4*r2**2* &
      r3
f(625) = a1**2*a3*r2*r3**2 - a1*a3**2*r1**2*r4 + a2**2*a4*r1*r4**2 - a2* &
      a4**2*r2**2*r3
f(626) = -a1**2*a3*r1**2*r4 + a1*a3**2*r2*r3**2 - a2**2*a4*r2**2*r3 + a2 &
      *a4**2*r1*r4**2
f(627) = -a1*a3*r1**2*r4**2 + a1*a3*r2**2*r3**2 + a2*a4*r1**2*r4**2 - a2 &
      *a4*r2**2*r3**2
f(628) = a1*a3*r1**3*r4 - a1*a3*r2*r3**3 - a2*a4*r1*r4**3 + a2*a4*r2**3* &
      r3
f(629) = -a1**2*a4*r1*r4 + a1*a4**2*r1*r4 - a2**2*a3*r2*r3 + a2*a3**2*r2 &
      *r3
f(630) = -a1**3*a4*r1*r4 + a1*a4**3*r1*r4 - a2**3*a3*r2*r3 + a2*a3**3*r2 &
      *r3
f(631) = a1*a4*r1**2*r4 - a1*a4*r1*r4**2 + a2*a3*r2**2*r3 - a2*a3*r2*r3 &
      **2
f(632) = -a1**2*a4*r1**2*r4 + a1*a4**2*r1*r4**2 - a2**2*a3*r2**2*r3 + a2 &
      *a3**2*r2*r3**2
f(633) = -a1**2*a4*r1*r4**2 + a1*a4**2*r1**2*r4 - a2**2*a3*r2*r3**2 + a2 &
      *a3**2*r2**2*r3
f(634) = a1*a4*r1**3*r4 - a1*a4*r1*r4**3 + a2*a3*r2**3*r3 - a2*a3*r2*r3 &
      **3
f(635) = -a1*a4*r1**2*r4 + a1*a4*r1*r4**2 - a2*a3*r2**2*r3 + a2*a3*r2*r3 &
      **2
f(636) = a1*b1**2*r1*r4 + a2*b1**2*r2*r3 - a3*b2**2*r2*r3 - a4*b2**2*r1* &
      r4
f(637) = a1**2*b1**2*r1*r4 + a2**2*b1**2*r2*r3 - a3**2*b2**2*r2*r3 - a4 &
      **2*b2**2*r1*r4
f(638) = a1*b1**2*r1*r4**2 + a2*b1**2*r2*r3**2 - a3*b2**2*r2**2*r3 - a4* &
      b2**2*r1**2*r4
f(639) = a1*b1**2*r1**2*r4 + a2*b1**2*r2**2*r3 - a3*b2**2*r2*r3**2 - a4* &
      b2**2*r1*r4**2
f(640) = a1*b2**2*r1*r4 + a2*b2**2*r2*r3 - a3*b1**2*r2*r3 - a4*b1**2*r1* &
      r4
f(641) = a1**2*b2**2*r1*r4 + a2**2*b2**2*r2*r3 - a3**2*b1**2*r2*r3 - a4 &
      **2*b1**2*r1*r4
f(642) = -a1*b2**2*r1*r4**2 - a2*b2**2*r2*r3**2 + a3*b1**2*r2**2*r3 + a4 &
      *b1**2*r1**2*r4
f(643) = -a1*b2**2*r1**2*r4 - a2*b2**2*r2**2*r3 + a3*b1**2*r2*r3**2 + a4 &
      *b1**2*r1*r4**2
f(644) = dtau**2*(-a1*r1*r4 - a2*r2*r3 + a3*r2*r3 + a4*r1*r4)
f(645) = dtau**2*(a1**2*r1*r4 + a2**2*r2*r3 - a3**2*r2*r3 - a4**2*r1*r4)
f(646) = dtau**2*(a1*r1*r4**2 + a2*r2*r3**2 - a3*r2**2*r3 - a4*r1**2*r4)
f(647) = dtau**2*(a1*r1**2*r4 + a2*r2**2*r3 - a3*r2*r3**2 - a4*r1*r4**2)
f(648) = -a1**2*a4*r2*r3 + a1*a4**2*r2*r3 - a2**2*a3*r1*r4 + a2*a3**2*r1 &
      *r4
f(649) = a1**3*a4*r2*r3 - a1*a4**3*r2*r3 + a2**3*a3*r1*r4 - a2*a3**3*r1* &
      r4
f(650) = -a1*a4*r2**2*r3 + a1*a4*r2*r3**2 - a2*a3*r1**2*r4 + a2*a3*r1*r4 &
      **2
f(651) = a1**2*a4*r2**2*r3 - a1*a4**2*r2*r3**2 + a2**2*a3*r1**2*r4 - a2* &
      a3**2*r1*r4**2
f(652) = a1**2*a4*r2*r3**2 - a1*a4**2*r2**2*r3 + a2**2*a3*r1*r4**2 - a2* &
      a3**2*r1**2*r4
f(653) = -a1*a4*r2**3*r3 + a1*a4*r2*r3**3 - a2*a3*r1**3*r4 + a2*a3*r1*r4 &
      **3
f(654) = a1*a4*r2**2*r3 - a1*a4*r2*r3**2 + a2*a3*r1**2*r4 - a2*a3*r1*r4 &
      **2
f(655) = a1*b1**2*r2*r3 + a2*b1**2*r1*r4 - a3*b2**2*r1*r4 - a4*b2**2*r2* &
      r3
f(656) = a1**2*b1**2*r2*r3 + a2**2*b1**2*r1*r4 - a3**2*b2**2*r1*r4 - a4 &
      **2*b2**2*r2*r3
f(657) = -a1*b1**2*r2*r3**2 - a2*b1**2*r1*r4**2 + a3*b2**2*r1**2*r4 + a4 &
      *b2**2*r2**2*r3
f(658) = -a1*b1**2*r2**2*r3 - a2*b1**2*r1**2*r4 + a3*b2**2*r1*r4**2 + a4 &
      *b2**2*r2*r3**2
f(659) = a1*b2**2*r2*r3 + a2*b2**2*r1*r4 - a3*b1**2*r1*r4 - a4*b1**2*r2* &
      r3
f(660) = a1**2*b2**2*r2*r3 + a2**2*b2**2*r1*r4 - a3**2*b1**2*r1*r4 - a4 &
      **2*b1**2*r2*r3
f(661) = a1*b2**2*r2*r3**2 + a2*b2**2*r1*r4**2 - a3*b1**2*r1**2*r4 - a4* &
      b1**2*r2**2*r3
f(662) = a1*b2**2*r2**2*r3 + a2*b2**2*r1**2*r4 - a3*b1**2*r1*r4**2 - a4* &
      b1**2*r2*r3**2
f(663) = dtau**2*(-a1*r2*r3 - a2*r1*r4 + a3*r1*r4 + a4*r2*r3)
f(664) = dtau**2*(-a1**2*r2*r3 - a2**2*r1*r4 + a3**2*r1*r4 + a4**2*r2*r3 &
      )
f(665) = dtau**2*(a1*r2*r3**2 + a2*r1*r4**2 - a3*r1**2*r4 - a4*r2**2*r3)
f(666) = dtau**2*(a1*r2**2*r3 + a2*r1**2*r4 - a3*r1*r4**2 - a4*r2*r3**2)
f(667) = b1*b2*(b1**2*r1*r4 + b1**2*r2*r3 - b2**2*r1*r4 - b2**2*r2*r3)
f(668) = b1*b2*(-r1**2*r4 + r1*r4**2 - r2**2*r3 + r2*r3**2)
f(669) = b1*b2*(r1**3*r4 - r1*r4**3 + r2**3*r3 - r2*r3**3)
f(670) = b1*b2*(r1**2*r4 - r1*r4**2 + r2**2*r3 - r2*r3**2)
f(671) = dtau*(b1*r1*r4 - b1*r2*r3 + b2*r1*r4 - b2*r2*r3)
f(672) = dtau**3*(-b1*r1*r4 + b1*r2*r3 - b2*r1*r4 + b2*r2*r3)
f(673) = dtau**2*(b1**2*r1*r4 + b1**2*r2*r3 - b2**2*r1*r4 - b2**2*r2*r3)
f(674) = dtau*(b1**3*r1*r4 - b1**3*r2*r3 + b2**3*r1*r4 - b2**3*r2*r3)
f(675) = dtau*(b1*r1*r4**2 - b1*r2*r3**2 + b2*r1**2*r4 - b2*r2**2*r3)
f(676) = dtau*(-b1*r1*r4**3 + b1*r2*r3**3 - b2*r1**3*r4 + b2*r2**3*r3)
f(677) = dtau*(b1*r1**2*r4 - b1*r2**2*r3 + b2*r1*r4**2 - b2*r2*r3**2)
f(678) = dtau*(-b1*r1**2*r4**2 + b1*r2**2*r3**2 - b2*r1**2*r4**2 + b2*r2 &
      **2*r3**2)
f(679) = dtau*(b1*r1**3*r4 - b1*r2**3*r3 + b2*r1*r4**3 - b2*r2*r3**3)
f(680) = -a1*a2*a3*r1 - a1*a2*a4*r2 + a1*a3*a4*r3 + a2*a3*a4*r4
f(681) = -a1**2*a3*a4*r3 + a1*a2*a3**2*r1 + a1*a2*a4**2*r2 - a2**2*a3*a4 &
      *r4
f(682) = a1**3*a3*a4*r3 - a1*a2*a3**3*r1 - a1*a2*a4**3*r2 + a2**3*a3*a4* &
      r4
f(683) = a1**2*a2*a4*r2 + a1*a2**2*a3*r1 - a1*a3*a4**2*r3 - a2*a3**2*a4* &
      r4
f(684) = a1**2*a2*a4**2*r2 - a1**2*a3*a4**2*r3 + a1*a2**2*a3**2*r1 - a2 &
      **2*a3**2*a4*r4
f(685) = a1**3*a2*a4*r2 + a1*a2**3*a3*r1 - a1*a3*a4**3*r3 - a2*a3**3*a4* &
      r4
f(686) = a1**2*a2*a3*r1 + a1*a2**2*a4*r2 - a1*a3**2*a4*r3 - a2*a3*a4**2* &
      r4
f(687) = a1**2*a2*a3**2*r1 - a1**2*a3**2*a4*r3 + a1*a2**2*a4**2*r2 - a2 &
      **2*a3*a4**2*r4
f(688) = -a1**2*a2**2*a3*r1 - a1**2*a2**2*a4*r2 + a1*a3**2*a4**2*r3 + a2 &
      *a3**2*a4**2*r4
f(689) = a1**3*a2*a3*r1 + a1*a2**3*a4*r2 - a1*a3**3*a4*r3 - a2*a3*a4**3* &
      r4
f(690) = a1*a2*a3*r1**2 + a1*a2*a4*r2**2 - a1*a3*a4*r3**2 - a2*a3*a4*r4 &
      **2
f(691) = -a1**2*a3*a4*r3**2 + a1*a2*a3**2*r1**2 + a1*a2*a4**2*r2**2 - a2 &
      **2*a3*a4*r4**2
f(692) = a1**2*a2*a4*r2**2 + a1*a2**2*a3*r1**2 - a1*a3*a4**2*r3**2 - a2* &
      a3**2*a4*r4**2
f(693) = -a1**2*a2*a3*r1**2 - a1*a2**2*a4*r2**2 + a1*a3**2*a4*r3**2 + a2 &
      *a3*a4**2*r4**2
f(694) = -a1*a2*a3*r1**3 - a1*a2*a4*r2**3 + a1*a3*a4*r3**3 + a2*a3*a4*r4 &
      **3
f(695) = -a1*a2*a3*r2 - a1*a2*a4*r1 + a1*a3*a4*r4 + a2*a3*a4*r3
f(696) = -a1**2*a3*a4*r4 + a1*a2*a3**2*r2 + a1*a2*a4**2*r1 - a2**2*a3*a4 &
      *r3
f(697) = -a1**3*a3*a4*r4 + a1*a2*a3**3*r2 + a1*a2*a4**3*r1 - a2**3*a3*a4 &
      *r3
f(698) = a1**2*a2*a3*r2 + a1*a2**2*a4*r1 - a1*a3**2*a4*r4 - a2*a3*a4**2* &
      r3
f(699) = -a1**2*a2*a3**2*r2 + a1**2*a3**2*a4*r4 - a1*a2**2*a4**2*r1 + a2 &
      **2*a3*a4**2*r3
f(700) = -a1**3*a2*a3*r2 - a1*a2**3*a4*r1 + a1*a3**3*a4*r4 + a2*a3*a4**3 &
      *r3
f(701) = -a1**2*a2*a4*r1 - a1*a2**2*a3*r2 + a1*a3*a4**2*r4 + a2*a3**2*a4 &
      *r3
f(702) = a1**2*a2*a4**2*r1 - a1**2*a3*a4**2*r4 + a1*a2**2*a3**2*r2 - a2 &
      **2*a3**2*a4*r3
f(703) = -a1**2*a2**2*a3*r2 - a1**2*a2**2*a4*r1 + a1*a3**2*a4**2*r4 + a2 &
      *a3**2*a4**2*r3
f(704) = a1**3*a2*a4*r1 + a1*a2**3*a3*r2 - a1*a3*a4**3*r4 - a2*a3**3*a4* &
      r3
f(705) = a1*a2*a3*r2**2 + a1*a2*a4*r1**2 - a1*a3*a4*r4**2 - a2*a3*a4*r3 &
      **2
f(706) = -a1**2*a3*a4*r4**2 + a1*a2*a3**2*r2**2 + a1*a2*a4**2*r1**2 - a2 &
      **2*a3*a4*r3**2
f(707) = a1**2*a2*a3*r2**2 + a1*a2**2*a4*r1**2 - a1*a3**2*a4*r4**2 - a2* &
      a3*a4**2*r3**2
f(708) = -a1**2*a2*a4*r1**2 - a1*a2**2*a3*r2**2 + a1*a3*a4**2*r4**2 + a2 &
      *a3**2*a4*r3**2
f(709) = a1*a2*a3*r2**3 + a1*a2*a4*r1**3 - a1*a3*a4*r4**3 - a2*a3*a4*r3 &
      **3
f(710) = a1*a2*b1**2*r1 + a1*a2*b1**2*r2 - a3*a4*b2**2*r3 - a3*a4*b2**2* &
      r4
f(711) = a1**2*a2*b1**2*r2 + a1*a2**2*b1**2*r1 - a3**2*a4*b2**2*r4 - a3* &
      a4**2*b2**2*r3
f(712) = a1**2*a2*b1**2*r1 + a1*a2**2*b1**2*r2 - a3**2*a4*b2**2*r3 - a3* &
      a4**2*b2**2*r4
f(713) = a1*a2*b1**2*r1**2 + a1*a2*b1**2*r2**2 - a3*a4*b2**2*r3**2 - a3* &
      a4*b2**2*r4**2
f(714) = -a1*a2*b2**2*r1 - a1*a2*b2**2*r2 + a3*a4*b1**2*r3 + a3*a4*b1**2 &
      *r4
f(715) = a1**2*a2*b2**2*r2 + a1*a2**2*b2**2*r1 - a3**2*a4*b1**2*r4 - a3* &
      a4**2*b1**2*r3
f(716) = a1**2*a2*b2**2*r1 + a1*a2**2*b2**2*r2 - a3**2*a4*b1**2*r3 - a3* &
      a4**2*b1**2*r4
f(717) = a1*a2*b2**2*r1**2 + a1*a2*b2**2*r2**2 - a3*a4*b1**2*r3**2 - a3* &
      a4*b1**2*r4**2
f(718) = dtau**2*(a1*a2*r1 + a1*a2*r2 - a3*a4*r3 - a3*a4*r4)
f(719) = dtau**2*(-a1**2*a2*r2 - a1*a2**2*r1 + a3**2*a4*r4 + a3*a4**2*r3 &
      )
f(720) = dtau**2*(a1**2*a2*r1 + a1*a2**2*r2 - a3**2*a4*r3 - a3*a4**2*r4)
f(721) = dtau**2*(-a1*a2*r1**2 - a1*a2*r2**2 + a3*a4*r3**2 + a3*a4*r4**2 &
      )
f(722) = a1*a2*a3*r3 + a1*a2*a4*r4 - a1*a3*a4*r1 - a2*a3*a4*r2
f(723) = a1**2*a2*a4*r4 + a1*a2**2*a3*r3 - a1*a3*a4**2*r1 - a2*a3**2*a4* &
      r2
f(724) = a1**3*a2*a4*r4 + a1*a2**3*a3*r3 - a1*a3*a4**3*r1 - a2*a3**3*a4* &
      r2
f(725) = a1**2*a2*a3*r3 + a1*a2**2*a4*r4 - a1*a3**2*a4*r1 - a2*a3*a4**2* &
      r2
f(726) = -a1**2*a2**2*a3*r3 - a1**2*a2**2*a4*r4 + a1*a3**2*a4**2*r1 + a2 &
      *a3**2*a4**2*r2
f(727) = -a1**3*a2*a3*r3 - a1*a2**3*a4*r4 + a1*a3**3*a4*r1 + a2*a3*a4**3 &
      *r2
f(728) = a1**2*a3*a4*r1 - a1*a2*a3**2*r3 - a1*a2*a4**2*r4 + a2**2*a3*a4* &
      r2
f(729) = a1**2*a2*a4**2*r4 - a1**2*a3*a4**2*r1 + a1*a2**2*a3**2*r3 - a2 &
      **2*a3**2*a4*r2
f(730) = a1**2*a2*a3**2*r3 - a1**2*a3**2*a4*r1 + a1*a2**2*a4**2*r4 - a2 &
      **2*a3*a4**2*r2
f(731) = -a1**3*a3*a4*r1 + a1*a2*a3**3*r3 + a1*a2*a4**3*r4 - a2**3*a3*a4 &
      *r2
f(732) = -a1*a2*a3*r3**2 - a1*a2*a4*r4**2 + a1*a3*a4*r1**2 + a2*a3*a4*r2 &
      **2
f(733) = a1**2*a2*a4*r4**2 + a1*a2**2*a3*r3**2 - a1*a3*a4**2*r1**2 - a2* &
      a3**2*a4*r2**2
f(734) = a1**2*a2*a3*r3**2 + a1*a2**2*a4*r4**2 - a1*a3**2*a4*r1**2 - a2* &
      a3*a4**2*r2**2
f(735) = -a1**2*a3*a4*r1**2 + a1*a2*a3**2*r3**2 + a1*a2*a4**2*r4**2 - a2 &
      **2*a3*a4*r2**2
f(736) = a1*a2*a3*r3**3 + a1*a2*a4*r4**3 - a1*a3*a4*r1**3 - a2*a3*a4*r2 &
      **3
f(737) = a1*a3*b1**2*r1 - a1*a3*b2**2*r3 + a2*a4*b1**2*r2 - a2*a4*b2**2* &
      r4
f(738) = -a1**2*a3*b2**2*r3 + a1*a3**2*b1**2*r1 - a2**2*a4*b2**2*r4 + a2 &
      *a4**2*b1**2*r2
f(739) = a1**2*a3*b1**2*r1 - a1*a3**2*b2**2*r3 + a2**2*a4*b1**2*r2 - a2* &
      a4**2*b2**2*r4
f(740) = -a1*a3*b1**2*r1**2 + a1*a3*b2**2*r3**2 - a2*a4*b1**2*r2**2 + a2 &
      *a4*b2**2*r4**2
f(741) = -a1*a3*b1**2*r3 + a1*a3*b2**2*r1 - a2*a4*b1**2*r4 + a2*a4*b2**2 &
      *r2
f(742) = -a1**2*a3*b1**2*r3 + a1*a3**2*b2**2*r1 - a2**2*a4*b1**2*r4 + a2 &
      *a4**2*b2**2*r2
f(743) = a1**2*a3*b2**2*r1 - a1*a3**2*b1**2*r3 + a2**2*a4*b2**2*r2 - a2* &
      a4**2*b1**2*r4
f(744) = a1*a3*b1**2*r3**2 - a1*a3*b2**2*r1**2 + a2*a4*b1**2*r4**2 - a2* &
      a4*b2**2*r2**2
f(745) = dtau**2*(a1*a3*r1 - a1*a3*r3 + a2*a4*r2 - a2*a4*r4)
f(746) = dtau**2*(-a1**2*a3*r3 + a1*a3**2*r1 - a2**2*a4*r4 + a2*a4**2*r2 &
      )
f(747) = dtau**2*(-a1**2*a3*r1 + a1*a3**2*r3 - a2**2*a4*r2 + a2*a4**2*r4 &
      )
f(748) = dtau**2*(a1*a3*r1**2 - a1*a3*r3**2 + a2*a4*r2**2 - a2*a4*r4**2)
f(749) = a1*a4*b1**2*r1 - a1*a4*b2**2*r4 + a2*a3*b1**2*r2 - a2*a3*b2**2* &
      r3
f(750) = a1**2*a4*b2**2*r4 - a1*a4**2*b1**2*r1 + a2**2*a3*b2**2*r3 - a2* &
      a3**2*b1**2*r2
f(751) = -a1**2*a4*b1**2*r1 + a1*a4**2*b2**2*r4 - a2**2*a3*b1**2*r2 + a2 &
      *a3**2*b2**2*r3
f(752) = -a1*a4*b1**2*r1**2 + a1*a4*b2**2*r4**2 - a2*a3*b1**2*r2**2 + a2 &
      *a3*b2**2*r3**2
f(753) = -a1*a4*b1**2*r4 + a1*a4*b2**2*r1 - a2*a3*b1**2*r3 + a2*a3*b2**2 &
      *r2
f(754) = -a1**2*a4*b1**2*r4 + a1*a4**2*b2**2*r1 - a2**2*a3*b1**2*r3 + a2 &
      *a3**2*b2**2*r2
f(755) = -a1**2*a4*b2**2*r1 + a1*a4**2*b1**2*r4 - a2**2*a3*b2**2*r2 + a2 &
      *a3**2*b1**2*r3
f(756) = -a1*a4*b1**2*r4**2 + a1*a4*b2**2*r1**2 - a2*a3*b1**2*r3**2 + a2 &
      *a3*b2**2*r2**2
f(757) = dtau**2*(a1*a4*r1 - a1*a4*r4 + a2*a3*r2 - a2*a3*r3)
f(758) = dtau**2*(-a1**2*a4*r4 + a1*a4**2*r1 - a2**2*a3*r3 + a2*a3**2*r2 &
      )
f(759) = dtau**2*(-a1**2*a4*r1 + a1*a4**2*r4 - a2**2*a3*r2 + a2*a3**2*r3 &
      )
f(760) = dtau**2*(a1*a4*r1**2 - a1*a4*r4**2 + a2*a3*r2**2 - a2*a3*r3**2)
f(761) = b1*b2*(a1*r1 + a2*r2 - a3*r3 - a4*r4)
f(762) = b1*b2*(a1*b2**2*r1 + a2*b2**2*r2 - a3*b1**2*r3 - a4*b1**2*r4)
f(763) = b1**2*b2**2*(a1*r1 + a2*r2 - a3*r3 - a4*r4)
f(764) = b1*b2*(-a1*b1**2*r1 - a2*b1**2*r2 + a3*b2**2*r3 + a4*b2**2*r4)
f(765) = b1*b2*(-a1**2*r1 - a2**2*r2 + a3**2*r3 + a4**2*r4)
f(766) = b1*b2*(-a1**3*r1 - a2**3*r2 + a3**3*r3 + a4**3*r4)
f(767) = b1*b2*(a1*r1**2 + a2*r2**2 - a3*r3**2 - a4*r4**2)
f(768) = b1*b2*(-a1**2*r1**2 - a2**2*r2**2 + a3**2*r3**2 + a4**2*r4**2)
f(769) = b1*b2*(-a1*r1**3 - a2*r2**3 + a3*r3**3 + a4*r4**3)
f(770) = dtau*(-a1*b1*r1 + a2*b1*r2 + a3*b2*r3 - a4*b2*r4)
f(771) = dtau**3*(a1*b1*r1 - a2*b1*r2 - a3*b2*r3 + a4*b2*r4)
f(772) = dtau**2*(a1*b1**2*r1 + a2*b1**2*r2 - a3*b2**2*r3 - a4*b2**2*r4)
f(773) = dtau*(a1*b1**3*r1 - a2*b1**3*r2 - a3*b2**3*r3 + a4*b2**3*r4)
f(774) = dtau*(a1**2*b1*r1 - a2**2*b1*r2 - a3**2*b2*r3 + a4**2*b2*r4)
f(775) = dtau*(a1**3*b1*r1 - a2**3*b1*r2 - a3**3*b2*r3 + a4**3*b2*r4)
f(776) = dtau*(a1*b1*r1**2 - a2*b1*r2**2 - a3*b2*r3**2 + a4*b2*r4**2)
f(777) = dtau*(a1**2*b1*r1**2 - a2**2*b1*r2**2 - a3**2*b2*r3**2 + a4**2* &
      b2*r4**2)
f(778) = dtau*(a1*b1*r1**3 - a2*b1*r2**3 - a3*b2*r3**3 + a4*b2*r4**3)
f(779) = dtau*(-a1*b2*r1 + a2*b2*r2 + a3*b1*r3 - a4*b1*r4)
f(780) = dtau**3*(-a1*b2*r1 + a2*b2*r2 + a3*b1*r3 - a4*b1*r4)
f(781) = dtau**2*(a1*b2**2*r1 + a2*b2**2*r2 - a3*b1**2*r3 - a4*b1**2*r4)
f(782) = dtau*(-a1*b2**3*r1 + a2*b2**3*r2 + a3*b1**3*r3 - a4*b1**3*r4)
f(783) = dtau*(-a1**2*b2*r1 + a2**2*b2*r2 + a3**2*b1*r3 - a4**2*b1*r4)
f(784) = dtau*(a1**3*b2*r1 - a2**3*b2*r2 - a3**3*b1*r3 + a4**3*b1*r4)
f(785) = dtau*(a1*b2*r1**2 - a2*b2*r2**2 - a3*b1*r3**2 + a4*b1*r4**2)
f(786) = dtau*(-a1**2*b2*r1**2 + a2**2*b2*r2**2 + a3**2*b1*r3**2 - a4**2 &
      *b1*r4**2)
f(787) = dtau*(-a1*b2*r1**3 + a2*b2*r2**3 + a3*b1*r3**3 - a4*b1*r4**3)
f(788) = a1*a2*a3*r4 + a1*a2*a4*r3 - a1*a3*a4*r2 - a2*a3*a4*r1
f(789) = -a1**2*a2*a3*r4 - a1*a2**2*a4*r3 + a1*a3**2*a4*r2 + a2*a3*a4**2 &
      *r1
f(790) = -a1**3*a2*a3*r4 - a1*a2**3*a4*r3 + a1*a3**3*a4*r2 + a2*a3*a4**3 &
      *r1
f(791) = -a1**2*a2*a4*r3 - a1*a2**2*a3*r4 + a1*a3*a4**2*r2 + a2*a3**2*a4 &
      *r1
f(792) = -a1**2*a2**2*a3*r4 - a1**2*a2**2*a4*r3 + a1*a3**2*a4**2*r2 + a2 &
      *a3**2*a4**2*r1
f(793) = -a1**3*a2*a4*r3 - a1*a2**3*a3*r4 + a1*a3*a4**3*r2 + a2*a3**3*a4 &
      *r1
f(794) = -a1**2*a3*a4*r2 + a1*a2*a3**2*r4 + a1*a2*a4**2*r3 - a2**2*a3*a4 &
      *r1
f(795) = a1**2*a2*a3**2*r4 - a1**2*a3**2*a4*r2 + a1*a2**2*a4**2*r3 - a2 &
      **2*a3*a4**2*r1
f(796) = a1**2*a2*a4**2*r3 - a1**2*a3*a4**2*r2 + a1*a2**2*a3**2*r4 - a2 &
      **2*a3**2*a4*r1
f(797) = a1**3*a3*a4*r2 - a1*a2*a3**3*r4 - a1*a2*a4**3*r3 + a2**3*a3*a4* &
      r1
f(798) = a1*a2*a3*r4**2 + a1*a2*a4*r3**2 - a1*a3*a4*r2**2 - a2*a3*a4*r1 &
      **2
f(799) = a1**2*a2*a3*r4**2 + a1*a2**2*a4*r3**2 - a1*a3**2*a4*r2**2 - a2* &
      a3*a4**2*r1**2
f(800) = -a1**2*a2*a4*r3**2 - a1*a2**2*a3*r4**2 + a1*a3*a4**2*r2**2 + a2 &
      *a3**2*a4*r1**2
f(801) = -a1**2*a3*a4*r2**2 + a1*a2*a3**2*r4**2 + a1*a2*a4**2*r3**2 - a2 &
      **2*a3*a4*r1**2
f(802) = -a1*a2*a3*r4**3 - a1*a2*a4*r3**3 + a1*a3*a4*r2**3 + a2*a3*a4*r1 &
      **3
f(803) = -a1*a4*b1**2*r2 + a1*a4*b2**2*r3 - a2*a3*b1**2*r1 + a2*a3*b2**2 &
      *r4
f(804) = -a1**2*a4*b2**2*r3 + a1*a4**2*b1**2*r2 - a2**2*a3*b2**2*r4 + a2 &
      *a3**2*b1**2*r1
f(805) = a1**2*a4*b1**2*r2 - a1*a4**2*b2**2*r3 + a2**2*a3*b1**2*r1 - a2* &
      a3**2*b2**2*r4
f(806) = a1*a4*b1**2*r2**2 - a1*a4*b2**2*r3**2 + a2*a3*b1**2*r1**2 - a2* &
      a3*b2**2*r4**2
f(807) = -a1*a4*b1**2*r3 + a1*a4*b2**2*r2 - a2*a3*b1**2*r4 + a2*a3*b2**2 &
      *r1
f(808) = -a1**2*a4*b1**2*r3 + a1*a4**2*b2**2*r2 - a2**2*a3*b1**2*r4 + a2 &
      *a3**2*b2**2*r1
f(809) = -a1**2*a4*b2**2*r2 + a1*a4**2*b1**2*r3 - a2**2*a3*b2**2*r1 + a2 &
      *a3**2*b1**2*r4
f(810) = a1*a4*b1**2*r3**2 - a1*a4*b2**2*r2**2 + a2*a3*b1**2*r4**2 - a2* &
      a3*b2**2*r1**2
f(811) = dtau**2*(a1*a4*r2 - a1*a4*r3 + a2*a3*r1 - a2*a3*r4)
f(812) = dtau**2*(-a1**2*a4*r3 + a1*a4**2*r2 - a2**2*a3*r4 + a2*a3**2*r1 &
      )
f(813) = dtau**2*(a1**2*a4*r2 - a1*a4**2*r3 + a2**2*a3*r1 - a2*a3**2*r4)
f(814) = dtau**2*(a1*a4*r2**2 - a1*a4*r3**2 + a2*a3*r1**2 - a2*a3*r4**2)
f(815) = -a1*a3*b1**2*r2 + a1*a3*b2**2*r4 - a2*a4*b1**2*r1 + a2*a4*b2**2 &
      *r3
f(816) = a1**2*a3*b2**2*r4 - a1*a3**2*b1**2*r2 + a2**2*a4*b2**2*r3 - a2* &
      a4**2*b1**2*r1
f(817) = -a1**2*a3*b1**2*r2 + a1*a3**2*b2**2*r4 - a2**2*a4*b1**2*r1 + a2 &
      *a4**2*b2**2*r3
f(818) = -a1*a3*b1**2*r2**2 + a1*a3*b2**2*r4**2 - a2*a4*b1**2*r1**2 + a2 &
      *a4*b2**2*r3**2
f(819) = -a1*a3*b1**2*r4 + a1*a3*b2**2*r2 - a2*a4*b1**2*r3 + a2*a4*b2**2 &
      *r1
f(820) = -a1**2*a3*b1**2*r4 + a1*a3**2*b2**2*r2 - a2**2*a4*b1**2*r3 + a2 &
      *a4**2*b2**2*r1
f(821) = -a1**2*a3*b2**2*r2 + a1*a3**2*b1**2*r4 - a2**2*a4*b2**2*r1 + a2 &
      *a4**2*b1**2*r3
f(822) = -a1*a3*b1**2*r4**2 + a1*a3*b2**2*r2**2 - a2*a4*b1**2*r3**2 + a2 &
      *a4*b2**2*r1**2
f(823) = dtau**2*(a1*a3*r2 - a1*a3*r4 + a2*a4*r1 - a2*a4*r3)
f(824) = dtau**2*(-a1**2*a3*r4 + a1*a3**2*r2 - a2**2*a4*r3 + a2*a4**2*r1 &
      )
f(825) = dtau**2*(a1**2*a3*r2 - a1*a3**2*r4 + a2**2*a4*r1 - a2*a4**2*r3)
f(826) = dtau**2*(a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 - a2*a4*r3**2)
f(827) = b1*b2*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(828) = b1*b2*(a1*b2**2*r2 + a2*b2**2*r1 - a3*b1**2*r4 - a4*b1**2*r3)
f(829) = b1**2*b2**2*(a1*r2 + a2*r1 - a3*r4 - a4*r3)
f(830) = b1*b2*(a1*b1**2*r2 + a2*b1**2*r1 - a3*b2**2*r4 - a4*b2**2*r3)
f(831) = b1*b2*(a1**2*r2 + a2**2*r1 - a3**2*r4 - a4**2*r3)
f(832) = b1*b2*(-a1**3*r2 - a2**3*r1 + a3**3*r4 + a4**3*r3)
f(833) = b1*b2*(a1*r2**2 + a2*r1**2 - a3*r4**2 - a4*r3**2)
f(834) = b1*b2*(a1**2*r2**2 + a2**2*r1**2 - a3**2*r4**2 - a4**2*r3**2)
f(835) = b1*b2*(-a1*r2**3 - a2*r1**3 + a3*r4**3 + a4*r3**3)
f(836) = dtau*(a1*b1*r2 - a2*b1*r1 - a3*b2*r4 + a4*b2*r3)
f(837) = dtau**3*(-a1*b1*r2 + a2*b1*r1 + a3*b2*r4 - a4*b2*r3)
f(838) = dtau**2*(a1*b1**2*r2 + a2*b1**2*r1 - a3*b2**2*r4 - a4*b2**2*r3)
f(839) = dtau*(-a1*b1**3*r2 + a2*b1**3*r1 + a3*b2**3*r4 - a4*b2**3*r3)
f(840) = dtau*(-a1**2*b1*r2 + a2**2*b1*r1 + a3**2*b2*r4 - a4**2*b2*r3)
f(841) = dtau*(-a1**3*b1*r2 + a2**3*b1*r1 + a3**3*b2*r4 - a4**3*b2*r3)
f(842) = dtau*(-a1*b1*r2**2 + a2*b1*r1**2 + a3*b2*r4**2 - a4*b2*r3**2)
f(843) = dtau*(a1**2*b1*r2**2 - a2**2*b1*r1**2 - a3**2*b2*r4**2 + a4**2* &
      b2*r3**2)
f(844) = dtau*(-a1*b1*r2**3 + a2*b1*r1**3 + a3*b2*r4**3 - a4*b2*r3**3)
f(845) = dtau*(a1*b2*r2 - a2*b2*r1 - a3*b1*r4 + a4*b1*r3)
f(846) = dtau**3*(a1*b2*r2 - a2*b2*r1 - a3*b1*r4 + a4*b1*r3)
f(847) = dtau**2*(a1*b2**2*r2 + a2*b2**2*r1 - a3*b1**2*r4 - a4*b1**2*r3)
f(848) = dtau*(a1*b2**3*r2 - a2*b2**3*r1 - a3*b1**3*r4 + a4*b1**3*r3)
f(849) = dtau*(-a1**2*b2*r2 + a2**2*b2*r1 + a3**2*b1*r4 - a4**2*b1*r3)
f(850) = dtau*(a1**3*b2*r2 - a2**3*b2*r1 - a3**3*b1*r4 + a4**3*b1*r3)
f(851) = dtau*(-a1*b2*r2**2 + a2*b2*r1**2 + a3*b1*r4**2 - a4*b1*r3**2)
f(852) = dtau*(a1**2*b2*r2**2 - a2**2*b2*r1**2 - a3**2*b1*r4**2 + a4**2* &
      b1*r3**2)
f(853) = dtau*(a1*b2*r2**3 - a2*b2*r1**3 - a3*b1*r4**3 + a4*b1*r3**3)
f(854) = -a1*a2*b2**2*r3 - a1*a2*b2**2*r4 + a3*a4*b1**2*r1 + a3*a4*b1**2 &
      *r2
f(855) = -a1**2*a2*b2**2*r4 - a1*a2**2*b2**2*r3 + a3**2*a4*b1**2*r2 + a3 &
      *a4**2*b1**2*r1
f(856) = a1**2*a2*b2**2*r3 + a1*a2**2*b2**2*r4 - a3**2*a4*b1**2*r1 - a3* &
      a4**2*b1**2*r2
f(857) = a1*a2*b2**2*r3**2 + a1*a2*b2**2*r4**2 - a3*a4*b1**2*r1**2 - a3* &
      a4*b1**2*r2**2
f(858) = -a1*a2*b1**2*r3 - a1*a2*b1**2*r4 + a3*a4*b2**2*r1 + a3*a4*b2**2 &
      *r2
f(859) = a1**2*a2*b1**2*r4 + a1*a2**2*b1**2*r3 - a3**2*a4*b2**2*r2 - a3* &
      a4**2*b2**2*r1
f(860) = a1**2*a2*b1**2*r3 + a1*a2**2*b1**2*r4 - a3**2*a4*b2**2*r1 - a3* &
      a4**2*b2**2*r2
f(861) = a1*a2*b1**2*r3**2 + a1*a2*b1**2*r4**2 - a3*a4*b2**2*r1**2 - a3* &
      a4*b2**2*r2**2
f(862) = dtau**2*(-a1*a2*r3 - a1*a2*r4 + a3*a4*r1 + a3*a4*r2)
f(863) = dtau**2*(a1**2*a2*r4 + a1*a2**2*r3 - a3**2*a4*r2 - a3*a4**2*r1)
f(864) = dtau**2*(a1**2*a2*r3 + a1*a2**2*r4 - a3**2*a4*r1 - a3*a4**2*r2)
f(865) = dtau**2*(-a1*a2*r3**2 - a1*a2*r4**2 + a3*a4*r1**2 + a3*a4*r2**2 &
      )
f(866) = b1*b2*(a1*r3 + a2*r4 - a3*r1 - a4*r2)
f(867) = b1*b2*(a1*b1**2*r3 + a2*b1**2*r4 - a3*b2**2*r1 - a4*b2**2*r2)
f(868) = b1**2*b2**2*(-a1*r3 - a2*r4 + a3*r1 + a4*r2)
f(869) = b1*b2*(a1*b2**2*r3 + a2*b2**2*r4 - a3*b1**2*r1 - a4*b1**2*r2)
f(870) = b1*b2*(a1**2*r3 + a2**2*r4 - a3**2*r1 - a4**2*r2)
f(871) = b1*b2*(a1**3*r3 + a2**3*r4 - a3**3*r1 - a4**3*r2)
f(872) = b1*b2*(a1*r3**2 + a2*r4**2 - a3*r1**2 - a4*r2**2)
f(873) = b1*b2*(a1**2*r3**2 + a2**2*r4**2 - a3**2*r1**2 - a4**2*r2**2)
f(874) = b1*b2*(a1*r3**3 + a2*r4**3 - a3*r1**3 - a4*r2**3)
f(875) = dtau*(a1*b2*r3 - a2*b2*r4 - a3*b1*r1 + a4*b1*r2)
f(876) = dtau**3*(a1*b2*r3 - a2*b2*r4 - a3*b1*r1 + a4*b1*r2)
f(877) = dtau**2*(a1*b2**2*r3 + a2*b2**2*r4 - a3*b1**2*r1 - a4*b1**2*r2)
f(878) = dtau*(a1*b2**3*r3 - a2*b2**3*r4 - a3*b1**3*r1 + a4*b1**3*r2)
f(879) = dtau*(a1**2*b2*r3 - a2**2*b2*r4 - a3**2*b1*r1 + a4**2*b1*r2)
f(880) = dtau*(-a1**3*b2*r3 + a2**3*b2*r4 + a3**3*b1*r1 - a4**3*b1*r2)
f(881) = dtau*(-a1*b2*r3**2 + a2*b2*r4**2 + a3*b1*r1**2 - a4*b1*r2**2)
f(882) = dtau*(a1**2*b2*r3**2 - a2**2*b2*r4**2 - a3**2*b1*r1**2 + a4**2* &
      b1*r2**2)
f(883) = dtau*(a1*b2*r3**3 - a2*b2*r4**3 - a3*b1*r1**3 + a4*b1*r2**3)
f(884) = dtau*(a1*b1*r3 - a2*b1*r4 - a3*b2*r1 + a4*b2*r2)
f(885) = dtau**3*(-a1*b1*r3 + a2*b1*r4 + a3*b2*r1 - a4*b2*r2)
f(886) = dtau**2*(a1*b1**2*r3 + a2*b1**2*r4 - a3*b2**2*r1 - a4*b2**2*r2)
f(887) = dtau*(-a1*b1**3*r3 + a2*b1**3*r4 + a3*b2**3*r1 - a4*b2**3*r2)
f(888) = dtau*(-a1**2*b1*r3 + a2**2*b1*r4 + a3**2*b2*r1 - a4**2*b2*r2)
f(889) = dtau*(-a1**3*b1*r3 + a2**3*b1*r4 + a3**3*b2*r1 - a4**3*b2*r2)
f(890) = dtau*(-a1*b1*r3**2 + a2*b1*r4**2 + a3*b2*r1**2 - a4*b2*r2**2)
f(891) = dtau*(a1**2*b1*r3**2 - a2**2*b1*r4**2 - a3**2*b2*r1**2 + a4**2* &
      b2*r2**2)
f(892) = dtau*(-a1*b1*r3**3 + a2*b1*r4**3 + a3*b2*r1**3 - a4*b2*r2**3)
f(893) = b1*b2*(a1*r4 + a2*r3 - a3*r2 - a4*r1)
f(894) = b1*b2*(a1*b1**2*r4 + a2*b1**2*r3 - a3*b2**2*r2 - a4*b2**2*r1)
f(895) = b1**2*b2**2*(-a1*r4 - a2*r3 + a3*r2 + a4*r1)
f(896) = b1*b2*(a1*b2**2*r4 + a2*b2**2*r3 - a3*b1**2*r2 - a4*b1**2*r1)
f(897) = b1*b2*(a1**2*r4 + a2**2*r3 - a3**2*r2 - a4**2*r1)
f(898) = b1*b2*(a1**3*r4 + a2**3*r3 - a3**3*r2 - a4**3*r1)
f(899) = b1*b2*(a1*r4**2 + a2*r3**2 - a3*r2**2 - a4*r1**2)
f(900) = b1*b2*(a1**2*r4**2 + a2**2*r3**2 - a3**2*r2**2 - a4**2*r1**2)
f(901) = b1*b2*(a1*r4**3 + a2*r3**3 - a3*r2**3 - a4*r1**3)
f(902) = dtau*(-a1*b2*r4 + a2*b2*r3 + a3*b1*r2 - a4*b1*r1)
f(903) = dtau**3*(-a1*b2*r4 + a2*b2*r3 + a3*b1*r2 - a4*b1*r1)
f(904) = dtau**2*(a1*b2**2*r4 + a2*b2**2*r3 - a3*b1**2*r2 - a4*b1**2*r1)
f(905) = dtau*(-a1*b2**3*r4 + a2*b2**3*r3 + a3*b1**3*r2 - a4*b1**3*r1)
f(906) = dtau*(a1**2*b2*r4 - a2**2*b2*r3 - a3**2*b1*r2 + a4**2*b1*r1)
f(907) = dtau*(-a1**3*b2*r4 + a2**3*b2*r3 + a3**3*b1*r2 - a4**3*b1*r1)
f(908) = dtau*(a1*b2*r4**2 - a2*b2*r3**2 - a3*b1*r2**2 + a4*b1*r1**2)
f(909) = dtau*(a1**2*b2*r4**2 - a2**2*b2*r3**2 - a3**2*b1*r2**2 + a4**2* &
      b1*r1**2)
f(910) = dtau*(-a1*b2*r4**3 + a2*b2*r3**3 + a3*b1*r2**3 - a4*b1*r1**3)
f(911) = dtau*(-a1*b1*r4 + a2*b1*r3 + a3*b2*r2 - a4*b2*r1)
f(912) = dtau**3*(a1*b1*r4 - a2*b1*r3 - a3*b2*r2 + a4*b2*r1)
f(913) = dtau**2*(-a1*b1**2*r4 - a2*b1**2*r3 + a3*b2**2*r2 + a4*b2**2*r1 &
      )
f(914) = dtau*(a1*b1**3*r4 - a2*b1**3*r3 - a3*b2**3*r2 + a4*b2**3*r1)
f(915) = dtau*(a1**2*b1*r4 - a2**2*b1*r3 - a3**2*b2*r2 + a4**2*b2*r1)
f(916) = dtau*(a1**3*b1*r4 - a2**3*b1*r3 - a3**3*b2*r2 + a4**3*b2*r1)
f(917) = dtau*(a1*b1*r4**2 - a2*b1*r3**2 - a3*b2*r2**2 + a4*b2*r1**2)
f(918) = dtau*(-a1**2*b1*r4**2 + a2**2*b1*r3**2 + a3**2*b2*r2**2 - a4**2 &
      *b2*r1**2)
f(919) = dtau*(a1*b1*r4**3 - a2*b1*r3**3 - a3*b2*r2**3 + a4*b2*r1**3)
f(920) = b1*b2*dtau**2*(r1 + r2 - r3 - r4)
f(921) = b1*b2*dtau*(-b1*r3 + b1*r4 + b2*r1 - b2*r2)
f(922) = b1*b2*dtau*(b1*r1 - b1*r2 - b2*r3 + b2*r4)
f(923) = b1*b2*dtau**2*(r1**2 + r2**2 - r3**2 - r4**2)
f(924) = b1*b2*dtau*(-b1*r3**2 + b1*r4**2 + b2*r1**2 - b2*r2**2)
f(925) = b1*b2*dtau*(-b1*r1**2 + b1*r2**2 + b2*r3**2 - b2*r4**2)
f(926) = a1*a2*a3*a4*(-a1 - a2 + a3 + a4)
f(927) = a1*a2*a3*a4*(a1**2 + a2**2 - a3**2 - a4**2)
f(928) = a1*a2*a3*a4*(-a1*a2 + a3*a4)
f(929) = a1*a2*a3*b1**2 + a1*a2*a4*b1**2 - a1*a3*a4*b2**2 - a2*a3*a4*b2 &
      **2
f(930) = -a1**2*a3*a4*b2**2 + a1*a2*a3**2*b1**2 + a1*a2*a4**2*b1**2 - a2 &
      **2*a3*a4*b2**2
f(931) = -a1**2*a2*a4*b1**2 - a1*a2**2*a3*b1**2 + a1*a3*a4**2*b2**2 + a2 &
      *a3**2*a4*b2**2
f(932) = -a1**2*a2*a3*b1**2 - a1*a2**2*a4*b1**2 + a1*a3**2*a4*b2**2 + a2 &
      *a3*a4**2*b2**2
f(933) = -a1*a2*a3*b2**2 - a1*a2*a4*b2**2 + a1*a3*a4*b1**2 + a2*a3*a4*b1 &
      **2
f(934) = a1**2*a3*a4*b1**2 - a1*a2*a3**2*b2**2 - a1*a2*a4**2*b2**2 + a2 &
      **2*a3*a4*b1**2
f(935) = a1**2*a2*a4*b2**2 + a1*a2**2*a3*b2**2 - a1*a3*a4**2*b1**2 - a2* &
      a3**2*a4*b1**2
f(936) = -a1**2*a2*a3*b2**2 - a1*a2**2*a4*b2**2 + a1*a3**2*a4*b1**2 + a2 &
      *a3*a4**2*b1**2
f(937) = dtau**2*(a1*a2*a3 + a1*a2*a4 - a1*a3*a4 - a2*a3*a4)
f(938) = dtau**2*(a1**2*a3*a4 - a1*a2*a3**2 - a1*a2*a4**2 + a2**2*a3*a4)
f(939) = dtau**2*(-a1**2*a2*a4 - a1*a2**2*a3 + a1*a3*a4**2 + a2*a3**2*a4 &
      )
f(940) = dtau**2*(-a1**2*a2*a3 - a1*a2**2*a4 + a1*a3**2*a4 + a2*a3*a4**2 &
      )
f(941) = b1*b2*(-a1*a2 + a3*a4)
f(942) = b1*b2*(-a1*a2*b2**2 + a3*a4*b1**2)
f(943) = b1**2*b2**2*(a1*a2 - a3*a4)
f(944) = b1*b2*(-a1*a2*b1**2 + a3*a4*b2**2)
f(945) = b1*b2*(-a1**2*a2 - a1*a2**2 + a3**2*a4 + a3*a4**2)
f(946) = b1*b2*(-a1**3*a2 - a1*a2**3 + a3**3*a4 + a3*a4**3)
f(947) = b1*b2*(-a1**2*a2**2 + a3**2*a4**2)
f(948) = dtau**2*(-a1*a2*b1**2 + a3*a4*b2**2)
f(949) = dtau*(-a1**2*a2*b1 + a1*a2**2*b1 + a3**2*a4*b2 - a3*a4**2*b2)
f(950) = dtau*(a1**3*a2*b1 - a1*a2**3*b1 - a3**3*a4*b2 + a3*a4**3*b2)
f(951) = dtau**2*(-a1*a2*b2**2 + a3*a4*b1**2)
f(952) = dtau*(a1**2*a2*b2 - a1*a2**2*b2 - a3**2*a4*b1 + a3*a4**2*b1)
f(953) = dtau*(-a1**3*a2*b2 + a1*a2**3*b2 + a3**3*a4*b1 - a3*a4**3*b1)
f(954) = b1*b2*(-a1*a3*b1**2 + a1*a3*b2**2 - a2*a4*b1**2 + a2*a4*b2**2)
f(955) = b1*b2*(a1**2*a3 - a1*a3**2 + a2**2*a4 - a2*a4**2)
f(956) = b1*b2*(-a1**3*a3 + a1*a3**3 - a2**3*a4 + a2*a4**3)
f(957) = dtau*(-a1*a3*b1 + a1*a3*b2 + a2*a4*b1 - a2*a4*b2)
f(958) = dtau**3*(-a1*a3*b1 + a1*a3*b2 + a2*a4*b1 - a2*a4*b2)
f(959) = dtau**2*(-a1*a3*b1**2 + a1*a3*b2**2 - a2*a4*b1**2 + a2*a4*b2**2 &
      )
f(960) = dtau*(-a1*a3*b1**3 + a1*a3*b2**3 + a2*a4*b1**3 - a2*a4*b2**3)
f(961) = dtau*(-a1**2*a3*b2 + a1*a3**2*b1 + a2**2*a4*b2 - a2*a4**2*b1)
f(962) = dtau*(a1**3*a3*b2 - a1*a3**3*b1 - a2**3*a4*b2 + a2*a4**3*b1)
f(963) = dtau*(a1**2*a3*b1 - a1*a3**2*b2 - a2**2*a4*b1 + a2*a4**2*b2)
f(964) = dtau*(a1**2*a3**2*b1 - a1**2*a3**2*b2 - a2**2*a4**2*b1 + a2**2* &
      a4**2*b2)
f(965) = dtau*(-a1**3*a3*b1 + a1*a3**3*b2 + a2**3*a4*b1 - a2*a4**3*b2)
f(966) = b1*b2*(a1*a4*b1**2 - a1*a4*b2**2 + a2*a3*b1**2 - a2*a3*b2**2)
f(967) = b1*b2*(-a1**2*a4 + a1*a4**2 - a2**2*a3 + a2*a3**2)
f(968) = b1*b2*(a1**3*a4 - a1*a4**3 + a2**3*a3 - a2*a3**3)
f(969) = dtau*(a1*a4*b1 + a1*a4*b2 - a2*a3*b1 - a2*a3*b2)
f(970) = dtau**3*(a1*a4*b1 + a1*a4*b2 - a2*a3*b1 - a2*a3*b2)
f(971) = dtau**2*(a1*a4*b1**2 - a1*a4*b2**2 + a2*a3*b1**2 - a2*a3*b2**2)
f(972) = dtau*(a1*a4*b1**3 + a1*a4*b2**3 - a2*a3*b1**3 - a2*a3*b2**3)
f(973) = dtau*(-a1**2*a4*b2 - a1*a4**2*b1 + a2**2*a3*b2 + a2*a3**2*b1)
f(974) = dtau*(a1**3*a4*b2 + a1*a4**3*b1 - a2**3*a3*b2 - a2*a3**3*b1)
f(975) = dtau*(a1**2*a4*b1 + a1*a4**2*b2 - a2**2*a3*b1 - a2*a3**2*b2)
f(976) = dtau*(a1**2*a4**2*b1 + a1**2*a4**2*b2 - a2**2*a3**2*b1 - a2**2* &
      a3**2*b2)
f(977) = dtau*(-a1**3*a4*b1 - a1*a4**3*b2 + a2**3*a3*b1 + a2*a3**3*b2)
f(978) = b1*b2*dtau**2*(-a1 - a2 + a3 + a4)
f(979) = b1*b2*dtau*(a1*b2 - a2*b2 - a3*b1 + a4*b1)
f(980) = b1*b2*dtau*(a1*b1 - a2*b1 - a3*b2 + a4*b2)
f(981) = b1*b2*dtau**2*(-a1**2 - a2**2 + a3**2 + a4**2)
f(982) = b1*b2*dtau*(-a1**2*b2 + a2**2*b2 + a3**2*b1 - a4**2*b1)
f(983) = b1*b2*dtau*(a1**2*b1 - a2**2*b1 - a3**2*b2 + a4**2*b2)
v = sum(f*params)
end function c2h4_dipole_b1u_n4_d6_ADF


!###############################################################################


function c2h4_dipole_b2u_n1_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(12)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(12)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r1 - r2 + r3 - r4
f(2) = r1**2 - r2**2 + r3**2 - r4**2
f(3) = r1**3 - r2**3 + r3**3 - r4**3
f(4) = -r1**4 + r2**4 - r3**4 + r4**4
f(5) = -r1**5 + r2**5 - r3**5 + r4**5
f(6) = -r1**6 + r2**6 - r3**6 + r4**6
f(7) = a1 - a2 + a3 - a4
f(8) = -a1**2 + a2**2 - a3**2 + a4**2
f(9) = a1**3 - a2**3 + a3**3 - a4**3
f(10) = -a1**4 + a2**4 - a3**4 + a4**4
f(11) = a1**5 - a2**5 + a3**5 - a4**5
f(12) = a1**6 - a2**6 + a3**6 - a4**6
v = sum(f*params)
end function c2h4_dipole_b2u_n1_d6_ADF


!###############################################################################


function c2h4_dipole_b2u_n2_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(177)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(177)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(r1 - r2 + r3 - r4)
f(2) = r0*(r1**2 - r2**2 + r3**2 - r4**2)
f(3) = r0*(r1**3 - r2**3 + r3**3 - r4**3)
f(4) = r0*(-r1**4 + r2**4 - r3**4 + r4**4)
f(5) = r0*(-r1**5 + r2**5 - r3**5 + r4**5)
f(6) = r0**2*(r1 - r2 + r3 - r4)
f(7) = r0**2*(-r1**2 + r2**2 - r3**2 + r4**2)
f(8) = r0**2*(r1**3 - r2**3 + r3**3 - r4**3)
f(9) = r0**2*(-r1**4 + r2**4 - r3**4 + r4**4)
f(10) = r0**3*(-r1 + r2 - r3 + r4)
f(11) = r0**3*(r1**2 - r2**2 + r3**2 - r4**2)
f(12) = r0**3*(r1**3 - r2**3 + r3**3 - r4**3)
f(13) = r0**4*(-r1 + r2 - r3 + r4)
f(14) = r0**4*(-r1**2 + r2**2 - r3**2 + r4**2)
f(15) = r0**5*(-r1 + r2 - r3 + r4)
f(16) = r0*(-a1 + a2 - a3 + a4)
f(17) = r0*(a1**2 - a2**2 + a3**2 - a4**2)
f(18) = r0*(-a1**3 + a2**3 - a3**3 + a4**3)
f(19) = r0*(a1**4 - a2**4 + a3**4 - a4**4)
f(20) = r0*(-a1**5 + a2**5 - a3**5 + a4**5)
f(21) = r0**2*(a1 - a2 + a3 - a4)
f(22) = r0**2*(-a1**2 + a2**2 - a3**2 + a4**2)
f(23) = r0**2*(a1**3 - a2**3 + a3**3 - a4**3)
f(24) = r0**2*(-a1**4 + a2**4 - a3**4 + a4**4)
f(25) = r0**3*(-a1 + a2 - a3 + a4)
f(26) = r0**3*(a1**2 - a2**2 + a3**2 - a4**2)
f(27) = r0**3*(-a1**3 + a2**3 - a3**3 + a4**3)
f(28) = r0**4*(a1 - a2 + a3 - a4)
f(29) = r0**4*(-a1**2 + a2**2 - a3**2 + a4**2)
f(30) = r0**5*(a1 - a2 + a3 - a4)
f(31) = r1**2*r2 - r1*r2**2 + r3**2*r4 - r3*r4**2
f(32) = -r1**3*r2 + r1*r2**3 - r3**3*r4 + r3*r4**3
f(33) = r1**4*r2 - r1*r2**4 + r3**4*r4 - r3*r4**4
f(34) = r1**5*r2 - r1*r2**5 + r3**5*r4 - r3*r4**5
f(35) = -r1**2*r2 + r1*r2**2 - r3**2*r4 + r3*r4**2
f(36) = r1**3*r2**2 - r1**2*r2**3 + r3**3*r4**2 - r3**2*r4**3
f(37) = r1**4*r2**2 - r1**2*r2**4 + r3**4*r4**2 - r3**2*r4**4
f(38) = -r1**4*r2**2 + r1**2*r2**4 - r3**4*r4**2 + r3**2*r4**4
f(39) = r1*r3 - r2*r4
f(40) = -r1**2*r3 - r1*r3**2 + r2**2*r4 + r2*r4**2
f(41) = -r1**3*r3 - r1*r3**3 + r2**3*r4 + r2*r4**3
f(42) = r1**4*r3 + r1*r3**4 - r2**4*r4 - r2*r4**4
f(43) = -r1**5*r3 - r1*r3**5 + r2**5*r4 + r2*r4**5
f(44) = -r1**2*r3**2 + r2**2*r4**2
f(45) = -r1**3*r3**2 - r1**2*r3**3 + r2**3*r4**2 + r2**2*r4**3
f(46) = -r1**4*r3**2 - r1**2*r3**4 + r2**4*r4**2 + r2**2*r4**4
f(47) = r1**3*r3**3 - r2**3*r4**3
f(48) = r1**2*r4 - r1*r4**2 - r2**2*r3 + r2*r3**2
f(49) = r1**3*r4 - r1*r4**3 - r2**3*r3 + r2*r3**3
f(50) = -r1**4*r4 + r1*r4**4 + r2**4*r3 - r2*r3**4
f(51) = -r1**5*r4 + r1*r4**5 + r2**5*r3 - r2*r3**5
f(52) = -r1**3*r4**2 + r1**2*r4**3 + r2**3*r3**2 - r2**2*r3**3
f(53) = r1**4*r4**2 - r1**2*r4**4 - r2**4*r3**2 + r2**2*r3**4
f(54) = -a1*r1 + a2*r2 - a3*r3 + a4*r4
f(55) = -a1**2*r1 + a2**2*r2 - a3**2*r3 + a4**2*r4
f(56) = -a1**3*r1 + a2**3*r2 - a3**3*r3 + a4**3*r4
f(57) = a1**4*r1 - a2**4*r2 + a3**4*r3 - a4**4*r4
f(58) = a1**5*r1 - a2**5*r2 + a3**5*r3 - a4**5*r4
f(59) = -a1*r1**2 + a2*r2**2 - a3*r3**2 + a4*r4**2
f(60) = -a1**2*r1**2 + a2**2*r2**2 - a3**2*r3**2 + a4**2*r4**2
f(61) = a1**3*r1**2 - a2**3*r2**2 + a3**3*r3**2 - a4**3*r4**2
f(62) = a1**4*r1**2 - a2**4*r2**2 + a3**4*r3**2 - a4**4*r4**2
f(63) = a1*r1**3 - a2*r2**3 + a3*r3**3 - a4*r4**3
f(64) = -a1**2*r1**3 + a2**2*r2**3 - a3**2*r3**3 + a4**2*r4**3
f(65) = a1**3*r1**3 - a2**3*r2**3 + a3**3*r3**3 - a4**3*r4**3
f(66) = -a1*r1**4 + a2*r2**4 - a3*r3**4 + a4*r4**4
f(67) = -a1**2*r1**4 + a2**2*r2**4 - a3**2*r3**4 + a4**2*r4**4
f(68) = -a1*r1**5 + a2*r2**5 - a3*r3**5 + a4*r4**5
f(69) = -a1*r2 + a2*r1 - a3*r4 + a4*r3
f(70) = a1**2*r2 - a2**2*r1 + a3**2*r4 - a4**2*r3
f(71) = a1**3*r2 - a2**3*r1 + a3**3*r4 - a4**3*r3
f(72) = -a1**4*r2 + a2**4*r1 - a3**4*r4 + a4**4*r3
f(73) = -a1**5*r2 + a2**5*r1 - a3**5*r4 + a4**5*r3
f(74) = a1*r2**2 - a2*r1**2 + a3*r4**2 - a4*r3**2
f(75) = -a1**2*r2**2 + a2**2*r1**2 - a3**2*r4**2 + a4**2*r3**2
f(76) = a1**3*r2**2 - a2**3*r1**2 + a3**3*r4**2 - a4**3*r3**2
f(77) = a1**4*r2**2 - a2**4*r1**2 + a3**4*r4**2 - a4**4*r3**2
f(78) = a1*r2**3 - a2*r1**3 + a3*r4**3 - a4*r3**3
f(79) = -a1**2*r2**3 + a2**2*r1**3 - a3**2*r4**3 + a4**2*r3**3
f(80) = a1**3*r2**3 - a2**3*r1**3 + a3**3*r4**3 - a4**3*r3**3
f(81) = -a1*r2**4 + a2*r1**4 - a3*r4**4 + a4*r3**4
f(82) = a1**2*r2**4 - a2**2*r1**4 + a3**2*r4**4 - a4**2*r3**4
f(83) = -a1*r2**5 + a2*r1**5 - a3*r4**5 + a4*r3**5
f(84) = -a1*r3 + a2*r4 - a3*r1 + a4*r2
f(85) = a1**2*r3 - a2**2*r4 + a3**2*r1 - a4**2*r2
f(86) = a1**3*r3 - a2**3*r4 + a3**3*r1 - a4**3*r2
f(87) = a1**4*r3 - a2**4*r4 + a3**4*r1 - a4**4*r2
f(88) = -a1**5*r3 + a2**5*r4 - a3**5*r1 + a4**5*r2
f(89) = -a1*r3**2 + a2*r4**2 - a3*r1**2 + a4*r2**2
f(90) = a1**2*r3**2 - a2**2*r4**2 + a3**2*r1**2 - a4**2*r2**2
f(91) = -a1**3*r3**2 + a2**3*r4**2 - a3**3*r1**2 + a4**3*r2**2
f(92) = -a1**4*r3**2 + a2**4*r4**2 - a3**4*r1**2 + a4**4*r2**2
f(93) = -a1*r3**3 + a2*r4**3 - a3*r1**3 + a4*r2**3
f(94) = a1**2*r3**3 - a2**2*r4**3 + a3**2*r1**3 - a4**2*r2**3
f(95) = a1**3*r3**3 - a2**3*r4**3 + a3**3*r1**3 - a4**3*r2**3
f(96) = a1*r3**4 - a2*r4**4 + a3*r1**4 - a4*r2**4
f(97) = -a1**2*r3**4 + a2**2*r4**4 - a3**2*r1**4 + a4**2*r2**4
f(98) = a1*r3**5 - a2*r4**5 + a3*r1**5 - a4*r2**5
f(99) = a1*r4 - a2*r3 + a3*r2 - a4*r1
f(100) = -a1**2*r4 + a2**2*r3 - a3**2*r2 + a4**2*r1
f(101) = -a1**3*r4 + a2**3*r3 - a3**3*r2 + a4**3*r1
f(102) = a1**4*r4 - a2**4*r3 + a3**4*r2 - a4**4*r1
f(103) = -a1**5*r4 + a2**5*r3 - a3**5*r2 + a4**5*r1
f(104) = -a1*r4**2 + a2*r3**2 - a3*r2**2 + a4*r1**2
f(105) = a1**2*r4**2 - a2**2*r3**2 + a3**2*r2**2 - a4**2*r1**2
f(106) = -a1**3*r4**2 + a2**3*r3**2 - a3**3*r2**2 + a4**3*r1**2
f(107) = -a1**4*r4**2 + a2**4*r3**2 - a3**4*r2**2 + a4**4*r1**2
f(108) = -a1*r4**3 + a2*r3**3 - a3*r2**3 + a4*r1**3
f(109) = -a1**2*r4**3 + a2**2*r3**3 - a3**2*r2**3 + a4**2*r1**3
f(110) = -a1**3*r4**3 + a2**3*r3**3 - a3**3*r2**3 + a4**3*r1**3
f(111) = a1*r4**4 - a2*r3**4 + a3*r2**4 - a4*r1**4
f(112) = -a1**2*r4**4 + a2**2*r3**4 - a3**2*r2**4 + a4**2*r1**4
f(113) = -a1*r4**5 + a2*r3**5 - a3*r2**5 + a4*r1**5
f(114) = -b1**2*r1 + b1**2*r2 - b2**2*r3 + b2**2*r4
f(115) = b1**4*r1 - b1**4*r2 + b2**4*r3 - b2**4*r4
f(116) = -b1**2*r1**2 + b1**2*r2**2 - b2**2*r3**2 + b2**2*r4**2
f(117) = -b1**4*r1**2 + b1**4*r2**2 - b2**4*r3**2 + b2**4*r4**2
f(118) = -b1**2*r1**3 + b1**2*r2**3 - b2**2*r3**3 + b2**2*r4**3
f(119) = b1**2*r1**4 - b1**2*r2**4 + b2**2*r3**4 - b2**2*r4**4
f(120) = b1**2*r3 - b1**2*r4 + b2**2*r1 - b2**2*r2
f(121) = b1**4*r3 - b1**4*r4 + b2**4*r1 - b2**4*r2
f(122) = -b1**2*r3**2 + b1**2*r4**2 - b2**2*r1**2 + b2**2*r2**2
f(123) = b1**4*r3**2 - b1**4*r4**2 + b2**4*r1**2 - b2**4*r2**2
f(124) = -b1**2*r3**3 + b1**2*r4**3 - b2**2*r1**3 + b2**2*r2**3
f(125) = b1**2*r3**4 - b1**2*r4**4 + b2**2*r1**4 - b2**2*r2**4
f(126) = dtau**2*(r1 - r2 + r3 - r4)
f(127) = dtau**4*(-r1 + r2 - r3 + r4)
f(128) = dtau**2*(r1**2 - r2**2 + r3**2 - r4**2)
f(129) = dtau**4*(-r1**2 + r2**2 - r3**2 + r4**2)
f(130) = dtau**2*(r1**3 - r2**3 + r3**3 - r4**3)
f(131) = dtau**2*(-r1**4 + r2**4 - r3**4 + r4**4)
f(132) = a1**2*a2 - a1*a2**2 + a3**2*a4 - a3*a4**2
f(133) = a1**3*a2 - a1*a2**3 + a3**3*a4 - a3*a4**3
f(134) = -a1**4*a2 + a1*a2**4 - a3**4*a4 + a3*a4**4
f(135) = -a1**5*a2 + a1*a2**5 - a3**5*a4 + a3*a4**5
f(136) = -a1**3*a2**2 + a1**2*a2**3 - a3**3*a4**2 + a3**2*a4**3
f(137) = a1**4*a2**2 - a1**2*a2**4 + a3**4*a4**2 - a3**2*a4**4
f(138) = -a1*a3 + a2*a4
f(139) = a1**2*a3 + a1*a3**2 - a2**2*a4 - a2*a4**2
f(140) = a1**3*a3 + a1*a3**3 - a2**3*a4 - a2*a4**3
f(141) = -a1**4*a3 - a1*a3**4 + a2**4*a4 + a2*a4**4
f(142) = -a1**5*a3 - a1*a3**5 + a2**5*a4 + a2*a4**5
f(143) = -a1**2*a3**2 + a2**2*a4**2
f(144) = a1**3*a3**2 + a1**2*a3**3 - a2**3*a4**2 - a2**2*a4**3
f(145) = a1**4*a3**2 + a1**2*a3**4 - a2**4*a4**2 - a2**2*a4**4
f(146) = -a1**3*a3**2 - a1**2*a3**3 + a2**3*a4**2 + a2**2*a4**3
f(147) = a1**3*a3**3 - a2**3*a4**3
f(148) = -a1**2*a4 + a1*a4**2 + a2**2*a3 - a2*a3**2
f(149) = -a1**3*a4 + a1*a4**3 + a2**3*a3 - a2*a3**3
f(150) = a1**4*a4 - a1*a4**4 - a2**4*a3 + a2*a3**4
f(151) = a1**5*a4 - a1*a4**5 - a2**5*a3 + a2*a3**5
f(152) = a1**3*a4**2 - a1**2*a4**3 - a2**3*a3**2 + a2**2*a3**3
f(153) = a1**4*a4**2 - a1**2*a4**4 - a2**4*a3**2 + a2**2*a3**4
f(154) = -a1*b1**2 + a2*b1**2 - a3*b2**2 + a4*b2**2
f(155) = a1*b1**4 - a2*b1**4 + a3*b2**4 - a4*b2**4
f(156) = a1**2*b1**2 - a2**2*b1**2 + a3**2*b2**2 - a4**2*b2**2
f(157) = -a1**2*b1**4 + a2**2*b1**4 - a3**2*b2**4 + a4**2*b2**4
f(158) = -a1**3*b1**2 + a2**3*b1**2 - a3**3*b2**2 + a4**3*b2**2
f(159) = -a1**4*b1**2 + a2**4*b1**2 - a3**4*b2**2 + a4**4*b2**2
f(160) = a1*b2**2 - a2*b2**2 + a3*b1**2 - a4*b1**2
f(161) = a1*b2**4 - a2*b2**4 + a3*b1**4 - a4*b1**4
f(162) = -a1**2*b2**2 + a2**2*b2**2 - a3**2*b1**2 + a4**2*b1**2
f(163) = -a1**2*b2**4 + a2**2*b2**4 - a3**2*b1**4 + a4**2*b1**4
f(164) = -a1**3*b2**2 + a2**3*b2**2 - a3**3*b1**2 + a4**3*b1**2
f(165) = -a1**4*b2**2 + a2**4*b2**2 - a3**4*b1**2 + a4**4*b1**2
f(166) = dtau**2*(-a1 + a2 - a3 + a4)
f(167) = dtau**4*(-a1 + a2 - a3 + a4)
f(168) = dtau**2*(a1**2 - a2**2 + a3**2 - a4**2)
f(169) = dtau**4*(a1**2 - a2**2 + a3**2 - a4**2)
f(170) = dtau**2*(-a1**3 + a2**3 - a3**3 + a4**3)
f(171) = dtau**2*(a1**4 - a2**4 + a3**4 - a4**4)
f(172) = dtau*(b1 + b2)
f(173) = dtau**3*(b1 + b2)
f(174) = dtau**5*(b1 + b2)
f(175) = dtau*(b1**3 + b2**3)
f(176) = dtau**3*(b1**3 + b2**3)
f(177) = dtau*(b1**5 + b2**5)
v = sum(f*params)
end function c2h4_dipole_b2u_n2_d6_ADF


!###############################################################################


function c2h4_dipole_b2u_n3_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(680)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(680)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(-r1**2*r2 + r1*r2**2 - r3**2*r4 + r3*r4**2)
f(2) = r0*(r1**3*r2 - r1*r2**3 + r3**3*r4 - r3*r4**3)
f(3) = r0*(r1**4*r2 - r1*r2**4 + r3**4*r4 - r3*r4**4)
f(4) = r0*(r1**3*r2**2 - r1**2*r2**3 + r3**3*r4**2 - r3**2*r4**3)
f(5) = r0*(-r1**3*r2 + r1*r2**3 - r3**3*r4 + r3*r4**3)
f(6) = r0*(-r1**3*r2**2 + r1**2*r2**3 - r3**3*r4**2 + r3**2*r4**3)
f(7) = r0*(-r1**4*r2 + r1*r2**4 - r3**4*r4 + r3*r4**4)
f(8) = r0**2*(r1**2*r2 - r1*r2**2 + r3**2*r4 - r3*r4**2)
f(9) = r0**2*(r1**3*r2 - r1*r2**3 + r3**3*r4 - r3*r4**3)
f(10) = r0**2*(-r1**3*r2 + r1*r2**3 - r3**3*r4 + r3*r4**3)
f(11) = r0**3*(-r1**2*r2 + r1*r2**2 - r3**2*r4 + r3*r4**2)
f(12) = r0**3*(r1**2*r2 - r1*r2**2 + r3**2*r4 - r3*r4**2)
f(13) = r0*(-r1*r3 + r2*r4)
f(14) = r0*(r1**2*r3 + r1*r3**2 - r2**2*r4 - r2*r4**2)
f(15) = r0*(-r1**3*r3 - r1*r3**3 + r2**3*r4 + r2*r4**3)
f(16) = r0*(-r1**4*r3 - r1*r3**4 + r2**4*r4 + r2*r4**4)
f(17) = r0*(-r1**2*r3**2 + r2**2*r4**2)
f(18) = r0*(r1**3*r3**2 + r1**2*r3**3 - r2**3*r4**2 - r2**2*r4**3)
f(19) = r0*(-r1**3*r3**2 - r1**2*r3**3 + r2**3*r4**2 + r2**2*r4**3)
f(20) = r0**2*(r1*r3 - r2*r4)
f(21) = r0**2*(r1**2*r3 + r1*r3**2 - r2**2*r4 - r2*r4**2)
f(22) = r0**2*(-r1**3*r3 - r1*r3**3 + r2**3*r4 + r2*r4**3)
f(23) = r0**2*(-r1**2*r3**2 + r2**2*r4**2)
f(24) = r0**3*(-r1*r3 + r2*r4)
f(25) = r0**3*(r1**2*r3 + r1*r3**2 - r2**2*r4 - r2*r4**2)
f(26) = r0**4*(-r1*r3 + r2*r4)
f(27) = r0*(-r1**2*r4 + r1*r4**2 + r2**2*r3 - r2*r3**2)
f(28) = r0*(r1**3*r4 - r1*r4**3 - r2**3*r3 + r2*r3**3)
f(29) = r0*(r1**4*r4 - r1*r4**4 - r2**4*r3 + r2*r3**4)
f(30) = r0*(r1**3*r4**2 - r1**2*r4**3 - r2**3*r3**2 + r2**2*r3**3)
f(31) = r0*(-r1**3*r4**2 + r1**2*r4**3 + r2**3*r3**2 - r2**2*r3**3)
f(32) = r0**2*(r1**2*r4 - r1*r4**2 - r2**2*r3 + r2*r3**2)
f(33) = r0**2*(r1**3*r4 - r1*r4**3 - r2**3*r3 + r2*r3**3)
f(34) = r0**3*(r1**2*r4 - r1*r4**2 - r2**2*r3 + r2*r3**2)
f(35) = r0**3*(-r1**2*r4 + r1*r4**2 + r2**2*r3 - r2*r3**2)
f(36) = r0*(a1*r1 - a2*r2 + a3*r3 - a4*r4)
f(37) = r0*(-a1**2*r1 + a2**2*r2 - a3**2*r3 + a4**2*r4)
f(38) = r0*(a1**3*r1 - a2**3*r2 + a3**3*r3 - a4**3*r4)
f(39) = r0*(a1**4*r1 - a2**4*r2 + a3**4*r3 - a4**4*r4)
f(40) = r0*(a1*r1**2 - a2*r2**2 + a3*r3**2 - a4*r4**2)
f(41) = r0*(-a1**2*r1**2 + a2**2*r2**2 - a3**2*r3**2 + a4**2*r4**2)
f(42) = r0*(-a1**3*r1**2 + a2**3*r2**2 - a3**3*r3**2 + a4**3*r4**2)
f(43) = r0*(a1*r1**3 - a2*r2**3 + a3*r3**3 - a4*r4**3)
f(44) = r0*(-a1**2*r1**3 + a2**2*r2**3 - a3**2*r3**3 + a4**2*r4**3)
f(45) = r0*(-a1*r1**4 + a2*r2**4 - a3*r3**4 + a4*r4**4)
f(46) = r0**2*(-a1*r1 + a2*r2 - a3*r3 + a4*r4)
f(47) = r0**2*(a1**2*r1 - a2**2*r2 + a3**2*r3 - a4**2*r4)
f(48) = r0**2*(-a1**3*r1 + a2**3*r2 - a3**3*r3 + a4**3*r4)
f(49) = r0**2*(-a1*r1**2 + a2*r2**2 - a3*r3**2 + a4*r4**2)
f(50) = r0**2*(a1**2*r1**2 - a2**2*r2**2 + a3**2*r3**2 - a4**2*r4**2)
f(51) = r0**2*(a1*r1**3 - a2*r2**3 + a3*r3**3 - a4*r4**3)
f(52) = r0**3*(a1*r1 - a2*r2 + a3*r3 - a4*r4)
f(53) = r0**3*(-a1**2*r1 + a2**2*r2 - a3**2*r3 + a4**2*r4)
f(54) = r0**3*(a1*r1**2 - a2*r2**2 + a3*r3**2 - a4*r4**2)
f(55) = r0**4*(a1*r1 - a2*r2 + a3*r3 - a4*r4)
f(56) = r0*(-a1*r2 + a2*r1 - a3*r4 + a4*r3)
f(57) = r0*(a1**2*r2 - a2**2*r1 + a3**2*r4 - a4**2*r3)
f(58) = r0*(-a1**3*r2 + a2**3*r1 - a3**3*r4 + a4**3*r3)
f(59) = r0*(a1**4*r2 - a2**4*r1 + a3**4*r4 - a4**4*r3)
f(60) = r0*(a1*r2**2 - a2*r1**2 + a3*r4**2 - a4*r3**2)
f(61) = r0*(a1**2*r2**2 - a2**2*r1**2 + a3**2*r4**2 - a4**2*r3**2)
f(62) = r0*(-a1**3*r2**2 + a2**3*r1**2 - a3**3*r4**2 + a4**3*r3**2)
f(63) = r0*(-a1*r2**3 + a2*r1**3 - a3*r4**3 + a4*r3**3)
f(64) = r0*(-a1**2*r2**3 + a2**2*r1**3 - a3**2*r4**3 + a4**2*r3**3)
f(65) = r0*(a1*r2**4 - a2*r1**4 + a3*r4**4 - a4*r3**4)
f(66) = r0**2*(a1*r2 - a2*r1 + a3*r4 - a4*r3)
f(67) = r0**2*(a1**2*r2 - a2**2*r1 + a3**2*r4 - a4**2*r3)
f(68) = r0**2*(-a1**3*r2 + a2**3*r1 - a3**3*r4 + a4**3*r3)
f(69) = r0**2*(a1*r2**2 - a2*r1**2 + a3*r4**2 - a4*r3**2)
f(70) = r0**2*(a1**2*r2**2 - a2**2*r1**2 + a3**2*r4**2 - a4**2*r3**2)
f(71) = r0**2*(a1*r2**3 - a2*r1**3 + a3*r4**3 - a4*r3**3)
f(72) = r0**3*(a1*r2 - a2*r1 + a3*r4 - a4*r3)
f(73) = r0**3*(-a1**2*r2 + a2**2*r1 - a3**2*r4 + a4**2*r3)
f(74) = r0**3*(a1*r2**2 - a2*r1**2 + a3*r4**2 - a4*r3**2)
f(75) = r0**4*(a1*r2 - a2*r1 + a3*r4 - a4*r3)
f(76) = r0*(-a1*r3 + a2*r4 - a3*r1 + a4*r2)
f(77) = r0*(-a1**2*r3 + a2**2*r4 - a3**2*r1 + a4**2*r2)
f(78) = r0*(-a1**3*r3 + a2**3*r4 - a3**3*r1 + a4**3*r2)
f(79) = r0*(-a1**4*r3 + a2**4*r4 - a3**4*r1 + a4**4*r2)
f(80) = r0*(-a1*r3**2 + a2*r4**2 - a3*r1**2 + a4*r2**2)
f(81) = r0*(a1**2*r3**2 - a2**2*r4**2 + a3**2*r1**2 - a4**2*r2**2)
f(82) = r0*(a1**3*r3**2 - a2**3*r4**2 + a3**3*r1**2 - a4**3*r2**2)
f(83) = r0*(a1*r3**3 - a2*r4**3 + a3*r1**3 - a4*r2**3)
f(84) = r0*(a1**2*r3**3 - a2**2*r4**3 + a3**2*r1**3 - a4**2*r2**3)
f(85) = r0*(a1*r3**4 - a2*r4**4 + a3*r1**4 - a4*r2**4)
f(86) = r0**2*(a1*r3 - a2*r4 + a3*r1 - a4*r2)
f(87) = r0**2*(a1**2*r3 - a2**2*r4 + a3**2*r1 - a4**2*r2)
f(88) = r0**2*(a1**3*r3 - a2**3*r4 + a3**3*r1 - a4**3*r2)
f(89) = r0**2*(-a1*r3**2 + a2*r4**2 - a3*r1**2 + a4*r2**2)
f(90) = r0**2*(-a1**2*r3**2 + a2**2*r4**2 - a3**2*r1**2 + a4**2*r2**2)
f(91) = r0**2*(-a1*r3**3 + a2*r4**3 - a3*r1**3 + a4*r2**3)
f(92) = r0**3*(a1*r3 - a2*r4 + a3*r1 - a4*r2)
f(93) = r0**3*(-a1**2*r3 + a2**2*r4 - a3**2*r1 + a4**2*r2)
f(94) = r0**3*(-a1*r3**2 + a2*r4**2 - a3*r1**2 + a4*r2**2)
f(95) = r0**4*(a1*r3 - a2*r4 + a3*r1 - a4*r2)
f(96) = r0*(-a1*r4 + a2*r3 - a3*r2 + a4*r1)
f(97) = r0*(a1**2*r4 - a2**2*r3 + a3**2*r2 - a4**2*r1)
f(98) = r0*(a1**3*r4 - a2**3*r3 + a3**3*r2 - a4**3*r1)
f(99) = r0*(a1**4*r4 - a2**4*r3 + a3**4*r2 - a4**4*r1)
f(100) = r0*(-a1*r4**2 + a2*r3**2 - a3*r2**2 + a4*r1**2)
f(101) = r0*(-a1**2*r4**2 + a2**2*r3**2 - a3**2*r2**2 + a4**2*r1**2)
f(102) = r0*(-a1**3*r4**2 + a2**3*r3**2 - a3**3*r2**2 + a4**3*r1**2)
f(103) = r0*(-a1*r4**3 + a2*r3**3 - a3*r2**3 + a4*r1**3)
f(104) = r0*(a1**2*r4**3 - a2**2*r3**3 + a3**2*r2**3 - a4**2*r1**3)
f(105) = r0*(a1*r4**4 - a2*r3**4 + a3*r2**4 - a4*r1**4)
f(106) = r0**2*(-a1*r4 + a2*r3 - a3*r2 + a4*r1)
f(107) = r0**2*(a1**2*r4 - a2**2*r3 + a3**2*r2 - a4**2*r1)
f(108) = r0**2*(-a1**3*r4 + a2**3*r3 - a3**3*r2 + a4**3*r1)
f(109) = r0**2*(-a1*r4**2 + a2*r3**2 - a3*r2**2 + a4*r1**2)
f(110) = r0**2*(-a1**2*r4**2 + a2**2*r3**2 - a3**2*r2**2 + a4**2*r1**2)
f(111) = r0**2*(a1*r4**3 - a2*r3**3 + a3*r2**3 - a4*r1**3)
f(112) = r0**3*(-a1*r4 + a2*r3 - a3*r2 + a4*r1)
f(113) = r0**3*(a1**2*r4 - a2**2*r3 + a3**2*r2 - a4**2*r1)
f(114) = r0**3*(-a1*r4**2 + a2*r3**2 - a3*r2**2 + a4*r1**2)
f(115) = r0**4*(-a1*r4 + a2*r3 - a3*r2 + a4*r1)
f(116) = r0*(-b1**2*r1 + b1**2*r2 - b2**2*r3 + b2**2*r4)
f(117) = r0*(-b1**4*r1 + b1**4*r2 - b2**4*r3 + b2**4*r4)
f(118) = r0*(b1**2*r1**2 - b1**2*r2**2 + b2**2*r3**2 - b2**2*r4**2)
f(119) = r0*(b1**2*r1**3 - b1**2*r2**3 + b2**2*r3**3 - b2**2*r4**3)
f(120) = r0**2*(b1**2*r1 - b1**2*r2 + b2**2*r3 - b2**2*r4)
f(121) = r0**2*(b1**2*r1**2 - b1**2*r2**2 + b2**2*r3**2 - b2**2*r4**2)
f(122) = r0**3*(b1**2*r1 - b1**2*r2 + b2**2*r3 - b2**2*r4)
f(123) = r0*(-b1**2*r3 + b1**2*r4 - b2**2*r1 + b2**2*r2)
f(124) = r0*(b1**4*r3 - b1**4*r4 + b2**4*r1 - b2**4*r2)
f(125) = r0*(b1**2*r3**2 - b1**2*r4**2 + b2**2*r1**2 - b2**2*r2**2)
f(126) = r0*(-b1**2*r3**3 + b1**2*r4**3 - b2**2*r1**3 + b2**2*r2**3)
f(127) = r0**2*(b1**2*r3 - b1**2*r4 + b2**2*r1 - b2**2*r2)
f(128) = r0**2*(b1**2*r3**2 - b1**2*r4**2 + b2**2*r1**2 - b2**2*r2**2)
f(129) = r0**3*(-b1**2*r3 + b1**2*r4 - b2**2*r1 + b2**2*r2)
f(130) = dtau**2*r0*(-r1 + r2 - r3 + r4)
f(131) = dtau**4*r0*(r1 - r2 + r3 - r4)
f(132) = dtau**2*r0*(-r1**2 + r2**2 - r3**2 + r4**2)
f(133) = dtau**2*r0*(-r1**3 + r2**3 - r3**3 + r4**3)
f(134) = dtau**2*r0**2*(-r1 + r2 - r3 + r4)
f(135) = dtau**2*r0**2*(-r1**2 + r2**2 - r3**2 + r4**2)
f(136) = dtau**2*r0**3*(r1 - r2 + r3 - r4)
f(137) = r0*(-a1**2*a2 + a1*a2**2 - a3**2*a4 + a3*a4**2)
f(138) = r0*(a1**3*a2 - a1*a2**3 + a3**3*a4 - a3*a4**3)
f(139) = r0*(-a1**4*a2 + a1*a2**4 - a3**4*a4 + a3*a4**4)
f(140) = r0*(-a1**3*a2**2 + a1**2*a2**3 - a3**3*a4**2 + a3**2*a4**3)
f(141) = r0**2*(a1**2*a2 - a1*a2**2 + a3**2*a4 - a3*a4**2)
f(142) = r0**2*(a1**3*a2 - a1*a2**3 + a3**3*a4 - a3*a4**3)
f(143) = r0**3*(a1**2*a2 - a1*a2**2 + a3**2*a4 - a3*a4**2)
f(144) = r0*(a1*a3 - a2*a4)
f(145) = r0*(a1**2*a3 + a1*a3**2 - a2**2*a4 - a2*a4**2)
f(146) = r0*(a1**3*a3 + a1*a3**3 - a2**3*a4 - a2*a4**3)
f(147) = r0*(a1**4*a3 + a1*a3**4 - a2**4*a4 - a2*a4**4)
f(148) = r0*(-a1**2*a3**2 + a2**2*a4**2)
f(149) = r0*(-a1**3*a3**2 - a1**2*a3**3 + a2**3*a4**2 + a2**2*a4**3)
f(150) = r0**2*(a1*a3 - a2*a4)
f(151) = r0**2*(a1**2*a3 + a1*a3**2 - a2**2*a4 - a2*a4**2)
f(152) = r0**2*(a1**3*a3 + a1*a3**3 - a2**3*a4 - a2*a4**3)
f(153) = r0**2*(a1**2*a3**2 - a2**2*a4**2)
f(154) = r0**3*(-a1*a3 + a2*a4)
f(155) = r0**3*(-a1**2*a3 - a1*a3**2 + a2**2*a4 + a2*a4**2)
f(156) = r0**4*(-a1*a3 + a2*a4)
f(157) = r0*(a1**2*a4 - a1*a4**2 - a2**2*a3 + a2*a3**2)
f(158) = r0*(-a1**3*a4 + a1*a4**3 + a2**3*a3 - a2*a3**3)
f(159) = r0*(-a1**4*a4 + a1*a4**4 + a2**4*a3 - a2*a3**4)
f(160) = r0*(a1**3*a4**2 - a1**2*a4**3 - a2**3*a3**2 + a2**2*a3**3)
f(161) = r0*(-a1**3*a4**2 + a1**2*a4**3 + a2**3*a3**2 - a2**2*a3**3)
f(162) = r0**2*(-a1**2*a4 + a1*a4**2 + a2**2*a3 - a2*a3**2)
f(163) = r0**2*(-a1**3*a4 + a1*a4**3 + a2**3*a3 - a2*a3**3)
f(164) = r0**3*(-a1**2*a4 + a1*a4**2 + a2**2*a3 - a2*a3**2)
f(165) = r0*(a1*b1**2 - a2*b1**2 + a3*b2**2 - a4*b2**2)
f(166) = r0*(a1*b1**4 - a2*b1**4 + a3*b2**4 - a4*b2**4)
f(167) = r0*(-a1**2*b1**2 + a2**2*b1**2 - a3**2*b2**2 + a4**2*b2**2)
f(168) = r0*(a1**3*b1**2 - a2**3*b1**2 + a3**3*b2**2 - a4**3*b2**2)
f(169) = r0**2*(-a1*b1**2 + a2*b1**2 - a3*b2**2 + a4*b2**2)
f(170) = r0**2*(a1**2*b1**2 - a2**2*b1**2 + a3**2*b2**2 - a4**2*b2**2)
f(171) = r0**3*(a1*b1**2 - a2*b1**2 + a3*b2**2 - a4*b2**2)
f(172) = r0*(-a1*b2**2 + a2*b2**2 - a3*b1**2 + a4*b1**2)
f(173) = r0*(a1*b2**4 - a2*b2**4 + a3*b1**4 - a4*b1**4)
f(174) = r0*(-a1**2*b2**2 + a2**2*b2**2 - a3**2*b1**2 + a4**2*b1**2)
f(175) = r0*(a1**3*b2**2 - a2**3*b2**2 + a3**3*b1**2 - a4**3*b1**2)
f(176) = r0**2*(a1*b2**2 - a2*b2**2 + a3*b1**2 - a4*b1**2)
f(177) = r0**2*(a1**2*b2**2 - a2**2*b2**2 + a3**2*b1**2 - a4**2*b1**2)
f(178) = r0**3*(a1*b2**2 - a2*b2**2 + a3*b1**2 - a4*b1**2)
f(179) = dtau**2*r0*(a1 - a2 + a3 - a4)
f(180) = dtau**4*r0*(-a1 + a2 - a3 + a4)
f(181) = dtau**2*r0*(a1**2 - a2**2 + a3**2 - a4**2)
f(182) = dtau**2*r0*(-a1**3 + a2**3 - a3**3 + a4**3)
f(183) = dtau**2*r0**2*(a1 - a2 + a3 - a4)
f(184) = dtau**2*r0**2*(-a1**2 + a2**2 - a3**2 + a4**2)
f(185) = dtau**2*r0**3*(-a1 + a2 - a3 + a4)
f(186) = dtau*r0*(b1 + b2)
f(187) = dtau**3*r0*(b1 + b2)
f(188) = dtau*r0*(b1**3 + b2**3)
f(189) = dtau*r0**2*(b1 + b2)
f(190) = dtau**3*r0**2*(b1 + b2)
f(191) = dtau*r0**2*(b1**3 + b2**3)
f(192) = dtau*r0**3*(b1 + b2)
f(193) = dtau*r0**4*(b1 + b2)
f(194) = r1*r2*r3 - r1*r2*r4 + r1*r3*r4 - r2*r3*r4
f(195) = -r1**2*r3*r4 - r1*r2*r3**2 + r1*r2*r4**2 + r2**2*r3*r4
f(196) = -r1**3*r3*r4 - r1*r2*r3**3 + r1*r2*r4**3 + r2**3*r3*r4
f(197) = -r1**4*r3*r4 - r1*r2*r3**4 + r1*r2*r4**4 + r2**4*r3*r4
f(198) = r1**2*r2*r4 - r1*r2**2*r3 - r1*r3*r4**2 + r2*r3**2*r4
f(199) = -r1**2*r2*r4**2 + r1**2*r3*r4**2 + r1*r2**2*r3**2 - r2**2*r3**2 &
      *r4
f(200) = r1**3*r3*r4**2 - r1**2*r2*r4**3 + r1*r2**2*r3**3 - r2**3*r3**2* &
      r4
f(201) = -r1**3*r2*r4 + r1*r2**3*r3 + r1*r3*r4**3 - r2*r3**3*r4
f(202) = -r1**3*r2*r4**2 + r1**2*r3*r4**3 + r1*r2**3*r3**2 - r2**2*r3**3 &
      *r4
f(203) = -r1**4*r2*r4 + r1*r2**4*r3 + r1*r3*r4**4 - r2*r3**4*r4
f(204) = -r1**2*r2*r3 + r1*r2**2*r4 - r1*r3**2*r4 + r2*r3*r4**2
f(205) = r1**2*r2*r3**2 + r1**2*r3**2*r4 - r1*r2**2*r4**2 - r2**2*r3*r4 &
      **2
f(206) = -r1**3*r3**2*r4 - r1**2*r2*r3**3 + r1*r2**2*r4**3 + r2**3*r3*r4 &
      **2
f(207) = r1**2*r2**2*r3 - r1**2*r2**2*r4 + r1*r3**2*r4**2 - r2*r3**2*r4 &
      **2
f(208) = r1**2*r2**2*r3**2 - r1**2*r2**2*r4**2 + r1**2*r3**2*r4**2 - r2 &
      **2*r3**2*r4**2
f(209) = -r1**3*r2**2*r4 + r1**2*r2**3*r3 + r1*r3**2*r4**3 - r2*r3**3*r4 &
      **2
f(210) = r1**3*r2*r3 - r1*r2**3*r4 + r1*r3**3*r4 - r2*r3*r4**3
f(211) = -r1**3*r2*r3**2 - r1**2*r3**3*r4 + r1*r2**3*r4**2 + r2**2*r3*r4 &
      **3
f(212) = -r1**3*r2**2*r3 + r1**2*r2**3*r4 - r1*r3**3*r4**2 + r2*r3**2*r4 &
      **3
f(213) = -r1**4*r2*r3 + r1*r2**4*r4 - r1*r3**4*r4 + r2*r3*r4**4
f(214) = -a1*r1*r2 + a2*r1*r2 - a3*r3*r4 + a4*r3*r4
f(215) = -a1**2*r1*r2 + a2**2*r1*r2 - a3**2*r3*r4 + a4**2*r3*r4
f(216) = a1**3*r1*r2 - a2**3*r1*r2 + a3**3*r3*r4 - a4**3*r3*r4
f(217) = a1**4*r1*r2 - a2**4*r1*r2 + a3**4*r3*r4 - a4**4*r3*r4
f(218) = -a1*r1*r2**2 + a2*r1**2*r2 - a3*r3*r4**2 + a4*r3**2*r4
f(219) = -a1**2*r1*r2**2 + a2**2*r1**2*r2 - a3**2*r3*r4**2 + a4**2*r3**2 &
      *r4
f(220) = a1**3*r1*r2**2 - a2**3*r1**2*r2 + a3**3*r3*r4**2 - a4**3*r3**2* &
      r4
f(221) = -a1*r1*r2**3 + a2*r1**3*r2 - a3*r3*r4**3 + a4*r3**3*r4
f(222) = a1**2*r1*r2**3 - a2**2*r1**3*r2 + a3**2*r3*r4**3 - a4**2*r3**3* &
      r4
f(223) = a1*r1*r2**4 - a2*r1**4*r2 + a3*r3*r4**4 - a4*r3**4*r4
f(224) = -a1*r1**2*r2 + a2*r1*r2**2 - a3*r3**2*r4 + a4*r3*r4**2
f(225) = -a1**2*r1**2*r2 + a2**2*r1*r2**2 - a3**2*r3**2*r4 + a4**2*r3*r4 &
      **2
f(226) = a1**3*r1**2*r2 - a2**3*r1*r2**2 + a3**3*r3**2*r4 - a4**3*r3*r4 &
      **2
f(227) = -a1*r1**2*r2**2 + a2*r1**2*r2**2 - a3*r3**2*r4**2 + a4*r3**2*r4 &
      **2
f(228) = -a1**2*r1**2*r2**2 + a2**2*r1**2*r2**2 - a3**2*r3**2*r4**2 + a4 &
      **2*r3**2*r4**2
f(229) = -a1*r1**2*r2**3 + a2*r1**3*r2**2 - a3*r3**2*r4**3 + a4*r3**3*r4 &
      **2
f(230) = -a1*r1**3*r2 + a2*r1*r2**3 - a3*r3**3*r4 + a4*r3*r4**3
f(231) = a1**2*r1**3*r2 - a2**2*r1*r2**3 + a3**2*r3**3*r4 - a4**2*r3*r4 &
      **3
f(232) = a1*r1**3*r2**2 - a2*r1**2*r2**3 + a3*r3**3*r4**2 - a4*r3**2*r4 &
      **3
f(233) = -a1*r1**4*r2 + a2*r1*r2**4 - a3*r3**4*r4 + a4*r3*r4**4
f(234) = -a1*r3*r4 + a2*r3*r4 - a3*r1*r2 + a4*r1*r2
f(235) = a1**2*r3*r4 - a2**2*r3*r4 + a3**2*r1*r2 - a4**2*r1*r2
f(236) = -a1**3*r3*r4 + a2**3*r3*r4 - a3**3*r1*r2 + a4**3*r1*r2
f(237) = -a1**4*r3*r4 + a2**4*r3*r4 - a3**4*r1*r2 + a4**4*r1*r2
f(238) = a1*r3*r4**2 - a2*r3**2*r4 + a3*r1*r2**2 - a4*r1**2*r2
f(239) = a1**2*r3*r4**2 - a2**2*r3**2*r4 + a3**2*r1*r2**2 - a4**2*r1**2* &
      r2
f(240) = -a1**3*r3*r4**2 + a2**3*r3**2*r4 - a3**3*r1*r2**2 + a4**3*r1**2 &
      *r2
f(241) = -a1*r3*r4**3 + a2*r3**3*r4 - a3*r1*r2**3 + a4*r1**3*r2
f(242) = a1**2*r3*r4**3 - a2**2*r3**3*r4 + a3**2*r1*r2**3 - a4**2*r1**3* &
      r2
f(243) = -a1*r3*r4**4 + a2*r3**4*r4 - a3*r1*r2**4 + a4*r1**4*r2
f(244) = a1*r3**2*r4 - a2*r3*r4**2 + a3*r1**2*r2 - a4*r1*r2**2
f(245) = a1**2*r3**2*r4 - a2**2*r3*r4**2 + a3**2*r1**2*r2 - a4**2*r1*r2 &
      **2
f(246) = -a1**3*r3**2*r4 + a2**3*r3*r4**2 - a3**3*r1**2*r2 + a4**3*r1*r2 &
      **2
f(247) = -a1*r3**2*r4**2 + a2*r3**2*r4**2 - a3*r1**2*r2**2 + a4*r1**2*r2 &
      **2
f(248) = -a1**2*r3**2*r4**2 + a2**2*r3**2*r4**2 - a3**2*r1**2*r2**2 + a4 &
      **2*r1**2*r2**2
f(249) = -a1*r3**2*r4**3 + a2*r3**3*r4**2 - a3*r1**2*r2**3 + a4*r1**3*r2 &
      **2
f(250) = -a1*r3**3*r4 + a2*r3*r4**3 - a3*r1**3*r2 + a4*r1*r2**3
f(251) = a1**2*r3**3*r4 - a2**2*r3*r4**3 + a3**2*r1**3*r2 - a4**2*r1*r2 &
      **3
f(252) = a1*r3**3*r4**2 - a2*r3**2*r4**3 + a3*r1**3*r2**2 - a4*r1**2*r2 &
      **3
f(253) = a1*r3**4*r4 - a2*r3*r4**4 + a3*r1**4*r2 - a4*r1*r2**4
f(254) = -b1**2*r1**2*r2 + b1**2*r1*r2**2 - b2**2*r3**2*r4 + b2**2*r3*r4 &
      **2
f(255) = b1**2*r1**3*r2 - b1**2*r1*r2**3 + b2**2*r3**3*r4 - b2**2*r3*r4 &
      **3
f(256) = b1**2*r1**2*r2 - b1**2*r1*r2**2 + b2**2*r3**2*r4 - b2**2*r3*r4 &
      **2
f(257) = -b1**2*r1**3*r2 + b1**2*r1*r2**3 - b2**2*r3**3*r4 + b2**2*r3*r4 &
      **3
f(258) = b1**2*r3**2*r4 - b1**2*r3*r4**2 + b2**2*r1**2*r2 - b2**2*r1*r2 &
      **2
f(259) = -b1**2*r3**3*r4 + b1**2*r3*r4**3 - b2**2*r1**3*r2 + b2**2*r1*r2 &
      **3
f(260) = -b1**2*r3**2*r4 + b1**2*r3*r4**2 - b2**2*r1**2*r2 + b2**2*r1*r2 &
      **2
f(261) = b1**2*r3**3*r4 - b1**2*r3*r4**3 + b2**2*r1**3*r2 - b2**2*r1*r2 &
      **3
f(262) = dtau**2*(-r1**2*r2 + r1*r2**2 - r3**2*r4 + r3*r4**2)
f(263) = dtau**2*(r1**3*r2 - r1*r2**3 + r3**3*r4 - r3*r4**3)
f(264) = dtau**2*(-r1**3*r2 + r1*r2**3 - r3**3*r4 + r3*r4**3)
f(265) = a1*r1*r3 - a2*r2*r4 + a3*r1*r3 - a4*r2*r4
f(266) = a1**2*r1*r3 - a2**2*r2*r4 + a3**2*r1*r3 - a4**2*r2*r4
f(267) = a1**3*r1*r3 - a2**3*r2*r4 + a3**3*r1*r3 - a4**3*r2*r4
f(268) = a1**4*r1*r3 - a2**4*r2*r4 + a3**4*r1*r3 - a4**4*r2*r4
f(269) = a1*r1*r3**2 - a2*r2*r4**2 + a3*r1**2*r3 - a4*r2**2*r4
f(270) = -a1**2*r1*r3**2 + a2**2*r2*r4**2 - a3**2*r1**2*r3 + a4**2*r2**2 &
      *r4
f(271) = -a1**3*r1*r3**2 + a2**3*r2*r4**2 - a3**3*r1**2*r3 + a4**3*r2**2 &
      *r4
f(272) = a1*r1*r3**3 - a2*r2*r4**3 + a3*r1**3*r3 - a4*r2**3*r4
f(273) = a1**2*r1*r3**3 - a2**2*r2*r4**3 + a3**2*r1**3*r3 - a4**2*r2**3* &
      r4
f(274) = -a1*r1*r3**4 + a2*r2*r4**4 - a3*r1**4*r3 + a4*r2**4*r4
f(275) = a1*r1**2*r3 - a2*r2**2*r4 + a3*r1*r3**2 - a4*r2*r4**2
f(276) = -a1**2*r1**2*r3 + a2**2*r2**2*r4 - a3**2*r1*r3**2 + a4**2*r2*r4 &
      **2
f(277) = -a1**3*r1**2*r3 + a2**3*r2**2*r4 - a3**3*r1*r3**2 + a4**3*r2*r4 &
      **2
f(278) = a1*r1**2*r3**2 - a2*r2**2*r4**2 + a3*r1**2*r3**2 - a4*r2**2*r4 &
      **2
f(279) = -a1**2*r1**2*r3**2 + a2**2*r2**2*r4**2 - a3**2*r1**2*r3**2 + a4 &
      **2*r2**2*r4**2
f(280) = -a1*r1**2*r3**3 + a2*r2**2*r4**3 - a3*r1**3*r3**2 + a4*r2**3*r4 &
      **2
f(281) = -a1*r1**3*r3 + a2*r2**3*r4 - a3*r1*r3**3 + a4*r2*r4**3
f(282) = -a1**2*r1**3*r3 + a2**2*r2**3*r4 - a3**2*r1*r3**3 + a4**2*r2*r4 &
      **3
f(283) = -a1*r1**3*r3**2 + a2*r2**3*r4**2 - a3*r1**2*r3**3 + a4*r2**2*r4 &
      **3
f(284) = -a1*r1**4*r3 + a2*r2**4*r4 - a3*r1*r3**4 + a4*r2*r4**4
f(285) = a1*r2*r4 - a2*r1*r3 + a3*r2*r4 - a4*r1*r3
f(286) = -a1**2*r2*r4 + a2**2*r1*r3 - a3**2*r2*r4 + a4**2*r1*r3
f(287) = -a1**3*r2*r4 + a2**3*r1*r3 - a3**3*r2*r4 + a4**3*r1*r3
f(288) = -a1**4*r2*r4 + a2**4*r1*r3 - a3**4*r2*r4 + a4**4*r1*r3
f(289) = -a1*r2*r4**2 + a2*r1*r3**2 - a3*r2**2*r4 + a4*r1**2*r3
f(290) = a1**2*r2*r4**2 - a2**2*r1*r3**2 + a3**2*r2**2*r4 - a4**2*r1**2* &
      r3
f(291) = a1**3*r2*r4**2 - a2**3*r1*r3**2 + a3**3*r2**2*r4 - a4**3*r1**2* &
      r3
f(292) = -a1*r2*r4**3 + a2*r1*r3**3 - a3*r2**3*r4 + a4*r1**3*r3
f(293) = -a1**2*r2*r4**3 + a2**2*r1*r3**3 - a3**2*r2**3*r4 + a4**2*r1**3 &
      *r3
f(294) = -a1*r2*r4**4 + a2*r1*r3**4 - a3*r2**4*r4 + a4*r1**4*r3
f(295) = -a1*r2**2*r4 + a2*r1**2*r3 - a3*r2*r4**2 + a4*r1*r3**2
f(296) = a1**2*r2**2*r4 - a2**2*r1**2*r3 + a3**2*r2*r4**2 - a4**2*r1*r3 &
      **2
f(297) = a1**3*r2**2*r4 - a2**3*r1**2*r3 + a3**3*r2*r4**2 - a4**3*r1*r3 &
      **2
f(298) = a1*r2**2*r4**2 - a2*r1**2*r3**2 + a3*r2**2*r4**2 - a4*r1**2*r3 &
      **2
f(299) = -a1**2*r2**2*r4**2 + a2**2*r1**2*r3**2 - a3**2*r2**2*r4**2 + a4 &
      **2*r1**2*r3**2
f(300) = -a1*r2**2*r4**3 + a2*r1**2*r3**3 - a3*r2**3*r4**2 + a4*r1**3*r3 &
      **2
f(301) = a1*r2**3*r4 - a2*r1**3*r3 + a3*r2*r4**3 - a4*r1*r3**3
f(302) = a1**2*r2**3*r4 - a2**2*r1**3*r3 + a3**2*r2*r4**3 - a4**2*r1*r3 &
      **3
f(303) = -a1*r2**3*r4**2 + a2*r1**3*r3**2 - a3*r2**2*r4**3 + a4*r1**2*r3 &
      **3
f(304) = -a1*r2**4*r4 + a2*r1**4*r3 - a3*r2*r4**4 + a4*r1*r3**4
f(305) = -b1**2*r1*r3 + b1**2*r2*r4 - b2**2*r1*r3 + b2**2*r2*r4
f(306) = b1**4*r1*r3 - b1**4*r2*r4 + b2**4*r1*r3 - b2**4*r2*r4
f(307) = -b1**2*r1*r3**2 + b1**2*r2*r4**2 - b2**2*r1**2*r3 + b2**2*r2**2 &
      *r4
f(308) = -b1**2*r1*r3**3 + b1**2*r2*r4**3 - b2**2*r1**3*r3 + b2**2*r2**3 &
      *r4
f(309) = -b1**2*r1**2*r3 + b1**2*r2**2*r4 - b2**2*r1*r3**2 + b2**2*r2*r4 &
      **2
f(310) = -b1**2*r1**2*r3**2 + b1**2*r2**2*r4**2 - b2**2*r1**2*r3**2 + b2 &
      **2*r2**2*r4**2
f(311) = b1**2*r1**3*r3 - b1**2*r2**3*r4 + b2**2*r1*r3**3 - b2**2*r2*r4 &
      **3
f(312) = dtau**2*(r1*r3 - r2*r4)
f(313) = dtau**4*(-r1*r3 + r2*r4)
f(314) = dtau**2*(r1**2*r3 + r1*r3**2 - r2**2*r4 - r2*r4**2)
f(315) = dtau**2*(r1**3*r3 + r1*r3**3 - r2**3*r4 - r2*r4**3)
f(316) = dtau**2*(-r1**2*r3**2 + r2**2*r4**2)
f(317) = -a1*r1*r4 + a2*r2*r3 - a3*r2*r3 + a4*r1*r4
f(318) = -a1**2*r1*r4 + a2**2*r2*r3 - a3**2*r2*r3 + a4**2*r1*r4
f(319) = -a1**3*r1*r4 + a2**3*r2*r3 - a3**3*r2*r3 + a4**3*r1*r4
f(320) = a1**4*r1*r4 - a2**4*r2*r3 + a3**4*r2*r3 - a4**4*r1*r4
f(321) = -a1*r1*r4**2 + a2*r2*r3**2 - a3*r2**2*r3 + a4*r1**2*r4
f(322) = -a1**2*r1*r4**2 + a2**2*r2*r3**2 - a3**2*r2**2*r3 + a4**2*r1**2 &
      *r4
f(323) = -a1**3*r1*r4**2 + a2**3*r2*r3**2 - a3**3*r2**2*r3 + a4**3*r1**2 &
      *r4
f(324) = -a1*r1*r4**3 + a2*r2*r3**3 - a3*r2**3*r3 + a4*r1**3*r4
f(325) = a1**2*r1*r4**3 - a2**2*r2*r3**3 + a3**2*r2**3*r3 - a4**2*r1**3* &
      r4
f(326) = a1*r1*r4**4 - a2*r2*r3**4 + a3*r2**4*r3 - a4*r1**4*r4
f(327) = -a1*r1**2*r4 + a2*r2**2*r3 - a3*r2*r3**2 + a4*r1*r4**2
f(328) = -a1**2*r1**2*r4 + a2**2*r2**2*r3 - a3**2*r2*r3**2 + a4**2*r1*r4 &
      **2
f(329) = -a1**3*r1**2*r4 + a2**3*r2**2*r3 - a3**3*r2*r3**2 + a4**3*r1*r4 &
      **2
f(330) = a1*r1**2*r4**2 - a2*r2**2*r3**2 + a3*r2**2*r3**2 - a4*r1**2*r4 &
      **2
f(331) = -a1**2*r1**2*r4**2 + a2**2*r2**2*r3**2 - a3**2*r2**2*r3**2 + a4 &
      **2*r1**2*r4**2
f(332) = -a1*r1**2*r4**3 + a2*r2**2*r3**3 - a3*r2**3*r3**2 + a4*r1**3*r4 &
      **2
f(333) = a1*r1**3*r4 - a2*r2**3*r3 + a3*r2*r3**3 - a4*r1*r4**3
f(334) = -a1**2*r1**3*r4 + a2**2*r2**3*r3 - a3**2*r2*r3**3 + a4**2*r1*r4 &
      **3
f(335) = -a1*r1**3*r4**2 + a2*r2**3*r3**2 - a3*r2**2*r3**3 + a4*r1**2*r4 &
      **3
f(336) = -a1*r1**4*r4 + a2*r2**4*r3 - a3*r2*r3**4 + a4*r1*r4**4
f(337) = -a1*r2*r3 + a2*r1*r4 - a3*r1*r4 + a4*r2*r3
f(338) = a1**2*r2*r3 - a2**2*r1*r4 + a3**2*r1*r4 - a4**2*r2*r3
f(339) = a1**3*r2*r3 - a2**3*r1*r4 + a3**3*r1*r4 - a4**3*r2*r3
f(340) = -a1**4*r2*r3 + a2**4*r1*r4 - a3**4*r1*r4 + a4**4*r2*r3
f(341) = a1*r2*r3**2 - a2*r1*r4**2 + a3*r1**2*r4 - a4*r2**2*r3
f(342) = a1**2*r2*r3**2 - a2**2*r1*r4**2 + a3**2*r1**2*r4 - a4**2*r2**2* &
      r3
f(343) = a1**3*r2*r3**2 - a2**3*r1*r4**2 + a3**3*r1**2*r4 - a4**3*r2**2* &
      r3
f(344) = a1*r2*r3**3 - a2*r1*r4**3 + a3*r1**3*r4 - a4*r2**3*r3
f(345) = -a1**2*r2*r3**3 + a2**2*r1*r4**3 - a3**2*r1**3*r4 + a4**2*r2**3 &
      *r3
f(346) = a1*r2*r3**4 - a2*r1*r4**4 + a3*r1**4*r4 - a4*r2**4*r3
f(347) = a1*r2**2*r3 - a2*r1**2*r4 + a3*r1*r4**2 - a4*r2*r3**2
f(348) = a1**2*r2**2*r3 - a2**2*r1**2*r4 + a3**2*r1*r4**2 - a4**2*r2*r3 &
      **2
f(349) = a1**3*r2**2*r3 - a2**3*r1**2*r4 + a3**3*r1*r4**2 - a4**3*r2*r3 &
      **2
f(350) = a1*r2**2*r3**2 - a2*r1**2*r4**2 + a3*r1**2*r4**2 - a4*r2**2*r3 &
      **2
f(351) = -a1**2*r2**2*r3**2 + a2**2*r1**2*r4**2 - a3**2*r1**2*r4**2 + a4 &
      **2*r2**2*r3**2
f(352) = -a1*r2**2*r3**3 + a2*r1**2*r4**3 - a3*r1**3*r4**2 + a4*r2**3*r3 &
      **2
f(353) = -a1*r2**3*r3 + a2*r1**3*r4 - a3*r1*r4**3 + a4*r2*r3**3
f(354) = a1**2*r2**3*r3 - a2**2*r1**3*r4 + a3**2*r1*r4**3 - a4**2*r2*r3 &
      **3
f(355) = -a1*r2**3*r3**2 + a2*r1**3*r4**2 - a3*r1**2*r4**3 + a4*r2**2*r3 &
      **3
f(356) = -a1*r2**4*r3 + a2*r1**4*r4 - a3*r1*r4**4 + a4*r2*r3**4
f(357) = b1**2*r1*r4 - b1**2*r2*r3 - b2**2*r1*r4 + b2**2*r2*r3
f(358) = -b1**4*r1*r4 + b1**4*r2*r3 + b2**4*r1*r4 - b2**4*r2*r3
f(359) = b1**2*r1*r4**2 - b1**2*r2*r3**2 - b2**2*r1**2*r4 + b2**2*r2**2* &
      r3
f(360) = b1**2*r1*r4**3 - b1**2*r2*r3**3 - b2**2*r1**3*r4 + b2**2*r2**3* &
      r3
f(361) = b1**2*r1**2*r4 - b1**2*r2**2*r3 - b2**2*r1*r4**2 + b2**2*r2*r3 &
      **2
f(362) = b1**2*r1**2*r4**2 - b1**2*r2**2*r3**2 - b2**2*r1**2*r4**2 + b2 &
      **2*r2**2*r3**2
f(363) = -b1**2*r1**3*r4 + b1**2*r2**3*r3 + b2**2*r1*r4**3 - b2**2*r2*r3 &
      **3
f(364) = dtau**2*(-r1**2*r4 + r1*r4**2 + r2**2*r3 - r2*r3**2)
f(365) = dtau**2*(r1**3*r4 - r1*r4**3 - r2**3*r3 + r2*r3**3)
f(366) = -a1*a2*r1 + a1*a2*r2 - a3*a4*r3 + a3*a4*r4
f(367) = -a1**2*a2*r2 + a1*a2**2*r1 - a3**2*a4*r4 + a3*a4**2*r3
f(368) = a1**3*a2*r2 - a1*a2**3*r1 + a3**3*a4*r4 - a3*a4**3*r3
f(369) = -a1**4*a2*r2 + a1*a2**4*r1 - a3**4*a4*r4 + a3*a4**4*r3
f(370) = a1**2*a2*r1 - a1*a2**2*r2 + a3**2*a4*r3 - a3*a4**2*r4
f(371) = a1**2*a2**2*r1 - a1**2*a2**2*r2 + a3**2*a4**2*r3 - a3**2*a4**2* &
      r4
f(372) = -a1**3*a2**2*r2 + a1**2*a2**3*r1 - a3**3*a4**2*r4 + a3**2*a4**3 &
      *r3
f(373) = -a1**3*a2*r1 + a1*a2**3*r2 - a3**3*a4*r3 + a3*a4**3*r4
f(374) = -a1**3*a2**2*r1 + a1**2*a2**3*r2 - a3**3*a4**2*r3 + a3**2*a4**3 &
      *r4
f(375) = -a1**4*a2*r1 + a1*a2**4*r2 - a3**4*a4*r3 + a3*a4**4*r4
f(376) = a1*a2*r1**2 - a1*a2*r2**2 + a3*a4*r3**2 - a3*a4*r4**2
f(377) = a1**2*a2*r2**2 - a1*a2**2*r1**2 + a3**2*a4*r4**2 - a3*a4**2*r3 &
      **2
f(378) = -a1**3*a2*r2**2 + a1*a2**3*r1**2 - a3**3*a4*r4**2 + a3*a4**3*r3 &
      **2
f(379) = a1**2*a2*r1**2 - a1*a2**2*r2**2 + a3**2*a4*r3**2 - a3*a4**2*r4 &
      **2
f(380) = -a1**2*a2**2*r1**2 + a1**2*a2**2*r2**2 - a3**2*a4**2*r3**2 + a3 &
      **2*a4**2*r4**2
f(381) = -a1**3*a2*r1**2 + a1*a2**3*r2**2 - a3**3*a4*r3**2 + a3*a4**3*r4 &
      **2
f(382) = a1*a2*r1**3 - a1*a2*r2**3 + a3*a4*r3**3 - a3*a4*r4**3
f(383) = -a1**2*a2*r2**3 + a1*a2**2*r1**3 - a3**2*a4*r4**3 + a3*a4**2*r3 &
      **3
f(384) = -a1**2*a2*r1**3 + a1*a2**2*r2**3 - a3**2*a4*r3**3 + a3*a4**2*r4 &
      **3
f(385) = a1*a2*r1**4 - a1*a2*r2**4 + a3*a4*r3**4 - a3*a4*r4**4
f(386) = -a1*a3*r1 - a1*a3*r3 + a2*a4*r2 + a2*a4*r4
f(387) = -a1**2*a3*r3 - a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2*r2
f(388) = a1**3*a3*r3 + a1*a3**3*r1 - a2**3*a4*r4 - a2*a4**3*r2
f(389) = -a1**4*a3*r3 - a1*a3**4*r1 + a2**4*a4*r4 + a2*a4**4*r2
f(390) = -a1**2*a3*r1 - a1*a3**2*r3 + a2**2*a4*r2 + a2*a4**2*r4
f(391) = a1**2*a3**2*r1 + a1**2*a3**2*r3 - a2**2*a4**2*r2 - a2**2*a4**2* &
      r4
f(392) = a1**3*a3**2*r3 + a1**2*a3**3*r1 - a2**3*a4**2*r4 - a2**2*a4**3* &
      r2
f(393) = -a1**3*a3*r1 - a1*a3**3*r3 + a2**3*a4*r2 + a2*a4**3*r4
f(394) = a1**3*a3**2*r1 + a1**2*a3**3*r3 - a2**3*a4**2*r2 - a2**2*a4**3* &
      r4
f(395) = -a1**4*a3*r1 - a1*a3**4*r3 + a2**4*a4*r2 + a2*a4**4*r4
f(396) = -a1*a3*r1**2 - a1*a3*r3**2 + a2*a4*r2**2 + a2*a4*r4**2
f(397) = a1**2*a3*r3**2 + a1*a3**2*r1**2 - a2**2*a4*r4**2 - a2*a4**2*r2 &
      **2
f(398) = -a1**3*a3*r3**2 - a1*a3**3*r1**2 + a2**3*a4*r4**2 + a2*a4**3*r2 &
      **2
f(399) = a1**2*a3*r1**2 + a1*a3**2*r3**2 - a2**2*a4*r2**2 - a2*a4**2*r4 &
      **2
f(400) = -a1**2*a3**2*r1**2 - a1**2*a3**2*r3**2 + a2**2*a4**2*r2**2 + a2 &
      **2*a4**2*r4**2
f(401) = a1**3*a3*r1**2 + a1*a3**3*r3**2 - a2**3*a4*r2**2 - a2*a4**3*r4 &
      **2
f(402) = a1*a3*r1**3 + a1*a3*r3**3 - a2*a4*r2**3 - a2*a4*r4**3
f(403) = a1**2*a3*r3**3 + a1*a3**2*r1**3 - a2**2*a4*r4**3 - a2*a4**2*r2 &
      **3
f(404) = a1**2*a3*r1**3 + a1*a3**2*r3**3 - a2**2*a4*r2**3 - a2*a4**2*r4 &
      **3
f(405) = a1*a3*r1**4 + a1*a3*r3**4 - a2*a4*r2**4 - a2*a4*r4**4
f(406) = -a1*a4*r1 + a1*a4*r4 + a2*a3*r2 - a2*a3*r3
f(407) = -a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 - a2*a3**2*r2
f(408) = a1**3*a4*r4 - a1*a4**3*r1 - a2**3*a3*r3 + a2*a3**3*r2
f(409) = -a1**4*a4*r4 + a1*a4**4*r1 + a2**4*a3*r3 - a2*a3**4*r2
f(410) = a1**2*a4*r1 - a1*a4**2*r4 - a2**2*a3*r2 + a2*a3**2*r3
f(411) = a1**2*a4**2*r1 - a1**2*a4**2*r4 - a2**2*a3**2*r2 + a2**2*a3**2* &
      r3
f(412) = -a1**3*a4**2*r4 + a1**2*a4**3*r1 + a2**3*a3**2*r3 - a2**2*a3**3 &
      *r2
f(413) = a1**3*a4*r1 - a1*a4**3*r4 - a2**3*a3*r2 + a2*a3**3*r3
f(414) = a1**3*a4**2*r1 - a1**2*a4**3*r4 - a2**3*a3**2*r2 + a2**2*a3**3* &
      r3
f(415) = -a1**4*a4*r1 + a1*a4**4*r4 + a2**4*a3*r2 - a2*a3**4*r3
f(416) = a1*a4*r1**2 - a1*a4*r4**2 - a2*a3*r2**2 + a2*a3*r3**2
f(417) = a1**2*a4*r4**2 - a1*a4**2*r1**2 - a2**2*a3*r3**2 + a2*a3**2*r2 &
      **2
f(418) = -a1**3*a4*r4**2 + a1*a4**3*r1**2 + a2**3*a3*r3**2 - a2*a3**3*r2 &
      **2
f(419) = a1**2*a4*r1**2 - a1*a4**2*r4**2 - a2**2*a3*r2**2 + a2*a3**2*r3 &
      **2
f(420) = -a1**2*a4**2*r1**2 + a1**2*a4**2*r4**2 + a2**2*a3**2*r2**2 - a2 &
      **2*a3**2*r3**2
f(421) = a1**3*a4*r1**2 - a1*a4**3*r4**2 - a2**3*a3*r2**2 + a2*a3**3*r3 &
      **2
f(422) = -a1*a4*r1**3 + a1*a4*r4**3 + a2*a3*r2**3 - a2*a3*r3**3
f(423) = a1**2*a4*r4**3 - a1*a4**2*r1**3 - a2**2*a3*r3**3 + a2*a3**2*r2 &
      **3
f(424) = a1**2*a4*r1**3 - a1*a4**2*r4**3 - a2**2*a3*r2**3 + a2*a3**2*r3 &
      **3
f(425) = a1*a4*r1**4 - a1*a4*r4**4 - a2*a3*r2**4 + a2*a3*r3**4
f(426) = a1*b1**2*r1 - a2*b1**2*r2 + a3*b2**2*r3 - a4*b2**2*r4
f(427) = a1*b1**4*r1 - a2*b1**4*r2 + a3*b2**4*r3 - a4*b2**4*r4
f(428) = -a1**2*b1**2*r1 + a2**2*b1**2*r2 - a3**2*b2**2*r3 + a4**2*b2**2 &
      *r4
f(429) = a1**3*b1**2*r1 - a2**3*b1**2*r2 + a3**3*b2**2*r3 - a4**3*b2**2* &
      r4
f(430) = a1*b1**2*r1**2 - a2*b1**2*r2**2 + a3*b2**2*r3**2 - a4*b2**2*r4 &
      **2
f(431) = -a1**2*b1**2*r1**2 + a2**2*b1**2*r2**2 - a3**2*b2**2*r3**2 + a4 &
      **2*b2**2*r4**2
f(432) = -a1*b1**2*r1**3 + a2*b1**2*r2**3 - a3*b2**2*r3**3 + a4*b2**2*r4 &
      **3
f(433) = a1*b2**2*r1 - a2*b2**2*r2 + a3*b1**2*r3 - a4*b1**2*r4
f(434) = a1*b2**4*r1 - a2*b2**4*r2 + a3*b1**4*r3 - a4*b1**4*r4
f(435) = -a1**2*b2**2*r1 + a2**2*b2**2*r2 - a3**2*b1**2*r3 + a4**2*b1**2 &
      *r4
f(436) = -a1**3*b2**2*r1 + a2**3*b2**2*r2 - a3**3*b1**2*r3 + a4**3*b1**2 &
      *r4
f(437) = a1*b2**2*r1**2 - a2*b2**2*r2**2 + a3*b1**2*r3**2 - a4*b1**2*r4 &
      **2
f(438) = -a1**2*b2**2*r1**2 + a2**2*b2**2*r2**2 - a3**2*b1**2*r3**2 + a4 &
      **2*b1**2*r4**2
f(439) = a1*b2**2*r1**3 - a2*b2**2*r2**3 + a3*b1**2*r3**3 - a4*b1**2*r4 &
      **3
f(440) = dtau**2*(-a1*r1 + a2*r2 - a3*r3 + a4*r4)
f(441) = dtau**4*(-a1*r1 + a2*r2 - a3*r3 + a4*r4)
f(442) = dtau**2*(a1**2*r1 - a2**2*r2 + a3**2*r3 - a4**2*r4)
f(443) = dtau**2*(-a1**3*r1 + a2**3*r2 - a3**3*r3 + a4**3*r4)
f(444) = dtau**2*(a1*r1**2 - a2*r2**2 + a3*r3**2 - a4*r4**2)
f(445) = dtau**2*(a1**2*r1**2 - a2**2*r2**2 + a3**2*r3**2 - a4**2*r4**2)
f(446) = dtau**2*(a1*r1**3 - a2*r2**3 + a3*r3**3 - a4*r4**3)
f(447) = -a1*a4*r2 + a1*a4*r3 + a2*a3*r1 - a2*a3*r4
f(448) = -a1**2*a4*r3 + a1*a4**2*r2 + a2**2*a3*r4 - a2*a3**2*r1
f(449) = a1**3*a4*r3 - a1*a4**3*r2 - a2**3*a3*r4 + a2*a3**3*r1
f(450) = a1**4*a4*r3 - a1*a4**4*r2 - a2**4*a3*r4 + a2*a3**4*r1
f(451) = a1**2*a4*r2 - a1*a4**2*r3 - a2**2*a3*r1 + a2*a3**2*r4
f(452) = -a1**2*a4**2*r2 + a1**2*a4**2*r3 + a2**2*a3**2*r1 - a2**2*a3**2 &
      *r4
f(453) = a1**3*a4**2*r3 - a1**2*a4**3*r2 - a2**3*a3**2*r4 + a2**2*a3**3* &
      r1
f(454) = a1**3*a4*r2 - a1*a4**3*r3 - a2**3*a3*r1 + a2*a3**3*r4
f(455) = -a1**3*a4**2*r2 + a1**2*a4**3*r3 + a2**3*a3**2*r1 - a2**2*a3**3 &
      *r4
f(456) = a1**4*a4*r2 - a1*a4**4*r3 - a2**4*a3*r1 + a2*a3**4*r4
f(457) = a1*a4*r2**2 - a1*a4*r3**2 - a2*a3*r1**2 + a2*a3*r4**2
f(458) = -a1**2*a4*r3**2 + a1*a4**2*r2**2 + a2**2*a3*r4**2 - a2*a3**2*r1 &
      **2
f(459) = a1**3*a4*r3**2 - a1*a4**3*r2**2 - a2**3*a3*r4**2 + a2*a3**3*r1 &
      **2
f(460) = -a1**2*a4*r2**2 + a1*a4**2*r3**2 + a2**2*a3*r1**2 - a2*a3**2*r4 &
      **2
f(461) = -a1**2*a4**2*r2**2 + a1**2*a4**2*r3**2 + a2**2*a3**2*r1**2 - a2 &
      **2*a3**2*r4**2
f(462) = -a1**3*a4*r2**2 + a1*a4**3*r3**2 + a2**3*a3*r1**2 - a2*a3**3*r4 &
      **2
f(463) = -a1*a4*r2**3 + a1*a4*r3**3 + a2*a3*r1**3 - a2*a3*r4**3
f(464) = -a1**2*a4*r3**3 + a1*a4**2*r2**3 + a2**2*a3*r4**3 - a2*a3**2*r1 &
      **3
f(465) = -a1**2*a4*r2**3 + a1*a4**2*r3**3 + a2**2*a3*r1**3 - a2*a3**2*r4 &
      **3
f(466) = a1*a4*r2**4 - a1*a4*r3**4 - a2*a3*r1**4 + a2*a3*r4**4
f(467) = -a1*a3*r2 - a1*a3*r4 + a2*a4*r1 + a2*a4*r3
f(468) = -a1**2*a3*r4 - a1*a3**2*r2 + a2**2*a4*r3 + a2*a4**2*r1
f(469) = a1**3*a3*r4 + a1*a3**3*r2 - a2**3*a4*r3 - a2*a4**3*r1
f(470) = a1**4*a3*r4 + a1*a3**4*r2 - a2**4*a4*r3 - a2*a4**4*r1
f(471) = -a1**2*a3*r2 - a1*a3**2*r4 + a2**2*a4*r1 + a2*a4**2*r3
f(472) = -a1**2*a3**2*r2 - a1**2*a3**2*r4 + a2**2*a4**2*r1 + a2**2*a4**2 &
      *r3
f(473) = -a1**3*a3**2*r4 - a1**2*a3**3*r2 + a2**3*a4**2*r3 + a2**2*a4**3 &
      *r1
f(474) = -a1**3*a3*r2 - a1*a3**3*r4 + a2**3*a4*r1 + a2*a4**3*r3
f(475) = -a1**3*a3**2*r2 - a1**2*a3**3*r4 + a2**3*a4**2*r1 + a2**2*a4**3 &
      *r3
f(476) = a1**4*a3*r2 + a1*a3**4*r4 - a2**4*a4*r1 - a2*a4**4*r3
f(477) = -a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3**2
f(478) = -a1**2*a3*r4**2 - a1*a3**2*r2**2 + a2**2*a4*r3**2 + a2*a4**2*r1 &
      **2
f(479) = a1**3*a3*r4**2 + a1*a3**3*r2**2 - a2**3*a4*r3**2 - a2*a4**3*r1 &
      **2
f(480) = -a1**2*a3*r2**2 - a1*a3**2*r4**2 + a2**2*a4*r1**2 + a2*a4**2*r3 &
      **2
f(481) = -a1**2*a3**2*r2**2 - a1**2*a3**2*r4**2 + a2**2*a4**2*r1**2 + a2 &
      **2*a4**2*r3**2
f(482) = -a1**3*a3*r2**2 - a1*a3**3*r4**2 + a2**3*a4*r1**2 + a2*a4**3*r3 &
      **2
f(483) = a1*a3*r2**3 + a1*a3*r4**3 - a2*a4*r1**3 - a2*a4*r3**3
f(484) = -a1**2*a3*r4**3 - a1*a3**2*r2**3 + a2**2*a4*r3**3 + a2*a4**2*r1 &
      **3
f(485) = -a1**2*a3*r2**3 - a1*a3**2*r4**3 + a2**2*a4*r1**3 + a2*a4**2*r3 &
      **3
f(486) = a1*a3*r2**4 + a1*a3*r4**4 - a2*a4*r1**4 - a2*a4*r3**4
f(487) = a1*b1**2*r2 - a2*b1**2*r1 + a3*b2**2*r4 - a4*b2**2*r3
f(488) = a1*b1**4*r2 - a2*b1**4*r1 + a3*b2**4*r4 - a4*b2**4*r3
f(489) = -a1**2*b1**2*r2 + a2**2*b1**2*r1 - a3**2*b2**2*r4 + a4**2*b2**2 &
      *r3
f(490) = a1**3*b1**2*r2 - a2**3*b1**2*r1 + a3**3*b2**2*r4 - a4**3*b2**2* &
      r3
f(491) = a1*b1**2*r2**2 - a2*b1**2*r1**2 + a3*b2**2*r4**2 - a4*b2**2*r3 &
      **2
f(492) = -a1**2*b1**2*r2**2 + a2**2*b1**2*r1**2 - a3**2*b2**2*r4**2 + a4 &
      **2*b2**2*r3**2
f(493) = a1*b1**2*r2**3 - a2*b1**2*r1**3 + a3*b2**2*r4**3 - a4*b2**2*r3 &
      **3
f(494) = a1*b2**2*r2 - a2*b2**2*r1 + a3*b1**2*r4 - a4*b1**2*r3
f(495) = a1*b2**4*r2 - a2*b2**4*r1 + a3*b1**4*r4 - a4*b1**4*r3
f(496) = a1**2*b2**2*r2 - a2**2*b2**2*r1 + a3**2*b1**2*r4 - a4**2*b1**2* &
      r3
f(497) = a1**3*b2**2*r2 - a2**3*b2**2*r1 + a3**3*b1**2*r4 - a4**3*b1**2* &
      r3
f(498) = -a1*b2**2*r2**2 + a2*b2**2*r1**2 - a3*b1**2*r4**2 + a4*b1**2*r3 &
      **2
f(499) = a1**2*b2**2*r2**2 - a2**2*b2**2*r1**2 + a3**2*b1**2*r4**2 - a4 &
      **2*b1**2*r3**2
f(500) = a1*b2**2*r2**3 - a2*b2**2*r1**3 + a3*b1**2*r4**3 - a4*b1**2*r3 &
      **3
f(501) = dtau**2*(a1*r2 - a2*r1 + a3*r4 - a4*r3)
f(502) = dtau**4*(a1*r2 - a2*r1 + a3*r4 - a4*r3)
f(503) = dtau**2*(a1**2*r2 - a2**2*r1 + a3**2*r4 - a4**2*r3)
f(504) = dtau**2*(-a1**3*r2 + a2**3*r1 - a3**3*r4 + a4**3*r3)
f(505) = dtau**2*(a1*r2**2 - a2*r1**2 + a3*r4**2 - a4*r3**2)
f(506) = dtau**2*(a1**2*r2**2 - a2**2*r1**2 + a3**2*r4**2 - a4**2*r3**2)
f(507) = dtau**2*(-a1*r2**3 + a2*r1**3 - a3*r4**3 + a4*r3**3)
f(508) = a1*a2*r3 - a1*a2*r4 + a3*a4*r1 - a3*a4*r2
f(509) = a1**2*a2*r4 - a1*a2**2*r3 + a3**2*a4*r2 - a3*a4**2*r1
f(510) = a1**3*a2*r4 - a1*a2**3*r3 + a3**3*a4*r2 - a3*a4**3*r1
f(511) = a1**4*a2*r4 - a1*a2**4*r3 + a3**4*a4*r2 - a3*a4**4*r1
f(512) = -a1**2*a2*r3 + a1*a2**2*r4 - a3**2*a4*r1 + a3*a4**2*r2
f(513) = -a1**2*a2**2*r3 + a1**2*a2**2*r4 - a3**2*a4**2*r1 + a3**2*a4**2 &
      *r2
f(514) = -a1**3*a2**2*r4 + a1**2*a2**3*r3 - a3**3*a4**2*r2 + a3**2*a4**3 &
      *r1
f(515) = -a1**3*a2*r3 + a1*a2**3*r4 - a3**3*a4*r1 + a3*a4**3*r2
f(516) = -a1**3*a2**2*r3 + a1**2*a2**3*r4 - a3**3*a4**2*r1 + a3**2*a4**3 &
      *r2
f(517) = a1**4*a2*r3 - a1*a2**4*r4 + a3**4*a4*r1 - a3*a4**4*r2
f(518) = -a1*a2*r3**2 + a1*a2*r4**2 - a3*a4*r1**2 + a3*a4*r2**2
f(519) = -a1**2*a2*r4**2 + a1*a2**2*r3**2 - a3**2*a4*r2**2 + a3*a4**2*r1 &
      **2
f(520) = -a1**3*a2*r4**2 + a1*a2**3*r3**2 - a3**3*a4*r2**2 + a3*a4**3*r1 &
      **2
f(521) = -a1**2*a2*r3**2 + a1*a2**2*r4**2 - a3**2*a4*r1**2 + a3*a4**2*r2 &
      **2
f(522) = a1**2*a2**2*r3**2 - a1**2*a2**2*r4**2 + a3**2*a4**2*r1**2 - a3 &
      **2*a4**2*r2**2
f(523) = -a1**3*a2*r3**2 + a1*a2**3*r4**2 - a3**3*a4*r1**2 + a3*a4**3*r2 &
      **2
f(524) = -a1*a2*r3**3 + a1*a2*r4**3 - a3*a4*r1**3 + a3*a4*r2**3
f(525) = a1**2*a2*r4**3 - a1*a2**2*r3**3 + a3**2*a4*r2**3 - a3*a4**2*r1 &
      **3
f(526) = a1**2*a2*r3**3 - a1*a2**2*r4**3 + a3**2*a4*r1**3 - a3*a4**2*r2 &
      **3
f(527) = -a1*a2*r3**4 + a1*a2*r4**4 - a3*a4*r1**4 + a3*a4*r2**4
f(528) = a1*b2**2*r3 - a2*b2**2*r4 + a3*b1**2*r1 - a4*b1**2*r2
f(529) = -a1*b2**4*r3 + a2*b2**4*r4 - a3*b1**4*r1 + a4*b1**4*r2
f(530) = a1**2*b2**2*r3 - a2**2*b2**2*r4 + a3**2*b1**2*r1 - a4**2*b1**2* &
      r2
f(531) = -a1**3*b2**2*r3 + a2**3*b2**2*r4 - a3**3*b1**2*r1 + a4**3*b1**2 &
      *r2
f(532) = a1*b2**2*r3**2 - a2*b2**2*r4**2 + a3*b1**2*r1**2 - a4*b1**2*r2 &
      **2
f(533) = -a1**2*b2**2*r3**2 + a2**2*b2**2*r4**2 - a3**2*b1**2*r1**2 + a4 &
      **2*b1**2*r2**2
f(534) = -a1*b2**2*r3**3 + a2*b2**2*r4**3 - a3*b1**2*r1**3 + a4*b1**2*r2 &
      **3
f(535) = a1*b1**2*r3 - a2*b1**2*r4 + a3*b2**2*r1 - a4*b2**2*r2
f(536) = -a1*b1**4*r3 + a2*b1**4*r4 - a3*b2**4*r1 + a4*b2**4*r2
f(537) = -a1**2*b1**2*r3 + a2**2*b1**2*r4 - a3**2*b2**2*r1 + a4**2*b2**2 &
      *r2
f(538) = a1**3*b1**2*r3 - a2**3*b1**2*r4 + a3**3*b2**2*r1 - a4**3*b2**2* &
      r2
f(539) = a1*b1**2*r3**2 - a2*b1**2*r4**2 + a3*b2**2*r1**2 - a4*b2**2*r2 &
      **2
f(540) = a1**2*b1**2*r3**2 - a2**2*b1**2*r4**2 + a3**2*b2**2*r1**2 - a4 &
      **2*b2**2*r2**2
f(541) = a1*b1**2*r3**3 - a2*b1**2*r4**3 + a3*b2**2*r1**3 - a4*b2**2*r2 &
      **3
f(542) = dtau**2*(-a1*r3 + a2*r4 - a3*r1 + a4*r2)
f(543) = dtau**4*(a1*r3 - a2*r4 + a3*r1 - a4*r2)
f(544) = dtau**2*(-a1**2*r3 + a2**2*r4 - a3**2*r1 + a4**2*r2)
f(545) = dtau**2*(a1**3*r3 - a2**3*r4 + a3**3*r1 - a4**3*r2)
f(546) = dtau**2*(-a1*r3**2 + a2*r4**2 - a3*r1**2 + a4*r2**2)
f(547) = dtau**2*(a1**2*r3**2 - a2**2*r4**2 + a3**2*r1**2 - a4**2*r2**2)
f(548) = dtau**2*(a1*r3**3 - a2*r4**3 + a3*r1**3 - a4*r2**3)
f(549) = a1*b2**2*r4 - a2*b2**2*r3 + a3*b1**2*r2 - a4*b1**2*r1
f(550) = -a1*b2**4*r4 + a2*b2**4*r3 - a3*b1**4*r2 + a4*b1**4*r1
f(551) = a1**2*b2**2*r4 - a2**2*b2**2*r3 + a3**2*b1**2*r2 - a4**2*b1**2* &
      r1
f(552) = -a1**3*b2**2*r4 + a2**3*b2**2*r3 - a3**3*b1**2*r2 + a4**3*b1**2 &
      *r1
f(553) = a1*b2**2*r4**2 - a2*b2**2*r3**2 + a3*b1**2*r2**2 - a4*b1**2*r1 &
      **2
f(554) = a1**2*b2**2*r4**2 - a2**2*b2**2*r3**2 + a3**2*b1**2*r2**2 - a4 &
      **2*b1**2*r1**2
f(555) = -a1*b2**2*r4**3 + a2*b2**2*r3**3 - a3*b1**2*r2**3 + a4*b1**2*r1 &
      **3
f(556) = -a1*b1**2*r4 + a2*b1**2*r3 - a3*b2**2*r2 + a4*b2**2*r1
f(557) = a1*b1**4*r4 - a2*b1**4*r3 + a3*b2**4*r2 - a4*b2**4*r1
f(558) = -a1**2*b1**2*r4 + a2**2*b1**2*r3 - a3**2*b2**2*r2 + a4**2*b2**2 &
      *r1
f(559) = a1**3*b1**2*r4 - a2**3*b1**2*r3 + a3**3*b2**2*r2 - a4**3*b2**2* &
      r1
f(560) = -a1*b1**2*r4**2 + a2*b1**2*r3**2 - a3*b2**2*r2**2 + a4*b2**2*r1 &
      **2
f(561) = a1**2*b1**2*r4**2 - a2**2*b1**2*r3**2 + a3**2*b2**2*r2**2 - a4 &
      **2*b2**2*r1**2
f(562) = -a1*b1**2*r4**3 + a2*b1**2*r3**3 - a3*b2**2*r2**3 + a4*b2**2*r1 &
      **3
f(563) = dtau**2*(-a1*r4 + a2*r3 - a3*r2 + a4*r1)
f(564) = dtau**4*(a1*r4 - a2*r3 + a3*r2 - a4*r1)
f(565) = dtau**2*(-a1**2*r4 + a2**2*r3 - a3**2*r2 + a4**2*r1)
f(566) = dtau**2*(a1**3*r4 - a2**3*r3 + a3**3*r2 - a4**3*r1)
f(567) = dtau**2*(-a1*r4**2 + a2*r3**2 - a3*r2**2 + a4*r1**2)
f(568) = dtau**2*(-a1**2*r4**2 + a2**2*r3**2 - a3**2*r2**2 + a4**2*r1**2 &
      )
f(569) = dtau**2*(-a1*r4**3 + a2*r3**3 - a3*r2**3 + a4*r1**3)
f(570) = b1*b2*(r1 - r2 + r3 - r4)
f(571) = b1*b2*(b1**2*r3 - b1**2*r4 + b2**2*r1 - b2**2*r2)
f(572) = b1**2*b2**2*(-r1 + r2 - r3 + r4)
f(573) = b1*b2*(b1**2*r1 - b1**2*r2 + b2**2*r3 - b2**2*r4)
f(574) = b1*b2*(r1**2 - r2**2 + r3**2 - r4**2)
f(575) = b1*b2*(-b1**2*r3**2 + b1**2*r4**2 - b2**2*r1**2 + b2**2*r2**2)
f(576) = b1**2*b2**2*(r1**2 - r2**2 + r3**2 - r4**2)
f(577) = b1*b2*(b1**2*r1**2 - b1**2*r2**2 + b2**2*r3**2 - b2**2*r4**2)
f(578) = b1*b2*(-r1**3 + r2**3 - r3**3 + r4**3)
f(579) = b1*b2*(r1**4 - r2**4 + r3**4 - r4**4)
f(580) = dtau*(b1*r1 + b1*r2 + b2*r3 + b2*r4)
f(581) = dtau**3*(b1*r1 + b1*r2 + b2*r3 + b2*r4)
f(582) = dtau**2*(b1**2*r1 - b1**2*r2 + b2**2*r3 - b2**2*r4)
f(583) = dtau*(b1**3*r1 + b1**3*r2 + b2**3*r3 + b2**3*r4)
f(584) = dtau*(b1*r1**2 + b1*r2**2 + b2*r3**2 + b2*r4**2)
f(585) = dtau**3*(b1*r1**2 + b1*r2**2 + b2*r3**2 + b2*r4**2)
f(586) = dtau**2*(b1**2*r1**2 - b1**2*r2**2 + b2**2*r3**2 - b2**2*r4**2)
f(587) = dtau*(b1**3*r1**2 + b1**3*r2**2 + b2**3*r3**2 + b2**3*r4**2)
f(588) = dtau*(b1*r1**3 + b1*r2**3 + b2*r3**3 + b2*r4**3)
f(589) = dtau*(b1*r1**4 + b1*r2**4 + b2*r3**4 + b2*r4**4)
f(590) = dtau*(b1*r3 + b1*r4 + b2*r1 + b2*r2)
f(591) = dtau**3*(b1*r3 + b1*r4 + b2*r1 + b2*r2)
f(592) = dtau**2*(b1**2*r3 - b1**2*r4 + b2**2*r1 - b2**2*r2)
f(593) = dtau*(b1**3*r3 + b1**3*r4 + b2**3*r1 + b2**3*r2)
f(594) = dtau*(b1*r3**2 + b1*r4**2 + b2*r1**2 + b2*r2**2)
f(595) = dtau**3*(b1*r3**2 + b1*r4**2 + b2*r1**2 + b2*r2**2)
f(596) = dtau**2*(b1**2*r3**2 - b1**2*r4**2 + b2**2*r1**2 - b2**2*r2**2)
f(597) = dtau*(b1**3*r3**2 + b1**3*r4**2 + b2**3*r1**2 + b2**3*r2**2)
f(598) = dtau*(b1*r3**3 + b1*r4**3 + b2*r1**3 + b2*r2**3)
f(599) = dtau*(b1*r3**4 + b1*r4**4 + b2*r1**4 + b2*r2**4)
f(600) = -a1*a2*a3 + a1*a2*a4 - a1*a3*a4 + a2*a3*a4
f(601) = a1**2*a3*a4 + a1*a2*a3**2 - a1*a2*a4**2 - a2**2*a3*a4
f(602) = -a1**3*a3*a4 - a1*a2*a3**3 + a1*a2*a4**3 + a2**3*a3*a4
f(603) = -a1**4*a3*a4 - a1*a2*a3**4 + a1*a2*a4**4 + a2**4*a3*a4
f(604) = -a1**2*a2*a4 + a1*a2**2*a3 + a1*a3*a4**2 - a2*a3**2*a4
f(605) = a1**2*a2*a4**2 - a1**2*a3*a4**2 - a1*a2**2*a3**2 + a2**2*a3**2* &
      a4
f(606) = -a1**3*a3*a4**2 + a1**2*a2*a4**3 - a1*a2**2*a3**3 + a2**3*a3**2 &
      *a4
f(607) = -a1**3*a2*a4 + a1*a2**3*a3 + a1*a3*a4**3 - a2*a3**3*a4
f(608) = a1**3*a2*a4**2 - a1**2*a3*a4**3 - a1*a2**3*a3**2 + a2**2*a3**3* &
      a4
f(609) = -a1**4*a2*a4 + a1*a2**4*a3 + a1*a3*a4**4 - a2*a3**4*a4
f(610) = a1**2*a2*a3 - a1*a2**2*a4 + a1*a3**2*a4 - a2*a3*a4**2
f(611) = a1**2*a2*a3**2 + a1**2*a3**2*a4 - a1*a2**2*a4**2 - a2**2*a3*a4 &
      **2
f(612) = a1**3*a3**2*a4 + a1**2*a2*a3**3 - a1*a2**2*a4**3 - a2**3*a3*a4 &
      **2
f(613) = -a1**2*a2**2*a3 + a1**2*a2**2*a4 - a1*a3**2*a4**2 + a2*a3**2*a4 &
      **2
f(614) = -a1**2*a2**2*a3**2 + a1**2*a2**2*a4**2 - a1**2*a3**2*a4**2 + a2 &
      **2*a3**2*a4**2
f(615) = a1**3*a2**2*a4 - a1**2*a2**3*a3 - a1*a3**2*a4**3 + a2*a3**3*a4 &
      **2
f(616) = a1**3*a2*a3 - a1*a2**3*a4 + a1*a3**3*a4 - a2*a3*a4**3
f(617) = a1**3*a2*a3**2 + a1**2*a3**3*a4 - a1*a2**3*a4**2 - a2**2*a3*a4 &
      **3
f(618) = a1**3*a2**2*a3 - a1**2*a2**3*a4 + a1*a3**3*a4**2 - a2*a3**2*a4 &
      **3
f(619) = -a1**4*a2*a3 + a1*a2**4*a4 - a1*a3**4*a4 + a2*a3*a4**4
f(620) = a1**2*a2*b1**2 - a1*a2**2*b1**2 + a3**2*a4*b2**2 - a3*a4**2*b2 &
      **2
f(621) = -a1**3*a2*b1**2 + a1*a2**3*b1**2 - a3**3*a4*b2**2 + a3*a4**3*b2 &
      **2
f(622) = -a1**2*a2*b2**2 + a1*a2**2*b2**2 - a3**2*a4*b1**2 + a3*a4**2*b1 &
      **2
f(623) = -a1**3*a2*b2**2 + a1*a2**3*b2**2 - a3**3*a4*b1**2 + a3*a4**3*b1 &
      **2
f(624) = dtau**2*(a1**2*a2 - a1*a2**2 + a3**2*a4 - a3*a4**2)
f(625) = dtau**2*(a1**3*a2 - a1*a2**3 + a3**3*a4 - a3*a4**3)
f(626) = -a1*a3*b1**2 - a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2**2
f(627) = a1*a3*b1**4 + a1*a3*b2**4 - a2*a4*b1**4 - a2*a4*b2**4
f(628) = a1**2*a3*b2**2 + a1*a3**2*b1**2 - a2**2*a4*b2**2 - a2*a4**2*b1 &
      **2
f(629) = a1**3*a3*b2**2 + a1*a3**3*b1**2 - a2**3*a4*b2**2 - a2*a4**3*b1 &
      **2
f(630) = a1**2*a3*b1**2 + a1*a3**2*b2**2 - a2**2*a4*b1**2 - a2*a4**2*b2 &
      **2
f(631) = -a1**2*a3**2*b1**2 - a1**2*a3**2*b2**2 + a2**2*a4**2*b1**2 + a2 &
      **2*a4**2*b2**2
f(632) = a1**3*a3*b1**2 + a1*a3**3*b2**2 - a2**3*a4*b1**2 - a2*a4**3*b2 &
      **2
f(633) = dtau**2*(a1*a3 - a2*a4)
f(634) = dtau**4*(-a1*a3 + a2*a4)
f(635) = dtau**2*(-a1**2*a3 - a1*a3**2 + a2**2*a4 + a2*a4**2)
f(636) = dtau**2*(a1**3*a3 + a1*a3**3 - a2**3*a4 - a2*a4**3)
f(637) = dtau**2*(a1**2*a3**2 - a2**2*a4**2)
f(638) = a1*a4*b1**2 - a1*a4*b2**2 - a2*a3*b1**2 + a2*a3*b2**2
f(639) = -a1*a4*b1**4 + a1*a4*b2**4 + a2*a3*b1**4 - a2*a3*b2**4
f(640) = -a1**2*a4*b2**2 + a1*a4**2*b1**2 + a2**2*a3*b2**2 - a2*a3**2*b1 &
      **2
f(641) = -a1**3*a4*b2**2 + a1*a4**3*b1**2 + a2**3*a3*b2**2 - a2*a3**3*b1 &
      **2
f(642) = -a1**2*a4*b1**2 + a1*a4**2*b2**2 + a2**2*a3*b1**2 - a2*a3**2*b2 &
      **2
f(643) = -a1**2*a4**2*b1**2 + a1**2*a4**2*b2**2 + a2**2*a3**2*b1**2 - a2 &
      **2*a3**2*b2**2
f(644) = -a1**3*a4*b1**2 + a1*a4**3*b2**2 + a2**3*a3*b1**2 - a2*a3**3*b2 &
      **2
f(645) = dtau**2*(a1**2*a4 - a1*a4**2 - a2**2*a3 + a2*a3**2)
f(646) = dtau**2*(-a1**3*a4 + a1*a4**3 + a2**3*a3 - a2*a3**3)
f(647) = b1*b2*(-a1 + a2 - a3 + a4)
f(648) = b1*b2*(-a1*b2**2 + a2*b2**2 - a3*b1**2 + a4*b1**2)
f(649) = b1**2*b2**2*(a1 - a2 + a3 - a4)
f(650) = b1*b2*(-a1*b1**2 + a2*b1**2 - a3*b2**2 + a4*b2**2)
f(651) = b1*b2*(-a1**2 + a2**2 - a3**2 + a4**2)
f(652) = b1*b2*(a1**2*b2**2 - a2**2*b2**2 + a3**2*b1**2 - a4**2*b1**2)
f(653) = b1**2*b2**2*(-a1**2 + a2**2 - a3**2 + a4**2)
f(654) = b1*b2*(a1**2*b1**2 - a2**2*b1**2 + a3**2*b2**2 - a4**2*b2**2)
f(655) = b1*b2*(a1**3 - a2**3 + a3**3 - a4**3)
f(656) = b1*b2*(-a1**4 + a2**4 - a3**4 + a4**4)
f(657) = dtau*(a1*b1 + a2*b1 + a3*b2 + a4*b2)
f(658) = dtau**3*(a1*b1 + a2*b1 + a3*b2 + a4*b2)
f(659) = dtau**2*(a1*b1**2 - a2*b1**2 + a3*b2**2 - a4*b2**2)
f(660) = dtau*(a1*b1**3 + a2*b1**3 + a3*b2**3 + a4*b2**3)
f(661) = dtau*(a1**2*b1 + a2**2*b1 + a3**2*b2 + a4**2*b2)
f(662) = dtau**3*(a1**2*b1 + a2**2*b1 + a3**2*b2 + a4**2*b2)
f(663) = dtau**2*(a1**2*b1**2 - a2**2*b1**2 + a3**2*b2**2 - a4**2*b2**2)
f(664) = dtau*(a1**2*b1**3 + a2**2*b1**3 + a3**2*b2**3 + a4**2*b2**3)
f(665) = dtau*(a1**3*b1 + a2**3*b1 + a3**3*b2 + a4**3*b2)
f(666) = dtau*(a1**4*b1 + a2**4*b1 + a3**4*b2 + a4**4*b2)
f(667) = dtau*(a1*b2 + a2*b2 + a3*b1 + a4*b1)
f(668) = dtau**3*(a1*b2 + a2*b2 + a3*b1 + a4*b1)
f(669) = dtau**2*(-a1*b2**2 + a2*b2**2 - a3*b1**2 + a4*b1**2)
f(670) = dtau*(a1*b2**3 + a2*b2**3 + a3*b1**3 + a4*b1**3)
f(671) = dtau*(a1**2*b2 + a2**2*b2 + a3**2*b1 + a4**2*b1)
f(672) = dtau**3*(a1**2*b2 + a2**2*b2 + a3**2*b1 + a4**2*b1)
f(673) = dtau**2*(a1**2*b2**2 - a2**2*b2**2 + a3**2*b1**2 - a4**2*b1**2)
f(674) = dtau*(a1**2*b2**3 + a2**2*b2**3 + a3**2*b1**3 + a4**2*b1**3)
f(675) = dtau*(a1**3*b2 + a2**3*b2 + a3**3*b1 + a4**3*b1)
f(676) = dtau*(a1**4*b2 + a2**4*b2 + a3**4*b1 + a4**4*b1)
f(677) = b1*b2*dtau*(b1 + b2)
f(678) = b1*b2*dtau**3*(b1 + b2)
f(679) = b1*b2*dtau*(b1**3 + b2**3)
f(680) = b1**2*b2**2*dtau*(b1 + b2)
v = sum(f*params)
end function c2h4_dipole_b2u_n3_d6_ADF


!###############################################################################


function c2h4_dipole_b2u_n4_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(985)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(985)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(r1*r2*r3 - r1*r2*r4 + r1*r3*r4 - r2*r3*r4)
f(2) = r0*(-r1**2*r3*r4 - r1*r2*r3**2 + r1*r2*r4**2 + r2**2*r3*r4)
f(3) = r0*(-r1**3*r3*r4 - r1*r2*r3**3 + r1*r2*r4**3 + r2**3*r3*r4)
f(4) = r0*(r1**2*r2*r4 - r1*r2**2*r3 - r1*r3*r4**2 + r2*r3**2*r4)
f(5) = r0*(-r1**2*r2*r4**2 + r1**2*r3*r4**2 + r1*r2**2*r3**2 - r2**2*r3 &
      **2*r4)
f(6) = r0*(-r1**3*r2*r4 + r1*r2**3*r3 + r1*r3*r4**3 - r2*r3**3*r4)
f(7) = r0*(-r1**2*r2*r3 + r1*r2**2*r4 - r1*r3**2*r4 + r2*r3*r4**2)
f(8) = r0*(r1**2*r2*r3**2 + r1**2*r3**2*r4 - r1*r2**2*r4**2 - r2**2*r3* &
      r4**2)
f(9) = r0*(r1**2*r2**2*r3 - r1**2*r2**2*r4 + r1*r3**2*r4**2 - r2*r3**2* &
      r4**2)
f(10) = r0*(r1**3*r2*r3 - r1*r2**3*r4 + r1*r3**3*r4 - r2*r3*r4**3)
f(11) = r0**2*(-r1*r2*r3 + r1*r2*r4 - r1*r3*r4 + r2*r3*r4)
f(12) = r0**2*(r1**2*r3*r4 + r1*r2*r3**2 - r1*r2*r4**2 - r2**2*r3*r4)
f(13) = r0**2*(r1**2*r2*r4 - r1*r2**2*r3 - r1*r3*r4**2 + r2*r3**2*r4)
f(14) = r0**2*(-r1**2*r2*r3 + r1*r2**2*r4 - r1*r3**2*r4 + r2*r3*r4**2)
f(15) = r0**3*(r1*r2*r3 - r1*r2*r4 + r1*r3*r4 - r2*r3*r4)
f(16) = r0*(-a1*r1*r2 + a2*r1*r2 - a3*r3*r4 + a4*r3*r4)
f(17) = r0*(-a1**2*r1*r2 + a2**2*r1*r2 - a3**2*r3*r4 + a4**2*r3*r4)
f(18) = r0*(a1**3*r1*r2 - a2**3*r1*r2 + a3**3*r3*r4 - a4**3*r3*r4)
f(19) = r0*(-a1*r1*r2**2 + a2*r1**2*r2 - a3*r3*r4**2 + a4*r3**2*r4)
f(20) = r0*(-a1**2*r1*r2**2 + a2**2*r1**2*r2 - a3**2*r3*r4**2 + a4**2*r3 &
      **2*r4)
f(21) = r0*(a1*r1*r2**3 - a2*r1**3*r2 + a3*r3*r4**3 - a4*r3**3*r4)
f(22) = r0*(-a1*r1**2*r2 + a2*r1*r2**2 - a3*r3**2*r4 + a4*r3*r4**2)
f(23) = r0*(-a1**2*r1**2*r2 + a2**2*r1*r2**2 - a3**2*r3**2*r4 + a4**2*r3 &
      *r4**2)
f(24) = r0*(-a1*r1**2*r2**2 + a2*r1**2*r2**2 - a3*r3**2*r4**2 + a4*r3**2 &
      *r4**2)
f(25) = r0*(a1*r1**3*r2 - a2*r1*r2**3 + a3*r3**3*r4 - a4*r3*r4**3)
f(26) = r0**2*(-a1*r1*r2 + a2*r1*r2 - a3*r3*r4 + a4*r3*r4)
f(27) = r0**2*(a1**2*r1*r2 - a2**2*r1*r2 + a3**2*r3*r4 - a4**2*r3*r4)
f(28) = r0**2*(-a1*r1*r2**2 + a2*r1**2*r2 - a3*r3*r4**2 + a4*r3**2*r4)
f(29) = r0**2*(-a1*r1**2*r2 + a2*r1*r2**2 - a3*r3**2*r4 + a4*r3*r4**2)
f(30) = r0**3*(a1*r1*r2 - a2*r1*r2 + a3*r3*r4 - a4*r3*r4)
f(31) = r0*(-a1*r3*r4 + a2*r3*r4 - a3*r1*r2 + a4*r1*r2)
f(32) = r0*(a1**2*r3*r4 - a2**2*r3*r4 + a3**2*r1*r2 - a4**2*r1*r2)
f(33) = r0*(-a1**3*r3*r4 + a2**3*r3*r4 - a3**3*r1*r2 + a4**3*r1*r2)
f(34) = r0*(-a1*r3*r4**2 + a2*r3**2*r4 - a3*r1*r2**2 + a4*r1**2*r2)
f(35) = r0*(a1**2*r3*r4**2 - a2**2*r3**2*r4 + a3**2*r1*r2**2 - a4**2*r1 &
      **2*r2)
f(36) = r0*(a1*r3*r4**3 - a2*r3**3*r4 + a3*r1*r2**3 - a4*r1**3*r2)
f(37) = r0*(-a1*r3**2*r4 + a2*r3*r4**2 - a3*r1**2*r2 + a4*r1*r2**2)
f(38) = r0*(a1**2*r3**2*r4 - a2**2*r3*r4**2 + a3**2*r1**2*r2 - a4**2*r1* &
      r2**2)
f(39) = r0*(-a1*r3**2*r4**2 + a2*r3**2*r4**2 - a3*r1**2*r2**2 + a4*r1**2 &
      *r2**2)
f(40) = r0*(a1*r3**3*r4 - a2*r3*r4**3 + a3*r1**3*r2 - a4*r1*r2**3)
f(41) = r0**2*(-a1*r3*r4 + a2*r3*r4 - a3*r1*r2 + a4*r1*r2)
f(42) = r0**2*(a1**2*r3*r4 - a2**2*r3*r4 + a3**2*r1*r2 - a4**2*r1*r2)
f(43) = r0**2*(-a1*r3*r4**2 + a2*r3**2*r4 - a3*r1*r2**2 + a4*r1**2*r2)
f(44) = r0**2*(-a1*r3**2*r4 + a2*r3*r4**2 - a3*r1**2*r2 + a4*r1*r2**2)
f(45) = r0**3*(a1*r3*r4 - a2*r3*r4 + a3*r1*r2 - a4*r1*r2)
f(46) = r0*(-b1**2*r1**2*r2 + b1**2*r1*r2**2 - b2**2*r3**2*r4 + b2**2*r3 &
      *r4**2)
f(47) = r0*(b1**2*r1**2*r2 - b1**2*r1*r2**2 + b2**2*r3**2*r4 - b2**2*r3* &
      r4**2)
f(48) = r0*(b1**2*r3**2*r4 - b1**2*r3*r4**2 + b2**2*r1**2*r2 - b2**2*r1* &
      r2**2)
f(49) = r0*(-b1**2*r3**2*r4 + b1**2*r3*r4**2 - b2**2*r1**2*r2 + b2**2*r1 &
      *r2**2)
f(50) = dtau**2*r0*(r1**2*r2 - r1*r2**2 + r3**2*r4 - r3*r4**2)
f(51) = dtau**2*r0*(-r1**2*r2 + r1*r2**2 - r3**2*r4 + r3*r4**2)
f(52) = r0*(a1*r1*r3 - a2*r2*r4 + a3*r1*r3 - a4*r2*r4)
f(53) = r0*(a1**2*r1*r3 - a2**2*r2*r4 + a3**2*r1*r3 - a4**2*r2*r4)
f(54) = r0*(a1**3*r1*r3 - a2**3*r2*r4 + a3**3*r1*r3 - a4**3*r2*r4)
f(55) = r0*(a1*r1*r3**2 - a2*r2*r4**2 + a3*r1**2*r3 - a4*r2**2*r4)
f(56) = r0*(-a1**2*r1*r3**2 + a2**2*r2*r4**2 - a3**2*r1**2*r3 + a4**2*r2 &
      **2*r4)
f(57) = r0*(a1*r1*r3**3 - a2*r2*r4**3 + a3*r1**3*r3 - a4*r2**3*r4)
f(58) = r0*(a1*r1**2*r3 - a2*r2**2*r4 + a3*r1*r3**2 - a4*r2*r4**2)
f(59) = r0*(-a1**2*r1**2*r3 + a2**2*r2**2*r4 - a3**2*r1*r3**2 + a4**2*r2 &
      *r4**2)
f(60) = r0*(a1*r1**2*r3**2 - a2*r2**2*r4**2 + a3*r1**2*r3**2 - a4*r2**2* &
      r4**2)
f(61) = r0*(-a1*r1**3*r3 + a2*r2**3*r4 - a3*r1*r3**3 + a4*r2*r4**3)
f(62) = r0**2*(a1*r1*r3 - a2*r2*r4 + a3*r1*r3 - a4*r2*r4)
f(63) = r0**2*(a1**2*r1*r3 - a2**2*r2*r4 + a3**2*r1*r3 - a4**2*r2*r4)
f(64) = r0**2*(-a1*r1*r3**2 + a2*r2*r4**2 - a3*r1**2*r3 + a4*r2**2*r4)
f(65) = r0**2*(-a1*r1**2*r3 + a2*r2**2*r4 - a3*r1*r3**2 + a4*r2*r4**2)
f(66) = r0**3*(a1*r1*r3 - a2*r2*r4 + a3*r1*r3 - a4*r2*r4)
f(67) = r0*(a1*r2*r4 - a2*r1*r3 + a3*r2*r4 - a4*r1*r3)
f(68) = r0*(-a1**2*r2*r4 + a2**2*r1*r3 - a3**2*r2*r4 + a4**2*r1*r3)
f(69) = r0*(-a1**3*r2*r4 + a2**3*r1*r3 - a3**3*r2*r4 + a4**3*r1*r3)
f(70) = r0*(-a1*r2*r4**2 + a2*r1*r3**2 - a3*r2**2*r4 + a4*r1**2*r3)
f(71) = r0*(a1**2*r2*r4**2 - a2**2*r1*r3**2 + a3**2*r2**2*r4 - a4**2*r1 &
      **2*r3)
f(72) = r0*(-a1*r2*r4**3 + a2*r1*r3**3 - a3*r2**3*r4 + a4*r1**3*r3)
f(73) = r0*(-a1*r2**2*r4 + a2*r1**2*r3 - a3*r2*r4**2 + a4*r1*r3**2)
f(74) = r0*(a1**2*r2**2*r4 - a2**2*r1**2*r3 + a3**2*r2*r4**2 - a4**2*r1* &
      r3**2)
f(75) = r0*(a1*r2**2*r4**2 - a2*r1**2*r3**2 + a3*r2**2*r4**2 - a4*r1**2* &
      r3**2)
f(76) = r0*(a1*r2**3*r4 - a2*r1**3*r3 + a3*r2*r4**3 - a4*r1*r3**3)
f(77) = r0**2*(-a1*r2*r4 + a2*r1*r3 - a3*r2*r4 + a4*r1*r3)
f(78) = r0**2*(a1**2*r2*r4 - a2**2*r1*r3 + a3**2*r2*r4 - a4**2*r1*r3)
f(79) = r0**2*(a1*r2*r4**2 - a2*r1*r3**2 + a3*r2**2*r4 - a4*r1**2*r3)
f(80) = r0**2*(a1*r2**2*r4 - a2*r1**2*r3 + a3*r2*r4**2 - a4*r1*r3**2)
f(81) = r0**3*(-a1*r2*r4 + a2*r1*r3 - a3*r2*r4 + a4*r1*r3)
f(82) = r0*(-b1**2*r1*r3 + b1**2*r2*r4 - b2**2*r1*r3 + b2**2*r2*r4)
f(83) = r0*(-b1**2*r1*r3**2 + b1**2*r2*r4**2 - b2**2*r1**2*r3 + b2**2*r2 &
      **2*r4)
f(84) = r0*(-b1**2*r1**2*r3 + b1**2*r2**2*r4 - b2**2*r1*r3**2 + b2**2*r2 &
      *r4**2)
f(85) = r0**2*(-b1**2*r1*r3 + b1**2*r2*r4 - b2**2*r1*r3 + b2**2*r2*r4)
f(86) = dtau**2*r0*(r1*r3 - r2*r4)
f(87) = dtau**2*r0*(r1**2*r3 + r1*r3**2 - r2**2*r4 - r2*r4**2)
f(88) = dtau**2*r0**2*(r1*r3 - r2*r4)
f(89) = r0*(-a1*r1*r4 + a2*r2*r3 - a3*r2*r3 + a4*r1*r4)
f(90) = r0*(-a1**2*r1*r4 + a2**2*r2*r3 - a3**2*r2*r3 + a4**2*r1*r4)
f(91) = r0*(-a1**3*r1*r4 + a2**3*r2*r3 - a3**3*r2*r3 + a4**3*r1*r4)
f(92) = r0*(-a1*r1*r4**2 + a2*r2*r3**2 - a3*r2**2*r3 + a4*r1**2*r4)
f(93) = r0*(-a1**2*r1*r4**2 + a2**2*r2*r3**2 - a3**2*r2**2*r3 + a4**2*r1 &
      **2*r4)
f(94) = r0*(-a1*r1*r4**3 + a2*r2*r3**3 - a3*r2**3*r3 + a4*r1**3*r4)
f(95) = r0*(-a1*r1**2*r4 + a2*r2**2*r3 - a3*r2*r3**2 + a4*r1*r4**2)
f(96) = r0*(-a1**2*r1**2*r4 + a2**2*r2**2*r3 - a3**2*r2*r3**2 + a4**2*r1 &
      *r4**2)
f(97) = r0*(-a1*r1**2*r4**2 + a2*r2**2*r3**2 - a3*r2**2*r3**2 + a4*r1**2 &
      *r4**2)
f(98) = r0*(a1*r1**3*r4 - a2*r2**3*r3 + a3*r2*r3**3 - a4*r1*r4**3)
f(99) = r0**2*(-a1*r1*r4 + a2*r2*r3 - a3*r2*r3 + a4*r1*r4)
f(100) = r0**2*(-a1**2*r1*r4 + a2**2*r2*r3 - a3**2*r2*r3 + a4**2*r1*r4)
f(101) = r0**2*(-a1*r1*r4**2 + a2*r2*r3**2 - a3*r2**2*r3 + a4*r1**2*r4)
f(102) = r0**2*(-a1*r1**2*r4 + a2*r2**2*r3 - a3*r2*r3**2 + a4*r1*r4**2)
f(103) = r0**3*(a1*r1*r4 - a2*r2*r3 + a3*r2*r3 - a4*r1*r4)
f(104) = r0*(-a1*r2*r3 + a2*r1*r4 - a3*r1*r4 + a4*r2*r3)
f(105) = r0*(a1**2*r2*r3 - a2**2*r1*r4 + a3**2*r1*r4 - a4**2*r2*r3)
f(106) = r0*(a1**3*r2*r3 - a2**3*r1*r4 + a3**3*r1*r4 - a4**3*r2*r3)
f(107) = r0*(a1*r2*r3**2 - a2*r1*r4**2 + a3*r1**2*r4 - a4*r2**2*r3)
f(108) = r0*(a1**2*r2*r3**2 - a2**2*r1*r4**2 + a3**2*r1**2*r4 - a4**2*r2 &
      **2*r3)
f(109) = r0*(a1*r2*r3**3 - a2*r1*r4**3 + a3*r1**3*r4 - a4*r2**3*r3)
f(110) = r0*(a1*r2**2*r3 - a2*r1**2*r4 + a3*r1*r4**2 - a4*r2*r3**2)
f(111) = r0*(a1**2*r2**2*r3 - a2**2*r1**2*r4 + a3**2*r1*r4**2 - a4**2*r2 &
      *r3**2)
f(112) = r0*(-a1*r2**2*r3**2 + a2*r1**2*r4**2 - a3*r1**2*r4**2 + a4*r2** &
      2*r3**2)
f(113) = r0*(-a1*r2**3*r3 + a2*r1**3*r4 - a3*r1*r4**3 + a4*r2*r3**3)
f(114) = r0**2*(a1*r2*r3 - a2*r1*r4 + a3*r1*r4 - a4*r2*r3)
f(115) = r0**2*(a1**2*r2*r3 - a2**2*r1*r4 + a3**2*r1*r4 - a4**2*r2*r3)
f(116) = r0**2*(a1*r2*r3**2 - a2*r1*r4**2 + a3*r1**2*r4 - a4*r2**2*r3)
f(117) = r0**2*(a1*r2**2*r3 - a2*r1**2*r4 + a3*r1*r4**2 - a4*r2*r3**2)
f(118) = r0**3*(-a1*r2*r3 + a2*r1*r4 - a3*r1*r4 + a4*r2*r3)
f(119) = r0*(b1**2*r1*r4 - b1**2*r2*r3 - b2**2*r1*r4 + b2**2*r2*r3)
f(120) = r0*(b1**2*r1*r4**2 - b1**2*r2*r3**2 - b2**2*r1**2*r4 + b2**2*r2 &
      **2*r3)
f(121) = r0*(b1**2*r1**2*r4 - b1**2*r2**2*r3 - b2**2*r1*r4**2 + b2**2*r2 &
      *r3**2)
f(122) = r0**2*(b1**2*r1*r4 - b1**2*r2*r3 - b2**2*r1*r4 + b2**2*r2*r3)
f(123) = dtau**2*r0*(r1**2*r4 - r1*r4**2 - r2**2*r3 + r2*r3**2)
f(124) = dtau**2*r0*(-r1**2*r4 + r1*r4**2 + r2**2*r3 - r2*r3**2)
f(125) = r0*(-a1*a2*r1 + a1*a2*r2 - a3*a4*r3 + a3*a4*r4)
f(126) = r0*(-a1**2*a2*r2 + a1*a2**2*r1 - a3**2*a4*r4 + a3*a4**2*r3)
f(127) = r0*(a1**3*a2*r2 - a1*a2**3*r1 + a3**3*a4*r4 - a3*a4**3*r3)
f(128) = r0*(a1**2*a2*r1 - a1*a2**2*r2 + a3**2*a4*r3 - a3*a4**2*r4)
f(129) = r0*(a1**2*a2**2*r1 - a1**2*a2**2*r2 + a3**2*a4**2*r3 - a3**2*a4 &
      **2*r4)
f(130) = r0*(-a1**3*a2*r1 + a1*a2**3*r2 - a3**3*a4*r3 + a3*a4**3*r4)
f(131) = r0*(a1*a2*r1**2 - a1*a2*r2**2 + a3*a4*r3**2 - a3*a4*r4**2)
f(132) = r0*(a1**2*a2*r2**2 - a1*a2**2*r1**2 + a3**2*a4*r4**2 - a3*a4**2 &
      *r3**2)
f(133) = r0*(a1**2*a2*r1**2 - a1*a2**2*r2**2 + a3**2*a4*r3**2 - a3*a4**2 &
      *r4**2)
f(134) = r0*(a1*a2*r1**3 - a1*a2*r2**3 + a3*a4*r3**3 - a3*a4*r4**3)
f(135) = r0**2*(a1*a2*r1 - a1*a2*r2 + a3*a4*r3 - a3*a4*r4)
f(136) = r0**2*(a1**2*a2*r2 - a1*a2**2*r1 + a3**2*a4*r4 - a3*a4**2*r3)
f(137) = r0**2*(a1**2*a2*r1 - a1*a2**2*r2 + a3**2*a4*r3 - a3*a4**2*r4)
f(138) = r0**2*(a1*a2*r1**2 - a1*a2*r2**2 + a3*a4*r3**2 - a3*a4*r4**2)
f(139) = r0**3*(-a1*a2*r1 + a1*a2*r2 - a3*a4*r3 + a3*a4*r4)
f(140) = r0*(a1*a3*r1 + a1*a3*r3 - a2*a4*r2 - a2*a4*r4)
f(141) = r0*(-a1**2*a3*r3 - a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2*r2)
f(142) = r0*(a1**3*a3*r3 + a1*a3**3*r1 - a2**3*a4*r4 - a2*a4**3*r2)
f(143) = r0*(-a1**2*a3*r1 - a1*a3**2*r3 + a2**2*a4*r2 + a2*a4**2*r4)
f(144) = r0*(a1**2*a3**2*r1 + a1**2*a3**2*r3 - a2**2*a4**2*r2 - a2**2*a4 &
      **2*r4)
f(145) = r0*(a1**3*a3*r1 + a1*a3**3*r3 - a2**3*a4*r2 - a2*a4**3*r4)
f(146) = r0*(a1*a3*r1**2 + a1*a3*r3**2 - a2*a4*r2**2 - a2*a4*r4**2)
f(147) = r0*(-a1**2*a3*r3**2 - a1*a3**2*r1**2 + a2**2*a4*r4**2 + a2*a4** &
      2*r2**2)
f(148) = r0*(a1**2*a3*r1**2 + a1*a3**2*r3**2 - a2**2*a4*r2**2 - a2*a4**2 &
      *r4**2)
f(149) = r0*(a1*a3*r1**3 + a1*a3*r3**3 - a2*a4*r2**3 - a2*a4*r4**3)
f(150) = r0**2*(a1*a3*r1 + a1*a3*r3 - a2*a4*r2 - a2*a4*r4)
f(151) = r0**2*(-a1**2*a3*r3 - a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2*r2)
f(152) = r0**2*(a1**2*a3*r1 + a1*a3**2*r3 - a2**2*a4*r2 - a2*a4**2*r4)
f(153) = r0**2*(a1*a3*r1**2 + a1*a3*r3**2 - a2*a4*r2**2 - a2*a4*r4**2)
f(154) = r0**3*(-a1*a3*r1 - a1*a3*r3 + a2*a4*r2 + a2*a4*r4)
f(155) = r0*(a1*a4*r1 - a1*a4*r4 - a2*a3*r2 + a2*a3*r3)
f(156) = r0*(-a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 - a2*a3**2*r2)
f(157) = r0*(a1**3*a4*r4 - a1*a4**3*r1 - a2**3*a3*r3 + a2*a3**3*r2)
f(158) = r0*(a1**2*a4*r1 - a1*a4**2*r4 - a2**2*a3*r2 + a2*a3**2*r3)
f(159) = r0*(a1**2*a4**2*r1 - a1**2*a4**2*r4 - a2**2*a3**2*r2 + a2**2*a3 &
      **2*r3)
f(160) = r0*(-a1**3*a4*r1 + a1*a4**3*r4 + a2**3*a3*r2 - a2*a3**3*r3)
f(161) = r0*(a1*a4*r1**2 - a1*a4*r4**2 - a2*a3*r2**2 + a2*a3*r3**2)
f(162) = r0*(-a1**2*a4*r4**2 + a1*a4**2*r1**2 + a2**2*a3*r3**2 - a2*a3** &
      2*r2**2)
f(163) = r0*(a1**2*a4*r1**2 - a1*a4**2*r4**2 - a2**2*a3*r2**2 + a2*a3**2 &
      *r3**2)
f(164) = r0*(a1*a4*r1**3 - a1*a4*r4**3 - a2*a3*r2**3 + a2*a3*r3**3)
f(165) = r0**2*(a1*a4*r1 - a1*a4*r4 - a2*a3*r2 + a2*a3*r3)
f(166) = r0**2*(-a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 - a2*a3**2*r2)
f(167) = r0**2*(a1**2*a4*r1 - a1*a4**2*r4 - a2**2*a3*r2 + a2*a3**2*r3)
f(168) = r0**2*(a1*a4*r1**2 - a1*a4*r4**2 - a2*a3*r2**2 + a2*a3*r3**2)
f(169) = r0**3*(-a1*a4*r1 + a1*a4*r4 + a2*a3*r2 - a2*a3*r3)
f(170) = r0*(a1*b1**2*r1 - a2*b1**2*r2 + a3*b2**2*r3 - a4*b2**2*r4)
f(171) = r0*(-a1**2*b1**2*r1 + a2**2*b1**2*r2 - a3**2*b2**2*r3 + a4**2* &
      b2**2*r4)
f(172) = r0*(-a1*b1**2*r1**2 + a2*b1**2*r2**2 - a3*b2**2*r3**2 + a4*b2** &
      2*r4**2)
f(173) = r0**2*(-a1*b1**2*r1 + a2*b1**2*r2 - a3*b2**2*r3 + a4*b2**2*r4)
f(174) = r0*(a1*b2**2*r1 - a2*b2**2*r2 + a3*b1**2*r3 - a4*b1**2*r4)
f(175) = r0*(-a1**2*b2**2*r1 + a2**2*b2**2*r2 - a3**2*b1**2*r3 + a4**2* &
      b1**2*r4)
f(176) = r0*(-a1*b2**2*r1**2 + a2*b2**2*r2**2 - a3*b1**2*r3**2 + a4*b1** &
      2*r4**2)
f(177) = r0**2*(-a1*b2**2*r1 + a2*b2**2*r2 - a3*b1**2*r3 + a4*b1**2*r4)
f(178) = dtau**2*r0*(-a1*r1 + a2*r2 - a3*r3 + a4*r4)
f(179) = dtau**2*r0*(a1**2*r1 - a2**2*r2 + a3**2*r3 - a4**2*r4)
f(180) = dtau**2*r0*(-a1*r1**2 + a2*r2**2 - a3*r3**2 + a4*r4**2)
f(181) = dtau**2*r0**2*(a1*r1 - a2*r2 + a3*r3 - a4*r4)
f(182) = r0*(-a1*a4*r2 + a1*a4*r3 + a2*a3*r1 - a2*a3*r4)
f(183) = r0*(-a1**2*a4*r3 + a1*a4**2*r2 + a2**2*a3*r4 - a2*a3**2*r1)
f(184) = r0*(a1**3*a4*r3 - a1*a4**3*r2 - a2**3*a3*r4 + a2*a3**3*r1)
f(185) = r0*(-a1**2*a4*r2 + a1*a4**2*r3 + a2**2*a3*r1 - a2*a3**2*r4)
f(186) = r0*(-a1**2*a4**2*r2 + a1**2*a4**2*r3 + a2**2*a3**2*r1 - a2**2* &
      a3**2*r4)
f(187) = r0*(a1**3*a4*r2 - a1*a4**3*r3 - a2**3*a3*r1 + a2*a3**3*r4)
f(188) = r0*(-a1*a4*r2**2 + a1*a4*r3**2 + a2*a3*r1**2 - a2*a3*r4**2)
f(189) = r0*(a1**2*a4*r3**2 - a1*a4**2*r2**2 - a2**2*a3*r4**2 + a2*a3**2 &
      *r1**2)
f(190) = r0*(-a1**2*a4*r2**2 + a1*a4**2*r3**2 + a2**2*a3*r1**2 - a2*a3** &
      2*r4**2)
f(191) = r0*(a1*a4*r2**3 - a1*a4*r3**3 - a2*a3*r1**3 + a2*a3*r4**3)
f(192) = r0**2*(-a1*a4*r2 + a1*a4*r3 + a2*a3*r1 - a2*a3*r4)
f(193) = r0**2*(a1**2*a4*r3 - a1*a4**2*r2 - a2**2*a3*r4 + a2*a3**2*r1)
f(194) = r0**2*(-a1**2*a4*r2 + a1*a4**2*r3 + a2**2*a3*r1 - a2*a3**2*r4)
f(195) = r0**2*(-a1*a4*r2**2 + a1*a4*r3**2 + a2*a3*r1**2 - a2*a3*r4**2)
f(196) = r0**3*(a1*a4*r2 - a1*a4*r3 - a2*a3*r1 + a2*a3*r4)
f(197) = r0*(-a1*a3*r2 - a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(198) = r0*(-a1**2*a3*r4 - a1*a3**2*r2 + a2**2*a4*r3 + a2*a4**2*r1)
f(199) = r0*(a1**3*a3*r4 + a1*a3**3*r2 - a2**3*a4*r3 - a2*a4**3*r1)
f(200) = r0*(a1**2*a3*r2 + a1*a3**2*r4 - a2**2*a4*r1 - a2*a4**2*r3)
f(201) = r0*(-a1**2*a3**2*r2 - a1**2*a3**2*r4 + a2**2*a4**2*r1 + a2**2* &
      a4**2*r3)
f(202) = r0*(-a1**3*a3*r2 - a1*a3**3*r4 + a2**3*a4*r1 + a2*a4**3*r3)
f(203) = r0*(-a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3**2)
f(204) = r0*(a1**2*a3*r4**2 + a1*a3**2*r2**2 - a2**2*a4*r3**2 - a2*a4**2 &
      *r1**2)
f(205) = r0*(-a1**2*a3*r2**2 - a1*a3**2*r4**2 + a2**2*a4*r1**2 + a2*a4** &
      2*r3**2)
f(206) = r0*(-a1*a3*r2**3 - a1*a3*r4**3 + a2*a4*r1**3 + a2*a4*r3**3)
f(207) = r0**2*(-a1*a3*r2 - a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(208) = r0**2*(a1**2*a3*r4 + a1*a3**2*r2 - a2**2*a4*r3 - a2*a4**2*r1)
f(209) = r0**2*(-a1**2*a3*r2 - a1*a3**2*r4 + a2**2*a4*r1 + a2*a4**2*r3)
f(210) = r0**2*(-a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3**2)
f(211) = r0**3*(-a1*a3*r2 - a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(212) = r0*(-a1*b1**2*r2 + a2*b1**2*r1 - a3*b2**2*r4 + a4*b2**2*r3)
f(213) = r0*(-a1**2*b1**2*r2 + a2**2*b1**2*r1 - a3**2*b2**2*r4 + a4**2* &
      b2**2*r3)
f(214) = r0*(a1*b1**2*r2**2 - a2*b1**2*r1**2 + a3*b2**2*r4**2 - a4*b2**2 &
      *r3**2)
f(215) = r0**2*(a1*b1**2*r2 - a2*b1**2*r1 + a3*b2**2*r4 - a4*b2**2*r3)
f(216) = r0*(a1*b2**2*r2 - a2*b2**2*r1 + a3*b1**2*r4 - a4*b1**2*r3)
f(217) = r0*(a1**2*b2**2*r2 - a2**2*b2**2*r1 + a3**2*b1**2*r4 - a4**2*b1 &
      **2*r3)
f(218) = r0*(a1*b2**2*r2**2 - a2*b2**2*r1**2 + a3*b1**2*r4**2 - a4*b1**2 &
      *r3**2)
f(219) = r0**2*(a1*b2**2*r2 - a2*b2**2*r1 + a3*b1**2*r4 - a4*b1**2*r3)
f(220) = dtau**2*r0*(a1*r2 - a2*r1 + a3*r4 - a4*r3)
f(221) = dtau**2*r0*(a1**2*r2 - a2**2*r1 + a3**2*r4 - a4**2*r3)
f(222) = dtau**2*r0*(-a1*r2**2 + a2*r1**2 - a3*r4**2 + a4*r3**2)
f(223) = dtau**2*r0**2*(-a1*r2 + a2*r1 - a3*r4 + a4*r3)
f(224) = r0*(a1*a2*r3 - a1*a2*r4 + a3*a4*r1 - a3*a4*r2)
f(225) = r0*(a1**2*a2*r4 - a1*a2**2*r3 + a3**2*a4*r2 - a3*a4**2*r1)
f(226) = r0*(a1**3*a2*r4 - a1*a2**3*r3 + a3**3*a4*r2 - a3*a4**3*r1)
f(227) = r0*(-a1**2*a2*r3 + a1*a2**2*r4 - a3**2*a4*r1 + a3*a4**2*r2)
f(228) = r0*(-a1**2*a2**2*r3 + a1**2*a2**2*r4 - a3**2*a4**2*r1 + a3**2* &
      a4**2*r2)
f(229) = r0*(-a1**3*a2*r3 + a1*a2**3*r4 - a3**3*a4*r1 + a3*a4**3*r2)
f(230) = r0*(-a1*a2*r3**2 + a1*a2*r4**2 - a3*a4*r1**2 + a3*a4*r2**2)
f(231) = r0*(-a1**2*a2*r4**2 + a1*a2**2*r3**2 - a3**2*a4*r2**2 + a3*a4** &
      2*r1**2)
f(232) = r0*(a1**2*a2*r3**2 - a1*a2**2*r4**2 + a3**2*a4*r1**2 - a3*a4**2 &
      *r2**2)
f(233) = r0*(-a1*a2*r3**3 + a1*a2*r4**3 - a3*a4*r1**3 + a3*a4*r2**3)
f(234) = r0**2*(-a1*a2*r3 + a1*a2*r4 - a3*a4*r1 + a3*a4*r2)
f(235) = r0**2*(-a1**2*a2*r4 + a1*a2**2*r3 - a3**2*a4*r2 + a3*a4**2*r1)
f(236) = r0**2*(a1**2*a2*r3 - a1*a2**2*r4 + a3**2*a4*r1 - a3*a4**2*r2)
f(237) = r0**2*(-a1*a2*r3**2 + a1*a2*r4**2 - a3*a4*r1**2 + a3*a4*r2**2)
f(238) = r0**3*(-a1*a2*r3 + a1*a2*r4 - a3*a4*r1 + a3*a4*r2)
f(239) = r0*(-a1*b2**2*r3 + a2*b2**2*r4 - a3*b1**2*r1 + a4*b1**2*r2)
f(240) = r0*(a1**2*b2**2*r3 - a2**2*b2**2*r4 + a3**2*b1**2*r1 - a4**2*b1 &
      **2*r2)
f(241) = r0*(a1*b2**2*r3**2 - a2*b2**2*r4**2 + a3*b1**2*r1**2 - a4*b1**2 &
      *r2**2)
f(242) = r0**2*(a1*b2**2*r3 - a2*b2**2*r4 + a3*b1**2*r1 - a4*b1**2*r2)
f(243) = r0*(a1*b1**2*r3 - a2*b1**2*r4 + a3*b2**2*r1 - a4*b2**2*r2)
f(244) = r0*(-a1**2*b1**2*r3 + a2**2*b1**2*r4 - a3**2*b2**2*r1 + a4**2* &
      b2**2*r2)
f(245) = r0*(a1*b1**2*r3**2 - a2*b1**2*r4**2 + a3*b2**2*r1**2 - a4*b2**2 &
      *r2**2)
f(246) = r0**2*(a1*b1**2*r3 - a2*b1**2*r4 + a3*b2**2*r1 - a4*b2**2*r2)
f(247) = dtau**2*r0*(-a1*r3 + a2*r4 - a3*r1 + a4*r2)
f(248) = dtau**2*r0*(a1**2*r3 - a2**2*r4 + a3**2*r1 - a4**2*r2)
f(249) = dtau**2*r0*(-a1*r3**2 + a2*r4**2 - a3*r1**2 + a4*r2**2)
f(250) = dtau**2*r0**2*(-a1*r3 + a2*r4 - a3*r1 + a4*r2)
f(251) = r0*(-a1*b2**2*r4 + a2*b2**2*r3 - a3*b1**2*r2 + a4*b1**2*r1)
f(252) = r0*(a1**2*b2**2*r4 - a2**2*b2**2*r3 + a3**2*b1**2*r2 - a4**2*b1 &
      **2*r1)
f(253) = r0*(a1*b2**2*r4**2 - a2*b2**2*r3**2 + a3*b1**2*r2**2 - a4*b1**2 &
      *r1**2)
f(254) = r0**2*(a1*b2**2*r4 - a2*b2**2*r3 + a3*b1**2*r2 - a4*b1**2*r1)
f(255) = r0*(-a1*b1**2*r4 + a2*b1**2*r3 - a3*b2**2*r2 + a4*b2**2*r1)
f(256) = r0*(-a1**2*b1**2*r4 + a2**2*b1**2*r3 - a3**2*b2**2*r2 + a4**2* &
      b2**2*r1)
f(257) = r0*(-a1*b1**2*r4**2 + a2*b1**2*r3**2 - a3*b2**2*r2**2 + a4*b2** &
      2*r1**2)
f(258) = r0**2*(-a1*b1**2*r4 + a2*b1**2*r3 - a3*b2**2*r2 + a4*b2**2*r1)
f(259) = dtau**2*r0*(-a1*r4 + a2*r3 - a3*r2 + a4*r1)
f(260) = dtau**2*r0*(-a1**2*r4 + a2**2*r3 - a3**2*r2 + a4**2*r1)
f(261) = dtau**2*r0*(a1*r4**2 - a2*r3**2 + a3*r2**2 - a4*r1**2)
f(262) = dtau**2*r0**2*(a1*r4 - a2*r3 + a3*r2 - a4*r1)
f(263) = b1*b2*r0*(r1 - r2 + r3 - r4)
f(264) = b1*b2*r0*(-b1**2*r3 + b1**2*r4 - b2**2*r1 + b2**2*r2)
f(265) = b1**2*b2**2*r0*(-r1 + r2 - r3 + r4)
f(266) = b1*b2*r0*(b1**2*r1 - b1**2*r2 + b2**2*r3 - b2**2*r4)
f(267) = b1*b2*r0*(r1**2 - r2**2 + r3**2 - r4**2)
f(268) = b1*b2*r0*(r1**3 - r2**3 + r3**3 - r4**3)
f(269) = b1*b2*r0**2*(r1 - r2 + r3 - r4)
f(270) = b1*b2*r0**2*(r1**2 - r2**2 + r3**2 - r4**2)
f(271) = b1*b2*r0**3*(-r1 + r2 - r3 + r4)
f(272) = dtau*r0*(b1*r1 + b1*r2 + b2*r3 + b2*r4)
f(273) = dtau**3*r0*(b1*r1 + b1*r2 + b2*r3 + b2*r4)
f(274) = dtau**2*r0*(b1**2*r1 - b1**2*r2 + b2**2*r3 - b2**2*r4)
f(275) = dtau*r0*(b1**3*r1 + b1**3*r2 + b2**3*r3 + b2**3*r4)
f(276) = dtau*r0*(b1*r1**2 + b1*r2**2 + b2*r3**2 + b2*r4**2)
f(277) = dtau*r0*(b1*r1**3 + b1*r2**3 + b2*r3**3 + b2*r4**3)
f(278) = dtau*r0**2*(b1*r1 + b1*r2 + b2*r3 + b2*r4)
f(279) = dtau*r0**2*(b1*r1**2 + b1*r2**2 + b2*r3**2 + b2*r4**2)
f(280) = dtau*r0**3*(b1*r1 + b1*r2 + b2*r3 + b2*r4)
f(281) = dtau*r0*(b1*r3 + b1*r4 + b2*r1 + b2*r2)
f(282) = dtau**3*r0*(b1*r3 + b1*r4 + b2*r1 + b2*r2)
f(283) = dtau**2*r0*(b1**2*r3 - b1**2*r4 + b2**2*r1 - b2**2*r2)
f(284) = dtau*r0*(b1**3*r3 + b1**3*r4 + b2**3*r1 + b2**3*r2)
f(285) = dtau*r0*(b1*r3**2 + b1*r4**2 + b2*r1**2 + b2*r2**2)
f(286) = dtau*r0*(b1*r3**3 + b1*r4**3 + b2*r1**3 + b2*r2**3)
f(287) = dtau*r0**2*(b1*r3 + b1*r4 + b2*r1 + b2*r2)
f(288) = dtau*r0**2*(b1*r3**2 + b1*r4**2 + b2*r1**2 + b2*r2**2)
f(289) = dtau*r0**3*(b1*r3 + b1*r4 + b2*r1 + b2*r2)
f(290) = r0*(-a1*a2*a3 + a1*a2*a4 - a1*a3*a4 + a2*a3*a4)
f(291) = r0*(a1**2*a3*a4 + a1*a2*a3**2 - a1*a2*a4**2 - a2**2*a3*a4)
f(292) = r0*(-a1**3*a3*a4 - a1*a2*a3**3 + a1*a2*a4**3 + a2**3*a3*a4)
f(293) = r0*(-a1**2*a2*a4 + a1*a2**2*a3 + a1*a3*a4**2 - a2*a3**2*a4)
f(294) = r0*(a1**2*a2*a4**2 - a1**2*a3*a4**2 - a1*a2**2*a3**2 + a2**2*a3 &
      **2*a4)
f(295) = r0*(-a1**3*a2*a4 + a1*a2**3*a3 + a1*a3*a4**3 - a2*a3**3*a4)
f(296) = r0*(-a1**2*a2*a3 + a1*a2**2*a4 - a1*a3**2*a4 + a2*a3*a4**2)
f(297) = r0*(-a1**2*a2*a3**2 - a1**2*a3**2*a4 + a1*a2**2*a4**2 + a2**2* &
      a3*a4**2)
f(298) = r0*(-a1**2*a2**2*a3 + a1**2*a2**2*a4 - a1*a3**2*a4**2 + a2*a3** &
      2*a4**2)
f(299) = r0*(a1**3*a2*a3 - a1*a2**3*a4 + a1*a3**3*a4 - a2*a3*a4**3)
f(300) = r0**2*(a1*a2*a3 - a1*a2*a4 + a1*a3*a4 - a2*a3*a4)
f(301) = r0**2*(a1**2*a3*a4 + a1*a2*a3**2 - a1*a2*a4**2 - a2**2*a3*a4)
f(302) = r0**2*(a1**2*a2*a4 - a1*a2**2*a3 - a1*a3*a4**2 + a2*a3**2*a4)
f(303) = r0**2*(a1**2*a2*a3 - a1*a2**2*a4 + a1*a3**2*a4 - a2*a3*a4**2)
f(304) = r0**3*(a1*a2*a3 - a1*a2*a4 + a1*a3*a4 - a2*a3*a4)
f(305) = r0*(a1**2*a2*b1**2 - a1*a2**2*b1**2 + a3**2*a4*b2**2 - a3*a4**2 &
      *b2**2)
f(306) = r0*(-a1**2*a2*b2**2 + a1*a2**2*b2**2 - a3**2*a4*b1**2 + a3*a4** &
      2*b1**2)
f(307) = dtau**2*r0*(-a1**2*a2 + a1*a2**2 - a3**2*a4 + a3*a4**2)
f(308) = r0*(-a1*a3*b1**2 - a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2**2)
f(309) = r0*(a1**2*a3*b2**2 + a1*a3**2*b1**2 - a2**2*a4*b2**2 - a2*a4**2 &
      *b1**2)
f(310) = r0*(-a1**2*a3*b1**2 - a1*a3**2*b2**2 + a2**2*a4*b1**2 + a2*a4** &
      2*b2**2)
f(311) = r0**2*(-a1*a3*b1**2 - a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2**2)
f(312) = dtau**2*r0*(a1*a3 - a2*a4)
f(313) = dtau**2*r0*(-a1**2*a3 - a1*a3**2 + a2**2*a4 + a2*a4**2)
f(314) = dtau**2*r0**2*(-a1*a3 + a2*a4)
f(315) = r0*(a1*a4*b1**2 - a1*a4*b2**2 - a2*a3*b1**2 + a2*a3*b2**2)
f(316) = r0*(a1**2*a4*b2**2 - a1*a4**2*b1**2 - a2**2*a3*b2**2 + a2*a3**2 &
      *b1**2)
f(317) = r0*(a1**2*a4*b1**2 - a1*a4**2*b2**2 - a2**2*a3*b1**2 + a2*a3**2 &
      *b2**2)
f(318) = r0**2*(-a1*a4*b1**2 + a1*a4*b2**2 + a2*a3*b1**2 - a2*a3*b2**2)
f(319) = dtau**2*r0*(-a1**2*a4 + a1*a4**2 + a2**2*a3 - a2*a3**2)
f(320) = b1*b2*r0*(-a1 + a2 - a3 + a4)
f(321) = b1*b2*r0*(a1*b2**2 - a2*b2**2 + a3*b1**2 - a4*b1**2)
f(322) = b1**2*b2**2*r0*(-a1 + a2 - a3 + a4)
f(323) = b1*b2*r0*(a1*b1**2 - a2*b1**2 + a3*b2**2 - a4*b2**2)
f(324) = b1*b2*r0*(-a1**2 + a2**2 - a3**2 + a4**2)
f(325) = b1*b2*r0*(a1**3 - a2**3 + a3**3 - a4**3)
f(326) = b1*b2*r0**2*(-a1 + a2 - a3 + a4)
f(327) = b1*b2*r0**2*(a1**2 - a2**2 + a3**2 - a4**2)
f(328) = b1*b2*r0**3*(a1 - a2 + a3 - a4)
f(329) = dtau*r0*(a1*b1 + a2*b1 + a3*b2 + a4*b2)
f(330) = dtau**3*r0*(a1*b1 + a2*b1 + a3*b2 + a4*b2)
f(331) = dtau**2*r0*(a1*b1**2 - a2*b1**2 + a3*b2**2 - a4*b2**2)
f(332) = dtau*r0*(a1*b1**3 + a2*b1**3 + a3*b2**3 + a4*b2**3)
f(333) = dtau*r0*(a1**2*b1 + a2**2*b1 + a3**2*b2 + a4**2*b2)
f(334) = dtau*r0*(a1**3*b1 + a2**3*b1 + a3**3*b2 + a4**3*b2)
f(335) = dtau*r0**2*(a1*b1 + a2*b1 + a3*b2 + a4*b2)
f(336) = dtau*r0**2*(a1**2*b1 + a2**2*b1 + a3**2*b2 + a4**2*b2)
f(337) = dtau*r0**3*(a1*b1 + a2*b1 + a3*b2 + a4*b2)
f(338) = dtau*r0*(a1*b2 + a2*b2 + a3*b1 + a4*b1)
f(339) = dtau**3*r0*(a1*b2 + a2*b2 + a3*b1 + a4*b1)
f(340) = dtau**2*r0*(-a1*b2**2 + a2*b2**2 - a3*b1**2 + a4*b1**2)
f(341) = dtau*r0*(a1*b2**3 + a2*b2**3 + a3*b1**3 + a4*b1**3)
f(342) = dtau*r0*(a1**2*b2 + a2**2*b2 + a3**2*b1 + a4**2*b1)
f(343) = dtau*r0*(a1**3*b2 + a2**3*b2 + a3**3*b1 + a4**3*b1)
f(344) = dtau*r0**2*(a1*b2 + a2*b2 + a3*b1 + a4*b1)
f(345) = dtau*r0**2*(a1**2*b2 + a2**2*b2 + a3**2*b1 + a4**2*b1)
f(346) = dtau*r0**3*(a1*b2 + a2*b2 + a3*b1 + a4*b1)
f(347) = b1*b2*dtau*r0*(b1 + b2)
f(348) = b1*b2*dtau*r0**2*(b1 + b2)
f(349) = r1*r2*r3*r4*(r1 - r2 + r3 - r4)
f(350) = r1*r2*r3*r4*(-r1**2 + r2**2 - r3**2 + r4**2)
f(351) = r1*r2*r3*r4*(-r1*r3 + r2*r4)
f(352) = -a1*r1*r2*r3 + a2*r1*r2*r4 - a3*r1*r3*r4 + a4*r2*r3*r4
f(353) = -a1**2*r1*r2*r3 + a2**2*r1*r2*r4 - a3**2*r1*r3*r4 + a4**2*r2*r3 &
      *r4
f(354) = -a1**3*r1*r2*r3 + a2**3*r1*r2*r4 - a3**3*r1*r3*r4 + a4**3*r2*r3 &
      *r4
f(355) = a1*r1*r2*r3**2 - a2*r1*r2*r4**2 + a3*r1**2*r3*r4 - a4*r2**2*r3* &
      r4
f(356) = -a1**2*r1*r2*r3**2 + a2**2*r1*r2*r4**2 - a3**2*r1**2*r3*r4 + a4 &
      **2*r2**2*r3*r4
f(357) = -a1*r1*r2*r3**3 + a2*r1*r2*r4**3 - a3*r1**3*r3*r4 + a4*r2**3*r3 &
      *r4
f(358) = a1*r1*r2**2*r3 - a2*r1**2*r2*r4 + a3*r1*r3*r4**2 - a4*r2*r3**2* &
      r4
f(359) = a1**2*r1*r2**2*r3 - a2**2*r1**2*r2*r4 + a3**2*r1*r3*r4**2 - a4 &
      **2*r2*r3**2*r4
f(360) = -a1*r1*r2**2*r3**2 + a2*r1**2*r2*r4**2 - a3*r1**2*r3*r4**2 + a4 &
      *r2**2*r3**2*r4
f(361) = -a1*r1*r2**3*r3 + a2*r1**3*r2*r4 - a3*r1*r3*r4**3 + a4*r2*r3**3 &
      *r4
f(362) = a1*r1**2*r2*r3 - a2*r1*r2**2*r4 + a3*r1*r3**2*r4 - a4*r2*r3*r4 &
      **2
f(363) = -a1**2*r1**2*r2*r3 + a2**2*r1*r2**2*r4 - a3**2*r1*r3**2*r4 + a4 &
      **2*r2*r3*r4**2
f(364) = -a1*r1**2*r2*r3**2 + a2*r1*r2**2*r4**2 - a3*r1**2*r3**2*r4 + a4 &
      *r2**2*r3*r4**2
f(365) = -a1*r1**2*r2**2*r3 + a2*r1**2*r2**2*r4 - a3*r1*r3**2*r4**2 + a4 &
      *r2*r3**2*r4**2
f(366) = -a1*r1**3*r2*r3 + a2*r1*r2**3*r4 - a3*r1*r3**3*r4 + a4*r2*r3*r4 &
      **3
f(367) = a1*r1*r2*r4 - a2*r1*r2*r3 + a3*r2*r3*r4 - a4*r1*r3*r4
f(368) = a1**2*r1*r2*r4 - a2**2*r1*r2*r3 + a3**2*r2*r3*r4 - a4**2*r1*r3* &
      r4
f(369) = a1**3*r1*r2*r4 - a2**3*r1*r2*r3 + a3**3*r2*r3*r4 - a4**3*r1*r3* &
      r4
f(370) = -a1*r1*r2*r4**2 + a2*r1*r2*r3**2 - a3*r2**2*r3*r4 + a4*r1**2*r3 &
      *r4
f(371) = a1**2*r1*r2*r4**2 - a2**2*r1*r2*r3**2 + a3**2*r2**2*r3*r4 - a4 &
      **2*r1**2*r3*r4
f(372) = a1*r1*r2*r4**3 - a2*r1*r2*r3**3 + a3*r2**3*r3*r4 - a4*r1**3*r3* &
      r4
f(373) = -a1*r1**2*r2*r4 + a2*r1*r2**2*r3 - a3*r2*r3**2*r4 + a4*r1*r3*r4 &
      **2
f(374) = a1**2*r1**2*r2*r4 - a2**2*r1*r2**2*r3 + a3**2*r2*r3**2*r4 - a4 &
      **2*r1*r3*r4**2
f(375) = a1*r1**2*r2*r4**2 - a2*r1*r2**2*r3**2 + a3*r2**2*r3**2*r4 - a4* &
      r1**2*r3*r4**2
f(376) = a1*r1**3*r2*r4 - a2*r1*r2**3*r3 + a3*r2*r3**3*r4 - a4*r1*r3*r4 &
      **3
f(377) = -a1*r1*r2**2*r4 + a2*r1**2*r2*r3 - a3*r2*r3*r4**2 + a4*r1*r3**2 &
      *r4
f(378) = a1**2*r1*r2**2*r4 - a2**2*r1**2*r2*r3 + a3**2*r2*r3*r4**2 - a4 &
      **2*r1*r3**2*r4
f(379) = a1*r1*r2**2*r4**2 - a2*r1**2*r2*r3**2 + a3*r2**2*r3*r4**2 - a4* &
      r1**2*r3**2*r4
f(380) = -a1*r1**2*r2**2*r4 + a2*r1**2*r2**2*r3 - a3*r2*r3**2*r4**2 + a4 &
      *r1*r3**2*r4**2
f(381) = a1*r1*r2**3*r4 - a2*r1**3*r2*r3 + a3*r2*r3*r4**3 - a4*r1*r3**3* &
      r4
f(382) = a1*r1*r3*r4 - a2*r2*r3*r4 + a3*r1*r2*r3 - a4*r1*r2*r4
f(383) = a1**2*r1*r3*r4 - a2**2*r2*r3*r4 + a3**2*r1*r2*r3 - a4**2*r1*r2* &
      r4
f(384) = a1**3*r1*r3*r4 - a2**3*r2*r3*r4 + a3**3*r1*r2*r3 - a4**3*r1*r2* &
      r4
f(385) = -a1*r1**2*r3*r4 + a2*r2**2*r3*r4 - a3*r1*r2*r3**2 + a4*r1*r2*r4 &
      **2
f(386) = a1**2*r1**2*r3*r4 - a2**2*r2**2*r3*r4 + a3**2*r1*r2*r3**2 - a4 &
      **2*r1*r2*r4**2
f(387) = a1*r1**3*r3*r4 - a2*r2**3*r3*r4 + a3*r1*r2*r3**3 - a4*r1*r2*r4 &
      **3
f(388) = -a1*r1*r3*r4**2 + a2*r2*r3**2*r4 - a3*r1*r2**2*r3 + a4*r1**2*r2 &
      *r4
f(389) = -a1**2*r1*r3*r4**2 + a2**2*r2*r3**2*r4 - a3**2*r1*r2**2*r3 + a4 &
      **2*r1**2*r2*r4
f(390) = a1*r1**2*r3*r4**2 - a2*r2**2*r3**2*r4 + a3*r1*r2**2*r3**2 - a4* &
      r1**2*r2*r4**2
f(391) = a1*r1*r3*r4**3 - a2*r2*r3**3*r4 + a3*r1*r2**3*r3 - a4*r1**3*r2* &
      r4
f(392) = -a1*r1*r3**2*r4 + a2*r2*r3*r4**2 - a3*r1**2*r2*r3 + a4*r1*r2**2 &
      *r4
f(393) = -a1**2*r1*r3**2*r4 + a2**2*r2*r3*r4**2 - a3**2*r1**2*r2*r3 + a4 &
      **2*r1*r2**2*r4
f(394) = a1*r1**2*r3**2*r4 - a2*r2**2*r3*r4**2 + a3*r1**2*r2*r3**2 - a4* &
      r1*r2**2*r4**2
f(395) = a1*r1*r3**2*r4**2 - a2*r2*r3**2*r4**2 + a3*r1**2*r2**2*r3 - a4* &
      r1**2*r2**2*r4
f(396) = -a1*r1*r3**3*r4 + a2*r2*r3*r4**3 - a3*r1**3*r2*r3 + a4*r1*r2**3 &
      *r4
f(397) = -a1*r2*r3*r4 + a2*r1*r3*r4 - a3*r1*r2*r4 + a4*r1*r2*r3
f(398) = -a1**2*r2*r3*r4 + a2**2*r1*r3*r4 - a3**2*r1*r2*r4 + a4**2*r1*r2 &
      *r3
f(399) = -a1**3*r2*r3*r4 + a2**3*r1*r3*r4 - a3**3*r1*r2*r4 + a4**3*r1*r2 &
      *r3
f(400) = -a1*r2**2*r3*r4 + a2*r1**2*r3*r4 - a3*r1*r2*r4**2 + a4*r1*r2*r3 &
      **2
f(401) = a1**2*r2**2*r3*r4 - a2**2*r1**2*r3*r4 + a3**2*r1*r2*r4**2 - a4 &
      **2*r1*r2*r3**2
f(402) = -a1*r2**3*r3*r4 + a2*r1**3*r3*r4 - a3*r1*r2*r4**3 + a4*r1*r2*r3 &
      **3
f(403) = a1*r2*r3**2*r4 - a2*r1*r3*r4**2 + a3*r1**2*r2*r4 - a4*r1*r2**2* &
      r3
f(404) = -a1**2*r2*r3**2*r4 + a2**2*r1*r3*r4**2 - a3**2*r1**2*r2*r4 + a4 &
      **2*r1*r2**2*r3
f(405) = a1*r2**2*r3**2*r4 - a2*r1**2*r3*r4**2 + a3*r1**2*r2*r4**2 - a4* &
      r1*r2**2*r3**2
f(406) = -a1*r2*r3**3*r4 + a2*r1*r3*r4**3 - a3*r1**3*r2*r4 + a4*r1*r2**3 &
      *r3
f(407) = -a1*r2*r3*r4**2 + a2*r1*r3**2*r4 - a3*r1*r2**2*r4 + a4*r1**2*r2 &
      *r3
f(408) = -a1**2*r2*r3*r4**2 + a2**2*r1*r3**2*r4 - a3**2*r1*r2**2*r4 + a4 &
      **2*r1**2*r2*r3
f(409) = -a1*r2**2*r3*r4**2 + a2*r1**2*r3**2*r4 - a3*r1*r2**2*r4**2 + a4 &
      *r1**2*r2*r3**2
f(410) = -a1*r2*r3**2*r4**2 + a2*r1*r3**2*r4**2 - a3*r1**2*r2**2*r4 + a4 &
      *r1**2*r2**2*r3
f(411) = -a1*r2*r3*r4**3 + a2*r1*r3**3*r4 - a3*r1*r2**3*r4 + a4*r1**3*r2 &
      *r3
f(412) = b1**2*r1*r2*r3 - b1**2*r1*r2*r4 + b2**2*r1*r3*r4 - b2**2*r2*r3* &
      r4
f(413) = b1**2*r1*r2*r3**2 - b1**2*r1*r2*r4**2 + b2**2*r1**2*r3*r4 - b2 &
      **2*r2**2*r3*r4
f(414) = b1**2*r1**2*r2*r4 - b1**2*r1*r2**2*r3 - b2**2*r1*r3*r4**2 + b2 &
      **2*r2*r3**2*r4
f(415) = b1**2*r1**2*r2*r3 - b1**2*r1*r2**2*r4 + b2**2*r1*r3**2*r4 - b2 &
      **2*r2*r3*r4**2
f(416) = -b1**2*r1*r3*r4 + b1**2*r2*r3*r4 - b2**2*r1*r2*r3 + b2**2*r1*r2 &
      *r4
f(417) = -b1**2*r1**2*r3*r4 + b1**2*r2**2*r3*r4 - b2**2*r1*r2*r3**2 + b2 &
      **2*r1*r2*r4**2
f(418) = b1**2*r1*r3*r4**2 - b1**2*r2*r3**2*r4 - b2**2*r1**2*r2*r4 + b2 &
      **2*r1*r2**2*r3
f(419) = b1**2*r1*r3**2*r4 - b1**2*r2*r3*r4**2 + b2**2*r1**2*r2*r3 - b2 &
      **2*r1*r2**2*r4
f(420) = dtau**2*(-r1*r2*r3 + r1*r2*r4 - r1*r3*r4 + r2*r3*r4)
f(421) = dtau**2*(r1**2*r3*r4 + r1*r2*r3**2 - r1*r2*r4**2 - r2**2*r3*r4)
f(422) = dtau**2*(r1**2*r2*r4 - r1*r2**2*r3 - r1*r3*r4**2 + r2*r3**2*r4)
f(423) = dtau**2*(-r1**2*r2*r3 + r1*r2**2*r4 - r1*r3**2*r4 + r2*r3*r4**2 &
      )
f(424) = -a1**2*a2*r1*r2 + a1*a2**2*r1*r2 - a3**2*a4*r3*r4 + a3*a4**2*r3 &
      *r4
f(425) = -a1**3*a2*r1*r2 + a1*a2**3*r1*r2 - a3**3*a4*r3*r4 + a3*a4**3*r3 &
      *r4
f(426) = -a1*a2*r1**2*r2 + a1*a2*r1*r2**2 - a3*a4*r3**2*r4 + a3*a4*r3*r4 &
      **2
f(427) = a1**2*a2*r1**2*r2 - a1*a2**2*r1*r2**2 + a3**2*a4*r3**2*r4 - a3* &
      a4**2*r3*r4**2
f(428) = a1**2*a2*r1*r2**2 - a1*a2**2*r1**2*r2 + a3**2*a4*r3*r4**2 - a3* &
      a4**2*r3**2*r4
f(429) = -a1*a2*r1**3*r2 + a1*a2*r1*r2**3 - a3*a4*r3**3*r4 + a3*a4*r3*r4 &
      **3
f(430) = a1*a2*r1**2*r2 - a1*a2*r1*r2**2 + a3*a4*r3**2*r4 - a3*a4*r3*r4 &
      **2
f(431) = a1*a2*r1**3*r2 - a1*a2*r1*r2**3 + a3*a4*r3**3*r4 - a3*a4*r3*r4 &
      **3
f(432) = -a1*a3*r1*r2 - a1*a3*r3*r4 + a2*a4*r1*r2 + a2*a4*r3*r4
f(433) = a1**2*a3*r3*r4 + a1*a3**2*r1*r2 - a2**2*a4*r3*r4 - a2*a4**2*r1* &
      r2
f(434) = -a1**3*a3*r3*r4 - a1*a3**3*r1*r2 + a2**3*a4*r3*r4 + a2*a4**3*r1 &
      *r2
f(435) = a1**2*a3*r1*r2 + a1*a3**2*r3*r4 - a2**2*a4*r1*r2 - a2*a4**2*r3* &
      r4
f(436) = a1**2*a3**2*r1*r2 + a1**2*a3**2*r3*r4 - a2**2*a4**2*r1*r2 - a2 &
      **2*a4**2*r3*r4
f(437) = -a1**3*a3*r1*r2 - a1*a3**3*r3*r4 + a2**3*a4*r1*r2 + a2*a4**3*r3 &
      *r4
f(438) = a1*a3*r1*r2**2 + a1*a3*r3*r4**2 - a2*a4*r1**2*r2 - a2*a4*r3**2* &
      r4
f(439) = -a1**2*a3*r3*r4**2 - a1*a3**2*r1*r2**2 + a2**2*a4*r3**2*r4 + a2 &
      *a4**2*r1**2*r2
f(440) = a1**2*a3*r1*r2**2 + a1*a3**2*r3*r4**2 - a2**2*a4*r1**2*r2 - a2* &
      a4**2*r3**2*r4
f(441) = -a1*a3*r1*r2**3 - a1*a3*r3*r4**3 + a2*a4*r1**3*r2 + a2*a4*r3**3 &
      *r4
f(442) = a1*a3*r1**2*r2 + a1*a3*r3**2*r4 - a2*a4*r1*r2**2 - a2*a4*r3*r4 &
      **2
f(443) = -a1**2*a3*r3**2*r4 - a1*a3**2*r1**2*r2 + a2**2*a4*r3*r4**2 + a2 &
      *a4**2*r1*r2**2
f(444) = a1**2*a3*r1**2*r2 + a1*a3**2*r3**2*r4 - a2**2*a4*r1*r2**2 - a2* &
      a4**2*r3*r4**2
f(445) = -a1*a3*r1**2*r2**2 - a1*a3*r3**2*r4**2 + a2*a4*r1**2*r2**2 + a2 &
      *a4*r3**2*r4**2
f(446) = -a1*a3*r1**3*r2 - a1*a3*r3**3*r4 + a2*a4*r1*r2**3 + a2*a4*r3*r4 &
      **3
f(447) = -a1*a4*r1*r2 + a1*a4*r3*r4 + a2*a3*r1*r2 - a2*a3*r3*r4
f(448) = a1**2*a4*r3*r4 - a1*a4**2*r1*r2 - a2**2*a3*r3*r4 + a2*a3**2*r1* &
      r2
f(449) = -a1**3*a4*r3*r4 + a1*a4**3*r1*r2 + a2**3*a3*r3*r4 - a2*a3**3*r1 &
      *r2
f(450) = a1**2*a4*r1*r2 - a1*a4**2*r3*r4 - a2**2*a3*r1*r2 + a2*a3**2*r3* &
      r4
f(451) = -a1**2*a4**2*r1*r2 + a1**2*a4**2*r3*r4 + a2**2*a3**2*r1*r2 - a2 &
      **2*a3**2*r3*r4
f(452) = a1**3*a4*r1*r2 - a1*a4**3*r3*r4 - a2**3*a3*r1*r2 + a2*a3**3*r3* &
      r4
f(453) = a1*a4*r1*r2**2 - a1*a4*r3**2*r4 - a2*a3*r1**2*r2 + a2*a3*r3*r4 &
      **2
f(454) = -a1**2*a4*r3**2*r4 + a1*a4**2*r1*r2**2 + a2**2*a3*r3*r4**2 - a2 &
      *a3**2*r1**2*r2
f(455) = -a1**2*a4*r1*r2**2 + a1*a4**2*r3**2*r4 + a2**2*a3*r1**2*r2 - a2 &
      *a3**2*r3*r4**2
f(456) = a1*a4*r1*r2**3 - a1*a4*r3**3*r4 - a2*a3*r1**3*r2 + a2*a3*r3*r4 &
      **3
f(457) = a1*a4*r1**2*r2 - a1*a4*r3*r4**2 - a2*a3*r1*r2**2 + a2*a3*r3**2* &
      r4
f(458) = -a1**2*a4*r3*r4**2 + a1*a4**2*r1**2*r2 + a2**2*a3*r3**2*r4 - a2 &
      *a3**2*r1*r2**2
f(459) = -a1**2*a4*r1**2*r2 + a1*a4**2*r3*r4**2 + a2**2*a3*r1*r2**2 - a2 &
      *a3**2*r3**2*r4
f(460) = a1*a4*r1**2*r2**2 - a1*a4*r3**2*r4**2 - a2*a3*r1**2*r2**2 + a2* &
      a3*r3**2*r4**2
f(461) = a1*a4*r1**3*r2 - a1*a4*r3*r4**3 - a2*a3*r1*r2**3 + a2*a3*r3**3* &
      r4
f(462) = -a1*b1**2*r1*r2 + a2*b1**2*r1*r2 - a3*b2**2*r3*r4 + a4*b2**2*r3 &
      *r4
f(463) = -a1**2*b1**2*r1*r2 + a2**2*b1**2*r1*r2 - a3**2*b2**2*r3*r4 + a4 &
      **2*b2**2*r3*r4
f(464) = a1*b1**2*r1*r2**2 - a2*b1**2*r1**2*r2 + a3*b2**2*r3*r4**2 - a4* &
      b2**2*r3**2*r4
f(465) = a1*b1**2*r1**2*r2 - a2*b1**2*r1*r2**2 + a3*b2**2*r3**2*r4 - a4* &
      b2**2*r3*r4**2
f(466) = -a1*b2**2*r1*r2 + a2*b2**2*r1*r2 - a3*b1**2*r3*r4 + a4*b1**2*r3 &
      *r4
f(467) = a1**2*b2**2*r1*r2 - a2**2*b2**2*r1*r2 + a3**2*b1**2*r3*r4 - a4 &
      **2*b1**2*r3*r4
f(468) = -a1*b2**2*r1*r2**2 + a2*b2**2*r1**2*r2 - a3*b1**2*r3*r4**2 + a4 &
      *b1**2*r3**2*r4
f(469) = -a1*b2**2*r1**2*r2 + a2*b2**2*r1*r2**2 - a3*b1**2*r3**2*r4 + a4 &
      *b1**2*r3*r4**2
f(470) = dtau**2*(-a1*r1*r2 + a2*r1*r2 - a3*r3*r4 + a4*r3*r4)
f(471) = dtau**2*(a1**2*r1*r2 - a2**2*r1*r2 + a3**2*r3*r4 - a4**2*r3*r4)
f(472) = dtau**2*(-a1*r1*r2**2 + a2*r1**2*r2 - a3*r3*r4**2 + a4*r3**2*r4 &
      )
f(473) = dtau**2*(-a1*r1**2*r2 + a2*r1*r2**2 - a3*r3**2*r4 + a4*r3*r4**2 &
      )
f(474) = -a1**2*a2*r3*r4 + a1*a2**2*r3*r4 - a3**2*a4*r1*r2 + a3*a4**2*r1 &
      *r2
f(475) = -a1**3*a2*r3*r4 + a1*a2**3*r3*r4 - a3**3*a4*r1*r2 + a3*a4**3*r1 &
      *r2
f(476) = a1*a2*r3**2*r4 - a1*a2*r3*r4**2 + a3*a4*r1**2*r2 - a3*a4*r1*r2 &
      **2
f(477) = a1**2*a2*r3**2*r4 - a1*a2**2*r3*r4**2 + a3**2*a4*r1**2*r2 - a3* &
      a4**2*r1*r2**2
f(478) = a1**2*a2*r3*r4**2 - a1*a2**2*r3**2*r4 + a3**2*a4*r1*r2**2 - a3* &
      a4**2*r1**2*r2
f(479) = a1*a2*r3**3*r4 - a1*a2*r3*r4**3 + a3*a4*r1**3*r2 - a3*a4*r1*r2 &
      **3
f(480) = -a1*a2*r3**2*r4 + a1*a2*r3*r4**2 - a3*a4*r1**2*r2 + a3*a4*r1*r2 &
      **2
f(481) = a1*b2**2*r3*r4 - a2*b2**2*r3*r4 + a3*b1**2*r1*r2 - a4*b1**2*r1* &
      r2
f(482) = a1**2*b2**2*r3*r4 - a2**2*b2**2*r3*r4 + a3**2*b1**2*r1*r2 - a4 &
      **2*b1**2*r1*r2
f(483) = a1*b2**2*r3*r4**2 - a2*b2**2*r3**2*r4 + a3*b1**2*r1*r2**2 - a4* &
      b1**2*r1**2*r2
f(484) = a1*b2**2*r3**2*r4 - a2*b2**2*r3*r4**2 + a3*b1**2*r1**2*r2 - a4* &
      b1**2*r1*r2**2
f(485) = a1*b1**2*r3*r4 - a2*b1**2*r3*r4 + a3*b2**2*r1*r2 - a4*b2**2*r1* &
      r2
f(486) = -a1**2*b1**2*r3*r4 + a2**2*b1**2*r3*r4 - a3**2*b2**2*r1*r2 + a4 &
      **2*b2**2*r1*r2
f(487) = a1*b1**2*r3*r4**2 - a2*b1**2*r3**2*r4 + a3*b2**2*r1*r2**2 - a4* &
      b2**2*r1**2*r2
f(488) = a1*b1**2*r3**2*r4 - a2*b1**2*r3*r4**2 + a3*b2**2*r1**2*r2 - a4* &
      b2**2*r1*r2**2
f(489) = dtau**2*(-a1*r3*r4 + a2*r3*r4 - a3*r1*r2 + a4*r1*r2)
f(490) = dtau**2*(a1**2*r3*r4 - a2**2*r3*r4 + a3**2*r1*r2 - a4**2*r1*r2)
f(491) = dtau**2*(-a1*r3*r4**2 + a2*r3**2*r4 - a3*r1*r2**2 + a4*r1**2*r2 &
      )
f(492) = dtau**2*(-a1*r3**2*r4 + a2*r3*r4**2 - a3*r1**2*r2 + a4*r1*r2**2 &
      )
f(493) = b1*b2*(-r1**2*r2 + r1*r2**2 - r3**2*r4 + r3*r4**2)
f(494) = b1*b2*(-r1**3*r2 + r1*r2**3 - r3**3*r4 + r3*r4**3)
f(495) = b1*b2*(r1**2*r2 - r1*r2**2 + r3**2*r4 - r3*r4**2)
f(496) = b1*b2*(r1**3*r2 - r1*r2**3 + r3**3*r4 - r3*r4**3)
f(497) = dtau*(b1*r1*r2 + b2*r3*r4)
f(498) = dtau**3*(b1*r1*r2 + b2*r3*r4)
f(499) = dtau*(b1**3*r1*r2 + b2**3*r3*r4)
f(500) = dtau*(b1*r1**2*r2 + b1*r1*r2**2 + b2*r3**2*r4 + b2*r3*r4**2)
f(501) = dtau*(b1*r1**3*r2 + b1*r1*r2**3 + b2*r3**3*r4 + b2*r3*r4**3)
f(502) = dtau*(b1*r1**2*r2**2 + b2*r3**2*r4**2)
f(503) = dtau*(b1*r3*r4 + b2*r1*r2)
f(504) = dtau**3*(b1*r3*r4 + b2*r1*r2)
f(505) = dtau*(b1**3*r3*r4 + b2**3*r1*r2)
f(506) = dtau*(b1*r3**2*r4 + b1*r3*r4**2 + b2*r1**2*r2 + b2*r1*r2**2)
f(507) = dtau*(b1*r3**3*r4 + b1*r3*r4**3 + b2*r1**3*r2 + b2*r1*r2**3)
f(508) = dtau*(b1*r3**2*r4**2 + b2*r1**2*r2**2)
f(509) = a1*a2*r1*r3 - a1*a2*r2*r4 + a3*a4*r1*r3 - a3*a4*r2*r4
f(510) = -a1**2*a2*r2*r4 + a1*a2**2*r1*r3 - a3**2*a4*r2*r4 + a3*a4**2*r1 &
      *r3
f(511) = a1**3*a2*r2*r4 - a1*a2**3*r1*r3 + a3**3*a4*r2*r4 - a3*a4**3*r1* &
      r3
f(512) = -a1**2*a2*r1*r3 + a1*a2**2*r2*r4 - a3**2*a4*r1*r3 + a3*a4**2*r2 &
      *r4
f(513) = a1**2*a2**2*r1*r3 - a1**2*a2**2*r2*r4 + a3**2*a4**2*r1*r3 - a3 &
      **2*a4**2*r2*r4
f(514) = a1**3*a2*r1*r3 - a1*a2**3*r2*r4 + a3**3*a4*r1*r3 - a3*a4**3*r2* &
      r4
f(515) = -a1*a2*r1*r3**2 + a1*a2*r2*r4**2 - a3*a4*r1**2*r3 + a3*a4*r2**2 &
      *r4
f(516) = -a1**2*a2*r2*r4**2 + a1*a2**2*r1*r3**2 - a3**2*a4*r2**2*r4 + a3 &
      *a4**2*r1**2*r3
f(517) = a1**2*a2*r1*r3**2 - a1*a2**2*r2*r4**2 + a3**2*a4*r1**2*r3 - a3* &
      a4**2*r2**2*r4
f(518) = -a1*a2*r1*r3**3 + a1*a2*r2*r4**3 - a3*a4*r1**3*r3 + a3*a4*r2**3 &
      *r4
f(519) = -a1*a2*r1**2*r3 + a1*a2*r2**2*r4 - a3*a4*r1*r3**2 + a3*a4*r2*r4 &
      **2
f(520) = -a1**2*a2*r2**2*r4 + a1*a2**2*r1**2*r3 - a3**2*a4*r2*r4**2 + a3 &
      *a4**2*r1*r3**2
f(521) = a1**2*a2*r1**2*r3 - a1*a2**2*r2**2*r4 + a3**2*a4*r1*r3**2 - a3* &
      a4**2*r2*r4**2
f(522) = -a1*a2*r1**2*r3**2 + a1*a2*r2**2*r4**2 - a3*a4*r1**2*r3**2 + a3 &
      *a4*r2**2*r4**2
f(523) = a1*a2*r1**3*r3 - a1*a2*r2**3*r4 + a3*a4*r1*r3**3 - a3*a4*r2*r4 &
      **3
f(524) = -a1*a3*r1*r3 + a2*a4*r2*r4
f(525) = a1**2*a3*r1*r3 + a1*a3**2*r1*r3 - a2**2*a4*r2*r4 - a2*a4**2*r2* &
      r4
f(526) = -a1**3*a3*r1*r3 - a1*a3**3*r1*r3 + a2**3*a4*r2*r4 + a2*a4**3*r2 &
      *r4
f(527) = -a1**2*a3**2*r1*r3 + a2**2*a4**2*r2*r4
f(528) = a1*a3*r1**2*r3 + a1*a3*r1*r3**2 - a2*a4*r2**2*r4 - a2*a4*r2*r4 &
      **2
f(529) = -a1**2*a3*r1**2*r3 - a1*a3**2*r1*r3**2 + a2**2*a4*r2**2*r4 + a2 &
      *a4**2*r2*r4**2
f(530) = -a1**2*a3*r1*r3**2 - a1*a3**2*r1**2*r3 + a2**2*a4*r2*r4**2 + a2 &
      *a4**2*r2**2*r4
f(531) = -a1*a3*r1**3*r3 - a1*a3*r1*r3**3 + a2*a4*r2**3*r4 + a2*a4*r2*r4 &
      **3
f(532) = -a1*a3*r1**2*r3**2 + a2*a4*r2**2*r4**2
f(533) = -a1*a4*r1*r3 + a1*a4*r2*r4 - a2*a3*r1*r3 + a2*a3*r2*r4
f(534) = a1**2*a4*r2*r4 - a1*a4**2*r1*r3 - a2**2*a3*r1*r3 + a2*a3**2*r2* &
      r4
f(535) = -a1**3*a4*r2*r4 + a1*a4**3*r1*r3 + a2**3*a3*r1*r3 - a2*a3**3*r2 &
      *r4
f(536) = -a1**2*a4*r1*r3 + a1*a4**2*r2*r4 + a2**2*a3*r2*r4 - a2*a3**2*r1 &
      *r3
f(537) = -a1**2*a4**2*r1*r3 + a1**2*a4**2*r2*r4 - a2**2*a3**2*r1*r3 + a2 &
      **2*a3**2*r2*r4
f(538) = -a1**3*a4*r1*r3 + a1*a4**3*r2*r4 + a2**3*a3*r2*r4 - a2*a3**3*r1 &
      *r3
f(539) = -a1*a4*r1*r3**2 + a1*a4*r2**2*r4 - a2*a3*r1**2*r3 + a2*a3*r2*r4 &
      **2
f(540) = -a1**2*a4*r2**2*r4 + a1*a4**2*r1*r3**2 + a2**2*a3*r1**2*r3 - a2 &
      *a3**2*r2*r4**2
f(541) = -a1**2*a4*r1*r3**2 + a1*a4**2*r2**2*r4 + a2**2*a3*r2*r4**2 - a2 &
      *a3**2*r1**2*r3
f(542) = -a1*a4*r1*r3**3 + a1*a4*r2**3*r4 - a2*a3*r1**3*r3 + a2*a3*r2*r4 &
      **3
f(543) = -a1*a4*r1**2*r3 + a1*a4*r2*r4**2 - a2*a3*r1*r3**2 + a2*a3*r2**2 &
      *r4
f(544) = -a1**2*a4*r2*r4**2 + a1*a4**2*r1**2*r3 + a2**2*a3*r1*r3**2 - a2 &
      *a3**2*r2**2*r4
f(545) = -a1**2*a4*r1**2*r3 + a1*a4**2*r2*r4**2 + a2**2*a3*r2**2*r4 - a2 &
      *a3**2*r1*r3**2
f(546) = -a1*a4*r1**2*r3**2 + a1*a4*r2**2*r4**2 - a2*a3*r1**2*r3**2 + a2 &
      *a3*r2**2*r4**2
f(547) = -a1*a4*r1**3*r3 + a1*a4*r2*r4**3 - a2*a3*r1*r3**3 + a2*a3*r2**3 &
      *r4
f(548) = -a1*b1**2*r1*r3 + a2*b1**2*r2*r4 - a3*b2**2*r1*r3 + a4*b2**2*r2 &
      *r4
f(549) = -a1**2*b1**2*r1*r3 + a2**2*b1**2*r2*r4 - a3**2*b2**2*r1*r3 + a4 &
      **2*b2**2*r2*r4
f(550) = a1*b1**2*r1*r3**2 - a2*b1**2*r2*r4**2 + a3*b2**2*r1**2*r3 - a4* &
      b2**2*r2**2*r4
f(551) = a1*b1**2*r1**2*r3 - a2*b1**2*r2**2*r4 + a3*b2**2*r1*r3**2 - a4* &
      b2**2*r2*r4**2
f(552) = a1*b2**2*r1*r3 - a2*b2**2*r2*r4 + a3*b1**2*r1*r3 - a4*b1**2*r2* &
      r4
f(553) = -a1**2*b2**2*r1*r3 + a2**2*b2**2*r2*r4 - a3**2*b1**2*r1*r3 + a4 &
      **2*b1**2*r2*r4
f(554) = a1*b2**2*r1*r3**2 - a2*b2**2*r2*r4**2 + a3*b1**2*r1**2*r3 - a4* &
      b1**2*r2**2*r4
f(555) = a1*b2**2*r1**2*r3 - a2*b2**2*r2**2*r4 + a3*b1**2*r1*r3**2 - a4* &
      b1**2*r2*r4**2
f(556) = dtau**2*(-a1*r1*r3 + a2*r2*r4 - a3*r1*r3 + a4*r2*r4)
f(557) = dtau**2*(-a1**2*r1*r3 + a2**2*r2*r4 - a3**2*r1*r3 + a4**2*r2*r4 &
      )
f(558) = dtau**2*(a1*r1*r3**2 - a2*r2*r4**2 + a3*r1**2*r3 - a4*r2**2*r4)
f(559) = dtau**2*(a1*r1**2*r3 - a2*r2**2*r4 + a3*r1*r3**2 - a4*r2*r4**2)
f(560) = -a1*a3*r2*r4 + a2*a4*r1*r3
f(561) = a1**2*a3*r2*r4 + a1*a3**2*r2*r4 - a2**2*a4*r1*r3 - a2*a4**2*r1* &
      r3
f(562) = a1**3*a3*r2*r4 + a1*a3**3*r2*r4 - a2**3*a4*r1*r3 - a2*a4**3*r1* &
      r3
f(563) = a1**2*a3**2*r2*r4 - a2**2*a4**2*r1*r3
f(564) = a1*a3*r2**2*r4 + a1*a3*r2*r4**2 - a2*a4*r1**2*r3 - a2*a4*r1*r3 &
      **2
f(565) = -a1**2*a3*r2**2*r4 - a1*a3**2*r2*r4**2 + a2**2*a4*r1**2*r3 + a2 &
      *a4**2*r1*r3**2
f(566) = -a1**2*a3*r2*r4**2 - a1*a3**2*r2**2*r4 + a2**2*a4*r1*r3**2 + a2 &
      *a4**2*r1**2*r3
f(567) = -a1*a3*r2**3*r4 - a1*a3*r2*r4**3 + a2*a4*r1**3*r3 + a2*a4*r1*r3 &
      **3
f(568) = -a1*a3*r2**2*r4**2 + a2*a4*r1**2*r3**2
f(569) = a1*b1**2*r2*r4 - a2*b1**2*r1*r3 + a3*b2**2*r2*r4 - a4*b2**2*r1* &
      r3
f(570) = a1**2*b1**2*r2*r4 - a2**2*b1**2*r1*r3 + a3**2*b2**2*r2*r4 - a4 &
      **2*b2**2*r1*r3
f(571) = -a1*b1**2*r2*r4**2 + a2*b1**2*r1*r3**2 - a3*b2**2*r2**2*r4 + a4 &
      *b2**2*r1**2*r3
f(572) = -a1*b1**2*r2**2*r4 + a2*b1**2*r1**2*r3 - a3*b2**2*r2*r4**2 + a4 &
      *b2**2*r1*r3**2
f(573) = -a1*b2**2*r2*r4 + a2*b2**2*r1*r3 - a3*b1**2*r2*r4 + a4*b1**2*r1 &
      *r3
f(574) = -a1**2*b2**2*r2*r4 + a2**2*b2**2*r1*r3 - a3**2*b1**2*r2*r4 + a4 &
      **2*b1**2*r1*r3
f(575) = a1*b2**2*r2*r4**2 - a2*b2**2*r1*r3**2 + a3*b1**2*r2**2*r4 - a4* &
      b1**2*r1**2*r3
f(576) = -a1*b2**2*r2**2*r4 + a2*b2**2*r1**2*r3 - a3*b1**2*r2*r4**2 + a4 &
      *b1**2*r1*r3**2
f(577) = dtau**2*(a1*r2*r4 - a2*r1*r3 + a3*r2*r4 - a4*r1*r3)
f(578) = dtau**2*(a1**2*r2*r4 - a2**2*r1*r3 + a3**2*r2*r4 - a4**2*r1*r3)
f(579) = dtau**2*(-a1*r2*r4**2 + a2*r1*r3**2 - a3*r2**2*r4 + a4*r1**2*r3 &
      )
f(580) = dtau**2*(-a1*r2**2*r4 + a2*r1**2*r3 - a3*r2*r4**2 + a4*r1*r3**2 &
      )
f(581) = b1*b2*(-r1*r3 + r2*r4)
f(582) = b1*b2*(-b1**2*r1*r3 + b1**2*r2*r4 - b2**2*r1*r3 + b2**2*r2*r4)
f(583) = b1**2*b2**2*(r1*r3 - r2*r4)
f(584) = b1*b2*(-r1**2*r3 - r1*r3**2 + r2**2*r4 + r2*r4**2)
f(585) = b1*b2*(-r1**3*r3 - r1*r3**3 + r2**3*r4 + r2*r4**3)
f(586) = b1*b2*(r1**2*r3 + r1*r3**2 - r2**2*r4 - r2*r4**2)
f(587) = b1*b2*(-r1**2*r3**2 + r2**2*r4**2)
f(588) = dtau*(b1*r1*r3 + b1*r2*r4 + b2*r1*r3 + b2*r2*r4)
f(589) = dtau**3*(b1*r1*r3 + b1*r2*r4 + b2*r1*r3 + b2*r2*r4)
f(590) = dtau**2*(-b1**2*r1*r3 + b1**2*r2*r4 - b2**2*r1*r3 + b2**2*r2*r4 &
      )
f(591) = dtau*(b1**3*r1*r3 + b1**3*r2*r4 + b2**3*r1*r3 + b2**3*r2*r4)
f(592) = dtau*(b1*r1*r3**2 + b1*r2*r4**2 + b2*r1**2*r3 + b2*r2**2*r4)
f(593) = dtau*(b1*r1*r3**3 + b1*r2*r4**3 + b2*r1**3*r3 + b2*r2**3*r4)
f(594) = dtau*(b1*r1**2*r3 + b1*r2**2*r4 + b2*r1*r3**2 + b2*r2*r4**2)
f(595) = dtau*(b1*r1**2*r3**2 + b1*r2**2*r4**2 + b2*r1**2*r3**2 + b2*r2 &
      **2*r4**2)
f(596) = dtau*(b1*r1**3*r3 + b1*r2**3*r4 + b2*r1*r3**3 + b2*r2*r4**3)
f(597) = -a1*a2*r1*r4 + a1*a2*r2*r3 + a3*a4*r1*r4 - a3*a4*r2*r3
f(598) = a1**2*a2*r2*r3 - a1*a2**2*r1*r4 + a3**2*a4*r1*r4 - a3*a4**2*r2* &
      r3
f(599) = -a1**3*a2*r2*r3 + a1*a2**3*r1*r4 - a3**3*a4*r1*r4 + a3*a4**3*r2 &
      *r3
f(600) = a1**2*a2*r1*r4 - a1*a2**2*r2*r3 + a3**2*a4*r2*r3 - a3*a4**2*r1* &
      r4
f(601) = -a1**2*a2**2*r1*r4 + a1**2*a2**2*r2*r3 + a3**2*a4**2*r1*r4 - a3 &
      **2*a4**2*r2*r3
f(602) = -a1**3*a2*r1*r4 + a1*a2**3*r2*r3 - a3**3*a4*r2*r3 + a3*a4**3*r1 &
      *r4
f(603) = a1*a2*r1*r4**2 - a1*a2*r2*r3**2 - a3*a4*r1**2*r4 + a3*a4*r2**2* &
      r3
f(604) = -a1**2*a2*r2*r3**2 + a1*a2**2*r1*r4**2 - a3**2*a4*r1**2*r4 + a3 &
      *a4**2*r2**2*r3
f(605) = -a1**2*a2*r1*r4**2 + a1*a2**2*r2*r3**2 - a3**2*a4*r2**2*r3 + a3 &
      *a4**2*r1**2*r4
f(606) = a1*a2*r1*r4**3 - a1*a2*r2*r3**3 - a3*a4*r1**3*r4 + a3*a4*r2**3* &
      r3
f(607) = -a1*a2*r1**2*r4 + a1*a2*r2**2*r3 + a3*a4*r1*r4**2 - a3*a4*r2*r3 &
      **2
f(608) = -a1**2*a2*r2**2*r3 + a1*a2**2*r1**2*r4 - a3**2*a4*r1*r4**2 + a3 &
      *a4**2*r2*r3**2
f(609) = -a1**2*a2*r1**2*r4 + a1*a2**2*r2**2*r3 - a3**2*a4*r2*r3**2 + a3 &
      *a4**2*r1*r4**2
f(610) = a1*a2*r1**2*r4**2 - a1*a2*r2**2*r3**2 - a3*a4*r1**2*r4**2 + a3* &
      a4*r2**2*r3**2
f(611) = -a1*a2*r1**3*r4 + a1*a2*r2**3*r3 + a3*a4*r1*r4**3 - a3*a4*r2*r3 &
      **3
f(612) = a1*a3*r1*r4 + a1*a3*r2*r3 - a2*a4*r1*r4 - a2*a4*r2*r3
f(613) = -a1**2*a3*r2*r3 - a1*a3**2*r1*r4 + a2**2*a4*r1*r4 + a2*a4**2*r2 &
      *r3
f(614) = a1**3*a3*r2*r3 + a1*a3**3*r1*r4 - a2**3*a4*r1*r4 - a2*a4**3*r2* &
      r3
f(615) = a1**2*a3*r1*r4 + a1*a3**2*r2*r3 - a2**2*a4*r2*r3 - a2*a4**2*r1* &
      r4
f(616) = -a1**2*a3**2*r1*r4 - a1**2*a3**2*r2*r3 + a2**2*a4**2*r1*r4 + a2 &
      **2*a4**2*r2*r3
f(617) = a1**3*a3*r1*r4 + a1*a3**3*r2*r3 - a2**3*a4*r2*r3 - a2*a4**3*r1* &
      r4
f(618) = -a1*a3*r1*r4**2 - a1*a3*r2**2*r3 + a2*a4*r1**2*r4 + a2*a4*r2*r3 &
      **2
f(619) = a1**2*a3*r2**2*r3 + a1*a3**2*r1*r4**2 - a2**2*a4*r1**2*r4 - a2* &
      a4**2*r2*r3**2
f(620) = a1**2*a3*r1*r4**2 + a1*a3**2*r2**2*r3 - a2**2*a4*r2*r3**2 - a2* &
      a4**2*r1**2*r4
f(621) = a1*a3*r1*r4**3 + a1*a3*r2**3*r3 - a2*a4*r1**3*r4 - a2*a4*r2*r3 &
      **3
f(622) = -a1*a3*r1**2*r4 - a1*a3*r2*r3**2 + a2*a4*r1*r4**2 + a2*a4*r2**2 &
      *r3
f(623) = a1**2*a3*r2*r3**2 + a1*a3**2*r1**2*r4 - a2**2*a4*r1*r4**2 - a2* &
      a4**2*r2**2*r3
f(624) = a1**2*a3*r1**2*r4 + a1*a3**2*r2*r3**2 - a2**2*a4*r2**2*r3 - a2* &
      a4**2*r1*r4**2
f(625) = -a1*a3*r1**2*r4**2 - a1*a3*r2**2*r3**2 + a2*a4*r1**2*r4**2 + a2 &
      *a4*r2**2*r3**2
f(626) = a1*a3*r1**3*r4 + a1*a3*r2*r3**3 - a2*a4*r1*r4**3 - a2*a4*r2**3* &
      r3
f(627) = -a1**2*a4*r1*r4 + a1*a4**2*r1*r4 + a2**2*a3*r2*r3 - a2*a3**2*r2 &
      *r3
f(628) = a1**3*a4*r1*r4 - a1*a4**3*r1*r4 - a2**3*a3*r2*r3 + a2*a3**3*r2* &
      r3
f(629) = a1*a4*r1**2*r4 - a1*a4*r1*r4**2 - a2*a3*r2**2*r3 + a2*a3*r2*r3 &
      **2
f(630) = a1**2*a4*r1**2*r4 - a1*a4**2*r1*r4**2 - a2**2*a3*r2**2*r3 + a2* &
      a3**2*r2*r3**2
f(631) = a1**2*a4*r1*r4**2 - a1*a4**2*r1**2*r4 - a2**2*a3*r2*r3**2 + a2* &
      a3**2*r2**2*r3
f(632) = -a1*a4*r1**3*r4 + a1*a4*r1*r4**3 + a2*a3*r2**3*r3 - a2*a3*r2*r3 &
      **3
f(633) = -a1*a4*r1**2*r4 + a1*a4*r1*r4**2 + a2*a3*r2**2*r3 - a2*a3*r2*r3 &
      **2
f(634) = -a1*b1**2*r1*r4 + a2*b1**2*r2*r3 - a3*b2**2*r2*r3 + a4*b2**2*r1 &
      *r4
f(635) = -a1**2*b1**2*r1*r4 + a2**2*b1**2*r2*r3 - a3**2*b2**2*r2*r3 + a4 &
      **2*b2**2*r1*r4
f(636) = a1*b1**2*r1*r4**2 - a2*b1**2*r2*r3**2 + a3*b2**2*r2**2*r3 - a4* &
      b2**2*r1**2*r4
f(637) = a1*b1**2*r1**2*r4 - a2*b1**2*r2**2*r3 + a3*b2**2*r2*r3**2 - a4* &
      b2**2*r1*r4**2
f(638) = a1*b2**2*r1*r4 - a2*b2**2*r2*r3 + a3*b1**2*r2*r3 - a4*b1**2*r1* &
      r4
f(639) = a1**2*b2**2*r1*r4 - a2**2*b2**2*r2*r3 + a3**2*b1**2*r2*r3 - a4 &
      **2*b1**2*r1*r4
f(640) = a1*b2**2*r1*r4**2 - a2*b2**2*r2*r3**2 + a3*b1**2*r2**2*r3 - a4* &
      b1**2*r1**2*r4
f(641) = a1*b2**2*r1**2*r4 - a2*b2**2*r2**2*r3 + a3*b1**2*r2*r3**2 - a4* &
      b1**2*r1*r4**2
f(642) = dtau**2*(-a1*r1*r4 + a2*r2*r3 - a3*r2*r3 + a4*r1*r4)
f(643) = dtau**2*(-a1**2*r1*r4 + a2**2*r2*r3 - a3**2*r2*r3 + a4**2*r1*r4 &
      )
f(644) = dtau**2*(-a1*r1*r4**2 + a2*r2*r3**2 - a3*r2**2*r3 + a4*r1**2*r4 &
      )
f(645) = dtau**2*(-a1*r1**2*r4 + a2*r2**2*r3 - a3*r2*r3**2 + a4*r1*r4**2 &
      )
f(646) = -a1**2*a4*r2*r3 + a1*a4**2*r2*r3 + a2**2*a3*r1*r4 - a2*a3**2*r1 &
      *r4
f(647) = -a1**3*a4*r2*r3 + a1*a4**3*r2*r3 + a2**3*a3*r1*r4 - a2*a3**3*r1 &
      *r4
f(648) = -a1*a4*r2**2*r3 + a1*a4*r2*r3**2 + a2*a3*r1**2*r4 - a2*a3*r1*r4 &
      **2
f(649) = a1**2*a4*r2**2*r3 - a1*a4**2*r2*r3**2 - a2**2*a3*r1**2*r4 + a2* &
      a3**2*r1*r4**2
f(650) = a1**2*a4*r2*r3**2 - a1*a4**2*r2**2*r3 - a2**2*a3*r1*r4**2 + a2* &
      a3**2*r1**2*r4
f(651) = a1*a4*r2**3*r3 - a1*a4*r2*r3**3 - a2*a3*r1**3*r4 + a2*a3*r1*r4 &
      **3
f(652) = a1*a4*r2**2*r3 - a1*a4*r2*r3**2 - a2*a3*r1**2*r4 + a2*a3*r1*r4 &
      **2
f(653) = a1*b1**2*r2*r3 - a2*b1**2*r1*r4 + a3*b2**2*r1*r4 - a4*b2**2*r2* &
      r3
f(654) = a1**2*b1**2*r2*r3 - a2**2*b1**2*r1*r4 + a3**2*b2**2*r1*r4 - a4 &
      **2*b2**2*r2*r3
f(655) = -a1*b1**2*r2*r3**2 + a2*b1**2*r1*r4**2 - a3*b2**2*r1**2*r4 + a4 &
      *b2**2*r2**2*r3
f(656) = -a1*b1**2*r2**2*r3 + a2*b1**2*r1**2*r4 - a3*b2**2*r1*r4**2 + a4 &
      *b2**2*r2*r3**2
f(657) = -a1*b2**2*r2*r3 + a2*b2**2*r1*r4 - a3*b1**2*r1*r4 + a4*b1**2*r2 &
      *r3
f(658) = -a1**2*b2**2*r2*r3 + a2**2*b2**2*r1*r4 - a3**2*b1**2*r1*r4 + a4 &
      **2*b1**2*r2*r3
f(659) = -a1*b2**2*r2*r3**2 + a2*b2**2*r1*r4**2 - a3*b1**2*r1**2*r4 + a4 &
      *b1**2*r2**2*r3
f(660) = -a1*b2**2*r2**2*r3 + a2*b2**2*r1**2*r4 - a3*b1**2*r1*r4**2 + a4 &
      *b1**2*r2*r3**2
f(661) = dtau**2*(a1*r2*r3 - a2*r1*r4 + a3*r1*r4 - a4*r2*r3)
f(662) = dtau**2*(-a1**2*r2*r3 + a2**2*r1*r4 - a3**2*r1*r4 + a4**2*r2*r3 &
      )
f(663) = dtau**2*(a1*r2*r3**2 - a2*r1*r4**2 + a3*r1**2*r4 - a4*r2**2*r3)
f(664) = dtau**2*(a1*r2**2*r3 - a2*r1**2*r4 + a3*r1*r4**2 - a4*r2*r3**2)
f(665) = b1*b2*(-b1**2*r1*r4 + b1**2*r2*r3 + b2**2*r1*r4 - b2**2*r2*r3)
f(666) = b1*b2*(r1**2*r4 - r1*r4**2 - r2**2*r3 + r2*r3**2)
f(667) = b1*b2*(r1**3*r4 - r1*r4**3 - r2**3*r3 + r2*r3**3)
f(668) = b1*b2*(-r1**2*r4 + r1*r4**2 + r2**2*r3 - r2*r3**2)
f(669) = dtau*(b1*r1*r4 + b1*r2*r3 + b2*r1*r4 + b2*r2*r3)
f(670) = dtau**3*(b1*r1*r4 + b1*r2*r3 + b2*r1*r4 + b2*r2*r3)
f(671) = dtau**2*(-b1**2*r1*r4 + b1**2*r2*r3 + b2**2*r1*r4 - b2**2*r2*r3 &
      )
f(672) = dtau*(b1**3*r1*r4 + b1**3*r2*r3 + b2**3*r1*r4 + b2**3*r2*r3)
f(673) = dtau*(b1*r1*r4**2 + b1*r2*r3**2 + b2*r1**2*r4 + b2*r2**2*r3)
f(674) = dtau*(b1*r1*r4**3 + b1*r2*r3**3 + b2*r1**3*r4 + b2*r2**3*r3)
f(675) = dtau*(b1*r1**2*r4 + b1*r2**2*r3 + b2*r1*r4**2 + b2*r2*r3**2)
f(676) = dtau*(b1*r1**2*r4**2 + b1*r2**2*r3**2 + b2*r1**2*r4**2 + b2*r2 &
      **2*r3**2)
f(677) = dtau*(b1*r1**3*r4 + b1*r2**3*r3 + b2*r1*r4**3 + b2*r2*r3**3)
f(678) = a1*a2*a3*r1 - a1*a2*a4*r2 + a1*a3*a4*r3 - a2*a3*a4*r4
f(679) = a1**2*a3*a4*r3 + a1*a2*a3**2*r1 - a1*a2*a4**2*r2 - a2**2*a3*a4* &
      r4
f(680) = a1**3*a3*a4*r3 + a1*a2*a3**3*r1 - a1*a2*a4**3*r2 - a2**3*a3*a4* &
      r4
f(681) = -a1**2*a2*a4*r2 + a1*a2**2*a3*r1 + a1*a3*a4**2*r3 - a2*a3**2*a4 &
      *r4
f(682) = a1**2*a2*a4**2*r2 - a1**2*a3*a4**2*r3 - a1*a2**2*a3**2*r1 + a2 &
      **2*a3**2*a4*r4
f(683) = -a1**3*a2*a4*r2 + a1*a2**3*a3*r1 + a1*a3*a4**3*r3 - a2*a3**3*a4 &
      *r4
f(684) = -a1**2*a2*a3*r1 + a1*a2**2*a4*r2 - a1*a3**2*a4*r3 + a2*a3*a4**2 &
      *r4
f(685) = a1**2*a2*a3**2*r1 + a1**2*a3**2*a4*r3 - a1*a2**2*a4**2*r2 - a2 &
      **2*a3*a4**2*r4
f(686) = -a1**2*a2**2*a3*r1 + a1**2*a2**2*a4*r2 - a1*a3**2*a4**2*r3 + a2 &
      *a3**2*a4**2*r4
f(687) = a1**3*a2*a3*r1 - a1*a2**3*a4*r2 + a1*a3**3*a4*r3 - a2*a3*a4**3* &
      r4
f(688) = a1*a2*a3*r1**2 - a1*a2*a4*r2**2 + a1*a3*a4*r3**2 - a2*a3*a4*r4 &
      **2
f(689) = a1**2*a3*a4*r3**2 + a1*a2*a3**2*r1**2 - a1*a2*a4**2*r2**2 - a2 &
      **2*a3*a4*r4**2
f(690) = a1**2*a2*a4*r2**2 - a1*a2**2*a3*r1**2 - a1*a3*a4**2*r3**2 + a2* &
      a3**2*a4*r4**2
f(691) = -a1**2*a2*a3*r1**2 + a1*a2**2*a4*r2**2 - a1*a3**2*a4*r3**2 + a2 &
      *a3*a4**2*r4**2
f(692) = -a1*a2*a3*r1**3 + a1*a2*a4*r2**3 - a1*a3*a4*r3**3 + a2*a3*a4*r4 &
      **3
f(693) = -a1*a2*a3*r2 + a1*a2*a4*r1 - a1*a3*a4*r4 + a2*a3*a4*r3
f(694) = -a1**2*a3*a4*r4 - a1*a2*a3**2*r2 + a1*a2*a4**2*r1 + a2**2*a3*a4 &
      *r3
f(695) = -a1**3*a3*a4*r4 - a1*a2*a3**3*r2 + a1*a2*a4**3*r1 + a2**3*a3*a4 &
      *r3
f(696) = a1**2*a2*a3*r2 - a1*a2**2*a4*r1 + a1*a3**2*a4*r4 - a2*a3*a4**2* &
      r3
f(697) = -a1**2*a2*a3**2*r2 - a1**2*a3**2*a4*r4 + a1*a2**2*a4**2*r1 + a2 &
      **2*a3*a4**2*r3
f(698) = a1**3*a2*a3*r2 - a1*a2**3*a4*r1 + a1*a3**3*a4*r4 - a2*a3*a4**3* &
      r3
f(699) = -a1**2*a2*a4*r1 + a1*a2**2*a3*r2 + a1*a3*a4**2*r4 - a2*a3**2*a4 &
      *r3
f(700) = -a1**2*a2*a4**2*r1 + a1**2*a3*a4**2*r4 + a1*a2**2*a3**2*r2 - a2 &
      **2*a3**2*a4*r3
f(701) = -a1**2*a2**2*a3*r2 + a1**2*a2**2*a4*r1 - a1*a3**2*a4**2*r4 + a2 &
      *a3**2*a4**2*r3
f(702) = -a1**3*a2*a4*r1 + a1*a2**3*a3*r2 + a1*a3*a4**3*r4 - a2*a3**3*a4 &
      *r3
f(703) = -a1*a2*a3*r2**2 + a1*a2*a4*r1**2 - a1*a3*a4*r4**2 + a2*a3*a4*r3 &
      **2
f(704) = a1**2*a3*a4*r4**2 + a1*a2*a3**2*r2**2 - a1*a2*a4**2*r1**2 - a2 &
      **2*a3*a4*r3**2
f(705) = a1**2*a2*a3*r2**2 - a1*a2**2*a4*r1**2 + a1*a3**2*a4*r4**2 - a2* &
      a3*a4**2*r3**2
f(706) = -a1**2*a2*a4*r1**2 + a1*a2**2*a3*r2**2 + a1*a3*a4**2*r4**2 - a2 &
      *a3**2*a4*r3**2
f(707) = -a1*a2*a3*r2**3 + a1*a2*a4*r1**3 - a1*a3*a4*r4**3 + a2*a3*a4*r3 &
      **3
f(708) = a1*a2*b1**2*r1 - a1*a2*b1**2*r2 + a3*a4*b2**2*r3 - a3*a4*b2**2* &
      r4
f(709) = a1**2*a2*b1**2*r2 - a1*a2**2*b1**2*r1 + a3**2*a4*b2**2*r4 - a3* &
      a4**2*b2**2*r3
f(710) = -a1**2*a2*b1**2*r1 + a1*a2**2*b1**2*r2 - a3**2*a4*b2**2*r3 + a3 &
      *a4**2*b2**2*r4
f(711) = -a1*a2*b1**2*r1**2 + a1*a2*b1**2*r2**2 - a3*a4*b2**2*r3**2 + a3 &
      *a4*b2**2*r4**2
f(712) = a1*a2*b2**2*r1 - a1*a2*b2**2*r2 + a3*a4*b1**2*r3 - a3*a4*b1**2* &
      r4
f(713) = -a1**2*a2*b2**2*r2 + a1*a2**2*b2**2*r1 - a3**2*a4*b1**2*r4 + a3 &
      *a4**2*b1**2*r3
f(714) = a1**2*a2*b2**2*r1 - a1*a2**2*b2**2*r2 + a3**2*a4*b1**2*r3 - a3* &
      a4**2*b1**2*r4
f(715) = -a1*a2*b2**2*r1**2 + a1*a2*b2**2*r2**2 - a3*a4*b1**2*r3**2 + a3 &
      *a4*b1**2*r4**2
f(716) = dtau**2*(a1*a2*r1 - a1*a2*r2 + a3*a4*r3 - a3*a4*r4)
f(717) = dtau**2*(a1**2*a2*r2 - a1*a2**2*r1 + a3**2*a4*r4 - a3*a4**2*r3)
f(718) = dtau**2*(-a1**2*a2*r1 + a1*a2**2*r2 - a3**2*a4*r3 + a3*a4**2*r4 &
      )
f(719) = dtau**2*(a1*a2*r1**2 - a1*a2*r2**2 + a3*a4*r3**2 - a3*a4*r4**2)
f(720) = a1*a2*a3*r3 - a1*a2*a4*r4 + a1*a3*a4*r1 - a2*a3*a4*r2
f(721) = a1**2*a2*a4*r4 - a1*a2**2*a3*r3 - a1*a3*a4**2*r1 + a2*a3**2*a4* &
      r2
f(722) = -a1**3*a2*a4*r4 + a1*a2**3*a3*r3 + a1*a3*a4**3*r1 - a2*a3**3*a4 &
      *r2
f(723) = a1**2*a2*a3*r3 - a1*a2**2*a4*r4 + a1*a3**2*a4*r1 - a2*a3*a4**2* &
      r2
f(724) = a1**2*a2**2*a3*r3 - a1**2*a2**2*a4*r4 + a1*a3**2*a4**2*r1 - a2* &
      a3**2*a4**2*r2
f(725) = a1**3*a2*a3*r3 - a1*a2**3*a4*r4 + a1*a3**3*a4*r1 - a2*a3*a4**3* &
      r2
f(726) = a1**2*a3*a4*r1 + a1*a2*a3**2*r3 - a1*a2*a4**2*r4 - a2**2*a3*a4* &
      r2
f(727) = -a1**2*a2*a4**2*r4 + a1**2*a3*a4**2*r1 + a1*a2**2*a3**2*r3 - a2 &
      **2*a3**2*a4*r2
f(728) = -a1**2*a2*a3**2*r3 - a1**2*a3**2*a4*r1 + a1*a2**2*a4**2*r4 + a2 &
      **2*a3*a4**2*r2
f(729) = a1**3*a3*a4*r1 + a1*a2*a3**3*r3 - a1*a2*a4**3*r4 - a2**3*a3*a4* &
      r2
f(730) = a1*a2*a3*r3**2 - a1*a2*a4*r4**2 + a1*a3*a4*r1**2 - a2*a3*a4*r2 &
      **2
f(731) = -a1**2*a2*a4*r4**2 + a1*a2**2*a3*r3**2 + a1*a3*a4**2*r1**2 - a2 &
      *a3**2*a4*r2**2
f(732) = -a1**2*a2*a3*r3**2 + a1*a2**2*a4*r4**2 - a1*a3**2*a4*r1**2 + a2 &
      *a3*a4**2*r2**2
f(733) = -a1**2*a3*a4*r1**2 - a1*a2*a3**2*r3**2 + a1*a2*a4**2*r4**2 + a2 &
      **2*a3*a4*r2**2
f(734) = a1*a2*a3*r3**3 - a1*a2*a4*r4**3 + a1*a3*a4*r1**3 - a2*a3*a4*r2 &
      **3
f(735) = a1*a3*b1**2*r1 + a1*a3*b2**2*r3 - a2*a4*b1**2*r2 - a2*a4*b2**2* &
      r4
f(736) = a1**2*a3*b2**2*r3 + a1*a3**2*b1**2*r1 - a2**2*a4*b2**2*r4 - a2* &
      a4**2*b1**2*r2
f(737) = a1**2*a3*b1**2*r1 + a1*a3**2*b2**2*r3 - a2**2*a4*b1**2*r2 - a2* &
      a4**2*b2**2*r4
f(738) = -a1*a3*b1**2*r1**2 - a1*a3*b2**2*r3**2 + a2*a4*b1**2*r2**2 + a2 &
      *a4*b2**2*r4**2
f(739) = -a1*a3*b1**2*r3 - a1*a3*b2**2*r1 + a2*a4*b1**2*r4 + a2*a4*b2**2 &
      *r2
f(740) = a1**2*a3*b1**2*r3 + a1*a3**2*b2**2*r1 - a2**2*a4*b1**2*r4 - a2* &
      a4**2*b2**2*r2
f(741) = -a1**2*a3*b2**2*r1 - a1*a3**2*b1**2*r3 + a2**2*a4*b2**2*r2 + a2 &
      *a4**2*b1**2*r4
f(742) = a1*a3*b1**2*r3**2 + a1*a3*b2**2*r1**2 - a2*a4*b1**2*r4**2 - a2* &
      a4*b2**2*r2**2
f(743) = dtau**2*(a1*a3*r1 + a1*a3*r3 - a2*a4*r2 - a2*a4*r4)
f(744) = dtau**2*(-a1**2*a3*r3 - a1*a3**2*r1 + a2**2*a4*r4 + a2*a4**2*r2 &
      )
f(745) = dtau**2*(-a1**2*a3*r1 - a1*a3**2*r3 + a2**2*a4*r2 + a2*a4**2*r4 &
      )
f(746) = dtau**2*(a1*a3*r1**2 + a1*a3*r3**2 - a2*a4*r2**2 - a2*a4*r4**2)
f(747) = a1*a4*b1**2*r1 - a1*a4*b2**2*r4 - a2*a3*b1**2*r2 + a2*a3*b2**2* &
      r3
f(748) = -a1**2*a4*b2**2*r4 + a1*a4**2*b1**2*r1 + a2**2*a3*b2**2*r3 - a2 &
      *a3**2*b1**2*r2
f(749) = a1**2*a4*b1**2*r1 - a1*a4**2*b2**2*r4 - a2**2*a3*b1**2*r2 + a2* &
      a3**2*b2**2*r3
f(750) = a1*a4*b1**2*r1**2 - a1*a4*b2**2*r4**2 - a2*a3*b1**2*r2**2 + a2* &
      a3*b2**2*r3**2
f(751) = a1*a4*b1**2*r4 - a1*a4*b2**2*r1 - a2*a3*b1**2*r3 + a2*a3*b2**2* &
      r2
f(752) = -a1**2*a4*b1**2*r4 + a1*a4**2*b2**2*r1 + a2**2*a3*b1**2*r3 - a2 &
      *a3**2*b2**2*r2
f(753) = -a1**2*a4*b2**2*r1 + a1*a4**2*b1**2*r4 + a2**2*a3*b2**2*r2 - a2 &
      *a3**2*b1**2*r3
f(754) = a1*a4*b1**2*r4**2 - a1*a4*b2**2*r1**2 - a2*a3*b1**2*r3**2 + a2* &
      a3*b2**2*r2**2
f(755) = dtau**2*(a1*a4*r1 - a1*a4*r4 - a2*a3*r2 + a2*a3*r3)
f(756) = dtau**2*(-a1**2*a4*r4 + a1*a4**2*r1 + a2**2*a3*r3 - a2*a3**2*r2 &
      )
f(757) = dtau**2*(-a1**2*a4*r1 + a1*a4**2*r4 + a2**2*a3*r2 - a2*a3**2*r3 &
      )
f(758) = dtau**2*(a1*a4*r1**2 - a1*a4*r4**2 - a2*a3*r2**2 + a2*a3*r3**2)
f(759) = b1*b2*(a1*r1 - a2*r2 + a3*r3 - a4*r4)
f(760) = b1*b2*(-a1*b2**2*r1 + a2*b2**2*r2 - a3*b1**2*r3 + a4*b1**2*r4)
f(761) = b1**2*b2**2*(a1*r1 - a2*r2 + a3*r3 - a4*r4)
f(762) = b1*b2*(a1*b1**2*r1 - a2*b1**2*r2 + a3*b2**2*r3 - a4*b2**2*r4)
f(763) = b1*b2*(a1**2*r1 - a2**2*r2 + a3**2*r3 - a4**2*r4)
f(764) = b1*b2*(-a1**3*r1 + a2**3*r2 - a3**3*r3 + a4**3*r4)
f(765) = b1*b2*(-a1*r1**2 + a2*r2**2 - a3*r3**2 + a4*r4**2)
f(766) = b1*b2*(-a1**2*r1**2 + a2**2*r2**2 - a3**2*r3**2 + a4**2*r4**2)
f(767) = b1*b2*(a1*r1**3 - a2*r2**3 + a3*r3**3 - a4*r4**3)
f(768) = dtau*(a1*b1*r1 + a2*b1*r2 + a3*b2*r3 + a4*b2*r4)
f(769) = dtau**3*(a1*b1*r1 + a2*b1*r2 + a3*b2*r3 + a4*b2*r4)
f(770) = dtau**2*(-a1*b1**2*r1 + a2*b1**2*r2 - a3*b2**2*r3 + a4*b2**2*r4 &
      )
f(771) = dtau*(a1*b1**3*r1 + a2*b1**3*r2 + a3*b2**3*r3 + a4*b2**3*r4)
f(772) = dtau*(a1**2*b1*r1 + a2**2*b1*r2 + a3**2*b2*r3 + a4**2*b2*r4)
f(773) = dtau*(a1**3*b1*r1 + a2**3*b1*r2 + a3**3*b2*r3 + a4**3*b2*r4)
f(774) = dtau*(a1*b1*r1**2 + a2*b1*r2**2 + a3*b2*r3**2 + a4*b2*r4**2)
f(775) = dtau*(a1**2*b1*r1**2 + a2**2*b1*r2**2 + a3**2*b2*r3**2 + a4**2* &
      b2*r4**2)
f(776) = dtau*(a1*b1*r1**3 + a2*b1*r2**3 + a3*b2*r3**3 + a4*b2*r4**3)
f(777) = dtau*(a1*b2*r1 + a2*b2*r2 + a3*b1*r3 + a4*b1*r4)
f(778) = dtau**3*(a1*b2*r1 + a2*b2*r2 + a3*b1*r3 + a4*b1*r4)
f(779) = dtau**2*(-a1*b2**2*r1 + a2*b2**2*r2 - a3*b1**2*r3 + a4*b1**2*r4 &
      )
f(780) = dtau*(a1*b2**3*r1 + a2*b2**3*r2 + a3*b1**3*r3 + a4*b1**3*r4)
f(781) = dtau*(a1**2*b2*r1 + a2**2*b2*r2 + a3**2*b1*r3 + a4**2*b1*r4)
f(782) = dtau*(a1**3*b2*r1 + a2**3*b2*r2 + a3**3*b1*r3 + a4**3*b1*r4)
f(783) = dtau*(a1*b2*r1**2 + a2*b2*r2**2 + a3*b1*r3**2 + a4*b1*r4**2)
f(784) = dtau*(a1**2*b2*r1**2 + a2**2*b2*r2**2 + a3**2*b1*r3**2 + a4**2* &
      b1*r4**2)
f(785) = dtau*(a1*b2*r1**3 + a2*b2*r2**3 + a3*b1*r3**3 + a4*b1*r4**3)
f(786) = a1*a2*a3*r4 - a1*a2*a4*r3 + a1*a3*a4*r2 - a2*a3*a4*r1
f(787) = -a1**2*a2*a3*r4 + a1*a2**2*a4*r3 - a1*a3**2*a4*r2 + a2*a3*a4**2 &
      *r1
f(788) = -a1**3*a2*a3*r4 + a1*a2**3*a4*r3 - a1*a3**3*a4*r2 + a2*a3*a4**3 &
      *r1
f(789) = -a1**2*a2*a4*r3 + a1*a2**2*a3*r4 + a1*a3*a4**2*r2 - a2*a3**2*a4 &
      *r1
f(790) = -a1**2*a2**2*a3*r4 + a1**2*a2**2*a4*r3 - a1*a3**2*a4**2*r2 + a2 &
      *a3**2*a4**2*r1
f(791) = a1**3*a2*a4*r3 - a1*a2**3*a3*r4 - a1*a3*a4**3*r2 + a2*a3**3*a4* &
      r1
f(792) = a1**2*a3*a4*r2 + a1*a2*a3**2*r4 - a1*a2*a4**2*r3 - a2**2*a3*a4* &
      r1
f(793) = a1**2*a2*a3**2*r4 + a1**2*a3**2*a4*r2 - a1*a2**2*a4**2*r3 - a2 &
      **2*a3*a4**2*r1
f(794) = a1**2*a2*a4**2*r3 - a1**2*a3*a4**2*r2 - a1*a2**2*a3**2*r4 + a2 &
      **2*a3**2*a4*r1
f(795) = a1**3*a3*a4*r2 + a1*a2*a3**3*r4 - a1*a2*a4**3*r3 - a2**3*a3*a4* &
      r1
f(796) = -a1*a2*a3*r4**2 + a1*a2*a4*r3**2 - a1*a3*a4*r2**2 + a2*a3*a4*r1 &
      **2
f(797) = a1**2*a2*a3*r4**2 - a1*a2**2*a4*r3**2 + a1*a3**2*a4*r2**2 - a2* &
      a3*a4**2*r1**2
f(798) = a1**2*a2*a4*r3**2 - a1*a2**2*a3*r4**2 - a1*a3*a4**2*r2**2 + a2* &
      a3**2*a4*r1**2
f(799) = -a1**2*a3*a4*r2**2 - a1*a2*a3**2*r4**2 + a1*a2*a4**2*r3**2 + a2 &
      **2*a3*a4*r1**2
f(800) = -a1*a2*a3*r4**3 + a1*a2*a4*r3**3 - a1*a3*a4*r2**3 + a2*a3*a4*r1 &
      **3
f(801) = -a1*a4*b1**2*r2 + a1*a4*b2**2*r3 + a2*a3*b1**2*r1 - a2*a3*b2**2 &
      *r4
f(802) = -a1**2*a4*b2**2*r3 + a1*a4**2*b1**2*r2 + a2**2*a3*b2**2*r4 - a2 &
      *a3**2*b1**2*r1
f(803) = -a1**2*a4*b1**2*r2 + a1*a4**2*b2**2*r3 + a2**2*a3*b1**2*r1 - a2 &
      *a3**2*b2**2*r4
f(804) = a1*a4*b1**2*r2**2 - a1*a4*b2**2*r3**2 - a2*a3*b1**2*r1**2 + a2* &
      a3*b2**2*r4**2
f(805) = -a1*a4*b1**2*r3 + a1*a4*b2**2*r2 + a2*a3*b1**2*r4 - a2*a3*b2**2 &
      *r1
f(806) = -a1**2*a4*b1**2*r3 + a1*a4**2*b2**2*r2 + a2**2*a3*b1**2*r4 - a2 &
      *a3**2*b2**2*r1
f(807) = -a1**2*a4*b2**2*r2 + a1*a4**2*b1**2*r3 + a2**2*a3*b2**2*r1 - a2 &
      *a3**2*b1**2*r4
f(808) = -a1*a4*b1**2*r3**2 + a1*a4*b2**2*r2**2 + a2*a3*b1**2*r4**2 - a2 &
      *a3*b2**2*r1**2
f(809) = dtau**2*(-a1*a4*r2 + a1*a4*r3 + a2*a3*r1 - a2*a3*r4)
f(810) = dtau**2*(a1**2*a4*r3 - a1*a4**2*r2 - a2**2*a3*r4 + a2*a3**2*r1)
f(811) = dtau**2*(-a1**2*a4*r2 + a1*a4**2*r3 + a2**2*a3*r1 - a2*a3**2*r4 &
      )
f(812) = dtau**2*(-a1*a4*r2**2 + a1*a4*r3**2 + a2*a3*r1**2 - a2*a3*r4**2 &
      )
f(813) = a1*a3*b1**2*r2 + a1*a3*b2**2*r4 - a2*a4*b1**2*r1 - a2*a4*b2**2* &
      r3
f(814) = a1**2*a3*b2**2*r4 + a1*a3**2*b1**2*r2 - a2**2*a4*b2**2*r3 - a2* &
      a4**2*b1**2*r1
f(815) = a1**2*a3*b1**2*r2 + a1*a3**2*b2**2*r4 - a2**2*a4*b1**2*r1 - a2* &
      a4**2*b2**2*r3
f(816) = a1*a3*b1**2*r2**2 + a1*a3*b2**2*r4**2 - a2*a4*b1**2*r1**2 - a2* &
      a4*b2**2*r3**2
f(817) = a1*a3*b1**2*r4 + a1*a3*b2**2*r2 - a2*a4*b1**2*r3 - a2*a4*b2**2* &
      r1
f(818) = a1**2*a3*b1**2*r4 + a1*a3**2*b2**2*r2 - a2**2*a4*b1**2*r3 - a2* &
      a4**2*b2**2*r1
f(819) = -a1**2*a3*b2**2*r2 - a1*a3**2*b1**2*r4 + a2**2*a4*b2**2*r1 + a2 &
      *a4**2*b1**2*r3
f(820) = a1*a3*b1**2*r4**2 + a1*a3*b2**2*r2**2 - a2*a4*b1**2*r3**2 - a2* &
      a4*b2**2*r1**2
f(821) = dtau**2*(-a1*a3*r2 - a1*a3*r4 + a2*a4*r1 + a2*a4*r3)
f(822) = dtau**2*(a1**2*a3*r4 + a1*a3**2*r2 - a2**2*a4*r3 - a2*a4**2*r1)
f(823) = dtau**2*(-a1**2*a3*r2 - a1*a3**2*r4 + a2**2*a4*r1 + a2*a4**2*r3 &
      )
f(824) = dtau**2*(-a1*a3*r2**2 - a1*a3*r4**2 + a2*a4*r1**2 + a2*a4*r3**2 &
      )
f(825) = b1*b2*(-a1*r2 + a2*r1 - a3*r4 + a4*r3)
f(826) = b1*b2*(-a1*b2**2*r2 + a2*b2**2*r1 - a3*b1**2*r4 + a4*b1**2*r3)
f(827) = b1**2*b2**2*(-a1*r2 + a2*r1 - a3*r4 + a4*r3)
f(828) = b1*b2*(-a1*b1**2*r2 + a2*b1**2*r1 - a3*b2**2*r4 + a4*b2**2*r3)
f(829) = b1*b2*(-a1**2*r2 + a2**2*r1 - a3**2*r4 + a4**2*r3)
f(830) = b1*b2*(a1**3*r2 - a2**3*r1 + a3**3*r4 - a4**3*r3)
f(831) = b1*b2*(a1*r2**2 - a2*r1**2 + a3*r4**2 - a4*r3**2)
f(832) = b1*b2*(-a1**2*r2**2 + a2**2*r1**2 - a3**2*r4**2 + a4**2*r3**2)
f(833) = b1*b2*(-a1*r2**3 + a2*r1**3 - a3*r4**3 + a4*r3**3)
f(834) = dtau*(a1*b1*r2 + a2*b1*r1 + a3*b2*r4 + a4*b2*r3)
f(835) = dtau**3*(a1*b1*r2 + a2*b1*r1 + a3*b2*r4 + a4*b2*r3)
f(836) = dtau**2*(a1*b1**2*r2 - a2*b1**2*r1 + a3*b2**2*r4 - a4*b2**2*r3)
f(837) = dtau*(a1*b1**3*r2 + a2*b1**3*r1 + a3*b2**3*r4 + a4*b2**3*r3)
f(838) = dtau*(a1**2*b1*r2 + a2**2*b1*r1 + a3**2*b2*r4 + a4**2*b2*r3)
f(839) = dtau*(a1**3*b1*r2 + a2**3*b1*r1 + a3**3*b2*r4 + a4**3*b2*r3)
f(840) = dtau*(a1*b1*r2**2 + a2*b1*r1**2 + a3*b2*r4**2 + a4*b2*r3**2)
f(841) = dtau*(a1**2*b1*r2**2 + a2**2*b1*r1**2 + a3**2*b2*r4**2 + a4**2* &
      b2*r3**2)
f(842) = dtau*(a1*b1*r2**3 + a2*b1*r1**3 + a3*b2*r4**3 + a4*b2*r3**3)
f(843) = dtau*(a1*b2*r2 + a2*b2*r1 + a3*b1*r4 + a4*b1*r3)
f(844) = dtau**3*(a1*b2*r2 + a2*b2*r1 + a3*b1*r4 + a4*b1*r3)
f(845) = dtau**2*(a1*b2**2*r2 - a2*b2**2*r1 + a3*b1**2*r4 - a4*b1**2*r3)
f(846) = dtau*(a1*b2**3*r2 + a2*b2**3*r1 + a3*b1**3*r4 + a4*b1**3*r3)
f(847) = dtau*(a1**2*b2*r2 + a2**2*b2*r1 + a3**2*b1*r4 + a4**2*b1*r3)
f(848) = dtau*(a1**3*b2*r2 + a2**3*b2*r1 + a3**3*b1*r4 + a4**3*b1*r3)
f(849) = dtau*(a1*b2*r2**2 + a2*b2*r1**2 + a3*b1*r4**2 + a4*b1*r3**2)
f(850) = dtau*(a1**2*b2*r2**2 + a2**2*b2*r1**2 + a3**2*b1*r4**2 + a4**2* &
      b1*r3**2)
f(851) = dtau*(a1*b2*r2**3 + a2*b2*r1**3 + a3*b1*r4**3 + a4*b1*r3**3)
f(852) = -a1*a2*b2**2*r3 + a1*a2*b2**2*r4 - a3*a4*b1**2*r1 + a3*a4*b1**2 &
      *r2
f(853) = a1**2*a2*b2**2*r4 - a1*a2**2*b2**2*r3 + a3**2*a4*b1**2*r2 - a3* &
      a4**2*b1**2*r1
f(854) = -a1**2*a2*b2**2*r3 + a1*a2**2*b2**2*r4 - a3**2*a4*b1**2*r1 + a3 &
      *a4**2*b1**2*r2
f(855) = a1*a2*b2**2*r3**2 - a1*a2*b2**2*r4**2 + a3*a4*b1**2*r1**2 - a3* &
      a4*b1**2*r2**2
f(856) = a1*a2*b1**2*r3 - a1*a2*b1**2*r4 + a3*a4*b2**2*r1 - a3*a4*b2**2* &
      r2
f(857) = a1**2*a2*b1**2*r4 - a1*a2**2*b1**2*r3 + a3**2*a4*b2**2*r2 - a3* &
      a4**2*b2**2*r1
f(858) = -a1**2*a2*b1**2*r3 + a1*a2**2*b1**2*r4 - a3**2*a4*b2**2*r1 + a3 &
      *a4**2*b2**2*r2
f(859) = -a1*a2*b1**2*r3**2 + a1*a2*b1**2*r4**2 - a3*a4*b2**2*r1**2 + a3 &
      *a4*b2**2*r2**2
f(860) = dtau**2*(-a1*a2*r3 + a1*a2*r4 - a3*a4*r1 + a3*a4*r2)
f(861) = dtau**2*(-a1**2*a2*r4 + a1*a2**2*r3 - a3**2*a4*r2 + a3*a4**2*r1 &
      )
f(862) = dtau**2*(a1**2*a2*r3 - a1*a2**2*r4 + a3**2*a4*r1 - a3*a4**2*r2)
f(863) = dtau**2*(-a1*a2*r3**2 + a1*a2*r4**2 - a3*a4*r1**2 + a3*a4*r2**2 &
      )
f(864) = b1*b2*(-a1*r3 + a2*r4 - a3*r1 + a4*r2)
f(865) = b1*b2*(-a1*b1**2*r3 + a2*b1**2*r4 - a3*b2**2*r1 + a4*b2**2*r2)
f(866) = b1**2*b2**2*(-a1*r3 + a2*r4 - a3*r1 + a4*r2)
f(867) = b1*b2*(-a1*b2**2*r3 + a2*b2**2*r4 - a3*b1**2*r1 + a4*b1**2*r2)
f(868) = b1*b2*(a1**2*r3 - a2**2*r4 + a3**2*r1 - a4**2*r2)
f(869) = b1*b2*(-a1**3*r3 + a2**3*r4 - a3**3*r1 + a4**3*r2)
f(870) = b1*b2*(a1*r3**2 - a2*r4**2 + a3*r1**2 - a4*r2**2)
f(871) = b1*b2*(-a1**2*r3**2 + a2**2*r4**2 - a3**2*r1**2 + a4**2*r2**2)
f(872) = b1*b2*(a1*r3**3 - a2*r4**3 + a3*r1**3 - a4*r2**3)
f(873) = dtau*(a1*b2*r3 + a2*b2*r4 + a3*b1*r1 + a4*b1*r2)
f(874) = dtau**3*(a1*b2*r3 + a2*b2*r4 + a3*b1*r1 + a4*b1*r2)
f(875) = dtau**2*(a1*b2**2*r3 - a2*b2**2*r4 + a3*b1**2*r1 - a4*b1**2*r2)
f(876) = dtau*(a1*b2**3*r3 + a2*b2**3*r4 + a3*b1**3*r1 + a4*b1**3*r2)
f(877) = dtau*(a1**2*b2*r3 + a2**2*b2*r4 + a3**2*b1*r1 + a4**2*b1*r2)
f(878) = dtau*(a1**3*b2*r3 + a2**3*b2*r4 + a3**3*b1*r1 + a4**3*b1*r2)
f(879) = dtau*(a1*b2*r3**2 + a2*b2*r4**2 + a3*b1*r1**2 + a4*b1*r2**2)
f(880) = dtau*(a1**2*b2*r3**2 + a2**2*b2*r4**2 + a3**2*b1*r1**2 + a4**2* &
      b1*r2**2)
f(881) = dtau*(a1*b2*r3**3 + a2*b2*r4**3 + a3*b1*r1**3 + a4*b1*r2**3)
f(882) = dtau*(a1*b1*r3 + a2*b1*r4 + a3*b2*r1 + a4*b2*r2)
f(883) = dtau**3*(a1*b1*r3 + a2*b1*r4 + a3*b2*r1 + a4*b2*r2)
f(884) = dtau**2*(-a1*b1**2*r3 + a2*b1**2*r4 - a3*b2**2*r1 + a4*b2**2*r2 &
      )
f(885) = dtau*(a1*b1**3*r3 + a2*b1**3*r4 + a3*b2**3*r1 + a4*b2**3*r2)
f(886) = dtau*(a1**2*b1*r3 + a2**2*b1*r4 + a3**2*b2*r1 + a4**2*b2*r2)
f(887) = dtau*(a1**3*b1*r3 + a2**3*b1*r4 + a3**3*b2*r1 + a4**3*b2*r2)
f(888) = dtau*(a1*b1*r3**2 + a2*b1*r4**2 + a3*b2*r1**2 + a4*b2*r2**2)
f(889) = dtau*(a1**2*b1*r3**2 + a2**2*b1*r4**2 + a3**2*b2*r1**2 + a4**2* &
      b2*r2**2)
f(890) = dtau*(a1*b1*r3**3 + a2*b1*r4**3 + a3*b2*r1**3 + a4*b2*r2**3)
f(891) = b1*b2*(a1*r4 - a2*r3 + a3*r2 - a4*r1)
f(892) = b1*b2*(a1*b1**2*r4 - a2*b1**2*r3 + a3*b2**2*r2 - a4*b2**2*r1)
f(893) = b1**2*b2**2*(a1*r4 - a2*r3 + a3*r2 - a4*r1)
f(894) = b1*b2*(a1*b2**2*r4 - a2*b2**2*r3 + a3*b1**2*r2 - a4*b1**2*r1)
f(895) = b1*b2*(a1**2*r4 - a2**2*r3 + a3**2*r2 - a4**2*r1)
f(896) = b1*b2*(a1**3*r4 - a2**3*r3 + a3**3*r2 - a4**3*r1)
f(897) = b1*b2*(-a1*r4**2 + a2*r3**2 - a3*r2**2 + a4*r1**2)
f(898) = b1*b2*(a1**2*r4**2 - a2**2*r3**2 + a3**2*r2**2 - a4**2*r1**2)
f(899) = b1*b2*(a1*r4**3 - a2*r3**3 + a3*r2**3 - a4*r1**3)
f(900) = dtau*(a1*b2*r4 + a2*b2*r3 + a3*b1*r2 + a4*b1*r1)
f(901) = dtau**3*(a1*b2*r4 + a2*b2*r3 + a3*b1*r2 + a4*b1*r1)
f(902) = dtau**2*(-a1*b2**2*r4 + a2*b2**2*r3 - a3*b1**2*r2 + a4*b1**2*r1 &
      )
f(903) = dtau*(a1*b2**3*r4 + a2*b2**3*r3 + a3*b1**3*r2 + a4*b1**3*r1)
f(904) = dtau*(a1**2*b2*r4 + a2**2*b2*r3 + a3**2*b1*r2 + a4**2*b1*r1)
f(905) = dtau*(a1**3*b2*r4 + a2**3*b2*r3 + a3**3*b1*r2 + a4**3*b1*r1)
f(906) = dtau*(a1*b2*r4**2 + a2*b2*r3**2 + a3*b1*r2**2 + a4*b1*r1**2)
f(907) = dtau*(a1**2*b2*r4**2 + a2**2*b2*r3**2 + a3**2*b1*r2**2 + a4**2* &
      b1*r1**2)
f(908) = dtau*(a1*b2*r4**3 + a2*b2*r3**3 + a3*b1*r2**3 + a4*b1*r1**3)
f(909) = dtau*(a1*b1*r4 + a2*b1*r3 + a3*b2*r2 + a4*b2*r1)
f(910) = dtau**3*(a1*b1*r4 + a2*b1*r3 + a3*b2*r2 + a4*b2*r1)
f(911) = dtau**2*(-a1*b1**2*r4 + a2*b1**2*r3 - a3*b2**2*r2 + a4*b2**2*r1 &
      )
f(912) = dtau*(a1*b1**3*r4 + a2*b1**3*r3 + a3*b2**3*r2 + a4*b2**3*r1)
f(913) = dtau*(a1**2*b1*r4 + a2**2*b1*r3 + a3**2*b2*r2 + a4**2*b2*r1)
f(914) = dtau*(a1**3*b1*r4 + a2**3*b1*r3 + a3**3*b2*r2 + a4**3*b2*r1)
f(915) = dtau*(a1*b1*r4**2 + a2*b1*r3**2 + a3*b2*r2**2 + a4*b2*r1**2)
f(916) = dtau*(a1**2*b1*r4**2 + a2**2*b1*r3**2 + a3**2*b2*r2**2 + a4**2* &
      b2*r1**2)
f(917) = dtau*(a1*b1*r4**3 + a2*b1*r3**3 + a3*b2*r2**3 + a4*b2*r1**3)
f(918) = b1*b2*dtau**2*(r1 - r2 + r3 - r4)
f(919) = b1*b2*dtau*(b1*r3 + b1*r4 + b2*r1 + b2*r2)
f(920) = b1*b2*dtau*(b1*r1 + b1*r2 + b2*r3 + b2*r4)
f(921) = b1*b2*dtau**2*(r1**2 - r2**2 + r3**2 - r4**2)
f(922) = b1*b2*dtau*(b1*r3**2 + b1*r4**2 + b2*r1**2 + b2*r2**2)
f(923) = b1*b2*dtau*(b1*r1**2 + b1*r2**2 + b2*r3**2 + b2*r4**2)
f(924) = a1*a2*a3*a4*(-a1 + a2 - a3 + a4)
f(925) = a1*a2*a3*a4*(a1**2 - a2**2 + a3**2 - a4**2)
f(926) = a1*a2*a3*a4*(a1*a3 - a2*a4)
f(927) = a1*a2*a3*b1**2 - a1*a2*a4*b1**2 + a1*a3*a4*b2**2 - a2*a3*a4*b2 &
      **2
f(928) = a1**2*a3*a4*b2**2 + a1*a2*a3**2*b1**2 - a1*a2*a4**2*b1**2 - a2 &
      **2*a3*a4*b2**2
f(929) = -a1**2*a2*a4*b1**2 + a1*a2**2*a3*b1**2 + a1*a3*a4**2*b2**2 - a2 &
      *a3**2*a4*b2**2
f(930) = -a1**2*a2*a3*b1**2 + a1*a2**2*a4*b1**2 - a1*a3**2*a4*b2**2 + a2 &
      *a3*a4**2*b2**2
f(931) = a1*a2*a3*b2**2 - a1*a2*a4*b2**2 + a1*a3*a4*b1**2 - a2*a3*a4*b1 &
      **2
f(932) = a1**2*a3*a4*b1**2 + a1*a2*a3**2*b2**2 - a1*a2*a4**2*b2**2 - a2 &
      **2*a3*a4*b1**2
f(933) = a1**2*a2*a4*b2**2 - a1*a2**2*a3*b2**2 - a1*a3*a4**2*b1**2 + a2* &
      a3**2*a4*b1**2
f(934) = a1**2*a2*a3*b2**2 - a1*a2**2*a4*b2**2 + a1*a3**2*a4*b1**2 - a2* &
      a3*a4**2*b1**2
f(935) = dtau**2*(-a1*a2*a3 + a1*a2*a4 - a1*a3*a4 + a2*a3*a4)
f(936) = dtau**2*(-a1**2*a3*a4 - a1*a2*a3**2 + a1*a2*a4**2 + a2**2*a3*a4 &
      )
f(937) = dtau**2*(a1**2*a2*a4 - a1*a2**2*a3 - a1*a3*a4**2 + a2*a3**2*a4)
f(938) = dtau**2*(a1**2*a2*a3 - a1*a2**2*a4 + a1*a3**2*a4 - a2*a3*a4**2)
f(939) = b1*b2*(-a1**2*a2 + a1*a2**2 - a3**2*a4 + a3*a4**2)
f(940) = b1*b2*(a1**3*a2 - a1*a2**3 + a3**3*a4 - a3*a4**3)
f(941) = dtau*(a1*a2*b1 + a3*a4*b2)
f(942) = dtau**3*(a1*a2*b1 + a3*a4*b2)
f(943) = dtau*(a1*a2*b1**3 + a3*a4*b2**3)
f(944) = dtau*(a1**2*a2*b1 + a1*a2**2*b1 + a3**2*a4*b2 + a3*a4**2*b2)
f(945) = dtau*(a1**3*a2*b1 + a1*a2**3*b1 + a3**3*a4*b2 + a3*a4**3*b2)
f(946) = dtau*(a1**2*a2**2*b1 + a3**2*a4**2*b2)
f(947) = dtau*(a1*a2*b2 + a3*a4*b1)
f(948) = dtau**3*(a1*a2*b2 + a3*a4*b1)
f(949) = dtau*(a1*a2*b2**3 + a3*a4*b1**3)
f(950) = dtau*(a1**2*a2*b2 + a1*a2**2*b2 + a3**2*a4*b1 + a3*a4**2*b1)
f(951) = dtau*(a1**3*a2*b2 + a1*a2**3*b2 + a3**3*a4*b1 + a3*a4**3*b1)
f(952) = dtau*(a1**2*a2**2*b2 + a3**2*a4**2*b1)
f(953) = b1*b2*(-a1*a3 + a2*a4)
f(954) = b1*b2*(-a1*a3*b1**2 - a1*a3*b2**2 + a2*a4*b1**2 + a2*a4*b2**2)
f(955) = b1**2*b2**2*(-a1*a3 + a2*a4)
f(956) = b1*b2*(-a1**2*a3 - a1*a3**2 + a2**2*a4 + a2*a4**2)
f(957) = b1*b2*(-a1**3*a3 - a1*a3**3 + a2**3*a4 + a2*a4**3)
f(958) = b1*b2*(-a1**2*a3**2 + a2**2*a4**2)
f(959) = dtau*(a1*a3*b1 + a1*a3*b2 + a2*a4*b1 + a2*a4*b2)
f(960) = dtau**3*(a1*a3*b1 + a1*a3*b2 + a2*a4*b1 + a2*a4*b2)
f(961) = dtau**2*(a1*a3*b1**2 + a1*a3*b2**2 - a2*a4*b1**2 - a2*a4*b2**2)
f(962) = dtau*(a1*a3*b1**3 + a1*a3*b2**3 + a2*a4*b1**3 + a2*a4*b2**3)
f(963) = dtau*(a1**2*a3*b2 + a1*a3**2*b1 + a2**2*a4*b2 + a2*a4**2*b1)
f(964) = dtau*(a1**3*a3*b2 + a1*a3**3*b1 + a2**3*a4*b2 + a2*a4**3*b1)
f(965) = dtau*(a1**2*a3*b1 + a1*a3**2*b2 + a2**2*a4*b1 + a2*a4**2*b2)
f(966) = dtau*(a1**2*a3**2*b1 + a1**2*a3**2*b2 + a2**2*a4**2*b1 + a2**2* &
      a4**2*b2)
f(967) = dtau*(a1**3*a3*b1 + a1*a3**3*b2 + a2**3*a4*b1 + a2*a4**3*b2)
f(968) = b1*b2*(a1*a4*b1**2 - a1*a4*b2**2 - a2*a3*b1**2 + a2*a3*b2**2)
f(969) = b1*b2*(a1**2*a4 - a1*a4**2 - a2**2*a3 + a2*a3**2)
f(970) = b1*b2*(a1**3*a4 - a1*a4**3 - a2**3*a3 + a2*a3**3)
f(971) = dtau*(a1*a4*b1 + a1*a4*b2 + a2*a3*b1 + a2*a3*b2)
f(972) = dtau**3*(a1*a4*b1 + a1*a4*b2 + a2*a3*b1 + a2*a3*b2)
f(973) = dtau**2*(a1*a4*b1**2 - a1*a4*b2**2 - a2*a3*b1**2 + a2*a3*b2**2)
f(974) = dtau*(a1*a4*b1**3 + a1*a4*b2**3 + a2*a3*b1**3 + a2*a3*b2**3)
f(975) = dtau*(a1**2*a4*b2 + a1*a4**2*b1 + a2**2*a3*b2 + a2*a3**2*b1)
f(976) = dtau*(a1**3*a4*b2 + a1*a4**3*b1 + a2**3*a3*b2 + a2*a3**3*b1)
f(977) = dtau*(a1**2*a4*b1 + a1*a4**2*b2 + a2**2*a3*b1 + a2*a3**2*b2)
f(978) = dtau*(a1**2*a4**2*b1 + a1**2*a4**2*b2 + a2**2*a3**2*b1 + a2**2* &
      a3**2*b2)
f(979) = dtau*(a1**3*a4*b1 + a1*a4**3*b2 + a2**3*a3*b1 + a2*a3**3*b2)
f(980) = b1*b2*dtau**2*(-a1 + a2 - a3 + a4)
f(981) = b1*b2*dtau*(a1*b2 + a2*b2 + a3*b1 + a4*b1)
f(982) = b1*b2*dtau*(a1*b1 + a2*b1 + a3*b2 + a4*b2)
f(983) = b1*b2*dtau**2*(a1**2 - a2**2 + a3**2 - a4**2)
f(984) = b1*b2*dtau*(a1**2*b2 + a2**2*b2 + a3**2*b1 + a4**2*b1)
f(985) = b1*b2*dtau*(a1**2*b1 + a2**2*b1 + a3**2*b2 + a4**2*b2)
v = sum(f*params)
end function c2h4_dipole_b2u_n4_d6_ADF


!###############################################################################


function c2h4_dipole_b3u_n1_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(3)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(3)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = -b1 + b2
f(2) = b1**3 - b2**3
f(3) = b1**5 - b2**5
v = sum(f*params)
end function c2h4_dipole_b3u_n1_d6_ADF


!###############################################################################


function c2h4_dipole_b3u_n2_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(69)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(69)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(b1 - b2)
f(2) = r0*(b1**3 - b2**3)
f(3) = r0*(-b1**5 + b2**5)
f(4) = r0**2*(-b1 + b2)
f(5) = r0**2*(b1**3 - b2**3)
f(6) = r0**3*(-b1 + b2)
f(7) = r0**3*(-b1**3 + b2**3)
f(8) = r0**4*(-b1 + b2)
f(9) = r0**5*(b1 - b2)
f(10) = -b1*r1 - b1*r2 + b2*r3 + b2*r4
f(11) = -b1**3*r1 - b1**3*r2 + b2**3*r3 + b2**3*r4
f(12) = b1**5*r1 + b1**5*r2 - b2**5*r3 - b2**5*r4
f(13) = b1*r1**2 + b1*r2**2 - b2*r3**2 - b2*r4**2
f(14) = -b1**3*r1**2 - b1**3*r2**2 + b2**3*r3**2 + b2**3*r4**2
f(15) = b1*r1**3 + b1*r2**3 - b2*r3**3 - b2*r4**3
f(16) = b1**3*r1**3 + b1**3*r2**3 - b2**3*r3**3 - b2**3*r4**3
f(17) = b1*r1**4 + b1*r2**4 - b2*r3**4 - b2*r4**4
f(18) = -b1*r1**5 - b1*r2**5 + b2*r3**5 + b2*r4**5
f(19) = b1*r3 + b1*r4 - b2*r1 - b2*r2
f(20) = b1**3*r3 + b1**3*r4 - b2**3*r1 - b2**3*r2
f(21) = -b1**5*r3 - b1**5*r4 + b2**5*r1 + b2**5*r2
f(22) = -b1*r3**2 - b1*r4**2 + b2*r1**2 + b2*r2**2
f(23) = b1**3*r3**2 + b1**3*r4**2 - b2**3*r1**2 - b2**3*r2**2
f(24) = b1*r3**3 + b1*r4**3 - b2*r1**3 - b2*r2**3
f(25) = -b1**3*r3**3 - b1**3*r4**3 + b2**3*r1**3 + b2**3*r2**3
f(26) = b1*r3**4 + b1*r4**4 - b2*r1**4 - b2*r2**4
f(27) = b1*r3**5 + b1*r4**5 - b2*r1**5 - b2*r2**5
f(28) = dtau*(-r1 + r2 + r3 - r4)
f(29) = dtau**3*(r1 - r2 - r3 + r4)
f(30) = dtau**5*(-r1 + r2 + r3 - r4)
f(31) = dtau*(r1**2 - r2**2 - r3**2 + r4**2)
f(32) = dtau**3*(r1**2 - r2**2 - r3**2 + r4**2)
f(33) = dtau*(-r1**3 + r2**3 + r3**3 - r4**3)
f(34) = dtau**3*(r1**3 - r2**3 - r3**3 + r4**3)
f(35) = dtau*(-r1**4 + r2**4 + r3**4 - r4**4)
f(36) = dtau*(r1**5 - r2**5 - r3**5 + r4**5)
f(37) = a1*b1 + a2*b1 - a3*b2 - a4*b2
f(38) = a1*b1**3 + a2*b1**3 - a3*b2**3 - a4*b2**3
f(39) = a1*b1**5 + a2*b1**5 - a3*b2**5 - a4*b2**5
f(40) = a1**2*b1 + a2**2*b1 - a3**2*b2 - a4**2*b2
f(41) = -a1**2*b1**3 - a2**2*b1**3 + a3**2*b2**3 + a4**2*b2**3
f(42) = a1**3*b1 + a2**3*b1 - a3**3*b2 - a4**3*b2
f(43) = a1**3*b1**3 + a2**3*b1**3 - a3**3*b2**3 - a4**3*b2**3
f(44) = a1**4*b1 + a2**4*b1 - a3**4*b2 - a4**4*b2
f(45) = a1**5*b1 + a2**5*b1 - a3**5*b2 - a4**5*b2
f(46) = a1*b2 + a2*b2 - a3*b1 - a4*b1
f(47) = a1*b2**3 + a2*b2**3 - a3*b1**3 - a4*b1**3
f(48) = a1*b2**5 + a2*b2**5 - a3*b1**5 - a4*b1**5
f(49) = -a1**2*b2 - a2**2*b2 + a3**2*b1 + a4**2*b1
f(50) = a1**2*b2**3 + a2**2*b2**3 - a3**2*b1**3 - a4**2*b1**3
f(51) = -a1**3*b2 - a2**3*b2 + a3**3*b1 + a4**3*b1
f(52) = -a1**3*b2**3 - a2**3*b2**3 + a3**3*b1**3 + a4**3*b1**3
f(53) = -a1**4*b2 - a2**4*b2 + a3**4*b1 + a4**4*b1
f(54) = -a1**5*b2 - a2**5*b2 + a3**5*b1 + a4**5*b1
f(55) = dtau*(-a1 + a2 + a3 - a4)
f(56) = dtau**3*(a1 - a2 - a3 + a4)
f(57) = dtau**5*(-a1 + a2 + a3 - a4)
f(58) = dtau*(-a1**2 + a2**2 + a3**2 - a4**2)
f(59) = dtau**3*(-a1**2 + a2**2 + a3**2 - a4**2)
f(60) = dtau*(a1**3 - a2**3 - a3**3 + a4**3)
f(61) = dtau**3*(a1**3 - a2**3 - a3**3 + a4**3)
f(62) = dtau*(a1**4 - a2**4 - a3**4 + a4**4)
f(63) = dtau*(-a1**5 + a2**5 + a3**5 - a4**5)
f(64) = b1*b2*(b1 - b2)
f(65) = b1*b2*(b1**3 - b2**3)
f(66) = b1**2*b2**2*(b1 - b2)
f(67) = dtau**2*(b1 - b2)
f(68) = dtau**4*(-b1 + b2)
f(69) = dtau**2*(b1**3 - b2**3)
v = sum(f*params)
end function c2h4_dipole_b3u_n2_d6_ADF


!###############################################################################


function c2h4_dipole_b3u_n3_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(433)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(433)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(2) = r0*(b1**3*r1 + b1**3*r2 - b2**3*r3 - b2**3*r4)
f(3) = r0*(b1*r1**2 + b1*r2**2 - b2*r3**2 - b2*r4**2)
f(4) = r0*(-b1**3*r1**2 - b1**3*r2**2 + b2**3*r3**2 + b2**3*r4**2)
f(5) = r0*(b1*r1**3 + b1*r2**3 - b2*r3**3 - b2*r4**3)
f(6) = r0*(-b1*r1**4 - b1*r2**4 + b2*r3**4 + b2*r4**4)
f(7) = r0**2*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(8) = r0**2*(-b1**3*r1 - b1**3*r2 + b2**3*r3 + b2**3*r4)
f(9) = r0**2*(b1*r1**2 + b1*r2**2 - b2*r3**2 - b2*r4**2)
f(10) = r0**2*(-b1*r1**3 - b1*r2**3 + b2*r3**3 + b2*r4**3)
f(11) = r0**3*(-b1*r1 - b1*r2 + b2*r3 + b2*r4)
f(12) = r0**3*(b1*r1**2 + b1*r2**2 - b2*r3**2 - b2*r4**2)
f(13) = r0**4*(-b1*r1 - b1*r2 + b2*r3 + b2*r4)
f(14) = r0*(b1*r3 + b1*r4 - b2*r1 - b2*r2)
f(15) = r0*(-b1**3*r3 - b1**3*r4 + b2**3*r1 + b2**3*r2)
f(16) = r0*(b1*r3**2 + b1*r4**2 - b2*r1**2 - b2*r2**2)
f(17) = r0*(b1**3*r3**2 + b1**3*r4**2 - b2**3*r1**2 - b2**3*r2**2)
f(18) = r0*(-b1*r3**3 - b1*r4**3 + b2*r1**3 + b2*r2**3)
f(19) = r0*(b1*r3**4 + b1*r4**4 - b2*r1**4 - b2*r2**4)
f(20) = r0**2*(b1*r3 + b1*r4 - b2*r1 - b2*r2)
f(21) = r0**2*(b1**3*r3 + b1**3*r4 - b2**3*r1 - b2**3*r2)
f(22) = r0**2*(b1*r3**2 + b1*r4**2 - b2*r1**2 - b2*r2**2)
f(23) = r0**2*(-b1*r3**3 - b1*r4**3 + b2*r1**3 + b2*r2**3)
f(24) = r0**3*(-b1*r3 - b1*r4 + b2*r1 + b2*r2)
f(25) = r0**3*(b1*r3**2 + b1*r4**2 - b2*r1**2 - b2*r2**2)
f(26) = r0**4*(-b1*r3 - b1*r4 + b2*r1 + b2*r2)
f(27) = dtau*r0*(r1 - r2 - r3 + r4)
f(28) = dtau**3*r0*(r1 - r2 - r3 + r4)
f(29) = dtau*r0*(-r1**2 + r2**2 + r3**2 - r4**2)
f(30) = dtau**3*r0*(r1**2 - r2**2 - r3**2 + r4**2)
f(31) = dtau*r0*(-r1**3 + r2**3 + r3**3 - r4**3)
f(32) = dtau*r0*(-r1**4 + r2**4 + r3**4 - r4**4)
f(33) = dtau*r0**2*(-r1 + r2 + r3 - r4)
f(34) = dtau**3*r0**2*(r1 - r2 - r3 + r4)
f(35) = dtau*r0**2*(r1**2 - r2**2 - r3**2 + r4**2)
f(36) = dtau*r0**2*(-r1**3 + r2**3 + r3**3 - r4**3)
f(37) = dtau*r0**3*(r1 - r2 - r3 + r4)
f(38) = dtau*r0**3*(-r1**2 + r2**2 + r3**2 - r4**2)
f(39) = dtau*r0**4*(r1 - r2 - r3 + r4)
f(40) = r0*(-a1*b1 - a2*b1 + a3*b2 + a4*b2)
f(41) = r0*(-a1*b1**3 - a2*b1**3 + a3*b2**3 + a4*b2**3)
f(42) = r0*(-a1**2*b1 - a2**2*b1 + a3**2*b2 + a4**2*b2)
f(43) = r0*(a1**2*b1**3 + a2**2*b1**3 - a3**2*b2**3 - a4**2*b2**3)
f(44) = r0*(-a1**3*b1 - a2**3*b1 + a3**3*b2 + a4**3*b2)
f(45) = r0*(-a1**4*b1 - a2**4*b1 + a3**4*b2 + a4**4*b2)
f(46) = r0**2*(-a1*b1 - a2*b1 + a3*b2 + a4*b2)
f(47) = r0**2*(a1*b1**3 + a2*b1**3 - a3*b2**3 - a4*b2**3)
f(48) = r0**2*(a1**2*b1 + a2**2*b1 - a3**2*b2 - a4**2*b2)
f(49) = r0**2*(a1**3*b1 + a2**3*b1 - a3**3*b2 - a4**3*b2)
f(50) = r0**3*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(51) = r0**3*(-a1**2*b1 - a2**2*b1 + a3**2*b2 + a4**2*b2)
f(52) = r0**4*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(53) = r0*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(54) = r0*(a1*b2**3 + a2*b2**3 - a3*b1**3 - a4*b1**3)
f(55) = r0*(a1**2*b2 + a2**2*b2 - a3**2*b1 - a4**2*b1)
f(56) = r0*(-a1**2*b2**3 - a2**2*b2**3 + a3**2*b1**3 + a4**2*b1**3)
f(57) = r0*(a1**3*b2 + a2**3*b2 - a3**3*b1 - a4**3*b1)
f(58) = r0*(a1**4*b2 + a2**4*b2 - a3**4*b1 - a4**4*b1)
f(59) = r0**2*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(60) = r0**2*(a1*b2**3 + a2*b2**3 - a3*b1**3 - a4*b1**3)
f(61) = r0**2*(-a1**2*b2 - a2**2*b2 + a3**2*b1 + a4**2*b1)
f(62) = r0**2*(-a1**3*b2 - a2**3*b2 + a3**3*b1 + a4**3*b1)
f(63) = r0**3*(-a1*b2 - a2*b2 + a3*b1 + a4*b1)
f(64) = r0**3*(a1**2*b2 + a2**2*b2 - a3**2*b1 - a4**2*b1)
f(65) = r0**4*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(66) = dtau*r0*(a1 - a2 - a3 + a4)
f(67) = dtau**3*r0*(a1 - a2 - a3 + a4)
f(68) = dtau*r0*(a1**2 - a2**2 - a3**2 + a4**2)
f(69) = dtau**3*r0*(-a1**2 + a2**2 + a3**2 - a4**2)
f(70) = dtau*r0*(-a1**3 + a2**3 + a3**3 - a4**3)
f(71) = dtau*r0*(-a1**4 + a2**4 + a3**4 - a4**4)
f(72) = dtau*r0**2*(a1 - a2 - a3 + a4)
f(73) = dtau**3*r0**2*(-a1 + a2 + a3 - a4)
f(74) = dtau*r0**2*(a1**2 - a2**2 - a3**2 + a4**2)
f(75) = dtau*r0**2*(-a1**3 + a2**3 + a3**3 - a4**3)
f(76) = dtau*r0**3*(a1 - a2 - a3 + a4)
f(77) = dtau*r0**3*(-a1**2 + a2**2 + a3**2 - a4**2)
f(78) = dtau*r0**4*(a1 - a2 - a3 + a4)
f(79) = b1*b2*r0*(-b1 + b2)
f(80) = b1*b2*r0*(b1**3 - b2**3)
f(81) = b1**2*b2**2*r0*(b1 - b2)
f(82) = b1*b2*r0**2*(-b1 + b2)
f(83) = b1*b2*r0**3*(b1 - b2)
f(84) = dtau**2*r0*(-b1 + b2)
f(85) = dtau**4*r0*(b1 - b2)
f(86) = dtau**2*r0*(-b1**3 + b2**3)
f(87) = dtau**2*r0**2*(-b1 + b2)
f(88) = dtau**2*r0**3*(b1 - b2)
f(89) = b1*r1*r2 - b2*r3*r4
f(90) = b1**3*r1*r2 - b2**3*r3*r4
f(91) = -b1*r1**2*r2 - b1*r1*r2**2 + b2*r3**2*r4 + b2*r3*r4**2
f(92) = -b1**3*r1**2*r2 - b1**3*r1*r2**2 + b2**3*r3**2*r4 + b2**3*r3*r4 &
      **2
f(93) = b1*r1**3*r2 + b1*r1*r2**3 - b2*r3**3*r4 - b2*r3*r4**3
f(94) = b1*r1**4*r2 + b1*r1*r2**4 - b2*r3**4*r4 - b2*r3*r4**4
f(95) = b1*r1**2*r2**2 - b2*r3**2*r4**2
f(96) = -b1*r1**3*r2**2 - b1*r1**2*r2**3 + b2*r3**3*r4**2 + b2*r3**2*r4 &
      **3
f(97) = b1*r3*r4 - b2*r1*r2
f(98) = b1**3*r3*r4 - b2**3*r1*r2
f(99) = b1*r3**2*r4 + b1*r3*r4**2 - b2*r1**2*r2 - b2*r1*r2**2
f(100) = -b1**3*r3**2*r4 - b1**3*r3*r4**2 + b2**3*r1**2*r2 + b2**3*r1*r2 &
      **2
f(101) = b1*r3**3*r4 + b1*r3*r4**3 - b2*r1**3*r2 - b2*r1*r2**3
f(102) = b1*r3**4*r4 + b1*r3*r4**4 - b2*r1**4*r2 - b2*r1*r2**4
f(103) = b1*r3**2*r4**2 - b2*r1**2*r2**2
f(104) = -b1*r3**3*r4**2 - b1*r3**2*r4**3 + b2*r1**3*r2**2 + b2*r1**2*r2 &
      **3
f(105) = dtau*(r1**2*r2 - r1*r2**2 - r3**2*r4 + r3*r4**2)
f(106) = dtau**3*(-r1**2*r2 + r1*r2**2 + r3**2*r4 - r3*r4**2)
f(107) = dtau*(r1**3*r2 - r1*r2**3 - r3**3*r4 + r3*r4**3)
f(108) = dtau*(r1**4*r2 - r1*r2**4 - r3**4*r4 + r3*r4**4)
f(109) = dtau**3*(r1**2*r2 - r1*r2**2 - r3**2*r4 + r3*r4**2)
f(110) = dtau*(-r1**3*r2**2 + r1**2*r2**3 + r3**3*r4**2 - r3**2*r4**3)
f(111) = dtau*(r1**3*r2**2 - r1**2*r2**3 - r3**3*r4**2 + r3**2*r4**3)
f(112) = -b1*r1*r3 - b1*r2*r4 + b2*r1*r3 + b2*r2*r4
f(113) = b1**3*r1*r3 + b1**3*r2*r4 - b2**3*r1*r3 - b2**3*r2*r4
f(114) = -b1*r1*r3**2 - b1*r2*r4**2 + b2*r1**2*r3 + b2*r2**2*r4
f(115) = -b1**3*r1*r3**2 - b1**3*r2*r4**2 + b2**3*r1**2*r3 + b2**3*r2**2 &
      *r4
f(116) = b1*r1*r3**3 + b1*r2*r4**3 - b2*r1**3*r3 - b2*r2**3*r4
f(117) = b1*r1*r3**4 + b1*r2*r4**4 - b2*r1**4*r3 - b2*r2**4*r4
f(118) = -b1*r1**2*r3 - b1*r2**2*r4 + b2*r1*r3**2 + b2*r2*r4**2
f(119) = -b1**3*r1**2*r3 - b1**3*r2**2*r4 + b2**3*r1*r3**2 + b2**3*r2*r4 &
      **2
f(120) = -b1*r1**2*r3**2 - b1*r2**2*r4**2 + b2*r1**2*r3**2 + b2*r2**2*r4 &
      **2
f(121) = b1*r1**2*r3**3 + b1*r2**2*r4**3 - b2*r1**3*r3**2 - b2*r2**3*r4 &
      **2
f(122) = b1*r1**3*r3 + b1*r2**3*r4 - b2*r1*r3**3 - b2*r2*r4**3
f(123) = b1*r1**3*r3**2 + b1*r2**3*r4**2 - b2*r1**2*r3**3 - b2*r2**2*r4 &
      **3
f(124) = b1*r1**4*r3 + b1*r2**4*r4 - b2*r1*r3**4 - b2*r2*r4**4
f(125) = dtau*(-r1**2*r3 + r1*r3**2 + r2**2*r4 - r2*r4**2)
f(126) = dtau**3*(-r1**2*r3 + r1*r3**2 + r2**2*r4 - r2*r4**2)
f(127) = dtau*(-r1**3*r3 + r1*r3**3 + r2**3*r4 - r2*r4**3)
f(128) = dtau*(-r1**4*r3 + r1*r3**4 + r2**4*r4 - r2*r4**4)
f(129) = dtau**3*(r1**2*r3 - r1*r3**2 - r2**2*r4 + r2*r4**2)
f(130) = dtau*(-r1**3*r3**2 + r1**2*r3**3 + r2**3*r4**2 - r2**2*r4**3)
f(131) = dtau*(r1**3*r3**2 - r1**2*r3**3 - r2**3*r4**2 + r2**2*r4**3)
f(132) = b1*r1*r4 + b1*r2*r3 - b2*r1*r4 - b2*r2*r3
f(133) = -b1**3*r1*r4 - b1**3*r2*r3 + b2**3*r1*r4 + b2**3*r2*r3
f(134) = b1*r1*r4**2 + b1*r2*r3**2 - b2*r1**2*r4 - b2*r2**2*r3
f(135) = b1**3*r1*r4**2 + b1**3*r2*r3**2 - b2**3*r1**2*r4 - b2**3*r2**2* &
      r3
f(136) = -b1*r1*r4**3 - b1*r2*r3**3 + b2*r1**3*r4 + b2*r2**3*r3
f(137) = b1*r1*r4**4 + b1*r2*r3**4 - b2*r1**4*r4 - b2*r2**4*r3
f(138) = b1*r1**2*r4 + b1*r2**2*r3 - b2*r1*r4**2 - b2*r2*r3**2
f(139) = b1**3*r1**2*r4 + b1**3*r2**2*r3 - b2**3*r1*r4**2 - b2**3*r2*r3 &
      **2
f(140) = -b1*r1**2*r4**2 - b1*r2**2*r3**2 + b2*r1**2*r4**2 + b2*r2**2*r3 &
      **2
f(141) = b1*r1**2*r4**3 + b1*r2**2*r3**3 - b2*r1**3*r4**2 - b2*r2**3*r3 &
      **2
f(142) = -b1*r1**3*r4 - b1*r2**3*r3 + b2*r1*r4**3 + b2*r2*r3**3
f(143) = b1*r1**3*r4**2 + b1*r2**3*r3**2 - b2*r1**2*r4**3 - b2*r2**2*r3 &
      **3
f(144) = -b1*r1**4*r4 - b1*r2**4*r3 + b2*r1*r4**4 + b2*r2*r3**4
f(145) = dtau*(-r1*r4 + r2*r3)
f(146) = dtau**3*(r1*r4 - r2*r3)
f(147) = dtau*(r1**2*r4 + r1*r4**2 - r2**2*r3 - r2*r3**2)
f(148) = dtau**3*(-r1**2*r4 - r1*r4**2 + r2**2*r3 + r2*r3**2)
f(149) = dtau*(r1**3*r4 + r1*r4**3 - r2**3*r3 - r2*r3**3)
f(150) = dtau*(r1**4*r4 + r1*r4**4 - r2**4*r3 - r2*r3**4)
f(151) = dtau**3*(r1**2*r4 + r1*r4**2 - r2**2*r3 - r2*r3**2)
f(152) = dtau*(-r1**2*r4**2 + r2**2*r3**2)
f(153) = dtau*(r1**3*r4**2 + r1**2*r4**3 - r2**3*r3**2 - r2**2*r3**3)
f(154) = dtau*(-r1**3*r4 - r1*r4**3 + r2**3*r3 + r2*r3**3)
f(155) = -a1*b1*r1 - a2*b1*r2 + a3*b2*r3 + a4*b2*r4
f(156) = -a1*b1**3*r1 - a2*b1**3*r2 + a3*b2**3*r3 + a4*b2**3*r4
f(157) = a1**2*b1*r1 + a2**2*b1*r2 - a3**2*b2*r3 - a4**2*b2*r4
f(158) = a1**2*b1**3*r1 + a2**2*b1**3*r2 - a3**2*b2**3*r3 - a4**2*b2**3* &
      r4
f(159) = -a1**3*b1*r1 - a2**3*b1*r2 + a3**3*b2*r3 + a4**3*b2*r4
f(160) = -a1**4*b1*r1 - a2**4*b1*r2 + a3**4*b2*r3 + a4**4*b2*r4
f(161) = a1*b1*r1**2 + a2*b1*r2**2 - a3*b2*r3**2 - a4*b2*r4**2
f(162) = a1*b1**3*r1**2 + a2*b1**3*r2**2 - a3*b2**3*r3**2 - a4*b2**3*r4 &
      **2
f(163) = -a1**2*b1*r1**2 - a2**2*b1*r2**2 + a3**2*b2*r3**2 + a4**2*b2*r4 &
      **2
f(164) = a1**3*b1*r1**2 + a2**3*b1*r2**2 - a3**3*b2*r3**2 - a4**3*b2*r4 &
      **2
f(165) = -a1*b1*r1**3 - a2*b1*r2**3 + a3*b2*r3**3 + a4*b2*r4**3
f(166) = a1**2*b1*r1**3 + a2**2*b1*r2**3 - a3**2*b2*r3**3 - a4**2*b2*r4 &
      **3
f(167) = -a1*b1*r1**4 - a2*b1*r2**4 + a3*b2*r3**4 + a4*b2*r4**4
f(168) = a1*b2*r1 + a2*b2*r2 - a3*b1*r3 - a4*b1*r4
f(169) = a1*b2**3*r1 + a2*b2**3*r2 - a3*b1**3*r3 - a4*b1**3*r4
f(170) = -a1**2*b2*r1 - a2**2*b2*r2 + a3**2*b1*r3 + a4**2*b1*r4
f(171) = -a1**2*b2**3*r1 - a2**2*b2**3*r2 + a3**2*b1**3*r3 + a4**2*b1**3 &
      *r4
f(172) = a1**3*b2*r1 + a2**3*b2*r2 - a3**3*b1*r3 - a4**3*b1*r4
f(173) = a1**4*b2*r1 + a2**4*b2*r2 - a3**4*b1*r3 - a4**4*b1*r4
f(174) = -a1*b2*r1**2 - a2*b2*r2**2 + a3*b1*r3**2 + a4*b1*r4**2
f(175) = -a1*b2**3*r1**2 - a2*b2**3*r2**2 + a3*b1**3*r3**2 + a4*b1**3*r4 &
      **2
f(176) = a1**2*b2*r1**2 + a2**2*b2*r2**2 - a3**2*b1*r3**2 - a4**2*b1*r4 &
      **2
f(177) = -a1**3*b2*r1**2 - a2**3*b2*r2**2 + a3**3*b1*r3**2 + a4**3*b1*r4 &
      **2
f(178) = a1*b2*r1**3 + a2*b2*r2**3 - a3*b1*r3**3 - a4*b1*r4**3
f(179) = -a1**2*b2*r1**3 - a2**2*b2*r2**3 + a3**2*b1*r3**3 + a4**2*b1*r4 &
      **3
f(180) = a1*b2*r1**4 + a2*b2*r2**4 - a3*b1*r3**4 - a4*b1*r4**4
f(181) = dtau*(-a1*r1 + a2*r2 + a3*r3 - a4*r4)
f(182) = dtau**3*(a1*r1 - a2*r2 - a3*r3 + a4*r4)
f(183) = dtau*(a1**2*r1 - a2**2*r2 - a3**2*r3 + a4**2*r4)
f(184) = dtau**3*(-a1**2*r1 + a2**2*r2 + a3**2*r3 - a4**2*r4)
f(185) = dtau*(a1**3*r1 - a2**3*r2 - a3**3*r3 + a4**3*r4)
f(186) = dtau*(-a1**4*r1 + a2**4*r2 + a3**4*r3 - a4**4*r4)
f(187) = dtau*(-a1*r1**2 + a2*r2**2 + a3*r3**2 - a4*r4**2)
f(188) = dtau**3*(a1*r1**2 - a2*r2**2 - a3*r3**2 + a4*r4**2)
f(189) = dtau*(a1**2*r1**2 - a2**2*r2**2 - a3**2*r3**2 + a4**2*r4**2)
f(190) = dtau*(-a1**3*r1**2 + a2**3*r2**2 + a3**3*r3**2 - a4**3*r4**2)
f(191) = dtau*(a1*r1**3 - a2*r2**3 - a3*r3**3 + a4*r4**3)
f(192) = dtau*(a1**2*r1**3 - a2**2*r2**3 - a3**2*r3**3 + a4**2*r4**3)
f(193) = dtau*(a1*r1**4 - a2*r2**4 - a3*r3**4 + a4*r4**4)
f(194) = a1*b1*r2 + a2*b1*r1 - a3*b2*r4 - a4*b2*r3
f(195) = a1*b1**3*r2 + a2*b1**3*r1 - a3*b2**3*r4 - a4*b2**3*r3
f(196) = -a1**2*b1*r2 - a2**2*b1*r1 + a3**2*b2*r4 + a4**2*b2*r3
f(197) = -a1**2*b1**3*r2 - a2**2*b1**3*r1 + a3**2*b2**3*r4 + a4**2*b2**3 &
      *r3
f(198) = a1**3*b1*r2 + a2**3*b1*r1 - a3**3*b2*r4 - a4**3*b2*r3
f(199) = a1**4*b1*r2 + a2**4*b1*r1 - a3**4*b2*r4 - a4**4*b2*r3
f(200) = -a1*b1*r2**2 - a2*b1*r1**2 + a3*b2*r4**2 + a4*b2*r3**2
f(201) = -a1*b1**3*r2**2 - a2*b1**3*r1**2 + a3*b2**3*r4**2 + a4*b2**3*r3 &
      **2
f(202) = a1**2*b1*r2**2 + a2**2*b1*r1**2 - a3**2*b2*r4**2 - a4**2*b2*r3 &
      **2
f(203) = -a1**3*b1*r2**2 - a2**3*b1*r1**2 + a3**3*b2*r4**2 + a4**3*b2*r3 &
      **2
f(204) = a1*b1*r2**3 + a2*b1*r1**3 - a3*b2*r4**3 - a4*b2*r3**3
f(205) = -a1**2*b1*r2**3 - a2**2*b1*r1**3 + a3**2*b2*r4**3 + a4**2*b2*r3 &
      **3
f(206) = a1*b1*r2**4 + a2*b1*r1**4 - a3*b2*r4**4 - a4*b2*r3**4
f(207) = -a1*b2*r2 - a2*b2*r1 + a3*b1*r4 + a4*b1*r3
f(208) = -a1*b2**3*r2 - a2*b2**3*r1 + a3*b1**3*r4 + a4*b1**3*r3
f(209) = a1**2*b2*r2 + a2**2*b2*r1 - a3**2*b1*r4 - a4**2*b1*r3
f(210) = a1**2*b2**3*r2 + a2**2*b2**3*r1 - a3**2*b1**3*r4 - a4**2*b1**3* &
      r3
f(211) = -a1**3*b2*r2 - a2**3*b2*r1 + a3**3*b1*r4 + a4**3*b1*r3
f(212) = -a1**4*b2*r2 - a2**4*b2*r1 + a3**4*b1*r4 + a4**4*b1*r3
f(213) = a1*b2*r2**2 + a2*b2*r1**2 - a3*b1*r4**2 - a4*b1*r3**2
f(214) = a1*b2**3*r2**2 + a2*b2**3*r1**2 - a3*b1**3*r4**2 - a4*b1**3*r3 &
      **2
f(215) = -a1**2*b2*r2**2 - a2**2*b2*r1**2 + a3**2*b1*r4**2 + a4**2*b1*r3 &
      **2
f(216) = a1**3*b2*r2**2 + a2**3*b2*r1**2 - a3**3*b1*r4**2 - a4**3*b1*r3 &
      **2
f(217) = a1*b2*r2**3 + a2*b2*r1**3 - a3*b1*r4**3 - a4*b1*r3**3
f(218) = a1**2*b2*r2**3 + a2**2*b2*r1**3 - a3**2*b1*r4**3 - a4**2*b1*r3 &
      **3
f(219) = -a1*b2*r2**4 - a2*b2*r1**4 + a3*b1*r4**4 + a4*b1*r3**4
f(220) = dtau*(-a1*r2 + a2*r1 + a3*r4 - a4*r3)
f(221) = dtau**3*(a1*r2 - a2*r1 - a3*r4 + a4*r3)
f(222) = dtau*(-a1**2*r2 + a2**2*r1 + a3**2*r4 - a4**2*r3)
f(223) = dtau**3*(-a1**2*r2 + a2**2*r1 + a3**2*r4 - a4**2*r3)
f(224) = dtau*(-a1**3*r2 + a2**3*r1 + a3**3*r4 - a4**3*r3)
f(225) = dtau*(a1**4*r2 - a2**4*r1 - a3**4*r4 + a4**4*r3)
f(226) = dtau*(-a1*r2**2 + a2*r1**2 + a3*r4**2 - a4*r3**2)
f(227) = dtau**3*(a1*r2**2 - a2*r1**2 - a3*r4**2 + a4*r3**2)
f(228) = dtau*(-a1**2*r2**2 + a2**2*r1**2 + a3**2*r4**2 - a4**2*r3**2)
f(229) = dtau*(a1**3*r2**2 - a2**3*r1**2 - a3**3*r4**2 + a4**3*r3**2)
f(230) = dtau*(a1*r2**3 - a2*r1**3 - a3*r4**3 + a4*r3**3)
f(231) = dtau*(a1**2*r2**3 - a2**2*r1**3 - a3**2*r4**3 + a4**2*r3**3)
f(232) = dtau*(-a1*r2**4 + a2*r1**4 + a3*r4**4 - a4*r3**4)
f(233) = -a1*b2*r3 - a2*b2*r4 + a3*b1*r1 + a4*b1*r2
f(234) = -a1*b2**3*r3 - a2*b2**3*r4 + a3*b1**3*r1 + a4*b1**3*r2
f(235) = a1**2*b2*r3 + a2**2*b2*r4 - a3**2*b1*r1 - a4**2*b1*r2
f(236) = a1**2*b2**3*r3 + a2**2*b2**3*r4 - a3**2*b1**3*r1 - a4**2*b1**3* &
      r2
f(237) = -a1**3*b2*r3 - a2**3*b2*r4 + a3**3*b1*r1 + a4**3*b1*r2
f(238) = -a1**4*b2*r3 - a2**4*b2*r4 + a3**4*b1*r1 + a4**4*b1*r2
f(239) = a1*b2*r3**2 + a2*b2*r4**2 - a3*b1*r1**2 - a4*b1*r2**2
f(240) = a1*b2**3*r3**2 + a2*b2**3*r4**2 - a3*b1**3*r1**2 - a4*b1**3*r2 &
      **2
f(241) = -a1**2*b2*r3**2 - a2**2*b2*r4**2 + a3**2*b1*r1**2 + a4**2*b1*r2 &
      **2
f(242) = a1**3*b2*r3**2 + a2**3*b2*r4**2 - a3**3*b1*r1**2 - a4**3*b1*r2 &
      **2
f(243) = -a1*b2*r3**3 - a2*b2*r4**3 + a3*b1*r1**3 + a4*b1*r2**3
f(244) = a1**2*b2*r3**3 + a2**2*b2*r4**3 - a3**2*b1*r1**3 - a4**2*b1*r2 &
      **3
f(245) = -a1*b2*r3**4 - a2*b2*r4**4 + a3*b1*r1**4 + a4*b1*r2**4
f(246) = a1*b1*r3 + a2*b1*r4 - a3*b2*r1 - a4*b2*r2
f(247) = a1*b1**3*r3 + a2*b1**3*r4 - a3*b2**3*r1 - a4*b2**3*r2
f(248) = -a1**2*b1*r3 - a2**2*b1*r4 + a3**2*b2*r1 + a4**2*b2*r2
f(249) = -a1**2*b1**3*r3 - a2**2*b1**3*r4 + a3**2*b2**3*r1 + a4**2*b2**3 &
      *r2
f(250) = a1**3*b1*r3 + a2**3*b1*r4 - a3**3*b2*r1 - a4**3*b2*r2
f(251) = a1**4*b1*r3 + a2**4*b1*r4 - a3**4*b2*r1 - a4**4*b2*r2
f(252) = -a1*b1*r3**2 - a2*b1*r4**2 + a3*b2*r1**2 + a4*b2*r2**2
f(253) = -a1*b1**3*r3**2 - a2*b1**3*r4**2 + a3*b2**3*r1**2 + a4*b2**3*r2 &
      **2
f(254) = a1**2*b1*r3**2 + a2**2*b1*r4**2 - a3**2*b2*r1**2 - a4**2*b2*r2 &
      **2
f(255) = -a1**3*b1*r3**2 - a2**3*b1*r4**2 + a3**3*b2*r1**2 + a4**3*b2*r2 &
      **2
f(256) = a1*b1*r3**3 + a2*b1*r4**3 - a3*b2*r1**3 - a4*b2*r2**3
f(257) = -a1**2*b1*r3**3 - a2**2*b1*r4**3 + a3**2*b2*r1**3 + a4**2*b2*r2 &
      **3
f(258) = a1*b1*r3**4 + a2*b1*r4**4 - a3*b2*r1**4 - a4*b2*r2**4
f(259) = dtau*(a1*r3 - a2*r4 - a3*r1 + a4*r2)
f(260) = dtau**3*(-a1*r3 + a2*r4 + a3*r1 - a4*r2)
f(261) = dtau*(a1**2*r3 - a2**2*r4 - a3**2*r1 + a4**2*r2)
f(262) = dtau**3*(-a1**2*r3 + a2**2*r4 + a3**2*r1 - a4**2*r2)
f(263) = dtau*(a1**3*r3 - a2**3*r4 - a3**3*r1 + a4**3*r2)
f(264) = dtau*(a1**4*r3 - a2**4*r4 - a3**4*r1 + a4**4*r2)
f(265) = dtau*(-a1*r3**2 + a2*r4**2 + a3*r1**2 - a4*r2**2)
f(266) = dtau**3*(a1*r3**2 - a2*r4**2 - a3*r1**2 + a4*r2**2)
f(267) = dtau*(-a1**2*r3**2 + a2**2*r4**2 + a3**2*r1**2 - a4**2*r2**2)
f(268) = dtau*(-a1**3*r3**2 + a2**3*r4**2 + a3**3*r1**2 - a4**3*r2**2)
f(269) = dtau*(-a1*r3**3 + a2*r4**3 + a3*r1**3 - a4*r2**3)
f(270) = dtau*(-a1**2*r3**3 + a2**2*r4**3 + a3**2*r1**3 - a4**2*r2**3)
f(271) = dtau*(-a1*r3**4 + a2*r4**4 + a3*r1**4 - a4*r2**4)
f(272) = a1*b2*r4 + a2*b2*r3 - a3*b1*r2 - a4*b1*r1
f(273) = a1*b2**3*r4 + a2*b2**3*r3 - a3*b1**3*r2 - a4*b1**3*r1
f(274) = -a1**2*b2*r4 - a2**2*b2*r3 + a3**2*b1*r2 + a4**2*b1*r1
f(275) = -a1**2*b2**3*r4 - a2**2*b2**3*r3 + a3**2*b1**3*r2 + a4**2*b1**3 &
      *r1
f(276) = a1**3*b2*r4 + a2**3*b2*r3 - a3**3*b1*r2 - a4**3*b1*r1
f(277) = a1**4*b2*r4 + a2**4*b2*r3 - a3**4*b1*r2 - a4**4*b1*r1
f(278) = -a1*b2*r4**2 - a2*b2*r3**2 + a3*b1*r2**2 + a4*b1*r1**2
f(279) = -a1*b2**3*r4**2 - a2*b2**3*r3**2 + a3*b1**3*r2**2 + a4*b1**3*r1 &
      **2
f(280) = a1**2*b2*r4**2 + a2**2*b2*r3**2 - a3**2*b1*r2**2 - a4**2*b1*r1 &
      **2
f(281) = -a1**3*b2*r4**2 - a2**3*b2*r3**2 + a3**3*b1*r2**2 + a4**3*b1*r1 &
      **2
f(282) = a1*b2*r4**3 + a2*b2*r3**3 - a3*b1*r2**3 - a4*b1*r1**3
f(283) = -a1**2*b2*r4**3 - a2**2*b2*r3**3 + a3**2*b1*r2**3 + a4**2*b1*r1 &
      **3
f(284) = a1*b2*r4**4 + a2*b2*r3**4 - a3*b1*r2**4 - a4*b1*r1**4
f(285) = -a1*b1*r4 - a2*b1*r3 + a3*b2*r2 + a4*b2*r1
f(286) = -a1*b1**3*r4 - a2*b1**3*r3 + a3*b2**3*r2 + a4*b2**3*r1
f(287) = a1**2*b1*r4 + a2**2*b1*r3 - a3**2*b2*r2 - a4**2*b2*r1
f(288) = a1**2*b1**3*r4 + a2**2*b1**3*r3 - a3**2*b2**3*r2 - a4**2*b2**3* &
      r1
f(289) = -a1**3*b1*r4 - a2**3*b1*r3 + a3**3*b2*r2 + a4**3*b2*r1
f(290) = -a1**4*b1*r4 - a2**4*b1*r3 + a3**4*b2*r2 + a4**4*b2*r1
f(291) = a1*b1*r4**2 + a2*b1*r3**2 - a3*b2*r2**2 - a4*b2*r1**2
f(292) = a1*b1**3*r4**2 + a2*b1**3*r3**2 - a3*b2**3*r2**2 - a4*b2**3*r1 &
      **2
f(293) = -a1**2*b1*r4**2 - a2**2*b1*r3**2 + a3**2*b2*r2**2 + a4**2*b2*r1 &
      **2
f(294) = a1**3*b1*r4**2 + a2**3*b1*r3**2 - a3**3*b2*r2**2 - a4**3*b2*r1 &
      **2
f(295) = -a1*b1*r4**3 - a2*b1*r3**3 + a3*b2*r2**3 + a4*b2*r1**3
f(296) = a1**2*b1*r4**3 + a2**2*b1*r3**3 - a3**2*b2*r2**3 - a4**2*b2*r1 &
      **3
f(297) = -a1*b1*r4**4 - a2*b1*r3**4 + a3*b2*r2**4 + a4*b2*r1**4
f(298) = dtau*(-a1*r4 + a2*r3 + a3*r2 - a4*r1)
f(299) = dtau**3*(-a1*r4 + a2*r3 + a3*r2 - a4*r1)
f(300) = dtau*(-a1**2*r4 + a2**2*r3 + a3**2*r2 - a4**2*r1)
f(301) = dtau**3*(a1**2*r4 - a2**2*r3 - a3**2*r2 + a4**2*r1)
f(302) = dtau*(a1**3*r4 - a2**3*r3 - a3**3*r2 + a4**3*r1)
f(303) = dtau*(a1**4*r4 - a2**4*r3 - a3**4*r2 + a4**4*r1)
f(304) = dtau*(-a1*r4**2 + a2*r3**2 + a3*r2**2 - a4*r1**2)
f(305) = dtau**3*(a1*r4**2 - a2*r3**2 - a3*r2**2 + a4*r1**2)
f(306) = dtau*(-a1**2*r4**2 + a2**2*r3**2 + a3**2*r2**2 - a4**2*r1**2)
f(307) = dtau*(-a1**3*r4**2 + a2**3*r3**2 + a3**3*r2**2 - a4**3*r1**2)
f(308) = dtau*(a1*r4**3 - a2*r3**3 - a3*r2**3 + a4*r1**3)
f(309) = dtau*(-a1**2*r4**3 + a2**2*r3**3 + a3**2*r2**3 - a4**2*r1**3)
f(310) = dtau*(a1*r4**4 - a2*r3**4 - a3*r2**4 + a4*r1**4)
f(311) = b1*b2*(-b1*r3 - b1*r4 + b2*r1 + b2*r2)
f(312) = b1*b2*(-b1**3*r3 - b1**3*r4 + b2**3*r1 + b2**3*r2)
f(313) = b1*b2*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(314) = b1**2*b2**2*(-b1*r3 - b1*r4 + b2*r1 + b2*r2)
f(315) = b1**2*b2**2*(-b1*r1 - b1*r2 + b2*r3 + b2*r4)
f(316) = b1*b2*(b1**3*r1 + b1**3*r2 - b2**3*r3 - b2**3*r4)
f(317) = b1*b2*(-b1*r3**2 - b1*r4**2 + b2*r1**2 + b2*r2**2)
f(318) = b1*b2*(b1*r1**2 + b1*r2**2 - b2*r3**2 - b2*r4**2)
f(319) = b1*b2*(-b1*r3**3 - b1*r4**3 + b2*r1**3 + b2*r2**3)
f(320) = b1*b2*(b1*r1**3 + b1*r2**3 - b2*r3**3 - b2*r4**3)
f(321) = dtau**2*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(322) = dtau**4*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(323) = dtau*(-b1**2*r1 + b1**2*r2 + b2**2*r3 - b2**2*r4)
f(324) = dtau**3*(b1**2*r1 - b1**2*r2 - b2**2*r3 + b2**2*r4)
f(325) = dtau**2*(b1**3*r1 + b1**3*r2 - b2**3*r3 - b2**3*r4)
f(326) = dtau*(-b1**4*r1 + b1**4*r2 + b2**4*r3 - b2**4*r4)
f(327) = dtau**2*(b1*r1**2 + b1*r2**2 - b2*r3**2 - b2*r4**2)
f(328) = dtau*(b1**2*r1**2 - b1**2*r2**2 - b2**2*r3**2 + b2**2*r4**2)
f(329) = dtau**2*(b1*r1**3 + b1*r2**3 - b2*r3**3 - b2*r4**3)
f(330) = dtau*(b1**2*r1**3 - b1**2*r2**3 - b2**2*r3**3 + b2**2*r4**3)
f(331) = dtau**2*(b1*r3 + b1*r4 - b2*r1 - b2*r2)
f(332) = dtau**4*(-b1*r3 - b1*r4 + b2*r1 + b2*r2)
f(333) = dtau*(b1**2*r3 - b1**2*r4 - b2**2*r1 + b2**2*r2)
f(334) = dtau**3*(-b1**2*r3 + b1**2*r4 + b2**2*r1 - b2**2*r2)
f(335) = dtau**2*(b1**3*r3 + b1**3*r4 - b2**3*r1 - b2**3*r2)
f(336) = dtau*(b1**4*r3 - b1**4*r4 - b2**4*r1 + b2**4*r2)
f(337) = dtau**2*(b1*r3**2 + b1*r4**2 - b2*r1**2 - b2*r2**2)
f(338) = dtau*(-b1**2*r3**2 + b1**2*r4**2 + b2**2*r1**2 - b2**2*r2**2)
f(339) = dtau**2*(-b1*r3**3 - b1*r4**3 + b2*r1**3 + b2*r2**3)
f(340) = dtau*(-b1**2*r3**3 + b1**2*r4**3 + b2**2*r1**3 - b2**2*r2**3)
f(341) = -a1*a2*b1 + a3*a4*b2
f(342) = a1*a2*b1**3 - a3*a4*b2**3
f(343) = -a1**2*a2*b1 - a1*a2**2*b1 + a3**2*a4*b2 + a3*a4**2*b2
f(344) = -a1**2*a2*b1**3 - a1*a2**2*b1**3 + a3**2*a4*b2**3 + a3*a4**2*b2 &
      **3
f(345) = -a1**3*a2*b1 - a1*a2**3*b1 + a3**3*a4*b2 + a3*a4**3*b2
f(346) = a1**4*a2*b1 + a1*a2**4*b1 - a3**4*a4*b2 - a3*a4**4*b2
f(347) = a1**2*a2**2*b1 - a3**2*a4**2*b2
f(348) = -a1**3*a2**2*b1 - a1**2*a2**3*b1 + a3**3*a4**2*b2 + a3**2*a4**3 &
      *b2
f(349) = -a1*a2*b2 + a3*a4*b1
f(350) = -a1*a2*b2**3 + a3*a4*b1**3
f(351) = -a1**2*a2*b2 - a1*a2**2*b2 + a3**2*a4*b1 + a3*a4**2*b1
f(352) = a1**2*a2*b2**3 + a1*a2**2*b2**3 - a3**2*a4*b1**3 - a3*a4**2*b1 &
      **3
f(353) = -a1**3*a2*b2 - a1*a2**3*b2 + a3**3*a4*b1 + a3*a4**3*b1
f(354) = -a1**4*a2*b2 - a1*a2**4*b2 + a3**4*a4*b1 + a3*a4**4*b1
f(355) = -a1**2*a2**2*b2 + a3**2*a4**2*b1
f(356) = a1**3*a2**2*b2 + a1**2*a2**3*b2 - a3**3*a4**2*b1 - a3**2*a4**3* &
      b1
f(357) = dtau*(a1**2*a2 - a1*a2**2 - a3**2*a4 + a3*a4**2)
f(358) = dtau**3*(a1**2*a2 - a1*a2**2 - a3**2*a4 + a3*a4**2)
f(359) = dtau*(-a1**3*a2 + a1*a2**3 + a3**3*a4 - a3*a4**3)
f(360) = dtau*(-a1**4*a2 + a1*a2**4 + a3**4*a4 - a3*a4**4)
f(361) = dtau*(-a1**3*a2**2 + a1**2*a2**3 + a3**3*a4**2 - a3**2*a4**3)
f(362) = dtau*(a1**3*a2**2 - a1**2*a2**3 - a3**3*a4**2 + a3**2*a4**3)
f(363) = a1*a3*b1 - a1*a3*b2 + a2*a4*b1 - a2*a4*b2
f(364) = a1*a3*b1**3 - a1*a3*b2**3 + a2*a4*b1**3 - a2*a4*b2**3
f(365) = a1**2*a3*b2 - a1*a3**2*b1 + a2**2*a4*b2 - a2*a4**2*b1
f(366) = a1**2*a3*b2**3 - a1*a3**2*b1**3 + a2**2*a4*b2**3 - a2*a4**2*b1 &
      **3
f(367) = -a1**3*a3*b2 + a1*a3**3*b1 - a2**3*a4*b2 + a2*a4**3*b1
f(368) = -a1**4*a3*b2 + a1*a3**4*b1 - a2**4*a4*b2 + a2*a4**4*b1
f(369) = -a1**2*a3*b1 + a1*a3**2*b2 - a2**2*a4*b1 + a2*a4**2*b2
f(370) = a1**2*a3*b1**3 - a1*a3**2*b2**3 + a2**2*a4*b1**3 - a2*a4**2*b2 &
      **3
f(371) = a1**2*a3**2*b1 - a1**2*a3**2*b2 + a2**2*a4**2*b1 - a2**2*a4**2* &
      b2
f(372) = -a1**3*a3**2*b2 + a1**2*a3**3*b1 - a2**3*a4**2*b2 + a2**2*a4**3 &
      *b1
f(373) = a1**3*a3*b1 - a1*a3**3*b2 + a2**3*a4*b1 - a2*a4**3*b2
f(374) = a1**3*a3**2*b1 - a1**2*a3**3*b2 + a2**3*a4**2*b1 - a2**2*a4**3* &
      b2
f(375) = -a1**4*a3*b1 + a1*a3**4*b2 - a2**4*a4*b1 + a2*a4**4*b2
f(376) = dtau*(-a1**2*a3 + a1*a3**2 + a2**2*a4 - a2*a4**2)
f(377) = dtau**3*(a1**2*a3 - a1*a3**2 - a2**2*a4 + a2*a4**2)
f(378) = dtau*(a1**3*a3 - a1*a3**3 - a2**3*a4 + a2*a4**3)
f(379) = dtau*(-a1**4*a3 + a1*a3**4 + a2**4*a4 - a2*a4**4)
f(380) = dtau*(a1**3*a3**2 - a1**2*a3**3 - a2**3*a4**2 + a2**2*a4**3)
f(381) = -a1*a4*b1 + a1*a4*b2 - a2*a3*b1 + a2*a3*b2
f(382) = -a1*a4*b1**3 + a1*a4*b2**3 - a2*a3*b1**3 + a2*a3*b2**3
f(383) = -a1**2*a4*b2 + a1*a4**2*b1 - a2**2*a3*b2 + a2*a3**2*b1
f(384) = a1**2*a4*b2**3 - a1*a4**2*b1**3 + a2**2*a3*b2**3 - a2*a3**2*b1 &
      **3
f(385) = a1**3*a4*b2 - a1*a4**3*b1 + a2**3*a3*b2 - a2*a3**3*b1
f(386) = a1**4*a4*b2 - a1*a4**4*b1 + a2**4*a3*b2 - a2*a3**4*b1
f(387) = a1**2*a4*b1 - a1*a4**2*b2 + a2**2*a3*b1 - a2*a3**2*b2
f(388) = a1**2*a4*b1**3 - a1*a4**2*b2**3 + a2**2*a3*b1**3 - a2*a3**2*b2 &
      **3
f(389) = a1**2*a4**2*b1 - a1**2*a4**2*b2 + a2**2*a3**2*b1 - a2**2*a3**2* &
      b2
f(390) = -a1**3*a4**2*b2 + a1**2*a4**3*b1 - a2**3*a3**2*b2 + a2**2*a3**3 &
      *b1
f(391) = -a1**3*a4*b1 + a1*a4**3*b2 - a2**3*a3*b1 + a2*a3**3*b2
f(392) = a1**3*a4**2*b1 - a1**2*a4**3*b2 + a2**3*a3**2*b1 - a2**2*a3**3* &
      b2
f(393) = a1**4*a4*b1 - a1*a4**4*b2 + a2**4*a3*b1 - a2*a3**4*b2
f(394) = dtau*(a1*a4 - a2*a3)
f(395) = dtau**3*(a1*a4 - a2*a3)
f(396) = dtau*(-a1**2*a4 - a1*a4**2 + a2**2*a3 + a2*a3**2)
f(397) = dtau**3*(a1**2*a4 + a1*a4**2 - a2**2*a3 - a2*a3**2)
f(398) = dtau*(a1**3*a4 + a1*a4**3 - a2**3*a3 - a2*a3**3)
f(399) = dtau*(a1**4*a4 + a1*a4**4 - a2**4*a3 - a2*a3**4)
f(400) = dtau*(-a1**2*a4**2 + a2**2*a3**2)
f(401) = dtau*(a1**3*a4**2 + a1**2*a4**3 - a2**3*a3**2 - a2**2*a3**3)
f(402) = dtau*(-a1**3*a4**2 - a1**2*a4**3 + a2**3*a3**2 + a2**2*a3**3)
f(403) = b1*b2*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(404) = b1*b2*(-a1*b2**3 - a2*b2**3 + a3*b1**3 + a4*b1**3)
f(405) = b1*b2*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(406) = b1**2*b2**2*(-a1*b2 - a2*b2 + a3*b1 + a4*b1)
f(407) = b1**2*b2**2*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(408) = b1*b2*(a1*b1**3 + a2*b1**3 - a3*b2**3 - a4*b2**3)
f(409) = b1*b2*(-a1**2*b2 - a2**2*b2 + a3**2*b1 + a4**2*b1)
f(410) = b1*b2*(-a1**2*b1 - a2**2*b1 + a3**2*b2 + a4**2*b2)
f(411) = b1*b2*(-a1**3*b2 - a2**3*b2 + a3**3*b1 + a4**3*b1)
f(412) = b1*b2*(-a1**3*b1 - a2**3*b1 + a3**3*b2 + a4**3*b2)
f(413) = dtau**2*(-a1*b1 - a2*b1 + a3*b2 + a4*b2)
f(414) = dtau**4*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(415) = dtau*(a1*b1**2 - a2*b1**2 - a3*b2**2 + a4*b2**2)
f(416) = dtau**3*(a1*b1**2 - a2*b1**2 - a3*b2**2 + a4*b2**2)
f(417) = dtau**2*(a1*b1**3 + a2*b1**3 - a3*b2**3 - a4*b2**3)
f(418) = dtau*(a1*b1**4 - a2*b1**4 - a3*b2**4 + a4*b2**4)
f(419) = dtau**2*(-a1**2*b1 - a2**2*b1 + a3**2*b2 + a4**2*b2)
f(420) = dtau*(-a1**2*b1**2 + a2**2*b1**2 + a3**2*b2**2 - a4**2*b2**2)
f(421) = dtau**2*(-a1**3*b1 - a2**3*b1 + a3**3*b2 + a4**3*b2)
f(422) = dtau*(a1**3*b1**2 - a2**3*b1**2 - a3**3*b2**2 + a4**3*b2**2)
f(423) = dtau**2*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(424) = dtau**4*(-a1*b2 - a2*b2 + a3*b1 + a4*b1)
f(425) = dtau*(-a1*b2**2 + a2*b2**2 + a3*b1**2 - a4*b1**2)
f(426) = dtau**3*(-a1*b2**2 + a2*b2**2 + a3*b1**2 - a4*b1**2)
f(427) = dtau**2*(a1*b2**3 + a2*b2**3 - a3*b1**3 - a4*b1**3)
f(428) = dtau*(a1*b2**4 - a2*b2**4 - a3*b1**4 + a4*b1**4)
f(429) = dtau**2*(-a1**2*b2 - a2**2*b2 + a3**2*b1 + a4**2*b1)
f(430) = dtau*(a1**2*b2**2 - a2**2*b2**2 - a3**2*b1**2 + a4**2*b1**2)
f(431) = dtau**2*(-a1**3*b2 - a2**3*b2 + a3**3*b1 + a4**3*b1)
f(432) = dtau*(-a1**3*b2**2 + a2**3*b2**2 + a3**3*b1**2 - a4**3*b1**2)
f(433) = b1*b2*dtau**2*(-b1 + b2)
v = sum(f*params)
end function c2h4_dipole_b3u_n3_d6_ADF


!###############################################################################


function c2h4_dipole_b3u_n4_d6_ADF(coords, params) result(v)
use adf
implicit none
type(adf_realq), intent(in) :: coords(12)
real(ark), intent(in) :: params(887)
type(adf_realq) :: v
type(adf_realq) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
type(adf_realq) :: f(887)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = r0*(b1*r1*r2 - b2*r3*r4)
f(2) = r0*(b1**3*r1*r2 - b2**3*r3*r4)
f(3) = r0*(b1*r1**2*r2 + b1*r1*r2**2 - b2*r3**2*r4 - b2*r3*r4**2)
f(4) = r0*(b1*r1**3*r2 + b1*r1*r2**3 - b2*r3**3*r4 - b2*r3*r4**3)
f(5) = r0*(b1*r1**2*r2**2 - b2*r3**2*r4**2)
f(6) = r0*(-b1*r1**3*r2 - b1*r1*r2**3 + b2*r3**3*r4 + b2*r3*r4**3)
f(7) = r0**2*(b1*r1*r2 - b2*r3*r4)
f(8) = r0**2*(-b1*r1**2*r2 - b1*r1*r2**2 + b2*r3**2*r4 + b2*r3*r4**2)
f(9) = r0**3*(-b1*r1*r2 + b2*r3*r4)
f(10) = r0*(b1*r3*r4 - b2*r1*r2)
f(11) = r0*(b1**3*r3*r4 - b2**3*r1*r2)
f(12) = r0*(b1*r3**2*r4 + b1*r3*r4**2 - b2*r1**2*r2 - b2*r1*r2**2)
f(13) = r0*(-b1*r3**3*r4 - b1*r3*r4**3 + b2*r1**3*r2 + b2*r1*r2**3)
f(14) = r0*(b1*r3**2*r4**2 - b2*r1**2*r2**2)
f(15) = r0**2*(b1*r3*r4 - b2*r1*r2)
f(16) = r0**2*(-b1*r3**2*r4 - b1*r3*r4**2 + b2*r1**2*r2 + b2*r1*r2**2)
f(17) = r0**3*(-b1*r3*r4 + b2*r1*r2)
f(18) = dtau*r0*(-r1**2*r2 + r1*r2**2 + r3**2*r4 - r3*r4**2)
f(19) = dtau*r0*(-r1**3*r2 + r1*r2**3 + r3**3*r4 - r3*r4**3)
f(20) = dtau*r0*(r1**3*r2 - r1*r2**3 - r3**3*r4 + r3*r4**3)
f(21) = dtau*r0**2*(r1**2*r2 - r1*r2**2 - r3**2*r4 + r3*r4**2)
f(22) = r0*(-b1*r1*r3 - b1*r2*r4 + b2*r1*r3 + b2*r2*r4)
f(23) = r0*(b1**3*r1*r3 + b1**3*r2*r4 - b2**3*r1*r3 - b2**3*r2*r4)
f(24) = r0*(-b1*r1*r3**2 - b1*r2*r4**2 + b2*r1**2*r3 + b2*r2**2*r4)
f(25) = r0*(b1*r1*r3**3 + b1*r2*r4**3 - b2*r1**3*r3 - b2*r2**3*r4)
f(26) = r0*(-b1*r1**2*r3 - b1*r2**2*r4 + b2*r1*r3**2 + b2*r2*r4**2)
f(27) = r0*(-b1*r1**2*r3**2 - b1*r2**2*r4**2 + b2*r1**2*r3**2 + b2*r2**2 &
      *r4**2)
f(28) = r0*(b1*r1**3*r3 + b1*r2**3*r4 - b2*r1*r3**3 - b2*r2*r4**3)
f(29) = r0**2*(-b1*r1*r3 - b1*r2*r4 + b2*r1*r3 + b2*r2*r4)
f(30) = r0**2*(-b1*r1*r3**2 - b1*r2*r4**2 + b2*r1**2*r3 + b2*r2**2*r4)
f(31) = r0**2*(-b1*r1**2*r3 - b1*r2**2*r4 + b2*r1*r3**2 + b2*r2*r4**2)
f(32) = r0**3*(b1*r1*r3 + b1*r2*r4 - b2*r1*r3 - b2*r2*r4)
f(33) = dtau*r0*(-r1**2*r3 + r1*r3**2 + r2**2*r4 - r2*r4**2)
f(34) = dtau*r0*(r1**3*r3 - r1*r3**3 - r2**3*r4 + r2*r4**3)
f(35) = dtau*r0*(-r1**3*r3 + r1*r3**3 + r2**3*r4 - r2*r4**3)
f(36) = dtau*r0**2*(r1**2*r3 - r1*r3**2 - r2**2*r4 + r2*r4**2)
f(37) = dtau*r0**2*(-r1**2*r3 + r1*r3**2 + r2**2*r4 - r2*r4**2)
f(38) = r0*(b1*r1*r4 + b1*r2*r3 - b2*r1*r4 - b2*r2*r3)
f(39) = r0*(-b1**3*r1*r4 - b1**3*r2*r3 + b2**3*r1*r4 + b2**3*r2*r3)
f(40) = r0*(b1*r1*r4**2 + b1*r2*r3**2 - b2*r1**2*r4 - b2*r2**2*r3)
f(41) = r0*(-b1*r1*r4**3 - b1*r2*r3**3 + b2*r1**3*r4 + b2*r2**3*r3)
f(42) = r0*(b1*r1**2*r4 + b1*r2**2*r3 - b2*r1*r4**2 - b2*r2*r3**2)
f(43) = r0*(-b1*r1**2*r4**2 - b1*r2**2*r3**2 + b2*r1**2*r4**2 + b2*r2**2 &
      *r3**2)
f(44) = r0*(-b1*r1**3*r4 - b1*r2**3*r3 + b2*r1*r4**3 + b2*r2*r3**3)
f(45) = r0**2*(b1*r1*r4 + b1*r2*r3 - b2*r1*r4 - b2*r2*r3)
f(46) = r0**2*(b1*r1*r4**2 + b1*r2*r3**2 - b2*r1**2*r4 - b2*r2**2*r3)
f(47) = r0**2*(b1*r1**2*r4 + b1*r2**2*r3 - b2*r1*r4**2 - b2*r2*r3**2)
f(48) = r0**3*(-b1*r1*r4 - b1*r2*r3 + b2*r1*r4 + b2*r2*r3)
f(49) = dtau*r0*(-r1*r4 + r2*r3)
f(50) = dtau**3*r0*(r1*r4 - r2*r3)
f(51) = dtau*r0*(-r1**2*r4 - r1*r4**2 + r2**2*r3 + r2*r3**2)
f(52) = dtau*r0*(r1**3*r4 + r1*r4**3 - r2**3*r3 - r2*r3**3)
f(53) = dtau*r0*(-r1**2*r4**2 + r2**2*r3**2)
f(54) = dtau*r0*(-r1**3*r4 - r1*r4**3 + r2**3*r3 + r2*r3**3)
f(55) = dtau*r0**2*(r1*r4 - r2*r3)
f(56) = dtau*r0**2*(r1**2*r4 + r1*r4**2 - r2**2*r3 - r2*r3**2)
f(57) = dtau*r0**3*(-r1*r4 + r2*r3)
f(58) = r0*(-a1*b1*r1 - a2*b1*r2 + a3*b2*r3 + a4*b2*r4)
f(59) = r0*(a1*b1**3*r1 + a2*b1**3*r2 - a3*b2**3*r3 - a4*b2**3*r4)
f(60) = r0*(a1**2*b1*r1 + a2**2*b1*r2 - a3**2*b2*r3 - a4**2*b2*r4)
f(61) = r0*(-a1**3*b1*r1 - a2**3*b1*r2 + a3**3*b2*r3 + a4**3*b2*r4)
f(62) = r0*(a1*b1*r1**2 + a2*b1*r2**2 - a3*b2*r3**2 - a4*b2*r4**2)
f(63) = r0*(a1**2*b1*r1**2 + a2**2*b1*r2**2 - a3**2*b2*r3**2 - a4**2*b2* &
      r4**2)
f(64) = r0*(-a1*b1*r1**3 - a2*b1*r2**3 + a3*b2*r3**3 + a4*b2*r4**3)
f(65) = r0**2*(a1*b1*r1 + a2*b1*r2 - a3*b2*r3 - a4*b2*r4)
f(66) = r0**2*(a1**2*b1*r1 + a2**2*b1*r2 - a3**2*b2*r3 - a4**2*b2*r4)
f(67) = r0**2*(a1*b1*r1**2 + a2*b1*r2**2 - a3*b2*r3**2 - a4*b2*r4**2)
f(68) = r0**3*(-a1*b1*r1 - a2*b1*r2 + a3*b2*r3 + a4*b2*r4)
f(69) = r0*(-a1*b2*r1 - a2*b2*r2 + a3*b1*r3 + a4*b1*r4)
f(70) = r0*(-a1*b2**3*r1 - a2*b2**3*r2 + a3*b1**3*r3 + a4*b1**3*r4)
f(71) = r0*(-a1**2*b2*r1 - a2**2*b2*r2 + a3**2*b1*r3 + a4**2*b1*r4)
f(72) = r0*(a1**3*b2*r1 + a2**3*b2*r2 - a3**3*b1*r3 - a4**3*b1*r4)
f(73) = r0*(-a1*b2*r1**2 - a2*b2*r2**2 + a3*b1*r3**2 + a4*b1*r4**2)
f(74) = r0*(-a1**2*b2*r1**2 - a2**2*b2*r2**2 + a3**2*b1*r3**2 + a4**2*b1 &
      *r4**2)
f(75) = r0*(a1*b2*r1**3 + a2*b2*r2**3 - a3*b1*r3**3 - a4*b1*r4**3)
f(76) = r0**2*(-a1*b2*r1 - a2*b2*r2 + a3*b1*r3 + a4*b1*r4)
f(77) = r0**2*(-a1**2*b2*r1 - a2**2*b2*r2 + a3**2*b1*r3 + a4**2*b1*r4)
f(78) = r0**2*(-a1*b2*r1**2 - a2*b2*r2**2 + a3*b1*r3**2 + a4*b1*r4**2)
f(79) = r0**3*(a1*b2*r1 + a2*b2*r2 - a3*b1*r3 - a4*b1*r4)
f(80) = dtau*r0*(-a1*r1 + a2*r2 + a3*r3 - a4*r4)
f(81) = dtau**3*r0*(a1*r1 - a2*r2 - a3*r3 + a4*r4)
f(82) = dtau*r0*(a1**2*r1 - a2**2*r2 - a3**2*r3 + a4**2*r4)
f(83) = dtau*r0*(a1**3*r1 - a2**3*r2 - a3**3*r3 + a4**3*r4)
f(84) = dtau*r0*(a1*r1**2 - a2*r2**2 - a3*r3**2 + a4*r4**2)
f(85) = dtau*r0*(a1**2*r1**2 - a2**2*r2**2 - a3**2*r3**2 + a4**2*r4**2)
f(86) = dtau*r0*(a1*r1**3 - a2*r2**3 - a3*r3**3 + a4*r4**3)
f(87) = dtau*r0**2*(a1*r1 - a2*r2 - a3*r3 + a4*r4)
f(88) = dtau*r0**2*(-a1**2*r1 + a2**2*r2 + a3**2*r3 - a4**2*r4)
f(89) = dtau*r0**2*(a1*r1**2 - a2*r2**2 - a3*r3**2 + a4*r4**2)
f(90) = dtau*r0**3*(a1*r1 - a2*r2 - a3*r3 + a4*r4)
f(91) = r0*(-a1*b1*r2 - a2*b1*r1 + a3*b2*r4 + a4*b2*r3)
f(92) = r0*(a1*b1**3*r2 + a2*b1**3*r1 - a3*b2**3*r4 - a4*b2**3*r3)
f(93) = r0*(-a1**2*b1*r2 - a2**2*b1*r1 + a3**2*b2*r4 + a4**2*b2*r3)
f(94) = r0*(a1**3*b1*r2 + a2**3*b1*r1 - a3**3*b2*r4 - a4**3*b2*r3)
f(95) = r0*(a1*b1*r2**2 + a2*b1*r1**2 - a3*b2*r4**2 - a4*b2*r3**2)
f(96) = r0*(-a1**2*b1*r2**2 - a2**2*b1*r1**2 + a3**2*b2*r4**2 + a4**2*b2 &
      *r3**2)
f(97) = r0*(-a1*b1*r2**3 - a2*b1*r1**3 + a3*b2*r4**3 + a4*b2*r3**3)
f(98) = r0**2*(a1*b1*r2 + a2*b1*r1 - a3*b2*r4 - a4*b2*r3)
f(99) = r0**2*(-a1**2*b1*r2 - a2**2*b1*r1 + a3**2*b2*r4 + a4**2*b2*r3)
f(100) = r0**2*(a1*b1*r2**2 + a2*b1*r1**2 - a3*b2*r4**2 - a4*b2*r3**2)
f(101) = r0**3*(a1*b1*r2 + a2*b1*r1 - a3*b2*r4 - a4*b2*r3)
f(102) = r0*(a1*b2*r2 + a2*b2*r1 - a3*b1*r4 - a4*b1*r3)
f(103) = r0*(-a1*b2**3*r2 - a2*b2**3*r1 + a3*b1**3*r4 + a4*b1**3*r3)
f(104) = r0*(a1**2*b2*r2 + a2**2*b2*r1 - a3**2*b1*r4 - a4**2*b1*r3)
f(105) = r0*(-a1**3*b2*r2 - a2**3*b2*r1 + a3**3*b1*r4 + a4**3*b1*r3)
f(106) = r0*(a1*b2*r2**2 + a2*b2*r1**2 - a3*b1*r4**2 - a4*b1*r3**2)
f(107) = r0*(a1**2*b2*r2**2 + a2**2*b2*r1**2 - a3**2*b1*r4**2 - a4**2*b1 &
      *r3**2)
f(108) = r0*(a1*b2*r2**3 + a2*b2*r1**3 - a3*b1*r4**3 - a4*b1*r3**3)
f(109) = r0**2*(a1*b2*r2 + a2*b2*r1 - a3*b1*r4 - a4*b1*r3)
f(110) = r0**2*(a1**2*b2*r2 + a2**2*b2*r1 - a3**2*b1*r4 - a4**2*b1*r3)
f(111) = r0**2*(-a1*b2*r2**2 - a2*b2*r1**2 + a3*b1*r4**2 + a4*b1*r3**2)
f(112) = r0**3*(a1*b2*r2 + a2*b2*r1 - a3*b1*r4 - a4*b1*r3)
f(113) = dtau*r0*(-a1*r2 + a2*r1 + a3*r4 - a4*r3)
f(114) = dtau**3*r0*(a1*r2 - a2*r1 - a3*r4 + a4*r3)
f(115) = dtau*r0*(-a1**2*r2 + a2**2*r1 + a3**2*r4 - a4**2*r3)
f(116) = dtau*r0*(-a1**3*r2 + a2**3*r1 + a3**3*r4 - a4**3*r3)
f(117) = dtau*r0*(-a1*r2**2 + a2*r1**2 + a3*r4**2 - a4*r3**2)
f(118) = dtau*r0*(-a1**2*r2**2 + a2**2*r1**2 + a3**2*r4**2 - a4**2*r3**2 &
      )
f(119) = dtau*r0*(-a1*r2**3 + a2*r1**3 + a3*r4**3 - a4*r3**3)
f(120) = dtau*r0**2*(-a1*r2 + a2*r1 + a3*r4 - a4*r3)
f(121) = dtau*r0**2*(-a1**2*r2 + a2**2*r1 + a3**2*r4 - a4**2*r3)
f(122) = dtau*r0**2*(-a1*r2**2 + a2*r1**2 + a3*r4**2 - a4*r3**2)
f(123) = dtau*r0**3*(a1*r2 - a2*r1 - a3*r4 + a4*r3)
f(124) = r0*(-a1*b2*r3 - a2*b2*r4 + a3*b1*r1 + a4*b1*r2)
f(125) = r0*(a1*b2**3*r3 + a2*b2**3*r4 - a3*b1**3*r1 - a4*b1**3*r2)
f(126) = r0*(-a1**2*b2*r3 - a2**2*b2*r4 + a3**2*b1*r1 + a4**2*b1*r2)
f(127) = r0*(a1**3*b2*r3 + a2**3*b2*r4 - a3**3*b1*r1 - a4**3*b1*r2)
f(128) = r0*(-a1*b2*r3**2 - a2*b2*r4**2 + a3*b1*r1**2 + a4*b1*r2**2)
f(129) = r0*(a1**2*b2*r3**2 + a2**2*b2*r4**2 - a3**2*b1*r1**2 - a4**2*b1 &
      *r2**2)
f(130) = r0*(-a1*b2*r3**3 - a2*b2*r4**3 + a3*b1*r1**3 + a4*b1*r2**3)
f(131) = r0**2*(-a1*b2*r3 - a2*b2*r4 + a3*b1*r1 + a4*b1*r2)
f(132) = r0**2*(a1**2*b2*r3 + a2**2*b2*r4 - a3**2*b1*r1 - a4**2*b1*r2)
f(133) = r0**2*(a1*b2*r3**2 + a2*b2*r4**2 - a3*b1*r1**2 - a4*b1*r2**2)
f(134) = r0**3*(-a1*b2*r3 - a2*b2*r4 + a3*b1*r1 + a4*b1*r2)
f(135) = r0*(a1*b1*r3 + a2*b1*r4 - a3*b2*r1 - a4*b2*r2)
f(136) = r0*(-a1*b1**3*r3 - a2*b1**3*r4 + a3*b2**3*r1 + a4*b2**3*r2)
f(137) = r0*(a1**2*b1*r3 + a2**2*b1*r4 - a3**2*b2*r1 - a4**2*b2*r2)
f(138) = r0*(-a1**3*b1*r3 - a2**3*b1*r4 + a3**3*b2*r1 + a4**3*b2*r2)
f(139) = r0*(-a1*b1*r3**2 - a2*b1*r4**2 + a3*b2*r1**2 + a4*b2*r2**2)
f(140) = r0*(-a1**2*b1*r3**2 - a2**2*b1*r4**2 + a3**2*b2*r1**2 + a4**2* &
      b2*r2**2)
f(141) = r0*(a1*b1*r3**3 + a2*b1*r4**3 - a3*b2*r1**3 - a4*b2*r2**3)
f(142) = r0**2*(-a1*b1*r3 - a2*b1*r4 + a3*b2*r1 + a4*b2*r2)
f(143) = r0**2*(-a1**2*b1*r3 - a2**2*b1*r4 + a3**2*b2*r1 + a4**2*b2*r2)
f(144) = r0**2*(-a1*b1*r3**2 - a2*b1*r4**2 + a3*b2*r1**2 + a4*b2*r2**2)
f(145) = r0**3*(-a1*b1*r3 - a2*b1*r4 + a3*b2*r1 + a4*b2*r2)
f(146) = dtau*r0*(a1*r3 - a2*r4 - a3*r1 + a4*r2)
f(147) = dtau**3*r0*(a1*r3 - a2*r4 - a3*r1 + a4*r2)
f(148) = dtau*r0*(-a1**2*r3 + a2**2*r4 + a3**2*r1 - a4**2*r2)
f(149) = dtau*r0*(a1**3*r3 - a2**3*r4 - a3**3*r1 + a4**3*r2)
f(150) = dtau*r0*(a1*r3**2 - a2*r4**2 - a3*r1**2 + a4*r2**2)
f(151) = dtau*r0*(a1**2*r3**2 - a2**2*r4**2 - a3**2*r1**2 + a4**2*r2**2)
f(152) = dtau*r0*(-a1*r3**3 + a2*r4**3 + a3*r1**3 - a4*r2**3)
f(153) = dtau*r0**2*(a1*r3 - a2*r4 - a3*r1 + a4*r2)
f(154) = dtau*r0**2*(-a1**2*r3 + a2**2*r4 + a3**2*r1 - a4**2*r2)
f(155) = dtau*r0**2*(a1*r3**2 - a2*r4**2 - a3*r1**2 + a4*r2**2)
f(156) = dtau*r0**3*(-a1*r3 + a2*r4 + a3*r1 - a4*r2)
f(157) = r0*(a1*b2*r4 + a2*b2*r3 - a3*b1*r2 - a4*b1*r1)
f(158) = r0*(a1*b2**3*r4 + a2*b2**3*r3 - a3*b1**3*r2 - a4*b1**3*r1)
f(159) = r0*(a1**2*b2*r4 + a2**2*b2*r3 - a3**2*b1*r2 - a4**2*b1*r1)
f(160) = r0*(-a1**3*b2*r4 - a2**3*b2*r3 + a3**3*b1*r2 + a4**3*b1*r1)
f(161) = r0*(a1*b2*r4**2 + a2*b2*r3**2 - a3*b1*r2**2 - a4*b1*r1**2)
f(162) = r0*(a1**2*b2*r4**2 + a2**2*b2*r3**2 - a3**2*b1*r2**2 - a4**2*b1 &
      *r1**2)
f(163) = r0*(-a1*b2*r4**3 - a2*b2*r3**3 + a3*b1*r2**3 + a4*b1*r1**3)
f(164) = r0**2*(a1*b2*r4 + a2*b2*r3 - a3*b1*r2 - a4*b1*r1)
f(165) = r0**2*(a1**2*b2*r4 + a2**2*b2*r3 - a3**2*b1*r2 - a4**2*b1*r1)
f(166) = r0**2*(a1*b2*r4**2 + a2*b2*r3**2 - a3*b1*r2**2 - a4*b1*r1**2)
f(167) = r0**3*(-a1*b2*r4 - a2*b2*r3 + a3*b1*r2 + a4*b1*r1)
f(168) = r0*(a1*b1*r4 + a2*b1*r3 - a3*b2*r2 - a4*b2*r1)
f(169) = r0*(-a1*b1**3*r4 - a2*b1**3*r3 + a3*b2**3*r2 + a4*b2**3*r1)
f(170) = r0*(-a1**2*b1*r4 - a2**2*b1*r3 + a3**2*b2*r2 + a4**2*b2*r1)
f(171) = r0*(a1**3*b1*r4 + a2**3*b1*r3 - a3**3*b2*r2 - a4**3*b2*r1)
f(172) = r0*(-a1*b1*r4**2 - a2*b1*r3**2 + a3*b2*r2**2 + a4*b2*r1**2)
f(173) = r0*(-a1**2*b1*r4**2 - a2**2*b1*r3**2 + a3**2*b2*r2**2 + a4**2* &
      b2*r1**2)
f(174) = r0*(a1*b1*r4**3 + a2*b1*r3**3 - a3*b2*r2**3 - a4*b2*r1**3)
f(175) = r0**2*(-a1*b1*r4 - a2*b1*r3 + a3*b2*r2 + a4*b2*r1)
f(176) = r0**2*(-a1**2*b1*r4 - a2**2*b1*r3 + a3**2*b2*r2 + a4**2*b2*r1)
f(177) = r0**2*(-a1*b1*r4**2 - a2*b1*r3**2 + a3*b2*r2**2 + a4*b2*r1**2)
f(178) = r0**3*(a1*b1*r4 + a2*b1*r3 - a3*b2*r2 - a4*b2*r1)
f(179) = dtau*r0*(-a1*r4 + a2*r3 + a3*r2 - a4*r1)
f(180) = dtau**3*r0*(-a1*r4 + a2*r3 + a3*r2 - a4*r1)
f(181) = dtau*r0*(-a1**2*r4 + a2**2*r3 + a3**2*r2 - a4**2*r1)
f(182) = dtau*r0*(-a1**3*r4 + a2**3*r3 + a3**3*r2 - a4**3*r1)
f(183) = dtau*r0*(-a1*r4**2 + a2*r3**2 + a3*r2**2 - a4*r1**2)
f(184) = dtau*r0*(-a1**2*r4**2 + a2**2*r3**2 + a3**2*r2**2 - a4**2*r1**2 &
      )
f(185) = dtau*r0*(a1*r4**3 - a2*r3**3 - a3*r2**3 + a4*r1**3)
f(186) = dtau*r0**2*(-a1*r4 + a2*r3 + a3*r2 - a4*r1)
f(187) = dtau*r0**2*(a1**2*r4 - a2**2*r3 - a3**2*r2 + a4**2*r1)
f(188) = dtau*r0**2*(-a1*r4**2 + a2*r3**2 + a3*r2**2 - a4*r1**2)
f(189) = dtau*r0**3*(a1*r4 - a2*r3 - a3*r2 + a4*r1)
f(190) = b1*b2*r0*(-b1*r3 - b1*r4 + b2*r1 + b2*r2)
f(191) = b1*b2*r0*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(192) = b1*b2*r0*(b1*r3**2 + b1*r4**2 - b2*r1**2 - b2*r2**2)
f(193) = b1*b2*r0*(-b1*r1**2 - b1*r2**2 + b2*r3**2 + b2*r4**2)
f(194) = b1*b2*r0**2*(b1*r3 + b1*r4 - b2*r1 - b2*r2)
f(195) = b1*b2*r0**2*(-b1*r1 - b1*r2 + b2*r3 + b2*r4)
f(196) = dtau**2*r0*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(197) = dtau*r0*(-b1**2*r1 + b1**2*r2 + b2**2*r3 - b2**2*r4)
f(198) = dtau**2*r0*(b1*r1**2 + b1*r2**2 - b2*r3**2 - b2*r4**2)
f(199) = dtau*r0*(b1**2*r1**2 - b1**2*r2**2 - b2**2*r3**2 + b2**2*r4**2)
f(200) = dtau**2*r0**2*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(201) = dtau*r0**2*(b1**2*r1 - b1**2*r2 - b2**2*r3 + b2**2*r4)
f(202) = dtau**2*r0*(b1*r3 + b1*r4 - b2*r1 - b2*r2)
f(203) = dtau*r0*(b1**2*r3 - b1**2*r4 - b2**2*r1 + b2**2*r2)
f(204) = dtau**2*r0*(b1*r3**2 + b1*r4**2 - b2*r1**2 - b2*r2**2)
f(205) = dtau*r0*(-b1**2*r3**2 + b1**2*r4**2 + b2**2*r1**2 - b2**2*r2**2 &
      )
f(206) = dtau**2*r0**2*(b1*r3 + b1*r4 - b2*r1 - b2*r2)
f(207) = dtau*r0**2*(-b1**2*r3 + b1**2*r4 + b2**2*r1 - b2**2*r2)
f(208) = r0*(a1*a2*b1 - a3*a4*b2)
f(209) = r0*(a1*a2*b1**3 - a3*a4*b2**3)
f(210) = r0*(-a1**2*a2*b1 - a1*a2**2*b1 + a3**2*a4*b2 + a3*a4**2*b2)
f(211) = r0*(-a1**3*a2*b1 - a1*a2**3*b1 + a3**3*a4*b2 + a3*a4**3*b2)
f(212) = r0*(a1**2*a2**2*b1 - a3**2*a4**2*b2)
f(213) = r0**2*(-a1*a2*b1 + a3*a4*b2)
f(214) = r0**2*(-a1**2*a2*b1 - a1*a2**2*b1 + a3**2*a4*b2 + a3*a4**2*b2)
f(215) = r0**3*(a1*a2*b1 - a3*a4*b2)
f(216) = r0*(-a1*a2*b2 + a3*a4*b1)
f(217) = r0*(-a1*a2*b2**3 + a3*a4*b1**3)
f(218) = r0*(a1**2*a2*b2 + a1*a2**2*b2 - a3**2*a4*b1 - a3*a4**2*b1)
f(219) = r0*(a1**3*a2*b2 + a1*a2**3*b2 - a3**3*a4*b1 - a3*a4**3*b1)
f(220) = r0*(-a1**2*a2**2*b2 + a3**2*a4**2*b1)
f(221) = r0**2*(a1*a2*b2 - a3*a4*b1)
f(222) = r0**2*(-a1**2*a2*b2 - a1*a2**2*b2 + a3**2*a4*b1 + a3*a4**2*b1)
f(223) = r0**3*(-a1*a2*b2 + a3*a4*b1)
f(224) = dtau*r0*(-a1**2*a2 + a1*a2**2 + a3**2*a4 - a3*a4**2)
f(225) = dtau*r0*(a1**3*a2 - a1*a2**3 - a3**3*a4 + a3*a4**3)
f(226) = dtau*r0**2*(a1**2*a2 - a1*a2**2 - a3**2*a4 + a3*a4**2)
f(227) = r0*(a1*a3*b1 - a1*a3*b2 + a2*a4*b1 - a2*a4*b2)
f(228) = r0*(a1*a3*b1**3 - a1*a3*b2**3 + a2*a4*b1**3 - a2*a4*b2**3)
f(229) = r0*(-a1**2*a3*b2 + a1*a3**2*b1 - a2**2*a4*b2 + a2*a4**2*b1)
f(230) = r0*(-a1**3*a3*b2 + a1*a3**3*b1 - a2**3*a4*b2 + a2*a4**3*b1)
f(231) = r0*(a1**2*a3*b1 - a1*a3**2*b2 + a2**2*a4*b1 - a2*a4**2*b2)
f(232) = r0*(a1**2*a3**2*b1 - a1**2*a3**2*b2 + a2**2*a4**2*b1 - a2**2*a4 &
      **2*b2)
f(233) = r0*(a1**3*a3*b1 - a1*a3**3*b2 + a2**3*a4*b1 - a2*a4**3*b2)
f(234) = r0**2*(a1*a3*b1 - a1*a3*b2 + a2*a4*b1 - a2*a4*b2)
f(235) = r0**2*(-a1**2*a3*b2 + a1*a3**2*b1 - a2**2*a4*b2 + a2*a4**2*b1)
f(236) = r0**2*(-a1**2*a3*b1 + a1*a3**2*b2 - a2**2*a4*b1 + a2*a4**2*b2)
f(237) = r0**3*(-a1*a3*b1 + a1*a3*b2 - a2*a4*b1 + a2*a4*b2)
f(238) = dtau*r0*(-a1**2*a3 + a1*a3**2 + a2**2*a4 - a2*a4**2)
f(239) = dtau*r0*(a1**3*a3 - a1*a3**3 - a2**3*a4 + a2*a4**3)
f(240) = dtau*r0**2*(-a1**2*a3 + a1*a3**2 + a2**2*a4 - a2*a4**2)
f(241) = r0*(a1*a4*b1 - a1*a4*b2 + a2*a3*b1 - a2*a3*b2)
f(242) = r0*(a1*a4*b1**3 - a1*a4*b2**3 + a2*a3*b1**3 - a2*a3*b2**3)
f(243) = r0*(-a1**2*a4*b2 + a1*a4**2*b1 - a2**2*a3*b2 + a2*a3**2*b1)
f(244) = r0*(-a1**3*a4*b2 + a1*a4**3*b1 - a2**3*a3*b2 + a2*a3**3*b1)
f(245) = r0*(a1**2*a4*b1 - a1*a4**2*b2 + a2**2*a3*b1 - a2*a3**2*b2)
f(246) = r0*(a1**2*a4**2*b1 - a1**2*a4**2*b2 + a2**2*a3**2*b1 - a2**2*a3 &
      **2*b2)
f(247) = r0*(a1**3*a4*b1 - a1*a4**3*b2 + a2**3*a3*b1 - a2*a3**3*b2)
f(248) = r0**2*(a1*a4*b1 - a1*a4*b2 + a2*a3*b1 - a2*a3*b2)
f(249) = r0**2*(-a1**2*a4*b2 + a1*a4**2*b1 - a2**2*a3*b2 + a2*a3**2*b1)
f(250) = r0**2*(a1**2*a4*b1 - a1*a4**2*b2 + a2**2*a3*b1 - a2*a3**2*b2)
f(251) = r0**3*(-a1*a4*b1 + a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
f(252) = dtau*r0*(a1*a4 - a2*a3)
f(253) = dtau**3*r0*(a1*a4 - a2*a3)
f(254) = dtau*r0*(-a1**2*a4 - a1*a4**2 + a2**2*a3 + a2*a3**2)
f(255) = dtau*r0*(a1**3*a4 + a1*a4**3 - a2**3*a3 - a2*a3**3)
f(256) = dtau*r0*(-a1**2*a4**2 + a2**2*a3**2)
f(257) = dtau*r0**2*(-a1*a4 + a2*a3)
f(258) = dtau*r0**2*(-a1**2*a4 - a1*a4**2 + a2**2*a3 + a2*a3**2)
f(259) = dtau*r0**3*(a1*a4 - a2*a3)
f(260) = b1*b2*r0*(-a1*b2 - a2*b2 + a3*b1 + a4*b1)
f(261) = b1*b2*r0*(-a1*b1 - a2*b1 + a3*b2 + a4*b2)
f(262) = b1*b2*r0*(a1**2*b2 + a2**2*b2 - a3**2*b1 - a4**2*b1)
f(263) = b1*b2*r0*(-a1**2*b1 - a2**2*b1 + a3**2*b2 + a4**2*b2)
f(264) = b1*b2*r0**2*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(265) = b1*b2*r0**2*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(266) = dtau**2*r0*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(267) = dtau*r0*(a1*b1**2 - a2*b1**2 - a3*b2**2 + a4*b2**2)
f(268) = dtau**2*r0*(a1**2*b1 + a2**2*b1 - a3**2*b2 - a4**2*b2)
f(269) = dtau*r0*(-a1**2*b1**2 + a2**2*b1**2 + a3**2*b2**2 - a4**2*b2**2 &
      )
f(270) = dtau**2*r0**2*(a1*b1 + a2*b1 - a3*b2 - a4*b2)
f(271) = dtau*r0**2*(-a1*b1**2 + a2*b1**2 + a3*b2**2 - a4*b2**2)
f(272) = dtau**2*r0*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(273) = dtau*r0*(a1*b2**2 - a2*b2**2 - a3*b1**2 + a4*b1**2)
f(274) = dtau**2*r0*(-a1**2*b2 - a2**2*b2 + a3**2*b1 + a4**2*b1)
f(275) = dtau*r0*(a1**2*b2**2 - a2**2*b2**2 - a3**2*b1**2 + a4**2*b1**2)
f(276) = dtau**2*r0**2*(a1*b2 + a2*b2 - a3*b1 - a4*b1)
f(277) = dtau*r0**2*(a1*b2**2 - a2*b2**2 - a3*b1**2 + a4*b1**2)
f(278) = b1*b2*dtau**2*r0*(-b1 + b2)
f(279) = b1*r1*r2*r3 + b1*r1*r2*r4 - b2*r1*r3*r4 - b2*r2*r3*r4
f(280) = -b1**3*r1*r2*r3 - b1**3*r1*r2*r4 + b2**3*r1*r3*r4 + b2**3*r2*r3 &
      *r4
f(281) = -b1*r1*r2*r3**2 - b1*r1*r2*r4**2 + b2*r1**2*r3*r4 + b2*r2**2*r3 &
      *r4
f(282) = b1*r1*r2*r3**3 + b1*r1*r2*r4**3 - b2*r1**3*r3*r4 - b2*r2**3*r3* &
      r4
f(283) = -b1*r1**2*r2*r4 - b1*r1*r2**2*r3 + b2*r1*r3*r4**2 + b2*r2*r3**2 &
      *r4
f(284) = b1*r1**2*r2*r4**2 + b1*r1*r2**2*r3**2 - b2*r1**2*r3*r4**2 - b2* &
      r2**2*r3**2*r4
f(285) = -b1*r1**3*r2*r4 - b1*r1*r2**3*r3 + b2*r1*r3*r4**3 + b2*r2*r3**3 &
      *r4
f(286) = -b1*r1**2*r2*r3 - b1*r1*r2**2*r4 + b2*r1*r3**2*r4 + b2*r2*r3*r4 &
      **2
f(287) = b1*r1**2*r2*r3**2 + b1*r1*r2**2*r4**2 - b2*r1**2*r3**2*r4 - b2* &
      r2**2*r3*r4**2
f(288) = b1*r1**2*r2**2*r3 + b1*r1**2*r2**2*r4 - b2*r1*r3**2*r4**2 - b2* &
      r2*r3**2*r4**2
f(289) = -b1*r1**3*r2*r3 - b1*r1*r2**3*r4 + b2*r1*r3**3*r4 + b2*r2*r3*r4 &
      **3
f(290) = b1*r1*r3*r4 + b1*r2*r3*r4 - b2*r1*r2*r3 - b2*r1*r2*r4
f(291) = -b1**3*r1*r3*r4 - b1**3*r2*r3*r4 + b2**3*r1*r2*r3 + b2**3*r1*r2 &
      *r4
f(292) = -b1*r1**2*r3*r4 - b1*r2**2*r3*r4 + b2*r1*r2*r3**2 + b2*r1*r2*r4 &
      **2
f(293) = b1*r1**3*r3*r4 + b1*r2**3*r3*r4 - b2*r1*r2*r3**3 - b2*r1*r2*r4 &
      **3
f(294) = -b1*r1*r3*r4**2 - b1*r2*r3**2*r4 + b2*r1**2*r2*r4 + b2*r1*r2**2 &
      *r3
f(295) = b1*r1**2*r3*r4**2 + b1*r2**2*r3**2*r4 - b2*r1**2*r2*r4**2 - b2* &
      r1*r2**2*r3**2
f(296) = -b1*r1*r3*r4**3 - b1*r2*r3**3*r4 + b2*r1**3*r2*r4 + b2*r1*r2**3 &
      *r3
f(297) = -b1*r1*r3**2*r4 - b1*r2*r3*r4**2 + b2*r1**2*r2*r3 + b2*r1*r2**2 &
      *r4
f(298) = b1*r1**2*r3**2*r4 + b1*r2**2*r3*r4**2 - b2*r1**2*r2*r3**2 - b2* &
      r1*r2**2*r4**2
f(299) = b1*r1*r3**2*r4**2 + b1*r2*r3**2*r4**2 - b2*r1**2*r2**2*r3 - b2* &
      r1**2*r2**2*r4
f(300) = -b1*r1*r3**3*r4 - b1*r2*r3*r4**3 + b2*r1**3*r2*r3 + b2*r1*r2**3 &
      *r4
f(301) = dtau*(r1*r2*r3 - r1*r2*r4 - r1*r3*r4 + r2*r3*r4)
f(302) = dtau**3*(r1*r2*r3 - r1*r2*r4 - r1*r3*r4 + r2*r3*r4)
f(303) = dtau*(-r1**2*r3*r4 + r1*r2*r3**2 - r1*r2*r4**2 + r2**2*r3*r4)
f(304) = dtau*(-r1**3*r3*r4 + r1*r2*r3**3 - r1*r2*r4**3 + r2**3*r3*r4)
f(305) = dtau*(-r1**2*r2*r4 + r1*r2**2*r3 - r1*r3*r4**2 + r2*r3**2*r4)
f(306) = dtau*(-r1**2*r2*r4**2 - r1**2*r3*r4**2 + r1*r2**2*r3**2 + r2**2 &
      *r3**2*r4)
f(307) = dtau*(r1**3*r2*r4 - r1*r2**3*r3 + r1*r3*r4**3 - r2*r3**3*r4)
f(308) = dtau*(r1**2*r2*r3 - r1*r2**2*r4 - r1*r3**2*r4 + r2*r3*r4**2)
f(309) = dtau*(r1**2*r2*r3**2 - r1**2*r3**2*r4 - r1*r2**2*r4**2 + r2**2* &
      r3*r4**2)
f(310) = dtau*(r1**2*r2**2*r3 - r1**2*r2**2*r4 - r1*r3**2*r4**2 + r2*r3 &
      **2*r4**2)
f(311) = dtau*(-r1**3*r2*r3 + r1*r2**3*r4 + r1*r3**3*r4 - r2*r3*r4**3)
f(312) = -a1*b1*r1*r2 - a2*b1*r1*r2 + a3*b2*r3*r4 + a4*b2*r3*r4
f(313) = -a1*b1**3*r1*r2 - a2*b1**3*r1*r2 + a3*b2**3*r3*r4 + a4*b2**3*r3 &
      *r4
f(314) = -a1**2*b1*r1*r2 - a2**2*b1*r1*r2 + a3**2*b2*r3*r4 + a4**2*b2*r3 &
      *r4
f(315) = -a1**3*b1*r1*r2 - a2**3*b1*r1*r2 + a3**3*b2*r3*r4 + a4**3*b2*r3 &
      *r4
f(316) = -a1*b1*r1*r2**2 - a2*b1*r1**2*r2 + a3*b2*r3*r4**2 + a4*b2*r3**2 &
      *r4
f(317) = -a1**2*b1*r1*r2**2 - a2**2*b1*r1**2*r2 + a3**2*b2*r3*r4**2 + a4 &
      **2*b2*r3**2*r4
f(318) = a1*b1*r1*r2**3 + a2*b1*r1**3*r2 - a3*b2*r3*r4**3 - a4*b2*r3**3* &
      r4
f(319) = -a1*b1*r1**2*r2 - a2*b1*r1*r2**2 + a3*b2*r3**2*r4 + a4*b2*r3*r4 &
      **2
f(320) = -a1**2*b1*r1**2*r2 - a2**2*b1*r1*r2**2 + a3**2*b2*r3**2*r4 + a4 &
      **2*b2*r3*r4**2
f(321) = a1*b1*r1**2*r2**2 + a2*b1*r1**2*r2**2 - a3*b2*r3**2*r4**2 - a4* &
      b2*r3**2*r4**2
f(322) = a1*b1*r1**3*r2 + a2*b1*r1*r2**3 - a3*b2*r3**3*r4 - a4*b2*r3*r4 &
      **3
f(323) = a1*b2*r1*r2 + a2*b2*r1*r2 - a3*b1*r3*r4 - a4*b1*r3*r4
f(324) = a1*b2**3*r1*r2 + a2*b2**3*r1*r2 - a3*b1**3*r3*r4 - a4*b1**3*r3* &
      r4
f(325) = a1**2*b2*r1*r2 + a2**2*b2*r1*r2 - a3**2*b1*r3*r4 - a4**2*b1*r3* &
      r4
f(326) = -a1**3*b2*r1*r2 - a2**3*b2*r1*r2 + a3**3*b1*r3*r4 + a4**3*b1*r3 &
      *r4
f(327) = a1*b2*r1*r2**2 + a2*b2*r1**2*r2 - a3*b1*r3*r4**2 - a4*b1*r3**2* &
      r4
f(328) = a1**2*b2*r1*r2**2 + a2**2*b2*r1**2*r2 - a3**2*b1*r3*r4**2 - a4 &
      **2*b1*r3**2*r4
f(329) = a1*b2*r1*r2**3 + a2*b2*r1**3*r2 - a3*b1*r3*r4**3 - a4*b1*r3**3* &
      r4
f(330) = a1*b2*r1**2*r2 + a2*b2*r1*r2**2 - a3*b1*r3**2*r4 - a4*b1*r3*r4 &
      **2
f(331) = a1**2*b2*r1**2*r2 + a2**2*b2*r1*r2**2 - a3**2*b1*r3**2*r4 - a4 &
      **2*b1*r3*r4**2
f(332) = -a1*b2*r1**2*r2**2 - a2*b2*r1**2*r2**2 + a3*b1*r3**2*r4**2 + a4 &
      *b1*r3**2*r4**2
f(333) = a1*b2*r1**3*r2 + a2*b2*r1*r2**3 - a3*b1*r3**3*r4 - a4*b1*r3*r4 &
      **3
f(334) = dtau*(a1*r1*r2 - a2*r1*r2 - a3*r3*r4 + a4*r3*r4)
f(335) = dtau**3*(-a1*r1*r2 + a2*r1*r2 + a3*r3*r4 - a4*r3*r4)
f(336) = dtau*(-a1**2*r1*r2 + a2**2*r1*r2 + a3**2*r3*r4 - a4**2*r3*r4)
f(337) = dtau*(a1**3*r1*r2 - a2**3*r1*r2 - a3**3*r3*r4 + a4**3*r3*r4)
f(338) = dtau*(-a1*r1*r2**2 + a2*r1**2*r2 + a3*r3*r4**2 - a4*r3**2*r4)
f(339) = dtau*(a1**2*r1*r2**2 - a2**2*r1**2*r2 - a3**2*r3*r4**2 + a4**2* &
      r3**2*r4)
f(340) = dtau*(-a1*r1*r2**3 + a2*r1**3*r2 + a3*r3*r4**3 - a4*r3**3*r4)
f(341) = dtau*(-a1*r1**2*r2 + a2*r1*r2**2 + a3*r3**2*r4 - a4*r3*r4**2)
f(342) = dtau*(a1**2*r1**2*r2 - a2**2*r1*r2**2 - a3**2*r3**2*r4 + a4**2* &
      r3*r4**2)
f(343) = dtau*(a1*r1**2*r2**2 - a2*r1**2*r2**2 - a3*r3**2*r4**2 + a4*r3 &
      **2*r4**2)
f(344) = dtau*(-a1*r1**3*r2 + a2*r1*r2**3 + a3*r3**3*r4 - a4*r3*r4**3)
f(345) = -a1*b2*r3*r4 - a2*b2*r3*r4 + a3*b1*r1*r2 + a4*b1*r1*r2
f(346) = -a1*b2**3*r3*r4 - a2*b2**3*r3*r4 + a3*b1**3*r1*r2 + a4*b1**3*r1 &
      *r2
f(347) = -a1**2*b2*r3*r4 - a2**2*b2*r3*r4 + a3**2*b1*r1*r2 + a4**2*b1*r1 &
      *r2
f(348) = -a1**3*b2*r3*r4 - a2**3*b2*r3*r4 + a3**3*b1*r1*r2 + a4**3*b1*r1 &
      *r2
f(349) = -a1*b2*r3*r4**2 - a2*b2*r3**2*r4 + a3*b1*r1*r2**2 + a4*b1*r1**2 &
      *r2
f(350) = -a1**2*b2*r3*r4**2 - a2**2*b2*r3**2*r4 + a3**2*b1*r1*r2**2 + a4 &
      **2*b1*r1**2*r2
f(351) = a1*b2*r3*r4**3 + a2*b2*r3**3*r4 - a3*b1*r1*r2**3 - a4*b1*r1**3* &
      r2
f(352) = -a1*b2*r3**2*r4 - a2*b2*r3*r4**2 + a3*b1*r1**2*r2 + a4*b1*r1*r2 &
      **2
f(353) = -a1**2*b2*r3**2*r4 - a2**2*b2*r3*r4**2 + a3**2*b1*r1**2*r2 + a4 &
      **2*b1*r1*r2**2
f(354) = a1*b2*r3**2*r4**2 + a2*b2*r3**2*r4**2 - a3*b1*r1**2*r2**2 - a4* &
      b1*r1**2*r2**2
f(355) = a1*b2*r3**3*r4 + a2*b2*r3*r4**3 - a3*b1*r1**3*r2 - a4*b1*r1*r2 &
      **3
f(356) = a1*b1*r3*r4 + a2*b1*r3*r4 - a3*b2*r1*r2 - a4*b2*r1*r2
f(357) = a1*b1**3*r3*r4 + a2*b1**3*r3*r4 - a3*b2**3*r1*r2 - a4*b2**3*r1* &
      r2
f(358) = a1**2*b1*r3*r4 + a2**2*b1*r3*r4 - a3**2*b2*r1*r2 - a4**2*b2*r1* &
      r2
f(359) = a1**3*b1*r3*r4 + a2**3*b1*r3*r4 - a3**3*b2*r1*r2 - a4**3*b2*r1* &
      r2
f(360) = a1*b1*r3*r4**2 + a2*b1*r3**2*r4 - a3*b2*r1*r2**2 - a4*b2*r1**2* &
      r2
f(361) = a1**2*b1*r3*r4**2 + a2**2*b1*r3**2*r4 - a3**2*b2*r1*r2**2 - a4 &
      **2*b2*r1**2*r2
f(362) = a1*b1*r3*r4**3 + a2*b1*r3**3*r4 - a3*b2*r1*r2**3 - a4*b2*r1**3* &
      r2
f(363) = a1*b1*r3**2*r4 + a2*b1*r3*r4**2 - a3*b2*r1**2*r2 - a4*b2*r1*r2 &
      **2
f(364) = a1**2*b1*r3**2*r4 + a2**2*b1*r3*r4**2 - a3**2*b2*r1**2*r2 - a4 &
      **2*b2*r1*r2**2
f(365) = -a1*b1*r3**2*r4**2 - a2*b1*r3**2*r4**2 + a3*b2*r1**2*r2**2 + a4 &
      *b2*r1**2*r2**2
f(366) = a1*b1*r3**3*r4 + a2*b1*r3*r4**3 - a3*b2*r1**3*r2 - a4*b2*r1*r2 &
      **3
f(367) = dtau*(a1*r3*r4 - a2*r3*r4 - a3*r1*r2 + a4*r1*r2)
f(368) = dtau**3*(-a1*r3*r4 + a2*r3*r4 + a3*r1*r2 - a4*r1*r2)
f(369) = dtau*(a1**2*r3*r4 - a2**2*r3*r4 - a3**2*r1*r2 + a4**2*r1*r2)
f(370) = dtau*(-a1**3*r3*r4 + a2**3*r3*r4 + a3**3*r1*r2 - a4**3*r1*r2)
f(371) = dtau*(-a1*r3*r4**2 + a2*r3**2*r4 + a3*r1*r2**2 - a4*r1**2*r2)
f(372) = dtau*(-a1**2*r3*r4**2 + a2**2*r3**2*r4 + a3**2*r1*r2**2 - a4**2 &
      *r1**2*r2)
f(373) = dtau*(-a1*r3*r4**3 + a2*r3**3*r4 + a3*r1*r2**3 - a4*r1**3*r2)
f(374) = dtau*(-a1*r3**2*r4 + a2*r3*r4**2 + a3*r1**2*r2 - a4*r1*r2**2)
f(375) = dtau*(-a1**2*r3**2*r4 + a2**2*r3*r4**2 + a3**2*r1**2*r2 - a4**2 &
      *r1*r2**2)
f(376) = dtau*(a1*r3**2*r4**2 - a2*r3**2*r4**2 - a3*r1**2*r2**2 + a4*r1 &
      **2*r2**2)
f(377) = dtau*(-a1*r3**3*r4 + a2*r3*r4**3 + a3*r1**3*r2 - a4*r1*r2**3)
f(378) = b1*b2*(b1*r3*r4 - b2*r1*r2)
f(379) = b1*b2*(-b1*r1*r2 + b2*r3*r4)
f(380) = b1*b2*(-b1*r3**2*r4 - b1*r3*r4**2 + b2*r1**2*r2 + b2*r1*r2**2)
f(381) = b1*b2*(-b1*r1**2*r2 - b1*r1*r2**2 + b2*r3**2*r4 + b2*r3*r4**2)
f(382) = dtau**2*(b1*r1*r2 - b2*r3*r4)
f(383) = dtau**2*(b1*r1**2*r2 + b1*r1*r2**2 - b2*r3**2*r4 - b2*r3*r4**2)
f(384) = dtau*(-b1**2*r1**2*r2 + b1**2*r1*r2**2 + b2**2*r3**2*r4 - b2**2 &
      *r3*r4**2)
f(385) = dtau*(b1**2*r1**2*r2 - b1**2*r1*r2**2 - b2**2*r3**2*r4 + b2**2* &
      r3*r4**2)
f(386) = dtau**2*(b1*r3*r4 - b2*r1*r2)
f(387) = dtau**2*(b1*r3**2*r4 + b1*r3*r4**2 - b2*r1**2*r2 - b2*r1*r2**2)
f(388) = dtau*(-b1**2*r3**2*r4 + b1**2*r3*r4**2 + b2**2*r1**2*r2 - b2**2 &
      *r1*r2**2)
f(389) = a1*b1*r1*r3 + a2*b1*r2*r4 - a3*b2*r1*r3 - a4*b2*r2*r4
f(390) = -a1*b1**3*r1*r3 - a2*b1**3*r2*r4 + a3*b2**3*r1*r3 + a4*b2**3*r2 &
      *r4
f(391) = -a1**2*b1*r1*r3 - a2**2*b1*r2*r4 + a3**2*b2*r1*r3 + a4**2*b2*r2 &
      *r4
f(392) = a1**3*b1*r1*r3 + a2**3*b1*r2*r4 - a3**3*b2*r1*r3 - a4**3*b2*r2* &
      r4
f(393) = a1*b1*r1*r3**2 + a2*b1*r2*r4**2 - a3*b2*r1**2*r3 - a4*b2*r2**2* &
      r4
f(394) = a1**2*b1*r1*r3**2 + a2**2*b1*r2*r4**2 - a3**2*b2*r1**2*r3 - a4 &
      **2*b2*r2**2*r4
f(395) = a1*b1*r1*r3**3 + a2*b1*r2*r4**3 - a3*b2*r1**3*r3 - a4*b2*r2**3* &
      r4
f(396) = a1*b1*r1**2*r3 + a2*b1*r2**2*r4 - a3*b2*r1*r3**2 - a4*b2*r2*r4 &
      **2
f(397) = a1**2*b1*r1**2*r3 + a2**2*b1*r2**2*r4 - a3**2*b2*r1*r3**2 - a4 &
      **2*b2*r2*r4**2
f(398) = a1*b1*r1**2*r3**2 + a2*b1*r2**2*r4**2 - a3*b2*r1**2*r3**2 - a4* &
      b2*r2**2*r4**2
f(399) = a1*b1*r1**3*r3 + a2*b1*r2**3*r4 - a3*b2*r1*r3**3 - a4*b2*r2*r4 &
      **3
f(400) = a1*b2*r1*r3 + a2*b2*r2*r4 - a3*b1*r1*r3 - a4*b1*r2*r4
f(401) = -a1*b2**3*r1*r3 - a2*b2**3*r2*r4 + a3*b1**3*r1*r3 + a4*b1**3*r2 &
      *r4
f(402) = a1**2*b2*r1*r3 + a2**2*b2*r2*r4 - a3**2*b1*r1*r3 - a4**2*b1*r2* &
      r4
f(403) = -a1**3*b2*r1*r3 - a2**3*b2*r2*r4 + a3**3*b1*r1*r3 + a4**3*b1*r2 &
      *r4
f(404) = a1*b2*r1*r3**2 + a2*b2*r2*r4**2 - a3*b1*r1**2*r3 - a4*b1*r2**2* &
      r4
f(405) = -a1**2*b2*r1*r3**2 - a2**2*b2*r2*r4**2 + a3**2*b1*r1**2*r3 + a4 &
      **2*b1*r2**2*r4
f(406) = a1*b2*r1*r3**3 + a2*b2*r2*r4**3 - a3*b1*r1**3*r3 - a4*b1*r2**3* &
      r4
f(407) = a1*b2*r1**2*r3 + a2*b2*r2**2*r4 - a3*b1*r1*r3**2 - a4*b1*r2*r4 &
      **2
f(408) = -a1**2*b2*r1**2*r3 - a2**2*b2*r2**2*r4 + a3**2*b1*r1*r3**2 + a4 &
      **2*b1*r2*r4**2
f(409) = a1*b2*r1**2*r3**2 + a2*b2*r2**2*r4**2 - a3*b1*r1**2*r3**2 - a4* &
      b1*r2**2*r4**2
f(410) = a1*b2*r1**3*r3 + a2*b2*r2**3*r4 - a3*b1*r1*r3**3 - a4*b1*r2*r4 &
      **3
f(411) = dtau*(a1*r1*r3 - a2*r2*r4 - a3*r1*r3 + a4*r2*r4)
f(412) = dtau**3*(-a1*r1*r3 + a2*r2*r4 + a3*r1*r3 - a4*r2*r4)
f(413) = dtau*(-a1**2*r1*r3 + a2**2*r2*r4 + a3**2*r1*r3 - a4**2*r2*r4)
f(414) = dtau*(-a1**3*r1*r3 + a2**3*r2*r4 + a3**3*r1*r3 - a4**3*r2*r4)
f(415) = dtau*(-a1*r1*r3**2 + a2*r2*r4**2 + a3*r1**2*r3 - a4*r2**2*r4)
f(416) = dtau*(-a1**2*r1*r3**2 + a2**2*r2*r4**2 + a3**2*r1**2*r3 - a4**2 &
      *r2**2*r4)
f(417) = dtau*(-a1*r1*r3**3 + a2*r2*r4**3 + a3*r1**3*r3 - a4*r2**3*r4)
f(418) = dtau*(-a1*r1**2*r3 + a2*r2**2*r4 + a3*r1*r3**2 - a4*r2*r4**2)
f(419) = dtau*(-a1**2*r1**2*r3 + a2**2*r2**2*r4 + a3**2*r1*r3**2 - a4**2 &
      *r2*r4**2)
f(420) = dtau*(a1*r1**2*r3**2 - a2*r2**2*r4**2 - a3*r1**2*r3**2 + a4*r2 &
      **2*r4**2)
f(421) = dtau*(a1*r1**3*r3 - a2*r2**3*r4 - a3*r1*r3**3 + a4*r2*r4**3)
f(422) = a1*b1*r2*r4 + a2*b1*r1*r3 - a3*b2*r2*r4 - a4*b2*r1*r3
f(423) = a1*b1**3*r2*r4 + a2*b1**3*r1*r3 - a3*b2**3*r2*r4 - a4*b2**3*r1* &
      r3
f(424) = a1**2*b1*r2*r4 + a2**2*b1*r1*r3 - a3**2*b2*r2*r4 - a4**2*b2*r1* &
      r3
f(425) = -a1**3*b1*r2*r4 - a2**3*b1*r1*r3 + a3**3*b2*r2*r4 + a4**3*b2*r1 &
      *r3
f(426) = -a1*b1*r2*r4**2 - a2*b1*r1*r3**2 + a3*b2*r2**2*r4 + a4*b2*r1**2 &
      *r3
f(427) = -a1**2*b1*r2*r4**2 - a2**2*b1*r1*r3**2 + a3**2*b2*r2**2*r4 + a4 &
      **2*b2*r1**2*r3
f(428) = a1*b1*r2*r4**3 + a2*b1*r1*r3**3 - a3*b2*r2**3*r4 - a4*b2*r1**3* &
      r3
f(429) = -a1*b1*r2**2*r4 - a2*b1*r1**2*r3 + a3*b2*r2*r4**2 + a4*b2*r1*r3 &
      **2
f(430) = -a1**2*b1*r2**2*r4 - a2**2*b1*r1**2*r3 + a3**2*b2*r2*r4**2 + a4 &
      **2*b2*r1*r3**2
f(431) = a1*b1*r2**2*r4**2 + a2*b1*r1**2*r3**2 - a3*b2*r2**2*r4**2 - a4* &
      b2*r1**2*r3**2
f(432) = a1*b1*r2**3*r4 + a2*b1*r1**3*r3 - a3*b2*r2*r4**3 - a4*b2*r1*r3 &
      **3
f(433) = a1*b2*r2*r4 + a2*b2*r1*r3 - a3*b1*r2*r4 - a4*b1*r1*r3
f(434) = a1*b2**3*r2*r4 + a2*b2**3*r1*r3 - a3*b1**3*r2*r4 - a4*b1**3*r1* &
      r3
f(435) = -a1**2*b2*r2*r4 - a2**2*b2*r1*r3 + a3**2*b1*r2*r4 + a4**2*b1*r1 &
      *r3
f(436) = a1**3*b2*r2*r4 + a2**3*b2*r1*r3 - a3**3*b1*r2*r4 - a4**3*b1*r1* &
      r3
f(437) = -a1*b2*r2*r4**2 - a2*b2*r1*r3**2 + a3*b1*r2**2*r4 + a4*b1*r1**2 &
      *r3
f(438) = a1**2*b2*r2*r4**2 + a2**2*b2*r1*r3**2 - a3**2*b1*r2**2*r4 - a4 &
      **2*b1*r1**2*r3
f(439) = -a1*b2*r2*r4**3 - a2*b2*r1*r3**3 + a3*b1*r2**3*r4 + a4*b1*r1**3 &
      *r3
f(440) = -a1*b2*r2**2*r4 - a2*b2*r1**2*r3 + a3*b1*r2*r4**2 + a4*b1*r1*r3 &
      **2
f(441) = a1**2*b2*r2**2*r4 + a2**2*b2*r1**2*r3 - a3**2*b1*r2*r4**2 - a4 &
      **2*b1*r1*r3**2
f(442) = a1*b2*r2**2*r4**2 + a2*b2*r1**2*r3**2 - a3*b1*r2**2*r4**2 - a4* &
      b1*r1**2*r3**2
f(443) = a1*b2*r2**3*r4 + a2*b2*r1**3*r3 - a3*b1*r2*r4**3 - a4*b1*r1*r3 &
      **3
f(444) = dtau*(a1*r2*r4 - a2*r1*r3 - a3*r2*r4 + a4*r1*r3)
f(445) = dtau**3*(a1*r2*r4 - a2*r1*r3 - a3*r2*r4 + a4*r1*r3)
f(446) = dtau*(a1**2*r2*r4 - a2**2*r1*r3 - a3**2*r2*r4 + a4**2*r1*r3)
f(447) = dtau*(a1**3*r2*r4 - a2**3*r1*r3 - a3**3*r2*r4 + a4**3*r1*r3)
f(448) = dtau*(a1*r2*r4**2 - a2*r1*r3**2 - a3*r2**2*r4 + a4*r1**2*r3)
f(449) = dtau*(a1**2*r2*r4**2 - a2**2*r1*r3**2 - a3**2*r2**2*r4 + a4**2* &
      r1**2*r3)
f(450) = dtau*(a1*r2*r4**3 - a2*r1*r3**3 - a3*r2**3*r4 + a4*r1**3*r3)
f(451) = dtau*(a1*r2**2*r4 - a2*r1**2*r3 - a3*r2*r4**2 + a4*r1*r3**2)
f(452) = dtau*(a1**2*r2**2*r4 - a2**2*r1**2*r3 - a3**2*r2*r4**2 + a4**2* &
      r1*r3**2)
f(453) = dtau*(a1*r2**2*r4**2 - a2*r1**2*r3**2 - a3*r2**2*r4**2 + a4*r1 &
      **2*r3**2)
f(454) = dtau*(-a1*r2**3*r4 + a2*r1**3*r3 + a3*r2*r4**3 - a4*r1*r3**3)
f(455) = b1*b2*(b1*r1*r3 + b1*r2*r4 - b2*r1*r3 - b2*r2*r4)
f(456) = b1*b2*(-b1*r1**2*r3 - b1*r2**2*r4 + b2*r1*r3**2 + b2*r2*r4**2)
f(457) = b1*b2*(-b1*r1*r3**2 - b1*r2*r4**2 + b2*r1**2*r3 + b2*r2**2*r4)
f(458) = dtau**2*(-b1*r1*r3 - b1*r2*r4 + b2*r1*r3 + b2*r2*r4)
f(459) = dtau*(b1**2*r1*r3 - b1**2*r2*r4 - b2**2*r1*r3 + b2**2*r2*r4)
f(460) = dtau**2*(-b1*r1*r3**2 - b1*r2*r4**2 + b2*r1**2*r3 + b2*r2**2*r4 &
      )
f(461) = dtau*(b1**2*r1*r3**2 - b1**2*r2*r4**2 - b2**2*r1**2*r3 + b2**2* &
      r2**2*r4)
f(462) = dtau**2*(-b1*r1**2*r3 - b1*r2**2*r4 + b2*r1*r3**2 + b2*r2*r4**2 &
      )
f(463) = dtau*(b1**2*r1**2*r3 - b1**2*r2**2*r4 - b2**2*r1*r3**2 + b2**2* &
      r2*r4**2)
f(464) = a1*b1*r1*r4 + a2*b1*r2*r3 - a3*b2*r2*r3 - a4*b2*r1*r4
f(465) = a1*b1**3*r1*r4 + a2*b1**3*r2*r3 - a3*b2**3*r2*r3 - a4*b2**3*r1* &
      r4
f(466) = a1**2*b1*r1*r4 + a2**2*b1*r2*r3 - a3**2*b2*r2*r3 - a4**2*b2*r1* &
      r4
f(467) = a1**3*b1*r1*r4 + a2**3*b1*r2*r3 - a3**3*b2*r2*r3 - a4**3*b2*r1* &
      r4
f(468) = a1*b1*r1*r4**2 + a2*b1*r2*r3**2 - a3*b2*r2**2*r3 - a4*b2*r1**2* &
      r4
f(469) = -a1**2*b1*r1*r4**2 - a2**2*b1*r2*r3**2 + a3**2*b2*r2**2*r3 + a4 &
      **2*b2*r1**2*r4
f(470) = -a1*b1*r1*r4**3 - a2*b1*r2*r3**3 + a3*b2*r2**3*r3 + a4*b2*r1**3 &
      *r4
f(471) = a1*b1*r1**2*r4 + a2*b1*r2**2*r3 - a3*b2*r2*r3**2 - a4*b2*r1*r4 &
      **2
f(472) = -a1**2*b1*r1**2*r4 - a2**2*b1*r2**2*r3 + a3**2*b2*r2*r3**2 + a4 &
      **2*b2*r1*r4**2
f(473) = a1*b1*r1**2*r4**2 + a2*b1*r2**2*r3**2 - a3*b2*r2**2*r3**2 - a4* &
      b2*r1**2*r4**2
f(474) = a1*b1*r1**3*r4 + a2*b1*r2**3*r3 - a3*b2*r2*r3**3 - a4*b2*r1*r4 &
      **3
f(475) = a1*b2*r1*r4 + a2*b2*r2*r3 - a3*b1*r2*r3 - a4*b1*r1*r4
f(476) = a1*b2**3*r1*r4 + a2*b2**3*r2*r3 - a3*b1**3*r2*r3 - a4*b1**3*r1* &
      r4
f(477) = a1**2*b2*r1*r4 + a2**2*b2*r2*r3 - a3**2*b1*r2*r3 - a4**2*b1*r1* &
      r4
f(478) = -a1**3*b2*r1*r4 - a2**3*b2*r2*r3 + a3**3*b1*r2*r3 + a4**3*b1*r1 &
      *r4
f(479) = a1*b2*r1*r4**2 + a2*b2*r2*r3**2 - a3*b1*r2**2*r3 - a4*b1*r1**2* &
      r4
f(480) = a1**2*b2*r1*r4**2 + a2**2*b2*r2*r3**2 - a3**2*b1*r2**2*r3 - a4 &
      **2*b1*r1**2*r4
f(481) = -a1*b2*r1*r4**3 - a2*b2*r2*r3**3 + a3*b1*r2**3*r3 + a4*b1*r1**3 &
      *r4
f(482) = a1*b2*r1**2*r4 + a2*b2*r2**2*r3 - a3*b1*r2*r3**2 - a4*b1*r1*r4 &
      **2
f(483) = a1**2*b2*r1**2*r4 + a2**2*b2*r2**2*r3 - a3**2*b1*r2*r3**2 - a4 &
      **2*b1*r1*r4**2
f(484) = a1*b2*r1**2*r4**2 + a2*b2*r2**2*r3**2 - a3*b1*r2**2*r3**2 - a4* &
      b1*r1**2*r4**2
f(485) = a1*b2*r1**3*r4 + a2*b2*r2**3*r3 - a3*b1*r2*r3**3 - a4*b1*r1*r4 &
      **3
f(486) = dtau*(a1*r1*r4 - a2*r2*r3 - a3*r2*r3 + a4*r1*r4)
f(487) = dtau**3*(-a1*r1*r4 + a2*r2*r3 + a3*r2*r3 - a4*r1*r4)
f(488) = dtau*(a1**2*r1*r4 - a2**2*r2*r3 - a3**2*r2*r3 + a4**2*r1*r4)
f(489) = dtau*(a1**3*r1*r4 - a2**3*r2*r3 - a3**3*r2*r3 + a4**3*r1*r4)
f(490) = dtau*(a1*r1*r4**2 - a2*r2*r3**2 - a3*r2**2*r3 + a4*r1**2*r4)
f(491) = dtau*(-a1**2*r1*r4**2 + a2**2*r2*r3**2 + a3**2*r2**2*r3 - a4**2 &
      *r1**2*r4)
f(492) = dtau*(a1*r1*r4**3 - a2*r2*r3**3 - a3*r2**3*r3 + a4*r1**3*r4)
f(493) = dtau*(a1*r1**2*r4 - a2*r2**2*r3 - a3*r2*r3**2 + a4*r1*r4**2)
f(494) = dtau*(-a1**2*r1**2*r4 + a2**2*r2**2*r3 + a3**2*r2*r3**2 - a4**2 &
      *r1*r4**2)
f(495) = dtau*(-a1*r1**2*r4**2 + a2*r2**2*r3**2 + a3*r2**2*r3**2 - a4*r1 &
      **2*r4**2)
f(496) = dtau*(-a1*r1**3*r4 + a2*r2**3*r3 + a3*r2*r3**3 - a4*r1*r4**3)
f(497) = a1*b1*r2*r3 + a2*b1*r1*r4 - a3*b2*r1*r4 - a4*b2*r2*r3
f(498) = -a1*b1**3*r2*r3 - a2*b1**3*r1*r4 + a3*b2**3*r1*r4 + a4*b2**3*r2 &
      *r3
f(499) = a1**2*b1*r2*r3 + a2**2*b1*r1*r4 - a3**2*b2*r1*r4 - a4**2*b2*r2* &
      r3
f(500) = -a1**3*b1*r2*r3 - a2**3*b1*r1*r4 + a3**3*b2*r1*r4 + a4**3*b2*r2 &
      *r3
f(501) = a1*b1*r2*r3**2 + a2*b1*r1*r4**2 - a3*b2*r1**2*r4 - a4*b2*r2**2* &
      r3
f(502) = a1**2*b1*r2*r3**2 + a2**2*b1*r1*r4**2 - a3**2*b2*r1**2*r4 - a4 &
      **2*b2*r2**2*r3
f(503) = a1*b1*r2*r3**3 + a2*b1*r1*r4**3 - a3*b2*r1**3*r4 - a4*b2*r2**3* &
      r3
f(504) = a1*b1*r2**2*r3 + a2*b1*r1**2*r4 - a3*b2*r1*r4**2 - a4*b2*r2*r3 &
      **2
f(505) = a1**2*b1*r2**2*r3 + a2**2*b1*r1**2*r4 - a3**2*b2*r1*r4**2 - a4 &
      **2*b2*r2*r3**2
f(506) = a1*b1*r2**2*r3**2 + a2*b1*r1**2*r4**2 - a3*b2*r1**2*r4**2 - a4* &
      b2*r2**2*r3**2
f(507) = -a1*b1*r2**3*r3 - a2*b1*r1**3*r4 + a3*b2*r1*r4**3 + a4*b2*r2*r3 &
      **3
f(508) = a1*b2*r2*r3 + a2*b2*r1*r4 - a3*b1*r1*r4 - a4*b1*r2*r3
f(509) = a1*b2**3*r2*r3 + a2*b2**3*r1*r4 - a3*b1**3*r1*r4 - a4*b1**3*r2* &
      r3
f(510) = -a1**2*b2*r2*r3 - a2**2*b2*r1*r4 + a3**2*b1*r1*r4 + a4**2*b1*r2 &
      *r3
f(511) = a1**3*b2*r2*r3 + a2**3*b2*r1*r4 - a3**3*b1*r1*r4 - a4**3*b1*r2* &
      r3
f(512) = a1*b2*r2*r3**2 + a2*b2*r1*r4**2 - a3*b1*r1**2*r4 - a4*b1*r2**2* &
      r3
f(513) = -a1**2*b2*r2*r3**2 - a2**2*b2*r1*r4**2 + a3**2*b1*r1**2*r4 + a4 &
      **2*b1*r2**2*r3
f(514) = a1*b2*r2*r3**3 + a2*b2*r1*r4**3 - a3*b1*r1**3*r4 - a4*b1*r2**3* &
      r3
f(515) = a1*b2*r2**2*r3 + a2*b2*r1**2*r4 - a3*b1*r1*r4**2 - a4*b1*r2*r3 &
      **2
f(516) = -a1**2*b2*r2**2*r3 - a2**2*b2*r1**2*r4 + a3**2*b1*r1*r4**2 + a4 &
      **2*b1*r2*r3**2
f(517) = a1*b2*r2**2*r3**2 + a2*b2*r1**2*r4**2 - a3*b1*r1**2*r4**2 - a4* &
      b1*r2**2*r3**2
f(518) = -a1*b2*r2**3*r3 - a2*b2*r1**3*r4 + a3*b1*r1*r4**3 + a4*b1*r2*r3 &
      **3
f(519) = dtau*(a1*r2*r3 - a2*r1*r4 - a3*r1*r4 + a4*r2*r3)
f(520) = dtau**3*(a1*r2*r3 - a2*r1*r4 - a3*r1*r4 + a4*r2*r3)
f(521) = dtau*(-a1**2*r2*r3 + a2**2*r1*r4 + a3**2*r1*r4 - a4**2*r2*r3)
f(522) = dtau*(-a1**3*r2*r3 + a2**3*r1*r4 + a3**3*r1*r4 - a4**3*r2*r3)
f(523) = dtau*(-a1*r2*r3**2 + a2*r1*r4**2 + a3*r1**2*r4 - a4*r2**2*r3)
f(524) = dtau*(a1**2*r2*r3**2 - a2**2*r1*r4**2 - a3**2*r1**2*r4 + a4**2* &
      r2**2*r3)
f(525) = dtau*(-a1*r2*r3**3 + a2*r1*r4**3 + a3*r1**3*r4 - a4*r2**3*r3)
f(526) = dtau*(-a1*r2**2*r3 + a2*r1**2*r4 + a3*r1*r4**2 - a4*r2*r3**2)
f(527) = dtau*(a1**2*r2**2*r3 - a2**2*r1**2*r4 - a3**2*r1*r4**2 + a4**2* &
      r2*r3**2)
f(528) = dtau*(-a1*r2**2*r3**2 + a2*r1**2*r4**2 + a3*r1**2*r4**2 - a4*r2 &
      **2*r3**2)
f(529) = dtau*(a1*r2**3*r3 - a2*r1**3*r4 - a3*r1*r4**3 + a4*r2*r3**3)
f(530) = b1*b2*(-b1*r1*r4 - b1*r2*r3 + b2*r1*r4 + b2*r2*r3)
f(531) = b1*b2*(b1*r1**2*r4 + b1*r2**2*r3 - b2*r1*r4**2 - b2*r2*r3**2)
f(532) = b1*b2*(b1*r1*r4**2 + b1*r2*r3**2 - b2*r1**2*r4 - b2*r2**2*r3)
f(533) = dtau**2*(b1*r1*r4 + b1*r2*r3 - b2*r1*r4 - b2*r2*r3)
f(534) = dtau*(-b1**2*r1*r4 + b1**2*r2*r3 - b2**2*r1*r4 + b2**2*r2*r3)
f(535) = dtau**2*(b1*r1*r4**2 + b1*r2*r3**2 - b2*r1**2*r4 - b2*r2**2*r3)
f(536) = dtau*(-b1**2*r1*r4**2 + b1**2*r2*r3**2 - b2**2*r1**2*r4 + b2**2 &
      *r2**2*r3)
f(537) = dtau**2*(b1*r1**2*r4 + b1*r2**2*r3 - b2*r1*r4**2 - b2*r2*r3**2)
f(538) = dtau*(-b1**2*r1**2*r4 + b1**2*r2**2*r3 - b2**2*r1*r4**2 + b2**2 &
      *r2*r3**2)
f(539) = -a1*a2*b1*r1 - a1*a2*b1*r2 + a3*a4*b2*r3 + a3*a4*b2*r4
f(540) = a1*a2*b1**3*r1 + a1*a2*b1**3*r2 - a3*a4*b2**3*r3 - a3*a4*b2**3* &
      r4
f(541) = a1**2*a2*b1*r2 + a1*a2**2*b1*r1 - a3**2*a4*b2*r4 - a3*a4**2*b2* &
      r3
f(542) = a1**3*a2*b1*r2 + a1*a2**3*b1*r1 - a3**3*a4*b2*r4 - a3*a4**3*b2* &
      r3
f(543) = a1**2*a2*b1*r1 + a1*a2**2*b1*r2 - a3**2*a4*b2*r3 - a3*a4**2*b2* &
      r4
f(544) = a1**2*a2**2*b1*r1 + a1**2*a2**2*b1*r2 - a3**2*a4**2*b2*r3 - a3 &
      **2*a4**2*b2*r4
f(545) = a1**3*a2*b1*r1 + a1*a2**3*b1*r2 - a3**3*a4*b2*r3 - a3*a4**3*b2* &
      r4
f(546) = a1*a2*b1*r1**2 + a1*a2*b1*r2**2 - a3*a4*b2*r3**2 - a3*a4*b2*r4 &
      **2
f(547) = a1**2*a2*b1*r2**2 + a1*a2**2*b1*r1**2 - a3**2*a4*b2*r4**2 - a3* &
      a4**2*b2*r3**2
f(548) = a1**2*a2*b1*r1**2 + a1*a2**2*b1*r2**2 - a3**2*a4*b2*r3**2 - a3* &
      a4**2*b2*r4**2
f(549) = a1*a2*b1*r1**3 + a1*a2*b1*r2**3 - a3*a4*b2*r3**3 - a3*a4*b2*r4 &
      **3
f(550) = -a1*a2*b2*r1 - a1*a2*b2*r2 + a3*a4*b1*r3 + a3*a4*b1*r4
f(551) = a1*a2*b2**3*r1 + a1*a2*b2**3*r2 - a3*a4*b1**3*r3 - a3*a4*b1**3* &
      r4
f(552) = -a1**2*a2*b2*r2 - a1*a2**2*b2*r1 + a3**2*a4*b1*r4 + a3*a4**2*b1 &
      *r3
f(553) = a1**3*a2*b2*r2 + a1*a2**3*b2*r1 - a3**3*a4*b1*r4 - a3*a4**3*b1* &
      r3
f(554) = a1**2*a2*b2*r1 + a1*a2**2*b2*r2 - a3**2*a4*b1*r3 - a3*a4**2*b1* &
      r4
f(555) = -a1**2*a2**2*b2*r1 - a1**2*a2**2*b2*r2 + a3**2*a4**2*b1*r3 + a3 &
      **2*a4**2*b1*r4
f(556) = -a1**3*a2*b2*r1 - a1*a2**3*b2*r2 + a3**3*a4*b1*r3 + a3*a4**3*b1 &
      *r4
f(557) = -a1*a2*b2*r1**2 - a1*a2*b2*r2**2 + a3*a4*b1*r3**2 + a3*a4*b1*r4 &
      **2
f(558) = a1**2*a2*b2*r2**2 + a1*a2**2*b2*r1**2 - a3**2*a4*b1*r4**2 - a3* &
      a4**2*b1*r3**2
f(559) = a1**2*a2*b2*r1**2 + a1*a2**2*b2*r2**2 - a3**2*a4*b1*r3**2 - a3* &
      a4**2*b1*r4**2
f(560) = a1*a2*b2*r1**3 + a1*a2*b2*r2**3 - a3*a4*b1*r3**3 - a3*a4*b1*r4 &
      **3
f(561) = dtau*(a1*a2*r1 - a1*a2*r2 - a3*a4*r3 + a3*a4*r4)
f(562) = dtau**3*(-a1*a2*r1 + a1*a2*r2 + a3*a4*r3 - a3*a4*r4)
f(563) = dtau*(-a1**2*a2*r2 + a1*a2**2*r1 + a3**2*a4*r4 - a3*a4**2*r3)
f(564) = dtau*(-a1**3*a2*r2 + a1*a2**3*r1 + a3**3*a4*r4 - a3*a4**3*r3)
f(565) = dtau*(-a1**2*a2*r1 + a1*a2**2*r2 + a3**2*a4*r3 - a3*a4**2*r4)
f(566) = dtau*(-a1**2*a2**2*r1 + a1**2*a2**2*r2 + a3**2*a4**2*r3 - a3**2 &
      *a4**2*r4)
f(567) = dtau*(-a1**3*a2*r1 + a1*a2**3*r2 + a3**3*a4*r3 - a3*a4**3*r4)
f(568) = dtau*(a1*a2*r1**2 - a1*a2*r2**2 - a3*a4*r3**2 + a3*a4*r4**2)
f(569) = dtau*(-a1**2*a2*r2**2 + a1*a2**2*r1**2 + a3**2*a4*r4**2 - a3*a4 &
      **2*r3**2)
f(570) = dtau*(-a1**2*a2*r1**2 + a1*a2**2*r2**2 + a3**2*a4*r3**2 - a3*a4 &
      **2*r4**2)
f(571) = dtau*(a1*a2*r1**3 - a1*a2*r2**3 - a3*a4*r3**3 + a3*a4*r4**3)
f(572) = -a1*a3*b1*r1 + a1*a3*b2*r3 - a2*a4*b1*r2 + a2*a4*b2*r4
f(573) = -a1*a3*b1**3*r1 + a1*a3*b2**3*r3 - a2*a4*b1**3*r2 + a2*a4*b2**3 &
      *r4
f(574) = -a1**2*a3*b2*r3 + a1*a3**2*b1*r1 - a2**2*a4*b2*r4 + a2*a4**2*b1 &
      *r2
f(575) = -a1**3*a3*b2*r3 + a1*a3**3*b1*r1 - a2**3*a4*b2*r4 + a2*a4**3*b1 &
      *r2
f(576) = -a1**2*a3*b1*r1 + a1*a3**2*b2*r3 - a2**2*a4*b1*r2 + a2*a4**2*b2 &
      *r4
f(577) = -a1**2*a3**2*b1*r1 + a1**2*a3**2*b2*r3 - a2**2*a4**2*b1*r2 + a2 &
      **2*a4**2*b2*r4
f(578) = -a1**3*a3*b1*r1 + a1*a3**3*b2*r3 - a2**3*a4*b1*r2 + a2*a4**3*b2 &
      *r4
f(579) = a1*a3*b1*r1**2 - a1*a3*b2*r3**2 + a2*a4*b1*r2**2 - a2*a4*b2*r4 &
      **2
f(580) = -a1**2*a3*b2*r3**2 + a1*a3**2*b1*r1**2 - a2**2*a4*b2*r4**2 + a2 &
      *a4**2*b1*r2**2
f(581) = a1**2*a3*b1*r1**2 - a1*a3**2*b2*r3**2 + a2**2*a4*b1*r2**2 - a2* &
      a4**2*b2*r4**2
f(582) = a1*a3*b1*r1**3 - a1*a3*b2*r3**3 + a2*a4*b1*r2**3 - a2*a4*b2*r4 &
      **3
f(583) = a1*a3*b1*r3 - a1*a3*b2*r1 + a2*a4*b1*r4 - a2*a4*b2*r2
f(584) = -a1*a3*b1**3*r3 + a1*a3*b2**3*r1 - a2*a4*b1**3*r4 + a2*a4*b2**3 &
      *r2
f(585) = -a1**2*a3*b1*r3 + a1*a3**2*b2*r1 - a2**2*a4*b1*r4 + a2*a4**2*b2 &
      *r2
f(586) = -a1**3*a3*b1*r3 + a1*a3**3*b2*r1 - a2**3*a4*b1*r4 + a2*a4**3*b2 &
      *r2
f(587) = -a1**2*a3*b2*r1 + a1*a3**2*b1*r3 - a2**2*a4*b2*r2 + a2*a4**2*b1 &
      *r4
f(588) = -a1**2*a3**2*b1*r3 + a1**2*a3**2*b2*r1 - a2**2*a4**2*b1*r4 + a2 &
      **2*a4**2*b2*r2
f(589) = a1**3*a3*b2*r1 - a1*a3**3*b1*r3 + a2**3*a4*b2*r2 - a2*a4**3*b1* &
      r4
f(590) = -a1*a3*b1*r3**2 + a1*a3*b2*r1**2 - a2*a4*b1*r4**2 + a2*a4*b2*r2 &
      **2
f(591) = -a1**2*a3*b1*r3**2 + a1*a3**2*b2*r1**2 - a2**2*a4*b1*r4**2 + a2 &
      *a4**2*b2*r2**2
f(592) = a1**2*a3*b2*r1**2 - a1*a3**2*b1*r3**2 + a2**2*a4*b2*r2**2 - a2* &
      a4**2*b1*r4**2
f(593) = a1*a3*b1*r3**3 - a1*a3*b2*r1**3 + a2*a4*b1*r4**3 - a2*a4*b2*r2 &
      **3
f(594) = dtau*(-a1*a3*r1 + a1*a3*r3 + a2*a4*r2 - a2*a4*r4)
f(595) = dtau**3*(a1*a3*r1 - a1*a3*r3 - a2*a4*r2 + a2*a4*r4)
f(596) = dtau*(-a1**2*a3*r3 + a1*a3**2*r1 + a2**2*a4*r4 - a2*a4**2*r2)
f(597) = dtau*(a1**3*a3*r3 - a1*a3**3*r1 - a2**3*a4*r4 + a2*a4**3*r2)
f(598) = dtau*(-a1**2*a3*r1 + a1*a3**2*r3 + a2**2*a4*r2 - a2*a4**2*r4)
f(599) = dtau*(a1**2*a3**2*r1 - a1**2*a3**2*r3 - a2**2*a4**2*r2 + a2**2* &
      a4**2*r4)
f(600) = dtau*(-a1**3*a3*r1 + a1*a3**3*r3 + a2**3*a4*r2 - a2*a4**3*r4)
f(601) = dtau*(a1*a3*r1**2 - a1*a3*r3**2 - a2*a4*r2**2 + a2*a4*r4**2)
f(602) = dtau*(a1**2*a3*r3**2 - a1*a3**2*r1**2 - a2**2*a4*r4**2 + a2*a4 &
      **2*r2**2)
f(603) = dtau*(a1**2*a3*r1**2 - a1*a3**2*r3**2 - a2**2*a4*r2**2 + a2*a4 &
      **2*r4**2)
f(604) = dtau*(-a1*a3*r1**3 + a1*a3*r3**3 + a2*a4*r2**3 - a2*a4*r4**3)
f(605) = -a1*a4*b1*r1 + a1*a4*b2*r4 - a2*a3*b1*r2 + a2*a3*b2*r3
f(606) = -a1*a4*b1**3*r1 + a1*a4*b2**3*r4 - a2*a3*b1**3*r2 + a2*a3*b2**3 &
      *r3
f(607) = -a1**2*a4*b2*r4 + a1*a4**2*b1*r1 - a2**2*a3*b2*r3 + a2*a3**2*b1 &
      *r2
f(608) = -a1**3*a4*b2*r4 + a1*a4**3*b1*r1 - a2**3*a3*b2*r3 + a2*a3**3*b1 &
      *r2
f(609) = -a1**2*a4*b1*r1 + a1*a4**2*b2*r4 - a2**2*a3*b1*r2 + a2*a3**2*b2 &
      *r3
f(610) = a1**2*a4**2*b1*r1 - a1**2*a4**2*b2*r4 + a2**2*a3**2*b1*r2 - a2 &
      **2*a3**2*b2*r3
f(611) = -a1**3*a4*b1*r1 + a1*a4**3*b2*r4 - a2**3*a3*b1*r2 + a2*a3**3*b2 &
      *r3
f(612) = a1*a4*b1*r1**2 - a1*a4*b2*r4**2 + a2*a3*b1*r2**2 - a2*a3*b2*r3 &
      **2
f(613) = a1**2*a4*b2*r4**2 - a1*a4**2*b1*r1**2 + a2**2*a3*b2*r3**2 - a2* &
      a3**2*b1*r2**2
f(614) = -a1**2*a4*b1*r1**2 + a1*a4**2*b2*r4**2 - a2**2*a3*b1*r2**2 + a2 &
      *a3**2*b2*r3**2
f(615) = a1*a4*b1*r1**3 - a1*a4*b2*r4**3 + a2*a3*b1*r2**3 - a2*a3*b2*r3 &
      **3
f(616) = -a1*a4*b1*r4 + a1*a4*b2*r1 - a2*a3*b1*r3 + a2*a3*b2*r2
f(617) = -a1*a4*b1**3*r4 + a1*a4*b2**3*r1 - a2*a3*b1**3*r3 + a2*a3*b2**3 &
      *r2
f(618) = -a1**2*a4*b1*r4 + a1*a4**2*b2*r1 - a2**2*a3*b1*r3 + a2*a3**2*b2 &
      *r2
f(619) = -a1**3*a4*b1*r4 + a1*a4**3*b2*r1 - a2**3*a3*b1*r3 + a2*a3**3*b2 &
      *r2
f(620) = -a1**2*a4*b2*r1 + a1*a4**2*b1*r4 - a2**2*a3*b2*r2 + a2*a3**2*b1 &
      *r3
f(621) = a1**2*a4**2*b1*r4 - a1**2*a4**2*b2*r1 + a2**2*a3**2*b1*r3 - a2 &
      **2*a3**2*b2*r2
f(622) = a1**3*a4*b2*r1 - a1*a4**3*b1*r4 + a2**3*a3*b2*r2 - a2*a3**3*b1* &
      r3
f(623) = -a1*a4*b1*r4**2 + a1*a4*b2*r1**2 - a2*a3*b1*r3**2 + a2*a3*b2*r2 &
      **2
f(624) = -a1**2*a4*b1*r4**2 + a1*a4**2*b2*r1**2 - a2**2*a3*b1*r3**2 + a2 &
      *a3**2*b2*r2**2
f(625) = -a1**2*a4*b2*r1**2 + a1*a4**2*b1*r4**2 - a2**2*a3*b2*r2**2 + a2 &
      *a3**2*b1*r3**2
f(626) = -a1*a4*b1*r4**3 + a1*a4*b2*r1**3 - a2*a3*b1*r3**3 + a2*a3*b2*r2 &
      **3
f(627) = dtau*(-a1*a4*r1 - a1*a4*r4 + a2*a3*r2 + a2*a3*r3)
f(628) = dtau**3*(a1*a4*r1 + a1*a4*r4 - a2*a3*r2 - a2*a3*r3)
f(629) = dtau*(-a1**2*a4*r4 - a1*a4**2*r1 + a2**2*a3*r3 + a2*a3**2*r2)
f(630) = dtau*(a1**3*a4*r4 + a1*a4**3*r1 - a2**3*a3*r3 - a2*a3**3*r2)
f(631) = dtau*(a1**2*a4*r1 + a1*a4**2*r4 - a2**2*a3*r2 - a2*a3**2*r3)
f(632) = dtau*(a1**2*a4**2*r1 + a1**2*a4**2*r4 - a2**2*a3**2*r2 - a2**2* &
      a3**2*r3)
f(633) = dtau*(a1**3*a4*r1 + a1*a4**3*r4 - a2**3*a3*r2 - a2*a3**3*r3)
f(634) = dtau*(a1*a4*r1**2 + a1*a4*r4**2 - a2*a3*r2**2 - a2*a3*r3**2)
f(635) = dtau*(a1**2*a4*r4**2 + a1*a4**2*r1**2 - a2**2*a3*r3**2 - a2*a3 &
      **2*r2**2)
f(636) = dtau*(a1**2*a4*r1**2 + a1*a4**2*r4**2 - a2**2*a3*r2**2 - a2*a3 &
      **2*r3**2)
f(637) = dtau*(a1*a4*r1**3 + a1*a4*r4**3 - a2*a3*r2**3 - a2*a3*r3**3)
f(638) = b1*b2*(a1*b2*r1 + a2*b2*r2 - a3*b1*r3 - a4*b1*r4)
f(639) = b1*b2*(a1*b1*r1 + a2*b1*r2 - a3*b2*r3 - a4*b2*r4)
f(640) = b1*b2*(-a1**2*b2*r1 - a2**2*b2*r2 + a3**2*b1*r3 + a4**2*b1*r4)
f(641) = b1*b2*(a1**2*b1*r1 + a2**2*b1*r2 - a3**2*b2*r3 - a4**2*b2*r4)
f(642) = b1*b2*(-a1*b2*r1**2 - a2*b2*r2**2 + a3*b1*r3**2 + a4*b1*r4**2)
f(643) = b1*b2*(a1*b1*r1**2 + a2*b1*r2**2 - a3*b2*r3**2 - a4*b2*r4**2)
f(644) = dtau**2*(a1*b1*r1 + a2*b1*r2 - a3*b2*r3 - a4*b2*r4)
f(645) = dtau*(a1*b1**2*r1 - a2*b1**2*r2 - a3*b2**2*r3 + a4*b2**2*r4)
f(646) = dtau**2*(a1**2*b1*r1 + a2**2*b1*r2 - a3**2*b2*r3 - a4**2*b2*r4)
f(647) = dtau*(a1**2*b1**2*r1 - a2**2*b1**2*r2 - a3**2*b2**2*r3 + a4**2* &
      b2**2*r4)
f(648) = dtau**2*(a1*b1*r1**2 + a2*b1*r2**2 - a3*b2*r3**2 - a4*b2*r4**2)
f(649) = dtau*(-a1*b1**2*r1**2 + a2*b1**2*r2**2 + a3*b2**2*r3**2 - a4*b2 &
      **2*r4**2)
f(650) = dtau**2*(-a1*b2*r1 - a2*b2*r2 + a3*b1*r3 + a4*b1*r4)
f(651) = dtau*(a1*b2**2*r1 - a2*b2**2*r2 - a3*b1**2*r3 + a4*b1**2*r4)
f(652) = dtau**2*(-a1**2*b2*r1 - a2**2*b2*r2 + a3**2*b1*r3 + a4**2*b1*r4 &
      )
f(653) = dtau*(-a1**2*b2**2*r1 + a2**2*b2**2*r2 + a3**2*b1**2*r3 - a4**2 &
      *b1**2*r4)
f(654) = dtau**2*(-a1*b2*r1**2 - a2*b2*r2**2 + a3*b1*r3**2 + a4*b1*r4**2 &
      )
f(655) = dtau*(-a1*b2**2*r1**2 + a2*b2**2*r2**2 + a3*b1**2*r3**2 - a4*b1 &
      **2*r4**2)
f(656) = -a1*a4*b1*r2 + a1*a4*b2*r3 - a2*a3*b1*r1 + a2*a3*b2*r4
f(657) = a1*a4*b1**3*r2 - a1*a4*b2**3*r3 + a2*a3*b1**3*r1 - a2*a3*b2**3* &
      r4
f(658) = a1**2*a4*b2*r3 - a1*a4**2*b1*r2 + a2**2*a3*b2*r4 - a2*a3**2*b1* &
      r1
f(659) = -a1**3*a4*b2*r3 + a1*a4**3*b1*r2 - a2**3*a3*b2*r4 + a2*a3**3*b1 &
      *r1
f(660) = -a1**2*a4*b1*r2 + a1*a4**2*b2*r3 - a2**2*a3*b1*r1 + a2*a3**2*b2 &
      *r4
f(661) = a1**2*a4**2*b1*r2 - a1**2*a4**2*b2*r3 + a2**2*a3**2*b1*r1 - a2 &
      **2*a3**2*b2*r4
f(662) = -a1**3*a4*b1*r2 + a1*a4**3*b2*r3 - a2**3*a3*b1*r1 + a2*a3**3*b2 &
      *r4
f(663) = -a1*a4*b1*r2**2 + a1*a4*b2*r3**2 - a2*a3*b1*r1**2 + a2*a3*b2*r4 &
      **2
f(664) = -a1**2*a4*b2*r3**2 + a1*a4**2*b1*r2**2 - a2**2*a3*b2*r4**2 + a2 &
      *a3**2*b1*r1**2
f(665) = -a1**2*a4*b1*r2**2 + a1*a4**2*b2*r3**2 - a2**2*a3*b1*r1**2 + a2 &
      *a3**2*b2*r4**2
f(666) = a1*a4*b1*r2**3 - a1*a4*b2*r3**3 + a2*a3*b1*r1**3 - a2*a3*b2*r4 &
      **3
f(667) = a1*a4*b1*r3 - a1*a4*b2*r2 + a2*a3*b1*r4 - a2*a3*b2*r1
f(668) = -a1*a4*b1**3*r3 + a1*a4*b2**3*r2 - a2*a3*b1**3*r4 + a2*a3*b2**3 &
      *r1
f(669) = a1**2*a4*b1*r3 - a1*a4**2*b2*r2 + a2**2*a3*b1*r4 - a2*a3**2*b2* &
      r1
f(670) = -a1**3*a4*b1*r3 + a1*a4**3*b2*r2 - a2**3*a3*b1*r4 + a2*a3**3*b2 &
      *r1
f(671) = a1**2*a4*b2*r2 - a1*a4**2*b1*r3 + a2**2*a3*b2*r1 - a2*a3**2*b1* &
      r4
f(672) = a1**2*a4**2*b1*r3 - a1**2*a4**2*b2*r2 + a2**2*a3**2*b1*r4 - a2 &
      **2*a3**2*b2*r1
f(673) = -a1**3*a4*b2*r2 + a1*a4**3*b1*r3 - a2**3*a3*b2*r1 + a2*a3**3*b1 &
      *r4
f(674) = -a1*a4*b1*r3**2 + a1*a4*b2*r2**2 - a2*a3*b1*r4**2 + a2*a3*b2*r1 &
      **2
f(675) = -a1**2*a4*b1*r3**2 + a1*a4**2*b2*r2**2 - a2**2*a3*b1*r4**2 + a2 &
      *a3**2*b2*r1**2
f(676) = -a1**2*a4*b2*r2**2 + a1*a4**2*b1*r3**2 - a2**2*a3*b2*r1**2 + a2 &
      *a3**2*b1*r4**2
f(677) = a1*a4*b1*r3**3 - a1*a4*b2*r2**3 + a2*a3*b1*r4**3 - a2*a3*b2*r1 &
      **3
f(678) = dtau*(a1*a4*r2 + a1*a4*r3 - a2*a3*r1 - a2*a3*r4)
f(679) = dtau**3*(-a1*a4*r2 - a1*a4*r3 + a2*a3*r1 + a2*a3*r4)
f(680) = dtau*(-a1**2*a4*r3 - a1*a4**2*r2 + a2**2*a3*r4 + a2*a3**2*r1)
f(681) = dtau*(a1**3*a4*r3 + a1*a4**3*r2 - a2**3*a3*r4 - a2*a3**3*r1)
f(682) = dtau*(a1**2*a4*r2 + a1*a4**2*r3 - a2**2*a3*r1 - a2*a3**2*r4)
f(683) = dtau*(-a1**2*a4**2*r2 - a1**2*a4**2*r3 + a2**2*a3**2*r1 + a2**2 &
      *a3**2*r4)
f(684) = dtau*(a1**3*a4*r2 + a1*a4**3*r3 - a2**3*a3*r1 - a2*a3**3*r4)
f(685) = dtau*(-a1*a4*r2**2 - a1*a4*r3**2 + a2*a3*r1**2 + a2*a3*r4**2)
f(686) = dtau*(-a1**2*a4*r3**2 - a1*a4**2*r2**2 + a2**2*a3*r4**2 + a2*a3 &
      **2*r1**2)
f(687) = dtau*(-a1**2*a4*r2**2 - a1*a4**2*r3**2 + a2**2*a3*r1**2 + a2*a3 &
      **2*r4**2)
f(688) = dtau*(-a1*a4*r2**3 - a1*a4*r3**3 + a2*a3*r1**3 + a2*a3*r4**3)
f(689) = -a1*a3*b1*r2 + a1*a3*b2*r4 - a2*a4*b1*r1 + a2*a4*b2*r3
f(690) = -a1*a3*b1**3*r2 + a1*a3*b2**3*r4 - a2*a4*b1**3*r1 + a2*a4*b2**3 &
      *r3
f(691) = -a1**2*a3*b2*r4 + a1*a3**2*b1*r2 - a2**2*a4*b2*r3 + a2*a4**2*b1 &
      *r1
f(692) = a1**3*a3*b2*r4 - a1*a3**3*b1*r2 + a2**3*a4*b2*r3 - a2*a4**3*b1* &
      r1
f(693) = -a1**2*a3*b1*r2 + a1*a3**2*b2*r4 - a2**2*a4*b1*r1 + a2*a4**2*b2 &
      *r3
f(694) = -a1**2*a3**2*b1*r2 + a1**2*a3**2*b2*r4 - a2**2*a4**2*b1*r1 + a2 &
      **2*a4**2*b2*r3
f(695) = -a1**3*a3*b1*r2 + a1*a3**3*b2*r4 - a2**3*a4*b1*r1 + a2*a4**3*b2 &
      *r3
f(696) = -a1*a3*b1*r2**2 + a1*a3*b2*r4**2 - a2*a4*b1*r1**2 + a2*a4*b2*r3 &
      **2
f(697) = a1**2*a3*b2*r4**2 - a1*a3**2*b1*r2**2 + a2**2*a4*b2*r3**2 - a2* &
      a4**2*b1*r1**2
f(698) = -a1**2*a3*b1*r2**2 + a1*a3**2*b2*r4**2 - a2**2*a4*b1*r1**2 + a2 &
      *a4**2*b2*r3**2
f(699) = a1*a3*b1*r2**3 - a1*a3*b2*r4**3 + a2*a4*b1*r1**3 - a2*a4*b2*r3 &
      **3
f(700) = -a1*a3*b1*r4 + a1*a3*b2*r2 - a2*a4*b1*r3 + a2*a4*b2*r1
f(701) = -a1*a3*b1**3*r4 + a1*a3*b2**3*r2 - a2*a4*b1**3*r3 + a2*a4*b2**3 &
      *r1
f(702) = -a1**2*a3*b1*r4 + a1*a3**2*b2*r2 - a2**2*a4*b1*r3 + a2*a4**2*b2 &
      *r1
f(703) = a1**3*a3*b1*r4 - a1*a3**3*b2*r2 + a2**3*a4*b1*r3 - a2*a4**3*b2* &
      r1
f(704) = a1**2*a3*b2*r2 - a1*a3**2*b1*r4 + a2**2*a4*b2*r1 - a2*a4**2*b1* &
      r3
f(705) = -a1**2*a3**2*b1*r4 + a1**2*a3**2*b2*r2 - a2**2*a4**2*b1*r3 + a2 &
      **2*a4**2*b2*r1
f(706) = -a1**3*a3*b2*r2 + a1*a3**3*b1*r4 - a2**3*a4*b2*r1 + a2*a4**3*b1 &
      *r3
f(707) = -a1*a3*b1*r4**2 + a1*a3*b2*r2**2 - a2*a4*b1*r3**2 + a2*a4*b2*r1 &
      **2
f(708) = -a1**2*a3*b1*r4**2 + a1*a3**2*b2*r2**2 - a2**2*a4*b1*r3**2 + a2 &
      *a4**2*b2*r1**2
f(709) = -a1**2*a3*b2*r2**2 + a1*a3**2*b1*r4**2 - a2**2*a4*b2*r1**2 + a2 &
      *a4**2*b1*r3**2
f(710) = a1*a3*b1*r4**3 - a1*a3*b2*r2**3 + a2*a4*b1*r3**3 - a2*a4*b2*r1 &
      **3
f(711) = dtau*(a1*a3*r2 - a1*a3*r4 - a2*a4*r1 + a2*a4*r3)
f(712) = dtau**3*(-a1*a3*r2 + a1*a3*r4 + a2*a4*r1 - a2*a4*r3)
f(713) = dtau*(-a1**2*a3*r4 + a1*a3**2*r2 + a2**2*a4*r3 - a2*a4**2*r1)
f(714) = dtau*(a1**3*a3*r4 - a1*a3**3*r2 - a2**3*a4*r3 + a2*a4**3*r1)
f(715) = dtau*(-a1**2*a3*r2 + a1*a3**2*r4 + a2**2*a4*r1 - a2*a4**2*r3)
f(716) = dtau*(-a1**2*a3**2*r2 + a1**2*a3**2*r4 + a2**2*a4**2*r1 - a2**2 &
      *a4**2*r3)
f(717) = dtau*(-a1**3*a3*r2 + a1*a3**3*r4 + a2**3*a4*r1 - a2*a4**3*r3)
f(718) = dtau*(-a1*a3*r2**2 + a1*a3*r4**2 + a2*a4*r1**2 - a2*a4*r3**2)
f(719) = dtau*(-a1**2*a3*r4**2 + a1*a3**2*r2**2 + a2**2*a4*r3**2 - a2*a4 &
      **2*r1**2)
f(720) = dtau*(-a1**2*a3*r2**2 + a1*a3**2*r4**2 + a2**2*a4*r1**2 - a2*a4 &
      **2*r3**2)
f(721) = dtau*(-a1*a3*r2**3 + a1*a3*r4**3 + a2*a4*r1**3 - a2*a4*r3**3)
f(722) = b1*b2*(a1*b2*r2 + a2*b2*r1 - a3*b1*r4 - a4*b1*r3)
f(723) = b1*b2*(a1*b1*r2 + a2*b1*r1 - a3*b2*r4 - a4*b2*r3)
f(724) = b1*b2*(a1**2*b2*r2 + a2**2*b2*r1 - a3**2*b1*r4 - a4**2*b1*r3)
f(725) = b1*b2*(-a1**2*b1*r2 - a2**2*b1*r1 + a3**2*b2*r4 + a4**2*b2*r3)
f(726) = b1*b2*(-a1*b2*r2**2 - a2*b2*r1**2 + a3*b1*r4**2 + a4*b1*r3**2)
f(727) = b1*b2*(a1*b1*r2**2 + a2*b1*r1**2 - a3*b2*r4**2 - a4*b2*r3**2)
f(728) = dtau**2*(a1*b1*r2 + a2*b1*r1 - a3*b2*r4 - a4*b2*r3)
f(729) = dtau*(-a1*b1**2*r2 + a2*b1**2*r1 + a3*b2**2*r4 - a4*b2**2*r3)
f(730) = dtau**2*(-a1**2*b1*r2 - a2**2*b1*r1 + a3**2*b2*r4 + a4**2*b2*r3 &
      )
f(731) = dtau*(a1**2*b1**2*r2 - a2**2*b1**2*r1 - a3**2*b2**2*r4 + a4**2* &
      b2**2*r3)
f(732) = dtau**2*(a1*b1*r2**2 + a2*b1*r1**2 - a3*b2*r4**2 - a4*b2*r3**2)
f(733) = dtau*(a1*b1**2*r2**2 - a2*b1**2*r1**2 - a3*b2**2*r4**2 + a4*b2 &
      **2*r3**2)
f(734) = dtau**2*(-a1*b2*r2 - a2*b2*r1 + a3*b1*r4 + a4*b1*r3)
f(735) = dtau*(-a1*b2**2*r2 + a2*b2**2*r1 + a3*b1**2*r4 - a4*b1**2*r3)
f(736) = dtau**2*(a1**2*b2*r2 + a2**2*b2*r1 - a3**2*b1*r4 - a4**2*b1*r3)
f(737) = dtau*(a1**2*b2**2*r2 - a2**2*b2**2*r1 - a3**2*b1**2*r4 + a4**2* &
      b1**2*r3)
f(738) = dtau**2*(a1*b2*r2**2 + a2*b2*r1**2 - a3*b1*r4**2 - a4*b1*r3**2)
f(739) = dtau*(a1*b2**2*r2**2 - a2*b2**2*r1**2 - a3*b1**2*r4**2 + a4*b1 &
      **2*r3**2)
f(740) = -a1*a2*b2*r3 - a1*a2*b2*r4 + a3*a4*b1*r1 + a3*a4*b1*r2
f(741) = a1*a2*b2**3*r3 + a1*a2*b2**3*r4 - a3*a4*b1**3*r1 - a3*a4*b1**3* &
      r2
f(742) = a1**2*a2*b2*r4 + a1*a2**2*b2*r3 - a3**2*a4*b1*r2 - a3*a4**2*b1* &
      r1
f(743) = a1**3*a2*b2*r4 + a1*a2**3*b2*r3 - a3**3*a4*b1*r2 - a3*a4**3*b1* &
      r1
f(744) = -a1**2*a2*b2*r3 - a1*a2**2*b2*r4 + a3**2*a4*b1*r1 + a3*a4**2*b1 &
      *r2
f(745) = -a1**2*a2**2*b2*r3 - a1**2*a2**2*b2*r4 + a3**2*a4**2*b1*r1 + a3 &
      **2*a4**2*b1*r2
f(746) = a1**3*a2*b2*r3 + a1*a2**3*b2*r4 - a3**3*a4*b1*r1 - a3*a4**3*b1* &
      r2
f(747) = -a1*a2*b2*r3**2 - a1*a2*b2*r4**2 + a3*a4*b1*r1**2 + a3*a4*b1*r2 &
      **2
f(748) = -a1**2*a2*b2*r4**2 - a1*a2**2*b2*r3**2 + a3**2*a4*b1*r2**2 + a3 &
      *a4**2*b1*r1**2
f(749) = a1**2*a2*b2*r3**2 + a1*a2**2*b2*r4**2 - a3**2*a4*b1*r1**2 - a3* &
      a4**2*b1*r2**2
f(750) = -a1*a2*b2*r3**3 - a1*a2*b2*r4**3 + a3*a4*b1*r1**3 + a3*a4*b1*r2 &
      **3
f(751) = a1*a2*b1*r3 + a1*a2*b1*r4 - a3*a4*b2*r1 - a3*a4*b2*r2
f(752) = -a1*a2*b1**3*r3 - a1*a2*b1**3*r4 + a3*a4*b2**3*r1 + a3*a4*b2**3 &
      *r2
f(753) = a1**2*a2*b1*r4 + a1*a2**2*b1*r3 - a3**2*a4*b2*r2 - a3*a4**2*b2* &
      r1
f(754) = a1**3*a2*b1*r4 + a1*a2**3*b1*r3 - a3**3*a4*b2*r2 - a3*a4**3*b2* &
      r1
f(755) = -a1**2*a2*b1*r3 - a1*a2**2*b1*r4 + a3**2*a4*b2*r1 + a3*a4**2*b2 &
      *r2
f(756) = -a1**2*a2**2*b1*r3 - a1**2*a2**2*b1*r4 + a3**2*a4**2*b2*r1 + a3 &
      **2*a4**2*b2*r2
f(757) = a1**3*a2*b1*r3 + a1*a2**3*b1*r4 - a3**3*a4*b2*r1 - a3*a4**3*b2* &
      r2
f(758) = -a1*a2*b1*r3**2 - a1*a2*b1*r4**2 + a3*a4*b2*r1**2 + a3*a4*b2*r2 &
      **2
f(759) = a1**2*a2*b1*r4**2 + a1*a2**2*b1*r3**2 - a3**2*a4*b2*r2**2 - a3* &
      a4**2*b2*r1**2
f(760) = a1**2*a2*b1*r3**2 + a1*a2**2*b1*r4**2 - a3**2*a4*b2*r1**2 - a3* &
      a4**2*b2*r2**2
f(761) = a1*a2*b1*r3**3 + a1*a2*b1*r4**3 - a3*a4*b2*r1**3 - a3*a4*b2*r2 &
      **3
f(762) = dtau*(a1*a2*r3 - a1*a2*r4 - a3*a4*r1 + a3*a4*r2)
f(763) = dtau**3*(a1*a2*r3 - a1*a2*r4 - a3*a4*r1 + a3*a4*r2)
f(764) = dtau*(-a1**2*a2*r4 + a1*a2**2*r3 + a3**2*a4*r2 - a3*a4**2*r1)
f(765) = dtau*(-a1**3*a2*r4 + a1*a2**3*r3 + a3**3*a4*r2 - a3*a4**3*r1)
f(766) = dtau*(a1**2*a2*r3 - a1*a2**2*r4 - a3**2*a4*r1 + a3*a4**2*r2)
f(767) = dtau*(a1**2*a2**2*r3 - a1**2*a2**2*r4 - a3**2*a4**2*r1 + a3**2* &
      a4**2*r2)
f(768) = dtau*(a1**3*a2*r3 - a1*a2**3*r4 - a3**3*a4*r1 + a3*a4**3*r2)
f(769) = dtau*(a1*a2*r3**2 - a1*a2*r4**2 - a3*a4*r1**2 + a3*a4*r2**2)
f(770) = dtau*(-a1**2*a2*r4**2 + a1*a2**2*r3**2 + a3**2*a4*r2**2 - a3*a4 &
      **2*r1**2)
f(771) = dtau*(a1**2*a2*r3**2 - a1*a2**2*r4**2 - a3**2*a4*r1**2 + a3*a4 &
      **2*r2**2)
f(772) = dtau*(a1*a2*r3**3 - a1*a2*r4**3 - a3*a4*r1**3 + a3*a4*r2**3)
f(773) = b1*b2*(a1*b1*r3 + a2*b1*r4 - a3*b2*r1 - a4*b2*r2)
f(774) = b1*b2*(a1*b2*r3 + a2*b2*r4 - a3*b1*r1 - a4*b1*r2)
f(775) = b1*b2*(-a1**2*b1*r3 - a2**2*b1*r4 + a3**2*b2*r1 + a4**2*b2*r2)
f(776) = b1*b2*(a1**2*b2*r3 + a2**2*b2*r4 - a3**2*b1*r1 - a4**2*b1*r2)
f(777) = b1*b2*(a1*b1*r3**2 + a2*b1*r4**2 - a3*b2*r1**2 - a4*b2*r2**2)
f(778) = b1*b2*(a1*b2*r3**2 + a2*b2*r4**2 - a3*b1*r1**2 - a4*b1*r2**2)
f(779) = dtau**2*(a1*b2*r3 + a2*b2*r4 - a3*b1*r1 - a4*b1*r2)
f(780) = dtau*(-a1*b2**2*r3 + a2*b2**2*r4 + a3*b1**2*r1 - a4*b1**2*r2)
f(781) = dtau**2*(-a1**2*b2*r3 - a2**2*b2*r4 + a3**2*b1*r1 + a4**2*b1*r2 &
      )
f(782) = dtau*(-a1**2*b2**2*r3 + a2**2*b2**2*r4 + a3**2*b1**2*r1 - a4**2 &
      *b1**2*r2)
f(783) = dtau**2*(-a1*b2*r3**2 - a2*b2*r4**2 + a3*b1*r1**2 + a4*b1*r2**2 &
      )
f(784) = dtau*(a1*b2**2*r3**2 - a2*b2**2*r4**2 - a3*b1**2*r1**2 + a4*b1 &
      **2*r2**2)
f(785) = dtau**2*(-a1*b1*r3 - a2*b1*r4 + a3*b2*r1 + a4*b2*r2)
f(786) = dtau*(-a1*b1**2*r3 + a2*b1**2*r4 + a3*b2**2*r1 - a4*b2**2*r2)
f(787) = dtau**2*(a1**2*b1*r3 + a2**2*b1*r4 - a3**2*b2*r1 - a4**2*b2*r2)
f(788) = dtau*(-a1**2*b1**2*r3 + a2**2*b1**2*r4 + a3**2*b2**2*r1 - a4**2 &
      *b2**2*r2)
f(789) = dtau**2*(-a1*b1*r3**2 - a2*b1*r4**2 + a3*b2*r1**2 + a4*b2*r2**2 &
      )
f(790) = dtau*(-a1*b1**2*r3**2 + a2*b1**2*r4**2 + a3*b2**2*r1**2 - a4*b2 &
      **2*r2**2)
f(791) = b1*b2*(-a1*b1*r4 - a2*b1*r3 + a3*b2*r2 + a4*b2*r1)
f(792) = b1*b2*(-a1*b2*r4 - a2*b2*r3 + a3*b1*r2 + a4*b1*r1)
f(793) = b1*b2*(a1**2*b1*r4 + a2**2*b1*r3 - a3**2*b2*r2 - a4**2*b2*r1)
f(794) = b1*b2*(a1**2*b2*r4 + a2**2*b2*r3 - a3**2*b1*r2 - a4**2*b1*r1)
f(795) = b1*b2*(-a1*b1*r4**2 - a2*b1*r3**2 + a3*b2*r2**2 + a4*b2*r1**2)
f(796) = b1*b2*(a1*b2*r4**2 + a2*b2*r3**2 - a3*b1*r2**2 - a4*b1*r1**2)
f(797) = dtau**2*(a1*b2*r4 + a2*b2*r3 - a3*b1*r2 - a4*b1*r1)
f(798) = dtau*(a1*b2**2*r4 - a2*b2**2*r3 - a3*b1**2*r2 + a4*b1**2*r1)
f(799) = dtau**2*(a1**2*b2*r4 + a2**2*b2*r3 - a3**2*b1*r2 - a4**2*b1*r1)
f(800) = dtau*(a1**2*b2**2*r4 - a2**2*b2**2*r3 - a3**2*b1**2*r2 + a4**2* &
      b1**2*r1)
f(801) = dtau**2*(a1*b2*r4**2 + a2*b2*r3**2 - a3*b1*r2**2 - a4*b1*r1**2)
f(802) = dtau*(-a1*b2**2*r4**2 + a2*b2**2*r3**2 + a3*b1**2*r2**2 - a4*b1 &
      **2*r1**2)
f(803) = dtau**2*(-a1*b1*r4 - a2*b1*r3 + a3*b2*r2 + a4*b2*r1)
f(804) = dtau*(a1*b1**2*r4 - a2*b1**2*r3 - a3*b2**2*r2 + a4*b2**2*r1)
f(805) = dtau**2*(a1**2*b1*r4 + a2**2*b1*r3 - a3**2*b2*r2 - a4**2*b2*r1)
f(806) = dtau*(a1**2*b1**2*r4 - a2**2*b1**2*r3 - a3**2*b2**2*r2 + a4**2* &
      b2**2*r1)
f(807) = dtau**2*(-a1*b1*r4**2 - a2*b1*r3**2 + a3*b2*r2**2 + a4*b2*r1**2 &
      )
f(808) = dtau*(a1*b1**2*r4**2 - a2*b1**2*r3**2 - a3*b2**2*r2**2 + a4*b2 &
      **2*r1**2)
f(809) = b1*b2*dtau*(-r1 + r2 + r3 - r4)
f(810) = b1*b2*dtau**3*(r1 - r2 - r3 + r4)
f(811) = b1*b2*dtau**2*(-b1*r3 - b1*r4 + b2*r1 + b2*r2)
f(812) = b1*b2*dtau*(-b1**2*r3 + b1**2*r4 + b2**2*r1 - b2**2*r2)
f(813) = b1*b2*dtau**2*(b1*r1 + b1*r2 - b2*r3 - b2*r4)
f(814) = b1**2*b2**2*dtau*(-r1 + r2 + r3 - r4)
f(815) = b1*b2*dtau*(b1**2*r1 - b1**2*r2 - b2**2*r3 + b2**2*r4)
f(816) = b1*b2*dtau*(r1**2 - r2**2 - r3**2 + r4**2)
f(817) = b1*b2*dtau*(r1**3 - r2**3 - r3**3 + r4**3)
f(818) = a1*a2*a3*b1 + a1*a2*a4*b1 - a1*a3*a4*b2 - a2*a3*a4*b2
f(819) = a1*a2*a3*b1**3 + a1*a2*a4*b1**3 - a1*a3*a4*b2**3 - a2*a3*a4*b2 &
      **3
f(820) = -a1**2*a3*a4*b2 + a1*a2*a3**2*b1 + a1*a2*a4**2*b1 - a2**2*a3*a4 &
      *b2
f(821) = -a1**3*a3*a4*b2 + a1*a2*a3**3*b1 + a1*a2*a4**3*b1 - a2**3*a3*a4 &
      *b2
f(822) = a1**2*a2*a4*b1 + a1*a2**2*a3*b1 - a1*a3*a4**2*b2 - a2*a3**2*a4* &
      b2
f(823) = a1**2*a2*a4**2*b1 - a1**2*a3*a4**2*b2 + a1*a2**2*a3**2*b1 - a2 &
      **2*a3**2*a4*b2
f(824) = a1**3*a2*a4*b1 + a1*a2**3*a3*b1 - a1*a3*a4**3*b2 - a2*a3**3*a4* &
      b2
f(825) = a1**2*a2*a3*b1 + a1*a2**2*a4*b1 - a1*a3**2*a4*b2 - a2*a3*a4**2* &
      b2
f(826) = a1**2*a2*a3**2*b1 - a1**2*a3**2*a4*b2 + a1*a2**2*a4**2*b1 - a2 &
      **2*a3*a4**2*b2
f(827) = a1**2*a2**2*a3*b1 + a1**2*a2**2*a4*b1 - a1*a3**2*a4**2*b2 - a2* &
      a3**2*a4**2*b2
f(828) = a1**3*a2*a3*b1 + a1*a2**3*a4*b1 - a1*a3**3*a4*b2 - a2*a3*a4**3* &
      b2
f(829) = -a1*a2*a3*b2 - a1*a2*a4*b2 + a1*a3*a4*b1 + a2*a3*a4*b1
f(830) = -a1*a2*a3*b2**3 - a1*a2*a4*b2**3 + a1*a3*a4*b1**3 + a2*a3*a4*b1 &
      **3
f(831) = a1**2*a3*a4*b1 - a1*a2*a3**2*b2 - a1*a2*a4**2*b2 + a2**2*a3*a4* &
      b1
f(832) = a1**3*a3*a4*b1 - a1*a2*a3**3*b2 - a1*a2*a4**3*b2 + a2**3*a3*a4* &
      b1
f(833) = -a1**2*a2*a4*b2 - a1*a2**2*a3*b2 + a1*a3*a4**2*b1 + a2*a3**2*a4 &
      *b1
f(834) = a1**2*a2*a4**2*b2 - a1**2*a3*a4**2*b1 + a1*a2**2*a3**2*b2 - a2 &
      **2*a3**2*a4*b1
f(835) = a1**3*a2*a4*b2 + a1*a2**3*a3*b2 - a1*a3*a4**3*b1 - a2*a3**3*a4* &
      b1
f(836) = a1**2*a2*a3*b2 + a1*a2**2*a4*b2 - a1*a3**2*a4*b1 - a2*a3*a4**2* &
      b1
f(837) = a1**2*a2*a3**2*b2 - a1**2*a3**2*a4*b1 + a1*a2**2*a4**2*b2 - a2 &
      **2*a3*a4**2*b1
f(838) = -a1**2*a2**2*a3*b2 - a1**2*a2**2*a4*b2 + a1*a3**2*a4**2*b1 + a2 &
      *a3**2*a4**2*b1
f(839) = -a1**3*a2*a3*b2 - a1*a2**3*a4*b2 + a1*a3**3*a4*b1 + a2*a3*a4**3 &
      *b1
f(840) = dtau*(a1*a2*a3 - a1*a2*a4 - a1*a3*a4 + a2*a3*a4)
f(841) = dtau**3*(a1*a2*a3 - a1*a2*a4 - a1*a3*a4 + a2*a3*a4)
f(842) = dtau*(a1**2*a3*a4 - a1*a2*a3**2 + a1*a2*a4**2 - a2**2*a3*a4)
f(843) = dtau*(-a1**3*a3*a4 + a1*a2*a3**3 - a1*a2*a4**3 + a2**3*a3*a4)
f(844) = dtau*(-a1**2*a2*a4 + a1*a2**2*a3 - a1*a3*a4**2 + a2*a3**2*a4)
f(845) = dtau*(-a1**2*a2*a4**2 - a1**2*a3*a4**2 + a1*a2**2*a3**2 + a2**2 &
      *a3**2*a4)
f(846) = dtau*(a1**3*a2*a4 - a1*a2**3*a3 + a1*a3*a4**3 - a2*a3**3*a4)
f(847) = dtau*(-a1**2*a2*a3 + a1*a2**2*a4 + a1*a3**2*a4 - a2*a3*a4**2)
f(848) = dtau*(a1**2*a2*a3**2 - a1**2*a3**2*a4 - a1*a2**2*a4**2 + a2**2* &
      a3*a4**2)
f(849) = dtau*(-a1**2*a2**2*a3 + a1**2*a2**2*a4 + a1*a3**2*a4**2 - a2*a3 &
      **2*a4**2)
f(850) = dtau*(-a1**3*a2*a3 + a1*a2**3*a4 + a1*a3**3*a4 - a2*a3*a4**3)
f(851) = b1*b2*(-a1*a2*b2 + a3*a4*b1)
f(852) = b1*b2*(-a1*a2*b1 + a3*a4*b2)
f(853) = b1*b2*(a1**2*a2*b2 + a1*a2**2*b2 - a3**2*a4*b1 - a3*a4**2*b1)
f(854) = b1*b2*(a1**2*a2*b1 + a1*a2**2*b1 - a3**2*a4*b2 - a3*a4**2*b2)
f(855) = dtau**2*(a1*a2*b1 - a3*a4*b2)
f(856) = dtau**2*(-a1**2*a2*b1 - a1*a2**2*b1 + a3**2*a4*b2 + a3*a4**2*b2 &
      )
f(857) = dtau*(-a1**2*a2*b1**2 + a1*a2**2*b1**2 + a3**2*a4*b2**2 - a3*a4 &
      **2*b2**2)
f(858) = dtau**2*(-a1*a2*b2 + a3*a4*b1)
f(859) = dtau**2*(a1**2*a2*b2 + a1*a2**2*b2 - a3**2*a4*b1 - a3*a4**2*b1)
f(860) = dtau*(-a1**2*a2*b2**2 + a1*a2**2*b2**2 + a3**2*a4*b1**2 - a3*a4 &
      **2*b1**2)
f(861) = b1*b2*(a1*a3*b1 - a1*a3*b2 + a2*a4*b1 - a2*a4*b2)
f(862) = b1*b2*(-a1**2*a3*b1 + a1*a3**2*b2 - a2**2*a4*b1 + a2*a4**2*b2)
f(863) = b1*b2*(-a1**2*a3*b2 + a1*a3**2*b1 - a2**2*a4*b2 + a2*a4**2*b1)
f(864) = dtau**2*(a1*a3*b1 - a1*a3*b2 + a2*a4*b1 - a2*a4*b2)
f(865) = dtau*(a1*a3*b1**2 - a1*a3*b2**2 - a2*a4*b1**2 + a2*a4*b2**2)
f(866) = dtau**2*(-a1**2*a3*b2 + a1*a3**2*b1 - a2**2*a4*b2 + a2*a4**2*b1 &
      )
f(867) = dtau*(-a1**2*a3*b2**2 + a1*a3**2*b1**2 + a2**2*a4*b2**2 - a2*a4 &
      **2*b1**2)
f(868) = dtau**2*(a1**2*a3*b1 - a1*a3**2*b2 + a2**2*a4*b1 - a2*a4**2*b2)
f(869) = dtau*(-a1**2*a3*b1**2 + a1*a3**2*b2**2 + a2**2*a4*b1**2 - a2*a4 &
      **2*b2**2)
f(870) = b1*b2*(-a1*a4*b1 + a1*a4*b2 - a2*a3*b1 + a2*a3*b2)
f(871) = b1*b2*(-a1**2*a4*b1 + a1*a4**2*b2 - a2**2*a3*b1 + a2*a3**2*b2)
f(872) = b1*b2*(-a1**2*a4*b2 + a1*a4**2*b1 - a2**2*a3*b2 + a2*a3**2*b1)
f(873) = dtau**2*(a1*a4*b1 - a1*a4*b2 + a2*a3*b1 - a2*a3*b2)
f(874) = dtau*(-a1*a4*b1**2 - a1*a4*b2**2 + a2*a3*b1**2 + a2*a3*b2**2)
f(875) = dtau**2*(-a1**2*a4*b2 + a1*a4**2*b1 - a2**2*a3*b2 + a2*a3**2*b1 &
      )
f(876) = dtau*(-a1**2*a4*b2**2 - a1*a4**2*b1**2 + a2**2*a3*b2**2 + a2*a3 &
      **2*b1**2)
f(877) = dtau**2*(a1**2*a4*b1 - a1*a4**2*b2 + a2**2*a3*b1 - a2*a3**2*b2)
f(878) = dtau*(-a1**2*a4*b1**2 - a1*a4**2*b2**2 + a2**2*a3*b1**2 + a2*a3 &
      **2*b2**2)
f(879) = b1*b2*dtau*(a1 - a2 - a3 + a4)
f(880) = b1*b2*dtau**3*(-a1 + a2 + a3 - a4)
f(881) = b1*b2*dtau**2*(-a1*b2 - a2*b2 + a3*b1 + a4*b1)
f(882) = b1*b2*dtau*(-a1*b2**2 + a2*b2**2 + a3*b1**2 - a4*b1**2)
f(883) = b1*b2*dtau**2*(-a1*b1 - a2*b1 + a3*b2 + a4*b2)
f(884) = b1**2*b2**2*dtau*(-a1 + a2 + a3 - a4)
f(885) = b1*b2*dtau*(-a1*b1**2 + a2*b1**2 + a3*b2**2 - a4*b2**2)
f(886) = b1*b2*dtau*(-a1**2 + a2**2 + a3**2 - a4**2)
f(887) = b1*b2*dtau*(a1**3 - a2**3 - a3**3 + a4**3)
v = sum(f*params)
end function c2h4_dipole_b3u_n4_d6_ADF
