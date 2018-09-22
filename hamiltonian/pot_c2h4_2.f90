subroutine poten_c2h4_886666_ADF_noreexp(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: i, irank, j, nc
  real(ark) :: rad_, re0, re1, alphae, betae, taue, amorse1, amorse2, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7), &!
               gamma1, gamma2, delta1, delta2, delta3, delta4, delta5, coefs_sym(maxval(func%nparams))
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

  coefs_sym = 0
  j = 0
  do i=26, func%nparams(irank)
    j = j + 1
    coefs_sym(j) = func%params(i,irank)
  enddo
  nc = j


  vshort = func%params(25,irank) &!
         + pot_c2h4_init_poten(molec%nmodes, (/y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),beta1,beta2,y(12)/), coefs_sym(1:nc))

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

end subroutine poten_c2h4_886666_ADF_noreexp


!################################################################################


function pot_c2h4_init_poten(nmodes, local, f) result(res)
  use adf
  implicit none

  integer(ik), intent(in) :: nmodes
  type(adf_realq), intent(in) :: local(nmodes)
  real(ark), intent(in) :: f(:)
  type(adf_realq) :: res

  integer(ik), parameter :: maxnterms = 12000
  integer(ik) :: iterm, nterms, terms(12,maxnterms), imode
  real(ark) :: coefs(maxnterms)
  type(adf_realq) :: prod

  if (nmodes/=12) then
    write(out, '(/a,1x,i3,1x,a)') 'pot_c2h4_init_poten error: number of internal coordinates =', nmodes, '!= 12'
    stop
  endif

  ! compute all symmetry-equivalent expansion coefficients

  nterms = 0
  terms = 0
  coefs = 0.0

#include 'pot_c2h4_2_part2.f90'

  if (nterms==0) then
    write(out, '(/a)') 'pot_c2h4_init_poten error: part 2 of the file pot_c2h4_2.f90 was not included at compilation, recompile module "hamiltonian" using -fpp or -cpp flag'
    stop
  endif

  if (nterms>maxnterms) then
    write(out, '(/a,1x,i6,1x,a,1x,i6)') 'pot_c2h4_init_poten error: number of terms =', nterms, 'exceeds max number =', maxnterms
    stop
  endif

  ! compute potential

  res = 0.0_ark

  do iterm=1, nterms
    prod = 1.0_ark
    do imode=1, nmodes
      if (terms(imode,iterm)==1) then
        prod = prod * local(imode)
      elseif (terms(imode,iterm)>1) then
        prod = prod * local(imode)**terms(imode,iterm)
      endif
    enddo
    res = res + coefs(iterm) * prod
  enddo

end function pot_c2h4_init_poten
