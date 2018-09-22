! This unit defines all specific routines for a six-atomic molecule of ethylene-type

module pot_c2h4
use accuracy
use moltype

implicit none

public MLpoten_c2h4_88, MLpoten_c2h4_886, MLpoten_c2h4_886666, MLpoten_c2h4_886666_nolong, MLpoten_c2h4_lee, MLalpha0_c2h4

private

integer(ik), parameter :: verbose = 4 ! Verbosity level


contains


!---------------------------------------------------------------------


recursive subroutine MLalpha0_c2h4(rank,ncoords,natoms,r,xyz,f)

  integer(ik),intent(in) ::  rank,ncoords,natoms
  real(ark),intent(in)   ::  r(ncoords),xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iterm, imu

  iterm=1
  do imu=1, 6
    f(imu) = extF%coef(iterm, imu)
  enddo

end subroutine MLalpha0_c2h4


!---------------------------------------------------------------------



function MLpoten_c2h4_88(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: rad, y(12), a1, a2, a3, r1, r2, b1_zmat, b2_zmat, tau1_zmat, tau1, tau2, dtau, db1, db2
  real(ark) :: vshort, beta1, beta2, beta12, vlong, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7)
  real(ark) :: vdamp, gamma1, gamma2, delta1, delta2, delta3, delta4, delta5

  rad = pi/180.0_ark

  ! expansion functions

  r1 = force(1)
  r2 = force(2)
  a1 = force(3)*rad
  a2 = force(4)*rad
  a3 = force(5)*rad
  beta1 = force(6)
  beta2 = force(7)

  ! long-range function parameters

  gamma1 = force(8)  ! cc stretch
  gamma2 = force(9)  ! ch stretch

  de_str2(1) = force(10)  ! cc stretch
  de_str4(1) = force(11)  ! cc stretch

  de_str2(2:5) = force(12)  ! ch stretch
  de_str4(2:5) = force(13)  ! ch stretch

  de_bnd2(1:4) = force(14)  ! hcc bend
  de_bnd4(1:4) = force(15)  ! hcc bend

  de_bnd2(5:6) = force(16)  ! beta bend
  de_bnd4(5:6) = force(17)  ! beta bend

  de_bnd2(7) = force(18)  ! tau bend
  de_bnd4(7) = force(19)  ! tau bend

  ! damping-function parameters

  delta1 = force(20)
  delta2 = force(21)
  delta3 = force(22)
  delta4 = force(23)
  delta5 = force(24)

  ! expansion functions

  y(1) = 1.0_ark-exp(-beta1*(local(1)-r1))

  y(2) = 1.0_ark-exp(-beta2*(local(2)-r2))
  y(3) = 1.0_ark-exp(-beta2*(local(3)-r2))
  y(4) = 1.0_ark-exp(-beta2*(local(4)-r2))
  y(5) = 1.0_ark-exp(-beta2*(local(5)-r2))

  y(6) = local(6) - a1
  y(7) = local(7) - a1
  y(8) = local(8) - a1
  y(9) = local(9) - a1

  b1_zmat   = local(10)
  tau1_zmat = local(11)
  b2_zmat   = local(12)

  db1 = b1_zmat - pi
  db2 = b2_zmat - pi

  if (tau1_zmat>pi) then
    tau1 = 2.0_ark*pi - tau1_zmat
  elseif (tau1_zmat<pi) then
    tau1 = -tau1_zmat
  endif

  tau2 = db1 + db2 + tau1
  dtau = tau1 + tau2

  y(10:12) = (/db1, db2, dtau/)

  ! short-range potential

  vshort = force(25) &!
         + c2h4_poten_n1_d8( y, force(26:57) ) &!
         + c2h4_poten_n2_d8( y, force(58:433) )

  ! long-range potential

  vlong = sum(de_str2(1:5)*y(1:5)**2 + de_str4(1:5)*y(1:5)**4) &!
        + exp( -gamma1*(local(1)-r1)**2 -gamma2*sum((local(2:5)-r2)**2) ) &!
                * sum( de_bnd2(1:7)*y(6:12)**2 + de_bnd4(1:7)*y(6:12)**4 )

  ! damping function

  vdamp = exp( -delta1*(local(1)-r1)**2        -2.0_ark*delta1*(local(1)-r1)**4 ) &!
        * exp( -delta2*sum((local(2:5)-r2)**2) -2.0_ark*delta2*sum((local(2:5)-r2)**4) ) &!
        * exp( -delta3*sum(y(6:9)**2)           -2.0_ark*delta3*sum(y(6:9)**4) ) &!
        * exp( -delta4*sum(y(10:11)**2)         -2.0_ark*delta4*sum(y(10:11)**4) ) &!
        * exp( -delta5*sum(y(12:12)**2)         -2.0_ark*delta5*sum(y(12:12)**4) )

  f = vshort*vdamp + vlong

end function MLpoten_c2h4_88


!---------------------------------------------------------------------


function MLpoten_c2h4_886(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: rad, y(12), a1, a2, a3, r1, r2, b1_zmat, b2_zmat, tau1_zmat, tau1, tau2, dtau, db1, db2
  real(ark) :: vshort, beta1, beta2, beta12, vlong, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7)
  real(ark) :: vdamp, gamma1, gamma2, delta1, delta2, delta3, delta4, delta5

  rad = pi/180.0_ark

  ! expansion functions

  r1 = force(1)
  r2 = force(2)
  a1 = force(3)*rad
  a2 = force(4)*rad
  a3 = force(5)*rad
  beta1 = force(6)
  beta2 = force(7)

  ! long-range function parameters

  gamma1 = force(8)  ! cc stretch
  gamma2 = force(9)  ! ch stretch

  de_str2(1) = force(10)  ! cc stretch
  de_str4(1) = force(11)  ! cc stretch

  de_str2(2:5) = force(12)  ! ch stretch
  de_str4(2:5) = force(13)  ! ch stretch

  de_bnd2(1:4) = force(14)  ! hcc bend
  de_bnd4(1:4) = force(15)  ! hcc bend

  de_bnd2(5:6) = force(16)  ! beta bend
  de_bnd4(5:6) = force(17)  ! beta bend

  de_bnd2(7) = force(18)  ! tau bend
  de_bnd4(7) = force(19)  ! tau bend

  ! damping-function parameters

  delta1 = force(20)
  delta2 = force(21)
  delta3 = force(22)
  delta4 = force(23)
  delta5 = force(24)

  ! expansion functions

  y(1) = 1.0_ark-exp(-beta1*(local(1)-r1))

  y(2) = 1.0_ark-exp(-beta2*(local(2)-r2))
  y(3) = 1.0_ark-exp(-beta2*(local(3)-r2))
  y(4) = 1.0_ark-exp(-beta2*(local(4)-r2))
  y(5) = 1.0_ark-exp(-beta2*(local(5)-r2))

  y(6) = local(6) - a1
  y(7) = local(7) - a1
  y(8) = local(8) - a1
  y(9) = local(9) - a1

  b1_zmat   = local(10)
  tau1_zmat = local(11)
  b2_zmat   = local(12)

  db1 = b1_zmat - pi
  db2 = b2_zmat - pi

  if (tau1_zmat>pi) then
    tau1 = 2.0_ark*pi - tau1_zmat
  elseif (tau1_zmat<pi) then
    tau1 = -tau1_zmat
  endif

  tau2 = db1 + db2 + tau1
  dtau = tau1 + tau2

  y(10:12) = (/db1, db2, dtau/)

  ! short-range potential

  vshort = force(25) &
         + c2h4_poten_n1_d8( y, force(26:57) ) &
         + c2h4_poten_n2_d8( y, force(58:433) ) &
         + c2h4_poten_n3_d6( y, force(434:1149) )

  ! long-range potential

  vlong = sum(de_str2(1:5)*y(1:5)**2 + de_str4(1:5)*y(1:5)**4) &!
        + exp( -gamma1*(local(1)-r1)**2 -gamma2*sum((local(2:5)-r2)**2) ) &!
                * sum( de_bnd2(1:7)*y(6:12)**2 + de_bnd4(1:7)*y(6:12)**4 )

  ! damping function

  vdamp = exp( -delta1*(local(1)-r1)**2        -2.0_ark*delta1*(local(1)-r1)**4 ) &!
        * exp( -delta2*sum((local(2:5)-r2)**2) -2.0_ark*delta2*sum((local(2:5)-r2)**4) ) &!
        * exp( -delta3*sum(y(6:9)**2)           -2.0_ark*delta3*sum(y(6:9)**4) ) &!
        * exp( -delta4*sum(y(10:11)**2)         -2.0_ark*delta4*sum(y(10:11)**4) ) &!
        * exp( -delta5*sum(y(12:12)**2)         -2.0_ark*delta5*sum(y(12:12)**4) )

  f = vshort*vdamp + vlong

end function MLpoten_c2h4_886



!---------------------------------------------------------------------


function MLpoten_c2h4_886666(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: rad, y(12), a1, a2, a3, r1, r2, b1_zmat, b2_zmat, tau1_zmat, tau1, tau2, dtau, db1, db2
  real(ark) :: vshort, beta1, beta2, beta12, vlong, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7)
  real(ark) :: vdamp, gamma1, gamma2, delta1, delta2, delta3, delta4, delta5

  rad = pi/180.0_ark

  ! expansion functions

  r1 = force(1)
  r2 = force(2)
  a1 = force(3)*rad
  a2 = force(4)*rad
  a3 = force(5)*rad
  beta1 = force(6)
  beta2 = force(7)

  ! long-range function parameters

  gamma1 = force(8)  ! cc stretch
  gamma2 = force(9)  ! ch stretch

  de_str2(1) = force(10)  ! cc stretch
  de_str4(1) = force(11)  ! cc stretch

  de_str2(2:5) = force(12)  ! ch stretch
  de_str4(2:5) = force(13)  ! ch stretch

  de_bnd2(1:4) = force(14)  ! hcc bend
  de_bnd4(1:4) = force(15)  ! hcc bend

  de_bnd2(5:6) = force(16)  ! beta bend
  de_bnd4(5:6) = force(17)  ! beta bend

  de_bnd2(7) = force(18)  ! tau bend
  de_bnd4(7) = force(19)  ! tau bend

  ! damping-function parameters

  delta1 = force(20)
  delta2 = force(21)
  delta3 = force(22)
  delta4 = force(23)
  delta5 = force(24)

  ! expansion functions

  y(1) = 1.0_ark-exp(-beta1*(local(1)-r1))

  y(2) = 1.0_ark-exp(-beta2*(local(2)-r2))
  y(3) = 1.0_ark-exp(-beta2*(local(3)-r2))
  y(4) = 1.0_ark-exp(-beta2*(local(4)-r2))
  y(5) = 1.0_ark-exp(-beta2*(local(5)-r2))

  y(6) = local(6) - a1
  y(7) = local(7) - a1
  y(8) = local(8) - a1
  y(9) = local(9) - a1

  b1_zmat   = local(10)
  tau1_zmat = local(11)
  b2_zmat   = local(12)

  db1 = b1_zmat - pi
  db2 = b2_zmat - pi

  if (tau1_zmat>pi) then
    tau1 = 2.0_ark*pi - tau1_zmat
  elseif (tau1_zmat<pi) then
    tau1 = -tau1_zmat
  endif

  tau2 = db1 + db2 + tau1
  dtau = tau1 + tau2

  y(10:12) = (/db1, db2, dtau/)

  ! short-range potential

  vshort = force(25) &
         + c2h4_poten_n1_d8( y, force(26:57) ) &
         + c2h4_poten_n2_d8( y, force(58:433) ) &
         + c2h4_poten_n3_d6( y, force(434:1149) ) &
         + c2h4_poten_n4_d6( y, force(1150:2167) ) &
         + c2h4_poten_n5_d6( y, force(2168:2755) ) &
         + c2h4_poten_n6_d6( y, force(2756:2875) )

  ! long-range potential

  vlong = sum(de_str2(1:5)*y(1:5)**2 + de_str4(1:5)*y(1:5)**4) &!
        + exp( -gamma1*(local(1)-r1)**2 -gamma2*sum((local(2:5)-r2)**2) ) &!
                * sum( de_bnd2(1:7)*y(6:12)**2 + de_bnd4(1:7)*y(6:12)**4 )

  ! damping function

  vdamp = exp( -delta1*(local(1)-r1)**2        -2.0_ark*delta1*(local(1)-r1)**4 ) &!
        * exp( -delta2*sum((local(2:5)-r2)**2) -2.0_ark*delta2*sum((local(2:5)-r2)**4) ) &!
        * exp( -delta3*sum(y(6:9)**2)           -2.0_ark*delta3*sum(y(6:9)**4) ) &!
        * exp( -delta4*sum(y(10:11)**2)         -2.0_ark*delta4*sum(y(10:11)**4) ) &!
        * exp( -delta5*sum(y(12:12)**2)         -2.0_ark*delta5*sum(y(12:12)**4) )

  f = vshort*vdamp + vlong

end function MLpoten_c2h4_886666


!---------------------------------------------------------------------


function MLpoten_c2h4_886666_nolong(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: rad, y(12), a1, a2, a3, r1, r2, b1_zmat, b2_zmat, tau1_zmat, tau1, tau2, dtau, db1, db2
  real(ark) :: vshort, beta1, beta2, beta12, vlong, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7)
  real(ark) :: vdamp, gamma1, gamma2, delta1, delta2, delta3, delta4, delta5

  rad = pi/180.0_ark

  ! expansion functions

  r1 = force(1)
  r2 = force(2)
  a1 = force(3)*rad
  a2 = force(4)*rad
  a3 = force(5)*rad
  beta1 = force(6)
  beta2 = force(7)

  ! long-range function parameters

  gamma1 = force(8)  ! cc stretch
  gamma2 = force(9)  ! ch stretch

  de_str2(1) = force(10)  ! cc stretch
  de_str4(1) = force(11)  ! cc stretch

  de_str2(2:5) = force(12)  ! ch stretch
  de_str4(2:5) = force(13)  ! ch stretch

  de_bnd2(1:4) = force(14)  ! hcc bend
  de_bnd4(1:4) = force(15)  ! hcc bend

  de_bnd2(5:6) = force(16)  ! beta bend
  de_bnd4(5:6) = force(17)  ! beta bend

  de_bnd2(7) = force(18)  ! tau bend
  de_bnd4(7) = force(19)  ! tau bend

  ! damping-function parameters

  delta1 = force(20)
  delta2 = force(21)
  delta3 = force(22)
  delta4 = force(23)
  delta5 = force(24)

  ! expansion functions

  y(1) = 1.0_ark-exp(-beta1*(local(1)-r1))

  y(2) = 1.0_ark-exp(-beta2*(local(2)-r2))
  y(3) = 1.0_ark-exp(-beta2*(local(3)-r2))
  y(4) = 1.0_ark-exp(-beta2*(local(4)-r2))
  y(5) = 1.0_ark-exp(-beta2*(local(5)-r2))

  y(6) = local(6) - a1
  y(7) = local(7) - a1
  y(8) = local(8) - a1
  y(9) = local(9) - a1

  b1_zmat   = local(10)
  tau1_zmat = local(11)
  b2_zmat   = local(12)

  db1 = b1_zmat - pi
  db2 = b2_zmat - pi

  if (tau1_zmat>pi) then
    tau1 = 2.0_ark*pi - tau1_zmat
  elseif (tau1_zmat<pi) then
    tau1 = -tau1_zmat
  endif

  tau2 = db1 + db2 + tau1
  dtau = tau1 + tau2

  y(10:12) = (/db1, db2, dtau/)

  ! short-range potential

  vshort = force(25) &
         + c2h4_poten_n1_d8( y, force(26:57) ) &
         + c2h4_poten_n2_d8( y, force(58:433) ) &
         + c2h4_poten_n3_d6( y, force(434:1149) ) &
         + c2h4_poten_n4_d6( y, force(1150:2167) ) &
         + c2h4_poten_n5_d6( y, force(2168:2755) ) &
         + c2h4_poten_n6_d6( y, force(2756:2875) )

  f = vshort

end function MLpoten_c2h4_886666_nolong


!---------------------------------------------------------------------


function c2h4_poten_n1_d8(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(32)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(32)
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
end function c2h4_poten_n1_d8


!---------------------------------------------------------------------


function c2h4_poten_n2_d8(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(376)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(376)
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
end function c2h4_poten_n2_d8


!---------------------------------------------------------------------


function c2h4_poten_n3_d6(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(716)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(716)
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
end function c2h4_poten_n3_d6


!---------------------------------------------------------------------


function c2h4_poten_n4_d6(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(1018)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(1018)
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
end function c2h4_poten_n4_d6


!---------------------------------------------------------------------


function c2h4_poten_n5_d6(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(588)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(588)
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
end function c2h4_poten_n5_d6


!---------------------------------------------------------------------


function c2h4_poten_n6_d6(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(120)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(120)
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
end function c2h4_poten_n6_d6


!---------------------------------------------------------------------


function MLpoten_c2h4_lee(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f,W(12,12)

  real(ark) :: y(12)
  real(ark) :: tau4213,tau5124,tau5126,tau6213
  integer(ik) :: k1,k2


  tau4213 =-local(10)
  tau4213 = mod(tau4213+2.0_ark*pi,2.0_ark*pi)
  if (tau4213>pi) tau4213 = tau4213 - 2.0_ark*pi
  !
  tau5124 = local(11)
  tau5124 = mod(tau5124+2.0_ark*pi,2.0_ark*pi)
  !
  tau5126 =-local(12)
  tau5126 = mod(tau5126+2.0_ark*pi,2.0_ark*pi)
  if (tau5126>pi) tau5126 = tau5126 - 2.0_ark*pi
  !
  tau6213 = 2.0_ark*pi-tau5124-tau4213-tau5126
  tau6213 = mod(tau6213+2.0_ark*pi,2.0_ark*pi)

  y(1) = 0.5_ark*sum(local(2:5))-molec%req(2)*2.0_ark
  y(2) = local(1)-molec%req(1)
  y(3) = 0.5_ark*sum(local(6:9))-molec%alphaeq(1)*2.0_ark
  y(4) = (tau5126+tau4213-molec%taueq(1)*2.0_ark)/sqrt(2.0_ark)
  y(5) = 0.5_ark*( local(2)-local(3)-local(4)+local(5) )
  y(6) = 0.5_ark*( local(6)-local(7)-local(8)+local(9) )
  y(7) = (tau6213-tau5124)/sqrt(2.0_ark)
  y(8) = (tau5126-tau4213)/sqrt(2.0_ark)
  y(9) = 0.5_ark*(-local(2)-local(3)+local(4)+local(5) )
  y(10)= 0.5_ark*(-local(6)-local(7)+local(8)+local(9) )
  y(11)= 0.5_ark*(-local(2)+local(3)-local(4)+local(5) )
  y(12)= 0.5_ark*(-local(6)+local(7)-local(8)+local(9) )
  !
  W = 0
  !
  W( 1, 1)= force( 1)
  W( 2, 1)= force( 2)
  W( 3, 1)= force( 3)
  W( 2, 2)= force( 4)
  W( 3, 2)= force( 5)
  W( 3, 3)= force( 6)
  W( 4, 4)= force( 7)
  W( 5, 5)= force( 8)
  W( 6, 5)= force( 9)
  W( 6, 6)= force(10)
  W( 7, 7)= force(11)
  W( 8, 8)= force(12)
  W( 9, 9)= force(13)
  W(10, 9)= force(14)
  W(10,10)= force(15)
  W(11,11)= force(16)
  W(12,11)= force(17)
  W(12,12)= force(18)
  !
  f = 0
  !
  do k1 = 1,12
      f = f + w(k1,k1)*y(k1)**2*0.5_ark
  enddo
  !
  do k1 = 1,12
    do k2 = k1+1,12
      !
      W(k1,k2) = W(k2,k1)
      !
      f = f + w(k2,k1)*y(k1)*y(k2)
      !
    enddo
  enddo
  !
  f = f*1.0e-11/planck/vellgt
  !
  !f = dot_product(y(:),matmul(W,y(:)))*1.0e-11/planck/vellgt
  !
end function MLpoten_c2h4_lee



end module pot_c2h4
