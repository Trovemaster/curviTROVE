!
!  This unit is for a user defined potential 
!
module pot_user
  use accuracy
  use moltype

  implicit none

  public MLdipole,MLpoten,ML_MEP

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
 contains
 !
 !
 function ML_MEP(dim,rho)  result(f)

  integer(ik),intent(in) ::  dim
  real(ark),intent(in)   ::  rho
  real(ark)              ::  f(dim)
  !
  if (dim/=3) stop 'Illegal size of the function - must be 3'
  !
  f(:) = molec%local_eq(:)
  f(molec%Ncoords) = rho

 end function ML_MEP
 !
 !
 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
   !
   f = 0.0
   !
 end subroutine MLdipole
 !
 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   f = MLpoten_c2h4_886666(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !
 !
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



end module pot_user
