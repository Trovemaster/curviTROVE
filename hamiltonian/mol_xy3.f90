! Internal-to-Cartesian coordinate transformation for NH3 molecule.
! Internal coordinates: r1, r2, r3, s4=(2*beta1-beta2-beta3)/sqrt(6), s5=(beta2-beta3)/sqrt(2), rho,
!                       where r1, r2, r3 - NH1, NH2, NH3 bond stretches,
!                       beta1, beta2, beta3 - H2NH3, H1NH3, H1NH2 angle bends projected on a plane perpendicular to zN,
!                       and rho is zNH1=zNH2=zNH3 angle bend.
! Cartesian coordinates: x(N) =  0
!                        y(N) =  0
!                        z(N) =  0
!                        x(H1) =  0
!                        y(H1) = -r1*sin(rho)
!                        z(H1) =  r1*cos(rho)
!                        x(H2) =  r2*sin(rho)*sin(beta3)
!                        y(H2) = -r2*sin(rho)*cos(beta3)
!                        z(H2) =  r2*cos(rho)
!                        x(H3) = -r3*sin(rho)*sin(beta2)
!                        y(H3) = -r3*sin(rho)*cos(beta2)
!                        z(H3) =  r3*cos(rho)
! Relations between rho, beta1, beta2, beta3 and bending angles alpha1=H2NH3, alpha2=H1NH3, alpha3=H1NH2:
!   cos(rho)**2 = cos(alpha2) - sin(rho)**2*cos(beta2)
!   cos(rho)**2 = cos(alpha3) - sin(rho)**2*cos(beta3)
!   cos(rho)**2 = cos(alpha1) - sin(rho)**2*cos(beta2+beta3)

subroutine internal_to_cartesian_xy3_symbeta_tau(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)   :: internal(molec%nmodes)
  real(ark), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, d3, beta1, beta2, beta3, cm(3), rho, s4, s5

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s4 = internal(4)
  s5 = internal(5)
  rho = internal(6) + real(pi,ark)*0.5_ark
  !rho = internal(6)

  beta1 = sqrt(6.0_ark)/3.0_ark*s4 + 2.0_ark*real(pi,ark)/3.0_ark
  beta2 = -1.0_ark/sqrt(6.0_ark)*s4 + 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark
  beta3 = -1.0_ark/sqrt(6.0_ark)*s4 - 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark

  cartesian(1,1) =  0.0
  cartesian(1,2) =  0.0
  cartesian(1,3) =  0.0

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta3)
  cartesian(3,2) =  d2*sin(rho)*sin(beta3)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta2)
  cartesian(4,2) = -d3*sin(rho)*sin(beta2)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_symbeta_tau


!###############################################################################


! ADF-adapted version of "internal_to_cartesian_xy3_symbeta_tau" subroutine

subroutine internal_to_cartesian_xy3_symbeta_tau_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)  :: internal(molec%nmodes)
  type(adf_realq), intent(out) :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, d3, beta1, beta2, beta3, cm(3), rho, s4, s5

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s4 = internal(4)
  s5 = internal(5)
  rho = internal(6) + real(pi,ark)*0.5_ark
  !rho = internal(6)

  beta1 = sqrt(6.0_ark)/3.0_ark*s4 + 2.0_ark*real(pi,ark)/3.0_ark
  beta2 = -1.0_ark/sqrt(6.0_ark)*s4 + 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark
  beta3 = -1.0_ark/sqrt(6.0_ark)*s4 - 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta3)
  cartesian(3,2) =  d2*sin(rho)*sin(beta3)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta2)
  cartesian(4,2) = -d3*sin(rho)*sin(beta2)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_symbeta_tau_ADF


!###############################################################################


subroutine internal_to_cartesian_xy3_symbeta_sintau(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)   :: internal(molec%nmodes)
  real(ark), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, d3, beta1, beta2, beta3, cm(3), rho, s4, s5, tau

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s4 = internal(4)
  s5 = internal(5)
  tau = asin(internal(6))
  rho = tau + real(pi,ark)*0.5_ark

  beta1 = sqrt(6.0_ark)/3.0_ark*s4 + 2.0_ark*real(pi,ark)/3.0_ark
  beta2 = -1.0_ark/sqrt(6.0_ark)*s4 + 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark
  beta3 = -1.0_ark/sqrt(6.0_ark)*s4 - 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark

  cartesian(1,1) =  0.0
  cartesian(1,2) =  0.0
  cartesian(1,3) =  0.0

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta3)
  cartesian(3,2) =  d2*sin(rho)*sin(beta3)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta2)
  cartesian(4,2) = -d3*sin(rho)*sin(beta2)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_symbeta_sintau


!###############################################################################


! ADF-adapted version of "internal_to_cartesian_xy3_symbeta_sintau" subroutine

subroutine internal_to_cartesian_xy3_symbeta_sintau_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)  :: internal(molec%nmodes)
  type(adf_realq), intent(out) :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, d3, beta1, beta2, beta3, cm(3), rho, s4, s5, tau

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s4 = internal(4)
  s5 = internal(5)
  tau = asin(internal(6))
  rho = tau + real(pi,ark)*0.5_ark

  beta1 = sqrt(6.0_ark)/3.0_ark*s4 + 2.0_ark*real(pi,ark)/3.0_ark
  beta2 = -1.0_ark/sqrt(6.0_ark)*s4 + 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark
  beta3 = -1.0_ark/sqrt(6.0_ark)*s4 - 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta3)
  cartesian(3,2) =  d2*sin(rho)*sin(beta3)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta2)
  cartesian(4,2) = -d3*sin(rho)*sin(beta2)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_symbeta_sintau_ADF


!###############################################################################


subroutine internal_to_cartesian_xy3_symbeta_rho(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)   :: internal(molec%nmodes)
  real(ark), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, d3, beta1, beta2, beta3, cm(3), rho, s4, s5

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s4 = internal(4)
  s5 = internal(5)
  rho = internal(6)

  beta1 = sqrt(6.0_ark)/3.0_ark*s4 + 2.0_ark*real(pi,ark)/3.0_ark
  beta2 = -1.0_ark/sqrt(6.0_ark)*s4 + 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark
  beta3 = -1.0_ark/sqrt(6.0_ark)*s4 - 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark

  cartesian(1,1) =  0.0
  cartesian(1,2) =  0.0
  cartesian(1,3) =  0.0

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta3)
  cartesian(3,2) =  d2*sin(rho)*sin(beta3)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta2)
  cartesian(4,2) = -d3*sin(rho)*sin(beta2)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_symbeta_rho


!###############################################################################


! ADF-adapted version of "internal_to_cartesian_xy3_symbeta_rho" subroutine

subroutine internal_to_cartesian_xy3_symbeta_rho_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)  :: internal(molec%nmodes)
  type(adf_realq), intent(out) :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, d3, beta1, beta2, beta3, cm(3), rho, s4, s5

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s4 = internal(4)
  s5 = internal(5)
  rho = internal(6)

  beta1 = sqrt(6.0_ark)/3.0_ark*s4 + 2.0_ark*real(pi,ark)/3.0_ark
  beta2 = -1.0_ark/sqrt(6.0_ark)*s4 + 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark
  beta3 = -1.0_ark/sqrt(6.0_ark)*s4 - 1.0_ark/sqrt(2.0_ark)*s5 + 2.0_ark*real(pi,ark)/3.0_ark

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta3)
  cartesian(3,2) =  d2*sin(rho)*sin(beta3)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta2)
  cartesian(4,2) = -d3*sin(rho)*sin(beta2)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_symbeta_rho_ADF


!###############################################################################


subroutine internal_to_cartesian_xy3_rsymalpha_tau_pas(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)   :: internal(molec%nmodes)
  real(ark), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, d3, alpha21, alpha32, alpha31, cm(3), delta, rho, beta21, beta31, s1, s2

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s1 = internal(4)
  s2 = internal(5)
  delta = internal(6)

  call find_alpha_from_sindelta(s1,s2,sin(delta),alpha32,alpha31,alpha21)

  rho = real(pi,ark)*0.5_ark + delta

  beta21 = acos((cos(alpha21)-cos(rho)**2)/sin(rho)**2)
  beta31 = acos((cos(alpha31)-cos(rho)**2)/sin(rho)**2)

  cartesian(1,1) =  0.0
  cartesian(1,2) =  0.0
  cartesian(1,3) =  0.0

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta21)
  cartesian(3,2) =  d2*sin(rho)*sin(beta21)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta31)
  cartesian(4,2) = -d3*sin(rho)*sin(beta31)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_rsymalpha_tau_pas


!###############################################################################


! ADF-adapted version of "internal_to_cartesian_xy3_rsymalpha_tau_pas" subroutine

subroutine internal_to_cartesian_xy3_rsymalpha_tau_pas_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)   :: internal(molec%nmodes)
  type(adf_realq), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, d3, alpha21, alpha32, alpha31, cm(3), delta, rho, beta21, beta31, s1, s2

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  s1 = internal(4)
  s2 = internal(5)
  delta = internal(6)

  call find_alpha_from_sindelta_ADF(s1,s2,sin(delta),alpha32,alpha31,alpha21)

  rho = real(pi,ark)*0.5_ark + delta

  beta21 = acos((cos(alpha21)-cos(rho)**2)/sin(rho)**2)
  beta31 = acos((cos(alpha31)-cos(rho)**2)/sin(rho)**2)

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta21)
  cartesian(3,2) =  d2*sin(rho)*sin(beta21)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta31)
  cartesian(4,2) = -d3*sin(rho)*sin(beta31)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_rsymalpha_tau_pas_ADF


!###############################################################################


subroutine find_alpha_from_sindelta(s4,s5,sindelta,alpha1,alpha2,alpha3)
  use adf
  implicit none

  real(ark),intent(in)  :: s4,s5,sindelta
  real(ark),intent(out) :: alpha1,alpha2,alpha3

  real(ark) :: eps,f,rjacob,dx,s6,dx0
  real(ark) :: stadev,ssq,stadev_best,h
  integer(ik) :: iter,itmax,i

  h = 1e-3
  iter = 0
  stadev = 1.e10

  stadev_best = sqrt(epsilon(1.0_rk))*10.0_ark
  itmax = 60

  ! initial values for alpha1 and s6
  alpha1 = 2.0_ark*real(pi,rk)/3.0_ark
  s6 = alpha1*sqrt(3.0_ark)

  do while( iter<itmax .and. stadev>stadev_best )

    iter = iter + 1
    ssq = 0

    ! calculate the function

    f = calc_s2sindelta2(s4,s5,s6)

    eps = f - sindelta**2

    ssq = abs(eps)

    ! calculate derivatives

    rjacob  = ( calc_s2sindelta2(s4,s5,s6+h)-calc_s2sindelta2(s4,s5,s6-h) )/h*0.5_ark
    dx = eps/rjacob
    stadev = sqrt(ssq)

    do i =1,10
      f = calc_s2sindelta2(s4,s5,s6-dx*real(i,ark)/10.0_ark)
      if (abs(f - sindelta**2)>ssq) exit
      ssq = abs(f - sindelta**2)
      dx0 = dx*real(i,ark)/10.0_ark
    enddo

    s6 = s6 - dx0

  enddo

  if (iter==itmax) then
    write(out, '(/a,1x,i3,1x,a)') 'find_alpha_from_sindelta error: could not find solution after', itmax, 'iterations'
    stop
  endif

  alpha1 =(sqrt(2.0_ark)*s6+2.0_ark*s4                 )/sqrt(6.0_ark)
  alpha2 =(sqrt(2.0_ark)*s6-        s4+sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
  alpha3 =(sqrt(2.0_ark)*s6-        s4-sqrt(3.0_ark)*s5)/sqrt(6.0_ark)

  contains

  function calc_s2sindelta2(s4,s5,s6) result(sindelta_2)

    real(ark),intent(in)  :: S4,S5,S6
    real(ark)             :: alpha1,alpha2,alpha3,sindelta_2,tau_2,norm_2

     alpha1 =(sqrt(2.0_ark)*s6+2.0_ark*s4                 )/sqrt(6.0_ark)
     alpha2 =(sqrt(2.0_ark)*s6-        s4+sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
     alpha3 =(sqrt(2.0_ark)*s6-        s4-sqrt(3.0_ark)*s5)/sqrt(6.0_ark)

     tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 &
            +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)

     norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
              2._ark*cos(alpha3)*cos(alpha1)-2._ark*cos(alpha2)+ &
              2._ark*cos(alpha2)*cos(alpha3)-2._ark*cos(alpha1)+ &
              2._ark*cos(alpha2)*cos(alpha1)-2._ark*cos(alpha3)

     sindelta_2 = tau_2/norm_2

  end function calc_s2sindelta2

end subroutine find_alpha_from_sindelta


!###############################################################################


subroutine find_alpha_from_sindelta_ADF(s4,s5,sindelta,alpha1,alpha2,alpha3)
  use adf
  implicit none

  type(adf_realq),intent(in)  :: s4,s5,sindelta
  type(adf_realq),intent(out) :: alpha1,alpha2,alpha3

  type(adf_realq) :: eps,f,rjacob,dx,s6,dx0,h
  real(ark) :: stadev(1:adf_nterms),ssq,stadev_best
  integer(ik) :: iter,itmax,i

  h = real(1e-3,ark)
  iter = 0
  itmax = 10
  stadev = 1.e10
  stadev_best = sqrt(epsilon(1.0_rk))*10.0_ark

  ! initial values for alpha1 and s6
  alpha1 = 2.0_ark*real(pi,rk)/3.0_ark
  s6 = alpha1*sqrt(3.0_ark)

  do while( iter<itmax )!.and. all(stadev(:)>stadev_best) )

    iter = iter + 1
    ssq = 0

    ! calculate the function

    f = calc_s2sindelta2(s4,s5,s6)

    eps = f - sindelta**2

    ssq = abs(eps%v)

    ! calculate derivatives

    rjacob = ( calc_s2sindelta2(s4,s5,s6+h)-calc_s2sindelta2(s4,s5,s6-h) )/h*0.5_ark
    dx = eps / rjacob
    stadev(1:adf_nterms) = sqrt(abs(eps%d(1:adf_nterms)))

    !do i =1,10
    !  f = calc_s2sindelta2(s4,s5,s6-dx*real(i,ark)/10.0_ark)
    !  if (abs(f%v-sindelta%v**2)>ssq) exit
    !  ssq = abs(f%v-sindelta%v**2)
    !  dx0 = dx*real(i,ark)/10.0_ark
    !enddo
    dx0 = dx

    s6 = s6 - dx0

  enddo

  !if (iter==itmax) then
  !  write(out, '(/a,1x,i3,1x,a)') 'find_alpha_from_sindelta_ADF error: could not find solution after', itmax, 'iterations'
  !  stop
  !endif

  alpha1 =(sqrt(2.0_ark)*s6+2.0_ark*s4                 )/sqrt(6.0_ark)
  alpha2 =(sqrt(2.0_ark)*s6-        s4+sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
  alpha3 =(sqrt(2.0_ark)*s6-        s4-sqrt(3.0_ark)*s5)/sqrt(6.0_ark)

  contains

  function calc_s2sindelta2(s4,s5,s6) result(sindelta_2)

    type(adf_realq),intent(in)  :: s4,s5,s6
    type(adf_realq)             :: alpha1,alpha2,alpha3,sindelta_2,tau_2,norm_2

     alpha1 =(sqrt(2.0_ark)*s6+2.0_ark*s4                 )/sqrt(6.0_ark)
     alpha2 =(sqrt(2.0_ark)*s6-        s4+sqrt(3.0_ark)*s5)/sqrt(6.0_ark)
     alpha3 =(sqrt(2.0_ark)*s6-        s4-sqrt(3.0_ark)*s5)/sqrt(6.0_ark)

     tau_2 = 1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 &
           + 2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3)

     norm_2 = sin(alpha3)**2+sin(alpha2)**2+sin(alpha1)**2+&
              2._ark*cos(alpha3)*cos(alpha1)-2._ark*cos(alpha2)+ &
              2._ark*cos(alpha2)*cos(alpha3)-2._ark*cos(alpha1)+ &
              2._ark*cos(alpha2)*cos(alpha1)-2._ark*cos(alpha3)

     sindelta_2 = tau_2/norm_2

  end function calc_s2sindelta2

end subroutine find_alpha_from_sindelta_ADF


!###############################################################################


! Internal-to-Cartesian coordinate transformation for NH3 molecule.
! Internal coordinates: r1, r2, r3, alpha1, alpha2, alpha3,
!                       where r1, r2, r3 - NH1, NH2, NH3 bond stretches,
!                       and alpha1, alpha2, alpha3 - H2NH3, H1NH3, H1NH2 angle bends.

subroutine internal_to_cartesian_xy3_ralpha_zmat(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)   :: internal(molec%nmodes)
  real(ark), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, d3, alpha21, alpha32, alpha31, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  alpha21 = internal(6)
  alpha31 = internal(5)
  alpha32 = internal(4)

  cartesian(1,1) =  0.0
  cartesian(1,2) =  0.0
  cartesian(1,3) =  0.0

  cartesian(2,1) =  0.0
  cartesian(2,2) =  0.0
  cartesian(2,3) = -d1

  cartesian(3,1) =  d2*sin(alpha21)
  cartesian(3,2) =  0.0
  cartesian(3,3) = -d2*cos(alpha21)

  cartesian(4,1) =  d3*(-cos(alpha31)*cos(alpha21)+cos(alpha32))/sin(alpha21)
  cartesian(4,2) =  d3*sqrt(1.0_ark-cos(alpha21)**2-cos(alpha31)**2+2.0_ark*cos(alpha31)*cos(alpha21)*cos(alpha32)-cos(alpha32)**2)/sin(alpha21)
  cartesian(4,3) = -d3*cos(alpha31)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_ralpha_zmat


!###############################################################################


! ADF-adapted version of "internal_to_cartesian_xy3_ralpha_zmat" subroutine

subroutine internal_to_cartesian_xy3_ralpha_zmat_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)  :: internal(molec%nmodes)
  type(adf_realq), intent(out) :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, d3, alpha21, alpha32, alpha31, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  alpha21 = internal(6)
  alpha31 = internal(5)
  alpha32 = internal(4)

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  0.0_ark
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) = -d1

  cartesian(3,1) =  d2*sin(alpha21)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) = -d2*cos(alpha21)

  cartesian(4,1) =  d3*(-cos(alpha31)*cos(alpha21)+cos(alpha32))/sin(alpha21)
  cartesian(4,2) =  d3*sqrt(1.0_ark-cos(alpha21)**2-cos(alpha31)**2+2.0_ark*cos(alpha31)*cos(alpha21)*cos(alpha32)-cos(alpha32)**2)/sin(alpha21)
  cartesian(4,3) = -d3*cos(alpha31)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product(cartesian(1:natoms,ix), molec%masses(1:natoms)) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_ralpha_zmat_ADF


!###############################################################################


! Internal-to-Cartesian coordinate transformation for NH3 molecule.
! Internal coordinates: r1, r2, r3, alpha1, alpha2, alpha3,
!                       where r1, r2, r3 - NH1, NH2, NH3 bond stretches,
!                       and alpha1, alpha2, alpha3 - H2NH3, H1NH3, H1NH2 angle bends.
!
! While the coordinates are same as used in "internal_to_cartesian_xy3_ralpha_zmat" subroutine, the orientation of axes system is different

subroutine internal_to_cartesian_xy3_ralpha_pas(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)   :: internal(molec%nmodes)
  real(ark), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, d3, alpha21, alpha32, alpha31, cm(3), tau_2, norm_2, sinrho, rho, beta21, beta31

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  alpha21 = internal(6)
  alpha31 = internal(5)
  alpha32 = internal(4)

  tau_2 = 1.0_ark-cos(alpha21)**2-cos(alpha31)**2-cos(alpha32)**2+2.0_ark*cos(alpha21)*cos(alpha31)*cos(alpha32)
  norm_2 = sin(alpha21)**2+sin(alpha31)**2+sin(alpha32)**2+&
               2.0_ark*cos(alpha21)*cos(alpha32)-2.0_ark*cos(alpha31)+ &
               2.0_ark*cos(alpha31)*cos(alpha21)-2.0_ark*cos(alpha32)+ &
               2.0_ark*cos(alpha31)*cos(alpha32)-2.0_ark*cos(alpha21)

  sinrho = sqrt(tau_2)/sqrt(norm_2)
  rho = real(pi,ark)*0.5_ark - asin(sinrho)

  beta21 = acos((cos(alpha21)-cos(rho)**2)/sin(rho)**2)
  beta31 = acos((cos(alpha31)-cos(rho)**2)/sin(rho)**2)

  cartesian(1,1) =  0.0
  cartesian(1,2) =  0.0
  cartesian(1,3) =  0.0

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta21)
  cartesian(3,2) =  d2*sin(rho)*sin(beta21)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta31)
  cartesian(4,2) = -d3*sin(rho)*sin(beta31)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_ralpha_pas


!###############################################################################


! ADF-adapted version of "internal_to_cartesian_xy3_ralpha_pas" subroutine

subroutine internal_to_cartesian_xy3_ralpha_pas_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)  :: internal(molec%nmodes)
  type(adf_realq), intent(out) :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, d3, alpha21, alpha32, alpha31, cm(3), tau_2, norm_2, sinrho, rho, beta21, beta31

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  d3 = internal(3)
  alpha21 = internal(6)
  alpha31 = internal(5)
  alpha32 = internal(4)

  tau_2 = 1.0_ark-cos(alpha21)**2-cos(alpha31)**2-cos(alpha32)**2+2.0_ark*cos(alpha21)*cos(alpha31)*cos(alpha32)
  norm_2 = sin(alpha21)**2+sin(alpha31)**2+sin(alpha32)**2+&
               2.0_ark*cos(alpha21)*cos(alpha32)-2.0_ark*cos(alpha31)+ &
               2.0_ark*cos(alpha31)*cos(alpha21)-2.0_ark*cos(alpha32)+ &
               2.0_ark*cos(alpha31)*cos(alpha32)-2.0_ark*cos(alpha21)

  sinrho = sqrt(tau_2)/sqrt(norm_2)
  rho = real(pi,ark)*0.5_ark - asin(sinrho)

  beta21 = acos((cos(alpha21)-cos(rho)**2)/sin(rho)**2)
  beta31 = acos((cos(alpha31)-cos(rho)**2)/sin(rho)**2)

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)*cos(beta21)
  cartesian(3,2) =  d2*sin(rho)*sin(beta21)
  cartesian(3,3) =  d2*cos(rho)

  cartesian(4,1) =  d3*sin(rho)*cos(beta31)
  cartesian(4,2) = -d3*sin(rho)*sin(beta31)
  cartesian(4,3) =  d3*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product(cartesian(1:natoms,ix), molec%masses(1:natoms)) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy3_ralpha_pas_ADF
