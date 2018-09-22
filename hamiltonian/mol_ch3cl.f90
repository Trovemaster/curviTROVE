! Internal-to-Cartesian coordinate transformation for CH3Cl molecule.
! Internal coordinates: r0, r1, r2, r3, beta1, beta2, beta3, stau1=(2*tau1-tau2-tau3)/sqrt(6), stau2=(tau2-tau3)/sqrt(2),
!                       where r0, r1, r2, r3 - CCl, CH1, CH2, CH3 bond stretches,
!                       beta1, beta2, beta3 - H1CCl, H2CCl, H3CCl angle bends,
!                       and tau1, tau2, tau3 - H2CH3, H1CH3, H2CH3 angle bends projected on a plane perpendicular to CCl.

subroutine internal_to_cartesian_ch3cl_beta_symtau(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  real(ark) :: r0, r1, r2, r3, beta1, beta2, beta3, dbeta1, dbeta2, dbeta3, stau1, stau2, tau1, tau2, tau3, cm(3)

  natoms = molec%natoms

  r0 = internal(1)     ! r_{CCl}
  r1 = internal(2)     ! r_{CH1}
  r2 = internal(3)     ! r_{CH2}
  r3 = internal(4)     ! r_{CH3}
  beta1 = internal(5)  ! beta_{H1CCl}
  beta2 = internal(6)  ! beta_{H2CCl}
  beta3 = internal(7)  ! beta_{H3CCl}
  stau1 = internal(8)  ! stau1=(2*tau1-tau2-tau3)/sqrt(6)
  stau2 = internal(9)  ! stau2=(tau2-tau3)/sqrt(2)

  dbeta1 = real(pi,ark) - beta1
  dbeta2 = real(pi,ark) - beta2
  dbeta3 = real(pi,ark) - beta3

  tau1 = sqrt(6.0_ark)/3.0_ark*stau1 + 2.0_ark*real(pi,ark)/3.0_ark
  tau2 = -1.0_ark/sqrt(6.0_ark)*stau1 + 1.0_ark/sqrt(2.0_ark)*stau2 + 2.0_ark*real(pi,ark)/3.0_ark
  tau3 = -1.0_ark/sqrt(6.0_ark)*stau1 - 1.0_ark/sqrt(2.0_ark)*stau2 + 2.0_ark*real(pi,ark)/3.0_ark

  ! C
  cartesian(1,1) =  0.0
  cartesian(1,2) =  0.0
  cartesian(1,3) =  0.0

  ! Cl
  cartesian(2,1) =  0.0
  cartesian(2,2) =  0.0
  cartesian(2,3) =  r0

  ! H1
  cartesian(3,1) =  r1*sin(dbeta1)
  cartesian(3,2) =  0.0
  cartesian(3,3) = -r1*cos(dbeta1)

  ! H2
  cartesian(4,1) =  r2*sin(dbeta2)*cos(tau3)
  cartesian(4,2) =  r2*sin(dbeta2)*sin(tau3)
  cartesian(4,3) = -r2*cos(dbeta2)

  ! H3
  cartesian(5,1) =  r3*sin(dbeta3)*cos(tau2)
  cartesian(5,2) = -r3*sin(dbeta3)*sin(tau2)
  cartesian(5,3) = -r3*cos(dbeta3)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_ch3cl_beta_symtau


!################################################################################


! ADF-adapted version of "internal_to_cartesian_ch3cl_beta_symtau" subroutine.

subroutine internal_to_cartesian_ch3cl_beta_symtau_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: r0, r1, r2, r3, beta1, beta2, beta3, dbeta1, dbeta2, dbeta3, stau1, stau2, tau1, tau2, tau3, cm(3)

  natoms = molec%natoms

  r0 = internal(1)     ! r_{CCl}
  r1 = internal(2)     ! r_{CH1}
  r2 = internal(3)     ! r_{CH2}
  r3 = internal(4)     ! r_{CH3}
  beta1 = internal(5)  ! beta_{H1CCl}
  beta2 = internal(6)  ! beta_{H2CCl}
  beta3 = internal(7)  ! beta_{H3CCl}
  stau1 = internal(8)  ! stau1=(2*tau1-tau2-tau3)/sqrt(6)
  stau2 = internal(9)  ! stau2=(tau2-tau3)/sqrt(2)

  dbeta1 = real(pi,ark) - beta1
  dbeta2 = real(pi,ark) - beta2
  dbeta3 = real(pi,ark) - beta3

  tau1 = sqrt(6.0_ark)/3.0_ark*stau1 + 2.0_ark*real(pi,ark)/3.0_ark
  tau2 = -1.0_ark/sqrt(6.0_ark)*stau1 + 1.0_ark/sqrt(2.0_ark)*stau2 + 2.0_ark*real(pi,ark)/3.0_ark
  tau3 = -1.0_ark/sqrt(6.0_ark)*stau1 - 1.0_ark/sqrt(2.0_ark)*stau2 + 2.0_ark*real(pi,ark)/3.0_ark

  ! C
  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  ! Cl
  cartesian(2,1) =  0.0_ark
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  r0

  ! H1
  cartesian(3,1) =  r1*sin(dbeta1)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) = -r1*cos(dbeta1)

  ! H2
  cartesian(4,1) =  r2*sin(dbeta2)*cos(tau3)
  cartesian(4,2) =  r2*sin(dbeta2)*sin(tau3)
  cartesian(4,3) = -r2*cos(dbeta2)

  ! H3
  cartesian(5,1) =  r3*sin(dbeta3)*cos(tau2)
  cartesian(5,2) = -r3*sin(dbeta3)*sin(tau2)
  cartesian(5,3) = -r3*cos(dbeta3)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_ch3cl_beta_symtau_ADF
