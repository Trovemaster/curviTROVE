! Internal-to-Cartesian coordinate transformation for ethylene-type molecule.
!
!  H1       H3
!   \      /
!    C1 - C2
!   /      \
!  H2       H4
!
! Internal coordinates:
!   bond lengths r_{C_1C_2}, r_{C_1H_1}, r_{C_1H_2}, r_{C_2H_3}, r_{C_2H_4}
!   bond angles alpha_{H_1C_1C_2}, alpha_{H_2C_1C_2}, alpha_{H_3C_2C_1}, alpha_{H_4C_2C_1}
!   HCH book-type dihedral angles: beta_{H_1C_1C_2H_2}, beta_{H_3C_2C_1H_4} (0 < {beta_{H_1C_1C_2H_2}, beta_{H_3C_2C_1H_4}} < 360 deg)
!   symmetrized HCH-HCH internal rotation angle: delta=2*tau_{H_1C_1C_2H_3}-beta_{H_1C_1C_2H_2}+beta_{H_3C_2C_1H_4}

subroutine internal_to_cartesian_c2h4_2beta_1tau(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)   :: internal(molec%nmodes)
  real(ark), intent(out)  :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  real(ark) :: r0, r1, r2, r3, r4, alpha1, alpha2, alpha3, alpha4, beta1, beta2, tau, delta, cm(3)

  natoms = molec%natoms

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

  tau = 0.5_ark*(delta+beta1-beta2)

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) = -r0*0.5_ark

  cartesian(2,1) =  0.0_ark
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  r0*0.5_ark

  cartesian(3,1) =  r1*sin(alpha1)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) =  r1*cos(alpha1)-r0*0.5_ark

  cartesian(4,1) =  r2*sin(alpha2)*cos(beta1)
  cartesian(4,2) =  r2*sin(alpha2)*sin(beta1)
  cartesian(4,3) =  r2*cos(alpha2)-r0*0.5_ark

  cartesian(5,1) =  r3*sin(alpha3)*cos(tau)
  cartesian(5,2) =  r3*sin(alpha3)*sin(tau)
  cartesian(5,3) = -r3*cos(alpha3)+r0*0.5_ark

  cartesian(6,1) =  r4*sin(alpha4)*(cos(tau)*cos(beta2)-sin(tau)*sin(beta2))
  cartesian(6,2) =  r4*sin(alpha4)*(cos(tau)*sin(beta2)+sin(tau)*cos(beta2))
  cartesian(6,3) = -r4*cos(alpha4)+r0*0.5_ark

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h4_2beta_1tau


!###############################################################################


! ADF-adapted version of "internal_to_cartesian_c2h4_2beta_1tau" subroutine.

subroutine internal_to_cartesian_c2h4_2beta_1tau_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)  :: internal(molec%nmodes)
  type(adf_realq), intent(out) :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: r0, r1, r2, r3, r4, alpha1, alpha2, alpha3, alpha4, beta1, beta2, tau, delta, cm(3)

  natoms = molec%natoms

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

  tau = 0.5_ark*(delta+beta1-beta2)

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) = -r0*0.5_ark

  cartesian(2,1) =  0.0_ark
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  r0*0.5_ark

  cartesian(3,1) =  r1*sin(alpha1)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) =  r1*cos(alpha1)-r0*0.5_ark

  cartesian(4,1) =  r2*sin(alpha2)*cos(beta1)
  cartesian(4,2) =  r2*sin(alpha2)*sin(beta1)
  cartesian(4,3) =  r2*cos(alpha2)-r0*0.5_ark

  cartesian(5,1) =  r3*sin(alpha3)*cos(tau)
  cartesian(5,2) =  r3*sin(alpha3)*sin(tau)
  cartesian(5,3) = -r3*cos(alpha3)+r0*0.5_ark

  cartesian(6,1) =  r4*sin(alpha4)*(cos(tau)*cos(beta2)-sin(tau)*sin(beta2))
  cartesian(6,2) =  r4*sin(alpha4)*(cos(tau)*sin(beta2)+sin(tau)*cos(beta2))
  cartesian(6,3) = -r4*cos(alpha4)+r0*0.5_ark

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h4_2beta_1tau_ADF
