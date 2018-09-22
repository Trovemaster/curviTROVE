subroutine internal_to_cartesian_h2o2_ralphatau(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, R, r1, r2, alpha1, alpha2, tau, cm(3)

  natoms = molec%natoms

  R      = internal(1)
  r1     = internal(2)
  r2     = internal(3)
  alpha1 = internal(4)
  alpha2 = internal(5)
  tau    = internal(6)

  ! C1
  cartesian(1,1) = 0.0
  cartesian(1,2) = 0.0
  cartesian(1,3) = 0.0

  ! C2
  cartesian(2,1) = 0.0
  cartesian(2,2) = 0.0
  cartesian(2,3) = R

  ! H1
  cartesian(3,1) = r1*sin(alpha1)
  cartesian(3,2) = 0.0
  cartesian(3,3) = r1*cos(alpha1)

  ! H2
  cartesian(4,1) =  r2*sin(alpha2)*cos(tau)
  cartesian(4,2) =  r2*sin(alpha2)*sin(tau)
  cartesian(4,3) =  R-r2*cos(alpha2)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_h2o2_ralphatau


!################################################################################


! ADF-adapted version of "internal_to_cartesian_h2o2_ralphatau" subroutine.

subroutine internal_to_cartesian_h2o2_ralphatau_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: R, r1, r2, alpha1, alpha2, tau, cm(3)

  natoms = molec%natoms

  R      = internal(1)
  r1     = internal(2)
  r2     = internal(3)
  alpha1 = internal(4)
  alpha2 = internal(5)
  tau    = internal(6)

  ! C1
  cartesian(1,1) = 0.0_ark
  cartesian(1,2) = 0.0_ark
  cartesian(1,3) = 0.0_ark

  ! C2
  cartesian(2,1) = 0.0_ark
  cartesian(2,2) = 0.0_ark
  cartesian(2,3) = R

  ! H1
  cartesian(3,1) = r1*sin(alpha1)
  cartesian(3,2) = 0.0_ark
  cartesian(3,3) = r1*cos(alpha1)

  ! H2
  cartesian(4,1) =  r2*sin(alpha2)*cos(tau)
  cartesian(4,2) =  r2*sin(alpha2)*sin(tau)
  cartesian(4,3) =  R-r2*cos(alpha2)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_h2o2_ralphatau_ADF
