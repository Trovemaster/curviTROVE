subroutine internal_to_cartesian_xy2_rrho(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, rho, alpha, alpha_half, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  rho = internal(3)
  alpha = real(pi,ark) - rho
  alpha_half = 0.5_ark*alpha

  cartesian(1,1) = 0.0
  cartesian(1,2) = 0.0
  cartesian(1,3) = 0.0

  cartesian(3,1) = d1*sin(alpha_half)
  cartesian(3,2) = 0.0
  cartesian(3,3) = d1*cos(alpha_half)

  cartesian(2,1) = -d2*sin(alpha_half)
  cartesian(2,2) =  0.0
  cartesian(2,3) =  d2*cos(alpha_half)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy2_rrho


!################################################################################


! ADF-adapted version of "internal_to_cartesian_xy2_rrho" subroutine.

subroutine internal_to_cartesian_xy2_rrho_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, rho, alpha, alpha_half, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  rho = internal(3)
  alpha = real(pi,ark) - rho
  alpha_half = 0.5_ark*alpha

  cartesian(1,1) = 0.0_ark
  cartesian(1,2) = 0.0_ark
  cartesian(1,3) = 0.0_ark

  cartesian(3,1) = d1*sin(alpha_half)
  cartesian(3,2) = 0.0_ark
  cartesian(3,3) = d1*cos(alpha_half)

  cartesian(2,1) = -d2*sin(alpha_half)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) =  d2*cos(alpha_half)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy2_rrho_ADF


!################################################################################


subroutine internal_to_cartesian_xy2_rrho_zmat(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, rho, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  rho = internal(3)*0.5_ark

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) = -d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) =  d2*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy2_rrho_zmat


!################################################################################


! ADF-adapted version of "internal_to_cartesian_xy2_rrho_zmat" subroutine.

subroutine internal_to_cartesian_xy2_rrho_zmat_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, rho, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  rho = internal(3)*0.5_ark

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  d1*sin(rho)
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) = -d1*cos(rho)

  cartesian(3,1) =  d2*sin(rho)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) =  d2*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy2_rrho_zmat_ADF


!################################################################################


subroutine internal_to_cartesian_xy2_rrho_testsing(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, d1, d2, rho, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  rho = internal(3)

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  0.0_ark
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) = -d1

  cartesian(3,1) =  d2*sin(rho)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) =  d2*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy2_rrho_testsing


!################################################################################


! ADF-adapted version of "internal_to_cartesian_xy2_rrho_testsing" subroutine.

subroutine internal_to_cartesian_xy2_rrho_testsing_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: d1, d2, rho, cm(3)

  natoms = molec%natoms

  d1 = internal(1)
  d2 = internal(2)
  rho = internal(3)

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  0.0_ark
  cartesian(2,2) =  0.0_ark
  cartesian(2,3) = -d1

  cartesian(3,1) =  d2*sin(rho)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) =  d2*cos(rho)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy2_rrho_testsing_ADF
