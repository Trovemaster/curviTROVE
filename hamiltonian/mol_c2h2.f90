

!################################################################################

subroutine internal_to_cartesian_c2h2_rqxqy(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: pm=1, ix, natoms, iatom
  real(ark) :: r0, r1, r2, qx1, qx2, qy1, qy2, fm, cm(3), wx1_new, wx2_new, wy1_new, wy2_new

  natoms = molec%natoms

  !call fromlocal2cartesian(molec, pm, internal(1:7), cartesian)

  r0 = internal(1)
  r1 = internal(2)
  r2 = internal(3)

  wx1_new = internal(4)
  wy1_new = internal(5)
  wx2_new = internal(6)
  wy2_new = internal(7)

  qx1 = wy1_new*r1
  qy1 = -wx1_new*r1
  qx2 = wy2_new*r2
  qy2 = -wx2_new*r2

  ! C1
  cartesian(1,1) = 0.0
  cartesian(1,2) = 0.0
  cartesian(1,3) = 0.0

  ! C2
  cartesian(2,1) = 0.0
  cartesian(2,2) = 0.0
  cartesian(2,3) = r0

  ! H1
  cartesian(3,1) =   qx1
  cartesian(3,2) =   qy1
  cartesian(3,3) =  -sqrt(r1**2-qx1**2-qy1**2)

  ! H2
  cartesian(4,1) =  qx2
  cartesian(4,2) =  qy2
  cartesian(4,3) =  sqrt(r2**2-qx2**2-qy2**2)+r0

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_rqxqy


!################################################################################


! ADF-adapted version of "internal_to_cartesian_c2h2_rqxqy" subroutine.

subroutine internal_to_cartesian_c2h2_rqxqy_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: pm=1, ix, natoms, iatom
  real(ark) :: fm
  type(adf_realq) :: r0, r1, r2, qx1, qx2, qy1, qy2, cm(3), wx1_new, wx2_new, wy1_new, wy2_new

  natoms = molec%natoms

  !call fromlocal2cartesian_ADF(molec, pm, internal(1:7), cartesian)

  r0 = internal(1)
  r1 = internal(2)
  r2 = internal(3)
  wx1_new = internal(4)
  wy1_new = internal(5)
  wx2_new = internal(6)
  wy2_new = internal(7)

  qx1 = wy1_new*r1
  qy1 = -wx1_new*r1
  qx2 = wy2_new*r2
  qy2 = -wx2_new*r2

  ! C1
  cartesian(1,1) = 0.0_ark
  cartesian(1,2) = 0.0_ark
  cartesian(1,3) = 0.0_ark

  ! C2
  cartesian(2,1) = 0.0_ark
  cartesian(2,2) = 0.0_ark
  cartesian(2,3) = r0

  ! H1
  cartesian(3,1) =   qx1
  cartesian(3,2) =   qy1
  cartesian(3,3) =  -sqrt(r1**2-qx1**2-qy1**2)

  ! H2
  cartesian(4,1) =  qx2
  cartesian(4,2) =  qy2
  cartesian(4,3) =  sqrt(r2**2-qx2**2-qy2**2)+r0

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_rqxqy_ADF


!################################################################################


subroutine internal_to_cartesian_c2h2_rxyz(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, R, x1, z1, x2, y2, y1, z2, cm(3), x1s, x2s, y

  natoms = molec%natoms

  R   = internal(1)
  x1s = internal(2)
  x2s = internal(3)
  z1  = internal(4)
  z2  = internal(5)
  y   = internal(6)

  x1 = x1s + x2s
  x2 = x1s - x2s
  y1 =  y
  y2 = -y

  ! C1
  cartesian(1,1) = 0.0
  cartesian(1,2) = 0.0
  cartesian(1,3) = 0.0

  ! C2
  cartesian(2,1) = 0.0
  cartesian(2,2) = 0.0
  cartesian(2,3) = R

  ! H1
  cartesian(3,1) =  x1
  cartesian(3,2) =  y1
  cartesian(3,3) = -z1

  ! H2
  cartesian(4,1) =  x2
  cartesian(4,2) =  y2
  cartesian(4,3) =  z2+R

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_rxyz


!################################################################################


! ADF-adapted version of "internal_to_cartesian_c2h2_rxyz" subroutine.

subroutine internal_to_cartesian_c2h2_rxyz_ADF(molec, internal, cartesian)
  use adf

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: R, x1, z1, x2, y2, y1, z2, cm(3), x1s, x2s, y

  natoms = molec%natoms

  R   = internal(1)
  x1s = internal(2)
  x2s = internal(3)
  z1  = internal(4)
  z2  = internal(5)
  y   = internal(6)

  x1 = x1s + x2s
  x2 = x1s - x2s
  y1 =  y
  y2 = -y

  ! C1
  cartesian(1,1) = 0.0_ark
  cartesian(1,2) = 0.0_ark
  cartesian(1,3) = 0.0_ark

  ! C2
  cartesian(2,1) = 0.0_ark
  cartesian(2,2) = 0.0_ark
  cartesian(2,3) = R

  ! H1
  cartesian(3,1) =  x1
  cartesian(3,2) =  y1
  cartesian(3,3) = -z1

  ! H2
  cartesian(4,1) =  x2
  cartesian(4,2) =  y2
  cartesian(4,3) =  z2+R

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_rxyz_ADF


!################################################################################


subroutine internal_to_cartesian_c2h2_rsymalphatau(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, rho, phi, R, r1, r2, alpha1, alpha2, tau, cm(3), theta

  natoms = molec%natoms

  R      = internal(1)
  r1     = internal(2)
  r2     = internal(3)
  rho    = internal(4)
  phi    = internal(5)
  tau    = internal(6)

  alpha1 = pi - rho!rho * cos(phi)
  alpha2 = pi - phi!rho * sin(phi)

  ! C1
  cartesian(1,1) = 0.0
  cartesian(1,2) = 0.0
  cartesian(1,3) = 0.0

  ! C2
  cartesian(2,1) = 0.0
  cartesian(2,2) = 0.0
  cartesian(2,3) = R

  ! H1
  cartesian(3,1) =  r1*sin(alpha1)
  cartesian(3,2) =  0.0
  cartesian(3,3) =  r1*cos(alpha1)

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

end subroutine internal_to_cartesian_c2h2_rsymalphatau


!################################################################################


! ADF-adapted version of "internal_to_cartesian_c2h2_rsymalphatau" subroutine.

subroutine internal_to_cartesian_c2h2_rsymalphatau_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: R, r1, r2, alpha1, alpha2, tau, cm(3), rho, phi, theta

  natoms = molec%natoms

  R      = internal(1)
  r1     = internal(2)
  r2     = internal(3)
  rho    = internal(4)
  phi    = internal(5)
  tau    = internal(6)

  alpha1 = pi - rho!rho * cos(phi)
  alpha2 = pi - phi!rho * sin(phi)

  ! C1
  cartesian(1,1) = 0.0_ark
  cartesian(1,2) = 0.0_ark
  cartesian(1,3) = 0.0_ark

  ! C2
  cartesian(2,1) = 0.0_ark
  cartesian(2,2) = 0.0_ark
  cartesian(2,3) = R

  ! H1
  cartesian(3,1) =  r1*sin(alpha1)
  cartesian(3,2) =  0.0_ark
  cartesian(3,3) =  r1*cos(alpha1)

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

end subroutine internal_to_cartesian_c2h2_rsymalphatau_ADF


!################################################################################


subroutine internal_to_cartesian_c2h2_ralpha(molec, internal, cartesian)

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
  cartesian(4,1) =  r2*sin(alpha2)
  cartesian(4,2) =  r2*cos(alpha2)*sin(tau)
  cartesian(4,3) =  R-r2*cos(alpha2)*cos(tau)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_ralpha


!################################################################################


! ADF-adapted version of "internal_to_cartesian_c2h2_ralpha" subroutine.

subroutine internal_to_cartesian_c2h2_ralpha_ADF(molec, internal, cartesian)
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
  cartesian(4,1) =  r2*sin(alpha2)
  cartesian(4,2) =  r2*cos(alpha2)*sin(tau)
  cartesian(4,3) =  R-r2*cos(alpha2)*cos(tau)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_ralpha_ADF


!################################################################################


subroutine internal_to_cartesian_c2h2_r4alpha(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm, R, r1, r2, alpha1, alpha2, alpha3, alpha4, cm(3)

  natoms = molec%natoms

  R      = internal(1)
  r1     = internal(2)
  r2     = internal(3)
  alpha1 = internal(4)
  alpha2 = internal(5)
  alpha3 = internal(6)
  alpha4 = internal(7)

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
  cartesian(3,2) = -r1*cos(alpha1)*sin(alpha2)
  cartesian(3,3) = -r1*cos(alpha1)*cos(alpha2)

  ! H2
  cartesian(4,1) =  r2*sin(alpha3)
  cartesian(4,2) =  -r2*cos(alpha3)*sin(alpha4)
  cartesian(4,3) =  R+r2*cos(alpha3)*cos(alpha4)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_r4alpha


!################################################################################


! ADF-adapted version of "internal_to_cartesian_c2h2_r4alpha" subroutine.

subroutine internal_to_cartesian_c2h2_r4alpha_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) :: R, r1, r2, alpha1, alpha2, alpha3, alpha4, cm(3)

  natoms = molec%natoms

  R      = internal(1)
  r1     = internal(2)
  r2     = internal(3)
  alpha1 = internal(4)
  alpha2 = internal(5)
  alpha3 = internal(6)
  alpha4 = internal(7)

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
  cartesian(3,2) = -r1*cos(alpha1)*sin(alpha2)
  cartesian(3,3) = -r1*cos(alpha1)*cos(alpha2)

  ! H2
  cartesian(4,1) =  r2*sin(alpha3)
  cartesian(4,2) =  -r2*cos(alpha3)*sin(alpha4)
  cartesian(4,3) =  R+r2*cos(alpha3)*cos(alpha4)

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = dot_product( cartesian(1:natoms,ix), molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_c2h2_r4alpha_ADF
