! Same as "poten_xy4_alpha_ADF" in pot_xy4.f90, except it uses the de-symmetrized representation of the PES (defined in "pot_xy4_diff_V")

subroutine poten_xy4_alpha_ADF_noreexp(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: i, k, irank
  real(ark) :: coefs_sym(289), amorse1, r1eq
  type(adf_realq) :: y1, y2, y3, y4, y5, y6, y7, y8, y9

  irank = 1

  r1eq = func%params(1,irank)
  amorse1 = func%params(2,irank)

  coefs_sym = 0
  do i=3, func%nparams(irank)
    k = func%iparams(1,i,irank)
    coefs_sym(k) = func%params(i,irank)
  enddo

  if (trim(molec%coord_transform)=='XY4_RSYMALPHA') then

    y1 = 1.0_ark-exp(-amorse1*(internal(1)-r1eq))
    y2 = 1.0_ark-exp(-amorse1*(internal(2)-r1eq))
    y3 = 1.0_ark-exp(-amorse1*(internal(3)-r1eq))
    y4 = 1.0_ark-exp(-amorse1*(internal(4)-r1eq))

    y5 = internal(5)
    y6 = internal(6)
    y7 = internal(7)
    y8 = internal(8)
    y9 = internal(9)

  else

    write(out, '(/a,a,a)') 'poten_xy4_alpha_ADF_noreexp error: coordinate type ="',trim(molec%coord_transform), '" is not supported'
    stop

  endif

  f(irank) = pot_xy4_diff_V_init_poten(molec%nmodes, (/y1,y2,y3,y4,y5,y6,y7,y8,y9/), coefs_sym)

end subroutine poten_xy4_alpha_ADF_noreexp


!################################################################################


function pot_xy4_diff_V_init_poten(nmodes, local, f) result(res)
  use adf
  implicit none

  integer(ik), intent(in) :: nmodes
  type(adf_realq), intent(in) :: local(nmodes)
  real(ark), intent(in) :: f(:)
  type(adf_realq) :: res

  integer(ik), parameter :: maxnterms = 4000
  integer(ik) :: iterm, nterms, terms(9,maxnterms), imode
  real(ark) :: coefs(maxnterms)
  type(adf_realq) :: prod

  if (nmodes/=9) then
    write(out, '(/a,1x,i3,1x,a)') 'pot_xy4_diff_V_init_poten error: number of internal coordinates =', nmodes, '!= 9'
    stop
  endif

  ! compute all symmetry-equivalent expansion coefficients

  nterms = 0
  terms = 0
  coefs = 0.0

#include 'pot_xy4_2_part2.f90'

  if (nterms==0) then
    write(out, '(/a)') 'pot_xy4_diff_V_init_poten error: part 2 of the file pot_xy4_2.f90 was not included at compilation, recompile module "hamiltonian" using -fpp or -cpp flag'
    stop
  endif

  if (nterms>maxnterms) then
    write(out, '(/a,1x,i6,1x,a,1x,i6)') 'pot_xy4_diff_V_init_poten error: number of terms =', nterms, 'exceeds max number =', maxnterms
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

end function pot_xy4_diff_V_init_poten


!################################################################################


! Same as "dipole_xy4_alpha_ADF" in pot_xy4.f90, except it uses the de-symmetrized representation of the dipole surfaces (defined in "dip_xy4_diff_mu")

subroutine dipole_xy4_alpha_ADF_noreexp(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: i, j, k, imu, nmodes, iterm
  type(adf_realq) :: r1, r2, r3, r4, amat_len, tmat(4,3), amat(3,3), ainv(3,3), dip(3), y1, y2, y3, y4, y5, y6, y7, y8, y9
  real(ark) :: r1eq, coefs_amat(3,3,adf_nterms), coefs_ainv(3,3,adf_nterms), coefs_sym(681)

  nmodes = molec%nmodes

  if (.not.present(cart)) then
    write(out, '(/a)') 'dipole_xy4_alpha_ADF_noreexp error: Cartesian-coordinate derivatives are not present'
    stop
  endif

  r1eq = func%params(1,1)

  if (trim(molec%coord_transform)=='XY4_RSYMALPHA') then

    r1 = internal(1)
    r2 = internal(2)
    r3 = internal(3)
    r4 = internal(4)

    y1 = (r1-r1eq) *exp(-1.0_ark*(r1-r1eq)**2)
    y2 = (r2-r1eq) *exp(-1.0_ark*(r2-r1eq)**2)
    y3 = (r3-r1eq) *exp(-1.0_ark*(r3-r1eq)**2)
    y4 = (r4-r1eq) *exp(-1.0_ark*(r4-r1eq)**2)

    y5 = internal(5)
    y6 = internal(6)
    y7 = internal(7)
    y8 = internal(8)
    y9 = internal(9)

  else

    write(out, '(/a,a,a)') 'dipole_xy4_alpha_ADF_noreexp error: coordinate type ="',trim(molec%coord_transform), '" is not supported'
    stop

  endif

  ! compute symmetry-adapted dipole moment projections DIP_{S1,S2,S3}

  imu = 1
  coefs_sym = 0
  do i=2, func%nparams(imu)
    k = func%iparams(1,i,imu)
    coefs_sym(k) = func%params(i,imu)
  enddo

  do imu=1, 3
    dip(imu) = dip_xy4_diff_mu_init_dipole(imu, nmodes, (/y1,y2,y3,y4,y5,y6,y7,y8,y9/), coefs_sym)
  enddo

  ! compute AMAT: Cartesian to symmetry-adapted unit vector transformation

  tmat(1,1:3) = (cart(2,1:3)-cart(1,1:3))/r1 ! X--Y1
  tmat(2,1:3) = (cart(3,1:3)-cart(1,1:3))/r2 ! X--Y2
  tmat(3,1:3) = (cart(4,1:3)-cart(1,1:3))/r3 ! X--Y3
  tmat(4,1:3) = (cart(5,1:3)-cart(1,1:3))/r4 ! X--Y4

  amat(1,:) = (tmat(1,:)-tmat(2,:)+tmat(3,:)-tmat(4,:))*0.5_ark
  amat(2,:) = (tmat(1,:)-tmat(2,:)-tmat(3,:)+tmat(4,:))*0.5_ark
  amat(3,:) = (tmat(1,:)+tmat(2,:)-tmat(3,:)-tmat(4,:))*0.5_ark

  do i=1, 3
    amat_len = sqrt( amat(i,1)*amat(i,1) + amat(i,2)*amat(i,2) + amat(i,3)*amat(i,3) )
    amat(i,:) = amat(i,:)/amat_len
  enddo

  ! find AMAT^{-1} using algebraic approach

  do iterm=1, adf_nterms
    do i=1, 3
      do j=1, 3
        coefs_amat(j,i,iterm) = amat(j,i)%d(iterm)
      enddo
    enddo
  enddo

  call deriv_invmat(nmodes, adf_nterms, adf_terms, 3, coefs_amat, coefs_ainv)

  ! switch to ADF approach

  do i=1, 3
    do j=1, 3
      call adf_set_var(ainv(j,i), coefs_ainv(j,i,1:adf_nterms))
    enddo
  enddo

  ! compute DIP_{X,Y,Z} = AMAT^{-1} * DIP_{S1,S2,S3}

  do imu=1, 3
    f(imu) = dot_product(ainv(imu,1:3), dip)
  enddo

end subroutine dipole_xy4_alpha_ADF_noreexp


!################################################################################


function dip_xy4_diff_mu_init_dipole(imu, nmodes, local, f) result(res)
  use adf
  implicit none

  integer(ik), intent(in) :: imu, nmodes
  type(adf_realq), intent(in) :: local(nmodes)
  real(ark), intent(in) :: f(:)
  type(adf_realq) :: res

  integer(ik), parameter :: maxnterms = 5000
  integer(ik) :: iterm, nterms(3), terms(9,maxnterms,3), imode
  real(ark) :: coefs(maxnterms,3)
  type(adf_realq) :: prod

  if (nmodes/=9) then
    write(out, '(/a,1x,i3,1x,a)') 'dip_xy4_diff_mu_init_dipole error: number of internal coordinates =', nmodes, '!= 9'
    stop
  endif

  if (imu>3 .or. imu<1) then
    write(out, '(/a,1x,i3,1x,a)') 'dip_xy4_diff_mu_init_dipole error: dipole moment component =', imu, '!= 1, 2, or 3'
    stop
  endif

  ! compute all symmetry-equivalent expansion coefficients

  nterms = 0
  terms = 0
  coefs = 0.0

#include 'pot_xy4_2_part3.f90'

  if (any(nterms==0)) then
    write(out, '(/a)') 'dip_xy4_diff_mu_init_dipole error: part 3 of the file pot_xy4_3.f90 was not included at compilation, recompile module "hamiltonian" using -fpp or -cpp flag'
    stop
  endif

  if (any(nterms>maxnterms)) then
    write(out, '(/a,1x,i6,1x,a,1x,i6)') 'dip_xy4_diff_mu_init_dipole error: number of terms =', maxval(nterms), 'exceeds max number =', maxnterms
    stop
  endif

  ! compute dipoles

  res = 0.0_ark
  do iterm=1, nterms(imu)
    prod = 1.0_ark
    do imode=1, nmodes
      if (terms(imode,iterm,imu)==1) then
        prod = prod * local(imode)
      elseif (terms(imode,iterm,imu)>1) then
        prod = prod * local(imode)**terms(imode,iterm,imu)
      endif
    enddo
    res = res + coefs(iterm,imu) * prod
  enddo

end function dip_xy4_diff_mu_init_dipole
