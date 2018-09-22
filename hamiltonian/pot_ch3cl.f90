! C3v-symmetry-adapted potential energy function for CH3Cl-type molecule
! for details see CH3Cl/CH3D paper of Alec Owens et al.

subroutine poten_ch3cl_sym_ADF(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: ioper, i, j, icoord, irank
  type(adf_realq) :: r0, r1, r2, r3, beta1, beta2, beta3, stau1, stau2, y(9), xi(9,6), term, prod
  real(ark) :: r0eq, r1eq, betaeq, amorse0, amorse1

  irank = 1

  r0eq = func%params(1,irank)
  r1eq = func%params(2,irank)
  betaeq = func%params(3,irank)*real(pi,ark)/180.0_ark
  amorse0 = func%params(4,irank)
  amorse1 = func%params(5,irank)

  ! define internal coordinates

  if (trim(molec%coord_transform)=='ZXY3_BETA_SYMTAU') then

    r0 = internal(1)
    r1 = internal(2)
    r2 = internal(3)
    r3 = internal(4)
    beta1 = internal(5)
    beta2 = internal(6)
    beta3 = internal(7)
    stau1 = internal(8)
    stau2 = internal(9)

    y(1) = 1.0_ark-exp(-amorse0*(r0-r0eq))
    y(2) = 1.0_ark-exp(-amorse1*(r1-r1eq))
    y(3) = 1.0_ark-exp(-amorse1*(r2-r1eq))
    y(4) = 1.0_ark-exp(-amorse1*(r3-r1eq))

    y(5) = beta1 - betaeq
    y(6) = beta2 - betaeq
    y(7) = beta3 - betaeq

    y(8) = stau1
    y(9) = stau2

  else

    write(out, '(/a,a,a)') 'poten_ch3cl_sym_ADF error: coordinate type = "', trim(molec%coord_transform), '" is not supported'
    stop

  endif

  ! define all symmetric permutations of internal coordinates

  do ioper=1, 6
    call ch3cl_symmetry_transformation_local_ADF(molec%coord_transform, ioper, y, xi(:,ioper))
  enddo

  ! compute symmetry-adapted potential

  f = 0.0_ark

  do i=6, func%nparams(irank)
    term = 0.0_ark
    do ioper=1, 6
      prod = 1.0_ark
      do icoord=1, 9
        if (func%iparams(icoord,i,irank)==1) then
          prod = prod * xi(icoord,ioper)
        elseif (func%iparams(icoord,i,irank)>1) then
          prod = prod * xi(icoord,ioper)**func%iparams(icoord,i,irank)
        endif
      enddo
      term = term + prod
    enddo
    f(irank) = f(irank) + term * func%params(i,irank)
  enddo

end subroutine poten_ch3cl_sym_ADF


!################################################################################


! Transforms internal coordinates of CH3CL molecule according to "ioper" symmetry operation of C3v point group.

subroutine ch3cl_symmetry_transformation_local_ADF(coord_transform, ioper, src, dst)
  use adf
  implicit none

  character(*), intent(in) :: coord_transform
  integer, intent(in) :: ioper
  type(adf_realq), intent(in)  :: src(1:9)
  type(adf_realq), intent(out) :: dst(1:9)

  real(ark) :: repres(12,2,2),a,b,e,o

  dst(:) = src(:)

  if (trim(coord_transform)=='ZXY3_BETA_SYMTAU') then

    select case(ioper)

    case (1) ! identity
      dst = src

    case (3) ! (123)
      dst(2) = src(3)
      dst(3) = src(4)
      dst(4) = src(2)
      dst(5) = src(6)
      dst(6) = src(7)
      dst(7) = src(5)

    case (2) ! (321)
      dst(2) = src(4)
      dst(3) = src(2)
      dst(4) = src(3)
      dst(5) = src(7)
      dst(6) = src(5)
      dst(7) = src(6)

    case (6) ! (12)
      dst(2) = src(3)
      dst(3) = src(2)
      dst(4) = src(4)
      dst(5) = src(6)
      dst(6) = src(5)
      dst(7) = src(7)

    case (5) ! (13)
      dst(2) = src(4)
      dst(3) = src(3)
      dst(4) = src(2)
      dst(5) = src(7)
      dst(6) = src(6)
      dst(7) = src(5)

    case (4) ! (23)
      dst(2) = src(2)
      dst(3) = src(4)
      dst(4) = src(3)
      dst(5) = src(5)
      dst(6) = src(7)
      dst(7) = src(6)

    case default

      write(out, '(/a,1x,i3)') 'ch3cl_symmetry_transformation_local_ADF error: unknown symmetry operation number =', ioper
      stop

    end select

    a = 0.5_ark
    b = 0.5_ark*sqrt(3.0_ark)
    e = 1.0_ark
    o = 0.0_ark

    repres ( 1,:,:)= reshape((/ e, o,  &!
                                o, e/),(/2,2/))

    repres ( 3,:,:)= reshape((/-a,-b,  &!
                                b,-a/),(/2,2/))

    repres ( 2,:,:)= reshape((/-a, b,  &!
                               -b,-a/),(/2,2/))

    repres ( 4,:,:)= reshape((/ e, o,  &!
                                o,-e/),(/2,2/))

    repres ( 6,:,:)= reshape((/-a, b,  &!
                                b, a/),(/2,2/))

    repres ( 5,:,:)= reshape((/-a,-b,  &!
                               -b, a/),(/2,2/))

    dst(8) = repres(ioper,1,1)*src(8)+repres(ioper,1,2)*src(9)
    dst(9) = repres(ioper,2,1)*src(8)+repres(ioper,2,2)*src(9)

  else

    write(out, '(/a,a,a)') 'ch3cl_symmetry_transformation_local_ADF error: symmetry transformation rules for coordinate type = "', trim(coord_transform), '" are not defined'
    stop

  endif

end subroutine ch3cl_symmetry_transformation_local_ADF


!################################################################################


! C3v-symmetry-adapted dipole moment function for CH3Cl-type molecule
! for details see CH3Cl/CH3D paper of Alec Owens et al.

subroutine dipole_ch3cl_sym_ADF(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: ioper, i, j, icoord, imu, iterm, k_ind(9), nmodes
  type(adf_realq) :: r0, r1, r2, r3, beta1, beta2, beta3, stau1, stau2, y(9), xi(9,6), amat_len, tmat(4,3), amat(3,3), ainv(3,3), dip(3), mu(3), nu(3), term_dip(3)
  real(ark) :: r0eq, r1eq, betaeq, coefs_amat(3,3,adf_nterms), coefs_ainv(3,3,adf_nterms)

  nmodes = molec%nmodes

  if (.not.present(cart)) then
    write(out, '(/a)') 'dipole_ch3cl_sym_ADF error: Cartesian-coordinate derivatives are not present'
    stop
  endif

  r0eq = func%params(1,1)
  r1eq = func%params(2,1)
  betaeq = func%params(3,1)*real(pi,ark)/180.0_ark

  ! define internal coordinates

  if (trim(molec%coord_transform)=='ZXY3_BETA_SYMTAU') then

    r0 = internal(1)
    r1 = internal(2)
    r2 = internal(3)
    r3 = internal(4)
    beta1 = internal(5)
    beta2 = internal(6)
    beta3 = internal(7)
    stau1 = internal(8)
    stau2 = internal(9)

    y(1) = r0 - r0eq
    y(2) = r1 - r1eq
    y(3) = r2 - r1eq
    y(4) = r3 - r1eq

    y(5) = beta1 - betaeq
    y(6) = beta2 - betaeq
    y(7) = beta3 - betaeq

    y(8) = stau1
    y(9) = stau2

  else

    write(out, '(/a,a,a)') 'dipole_ch3cl_sym_ADF: coordinate type = "', trim(molec%coord_transform), '" is not supported'
    stop

  endif

  do ioper=1, 6
    call ch3cl_symmetry_transformation_local_ADF(molec%coord_transform, ioper, y, xi(:,ioper))
  enddo

  ! compute symmetry-adapted dipole moment functions

  dip(1:3) = 0.0_ark

  do imu=1, 3

    mu(1:3) = 0.0_ark

    do iterm=4, func%nparams(imu)

      k_ind(1:9) = func%iparams(1:9,iterm,imu)
      term_dip = 0.0_ark

      do ioper=1, 6

        mu(imu) = 1.0_ark
        do icoord=1, 9
          if (k_ind(icoord)==1) then
            mu(imu) = mu(imu) * xi(icoord,ioper)
          elseif (k_ind(icoord)>1) then
            mu(imu) = mu(imu) * xi(icoord,ioper)**k_ind(icoord)
          endif
        enddo

        call ch3cl_symmetry_transformation_dipole_ADF(ioper,mu,nu)

        term_dip(:) = term_dip(:) + nu(:)

      enddo

      dip(:) = dip(:) + term_dip(:) * func%params(iterm,imu)

    enddo ! iterm
  enddo ! imu

  ! form transformation from Cartesian directions to symmetry-adapted bond-projections
  ! AMAT * DIP_{X,Y,Z} = DIP_{S1,S2,S3}

  tmat(1,1:3) = (cart(2,1:3)-cart(1,1:3))/r0 ! C--Cl
  tmat(2,1:3) = (cart(3,1:3)-cart(1,1:3))/r1 ! C--H1
  tmat(3,1:3) = (cart(4,1:3)-cart(1,1:3))/r2 ! C--H2
  tmat(4,1:3) = (cart(5,1:3)-cart(1,1:3))/r3 ! C--H3

  amat(1,:) = (tmat(2,:)*2.0_ark-tmat(3,:)-tmat(4,:))/sqrt(6.0_ark)
  amat(2,:) = (                  tmat(3,:)-tmat(4,:))/sqrt(2.0_ark)
  amat(3,:) = tmat(1,:)

  do i=1, 3
    !amat_len = sqrt(sum(amat(i,:)**2))
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

end subroutine dipole_ch3cl_sym_ADF


!################################################################################


subroutine ch3cl_symmetry_transformation_dipole_ADF(ioper,src,dst)
    use adf
    implicit none
    !
    integer, intent(in) :: ioper  ! group operation
    type(adf_realq), intent(in)  :: src(1:3)
    type(adf_realq), intent(out) :: dst(1:3)
    !
    real(ark) :: repres(6,3,3),a,b,e,o
    !
    a = 0.5_ark
    b = 0.5_ark*sqrt(3.0_ark)
    e = 1.0_ark
    o = 0.0_ark
    !
    repres = 0
    !
    repres(:,3,3) = e !no operations on mu_z
    !
    repres(1,1,1) = e !identity
    repres(1,2,2) = e
    !
    repres(2,1,1) = -a !(321)
    repres(2,1,2) = -b
    repres(2,2,1) =  b
    repres(2,2,2) = -a
    !
    repres(3,1,1) = -a !(123)
    repres(3,1,2) =  b
    repres(3,2,1) = -b
    repres(3,2,2) = -a
    !
    repres(4,1,1) =  e !(23)
    repres(4,1,2) =  o
    repres(4,2,1) =  o
    repres(4,2,2) = -e
    !
    repres(5,1,1) = -a !(13)
    repres(5,1,2) = -b
    repres(5,2,1) = -b
    repres(5,2,2) =  a
    !
    repres(6,1,1) = -a !(12)
    repres(6,1,2) =  b
    repres(6,2,1) =  b
    repres(6,2,2) =  a
    !
    if (ioper<0.or.ioper>6) then
      write(out, '(/a,1x,i3,1x,a)') 'ch3cl_symmetry_transformation_diple_ADF error: operation =', ioper, 'is unknown'
      stop
    endif
    !
    dst = matmul(transpose(repres(ioper,:,:)),src)
    !
end subroutine ch3cl_symmetry_transformation_dipole_ADF
