subroutine shift_cm(molec, cart, cart_cm)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: cart(molec%natoms,3)
  real(ark), intent(out)          :: cart_cm(molec%natoms,3)

  integer(ik) :: i, iatom, natoms
  real(ark) :: cm(3), fm

  natoms = molec%natoms

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do i=1, 3
    cm(i) = sum( cart(1:natoms,i) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cart_cm(iatom,1:3) = cart(iatom,1:3) - cm(1:3)
  enddo

end subroutine shift_cm


!################################################################################


subroutine rotate_I(molec, cart, cart_I, rotmat)

  type(HM_molec_type), intent(in)  :: molec
  real(ark), intent(in)            :: cart(molec%natoms,3)
  real(ark), intent(out)           :: cart_I(molec%natoms,3)
  real(ark), intent(out), optional :: rotmat(3,3)

  integer(ik) :: i

  call shift_cm(molec, cart, cart_I)
  if (present(rotmat)) then
    rotmat = 0.0
    forall(i=1:3) rotmat(i,i) = 1.0_ark
  endif

end subroutine rotate_I


!################################################################################


! Rotates Cartesian coordinates of atoms to principal axes system (PAS).
! Rotation matrix is obtained by diagonalizing the inertia tensor.

subroutine rotate_pas_diag(molec, cart, cart_pas, rotmat)

  type(HM_molec_type), intent(in)  :: molec
  real(ark), intent(in)            :: cart(molec%natoms,3)
  real(ark), intent(out)           :: cart_pas(molec%natoms,3)
  real(ark), intent(out), optional :: rotmat(3,3)

  integer(ik) :: i, j, lwork, info, natoms
  real(ark) :: imat(3,3)
  double precision :: imatd(3,3), diagd(3), work(128*3)

  natoms = molec%natoms

  call shift_cm(molec, cart, cart_pas)

  ! compute inertia matrix
  do i=1, 3
    do j=i+1, 3
      imat(i,j) = -sum( cart_pas(1:natoms,i) * cart_pas(1:natoms,j) * molec%masses(1:natoms) )
      imat(j,i) = imat(i,j)
    enddo
  enddo
  imat(1,1) = sum( (cart_pas(1:natoms,2)**2 + cart_pas(1:natoms,3)**2) * molec%masses(1:natoms) )
  imat(2,2) = sum( (cart_pas(1:natoms,1)**2 + cart_pas(1:natoms,3)**2) * molec%masses(1:natoms) )
  imat(3,3) = sum( (cart_pas(1:natoms,1)**2 + cart_pas(1:natoms,2)**2) * molec%masses(1:natoms) )

  ! diagonalize inertia matrix
  imatd = dble(imat)
  lwork = size(work)
  call dsyev('V', 'L', 3, imatd, 3, diagd, work, lwork, info)
  if (info/=0) then
    write(out, '(/a,1x,i6)') 'rotate_pas_diag error: diagonalization of inertia matrix failed, lapack info =', info
    stop
  endif
  imat = real(imatd, kind=ark)

  ! rotate coordinates
  cart_pas = matmul(cart_pas, imat)
  if (present(rotmat)) rotmat = transpose(imat)

end subroutine rotate_pas_diag


!################################################################################


! Rotates Cartesian coordinates of atoms to principal axes system (PAS).
!
! Coordinate rotation:
!   cartesian^{pas} = exp(-kappa) * cartesian^{ref},
!   where kappa = -kappa^T is skew-symmetric rotation matrix,
!   and "cartesian" are Cartesian coordinates of atoms in reference axes system (ref) and PAS.
!
! PAS equations:
!   exp(-kappa)*u*exp(kappa) = diagonal-matrix,
!   where u_{ix,iy} = \sum_{iatom} mass_{iatom} * cartesian_{iatom,ix}^{ref} * cartesian_{iatom,iy}^{ref},
!   ix and iy = x, y, or z, and iatom = 1..number-of-atoms.
!
! Iterative solution:
!   (lambda-kappa)*u*(lambda^T+kappa) = diagonal matrix,
!   i.e. -kappa*u*lambda^T  + lambda*u*kappa = kappa*u*kappa - lambda*u*lambda^T (for off-diagonal elements),
!   where lambda = exp(-kappa) + kappa is subsequently calculated from kappa on previous iteration.

subroutine rotate_pas(molec, cart, cart_pas, rotmat)

  type(HM_molec_type), intent(in)  :: molec
  real(ark), intent(in)            :: cart(molec%natoms,3)
  real(ark), intent(out)           :: cart_pas(molec%natoms,3)
  real(ark), intent(out), optional :: rotmat(3,3)

  integer(ik) :: iter, i, j, info, ielem, lwork, rank, natoms
  real(ark) :: kappa(3,3), rot_mat(3,3), kappa_p(3,3), exp_kappa(3,3), lambda(3,3), umat(3,3), mat(3,3), vmat(3,3), kappa0(3,3), tmat(3,3), vec(3), normtol
  double precision :: matd(3,3), vecd(3,1), sv(3), work(128*3), rcond

  natoms = molec%natoms

  ! initial guess for rotation matrices
  rot_mat = 0.0
  forall(i=1:3) rot_mat(i,i) = 1.0_ark
  kappa = 0.0

  ! shift coordinates to the centre-of-mass
  call shift_cm(molec, cart, cart_pas)

  do iter=1, molec%pas_maxniter

    ! compute exp(-kappa) using scaled Taylor series expansion
    normtol = 0.01_ark
    call matexp_taylor(3, kappa, molec%matexp_deg2, normtol, molec%matexp_maxntaylor, molec%matexp_accur, exp_kappa)
    exp_kappa = transpose(exp_kappa)

    ! cumulative rotation matrix
    rot_mat = matmul(exp_kappa, rot_mat)

    ! lambda = exp(-kappa) + kappa
    lambda = exp_kappa + kappa

    ! rotate coordinates
    cart_pas = matmul(cart_pas, transpose(exp_kappa))

    ! compute system of linear equations for three elements of kappa:
    ! -kappa*u*lambda^T  + lambda*u*kappa = kappa*u*kappa - lambda*u*lambda^T

    ! u matrix
    umat = 0.0
    do i=1, 3
      do j=1, 3
        umat(i,j) = sum(molec%masses(:) * cart_pas(:,i) * cart_pas(:,j))
      enddo
    enddo

    ! right-hand side vector: kappa*u*kappa - lambda*u*lambda^T,
    vmat = matmul(kappa, matmul(umat, kappa)) - matmul(lambda, matmul(umat, transpose(lambda)))
    vec = (/vmat(1,2), vmat(1,3), vmat(2,3)/)

    ! matrix: -kappa*u*lambda^T  + lambda*u*kappa
    mat = 0.0
    ielem = 0
    do i=1, 3
      do j=i+1, 3
        kappa0 = 0.0
        kappa0(i,j) =  1.0_ark
        kappa0(j,i) = -1.0_ark
        ielem = ielem + 1
        tmat = matmul(lambda, matmul(umat, kappa0)) - matmul(kappa0, matmul(umat, transpose(lambda)))
        mat(1,ielem) = tmat(1,2)
        mat(2,ielem) = tmat(1,3)
        mat(3,ielem) = tmat(2,3)
      enddo
    enddo

    ! solve system of linear equations for three elements of kappa (double precision lapack)
    matd = dble(mat)
    vecd(1:3,1) = dble(vec)
    lwork = size(work)
    rcond = -1.0d-16
    call dgelss(3, 3, 1, matd, 3, vecd, 3, sv, rcond, rank, work, lwork, info)
    if (info/=0) then
      write(out, '(/a,1x,i6)') 'rotate_pas error: failed to solve system of linear equations, lapack info =', info
      stop
    endif
    vec = real(vecd(:,1), kind=ark)

    ! solve system of linear equations for three elements of kappa (quadrupole precision)
    !call linsolve3_cramer(mat, vec)

    ! new kappa matrix
    kappa(1,1:3) = (/ 0.0_ark,   vec(1),  vec(2)/)
    kappa(2,1:3) = (/ -vec(1),  0.0_ark,  vec(3)/)
    kappa(3,1:3) = (/ -vec(2),  -vec(3), 0.0_ark/)

    ! check convergence
    if (all(abs(kappa)<=molec%pas_accur)) exit

    if (iter==molec%pas_maxniter) then
      write(out, '(/a,1x,i3,1x,a,1x,es16.8)') 'rotate_pas error: failed to converge PAS equations in', molec%pas_maxniter, 'iterations to accuracy =', molec%pas_accur
      write(out, '(a,1x,es16.8,1x,a,3(es16.8),1x,a)') 'Delta kappa (>', molec%pas_accur, ') = [', vec, ']'
      write(out, '(a)') 'most likely the reference axes system is ill-defined or very far from the principal axes system'
      stop
    endif

  enddo ! iter

  if (present(rotmat)) rotmat = rot_mat

end subroutine rotate_pas


!################################################################################


! Rotates Cartesian coordinates of atoms to Eckart axes system.
!
! Coordinate rotation:
!   cartesian^{eckart} = exp(-kappa) * cartesian^{ref},
!   where kappa = -kappa^T is skew-symmetric rotation matrix,
!   and "cartesian" are Cartesian coordinates of atoms in reference axes system (ref) and Eckart axes system.
!
! Eckart equations:
!   u*exp(kappa) - exp(-kappa)*u^T = 0,
!   where u_{ix,iy} = \sum_{iatom} mass_{iatom} * cartesian_{iatom,ix}^{eq,pas} * cartesian_{iatom,iy}^{ref},
!   and cartesian^{eq,pas} are equilibrium Cartesian coordinates in principal axes system.
!
! Iterative solution:
!   u*(lambda^T+kappa) - (lambda-kappa)*u^T = 0,
!   where lambda = exp(-kappa) + kappa is subsequently calculated from kappa on previous iteration.

subroutine rotate_eckart(molec, cart, cart_eq_pas, cart_eckart, rotmat)

  type(HM_molec_type), intent(in)  :: molec
  real(ark), intent(in)            :: cart(molec%natoms,3)
  real(ark), intent(in)            :: cart_eq_pas(molec%natoms,3)
  real(ark), intent(out)           :: cart_eckart(molec%natoms,3)
  real(ark), intent(out), optional :: rotmat(3,3)

  integer(ik) :: iter, i, j, info, ielem, lwork, rank, iatom, natoms
  real(ark) :: kappa(3,3), rot_mat(3,3), kappa_p(3,3), exp_kappa(3,3), lambda(3,3), umat(3,3), mat(3,3), vmat(3,3), kappa0(3,3), tmat(3,3), vec(3), normtol, euler(2,3)
  double precision :: matd(3,3), vecd(3,1), sv(3), work(128*3), rcond

  natoms = molec%natoms

  write(out, '(/a)') 'rotate_eckart: solve Eckart equations'
  write(out, '(1x,a)') 'Cartesian coordinates:'
  do iatom=1, natoms
    write(out, '(1x,i3,3(1x,f30.16))') iatom, cart(iatom,1:3)
  enddo
  write(out, '(1x,a)') 'reference Cartesian coordinates:'
  do iatom=1, natoms
    write(out, '(1x,i3,3(1x,f30.16))') iatom, cart_eq_pas(iatom,1:3)
  enddo

  ! initial guess for Eckart system - PAS coordinates
  !call rotate_pas(natoms, masses, cart, maxit, accur, cart_eckart, rot_mat)

  ! initial guess for Eckart coordinates - centre-of-mass coordinates
  call shift_cm(molec, cart, cart_eckart)

  ! initial guess for rotation matrices
  rot_mat = 0.0
  forall(i=1:3) rot_mat(i,i) = 1.0_ark
  kappa = 0.0

  write(out, '(1x,a)') 'start iterations'

  do iter=1, molec%eckart_maxniter

    ! compute exp(-kappa) using scaled Taylor series expansion
    normtol = 0.01_ark
    call matexp_taylor(3, kappa, molec%matexp_deg2, normtol, molec%matexp_maxntaylor, molec%matexp_accur, exp_kappa)
    exp_kappa = transpose(exp_kappa)

    ! cumulative rotation matrix
    rot_mat = matmul(exp_kappa, rot_mat)

    ! lambda = exp(-kappa) + kappa
    lambda = exp_kappa + kappa

    ! rotate coordinates
    cart_eckart = matmul(cart_eckart, transpose(exp_kappa))

    ! compute system of linear equations for three elements of kappa:
    ! u*kappa + kappa*u^T = lambda*u^T - u*lambda^T,
    ! where lambda = exp(-kappa) + kappa

    ! u matrix
    umat = 0.0
    do i=1, 3
      do j=1, 3
        umat(i,j) = sum(molec%masses(:) * cart_eq_pas(:,i) * cart_eckart(:,j))
      enddo
    enddo

    ! right hand side vector: lambda*u^T - u*lambda^T
    vmat = matmul(lambda, transpose(umat)) - matmul(umat, transpose(lambda))
    vec = (/vmat(1,2), vmat(1,3), vmat(2,3)/)

    ! matrix: u*kappa + kappa*u^T
    mat = 0.0
    ielem = 0
    do i=1, 3
      do j=i+1, 3
        kappa0 = 0.0
        kappa0(i,j) =  1.0_ark
        kappa0(j,i) = -1.0_ark
        ielem = ielem + 1
        tmat = matmul(umat, kappa0) + matmul(kappa0, transpose(umat))
        mat(1,ielem) = tmat(1,2)
        mat(2,ielem) = tmat(1,3)
        mat(3,ielem) = tmat(2,3)
      enddo
    enddo

    ! solve system of linear equations for three elements of kappa (double precision lapack)
    matd = dble(mat)
    vecd(1:3,1) = dble(vec)
    lwork = size(work)
    rcond = -1.0d-16
    call dgelss(3, 3, 1, matd, 3, vecd, 3, sv, rcond, rank, work, lwork, info)
    if (info/=0) then
      write(out, '(/a,1x,i6)') 'rotate_eckart error: failed to solve system of linear equations, lapack info =', info
      stop
    endif
    vec = real(vecd(:,1), kind=ark)

    ! new kappa matrix
    kappa(1,1:3) = (/ 0.0_ark,   vec(1),  vec(2)/)
    kappa(2,1:3) = (/ -vec(1),  0.0_ark,  vec(3)/)
    kappa(3,1:3) = (/ -vec(2),  -vec(3), 0.0_ark/)

    ! check convergence
    if (all(abs(kappa)<=molec%eckart_accur)) exit

    if (iter==molec%eckart_maxniter) then
      write(out, '(/a,1x,i3,1x,a,1x,es16.8)') 'rotate_eckart error: failed to converge Eckart equations in', molec%eckart_maxniter, 'iterations to accuracy =', molec%eckart_accur
      write(out, '(a,1x,es16.8,1x,a,3(es16.8),1x,a)') 'Delta kappa (>', molec%eckart_accur, ') = [', vec, ']'
      write(out, '(a)') 'most likely the reference axes system is ill-defined or very far from the Eckart axes system'
      stop
    endif

    write(out, '(1x,i3,3(1x,f35.30))') iter, vec(1:3)

  enddo ! iter

  if (present(rotmat)) rotmat = rot_mat

  write(out, '(1x,a)') 'Cartesian coordinates in Eckart frame:'
  do iatom=1, natoms
    write(out, '(1x,i3,3(1x,f30.16))') iatom, cart_eckart(iatom,1:3)
  enddo

  write(out, '(1x,a,1x,a,1x,a,1x,a)') 'Euler angles (deg)', 'theta', 'psi', 'phi'
  call euler_from_rotmat(rot_mat, euler)
  do i=1, 2
    write(out, '(3(1x,es16.8))') euler(i,1:3)*180.0_ark/real(pi,ark)
  enddo

  write(out, '(a)') 'rotate_eckart: done'

end subroutine rotate_eckart


!################################################################################


subroutine rotate_eckart_lin(molec, cart, cart_eq_pas, cart_eckart, rotmat)

  type(HM_molec_type), intent(in)  :: molec
  real(ark), intent(in)            :: cart(molec%natoms,3)
  real(ark), intent(in)            :: cart_eq_pas(molec%natoms,3)
  real(ark), intent(out)           :: cart_eckart(molec%natoms,3)
  real(ark), intent(out), optional :: rotmat(3,3)

  integer(ik) :: i, j, info, ielem, lwork, rank, iatom, natoms
  real(ark) :: kappa(3,3), rot_mat(3,3), kappa_p(3,3), exp_kappa(3,3), lambda(3,3), umat(3,3), mat(3,3), vmat(3,3), kappa0(3,3), tmat(3,3), vec(3), normtol, euler(2,3)
  double precision :: matd(3,3), vecd(3,1), sv(3), work(128*3), rcond

  natoms = molec%natoms

  write(out, '(/a)') 'rotate_eckart_lin: solve quasi-Eckart equations'
  write(out, '(1x,a)') 'Cartesian coordinates:'
  do iatom=1, natoms
    write(out, '(1x,i3,3(1x,f30.16))') iatom, cart(iatom,1:3)
  enddo
  write(out, '(1x,a)') 'reference Cartesian coordinates:'
  do iatom=1, natoms
    write(out, '(1x,i3,3(1x,f30.16))') iatom, cart_eq_pas(iatom,1:3)
  enddo

  ! initial guess for Eckart system - PAS coordinates
  !call rotate_pas(natoms, masses, cart, maxit, accur, cart_eckart, rot_mat)

  ! initial guess for Eckart coordinates - centre-of-mass coordinates
  call shift_cm(molec, cart, cart_eckart)
  rot_mat = 0.0
  forall(i=1:3) rot_mat(i,i) = 1.0_ark

  kappa = 0.0
  lambda = 0.0
  forall(i=1:3) lambda(i,i) = 1.0_ark

  ! compute system of linear equations for three elements of kappa: u*kappa + kappa*u^T = lambda*u^T - u*lambda^T

  ! u matrix
  umat = 0.0
  do i=1, 3
    do j=1, 3
      umat(i,j) = sum(molec%masses(:) * cart_eq_pas(:,i) * cart_eckart(:,j))
    enddo
  enddo

  ! right hand side vector: lambda*u^T - u*lambda^T
  vmat = matmul(lambda, transpose(umat)) - matmul(umat, transpose(lambda))
  vec = (/vmat(1,2), vmat(1,3), vmat(2,3)/)

  ! matrix: u*kappa + kappa*u^T
  mat = 0.0
  ielem = 0
  do i=1, 3
    do j=i+1, 3
      kappa0 = 0.0
      kappa0(i,j) =  1.0_ark
      kappa0(j,i) = -1.0_ark
      ielem = ielem + 1
      tmat = matmul(umat, kappa0) + matmul(kappa0, transpose(umat))
      mat(1,ielem) = tmat(1,2)
      mat(2,ielem) = tmat(1,3)
      mat(3,ielem) = tmat(2,3)
    enddo
  enddo

  ! solve system of linear equations for three elements of kappa (double precision lapack)
  matd = dble(mat)
  vecd(1:3,1) = dble(vec)
  lwork = size(work)
  rcond = -1.0d-16
  call dgelss(3, 3, 1, matd, 3, vecd, 3, sv, rcond, rank, work, lwork, info)
  if (info/=0) then
    write(out, '(/a,1x,i6)') 'rotate_eckart_lin error: failed to solve system of linear equations, lapack info =', info
    stop
  endif
  vec = real(vecd(:,1), kind=ark)

  ! new kappa matrix
  kappa(1,1:3) = (/ 0.0_ark,   vec(1),  vec(2)/)
  kappa(2,1:3) = (/ -vec(1),  0.0_ark,  vec(3)/)
  kappa(3,1:3) = (/ -vec(2),  -vec(3), 0.0_ark/)

  ! compute exp(-kappa) using scaled Taylor series expansion
  normtol = 0.01_ark
  call matexp_taylor(3, kappa, molec%matexp_deg2, normtol, molec%matexp_maxntaylor, molec%matexp_accur, exp_kappa)
  exp_kappa = transpose(exp_kappa)

  ! cumulative rotation matrix
  rot_mat = matmul(exp_kappa, rot_mat)

  ! rotate coordinates
  cart_eckart = matmul(cart_eckart, transpose(exp_kappa))


  if (present(rotmat)) rotmat = rot_mat

  write(out, '(1x,a)') 'Cartesian coordinates in Eckart frame:'
  do iatom=1, natoms
    write(out, '(1x,i3,3(1x,f30.16))') iatom, cart_eckart(iatom,1:3)
  enddo

  write(out, '(1x,a,1x,a,1x,a,1x,a)') 'Euler angles (deg)', 'theta', 'psi', 'phi'
  call euler_from_rotmat(rot_mat, euler)
  do i=1, 2
    write(out, '(3(1x,es16.8))') euler(i,1:3)*180.0_ark/real(pi,ark)
  enddo

  write(out, '(a)') 'rotate_eckart_lin: done'

end subroutine rotate_eckart_lin


!################################################################################


! Determines orthogonal matrix of rotation between Cartesian coordinates "cart0" and "cart".

subroutine rotate_to_xyz(molec, cart0, cart, rotmat)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: cart0(molec%natoms,3)
  real(ark), intent(in)           :: cart(molec%natoms,3)
  real(ark), intent(out)          :: rotmat(3,3)

  integer(ik) :: iter, i, j, info, ielem, lwork, rank, natoms, irow, iatom
  real(ark) :: kappa(3,3), exp_kappa(3,3), lambda(3,3), tmat(molec%natoms,3), vec(3), normtol, cart0_cm(molec%natoms,3), cart_cm(molec%natoms,3), cart2(molec%natoms,3)
  double precision :: matd(molec%natoms*3,3), vecd(molec%natoms*3,1), sv(3), work(128*molec%natoms*3), rcond

  natoms = molec%natoms

  ! initial guess for rotation matrices
  rotmat = 0.0
  forall(i=1:3) rotmat(i,i) = 1.0_ark
  kappa = 0.0

  ! shift coordinates to the centre-of-mass
  call shift_cm(molec, cart0, cart0_cm)
  call shift_cm(molec, cart, cart_cm)

  cart2 = cart0_cm

  do iter=1, molec%pas_maxniter

    ! compute exp(-kappa) using scaled Taylor series expansion
    normtol = 0.001_ark
    call matexp_taylor(3, kappa, molec%matexp_deg2, normtol, molec%matexp_maxntaylor, molec%matexp_accur, exp_kappa)
    exp_kappa = transpose(exp_kappa)
    !call matexp_taylor(3, -kappa, molec%matexp_deg2, normtol, molec%matexp_maxntaylor, molec%matexp_accur, exp_kappa)

    ! cumulative rotation matrix
    rotmat = matmul(exp_kappa, rotmat)

    ! lambda = exp(-kappa) + kappa
    lambda = exp_kappa + kappa

    ! rotate coordinates
    cart2 = matmul(cart2, transpose(exp_kappa))

    ! compute system of linear equations for three elements of kappa

    matd = 0.0
    irow = 0
    do iatom=1, natoms
      irow = irow + 1
      matd(irow,1:3)   = dble( (/  cart2(iatom,2),  cart2(iatom,3),        0.0_ark /) )
      matd(irow+1,1:3) = dble( (/ -cart2(iatom,1),         0.0_ark,  cart2(iatom,3) /) )
      matd(irow+2,1:3) = dble( (/         0.0_ark, -cart2(iatom,1), -cart2(iatom,2) /) )
      irow = irow + 2
    enddo

    tmat = matmul(cart2, transpose(lambda)) - cart_cm
    vecd = 0.0
    irow = 0
    do iatom=1, natoms
      irow = irow + 1
      vecd(irow,1)   = dble( tmat(iatom,1) )
      vecd(irow+1,1) = dble( tmat(iatom,2) )
      vecd(irow+2,1) = dble( tmat(iatom,3) )
      irow = irow + 2
    enddo

    ! solve system of linear equations for three elements of kappa (double precision lapack)
    lwork = size(work)
    rcond = -1.0d-16
    call dgelss(natoms*3, 3, 1, matd, natoms*3, vecd, natoms*3, sv, rcond, rank, work, lwork, info)
    if (info/=0) then
      write(out, '(/a,1x,i6)') 'rotate_to_xyz error: failed to solve system of linear equations, lapack info =', info
      stop
    endif
    vec(1:3) = real(vecd(1:3,1), ark)

    ! new kappa matrix
    kappa(1,1:3) = (/ 0.0_ark,   vec(1),  vec(2)/)
    kappa(2,1:3) = (/ -vec(1),  0.0_ark,  vec(3)/)
    kappa(3,1:3) = (/ -vec(2),  -vec(3), 0.0_ark/)

    ! check convergence
    if (all(abs(kappa)<=molec%pas_accur)) exit

    if (iter==molec%pas_maxniter) then
      write(out, '(/a,1x,i3,1x,a,1x,es16.8)') 'rotate_to_xyz error: failed to converge in', molec%pas_maxniter, 'iterations to accuracy =', molec%pas_accur
      write(out, '(a,1x,es16.8,1x,a,3(es16.8),1x,a)') 'Delta kappa (>', molec%pas_accur, ') = [', vec, ']'
      stop
    endif

  enddo ! iter

end subroutine rotate_to_xyz


!################################################################################


! Simple algorithm to compute exponential of matrix using saled Taylor series expansion.

subroutine matexp_taylor(dimen, mat, maxdeg2, normtol, maxorder, accur, matexp)

  integer(ik), intent(in) :: dimen
  real(ark), intent(in)   :: mat(dimen,dimen)
  integer(ik), intent(in) :: maxdeg2, maxorder
  real(ark), intent(in)   :: normtol, accur
  real(ark), intent(out)  :: matexp(dimen,dimen)

  integer(ik) :: deg2, iorder, ideg, i
  real(ark) :: norm, mat_s(dimen,dimen), mat_p(dimen,dimen)

  ! scale matrix
  norm = maxval(abs(mat))
  deg2 = 0
  do
    if (deg2>=maxdeg2) then
      write(out, '(/a,1x,i3,1x,a)') 'matexp_taylor error: max power of 2 =', maxdeg2, 'is exceeded when computing scaling factor'
      exit
    endif
    if (abs(norm)<=normtol) exit
    norm = norm / real(2**deg2,kind=ark)
    deg2 = deg2 + 1
  end do
  mat_s = mat/real(2**deg2,kind=ark)

  ! compute exponential of scaled matrix
  matexp = 0.0; forall(i=1:dimen) matexp(i,i) = 1.0_ark
  mat_p = matexp
  do iorder=1, maxorder
    mat_p = matmul(mat_s, mat_p)/real(iorder,kind=ark)
    matexp = matexp + mat_p
    if (all(abs(mat_p)<=accur)) exit
    if (iorder==maxorder) then
      write(out, '(/a,1x,es16.8)') 'matexp_taylor error: failed to converge exponential of scaled matrix to accuracy =', accur
      write(out, '(a,i3,1x,a)') 'scaled_matrix **', maxorder, '= ['
      write(out, *) mat_p
      write(out, '(a)') ']'
      stop
    endif
  enddo

  ! compute exponential of original matrix
  do ideg=1, deg2
    matexp = matmul(matexp, matexp)
  enddo

end subroutine matexp_taylor


!################################################################################


! Solves 3x3 linear equations system using Cramer's approach

subroutine linsolve3_cramer(m, v)

  real(ark), intent(in) :: m(3,3)
  real(ark), intent(inout) :: v(3)

  real(ark) :: det0, det1, det2, det3

  det0 = m(1,3)*m(2,2)*m(3,1)-m(1,2)*m(2,3)*m(3,1)-m(1,3)*m(2,1)*m(3,2)+m(1,1)*m(2,3)*m(3,2)+m(1,2)*m(2,1)*m(3,3)-m(1,1)*m(2,2)*m(3,3)

  if (abs(det0)<1.0d-12) then

    write(out, '(/a,1x,es16.8)') 'linsolve3_cramer error: dependent system of linear equations, determinant =', det0
    write(out, '(a)') 'linear equations matrix:'
    write(out, *) m(:,:)
    stop

  else

    det1 = -(m(2,3)*m(3,2)*v(1))+m(2,2)*m(3,3)*v(1)+m(1,3)*m(3,2)*v(2)-m(1,2)*m(3,3)*v(2)-m(1,3)*m(2,2)*v(3)+m(1,2)*m(2,3)*v(3)
    det2 = m(2,3)*m(3,1)*v(1)-m(2,1)*m(3,3)*v(1)-m(1,3)*m(3,1)*v(2)+m(1,1)*m(3,3)*v(2)+m(1,3)*m(2,1)*v(3)-m(1,1)*m(2,3)*v(3)
    det3 = -(m(2,2)*m(3,1)*v(1))+m(2,1)*m(3,2)*v(1)+m(1,2)*m(3,1)*v(2)-m(1,1)*m(3,2)*v(2)-m(1,2)*m(2,1)*v(3)+m(1,1)*m(2,2)*v(3)

    v = (/det1/det0, det2/det0, det3/det0/)

  endif

end subroutine linsolve3_cramer


!################################################################################


! Computes Euler angles (2 sets) for given rotation matrix

subroutine euler_from_rotmat(rotmat, euler)

  real(ark), intent(in) :: rotmat(3,3)
  real(ark), intent(out) :: euler(2,3)

  real(ark) :: theta1, theta2, psi1, psi2, phi1, phi2

  if (abs(abs(rotmat(3,1))-1.0_ark)>epsilon(1.0_ark)) then

    theta1 = -asin(rotmat(3,1))
    theta2 = real(pi,ark) - theta1

    psi1 = atan2(rotmat(3,2)/cos(theta1),rotmat(3,3)/cos(theta1))
    psi2 = atan2(rotmat(3,2)/cos(theta2),rotmat(3,3)/cos(theta2))

    phi1 = atan2(rotmat(2,1)/cos(theta1),rotmat(1,1)/cos(theta1))
    phi2 = atan2(rotmat(2,1)/cos(theta2),rotmat(1,1)/cos(theta2))

  else

    phi1 = 0.0
    phi2 = 0.0

    if (abs(rotmat(3,1)+1.0_ark)<=epsilon(1.0_ark)) then

      theta1 = real(pi,ark)*0.5_ark
      theta2 = real(pi,ark)*0.5_ark

      psi1 = phi1 + atan2(rotmat(1,2),rotmat(1,3))
      psi2 = phi2 + atan2(rotmat(1,2),rotmat(1,3))

    else

      theta1 = -real(pi,ark)*0.5_ark
      theta2 = -real(pi,ark)*0.5_ark

      psi1 = -phi1 + atan2(-rotmat(1,2),-rotmat(1,3))
      psi2 = -phi2 + atan2(-rotmat(1,2),-rotmat(1,3))

    endif

  endif

  euler(1:2,1) = (/theta1, theta2/)
  euler(1:2,2) = (/psi1, psi2/)
  euler(1:2,3) = (/phi1, phi2/)

end subroutine euler_from_rotmat
