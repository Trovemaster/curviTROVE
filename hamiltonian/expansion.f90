! Computes derivatives of Cartesian coordinates with respect to internal coordinates using ADF.

subroutine deriv_cart_ADF( molec, x0_npoints, x0, order, nmodes_active, modes_active, modes_order, nterms, terms, coefs )
  use adf
  implicit none

  type(HM_molec_type), intent(in)       :: molec
  integer(ik), intent(in)               :: x0_npoints
  real(ark), intent(in)                 :: x0(molec%nmodes,x0_npoints)
  integer(ik), intent(in)               :: order
  integer(ik), intent(in)               :: nmodes_active
  integer(ik), intent(in)               :: modes_active(nmodes_active)
  integer(ik), intent(in)               :: modes_order(nmodes_active)
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)
  real(ark), allocatable, intent(out)   :: coefs(:,:,:,:)

  integer(ik) :: natoms, nmodes, imode, jmode, ndeg(molec%nmodes), iorder, iterm, iatom, ix, info, deg(order+1,molec%nmodes), ipoint
  type(adf_realq) :: x0_ADF(molec%nmodes), cartesian_ADF(molec%natoms,3)

  natoms = molec%natoms
  nmodes = molec%nmodes

  if (associated(internal_to_cartesian_ADF_ptr).eqv..false.) then
    write(out, '(/a)') 'deriv_cart_ADF error: function pointer "internal_to_cartesian_ADF_ptr" is not associated'
    stop
  endif

  ! generate derivative terms for expansion 0..order

  ndeg = 1
  deg = 0
  do imode=1, nmodes_active
    jmode = modes_active(imode)
    ndeg(jmode) = modes_order(imode)+1
    deg(1:ndeg(jmode), jmode) = (/(iorder, iorder=0, modes_order(imode))/)
  enddo
  call expansion_terms(nmodes, ndeg(1:nmodes), deg(1:maxval(ndeg),1:nmodes), 0, order, nterms, terms) ! "terms" is (re-)allocated inside "expansion_terms"

  write(out, '(3x,a,1x,i6)') 'nterms =', nterms

  ! allocate array to keep derivative values

  if (allocated(coefs)) deallocate(coefs)
  allocate(coefs(natoms,3,nterms,x0_npoints), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'deriv_cart_ADF error: failed to allocate coefs(natoms,3,nterms,x0_npoints)', 'natoms, nterms, x0_npoints =', natoms, nterms, x0_npoints
    stop
  endif
  coefs = 0.0

  ! compute derivatives

  ipoint = 1
  call adf_init(nmodes, x0_ADF, x0(:,ipoint), nterms, terms)
  call internal_to_cartesian_ADF_ptr(molec, x0_ADF, cartesian_ADF)
  do iterm=1, nterms
    do ix=1, 3
      do iatom=1, natoms
        coefs(iatom,ix,iterm,ipoint) = cartesian_ADF(iatom,ix)%d(iterm)
      enddo
    enddo
  enddo

  do ipoint=2, x0_npoints
    call adf_reset_x0(nmodes, x0_ADF, x0(:,ipoint))
    call internal_to_cartesian_ADF_ptr(molec, x0_ADF, cartesian_ADF)
    do iterm=1, nterms
      do ix=1, 3
        do iatom=1, natoms
          coefs(iatom,ix,iterm,ipoint) = cartesian_ADF(iatom,ix)%d(iterm)
        enddo
      enddo
    enddo
  enddo

  call adf_finalize

end subroutine deriv_cart_ADF


!################################################################################


! Computes derivatives of external function with respect to internal coordinates using ADF.

subroutine deriv_func_ADF( molec, func, x0_npoints, x0, order, nmodes_active, modes_active, modes_order, nterms_cart, terms_cart, coefs_cart, nterms, terms, coefs )
  use adf
  implicit none

  type(HM_molec_type), intent(in)       :: molec
  type(HM_func_type), intent(in)        :: func
  integer(ik), intent(in)               :: x0_npoints
  real(ark), intent(in)                 :: x0(molec%nmodes,x0_npoints)
  integer(ik), intent(in)               :: order
  integer(ik), intent(in)               :: nmodes_active
  integer(ik), intent(in)               :: modes_active(nmodes_active)
  integer(ik), intent(in)               :: modes_order(nmodes_active)
  integer(ik), intent(in)               :: nterms_cart
  integer(ik), intent(in)               :: terms_cart(molec%nmodes,nterms_cart)
  real(ark), intent(in)                 :: coefs_cart(molec%natoms,3,nterms_cart,x0_npoints)
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)
  real(ark), allocatable, intent(out)   :: coefs(:,:,:)

  integer(ik) :: natoms, nmodes, imode, jmode, ndeg(molec%nmodes), deg(order+1,molec%nmodes), iorder, iterm, iatom, ix, info, ipoint, irank, rank, ipos
  real(ark) :: fx0(molec%nmodes,x0_npoints), cosx0, amorse
  real(ark), allocatable :: coefs_cart_(:,:,:,:)
  type(adf_realq) :: x0_ADF(molec%nmodes), fx0_ADF(molec%nmodes), f(func%rank), cart(molec%natoms,3)

  natoms = molec%natoms
  nmodes = molec%nmodes
  rank = func%rank

  ! check if external function is assigned to pointer func_ADF_ptr

  if (associated(func_ADF_ptr).eqv..false.) then
    write(out, '(/a)') 'deriv_func_ADF error: function pointer "func_ADF_ptr" is not associated'
    stop
  endif

  ! generate derivative terms for expansion 0..order

  ndeg = 1
  deg = 0
  do imode=1, nmodes_active
    jmode = modes_active(imode)
    ndeg(jmode) = modes_order(imode)+1
    deg(1:ndeg(jmode), jmode) = (/(iorder, iorder=0, modes_order(imode))/)
  enddo
  call expansion_terms(nmodes, ndeg(1:nmodes), deg(1:maxval(ndeg),1:nmodes), 0, order, nterms, terms) ! "terms" is (re-)allocated inside "expansion_terms"

  ! allocate array to keep derivative values

  if (allocated(coefs)) deallocate(coefs)
  allocate(coefs(nterms,x0_npoints,rank), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'deriv_func_ADF error: failed to allocate coefs(nterms,x0_npoints,rank)', 'nterms, x0_npoints, rank =', nterms, x0_npoints, rank
    stop
  endif
  coefs = 0.0

  ! compute reference expansion functions, i.e. alpha_ref -> cos(alpha_ref)

  do imode=1, nmodes
    select case(trim(func%coord_type(imode)))
    case default
      write(out, '(/a,1x,a,1x,a,1x,i3)') 'deriv_func_ADF error: unknown type of expansion function = "', trim(func%coord_type(imode)), '" for mode #', imode
      stop
    case('LINEAR')
      fx0(imode,1:x0_npoints) = x0(imode,1:x0_npoints)
    case('MORSE')
      fx0(imode,:) = 0.0_ark
    case('COSRHO')
      fx0(imode,:) = 0.0_ark
    case('SINRHO')
      fx0(imode,:) = 0.0_ark
    end select
  enddo

  ! establish the correspondence between indices of Cartesian-coordinate derivatives terms in "terms" and "molec%terms_cart"

  if (func%rotate) then

    allocate(coefs_cart_(natoms,3,nterms,x0_npoints), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') 'deriv_func_ADF error: failed to allocate coefs_cart_(natoms,3,nterms,x0_npoints)', 'natoms, nterms, x0_npoints =', natoms, nterms, x0_npoints
      stop
    endif
    coefs_cart_ = 0.0

    do iterm=1, nterms
      ipos = index_iarr1(terms(1:nmodes,iterm), terms_cart(1:nmodes,1:nterms_cart))
      if (ipos<=0) then
        write(out, '(/a,<nmodes>(1x,i3),a)') 'deriv_func_ADF error: Cartesian-coordinate derivative term = (', terms(1:nmodes,iterm), ') is not initialized'
        stop
      else
        coefs_cart_(:,:,iterm,:) = coefs_cart(:,:,ipos,:)
      endif
    enddo

  endif

  ! compute derivatives for all reference points

  do ipoint=1, x0_npoints

    ! init expansion points and ADF structures

    if (ipoint==1) then
      call adf_init(nmodes, fx0_ADF, fx0(:,ipoint), nterms, terms)
    else
      call adf_reset_x0(nmodes, fx0_ADF, fx0(:,ipoint))
    endif

    ! obtain local coordinates from expansion functions, i.e. cos(alpha) -> alpha

    do imode=1, nmodes
      select case(trim(func%coord_type(imode)))
      case default
        write(out, '(/a,1x,a,1x,a,1x,i3)') 'deriv_func_ADF error: unknown type of expansion function = "', trim(func%coord_type(imode)), '" for mode #', imode
        stop
      case('LINEAR')
        x0_ADF(imode) = fx0_ADF(imode)
      case('MORSE')
        amorse = func%coord_params(imode,1)
        x0_ADF(imode) = -log(1.0_ark-fx0_ADF(imode))/amorse + x0(imode,ipoint)
      case('COSRHO')
        x0_ADF(imode) = acos(fx0_ADF(imode) + sin(x0(imode,ipoint)))
      case('SINRHO')
        x0_ADF(imode) = asin(fx0_ADF(imode) + sin(x0(imode,ipoint)))
      end select
    enddo

    ! init derivatives of Cartesian coordinates

    if (func%rotate) then
      do ix=1, 3
        do iatom=1, natoms
          call adf_set_var(cart(iatom,ix), coefs_cart_(iatom,ix,1:nterms,ipoint))
        enddo
      enddo
    endif

    ! compute external function

    if (func%rotate) then
      call func_ADF_ptr(molec, func, x0_ADF, f(1:rank), cart)
    else
      call func_ADF_ptr(molec, func, x0_ADF, f(1:rank))
    endif

    ! copy derivatives

    do irank=1, rank
      coefs(1:nterms,ipoint,irank) = f(irank)%d(1:nterms)
    enddo

  enddo ! ipoint

  call adf_finalize

  if (func%rotate) then
    deallocate(coefs_cart_)
  endif

end subroutine deriv_func_ADF


!################################################################################


! Computes derivatives of inverse matrix A^{-1} providing that derivatives of A are known.

subroutine deriv_invmat(nmodes, nterms, terms, dimen, coefs_amat, coefs_ainv)

  integer(ik), intent(in) :: nmodes
  integer(ik), intent(in) :: nterms
  integer(ik), intent(in) :: terms(nmodes,nterms)
  integer(ik), intent(in) :: dimen
  real(ark), intent(in)   :: coefs_amat(dimen,dimen,nterms)
  real(ark), intent(out)  :: coefs_ainv(dimen,dimen,nterms)

  integer(ik) :: iterm, ielem, iterm1, iterm2, info, ipos, imode, i
  integer(ik), allocatable :: ind_chr(:,:,:), nterms_chr(:)
  real(ark) :: mat(dimen,dimen), ainv(dimen,dimen)
  logical :: coefs_ainv_init(nterms)

  ! generate chain rule expressions for all required derivatives of a product of two variables

  allocate(nterms_chr(nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'deriv_invmat error: failed to allocate nterms_chr(nterms)', 'nterms =', nterms
    stop
  endif
  call chain_rule_product_ind(nmodes, nterms, terms, 2, nterms_chr, ind_chr)

  ! compute inverse of A

  ! locate zero-order derivative of A
  ipos = index_iarr1((/(0, imode=1,nmodes)/), terms)
  if (ipos<=0) then
    write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_invmat error: derivative of A = (', (/(0, imode=1,nmodes)/), ') is not initialized'
    stop
  endif

  !call matrix_inverse_ark(dimen, coefs_amat(:,:,ipos), ainv)
  call matrix_inverse2_ark(dimen, coefs_amat(:,:,ipos), ainv)

  ! start recursive calculations of derivatives A^{-1}

  coefs_ainv = 0.0
  coefs_ainv_init = .false.

  do iterm=1, nterms

    mat = 0.0
    if (all(terms(:,iterm)==0)) forall(i=1:dimen) mat(i,i) = 1.0_ark

    ! loop over elements of the chain-rule expression for terms(:,iterm)-derivative of product A * A^{-1}

    do ielem=1, nterms_chr(iterm)

      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)

      if (iterm2==iterm) then

        if (any(terms(:,iterm1)/=0)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_invmat error: corrupted chain-rule expression for derivative = (', terms(:,iterm), ')'
          stop
        endif

      else

        if (.not.coefs_ainv_init(iterm2)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_invmat error: derivative of A^{-1} = (', terms(:,iterm2), ') is not initialized'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of A^{-1} = (', terms(:,iterm), ')'
          write(out, '(a)') 'IMPORTANT: derivative indices must be sorted ascendingly wrt the total derivative order'
          stop
        endif

        mat = mat - matmul(coefs_amat(:,:,iterm1), coefs_ainv(:,:,iterm2))

      endif
    enddo

    coefs_ainv(:,:,iterm) = matmul(ainv, mat)
    coefs_ainv_init(iterm) = .true.

  enddo ! iterm


  if (allocated(nterms_chr)) deallocate(nterms_chr)
  if (allocated(ind_chr)) deallocate(ind_chr)

end subroutine deriv_invmat


!################################################################################


! Computes product of two polynomials.

subroutine poly_prod(nmodes, nterms, terms, coefs1, coefs2, coefs)

  integer(ik), intent(in)         :: nmodes, nterms
  integer(ik), intent(in)         :: terms(nmodes, nterms)
  real(ark), intent(in)           :: coefs1(nterms)
  real(ark), intent(in)           :: coefs2(nterms)
  real(ark), intent(out)          :: coefs(nterms)

  integer(ik) :: iterm, ielem, iterm1, iterm2, info
  integer(ik), allocatable, save :: nterms_chr(:), ind_chr(:,:,:)

  allocate(nterms_chr(nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'poly_prod error: failed to allocate nterms_chr(nterms)', 'nterms =', nterms
    stop
  endif
  call chain_rule_product_ind(nmodes, nterms, terms, 2, nterms_chr, ind_chr)

  do iterm=1, nterms
    coefs(iterm) = 0.0
    do ielem=1, nterms_chr(iterm)
      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)
      coefs(iterm) = coefs(iterm) + coefs1(iterm1) * coefs2(iterm2)
    enddo
  enddo

  deallocate(nterms_chr)
  deallocate(ind_chr)

end subroutine poly_prod


!################################################################################


! Computes derivatives of Cartesian coordinates of atoms wrt internal coordinates in "rotated" coordinate frame.

subroutine deriv_cart_rotated(molec, nterms, terms, nterms_chr, ind_chr, coefs_cart, coefs_rotmat, coefs)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: nterms
  integer(ik), intent(in)         :: terms(molec%nmodes,nterms)
  integer(ik), intent(in)         :: nterms_chr(nterms)
  integer(ik), intent(in)         :: ind_chr(2,maxval(nterms_chr),nterms)
  real(ark), intent(in)           :: coefs_cart(molec%natoms,3,nterms)
  real(ark), intent(in)           :: coefs_rotmat(3,3,nterms)
  real(ark), intent(out)          :: coefs(molec%natoms,3,nterms)

  integer(ik) :: iterm, ielem, iterm1, iterm2, info, natoms, nmodes

  natoms = molec%natoms
  nmodes = molec%nmodes

  do iterm=1, nterms
    ! sum up all elements in the chain rule expression for terms(:,iterm)-derivative of rotation_matrix*cartesian product
    coefs(:,:,iterm) = 0.0
    do ielem=1, nterms_chr(iterm)
      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)
      coefs(:,:,iterm) = coefs(:,:,iterm) + matmul(coefs_cart(:,:,iterm2), transpose(coefs_rotmat(:,:,iterm1)))
    enddo
  enddo

end subroutine deriv_cart_rotated


!################################################################################


! Computes derivatives of Eckart-frame rotation matrix wrt internal coordinates.

subroutine deriv_rotmat_eckart(molec, cart0_ref, nterms, terms, nterms_chr, ind_chr, coefs_cart, coefs_expkappa)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: cart0_ref(molec%natoms,3)
  integer(ik), intent(in)         :: nterms
  integer(ik), intent(in)         :: terms(molec%nmodes,nterms)
  integer(ik), intent(in)         :: nterms_chr(nterms)
  integer(ik), intent(in)         :: ind_chr(2,maxval(nterms_chr),nterms)
  real(ark), intent(in)           :: coefs_cart(molec%natoms,3,nterms)
  real(ark), intent(out)          :: coefs_expkappa(3,3,nterms)

  integer(ik) :: iterm, i, j, ielem, iterm1, iterm2, ideg, iter, natoms, nmodes, imode, ipos, info
  real(ark) :: coefs_umat(3,3,nterms), coefs_kappa(3,3,molec%matexp_maxntaylor,nterms), coefs_kappa_scaled(3,3,molec%matexp_maxntaylor,nterms), &!
               coefs_expkappa_scaled(3,3,molec%matexp_deg2,nterms), tmat(3,3), umat(3,3), uinv(3,3), vec0(3), kappa0(3,3), vecp(3,molec%eckart_maxniter), &!
               lambda(3,3), vec(3), mat(3,3)
  logical :: coefs_kappa_init(molec%matexp_maxntaylor,nterms), coefs_expkappa_init(nterms), coefs_umat_init(nterms), coefs_kappa_scaled_init(molec%matexp_maxntaylor,nterms), &!
             coefs_expkappa_scaled_init(molec%matexp_deg2,nterms)

  natoms = molec%natoms
  nmodes = molec%nmodes

  ! precompute derivatives of u-matrix

  coefs_umat_init(:) = .false.
  do iterm=1, nterms
    do i=1, 3
      do j=1, 3
        coefs_umat(i,j,iterm) = sum( molec%masses(1:natoms) * cart0_ref(1:natoms,i) * coefs_cart(1:natoms,j,iterm) )
        coefs_umat_init(iterm) = .true.
      enddo
    enddo
  enddo

  ! precompute linear system matrix for kappa-derivatives (same for all derivatives)

  ! locate zero-order derivative
  ipos = index_iarr1((/(0, imode=1,nmodes)/), terms)
  if (ipos<=0) then
    write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: derivative of u-matrix = (', (/(0, imode=1,nmodes)/), ') is not initialized'
    stop
  endif

  umat = coefs_umat(:,:,ipos)
  mat = 0.0
  ielem = 0
  do i=1, 3
    do j=i+1, 3
      ielem = ielem + 1
      kappa0 = 0.0
      kappa0(i,j) =  1.0_ark
      kappa0(j,i) = -1.0_ark
      tmat = matmul(umat, kappa0) + matmul(kappa0, transpose(umat))
      mat(1,ielem) = tmat(1,2)
      mat(2,ielem) = tmat(1,3)
      mat(3,ielem) = tmat(2,3)
    enddo
  enddo

  call matrix_inverse_ark(3, mat, uinv)

  ! start recursive calculations of derivatives of exp(-kappa)

  coefs_kappa(:,:,:,:) = 0.0
  coefs_kappa_init(:,:) = .false.
  coefs_expkappa(:,:,:) = 0.0
  coefs_expkappa_init(:) = .false.

  coefs_kappa_scaled(:,:,:,:) = 0.0
  coefs_kappa_scaled_init(:,:) = .false.
  coefs_expkappa_scaled(:,:,:,:) = 0.0
  coefs_expkappa_scaled_init(:,:) = .false.

  do iterm=1, nterms

    ! compute right-hand side vector for current derivative of exp(-kappa)

    tmat = 0.0

    do ielem=1, nterms_chr(iterm)

      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)

      if (.not.coefs_umat_init(iterm1)) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: derivative of u-matrix = (', terms(:,iterm1), ') is not initialized'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
        stop
      endif

      if (iterm2==iterm) then

        if (any(terms(:,iterm1)/=0)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: corrupted chain rule expression for derivative = (', terms(:,iterm), ')'
          stop
        endif

      else

        if (.not.coefs_expkappa_init(iterm2)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: derivative of exp(-kappa) = (', terms(:,iterm2), ') is not initialized'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
          write(out, '(a)') 'IMPORTANT: derivative indices must be sorted ascendingly wrt the total derivative order'
          stop
        endif

        tmat = tmat + matmul(coefs_umat(:,:,iterm1), transpose(coefs_expkappa(:,:,iterm2)))

      endif

    enddo ! ielem

    ! compute part of the right hand side vector that does not depend on current derivative of exp(-kappa)
    tmat = tmat - transpose(tmat)
    vec0 = (/tmat(1,2), tmat(1,3), tmat(2,3)/)

    ! solve system of non-linear equations for current derivative of exp(-kappa)

    coefs_kappa(:,:,:,iterm) = 0.0
    coefs_kappa_init(1,iterm) = .true.
    coefs_kappa_init(2:,iterm) = .false.
    coefs_expkappa_init(iterm) = .false.

    vecp = 0.0

    do iter=1, molec%eckart_maxniter

      ! compute current derivative of exp(-kappa) using scaled Taylor series expansion

      coefs_kappa_init(2:,iterm) = .false.
      coefs_kappa(:,:,2:,iterm) = 0.0

      coefs_expkappa_scaled_init(:,iterm) = .false.
      coefs_expkappa_scaled(:,:,:,iterm) = 0.0
      coefs_kappa_scaled_init(:,iterm) = .false.
      coefs_kappa_scaled(:,:,:,iterm) = 0.0

      coefs_kappa_scaled(:,:,1,iterm) = coefs_kappa(:,:,1,iterm)/real(2**molec%matexp_deg2,kind=ark)
      coefs_kappa_scaled_init(1,iterm) = .true.

      ! compute current derivative for all matrix powers (2..molec%matexp_maxntaylor) of scaled kappa (scaling step)

      do ideg=2, molec%matexp_maxntaylor
        tmat = 0.0
        do ielem=1, nterms_chr(iterm)
          iterm1 = ind_chr(1,ielem,iterm)
          iterm2 = ind_chr(2,ielem,iterm)
          if (.not.coefs_kappa_scaled_init(1,iterm1)) then
            write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: derivative of scaled kappa = (', terms(:,iterm1), ') is not initialized'
            write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
            write(out, '(a)') 'IMPORTANT: derivative indices must be sorted ascendingly wrt the total derivative order'
            stop
          endif
          if (.not.coefs_kappa_scaled_init(ideg-1,iterm2)) then
            write(out, '(/a,1x,i3,1x,a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: derivative of scaled kappa**', ideg-1, ' = (', terms(:,iterm2), ') is not initialized'
            write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
            write(out, '(a)') 'IMPORTANT: derivative indices must be sorted ascendingly wrt the total derivative order'
            stop
          endif
          tmat = tmat + matmul(coefs_kappa_scaled(:,:,1,iterm1), coefs_kappa_scaled(:,:,ideg-1,iterm2))
        enddo
        coefs_kappa_scaled(:,:,ideg,iterm) = tmat / real( ideg, kind=ark )
        coefs_kappa_scaled_init(ideg,iterm) = .true.
      enddo

      ! compute current derivative of exponential of scaled kappa (expansion step)

      tmat = 0.0
      if (all(terms(:,iterm)==0)) forall(i=1:3) tmat(i,i) = 1.0_ark
      tmat = tmat - coefs_kappa_scaled(:,:,1,iterm)
      do ideg=2, molec%matexp_maxntaylor
        if (mod(ideg,2)==0) then
          tmat = tmat + coefs_kappa_scaled(:,:,ideg,iterm)
        else
          tmat = tmat - coefs_kappa_scaled(:,:,ideg,iterm)
        endif
        if (all(abs(coefs_kappa_scaled(:,:,ideg,iterm))<molec%matexp_accur)) exit
        if (ideg==molec%matexp_maxntaylor) then
          write(out, '(/a,1x,es16.8)') 'deriv_rotmat_eckart error: failed to Taylor converge derivative of exp(-kappa_scaled) to accuracy =', molec%matexp_accur
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
          write(out, '(a,1x,i3,1x,a)') 'derivative of kappa_scaled**', molec%matexp_maxntaylor, '='
          write(out, *) coefs_kappa_scaled(:,:,molec%matexp_maxntaylor,iterm)
          stop
        endif
      enddo
      coefs_expkappa_scaled(:,:,1,iterm) = tmat
      coefs_expkappa_scaled_init(1,iterm) = .true.

      ! compute current derivative of exponential of kappa (squaring step)

      do ideg=1, molec%matexp_deg2
        tmat = 0.0
        do ielem=1, nterms_chr(iterm)
          iterm1 = ind_chr(1,ielem,iterm)
          iterm2 = ind_chr(2,ielem,iterm)
          if (.not.coefs_expkappa_scaled_init(ideg,iterm1)) then
            write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: derivative of exp(-kappa_scaled) = (', terms(:,iterm1), ') is not initialized'
            write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
            write(out, '(a)') 'IMPORTANT: derivative indices must be sorted ascendingly wrt the total derivative order'
            stop
          endif
          if (.not.coefs_expkappa_scaled_init(ideg,iterm2)) then
            write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_rotmat_eckart error: derivative of exp(-kappa_scaled) = (', terms(:,iterm2), ') is not initialized'
            write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
            write(out, '(a)') 'IMPORTANT: derivative indices must be sorted ascendingly wrt the total derivative order'
            stop
          endif
          tmat = tmat + matmul(coefs_expkappa_scaled(:,:,ideg,iterm1), coefs_expkappa_scaled(:,:,ideg,iterm2))
        enddo
        if (ideg==molec%matexp_deg2) then
          coefs_expkappa(:,:,iterm) = tmat
        else
          coefs_expkappa_scaled(:,:,ideg+1,iterm) = tmat
          coefs_expkappa_scaled_init(ideg+1,iterm) = .true.
        endif
      enddo

      ! current derivative of lambda = exp(-kappa) + kappa

      lambda = coefs_expkappa(:,:,iterm) + coefs_kappa(:,:,1,iterm)

      ! compute right-hand side vector

      tmat = matmul(lambda, transpose(umat)) - matmul(umat, transpose(lambda))
      vec = (/ tmat(1,2)-vec0(1), tmat(1,3)-vec0(2), tmat(2,3)-vec0(3) /)

      ! find new kappa

      vec = matmul(uinv, vec)

      ! check convergence

      !if (iter>1) then
      !  vec = vecp(:,iter-1) + eckart_damp*(vec-vecp(:,iter-1))
      !endif

      if (iter>3) then
        if (all(abs(vec-vecp(:,iter-1))<=molec%eckart_accur)) then
          vecp(:,iter) = vec
          exit
        endif
      endif

      vecp(:,iter) = vec

      if (iter==molec%eckart_maxniter) then
        write(out, '(/a,1x,es16.8,1x,a,1x,i3,1x,a)') 'deriv_rotmat_eckart error: failed to converge Eckart equations to accuracy =', molec%eckart_accur, 'in', molec%eckart_maxniter, 'iterations'
        write(out, '(a)') 'most likely the reference axes system is ill-defined or very far from the Eckart axes system'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
        write(out, '(1x,a,17x,a,15x,a,15x,a)') 'iter', 'kappa(1,2)', 'kappa(1,3)', 'kappa(2,3)'
        do i=1, iter
          write(out, '(1x,i3,3x,3(1x,es40.32))') i, vecp(1:3,i)
        enddo
        stop
      endif

      ! update kappa

      coefs_kappa(1,1:3,1,iterm) = (/ 0.0_ark,   vec(1),  vec(2)/)
      coefs_kappa(2,1:3,1,iterm) = (/ -vec(1),  0.0_ark,  vec(3)/)
      coefs_kappa(3,1:3,1,iterm) = (/ -vec(2),  -vec(3), 0.0_ark/)

    enddo ! iter

    coefs_expkappa_init(iterm) = .true.

  enddo ! iterm

end subroutine deriv_rotmat_eckart


!################################################################################


! Computes derivatives of s-vectors from derivatives of Cartesian coordinates.

subroutine deriv_svec(molec, nterms_cart, terms_cart, coefs_cart, nterms, terms, nterms_chr, ind_chr, coefs_svec)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: nterms_cart
  integer(ik), intent(in)         :: terms_cart(molec%nmodes,nterms_cart)
  real(ark), intent(in)           :: coefs_cart(molec%natoms,3,nterms_cart)
  integer(ik), intent(in)         :: nterms
  integer(ik), intent(in)         :: terms(molec%nmodes,nterms)
  integer(ik), intent(in)         :: nterms_chr(nterms)
  integer(ik), intent(in)         :: ind_chr(2,maxval(nterms_chr),nterms)
  real(ark), intent(out)          :: coefs_svec(molec%natoms*3,molec%natoms*3,nterms)

  integer(ik) :: ncoords, iterm, iterm1, iterm2, ielem, i, j, t(molec%nmodes), ipos, imode, ix, iatom, icoord, info, natoms, nmodes
  real(ark) :: tmat0(molec%natoms*3,molec%natoms*3), tmat0_inv(molec%natoms*3,molec%natoms*3), bvec(molec%natoms*3,molec%natoms*3), coefs_tvec(molec%natoms*3,molec%natoms*3), asym_tens(3,3,3)
  logical :: coefs_svec_init(nterms), tmat0_init, ifterm0
  real(ark) :: etamat(molec%natoms*3,molec%natoms*3), cs(molec%natoms*3,molec%natoms*3)

  natoms = molec%natoms
  nmodes = molec%nmodes

  ! Levi-Civita tensor

  asym_tens = 0.0
  asym_tens(1,2,3) = 1.0_ark
  asym_tens(1,3,2) =-1.0_ark
  asym_tens(2,1,3) =-1.0_ark
  asym_tens(2,3,1) = 1.0_ark
  asym_tens(3,1,2) = 1.0_ark
  asym_tens(3,2,1) =-1.0_ark

  ! start recursive calculations of s-vectors derivatives

  coefs_svec(:,:,:) = 0.0
  coefs_svec_init(:) = .false.

  tmat0 = 0.0
  tmat0_init = .false.

  do iterm=1, nterms

    ! set up system of linear equations for terms(:,iterm)-derivative of s-vectors

    bvec = 0.0
    ifterm0 = .false.

    ! for 0-order derivative right-hand side is delta symbol
    if (all(terms(1:nmodes,iterm)==0)) then
      forall(i=1:natoms*3) bvec(i,i) = 1.0_ark
      ifterm0 = .true.
    endif

    ! loop over elements of the chain-rule expression for terms(:,iterm)-derivative of t*s product

    do ielem=1, nterms_chr(iterm)

      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)

      coefs_tvec(:,:) = 0.0

      ! derivative of vibrational t-vector: t_{iatom,ix,imode}^{vib} = d cartesian_{iatom,ix} / d internal_{imode}

      do imode=1, nmodes
        t = terms(1:nmodes,iterm1)
        t(imode) = t(imode) + 1
        ipos = index_iarr1(t, terms_cart)
        if (ipos==0) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: Cartesian-coordinate derivative = (', t, ') is not found'
          write(out, '(a,1x,i3,1x,a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of vibrational t-vector for', imode, 'mode = (', terms(:,iterm1), ')'
          stop
        else
          icoord = 0
          do ix=1, 3
            do iatom=1, natoms
              icoord = icoord + 1
              coefs_tvec(icoord,imode) = coefs_cart(iatom,ix,ipos)
            enddo
          enddo
        endif
      enddo

      ! derivative of rotational t-vector: t_{iatom,ix,imode}^{rot} = \sum_{iy} asym_tensor_{ix,imode,iy} * cartesian_{iatom,iy}

      t = terms(1:nmodes,iterm1)
      ipos = index_iarr1(t, terms_cart)
      if (ipos==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: Cartesian-coordinate derivative = (', t, ') is not found'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of rotational t-vector = (', terms(:,iterm1), ')'
        stop
      else
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            coefs_tvec(icoord,nmodes+1) = dot_product(coefs_cart(iatom,1:3,ipos), asym_tens(ix,1,1:3))
            coefs_tvec(icoord,nmodes+2) = dot_product(coefs_cart(iatom,1:3,ipos), asym_tens(ix,2,1:3))
            coefs_tvec(icoord,nmodes+3) = dot_product(coefs_cart(iatom,1:3,ipos), asym_tens(ix,3,1:3))
          enddo
        enddo
      endif

      ! derivative of translational t-vector: t_{iatom,ix,imode}^{tran} = \delta_{ix,imode}

      if (all(terms(1:nmodes,iterm1)==0)) then
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            coefs_tvec(icoord,nmodes+3+ix) = 1.0_ark
          enddo
        enddo
      endif

      ! right-hand-side vector and linear system matrix

      if (iterm2==iterm) then
        if (any(terms(:,iterm1)/=0)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: corrupted chain rule expression for derivative = (', terms(:,iterm), ')'
          stop
        endif
        if (ifterm0) then
          tmat0(:,:) = transpose(coefs_tvec(:,:))
        endif
      else
        if (.not.coefs_svec_init(iterm2)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: s-vector derivative = (', terms(:,iterm2), ') is not initialized'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of s-vectors = (', terms(:,iterm), ')'
          write(out, '(a)') 'check if derivative indices are sorted ascendingly wrt the total derivative order'
          stop
        else
          bvec(:,:) = bvec(:,:) - matmul(coefs_svec(:,:,iterm2), coefs_tvec(:,:))
        endif
      endif

    enddo ! ielem

    ! invert t0 matrix (only once since the linear equations matrix is the same for all s-vector derivatives)

    if (ifterm0.and.(.not.tmat0_init)) then
      call matrix_inverse_ark(natoms*3, tmat0, tmat0_inv)
      tmat0_init = .true.
      !do i=1, natoms*3
      !  do j=1, natoms*3
      !    write(out, '(1x,i3,1x,i3,1x,f)') i,j, tmat0_inv(j,i)
      !  enddo
      !enddo
    endif

    ! compute derivative of s-vectors (solve system of linear equations)

    if (tmat0_init) then
      coefs_svec(:,:,iterm) = matmul(bvec, tmat0_inv)
      coefs_svec_init(iterm) = .true.
    else
      write(out, '(/a)') 'deriv_svec error: inverse of t-matrix is not initialized'
      stop
    endif

!================= test for near-singular point =====================!
!
!    etamat(:,:) = 0.0_ark
!    do icoord=1, natoms*3
!      etamat(icoord,icoord) = 0.0_ark
!    enddo
!    etamat(6,6) = 1.0_rk!(-0.01_ark)*(-0.01_ark)
!    etamat(9,9) = 1.0_rk!(-0.01_ark)*(-0.01_ark)
!    cs = matmul(etamat, coefs_svec(:,:,iterm))
!
!
!    write(out, '(/1x,a,<nmodes>(1x,i3),1x,a)') 's-vector derivative = (', terms(:,iterm), ')'
!
!    do i=1, natoms*3
!      do j=1, natoms*3
!        write(out, '(1x,i3,1x,i3,3x,f30.12)') i, j, cs(i,j)
!      enddo
!    enddo
!
!====================================================================!

  enddo ! iterm

end subroutine deriv_svec


!################################################################################


! Computes derivatives of s-vectors from derivatives of Cartesian coordinates (3N-5 vibrational coordinates case).

subroutine deriv_svec_3n5(molec, nterms_cart, terms_cart, coefs_cart, nterms, terms, nterms_chr, ind_chr, coefs_svec)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: nterms_cart
  integer(ik), intent(in)         :: terms_cart(molec%nmodes,nterms_cart)
  real(ark), intent(in)           :: coefs_cart(molec%natoms,3,nterms_cart)
  integer(ik), intent(in)         :: nterms
  integer(ik), intent(in)         :: terms(molec%nmodes,nterms)
  integer(ik), intent(in)         :: nterms_chr(nterms)
  integer(ik), intent(in)         :: ind_chr(2,maxval(nterms_chr),nterms)
  real(ark), intent(out)          :: coefs_svec(molec%natoms*3,molec%natoms*3,nterms)

  integer(ik) :: ncoords, iterm, iterm1, iterm2, ielem, i, j, t(molec%nmodes), ipos, imode, ix, iatom, icoord, info, natoms, nmodes
  real(ark) :: tmat0(molec%natoms*3,molec%natoms*3), tmat0_inv(molec%natoms*3,molec%natoms*3), bvec(molec%natoms*3,molec%natoms*3), coefs_tvec(molec%natoms*3,molec%natoms*3), asym_tens(3,3,3)
  logical :: coefs_svec_init(nterms), tmat0_init, ifterm0

  natoms = molec%natoms
  nmodes = molec%nmodes

  ! Levi-Civita tensor

  asym_tens = 0.0
  asym_tens(1,2,3) = 1.0_ark
  asym_tens(1,3,2) =-1.0_ark
  asym_tens(2,1,3) =-1.0_ark
  asym_tens(2,3,1) = 1.0_ark
  asym_tens(3,1,2) = 1.0_ark
  asym_tens(3,2,1) =-1.0_ark

  ! start recursive calculations of s-vectors derivatives

  coefs_svec(:,:,:) = 0.0
  coefs_svec_init(:) = .false.

  tmat0 = 0.0
  tmat0_init = .false.

  do iterm=1, nterms

    ! set up system of linear equations for terms(:,iterm)-derivative of s-vectors

    bvec = 0.0
    ifterm0 = .false.

    ! for 0-order derivative right-hand side is delta symbol
    if (all(terms(1:nmodes,iterm)==0)) then
      forall(i=1:natoms*3) bvec(i,i) = 1.0_ark
      ifterm0 = .true.
    endif

    ! loop over elements of the chain-rule expression for terms(:,iterm)-derivative of t*s product

    do ielem=1, nterms_chr(iterm)

      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)

      coefs_tvec(:,:) = 0.0

      ! derivative of vibrational t-vector: t_{iatom,ix,imode}^{vib} = d cartesian_{iatom,ix} / d internal_{imode}

      do imode=1, nmodes
        t = terms(1:nmodes,iterm1)
        t(imode) = t(imode) + 1
        ipos = index_iarr1(t, terms_cart)
        if (ipos==0) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: Cartesian-coordinate derivative = (', t, ') is not found'
          write(out, '(a,1x,i3,1x,a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of vibrational t-vector for', imode, 'mode = (', terms(:,iterm1), ')'
          stop
        else
          icoord = 0
          do ix=1, 3
            do iatom=1, natoms
              icoord = icoord + 1
              coefs_tvec(icoord,imode) = coefs_cart(iatom,ix,ipos)
            enddo
          enddo
        endif
      enddo

      ! derivative of rotational t-vector: t_{iatom,ix,imode}^{rot} = \sum_{iy} asym_tensor_{ix,imode,iy} * cartesian_{iatom,iy}

      t = terms(1:nmodes,iterm1)
      ipos = index_iarr1(t, terms_cart)
      if (ipos==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: Cartesian-coordinate derivative = (', t, ') is not found'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of rotational t-vector = (', terms(:,iterm1), ')'
        stop
      else
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            coefs_tvec(icoord,nmodes+1) = dot_product(coefs_cart(iatom,1:3,ipos), asym_tens(ix,1,1:3))
            coefs_tvec(icoord,nmodes+2) = dot_product(coefs_cart(iatom,1:3,ipos), asym_tens(ix,2,1:3))
          enddo
        enddo
      endif

      ! derivative of translational t-vector: t_{iatom,ix,imode}^{tran} = \delta_{ix,imode}

      if (all(terms(1:nmodes,iterm1)==0)) then
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            coefs_tvec(icoord,nmodes+2+ix) = 1.0_ark
          enddo
        enddo
      endif

      ! right-hand-side vector and linear system matrix

      if (iterm2==iterm) then
        if (any(terms(:,iterm1)/=0)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: corrupted chain rule expression for derivative = (', terms(:,iterm), ')'
          stop
        endif
        if (ifterm0) then
          tmat0(:,:) = transpose(coefs_tvec(:,:))
        endif
      else
        if (.not.coefs_svec_init(iterm2)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_svec error: s-vector derivative = (', terms(:,iterm2), ') is not initialized'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of s-vectors = (', terms(:,iterm), ')'
          write(out, '(a)') 'check if derivative indices are sorted ascendingly wrt the total derivative order'
          stop
        else
          bvec(:,:) = bvec(:,:) - matmul(coefs_svec(:,:,iterm2), coefs_tvec(:,:))
        endif
      endif

    enddo ! ielem

    ! invert t0 matrix (only once since the linear equations matrix is the same for all s-vector derivatives)

    if (ifterm0.and.(.not.tmat0_init)) then
      call matrix_inverse_ark(natoms*3, tmat0, tmat0_inv)
      tmat0_init = .true.
    endif

    ! compute derivative of s-vectors (solve system of linear equations)

    if (tmat0_init) then
      coefs_svec(:,:,iterm) = matmul(bvec, tmat0_inv)
      coefs_svec_init(iterm) = .true.
    else
      write(out, '(/a)') 'deriv_svec error: inverse of t-matrix is not initialized'
      stop
    endif

  enddo ! iterm

end subroutine deriv_svec_3n5


!################################################################################


! Computes derivative of t-matrix from derivatives of Cartesian coordinates.

subroutine deriv_tvec(molec, npoints, nterms_cart, terms_cart, coefs_cart, term, coefs_tvec)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: npoints
  integer(ik), intent(in)         :: nterms_cart
  integer(ik), intent(in)         :: terms_cart(molec%nmodes,nterms_cart)
  real(ark), intent(in)           :: coefs_cart(molec%natoms,3,nterms_cart,npoints)
  integer(ik), intent(in)         :: term(molec%nmodes)
  real(ark), intent(out)          :: coefs_tvec(molec%natoms*3,molec%natoms*3,npoints)

  integer(ik) :: t(molec%nmodes), ipos, imode, ix, iatom, icoord, natoms, nmodes, ipoint
  real(ark) :: asym_tens(3,3,3)

  natoms = molec%natoms
  nmodes = molec%nmodes

  ! Levi-Civita tensor

  asym_tens = 0.0
  asym_tens(1,2,3) = 1.0_ark
  asym_tens(1,3,2) =-1.0_ark
  asym_tens(2,1,3) =-1.0_ark
  asym_tens(2,3,1) = 1.0_ark
  asym_tens(3,1,2) = 1.0_ark
  asym_tens(3,2,1) =-1.0_ark

  coefs_tvec(:,:,:) = 0.0

  ! derivative of vibrational t-vector

  do imode=1, nmodes
    t = term(1:nmodes)
    t(imode) = t(imode) + 1
    ipos = index_iarr1(t, terms_cart)
    if (ipos==0) then
      write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_tvec error: Cartesian-coordinate derivative = (', t, ') is not initialized'
      write(out, '(a,1x,i3,1x,a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of vibrational t-vector for', imode, 'mode = (', term, ')'
      stop
    else
      do ipoint=1, npoints
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            coefs_tvec(icoord,imode,ipoint) = coefs_cart(iatom,ix,ipos,ipoint)
          enddo
        enddo
      enddo
    endif
  enddo

  ! derivative of rotational t-vector

  t = term(1:nmodes)
  ipos = index_iarr1(t, terms_cart)
  if (ipos==0) then
    write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_tvec error: Cartesian-coordinate derivative = (', t, ') is not initialized'
    write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of rotational t-vector = (', term, ')'
    stop
  else
    do ipoint=1, npoints
      icoord = 0
      do ix=1, 3
        do iatom=1, natoms
          icoord = icoord + 1
          coefs_tvec(icoord,nmodes+1,ipoint) = dot_product(coefs_cart(iatom,1:3,ipos,ipoint), asym_tens(ix,1,1:3))
          coefs_tvec(icoord,nmodes+2,ipoint) = dot_product(coefs_cart(iatom,1:3,ipos,ipoint), asym_tens(ix,2,1:3))
          coefs_tvec(icoord,nmodes+3,ipoint) = dot_product(coefs_cart(iatom,1:3,ipos,ipoint), asym_tens(ix,3,1:3))
        enddo
      enddo
    enddo
  endif

  ! derivative of translational t-vector

  if (all(term(1:nmodes)==0)) then
    do ipoint=1, npoints
      icoord = 0
      do ix=1, 3
        do iatom=1, natoms
          icoord = icoord + 1
          coefs_tvec(icoord,nmodes+3+ix,ipoint) = 1.0_ark
        enddo
      enddo
    enddo
  endif

end subroutine deriv_tvec


subroutine deriv_tvec_cmplx(molec, npoints, nterms_cart, terms_cart, coefs_cart, term, coefs_tvec)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: npoints
  integer(ik), intent(in)         :: nterms_cart
  integer(ik), intent(in)         :: terms_cart(molec%nmodes,nterms_cart)
  real(ark), intent(in)           :: coefs_cart(molec%natoms,3,nterms_cart,npoints)
  integer(ik), intent(in)         :: term(molec%nmodes)
  complex(ark), intent(out)       :: coefs_tvec(molec%natoms*3,molec%natoms*3,npoints)

  integer(ik) :: t(molec%nmodes), ipos, imode, ix, iatom, icoord, natoms, nmodes, ipoint
  real(ark) :: asym_tens(3,3,3)

  natoms = molec%natoms
  nmodes = molec%nmodes

  ! Levi-Civita tensor

  asym_tens = 0.0
  asym_tens(1,2,3) = 1.0_ark
  asym_tens(1,3,2) =-1.0_ark
  asym_tens(2,1,3) =-1.0_ark
  asym_tens(2,3,1) = 1.0_ark
  asym_tens(3,1,2) = 1.0_ark
  asym_tens(3,2,1) =-1.0_ark

  coefs_tvec(:,:,:) = 0.0

  ! derivative of vibrational t-vector

  do imode=1, nmodes
    t = term(1:nmodes)
    t(imode) = t(imode) + 1
    ipos = index_iarr1(t, terms_cart)
    if (ipos==0) then
      write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_tvec error: Cartesian-coordinate derivative = (', t, ') is not initialized'
      write(out, '(a,1x,i3,1x,a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of vibrational t-vector for', imode, 'mode = (', term, ')'
      stop
    else
      do ipoint=1, npoints
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            coefs_tvec(icoord,imode,ipoint) = coefs_cart(iatom,ix,ipos,ipoint)
            if (imode==6) then
              coefs_tvec(icoord,imode,ipoint) = coefs_tvec(icoord,imode,ipoint) * cmplx(0.0_ark,1.0_ark)
            endif
          enddo
        enddo
      enddo
    endif
  enddo

  ! derivative of rotational t-vector

  t = term(1:nmodes)
  ipos = index_iarr1(t, terms_cart)
  if (ipos==0) then
    write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_tvec error: Cartesian-coordinate derivative = (', t, ') is not initialized'
    write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of rotational t-vector = (', term, ')'
    stop
  else
    do ipoint=1, npoints
      icoord = 0
      do ix=1, 3
        do iatom=1, natoms
          icoord = icoord + 1
          coefs_tvec(icoord,nmodes+1,ipoint) = dot_product(coefs_cart(iatom,1:3,ipos,ipoint), asym_tens(ix,1,1:3))
          coefs_tvec(icoord,nmodes+2,ipoint) = dot_product(coefs_cart(iatom,1:3,ipos,ipoint), asym_tens(ix,2,1:3))
          coefs_tvec(icoord,nmodes+3,ipoint) = dot_product(coefs_cart(iatom,1:3,ipos,ipoint), asym_tens(ix,3,1:3))
          coefs_tvec(icoord,nmodes+3,ipoint) = coefs_tvec(icoord,nmodes+3,ipoint) * cmplx(0.0_ark,1.0_ark)
        enddo
      enddo
    enddo
  endif

  ! derivative of translational t-vector

  if (all(term(1:nmodes)==0)) then
    do ipoint=1, npoints
      icoord = 0
      do ix=1, 3
        do iatom=1, natoms
          icoord = icoord + 1
          coefs_tvec(icoord,nmodes+3+ix,ipoint) = 1.0_ark
        enddo
      enddo
    enddo
  endif

end subroutine deriv_tvec_cmplx


!################################################################################


! Computes derivatives of kinetic energy G-matrix from derivatives of s-vectors.

subroutine deriv_gmat(molec, nterms_svec, terms_svec, coefs_svec, nterms, terms, nterms_chr, ind_chr, coefs_gmat)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: nterms_svec
  integer(ik), intent(in)         :: terms_svec(molec%nmodes,nterms_svec)
  real(ark), intent(in)           :: coefs_svec(molec%natoms*3,molec%natoms*3,nterms_svec)
  integer(ik), intent(in)         :: nterms
  integer(ik), intent(in)         :: terms(molec%nmodes,nterms)
  integer(ik), intent(in)         :: nterms_chr(nterms)
  integer(ik), intent(in)         :: ind_chr(2,maxval(nterms_chr),nterms)
  real(ark), intent(out)          :: coefs_gmat(molec%natoms*3,molec%natoms*3,nterms)

  integer(ik) :: natoms, nmodes, iterm, ielem, iatom, iterm1, iterm2, ipos1, ipos2, icoord, jcoord, ix, verbose, imode
  real(ark) :: inv_masses(molec%natoms*3), coef

  natoms = molec%natoms
  nmodes = molec%nmodes
  verbose = molec%expkeo_verbose

  if (verbose>=6) then
    write(out, '(/a)') 'deriv_gmat/start: compute derivatives of G-matrix'
  endif

  ! precompute inverse masses
  inv_masses = 0.0
  icoord = 0
  do ix=1, 3
    do iatom=1, natoms
      icoord = icoord + 1
      inv_masses(icoord) = 1.0_ark / molec%masses(iatom)
    enddo
  enddo

  coefs_gmat = 0.0

  do iterm=1, nterms

    ! loop over elements of the chain-rule expression for derivative of s*s product defined by terms(1:nmodes,iterm) indexes

    do ielem=1, nterms_chr(iterm)

      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)

      ipos1 = index_iarr1(terms(1:nmodes,iterm1), terms_svec)
      ipos2 = index_iarr1(terms(1:nmodes,iterm2), terms_svec)
      if (ipos1==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_gmat error: s-vector derivative = (', terms(:,iterm1), ') is not initialized'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of G-matrix = (', terms(:,iterm), ')'
        stop
      endif
      if (ipos2==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_gmat error: s-vector derivative = (', terms(:,iterm2), ') is not initialized'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of G-matrix = (', terms(:,iterm), ')'
        stop
      endif

      do icoord=1, natoms*3
        do jcoord=1, natoms*3
          coef = sum( coefs_svec(jcoord,:,ipos1) * coefs_svec(icoord,:,ipos2) * inv_masses(:) )
          coefs_gmat(jcoord,icoord,iterm) = coefs_gmat(jcoord,icoord,iterm) + coef
        enddo
      enddo

    enddo ! ielem

  enddo ! iterm


  if (verbose>=6) then
    do icoord=1, natoms*3
      do jcoord=1, icoord
        write(out, '(/1x,a,1x,i3,1x,i3,1x,a)') 'derivatives of G(', jcoord, icoord, ')'
        do iterm=1, nterms
          if (abs(coefs_gmat(jcoord,icoord,iterm))>1.0d-08) then
            write(out, '(1x,i6,3x,a,<nmodes>(1x,i3),1x,a,3x,f)') iterm, '(', terms(1:nmodes,iterm), ')', coefs_gmat(jcoord,icoord,iterm)
          endif
        enddo
      enddo
    enddo
  endif

  if (verbose>=6) then
    write(out, '(/a)') 'deriv_gmat/done'
  endif

end subroutine deriv_gmat


!################################################################################




! Computes derivatives of pseudo-potential from derivatives of s-vectors.

subroutine deriv_pseudo(molec, nterms_svec, terms_svec, coefs_svec, nterms, terms, nterms_chr, ind_chr, coefs)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: nterms_svec
  integer(ik), intent(in)         :: terms_svec(molec%nmodes,nterms_svec)
  real(ark), intent(in)           :: coefs_svec(molec%natoms*3,molec%natoms*3,nterms_svec)
  integer(ik), intent(in)         :: nterms
  integer(ik), intent(in)         :: terms(molec%nmodes,nterms)
  integer(ik), intent(in)         :: nterms_chr(nterms)
  integer(ik), intent(in)         :: ind_chr(2,maxval(nterms_chr),nterms)
  real(ark), intent(out)          :: coefs(nterms)

  integer(ik) :: natoms, nmodes, iterm, ielem, iatom, iterm1, iterm2, ipos, imode, jmode, ix, iy, icoord, t(molec%nmodes), i, j, k, l, info, verbose
  real(ark) :: inv_masses(molec%natoms), pseudo1(molec%natoms), pseudo2(molec%natoms), pseudo3(molec%natoms), pseudo4(molec%natoms), &!
               s10(molec%natoms*3,molec%natoms,3), s20(molec%natoms*3,molec%natoms,3), s11(molec%natoms*3,molec%natoms,3,molec%nmodes), &!
               s21(molec%natoms*3,molec%natoms,3,molec%nmodes), s22(molec%nmodes,molec%natoms,3), s11_(3), s21_(3), tmat(molec%natoms*3,molec%natoms*3), &!
               asym_tens(3,3,3), asym_tens2(3,3,3,3)

  natoms = molec%natoms
  nmodes = molec%nmodes
  verbose = molec%expkeo_verbose

  if (verbose>=6) then
    write(out, '(/a)') 'deriv_pseudo/start: compute derivatives of pseudo-potential'
  endif

  ! precompute 3D levi-civita tensor

  asym_tens = 0.0
  asym_tens(1,2,3) = 1.0_ark
  asym_tens(1,3,2) =-1.0_ark
  asym_tens(2,1,3) =-1.0_ark
  asym_tens(2,3,1) = 1.0_ark
  asym_tens(3,1,2) = 1.0_ark
  asym_tens(3,2,1) =-1.0_ark
  forall(i=1:3,j=1:3,k=1:3,l=1:3) asym_tens2(i,j,k,l) = dot_product(asym_tens(:,i,j), asym_tens(:,k,l))


  ! precompute inverse masses

  inv_masses = 0.0
  do iatom=1, natoms
    inv_masses(iatom) = 1.0_ark / molec%masses(iatom)
  enddo


  ! compute derivatives of pseudo-potential

  coefs = 0.0

  do iterm=1, nterms

    pseudo1 = 0.0
    pseudo2 = 0.0
    pseudo3 = 0.0
    pseudo4 = 0.0

    do ielem=1, nterms_chr(iterm)

      s10 = 0.0
      s20 = 0.0
      s11 = 0.0
      s21 = 0.0
      s22 = 0.0

      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)

      ipos = index_iarr1(terms(1:nmodes,iterm1), terms_svec)
      if (ipos==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_pseudo error: s-vector derivative = (', terms(:,iterm1), ') is not found'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of pseudo-potential = (', terms(:,iterm), ')'
        stop
      endif
      icoord = 0
      do ix=1, 3
        do iatom=1, natoms
          icoord = icoord + 1
          s10(:,iatom,ix) = coefs_svec(:,icoord,ipos)
        enddo
      enddo

      ipos = index_iarr1(terms(1:nmodes,iterm2), terms_svec)
      if (ipos==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_pseudo error: s-vector derivative = (', terms(:,iterm2), ') is not found'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of pseudo-potential = (', terms(:,iterm), ')'
        stop
      endif
      icoord = 0
      do ix=1, 3
        do iatom=1, natoms
          icoord = icoord + 1
          s20(:,iatom,ix) = coefs_svec(:,icoord,ipos)
        enddo
      enddo

      do imode=1, nmodes

        t = terms(1:nmodes,iterm1)
        t(imode) = t(imode) + 1
        ipos = index_iarr1(t, terms_svec)
        if (ipos==0) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_pseudo error: s-vector derivative = (', t, ') is not found'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of pseudo-potential = (', terms(:,iterm), ')'
          stop
        endif
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            s11(:,iatom,ix,imode) = coefs_svec(:,icoord,ipos)
          enddo
        enddo

        t = terms(1:nmodes,iterm2)
        t(imode) = t(imode) + 1
        ipos = index_iarr1(t, terms_svec)
        if (ipos==0) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_pseudo error: s-vector derivative = (', t, ') is not found'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of pseudo-potential = (', terms(:,iterm), ')'
          stop
        endif
        icoord = 0
        do ix=1, 3
          do iatom=1, natoms
            icoord = icoord + 1
            s21(:,iatom,ix,imode) = coefs_svec(:,icoord,ipos)
          enddo
        enddo

        do jmode=1, nmodes
          t = terms(1:nmodes,iterm2)
          t(imode) = t(imode) + 1
          t(jmode) = t(jmode) + 1
          ipos = index_iarr1(t, terms_svec)
          if (ipos==0) then
            write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'deriv_pseudo error: s-vector derivative = (', t, ') is not found'
            write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of pseudo-potential = (', terms(:,iterm), ')'
            stop
          endif
          icoord = 0
          do ix=1, 3
            do iatom=1, natoms
              icoord = icoord + 1
              s22(imode,iatom,ix) = s22(imode,iatom,ix) + coefs_svec(jmode,icoord,ipos)
            enddo
          enddo
        enddo ! jmode

      enddo ! imode

      do iatom=1, natoms

        do iy=1, 3
          do ix=1, 3
            tmat(ix,iy) = sum( asym_tens2(ix,iy,1:3,1:3) * s20(nmodes+1:nmodes+3,iatom,1:3) )
          enddo
        enddo
        pseudo1(iatom) = pseudo1(iatom) + sum( tmat(1:3,1:3) * s10(nmodes+1:nmodes+3,iatom,1:3) )

        do ix=1, 3
          do imode=1, nmodes
            tmat(imode,ix) = sum( s21(nmodes+1:nmodes+3,iatom,1:3,imode) * transpose(asym_tens(1:3,1:3,ix)) )
          enddo
        enddo
        pseudo2(iatom) = pseudo2(iatom) + sum( tmat(1:nmodes,1:3) * s10(1:nmodes,iatom,1:3) )

        pseudo3(iatom) = pseudo3(iatom) + sum( s10(1:nmodes,iatom,1:3) * s22(1:nmodes,iatom,1:3) )

        s11_ = (/(sum( (/(s11(imode,iatom,ix,imode), imode=1, nmodes)/) ), ix=1, 3)/)
        s21_ = (/(sum( (/(s21(imode,iatom,ix,imode), imode=1, nmodes)/) ), ix=1, 3)/)
        pseudo4(iatom) = pseudo4(iatom) + sum( s11_ * s21_ )

      enddo ! iatom

    enddo ! ielem

    pseudo1 = pseudo1 * 0.5_ark
    pseudo4 = pseudo4 * 0.5_ark

    coefs(iterm) = sum((pseudo1 + pseudo2 - pseudo3 - pseudo4) * inv_masses) * 0.25_ark

  enddo ! iterm

  if (verbose>=6) then
    write(out, '(/1x,a)') 'derivatives of pseudo-potential'
    do iterm=1, nterms
      if (abs(coefs(iterm))>1.0d-08) then
        write(out, '(1x,i6,3x,a,<nmodes>(1x,i3),1x,a,3x,f)') iterm, '(', terms(1:nmodes,iterm), ')', coefs(iterm)
      endif
    enddo
  endif

  if (verbose>=6) then
    write(out, '(/a)') 'deriv_pseudo/done'
  endif

end subroutine deriv_pseudo



!################################################################################


! Computes inverse of matrix.

subroutine matrix_inverse_ark(nelem, mat, invmat, info)

  integer(ik), intent(in) :: nelem
  real(ark), intent(in) :: mat(nelem,nelem)
  real(ark), intent(out) :: invmat(nelem,nelem)
  integer(ik), intent(out), optional :: info

  integer(ik) lwork, info_, i, j
  double precision work1(1), matd(nelem,nelem), matu(nelem,nelem), matvt(nelem,nelem), invmatd(nelem,nelem), mat_d(nelem,nelem), sv(nelem), tmat(nelem,nelem), tol
  double precision, allocatable :: work(:)

  tol = epsilon(1.0d0) ! 1.0d-14
  matd = dble(mat)

  lwork = -1
  call dgesvd('A', 'A', nelem, nelem, matd, nelem, sv, matu, nelem, matvt, nelem, work1, lwork, info_)
  lwork = int(work1(1), kind=ik)

  allocate(work(lwork), stat=info_)
  if (info_/=0) then
    write(out, '(/a,1x,i6)') 'matrix_inverse_ark error: failed to allocate workspace array for diagonalization, size =', lwork
    stop
  endif

  call dgesvd('A', 'A', nelem, nelem, matd, nelem, sv, matu, nelem, matvt, nelem, work, lwork, info_)

  if (info_/=0) then
    write(out, '(/a,1x,i6)') 'matrix_inverse_ark error: SVD failed, info =', info_
    stop
  endif

  if (present(info)) info = 0

  mat_d = 0.0
  do i=1, nelem
    if (sv(i)>=tol) then
      mat_d(i,i) = 1.0d0/sv(i)
    else
      if (present(info)) then
        info = i
        mat_d(i,i) = 0.0
      else
        write(out, '(/a)') 'matrix_inverse_ark error: matrix is singular'
        write(out, '(a)') 'matrix'
        do j=1, nelem
          write(out, '(<nelem>(1x,f10.6))') mat(1:nelem,j)
        enddo
        write(out, '(a)') 'singular elements'
        do j=1, nelem
          write(out, '(1x,f)') sv(j)
        enddo
        stop
      endif
    endif
  enddo

  call dgemm('N', 'N', nelem, nelem, nelem, 1.0d0, mat_d, nelem, matvt, nelem, 0.0d0, tmat, nelem)
  call dgemm('N', 'N', nelem, nelem, nelem, 1.0d0, matu, nelem, tmat, nelem, 0.0d0, invmatd, nelem)

  invmat = real(invmatd,kind=ark)

  deallocate(work)

end subroutine matrix_inverse_ark


subroutine matrix_inverse2_ark(nelem, mat, invmat)

  integer(ik), intent(in) :: nelem
  real(ark), intent(in) :: mat(nelem,nelem)
  real(ark), intent(out) :: invmat(nelem,nelem)

  integer(ik) lwork, info, ipiv(nelem)
  double precision matd(nelem,nelem), work(nelem*2)

  matd = dble(mat)

  call dgetrf( nelem, nelem, matd, nelem, ipiv, info )

  if (info/=0) then
    write(out, '(/a,1x,i6,1x,a)') 'matrix_inverse2_ark error: LU factorization of matrix failed, info =', info, ', singular matrix?'
    stop
  endif

  lwork = size(work)

  call dgetri( nelem, matd, nelem, ipiv, work, lwork, info )

  if (info/=0) then
    write(out, '(/a,1x,i6)') 'matrix_inverse2_ark error: inversion of LU-factored matrix failed, info =', info
    stop
  endif

  invmat = real(matd,kind=ark)

end subroutine matrix_inverse2_ark


!################################################################################


subroutine expansion_terms_pseudo2cart(nmodes, nterms, terms)

  integer(ik), intent(in)                 :: nmodes
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 1000
  !integer(ik), parameter :: terms_incr = 7000
  integer(ik) :: nterms_aug, iterm, imode, t(nmodes), ipos, info, terms_size, iorder, jmode
  integer(ik), allocatable :: tmp(:,:)

  terms_size = nterms
  nterms_aug = nterms

  ! required derivatives of Cartesian coordinates:
  ! dX/dxi_n, d^2X/(dxi_n^2), d^3X/(dxi_m^2 dxi_n) for n,m=1..nmodes

  do iterm=1, nterms
    do imode=1, nmodes
      do iorder=1, 2
        do jmode=0, nmodes
          t = terms(:,iterm)
          t(imode) = t(imode) + iorder
          if (jmode>0.and.iorder==2) t(jmode) = t(jmode) + 1
          !ipos = index_iarr1(t, terms(:,1:nterms_aug))
          ipos = index_iarr1_2(nmodes, nterms_aug, t, terms(:,1:nterms_aug))
          if (ipos==0) then
            if (nterms_aug+1>terms_size) then
              terms_size = terms_size + terms_incr
              allocate(tmp(nmodes,terms_size), stat=info)
              if (info/=0) then
                write(out, '(/a/a,2(1x,i6))') 'expansion_terms_pseudo2cart error: failed to allocate tmp(nmodes,terms_size)', 'nmodes, terms_size =', nmodes, terms_size
                stop
              endif
              tmp(:,1:nterms_aug) = terms(:,1:nterms_aug)
              deallocate(terms)
              call move_alloc(tmp, terms)
            endif
            nterms_aug = nterms_aug + 1
            terms(:,nterms_aug) = t
          endif
        enddo
      enddo
    enddo
  enddo

  allocate(tmp(nmodes,nterms_aug), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') 'expansion_terms_pseudo2cart error: failed to allocate tmp(nvar,nmodes)', 'nmodes, nterms_aug =', nmodes, nterms_aug
    stop
  endif

  tmp = terms(:,1:nterms_aug)
  deallocate(terms)
  call move_alloc(tmp, terms)
  nterms = nterms_aug

  if (allocated(tmp)) deallocate(tmp)

end subroutine expansion_terms_pseudo2cart


!################################################################################


subroutine expansion_terms_svec2cart(nmodes, nterms, terms)

  integer(ik), intent(in)                 :: nmodes
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 1000
  !integer(ik), parameter :: terms_incr = 7000
  integer(ik) :: nterms_aug, iterm, imode, t(nmodes), ipos, info, terms_size
  integer(ik), allocatable :: tmp(:,:)

  terms_size = nterms
  nterms_aug = nterms

  do iterm=1, nterms
    do imode=1, nmodes
      t = terms(:,iterm)
      t(imode) = t(imode) + 1
      ipos = index_iarr1(t, terms(:,1:nterms_aug))
      if (ipos==0) then
        if (nterms_aug+1>terms_size) then
          terms_size = terms_size + terms_incr
          allocate(tmp(nmodes,terms_size), stat=info)
          if (info/=0) then
            write(out, '(/a/a,2(1x,i6))') 'expansion_terms_svec2cart error: failed to allocate tmp(nmodes,terms_size)', 'nmodes, terms_size =', nmodes, terms_size
            stop
          endif
          tmp(:,1:nterms_aug) = terms(:,1:nterms_aug)
          deallocate(terms)
          call move_alloc(tmp, terms)
        endif
        nterms_aug = nterms_aug + 1
        terms(:,nterms_aug) = t
      endif
    enddo
  enddo

  allocate(tmp(nmodes,nterms_aug), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &!
    'expansion_terms_svec2cart error: failed to allocate tmp(nvar,nmodes)', 'nmodes, nterms_aug =', nmodes, nterms_aug
    stop
  endif

  tmp = terms(:,1:nterms_aug)
  deallocate(terms)
  call move_alloc(tmp, terms)
  nterms = nterms_aug

  if (allocated(tmp)) deallocate(tmp)

end subroutine expansion_terms_svec2cart


!################################################################################


subroutine expansion_terms_dgmat2gmat(nmodes, nterms, terms)

  integer(ik), intent(in)                 :: nmodes
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 1000
  !integer(ik), parameter :: terms_incr = 7000
  integer(ik) :: nterms_aug, iterm, imode, t(nmodes), ipos, info, terms_size
  integer(ik), allocatable :: tmp(:,:)

  terms_size = nterms
  nterms_aug = nterms

  do iterm=1, nterms
    do imode=1, nmodes
      t = terms(:,iterm)
      t(imode) = t(imode) + 1
      ipos = index_iarr1(t, terms(:,1:nterms_aug))
      if (ipos==0) then
        if (nterms_aug+1>terms_size) then
          terms_size = terms_size + terms_incr
          allocate(tmp(nmodes,terms_size), stat=info)
          if (info/=0) then
            write(out, '(/a/a,2(1x,i6))') 'expansion_terms_dgmat2gmat error: failed to allocate tmp(nmodes,terms_size)', 'nmodes, terms_size =', nmodes, terms_size
            stop
          endif
          tmp(:,1:nterms_aug) = terms(:,1:nterms_aug)
          deallocate(terms)
          call move_alloc(tmp, terms)
        endif
        nterms_aug = nterms_aug + 1
        terms(:,nterms_aug) = t
      endif
    enddo
  enddo

  allocate(tmp(nmodes,nterms_aug), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &!
    'expansion_terms_dgmat2gmat error: failed to allocate tmp(nvar,nmodes)', 'nmodes, nterms_aug =', nmodes, nterms_aug
    stop
  endif

  tmp = terms(:,1:nterms_aug)
  deallocate(terms)
  call move_alloc(tmp, terms)
  nterms = nterms_aug

  if (allocated(tmp)) deallocate(tmp)

end subroutine expansion_terms_dgmat2gmat


!################################################################################


subroutine expansion_terms_pseudo2svec(nmodes, nterms, terms)

  integer(ik), intent(in)                 :: nmodes
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 1000
  !integer(ik), parameter :: terms_incr = 7000
  integer(ik) :: nterms_aug, iterm, imode, jmode, t(nmodes), ipos, info, terms_size
  integer(ik), allocatable :: tmp(:,:)

  terms_size = nterms
  nterms_aug = nterms

  do iterm=1, nterms
    do imode=1, nmodes
      do jmode=0, nmodes
        t = terms(:,iterm)
        t(imode) = t(imode) + 1
        if (jmode>0) t(jmode) = t(jmode) + 1
        ipos = index_iarr1(t, terms(:,1:nterms_aug))
        if (ipos==0) then
          if (nterms_aug+1>terms_size) then
            terms_size = terms_size + terms_incr
            allocate(tmp(nmodes,terms_size), stat=info)
            if (info/=0) then
              write(out, '(/a/a,2(1x,i6))') &
              'expansion_terms_pseudo2svec error: failed to allocate tmp(nmodes,terms_size)', 'nmodes, terms_size =', &!
              nmodes, terms_size
              stop
            endif
            tmp(:,1:nterms_aug) = terms(:,1:nterms_aug)
            deallocate(terms)
            call move_alloc(tmp, terms)
          endif
          nterms_aug = nterms_aug + 1
          terms(:,nterms_aug) = t
        endif
      enddo
    enddo
  enddo

  allocate(tmp(nmodes,nterms_aug), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &!
    'expansion_terms_pseudo2svec error: failed to allocate tmp(nvar,nmodes)', 'nmodes, nterms_aug =', nmodes, nterms_aug
    stop
  endif

  tmp = terms(:,1:nterms_aug)
  deallocate(terms)
  call move_alloc(tmp, terms)
  nterms = nterms_aug

  if (allocated(tmp)) deallocate(tmp)

end subroutine expansion_terms_pseudo2svec


!################################################################################


subroutine expansion_terms_complete(nvar, nterms, terms)

  integer(ik), intent(in)                 :: nvar
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 1000
  !integer(ik), parameter :: terms_incr = 7000
  integer(ik) :: nterms_aug, iterm, ipos, info, nelem, ielem, terms_size
  integer(ik), allocatable :: term2(:,:,:), tmp(:,:)

  terms_size = nterms
  nterms_aug = nterms

  do iterm=1, nterms
    call chain_rule_product(nvar, terms(1:nvar,iterm), 2, nelem, term2)
    do ielem=1, nelem
      ipos = index_iarr1(term2(:,1,ielem), terms(:,1:nterms_aug))
      if (ipos==0) then
        if (nterms_aug+1>terms_size) then
          terms_size = terms_size + terms_incr
          allocate(tmp(nvar,terms_size), stat=info)
          if (info/=0) then
            write(out, '(/a/a,2(1x,i6))') &!
            'expansion_terms_complete error: failed to allocate tmp(nvar,terms_size)', 'nvar, terms_size =', &!
            nvar, terms_size
            stop
          endif
          tmp(:,1:nterms_aug) = terms(:,1:nterms_aug)
          deallocate(terms)
          call move_alloc(tmp, terms)
        endif
        nterms_aug = nterms_aug + 1
        terms(:,nterms_aug) = term2(:,1,ielem)
      endif
    enddo
  enddo

  allocate(tmp(nvar,nterms_aug), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &!
    'expansion_terms_complete error: failed to allocate tmp(nvar,nterms_aug)', 'nvar, nterms_aug =', nvar, nterms_aug
    stop
  endif

  tmp = terms(:,1:nterms_aug)
  deallocate(terms)
  call move_alloc(tmp, terms)
  nterms = nterms_aug

  if (allocated(tmp)) deallocate(tmp)
  if (allocated(term2)) deallocate(term2)

end subroutine expansion_terms_complete


!################################################################################


function index_iarr1(ind, ind2) result(ipos)

  integer(ik), intent(in) :: ind(:), ind2(:,:)
  integer(ik) :: ipos, i, nelem

  nelem = size(ind2,dim=2)

  ipos = 0
  do i=1, nelem
    if (all(ind==ind2(:,i))) then
      ipos = i
      exit
    endif
  enddo

end function index_iarr1


!################################################################################


function index_iarr1_2(dimen, nelem, ind, ind2) result(ipos)

  integer(ik), intent(in) :: dimen, nelem, ind(dimen), ind2(dimen, nelem)
  integer(ik) :: ipos, i, j, nelem_, nelem0, ielem(nelem), ielem0(nelem)

  nelem0 = nelem
  forall(j=1:nelem0) ielem0(j) = j
  do i=1, dimen
    nelem_ = 0
    do j=1, nelem0
      if (ind2(i,ielem0(j))==ind(i)) then
        nelem_ = nelem_ + 1
        ielem(nelem_) = ielem0(j)
      endif
    enddo
    nelem0 = nelem_
    ielem0(1:nelem0) = ielem(1:nelem0)
  enddo

  if (nelem0>1) then
    write(out, '(/a)') 'index_iarr1_2: multiple occurrence of ind in ind2'
    stop
  endif

  if (nelem0==0) then
    ipos = 0
  else
    ipos = ielem0(1)
  endif

end function index_iarr1_2


!################################################################################


subroutine chain_rule_product(nvar, term, nelem_prod, nterms, term2)

  integer(ik), intent(in)                         :: nvar
  integer(ik), intent(in)                         :: term(nvar)
  integer(ik), intent(in)                         :: nelem_prod
  integer(ik), intent(out)                        :: nterms
  integer(ik), allocatable, intent(out), optional :: term2(:,:,:)

  integer(ik) :: ivar, ideriv, iterm, jterm, ielem, info
  integer(ik), allocatable :: term0(:,:,:)

  nterms = 1
  do ivar=1, nvar
    do ideriv=1, term(ivar)
      nterms = nterms*nelem_prod
    enddo
  enddo

  if (present(term2)) then

    if (allocated(term2)) deallocate(term2)
    allocate(term2(nvar,nelem_prod,nterms), term0(nvar,nelem_prod,nterms), stat=info)
    if (info/=0) then
      write(out, '(/a/a,3(1x,i6))') &!
      'chain_rule_product error: failed to allocate term0/term2(nvar,nelem_prod,nterms)', &!
      'nvar, nelem_prod, nterms =', nvar, nelem_prod, nterms
      stop
    endif

    term0 = 0
    term2 = 0
    nterms = 1
    do ivar=1, nvar
      do ideriv=1, term(ivar)
        jterm = 0
        do iterm=1, nterms
          do ielem=1, nelem_prod
            jterm = jterm + 1
            term2(1:nvar,1:nelem_prod,jterm) = term0(1:nvar,1:nelem_prod,iterm)
            term2(ivar,ielem,jterm) = term2(ivar,ielem,jterm) + 1
          enddo
        enddo
        nterms = jterm
        term0 = term2
      enddo
    enddo

    if (allocated(term0)) deallocate(term0)

  endif

end subroutine chain_rule_product


!################################################################################


subroutine chain_rule_product_ind(nvar, nterms, terms, nelem_prod, nterms_chr, ind_chr)

  integer(ik), intent(in)               :: nvar
  integer(ik), intent(in)               :: nterms
  integer(ik), intent(in)               :: terms(nvar,nterms)
  integer(ik), intent(in)               :: nelem_prod
  integer(ik), intent(out)              :: nterms_chr(nterms)
  integer(ik), allocatable, intent(out) :: ind_chr(:,:,:)

  integer(ik) :: iterm, max_nelem, info, jterm, iprod, ipos
  integer(ik), allocatable :: terms2(:,:,:)

  ! estimate max number of elements in the chain rule expression among all terms

  do iterm=1, nterms
    call chain_rule_product(nvar, terms(1:nvar,iterm), nelem_prod, nterms_chr(iterm))
  enddo
  max_nelem = maxval(nterms_chr(1:nterms))

  ! allocate arrays to keep chain rule expressions

  if (allocated(ind_chr)) deallocate(ind_chr)
  allocate(ind_chr(nelem_prod,max_nelem,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,3(1x,i8))') &!
    'chain_rule_product_ind error: failed to allocate ind_chr(nelem_prod,max_nelem,nterms)', &!
    'nelem_prod, max_nelem, nterms =', nelem_prod, max_nelem, nterms
    stop
  endif

  ! generate chain rule expressions

  do iterm=1, nterms

    call chain_rule_product(nvar, terms(1:nvar,iterm), nelem_prod, nterms_chr(iterm), terms2)

    do jterm=1, nterms_chr(iterm)
      do iprod=1, nelem_prod
        ipos = index_iarr1(terms2(:,iprod,jterm), terms)
        if (ipos==0) then
          write(out, '(/a,1x,<nvar>(1x,i3),1x,a)') &!
          'chain_rule_product_ind error: failed to index expansion term = (', terms2(:,iprod,jterm), ')'
          stop
        else
          ind_chr(iprod,jterm,iterm) = ipos
        endif
      enddo
    enddo

  enddo

  if (allocated(terms2)) deallocate(terms2)

end subroutine chain_rule_product_ind


!################################################################################


subroutine nmode_distribute(nmodes, nterms, terms, ncomb, comb, nterms_n, terms_n)

  integer(ik), intent(in) :: nmodes
  integer(ik), intent(in) :: nterms
  integer(ik), intent(in) :: terms(nmodes,nterms)
  integer(ik), intent(out) :: ncomb
  integer(ik), allocatable, intent(out) :: comb(:,:)
  integer(ik), allocatable, intent(out) :: nterms_n(:)
  integer(ik), allocatable, intent(out) :: terms_n(:,:,:)

  integer(ik), parameter :: ncomb_incr = 10, nterms_incr = 1000
  !integer(ik), parameter :: ncomb_incr = 10, nterms_incr = 1000
  integer(ik) :: iterm, imode, icomb, ind(nmodes), n, ipos, ncomb_, nterms_, info, maxnterms
  integer(ik), allocatable :: tmp_comb(:,:), tmp_nterms(:), tmp_terms(:,:,:)

  if (allocated(comb)) deallocate(comb)
  if (allocated(nterms_n)) deallocate(nterms_n)
  if (allocated(terms_n)) deallocate(terms_n)
  allocate(comb(nmodes,1), nterms_n(1), terms_n(nmodes,1,1), stat=info)
  if (info/=0) then
    write(out, '(/a/a,1(1x,i6))') &!
    'nmode_distribute error: failed to allocate comb(nmodes,1), nterms_n(1), terms_n(nmodes,1,1)', &!
    'nmodes =', nmodes
    stop
  endif
  comb = 0
  nterms_n = 0
  terms_n = 0

  ncomb = 0

  ! loop over derivative terms
  do iterm=1, nterms

    ! determine combination of modes - indices of modes wrt which term(1:nmodes,iterm) derivative is taken
    n = 0
    ind = 0
    do imode=1, nmodes
      if (terms(imode,iterm)>0) then
        n = n + 1
        ind(n) = imode
      endif
    enddo

    ! check if current combination of modes was already added
    icomb = index_iarr1(ind(1:n), comb(1:n,1:ncomb))

    ! if not - add new combination
    if (icomb==0) then
      ncomb = ncomb + 1
      ! increase sizes of arrays (by ncomb_incr) if max number of combinations is exceeded
      if (ncomb>size(comb,dim=2)) then
        ncomb_ = size(comb,dim=2) + ncomb_incr
        nterms_ = size(terms_n,dim=2)
        call move_alloc(comb, tmp_comb)
        call move_alloc(nterms_n, tmp_nterms)
        call move_alloc(terms_n, tmp_terms)
        allocate(comb(nmodes,ncomb_), nterms_n(ncomb_), terms_n(nmodes,nterms_,ncomb_), stat=info)
        if (info/=0) then
          write(out, '(/a/a,3(1x,i6))') &!
          'nmode_distribute error: failed to allocate comb(nmodes,ncomb_), nterms_n(ncomb_), terms_n(nmodes,nterms_,ncomb_)', &!
          'nmodes, ncomb_, nterms_ =', nmodes, ncomb_, nterms_
          stop
        endif
        comb = 0
        nterms_n = 0
        terms_n = 0
        ncomb_ = size(tmp_comb,dim=2)
        comb(1:nmodes,1:ncomb_) = tmp_comb
        nterms_n(1:ncomb_) = tmp_nterms
        terms_n(1:nmodes,1:nterms_,1:ncomb_) = tmp_terms
        deallocate(tmp_comb)
        deallocate(tmp_nterms)
        deallocate(tmp_terms)
      endif
      ! add new combination of modes
      comb(1:n,ncomb) = ind(1:n)
      icomb = ncomb
    endif

    ! check if current derivative term was already added
    ipos = index_iarr1(terms(1:nmodes,iterm), terms_n(1:nmodes,1:nterms_n(icomb),icomb))

    ! if not - add new derivative term for icomb combination
    if (ipos==0) then
      nterms_n(icomb) = nterms_n(icomb) + 1
      ! increase size of terms_n array (by nterms_incr) if max number of terms is exceeded
      if (nterms_n(icomb)>size(terms_n,dim=2)) then
        ncomb_ = size(terms_n,dim=3)
        nterms_ = size(terms_n,dim=2) + nterms_incr
        call move_alloc(terms_n, tmp_terms)
        allocate(terms_n(nmodes,nterms_,ncomb_), stat=info)
        if (info/=0) then
          write(out, '(/a/a,3(1x,i6))') &!
          'nmode_distribute error: failed to allocate terms_n(nmodes,nterms_,ncomb_)', &!
          'nmodes, nterms_, ncomb_', nmodes, nterms_, ncomb_
          stop
        endif
        terms_n = 0
        nterms_ = size(tmp_terms,dim=2)
        terms_n(1:nmodes,1:nterms_,1:ncomb_) = tmp_terms
        deallocate(tmp_terms)
      endif
      ! add new expansion term
      terms_n(1:nmodes,nterms_n(icomb),icomb) = terms(1:nmodes,iterm)
    endif
  enddo

  ! free unused parts of arrays
  maxnterms = maxval(nterms_n)
  if (size(terms_n,dim=2)>maxnterms) then
    call move_alloc(comb, tmp_comb)
    call move_alloc(nterms_n, tmp_nterms)
    call move_alloc(terms_n, tmp_terms)
    allocate(comb(nmodes,ncomb), nterms_n(ncomb), terms_n(nmodes,maxnterms,ncomb), stat=info)
    if (info/=0) then
      write(out, '(/a/a,3(1x,i6))') &
      'nmode_distribute error: failed to allocate comb(nmodes,ncomb), nterms_n(ncomb), terms_n(nmodes,maxnterms,ncomb)', &!
      'nmodes, maxnterms, ncomb', nmodes, maxnterms, ncomb
      stop
    endif
    comb(1:nmodes,1:ncomb) = tmp_comb(1:nmodes,1:ncomb)
    nterms_n(1:ncomb) = tmp_nterms(1:ncomb)
    terms_n(1:nmodes,1:maxnterms,1:ncomb) = tmp_terms(1:nmodes,1:maxnterms,1:ncomb)
    deallocate(tmp_comb)
    deallocate(tmp_nterms)
    deallocate(tmp_terms)
  endif

end subroutine nmode_distribute


!################################################################################


subroutine nmode_expansion(nmodes, nmax, ndeg, deg, mindeg, maxdeg, nterms, terms)

  integer(ik), intent(in)               :: nmodes
  integer(ik), intent(in)               :: nmax
  integer(ik), intent(in)               :: ndeg(nmodes)
  integer(ik), intent(in)               :: deg(maxval(ndeg),nmodes)
  integer(ik), intent(in)               :: mindeg(nmax)
  integer(ik), intent(in)               :: maxdeg(nmax)
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)

  integer(ik) :: imode, info, iterm1, iterm2, nt, n, ncomb, icomb, nd(nmodes), d(maxval(ndeg),nmodes)
  integer(ik), allocatable :: comb(:,:), t(:,:), tmp(:,:)

  if (allocated(terms)) deallocate(terms)

  ! 0-order derivative term
  nterms = 1
  allocate(terms(nmodes,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') 'nmode_expansion error: failed to allocate terms(nmodes,nterms)', &!
    'nmodes, nterms =', nmodes, nterms
    stop
  endif
  iterm1 = 1
  iterm2 = 1
  terms(:,1) = 0

  do n=1, nmax

    ! generate all possible combinations of modes for current N
    call combinations(nmodes, (/(imode,imode=1,nmodes)/), n, ncomb, comb)

    do icomb=1, ncomb

      ! generate list of expansion degrees for each mode in current combination
      nd = 1
      d = 0
      do imode=1, n
        nd(comb(imode,icomb)) = ndeg(comb(imode,icomb))
        d(:,comb(imode,icomb)) = deg(:,comb(imode,icomb))
      enddo

      ! generate expansion terms for current combination of modes
      call expansion_terms(nmodes, nd(1:nmodes), d(1:maxval(nd),1:nmodes), mindeg(n), maxdeg(n), nt, t)

      iterm1 = iterm2 + 1
      iterm2 = iterm1 + nt - 1

      ! increase size of terms array if max number of terms is exceeded
      if (iterm2>size(terms,dim=2)) then
        nt = size(terms,dim=2) + nt
        call move_alloc(terms, tmp)
        allocate(terms(nmodes,nt), stat=info)
        if (info/=0) then
          write(out, '(/a/a,2(1x,i6))') 'nmode_expansion error: failed to allocate terms(nmodes,nt)', &!
          'nmodes, nt =', nmodes, nt
          stop
        endif
        terms(1:nmodes,1:size(tmp,dim=2)) = tmp
        deallocate(tmp)
      endif

      ! add derivative terms
      terms(1:nmodes,iterm1:iterm2) = t(1:nmodes,1:nt)
      nterms = iterm2

    enddo !icomb
  enddo !n

  ! free unused part of terms array
  call move_alloc(terms, tmp)
  allocate(terms(nmodes,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') 'nmode_expansion error: failed to allocate terms(nmodes,nterms)', &!
    'nmodes, nterms =', nmodes, nterms
    stop
  endif
  terms(1:nmodes,1:nterms) = tmp(1:nmodes,1:nterms)
  deallocate(tmp)

  if (allocated(t)) deallocate(t)
  if (allocated(comb)) deallocate(comb)

end subroutine nmode_expansion


!################################################################################


subroutine delete_same_terms(nmodes, nterms, terms)

  integer(ik), intent(in)                 :: nmodes
  integer(ik), intent(inout)              :: nterms
  integer(ik), intent(inout), allocatable :: terms(:,:)

  integer(ik) :: iterm, ipos, info, nterms_uniq
  integer(ik), allocatable :: terms_uniq(:,:)

  allocate(terms_uniq(nmodes,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'delete_same_terms error: failed to allocate terms_uniq(nmodes,nterms)', &
    'nmodes, nterms =', nmodes, nterms
    stop
  endif
  terms_uniq = 0
  nterms_uniq = 0

  do iterm=1, nterms
    if (nterms_uniq==0) then
      ipos = 0
    else
      ipos = index_iarr1(terms(1:nmodes,iterm), terms_uniq(1:nmodes,1:nterms_uniq))
    endif
    if (ipos==0) then
      nterms_uniq = nterms_uniq + 1
      terms_uniq(1:nmodes,nterms_uniq) = terms(1:nmodes,iterm)
    endif
  enddo

  deallocate(terms)
  allocate(terms(nmodes,nterms_uniq), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'delete_same_terms error: failed to allocate terms(nmodes,nterms_uniq)', &
    'nmodes, nterms_uniq =', nmodes, nterms_uniq
    stop
  endif
  terms(1:nmodes,1:nterms_uniq) = terms_uniq(1:nmodes,1:nterms_uniq)
  nterms = nterms_uniq
  deallocate(terms_uniq)

end subroutine delete_same_terms


!################################################################################


subroutine sort_terms_order(nmodes, nterms, terms)

  integer(ik), intent(in)    :: nmodes
  integer(ik), intent(in)    :: nterms
  integer(ik), intent(inout) :: terms(nmodes,nterms)

  integer(ik) :: iterm, terms_tmp(nmodes,nterms), ind_term(nterms)
  real(rk) :: rorder(nterms)

  terms_tmp = terms

  rorder(1:nterms) = (/( real(sum(terms(1:nmodes,iterm)),rk), iterm=1, nterms)/)
  ind_term(1:nterms) = (/(iterm, iterm=1, nterms)/)

  call sort2(nterms, rorder, ind_term)

  do iterm=1, nterms
    terms(1:nmodes,iterm) = terms_tmp(1:nmodes,ind_term(iterm))
  enddo

end subroutine sort_terms_order


!################################################################################


subroutine expansion_terms(n, ndeg, deg, mindeg, maxdeg, nterms, terms)

  integer(ik), intent(in)               :: n
  integer(ik), intent(in)               :: ndeg(n)
  integer(ik), intent(in)               :: deg(maxval(ndeg),n)
  integer(ik), intent(in)               :: mindeg
  integer(ik), intent(in)               :: maxdeg
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)

  integer(ik) :: ii(n), info

  ! estimate number of terms
  call cartesian_product(n, ndeg(1:n), deg(1:maxval(ndeg),1:n), .false., mindeg, maxdeg, 1, ii, nterms)

  ! allocate array to keep terms
  if (allocated(terms)) deallocate(terms)
  allocate(terms(n,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') 'expansion_terms error: failed to allocate terms(n,nterms)', 'n, nterms =', n, nterms
    stop
  endif
  terms = 0

  ! generate terms
  call cartesian_product(n, ndeg(1:n), deg(1:maxval(ndeg),1:n), .false., mindeg, maxdeg, 1, ii, nterms, terms)

end subroutine expansion_terms


!################################################################################


subroutine combinations(nmodes, ind_modes, n, ncomb, comb)

  integer(ik), intent(in)               :: nmodes
  integer(ik), intent(in)               :: ind_modes(nmodes)
  integer(ik), intent(in)               :: n
  integer(ik), intent(out)              :: ncomb
  integer(ik), allocatable, intent(out) :: comb(:,:)

  integer(ik) :: i, imode, ind(nmodes,n), minsum, maxsum, ii(n), info

  ind = 0
  do imode=1, n
    ind(1:nmodes,imode) = ind_modes(1:nmodes)
  enddo
  minsum = 0
  maxsum = sum(ind(:,1))

  ! estimate number of different combinations
  call cartesian_product(n, (/(nmodes, i=1, n)/), ind, .true., minsum, maxsum, 1, ii, ncomb)

  ! allocate array to keep combination indices
  if (allocated(comb)) deallocate(comb)
  allocate(comb(n,ncomb), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') 'combinations error: failed to allocate comb(n,ncomb)', 'n, ncomb =', n, ncomb
    stop
  endif
  comb = 0

  ! generate combination indices
  call cartesian_product(n, (/(nmodes, i=1, n)/), ind, .true., minsum, maxsum, 1, ii, ncomb, comb)

end subroutine combinations


!################################################################################


recursive subroutine cartesian_product(nvec, nelem, elem, sym, minsum, maxsum, ivec, ind, nterms, cprod)

  integer(ik), intent(in)              :: nvec
  integer(ik), intent(in)              :: nelem(nvec)
  integer(ik), intent(in)              :: elem(maxval(nelem),nvec)
  logical, intent(in)                  :: sym
  integer(ik), intent(in)              :: minsum
  integer(ik), intent(in)              :: maxsum
  integer(ik), intent(in)              :: ivec
  integer(ik), intent(inout)           :: ind(nvec)
  integer(ik), intent(inout)           :: nterms
  integer(ik), intent(inout), optional :: cprod(:,:)

  integer(ik) :: ielem, i
  integer(ik) :: sumind

  if (ivec==1) nterms = 0

  do ielem=1, nelem(ivec)
    ind(ivec) = elem(ielem,ivec)
    sumind = sum(ind(1:ivec))
    if (sumind>maxsum) cycle
    if (ivec>1.and.sym.and..not.all((/(ind(i-1)<ind(i), i=2, ivec)/))) cycle
    if (ivec==nvec) then
      if (sumind<minsum) cycle
      nterms = nterms + 1
      if (present(cprod)) cprod(:,nterms) = ind
    else
      call cartesian_product(nvec, nelem, elem, sym, minsum, maxsum, ivec+1, ind, nterms, cprod)
    endif

  enddo

end subroutine cartesian_product


!###############################################################################


! Computes inverse (Moore-Penrose) of complex matrix

subroutine matrix_inverse_cmplx(nelem, mat, invmat, info)

  integer(ik), intent(in) :: nelem
  complex(rk), intent(in) :: mat(nelem,nelem)
  complex(rk), intent(out) :: invmat(nelem,nelem)
  integer(ik), intent(out), optional :: info

  integer(ik) lwork, info_, i, j
  complex(rk) :: work1(1), matt(nelem,nelem), matu(nelem,nelem), matvt(nelem,nelem), mat_d(nelem,nelem), tmat(nelem,nelem)
  complex(rk), allocatable :: work(:)
  double precision :: rwork(5*nelem), tol, sv(nelem)

  tol = 1.0d-08!epsilon(1.0_rk)
  matt = mat

  lwork = -1
  call zgesvd('A', 'A', nelem, nelem, matt, nelem, sv, matu, nelem, matvt, nelem, work1, lwork, rwork, info_)
  lwork = int(work1(1), kind=ik)

  allocate(work(lwork), stat=info_)
  if (info_/=0) then
    write(out, '(/a,1x,i6)') 'matrix_inverse_cmplx error: failed to allocate workspace array for diagonalization, size =', lwork
    stop
  endif

  call zgesvd('A', 'A', nelem, nelem, matt, nelem, sv, matu, nelem, matvt, nelem, work, lwork, rwork, info_)

  if (info_/=0) then
    write(out, '(/a,1x,i6)') 'matrix_inverse_cmplx error: SVD failed, info =', info_
    stop
  endif

  if (present(info)) info = 0

  mat_d = 0.0
  do i=1, nelem
    if (sv(i)>=tol) then
      mat_d(i,i) = 1.0_rk/sv(i)
    else
      if (present(info)) then
        info = i
        mat_d(i,i) = 0.0
      else
        write(out, '(/a)') 'matrix_inverse_cmplx error: matrix is singular'
        write(out, '(a)') 'singular elements:'
        do j=1, nelem
          write(out, '(1x,i3,1x,f)') j, sv(j)
        enddo
        stop
      endif
    endif
  enddo

  tmat = matmul(mat_d, conjg(transpose(matu)))
  invmat = matmul(conjg(transpose(matvt)), tmat)

  deallocate(work)

end subroutine matrix_inverse_cmplx


!###############################################################################


! General Z-matrix internal-to-Cartesian coordinate transformation (ADF version).

subroutine fromlocal2cartesian_ADF(molec, pm, r, cartesian)

  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: pm
  type(adf_realq), intent(in)     :: r(:)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  type(adf_realq) :: x(molec%natoms,3), n1(3), n2(3), n3(3), cm(3), alpha, phi, v12(3), v23(3), f_t, beta, alpha3, cosa3, cosphi, rbond
  integer(ik) :: iangle, idihedral, iatom, j, p0, p1, p2, p3, p4, ix, zeta, kappa(3), ikappa, idelta
  real(ark) :: fm

  if (abs(pm)/=1) then
    write(out, '(/a,1x,i3)') 'fromlocal2cartesian_ADF error: illegal pm =', pm
    stop
  endif

  x = 0.0_ark

  ! #1 - at the reference coordinate of atom #1

  x(1,:) = 0.0_ark

  ! #2 - on the vector between reference atoms 1 and 2

  n1(:) = (/0.0_ark, 0.0_ark, 1.0_ark/)

  do ix=1, 3
    x(2,ix) = x(1,ix) + n1(ix) * r(1)
  enddo

  if (molec%natoms>2) then
    iatom = 3
    p1 = molec%zmat_connect(1,iatom)
    p2 = molec%zmat_connect(2,iatom)
    p3 = molec%zmat_connect(3,iatom)
    alpha = r(molec%nbonds+1)
    n2(1) = sin(alpha)  ;  n2(2) = 0.0_ark  ;  n2(3) = pm*cos(alpha)
    do ix=1, 3
      x(3,ix) = x(p1,ix) + n2(ix) * r(2)
    enddo
  endif

  iangle = 1
  idihedral = 0

  do iatom=4, molec%natoms

    J = molec%zmat_connect(4,iatom)

    ! it is special when we work with the third atom
    ! it will be distinguished by "J"

    rbond = r(iatom-1)

    iangle = iangle + 1
    alpha = r(molec%nbonds+iangle)

    ! connections

    p0 = iatom
    p1 = molec%zmat_connect(1,iatom)
    p2 = molec%zmat_connect(2,iatom)
    p3 = molec%zmat_connect(3,iatom)

    select case (J)

    case(-1,0)

      iangle = iangle + 1
      beta = r(molec%nbonds+iangle)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p1,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product_ADF(v23,v12)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product_ADF(n2,n3)

      cosa3 = sum(n2(:)*v23(:))/sqrt(sum(v23(:)**2))
      alpha3 = acos(cosa3)
      cosphi =  ( cos(beta)-cos(alpha)*cos(alpha3) )/( sin(alpha)*sin(alpha3) )
      phi = acos(cosphi)

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       + sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(1)

      idihedral = idihedral + 1
      phi = r(molec%nbonds+molec%nangles+idihedral)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p1,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product_ADF(v12,v23)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product_ADF(n2,n3)

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       - sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(-2,2)

      idihedral = idihedral + 1
      phi = r(molec%nbonds+molec%nangles+idihedral)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p2,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product_ADF(v23,v12)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product_ADF(n2,n3)

      if (J<0) n3 = -n3

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       - sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(4:100)

      iangle = iangle + 1
      beta = r(molec%nbonds+iangle)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p1,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product_ADF(v12,v23)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product_ADF(n3,n2)

      cosa3 = sum(n2(:)*v23(:))/sqrt(sum(v23(:)**2))
      alpha3 = acos(cosa3)
      cosphi =  ( cos(beta)-cos(alpha)*cos(alpha3) )/( sin(alpha)*sin(alpha3) )
      phi = acos(cosphi)

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       + sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(101)

      idihedral = idihedral + 1
      zeta = molec%zmat_connect(3,iatom)

      idelta = 0
      kappa(3) = zeta
      do ikappa = 1,3
        if (ikappa==zeta) cycle
        idelta = idelta + 1
        kappa(idelta) = ikappa
      enddo

      n1 = 0.0_ark
      n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
      n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
      n1(kappa(3)) = sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
      do ix=1, 3
        x(iatom,ix) = x(p1,ix)+rbond*n1(ix)
      enddo

    case(103)

      idihedral = idihedral + 1
      zeta = molec%zmat_connect(3,iatom)

      idelta = 0
      kappa(3) = zeta
      do ikappa = 1,3
        if (ikappa==zeta) cycle
        idelta = idelta + 1
        kappa(idelta) = ikappa
      enddo

      n1 = 0.0_ark
      n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
      n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
      n1(kappa(3)) = sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
      do ix=1, 3
        x(iatom,ix) = x(p1,ix)+rbond*n1(ix)
      enddo
    end select

  enddo ! iatom

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:molec%natoms))

  do ix=1, 3
    cm(ix) = sum( x(1:molec%natoms,ix) * molec%masses(1:molec%natoms) ) * fm
  enddo
  do iatom=1, molec%natoms
    cartesian(iatom,1:3) = x(iatom,1:3) - cm(1:3)
  enddo

end subroutine fromlocal2cartesian_ADF


!###############################################################################


! General Z-matrix internal-to-Cartesian coordinate transformation.

subroutine fromlocal2cartesian(molec, pm, r, cartesian)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: pm
  real(ark), intent(in)     :: r(:)
  real(ark), intent(out)    :: cartesian(molec%natoms,3)

  real(ark) :: x(molec%natoms,3), n1(3), n2(3), n3(3), cm(3), alpha, phi, v12(3), v23(3), f_t, beta, alpha3, cosa3, cosphi, rbond
  integer(ik) :: iangle, idihedral, iatom, j, p0, p1, p2, p3, p4, ix, zeta, kappa(3), ikappa, idelta
  real(ark) :: fm

  if (abs(pm)/=1) then
    write(out, '(/a,1x,i3)') 'fromlocal2cartesian error: illegal pm =', pm
    stop
  endif

  x = 0.0_ark

  ! #1 - at the reference coordinate of atom #1

  x(1,:) = 0.0_ark

  ! #2 - on the vector between reference atoms 1 and 2

  n1(:) = (/0.0_ark, 0.0_ark, 1.0_ark/)

  do ix=1, 3
    x(2,ix) = x(1,ix) + n1(ix) * r(1)
  enddo

  if (molec%natoms>2) then
    iatom = 3
    p1 = molec%zmat_connect(1,iatom)
    p2 = molec%zmat_connect(2,iatom)
    p3 = molec%zmat_connect(3,iatom)
    alpha = r(molec%nbonds+1)
    n2(1) = sin(alpha)  ;  n2(2) = 0.0_ark  ;  n2(3) =pm*cos(alpha)
    do ix=1, 3
      x(3,ix) = x(p1,ix) + n2(ix) * r(2)
    enddo
  endif

  iangle = 1
  idihedral = 0

  do iatom=4, molec%natoms

    J = molec%zmat_connect(4,iatom)

    ! it is special when we work with the third atom
    ! it will be distinguished by "J"

    rbond = r(iatom-1)

    iangle = iangle + 1
    alpha = r(molec%nbonds+iangle)

    ! connections

    p0 = iatom
    p1 = molec%zmat_connect(1,iatom)
    p2 = molec%zmat_connect(2,iatom)
    p3 = molec%zmat_connect(3,iatom)

    select case (J)

    case(-1,0)

      iangle = iangle + 1
      beta = r(molec%nbonds+iangle)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p1,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product(v23,v12)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product(n2,n3)

      cosa3 = sum(n2(:)*v23(:))/sqrt(sum(v23(:)**2))
      alpha3 = acos(cosa3)
      cosphi =  ( cos(beta)-cos(alpha)*cos(alpha3) )/( sin(alpha)*sin(alpha3) )
      phi = acos(cosphi)

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       + sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(1)

      idihedral = idihedral + 1
      phi = r(molec%nbonds+molec%nangles+idihedral)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p1,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product(v12,v23)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product(n2,n3)

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       - sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(-2,2)

      idihedral = idihedral + 1
      phi = r(molec%nbonds+molec%nangles+idihedral)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p2,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product(v23,v12)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product(n2,n3)

      if (J<0) n3 = -n3

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       - sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(4:100)

      iangle = iangle + 1
      beta = r(molec%nbonds+iangle)

      v12 = x(p2,:) - x(p1,:)
      v23 = x(p3,:) - x(p1,:)
      n2 = v12/sqrt(sum(v12(:)**2))
      n3 = vector_product(v12,v23)

      f_t = sqrt(sum(n3(:)**2))
      n3 = n3/f_t
      n1 = vector_product(n3,n2)

      cosa3 = sum(n2(:)*v23(:))/sqrt(sum(v23(:)**2))
      alpha3 = acos(cosa3)
      cosphi =  ( cos(beta)-cos(alpha)*cos(alpha3) )/( sin(alpha)*sin(alpha3) )
      phi = acos(cosphi)

      do ix=1, 3
        x(iatom,ix) = x(p1,ix) + rbond*( cos(alpha)*n2(ix) &!
                                       + sin(alpha)*cos(phi)*n1(ix) &!
                                       + sin(alpha)*sin(phi)*n3(ix) )
      enddo

    case(101)

      idihedral = idihedral + 1
      zeta = molec%zmat_connect(3,iatom)

      idelta = 0
      kappa(3) = zeta
      do ikappa = 1,3
        if (ikappa==zeta) cycle
        idelta = idelta + 1
        kappa(idelta) = ikappa
      enddo

      n1 = 0
      n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
      n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral)
      n1(kappa(3)) = sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
      do ix=1, 3
        x(iatom,ix) = x(p1,ix)+rbond*n1(ix)
      enddo

    case(103)

      idihedral = idihedral + 1
      zeta = molec%zmat_connect(3,iatom)

      idelta = 0
      kappa(3) = zeta
      do ikappa = 1,3
        if (ikappa==zeta) cycle
        idelta = idelta + 1
        kappa(idelta) = ikappa
      enddo

      n1 = 0
      n1(kappa(1)) = r(molec%Nbonds+molec%Nangles+idihedral)
      n1(kappa(2)) = r(molec%Nbonds+molec%Nangles+idihedral+1)
      n1(kappa(3)) =sqrt(1.0_ark-(n1(kappa(1))**2+n1(kappa(2))**2))
      do ix=1, 3
        x(iatom,ix) = x(p1,ix)+rbond*n1(ix)
      enddo

    end select

  enddo ! iatom

  ! shift to nuclear center of mass

  fm = 1.0_ark / sum(molec%masses(1:molec%natoms))

  do ix=1, 3
    cm(ix) = sum( x(1:molec%natoms,ix) * molec%masses(1:molec%natoms) ) * fm
  enddo
  do iatom=1, molec%natoms
    cartesian(iatom,1:3) = x(iatom,1:3) - cm(1:3)
  enddo

end subroutine fromlocal2cartesian


!###############################################################################


function vector_product_ADF(v1,v2) result (v)
  use adf
  implicit none
  type(adf_realq), intent(in) :: v1(3), v2(3)
  type(adf_realq) :: v(3)
  v(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v(3) = v1(1)*v2(2)-v1(2)*v2(1)
end function vector_product_ADF


!###############################################################################


function vector_product(v1,v2) result (v)
  real(ark), intent(in) :: v1(3), v2(3)
  real(ark) :: v(3)
  v(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v(3) = v1(1)*v2(2)-v1(2)*v2(1)
end function vector_product


!###############################################################################
