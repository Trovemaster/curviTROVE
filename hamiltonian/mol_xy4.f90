
! Given equations for five symmetry-adapted coordinates S as functions of five bending angles ALPHA
! this subroutine computes "max_order"-order Taylor series expansions of ALPHA in terms of S.
! The resulting expansion terms and coefficients are stored in global variables "xy4_terms_alpha2sym" and "xy4_coefs_alpha2sym",
! which are required later at a stage of internal-to-Cartesian transformation (see, for example, "internal_to_cartesian_xy4_rsymalpha_ADF" function).

subroutine xy4_convert_symalpha_to_alpha(molec, max_order, s_eq)

  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in) :: max_order
  real(ark), intent(in) :: s_eq(5)

  integer(ik) :: nmodes, ndeg(5), deg(100,5), imode, nterms, info, term0(5), ipos, i, iorder, nterms1, iterm, verbose
  integer(ik), allocatable :: terms(:,:), terms1(:,:)
  real(ark) :: eq_guess(5), alpha_eq(5), alpha34_eq, fac
  real(ark), allocatable :: coefs_smat(:,:,:), coefs_alpha(:,:,:)
  type(adf_realq) :: alpha_adf(5), alpha12, alpha13, alpha14, alpha23, alpha24, cosbeta, beta312, beta412, cosa34, alpha34, s(5)

  nmodes = 5
  verbose = molec%expkeo_verbose

  if (verbose>=6) write(out, '(/a)') 'xy4_convert_symalpha_to_alpha/start'

  ! find equilibrium for ALPHA

  if (verbose>=6) then
    write(out, '(1x,a,<nmodes>(/1x,i3,1x,es16.8))') 'equilibrium values of S-coordinates:', (imode, s_eq(imode), imode=1, nmodes)
  endif

  eq_guess = 109.5_ark*real(pi,ark)/180.0_ark

  call from_sym2alphaII(molec, s_eq(1:nmodes), eq_guess, alpha_eq(1:nmodes), alpha34_eq)

  if (verbose>=6) then
    write(out, '(1x,a,<nmodes>(/1x,i3,1x,es16.8))') 'equilibrium values of ALPHA-coordinates:', (imode, alpha_eq(imode), imode=1, nmodes)
  endif

  ! initialize expansion terms

  ndeg = 1
  deg = 0
  do imode=1, nmodes
    ndeg(imode) = max_order + 1
    deg(1:ndeg(imode),imode) = (/(iorder, iorder=0, max_order)/)
  enddo
  call expansion_terms(nmodes, ndeg(1:nmodes), deg(1:maxval(ndeg),1:nmodes), 0, max_order, nterms, terms) ! array "terms" is (re-)allocated inside "expansion_terms"


  ! build expansions of S in terms of ALPHA at ALPHA=ALPHA_EQ using ADF

  call adf_init(nmodes, alpha_adf(1:nmodes), alpha_eq(1:nmodes), nterms, terms)

  alpha12 = alpha_adf(1)
  alpha13 = alpha_adf(2)
  alpha14 = alpha_adf(3)
  alpha23 = alpha_adf(4)
  alpha24 = alpha_adf(5)

  cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
  beta312 = acos(cosbeta)

  cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
  beta412 = acos(cosbeta)

  cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
  alpha34 = acos(cosa34)

  s(1) = (2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
  s(2) = (alpha13-alpha14-alpha23+alpha24)*0.5_ark
  s(3) = (alpha24-alpha13)/sqrt(2.0_ark)
  s(4) = (alpha23-alpha14)/sqrt(2.0_ark)
  s(5) = (alpha34-alpha12)/sqrt(2.0_ark)


  if (verbose>=6) then
    write(out, '(1x,a)') 'derivatives of S wrt ALPHA'
    do iterm=1, adf_nterms
      write(out, '(1x,i3,5x,<nmodes>(1x,i3),5x,<nmodes>(1x,es16.8))') iterm, adf_terms(1:nmodes,iterm), (s(i)%d(iterm), i=1, nmodes)
    enddo
  endif


  ! build matrix of first-order derivatives SMAT(k,l) = dS(k)/dALPHA(l) as expansions in terms of ALPHA coordinates

  allocate(coefs_smat(nmodes,nmodes,nterms), coefs_alpha(nmodes,nmodes,nterms), terms1(nmodes,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'xy4_convert_symalpha_to_alpha error: failed to allocate coefs_smat(nmodes,nmodes,nterms), coefs_alpha(nmodes,nmodes,nterms), terms1(nmodes,nterms)', &!
    'nmodes, nterms =', nmodes, nterms
    stop
  endif

  coefs_smat = 0.0
  nterms1 = 0
  terms1 = 0

  do iterm=1, nterms
    do imode=1, nmodes
      term0 = terms(1:nmodes,iterm)
      if (term0(imode)>0) then
        term0(imode) = term0(imode) - 1
        if (nterms1==0) then
          ipos = 0
        else
          ipos = index_iarr1(term0, terms1(1:nmodes,1:nterms1))
        endif
        if (ipos<=0) then
          nterms1 = nterms1 + 1
          terms1(:,nterms1) = term0
          ipos = nterms1
        endif
        coefs_smat(1:nmodes,imode,ipos) = s(1:nmodes)%d(iterm)
      endif
    enddo
  enddo

  call adf_finalize


  ! compute expansion of ALPHA' = (S')^{-1}

  coefs_alpha = 0.0

  call deriv_invmat(nmodes, nterms1, terms1(:,1:nterms1), nmodes, coefs_smat(:,:,1:nterms1), coefs_alpha(:,:,1:nterms1))


  ! store expansion coefficients in global arrays

  do iterm=1, nterms1
    fac = 1.0_ark
    do imode=1, nmodes
      if (imode>0) then
        !fac = fac / qfactorial(terms1(imode,iterm)+1)
        fac = fac / qfactorial(terms1(imode,iterm))
      endif
    enddo
   ! coefs_alpha(:,:,iterm) = coefs_alpha(:,:,iterm) * fac
  enddo

  if (allocated(xy4_terms_alpha2sym)) deallocate(xy4_terms_alpha2sym)
  if (allocated(xy4_coefs_alpha2sym)) deallocate(xy4_coefs_alpha2sym)
  if (allocated(xy4_coefs_alpha2sym0)) deallocate(xy4_coefs_alpha2sym0)
  allocate(xy4_terms_alpha2sym(nmodes,nterms1), xy4_coefs_alpha2sym(nmodes,nmodes,nterms1), xy4_coefs_alpha2sym0(nmodes), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'xy4_convert_symalpha_to_alpha error: failed to allocate xy4_terms_alpha2sym(nmodes,nterms1), xy4_coefs_alpha2sym(nmodes,nmodes,nterms1)', &!
    'nmodes, nterms1 =', nmodes, nterms1
    stop
  endif

  xy4_terms_alpha2sym(1:nmodes,1:nterms1) = terms1(1:nmodes,1:nterms1)
  xy4_coefs_alpha2sym(1:nmodes,1:nmodes,1:nterms1) = coefs_alpha(1:nmodes,1:nmodes,1:nterms1)
  xy4_coefs_alpha2sym0(1:nmodes) = alpha_eq(1:nmodes)


  deallocate(terms, terms1, coefs_smat, coefs_alpha)

  if (verbose>=6) write(out, '(a/)') 'xy4_convert_symalpha_to_alpha/done'

end subroutine xy4_convert_symalpha_to_alpha


!################################################################################


! Internal-to-Cartesian coordinate transformation for XY4 molecule.
! Internal coordinates: r1, r2, r3, r4, s1, s2, s3, s4, s5
!                       where r1, r2, r3, r4 - r(X-Y1), r(X-Y2), r(X-Y3), r(X-Y4) bond stretches,
!                       s1 = (2*alpha12-alpha13-alpha14-alpha23-alpha24+2*alpha34)/sqrt(12),
!                       s2 = (alpha13-alpha14-alpha23+alpha24)/2,
!                       s3 = (alpha24-alpha13)/sqrt(2),
!                       s4 = (alpha23-alpha14)/sqrt(2),
!                       s5 = (alpha34-alpha12)/sqrt(2),
!                       and alpha12, alpha13, alpha14, alpha23, alpha24, alpha34, are bond angles between respective Y atoms.

subroutine internal_to_cartesian_xy4_rsymalpha(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: pm=1
  real(ark) :: fm
  real(ark) ::  r1, r2, r3, r4, local(10), eq(5), cm(3), alpha(6), sqrt3, re14

  local(1) = internal(1) ! r_{XY1}
  local(2) = internal(2) ! r_{XY2}
  local(3) = internal(3) ! r_{XY3}
  local(4) = internal(4) ! r_{XY4}

  ! compute alpha-coordinates from symmetrized s-coordinates

  ! initial values for alphas
  eq(1:5) = 109.5_ark*real(pi,ark)/180.0_ark

  call from_sym2alphaII(molec,internal(5:9),eq(1:5),alpha(1:5),alpha(6))
  ! order of alphas on output: 1=alpha12, 2=alpha13, 3=alpha14, 4=alpha23, 5=alpha24, 6=alpha34

  local(5) = alpha(1)
  local(6) = alpha(2)
  local(7) = alpha(4)
  local(8) = alpha(3)
  local(9) = alpha(5)

  ! compute Z-matrix Cartesian coordinates

  call fromlocal2cartesian(molec, pm, local(1:9), cartesian)

end subroutine internal_to_cartesian_xy4_rsymalpha


!################################################################################


! ADF-adapted version of "internal_to_cartesian_xy4_rsymalpha" subroutine.

subroutine internal_to_cartesian_xy4_rsymalpha_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: pm=1
  real(ark) :: fm, eq(5)
  type(adf_realq) :: r1, r2, r3, r4, alpha(6), cm(3), local(10)

  local(1) = internal(1) ! r_{XY1}
  local(2) = internal(2) ! r_{XY2}
  local(3) = internal(3) ! r_{XY3}
  local(4) = internal(4) ! r_{XY4}

  ! compute alpha-coordinates from symmetrized s-coordinates

  ! initial values for alphas
  eq(1:5) = 109.5_ark*real(pi,ark)/180.0_ark

  call from_sym2alphaII_ADF(molec,internal(5:9),eq(1:5),alpha(1:5),alpha(6))
  ! order of alphas on output: 1=alpha12, 2=alpha13, 3=alpha14, 4=alpha23, 5=alpha24, 6=alpha34

  local(5) = alpha(1)
  local(6) = alpha(2)
  local(7) = alpha(4)
  local(8) = alpha(3)
  local(9) = alpha(5)

  ! compute Z-matrix Cartesian coordinates

  call fromlocal2cartesian_ADF(molec, pm, local(1:9), cartesian)

end subroutine internal_to_cartesian_xy4_rsymalpha_ADF


!################################################################################


subroutine internal_to_cartesian_xy4_equilibrium_td(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  real(ark) :: sqrt3, re14

  re14 = internal(1)
  sqrt3 = 1.0_ark/sqrt(3.0_ark)

  cartesian(1,1:3) = (/ 0.0_ark, 0.0_ark, 0.0_ark /)
  cartesian(2,1:3) = (/ -re14*sqrt3, re14*sqrt3, re14*sqrt3 /)
  cartesian(3,1:3) = (/ -re14*sqrt3, -re14*sqrt3, -re14*sqrt3 /)
  cartesian(4,1:3) = (/ re14*sqrt3, re14*sqrt3, -re14*sqrt3 /)
  cartesian(5,1:3) = (/ re14*sqrt3, -re14*sqrt3, re14*sqrt3 /)

end subroutine internal_to_cartesian_xy4_equilibrium_td


!################################################################################


! Returns six valence-angle bending coordinates of XY4 molecule for given five symmetry-adapted bending coordinates
! slightly modified version of mol_xy4.f90/from_sym2alphaII

subroutine from_sym2alphaII(molec,s,local_eq,local,alpha34)
  implicit none

  type(HM_molec_type), intent(in) :: molec
  real(ark),intent(in)  :: s(5)
  real(ark), intent(in) :: local_eq(5)
  real(ark),intent(out) :: local(5),alpha34

  real(ark) :: eps(5),rjacob(5,5),s_r(5),s_l(5),am(5,5),ainv(5,5),bm(5),cm(5),a(5,5),b(5,1),alpha34_
  real(ark) :: stadev,stadev_best,h
  integer(ik) :: iter,itmax,k,i,j,iterm

  ! initial values
  local(:) = local_eq(5)

  iter = 0
  itmax = 30
  stadev =  1.e10
  stadev_best = sqrt(epsilon(1.0_rk))*0.001_ark

  do while( iter<itmax .and. stadev>stadev_best )

    iter = iter + 1

    ! calculate function

    call calc_sym_from_alpha(local,s_r,alpha34)

    eps(:) = s(:)-s_r(:)

    ! calculate gradient

    do k = 1,5

      h = 0.001_ark!*abs(local(k))
      !if (h<1e-12) h = 1e-7

      local(k) = local(k) + h

      call calc_sym_from_alpha(local,s_r,alpha34_)

      local(k) = local(k) - h - h

      call calc_sym_from_alpha(local,s_l,alpha34_)

      rjacob(:,k) = (s_r(:)-s_l(:))/h*0.5_ark

      local(k) = local(k) + h

    enddo

    ! construct a set of linear equations A x = B

    ! form A matrix
    do i=1,5
      do j=1,i
        am(i,j)=sum(rjacob(:,j)*rjacob(:,i))
        am(j,i)=am(i,j)
      enddo
    enddo

    ! form B matrix
    do i=1,5
      bm(i)=sum(eps(:)*rjacob(:,i))
    enddo

    ! solve the system of linear equations

    !call matrix_inverse_ark(5, am, ainv)
    call matrix_inverse2_ark(5, am, ainv)
    do i=1, 5
      cm(i) = dot_product(ainv(i,:), bm)
    enddo

    local(:) = local(:) + cm(:)

    stadev = sqrt(sum(eps(:)**2))/sqrt(5.0_ark)

  enddo

  if (iter==itmax) then
    write(out, '(/a,1x,i3,1x,a)') 'from_sym2alphaII error: could not find solution after', itmax, 'iterations'
    stop
  endif

  contains

  ! Returns five symmetry-adapted bending coordinates in terms of six valence-angle bending coordinates

  subroutine calc_sym_from_alpha(src,dst,alpha34)

    real(ark),intent(in)  :: src(5)
    real(ark),intent(out) :: dst(5),alpha34
    real(ark)             :: alpha12,alpha13,alpha23,alpha14,alpha24,cosbeta,beta312,beta412,cosa34

    alpha12 = src(1)
    alpha13 = src(2)
    alpha14 = src(3)
    alpha23 = src(4)
    alpha24 = src(5)

    cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
    beta312 = acos(cosbeta)

    cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
    beta412 = acos(cosbeta)

    cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
    alpha34 = acos(cosa34)

    dst(1)=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
    dst(2)=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
    dst(3)=(alpha24-alpha13)/sqrt(2.0_ark)
    dst(4)=(alpha23-alpha14)/sqrt(2.0_ark)
    dst(5)=(alpha34-alpha12)/sqrt(2.0_ark)

  end subroutine calc_sym_from_alpha

end subroutine from_sym2alphaII


!################################################################################


! Returns six valence-angle bending coordinates of XY4 molecule for given five symmetry-adapted bending coordinates (ADF version)
! slightly modified version of mol_xy4.f90/from_sym2alphaII

subroutine from_sym2alphaII_ADF(molec,s,local_eq,local,alpha34)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq),intent(in)  :: s(5)
  real(ark), intent(in) :: local_eq(5)
  type(adf_realq),intent(out) :: local(5),alpha34

  type(adf_realq) :: eps(5),rjacob(5,5),s_r(5),s_l(5),am(5,5),ainv(5,5),bm(5),cm(5),a(5,5),b(5,1),alpha34_,h
  real(ark) :: stadev(adf_nterms),stadev_best
  real(ark) :: coefs_am(5,5,adf_nterms),coefs_ainv(5,5,adf_nterms)
  integer(ik) :: iter,itmax,k,i,j,iterm

  ! initial values
  local(:) = local_eq(5)

  iter = 0
  itmax = 10
  stadev = 1.e10
  stadev_best = sqrt(epsilon(1.0_rk))*0.001_ark

  do while( iter<itmax )! .and. all(stadev(:)>stadev_best) )

    iter = iter + 1

    ! calculate function

    call calc_sym_from_alpha(local, s_r, alpha34)

    eps(:) = s(:) - s_r(:)

    ! calculate gradient

    do k = 1,5

      h = 0.001_ark!*abs(local(k)%v)
      !if (h%v<1e-12) h = real(1e-7,ark)

      local(k) = local(k) + h

      call calc_sym_from_alpha(local, s_r, alpha34_)

      local(k) = local(k) - h - h

      call calc_sym_from_alpha(local, s_l, alpha34_)

      rjacob(:,k) = ( s_r(:) - s_l(:) ) / h * 0.5_ark

      local(k) = local(k) + h

    enddo

    ! construct a set of linear equations A x = B

    ! form A matrix
    do i=1,5
      do j=1,i
        am(i,j)=sum(rjacob(:,j)*rjacob(:,i))
        am(j,i)=am(i,j)
      enddo
    enddo

    ! form B matrix
    do i=1,5
      bm(i)=sum(eps(:)*rjacob(:,i))
    enddo

    ! solve the system of linear equations

    do iterm=1, adf_nterms
      do i=1, 5
        do j=1, 5
          coefs_am(j,i,iterm) = am(j,i)%d(iterm)
        enddo
      enddo
    enddo
    call deriv_invmat(molec%nmodes, adf_nterms, adf_terms, 5, coefs_am, coefs_ainv)
    do i=1, 5
      do j=1, 5
        call adf_set_var(ainv(j,i), coefs_ainv(j,i,1:adf_nterms))
      enddo
    enddo
    do i=1, 5
      cm(i) = dot_product(ainv(i,:), bm)
    enddo

    local(:) = local(:) + cm(:)

    do iterm=1, adf_nterms
      stadev(iterm) = sqrt(sum(eps(:)%d(iterm)**2))/sqrt(5.0_ark)
    enddo

  enddo

  !if (iter==itmax) then
  !  write(out, '(/a,1x,i3,1x,a)') 'from_sym2alphaII_ADF error: could not find solution after', itmax, 'iterations'
  !  stop
  !endif

  contains

  ! Returns five symmetry-adapted bending coordinates in terms of six valence-angle bending coordinates (ADF version)

  subroutine calc_sym_from_alpha(src,dst,alpha34)

    type(adf_realq),intent(in)  :: src(5)
    type(adf_realq),intent(out) :: dst(5),alpha34
    type(adf_realq)             :: alpha12,alpha13,alpha23,alpha14,alpha24,cosbeta,beta312,beta412,cosa34

    alpha12 = src(1)
    alpha13 = src(2)
    alpha14 = src(3)
    alpha23 = src(4)
    alpha24 = src(5)

    cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
    beta312 = acos(cosbeta)

    cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
    beta412 = acos(cosbeta)

    cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
    alpha34 = acos(cosa34)

    dst(1)=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
    dst(2)=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
    dst(3)=(alpha24-alpha13)/sqrt(2.0_ark)
    dst(4)=(alpha23-alpha14)/sqrt(2.0_ark)
    dst(5)=(alpha34-alpha12)/sqrt(2.0_ark)

  end subroutine calc_sym_from_alpha

end subroutine from_sym2alphaII_ADF


!################################################################################


subroutine internal_to_cartesian_xy4_symbeta_tau(molec, internal, cartesian)

  type(HM_molec_type), intent(in) :: molec
  real(ark), intent(in)           :: internal(molec%nmodes)
  real(ark), intent(out)          :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  real(ark) ::  r1, r2, r3, r4, s1, s2, s3, rho1, rho2, beta13, beta23, beta24, beta14, cm(3)

  natoms = molec%natoms

  r1 = internal(1)
  r2 = internal(2)
  r3 = internal(3)
  r4 = internal(4)
  s1 = internal(5)   ! s1 = (beta13-beta14-beta23+beta24)/2
  s2 = internal(6)   ! s2 = (beta24-beta13)/sqrt(2)
  s3 = internal(7)   ! s3 = (beta23-beta14)/sqrt(2)
  rho1 = real(pi,ark)*0.5_ark - internal(8)
  rho2 = real(pi,ark)*0.5_ark + internal(9)

  beta13 = ( real(pi,ark) + s1 - s2 * sqrt(2.0_ark) )*0.5_ark
  beta23 = ( real(pi,ark) - s1 + s3 * sqrt(2.0_ark) )*0.5_ark
  beta24 = ( real(pi,ark) + s1 + s2 * sqrt(2.0_ark) )*0.5_ark
  beta14 = ( real(pi,ark) - s1 - s3 * sqrt(2.0_ark) )*0.5_ark

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  0.0_ark
  cartesian(2,2) = -r1*sin(rho1)
  cartesian(2,3) =  r1*cos(rho1)

  cartesian(4,1) =  r3*sin(rho2)*sin(beta13)
  cartesian(4,2) = -r3*sin(rho2)*cos(beta13)
  cartesian(4,3) =  r3*cos(rho2)

  cartesian(3,1) =  r2*sin(rho1)*sin(beta13+beta23)
  cartesian(3,2) = -r2*sin(rho1)*cos(beta13+beta23)
  cartesian(3,3) =  r2*cos(rho1)

  cartesian(5,1) = -r4*sin(rho2)*sin(beta14)
  cartesian(5,2) = -r4*sin(rho2)*cos(beta14)
  cartesian(5,3) =  r4*cos(rho2)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy4_symbeta_tau


!################################################################################


! ADF-adapted version of "internal_to_cartesian_xy4_symbeta_tau" subroutine.

subroutine internal_to_cartesian_xy4_symbeta_tau_ADF(molec, internal, cartesian)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: cartesian(molec%natoms,3)

  integer(ik) :: ix, iatom, natoms
  real(ark) :: fm
  type(adf_realq) ::  r1, r2, r3, r4, s1, s2, s3, rho1, rho2, beta13, beta23, beta24, beta14, cm(3)

  natoms = molec%natoms

  r1 = internal(1)
  r2 = internal(2)
  r3 = internal(3)
  r4 = internal(4)
  s1 = internal(5)   ! s1 = (beta13-beta14-beta23+beta24)/2
  s2 = internal(6)   ! s2 = (beta24-beta13)/sqrt(2)
  s3 = internal(7)   ! s3 = (beta23-beta14)/sqrt(2)
  rho1 = real(pi,ark)*0.5_ark - internal(8)
  rho2 = real(pi,ark)*0.5_ark + internal(9)

  beta13 = ( real(pi,ark) + s1 - s2 * sqrt(2.0_ark) )*0.5_ark
  beta23 = ( real(pi,ark) - s1 + s3 * sqrt(2.0_ark) )*0.5_ark
  beta24 = ( real(pi,ark) + s1 + s2 * sqrt(2.0_ark) )*0.5_ark
  beta14 = ( real(pi,ark) - s1 - s3 * sqrt(2.0_ark) )*0.5_ark

  cartesian(1,1) =  0.0_ark
  cartesian(1,2) =  0.0_ark
  cartesian(1,3) =  0.0_ark

  cartesian(2,1) =  0.0_ark
  cartesian(2,2) = -r1*sin(rho1)
  cartesian(2,3) =  r1*cos(rho1)

  cartesian(4,1) =  r3*sin(rho2)*sin(beta13)
  cartesian(4,2) = -r3*sin(rho2)*cos(beta13)
  cartesian(4,3) =  r3*cos(rho2)

  cartesian(3,1) =  r2*sin(rho1)*sin(beta13+beta23)
  cartesian(3,2) = -r2*sin(rho1)*cos(beta13+beta23)
  cartesian(3,3) =  r2*cos(rho1)

  cartesian(5,1) = -r4*sin(rho2)*sin(beta14)
  cartesian(5,2) = -r4*sin(rho2)*cos(beta14)
  cartesian(5,3) =  r4*cos(rho2)

  ! shift to nuclear centre of mass

  fm = 1.0_ark / sum(molec%masses(1:natoms))

  do ix=1, 3
    cm(ix) = sum( cartesian(1:natoms,ix) * molec%masses(1:natoms) ) * fm
  enddo
  do iatom=1, natoms
    cartesian(iatom,1:3) = cartesian(iatom,1:3) - cm(1:3)
  enddo

end subroutine internal_to_cartesian_xy4_symbeta_tau_ADF
