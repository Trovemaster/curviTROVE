module hamiltonian
use accuracy
implicit none


private
public HM_molec_type, HM_func_type, read_kinetic, expand_kinetic, expand_func, pointwise_kinetic, pointwise_potential, expand_kinetic_1d, expansion_terms
public delete_same_terms, index_iarr1


type HM_molec_type
  integer(ik)              :: natoms
  integer(ik)              :: nmodes
  real(ark), allocatable   :: masses(:)
  integer(ik)              :: nbonds, nangles
  integer(ik), allocatable :: zmat_connect(:,:)
  !
  character(cl)            :: deriv_method = 'ADF'
  character(cl)            :: coord_transform = 'NONE'
  character(cl)            :: frame_rotation = 'ECKART'
  character(cl)            :: reference_frame = 'PAS'
  character(cl)            :: pas_function = 'EXPKAPPA'
  integer(ik)              :: matexp_maxntaylor = 30
  integer(ik)              :: matexp_deg2 = 6
  real(ark)                :: matexp_accur = 1.0d-60
  integer(ik)              :: eckart_maxniter = 10
  real(ark)                :: eckart_accur = 1.0d-30
  integer(ik)              :: pas_maxniter = 100
  real(ark)                :: pas_accur = 1.0d-30!1.0d-60!epsilon(1.0d0)
  integer(ik)              :: expkeo_verbose = 6
  integer(ik)              :: expfunc_verbose = 3
  integer(ik)              :: expcart_verbose = 3
  !
  logical                  :: read_keo = .false.
  character(cl)            :: fname_gmat = 'keo_gmat.dat'
  character(cl)            :: fname_pseudo = 'keo_pseudo.dat'
  !
  integer(ik)              :: nonrigid_imode = 0
  logical                  :: vib_3n5 = .false.
  !
  character(cl), allocatable              :: sing_type(:)
  integer(ik)                :: sing_nmodes = 0  ! number of coordinates that sample singular point
  integer(ik), allocatable   :: sing_imode(:)    ! indexes of coordinates sampling singular point
  real(ark), allocatable     :: sing_x(:)        ! values of coordinates at singular point
  character(cl), allocatable :: sing_fx(:)
  real(ark)                  :: sing_xtol = 1.0d-06
  real(ark)                  :: sing_coef_tol = 1.0d-08
  integer(ik)                :: sing_npoints = 10
  real(ark), allocatable     :: sing_points(:,:)  ! extrapolation grid around singular point
  real(ark)                  :: sing_thresh = 1.0d-08
  integer(ik)                :: sing_nterms       ! extrapolation polynomial
  integer(ik), allocatable   :: sing_terms(:,:)
  real(ark)                    :: sing_zero_thresh = 1.0d-08
  logical, allocatable       :: sing_tvec(:,:)
  !
  integer(ik)                :: imode_sing
  logical, allocatable       :: ivec_sing(:)
  real(ark)                  :: zero_sing = 1.0d-10
  !
end type HM_molec_type


integer(ik), allocatable :: xy4_terms_alpha2sym(:,:)
real(ark), allocatable :: xy4_coefs_alpha2sym(:,:,:), xy4_coefs_alpha2sym0(:)




type HM_func_type
  character(cl)              :: func_type
  character(cl)              :: ref_poten_type
  character(cl), allocatable :: coord_type(:)
  integer(ik)                :: coord_nparams
  real(ark), allocatable     :: coord_params(:,:)
  integer(ik)                :: rank
  integer(ik), allocatable   :: nparams(:)
  integer(ik), allocatable   :: iparams(:,:,:)
  integer(ik), allocatable   :: ifit(:,:)
  real(ark), allocatable     :: params(:,:)
  logical                    :: rotate = .false.
end type HM_func_type


procedure(internal_to_cartesian), pointer     :: internal_to_cartesian_ptr => null()
procedure(internal_to_cartesian_ADF), pointer :: internal_to_cartesian_ADF_ptr => null()
procedure(func_ADF), pointer                  :: func_ADF_ptr => null()


abstract interface
  subroutine internal_to_cartesian_ADF(molec, internal_ADF, cartesian_ADF)
    use accuracy
    use adf
    import HM_molec_type
    type(HM_molec_type), intent(in) :: molec
    type(adf_realq), intent(in)  :: internal_ADF(molec%nmodes)
    type(adf_realq), intent(out) :: cartesian_ADF(molec%natoms,3)
  end subroutine internal_to_cartesian_ADF
end interface


abstract interface
  subroutine internal_to_cartesian(molec, internal, cartesian)
    use accuracy
    import HM_molec_type
    type(HM_molec_type), intent(in) :: molec
    real(ark), intent(in)  :: internal(molec%nmodes)
    real(ark), intent(out) :: cartesian(molec%natoms,3)
  end subroutine internal_to_cartesian
end interface


abstract interface
  subroutine func_ADF(molec, func, internal_ADF, f_ADF, cart)
    use accuracy
    use adf
    import HM_molec_type
    import HM_func_type
    type(HM_molec_type), intent(in) :: molec
    type(HM_func_type), intent(in) :: func
    type(adf_realq), intent(in)  :: internal_ADF(molec%nmodes)
    type(adf_realq), intent(out) :: f_ADF(func%rank)
    type(adf_realq), intent(in), optional :: cart(molec%natoms,3)
  end subroutine func_ADF
end interface


contains


#ifdef _MOL_XY2_
#include 'mol_xy2.f90'
#include 'pot_xy2.f90'
#endif
#ifdef _MOL_XY3_
#include 'mol_xy3.f90'
#include 'pot_xy3.f90'
#endif
#ifdef _MOL_XY4_
#include 'mol_xy4.f90'
!#include 'pot_xy4.f90'
#include 'pot_xy4_2.f90'
#endif
#ifdef _MOL_ZXY3_
#include 'mol_ch3cl.f90'
#include 'pot_ch3cl.f90'
#endif
#ifdef _MOL_C2H4_
#include 'mol_c2h4.f90'
#include 'pot_c2h4.f90'
!#include 'pot_c2h4_2.f90'
#endif
#ifdef _MOL_C2H2_
#include 'mol_c2h2.f90'
#include 'pot_c2h2.f90'
#endif
#ifdef _MOL_H2O2_
#include 'mol_h2o2.f90'
#include 'pot_h2o2.f90'
#endif


#include 'potref.f90'
#include 'expansion.f90'
#include 'coord_transform.f90'
#include 'numrec.f90'


!################################################################################


subroutine init_coord_transform(molec)

  type(HM_molec_type), intent(in) :: molec

  select case(trim(molec%coord_transform))
#ifdef _MOL_XY2_
  case('XY2_RRHO')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy2_rrho
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy2_rrho_ADF
  case('XY2_RRHO_ZMAT')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy2_rrho_zmat
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy2_rrho_zmat_ADF
  case('XY2_RRHO_TESTSING')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy2_rrho_testsing
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy2_rrho_testsing_ADF
#endif
#ifdef _MOL_XY3_
  case('XY3_SYMBETA_TAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy3_symbeta_tau
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy3_symbeta_tau_ADF
  case('XY3_SYMBETA_SINTAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy3_symbeta_sintau
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy3_symbeta_sintau_ADF
  case('XY3_SYMBETA_RHO')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy3_symbeta_rho
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy3_symbeta_rho_ADF
  case('XY3_SYMALPHA_TAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy3_rsymalpha_tau_pas
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy3_rsymalpha_tau_pas_ADF
  case('XY3_RALPHA_ZMAT')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy3_ralpha_zmat
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy3_ralpha_zmat_ADF
  case('XY3_RALPHA_PAS')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy3_ralpha_pas
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy3_ralpha_pas_ADF
#endif
#ifdef _MOL_ZXY3_
  case('ZXY3_BETA_SYMTAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_ch3cl_beta_symtau
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_ch3cl_beta_symtau_ADF
#endif
#ifdef _MOL_C2H4_
  case('C2H4_2BETA_1TAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_c2h4_2beta_1tau
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_c2h4_2beta_1tau_ADF
#endif
#ifdef _MOL_C2H2_
  case('C2H2_RQXQY')
    internal_to_cartesian_ptr      => internal_to_cartesian_c2h2_rqxqy
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_c2h2_rqxqy_ADF
  case('C2H2_RXYZ')
    internal_to_cartesian_ptr      => internal_to_cartesian_c2h2_rxyz
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_c2h2_rxyz_ADF
  case('C2H2_RSYMALPHATAU','C2H2_RALPHATAU','C2H2_RDALPHATAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_c2h2_rsymalphatau
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_c2h2_rsymalphatau_ADF
  case('C2H2_RALPHA')
    internal_to_cartesian_ptr      => internal_to_cartesian_c2h2_ralpha
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_c2h2_ralpha_ADF
  case('C2H2_R4ALPHA')
    internal_to_cartesian_ptr      => internal_to_cartesian_c2h2_r4alpha
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_c2h2_r4alpha_ADF
#endif
#ifdef _MOL_H2O2_
  case('H2O2_RALPHATAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_h2o2_ralphatau
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_h2o2_ralphatau_ADF
#endif
#ifdef _MOL_XY4_
  case('XY4_RSYMALPHA')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy4_rsymalpha
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy4_rsymalpha_ADF
  case('XY4_RSYMBETA_TAU')
    internal_to_cartesian_ptr      => internal_to_cartesian_xy4_symbeta_tau
    internal_to_cartesian_ADF_ptr  => internal_to_cartesian_xy4_symbeta_tau_ADF
#endif
  case default
    write(out, '(/a,a,a)') 'init_coord_transform error: unknown internal-to-Cartesian coordinate transformation = "', trim(molec%coord_transform), '"'
    stop
  end select

end subroutine init_coord_transform


!################################################################################


subroutine init_func(func)

  type(HM_func_type), intent(in) :: func

  select case(trim(func%func_type))
#ifdef _MOL_XY2_
  case('POTEN_XY2_MORBID')
    func_ADF_ptr  => poten_xy2_morbid_ADF
#endif
#ifdef _MOL_XY3_
  case('POTEN_XY3_MORBID_10')
    func_ADF_ptr  => poten_xy3_morbid_ADF
  case('POTEN_XY3_MORBID_10_DBOC')
    func_ADF_ptr  => poten_xy3_morbid_DBOC_ADF
  case('DIPOLE_XY3_MB')
    func_ADF_ptr  => dipole_xy3_mb2xyz_ADF
  case('DIPOLE_XY3_SYMMB','POLARIZ_XY3_SYMMB')
    func_ADF_ptr  => dipole_xy3_symmb2xyz_ADF
  case('XY3_SYMFUNC_MORBID_10')
    func_ADF_ptr  => extfunc_xy3_morbid_ADF
#endif
#ifdef _MOL_ZXY3_
  case('POTEN_CH3CL_SYM','POTEN_ZXY3_SYM')
    func_ADF_ptr  => poten_ch3cl_sym_ADF
  case('ZXY3_SYM')
    func_ADF_ptr  => dipole_ch3cl_sym_ADF
#endif
#ifdef _MOL_C2H4_
!  case('POTEN_C2H4_88666678')
!    func_ADF_ptr  => poten_c2h4_88666678_ADF
  case('POTEN_C2H4_886666')
    func_ADF_ptr  => poten_c2h4_886666_ADF
!    func_ADF_ptr  => poten_c2h4_886666_ADF_noreexp
  case('DIPOLE_C2H4_4M')
    func_ADF_ptr  => dipole_c2h4_4m_ADF
#endif
#ifdef _MOL_C2H2_
  case('POTEN_C2H2_MORSE')
    func_ADF_ptr  => poten_c2h2_morse_ADF
  case('POTEN_C2H2_STREY_MILLS')
    func_ADF_ptr  => poten_c2h2_strey_mills_ADF
  case('POTEN_C2H2_7')
    func_ADF_ptr  => poten_c2h2_7coords_ADF
#endif
#ifdef _MOL_H2O2_
  case('POTEN_H2O2_KOPUT')
    func_ADF_ptr  => poten_h2o2_koput_ADF
#endif
#ifdef _MOL_XY4_
  case('POTEN_XY4_ALPHA')
  !  func_ADF_ptr  => poten_xy4_alpha_ADF
    func_ADF_ptr  => poten_xy4_alpha_ADF_noreexp
  case('DIPOLE_XY4_ALPHA')
  !  func_ADF_ptr  => dipole_xy4_alpha_ADF
    func_ADF_ptr  => dipole_xy4_alpha_ADF_noreexp
#endif
  case('POTENTIAL')
    func_ADF_ptr => potref_general
  case default
    write(out, '(/a,a,a)') 'init_func error: unknown function = "', trim(func%func_type), '"'
    stop
  end select

end subroutine init_func


!################################################################################


subroutine read_kinetic(molec, nmax_g, maxdeg_g, nmax_u, maxdeg_u, x0_npoints, x0, nterms_gmat, terms_gmat, coefs_gmat, nterms_u, terms_u, coefs_u)

  type(HM_molec_type), intent(in)         :: molec
  integer(ik), intent(in)                 :: nmax_g, nmax_u, maxdeg_g, maxdeg_u
  integer(ik), intent(in)                 :: x0_npoints
  real(ark), intent(in)                   :: x0(:,:)
  integer(ik), intent(out)                :: nterms_gmat
  integer(ik), allocatable, intent(out)   :: terms_gmat(:,:)
  real(ark), allocatable, intent(out)     :: coefs_gmat(:,:,:,:)
  integer(ik), intent(out)                :: nterms_u
  integer(ik), allocatable, intent(out)   :: terms_u(:,:)
  real(ark), allocatable, intent(out)     :: coefs_u(:,:)

  integer(ik), parameter :: nterms_max = 1000000
  integer(ik) :: nmodes, natoms3, iounit, info, iline, term(molec%nmodes), x0_npoints_, nmodes_, maxorder, ipoint, jpoint, lambda, mu, ipos, nterms
  integer(ik), allocatable :: terms(:,:), ind_ipoint(:)
  real(ark) :: gval, x0_(molec%nmodes), planck_factor
  character(300) :: sbuf, sbuf2, sbuf3

  write(out, '(/a)') 'read_kinetic/start: read KEO from file'

  planck_factor = real(planck,ark)*real(avogno,ark)*(1.0q+16)/(4.0_ark*real(pi,ark)*real(pi,ark)*real(vellgt,ark))

  nmodes = molec%nmodes
  natoms3 = molec%natoms*3
  iounit = 1000

  write(out, '(/1x,a,1x,i3)') 'max requested expansion order for G matrix:', maxdeg_g
  write(out, '(1x,a,1x,i3)') 'max requested expansion order for pseudo-potential:', maxdeg_u
  write(out, '(1x,a,1x,i3)') 'requested number of reference points:', x0_npoints

  if (molec%expkeo_verbose>=3) then
    write(out, '(1x,a,1x,i6,1x,a)') 'reference geometries (', x0_npoints, '):'
    do ipoint=1, x0_npoints
      write(out, '(1x,i6,3x,<nmodes>(1x,es16.8))') ipoint, x0(1:nmodes,ipoint)
    enddo
  endif


  !--------------------------------------------------------------------------------------
  ! read expansions for kinetic energy matrix G_{lambda,mu}, where lambda and mu = 1..3N
  !--------------------------------------------------------------------------------------

  write(out, '(/1x,a,1x,a,1x,a,1x,i4,1x,a)') 'read elements of G-matrix from file', trim(molec%fname_gmat), '(I/O unit =', iounit, ')'

  open(unit=iounit,form='formatted',action='read',position='rewind',status='old',file=molec%fname_gmat,iostat=info)
  if (info/=0) then
    write(out, '(/a,1x,a)') 'read_kinetic error: could not open file', trim(molec%fname_gmat)
    stop
  endif

  ! read number of expansion points, max expansion orders, and reference geometries

  rewind(iounit)

  read(iounit,'(a)') sbuf
  read(iounit,'(a)') sbuf2
  read(iounit,'(a)') sbuf3

  write(out, '(1x,a,1x,1(1h|),1x,a,2(/14x,1(1h|),1x,a))') 'file header:', trim(sbuf), trim(sbuf2), trim(sbuf3)

  read(iounit,*) sbuf, nmodes_
  read(iounit,*) sbuf, maxorder
  read(iounit,*) sbuf, x0_npoints_

  write(out, '(1x,a,1x,i4)') 'number of vibrational modes:', nmodes_
  write(out, '(1x,a,1x,i4)') 'max expansion order:', maxorder
  write(out, '(1x,a,1x,i4)') 'number of reference points:', x0_npoints_

  if (nmodes_/=nmodes) then
    write(out, '(/a,1x,i3,1x,a,1x,i3,1x,a)') 'read_kinetic error: number of vibrational modes =', nmodes_, 'differs from Nmodes =', nmodes, 'defined in the input'
    stop
  endif

  allocate(ind_ipoint(x0_npoints_), terms(nmodes,nterms_max), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'read_kinetic error: failed to allocate ind_ipoint(x0_npoints_), terms(nmodes,nterms_max)', &!
    'x0_npoints_, nmodes, nterms_max =', x0_npoints_, nmodes, nterms_max
    stop
  endif

  ! read reference geometries from the G-matrix file and match them with those in the input to this subroutine
  ind_ipoint = 0
  do ipoint=1, x0_npoints_
    read(iounit,*) x0_(1:nmodes)
    do jpoint=1, x0_npoints
      if (all(abs(x0(1:nmodes,jpoint)-x0_(1:nmodes))<=1.0d-14)) then
        ind_ipoint(ipoint) = jpoint
        exit
      endif
    enddo
  enddo

  ! check if all input reference geometries have been found in the G-matrix file
  do ipoint=1, x0_npoints
    if (.not.any(ind_ipoint(1:x0_npoints_)==ipoint)) then
      write(out, '(/a,<nmodes>(/1x,es16.8)/a,1x,a)') &!
      'read_kinetic error: could not match input geometry = (', x0(1:nmodes,ipoint), ') with any from the G-matrix file', trim(molec%fname_gmat)
      stop
    endif
  enddo

  ! estimate number of expansion terms

  iline = x0_npoints_+6
  nterms = 0
  terms = 0
  do
    read(iounit,*,iostat=info) ipoint, lambda, mu, term(1:nmodes)
    iline = iline + 1
    if (info>0) then
      write(out, '(/a,1x,a,1x,a,1x,i6)') 'read_kinetic error while reading file', trim(molec%fname_gmat), ', line =', iline
      stop
    elseif (info<0) then
      exit
    endif
    if (ind_ipoint(ipoint)==0) cycle
    if (count(term(1:nmodes)>0)>nmax_g) cycle
    if (sum(term(1:nmodes))>maxdeg_g) cycle
    if (nterms==0) then
      ipos = 0
    else
      ipos = index_iarr1(term(1:nmodes), terms(1:nmodes,1:nterms))
    endif
    if (ipos==0) then
      nterms = nterms + 1
      if (nterms>nterms_max) then
        write(out, '(/a,1x,i6,1x,a/a)') 'read_kinetic error: number of expansion terms =', nterms, 'exceeds size of array "terms"', &!
        'increase value of "nterms_max" in subroutine "read_kinetic"'
        stop
      endif
      terms(1:nmodes,nterms) = term(1:nmodes)
    endif
  enddo

  nterms_gmat = nterms

  write(out, '(1x,a,1x,i6)') 'number of expansion terms:', nterms_gmat

  ! allocate arrays to keep expansion terms and coefficients

  if (allocated(terms_gmat)) deallocate(terms_gmat)
  if (allocated(coefs_gmat)) deallocate(coefs_gmat)

  allocate(coefs_gmat(natoms3,natoms3,nterms_gmat,x0_npoints), terms_gmat(nmodes,nterms_gmat), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_gmat(natoms3,natoms3,nterms_gmat,x0_npoints), terms_gmat(nmodes,nterms_gmat)', &!
    'natoms3, nterms_gmat, x0_npoints =', natoms3, nterms_gmat, x0_npoints
    stop
  endif
  coefs_gmat = 0.0
  terms_gmat = 0

  ! read expansion terms and coefficients form file

  rewind(iounit)
  read(iounit,*) sbuf
  read(iounit,*) sbuf
  read(iounit,*) sbuf
  read(iounit,*) sbuf, nmodes_
  read(iounit,*) sbuf, maxorder
  read(iounit,*) sbuf, x0_npoints_
  do ipoint=1, x0_npoints_
    read(iounit,*) x0_(1)
  enddo

  iline = x0_npoints_+6
  nterms_gmat = 0
  do
    read(iounit,*,iostat=info) ipoint, lambda, mu, term(1:nmodes), gval
    iline = iline + 1
    if (info>0) then
      write(out, '(/a,1x,a,1x,a,1x,i4)') 'read_kinetic error while reading file', trim(molec%fname_gmat), ', line =', iline
      stop
    elseif (info<0) then
      exit
    endif
    if (ind_ipoint(ipoint)==0) cycle
    if (count(term(1:nmodes)>0)>nmax_g) cycle
    if (sum(term(1:nmodes))>maxdeg_g) cycle
    if (nterms_gmat==0) then
      ipos = 0
    else
      ipos = index_iarr1(term(1:nmodes), terms_gmat(1:nmodes,1:nterms_gmat))
    endif
    if (ipos==0) then
      nterms_gmat = nterms_gmat + 1
      terms_gmat(1:nmodes,nterms_gmat) = term(1:nmodes)
      ipos = nterms_gmat
    endif
    coefs_gmat(lambda,mu,ipos,ind_ipoint(ipoint)) = gval * planck_factor
  enddo
  close(iounit)

  deallocate(ind_ipoint, terms)


  !--------------------------------------
  ! read expansions for pseudo-potential
  !--------------------------------------

  write(out, '(/1x,a,1x,a,1x,a,1x,i4,1x,a)') 'read pseudo-potential from file', trim(molec%fname_pseudo), '(I/O unit =', iounit, ')'

  open(unit=iounit,form='formatted',action='read',position='rewind',status='old',file=molec%fname_pseudo,iostat=info)
  if (info/=0) then
    write(out, '(/a,1x,a)') 'read_kinetic error: could not open file', trim(molec%fname_pseudo)
    stop
  endif

  ! read number of expansion points, max expansion orders, and reference geometries

  rewind(iounit)

  read(iounit,'(a)') sbuf
  read(iounit,'(a)') sbuf2
  read(iounit,'(a)') sbuf3

  write(out, '(1x,a,1x,1(1h|),1x,a,2(/14x,1(1h|),1x,a))') 'file header:', trim(sbuf), trim(sbuf2), trim(sbuf3)

  read(iounit,*) sbuf, nmodes_
  read(iounit,*) sbuf, maxorder
  read(iounit,*) sbuf, x0_npoints_

  write(out, '(1x,a,1x,i4)') 'number of vibrational modes:', nmodes_
  write(out, '(1x,a,1x,i4)') 'max expansion order:', maxorder
  write(out, '(1x,a,1x,i4)') 'number of reference points:', x0_npoints_

  if (nmodes_/=nmodes) then
    write(out, '(/a,1x,i3,1x,a,1x,i3,1x,a)') 'read_kinetic error: number of vibrational modes =', nmodes_, 'differs from Nmodes =', nmodes, 'defined in the input'
    stop
  endif

  allocate(ind_ipoint(x0_npoints_), terms(nmodes,nterms_max), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'read_kinetic error: failed to allocate ind_ipoint(x0_npoints_), terms(nmodes,nterms_max)', 'x0_npoints_, nmodes, nterms_max =', x0_npoints_, nmodes, nterms_max
    stop
  endif

  ! read reference geometries from the pseudo-potential file and match them with those in the input to this subroutine
  ind_ipoint = 0
  do ipoint=1, x0_npoints_
    read(iounit,*) x0_(1:nmodes)
    do jpoint=1, x0_npoints
      if (all(abs(x0(1:nmodes,jpoint)-x0_(1:nmodes))<=1.0d-14)) then
        ind_ipoint(ipoint) = jpoint
        exit
      endif
    enddo
  enddo

  ! check if all input reference geometries have been found in the pseudo-potential file
  do ipoint=1, x0_npoints
    if (.not.any(ind_ipoint(1:x0_npoints_)==ipoint)) then
      write(out, '(/a,<nmodes>(/1x,es16.8)/a,1x,a)') &!
      'read_kinetic error: could not match input geometry = (', x0(1:nmodes,ipoint), ') with any from the pseudo-potential file', trim(molec%fname_pseudo)
      stop
    endif
  enddo

  ! estimate number of expansion terms

  iline = x0_npoints_+6
  nterms = 0
  terms = 0
  do
    read(iounit,*,iostat=info) ipoint, term(1:nmodes)
    iline = iline + 1
    if (info>0) then
      write(out, '(/a,1x,a,1x,a,1x,i6)') 'read_kinetic error while reading file', trim(molec%fname_pseudo), ', line =', iline
      stop
    elseif (info<0) then
      exit
    endif
    if (ind_ipoint(ipoint)==0) cycle
    if (count(term(1:nmodes)>0)>nmax_u) cycle
    if (sum(term(1:nmodes))>maxdeg_u) cycle
    if (nterms==0) then
      ipos = 0
    else
      ipos = index_iarr1(term(1:nmodes), terms(1:nmodes,1:nterms))
    endif
    if (ipos==0) then
      nterms = nterms + 1
      if (nterms>nterms_max) then
        write(out, '(/a,1x,i6,1x,a/a)') 'read_kinetic error: number of expansion terms =', nterms, 'exceeds size of array "terms"', 'increase value of "nterms_max" in subroutine "read_kinetic"'
        stop
      endif
      terms(1:nmodes,nterms) = term(1:nmodes)
    endif
  enddo

  nterms_u = nterms

  write(out, '(1x,a,1x,i6)') 'number of expansion terms:', nterms_gmat

  ! allocate arrays to keep expansion terms and coefficients

  if (allocated(terms_u)) deallocate(terms_u)
  if (allocated(coefs_u)) deallocate(coefs_U)

  allocate(coefs_u(nterms_u,x0_npoints), terms_u(nmodes,nterms_u), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_u(nterms_u,x0_npoints), terms_u(nmodes,nterms_u)', &!
    'natoms3, nterms_u, x0_npoints =', natoms3, nterms_u, x0_npoints
    stop
  endif
  coefs_u = 0.0
  terms_u = 0

  ! read expansion terms and coefficients form file

  rewind(iounit)
  read(iounit,*) sbuf
  read(iounit,*) sbuf
  read(iounit,*) sbuf
  read(iounit,*) sbuf, nmodes_
  read(iounit,*) sbuf, maxorder
  read(iounit,*) sbuf, x0_npoints_
  do ipoint=1, x0_npoints_
    read(iounit,*) x0_(1)
  enddo

  iline = x0_npoints_+6
  nterms_u = 0
  do
    read(iounit,*,iostat=info) ipoint, term(1:nmodes), gval
    iline = iline + 1
    if (info>0) then
      write(out, '(/a,1x,a,1x,a,1x,i4)') 'read_kinetic error while reading file', trim(molec%fname_pseudo), ', line =', iline
      stop
    elseif (info<0) then
      exit
    endif
    if (ind_ipoint(ipoint)==0) cycle
    if (count(term(1:nmodes)>0)>nmax_u) cycle
    if (sum(term(1:nmodes))>maxdeg_u) cycle
    if (nterms_u==0) then
      ipos = 0
    else
      ipos = index_iarr1(term(1:nmodes), terms_u(1:nmodes,1:nterms_u))
    endif
    if (ipos==0) then
      nterms_u = nterms_u + 1
      terms_u(1:nmodes,nterms_u) = term(1:nmodes)
      ipos = nterms_u
    endif
    coefs_u(ipos,ind_ipoint(ipoint)) = gval * planck_factor
  enddo
  close(iounit)

  deallocate(ind_ipoint, terms)


  write(out, '(/a)') 'read_kinetic/done'

end subroutine read_kinetic


!################################################################################


subroutine expand_kinetic(molec, x0_npoints,x0, nmax_g,nmax_u,ndeg,deg,mindeg_g,maxdeg_g,mindeg_u,maxdeg_u, nterms_gmat,terms_gmat,coefs_gmat, nterms_u,terms_u,coefs_u, coefs_dgmat)

  type(HM_molec_type), intent(in)               :: molec
  integer(ik), intent(in)                       :: x0_npoints
  real(ark), intent(in)                         :: x0(molec%nmodes,x0_npoints)
  integer(ik), intent(in)                       :: nmax_g
  integer(ik), intent(in)                       :: nmax_u
  integer(ik), intent(in)                       :: ndeg(molec%nmodes)
  integer(ik), intent(in)                       :: deg(maxval(ndeg),molec%nmodes)
  integer(ik), intent(in)                       :: mindeg_g(nmax_g)
  integer(ik), intent(in)                       :: maxdeg_g(nmax_g)
  integer(ik), intent(in)                       :: mindeg_u(nmax_u)
  integer(ik), intent(in)                       :: maxdeg_u(nmax_u)
  integer(ik), intent(out)                      :: nterms_gmat
  integer(ik), allocatable, intent(out)         :: terms_gmat(:,:)
  real(ark), allocatable, intent(out)           :: coefs_gmat(:,:,:,:)
  integer(ik), intent(out)                      :: nterms_u
  integer(ik), allocatable, intent(out)         :: terms_u(:,:)
  real(ark), allocatable, intent(out)           :: coefs_u(:,:)
  real(ark), allocatable, intent(out), optional :: coefs_dgmat(:,:,:,:,:)

  integer(ik) :: natoms, nmodes, info, nterms_svec, nterms_tmp, nterms_cart, icomb, ncomb_cart, n, max_order, nterms, iatom, iterm, imode, ipos, i, ipoint, &!
                 natoms3, max_order_modes(molec%nmodes), verbose, nterms_gmat_d, term0(molec%nmodes)
  integer(ik), allocatable :: terms_svec(:,:), terms_tmp(:,:), terms_cart(:,:), comb_cart(:,:), nterms_cart_split(:), &!
                              terms_cart_split(:,:,:), terms(:,:), nterms_chr(:), ind_chr(:,:,:), terms_gmat_d(:,:)
  real(ark) :: fac, planck_factor
  real(ark), allocatable :: cart0(:,:,:), cart0_ref(:,:,:), rotmat0_ref(:,:,:), coefs(:,:,:,:), coefs_rotmat(:,:,:,:), coefs_cart_ref(:,:,:,:), &!
                            coefs_cart(:,:,:,:), coefs_svec(:,:,:,:), coefs_gmat_d(:,:,:,:)
  logical, allocatable :: coefs_cart_init(:)


  planck_factor = real(planck,ark)*real(avogno,ark)*(1.0q+16)/(4.0_ark*real(pi,ark)*real(pi,ark)*real(vellgt,ark))

  natoms = molec%natoms
  nmodes = molec%nmodes
  natoms3 = natoms*3
  verbose = molec%expkeo_verbose


  if (verbose>=3) then
    write(out, '(/a)') 'expand_kinetic/start: power series expansion of kinetic energy matrix and pseudo-potential'
    write(out, '(/1x,a,<nmodes>(1x,i3))') 'min N-mode expansion orders for kinetic energy matrix:', mindeg_g(1:nmodes)
    write(out, '(1x,a,<nmodes>(1x,i3))') 'max N-mode expansion orders for kinetic energy matrix:', maxdeg_g(1:nmodes)
    write(out, '(1x,a,<nmodes>(1x,i3))') 'min N-mode expansion orders for pseudo-potential:', mindeg_u(1:nmodes)
    write(out, '(1x,a,<nmodes>(1x,i3))') 'max N-mode expansion orders for pseudo-potential:', maxdeg_u(1:nmodes)
    write(out, '(1x,a)') 'expansion orders for each mode:'
    do imode=1, nmodes
     write(out, '(1x,i3,5x,100(1x,i3))') imode, deg(1:ndeg(imode),imode)
    enddo
  endif

  if (verbose>=3) then
    write(out, '(/1x,a,1x,i6,1x,a)') 'reference internal coordinates (', x0_npoints, '):'
    do ipoint=1, x0_npoints
      write(out, '(1x,i6,3x,<nmodes>(1x,es16.8))') ipoint, x0(1:nmodes,ipoint)
    enddo
  endif


  !-----------------------------------------!
  ! compute reference Cartesian coordinates !
  !-----------------------------------------!

  ! allocate arrays to store reference geometries

  allocate(cart0(molec%natoms,3,x0_npoints), cart0_ref(molec%natoms,3,x0_npoints), rotmat0_ref(3,3,x0_npoints), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate cart0(molec%natoms,3,x0_npoints), cart0_ref(molec%natoms,3,x0_npoints), rotmat0_ref(3,3,x0_npoints)', &!
    'molec%nmodes, x0_npoints =', molec%nmodes, x0_npoints
  endif
  cart0 = 0.0
  cart0_ref = 0.0
  rotmat0_ref = 0.0

  ! initialize internal-to-Cartesian transformation function

  call init_coord_transform(molec)

  if (associated(internal_to_cartesian_ptr).eqv..false.) then
    write(out, '(/a)') 'expand_kinetic error: function pointer "internal_to_cartesian_ptr" is not associated'
    stop
  endif

  ! compute reference Cartesian coordinates in molec%coord_transform axes system

  do ipoint=1, x0_npoints
    call internal_to_cartesian_ptr(molec, x0(1:nmodes,ipoint), cart0(1:natoms,1:3,ipoint))
  enddo

  if (verbose>=3) then
    write(out, '(/1x,a,1x,a,1x,a)') 'Cartesian coordinates in', trim(molec%coord_transform), 'frame'
    do ipoint=1, x0_npoints
      do iatom=1, natoms
        write(out, '(1x,i6,3x,i3,3x,3(f20.16))') ipoint, iatom, cart0(iatom,1:3,ipoint)
      enddo
    enddo
  endif

  ! rotate reference Cartesian coordinates to molec%reference_frame axes system

  if(trim(molec%reference_frame)==trim(molec%coord_transform)) then
    cart0_ref = cart0
    rotmat0_ref = 0.0
    forall(i=1:3) rotmat0_ref(i,i,:) = 1.0_ark
  elseif(trim(molec%reference_frame)=='PAS') then
    if (trim(molec%pas_function)=='EXPKAPPA') then
      do ipoint=1, x0_npoints
        call rotate_pas(molec, cart0(1:natoms,1:3,ipoint), cart0_ref(1:natoms,1:3,ipoint), rotmat0_ref(1:3,1:3,ipoint))
      enddo
    elseif (trim(molec%pas_function)=='LAPACK') then
      do ipoint=1, x0_npoints
        call rotate_pas_diag(molec, cart0(1:natoms,1:3,ipoint), cart0_ref(1:natoms,1:3,ipoint), rotmat0_ref(1:3,1:3,ipoint))
      enddo
    else
      write(out, '(/a,1x,a)') 'expand_kinetic error: unknown principal-axes-system rotation function =', trim(molec%pas_function)
      stop
    endif
#ifdef _MOL_XY4_
  elseif(trim(molec%reference_frame)=='XY4_TD') then
    do ipoint=1, x0_npoints
      call internal_to_cartesian_xy4_equilibrium_td(molec, x0(1:nmodes,ipoint), cart0_ref(1:natoms,1:3,ipoint))
      call rotate_to_xyz(molec, cart0(1:natoms,1:3,ipoint), cart0_ref(1:natoms,1:3,ipoint), rotmat0_ref(1:3,1:3,ipoint))
      cart0_ref(1:natoms,1:3,ipoint) = matmul(cart0(1:natoms,1:3,ipoint),transpose(rotmat0_ref(1:3,1:3,ipoint)))
    enddo
#endif
  else
    write(out, '(/a,a,a)') 'expand_kinetic error: unknown reference frame = "', trim(molec%reference_frame), '"'
    stop
  endif

  if (verbose>=3) then
    if (trim(molec%reference_frame)/=trim(molec%coord_transform)) then
      write(out, '(/1x,a,1x,a,1x,a)') 'Cartesian coordinates in', trim(molec%reference_frame), 'frame'
      do ipoint=1, x0_npoints
        do iatom=1, natoms
          write(out, '(1x,i6,3x,i3,3x,3(f20.16))') ipoint, iatom, cart0_ref(iatom,1:3,ipoint)
        enddo
      enddo
    endif
  endif


  !---------------------------------!
  ! generate N-mode expansion terms !
  !---------------------------------!

  ! generate expansion terms for G-matrix

  if (verbose>=3) then
    write(out, '(/1x,a)') 'generate expansion terms for kinetic energy matrix'
  endif

  call nmode_expansion(nmodes, nmax_g, ndeg, deg, mindeg_g, maxdeg_g, nterms_gmat, terms_gmat)
  call delete_same_terms(nmodes, nterms_gmat, terms_gmat)

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6,1x,a)') 'expansion terms (', nterms_gmat, ')'
    if (verbose>=6) then
      do iterm=1, nterms_gmat
        write(out, '(1x,i6,5x,<nmodes>(1x,i3))') iterm, terms_gmat(1:nmodes,iterm)
      enddo
    endif
  endif

  nterms_gmat_d = nterms_gmat
  allocate(terms_gmat_d(nmodes,nterms_gmat_d), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate terms_gmat_d(nmodes,nterms_gmat_d)', 'nmodes, nterms_gmat_d =', nmodes, nterms_gmat_d
    stop
  endif
  terms_gmat_d = terms_gmat

  if (present(coefs_dgmat)) then
    call expansion_terms_dgmat2gmat(nmodes, nterms_gmat_d, terms_gmat_d)
    call delete_same_terms(nmodes, nterms_gmat_d, terms_gmat_d)
    !call sort_terms_order(nmodes, nterms_gmat_d, terms_gmat_d)
    if (verbose>=3) then
      write(out, '(1x,a,1x,i6,1x,a)') 'expansion terms required for first-order derivatives of kinetic energy matrix (', nterms_gmat_d, ')'
      if (verbose>=6) then
        do iterm=1, nterms_gmat_d
          write(out, '(1x,i6,5x,<nmodes>(1x,i3))') iterm, terms_gmat_d(1:nmodes,iterm)
        enddo
      endif
    endif
  endif

  ! generate expansion terms for pseudo-potential

  if (verbose>=3) then
    write(out, '(/1x,a)') 'generate expansion terms for pseudo-potential'
  endif

  call nmode_expansion(nmodes, nmax_u, ndeg, deg, mindeg_u, maxdeg_u, nterms_u, terms_u)
  call delete_same_terms(nmodes, nterms_u, terms_u)

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6,1x,a)') 'expansion terms (', nterms_u, ')'
    if (verbose>=6) then
      do iterm=1, nterms_u
        write(out, '(1x,i6,5x,<nmodes>(1x,i3))') iterm, terms_u(1:nmodes,iterm)
      enddo
    endif
  endif

  ! generate complete set of expansion terms required for recursive determination of all s-vector derivatives

  if (verbose>=3) then
    write(out, '(/1x,a)') 'generate expansion terms for s-vectors'
  endif

  nterms_tmp = nterms_u
  allocate(terms_tmp(nmodes,nterms_tmp), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate terms_tmp(nmodes,nterms_tmp)', 'nmodes, nterms_tmp =', nmodes, nterms_tmp
    stop
  endif
  terms_tmp = terms_u
  call expansion_terms_pseudo2svec(nmodes, nterms_tmp, terms_tmp)

  nterms_svec = nterms_tmp + nterms_gmat_d
  allocate(terms_svec(nmodes,nterms_svec), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate terms_svec(nmodes,nterms_svec)', 'nmodes, nterms_svec =', nmodes, nterms_svec
    stop
  endif
  terms_svec(:,1:nterms_tmp) = terms_tmp
  terms_svec(:,nterms_tmp+1:) = terms_gmat_d
  deallocate(terms_tmp)
  call delete_same_terms(nmodes, nterms_svec, terms_svec)
  !call expansion_terms_complete(nmodes, nterms_svec, terms_svec)
  call sort_terms_order(nmodes, nterms_svec, terms_svec)

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6,1x,a)') 'expansion terms (', nterms_svec, ')'
    if (verbose>=6) then
      do iterm=1, nterms_svec
        write(out, '(1x,i6,5x,<nmodes>(1x,i3))') iterm, terms_svec(1:nmodes,iterm)
      enddo
    endif
  endif

  ! generate expansion terms for Cartesian coordinates

  if (verbose>=3) then
    write(out, '(/1x,a)') 'generate expansion terms for Cartesian coordinates'
  endif

  nterms_cart = nterms_svec
  if (allocated(terms_cart)) deallocate(terms_cart)
  allocate(terms_cart(nmodes,nterms_cart), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate terms_cart(nmodes,nterms_cart)', 'nmodes, nterms_cart =', nmodes, nterms_cart
    stop
  endif
  terms_cart = terms_svec
  call expansion_terms_svec2cart(nmodes, nterms_cart, terms_cart)
  call delete_same_terms(nmodes, nterms_cart, terms_cart)
  !call expansion_terms_complete(nmodes, nterms_cart, terms_cart)
  call sort_terms_order(nmodes, nterms_cart, terms_cart)

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6,1x,a)') 'expansion terms (', nterms_cart, ')'
    if (verbose>=6) then
      do iterm=1, nterms_cart
        write(out, '(1x,i6,5x,<nmodes>(1x,i3))') iterm, terms_cart(1:nmodes,iterm)
      enddo
    endif
  endif

  ! split expansion terms for Cartesian coordinates into groups for different combinations of modes

  if (verbose>=3) then
    write(out, '(/1x,a)') 'cast expansion terms for Cartesian coordinates into the N-mode expansion form'
  endif

  call nmode_distribute(nmodes, nterms_cart, terms_cart, ncomb_cart, comb_cart, nterms_cart_split, terms_cart_split)

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6)') 'number of mode-combinations:', ncomb_cart
    do icomb=1, ncomb_cart
      n = count(comb_cart(:,icomb)>0)
      write(out, '(1x,a,1x,i6,3x,a,1x,i6,3x,a,100(1x,i3))') 'icomb:', icomb, 'nterms:', nterms_cart_split(icomb), 'imodes:', comb_cart(1:n,icomb)
    enddo
  endif


  !----------------------------------------------!
  ! compute derivatives of Cartesian coordinates !
  !----------------------------------------------!

  ! 1. Compute derivatives of Cartesian coordinates in a frame defined by molec%coord_transform function
  ! 2. Compute derivatives of frame rotation matrix (due to Eckart, PAS, etc. conditions)
  ! 3. Apply chain rules to compute derivatives of "rotated" Cartesian coordinates

  if (verbose>=3) then
    write(out, '(/1x,a)') 'compute derivatives of Cartesian coordinates'
  endif

  allocate(coefs_cart(natoms,3,nterms_cart,x0_npoints), coefs_cart_ref(natoms,3,nterms_cart,x0_npoints), coefs_rotmat(3,3,nterms_cart,x0_npoints), coefs_cart_init(nterms_cart), stat=info)
  if (info/=0) then
    write(out, '(/a,1x,a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_cart(natoms,3,nterms_cart,x0_npoints), coefs_cart_ref(..)', &!
    'coefs_rotmat(3,3,nterms_cart,x0_npoints), coefs_cart_init(nterms_cart)', 'natoms, nterms_cart, x0_npoints =', natoms, nterms_cart, x0_npoints
    stop
  endif

  ! (1) compute derivatives of Cartesian coordinates in COORD_TRANSFORM axes system

  if (verbose>=3) then
    write(out, '(1x,a,1x,a)') trim(molec%coord_transform), 'frame ...'
  endif

  coefs_cart_ref = 0.0
  coefs_cart_init = .false.

  !$omp parallel do private(icomb,max_order,n,max_order_modes,imode,nterms,terms,coefs,iterm,ipoint,ipos) schedule(dynamic)
  do icomb=1, ncomb_cart

    ! estimate total max derivative order, max derivative orders for each mode, and N for current combination of modes

    max_order = maxval(sum(terms_cart_split(1:nmodes,1:nterms_cart_split(icomb),icomb),dim=1))
    n = count(comb_cart(:,icomb)>0)
    max_order_modes = maxval((/( terms_cart_split(comb_cart(imode,icomb),1:nterms_cart_split(icomb),icomb), imode=1, n )/))

    if (verbose>=3) then
      write(out, '(1x,a,1x,i6,1x,a,100(1x,i3))') 'icomb:', icomb, 'max_order:', max_order_modes(1:n)
    endif

    ! compute derivatives using ADF

    call deriv_cart_ADF(molec, x0_npoints, x0, max_order, n, comb_cart(1:n,icomb), max_order_modes(1:n), nterms, terms, coefs)

    ! rotate to equilibrium-reference axes system

    do ipoint=1, x0_npoints
      do iterm=1, nterms
        coefs(1:natoms,1:3,iterm,ipoint) = matmul(coefs(1:natoms,1:3,iterm,ipoint), transpose(rotmat0_ref(1:3,1:3,ipoint)))
      enddo
    enddo

    ! keep derivatives for current icomb-combination of modes

    do iterm=1, nterms
      ipos = index_iarr1(terms(1:nmodes,iterm), terms_cart(1:nmodes,1:nterms_cart))
      if (ipos/=0) then
        coefs_cart_ref(1:natoms,1:3,ipos,1:x0_npoints) = coefs(1:natoms,1:3,iterm,1:x0_npoints)
        coefs_cart_init(ipos) = .true.
      endif
    enddo

    deallocate(terms)
    deallocate(coefs)

  enddo ! icomb
  !$omp end parallel do

  ! check if all Cartesian-coordinate derivatives have been initialized

  if (.not.all(coefs_cart_init)) then
    write(out, '(/a)') 'expand_kinetic error: following Cartesian-coordinate derivatives have not been initialized'
    do iterm=1, nterms_cart
      if (.not.coefs_cart_init(iterm)) then
        write(out, '(1x,i3,5x,<nmodes>(1x,i3))') iterm, terms_cart(1:nmodes,iterm)
      endif
    enddo
    stop
  endif

  ! (2) compute derivatives of frame rotation matrix

  if (verbose>=3) then
    write(out, '(1x,a,1x,a)') trim(molec%frame_rotation), 'rotation ...'
  endif

  ! initialize chain-rule expressions for product cartesian * rotmat

  allocate(nterms_chr(nterms_cart), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate nterms_chr(nterms_cart)', 'nterms_cart =', nterms_cart
    stop
  endif
  call chain_rule_product_ind(nmodes, nterms_cart, terms_cart, 2, nterms_chr, ind_chr)

  coefs_rotmat = 0.0

  if (trim(molec%frame_rotation)=='ECKART') then

    !$omp parallel do private(ipoint) schedule(dynamic)
    do ipoint=1, x0_npoints
      call deriv_rotmat_eckart(molec, cart0_ref(1:natoms,1:3,ipoint), nterms_cart, terms_cart, nterms_chr, ind_chr, &!
                               coefs_cart_ref(1:natoms,1:3,1:nterms_cart,ipoint), coefs_rotmat(1:3,1:3,1:nterms_cart,ipoint))
    enddo
    !$omp end parallel do

  elseif (trim(molec%frame_rotation)=='I') then

    coefs_rotmat = 0.0
    forall(iterm=1:nterms_cart,i=1:3,all(terms_cart(1:nmodes,iterm)==0)) coefs_rotmat(i,i,iterm,1:x0_npoints) = 1.0_ark

  else

    write(out, '(/a,a,a)') 'expand_kinetic error: unknown molecule-fixed frame embedding = "', trim(molec%frame_rotation), '"'
    stop

  endif

  ! (3) compute derivatives of Cartesian coordinates in FRAME_ROTATION axes system

  if (verbose>=3) then
    write(out, '(1x,a,1x,a)') trim(molec%frame_rotation), 'frame ...'
  endif

  if (trim(molec%frame_rotation)=='I') then

    coefs_cart = coefs_cart_ref

  else

    !$omp parallel do private(ipoint) schedule(dynamic)
    do ipoint=1, x0_npoints
      call deriv_cart_rotated(molec, nterms_cart, terms_cart, nterms_chr, ind_chr, coefs_cart_ref(1:natoms,1:3,1:nterms_cart,ipoint), &!
                              coefs_rotmat(1:3,1:3,1:nterms_cart,ipoint), coefs_cart(1:natoms,1:3,1:nterms_cart,ipoint))
    enddo
    !$omp end parallel do

  endif

  deallocate(nterms_chr, ind_chr)
  deallocate(coefs_rotmat, coefs_cart_init)


  !----------------------------------!
  ! compute derivatives of s-vectors !
  !----------------------------------!

  if (verbose>=3) then
    write(out, '(/1x,a)') 'compute derivatives of s-vectors'
  endif

  ! allocate arrays to store derivatives of s-vectors

  allocate(coefs_svec(natoms3,natoms3,nterms_svec,x0_npoints), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_svec(natoms3,natoms3,nterms_svec,x0_npoints)', &!
    'natoms3, nterms_svec, x0_npoints =', natoms3, nterms_svec, x0_npoints
    stop
  endif
  coefs_svec = 0.0

  ! initialize chain-rule expressions for product t-vector * s-vector

  allocate(nterms_chr(nterms_svec), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate nterms_chr(nterms_svec)', 'nterms_svec =', nterms_svec
    stop
  endif
  call chain_rule_product_ind(nmodes, nterms_svec, terms_svec, 2, nterms_chr, ind_chr)

  ! compute derivatives of s-vectors at each point

  if (molec%vib_3n5) then

    !$omp parallel do private(ipoint) schedule(dynamic)
    do ipoint=1, x0_npoints
      call deriv_svec_3n5(molec, nterms_cart, terms_cart, coefs_cart(1:natoms,1:3,1:nterms_cart,ipoint), nterms_svec, terms_svec, nterms_chr, ind_chr, &!
                          coefs_svec(1:natoms3,1:natoms3,1:nterms_svec,ipoint))
    enddo
    !$omp end parallel do

  else

    !$omp parallel do private(ipoint) schedule(dynamic)
    do ipoint=1, x0_npoints
      call deriv_svec(molec, nterms_cart, terms_cart, coefs_cart(1:natoms,1:3,1:nterms_cart,ipoint), nterms_svec, terms_svec, nterms_chr, ind_chr, &!
                      coefs_svec(1:natoms3,1:natoms3,1:nterms_svec,ipoint))
    enddo
    !$omp end parallel do

  endif

  ! deallocate chain rule expressions

  deallocate(nterms_chr, ind_chr)


  !---------------------------------!
  ! compute derivatives of G-matrix !
  !---------------------------------!

  if (verbose>=3) then
    write(out, '(/1x,a)') 'compute derivatives of kinetic energy matrix'
  endif

  ! allocate arrays to store G and G' expansion coefficients

  allocate(coefs_gmat_d(natoms3,natoms3,nterms_gmat_d,x0_npoints), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_gmat_d(natoms3,natoms3,nterms_gmat_d,x0_npoints)', &!
    'natoms3, nterms_gmat_d, x0_npoints =', natoms3, nterms_gmat_d, x0_npoints
    stop
  endif
  coefs_gmat_d = 0.0

  if (allocated(coefs_gmat)) deallocate(coefs_gmat)
  allocate(coefs_gmat(natoms3,natoms3,nterms_gmat,x0_npoints), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_gmat(natoms3,natoms3,nterms_gmat,x0_npoints)', &!
    'natoms3, nterms_gmat, x0_npoints =', natoms3, nterms_gmat, x0_npoints
    stop
  endif
  coefs_gmat = 0.0

  if (present(coefs_dgmat)) then
    if (allocated(coefs_dgmat)) deallocate(coefs_dgmat)
    allocate(coefs_dgmat(natoms3,natoms3,nterms_gmat,x0_npoints,nmodes), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_dgmat(natoms3,natoms3,nterms_gmat,x0_npoints,nmodes)', &!
      'natoms3, nterms_gmat, x0_npoints, nmodes =', natoms3, nterms_gmat, x0_npoints, nmodes
      stop
    endif
    coefs_dgmat = 0.0
  endif

  ! initialize chain-rule expressions for product s-vector * s-vector

  allocate(nterms_chr(nterms_gmat_d), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate nterms_chr(nterms_gmat_d)', 'nterms_gmat_d =', nterms_gmat_d
    stop
  endif
  call chain_rule_product_ind(nmodes, nterms_gmat_d, terms_gmat_d, 2, nterms_chr, ind_chr)

  ! compute derivatives of G-matrix at every reference point

  !$omp parallel do private(ipoint) schedule(dynamic)
  do ipoint=1, x0_npoints
    call deriv_gmat(molec, nterms_svec, terms_svec, coefs_svec(1:natoms3,1:natoms3,1:nterms_svec,ipoint), nterms_gmat_d, terms_gmat_d, nterms_chr, ind_chr, &!
                    coefs_gmat_d(1:natoms3,1:natoms3,1:nterms_gmat_d,ipoint))
  enddo
  !$omp end parallel do

  ! copy derivatives of G

  coefs_gmat = 0.0

  do iterm=1, nterms_gmat
    term0 = terms_gmat(1:nmodes,iterm)
    ipos = index_iarr1(term0, terms_gmat_d)
    if (ipos>0) then
      coefs_gmat(1:natoms3,1:natoms3,iterm,1:x0_npoints) = coefs_gmat_d(1:natoms3,1:natoms3,ipos,1:x0_npoints)
    else
      write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'expand_kinetic error: G-matrix derivative = (', term0(1:nmodes), ') is not initialized'
      stop
    endif
  enddo

  ! copy derivatives of G'

  if (present(coefs_dgmat)) then
    do iterm=1, nterms_gmat
      do imode=1, nmodes
        term0 = terms_gmat(1:nmodes,iterm)
        term0(imode) = term0(imode) + 1
        ipos = index_iarr1(term0, terms_gmat_d)
        if (ipos>0) then
          coefs_dgmat(1:natoms3,1:natoms3,iterm,1:x0_npoints,imode) = coefs_gmat_d(1:natoms3,1:natoms3,ipos,1:x0_npoints)
        else
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'expand_kinetic error: G-matrix derivative = (', term0(1:nmodes), ') is not initialized'
          stop
        endif
      enddo
    enddo
  endif

  ! deallocate chain rule expressions

  deallocate(nterms_chr, ind_chr)

  do iterm=1, nterms_gmat
    fac = 1.0_ark
    do imode=1, nmodes
      !fac = fac * real( exp(factln(terms_gmat(imode,iterm))), kind=ark)
      fac = fac * qfactorial(terms_gmat(imode,iterm))
    enddo
    fac = planck_factor / fac
    coefs_gmat(1:natoms3,1:natoms3,iterm,1:x0_npoints) = coefs_gmat(1:natoms3,1:natoms3,iterm,1:x0_npoints) * fac
    if (present(coefs_dgmat)) then
      coefs_dgmat(1:natoms3,1:natoms3,iterm,1:x0_npoints,1:nmodes) = coefs_dgmat(1:natoms3,1:natoms3,iterm,1:x0_npoints,1:nmodes) * fac
    endif
  enddo


  !-----------------------------------------!
  ! compute derivatives of pseudo-potential !
  !-----------------------------------------!

  if (verbose>=3) then
    write(out, '(/1x,a)') 'compute derivatives of pseudo-potential'
  endif

  ! allocate arrays to store expansion coefficients of pseudo-potential

  if (allocated(coefs_u)) deallocate(coefs_u)
  allocate(coefs_u(nterms_u,x0_npoints), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate coefs_u(nterms_u,x0_npoints)', 'nterms_u, x0_npoints =', nterms_u, x0_npoints
    stop
  endif
  coefs_u = 0.0

  ! initialize chain-rule expressions for product s-vector * s-vector

  allocate(nterms_chr(nterms_u), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate nterms_chr(nterms_u)', 'nterms_u =', nterms_u
    stop
  endif
  call chain_rule_product_ind(nmodes, nterms_u, terms_u, 2, nterms_chr, ind_chr)

  ! compute derivatives of pseudo-potential at every reference point

  if (molec%vib_3n5) then

    coefs_u = 0

  else

    !$omp parallel do private(ipoint) schedule(dynamic)
    do ipoint=1, x0_npoints
      call deriv_pseudo(molec, nterms_svec, terms_svec, coefs_svec(1:natoms3,1:natoms3,1:nterms_svec,ipoint), nterms_u, terms_u, nterms_chr, ind_chr, coefs_u(1:nterms_u,ipoint))
    enddo
    !$omp end parallel do

  endif

  ! deallocate chain rule expressions

  deallocate(nterms_chr, ind_chr)

  do iterm=1, nterms_u
    fac = 1.0_ark
    do imode=1, nmodes
      !fac = fac * real( exp(factln(terms_u(imode,iterm))), kind=ark)
      fac = fac * qfactorial(terms_u(imode,iterm))
    enddo
    fac = planck_factor / fac
    coefs_u(iterm,1:x0_npoints) = coefs_u(iterm,1:x0_npoints) * fac
  enddo


  if (allocated(coefs_cart)) deallocate(coefs_cart)
  if (allocated(terms_cart)) deallocate(terms_cart)

  deallocate(cart0, cart0_ref, rotmat0_ref)
  deallocate(comb_cart, nterms_cart_split, terms_cart_split)
  deallocate(terms_svec, coefs_svec)
  deallocate(terms_gmat_d, coefs_gmat_d)

  if (verbose>=3) then
    write(out, '(/a)') 'expand_kinetic/done'
  endif

end subroutine expand_kinetic


!################################################################################


subroutine expand_func(molec, func, x0_npoints, x0, nmax, ndeg, deg, mindeg, maxdeg, nterms_func, terms_func, coefs_func)

  type(HM_molec_type), intent(in)       :: molec
  type(HM_func_type), intent(in)        :: func
  integer(ik), intent(in)               :: x0_npoints
  real(ark), intent(in)                 :: x0(molec%nmodes,x0_npoints)
  integer(ik), intent(in)               :: nmax
  integer(ik), intent(in)               :: ndeg(molec%nmodes)
  integer(ik), intent(in)               :: deg(maxval(ndeg),molec%nmodes)
  integer(ik), intent(in)               :: mindeg(nmax)
  integer(ik), intent(in)               :: maxdeg(nmax)
  integer(ik), intent(out)              :: nterms_func
  integer(ik), intent(out), allocatable :: terms_func(:,:)
  real(ark), intent(out), allocatable   :: coefs_func(:,:,:)

  integer(ik) :: natoms, nmodes, verbose, imode, iatom, iterm, n, icomb, info, max_order, ncomb, nterms, nterms_cart, ipoint, ipos, max_order_modes(molec%nmodes), rank
  integer(ik), allocatable :: nterms_split(:), terms_split(:,:,:), comb(:,:), terms(:,:), terms_cart(:,:)
  real(ark) :: fac
  real(ark), allocatable :: coefs(:,:,:), coefs_cart(:,:,:,:)
  logical, allocatable :: coefs_func_init(:)

  natoms = molec%natoms
  nmodes = molec%nmodes
  rank = func%rank
  verbose = molec%expfunc_verbose


  if (verbose>=3) then
    write(out, '(/a)') 'expand_func/start: power series expansion of a function'

    write(out, '(/1x,a,<nmodes>(1x,i3))') 'min N-mode expansion orders:', mindeg(1:nmodes)
    write(out, '(1x,a,<nmodes>(1x,i3))') 'max N-mode expansion orders:', maxdeg(1:nmodes)
    write(out, '(1x,a)') 'expansion orders for each mode:'
    do imode=1, nmodes
     write(out, '(1x,i3,5x,100(1x,i3))') imode, deg(1:ndeg(imode),imode)
    enddo

    write(out, '(/1x,a,1x,i6,1x,a)') 'reference points (', x0_npoints, '):'
    do ipoint=1, x0_npoints
      write(out, '(1x,i6,3x,<nmodes>(1x,es16.8))') ipoint, x0(1:nmodes,ipoint)
    enddo
  endif


  call init_coord_transform(molec)
  call init_func(func)


  if (associated(internal_to_cartesian_ptr).eqv..false.) then
    write(out, '(/a)') 'expand_func error: function pointer "internal_to_cartesian_ptr" is not associated'
    stop
  endif


  !---------------------------------!
  ! generate N-mode expansion terms !
  !---------------------------------!

  ! generate derivative terms

  call nmode_expansion(nmodes, nmax, ndeg, deg, mindeg, maxdeg, nterms_func, terms_func)
  call delete_same_terms(nmodes, nterms_func, terms_func)

  if (verbose>=3) then
    write(out, '(/1x,a,1x,i6,1x,a)') 'expansion terms (', nterms_func, ')'
    !do iterm=1, nterms_func
    !  write(out, '(1x,i6,5x,<nmodes>(1x,i3))') iterm, terms_func(1:nmodes,iterm)
    !enddo
  endif

  ! split derivative terms into groups for different combinations of modes

  if (verbose>=3) then
    write(out, '(/1x,a)') 'cast expansion terms into the N-mode expansion form'
  endif

  call nmode_distribute(nmodes, nterms_func, terms_func, ncomb, comb, nterms_split, terms_split)

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6)') 'number of mode-combinations:', ncomb
    do icomb=1, ncomb
      n = count(comb(:,icomb)>0)
      write(out, '(1x,a,1x,i6,3x,a,1x,i6,3x,a,100(1x,i3))') 'icomb:', icomb, 'nterms:', nterms_split(icomb), 'imodes:', comb(1:n,icomb)
    enddo
  endif

  ! allocate array to keep derivative values

  allocate(coefs_func(nterms_func,x0_npoints,func%rank), coefs_func_init(nterms_func), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_func error: failed to allocate coefs_func(nterms_func,x0_npoints,func%rank), coefs_func_init(nterms_func)', &!
    'nterms_func, x0_npoints, func%rank', nterms_func, x0_npoints, func%rank
    stop
  endif


  ! compute derivatives of Cartesian coordinates

  nterms_cart = nterms_func
  allocate(terms_cart(nmodes,nterms_cart), coefs_cart(natoms,3,nterms_cart,x0_npoints), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_func error: failed to allocate terms_cart(nmodes,nterms_cart), coefs_cart(natoms,3,nterms_cart,x0_npoints)', &!
    'nmodes, nterms_cart, natoms, x0_npoints =', nmodes, nterms_cart, natoms, x0_npoints
    stop
  endif
  terms_cart = terms_func
  coefs_cart = 0.0

  if (func%rotate) then
    call expand_cartesian(molec, x0_npoints, x0, nterms_cart, terms_cart, coefs_cart)
  endif


  !--------------------------------!
  ! compute expansion coefficients !
  !--------------------------------!

  if (verbose>=3) then
    write(out, '(/1x,a)') 'compute derivatives'
  endif

  coefs_func = 0.0
  coefs_func_init = .false.

  !$omp parallel do private(icomb,max_order,n,imode,max_order_modes,nterms,terms,coefs,iterm,ipos) schedule(dynamic)
  do icomb=1, ncomb

    if (verbose>=3) then
      write(out, '(1x,a,1x,i6)') 'icomb:', icomb
    endif

    ! estimate total max derivative order, max derivative orders for each mode, and N for current combination of modes

    max_order = maxval(sum(terms_split(1:nmodes,1:nterms_split(icomb),icomb),dim=1))
    n = count(comb(:,icomb)>0)
    max_order_modes = maxval((/( terms_split(comb(imode,icomb),1:nterms_split(icomb),icomb), imode=1, n )/))

    ! estimate number of derivative terms for expansion 0..max_order

    call deriv_func_ADF(molec, func, x0_npoints, x0, max_order, n, comb(1:n,icomb), max_order_modes(1:n), nterms_cart, terms_cart, coefs_cart, nterms, terms, coefs)

    ! keep derivatives for current combination of modes

    do iterm=1, nterms
      ipos = index_iarr1(terms(1:nmodes,iterm), terms_func(1:nmodes,1:nterms_func))
      if (ipos/=0) then
        coefs_func(ipos,1:x0_npoints,1:func%rank) = coefs(iterm,1:x0_npoints,1:func%rank)
        coefs_func_init(ipos) = .true.
      endif
    enddo

    deallocate(terms)
    deallocate(coefs)

  enddo ! icomb
  !$omp end parallel do

  ! check if all derivatives have been calculated

  if (.not.all(coefs_func_init)) then
    write(out, '(/a)') 'expand_func error: following derivatives have not been initialized'
    do iterm=1, nterms_func
      if (.not.coefs_func_init(iterm)) then
        write(out, '(1x,i3,5x,100(1x,i3))') iterm, terms_func(1:nmodes,iterm)
      endif
    enddo
    stop
  endif


  if (verbose>=6) then
    do iterm=1, nterms_func
      write(out, '(/1x,a,1x,<nmodes>(1x,i3),1x,a)') '(', terms_func(1:nmodes,iterm), ') derivative of external function'
      do ipoint=1, x0_npoints
        write(out, '(4x,i6,1x,<nmodes>(1x,es16.8),1x,<rank>(3x,es16.8))') ipoint, x0(1:nmodes,ipoint), coefs_func(iterm,ipoint,1:rank)
      enddo
    enddo
  endif

  ! compute Taylor series expansion coefficients

  do iterm=1, nterms_func
    fac = 1.0_ark
    do imode=1, nmodes
      fac = fac * real( exp(factln(terms_func(imode,iterm))), kind=ark)
    enddo
    coefs_func(iterm,1:x0_npoints,1:func%rank) = coefs_func(iterm,1:x0_npoints,1:func%rank) / fac
  enddo


  deallocate(comb)
  deallocate(nterms_split)
  deallocate(terms_split)
  deallocate(coefs_func_init)
  deallocate(terms_cart, coefs_cart)


  if (verbose>=3) then
    write(out, '(/a)') 'expand_func/done'
  endif

end subroutine expand_func


!###############################################################################


subroutine expand_cartesian(molec, x0_npoints, x0, nterms_cart, terms_cart, coefs_cart)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: x0_npoints
  real(ark), intent(in)           :: x0(molec%nmodes,x0_npoints)
  integer(ik), intent(in)         :: nterms_cart
  integer(ik), intent(in)         :: terms_cart(:,:)
  real(ark), intent(out)          :: coefs_cart(:,:,:,:)

  integer(ik) :: natoms, nmodes, info, icomb, ncomb_cart, n, max_order, nterms, iatom, iterm, imode, ipos, i, ipoint, max_order_modes(molec%nmodes), verbose
  integer(ik), allocatable :: comb_cart(:,:), nterms_cart_split(:), terms_cart_split(:,:,:), terms(:,:), nterms_chr(:), ind_chr(:,:,:)
  real(ark) :: cart0(molec%natoms,3,x0_npoints), cart0_ref(molec%natoms,3,x0_npoints), rotmat0_ref(3,3,x0_npoints)
  real(ark), allocatable :: coefs(:,:,:,:), coefs_rotmat(:,:,:,:), coefs_cart_ref(:,:,:,:)
  logical, allocatable :: coefs_cart_init(:)

  natoms = molec%natoms
  nmodes = molec%nmodes
  verbose = molec%expcart_verbose


  if (verbose>=3) then
    write(out, '(/a)') 'expand_cartesian/start: power series expansion of Cartesian coordinates'
    write(out, '(1x,a,1x,i6,1x,a)') 'expansion terms (', nterms_cart, ')'
    !do iterm=1, nterms_cart
    !  write(out, '(1x,i6,5x,<nmodes>(1x,i3))') iterm, terms_cart(1:nmodes,iterm)
    !enddo
  endif

  call init_coord_transform(molec)


  ! compute reference Cartesian coordinates

  if (associated(internal_to_cartesian_ptr).eqv..false.) then
    write(out, '(/a)') 'expand_cartesian error: function pointer "internal_to_cartesian_ptr" is not associated'
    stop
  endif

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6,1x,a)') 'reference internal coordinates (', x0_npoints, '):'
    do ipoint=1, x0_npoints
      write(out, '(1x,i6,3x,<nmodes>(1x,es16.8))') ipoint, x0(1:nmodes,ipoint)
    enddo
  endif

  cart0_ref = 0.0
  rotmat0_ref = 0.0

  ! compute reference Cartesian coordinates in COORD_TRANSFORM axes system

  do ipoint=1, x0_npoints
    call internal_to_cartesian_ptr(molec, x0(1:nmodes,ipoint), cart0(1:natoms,1:3,ipoint))
  enddo

  if (verbose>=3) then
    write(out, '(1x,a,1x,a,1x,a)') 'Cartesian coordinates in', trim(molec%coord_transform), 'frame'
    do ipoint=1, x0_npoints
      do iatom=1, natoms
        write(out, '(1x,i6,3x,i3,3x,3(f20.16))') ipoint, iatom, cart0(iatom,1:3,ipoint)
      enddo
    enddo
  endif

  ! rotate reference Cartesian coordinates to REFERENCE_FRAME axes system

  if(trim(molec%reference_frame)==trim(molec%coord_transform)) then
    cart0_ref = cart0
    rotmat0_ref = 0.0
    forall(i=1:3) rotmat0_ref(i,i,:) = 1.0_ark
  elseif(trim(molec%reference_frame)=='PAS') then
    if (trim(molec%pas_function)=='EXPKAPPA') then
      do ipoint=1, x0_npoints
        call rotate_pas(molec, cart0(1:natoms,1:3,ipoint), cart0_ref(1:natoms,1:3,ipoint), rotmat0_ref(1:3,1:3,ipoint))
      enddo
    elseif (trim(molec%pas_function)=='LAPACK') then
      do ipoint=1, x0_npoints
        call rotate_pas_diag(molec, cart0(1:natoms,1:3,ipoint), cart0_ref(1:natoms,1:3,ipoint), rotmat0_ref(1:3,1:3,ipoint))
      enddo
    else
      write(out, '(/a,1x,a)') 'expand_cartesian error: unknown principal-axes-system rotation function =', trim(molec%pas_function)
      stop
    endif
#ifdef _MOL_XY4_
  elseif(trim(molec%reference_frame)=='XY4_TD') then
    do ipoint=1, x0_npoints
      call internal_to_cartesian_xy4_equilibrium_td(molec, x0(1:nmodes,ipoint), cart0_ref(1:natoms,1:3,ipoint))
      call rotate_to_xyz(molec, cart0(1:natoms,1:3,ipoint), cart0_ref(1:natoms,1:3,ipoint), rotmat0_ref(1:3,1:3,ipoint))
      cart0_ref(1:natoms,1:3,ipoint) = matmul(cart0(1:natoms,1:3,ipoint),transpose(rotmat0_ref(1:3,1:3,ipoint)))
    enddo
#endif
  else
    write(out, '(/a,a,a)') 'expand_cartesian error: unknown reference frame = "', trim(molec%reference_frame), '"'
    stop
  end if

  if (trim(molec%reference_frame)/=trim(molec%coord_transform)) then
    if (verbose>=3) then
      write(out, '(1x,a,1x,a,1x,a)') 'Cartesian coordinates in', trim(molec%reference_frame), 'frame'
      do ipoint=1, x0_npoints
        do iatom=1, natoms
          write(out, '(1x,i6,3x,i3,3x,3(f20.16))') ipoint, iatom, cart0_ref(iatom,1:3,ipoint)
        enddo
      enddo
    endif
  endif


  ! split expansion terms for Cartesian coordinates into groups for different combinations of modes

  if (verbose>=3) then
    write(out, '(/1x,a)') 'cast expansion terms for Cartesian coordinates into the N-mode expansion form'
  endif

  call nmode_distribute(nmodes, nterms_cart, terms_cart, ncomb_cart, comb_cart, nterms_cart_split, terms_cart_split)

  if (verbose>=3) then
    write(out, '(1x,a,1x,i6)') 'number of mode-combinations:', ncomb_cart
    do icomb=1, ncomb_cart
      n = count(comb_cart(:,icomb)>0)
      write(out, '(1x,a,1x,i6,3x,a,1x,i6,3x,a,100(1x,i3))') 'icomb:', icomb, 'nterms:', nterms_cart_split(icomb), 'imodes:', comb_cart(1:n,icomb)
    enddo
  endif


  ! allocate arrays to keep expansion coefficients for Cartesian coordinates, s-vectors, G-matrix, and pseudo-potential

  allocate(coefs_cart_ref(natoms,3,nterms_cart,x0_npoints), coefs_rotmat(3,3,nterms_cart,x0_npoints), coefs_cart_init(nterms_cart), stat=info)
  if (info/=0) then
    write(out, '(/a,1x,a/a,10(1x,i6))') 'expand_cartesian error: failed to allocate coefs_cart_ref(natoms,3,nterms_cart,x0_npoints)', &!
    'coefs_rotmat(3,3,nterms_cart,x0_npoints), coefs_cart_init(nterms_cart)', 'natoms, nterms_cart, x0_npoints =', natoms, nterms_cart, x0_npoints
    stop
  endif


  !-----------------------------------------------!
  ! compute derivatives for Cartesian coordinates !
  !-----------------------------------------------!

  ! 1. Compute derivatives of Cartesian coordinates in a frame defined by COORD_TRANSFORM function
  ! 2. Compute derivatives of frame rotation matrix (due to Eckart, PAS, etc. conditions)
  ! 3. Apply chain rules to compute derivatives of "rotated" Cartesian coordinates

  if (verbose>=3) then
    write(out, '(/1x,a)') 'compute derivatives of Cartesian coordinates'
  endif

  ! (1) compute derivatives of Cartesian coordinates in COORD_TRANSFORM axes system

  if (verbose>=3) then
    write(out, '(1x,a,1x,a)') trim(molec%coord_transform), 'frame ...'
  endif

  coefs_cart_ref = 0.0
  coefs_cart_init = .false.

  !$omp parallel do private(icomb,max_order,n,max_order_modes,imode,nterms,terms,coefs,iterm,ipoint,ipos) schedule(dynamic)
  do icomb=1, ncomb_cart

    if (verbose>=3) then
      write(out, '(1x,a,1x,i6)') 'icomb:', icomb
    endif

    ! estimate total max derivative order, max derivative orders for each mode, and N for current combination of modes

    max_order = maxval(sum(terms_cart_split(1:nmodes,1:nterms_cart_split(icomb),icomb),dim=1))
    n = count(comb_cart(:,icomb)>0)
    max_order_modes = maxval((/( terms_cart_split(comb_cart(imode,icomb),1:nterms_cart_split(icomb),icomb), imode=1, n )/))

    ! compute derivatives using ADF

    call deriv_cart_ADF(molec, x0_npoints, x0, max_order, n, comb_cart(1:n,icomb), max_order_modes(1:n), nterms, terms, coefs)

    ! rotate to equilibrium-reference axes system

    do ipoint=1, x0_npoints
      do iterm=1, nterms
        coefs(1:natoms,1:3,iterm,ipoint) = matmul(coefs(1:natoms,1:3,iterm,ipoint), transpose(rotmat0_ref(1:3,1:3,ipoint)))
      enddo
    enddo

    ! keep derivatives for current icomb-combination of modes

    do iterm=1, nterms
      ipos = index_iarr1(terms(1:nmodes,iterm), terms_cart(1:nmodes,1:nterms_cart))
      if (ipos/=0) then
        coefs_cart_ref(1:natoms,1:3,ipos,1:x0_npoints) = coefs(1:natoms,1:3,iterm,1:x0_npoints)
        coefs_cart_init(ipos) = .true.
      endif
    enddo

    deallocate(terms)
    deallocate(coefs)

  enddo ! icomb
  !$omp end parallel do

  ! check if all Cartesian-coordinate derivatives have been initialized

  if (.not.all(coefs_cart_init)) then
    write(out, '(/a)') 'expand_cartesian error: following Cartesian-coordinate derivatives have not been initialized'
    do iterm=1, nterms_cart
      if (.not.coefs_cart_init(iterm)) then
        write(out, '(1x,i3,5x,<nmodes>(1x,i3))') iterm, terms_cart(1:nmodes,iterm)
      endif
    enddo
    stop
  endif


  ! (2) compute derivatives of frame rotation matrix

  if (verbose>=3) then
    write(out, '(1x,a,1x,a)') trim(molec%frame_rotation), 'rotation ...'
  endif

  ! initialize chain-rule expressions for product cartesian * rotmat

  allocate(nterms_chr(nterms_cart), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'expand_kinetic error: failed to allocate nterms_chr(nterms_cart)', 'nterms_cart =', nterms_cart
    stop
  endif
  call chain_rule_product_ind(nmodes, nterms_cart, terms_cart, 2, nterms_chr, ind_chr)

  coefs_rotmat = 0.0

  if (trim(molec%frame_rotation)=='ECKART') then

    !$omp parallel do private(ipoint) schedule(dynamic)
    do ipoint=1, x0_npoints
      call deriv_rotmat_eckart(molec, cart0_ref(1:natoms,1:3,ipoint), nterms_cart, terms_cart, nterms_chr, ind_chr, &!
                               coefs_cart_ref(1:natoms,1:3,1:nterms_cart,ipoint), coefs_rotmat(1:3,1:3,1:nterms_cart,ipoint))
    enddo
    !$omp end parallel do

  elseif (trim(molec%frame_rotation)=='I') then

    coefs_rotmat = 0.0
    forall(iterm=1:nterms_cart,i=1:3,all(terms_cart(1:nmodes,iterm)==0)) coefs_rotmat(i,i,iterm,1:x0_npoints) = 1.0_ark

  else

    write(out, '(/a,a,a)') 'expand_kinetic error: unknown molecule-fixed frame embedding = "', trim(molec%frame_rotation), '"'
    stop

  endif

  ! (3) compute derivatives of Cartesian coordinates in FRAME_ROTATION axes system

  if (verbose>=3) then
    write(out, '(1x,a,1x,a)') trim(molec%frame_rotation), 'frame ...'
  endif

  if (trim(molec%frame_rotation)=='I') then

    coefs_cart = coefs_cart_ref

  else

    !$omp parallel do private(ipoint) schedule(dynamic)
    do ipoint=1, x0_npoints
      call deriv_cart_rotated(molec, nterms_cart, terms_cart, nterms_chr, ind_chr, coefs_cart_ref(1:natoms,1:3,1:nterms_cart,ipoint), &!
                              coefs_rotmat(1:3,1:3,1:nterms_cart,ipoint), coefs_cart(1:natoms,1:3,1:nterms_cart,ipoint))
    enddo
    !$omp end parallel do

  endif

  deallocate(nterms_chr, ind_chr)
  deallocate(coefs_rotmat, coefs_cart_init)

  deallocate(comb_cart)
  deallocate(nterms_cart_split)
  deallocate(terms_cart_split)

  if (verbose>=3) then
    write(out, '(/a)') 'expand_cartesian/done'
  endif

end subroutine expand_cartesian


!###############################################################################


subroutine pointwise_kinetic(molec, x0_npoints, x0, gmat, pseudo, dgmat)

  type(HM_molec_type), intent(in)  :: molec
  integer(ik), intent(in)          :: x0_npoints
  real(ark), intent(in)            :: x0(molec%nmodes,x0_npoints)
  real(ark), intent(out)           :: gmat(molec%natoms*3,molec%natoms*3,x0_npoints)
  real(ark), intent(out)           :: pseudo(x0_npoints)
  real(ark), intent(out), optional :: dgmat(molec%natoms*3,molec%natoms*3,x0_npoints,molec%nmodes)

  integer(ik) :: imode, nmodes, nmax_g, nmax_u, ndeg(molec%nmodes), deg(100,molec%nmodes), mindeg_g(molec%nmodes), maxdeg_g(molec%nmodes), &!
                 mindeg_u(molec%nmodes), maxdeg_u(molec%nmodes), nterms_gmat, nterms_u, nterms_gmat_sing, nterms_u_sing, term0(molec%nmodes), ipos
  integer(ik), allocatable :: terms_gmat(:,:), terms_u(:,:), terms_gmat_sing(:,:), terms_u_sing(:,:)
  real(ark), allocatable :: coefs_gmat(:,:,:,:), coefs_dgmat(:,:,:,:,:), coefs_u(:,:), coefs_gmat_sing(:,:,:,:), coefs_u_sing(:,:)

  nmodes = molec%nmodes

  nmax_g = 1
  nmax_u = 1
  ndeg(1:nmodes) = 1
  deg(:,:) = 0
  do imode=1, nmodes
    deg(1:1,imode) = (/0/)
  enddo
  mindeg_g(:) = 0
  mindeg_u(:) = 0
  maxdeg_g(:) = 0
  maxdeg_g(:) = 0
  maxdeg_u(:) = 0
  maxdeg_u(:) = 0

  if (.not.molec%read_keo) then

    if (present(dgmat)) then
      call expand_kinetic(molec, x0_npoints,x0, nmax_g,nmax_u,ndeg,deg(1:maxval(ndeg),1:nmodes),mindeg_g,maxdeg_g,mindeg_u,maxdeg_u, &!
                          nterms_gmat,terms_gmat,coefs_gmat, nterms_u,terms_u,coefs_u, coefs_dgmat)
    else
      call expand_kinetic(molec, x0_npoints,x0, nmax_g,nmax_u,ndeg,deg(1:maxval(ndeg),1:nmodes),mindeg_g,maxdeg_g,mindeg_u,maxdeg_u, &!
                          nterms_gmat,terms_gmat,coefs_gmat, nterms_u,terms_u,coefs_u)
    endif

  else

    if (present(dgmat)) then
      write(out, '(/a)') "pointwise_kinetic error: expansion of G' is not available if you read the KEO from file"
      stop
    endif

    call read_kinetic(molec, nmax_g, maxdeg_g(1), nmax_u, maxdeg_u(1), x0_npoints, x0, nterms_gmat, terms_gmat, coefs_gmat, nterms_u, terms_u, coefs_u)

  endif

  ! output for G-matrix

  term0 = 0
  ipos = index_iarr1(term0, terms_gmat)
  if (ipos<=0) then
    write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'pointwise_kinetic error: G-matrix derivative = (', term0, ') is not initialized'
    stop
  endif
  gmat(:,:,1:x0_npoints) = coefs_gmat(:,:,ipos,1:x0_npoints)
  if (present(dgmat)) then
    dgmat(:,:,1:x0_npoints,1:nmodes) = coefs_dgmat(:,:,ipos,1:x0_npoints,1:nmodes)
  endif

  ! output for pseudo-potential

  term0 = 0
  ipos = index_iarr1(term0, terms_u)
  if (ipos<=0) then
    write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'pointwise_kinetic error: pseudopotential derivative = (', term0, ') is not initialized'
    stop
  endif
  pseudo(1:x0_npoints) = coefs_u(ipos,1:x0_npoints)

  deallocate(terms_gmat, terms_u, coefs_gmat, coefs_u)
  if (allocated(coefs_dgmat)) deallocate(coefs_dgmat)

end subroutine pointwise_kinetic


!###############################################################################


subroutine pointwise_potential(molec, func, x0_npoints, x0, poten, poten2)

  type(HM_molec_type), intent(in)  :: molec
  type(HM_func_type), intent(in)   :: func
  integer(ik), intent(in)          :: x0_npoints
  real(ark), intent(in)            :: x0(molec%nmodes,x0_npoints)
  real(ark), intent(out)           :: poten(x0_npoints)
  real(ark), intent(out), optional :: poten2(molec%nmodes,x0_npoints)

  integer(ik) :: imode, nmodes, nmax, ndeg(molec%nmodes), deg(100,molec%nmodes), mindeg(molec%nmodes), maxdeg(molec%nmodes), nterms_pot, term0(molec%nmodes), ipos
  integer(ik), allocatable :: terms_pot(:,:)
  real(ark), allocatable :: coefs_pot(:,:,:)

  nmodes = molec%nmodes

  nmax = 0
  ndeg(1:nmodes) = 1
  deg = 0
  mindeg(:) = 0
  maxdeg(:) = 0

  if (present(poten2)) then
    nmax = 1
    ndeg(1:nmodes) = 2
    do imode=1, nmodes
      deg(1:2,imode) = (/0,2/)
    enddo
    maxdeg(1:nmax) = 2
  endif

  call expand_func(molec, func, x0_npoints, x0, nmax, ndeg, deg(1:maxval(ndeg),:), mindeg, maxdeg, nterms_pot, terms_pot, coefs_pot)

  term0 = 0
  ipos = index_iarr1(term0, terms_pot)
  if (ipos<=0) then
    write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'pointwise_potential error: potential energy derivative = (', term0, ') is not initialized'
    stop
  endif
  poten(1:x0_npoints) = coefs_pot(ipos,1:x0_npoints,1)

  if (present(poten2)) then
    do imode=1, nmodes
      term0 = 0
      term0(imode) = 2
      ipos = index_iarr1(term0, terms_pot)
      if (ipos<=0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'pointwise_potential error: potential energy derivative = (', term0, ') is not initialized'
        stop
      endif
      poten2(imode,1:x0_npoints) = coefs_pot(ipos,1:x0_npoints,1)
    enddo
  endif

  deallocate(terms_pot, coefs_pot)

end subroutine pointwise_potential


!###############################################################################


subroutine expand_kinetic_1d(molec, x0_npoints, x0, imode, ndeg_imode, deg_imode, gmat, pseudo)

  type(HM_molec_type), intent(in) :: molec
  integer(ik), intent(in)         :: x0_npoints
  real(ark), intent(in)           :: x0(molec%nmodes,x0_npoints)
  integer(ik), intent(in)         :: imode, ndeg_imode, deg_imode(ndeg_imode)
  real(ark), intent(out)          :: gmat(molec%natoms*3,molec%natoms*3,x0_npoints,ndeg_imode)
  real(ark), intent(out)          :: pseudo(x0_npoints,ndeg_imode)

  integer(ik) :: nmodes, nmax_g, nmax_u, ndeg(molec%nmodes), deg(100,molec%nmodes), mindeg_g(molec%nmodes), maxdeg_g(molec%nmodes), i, &!
                 mindeg_u(molec%nmodes), maxdeg_u(molec%nmodes), nterms_gmat, nterms_u, term0(molec%nmodes), ipos, ideg
  integer(ik), allocatable :: terms_gmat(:,:), terms_u(:,:)
  real(ark), allocatable :: coefs_gmat(:,:,:,:), coefs_u(:,:)

  nmodes = molec%nmodes

  ndeg(1:nmodes) = 1
  deg(:,:) = 0
  mindeg_g = 0
  maxdeg_g = 0
  mindeg_u = 0
  maxdeg_u = 0

  nmax_g = 1
  nmax_u = 1
  ndeg(imode) = maxval(deg_imode(1:ndeg_imode))+1
  deg(1:ndeg(imode),imode) = (/( i, i=0,ndeg(imode)-1 )/)
  mindeg_g(1) = minval(deg_imode(1:ndeg_imode))
  mindeg_u(1) = minval(deg_imode(1:ndeg_imode))
  maxdeg_g(1) = maxval(deg_imode(1:ndeg_imode))
  maxdeg_u(1) = maxval(deg_imode(1:ndeg_imode))

  call expand_kinetic(molec, x0_npoints, x0, nmax_g, nmax_u, ndeg, deg(1:maxval(ndeg),1:nmodes), mindeg_g, maxdeg_g, mindeg_u, maxdeg_u, &!
                      nterms_gmat, terms_gmat, coefs_gmat, nterms_u, terms_u, coefs_u)

  ! output for G-matrix

  do ideg=1, ndeg_imode
    term0 = 0
    term0(imode) = deg_imode(ideg)
    ipos = index_iarr1(term0, terms_gmat)
    if (ipos<=0) then
      write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'expand_kinetic_1d error: G-matrix derivative = (', term0, ') is not initialized'
      stop
    endif
    gmat(:,:,1:x0_npoints,ideg) = coefs_gmat(:,:,ipos,1:x0_npoints)
  enddo

  ! output for pseudopotential

  do ideg=1, ndeg_imode
    term0 = 0
    term0(imode) = deg_imode(ideg)
    ipos = index_iarr1(term0, terms_u)
    if (ipos<=0) then
      write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') 'expand_kinetic_1d error: pseudopotential derivative = (', term0, ') is not initialized'
      stop
    endif
    pseudo(1:x0_npoints,ideg) = coefs_u(ipos,1:x0_npoints)
  enddo

  deallocate(terms_gmat, terms_u, coefs_gmat, coefs_u)

end subroutine expand_kinetic_1d


!###############################################################################


end module hamiltonian
