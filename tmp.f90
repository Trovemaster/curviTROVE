module hamiltonian
use accuracy
use fields
include 'RAinclude.i90'
use higherOrderTensorUtil
implicit none


interface
  subroutine internal_to_cartesian_RA(natoms, nmodes, masses, x0_RA, cartesian_RA)
    use accuracy
    include 'RAinclude.i90'
    integer(ik), intent(in)   :: nmodes, natoms
    real(rk), intent(in)      :: masses(natoms)
    type(RARealD), intent(in) :: x0_RA(nmodes)
    type(RARealD)             :: cartesian_RA(natoms,3)
  end subroutine internal_to_cartesian_RA
end interface

interface
  subroutine internal_to_cartesian(natoms, nmodes, masses, x0, cartesian)
    use accuracy
    integer(ik), intent(in) :: nmodes, natoms
    real(rk), intent(in)    :: masses(natoms)
    real(rk), intent(in)    :: x0(nmodes)
    real(rk)                :: cartesian(natoms,3)
  end subroutine internal_to_cartesian
end interface

interface
  subroutine internal_to_cartesian_ark(natoms, nmodes, masses, x0, cartesian)
    use accuracy
    integer(ik), intent(in) :: nmodes, natoms
    real(ark), intent(in)   :: masses(natoms)
    real(ark), intent(in)   :: x0(nmodes)
    real(ark)               :: cartesian(natoms,3)
  end subroutine internal_to_cartesian_ark
end interface

procedure(internal_to_cartesian_RA), pointer :: internal_to_cartesian_RA_ptr
procedure(internal_to_cartesian), pointer :: internal_to_cartesian_ptr
procedure(internal_to_cartesian_ark), pointer :: internal_to_cartesian_ark_ptr


contains


include 'coord_transform.f90'
include 'mol_xy2.f90'



subroutine HAinitilize_Kinetic

  integer(ik), parameter :: max_order=60
  integer(ik) :: natoms, nmodes, kin_nmax, pot_nmax, ppot_nmax, fdf_npoints
  integer(ik) :: rot_maxit, rot_maxorder
  integer(ik), allocatable :: modes_order(:,:), modes_norder(:), kin_order(:), pot_order(:), ppot_order(:)
  real(rk), allocatable :: masses(:), x0(:), fdf_step(:)
  real(rk) :: rot_accur
  character(cl) :: deriv_method, coord_transform, frame_rotation

  integer(ik), parameter :: max_nwords=100
  integer(ik) :: nwords, info, iatom, imode, jmode, i, ndeg
  real(rk) :: fdf_step_
  character(500) buf500
  character(cl) :: word(max_nwords)

  natoms = 0
  nmodes = 0
  kin_nmax = 0
  pot_nmax = 0
  ppot_nmax = 0
  deriv_method = 'ADF-NUM'
  coord_transform = 'UNDEFINED'
  frame_rotation = 'ECKART'
  fdf_npoints = 3
  fdf_step_ = 0.001_rk
  rot_maxit = 30
  rot_maxorder = 30
  rot_accur = 1.0d-14

  rewind(inp)
  do
    read(inp,'(a)',iostat=info) buf500
    if (info/=0) exit
    call tokenize_str(buf500, nwords, word)
    select case(trim(to_upper(word(1))))
      !
    case('NATOMS') ! number of atoms
      !
      read(word(2),*) natoms
      !
    case('NMODES') ! number of modes
      !
      read(word(2),*) nmodes
      !
      kin_nmax = nmodes
      pot_nmax = nmodes
      ppot_nmax = nmodes
      !
    case('ZMAT') ! masses
      !
      if (natoms==0) then
        write(out, '(/a)') 'HMinitilize_Kinetic error: "NATOMS" keyword must appear before "ZMAT"'
        stop
      endif
      allocate(masses(natoms), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate masses(natoms)', 'natoms =', natoms
        stop
      endif
      masses = 0.0
      do iatom=1, natoms
        read(inp,'(a)') buf500
        call tokenize_str(buf500, nwords, word)
        if (trim(to_upper(word(nwords)))=='END') then
          write(out, '(/a)') 'HMinitialize_Kinetic error: number of rows in "ZMAT" section is less than the number of atoms'
          stop
        endif
        read(word(nwords),*) masses(iatom)
      enddo
      !
      if (any(masses<1.0)) then
        write(out, '(/a/a,100(1x,es16.8))') 'HMinitilize_Kinetic error: too small atomic masses', &
        'masses =', masses
        stop
      endif
      !
    case('EQUILIBRIUM') ! equilibrium values of internal coordinates
      !
      if (nmodes==0) then
        write(out, '(/a)') 'HMinitilize_Kinetic error: "NMODES" keyword must appear before "EQUILIBRIUM"'
        stop
      endif
      allocate(x0(nmodes), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate x0(nmodes)', 'nmodes =', nmodes
        stop
      endif
      x0 = 0.0
      do imode=1, nmodes
        read(inp,'(a)') buf500
        call tokenize_str(buf500, nwords, word)
        if (trim(to_upper(word(nwords)))=='DEG') then
          read(word(nwords-1),*) x0(imode)
          x0(imode) = x0(imode) * real(pi,rk)/180.0_rk
        elseif (trim(to_upper(word(nwords)))=='END') then
          write(out, '(/a)') 'HMinitialize_Kinetic error: number of rows in "EQUILIBRIUM" section is less than the number of modes'
          stop
        else
          read(word(nwords),*) x0(imode)
        endif
      enddo
      !
    case('TRANSFORM') ! internal-to-Cartesian coordinate transformation
      !
      coord_transform = trim(to_upper(word(2)))
      !
    case('EMBEDDING') ! Cartesian-coordinate frame rotation
      !
      frame_rotation = trim(to_upper(word(2)))
      !
      if (.not.any((/'ECKART','PAS','I'/)==trim(frame_rotation))) then
        write(out, '(/a,1x,a)') 'HMinitilize_Kinetic error: unknown frame rotation method =', trim(frame_rotation)
        stop
      endif
      !
    case('DERIV_METHOD') ! method for derivatives
      !
      deriv_method = trim(to_upper(word(2)))
      !
      if (.not.any((/'FDF8','FDF16','FDF8-NUM','FDF16-NUM','ADF-NUM'/)==trim(deriv_method))) then
        write(out, '(/a,1x,a)') 'HMinitilize_Kinetic error: unknown derivative method =', trim(deriv_method)
        stop
      endif
      !
    case('FDF_STEP')
      !
      if (nmodes==0) then
        write(out, '(/a)') 'HMinitilize_Kinetic error: "NMODES" keyword must appear before "FDFSTEP"'
        stop
      endif
      allocate(fdf_step(nmodes), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate fdf_step(nmodes)', 'nmodes =', nmodes
        stop
      endif
      fdf_step = fdf_step_
      do i=2, nwords
        read(word(i),*) fdf_step(i-1)
      enddo
      do i=min(nmodes,nwords)+1, nmodes+1
        fdf_step(i-1) = fdf_step(i-2)
      enddo
      !
      if (any(fdf_step<5.0d-4)) then
        write(out, '(/a/a,100(1x,f10.6))') 'HMinitilize_Kinetic error: input finite-difference step size is too small (<0.0005)', &
       'fdf_step =', fdf_step
        stop
      endif
      !
    case('FDF_NPOINTS')
      !
      read(word(2),*) fdf_npoints
      !
      if (.not.any((/3,5,7,9/)==fdf_npoints)) then
        write(out, '(/a,1x,i3,1x,a)') 'HMinitilize_Kinetic error: wrong number of fdf points =', fdf_npoints, '(only 3, 5, 7, or 9)'
        stop
      endif
      !
    case('KINNMAX')
      !
      read(word(2),*) kin_nmax
      kin_nmax = min(nmodes,kin_nmax)
      !
    case('KINORDER')
      !
      if (nmodes==0) then
        write(out, '(/a)') 'HMinitilize_Kinetic error: "NMODES" keyword must appear before "KINORDER"'
        stop
      endif
      allocate(kin_order(nmodes), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate kin_order(nmodes)', 'nmodes =', nmodes
        stop
      endif
      kin_order = 0
      do i=2, min(nmodes,nwords)
        read(word(i),*) kin_order(i-1)
      enddo
      do i=min(nmodes,nwords)+1, nmodes+1
        kin_order(i-1) = kin_order(i-2)
      enddo
      !
      if (any(kin_order<0).or.any(kin_order>20)) then
        write(out, '(/a/a,100(1x,i3))') 'HMinitilize_Kinetic error: negative or too large expansion orders for KEO', &
        '"KINORDER"=', kin_order
        stop
      endif
      !
    case('POTNMAX')
      !
      read(word(2),*) pot_nmax
      pot_nmax = min(nmodes,pot_nmax)
      !
    case('POTORDER')
      !
      if (nmodes==0) then
        write(out, '(/a)') 'HMinitilize_Kinetic error: "NMODES" keyword must appear before "POTORDER"'
        stop
      endif
      allocate(pot_order(nmodes), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate pot_order(nmodes)', 'nmodes =', nmodes
        stop
      endif
      pot_order = 0
      do i=2, min(nmodes,nwords)
        read(word(i),*) pot_order(i-1)
      enddo
      do i=min(nmodes,nwords)+1, nmodes+1
        pot_order(i-1) = pot_order(i-2)
      enddo
      !
      if (any(pot_order<0).or.any(pot_order>20)) then
        write(out, '(/a/a,100(1x,i3))') 'HMinitilize_Kinetic error: negative or too large expansion orders for PES', &
        '"POTORDER"=', pot_order
        stop
      endif
      !
    case('MODESORDER')
      !
      if (nmodes==0) then
        write(out, '(/a)') 'HMinitilize_Kinetic error: "NMODES" keyword must appear before "MODESORDER"'
        stop
      endif
      allocate(modes_order(max_order,nmodes), modes_norder(nmodes), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate modes_order(max_order,nmodes), modes_norder(nmodes)', &
        'max_order, nmodes =', max_order, nmodes
        stop
      endif
      modes_order = 0
      modes_norder = 1
      do imode=1, nmodes
        read(inp,'(a)') buf500
        call tokenize_str(buf500, nwords, word)
        if (trim(to_upper(word(1)))=='END') exit
        read(word(1),*) jmode
        do i=2, nwords
          read(word(i),*) modes_order(i-1,jmode)
        enddo
        modes_norder(jmode) = nwords-1
      enddo
      !
    end select
  enddo


  if (nmodes<=0) then
    write(out, '(/a)') 'HMinitilize_Kinetic error: wrong number of modes ("NMODES" keyword)'
    stop
  endif

  if (natoms<=0) then
    write(out, '(/a)') 'HMinitilize_Kinetic error: wrong number of atoms ("NATOMS" keyword)'
    stop
  endif

  if (.not.allocated(masses)) then
    write(out, '(/a)') 'HMinitilize_Kinetic error: atomic masses are not initialized ("ZMAT" keyword)'
    stop
  endif

  if (.not.allocated(x0)) then
    write(out, '(/a)') 'HMinitilize_Kinetic error: equilibrium coordinates are not initialized ("EQUILIBRIUM" keyword)'
    stop
  endif

  if (trim(coord_transform)=='UNDEFINED') then
    write(out, '(/a)') 'HMinitilize_Kinetic error: internal-to-Cartesian transformation is not initialized ("TRANSFORM" keyword)'
    stop
  endif

  if (.not.allocated(kin_order)) then
    write(out, '(/a)') 'HMinitilize_Kinetic error: KEO expansion is not initialized ("KINORDER" keyword)'
    stop
  endif

  if (.not.allocated(pot_order)) then
    write(out, '(/a)') 'HMinitilize_Kinetic error: PES expansion is not initialized ("POTORDER" keyword)'
    stop
  endif

  if (.not.allocated(ppot_order)) then
    allocate(ppot_order(nmodes), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate ppot_order(nmodes)', 'nmodes =', nmodes
      stop
    endif
    ppot_order = kin_order
  endif

  if (.not.allocated(modes_order)) then
    allocate(modes_order(max_order,nmodes), modes_norder(nmodes), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate modes_order(max_order,nmodes), modes_norder(nmodes)', &
      'max_order, nmodes =', max_order, nmodes
      stop
    endif
    modes_order = 0
    modes_norder = 1
    do imode=1, nmodes
      ndeg = max(maxval(pot_order), maxval(kin_order), maxval(ppot_order))
      do i=0, ndeg
        modes_order(i+1,imode) = i
      enddo
      modes_norder(imode) = ndeg+1
    enddo
  endif

  if (.not.allocated(fdf_step)) then
    allocate(fdf_step(nmodes), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') 'HMinitilize_Kinetic error: failed to allocate fdf_step(nmodes)', 'nmodes =', nmodes
      stop
    endif
    fdf_step = fdf_step_
  endif

  write(out, '(1x,a,1x,i3)') 'number of modes:', nmodes
  write(out, '(1x,a,1x,i3)') 'number of atoms:', natoms
  write(out, '(1x,a)') 'atomic masses:'
  do iatom=1, natoms
    write(out, '(1x,i3,1x,es16.8)') iatom, masses(iatom)
  enddo
  write(out, '(1x,a)') 'equilibrium coordinates:'
  do imode=1, nmodes
    write(out, '(1x,i3,1x,es16.8)') imode, x0(imode)
  enddo
  write(out, '(1x,a,1x,a)') 'internal-to-Cartesian transformation:', trim(coord_transform)
  write(out, '(1x,a,1x,a)') 'coordnate frame rotation:', trim(frame_rotation)
  write(out, '(1x,a,1x,a)') 'derivatives'' method:', trim(deriv_method)
  if (index(deriv_method,'FDF')/=0) then
    write(out, '(1x,a,1x,i3)') 'number of finite-difference points:', fdf_npoints
    write(out, '(1x,a,100(1x,f10.6))') 'finite-difference step size:', fdf_step(:)
  endif

  write(out, '(1x,a,100(1x,i3))') 'n-mode max expansion orders for KEO:', kin_order(1:kin_nmax)
  write(out, '(1x,a,100(1x,i3))') 'n-mode max expansion orders for PES:', pot_order(1:pot_nmax)
  write(out, '(1x,a,100(1x,i3))') 'expansion orders for each mode:'
  do imode=1, nmodes
    write(out, '(1x,i3,1x,a,100(1x,i3))') imode, ':', modes_order(1:modes_norder(imode),imode)
  enddo

  call expand_gmat(natoms, nmodes, masses, x0, kin_nmax, modes_norder, modes_order(1:maxval(modes_norder),1:nmodes), &
                   (/(0, i=1, kin_nmax)/), kin_order(1:kin_nmax), deriv_method, coord_transform, &
                   frame_rotation, fdf_npoints, fdf_step, rot_maxorder, rot_maxit, rot_accur)

end subroutine HAinitilize_Kinetic



subroutine expand_gmat(natoms, nmodes, masses, x0, nmax, ndeg, deg, mindeg, maxdeg, deriv_method, coord_transform, &
                       frame_rotation, fdf_npt, fdf_step, rot_maxorder, rot_maxit, rot_accur)

  integer(ik), intent(in)  :: natoms
  integer(ik), intent(in)  :: nmodes
  real(rk), intent(in)     :: masses(natoms)
  real(rk), intent(in)     :: x0(nmodes)
  integer(ik), intent(in)  :: nmax
  integer(ik), intent(in)  :: ndeg(nmodes)
  integer(ik), intent(in)  :: deg(maxval(ndeg),nmodes)
  integer(ik), intent(in)  :: mindeg(nmax)
  integer(ik), intent(in)  :: maxdeg(nmax)
  character(*), intent(in) :: deriv_method
  character(*), intent(in) :: coord_transform
  character(*), intent(in) :: frame_rotation
  integer(ik), intent(in)  :: fdf_npt
  real(rk), intent(in)     :: fdf_step(nmodes)
  integer(ik), intent(in)  :: rot_maxorder
  integer(ik), intent(in)  :: rot_maxit
  real(rk), intent(in)     :: rot_accur

  integer(ik) :: info, nterms_gmat, nterms_svec, nterms_cart, icomb, ncomb_cart, n, iorder, max_order, iterm1, iterm2, &!
                 nterms, iatom, iterm, imode, tot_nterms, nterms_cart_icomb, ipos, i
  integer(ik), allocatable :: terms_gmat(:,:), terms_svec(:,:), terms_cart(:,:), comb_cart(:,:), nterms_cart_split(:), &!
                              terms_cart_split(:,:,:), terms(:,:), tmp_terms(:,:), terms_cart_icomb(:,:)
  real(rk) :: cart0(natoms,3), cart0_pas(natoms,3), fac
  real(rk), allocatable :: coefs(:,:,:), tmp_coefs(:,:,:), coefs_cart_sf(:,:,:), coefs_rotmat(:,:,:), coefs_cart_icomb(:,:,:), &!
                           coefs_cart(:,:,:), coefs_svec(:,:,:), coefs_gmat(:,:,:)
  logical, allocatable :: coefs_cart_init(:)

  ! init internal-to-Cartesian transformation function
  select case(trim(coord_transform))
  case('XY2_RALPHA')
    internal_to_cartesian_RA_ptr => internal_to_cartesian_xy2_ralpha_RA
    internal_to_cartesian_ptr => internal_to_cartesian_xy2_ralpha
    internal_to_cartesian_ark_ptr => internal_to_cartesian_xy2_ralpha_ark
  case default
    write(out, '(/a,1x,a)') 'expand_gmat error: unknown internal-to-Cartesian coordinate transformation =', trim(coord_transform)
    stop
  end select

  ! equilibrium Cartesian coordinates of atoms
  call internal_to_cartesian_ptr(natoms, nmodes, masses, x0, cart0)
  call rotate_pas(natoms, masses, cart0, rot_maxorder, rot_maxit, rot_accur, cart0_pas)
  write(out, '(1x,a)') 'Cartesian coordinates of atoms:'
  do iatom=1, natoms
    write(out, '(1x,i3,3x,3(1x,f20.16))') iatom, cart0(iatom,1:3)
  enddo
  write(out, '(1x,a)') 'Cartesian coordinates of atoms in PAS:'
  do iatom=1, natoms
    write(out, '(1x,i3,3x,3(1x,f20.16))') iatom, cart0_pas(iatom,1:3)
  enddo

  ! generate G-matrix derivative terms
  write(out, '(/1x,a)') 'generate G-matrix derivative terms ...'
  call nmode_expansion(nmodes, nmax, ndeg, deg, mindeg, maxdeg, nterms_gmat, terms_gmat)
  write(out, '(1x,a,1x,i8)') 'number of terms:', nterms_gmat
  do iterm=1, nterms_gmat
    write(out, '(1x,i3,5x,100(1x,i3))') iterm, terms_gmat(1:nmodes,iterm)
  enddo

  ! generate complete set of derivative terms required for recursive determination of all s-vector derivatives
  write(out, '(/1x,a)') 'generate s-vector derivative terms ...'
  nterms_svec = nterms_gmat
  allocate(terms_svec(nmodes,nterms_svec), coefs_svec(natoms*3,natoms*3,nterms_svec), stat=info)
  if (info/=0) then
    write(out, '(/a/a,3(1x,i6))') &
    'expand_gmat error: failed to allocate terms_svec(nmodes,nterms_svec), coefs_svec(natoms*3,natoms*3,nterms_svec)', &!
    'nmodes, natoms, nterms_svec =', nmodes, natoms, nterms_svec
    stop
  endif
  coefs_svec = 0.0
  terms_svec = terms_gmat
  call expansion_terms_complete(nmodes, nterms_svec, terms_svec)
  write(out, '(1x,a,1x,i8)') 'number of terms:', nterms_svec
  do iterm=1, nterms_svec
    write(out, '(1x,i3,5x,100(1x,i3))') iterm, terms_svec(1:nmodes,iterm)
  enddo

  ! generate Cartesian-coordinate derivative terms
  write(out, '(/1x,a)') 'generate Cartesian-coordinate derivative terms ...'
  nterms_cart = nterms_svec
  allocate(terms_cart(nmodes,nterms_cart), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &
    'expand_gmat error: failed to allocate terms_cart(nmodes,nterms_cart)', &!
    'nmodes, nterms_cart =', nmodes, nterms_cart
    stop
  endif
  terms_cart = terms_svec
  call expansion_terms_svec2cart(nmodes, nterms_cart, terms_cart)
  write(out, '(1x,a,1x,i8)') 'number of terms:', nterms_cart
  do iterm=1, nterms_cart
    write(out, '(1x,i3,5x,100(1x,i3))') iterm, terms_cart(1:nmodes,iterm)
  enddo

  ! allocate array to keep Cartesian-coordinate derivatives
  allocate(coefs_cart(natoms,3,nterms_cart), coefs_cart_init(nterms_cart), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &
    'expand_gmat error: failed to allocate coefs_cart(natoms,3,nterms_cart), coefs_cart_init(nterms_cart)', &!
    'natoms, nterms_cart', natoms, nterms_cart
    stop
  endif
  coefs_cart = 0.0
  coefs_cart_init = .false.

  ! split Cartesian-coordinate derivative terms into groups for different combinations of modes;
  ! this is necessary to parallelize the computations of Cartesian-coordinate derivatives
  write(out, '(/1x,a)') 'cast Cartesian-coordinate derivative terms into n-mode expansion form ...'
  call nmode_distribute(nmodes, nterms_cart, terms_cart, ncomb_cart, comb_cart, nterms_cart_split, terms_cart_split)
  write(out, '(1x,a,1x,i3)') 'number of mode-combinations:', ncomb_cart
  do icomb=1, ncomb_cart
    n = count(comb_cart(:,icomb)>0)
    write(out, '(1x,a,1x,i3,3x,a,1x,i6,3x,a,100(1x,i3))') &
   'icomb:', icomb, 'nterms:', nterms_cart_split(icomb), 'imodes:', comb_cart(1:n,icomb)
  enddo

  write(out, '(/1x,a)') 'compute Cartesian-coordinate derivatives ...'

  ! compute derivatives of Cartesian coordinates independently for each combination of modes
  do icomb=1, ncomb_cart

    ! max derivative order and N for current combination of modes
    max_order = maxval(sum(terms_cart_split(1:nmodes,1:nterms_cart_split(icomb),icomb),dim=1))
    n = count(comb_cart(:,icomb)>0)

    ! estimate number of derivative terms for expansion 0..max_order
    tot_nterms = 1
    do iorder=1, max_order
      tot_nterms = tot_nterms + bico(n+iorder-1,iorder)
    enddo

    ! allocate arrays to keep derivatives of Cartesian coordinates,
    allocate(terms_cart_icomb(nmodes,tot_nterms), coefs_cart_icomb(natoms,3,tot_nterms), stat=info)
    if (info/=0) then
      write(out, '(/a/a,2(1x,i6))') &
      'expand_gmat error: failed to allocate terms_cart_icomb(nmodes,tot_nterms), coefs_cart_icomb(natoms,3,tot_nterms)', &
      'nmodes, tot_nterms =', nmodes, tot_nterms
      stop
    endif

    select case(trim(deriv_method))

    ! compute derivatives using finite-differences
    case('FDF8','FDF16')

      terms_cart_icomb = 0
      coefs_cart_icomb = 0.0

      iterm2 = 0
      do iorder=0, max_order
        select case(trim(deriv_method))
        case('FDF16')
          ! finite-difference method (quadrupole precision)
          call deriv_cart_fdf_ark(natoms, nmodes, masses, x0, cart0_pas, iorder, n, comb_cart(1:n,icomb), fdf_npt, &
                                  fdf_step, frame_rotation, rot_maxorder, rot_maxit, rot_accur, nterms, terms, coefs)
        case('FDF8')
          ! finite-difference method
          call deriv_cart_fdf(natoms, nmodes, masses, x0, cart0_pas, iorder, n, comb_cart(1:n,icomb), fdf_npt, &
                              fdf_step, frame_rotation, rot_maxorder, rot_maxit, rot_accur, nterms, terms, coefs)
        case default
          write(out, '(/a,1x,a)') 'expand_gmat error: unknown derivative method =', trim(deriv_method)
          stop
        end select
        iterm1 = iterm2 + 1
        iterm2 = iterm1 + nterms - 1
        if (iterm2>tot_nterms) then
          write(out, '(/a,1x,i6,1x,a,1x,<n>(1x,i3),1x,a)') 'expand_gmat error: estimated number of expansion terms =', &
          tot_nterms, 'for mode-combination = (', comb_cart(1:n,icomb), ') is exceeded'
          stop
        endif
        terms_cart_icomb(1:nmodes,iterm1:iterm2) = terms
        coefs_cart_icomb(1:natoms,1:3,iterm1:iterm2) = coefs
      enddo

      nterms_cart_icomb = iterm2

      if (allocated(terms)) deallocate(terms)
      if (allocated(coefs)) deallocate(coefs)

    ! compute derivatives using finite-differences/automatic differentiation
    ! combined with numerical approach for coordinate frame rotation matrix
    case('FDF8-NUM','FDF16-NUM','ADF-NUM')

      ! allocate arrays to keep derivatives of space-fixed Cartesian coordinates and coordinate frame rotation matrix
      allocate(coefs_cart_sf(natoms,3,tot_nterms), coefs_rotmat(3,3,tot_nterms), stat=info)
      if (info/=0) then
        write(out, '(/a/a,2(1x,i6))') &
        'expand_gmat error: failed to allocate coefs_cart_sf(natoms,3,tot_nterms), coefs_rotmat(3,3,tot_nterms)', &
        'nmodes, tot_nterms =', nmodes, tot_nterms
        stop
      endif

      ! compute derivatives of space-fixed Cartesian coordinates
      terms_cart_icomb = 0
      coefs_cart_sf = 0.0
      iterm2 = 0
      do iorder=0, max_order
        select case(deriv_method)
        case('ADF-NUM')
          ! automatic differentiation method
          call deriv_cart_rapsodia(natoms, nmodes, masses, x0, iorder, n, comb_cart(1:n,icomb), nterms, terms, coefs)
        case('FDF8-NUM')
          ! finite-difference method
          call deriv_cart_fdf(natoms, nmodes, masses, x0, cart0_pas, iorder, n, comb_cart(1:n,icomb), fdf_npt, &
                              fdf_step, 'I', rot_maxorder, rot_maxit, rot_accur, nterms, terms, coefs)
        case('FDF16-NUM')
          ! finite-difference method (quadrupole precision)
          call deriv_cart_fdf_ark(natoms, nmodes, masses, x0, cart0_pas, iorder, n, comb_cart(1:n,icomb), fdf_npt, &
                                  fdf_step, 'I', rot_maxorder, rot_maxit, rot_accur, nterms, terms, coefs)
        case default
          write(out, '(/a,1x,a)') 'expand_gmat error: unknown derivative method =', trim(deriv_method)
          stop
        end select
        iterm1 = iterm2 + 1
        iterm2 = iterm1 + nterms - 1
        if (iterm2>tot_nterms) then
          write(out, '(/a,1x,i6,1x,a,1x,<n>(1x,i3),1x,a)') 'expand_gmat error: estimated number of expansion terms =', &
          tot_nterms, 'for mode-combination = (', comb_cart(1:n,icomb), ') is exceeded'
          stop
        endif
        terms_cart_icomb(1:nmodes,iterm1:iterm2) = terms
        coefs_cart_sf(1:natoms,1:3,iterm1:iterm2) = coefs
      enddo

      nterms_cart_icomb = iterm2

      if (allocated(terms)) deallocate(terms)
      if (allocated(coefs)) deallocate(coefs)

      ! compute derivatives of coordinate frame rotation matrix
      coefs_rotmat = 0.0
      select case(trim(frame_rotation))
      case('ECKART')
        call deriv_rotmat_eckart(natoms, nmodes, masses, cart0_pas, nterms_cart_icomb, terms_cart_icomb, &
                                 coefs_cart_sf, rot_maxorder, rot_maxit, rot_accur, coefs_rotmat)
      case('I')
        coefs_rotmat = 0.0
        forall(iterm=1:nterms_cart_icomb,i=1:3,all(terms_cart_icomb(:,iterm)==0)) coefs_rotmat(i,i,iterm) = 1.0_rk
      case default
        write(out, '(/a,1x,a)') 'expand_gmat error: unknown coordinate frame rotation function =', trim(frame_rotation)
        stop
      end select

      ! compute derivatives of body-fixed Cartesian coordinates
      coefs_cart_icomb = 0.0
      call deriv_cart_rotated(natoms, nmodes, nterms_cart_icomb, terms_cart_icomb, coefs_cart_sf, coefs_rotmat, &
                              coefs_cart_icomb)

    case default
      write(out, '(/a,1x,a)') 'expand_gmat error: unknown derivative method =', trim(deriv_method)
      stop

    end select

    ! save derivatives for current combination of modes
    do iterm=1, nterms_cart_icomb
      ipos = index_iarr1(terms_cart_icomb(1:nmodes,iterm), terms_cart)
      if (ipos/=0) then
        coefs_cart(1:natoms,1:3,ipos) = coefs_cart_icomb(1:natoms,1:3,iterm)
        coefs_cart_init(ipos) = .true.
      endif
    enddo

    deallocate(coefs_cart_icomb)
    deallocate(terms_cart_icomb)
    if (allocated(coefs_cart_sf)) deallocate(coefs_cart_sf)
    if (allocated(coefs_rotmat)) deallocate(coefs_rotmat)

  enddo ! icomb

  if (allocated(comb_cart)) deallocate(comb_cart)
  if (allocated(nterms_cart_split)) deallocate(nterms_cart_split)
  if (allocated(terms_cart_split)) deallocate(terms_cart_split)

  ! check if all Cartesian-coordinate derivatives are initialized
  if (.not.all(coefs_cart_init)) then
    write(out, '(/a)') 'expand_gmat error: following Cartesian-coordinate derivatives have not been initialized'
    do iterm=1, nterms_cart
      if (.not.coefs_cart_init(iterm)) then
        write(out, '(1x,i3,5x,100(1x,i3))') iterm, terms_cart(:,iterm)
      endif
    enddo
    stop
  endif

  deallocate(coefs_cart_init)

  do iterm=1, nterms_cart
    write(out, '(100(1x,i3))') terms_cart(:,iterm)
    do iatom=1, natoms
      write(out, '(30x,100(1x,es16.8))') coefs_cart(iatom,1:3,iterm)
    enddo
  enddo
stop
  ! compute derivatives of s-vectors

  call deriv_svec(natoms, nmodes, nterms_cart, terms_cart, coefs_cart, nterms_svec, terms_svec, coefs_svec)

  ! compute derivatives of G-matrix

  allocate(coefs_gmat(natoms*3,natoms*3,nterms_gmat))

  call deriv_gmat(natoms, nmodes, masses, nterms_svec, terms_svec, coefs_svec, nterms_gmat, terms_gmat, coefs_gmat)

  do iterm=1, nterms_gmat
    !fac = product(  (/(exp(factln(terms_gmat(imode,iterm))), imode=1, nmodes)/) )
    write(out, '(100(1x,i3))') terms_gmat(:,iterm)
    do imode=1, natoms*3
      write(out, '(30x,100(1x,es16.8))') coefs_gmat(:,imode,iterm)!/fac
    enddo
  enddo

end subroutine expand_gmat


! Computes derivatives of space-fixed Cartesian coordinates of atoms with respect to internal coordinates
! using the Automatic Differentiation technique as implemented in Rapsodia library.
!
! Input:
!   natoms - number of atoms;
!   nmodes - number of internal coordinates (3*natoms-6);
!   masses(1:natoms) - atomic masses;
!   x0(1:nmodes) - equilibrium internal coordinates;
!   order - max derivative order, 0..order derivatives will be calculated;
!   nmodes_active, modes_active(1:nmodes_active) - number and indices of modes wrt which derivatives to be calculated;
!
! Output:
!   nterms, terms(1:nmodes,1:nterms) - number and indices of Cartesian-coordinate derivatives;
!   coefs(1:natoms,1:3,1:nterms) - Cartesian-coordinate derivative values;
subroutine deriv_cart_rapsodia(natoms, nmodes, masses, x0, order, nmodes_active, modes_active, nterms, terms, coefs)

  integer(ik), intent(in)               :: natoms
  integer(ik), intent(in)               :: nmodes
  real(rk), intent(in)                  :: masses(natoms)
  real(rk), intent(in)                  :: x0(nmodes)
  integer(ik), intent(in)               :: order
  integer(ik), intent(in)               :: nmodes_active
  integer(ik), intent(in)               :: modes_active(nmodes_active)
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)
  real(rk), allocatable, intent(out)    :: coefs(:,:,:)

  integer(ik) :: info, iterm, imode, ix, iatom, iorder
  type(higherOrderTensor) :: Tens
  integer(ik), allocatable :: SeedMatrix(:,:)
  real(RAdKind), allocatable :: TaylorCoefficients(:,:,:,:)
  real(RAdKind), allocatable :: CompressedTensor(:,:,:)
  type(RARealD) :: x0_RA(nmodes), cartesian_RA(natoms,3)

  ! estimate number of expansion terms

  if (order==0) then
    nterms = 1
  else
    call setNumberOfIndependents(Tens, nmodes_active)
    call setHighestDerivativeDegree(Tens, order)
    nterms = getDirectionCount(Tens)
  endif

  ! allocate arrays to keep derivative indices and values

  if (allocated(terms)) deallocate(terms)
  if (allocated(coefs)) deallocate(coefs)
  allocate(terms(nmodes,nterms), coefs(natoms,3,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,3(1x,i6))') &!
    'deriv_cart_rapsodia error: failed to allocate terms(nmodes,nterms), coefs(natoms,3,nterms)', &!
    'nmodes, nterms, natoms =', nmodes, nterms, natoms
    stop
  endif

  ! compute derivatives

  if (order==0) then

      x0_RA(1:nmodes) = x0(1:nmodes)
      do imode=1, nmodes_active
        call RAset(x0_RA(modes_active(imode)), 1, 1, real(imode, kind=RAdKind))
      enddo
      call internal_to_cartesian_RA_ptr(natoms, nmodes, masses, x0_RA, cartesian_RA)

      terms(1:nmodes,1) = 0
      coefs(1:natoms,1:3,1) = cartesian_RA%v

  else

    ! allocate array to keep derivative indices for active modes

    allocate(SeedMatrix(nmodes_active,nterms), stat=info)
    if (info/=0) then
      write(out, '(/a/a,2(1x,i6))') &!
      'deriv_cart_rapsodia error: failed to allocate SeedMatrix(nmodes_active,nterms)', &!
      'nmodes_active, nterms =', nmodes_active, nterms
      stop
    endif

    ! generate derivative indices (rays) for modes wrt which derivatives are taken

    call getSeedMatrix(Tens, SeedMatrix)

    ! compute cartesian(1:natoms,1:3) = f(internal)

    x0_RA(1:nmodes) = x0(1:nmodes)
    do iterm=1, nterms
      do imode=1, nmodes_active
        call RAset(x0_RA(modes_active(imode)), iterm, 1, real(SeedMatrix(imode, iterm), kind=RAdKind))
      enddo
    enddo
    call internal_to_cartesian_RA_ptr(natoms, nmodes, masses, x0_RA, cartesian_RA)

    ! retrieve Taylor expansion coefficients

    allocate(TaylorCoefficients(order,nterms,natoms,3), stat=info)
    if (info/=0) then
      write(out, '(/a/a,3(1x,i6))') &!
      'deriv_cart_rapsodia error: failed to allocate TaylorCoefficients(order,nterms,natoms,3)', &!
      'order, nterms, natoms =', order, nterms, natoms
      stop
    endif

    do ix=1, 3
      do iatom=1, natoms
        do iorder=1, order
          do iterm=1, nterms
            call RAget(cartesian_RA(iatom,ix), iterm, iorder, TaylorCoefficients(iorder,iterm,iatom,ix))
          enddo
        enddo
      enddo
    enddo

    ! retrieve compressed derivative tensor

    allocate(CompressedTensor(nterms,natoms,3), stat=info)
    if (info/=0) then
      write(out, '(/a/a,2(1x,i6))') &!
      'deriv_cart_rapsodia error: failed to allocate CompressedTensor(nterms,natoms,3)', &!
      'nterms, natoms =', nterms, natoms
      stop
    endif

    do ix=1, 3
      do iatom=1, natoms
        call setTaylorCoefficients(Tens, TaylorCoefficients(1:order,1:nterms,iatom,ix))
        call getCompressedTensor(Tens, order, CompressedTensor(1:nterms,iatom,ix))
      enddo
    enddo

    ! retrieve expansion coefficients

    do iterm=1, nterms
      terms(:,iterm) = 0
      do imode=1, nmodes_active
        terms(modes_active(imode),iterm) = SeedMatrix(imode,iterm)
      enddo
      coefs(1:natoms,1:3,iterm) = CompressedTensor(iterm,1:natoms,1:3)
    enddo

    deallocate(TaylorCoefficients)
    deallocate(CompressedTensor)
    deallocate(SeedMatrix)

  endif

end subroutine deriv_cart_rapsodia


! Computes derivatives of space-fixed Cartesian coordinates of atoms with respect to internal coordinates
! using the finite-difference technique.
!
! Input:
!   natoms - number of atoms;
!   nmodes - number of internal coordinates (3*natoms-6);
!   masses(1:natoms) - atomic masses;
!   x0(1:nmodes) - equilibrium internal coordinates;
!   cart0_pas(1:natoms,1:3) - equilibrium Cartesian coordinates in principal axes system;
!   order - max derivative order, 0..order derivatives will be calculated;
!   nmodes_active, modes_active(1:nmodes_active) - number and indices of modes wrt which derivatives to be calculated;
!   fdf_npt - number of points in finite-difference formula (3, 5, 7, or 9);
!   fdf_step(1:nmodes) - finite-difference step size for each mode;
!   frame_rotation - keyword that defines Cartesian coordinate frame rotation function ("eckart", "pas", "I" ...);
!   rot_maxorder - max truncation order in Taylor series expansion of matrix exponential in frame rotation function;
!   rot_maxit - max number of iterations in frame rotation function;
!   rot_accur - accuracy of solution in frame rotation function;
!
! Output:
!   nterms, terms(1:nmodes,1:nterms) - number and indices of Cartesian-coordinate derivatives;
!   coefs(1:natoms,1:3,1:nterms) - Cartesian-coordinate derivative values;
subroutine deriv_cart_fdf(natoms, nmodes, masses, x0, cart0_pas, order, nmodes_active, modes_active, fdf_npt, &
                          fdf_step, frame_rotation, rot_maxorder, rot_maxit, rot_accur, nterms, terms, coefs)

  integer(ik), intent(in)               :: natoms
  integer(ik), intent(in)               :: nmodes
  real(rk), intent(in)                  :: masses(natoms)
  real(rk), intent(in)                  :: x0(nmodes)
  real(rk), intent(in)                  :: cart0_pas(natoms,3)
  integer(ik), intent(in)               :: order
  integer(ik), intent(in)               :: nmodes_active
  integer(ik), intent(in)               :: modes_active(nmodes_active)
  integer(ik), intent(in)               :: fdf_npt
  real(rk), intent(in)                  :: fdf_step(nmodes)
  character(*), intent(in)              :: frame_rotation
  integer(ik), intent(in)               :: rot_maxorder
  integer(ik), intent(in)               :: rot_maxit
  real(rk), intent(in)                  :: rot_accur
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)
  real(rk), allocatable, intent(out)    :: coefs(:,:,:)

  integer(ik) :: imode, jmode, ndeg(nmodes), deg(order+1,nmodes), iorder, iterm, info
  real(rk) :: fdf_coefs(4,9), fdf_h(9)

  ! init 3, 5, 7, and 9-point finite-difference coefficients

  if (count((/3,5,7,9/)==fdf_npt)==0) then
    write(out, '(/a,1x,i3,1x,a)') &!
    'deriv_cart_fdf error: wrong number of difference points fdf_npt =', fdf_npt, ' (can be 3, 5, 7, or 9)'
    stop
  endif
  fdf_coefs = 0.0
  fdf_coefs(1:1,3) = (/1.0_rk/)
  fdf_coefs(1:2,5) = (/8.0_rk, -1.0_rk/)
  fdf_coefs(1:3,7) = (/45.0_rk, -9.0_rk, 1.0_rk/)
  fdf_coefs(1:4,9) = (/672.0_rk, -168.0_rk, 32.0_rk, -3.0_rk/)
  fdf_h(3) = 2.0_rk
  fdf_h(5) = 12.0_rk
  fdf_h(7) = 60.0_rk
  fdf_h(9) = 840.0_rk

  ! generate derivative terms for expansion 0..order

  ndeg = 1
  deg = 0
  do imode=1, nmodes_active
    jmode = modes_active(imode)
    ndeg(jmode) = order+1
    deg(1:ndeg(jmode), jmode) = (/(iorder, iorder=0, order)/)
  enddo
  call expansion_terms(nmodes, ndeg(1:nmodes), deg(1:maxval(ndeg),1:nmodes), order, order, nterms, terms)

  ! allocate array to keep derivatives

  if (allocated(coefs)) deallocate(coefs)
  allocate(coefs(natoms,3,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &!
    'deriv_cart_fdf error: failed to allocate coefs(natoms,3,nterms)', 'natoms, nterms =', natoms, nterms
    stop
  endif

  ! compute finite-difference derivatives

  do iterm=1, nterms
    coefs(1:natoms,1:3,iterm) = rec_fdf(x0, terms(1:nmodes,iterm))
  enddo

  contains

  recursive function rec_fdf(x, term) result(f)

    real(rk), intent(in)    :: x(nmodes)
    integer(ik), intent(in) :: term(nmodes)
    real(rk)                :: f(natoms,3)

    integer(ik) :: imode, t(nmodes), ipt
    real(rk) :: xp(nmodes), xm(nmodes), fp(natoms,3), fm(natoms,3), cart(natoms,3)

    do imode=1, nmodes
      if (term(imode)>0) then
        t = term
        t(imode) = t(imode) - 1
        f = 0.0
        do ipt=1, (fdf_npt-1)/2+1
          xp = x
          xm = x
          xp(imode) = xp(imode) + ipt*fdf_step(imode)
          xm(imode) = xm(imode) - ipt*fdf_step(imode)
          fp = rec_fdf(xp, t)
          fm = rec_fdf(xm, t)
          f = f + fdf_coefs(ipt,fdf_npt) * (fp - fm)
        enddo
        f = f / (fdf_h(fdf_npt)*fdf_step(imode))
        exit
      endif
    enddo
    if (all(term==0)) then
      call internal_to_cartesian_ptr(natoms, nmodes, masses, x, cart)
      select case(trim(frame_rotation))
      case('ECKART')
        call rotate_eckart(natoms, masses, cart, cart0_pas, rot_maxorder, rot_maxit, rot_accur, f)
      case('I')
        call rotate_I(natoms, masses, cart, f)
      case default
        write(out, '(/a,1x,a)') &!
        'deriv_cart_fdf error: unknown coordinate frame rotation function =', trim(frame_rotation)
        stop
      end select
    endif

  end function rec_fdf

end subroutine deriv_cart_fdf


! Computes derivatives of space-fixed Cartesian coordinates of atoms with respect to internal coordinates
! using the finite-difference technique (quadrupole precision).
!
! Input:
!   natoms - number of atoms;
!   nmodes - number of internal coordinates (3*natoms-6);
!   masses(1:natoms) - atomic masses;
!   x0(1:nmodes) - equilibrium internal coordinates;
!   cart0_pas(1:natoms,1:3) - equilibrium Cartesian coordinates in principal axes system;
!   order - max derivative order, 0..order derivatives will be calculated;
!   nmodes_active, modes_active(1:nmodes_active) - number and indices of modes wrt which derivatives to be calculated;
!   fdf_npt - number of points in finite-difference formula (3, 5, 7, or 9);
!   fdf_step(1:nmodes) - finite-difference step size for each mode;
!   frame_rotation - keyword that defines Cartesian coordinate frame rotation function ("eckart", "pas", "I" ...);
!   rot_maxorder - max truncation order in Taylor series expansion of matrix exponential in frame rotation function;
!   rot_maxit - max number of iterations in frame rotation function;
!   rot_accur - accuracy of solution in frame rotation function;
!
! Output:
!   nterms, terms(1:nmodes,1:nterms) - number and indices of Cartesian-coordinate derivatives;
!   coefs(1:natoms,1:3,1:nterms) - Cartesian-coordinate derivative values;
subroutine deriv_cart_fdf_ark(natoms, nmodes, masses, x0, cart0_pas, order, nmodes_active, modes_active, fdf_npt, &!
                              fdf_step, frame_rotation, rot_maxorder, rot_maxit, rot_accur, nterms, terms, coefs)

  integer(ik), intent(in)               :: natoms
  integer(ik), intent(in)               :: nmodes
  real(rk), intent(in)                  :: masses(natoms)
  real(rk), intent(in)                  :: x0(nmodes)
  real(rk), intent(in)                  :: cart0_pas(natoms,3)
  integer(ik), intent(in)               :: order
  integer(ik), intent(in)               :: nmodes_active
  integer(ik), intent(in)               :: modes_active(nmodes_active)
  integer(ik), intent(in)               :: fdf_npt
  real(rk), intent(in)                  :: fdf_step(nmodes)
  character(*), intent(in)              :: frame_rotation
  integer(ik), intent(in)               :: rot_maxorder
  integer(ik), intent(in)               :: rot_maxit
  real(rk), intent(in)                  :: rot_accur
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)
  real(rk), allocatable, intent(out)    :: coefs(:,:,:)

  integer(ik) :: imode, jmode, ndeg(nmodes), deg(order+1,nmodes), iorder, iterm, info, iatom
  real(ark) :: masses_ark(natoms), cart0_pas_ark(natoms,3), x0_ark(nmodes), fdf_step_ark(nmodes), rot_accur_ark, &!
               fdf_coefs(4,9), fdf_h(9), coefs_ark(natoms,3)
  character(1000) :: s

  ! convert masses, equilibrium coordinates, etc. into quadrupole precision

  do iatom=1, natoms
    s = ''
    write(s,*) masses(iatom)
    read(s,*) masses_ark(iatom)
    s = ''
    write(s,*) cart0_pas(iatom,1:3)
    read(s,*) cart0_pas_ark(iatom,1:3)
  enddo
  do imode=1, nmodes
    s = ''
    write(s,*) x0(imode)
    read(s,*) x0_ark(imode)
    s = ''
    write(s,*) fdf_step(imode)
    read(s,*) fdf_step_ark(imode)
  enddo
  s = ''
  write(s,*) rot_accur
  read(s,*) rot_accur_ark

  ! init 3, 5, 7, and 9-point finite-difference coefficients

  if (count((/3,5,7,9/)==fdf_npt)==0) then
    write(out, '(/a,1x,i3,1x,a)') &!
    'deriv_cart_fdf_ark error: wrong number of difference points fdf_npt =', fdf_npt, ' (can be 3, 5, 7, or 9)'
    stop
  endif
  fdf_coefs = 0.0
  fdf_coefs(1:1,3) = (/1.0_ark/)
  fdf_coefs(1:2,5) = (/8.0_ark, -1.0_ark/)
  fdf_coefs(1:3,7) = (/45.0_ark, -9.0_ark, 1.0_ark/)
  fdf_coefs(1:4,9) = (/672.0_ark, -168.0_ark, 32.0_ark, -3.0_ark/)
  fdf_h(3) = 2.0_ark
  fdf_h(5) = 12.0_ark
  fdf_h(7) = 60.0_ark
  fdf_h(9) = 840.0_ark

  ! generate derivative terms for expansion 0..order

  ndeg = 1
  deg = 0
  do imode=1, nmodes_active
    jmode = modes_active(imode)
    ndeg(jmode) = order+1
    deg(1:ndeg(jmode), jmode) = (/(iorder, iorder=0, order)/)
  enddo
  call expansion_terms(nmodes, ndeg(1:nmodes), deg(1:maxval(ndeg),1:nmodes), order, order, nterms, terms)

  ! allocate array to keep derivatives

  if (allocated(coefs)) deallocate(coefs)
  allocate(coefs(natoms,3,nterms), stat=info)
  if (info/=0) then
    write(out, '(/a/a,2(1x,i6))') &!
    'deriv_cart_fdf_ark error: failed to allocate coefs(natoms,3,nterms)', 'natoms, nterms =', natoms, nterms
    stop
  endif

  ! compute finite-difference derivatives

  do iterm=1, nterms
    coefs_ark(1:natoms,1:3) = rec_fdf(x0_ark, terms(1:nmodes,iterm))
    coefs(1:natoms,1:3,iterm) = real(coefs_ark(1:natoms,1:3), kind=rk)
  enddo

  contains

  recursive function rec_fdf(x, term) result(f)

    real(ark), intent(in)   :: x(nmodes)
    integer(ik), intent(in) :: term(nmodes)
    real(ark)               :: f(natoms,3)

    integer(ik) :: imode, t(nmodes), ipt
    real(ark) :: xp(nmodes), xm(nmodes), fp(natoms,3), fm(natoms,3), cart(natoms,3)

    do imode=1, nmodes
      if (term(imode)>0) then
        t = term
        t(imode) = t(imode) - 1
        f = 0.0
        do ipt=1, (fdf_npt-1)/2+1
          xp = x
          xm = x
          xp(imode) = xp(imode) + ipt*fdf_step_ark(imode)
          xm(imode) = xm(imode) - ipt*fdf_step_ark(imode)
          fp = rec_fdf(xp, t)
          fm = rec_fdf(xm, t)
          f = f + fdf_coefs(ipt,fdf_npt) * (fp - fm)
        enddo
        f = f / (fdf_h(fdf_npt)*fdf_step(imode))
        exit
      endif
    enddo
    if (all(term==0)) then
      call internal_to_cartesian_ark_ptr(natoms, nmodes, masses_ark, x, cart)
      select case(trim(frame_rotation))
      case('ECKART')
        call rotate_eckart_ark(natoms, masses_ark, cart, cart0_pas_ark, rot_maxorder, rot_maxit, rot_accur_ark, f)
      case('I')
        call rotate_I_ark(natoms, masses_ark, cart, f)
      case default
        write(out, '(/a,1x,a)') &!
        'deriv_cart_fdf_ark error: unknown coordinate frame rotation function =', trim(frame_rotation)
        stop
      end select
    endif

  end function rec_fdf

end subroutine deriv_cart_fdf_ark


! Computes derivatives of Cartesian coordinates of atoms in rotated coordinate frame wrt internal coordinates.
!
! Coordinate rotation:
!     cartesian_{rotated} = rotation_matrix * cartesian,
!     where cartesian are Cartesian coordinates of atoms.
!
! Input:
!   natoms - number of atoms;
!   nmodes - number of internal coordinates (3*natoms-6);
!   nterms, terms(1:nmodes,1:nterms) - number and indices of derivatives to be computed;
!   coefs_cart(1:natoms,1:3,1:nterms) - derivatives of space-fixed Cartesian coordinates of atoms;
!   coefs_rotmat(1:3,1:3,1:nterms) - derivatives of coordinate frame rotation matrix;
!
! Output:
!   coefs(1:natoms,1:3,1:nterms) - derivatives of Cartesian coordinates of atoms in rotated coordinate frame;
subroutine deriv_cart_rotated(natoms, nmodes, nterms, terms, coefs_cart, coefs_rotmat, coefs)

  integer(ik), intent(in) :: natoms
  integer(ik), intent(in) :: nmodes
  integer(ik), intent(in) :: nterms
  integer(ik), intent(in) :: terms(nmodes,nterms)
  real(rk), intent(in)    :: coefs_cart(natoms,3,nterms)
  real(rk), intent(in)    :: coefs_rotmat(3,3,nterms)
  real(rk), intent(out)   :: coefs(natoms,3,nterms)

  integer(ik) :: nterms_chr(nterms), iterm, ielem, iterm1, iterm2, info
  integer(ik), allocatable :: ind_chr(:,:,:)

  ! generate chain rule expressions for all derivatives of product rotmat*cartesian

  call chain_rule_product_ind(nmodes, nterms, terms, 2, nterms_chr, ind_chr)

  do iterm=1, nterms

    ! sum up all elements in the chain rule expression for terms(:,iterm)-derivative of rotmat*cartesian product

    coefs(:,:,iterm) = 0.0
    do ielem=1, nterms_chr(iterm)
      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)
      coefs(:,:,iterm) = coefs(:,:,iterm) + matmul(coefs_cart(:,:,iterm2), transpose(coefs_rotmat(:,:,iterm1)))
    enddo

  enddo

  if (allocated(ind_chr)) deallocate(ind_chr)

end subroutine deriv_cart_rotated


! Computes derivatives of Eckart coordinate rotation matrix.
!
! Coordinate rotation:
!     cartesian_{eckart} = exp(-kappa) * cartesian,
!     where kappa = -kappa^T is skew-symmetric rotation matrix,
!     and cartesian are Cartesian coordinates of atoms.
!
! Eckart equations:
!     u*exp(kappa) - exp(-kappa)*u^T = 0,
!     where u_{ix,iy} = \sum_{iatom} mass_{iatom} * cartesian_eq_pas_{iatom,ix} * cartesian_{iatom,iy},
!     and cartesian_eq_pas are equilibrium Cartesian coordinates in principal axes system.
!
! Eckart equations for p-th derivative of exp(-kappa):
!     u*exp(kappa)^{(p)} - exp(-kappa)^{(p)}*u^T + vec = 0,
!     where exp(kappa)^{(p)} is a p-th derivative of exp(kappa),
!     and vec includes the rest terms of the expression for p-th derivative of the Eckart equations.
!
! Iterative solution for p-th derivative of exp(-kappa):
!     u*kappa^{(p)} + kappa^{(p)}*u^T = lambda^{(p)}*u^T - u*(lambda^{(p)})^T - vec
!     where lambda^{(p)} = exp(-kappa)^{(p)} + kappa^{(p)} is calculated from kappa^{(p)} on previous iteration.
!
! Input:
!     natoms - number of atoms;
!     nmodes - number of internal coordinates (3*natoms-6);
!     masses(1:natoms) - atomic masses;
!     cart_eq_pas(1:natoms,1:3) - equilibrium Cartesian coordinates of atoms in principal axes system;
!     nterms, terms(1:nmodes,1:nterms) - number and indices of frame rotation-matrix derivatives to be calculated;
!     coefs_cart(1:natoms,1:3,1:nterms) - derivatives of space-fixed Cartesian coordinates of atoms;
!     maxdeg - max truncation order in Taylor series expansion of matrix exponential (exp(-kappa));
!     maxit - max number of iterations for Eckart solution;
!     accur - desired accuracy of Eckart solution;
!
! Output:
!     coefs_expkappa(1:3,1:3,1:nterms) - derivatives of Eckart coordinate rotation matrix (exp(-kappa));
subroutine deriv_rotmat_eckart(natoms, nmodes, masses, cart_eq_pas, nterms, terms, coefs_cart, maxdeg, maxit, accur, &!
                               coefs_expkappa)

  integer(ik), intent(in) :: natoms
  integer(ik), intent(in) :: nmodes
  real(rk), intent(in)    :: masses(natoms)
  real(rk), intent(in)    :: cart_eq_pas(natoms,3)
  integer(ik), intent(in) :: nterms
  integer(ik), intent(in) :: terms(nmodes,nterms)
  real(rk), intent(in)    :: coefs_cart(natoms,3,nterms)
  integer(ik), intent(in) :: maxdeg
  integer(ik), intent(in) :: maxit
  real(rk), intent(in)    :: accur
  real(rk), intent(out)   :: coefs_expkappa(3,3,nterms)

  integer(ik) :: nterms_chr(nterms), iterm, i, j, ielem, iterm1, iterm2, ideg, iter
  integer(ik), allocatable :: ind_chr(:,:,:)
  real(rk) :: coefs_umat(3,3,nterms), coefs_kappa(3,3,maxdeg,nterms), tmat(3,3), umat(3,3), vec0(3), &!
              mat(3,3), kappa0(3,3), vecp(3), lambda(3,3), vec(3)
  logical :: coefs_kappa_init(maxdeg,nterms), coefs_expkappa_init(nterms), coefs_umat_init(nterms)
  integer(ik) :: lwork, info, rank
  double precision :: matd(3,3), vecd(3,1), sv(3), work(128*3), rcond

  ! generate chain rule expressions for all derivatives of product u*exp(kappa)

  call chain_rule_product_ind(nmodes, nterms, terms, 2, nterms_chr, ind_chr)

  ! precompute derivatives of u matrix

  coefs_umat_init(:) = .false.
  do iterm=1, nterms
    do i=1, 3
      do j=1, 3
        coefs_umat(i,j,iterm) = sum( masses(1:natoms) * cart_eq_pas(1:natoms,i) * coefs_cart(1:natoms,j,iterm) )
        coefs_umat_init(iterm) = .true.
      enddo
    enddo
  enddo

  ! start recursive calculations of derivatives of exp(-kappa)

  coefs_kappa(:,:,:,:) = 0.0
  coefs_kappa_init(:,:) = .false.
  coefs_expkappa(:,:,:) = 0.0
  coefs_expkappa_init(:) = .false.

  do iterm=1, nterms

    ! compute system of non-linear equations for terms(:,iterm) derivative of exp(-kappa)
    ! u*kappa^{(iterm)} + kappa^{(iterm)}*u^T = lambda^{(iterm)}*u^T - u*(lambda^{(iterm)})^T - vec

    tmat = 0.0
    umat = 0.0
    do ielem=1, nterms_chr(iterm)
      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)
      if (.not.coefs_umat_init(iterm1)) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
        'deriv_rotmat_eckart error: derivative = (', terms(:,iterm1), ') of u-matrix is not initialized'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
        'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
        stop
      endif
      if (iterm2==iterm) then
        ! compute u
        umat = coefs_umat(:,:,iterm1)
      else
        if (.not.coefs_expkappa_init(iterm2)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
          'deriv_rotmat_eckart error: derivative = (', terms(:,iterm2), ') of exp(-kappa) is not initialized'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
          'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
          write(out, '(a)') &!
          'possible cause: derivative indices are not sorted ascendingly wrt total derivative order'
          stop
        endif
        ! compute vec
        tmat = tmat + matmul(coefs_umat(:,:,iterm1), transpose(coefs_expkappa(:,:,iterm2)))
      endif
    enddo

    ! compute vec - part of the right hand side vector that does not depend on terms(:,iterm) derivative of exp(-kappa)

    tmat = tmat - transpose(tmat)
    vec0 = (/tmat(1,2), tmat(1,3), tmat(2,3)/)

    ! compute matrix u*kappa^{(iterm)} + kappa^{(iterm)}*u^T

    mat = 0.0
    ielem = 0
    do i=1, 3
      do j=i+1, 3
        ielem = ielem + 1
        kappa0 = 0.0
        kappa0(i,j) = 1.0_rk
        kappa0(j,i) = -1.0_rk
        tmat = matmul(umat, kappa0) + matmul(kappa0, transpose(umat))
        mat(1,ielem) = tmat(1,2)
        mat(2,ielem) = tmat(1,3)
        mat(3,ielem) = tmat(2,3)
      enddo
    enddo

    ! solve system of non-linear equations for terms(:,iterm) derivative of exp(-kappa)

    coefs_kappa(:,:,:,iterm) = 0.0
    coefs_kappa_init(:,iterm) = .false.
    coefs_kappa_init(1,iterm) = .true.
    coefs_expkappa_init(iterm) = .false.
    vecp = 0.0

    do iter=1, maxit

      coefs_kappa_init(2:,iterm) = .false.
      coefs_kappa(:,:,2:,iterm) = 0.0

      ! compute term(:,iterm) derivative of all powers (2..taylor_maxorder) of kappa

      do ideg=2, maxdeg
        tmat = 0.0
        do ielem=1, nterms_chr(iterm)
          iterm1 = ind_chr(1,ielem,iterm)
          iterm2 = ind_chr(2,ielem,iterm)
          if (.not.coefs_kappa_init(1,iterm1)) then
            write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
            'deriv_rotmat_eckart error: derivative = (', terms(:,iterm1), ') of kappa**1 is not initialized'
            write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
            'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
            write(out, '(a)') &!
            'possible cause: derivative indices are not sorted ascendingly wrt total derivative order'
            stop
          endif
          if (.not.coefs_kappa_init(ideg-1,iterm2)) then
            write(out, '(/a,1x,<nmodes>(1x,i3),1x,a,i2,1x,a)') &!
            'deriv_rotmat_eckart error: derivative = (', terms(:,iterm2), ') of kappa**', ideg-1, 'is not initialized'
            write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
            'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
            write(out, '(a)') &!
            'possible cause: derivative indices are not sorted ascendingly wrt total derivative order'
            stop
          endif
          tmat = tmat + matmul(coefs_kappa(:,:,1,iterm1), coefs_kappa(:,:,ideg-1,iterm2))
        enddo
        coefs_kappa(:,:,ideg,iterm) = tmat/real(ideg, kind=rk)
        coefs_kappa_init(ideg,iterm) = .true.
      enddo

      ! compute term(:,iterm) derivative of exp(-kappa)

      tmat = 0.0
      if (all(terms(:,iterm)==0)) forall(i=1:3) tmat(i,i) = 1.0_rk
      tmat = tmat - coefs_kappa(:,:,1,iterm)
      do ideg=2, maxdeg
        if (mod(ideg,2)==0) then
          tmat = tmat + coefs_kappa(:,:,ideg,iterm)
        else
          tmat = tmat - coefs_kappa(:,:,ideg,iterm)
        endif
        if (all(abs(coefs_kappa(:,:,ideg,iterm))<accur)) exit
        if (ideg==maxdeg) then
          write(out, '(/a,1x,es16.8)') &!
          'deriv_rotmat_eckart error: failed to converge derivatives of exp(-kappa) to accuracy =', accur
          stop
        endif
      enddo
      coefs_expkappa(:,:,iterm) = tmat

      ! compute term(:,iterm) derivative of lambda = exp(-kappa) + kappa

      lambda = coefs_expkappa(:,:,iterm) + coefs_kappa(:,:,1,iterm)

      ! compute right hand side vector

      tmat = matmul(lambda, transpose(umat)) - matmul(umat, transpose(lambda))
      vec = (/ tmat(1,2)-vec0(1), tmat(1,3)-vec0(2), tmat(2,3)-vec0(3) /)

      ! solve system of linear equations

      matd = dble(mat)
      vecd(1:3,1) = dble(vec)
      lwork = size(work)
      rcond = -1.0d-14
      call dgelss(3, 3, 1, matd, 3, vecd, 3, sv, rcond, rank, work, lwork, info)

      if (info/=0) then
        write(out, '(/a,1x,i6)') 'deriv_rotmat_eckart error: failed to solve system of linear equations, lapack info =', info
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
        stop
      endif
      vec = real(vecd(:,1), kind=rk)

      ! check convergence

      if (all(abs(vec-vecp)<accur)) exit
      vecp = vec

      if (iter==maxit) then
        write(out, '(/a,1x,es16.8,1x,a,1x,i3,1x,a)') &!
        'deriv_rotmat_eckart error: failed to converge Eckart equations to accuracy =', accur, 'in', maxit, 'iterations'
         write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of exp(-kappa) = (', terms(:,iterm), ')'
        stop
      endif

      ! update kappa

      coefs_kappa(1,1:3,1,iterm) = (/ 0.0_rk,  vec(1), vec(2)/)
      coefs_kappa(2,1:3,1,iterm) = (/-vec(1),  0.0_rk, vec(3)/)
      coefs_kappa(3,1:3,1,iterm) = (/-vec(2), -vec(3), 0.0_rk/)

    enddo !iter

    coefs_expkappa_init(iterm) = .true.

  enddo !iterm

  if (allocated(ind_chr)) deallocate(ind_chr)

end subroutine deriv_rotmat_eckart


! Computes derivatives of s-vectors with respect to internal coordinates.
!
! Definition of s-vectors:
!     s_{imode,iatom,ix} = d internal_{imode} / d cartesian_{iatom,ix},
!     where imode=1..3*natoms (3N-6 vibrations, 3 Euler angles, and 3 translations), ix=x,y,z, and iatom=1..natoms.
!
! Linear equations for s vectors:
!     \sum_{iatom}\sum_{ix=x,y,z} t_{iatom,ix,imode} * s_{jmode,iatom,ix} = \delta_{imode,jmode},
!     where t_{iatom,ix,imode} = d cartesian_{iatom,ix} / d internal_{imode}.
!
! Linear equations for p-th derivative of s-vectors:
!     \sum_{iatom}\sum_{ix=x,y,z} t_{iatom,ix,imode} * s_{jmode,iatom,ix}^{(p)} + vec = 0,
!     where s^{(p)} is a p-th derivative of s-vector,
!     and vec includes the rest terms of the expression for p-th derivative of t*s product.
!
! Input:
!     natoms - number of atoms;
!     nmodes - number of internal coordinates (3*natoms-6);
!     nterms_cart, terms_cart(1:nmodes,1:nterms_cart), coefs_cart(1:natoms,1:3,1:nterms_cart) - number, indices
!       and values of Cartesian-coordinate derivatives;
!     nterms, terms(1:nmodes,1:nterms) - number and indices of desired s-vector derivatives;
!
! Output:
!     coefs_svec(ncoords,3*natoms,1:nterms) - s-vector derivative values, the first index runs through
!       3*natoms-6 vibrations, 3 Euler rotations, and 3 translations, and the second index runs through
!       3*natoms Cartesian coordinates of atoms;
subroutine deriv_svec(natoms, nmodes, nterms_cart, terms_cart, coefs_cart, nterms, terms, coefs_svec)

  integer(ik), intent(in) :: natoms
  integer(ik), intent(in) :: nmodes
  integer(ik), intent(in) :: nterms_cart
  integer(ik), intent(in) :: terms_cart(nmodes,nterms_cart)
  real(rk), intent(in)    :: coefs_cart(natoms,3,nterms_cart)
  integer(ik), intent(in) :: nterms
  integer(ik), intent(in) :: terms(nmodes,nterms)
  real(rk), intent(out)   :: coefs_svec(:,:,:)

  integer(ik) :: ncoords, nterms_chr(nterms), iterm, iterm1, iterm2, ielem, i, j, k, t(nmodes), ipos, imode, ix, iatom
  integer(ik), allocatable :: ind_chr(:,:,:)
  real(rk) :: tmat(natoms*3,natoms*3), bvec(natoms*3,natoms*3), coefs_tvec(natoms*3,natoms*3), asym_tens(3,3,3)
  logical :: coefs_svec_init(nterms)
  integer(ik) :: lwork, info, rank
  double precision :: matd(natoms*3,natoms*3), vecd(natoms*3,natoms*3), sv(natoms*3), work(128*natoms*3), rcond

  ncoords = natoms*3

  ! precompute 3D levi-civita tensor

  asym_tens = 0.0
  asym_tens(1,2,3) = 1.0_rk
  asym_tens(1,3,2) =-1.0_rk
  asym_tens(2,1,3) =-1.0_rk
  asym_tens(2,3,1) = 1.0_rk
  asym_tens(3,1,2) = 1.0_rk
  asym_tens(3,2,1) =-1.0_rk

  ! generate chain rule expressions for all derivatives of product t*s

  call chain_rule_product_ind(nmodes, nterms, terms, 2, nterms_chr, ind_chr)

  ! start recursive calculations of s-vectors derivatives

  coefs_svec(:,:,:) = 0.0
  coefs_svec_init(:) = .false.

  do iterm=1, nterms

    ! compute matrix and right hand side vectors for system of linear equations:
    ! d(t_{imode}*s_{jmode}) = d(\delta_{imode,jmode})

    tmat = 0.0
    bvec = 0.0
    if (all(terms(1:nmodes,iterm)==0)) forall(i=1:ncoords) bvec(i,i) = 1.0_rk

    do ielem=1, nterms_chr(iterm)
      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)

      ! compute terms(:,iterm1)-derivative of vibrational, rotational and translational t-vectors;
      ! details in S.N.Yurchenko, W.Thiel, and P.Jensen, J. Mol. Spectrosc. 245 (2007), 126-140,
      ! and G.O.Sorensen, M.J.S. Dewar (Ed.) et al., Topics in Current Chemistry, vol. 82, 1979, p. 99

      coefs_tvec(:,:) = 0.0

      ! derivative of vibrational t-vector:
      ! t_{iatom,ix,imode}^{vib} = d cartesian_{iatom,ix} / d internal_{imode}

      do imode=1, nmodes
        t = terms(1:nmodes,iterm1)
        t(imode) = t(imode) + 1
        ipos = index_iarr1(t, terms_cart)
        if (ipos==0) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
          'deriv_svec error: derivative = (', t, ') of Cartesian coordinates is not found'
          write(out, '(a,1x,i3,1x,a,1x,<nmodes>(1x,i3),1x,a)') &!
          'target derivative of vibrational t-vector for', imode, 'mode = (', terms(:,iterm1), ')'
          stop
        else
          coefs_tvec(:,imode) = reshape(coefs_cart(:,:,ipos), (/ncoords/))
        endif
      enddo

      ! derivative of rotational t-vector:
      ! t_{iatom,ix,imode}^{rot} = \sum_{iy} asym_tensor_{ix,imode,iy} * cartesian_{iatom,iy}

      t = terms(1:nmodes,iterm1)
      ipos = index_iarr1(t, terms_cart)
      if (ipos==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
        'deriv_svec error: derivative = (', t, ') of Cartesian coordinates is not found'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
        'target derivative of rotational t-vector = (', terms(:,iterm1), ')'
        stop
      else
        coefs_tvec(:,nmodes+1) = reshape( matmul(coefs_cart(:,:,ipos), transpose(asym_tens(:,1,:))), (/ncoords/) )
        coefs_tvec(:,nmodes+2) = reshape( matmul(coefs_cart(:,:,ipos), transpose(asym_tens(:,2,:))), (/ncoords/) )
        coefs_tvec(:,nmodes+3) = reshape( matmul(coefs_cart(:,:,ipos), transpose(asym_tens(:,3,:))), (/ncoords/) )
      endif

      ! derivative of translational t-vector:
      ! t_{iatom,ix,imode}^{tran} = \delta_{ix,imode}

      if (all(terms(1:nmodes,iterm1)==0)) then
        i = 0
        do ix=1, 3
          do iatom=1, natoms
            i = i + 1
            coefs_tvec(i,nmodes+3+ix) = 1.0_rk
          enddo
        enddo
      endif

      ! if current element of the chain rule expression (ielem) contains target derivative of s-vectors (iterm),
      ! then the corresponding derivative of t-vectors (coefs_tvec) forms the linear equations matrix,
      ! else it is multiplied by already calculated (on previous steps) derivative of s-vectors and result is added
      ! to the right hand side vectors of system of linear equations.

      if (iterm2==iterm) then
        tmat = transpose(coefs_tvec)
      else
        if (.not.coefs_svec_init(iterm2)) then
          write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
          'deriv_svec error: derivative = (', terms(:,iterm2), ') of s-vectors is not initialized'
          write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
          'target derivative of s-vectors = (', terms(:,iterm), ')'
          write(out, '(a)') &!
          'possible cause: derivative indices are not sorted ascendingly wrt total derivative order'
          stop
        else
          bvec = bvec - transpose( matmul(coefs_svec(:,:,iterm2), coefs_tvec) )
        endif
      endif

    enddo !ielem

    ! solve 3*natoms systems of linear equations for terms(:,iterm) derivative of s-vectors

    matd = dble(tmat)
    vecd = dble(bvec)
    lwork = size(work)
    rcond = -1.0d-14

    call dgelss(ncoords, ncoords, ncoords, matd, ncoords, vecd, ncoords, sv, rcond, rank, work, lwork, info)

    if (info/=0) then
      write(out, '(/a,1x,i6)') &!
      'deriv_svec error: failed to solve system of linear equations, lapack info =', info
      write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') 'target derivative of s-vectors = (', terms(:,iterm), ')'
      stop
    endif

    coefs_svec(:,:,iterm) = transpose( real(vecd, kind=rk) )
    coefs_svec_init(iterm) = .true.

  enddo !iterm

  if (allocated(ind_chr)) deallocate(ind_chr)

end subroutine deriv_svec


! Computes derivatives of kinetic energy matrix G.
!
! Definition of kinetic energy matrix:
!   G_{imode,jmode} = \sum_{iatom} 1/mass_{iatom} \sum_{ix=x,y,z} s_{imode,iatom,ix} * s_{jmode,iatom,ix},
!   where s_{imode,iatom,ix} = d internal_{imode} / d cartesian_{iatom,ix},
!   and imode,jmode=1..3*natoms (3N-6 vibrations, 3 Euler angles, and 3 translations), ix=x,y,z, and iatom=1..natoms.
!
! Input:
!   natoms - number of atoms;
!   nmodes - number of internal coordinates (modes);
!   masses(1:natoms) - atomic masses;
!   nterms_svec - number of s-vector derivatives;
!   terms_svec(1:nmodes,1:nterms_svec) - s-vector derivative indices;
!   coefs_svec(1:natoms*3,1:natoms*3,1:nterms_svec) - s-vector derivative values, where the first index runs through
!     3*natoms-6 vibrations, 3 Euler rotations, and 3 translations, and the second index runs through
!     3*natoms Cartesian coordinates;
!   nterms, terms(1:nmodes,1:nterms) - number and indices of G-matrix derivatives to be calculated;
!
! Output:
!   coefs(1:natoms*3,1:natoms*3) - G-matrix derivative values;
subroutine deriv_gmat(natoms, nmodes, masses, nterms_svec, terms_svec, coefs_svec, nterms, terms, coefs)

  integer(ik), intent(in) :: natoms
  integer(ik), intent(in) :: nmodes
  real(rk), intent(in)    :: masses(natoms)
  integer(ik), intent(in) :: nterms_svec
  integer(ik), intent(in) :: terms_svec(nmodes,nterms_svec)
  real(rk), intent(in)    :: coefs_svec(natoms*3,natoms*3,nterms_svec)
  integer(ik), intent(in) :: nterms
  integer(ik), intent(in) :: terms(nmodes,nterms)
  real(rk), intent(out)   :: coefs(natoms*3,natoms*3,nterms)

  integer(ik) :: nterms_chr(nterms), iterm, ielem, iatom, iterm1, iterm2, ipos1, ipos2, icoord, jcoord, ncoords
  integer(ik), allocatable :: ind_chr(:,:,:)
  real(rk) :: inv_masses(natoms), deriv(natoms*3,natoms*3,natoms), deriv1(natoms*3,natoms,3), deriv2(natoms*3,natoms,3)

  ncoords = natoms*3

  ! generate chain rule expressions for all derivatives of product s*s

  call chain_rule_product_ind(nmodes, nterms, terms, 2, nterms_chr, ind_chr)

  ! precompute inverse masses

  inv_masses = 0.0
  do iatom=1, natoms
    inv_masses(iatom) = 1.0_rk / masses(iatom)
  enddo

  do iterm=1, nterms

    deriv = 0.0
    do ielem=1, nterms_chr(iterm)
      iterm1 = ind_chr(1,ielem,iterm)
      iterm2 = ind_chr(2,ielem,iterm)
      ipos1 = index_iarr1(terms(:,iterm1), terms_svec)
      ipos2 = index_iarr1(terms(:,iterm2), terms_svec)
      if (ipos1==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
        'deriv_gmat error: derivative = (', terms(:,iterm1), ') of s-vectors is not found'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
        'target derivative of G-matrix = (', terms(:,iterm), ')'
        stop
      endif
      if (ipos2==0) then
        write(out, '(/a,1x,<nmodes>(1x,i3),1x,a)') &!
        'deriv_gmat error: derivative = (', terms(:,iterm2), ') of s-vectors is not found'
        write(out, '(a,1x,<nmodes>(1x,i3),1x,a)') &!
        'target derivative of G-matrix = (', terms(:,iterm), ')'
        stop
      endif
      ! compute intermediate: deriv_{imode,jmode,iatom} = \sum_{ix=x,y,z} s_{imode,iatom,ix} * s_{jmode,iatom,ix}
      do icoord=1, ncoords
        deriv1(icoord,1:natoms,1:3) = reshape(coefs_svec(icoord,:,ipos1), (/natoms,3/))
        deriv2(icoord,1:natoms,1:3) = reshape(coefs_svec(icoord,:,ipos2), (/natoms,3/))
      enddo
      do iatom=1, natoms
        deriv(:,:,iatom) = deriv(:,:,iatom) + matmul(deriv1(:,iatom,:), transpose(deriv2(:,iatom,:)))
      enddo
    enddo ! ielem

    do icoord=1, ncoords
      do jcoord=1, ncoords
        coefs(icoord,jcoord,iterm) = dot_product(deriv(icoord,jcoord,:), inv_masses(:))
      enddo
    enddo

  enddo ! iterm

  if (allocated(ind_chr)) deallocate(ind_chr)

end subroutine deriv_gmat


subroutine expansion_terms_svec2cart(nmodes, nterms, terms)

  integer(ik), intent(in)                 :: nmodes
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 100
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
            write(out, '(/a/a,2(1x,i6))') &!
            'expansion_terms_svec2cart error: failed to allocate tmp(nmodes,terms_size)', 'nmodes, terms_size =', &!
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


subroutine expansion_terms_pseudo2svec(nmodes, nterms, terms)

  integer(ik), intent(in)                 :: nmodes
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 100
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
              write(out, '(/a/a,2(1x,i6))') &!
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


subroutine expansion_terms_complete(nvar, nterms, terms)

  integer(ik), intent(in)                 :: nvar
  integer(ik), intent(inout)              :: nterms
  integer(ik), allocatable, intent(inout) :: terms(:,:)

  integer(ik), parameter :: terms_incr = 100
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


function index_iarr1(ind, ind2) result(ipos)

  integer(ik), intent(in) :: ind(:), ind2(:,:)
  integer(ik) :: ipos, i

  ipos = 0
  do i=1, size(ind2,dim=2)
    if (all(ind==ind2(:,i))) then
      ipos = i
      exit
    endif
  enddo

end function index_iarr1


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
    write(out, '(/a/a,3(1x,i6))') &!
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


subroutine nmode_distribute(nmodes, nterms, terms, ncomb, comb, nterms_n, terms_n)

  integer(ik), intent(in) :: nmodes
  integer(ik), intent(in) :: nterms
  integer(ik), intent(in) :: terms(nmodes,nterms)
  integer(ik), intent(out) :: ncomb
  integer(ik), allocatable, intent(out) :: comb(:,:)
  integer(ik), allocatable, intent(out) :: nterms_n(:)
  integer(ik), allocatable, intent(out) :: terms_n(:,:,:)

  integer(ik), parameter :: ncomb_incr = 10, nterms_incr = 100
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
      write(out, '(/a/a,3(1x,i6))') &!
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


subroutine nmode_expansion(nmodes, nmax, ndeg, deg, mindeg, maxdeg, nterms, terms)

  integer(ik), intent(in)               :: nmodes
  integer(ik), intent(in)               :: nmax
  integer(ik), intent(in)               :: ndeg(nmodes)
  integer(ik), intent(in)               :: deg(maxval(ndeg),nmodes)
  integer(ik), intent(in)               :: mindeg(nmax)
  integer(ik), intent(in)               :: maxdeg(nmax)
  integer(ik), intent(out)              :: nterms
  integer(ik), allocatable, intent(out) :: terms(:,:)

  integer(ik) :: imode, info, iterm, iterm1, iterm2, nt, n, ncomb, icomb, nd(nmodes), d(maxval(ndeg),nmodes)
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


! Computes binomial coefficient n choose k
function bico(n, k)

  integer(ik), intent(in) :: k, n
  real(rk) :: bico

  bico = nint(exp(factln(n)-factln(k)-factln(n-k)))

end function bico


! Returns ln(n!)
function factln(n)

  integer(ik), intent(in) :: n
  real(rk) :: factln

  integer(ik), parameter :: nmax = 100
  real(rk), save :: a(0:nmax) = -1.0_rk

  if (n<0) then
    write(out, '(/a,1x,i4)') 'factln error: negative factorial=', n
    stop
  endif

  if (n<=nmax) then
    if (a(n+1)<0.0_rk) then
      a(n+1) = gammln(real(n+1,kind=rk))
    endif
    factln = a(n+1)
  else
    factln = gammln(real(n+1,kind=rk))
  endif

end function factln


! Returns ln(gamma(x))); call C lgamma function
function gammln(x)

  use iso_c_binding
  real(rk), intent(in) :: x
  real(rk) :: gammln

  interface
    real(c_double) function lgamma (y) bind(c)
    use iso_c_binding
    real(c_double), value :: y
    end function lgamma
  end interface

  gammln = lgamma(x)

end function gammln


subroutine tokenize_str(str0, n, word)

  character(*), intent(in) :: str0
  integer(ik), intent(out) :: n
  character(cl), intent(out) :: word(:)

  integer(ik) :: pos1
  character(len(str0)) :: str, str_

  pos1 = index(str0, '(')
  if (pos1>0) then
    str = trim(adjustl(str0(:pos1-1)))
  else
    str = trim(adjustl(str0))
  endif

  n = 0
  do
    pos1 = index(trim(str), ' ')
    if (pos1==0) then
       n = n + 1
       word(n) = trim(adjustl(str))
       exit
    endif
    n = n + 1
    word(n) = str(:pos1-1)
    str_ = trim(adjustl(str(pos1+1:)))
    str = str_
 enddo

end subroutine tokenize_str


function to_upper(str_in) result(str_out)

  character(*), intent(in) :: str_in
  character(len(str_in)) :: str_out
  integer :: i, j

  do i=1, len(str_in)
    j = iachar(str_in(i:i))
    if (j>=iachar('a').and.j<=iachar('z') ) then
      str_out(i:i) = achar(iachar(str_in(i:i))-32)
    else
      str_out(i:i) = str_in(i:i)
    endif
  enddo

end function to_upper

end module hamiltonian

