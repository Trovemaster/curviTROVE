subroutine potref_general(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: nmodes, rank, irank, nparams, info, iparam
  type(HM_func_type), allocatable :: func0


  nmodes = molec%nmodes


  ! copy func into func0

  allocate(func0, stat=info)
  if (info/=0) then
    write(out, '(/a)') 'potref_general error: failed to allocate func0'
    stop
  endif

  func0%func_type = func%ref_poten_type
  func0%coord_nparams = func%coord_nparams

  allocate(func0%coord_type(nmodes), func0%coord_params(nmodes,func0%coord_nparams), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'potref_general error: failed to allocate func0%coord_type(nmodes), func0%coord_params(nmodes,func0%coord_nparams)', &!
    'nmodes, func0%coord_nparams =', nmodes, func0%coord_nparams
    stop
  endif

  func0%coord_type(1:nmodes) = func%coord_type(1:nmodes)
  func0%coord_params(1:nmodes,1) = func%coord_params(1:nmodes,1)


  nparams = func%rank
  rank = nparams

  allocate(func0%nparams(1), func0%iparams(nmodes,nparams,1), func0%params(nparams,1), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'potref_general error: failed to allocate func0%nparams(1), func0%iparams(nmodes,nparams,1), func0%params(nparams,1)', &!
    'nparams, nmodes =', nparams, nmodes
    stop
  endif

  func0%rank = 1
  func0%nparams(1) = nparams
  func0%iparams(1:nmodes,1:nparams,1) = func%iparams(1:nmodes,1,1:rank)
  func0%rotate = .false.
  do iparam=1, nparams
    irank = iparam
    if (func%ifit(1,irank)==0) then
      func0%params(iparam,1) = func%params(1,irank)
    else
      func0%params(iparam,1) = 0
    endif
  enddo


  ! compute derivatives of potential wrt expansion coefficients

  f = 0.0_ark

  do iparam=1, nparams

    irank = iparam

    if (func%ifit(1,irank)==0) cycle

    func0%params(iparam,1) = 1.0_ark

    select case(trim(func0%func_type))
#ifdef _MOL_ZXY3_
    case('POTEN_CH3CL_SYM','POTEN_ZXY3_SYM')
      call poten_ch3cl_sym_ADF(molec, func0, internal, f(iparam), cart)
#endif
#ifdef _MOL_XY4_
    case('POTEN_XY4_ALPHA')
      call poten_xy4_alpha_ADF_noreexp(molec, func0, internal, f(iparam), cart)
#endif
#ifdef _MOL_C2H4_
    case('POTEN_C2H4_886666')
      call  poten_c2h4_886666_ADF(molec, func0, internal, f(iparam), cart)
    !  call  poten_c2h4_886666_ADF_noreexp(molec, func0, internal, f(iparam), cart)
#endif
    case default
      write(out, '(/a,a,a)') 'potref_general error: unknown type of potential function = "', trim(func0%func_type), '"'
      stop
    end select

    func0%params(iparam,1) = 0

  enddo

  deallocate(func0)

end subroutine potref_general
