subroutine poten_h2o2_koput_ADF(molec, func, internal, f, cart)
  use adf
  implicit none

  type(HM_molec_type), intent(in) :: molec
  type(HM_func_type), intent(in)  :: func
  type(adf_realq), intent(in)     :: internal(molec%nmodes)
  type(adf_realq), intent(out)    :: f(func%rank)
  type(adf_realq), intent(in), optional :: cart(molec%natoms,3)

  integer(ik) :: irank, i, icoord, k_ind(6)
  real(ark) :: e1, e2, pd, e4, e6, a1, a2, tocm=219474.63067_ark
  type(adf_realq) :: r0, r1, r2, alpha1, alpha2, tau, q(6), y(6)

  irank = 1

  ! internal coordinates

  if (trim(molec%coord_transform)=='H2O2_RALPHATAU') then

    r0     = internal(1)
    r1     = internal(2)
    r2     = internal(3)
    alpha1 = internal(4)
    alpha2 = internal(5)
    tau    = internal(6)

  else

    write(out, '(/a,a,a)') 'poten_h2o2_koput_ADF error: coordinate type = "', trim(molec%coord_transform), '" is not supported'
    stop

  endif

  ! expansion functions

  pd = real(pi,ark)/180.0_ark
  e1 = func%params(1,irank)
  e2 = func%params(2,irank)
  e4 = func%params(3,irank)*pd
  e6 = func%params(4,irank)*pd

  q(1) = (r0-e1)/r0
  q(2) = (r1-e2)/r1
  q(3) = (r2-e2)/r2
  q(4) = (alpha1-e4)
  q(5) = (alpha2-e4)
  q(6) = tau

  f = 0.0_ark

  do i=5, func%nparams(irank)

    k_ind(1) = func%iparams(1,i,irank)
    k_ind(2) = func%iparams(2,i,irank)
    k_ind(3) = func%iparams(3,i,irank)
    k_ind(4) = func%iparams(4,i,irank)
    k_ind(5) = func%iparams(5,i,irank)
    k_ind(6) = func%iparams(6,i,irank)

    y(1:6) = 1.0_ark

    do icoord=1, 5
      if (k_ind(icoord)==1) then
        y(icoord) = q(icoord)
      elseif (k_ind(icoord)>1) then
        y(icoord) = q(icoord)**k_ind(icoord)
      endif
    enddo

    if (k_ind(6).ne.0) y(6) = cos(real(k_ind(6),ark)*q(6))

    f(irank) = f(irank) + func%params(i,irank) * product(y(1:6))

  enddo

  f = f * tocm

end subroutine poten_h2o2_koput_ADF
