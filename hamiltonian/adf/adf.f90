module adf
use accuracy!, only: inp, out, cl, ik, rk, ark=>rk, small_, pi, avogno, vellgt, planck
implicit none


!private
!public adf_init, adf_reset_x0, adf_finalize


integer(ik), parameter :: adf_maxnterms = 30000
integer(ik), parameter :: adf_maxnvar = 12


type chrule_element_type
  integer(ik)              :: nprod
  integer(ik), allocatable :: order(:)
  real(ark)                :: coef
  integer(ik), allocatable :: pow(:)
  integer(ik), allocatable :: deriv(:,:)
  integer(ik), allocatable :: ind(:)
end type chrule_element_type


type chrule_type
  integer(ik)                            :: nelem
  type(chrule_element_type), allocatable :: chrule(:)
end type chrule_type


real(ark) :: small_number_check
real(ark) :: huge_number_check
real(ark) :: small_number



integer(ik)                    :: adf_nterms = 0
integer(ik)                    :: adf_nvar = 0
integer(ik), allocatable       :: adf_order(:)
integer(ik), allocatable       :: adf_terms(:,:)
integer(ik), allocatable       :: adf_terms_ivar(:,:)
type(chrule_type), allocatable :: adf_prod(:)
type(chrule_type), allocatable :: adf_intrf(:)
integer(ik)                    :: adf_iterm0 = 0


type, public :: adf_realq
  real(ark)                     :: v
  real(ark)                     :: d(adf_maxnterms)
  logical                       :: ifd0(adf_maxnterms)
  integer(ik)                   :: ivar(adf_maxnvar)
  character(len=:), allocatable :: sexpr
  character(cl)                 :: except = ''
end type adf_realq


!$omp threadprivate(adf_nterms,adf_nvar,adf_order,adf_terms,adf_terms_ivar,adf_prod,adf_intrf,adf_iterm0,small_number_check,small_number,huge_number_check)


interface operator (+)
  module procedure add_xy, add_xa, add_ax, &
                   add_xy_11, add_xy_22, &
                   add_xa_11, add_xa_22, &
                   add_ax_11, add_ax_22
end interface
public operator (+)


interface operator (-)
  module procedure subtr_xy, subtr_xa, subtr_ax, subtr_0x, &
                   subtr_xy_11, subtr_xy_22, subtr_0x_1, subtr_0x_2, &
                   subtr_xa_11, subtr_xa_22, &
                   subtr_ax_11, subtr_ax_22
end interface
public operator (-)


interface operator (*)
  module procedure prod_xy, prod_xa, prod_ax, prod_xi, prod_ix, &
                   prod_xy_11, prod_xa_10, prod_xa_20, prod_xa_11, prod_xi_10, prod_xi_20, prod_ax_11
end interface
public operator (*)


interface operator (/)
  module procedure div_xy, div_xa, div_ax, &
                   div_xy_10, div_xa_10, div_xa_20
end interface
public operator (/)


interface operator (**)
  module procedure ipow_x, qpow_x, ipow_x_10, ipow_x_11
end interface
public operator (**)


interface assignment (=)
  module procedure assign_xy, assign_xa, assign_ax, &
                   assign_xy_10, assign_xy_20, assign_xy_11, assign_xy_22, &
                   assign_xa_10, assign_xa_20, assign_xa_11, assign_xa_22, &
                   assign_ax_10, assign_ax_20, assign_ax_11, assign_ax_22
end interface
public assignment (=)


interface operator (==)
  module procedure eq_xy, eq_xa, eq_ax, &
                   eq_xy_10, eq_xy_20, eq_xy_11, eq_xy_22, &
                   eq_xa_10, eq_xa_20, eq_xa_11, eq_xa_22, &
                   eq_ax_10, eq_ax_20, eq_ax_11, eq_ax_22
end interface
public operator (==)


interface operator (<)
  module procedure lt_xy, lt_xa, lt_ax, &
                   lt_xy_10, lt_xy_20, lt_xy_11, lt_xy_22, &
                   lt_xa_10, lt_xa_20, lt_xa_11, lt_xa_22, &
                   lt_ax_10, lt_ax_20, lt_ax_11, lt_ax_22
end interface
public operator (<)


interface operator (>)
  module procedure gt_xy, gt_xa, gt_ax, &
                   gt_xy_10, gt_xy_20, gt_xy_11, gt_xy_22, &
                   gt_xa_10, gt_xa_20, gt_xa_11, gt_xa_22, &
                   gt_ax_10, gt_ax_20, gt_ax_11, gt_ax_22
end interface
public operator (>)


interface operator (<=)
  module procedure le_xy, le_xa, le_ax, &
                   le_xy_10, le_xy_20, le_xy_11, le_xy_22, &
                   le_xa_10, le_xa_20, le_xa_11, le_xa_22, &
                   le_ax_10, le_ax_20, le_ax_11, le_ax_22
end interface
public operator (<=)


interface operator (>=)
  module procedure ge_xy, ge_xa, ge_ax, &
                   ge_xy_10, ge_xy_20, ge_xy_11, ge_xy_22, &
                   ge_xa_10, ge_xa_20, ge_xa_11, ge_xa_22, &
                   ge_ax_10, ge_ax_20, ge_ax_11, ge_ax_22
end interface
public operator (>=)


interface abs
  module procedure abs_x, abs_x1, abs_x2
end interface
public abs


interface sum
  module procedure sum_x1
end interface
public sum


interface product
  module procedure product_x
end interface
public product


interface maxval
  module procedure maxval_x1, maxval_x2
end interface
public maxval


interface sqrt
  module procedure sqrt_x
end interface
public sqrt


interface cos
  module procedure cos_x
end interface
public cos


interface sin
  module procedure sin_x
end interface
public sin


interface exp
  module procedure exp_x
end interface
public exp


interface log
  module procedure log_x
end interface
public log


interface acos
  module procedure acos_x
end interface
public acos


interface asin
  module procedure asin_x
end interface
public asin


interface dot_product
  module procedure dot_product_xy, dot_product_xa, dot_product_ax
end interface
public dot_product


interface matmul
  module procedure matmul_xy, matmul_xa, matmul_ax, matmul_ax_1
end interface
public matmul


interface transpose
  module procedure transpose_x
end interface
public transpose


interface check_nan_inf
  module procedure check_nan_inf_fx, check_nan_inf_xy, check_nan_inf_xa, check_nan_inf_ax, check_nan_inf_xi, check_nan_inf_ix, &
                   check_nan_inf_dfx, check_nan_inf_dxy, check_nan_inf_dxa, check_nan_inf_dax, check_nan_inf_dxi, check_nan_inf_dix
end interface


interface expr
  module procedure expr_fx, expr_xy, expr_xa, expr_ax, expr_xi, expr_ix
end interface


contains


#include 'add.f90'
#include 'subtr.f90'
#include 'assign.f90'
#include 'div.f90'
#include 'prod.f90'
#include 'pow.f90'
#include 'dot_product.f90'
#include 'matmul.f90'
#include 'transpose.f90'
#include 'eq.f90'
#include 'lt.f90'
#include 'gt.f90'
#include 'le.f90'
#include 'ge.f90'
#include 'abs.f90'
#include 'sum.f90'
#include 'product.f90'
#include 'maxval.f90'
#include 'sqrt.f90'
#include 'cos.f90'
#include 'sin.f90'
#include 'exp.f90'
#include 'log.f90'
#include 'acos.f90'
#include 'asin.f90'



subroutine adf_init(nvar, x, x0, nterms, terms)

  integer(ik), intent(in)      :: nvar
  type(adf_realq), intent(out) :: x(nvar)
  real(ark), intent(in)        :: x0(nvar)
  integer(ik), intent(in)      :: nterms
  integer(ik), intent(in)      :: terms(nvar,nterms)

  integer(ik) :: iterm, info, ielem, ipos, ivar, iprod, term(nvar)
  character(cl) :: svar, sval

  small_number_check = real(epsilon(1.0d0), kind=ark)
  huge_number_check = real(huge(1.0d0), kind=ark)
  !small_number = real(1.0d-12, kind=ark)
  small_number = small_number_check

  adf_nterms = nterms
  if (nterms>adf_maxnterms) then
    write(out, '(/a,1x,i6)') 'adf_init error: number of terms exceeds the maximum =', adf_maxnterms
    stop
  endif

  adf_nvar = nvar
  if (nvar>adf_maxnvar) then
    write(out, '(/a,1x,i6)') 'adf_init error: number of variables exceeds the maximum =', adf_maxnvar
    stop
  endif

  if (allocated(adf_order)) deallocate(adf_order)
  if (allocated(adf_terms)) deallocate(adf_terms)
  if (allocated(adf_terms_ivar)) deallocate(adf_terms_ivar)
  allocate( adf_order(adf_nterms), adf_terms(nvar,adf_nterms), adf_terms_ivar(nvar,adf_nterms), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'adf_init error: failed to allocate adf_order(adf_nterms), adf_terms(nvar,adf_nterms), adf_terms_ivar(nvar,adf_nterms)', &
    'nvar, adf_nterms =', nvar, adf_nterms
    stop
  endif
  do iterm=1, adf_nterms
    adf_order(iterm) = sum(terms(1:nvar,iterm))
    adf_terms(1:nvar,iterm) = terms(1:nvar,iterm)
    where(adf_terms(1:nvar,iterm)>0)
      adf_terms_ivar(1:nvar,iterm) = 1
    elsewhere
      adf_terms_ivar(1:nvar,iterm) = 0
    endwhere
  enddo

  term = 0
  ipos = index_iarr1(term(1:nvar), adf_terms(1:nvar,1:adf_nterms))
  if (ipos==0) then
    write(out, '(/a)') 'adf_init error: zero-order derivative term is not found'
    stop
  else
    adf_iterm0 = ipos
  endif

  ! initialize active variables

  do ivar=1, nvar

    x(ivar)%v = x0(ivar)
    x(ivar)%d(1:adf_nterms) = 0.0

    term = 0
    ipos = index_iarr1(term(1:nvar), adf_terms(1:nvar,1:adf_nterms))
    if (ipos==0) then
      write(out, '(/a,1x,i3)') 'adf_init error: zero-order derivative term is not found for variable #', ivar
      stop
    else
      x(ivar)%d(ipos) = x0(ivar)
    endif

    if (any(adf_terms(ivar,:)>0)) then
      term = 0
      term(ivar) = 1
      ipos = index_iarr1(term(1:nvar), adf_terms(1:nvar,1:adf_nterms))
      if (ipos==0) then
        write(out, '(/a,1x,i3)') 'adf_init error: first-order derivative term is not found for variable #', ivar
        stop
      else
        x(ivar)%d(ipos) = 1.0_ark
      endif
    endif

#ifdef _ADF_OPT2_
    x(ivar)%ivar(:) = 0
    x(ivar)%ivar(ivar) = 1
#endif

#ifdef _ADF_OPT1_
    where(abs(x(ivar)%d(1:adf_nterms))<=small_number)
      x(ivar)%ifd0(1:adf_nterms) = .true.
    elsewhere
      x(ivar)%ifd0(1:adf_nterms) = .false.
    endwhere
#endif

    write(svar, '(i6)') ivar
    write(sval, '(es16.8)') x(ivar)%v
    svar = 'x'//trim(adjustl(svar))//'[='//trim(adjustl(sval))//']'
    if (allocated(x(ivar)%sexpr)) deallocate(x(ivar)%sexpr)
    allocate(character(len=len(trim(svar))) :: x(ivar)%sexpr)
    x(ivar)%sexpr = trim(adjustl(svar))

  enddo ! ivar

  ! allocate arrays to keep chain rule expressions

  if (allocated(adf_intrf)) deallocate(adf_intrf)
  allocate( adf_intrf(adf_nterms), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'adf_init error: failed to allocate adf_intrf(adf_nterms)', 'adf_nterms =', adf_nterms
    stop
  endif

  if (allocated(adf_prod)) deallocate(adf_prod)
  allocate( adf_prod(adf_nterms), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'adf_init error: failed to allocate adf_prod(adf_nterms)', 'adf_nterms =', adf_nterms
    stop
  endif

  ! generate chain rule expressions

  do iterm=1, adf_nterms

    call chrule_fx(nvar, adf_terms(1:nvar,iterm), adf_intrf(iterm)%nelem, adf_intrf(iterm)%chrule)
    do ielem=1, adf_intrf(iterm)%nelem
      do iprod=1, adf_intrf(iterm)%chrule(ielem)%nprod
        ipos = index_iarr1(adf_intrf(iterm)%chrule(ielem)%deriv(1:nvar,iprod), adf_terms)
        if (ipos==0) then
          write(out, '(/a,1x,<nvar>(1x,i3),1x,a)') 'adf_init error: failed to index derivative term = (', &
          adf_intrf(iterm)%chrule(ielem)%deriv(1:nvar,iprod), ') in the chain rule expression for f(x)'
          stop
        else
          adf_intrf(iterm)%chrule(ielem)%ind(iprod) = ipos
        endif
      enddo
      deallocate(adf_intrf(iterm)%chrule(ielem)%deriv)
    enddo

!#ifdef _ADF_SAVE_CHRULES_
!!$omp critical
!    write(1000,*) adf_terms(1:nvar,iterm), adf_intrf(iterm)%nelem
!    do ielem=1, adf_intrf(iterm)%nelem
!      write(1000,*) adf_intrf(iterm)%chrule(ielem)%nprod, (adf_intrf(iterm)%chrule(ielem)%ind(iprod), iprod=1, adf_intrf(iterm)%chrule(ielem)%nprod)
!    enddo
!!$omp end critical
!#endif

    call chrule_prod(nvar, adf_terms(1:nvar,iterm), adf_prod(iterm)%nelem, adf_prod(iterm)%chrule)
    do ielem=1, adf_prod(iterm)%nelem
      do iprod=1, adf_prod(iterm)%chrule(ielem)%nprod
        ipos = index_iarr1(adf_prod(iterm)%chrule(ielem)%deriv(1:nvar,iprod), adf_terms)
        if (ipos==0) then
          write(out, '(/a,1x,<nvar>(1x,i3),1x,a)') 'adf_init error: failed to index derivative term = (', &
          adf_prod(iterm)%chrule(ielem)%deriv(1:nvar,iprod), ') in the chain rule expression for multiplication operator'
          stop
        else
          adf_prod(iterm)%chrule(ielem)%ind(iprod) = ipos
        endif
      enddo
      deallocate(adf_prod(iterm)%chrule(ielem)%deriv)
    enddo

  enddo ! iterm

end subroutine adf_init


subroutine adf_reset_x0(nvar, x, x0)

  integer(ik), intent(in)      :: nvar
  type(adf_realq), intent(out) :: x(nvar)
  real(ark), intent(in)        :: x0(nvar)

  integer(ik) :: iterm, info, ielem, ipos, ivar, iprod, term(nvar)
  character(cl) :: svar, sval

  ! initialize active variables

  do ivar=1, nvar

    x(ivar)%v = x0(ivar)
    x(ivar)%d(1:adf_nterms) = 0.0

    term = 0
    ipos = index_iarr1(term(1:nvar), adf_terms(1:nvar,1:adf_nterms))
    if (ipos==0) then
      write(out, '(/a,1x,i3)') 'adf_reset_x0 error: zero-order derivative term is not found for variable #', ivar
      stop
    else
      x(ivar)%d(ipos) = x0(ivar)
    endif

    if (any(adf_terms(ivar,:)>0)) then
      term = 0
      term(ivar) = 1
      ipos = index_iarr1(term(1:nvar), adf_terms(1:nvar,1:adf_nterms))
      if (ipos==0) then
        write(out, '(/a,1x,i3)') 'adf_reset_x0 error: first-order derivative term is not found for variable #', ivar
        stop
      else
        x(ivar)%d(ipos) = 1.0_ark
      endif
    endif

#ifdef _ADF_OPT2_
    x(ivar)%ivar(:) = 0
    x(ivar)%ivar(ivar) = 1
#endif

#ifdef _ADF_OPT1_
    where(abs(x(ivar)%d(1:adf_nterms))<=small_number)
      x(ivar)%ifd0(1:adf_nterms) = .true.
    elsewhere
      x(ivar)%ifd0(1:adf_nterms) = .false.
    endwhere
#endif

    write(svar, '(i6)') ivar
    write(sval, '(es16.8)') x(ivar)%v
    svar = 'x'//trim(adjustl(svar))//'[='//trim(adjustl(sval))//']'
    if (allocated(x(ivar)%sexpr)) deallocate(x(ivar)%sexpr)
    allocate(character(len=len(trim(svar))) :: x(ivar)%sexpr)
    x(ivar)%sexpr = trim(adjustl(svar))

  enddo ! ivar

end subroutine adf_reset_x0


subroutine adf_finalize

  if (allocated(adf_order)) deallocate(adf_order)
  if (allocated(adf_terms)) deallocate(adf_terms)
  if (allocated(adf_terms_ivar)) deallocate(adf_terms_ivar)
  if (allocated(adf_prod)) deallocate(adf_prod)
  if (allocated(adf_intrf)) deallocate(adf_intrf)

  adf_nterms = 0
  adf_nvar = 0
  adf_iterm0 = 0

end subroutine adf_finalize


subroutine adf_set_var(var, dcoefs)

  type(adf_realq), intent(inout) :: var
  real(ark), intent(in) :: dcoefs(:)

  integer(ik) :: iterm, ivar

  if (size(dcoefs)/=adf_nterms) then
    write(out, '(/a,1x,i6,1x,a,1x,i6)') 'adf_set_var error: size of array "dcoefs" =', size(dcoefs), 'is not equal to the number ADF terms =', adf_nterms
    stop
  endif

  var%v = dcoefs(adf_iterm0)
  do iterm=1, adf_nterms
    var%d(iterm) = dcoefs(iterm)
  enddo

#ifdef _ADF_OPT1_
    where(abs(var%d(1:adf_nterms))<=small_number)
      var%ifd0(1:adf_nterms) = .true.
    elsewhere
      var%ifd0(1:adf_nterms) = .false.
    endwhere
#endif

#ifdef _ADF_OPT2_
    do ivar=1, adf_nvar
      if (all(adf_terms(ivar,:)==0)) then
        var%ivar(ivar) = 0
      else
        var%ivar(ivar) = 1
      endif
    enddo
#endif

end subroutine adf_set_var


! Generates chain rule expression for derivative of a f(x)
!
subroutine chrule_fx(nvar, term, nelem, chrule)

  integer(ik), intent(in)                             :: nvar
  integer(ik), intent(in)                             :: term(nvar)
  integer(ik), intent(out)                            :: nelem
  type(chrule_element_type), allocatable, intent(out) :: chrule(:)

  integer(ik), parameter :: max_nelem = 1000000, max_nprod = 20
  integer(ik) :: ivar, iorder, ielem, jelem, ind(nvar), ipos, nelem0, iprod, info, np, jprod
  integer(ik) :: nprod(max_nelem), nprod0(max_nelem), pow(max_nprod,max_nelem), pow0(max_nprod,max_nelem), &
                 deriv(nvar,max_nprod,max_nelem), deriv0(nvar,max_nprod,max_nelem), imatch(max_nprod)
  real(ark) :: coef(max_nelem), coef0(max_nelem)

  nelem0 = 1
  nprod0(1) = 1
  coef0(1) = 1.0_ark
  pow0(1,1) = 1
  deriv0(:,1,1) = 0

  do ivar=1, nvar
    do iorder=1, term(ivar)

      jelem = 0
      do ielem=1, nelem0

        ! derivative of f(x): f(x)' = f'(x) * x'

        jelem = jelem + 1
        call check_nelem(jelem)
        nprod(jelem) = nprod0(ielem)
        pow(1:nprod(jelem),jelem) = pow0(1:nprod0(ielem),ielem)
        coef(jelem) = coef0(ielem)
        deriv(:,1:nprod(jelem),jelem) = deriv0(:,1:nprod0(ielem),ielem)

        deriv(ivar,1,jelem) = deriv(ivar,1,jelem) + 1
        ind = 0; ind(ivar)=1
        if (nprod(jelem)>1) then
          ipos = index_iarr1(ind, deriv(:,2:nprod(jelem),jelem))
        else
          ipos = 0
        endif
        if (ipos==0) then
          nprod(jelem) = nprod(jelem) + 1
          call check_nprod(nprod(jelem))
          pow(nprod(jelem),jelem) = 1
          deriv(:,nprod(jelem),jelem) = ind
        else
          pow(ipos+1,jelem) = pow(ipos+1,jelem) + 1
        endif

        ! derivative of product of x-derivatives

        do iprod=2, nprod0(ielem)

          jelem = jelem + 1
          call check_nelem(jelem)
          nprod(jelem) = nprod0(ielem)
          pow(1:nprod(jelem),jelem) = pow0(1:nprod0(ielem),ielem)
          coef(jelem) = coef0(ielem)
          deriv(:,1:nprod(jelem),jelem) = deriv0(:,1:nprod0(ielem),ielem)

          if (pow(iprod,jelem)==1) then
            deriv(ivar,iprod,jelem) = deriv(ivar,iprod,jelem) + 1
          elseif (pow(iprod,jelem)>1) then
            coef(jelem) = coef(jelem) * pow(iprod,jelem)
            pow(iprod,jelem) = pow(iprod,jelem) - 1
            ind = deriv(:,iprod,jelem)
            ind(ivar) = ind(ivar) + 1
            if (nprod(jelem)>1) then
              ipos = index_iarr1(ind, deriv(:,2:nprod(jelem),jelem))
            else
              ipos = 0
            endif
            if (ipos==0) then
              nprod(jelem) = nprod(jelem) + 1
              call check_nprod(nprod(jelem))
              pow(nprod(jelem),jelem) = 1
              deriv(:,nprod(jelem),jelem) = ind
            else
              pow(ipos+1,jelem) = pow(ipos+1,jelem) + 1
            endif
          else
            write(out, '(/a)') 'chrule_fx error: x**0 must not appear in the chain rule expression for derivative of f(x)'
            stop
          endif

        enddo ! iprod

      enddo ! ielem

      ! sum up same elements

      nelem = jelem
      nelem0 = 0
      do ielem=1, nelem
        ipos = 0
        do jelem=1, nelem0
          if (nprod(ielem)==nprod0(jelem).and.all(deriv(1:nvar,1,ielem)==deriv0(1:nvar,1,jelem)).and.pow(1,ielem)==pow0(1,jelem)) then
            imatch = 0
            imatch(1) = 1
            do iprod=2, nprod(ielem)
              do jprod=2, nprod(ielem)
                if (any(imatch(1:iprod-1)==jprod)) cycle
                if (all(deriv(1:nvar,iprod,ielem)==deriv0(1:nvar,jprod,jelem)).and.pow(iprod,ielem)==pow0(jprod,jelem)) then
                  imatch(iprod) = jprod
                  exit
                endif
              enddo
            enddo
            if (all(imatch(1:nprod(ielem))>0)) then
              ipos = jelem
              exit
            endif
          endif
        enddo
        if (ipos>0) then
          coef0(ipos) = coef0(ipos) + coef(ielem)
        else
          nelem0 = nelem0 + 1
          nprod0(nelem0) = nprod(ielem)
          coef0(nelem0) = coef(ielem)
          pow0(1:nprod0(nelem0),nelem0) = pow(1:nprod(ielem),ielem)
          deriv0(1:nvar,1:nprod0(nelem0),nelem0) = deriv(1:nvar,1:nprod(ielem),ielem)
        endif
      enddo

    enddo ! iorder
  enddo ! ivar

  ! copy the result into "chrule" variable

  nelem = nelem0
  if (allocated(chrule)) deallocate(chrule)
  allocate( chrule(nelem), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'chrule_fx error: failed to allocate chrule', 'nelem =', nelem
    stop
  endif
  do ielem=1, nelem
    chrule(ielem)%coef = coef0(ielem)
    chrule(ielem)%nprod = nprod0(ielem)
    np = chrule(ielem)%nprod
    allocate( chrule(ielem)%pow(np), chrule(ielem)%deriv(nvar,np), chrule(ielem)%ind(np), chrule(ielem)%order(np), stat=info )
    if (info/=0) then
      write(out, '(/a/a,10(1x,i8))') 'chrule_fx error: failed to allocate chrule%pow, chrule%deriv, chrule%ind', 'nelem, np =', nelem, np
      stop
    endif
    chrule(ielem)%pow(1:np) = pow0(1:np,ielem)
    chrule(ielem)%deriv(1:nvar,1:np) = deriv0(1:nvar,1:np,ielem)
    do iprod=1, np
      chrule(ielem)%order(iprod) = sum(deriv0(1:nvar,iprod,ielem))
    enddo
    chrule(ielem)%ind = 0
  enddo

  contains

  subroutine check_nelem(ielem)
    integer(ik) :: ielem
    if (ielem>max_nelem) then
      write(out, '(/a,1x,i8)') 'chrule_fx error: number of elements in the chain rule expression exceeds maximum =', max_nelem
      write(out, '(a)') 'increase max_nelem in adf.f90/chrule_fx and adf.f90/chrule_prod'
      stop
    endif
  end subroutine check_nelem

  subroutine check_nprod(iprod)
    integer(ik) :: iprod
    if (iprod>max_nprod) then
      write(out, '(/a,1x,i8)') 'chrule_fx error: number of product-elements in the chain rule expression exceeds maximum =', max_nprod
      write(out, '(a)') 'increase max_nprod in adf.f90/chrule_fx'
      stop
    endif
  end subroutine check_nprod

end subroutine chrule_fx



! Generates chain rule expression for derivative of f(x)*g(x)
!
subroutine chrule_prod(nvar, term, nelem, chrule)

  integer(ik), intent(in)                             :: nvar
  integer(ik), intent(in)                             :: term(nvar)
  integer(ik), intent(out)                            :: nelem
  type(chrule_element_type), allocatable, intent(out) :: chrule(:)


  integer(ik), parameter :: max_nelem = 1000000, max_nprod = 2
  integer(ik) :: ivar, iorder, ielem, jelem, ind(nvar), ipos, nelem0, iprod, info, np, jprod, nelem2
  integer(ik) :: nprod(max_nelem), nprod0(max_nelem), pow(max_nprod,max_nelem), pow0(max_nprod,max_nelem), &
                 deriv(nvar,max_nprod,max_nelem), deriv0(nvar,max_nprod,max_nelem), imatch(max_nprod)
  real(ark) :: coef(max_nelem), coef0(max_nelem)

  nelem0 = 1
  nprod0(1) = 2
  coef0(1) = 1.0_ark
  pow0(1:2,1) = 1
  deriv0(1:nvar,1:2,1) = 0

  do ivar=1, nvar
    do iorder=1, term(ivar)

      jelem = 0
      do ielem=1, nelem0
        do iprod=1, nprod0(ielem)
          jelem = jelem + 1
          call check_nelem(jelem)
          nprod(jelem) = nprod0(ielem)
          pow(1:nprod(jelem),jelem) = pow0(1:nprod0(ielem),ielem)
          coef(jelem) = coef0(ielem)
          deriv(1:nvar,1:nprod(jelem),jelem) = deriv0(1:nvar,1:nprod0(ielem),ielem)
          deriv(ivar,iprod,jelem) = deriv(ivar,iprod,jelem) + 1
        enddo
      enddo

      ! sum up same elements

      nelem2 = jelem
      nelem0 = 0
      do ielem=1, nelem2
        ipos = 0
        do jelem=1, nelem0
          if (nprod(ielem)==nprod0(jelem)) then
            if (all(deriv(1:nvar,1:nprod(ielem),ielem)==deriv0(1:nvar,1:nprod0(jelem),jelem))) then
              ipos = jelem
              exit
            endif
          endif
        enddo
        if (ipos>0) then
          coef0(ipos) = coef0(ipos) + coef(ielem)
        else
          nelem0 = nelem0 + 1
          nprod0(nelem0) = nprod(ielem)
          coef0(nelem0) = coef(ielem)
          pow0(1:nprod0(nelem0),nelem0) = pow(1:nprod(ielem),ielem)
          deriv0(1:nvar,1:nprod0(nelem0),nelem0) = deriv(1:nvar,1:nprod(ielem),ielem)
        endif
      enddo

    enddo ! iorder
  enddo ! ivar

  ! copy the result into "chrule" variable

  nelem = nelem0
  if (allocated(chrule)) deallocate(chrule)
  allocate( chrule(nelem), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'chrule_prod error: failed to allocate chrule', 'nelem =', nelem
    stop
  endif
  do ielem=1, nelem
    chrule(ielem)%coef = coef0(ielem)
    chrule(ielem)%nprod = nprod0(ielem)
    np = chrule(ielem)%nprod
    allocate( chrule(ielem)%pow(np), chrule(ielem)%deriv(nvar,np), chrule(ielem)%ind(np), chrule(ielem)%order(np), stat=info )
    if (info/=0) then
      write(out, '(/a/a,10(1x,i8))') 'chrule_prod error: failed to allocate chrule%pow, chrule%deriv, chrule%ind', 'nelem, np =', nelem, np
      stop
    endif
    chrule(ielem)%pow(1:np) = pow0(1:np,ielem)
    chrule(ielem)%deriv(1:nvar,1:np) = deriv0(1:nvar,1:np,ielem)
    do iprod=1, np
      chrule(ielem)%order(iprod) = sum(deriv0(1:nvar,iprod,ielem))
    enddo
    chrule(ielem)%ind = 0
  enddo

  contains

  subroutine check_nelem(ielem)
    integer(ik) :: ielem
    if (ielem>max_nelem) then
      write(out, '(/a,1x,i8)') 'chrule_prod error: number of elements in the chain rule expression exceeds maximum =', max_nelem
      write(out, '(a)') 'increase max_nelem in adf.f90/chrule_fx and adf.f90/chrule_prod'
      stop
    endif
  end subroutine check_nelem

end subroutine chrule_prod



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



subroutine check_nan_inf_fx(ftest, sfunc, x)

  real(ark), intent(in) :: ftest
  character(*), intent(in) :: sfunc
  type(adf_realq), intent(in) :: x

  if (ftest/=ftest) then
    write(out, '(/a,1x,a,a,1x,es16.8)') 'function', trim(sfunc), '(x) = NaN at x =', x%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,1x,a,a,1x,es16.8)') 'function', trim(sfunc), '(x) = Infinity at x =', x%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_fx



subroutine check_nan_inf_dfx(deriv, ftest, sfunc, x)

  integer(ik), intent(in) :: deriv(:)
  real(ark), intent(in) :: ftest
  character(*), intent(in) :: sfunc
  type(adf_realq), intent(in) :: x

  integer(ik) :: nvar

  nvar = size(deriv)

  if (ftest/=ftest) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,1x,a,a,1x,es16.8)') '(', deriv, ') derivative of function', trim(sfunc), '(x) = NaN at x =', x%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,1x,a,a,1x,es16.8)') '(', deriv, ') derivative of function', trim(sfunc), '(x) = Infinity at x =', x%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_dfx



subroutine check_nan_inf_xy(ftest, soper, x, y)

  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x, y

  if (ftest/=ftest) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,es16.8)') 'x', trim(soper), 'y = NaN at x =', x%v, 'and y =', y%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'y =', trim(y%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,es16.8)') 'x', trim(soper), 'y = Infinity at x =', x%v, 'and y =', y%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'y =', trim(y%sexpr)
    stop
  endif

end subroutine check_nan_inf_xy



subroutine check_nan_inf_dxy(deriv, ftest, soper, x, y)

  integer(ik), intent(in) :: deriv(:)
  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x, y

  integer(ik) :: nvar

  nvar = size(deriv)

  if (ftest/=ftest) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,es16.8)') '(', deriv, ') derivative of x', trim(soper), 'y = NaN at x =', x%v, 'and y =', y%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'y =', trim(y%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,es16.8)') '(', deriv, ') derivative of x', trim(soper), 'y = Infinity at x =', x%v, 'and y =', y%v
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    write(out, '(a,1x,a)') 'y =', trim(y%sexpr)
    stop
  endif

end subroutine check_nan_inf_dxy



subroutine check_nan_inf_xa(ftest, soper, x, a)

  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: a

  if (ftest/=ftest) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,es16.8)') 'x', trim(soper), 'a = NaN at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,es16.8)') 'x', trim(soper), 'a = Infinity at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_xa



subroutine check_nan_inf_dxa(deriv, ftest, soper, x, a)

  integer(ik), intent(in) :: deriv(:)
  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: a

  integer(ik) :: nvar

  nvar = size(deriv)

  if (ftest/=ftest) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,es16.8)') '(', deriv, ') derivative of x', trim(soper), 'a = NaN at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,es16.8)') '(', deriv, ') derivative of x', trim(soper), 'a = Infinity at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_dxa



subroutine check_nan_inf_ax(ftest, soper, a, x)

  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  real(ark), intent(in) :: a
  type(adf_realq), intent(in) :: x

  if (ftest/=ftest) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,es16.8)') 'a', trim(soper), 'x = NaN at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,es16.8)') 'a', trim(soper), 'x = Infinity at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_ax



subroutine check_nan_inf_dax(deriv, ftest, soper, a, x)

  integer(ik), intent(in) :: deriv(:)
  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  real(ark), intent(in) :: a
  type(adf_realq), intent(in) :: x

  integer(ik) :: nvar

  nvar = size(deriv)

  if (ftest/=ftest) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,es16.8)') '(', deriv, ') derivative of a', trim(soper), 'x = NaN at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,es16.8)') '(', deriv, ') derivative of a', trim(soper), 'x = Infinity at x =', x%v, 'and a =', a
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_dax



subroutine check_nan_inf_xi(ftest, soper, x, i)

  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  integer(ik), intent(in) :: i

  if (ftest/=ftest) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,i6)') 'x', trim(soper), 'i = NaN at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,i6)') 'x', trim(soper), 'i = Infinity at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_xi



subroutine check_nan_inf_dxi(deriv, ftest, soper, x, i)

  integer(ik), intent(in) :: deriv(:)
  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  integer(ik), intent(in) :: i

  integer(ik) :: nvar

  nvar = size(deriv)

  if (ftest/=ftest) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,i6)') '(', deriv, ') derivative of x', trim(soper), 'i = NaN at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,i6)') '(', deriv, ') derivative of x', trim(soper), 'i = Infinity at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_dxi



subroutine check_nan_inf_ix(ftest, soper, i, x)

  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  integer(ik), intent(in) :: i

  if (ftest/=ftest) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,i6)') 'i', trim(soper), 'x = NaN at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,a,a,1x,es16.8,1x,a,1x,i6)') 'i', trim(soper), 'x = Infinity at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_ix



subroutine check_nan_inf_dix(deriv, ftest, soper, i, x)

  integer(ik), intent(in) :: deriv(:)
  real(ark), intent(in) :: ftest
  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  integer(ik), intent(in) :: i

  integer(ik) :: nvar

  nvar = size(deriv)

  if (ftest/=ftest) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,i6)') '(', deriv, ') derivative of i', trim(soper), 'x = NaN at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

  if (abs(ftest)>=huge_number_check) then
    write(out, '(/a,<nvar>(1x,i3),1x,a,a,a,1x,es16.8,1x,a,1x,i6)') '(', deriv, ') derivative of i', trim(soper), 'x = Infinity at x =', x%v, 'and i =', i
    write(out, '(a,1x,a)') 'x =', trim(x%sexpr)
    stop
  endif

end subroutine check_nan_inf_dix



subroutine expr_fx(sfunc, x, f)

  character(*), intent(in) :: sfunc
  type(adf_realq), intent(in) :: x
  type(adf_realq), intent(inout) :: f

  integer(ik) :: l

  l = len(trim(x%sexpr)) + len(trim(adjustl(sfunc))) + 2
  if (allocated(f%sexpr)) deallocate(f%sexpr)
  allocate(character(len=l) :: f%sexpr)
  f%sexpr = trim(adjustl(sfunc))//'('//trim(x%sexpr)//')'

end subroutine expr_fx



subroutine expr_xy(soper, x, y, f)

  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x, y
  type(adf_realq), intent(inout) :: f

  integer(ik) :: l

  l = len(trim(x%sexpr)) + len(trim(adjustl(soper))) + len(trim(y%sexpr)) + 2
  if (allocated(f%sexpr)) deallocate(f%sexpr)
  allocate(character(len=l) :: f%sexpr)
  f%sexpr = '('//trim(x%sexpr)//trim(adjustl(soper))//trim(y%sexpr)//')'

end subroutine expr_xy



subroutine expr_xa(soper, x, a, f)

  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  real(ark), intent(in) :: a
  type(adf_realq), intent(inout) :: f

  integer(ik) :: l
  character(cl) :: sa

  write(sa, '(es16.8)') a
  l = len(trim(x%sexpr)) + len(trim(adjustl(soper))) + len(trim(adjustl(sa))) + 2
  if (allocated(f%sexpr)) deallocate(f%sexpr)
  allocate(character(len=l) :: f%sexpr)
  f%sexpr = '('//trim(x%sexpr)//trim(adjustl(soper))//trim(adjustl(sa))//')'

end subroutine expr_xa



subroutine expr_ax(soper, a, x, f)

  character(*), intent(in) :: soper
  real(ark), intent(in) :: a
  type(adf_realq), intent(in) :: x
  type(adf_realq), intent(inout) :: f

  integer(ik) :: l
  character(cl) :: sa

  write(sa, '(es16.8)') a
  l = len(trim(adjustl(sa))) + len(trim(adjustl(soper))) + len(trim(x%sexpr)) + 2
  if (allocated(f%sexpr)) deallocate(f%sexpr)
  allocate(character(len=l) :: f%sexpr)
  f%sexpr = '('//trim(adjustl(sa))//trim(adjustl(soper))//trim(x%sexpr)//')'

end subroutine expr_ax



subroutine expr_xi(soper, x, i, f)

  character(*), intent(in) :: soper
  type(adf_realq), intent(in) :: x
  integer(ik), intent(in) :: i
  type(adf_realq), intent(inout) :: f

  integer(ik) :: l
  character(cl) :: si

  write(si, '(i6)') i
  l = len(trim(x%sexpr)) + len(trim(adjustl(soper))) + len(trim(adjustl(si))) + 2
  if (allocated(f%sexpr)) deallocate(f%sexpr)
  allocate(character(len=l) :: f%sexpr)
  f%sexpr = '('//trim(x%sexpr)//trim(adjustl(soper))//trim(adjustl(si))//')'

end subroutine expr_xi



subroutine expr_ix(soper, i, x, f)

  character(*), intent(in) :: soper
  integer(ik), intent(in) :: i
  type(adf_realq), intent(in) :: x
  type(adf_realq), intent(inout) :: f

  integer(ik) :: l
  character(cl) :: si

  write(si, '(i6)') i
  l = len(trim(adjustl(si))) + len(trim(adjustl(soper))) + len(trim(x%sexpr)) + 2
  if (allocated(f%sexpr)) deallocate(f%sexpr)
  allocate(character(len=l) :: f%sexpr)
  f%sexpr = '('//trim(adjustl(si))//trim(adjustl(soper))//trim(x%sexpr)//')'

end subroutine expr_ix



end module adf
