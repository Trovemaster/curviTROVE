subroutine assign_xy(x, y)

  type(adf_realq), intent(inout) :: x
  type(adf_realq), intent(in) :: y

  integer(ik) :: iterm

  x%v = y%v

#ifdef _ADF_OPT2_
  where(y%ivar(1:adf_nvar)==1)
    x%ivar(1:adf_nvar) = 1
  elsewhere
    x%ivar(1:adf_nvar) = 0
  endwhere
#endif

  x%d(1:adf_nterms) = y%d(1:adf_nterms)

  x%except = y%except

  if (allocated(x%sexpr)) deallocate(x%sexpr)
  allocate(character(len=len(trim(y%sexpr))) :: x%sexpr)
  x%sexpr = y%sexpr

#ifdef _ADF_OPT1_
  where(abs(x%d(1:adf_nterms))<=small_number)
    x%ifd0(1:adf_nterms) = .true.
  elsewhere
    x%ifd0(1:adf_nterms) = .false.
  endwhere
#endif

end subroutine assign_xy



subroutine assign_xy_10(x, y)

  type(adf_realq), intent(inout) :: x(:)
  type(adf_realq), intent(in) :: y

  integer(ik) :: i

  do i=1, size(x)
    x(i) = y
  enddo

end subroutine assign_xy_10



subroutine assign_xy_20(x, y)

  type(adf_realq), intent(inout) :: x(:,:)
  type(adf_realq), intent(in) :: y

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      x(j,i) = y
    enddo
  enddo

end subroutine assign_xy_20



subroutine assign_xy_11(x, y)

  type(adf_realq), intent(inout) :: x(:)
  type(adf_realq), intent(in) :: y(:)

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'assign_xy_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    x(i) = y(i)
  enddo

end subroutine assign_xy_11



subroutine assign_xy_22(x, y)

  type(adf_realq), intent(inout) :: x(:,:)
  type(adf_realq), intent(in) :: y(:,:)

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'assign_xy_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      x(j,i) = y(j,i)
    enddo
  enddo

end subroutine assign_xy_22



subroutine assign_xa(x, y)

  type(adf_realq), intent(inout) :: x
  real(ark), intent(in) :: y

  integer(ik) :: iterm
  character(cl) :: sy

  x%v = y

#ifdef _ADF_OPT2_
  x%ivar(1:adf_nvar) = 0
#endif

  !omp parallel do private(iterm) schedule(dynamic)
  do iterm=1, adf_nterms
    if (adf_order(iterm)==0) then
      x%d(iterm) = x%v
    else
      x%d(iterm) = 0.0
    endif
  enddo
  !omp end parallel do

  write(sy, '(es16.8)') y
  if (allocated(x%sexpr)) deallocate(x%sexpr)
  allocate(character(len=len(trim(adjustl(sy)))) :: x%sexpr)
  x%sexpr = trim(adjustl(sy))

#ifdef _ADF_OPT1_
  where(abs(x%d(1:adf_nterms))<=small_number)
    x%ifd0(1:adf_nterms) = .true.
  elsewhere
    x%ifd0(1:adf_nterms) = .false.
  endwhere
#endif

end subroutine assign_xa



subroutine assign_xa_10(x, y)

  type(adf_realq), intent(inout) :: x(:)
  real(ark), intent(in) :: y

  integer(ik) :: i

  do i=1, size(x)
    x(i) = y
  enddo

end subroutine assign_xa_10



subroutine assign_xa_20(x, y)

  type(adf_realq), intent(inout) :: x(:,:)
  real(ark), intent(in) :: y

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      x(j,i) = y
    enddo
  enddo

end subroutine assign_xa_20



subroutine assign_xa_11(x, y)

  type(adf_realq), intent(inout) :: x(:)
  real(ark), intent(in) :: y(:)

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'assign_xa_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    x(i) = y(i)
  enddo

end subroutine assign_xa_11



subroutine assign_xa_22(x, y)

  type(adf_realq), intent(inout) :: x(:,:)
  real(ark), intent(in) :: y(:,:)

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'assign_xa_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      x(j,i) = y(j,i)
    enddo
  enddo

end subroutine assign_xa_22



subroutine assign_ax(x, y)

  real(ark), intent(inout) :: x
  type(adf_realq), intent(in) :: y

  integer(ik) :: iterm

  x = y%v

end subroutine assign_ax



subroutine assign_ax_10(x, y)

  real(ark), intent(inout) :: x(:)
  type(adf_realq), intent(in) :: y

  integer(ik) :: i

  do i=1, size(x)
    x(i) = y
  enddo

end subroutine assign_ax_10



subroutine assign_ax_20(x, y)

  real(ark), intent(inout) :: x(:,:)
  type(adf_realq), intent(in) :: y

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      x(j,i) = y
    enddo
  enddo

end subroutine assign_ax_20



subroutine assign_ax_11(x, y)

  real(ark), intent(inout) :: x(:)
  type(adf_realq), intent(in) :: y(:)

  integer(ik) :: i

  if (size(x)/=size(y)) then
    write(out, '(/a)') 'assign_ax_11 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x)
    x(i) = y(i)
  enddo

end subroutine assign_ax_11



subroutine assign_ax_22(x, y)

  real(ark), intent(inout) :: x(:,:)
  type(adf_realq), intent(in) :: y(:,:)

  integer(ik) :: i, j

  if (.not.all(shape(x)==shape(y))) then
    write(out, '(/a)') 'assign_ax_22 error: shapes of x and y do not conform'
    stop
  endif

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      x(j,i) = y(j,i)
    enddo
  enddo

end subroutine assign_ax_22
