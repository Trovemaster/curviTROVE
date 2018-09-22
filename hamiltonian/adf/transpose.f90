function transpose_x(x) result(f)

  type(adf_realq), intent(in) :: x(:,:)
  type(adf_realq) :: f(size(x,dim=2),size(x,dim=1))

  integer(ik) :: i, j

  do i=1, size(x,dim=2)
    do j=1, size(x,dim=1)
      f(i,j) = x(j,i)
    enddo
  enddo

end function transpose_x
