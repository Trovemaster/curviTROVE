
! Computes binomial coefficient n choose k

function bico(n, k)

  integer(ik), intent(in) :: k, n
  real(rk) :: bico

  bico = nint(exp(factln(n)-factln(k)-factln(n-k)))

end function bico


!################################################################################


! Returns ln(n!)

function factln(n)

  integer(ik), intent(in) :: n
  real(rk) :: factln

  if (n<0) then
    write(out, '(/a,1x,i4)') 'factln error: negative factorial =', n
    stop
  endif

  factln = gammln(real(n+1,kind=rk))

end function factln


!################################################################################


! Returns n!

function qfactorial(n)

  integer(ik), intent(in) :: n
  real(ark) :: qfactorial

  real(ark), parameter :: qfac(0:20) = (/1.0_ark, &!
                                         1.0_ark, &!
                                         2.0_ark, &!
                                         6.0_ark, &!
                                         24.0_ark, &!
                                         120.0_ark, &!
                                         720.0_ark, &!
                                         5040.0_ark, &!
                                         40320.0_ark, &!
                                         362880.0_ark, &!
                                         3628800.0_ark, &!
                                         39916800.0_ark, &!
                                         479001600.0_ark, &!
                                         6227020800.0_ark, &!
                                         87178291200.0_ark, &!
                                         1307674368000.0_ark, &!
                                         20922789888000.0_ark, &!
                                         355687428096000.0_ark, &!
                                         6402373705728000.0_ark, &!
                                         121645100408832000.0_ark, &!
                                         2432902008176640000.0_ark /)

  if (n<0) then
    write(out, '(/a,1x,i4)') 'factorial error: negative factorial =', n
    stop
  endif

  if (n<=20) then
    qfactorial = qfac(n)
  else
    qfactorial = real( exp(factln(n)), kind=ark )
  endif

end function qfactorial


!################################################################################


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


!################################################################################


subroutine tokenize_str(str0, n, word)

  character(*), intent(in) :: str0
  integer(ik), intent(out) :: n
  character(cl), intent(out) :: word(:)

  integer(ik) :: pos1
  character(len(str0)) :: str, str_

  word = ''

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


!################################################################################


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


!################################################################################


SUBROUTINE sort2(n,arr,brr)
  INTEGER(ik) n,M,NSTACK
  REAL(rk) arr(n)
  INTEGER(ik) brr(n)
  PARAMETER (M=7,NSTACK=50)
  !Sorts an array arr(1:n) into ascending order using Quicksort, while making the corresponding
  !rearrangement of the array brr(1:n).
  INTEGER(ik) i,ir,j,jstack,k,l,istack(NSTACK)
  REAL(rk) a,atemp
  INTEGER(ik) b,btemp
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
    do j=l+1,ir
      a=arr(j)
      b=brr(j)
      do i=j-1,l,-1
        if(arr(i).le.a)goto 2
        arr(i+1)=arr(i)
        brr(i+1)=brr(i)
      enddo
      i=l-1
2     arr(i+1)=a
      brr(i+1)=b
    enddo
    if(jstack.eq.0)return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
  else
    k=(l+ir)/2
    atemp=arr(k)
    arr(k)=arr(l+1)
    arr(l+1)=atemp
    btemp=brr(k)
    brr(k)=brr(l+1)
    brr(l+1)=btemp
    if(arr(l).gt.arr(ir))then
      atemp=arr(l)
      arr(l)=arr(ir)
      arr(ir)=atemp
      btemp=brr(l)
      brr(l)=brr(ir)
      brr(ir)=btemp
    endif
    if(arr(l+1).gt.arr(ir))then
      atemp=arr(l+1)
      arr(l+1)=arr(ir)
      arr(ir)=atemp
      btemp=brr(l+1)
      brr(l+1)=brr(ir)
      brr(ir)=btemp
    endif
    if(arr(l).gt.arr(l+1))then
      atemp=arr(l)
      arr(l)=arr(l+1)
      arr(l+1)=atemp
      btemp=brr(l)
      brr(l)=brr(l+1)
      brr(l+1)=btemp
    endif
    i=l+1
    j=ir
    a=arr(l+1)
    b=brr(l+1)
3   continue
    i=i+1
    if(arr(i).lt.a)goto 3
4   continue
    j=j-1
    if(arr(j).gt.a)goto 4
    if(j.lt.i)goto 5
    atemp=arr(i)
    arr(i)=arr(j)
    arr(j)=atemp
    btemp=brr(i)
    brr(i)=brr(j)
    brr(j)=btemp
    goto 3
5   arr(l+1)=arr(j)
    arr(j)=a
    brr(l+1)=brr(j)
    brr(j)=b
    jstack=jstack+2
    if(jstack.gt.NSTACK) then
      write(out, '(/a,1x,i3,1x,a)') 'sort2 error: nstack=', nstack, 'is too small'
     stop
    endif
    if(ir-i+1.ge.j-l)then
      istack(jstack)=ir
      istack(jstack-1)=i
      ir=j-1
    else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    endif
  endif
  goto 1
end subroutine sort2


!################################################################################


! Computes x/sin(x)

function x_div_sinx(x, thresh) result(f)

  real(ark), intent(in) :: x
  real(ark), intent(in), optional :: thresh
  real(ark) :: f

  integer(ik) :: i
  real(ark) :: thresh_

  ! x/sin(x) Taylor series expansion coefficients at x=0
  real(ark), parameter :: coefs(0:10) = (/ 1.0_ark, &!
                                           0.16666666666666666666666666666667_ark, &!
                                           0.019444444444444444444444444444444_ark, &!
                                           0.0020502645502645502645502645502646_ark, &!
                                           0.00020998677248677248677248677248677_ark, &!
                                           0.000021336045641601197156752712308268_ark, &!
                                           2.1633474427786597098766410935723_ark*10.0_ark**(-6), &!
                                           2.1923271344567640863937160233457_ark*10.0_ark**(-7), &!
                                           2.2213930853920414559485002477014_ark*10.0_ark**(-8), &!
                                           2.2507674795567867297317145265855_ark*10.0_ark**(-9), &!
                                           2.2805107707218211704626350426388_ark*10.0_ark**(-10) /)
  thresh_ = 0.01_ark
  if (present(thresh)) thresh_ = abs(thresh)

  if (abs(x)>thresh) then
    !if (abs(abs(x)-real(pi,ark))<=epsilon()) then
    !endif
    f = x/sin(x)
  else
    f = 0.0
    do i=0, 10
      f = f + coefs(i)*x**(2*i)
    enddo
  endif


end function x_div_sinx


!################################################################################


subroutine polint(xa,ya,n,x,y,dy)
  integer(ik) :: n, NMAX
  real(ark) :: dy,x,y,xa(n),ya(n)
  parameter (NMAX=1000)
  integer(ik) :: i,m,ns
  real(ark) :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=abs(x-xa(1))
  do i=1,n
    dift=abs(x-xa(i))
    if (dift.lt.dif) then
      ns=i
      dif=dift
    endif
    c(i)=ya(i)
    d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
    do i=1,n-m
      ho=xa(i)-x
      hp=xa(i+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
      if(abs(den)<=epsilon(1.0_ark)) then
        write(out, '(/a)') 'failure in polint'
        stop
      endif
      den=w/den
      d(i)=hp*den
      c(i)=ho*den
    enddo
    if (2*ns.lt.n-m)then
      dy=c(ns+1)
    else
      dy=d(ns)
      ns=ns-1
    endif
    y=y+dy
  enddo
end subroutine polint


!################################################################################


subroutine polint_fit(dimen, npoints, x, y, x0, nterms, terms, y0, rms)

  integer(ik), intent(in) :: dimen, npoints, nterms, terms(dimen,nterms)
  real(ark), intent(in) :: x(dimen,npoints), y(npoints), x0(dimen)
  real(ark), intent(out) :: y0, rms

  integer(ik) :: iterm, jterm, ipoint, nsv, lwork, info
  real(ark) :: deriv(npoints,nterms), coefs(nterms), mat, vec, f, df
  double precision :: matd(nterms,nterms), vecd(nterms,1), sv(nterms), work(nterms*nterms*32), tol

  ! compute derivatives

  deriv = 0.0
  do iterm=1, nterms
    coefs = 0
    coefs(iterm) = 1.0_ark
    do ipoint=1, npoints
      deriv(ipoint,iterm) = fitfunc(dimen, nterms, terms, coefs, x(1:dimen,ipoint), x0)
    enddo
  enddo

  ! compute system of linear equations

  do iterm=1, nterms
    do jterm=1, nterms
      mat = sum( deriv(1:npoints,iterm) * deriv(1:npoints,jterm) )
      matd(jterm,iterm) = dble( mat )
    enddo
    vec = sum( deriv(1:npoints,iterm) * y(1:npoints) )
    vecd(iterm,1)=dble( vec )
  enddo

  ! Lapack-solve system of linear equations

  tol = -1.0d-12
  lwork = size(work)
  call dgelss( nterms, nterms, 1, matd, nterms, vecd, nterms, sv, tol, nsv, work, lwork, info )
  if (info/=0) then
    write(out, '(/a)') 'polint_fit error: SVD failed'
    stop
  endif

  ! coefficients

  do iterm=1, nterms
    coefs(iterm) = real(vecd(iterm,1),kind=ark)
  enddo

  ! rms

  rms = 0.0
  do ipoint=1, npoints
    f = fitfunc(dimen, nterms, terms, coefs, x(1:dimen,ipoint), x0)
    df = y(ipoint) - f
    rms = rms + df**2
  enddo
  rms = sqrt(rms/real(npoints,kind=ark))

  ! y0

  y0 = fitfunc(dimen, nterms, terms, coefs, x0(1:dimen), x0)

  contains

  function fitfunc(dimen, nterms, terms, coefs, x, x0) result(f)
    integer(ik), intent(in) :: dimen, nterms, terms(dimen,nterms)
    real(ark), intent(in) :: coefs(nterms), x(dimen), x0(dimen)
    real(ark) :: f
    integer(ik) :: iterm
    f = 0.0
    do iterm=1, nterms
      f = f + coefs(iterm) * product( (x(1:dimen)-x0(1:dimen))**terms(1:dimen,iterm) )
    enddo
  end function fitfunc

end subroutine polint_fit


!################################################################################

   subroutine diag_ulen_ark(n,a,d,ve)
      !
      integer(ik)  ::  p,q,p1,q1,irot,i,n2,kp,kp1,ip1,j,iq1,kq1,ii,n
      real(ark)     ::  a(n,n),d(n),ve(n,n)
      real(ark)     ::  sm,tresh,g,h,s,c,t,tau,theta,ff1,ff2,ff3,ff4
      real(ark)     ::  err
      real(ark),allocatable  ::  b(:),z(:)
      !
      err = small_
      !
      allocate(b(n+10),z(n+10))


 101  format(5e14.5)
      do 10 p=1,n
      do 10 q=1,n
      ve(p,q)=0.0_ark
      if(p.eq.q) ve(p,q)=1.0_ark
  10  continue
      do 99 p=1,n
      z(p)=0.0_ark
      d(p)=a(p,p)
      b(p)=d(p)
 99   continue
      irot=0
      do 50 i=1,50
      sm=0.0_ark
      n2=n-1
      do 30 p=1,n2
      kp=p+1
      do 30 q=kp,n
      sm=sm+abs(a(p,q))
  30  continue
      if(sm.le.err) goto 50
      tresh=0.0_ark
      if(i-4) 3,4,4
  3   tresh=0.2_ark*sm/(n*n)
  4     do 33 p1=1,n2
        kp1=p1+1
      do 33 q1=kp1,n
      g=100*abs(a(p1,q1))
      ff1=abs(d(p1)+g)
      ff2=abs(d(p1))
      ff3=abs(d(q1)+g)
      ff4=abs(d(q1))
      if(i.gt.4.and.(ff1.eq.ff2).and.ff3.eq.ff4) goto 7
      ff1=abs(a(p1,q1))
        if(ff1.le.tresh) goto 33
      h=d(q1)-d(p1)
      ff1=abs(h)+g
      ff2=abs(h)
      if(ff1.ne.ff2) goto 13
      t=a(p1,q1)/h
      goto 6
  13    theta=0.5_ark*h/a(p1,q1)
        t=1._ark/(abs(theta)+sqrt(1._ark+theta*theta))
      if(theta) 5,6,6
  5     t=-t
  6     c=1._ark/sqrt(1._ark+t*t)
        s=t*c
      tau=s/(1._ark+c)
      h=t*a(p1,q1)
      z(p1)=z(p1)-h
      z(q1)=z(q1)+h
      d(p1)=d(p1)-h
      d(q1)=d(q1)+h
      a(p1,q1)=0.0_ark
      ip1=p1-1
        do 20 j=1,ip1
        g=a(j,p1)
        h=a(j,q1)
        a(j,p1)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  20      continue
        iq1=q1-1
        do 21 j=kp1,iq1
        g=a(p1,j)
        h=a(j,q1)
        a(p1,j)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  21    continue
        kq1=q1+1
        do 26 j=kq1,n
        g=a(p1,j)
        h=a(q1,j)
        a(p1,j)=g-s*(h+g*tau)
        a(q1,j)=h+s*(g-h*tau)
  26    continue
          do 29 j=1,n
        g=ve(j,p1)
        h=ve(j,q1)
        ve(j,p1)=g-s*(h+g*tau)
        ve(j,q1)=h+s*(g-h*tau)
  29      continue
        irot=irot+1
  7     a(p1,q1)=0.0_ark
  33    continue
        do 44 ii=1,n
      d(ii)=b(ii)+z(ii)
      b(ii)=d(ii)
      z(ii)=0.0d0
  44  continue
  50  continue



  deallocate(b,z)

  end subroutine diag_ulen_ark
