module plasma
 
!dec$ define plasma_ = 0

!
!  Simplistic type-agnostic PLASMA interface
!
  use accuracy
  use timer
  implicit none
  private verbose

  interface plasma_sytrdx
    module procedure plasma_diag_dsytrdx
  end interface ! plasma_sytrdx

  integer,parameter:: verbose = 4
  !
  contains


  subroutine plasma_diag_dsytrdx(n,h,e,nprocs)
    !
    !dec$ if (plasma_ > 0) 
      INCLUDE "plasmaf.h"
      integer(ik), parameter :: VEC = PlasmaVec
      integer(ik), parameter :: UPLO = PlasmaLower
      EXTERNAL PLASMA_DSYTRDX
      INTEGER PLASMA_DSYTRDX
    !dec$ end if
    !
    integer         , intent(in)    :: n
    6double precision, intent(inout) :: h(n,n)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(n)    ! Out: Eigenvalues
    integer         , intent(in)    :: nprocs
    integer          :: info,nprocs_,VBLKSIZ
    integer          :: nh,alloc,k,OMP_GET_NUM_THREADS,tid
    integer          :: msize,nb,ib,lda,ldq,corea,coreb,corec,i,j
    integer(hik)     :: hsize
    double precision,allocatable :: a(:,:),w(:),d(:)
    double precision :: mat_elem,dot_product
    !
    if (verbose>=4) call TimerStart('plasma_sytrdx: diagonalization')
    !
    !dec$ if (plasma_ == 0) 
       write(out,"('Plasma is not activated, in plasma.f90 please set plasma_ to 1')")
    !dec$ end if
    !
    NB = 320
    if (n>150000) NB = 384
    IB = 64
    LDA = N
    LDQ = N
    !
    !$omp parallel private(tid)
      if (tid==0) then
         nprocs_ = OMP_GET_NUM_THREADS()
      endif
    !$omp end parallel
    !
    if (verbose>=4) write(out,"('nthreads = ',i7)") nprocs_
    !
    COREA = nprocs_
    if (n>=64000)  corea = max( min(int( ( 632.0+6.0*N/1000.0 )/17.0,ik ),nprocs_),1 )
    if (n>200000)  corea = min( 96,corea )
    !
    !???????????????
    !!!!COREA = nprocs_/2
    !
    COREB = nprocs_
    COREC = nprocs_
    VBLKSIZ = 128
    !
    if (verbose>=4) then
      !
      write(out,"('COREA = ',i7)") COREA
      call print_openmp
      !
    endif
    !
    call OMP_SET_NUM_THREADS(1) 
    call mkl_set_num_threads(1)
    !
    ! check openmp environment variables
    !
    if (verbose>=4) then
      call print_openmp
    endif
    !
    if (verbose>=4) write(out,"(' allocation of a matrix of ',i,'x',i)") N, N
    !
    hsize = int(n,hik)*int(n,hik)
    !
    allocate (a(LDA,N), stat=alloc)
    call ArrayStart('plasma_sytrdx-a',alloc,1,kind(a),hsize)
    allocate (d(N), stat=alloc)
    call ArrayStart('plasma_sytrdx-d',alloc,size(d),kind(d))
    allocate ( w(N), stat=alloc)
    call ArrayStart('plasma_sytrdx-w',alloc,size(w),kind(w))
    !
    call estimate_plasma_memory(N,NB,IB,VBLKSIZ,hsize)
    call ArrayStart('plasma_sytrdx-all-plasma',alloc,1,kind(w),hsize)
    !
    if (verbose>=4) call MemoryReport
    !
    if (verbose>=3) write(out,"(/   'Diagonalization with PLASMA_dsytrdX...')")
    !
    info = -1
    !
    e = 0 
    !
    !dec$ if (plasma_ > 0)
        !
        INFO = PLASMA_DSYTRDX( VEC, UPLO, N, h, LDA, D, W, E, a, LDQ, COREA, COREB, COREC, NB, IB )
        !
        IF (INFO .NE. PLASMA_SUCCESS) THEN
           PRINT*,"PLASMA_DSYTRDX RETURNS ", INFO, " ABORTING!"
           STOP
        END IF
    !dec$ end if
    !
    call OMP_SET_NUM_THREADS(nprocs_) 
    call mkl_set_num_threads(nprocs_)
    !
    if (verbose>=3) write(out,"('   ...done!')") 
    !
    !if (verbose>=4) write(out,"('nthreads = ',i7)") OMP_GET_NUM_THREADS()
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,n
       h(:,k) = a(:,k)
    enddo
    !$omp end parallel do
    !
    !if (verbose>=5) then 
    !  write(out,"(//'olution:')") 
    !  do i=1,n
    !     write(out,"(i8,f19.8)") i,e(i)
    !     do j=1,n
    !       write(out,"(2i8,f19.8)") i,j,a(j,i)
    !     enddo
    !  enddo
    !  !
    !  call MemoryReport
    !  !
    !endif 
    !
    !
    if (verbose>=5) then 
      write(out,"(/'Eigenvalues:')")
      !
      !omp parallel do private(i,j,mat_elem) schedule(dynamic)
      do i=1,n
         !
         write(out,"(i8,f19.8)") i,e(i)-e(1)
         !
         !do j=i,n
         !  !
         !  mat_elem = dot_product(h(:,i),h(:,j))
         !  !
         !  write(out,"(2i8,g19.8)") i,j,mat_elem
         !  !
         !enddo
      enddo
      !omp end parallel do
      !
      call MemoryReport
      !
    endif 
    !
    deallocate(a,d,w)
    call ArrayStop('plasma_sytrdx-a')
    call ArrayStop('plasma_sytrdx-d')
    call ArrayStop('plasma_sytrdx-w')
    call ArrayStop('plasma_sytrdx-all-plasma')
    !
    if (verbose>=4) then
      call print_openmp
    endif
    !
    if (verbose>=4) call TimerStop('plasma_sytrdx: diagonalization')
    !
  end subroutine plasma_diag_dsytrdx

  subroutine print_openmp
    !
    implicit none
    integer(ik)      :: TID,PROCS,NTHREADS,MAXT
    logical          :: INPAR,DYNAMIC,NESTED,OMP_IN_PARALLEL,OMP_GET_DYNAMIC,OMP_GET_NESTED
    integer(ik)      :: OMP_GET_THREAD_NUM,OMP_GET_NUM_PROCS,OMP_GET_MAX_THREADS,OMP_GET_NUM_THREADS
    !
    !     Start parallel region
    !$OMP PARALLEL PRIVATE(NTHREADS, TID)
    !
    !  Obtain thread number
    TID = OMP_GET_THREAD_NUM()
    !
    !     Only master thread does this
    IF (TID .EQ. 0) THEN
        !
        PRINT *, 'Thread',tid,'getting environment information'
        !     Get environment information
        PROCS = OMP_GET_NUM_PROCS() 
        NTHREADS = OMP_GET_NUM_THREADS()
        MAXT = OMP_GET_MAX_THREADS()
        INPAR = OMP_IN_PARALLEL()
        DYNAMIC = OMP_GET_DYNAMIC()
        NESTED = OMP_GET_NESTED()
        !     Print environment information
        !
        PRINT *, 'Number of processors = ', PROCS
        PRINT *, 'Number of threads = ', NTHREADS
        PRINT *, 'Max threads = ', MAXT
        PRINT *, 'In parallel? = ', INPAR
        PRINT *, 'Dynamic threads enabled? = ', DYNAMIC
        PRINT *, 'Nested parallelism supported? = ', NESTED
        !
        call OMP_SET_NUM_THREADS(1) 
        !
      END IF
      !     Done
    !$OMP END PARALLEL
    !
  end subroutine print_openmp

  subroutine estimate_plasma_memory(N,NB,IB,VBLKSIZ,hsize)
   !
   integer,intent(in)  :: n,NB,ib,vblksiz
   integer(hik),intent(out) :: hsize
   integer :: m,bcnt
   real(rk) :: total
     !
     total = 0
     total = total + n*n
     total = total + (n+real(ib)/real(nb)+1.0)*(n+real(ib)/real(nb)+1.0)
     total = total + real((2*nb+6),rk)*real(n,rk)
     if (verbose>=4) print*,'Memory required for computing plasma eigenvalues'
     if (verbose>=4) print*,'only is ',total+6*n,' double words.'
     m = n/vblksiz + 1
     bcnt = m * (real(n)/real(nb) + 1) - (m*(m-1)/2)*real(vblksiz)/real(nb)
     total = total + vblksiz*bcnt*vblksiz
     total = total + bcnt*vblksiz
     total = total + (nb+vblksiz-1)*bcnt*vblksiz
     !
     if (verbose>=4) print*,'Memory required for both plasma eigenvalues and eigenvectors'
     if (verbose>=4) print*,'is ',total,' double words'
     !
     hsize = int(total,hik)
     !
   end subroutine estimate_plasma_memory


end module plasma
