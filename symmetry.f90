!
module symmetry
  use accuracy

  implicit none
  public SymmetryInitialize,sym,max_irreps


  type  ScIIT
     integer(ik)          :: Noper     ! Number of operations in the CII operator
     integer(ik)          :: Nzeta     ! Number of symmetric elements taking into account the degeneracies 
     real(ark),pointer    :: ioper(:)  ! the operation number in the MS group
     integer(ik),pointer  :: coeff(:)  ! coefficients of the CII operator
     integer(ik),pointer  :: izeta(:)  ! symmetry indentification as a eigenvalues of the CII operator
  end type ScIIT


  type  SymmetryT
     character(len=cl)    :: group = 'C' ! The symmetry group 
     integer(ik)          :: Nrepresen = 1  ! Number of irreduc. represent.
     integer(ik)          :: Noper     = 1  ! Number of operations
     integer(ik)          :: Nclasses  = 1  ! Number of classes
     integer(ik),pointer  :: Nelements(:)   ! Number of elements in a class
     real(rk),pointer     :: characters(:,:)! Character table
     type(SrepresT),pointer :: irr(:,:)     ! irreducible representaion 
     integer(ik),pointer  :: degen(:)       ! degeneracy
     character(len=3),pointer  :: label(:)  ! The symmetry label 
     integer(ik)          :: Maxdegen  = 1  ! Maximal degeneracy order
     integer(ik),pointer  :: igenerator(:)  ! address of the class generator in the sym%Ngroup list
     type(ScIIT)          :: CII            ! the elements of the CII operator 
     real(ark),pointer    :: euler(:,:)     ! rotational angles equivalent to the group operations
     integer(ik)          :: class_size_max = 8 ! current maximal class size 
     integer(ik)          :: N  = 1         ! The group order, currently desgined for Dnh where N is odd 
     integer(ik),pointer  :: lquant(:)      ! Store the value of the (vibrational) angular momentum 
     !
  end type SymmetryT

  type  SrepresT
     real(ark),pointer  :: repres(:,:)      ! matrix representation of the group 
  end type SrepresT


  type(SymmetryT) , save  :: sym
  integer(ik),parameter   :: max_irreps=100
  integer(ik),parameter   :: verbose_ = 3

contains 


  subroutine SymmetryInitialize(sym_group)
  character(len=cl),intent(in) :: sym_group
  integer(ik):: alloc,iclass,gamma,ioper,ielem,irepr,Nrot,irep,k,irot,N_Cn,ioper_
  real(ark)  :: a,b,e,o,p2,p3,p4,p23,p43,phi,phi_n,factor
  character(len=4) :: Kchar
  !   
  sym%group=sym_group
  !
  select case(trim(sym_group))

  case("C(M)","C")

    sym%Nrepresen=1
    sym%Noper=1
    sym%Nclasses=1
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters(1,1)=1.0_rk
    sym%degen=(/1/)
    sym%Nelements=(/1/)
    sym%label=(/'A'/)

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 


  case("CS(M)","CS")

    sym%Nrepresen=2
    sym%Noper=2
    sym%Nclasses=2
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                 (/ 1, 1,  & 
                                    1,-1  /),(/2,2/))
    sym%degen=(/1,1/)
    sym%Nelements=(/1,1/)
    sym%label=(/'A''','A"'/)

    call irr_allocation

  case("C2V(M)","C2V")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! A1
                                     1, 1,-1,-1, &   ! A2
                                     1,-1,-1, 1, &   ! B1
                                     1,-1, 1,-1 /),(/4,4/)) ! B2
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A1','A2','B1','B2'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o /)
    sym%euler( 3,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/p2,pi,p3/)
    !
    call irr_allocation


  case("C2H(M)","C2H")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! Ag
                                     1, 1,-1,-1, &   ! Au
                                     1,-1,-1, 1, &   ! Bg
                                     1,-1, 1,-1 /),(/4,4/)) ! Bu
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'Ag','Au','Bg','Bu'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 2,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    !
    ! 1324 and 1234 are working
    !
    call irr_allocation
    !
  case("CS(EM)")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2a  E* sab i
                                  (/ 1,  1,  1,  1, &   ! Ag 
                                     1,  1, -1, -1, &   ! Au 
                                     1, -1,  1, -1, &   ! Bg
                                     1, -1, -1,  1/),(/4 ,4/)) ! Bu
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'Ag','Au','Bg','Bu'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 4,:) = (/o,pi,o/)
    !
    call irr_allocation
    !
  case("G4(M)","G4")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &       ! E (12)(34)  E*  (12)(34)*
                                  (/ 1, 1, 1, 1, &   ! A+
                                     1, 1,-1,-1, &   ! A-
                                     1,-1,-1, 1, &   ! B+
                                     1,-1, 1,-1 /),(/4,4/)) ! B-
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A+','A-','B-','B+'/)
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p2,pi,p3/)
    sym%euler( 3,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    !
    call irr_allocation
    !
  case("G4(EM)")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E   a   b  ab   E' E'a E'b E'ab 
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ags
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! Aus
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! Bgs
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! Bus
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! Agd 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Aud
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! Bgd
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! Bud
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ags','Aus','Bgs','Bus','Agd','Aud','Bgd','Bud'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 3,:) = (/p2,pi,p3/)
    sym%euler( 2,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/pi,o,o/)
    sym%euler( 5,:) = (/pi,o,o/)
    sym%euler( 7,:) = (/o,pi,o/)
    sym%euler( 6,:) = (/p2,pi,p3/)
    sym%euler( 8,:) = 0


    !sym%euler( 1,:) = 0
    !sym%euler( 3,:) = (/p2,pi,p3/)
    !sym%euler( 4,:) = (/o,pi,o/)
    !sym%euler( 2,:) = (/pi,o,o/)
    !sym%euler( 7,:) = (/pi,o,o/)
    !sym%euler( 5,:) = (/o,pi,o/)
    !sym%euler( 6,:) = (/p2,pi,p3/)
    !sym%euler( 8,:) = 0


    !
    call irr_allocation

  case("D2H(M)")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2a  C2b C2c E* sab sac sbc i
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ag 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Au 
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! B1g
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! B1u
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! B2g 
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! B2u
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! B3g
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! B3u
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ag','Au','B1g','B1u','B2g','B2u','B3g','B3u'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/o,pi,o/) ! Ry
    sym%euler( 4,:) = (/p2,pi,p3/) ! Rx
    !
    sym%euler( 8,:) = (/p2,pi,p3/)
    sym%euler( 7,:) = (/o,pi,o/)
    sym%euler( 6,:) = (/pi,o,o/)
    sym%euler( 5,:) = 0
    !
    call irr_allocation
    !
  case("D2H")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E  C2z  C2y C2x i  sxy sxz syz  
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ag 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Au 
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! B1g
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! B1u
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! B2g 
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! B2u
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! B3g
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! B3u
                                     

    sym%characters = transpose(sym%characters)
    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ag','Au','B1g','B1u','B2g','B2u','B3g','B3u'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o/) ! Rz
    sym%euler( 3,:) = (/o,pi,o/) ! Ry
    sym%euler( 4,:) = (/p2,pi,p3/) ! Rx
    !
    sym%euler( 7,:) = (/p2,pi,p3/)
    sym%euler( 8,:) = (/o,pi,o/)
    sym%euler( 5,:) = (/pi,o,o/)
    sym%euler( 6,:) = 0
    !
    call irr_allocation    
    !
  case("C3V(M)","C3V")
    !
    sym%Nrepresen=3
    sym%Noper=6
    sym%Nclasses=3
    sym%CII%Noper = 0
    !
    call simple_arrays_allocation


    sym%characters= reshape( &      !A1 A2 E
                                  (/ 1, 1, 2, &  
                                     1, 1,-1, &  
                                     1,-1, 0 /),(/3,3/)) 
    sym%degen=(/1,1,2/)
    sym%Nelements=(/1,2,3/)
    sym%label=(/'A1','A2','E '/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 2.0_ark/3.0_ark*pi
    p4 = 4.0_ark/3.0_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p4,o ,o /)
    sym%euler( 3,:) = (/p3,o ,o /)
    sym%euler( 4,:) = (/o ,pi,o /)
    sym%euler( 5,:) = (/p3,pi,o /)
    sym%euler( 6,:) = (/p4,pi,o /)
    !
    sym%lquant(3) = 1
    !
    call irr_allocation
       !
       sym%irr(3,1)%repres = reshape((/1.0_ark,               0.0_ark, &
                                       0.0_ark,               1.0_ark/),(/2,2/))
       !
       sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                        0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
       !
       sym%irr(3,3)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
                                       -0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
       !
       sym%irr(3,4)%repres = reshape((/ 1.0_ark,               0.0_ark, &
                                        0.0_ark,              -1.0_ark/),(/2,2/))
       !
       sym%irr(3,5)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
                                        0.5_ark*sqrt(3.0_ark), 0.5_ark            /),(/2,2/))
       !
       sym%irr(3,6)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                       -0.5_ark*sqrt(3.0_ark), 0.5_ark             /),(/2,2/))
       !
  case("D3H(M)","D3H")

    sym%Nrepresen=6
    sym%Noper=12
    sym%Nclasses=6
    sym%CII%Noper = 0

    call simple_arrays_allocation


    sym%characters= reshape( &      !A1 A2 E  A1 A2 E
                                  (/ 1, 1, 2, 1, 1, 2, &  
                                     1, 1,-1, 1, 1,-1, &  
                                     1,-1, 0, 1,-1, 0, &          
                                     1, 1, 2,-1,-1,-2, &  
                                     1, 1,-1,-1,-1, 1, &  
                                     1,-1, 0,-1, 1, 0 /),(/6,6/)) 
    sym%degen=(/1,1,2,1,1,2/)
    sym%Nelements=(/1,2,3,1,2,3/)
    sym%label=(/'A1''','A2''','E'' ','A1"','A2"','E" '/)
    !
    o   = 0.0_ark
    p2  = 0.5_ark*pi
    p3  = pi/3.0_ark
    p23 = 2.0_ark/3.0_ark*pi
    p43 = 4.0_ark/3.0_ark*pi
    !
    sym%euler( 1,:) = 0               ! (E)
    sym%euler( 2,:) = (/p23,o ,o /)   ! (132)
    sym%euler( 3,:) = (/p43,o ,o /)   ! (123)
    sym%euler( 4,:) = (/ pi,pi,o /)   ! (23)
    sym%euler( 5,:) = (/ p3,pi,o /)   ! (12)
    sym%euler( 6,:) = (/-p3,pi,o /)   ! (13)
    sym%euler( 7,:) = (/ pi, o,o /)   ! (E)*
    sym%euler( 8,:) = (/-p3, o,o /)   ! (132)*
    sym%euler( 9,:) = (/ p3, o,o /)   ! (123)*
    sym%euler(10,:) = (/  o,pi,o /)   ! (23)*
    sym%euler(11,:) = (/p43,pi,o /)   ! (12)*
    sym%euler(12,:) = (/p23,pi,o /)   ! (13)*
    !
    call irr_allocation
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    sym%irr(3,1)%repres = reshape((/ e, o,  &
                                     o, e/),(/2,2/))
    !
    sym%irr(3,2)%repres = reshape((/-a, b,  &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(3,3)%repres = reshape((/-a,-b,  &
                                     b,-a/),(/2,2/))
    !
    sym%irr(3,4)%repres = reshape((/ e, o,  &
                                     o,-e/),(/2,2/))
    !
    sym%irr(3,5)%repres = reshape((/-a,-b,  &
                                    -b, a /),(/2,2/))
    !
    sym%irr(3,6)%repres = reshape((/-a, b,  &
                                     b, a  /),(/2,2/))
    !
    !
    sym%irr(3,7)%repres = reshape((/ e, o, &
                                     o, e/),(/2,2/))
    !
    sym%irr(3,8)%repres = reshape((/-a, b, &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(3,9)%repres = reshape((/-a,-b, &
                                     b,-a/),(/2,2/))
    !
    sym%irr(3,10)%repres= reshape((/ e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(3,11)%repres= reshape((/-a,-b, &
                                    -b, a /),(/2,2/))
    !
    sym%irr(3,12)%repres= reshape((/-a, b, &
                                     b, a  /),(/2,2/))
    !
    sym%irr(6,1)%repres = reshape((/e, o,  &
                                    o, e/),(/2,2/))
    !
    sym%irr(6,2)%repres = reshape((/-a, b, &
                                    -b,-a/),(/2,2/))
    !
    sym%irr(6,3)%repres = reshape((/-a,-b, &
                                     b,-a/),(/2,2/))
    !
    sym%irr(6,4)%repres = reshape((/-e, o, &
                                     o, e/),(/2,2/))
    !
    sym%irr(6,5)%repres = reshape((/ a, b, &
                                     b,-a /),(/2,2/))
    !
    sym%irr(6,6)%repres = reshape((/ a,-b, &
                                    -b,-a  /),(/2,2/))
    !
    sym%irr(6,7)%repres = reshape((/-e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(6,8)%repres = reshape((/ a,-b, &
                                     b, a/),(/2,2/))
    !
    sym%irr(6,9)%repres = reshape((/ a, b, &
                                    -b, a/),(/2,2/))
    !
    sym%irr(6,10)%repres= reshape((/ e, o, &
                                     o,-e/),(/2,2/))
    !
    sym%irr(6,11)%repres= reshape((/-a,-b, &
                                    -b, a /),(/2,2/))
    !
    sym%irr(6,12)%repres= reshape((/-a, b, &
                                     b, a  /),(/2,2/))
    !
    sym%lquant(1:2) = 0
    sym%lquant(4:5) = 0
    sym%lquant(3) = 1
    sym%lquant(6) = 1
    !
  case("TD(M)","TD")

    sym%Nrepresen=5
    sym%Noper=24
    sym%Nclasses=5
    sym%CII%Noper = 6

    call simple_arrays_allocation

    sym%CII%coeff = (/1.0_ark,1.0_ark,4.0_ark,10.0_ark,1.0_ark,1.0_ark/)
    sym%CII%ioper = (/19,20,21,22,23,24/)



    sym%characters= reshape( &      !A1 A2 E  F1 F2
                                      (/ 1, 1, 2, 3, 3, &  
                                         1, 1,-1, 0, 0, &  
                                         1, 1, 2,-1,-1, &
                                         1,-1, 0, 1,-1, &
                                         1,-1, 0,-1, 1  /),(/5,5/)) 
    sym%degen=(/1,1,2,3,3/)
    sym%Nelements=(/1,8,3,6,6/)
    sym%label=(/'A1','A2','E ','F1','F2'/)

    sym%CII%Nzeta = sum(sym%degen(:))
    allocate(sym%CII%izeta(sym%CII%Nzeta),stat=alloc)
    !
    sym%CII%izeta = (/18,-18,12,-12,-14,-8,4,14,-4,8/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/p3,p2,o /)
    sym%euler( 3,:) = (/pi,p2,p3/)
    sym%euler( 4,:) = (/ o,p2,p3/)
    sym%euler( 5,:) = (/p3,p2,pi/)
    sym%euler( 6,:) = (/p2,p2,o /)
    sym%euler( 7,:) = (/pi,p2,p2/)
    sym%euler( 8,:) = (/ o,p2,p2/)
    sym%euler( 9,:) = (/p2,p2,pi/)
    sym%euler(10,:) = (/p2,pi,p3/)
    sym%euler(11,:) = (/ o,pi,o /)
    sym%euler(12,:) = (/pi, o,o /)
    sym%euler(13,:) = (/ o,p3,o /)
    sym%euler(14,:) = (/ o,p2,o /)
    sym%euler(15,:) = (/p2, o,o /)
    sym%euler(16,:) = (/p3, o,o /)
    sym%euler(17,:) = (/p2,p3,p3/)
    sym%euler(18,:) = (/p2,p2,p3/)
    sym%euler(19,:) = (/p2,p2,p2/)
    sym%euler(20,:) = (/ o,p2,pi/)
    sym%euler(21,:) = (/ o,pi,p2/)
    sym%euler(22,:) = (/pi,pi,p2/)
    sym%euler(23,:) = (/pi,p2,o /)
    sym%euler(24,:) = (/p3,p2,p3/)
    !
    !sym%euler(:,1) = sym%euler(:,1)+pi*0.25_ark
    !sym%euler(:,2) = sym%euler(:,2)
    !sym%euler(:,3) = sym%euler(:,3)-pi*0.25_ark
    !
    call irr_allocation
    !
    sym%irr(3,1)%repres = reshape((/ 1.0_ark,               0.0_ark, &
                                     0.0_ark,               1.0_ark/),(/2,2/))
    !
    !sym%irr(3,2)%repres = reshape((/-0.5_ark,               0.5_ark*sqrt(3.0_ark),&
    !                                -0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
    ! working 
    sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
                                     0.5_ark*sqrt(3.0_ark),-0.5_ark             /),(/2,2/))
    !
    !sym%irr(3,2)%repres = reshape((/-0.5_ark,              -0.5_ark*sqrt(3.0_ark),&
    !                                 0.5_ark*sqrt(3.0_ark),-0.5_ark           /),(/2,2/))
    !
    ! 
    !sym%irr(3,21)%repres = reshape((/-0.5_ark,              0.5_ark*sqrt(3.0_ark),&
    !                                  0.5_ark*sqrt(3.0_ark),0.5_ark             /),(/2,2/))
    !working    
    sym%irr(3,21)%repres = reshape((/-0.5_ark,              0.5_ark*sqrt(3.0_ark),&
                                      0.5_ark*sqrt(3.0_ark),0.5_ark             /),(/2,2/))

    !
    sym%irr(4,1)%repres = reshape((/ 1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark, 1.0_ark, 0.0_ark, &
                                     0.0_ark, 0.0_ark, 1.0_ark/),(/3,3/))
    !
    !sym%irr(4,2)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    !!sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !!                                 1.0_ark, 0.0_ark, 0.0_ark, &
    !!                                 0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
    ! working ->
    !
    sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
                                     1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))

    !sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark, 1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark, &
    !                                 0.0_ark, 1.0_ark, 0.0_ark/),(/3,3/))


    !sym%irr(4,2)%repres = transpose(reshape((/ 0.0_ark, 0.5_ark*sqrt(2.0_ark), 0.5_ark*sqrt(2.0_ark), &
    !                                 0.5_ark*sqrt(2.0_ark), 0.5_ark,-0.5_ark, &
    !                                -0.5_ark*sqrt(2.0_ark), 0.5_ark,-0.5_ark/),(/3,3/)))
    !
    !
    !sym%irr(4,21)%repres = reshape((/-1.0_ark, 0.0_ark, 0.0_ark, &
    !                                  0.0_ark, 0.0_ark,-1.0_ark, &
    !                                  0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))

    ! working->
    sym%irr(4,21)%repres = reshape((/  0.0_ark, 1.0_ark, 0.0_ark, &
                                       1.0_ark, 0.0_ark, 0.0_ark, &
                                       0.0_ark, 0.0_ark,-1.0_ark/),(/3,3/))


    !
    ! roman
    !

     !  sym%irr(4,2)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
     !                                   1.0_ark, 0.0_ark, 0.0_ark, &
     !                                   0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
     !  !
     !  sym%irr(4,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
     !                                    0.0_ark,-1.0_ark, 0.0_ark, &
     !                                   -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))




    !
    !sym%irr(4,21)%repres = reshape((/-1.0_ark, 0.0_ark, 0.0_ark, &
    !                                  0.0_ark, 1.0_ark, 0.0_ark, &
    !                                  0.0_ark, 0.0_ark,-1.0_ark/),(/3,3/))
    
    !sym%irr(4,21)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 -1.0_ark, 0.0_ark, 0.0_ark, &
    !!                                  0.0_ark, 0.0_ark, -1.0_ark/),(/3,3/))
    !


    !
    !!!sym%irr(4,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !!!                                  0.0_ark,-1.0_ark, 0.0_ark, &
    !!!                                 -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,1)%repres = reshape((/ 1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark, 1.0_ark, 0.0_ark, &
                                     0.0_ark, 0.0_ark, 1.0_ark/),(/3,3/))
    ! working
    !sym%irr(5,2)%repres = reshape((/ 0.0_ark,-1.0_ark, 0.0_ark, &
    !                                 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                 1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,2)%repres = reshape((/ 0.0_ark, 0.0_ark, 1.0_ark, &
                                    -1.0_ark, 0.0_ark, 0.0_ark, &
                                     0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))
    ! working
    !sym%irr(5,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
    !                                  0.0_ark, 1.0_ark, 0.0_ark, &
    !                                 -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    sym%irr(5,21)%repres = reshape((/ 0.0_ark, 0.0_ark,-1.0_ark, &
                                      0.0_ark, 1.0_ark, 0.0_ark, &
                                     -1.0_ark, 0.0_ark, 0.0_ark/),(/3,3/))
    !
    !sym%irr(5,2)%repres = reshape((/ 0.0_ark, 0.0_ark,  1.0_ark, &
    !                                -1.0_ark, 0.0_ark, 0.0_ark, &
    !                                 0.0_ark,-1.0_ark, 0.0_ark/),(/3,3/))


    !
    ! (2)  = (123)
    ! (21) = (14)*
    do irepr=3,5
      sym%irr(irepr,3)%repres = matmul(sym%irr(irepr,2)%repres,sym%irr(irepr,2)%repres)   ! (132)
      sym%irr(irepr,13)%repres= matmul(sym%irr(irepr,21)%repres,sym%irr(irepr,2)%repres)  ! (1423)*
      sym%irr(irepr,8)%repres = matmul(sym%irr(irepr,13)%repres,sym%irr(irepr,21)%repres) ! (234)
      sym%irr(irepr,5)%repres = matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,3)%repres)   ! (134)
      sym%irr(irepr,4)%repres = matmul(sym%irr(irepr,5)%repres,sym%irr(irepr,5)%repres)   ! (143)
      sym%irr(irepr,6)%repres = matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,3)%repres)   ! (142)
      sym%irr(irepr,7)%repres = matmul(sym%irr(irepr,6)%repres,sym%irr(irepr,6)%repres)   ! (124)
      sym%irr(irepr,9)%repres = matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,8)%repres)   ! (243)
      sym%irr(irepr,10)%repres= matmul(sym%irr(irepr,7)%repres,sym%irr(irepr,2)%repres)   ! (13)(24)
      sym%irr(irepr,11)%repres= matmul(sym%irr(irepr,8)%repres,sym%irr(irepr,2)%repres)   ! (12)(34)
      sym%irr(irepr,12)%repres= matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,2)%repres)   ! (14)(23)
      sym%irr(irepr,14)%repres= matmul(sym%irr(irepr,3)%repres,sym%irr(irepr,21)%repres)  ! (1324)*
      sym%irr(irepr,15)%repres= matmul(sym%irr(irepr,11)%repres,sym%irr(irepr,21)%repres) ! (1243)*
      sym%irr(irepr,16)%repres= matmul(sym%irr(irepr,10)%repres,sym%irr(irepr,21)%repres) ! (1342)*
      sym%irr(irepr,17)%repres= matmul(sym%irr(irepr,9)%repres,sym%irr(irepr,21)%repres)  ! (1432)*
      sym%irr(irepr,18)%repres= matmul(sym%irr(irepr,2)%repres,sym%irr(irepr,21)%repres)  ! (1234)*
      sym%irr(irepr,19)%repres= matmul(sym%irr(irepr,5)%repres,sym%irr(irepr,21)%repres)  ! (13)* 
      sym%irr(irepr,20)%repres= matmul(sym%irr(irepr,7)%repres,sym%irr(irepr,21)%repres)  ! (12)*
      sym%irr(irepr,22)%repres= matmul(sym%irr(irepr,12)%repres,sym%irr(irepr,21)%repres) ! (23)*
      sym%irr(irepr,23)%repres= matmul(sym%irr(irepr,4)%repres,sym%irr(irepr,21)%repres)  ! (34)*
      sym%irr(irepr,24)%repres= matmul(sym%irr(irepr,6)%repres,sym%irr(irepr,21)%repres)  ! (24)*
    enddo
    sym%lquant(1:2) = 0
    sym%lquant(3) = 1
    sym%lquant(4:5) = 2
    !
  case("DNH(M)","DNH") ! D_infinity_H(M)
    !
    if (mod(sym%N,2)/=1) then
      write(out,"('symmetry: currently Dnh is only for an even ggroup order N ',i8)") sym%N
      stop 'symmetry: illegal order of Dnh group '
    endif
    !
    ! Number of rotations to test for < infinity 
    !
    Nrot = sym%N ! must be >=1
    !
    ! Number of Cn classes 
    N_Cn = sym%N/2
    !
    sym%Noper=2+4*N_Cn+2*Nrot
    sym%Nclasses=4+N_Cn*2
    sym%Nrepresen= 4+N_Cn*2
    sym%CII%Noper = 0
    !
    phi = 2.0_ark*pi/real(Nrot,ark)
    !
    call simple_arrays_allocation
    !
    ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
    !
    sym%characters(:,:) = 0
    !
    ! A1g,A1u,A2g,A2u:
    ! E
    sym%characters(1:4,1) = 1
    ! Cinf
    sym%characters(1:4,1+1:1+N_Cn) = 1
    ! sigmav
    sym%characters(1,1+N_Cn+1) = 1
    sym%characters(2,1+N_Cn+1) =-1
    sym%characters(3,1+N_Cn+1) = 1
    sym%characters(4,1+N_Cn+1) =-1
    ! i
    sym%characters(1,1+N_Cn+2) = 1
    sym%characters(2,1+N_Cn+2) = 1
    sym%characters(3,1+N_Cn+2) =-1
    sym%characters(4,1+N_Cn+2) =-1
    ! Sinf
    sym%characters(1,1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 1
    sym%characters(2,1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 1
    sym%characters(3,1+N_Cn+2+1:1+N_Cn+2+N_Cn) =-1
    sym%characters(4,1+N_Cn+2+1:1+N_Cn+2+N_Cn) =-1
    ! C'inf
    sym%characters(1,1+N_Cn+2+N_Cn+1) = 1
    sym%characters(2,1+N_Cn+2+N_Cn+1) =-1
    sym%characters(3,1+N_Cn+2+N_Cn+1) =-1
    sym%characters(4,1+N_Cn+2+N_Cn+1) = 1
    !
    !sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
    !
    sym%label(1:4)=(/'A1''','A2''','A1"','A2"'/)
    !
    !sym%label=(/'A1''','A2''','E'' ','A1"','A2"','E" '/)
    !
    ! E1' E1" E2' E2" E3' E3" ....
    !
    sym%lquant(1:4) = 0 
    !
    irep = 4
    do k = 1,(sym%Nrepresen-4)/2
      !
      irep = irep + 1
      !
      write(Kchar, '(i4)') K
      !
      sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//''''
      sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'"'
      !
      sym%characters(irep  ,1) = 2.0_ark
      sym%characters(irep+1,1) = 2.0_ark
      !
      sym%characters(irep  ,1+N_Cn+2) = 2.0_ark
      sym%characters(irep+1,1+N_Cn+2) =-2.0_ark
      !
      do irot = 1,N_Cn
        sym%characters(irep  ,1+irot)          = 2.0_ark*cos(phi*irot*k)
        sym%characters(irep+1,1+irot)          = 2.0_ark*cos(phi*irot*k)
        sym%characters(irep  ,1+N_Cn+2+irot)   = 2.0_ark*cos(phi*irot*k)
        sym%characters(irep+1,1+N_Cn+2+irot)   =-2.0_ark*cos(phi*irot*k)
        !
        sym%lquant(irep  ) = k
        sym%lquant(irep+1) = k
        !
      enddo
      !
      irep = irep + 1
      !
    enddo
    !
    sym%degen(:)   = 2
    sym%degen(1:4) = 1
    !
    sym%Nelements(1) = 1
    sym%Nelements(1+ 1:1+ N_Cn) = 2
    sym%Nelements(1+N_Cn+1) = Nrot
    sym%Nelements(1+N_Cn+2) = 1
    sym%Nelements(1+N_Cn+2+1:1+N_Cn+2+N_Cn) = 2
    sym%Nelements(1+N_Cn+2+N_Cn+1) = Nrot
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler(:,:) = 0
    !
    ioper = 1
    do irot = 1,N_Cn
      !
      sym%euler(1+ioper  ,:) = (/o, phi*irot,o/) ! Rz
      sym%euler(1+ioper+1,:) = (/o,-phi*irot,o/) ! Rz
      sym%euler(1+2*N_Cn+Nrot+1+ioper,:)   = (/o,pi+phi*irot,o/) ! Rz
      sym%euler(1+2*N_Cn+Nrot+1+ioper+1,:) = (/o,pi-phi*irot,o/) ! Rz
      !
      ioper = ioper + 2
    enddo
    !
    call irr_allocation
    !
    ! Generate irr-representations
    !
    irep = 4
    do k = 1,(sym%Nrepresen-4)/2
      !
      irep = irep + 1
      !
      ioper = 1
      do ioper = 1,sym%Noper
        !
        factor = 1.0_ark
        !
        if (ioper==1) then ! E 
          !
          !
          sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
          sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
          !
          sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
          sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
          !
          sym%irr(irep+1,ioper)%repres(1,1) = 1.0_ark
          sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
          !
          sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
          sym%irr(irep+1,ioper)%repres(2,2) = 1.0_ark
          !
        elseif (ioper<=1+2*N_Cn) then !  Cinf
          !
          ioper_ =ioper-1 
          irot = (ioper_+1)/2
          !
          phi_n = phi*irot*k
          !
          ! Second oper in a class is with negative phi
          if ( mod(ioper_,2)==0 ) phi_n = -phi_n
          !
          sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
          sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
          !
          sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
          sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
          sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
          sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
          !
        elseif (ioper<=1+2*N_Cn+Nrot) then !  sigmav
          !
          irot = ioper-(1+2*N_Cn)
          !
          phi_n = phi*irot*k
          !
          sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
          sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
          !
          sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
          sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(1,1) = cos(phi_n)
          sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(2,1) =-sin(phi_n)
          sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
          !
        elseif (ioper==1+2*N_Cn+Nrot+1) then ! i
          !
          phi_n = 0
          factor = -1.0_ark
          !
          sym%irr(irep,ioper)%repres(1,1) = 1.0_ark
          sym%irr(irep,ioper)%repres(1,2) = 0.0_ark
          !
          sym%irr(irep,ioper)%repres(2,1) = 0.0_ark
          sym%irr(irep,ioper)%repres(2,2) = 1.0_ark
          !
          sym%irr(irep+1,ioper)%repres(1,1) =-1.0_ark
          sym%irr(irep+1,ioper)%repres(1,2) = 0.0_ark
          !
          sym%irr(irep+1,ioper)%repres(2,1) = 0.0_ark
          sym%irr(irep+1,ioper)%repres(2,2) =-1.0_ark
          !
        elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn) then !  Sinf
          !
          ioper_ = ioper-(1+2*N_Cn+Nrot+1)
          !
          irot = (ioper_+1)/2
          !
          phi_n = phi*irot*k
          !
          ! Second oper in a class is with negative phi
          if (mod(ioper_,2)==0)  phi_n = -phi_n
          !
          sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
          sym%irr(irep,ioper)%repres(1,2) =-sin(phi_n)
          !
          sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
          sym%irr(irep,ioper)%repres(2,2) = cos(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(1,1) =-cos(phi_n)
          sym%irr(irep+1,ioper)%repres(1,2) =-sin(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
          sym%irr(irep+1,ioper)%repres(2,2) =-cos(phi_n)
          !
        elseif (ioper<=1+2*N_Cn+Nrot+1+2*N_Cn+Nrot) then !  C'2
          !
          irot = ioper-(1+2*N_Cn+Nrot+1+2*N_Cn)
          !
          phi_n  = phi*irot*k
          !
          sym%irr(irep,ioper)%repres(1,1) = cos(phi_n)
          sym%irr(irep,ioper)%repres(1,2) = sin(phi_n)
          !
          sym%irr(irep,ioper)%repres(2,1) = sin(phi_n)
          sym%irr(irep,ioper)%repres(2,2) =-cos(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(1,1) =-cos(phi_n)
          sym%irr(irep+1,ioper)%repres(1,2) = sin(phi_n)
          !
          sym%irr(irep+1,ioper)%repres(2,1) = sin(phi_n)
          sym%irr(irep+1,ioper)%repres(2,2) = cos(phi_n)
          !
        else
          !
          stop  'symmetry: illegal ioper'
          !
        endif
        !
      enddo
      !
      irep = irep + 1
      !
    enddo
    !
  case("DINFTYH(M)") ! D_infinity_H(M)

    ! Number of rotations to test for < infinity 
    !
    Nrot = 37 ! must be >=1
    !
    sym%Noper=6+2*Nrot
    sym%Nclasses=6
    sym%Nrepresen= max_irreps  ! must be even and >=4
    sym%CII%Noper = 0
    !
    phi = 2.0_ark*pi/real(Nrot,ark)
    !
    call simple_arrays_allocation
    !
    ! E nrotxCinf nrotxsigmav i  nrotxSinf nrotxC'2
    !
    sym%characters(:,:) = 0
    !
    ! A1g,A1u,A2g,A2u:
    ! E
    sym%characters(1:4,1) = 1
    ! Cinf
    sym%characters(1:4,2) = 1
    ! sigmav
    sym%characters(1,3) = 1
    sym%characters(2,3) = 1
    sym%characters(3,3) =-1
    sym%characters(4,3) =-1
    ! i
    sym%characters(1,4) = 1
    sym%characters(2,4) =-1
    sym%characters(3,4) = 1
    sym%characters(4,4) =-1
    ! Sinf
    sym%characters(1,5) = 1
    sym%characters(2,5) =-1
    sym%characters(3,5) = 1
    sym%characters(4,5) =-1
    ! C'inf
    sym%characters(1,6) = 1
    sym%characters(2,6) =-1
    sym%characters(3,6) =-1
    sym%characters(4,6) = 1
    !
    sym%label(1:4)=(/'A1g','A2g','A1u','A2u'/)
    !
    ! E_1g E_1u E_2g E_2u E_3g E_3u .... 
    !
    irep = 4
    do k = 1,(sym%Nrepresen-4)/2
      !
      irep = irep + 1
      !
      write(Kchar, '(i4)') K
      !
      sym%label(irep  ) = 'E'//trim(adjustl(Kchar))//'g'
      sym%label(irep+1) = 'E'//trim(adjustl(Kchar))//'u'
      !
      sym%characters(irep  ,1) = 2.0_ark
      sym%characters(irep+1,1) = 2.0_ark
      !
      sym%characters(irep  ,4) = 2.0_ark
      sym%characters(irep+1,4) =-2.0_ark
      !
      !do irot = 1,Nrot
      !  sym%characters(irep  ,2)          = 2.0_ark*cos(phi*irot*k)
      !  sym%characters(irep+1,2)          = 2.0_ark*cos(phi*irot*k)
      !  sym%characters(irep  ,5) = 2.0_ark*cos(phi*irot*k)*(-1)**k
      !  sym%characters(irep+1,5) = 2.0_ark*cos(phi*irot*k)*(-1)**(k+1)
      !enddo
      !
      sym%characters(irep  ,2)          = 2.0_ark*cos(phi*k)
      sym%characters(irep+1,2)          = 2.0_ark*cos(phi*k)
      sym%characters(irep  ,5) = 2.0_ark*cos(phi*k)*(-1)**k
      sym%characters(irep+1,5) = 2.0_ark*cos(phi*k)*(-1)**(k+1)
      !
      irep = irep + 1
      !
    enddo
    !
    sym%degen(:)   = 2
    sym%degen(1:4) = 1
    !
    sym%Nelements = (/1,2,Nrot,1,2,Nrot/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler(:,:) = 0
    !
    !do irot = 1,Nrot
      sym%euler(2,:) = (/o, phi,o/) ! Rz
      sym%euler(3,:) = (/o,-phi,o/) ! Rz
      sym%euler(3+Nrot+2,:) = (/o,pi+phi,o/) ! Rz
      sym%euler(3+Nrot+2,:) = (/o,pi-phi,o/) ! Rz
    !enddo
    !
    call irr_allocation        !
  case default

    write(out,"('symmetry: undefined symmetry group ',a)") trim(sym_group)
    stop 'symmetry: undefined symmetry group '

  end select
  !
  if (max_irreps<sym%Nrepresen) then 
    !
    write(out,"('symmetry: number of elements in _select_gamma_ is too small: ',i8)") 100 ! size(job%select_gamma)
    stop 'symmetry: size of _select_gamma_ is too small'
    !
  endif 
  !
  sym%maxdegen = maxval(sym%degen(:),dim=1)

  !
  ! store the address of the group generator from ioper = 1..Noper list 
  !
  ioper = 1
  !
  do iclass = 1,sym%Nclasses
    !
    sym%igenerator(iclass) = ioper
    ioper = ioper + sym%Nelements(iclass)
    !
  enddo
  
  !
  ! check 
  call check_characters_and_representation
  !

  contains

   subroutine simple_arrays_allocation

    integer(ik) :: alloc,nCII

    nCII = max(1,sym%CII%Noper)
    !
    allocate (sym%characters(sym%Nrepresen,sym%Nclasses),sym%irr(sym%Nrepresen,sym%Noper),&  
              sym%degen(sym%Nrepresen),sym%Nelements(sym%Nclasses),sym%label(sym%Nrepresen),&
              sym%igenerator(sym%Nclasses),&
              sym%CII%ioper(nCII),sym%CII%coeff(nCII),sym%euler(sym%Noper,3),sym%lquant(sym%Nrepresen),stat=alloc)

    if (alloc/=0) stop 'simple_arrays_allocation - out of memory'
    !
    sym%CII%coeff = 0
    sym%CII%ioper = 1
    sym%euler = 0
    sym%lquant = 0 
    !
   end subroutine simple_arrays_allocation



   subroutine irr_allocation

    integer(ik) :: gamma,ioper,iclass,ielem,alloc

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 


   if (alloc/=0) then
       write (out,"(' symmetryInitialize ',i9,' error trying to allocate symmetry')") alloc
      stop 'symmetryInitialize, symmetries - out of memory'
   end if


   end subroutine irr_allocation


 end subroutine symmetryInitialize



   subroutine check_characters_and_representation

    integer(ik) :: igamma,jgamma,ioper,ideg,iclass,ielem
    real(ark)   :: temp

    do igamma = 1,sym%Nrepresen
      do jgamma = igamma,sym%Nrepresen
        !
        temp = sum(real(sym%Nelements(:),ark)*sym%characters(igamma,:)*sym%characters(jgamma,:))
        !
        if (igamma/=jgamma.and.abs(temp)>sqrt(small_)) then 
          write (out,"(' check_characters_and_representation: not orhogonal for igamma,jgamma = ',2i4,' -> ',g18.6)") igamma,jgamma,temp
          stop 'check_characters_and_representation: not orhogonal'
        endif
        !
        if (igamma==jgamma.and.abs(temp-sym%Noper)>sqrt(small_)) then 
          write (out,"(' check_characters_and_representation: dot procut ',f16.2,' for isym = ',i4,' /= size of the group ',f16.0)") igamma,temp,sym%Noper
          stop 'check_characters_and_representation: not orhogonal'
        endif
        !
      enddo 
    enddo 

    !
    ! Check characters
    !
    if (verbose_>=5) write(out,"('Irrep matrices:')")
    !
    ioper = 0
    do iclass = 1,sym%Nclasses
      !
      do ielem = 1,sym%Nelements(iclass)
        !
        ioper = ioper + 1
        !
        do igamma = 1,sym%Nrepresen
          !
          temp = 0
          !
          do ideg = 1,sym%degen(igamma)
            !
            temp = temp + sym%irr(igamma,ioper)%repres(ideg,ideg)
            !
          enddo
          !
          if (verbose_>=5) then
            write(out,"('igamma,iclass,ioper = ',3i6)") igamma,iclass,ioper
            do ideg = 1,sym%degen(igamma)
              write(out,"(<sym%degen(igamma)>f18.8)") sym%irr(igamma,ioper)%repres(ideg,:)
            enddo
          endif
          !
          if (abs(temp-sym%characters(igamma,iclass))>sqrt(small_)) then 
            write (out,"(' symmetry: character and representation do not agree ',2f18.8,', igamma,iclass,ioper = ',3i5)") sym%characters(igamma,iclass),temp,igamma,iclass,ioper
            stop 'symmetry: character and representation do not agree'
          endif
          !
        enddo
        !
      enddo
      !
    enddo


   end subroutine check_characters_and_representation

end module symmetry



