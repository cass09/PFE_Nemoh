!--------------------------------------------------------------------------------------
!
! Resolution of interaction theory

 ! Attention cas profondeur infini ou pas
  ! Attention cas symétrique ou pas

  ! Question : boucle sur les fréquences ici ou avant appel fonction ?
  ! Question : création d'un autre exécutable ou pas nécessaire ?
  ! Nouvelles données input (nombre et positions des WEC; waves direction of the entire system
  ! and discrétisation dtheta et dz ? ) -> ajout dans Nemoh.cal ou création autre fichier spécifique ?

  ! PLAN : 
  ! 1 - recupérer Madd, Crad, Fex et potentiels pour un solide isolé :
      ! solve BEM recup que Potential scattering et radiation sur cylindrical mesh
      ! coeffs recup dans fichiers output ou avec fonction commune ReadTresults ou autre ? 
      ! mesh = 1 solide -> solve BEM que sur un puis interaction theory sur l'ensemble !
  ! 2 - calcul D, G et radiation coefs (fonction tranfert dans code python)
          ! - calcul a_s_scat à partir du flux phi_scat 
          ! - calcul a_s_rad à partir du flux phi_rad
          ! Besoin des fonctions de Bessel
          ! - calcul a_i
          ! - Resolution a_i * D = a_s_scat 
          ! - Resolution a_i * G = fex
          ! Besoin solveur LU, GMRES 
          ! - troncature et réduction ???
  ! 3 - calcul Madd, Crad, Fex pour l'ensemble du système (inspiration code python)
          ! Attention : transposée ou pas ??
          ! recup params sur la ferme
          ! - calcul matrice transformation T
          ! - troncature ???
          ! Besoin des fonctions de Bessel
          ! - calculer Fex à partir scatting (theory OK)
          ! - calculer Madd, Crad en passant par FR à partir radiation (theory OK)


!
!--------------------------------------------------------------------------------------
MODULE SOLVE_INTERACTION_THEORY

  USE Constants
  USE MEnvironment,       ONLY: TEnvironment
  USE M_SOLVER,           ONLY:GAUSSZ,LU_INVERS_MATRIX,GMRES_SOLVER, LU_SOLVER, &
                               ID_GAUSS,ID_GMRES,TSolver
  USE MBESSEL,            ONLY:fun_BESSJ
  USE Elementary_functions, ONLY: X0
  IMPLICIT NONE

  PUBLIC :: SOLVE_POTENTIAL_MATRIX


  TYPE Tcylsurface 
    INTEGER   :: Ntheta, Nz
    REAL      :: R
  END TYPE

  TYPE TInteractionTheory
    INTEGER       :: Nb !Nb bodies
    REAL, DIMENSION(:,:),     ALLOCATABLE :: Coord   ! coord bodies center Nbx2 (même altitude)
    REAL,DIMENSION(:),        ALLOCATABLE :: wave_dir
    REAL,DIMENSION(:,:,:),    ALLOCATABLE :: Madd, Crad
    COMPLEX,DIMENSION(:,:,:), ALLOCATABLE :: Fex
    COMPLEX,DIMENSION(:,:,:), ALLOCATABLE :: D, G
    TYPE(Tcylsurface)   :: Mcyl
  END TYPE TInteractionTheory

  PRIVATE :: ReadTInteractionTheory, CloseIT, &
   ReadIndex, ReadParamIso, ReadTwoPotential, CalculMatrix, SolveITproblem
  ! OneToTwoPotential, CalculParamIso
  
  ! Those variables will be conserved between calls of the subroutine.
  INTEGER :: Nw, Nbeta_iso, Nrad, Nint, Ndir, Mtronc
  REAL, DIMENSION(:), ALLOCATABLE              :: omega

CONTAINS

  SUBROUTINE SOLVE_POTENTIAL_MATRIX             &
  (Env,SolverOpt, wd)
  IMPLICIT NONE

  TYPE(TEnvironment),                   INTENT(IN) :: Env
  TYPE(TSolver),                        INTENT(IN) :: SolverOpt
  ! COMPLEX, DIMENSION(Nproblems, Mesh%Npoints),  INTENT(IN) :: Potential          ! Computed potential
  ! INTEGER,                              INTENT(IN) :: Nproblems
  ! INTEGER,DIMENSION(Nproblems),         INTENT(IN) :: Switch_Type
  CHARACTER(LEN=*),                     INTENT(IN) :: wd

  TYPE(TInteractionTheory)              :: ParamsIT
  REAL, DIMENSION(:),       ALLOCATABLE :: beta_iso
  REAL,DIMENSION(:,:,:),    ALLOCATABLE :: Madd_iso, Crad_iso
  COMPLEX,DIMENSION(:,:,:), ALLOCATABLE :: Fex_iso
  
  WRITE(*,*) "-------------------Début"
  CALL ReadTInteractionTheory(ParamsIT, wd)
  CALL ReadIndex(beta_iso, wd)
  
  ALLOCATE(Fex_iso(Nw, Nbeta_iso, Nint), Madd_iso(Nw, Nrad, Nint), Crad_iso(Nw, Nrad, Nint))
  CALL ReadParamIso(Fex_iso, Madd_iso, Crad_iso, wd)
  WRITE(*,*) "-----------------Params ok"
  
  ALLOCATE(ParamsIT%D(Nw, 2*Mtronc+1, 2*Mtronc+1), ParamsIT%G(Nw, 2*Mtronc+1, Nint))

  CALL CalculMatrix  &
  (Env, SolverOpt, ParamsIT, beta_iso, Fex_iso, wd)

  CALL SolveITproblem(ParamsIT, Madd_iso, Crad_iso)

  DEALLOCATE(beta_iso, omega)
  DEALLOCATE(Fex_iso, Crad_iso, Madd_iso)
  WRITE(*,*) "-------------------END"
  
  CALL CloseIT(ParamsIT)

  END SUBROUTINE SOLVE_POTENTIAL_MATRIX

  SUBROUTINE CloseIT(ParamsIT)

  TYPE(TInteractionTheory),       INTENT(INOUT) :: ParamsIT

  IF (ALLOCATED(ParamsIT%Coord)) DEALLOCATE(ParamsIT%Coord)
  IF (ALLOCATED(ParamsIT%wave_dir)) DEALLOCATE(ParamsIT%wave_dir)
  IF (ALLOCATED(ParamsIT%Fex)) DEALLOCATE(ParamsIT%Fex)
  IF (ALLOCATED(ParamsIT%Madd)) DEALLOCATE(ParamsIT%Madd)
  IF (ALLOCATED(ParamsIT%Crad)) DEALLOCATE(ParamsIT%Crad)
  IF (ALLOCATED(ParamsIT%D)) DEALLOCATE(ParamsIT%D)
  IF (ALLOCATED(ParamsIT%G)) DEALLOCATE(ParamsIT%G)

  END SUBROUTINE CloseIT

  ! SUBROUTINE ReadTInteractionTheory             &
  ! (InpNEMOHCAL, InputIT)

  ! IMPLICIT NONE
  ! TYPE(TNemCal),                :: InpNEMOHCAL
  ! TYPE(TInteractionTheory),    INTENT(OUT) :: InputIT

  ! INTEGER :: j, k
  ! REAL    :: dir_min, dir_max

  ! InputIT%Nb =InpNEMOHCAL%IntTheory%Nb
  ! ALLOCATE(InputIT%Coord(InputIT%Nb, 2))
  ! InputIT%Coord =InpNEMOHCAL%IntTheory%Bcoord
  ! Ndir =InpNEMOHCAL%IntTheory%NDir
  ! ALLOCATE(InputIT%wave_dir(Ndir))
  ! dir_max =InpNEMOHCAL%IntTheory%DirMax
  ! dir_min =InpNEMOHCAL%IntTheory%DirMin
  ! DO j=1,Ndir
  !   InputIT%wave_dir(j)=(dir_min+(dir_max-dir_min)*(j-1)/(Ndir-1))*PI/180.
  ! END DO

  ! END SUBROUTINE ReadTInteractionTheory

  SUBROUTINE ReadTInteractionTheory             &
  (InputIT, wd)

  IMPLICIT NONE
  TYPE(TInteractionTheory),       INTENT(OUT) :: InputIT
  CHARACTER(LEN=*),               INTENT(IN)  :: wd 

  INTEGER :: j, k
  INTEGER :: Nb_NEMOHcal
  REAL    :: dir_min, dir_max

  OPEN(15,FILE=TRIM(wd)//'/input_IT.dat')
    READ(15,*) 
    READ(15,*) 
    READ(15,*) InputIT%Nb
    ALLOCATE(InputIT%Coord(InputIT%Nb, 2))
    DO k=1,InputIT%Nb
        READ(15,*) InputIT%Coord(k,1),InputIT%Coord(k,2)
    END DO
    READ(15,*) Ndir, dir_min, dir_max
    READ(15,*) InputIT%Mcyl%R, InputIT%Mcyl%Ntheta, InputIT%Mcyl%Nz
  CLOSE(15)
  ALLOCATE(InputIT%wave_dir(Ndir))
  DO j=1,Ndir
    IF (Ndir > 1) THEN
        InputIT%wave_dir(j) = (dir_min + (dir_max - dir_min) * (j - 1) / (Ndir - 1)) * PI / 180.
    ELSE
        InputIT%wave_dir(j) = dir_min * PI / 180.  
    END IF
  END DO
  ALLOCATE(InputIT%Madd(Nw, InputIT%Nb*Nrad, InputIT%Nb*Nint))
  ALLOCATE(InputIT%Crad(Nw, InputIT%Nb*Nrad, InputIT%Nb*Nint))
  ALLOCATE(InputIT%Fex(Nw, Ndir, InputIT%Nb*Nint))
  END SUBROUTINE ReadTInteractionTheory


  SUBROUTINE ReadIndex &
  (beta_iso, wd)
  IMPLICIT NONE 

  REAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: beta_iso
  CHARACTER(LEN=*),                INTENT(IN)  :: wd 

  INTEGER :: k

  OPEN(17,FILE=TRIM(wd)//'/results/index.dat')
    READ(17,*) Nw,Nbeta_iso,Nrad,Nint
    READ(17,*)
    DO k=1,Nint
        READ(17,*) 
    END DO
    READ(17,*)
    DO k=1,Nrad
        READ(17,*)
    END DO
    ALLOCATE(beta_iso(Nbeta_iso), omega(Nw))
    READ(17,*) (beta_iso(k),k=1,Nbeta_iso)
    READ(17,*) (omega(k),k=1,Nw)
  CLOSE(17)
  Mtronc = INT(Nbeta_iso-1)/2

  END SUBROUTINE ReadIndex


  ! SUBROUTINE CalculParamIso   &
  ! (Env, PHI_S, PHI_R)
  ! IMPLICIT NONE

  ! TYPE(TEnvironment),          INTENT(IN) :: Env
  ! COMPLEX, DIMENSION(:,:,:),   INTENT(IN) :: PHI_S, PHI_R  

  ! ! calcul momentum avec rho, phi, normal   
  ! ! calcul Fex, Madd, Crad avec omega et momentum
  ! ! + dimensions
  ! ! Tout ok sauf normal à lire dans fichier mais boucle Nproblème

  ! END SUBROUTINE CalculParamIso

  SUBROUTINE ReadParamIso   &
  (Fex, Madd, Crad, wd)
  IMPLICIT NONE

  CHARACTER(LEN=*),                     INTENT(IN)  :: wd 
  REAL,DIMENSION(Nw,Nrad,Nint),         INTENT(OUT) :: Madd, Crad
  COMPLEX,DIMENSION(Nw,Nbeta_iso,Nint), INTENT(OUT) :: Fex

  REAL,DIMENSION(:),ALLOCATABLE :: line
  INTEGER :: i, j, k, c

  ALLOCATE(line(2*Nint))
  OPEN(18,FILE=TRIM(wd)//'/results/Forces.dat')
  READ(18,*)
    DO i=1,Nw
      DO j=1,Nbeta_iso
        READ(18,*) (line(c),c=1,2*Nint)
        DO k=1,Nint
            Fex(i,j,k)=line(2*k-1)*CEXP(CMPLX(0.,1.)*line(2*k))
        END DO
      END DO
      DO j=1,Nrad
        READ(18,*) (line(c),c=1,2*Nint)
        DO k=1,Nint
            Madd(i,j,k)=line(2*k-1)
            Crad(i,j,k)=line(2*k)
        END DO
      END DO
    END DO
  CLOSE(18)
  DEALLOCATE(line)

  ALLOCATE(line(2*Nint+1))
  OPEN(19,FILE=TRIM(wd)//'/results/FKForce.tec')
  READ(19,*)
  DO k=1,Nint
      READ(19,*)
  END DO
  DO j=1,Nbeta_iso
      READ(19,*)
      DO i=1,Nw
          READ(19,*) (line(k),k=1,1+2*Nint)
          c=2
          DO k=1,Nint
              ! Results%FroudeKrylovForce(i,j,k)=line(c)*CEXP(CMPLX(0.,1.)*line(c+1))
              Fex(i,j,k)=Fex(i,j,k)+line(c)*CEXP(CMPLX(0.,1.)*line(c+1))
              c=c+2
          END DO
      END DO
  END DO
  CLOSE(19)
  DEALLOCATE(line)
  END SUBROUTINE ReadParamIso

  ! SUBROUTINE OneToTwoPotential             &
  ! (Nmesh, Potential, Switch_type, PHI_S, PHI_R)
  ! IMPLICIT NONE
  ! ! Input/output
  ! INTEGER,                                        INTENT(IN)  :: Nmesh
  ! COMPLEX, DIMENSION(Nw*(Nbeta_iso+Nrad), Nmesh), INTENT(IN)  :: Potential          ! Computed potential
  ! INTEGER,DIMENSION(Nw*(Nbeta_iso+Nrad)),         INTENT(IN)  :: Switch_Type
  ! COMPLEX, DIMENSION(Nw,Nbeta_iso,Nmesh),       INTENT(INOUT) :: PHI_S
  ! COMPLEX, DIMENSION(Nw,Nrad,Nmesh),            INTENT(INOUT) :: PHI_R

  ! COMPLEX, DIMENSION(:,:), ALLOCATABLE         :: PHI_Sw, PHI_Rw  
  ! INTEGER :: i, j, indice_S, indice_R

  ! ALLOCATE(PHI_Rw(Nw*Nrad, Nmesh), PHI_Sw(Nw*Nbeta_iso, Nmesh))
  ! indice_S = 1
  ! indice_R = 1
  ! DO i=1, Nw*(Nbeta_iso+Nrad)
  !   IF (switch_type(i) == DIFFRACTION_PROBLEM) THEN
  !         PHI_Sw(indice_S,:) = Potential(i,:) 
  !         indice_S = indice_S+1
  !   ELSE IF (switch_type(i) == RADIATION_PROBLEM) THEN
  !         PHI_Rw(indice_R,:) = Potential(i,:)  
  !         indice_R=indice_R+1      
  !   END IF
  ! END DO
  ! indice_S = 1
  ! indice_R = 1
  ! DO i=1,Nw  ! boucle à étendre jusqu'à la fin ??
  !     DO j=1,Nbeta_iso
  !       PHI_S(i, j,:) = PHI_Sw(indice_S,:)
  !       indice_S = indice_S+1
  !     END DO
  !     DO j=1,Nrad
  !       PHI_R(i,j,:) = PHI_Rw(indice_R,:)
  !       indice_R=indice_R+1 
  !     END DO
  ! END DO
  ! DEALLOCATE(PHI_Rw, PHI_Sw)
  ! ! phi_S et phi_R OK (Nw, N_, Npanels)               
  ! END SUBROUTINE OneToTwoPotential

  SUBROUTINE ReadTwoPotential             &
  (PHI_S, PHI_R, Mcyl, G, wd)
  IMPLICIT NONE

  TYPE(Tcylsurface),                            INTENT(IN)  :: Mcyl
  REAL,                                         INTENT(IN)  :: G
  COMPLEX, DIMENSION(Nw,Nbeta_iso,Mcyl%Ntheta, Mcyl%Nz),INTENT(INOUT) :: PHI_S
  COMPLEX, DIMENSION(Nw,Nrad,Mcyl%Ntheta, Mcyl%Nz),     INTENT(INOUT) :: PHI_R
  CHARACTER(LEN=*),                             INTENT(IN)  :: wd 



  COMPLEX, DIMENSION(Nw*(Nrad+Nbeta_iso), Mcyl%Ntheta, Mcyl%Nz) :: ETAc 
  CHARACTER(LEN=30) :: filename 
  INTEGER :: Np, i, j, k, indice_p, unit
  REAL    :: x, y, z, eta_amp, eta_phase

  Np = Nw*(Nrad+Nbeta_iso)
  DO i=1, Np
    print*, "Unpack problem ", i, " / ", Np
    ! Construire le nom du fichier
    write(filename, '(A,I5.5,A)') TRIM(wd)//'/results/cylsurface.', i, '.dat'
    unit = 20 + i  
    open(unit=unit, file=trim(filename))
    READ(unit,*)
    READ(unit,*)
    DO j=1, Mcyl%Ntheta
      DO k=1, Mcyl%Nz
        READ(unit,*) x, y, z, eta_amp, eta_phase
        ETAc(i,j,k)=eta_amp*EXP(II*eta_phase)
      END DO
    END DO
    CLOSE(unit) 
  END DO     
  indice_p=1
  DO i=1, Nw
    ETAc(indice_p, :, :)=ETAc(indice_p, :, :)*G/(II*omega(i))
    DO j=1, Nbeta_iso
      PHI_S(i,j, :, :)=ETAc(indice_p, :, :)
      indice_p=indice_p+1
    END DO 
    DO j=1, Nrad
      PHI_R(i,j, :, :)=ETAc(indice_p, :, :)
      indice_p=indice_p+1
    END DO 
  END DO 

  ! freesurface
  ! ETA(j) = II*omega/Env%G*PHI

  ! python
  ! Scattering[ind] *= g/1j/fr
  ! Radiation[ind] *= g/1j/fr*(-1j*fr)

  END SUBROUTINE ReadTwoPotential



  SUBROUTINE CalculMatrix             &
  (Env, SolverOpt, ParamsIT, beta_iso, Fex_iso, wd)

  IMPLICIT NONE
  TYPE(TInteractionTheory),                      INTENT(INOUT) :: ParamsIT
  TYPE(TEnvironment),                             INTENT(IN)   :: Env
  TYPE(TSolver),                                  INTENT(IN)   :: SolverOpt
  REAL, DIMENSION(Nbeta_iso),                     INTENT(IN)   :: beta_iso
  COMPLEX,DIMENSION(Nw,Nbeta_iso,Nint),           INTENT(IN)   :: Fex_iso
  CHARACTER(LEN=*),                               INTENT(IN)   :: wd 


  COMPLEX, DIMENSION(Nw, Nrad, ParamsIT%Mcyl%Ntheta, ParamsIT%Mcyl%Nz)      :: PHI_R 
  COMPLEX, DIMENSION(Nw, Nbeta_iso, ParamsIT%Mcyl%Ntheta, ParamsIT%Mcyl%Nz) :: PHI_S 
  COMPLEX, DIMENSION(Nw, Nrad, 2*Mtronc+1)       :: A_RAD
  COMPLEX, DIMENSION(Nbeta_iso, 2*Mtronc+1)      :: A_SCAT, A_I
  INTEGER :: i, m, i_m
  REAL    :: coef, k
  COMPLEX :: Hankel_2
  COMPLEX, DIMENSION(:), ALLOCATABLE :: int_R, int_S

  CALL ReadTwoPotential(PHI_S, PHI_R, ParamsIT%Mcyl, Env%G, wd)
  ALLOCATE(int_R(Nrad), int_S(Nbeta_iso))
  DO i=1, Nw 
    ! calcul wave number
    IF (omega(i)**2*Env%depth/Env%g >= 20) THEN
      k = omega(i)**2/Env%g
    ELSE
      k = X0(omega(i)**2*Env%depth/Env%g)/Env%depth
      ! X0(y) returns the solution of y = x * tanh(x)
    END IF

    ! calcul coefficient a
    coef=2*COSH(k*Env%Depth)/(Env%Depth*(1+SINH(2*k*Env%Depth)/(2*k*Env%Depth)))
    coef=coef*(-omega(i)/(2*PI*Env%G)) ! signe opposé sur python
    DO i_m=1, 2*Mtronc+1
      m=i_m-Mtronc-1
      ! fonction de Hankel d'ordre 2
      ! Hm(2)= -iJm (exp(i*pi*m)-(-1)^m)/sin(pi*m) car J-m = (-1)^m Jm
      Hankel_2=-II*fun_BESSJ(m, k*ParamsIT%Mcyl%R)*(EXP(II*PI*m)-(-1)**m)/SIN(PI*m)
      int_R=CALCUL_INT_A(m, k, Env%Depth, PHI_R(i,:,:,:), ParamsIT%Mcyl)
      int_S=CALCUL_INT_A(m, k, Env%Depth, PHI_S(i,:,:,:), ParamsIT%Mcyl)
      A_RAD(i,:,i_m)=II*coef*int_R(:)/Hankel_2
      A_SCAT(:,i_m)=II*coef*int_S(:)/Hankel_2  
      A_I(:,i_m)=EXP(-II*m*(beta_iso(:)+PI/2)) ! theory
      ! A_I(:,i_m)=EXP(II*m*(-beta_iso(:)+PI/2)) ! convention Nemoh python -> le conjugué
    END DO
    
    ! Nmode = 2*Mtruc+1     = Nbeta_iso si impair !!
    ! LU_solver OK que pour une matrice carrée
    CALL LU_SOLVER(A_I,A_SCAT,ParamsIT%D(i,:,:),Nbeta_iso,2*Mtronc+1,Nbeta_iso)
    CALL LU_SOLVER(A_I,Fex_iso(i,:,:),ParamsIT%G(i,:,:),Nbeta_iso,2*Mtronc+1,Nint)

    ! réduction ?
  END DO
  WRITE(*,*) "-------- Calcul Matrix OK"

    ! - calcul a_s_scat à partir du flux phi_scat 
    ! - calcul a_s_rad à partir du flux phi_rad
    ! Besoin des fonctions de Bessel  -> PROBLEME !!
    ! - calcul a_i
    ! - Resolution a_i * D = a_s_scat 
    ! - Resolution a_i * G = fex
    ! Besoin solveur LU, GMRES 
    ! - troncature = modes pour les directions des vagues avec Hm
    ! - réduction pour enlever modes non significatifs

  END SUBROUTINE CalculMatrix

  FUNCTION CALCUL_INT_A(m, k, h, PHI, Mcyl) RESULT(INT_A)
  IMPLICIT NONE

  INTEGER                    :: m
  TYPE(Tcylsurface)          :: Mcyl
  REAL                       :: k, h
  COMPLEX, DIMENSION(:,:,:)  :: PHI     ! N x Ntheta x Nz
  COMPLEX, DIMENSION(SIZE(PHI, 1))   :: INT_A

  REAL, DIMENSION(Mcyl%Ntheta) :: theta
  REAL, DIMENSION(Mcyl%Nz)     :: z
  REAL    :: dth, dz
  INTEGER :: i, j
  COMPLEX, DIMENSION(SIZE(PHI, 1), SIZE(PHI, 2), SIZE(PHI, 3)) :: COEF
  COMPLEX, DIMENSION(SIZE(PHI, 1), SIZE(PHI, 3)) :: INT_TH
  COMPLEX, DIMENSION(SIZE(PHI, 1)) :: INT_Z

  ! à mettre avant boucle sur m ?
  DO i=1,Mcyl%Ntheta
    theta(i)=2*PI*(i-1)/(Mcyl%Ntheta-1)
  END DO
  dth=2*PI/Mcyl%Ntheta
  DO i=1,Mcyl%Nz
    z(i)=-h*(i-1)/(Mcyl%Nz-1)
  END DO
  dz=h/Mcyl%Nz

  DO j=1,Mcyl%Nz
    DO i=1,Mcyl%Ntheta
        COEF(:, i, j) = COSH(k*(h+z(j)))*EXP(-II*m*theta(j))*PHI(:, i, j)
    END DO
  END DO
  INT_TH=0
  INT_Z=0
  DO j=1,Mcyl%Nz
    DO i=1,Mcyl%Ntheta-1
        INT_TH(:,j)=INT_TH(:,j)+(COEF(:,i, j)+COEF(:,i,j+1))*dth/2
    END DO
  END DO
  DO j=1,Mcyl%Nz-1
    INT_Z(:)=INT_Z(:)+(INT_TH(:,j)+INT_TH(:,j+1))*dz/2
  END DO

  INT_A=INT_Z
  END FUNCTION CALCUL_INT_A


  SUBROUTINE SolveITproblem             &
  (ParamsIT, Madd_iso, Crad_iso)

  IMPLICIT NONE
  TYPE(TInteractionTheory),    INTENT(INOUT) :: ParamsIT
  REAL,DIMENSION(Nw,Nrad,Nint),  INTENT(IN)  :: Madd_iso, Crad_iso

    
  WRITE(*,*) "-------- Solve IT"
  ! 3 - calcul Madd, Crad, Fex pour l'ensemble du système (inspiration code python)
          ! Attention : transposée ou pas ??
          ! - calcul matrice transformation T
          ! - troncature ???
          ! Besoin des fonctions de Bessel
          ! - calculer Fex à partir scatting (theory OK)
          ! - calculer Madd, Crad en passant par FR à partir radiation (theory OK)
  END SUBROUTINE SolveITproblem

END MODULE
