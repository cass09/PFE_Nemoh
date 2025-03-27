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
      ! solve BEM recup que Potential scattering et radiation
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

  ! Attention code python adpatée à la géométrie du cylindre 
  ! -> ajustement à faire ? utiliser index panel et non coord cylindriques ?


!
!--------------------------------------------------------------------------------------
MODULE SOLVE_INTERACTION_THEORY

  USE Constants
  USE MMesh,              ONLY: TMesh
  USE MFace,              ONLY: TVFace
  USE MEnvironment,       ONLY: TEnvironment
  USE M_SOLVER,           ONLY:GAUSSZ,LU_INVERS_MATRIX,GMRES_SOLVER, &
                               ID_GAUSS,ID_GMRES,TSolver

  IMPLICIT NONE

  PUBLIC :: SOLVE_POTENTIAL_MATRIX, ReadTInteractionTheory, CloseIT, &
   ReadIndex, ReadParamIso, OneToTwoPotential, CalculParamIso, &
   CalculMatrix, SolveITproblem

  TYPE TInteractionTheory
    INTEGER       :: run_IT, Nb   !Nb bodies
    REAL, DIMENSION(:,:) , ALLOCATABLE :: Coord   ! coord bodies center Nbx2 (même altitude)
    REAL,DIMENSION(:)    , ALLOCATABLE :: wave_dir
    REAL,DIMENSION(:,:,:), ALLOCATABLE :: Madd, Crad
    COMPLEX,DIMENSION(:,:,:), ALLOCATABLE :: Fex
    REAL,DIMENSION(:,:,:), ALLOCATABLE :: D, G
  END TYPE TInteractionTheory

  PRIVATE
  ! Those variables will be conserved between calls of the subroutine.
  INTEGER :: Nw, Nbeta_iso, Nrad, Nint, Ndir
  REAL, DIMENSION(:), ALLOCATABLE              :: omega

CONTAINS

  SUBROUTINE SOLVE_POTENTIAL_MATRIX             &
  (VFace, Mesh, Env,SolverOpt, Nproblems, Potential, Switch_type, wd)
  IMPLICIT NONE

  TYPE(TVFace),                                   INTENT(IN) :: VFace
  TYPE(TMesh),                                    INTENT(IN) :: Mesh
  TYPE(TEnvironment),                             INTENT(IN) :: Env
  TYPE(TSolver),                                  INTENT(IN) :: SolverOpt
  COMPLEX, DIMENSION(Nproblems, Mesh%Npoints),    INTENT(IN) :: Potential          ! Computed potential
  INTEGER,                                        INTENT(IN) :: Nproblems
  INTEGER,DIMENSION(Nproblems),                   INTENT(IN) :: Switch_Type
  CHARACTER(LEN=*),                               INTENT(IN) :: wd

  TYPE(TInteractionTheory)                     :: ParamsIT
  REAL, DIMENSION(:), ALLOCATABLE              :: beta_iso
  REAL,DIMENSION(:,:,:), ALLOCATABLE :: Madd_iso, Crad_iso
  COMPLEX,DIMENSION(:,:,:), ALLOCATABLE :: Fex_iso
  !  use Nmesh au lieu de Mesh entier parfois
  
  CALL ReadTInteractionTheory(ParamsIT, wd)

  IF (ParamsIT%run_IT==1) THEN
    WRITE(*,*) "-------------------Début"
    CALL ReadIndex(beta_iso, wd)
   
    ALLOCATE(Fex_iso(Nw, Nbeta_iso, Nint), Madd_iso(Nw, Nrad, Nint), Crad_iso(Nw, Nrad, Nint))
    CALL ReadParamIso(Fex_iso, Madd_iso, Crad_iso, wd)
    WRITE(*,*) "----------------ReadParams Iso OK"

    CALL CalculMatrix  &
    (Env, Mesh%Npoints, Potential, Switch_Type, ParamsIT, beta_iso, Fex_iso)

    CALL SolveITproblem(ParamsIT, Madd_iso, Crad_iso)

    DEALLOCATE(beta_iso, omega)
    DEALLOCATE(Fex_iso, Crad_iso, Madd_iso)
    WRITE(*,*) "-------------------END"
  ELSE 
    WRITE(*,*) "--------- Interaction Theory not activated ------------"
    
  END IF
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

  SUBROUTINE ReadTInteractionTheory             &
  (InputIT, wd)

  IMPLICIT NONE
  TYPE(TInteractionTheory),       INTENT(OUT) :: InputIT
  CHARACTER(LEN=*),               INTENT(IN)  :: wd 

  INTEGER :: j, k
  INTEGER :: Nb_NEMOHcal
  REAL    :: dir_min, dir_max

  OPEN(15,FILE=TRIM(wd)//'/input_IT.dat')
    READ(15,*) InputIT%run_IT
    ! IF (ParamsIT%run_IT==0) THEN
    !   WRITE(*,*) "--------- Interaction Theory not activated ------------"
    !   CLOSE(15)
    !   RETURN
    ! END IF
    READ(15,*) InputIT%Nb
    ALLOCATE(InputIT%Coord(InputIT%Nb, 2))
    DO k=1,InputIT%Nb
        READ(15,*) InputIT%Coord(k,1),InputIT%Coord(k,2)
    END DO
    READ(15,*) Ndir, dir_min, dir_max
  CLOSE(15)
  ALLOCATE(InputIT%wave_dir(Ndir))
  DO j=1,Ndir
    InputIT%wave_dir(j)=(dir_min+(dir_max-dir_min)*(j-1)/(Ndir-1))*PI/180.
  END DO

  OPEN(16, FILE=TRIM(wd)//'/Nemoh.cal')
    READ(16,*) !--- Environment ----------------------------!
    READ(16,*) ! RHO
    READ(16,*) ! G
    READ(16,*) ! Depth
    READ(16,*) ! Xeff, Yeff
    READ(16,*) !--- Description of floating bodies-----------!
    READ(16,*) Nb_NEMOHcal
  CLOSE(16)
  ! Récupérer le nombre de body dans le fichier Nemoh.cal
  ! InteractionTheory appliquée que si BEM exécutée sur UN body 
  IF (Nb_NEMOHcal .NE. 1) THEN 
    InputIT%run_IT = 0
    WRITE(*,*) "BEM resolution must be on only ONE body to apply Interaction Theory"
  END IF

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
  
  END SUBROUTINE ReadIndex


  SUBROUTINE CalculParamIso   &
  (Env, PHI_S, PHI_R)
  IMPLICIT NONE

  TYPE(TEnvironment),           INTENT(IN)  :: Env
  COMPLEX, DIMENSION(:,:,:),    INTENT(IN)  :: PHI_S, PHI_R  

  ! calcul momentum avec rho, phi, normal   
  ! calcul Fex, Madd, Crad avec omega et momentum
  ! + dimensions
  ! Tout ok sauf normal à lire dans fichier mais boucle Nproblème

  END SUBROUTINE CalculParamIso

  SUBROUTINE ReadParamIso   &
  (Fex, Madd, Crad, wd)
  IMPLICIT NONE

  CHARACTER(LEN=*),             INTENT(IN)  :: wd 
  REAL,DIMENSION(Nw,Nrad,Nint),        INTENT(OUT) :: Madd, Crad
  COMPLEX,DIMENSION(Nw,Nbeta_iso,Nint),     INTENT(OUT) :: Fex

  REAL,DIMENSION(:),ALLOCATABLE :: line
  INTEGER :: i, j, k, c

  ! WRITE(*,*) "-----------------Param iso"
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
        ! WRITE(*,*) "-----------------Fex ok"
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
  ! WRITE(*,*) "-----------------Param iso ok"
  END SUBROUTINE ReadParamIso

  SUBROUTINE OneToTwoPotential             &
  (Nmesh, Potential, Switch_type, PHI_S, PHI_R)
  IMPLICIT NONE
  ! Input/output
  INTEGER,                                         INTENT(IN) :: Nmesh
  COMPLEX, DIMENSION(Nw*(Nbeta_iso+Nrad), Nmesh),  INTENT(IN) :: Potential          ! Computed potential
  INTEGER,DIMENSION(Nw*(Nbeta_iso+Nrad)),          INTENT(IN) :: Switch_Type
  COMPLEX, DIMENSION(Nw,Nbeta_iso,Nmesh),       INTENT(INOUT) :: PHI_S
  COMPLEX, DIMENSION(Nw,Nrad,Nmesh),            INTENT(INOUT) :: PHI_R


  COMPLEX, DIMENSION(:,:), ALLOCATABLE         :: PHI_Sw, PHI_Rw  
  INTEGER :: i, j, indice_S, indice_R

  ALLOCATE(PHI_Rw(Nw*Nrad, Nmesh), PHI_Sw(Nw*Nbeta_iso, Nmesh))
  ! WRITE(*,*) "allocate ok"
  indice_S = 1
  indice_R = 1
  DO i=1, Nw*(Nbeta_iso+Nrad)
    IF (switch_type(i) == DIFFRACTION_PROBLEM) THEN
          PHI_Sw(indice_S,:) = Potential(i,:) 
          indice_S = indice_S+1
    ELSE IF (switch_type(i) == RADIATION_PROBLEM) THEN
          PHI_Rw(indice_R,:) = Potential(i,:)  
          indice_R=indice_R+1      
    END IF
  END DO
  ! WRITE(*,*) "Boucle 1 ok"
  indice_S = 1
  indice_R = 1
  DO i=1,Nw  ! boucle à étendre jusqu'à la fin ??
      DO j=1,Nbeta_iso
        PHI_S(i, j,:) = PHI_Sw(indice_S,:)
        indice_S = indice_S+1
      END DO
      DO j=1,Nrad
        PHI_R(i,j,:) = PHI_Rw(indice_R,:)
        indice_R=indice_R+1 
      END DO
  END DO
  ! WRITE(*,*) "boucle 2 ok"
  DEALLOCATE(PHI_Rw, PHI_Sw)
  ! phi_S et phi_R OK (Nw, N_, Npanels)               


  END SUBROUTINE OneToTwoPotential





  SUBROUTINE CalculMatrix             &
  (Env, Nmesh, Potential, Switch_Type, ParamsIT, beta_iso, Fex_iso)

  IMPLICIT NONE
  TYPE(TInteractionTheory),                      INTENT(INOUT) :: ParamsIT
  INTEGER,                                        INTENT(IN)   :: Nmesh
  TYPE(TEnvironment),                             INTENT(IN)   :: Env
  COMPLEX, DIMENSION(Nw*(Nbeta_iso+Nrad), Nmesh), INTENT(IN)   :: Potential         
  INTEGER,DIMENSION(Nw*(Nbeta_iso+Nrad)),         INTENT(IN)   :: Switch_Type
  REAL, DIMENSION(Nbeta_iso),                     INTENT(IN)   :: beta_iso
  COMPLEX,DIMENSION(Nw,Nbeta_iso,Nint),           INTENT(IN)   :: Fex_iso

  COMPLEX, DIMENSION(:,:,:), ALLOCATABLE       :: PHI_S, PHI_R 

  ALLOCATE(PHI_R(Nw, Nrad, Nmesh), PHI_S(Nw,Nbeta_iso, Nmesh))
  CALL OneToTwoPotential(Nmesh, Potential, Switch_type, PHI_S, PHI_R)
  WRITE(*,*) "-------- Potential OK"

    ! - calcul a_s_scat à partir du flux phi_scat 
    ! - calcul a_s_rad à partir du flux phi_rad
    ! Besoin des fonctions de Bessel
    ! - calcul a_i
    ! - Resolution a_i * D = a_s_scat 
    ! - Resolution a_i * G = fex
    ! Besoin solveur LU, GMRES 
    ! - troncature et réduction ???


  DEALLOCATE(PHI_R, PHI_S)
  END SUBROUTINE CalculMatrix




  SUBROUTINE SolveITproblem             &
  (ParamsIT, Madd_iso, Crad_iso)

  IMPLICIT NONE
  TYPE(TInteractionTheory),    INTENT(INOUT) :: ParamsIT
  REAL,DIMENSION(Nw,Nrad,Nint),   INTENT(IN) :: Madd_iso, Crad_iso

    

  END SUBROUTINE SolveITproblem

END MODULE
