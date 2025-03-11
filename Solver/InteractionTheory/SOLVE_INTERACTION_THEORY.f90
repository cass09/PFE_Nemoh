!--------------------------------------------------------------------------------------
!
! Resolution of interaction theory
!
!--------------------------------------------------------------------------------------
MODULE SOLVE_INTERACTION_THEORY

  USE Constants
  USE MMesh,              ONLY: TMesh
  USE MFace,              ONLY: TVFace
  USE MEnvironment,       ONLY: TEnvironment

  ! Solver for linear problem
  USE M_SOLVER,           ONLY:GAUSSZ,LU_INVERS_MATRIX,GMRES_SOLVER, &
                               ID_GAUSS,ID_GMRES,TSolver

  IMPLICIT NONE

  PUBLIC :: SOLVE_POTENTIAL_MATRIX

  PRIVATE
  ! Those variables will be conserved between calls of the subroutine.
  REAL :: Omega_previous = -1.0


CONTAINS

  SUBROUTINE SOLVE_POTENTIAL_MATRIX             &
  ( VFace, Mesh, Env, omega, wavenumber, &
    NVel, S,V,Vinv,Potential,SolverOpt,wd)

  ! Input/output
  TYPE(TVFace),                                   INTENT(IN)    :: VFace
  TYPE(TMesh),                                    INTENT(IN)    :: Mesh
  TYPE(TEnvironment),                             INTENT(IN)    :: Env
  REAL,                                           INTENT(IN)    :: omega, wavenumber
  TYPE(TSolver),                                  INTENT(IN)    :: SolverOpt
  COMPLEX, DIMENSION(Mesh%Npanels*2**Mesh%Isym),  INTENT(IN)    :: NVel
  COMPLEX, DIMENSION(Mesh%Npanels,Mesh%Npanels,2**Mesh%Isym),                     &
                                                  INTENT(INOUT) :: V,S,Vinv ! Influence Coef
  COMPLEX, DIMENSION(Mesh%Npanels*2**Mesh%Isym),  INTENT(OUT)   :: Potential

  INTEGER :: I, J,FLAG_CAL
  CHARACTER(LEN=*),                               INTENT(IN)  :: wd

  FLAG_CAL=0
  
  ! Attention cas profondeur infini ou pas
  ! Attention cas symétrique ou pas

  ! Question : boucle sur les fréquences ici ou avant appel fonction ?
  ! Question : création d'un autre exécutable ou pas nécessaire ?

  ! PLAN : 
  ! 1 - recupérer Madd, Crad et Fex pour un solide isolé :
      ! avec fonction commune ReadTresults ou autre ? 
      ! solve BEM recup que Potential différencié entre scattering et radiation !!
      ! Question : mesh = 1 solide -> prise en compte organisation de la ferme ?
      !          ou tous les solides sont configurés dans quel cas appel solve BEM que pour un !
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
          ! - calcul matrice transformation T
          ! - troncature ???
          ! Besoin des fonctions de Bessel
          ! - calculer Fex à partir scatting (theory OK)
          ! - calculer Madd, Crad en passant par FR à partir radiation (theory to check)

  ! Attention code python adpatée à la géométrie du cylindre -> ajustement à faire ? utiliser index panel et non coord cylindriques ?

  IF (omega /= omega_previous) THEN
      ! Do not recompute if the same frequency is studied twice
      omega_previous = omega
      FLAG_CAL=1
  END IF 

  RETURN
  END SUBROUTINE

END MODULE
