!--------------------------------------------------------------------------------------
!
!    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!   Contributors list:
!   - G. Delhommeau
!   - J. Singh
!   - P. Guével
!   - J.C. Daubisse
!   - R. Kurnia (2020)
!--------------------------------------------------------------------------------------
!
! Resolution of the boundary elements problem
!
!--------------------------------------------------------------------------------------
MODULE SOLVE_BEM_DIRECT

  USE Constants
  USE MMesh,              ONLY: TMesh
  USE MFace,              ONLY: TVFace
  USE MEnvironment,       ONLY: TEnvironment

  ! Green functions
  USE M_INITIALIZE_GREEN, ONLY: TGREEN, LISC
  USE GREEN_2,            ONLY: VNSINFD, VNSFD

  ! Solver for linear problem
  USE M_SOLVER,           ONLY:GAUSSZ,LU_INVERS_MATRIX,GMRES_SOLVER, &
                               ID_GAUSS,ID_GMRES,TSolver

  IMPLICIT NONE

  PUBLIC :: SOLVE_POTENTIAL_DIRECT

  PRIVATE
  ! Those variables will be conserved between calls of the subroutine.
  REAL :: Omega_previous = -1.0


CONTAINS

  SUBROUTINE SOLVE_POTENTIAL_DIRECT              &
  ( VFace, Mesh, Env, omega, wavenumber, IGreen, &
    NVel, S,V,Vinv,ZIGB, ZIGS,Potential,SolverOpt,wd)

  ! Input/output
  TYPE(TVFace),                                   INTENT(IN)    :: VFace
  TYPE(TMesh),                                    INTENT(IN)    :: Mesh
  TYPE(TEnvironment),                             INTENT(IN)    :: Env
  REAL,                                           INTENT(IN)    :: omega, wavenumber
  TYPE(TGREEN),                                   INTENT(INOUT) :: IGreen
  TYPE(TSolver),                                  INTENT(IN)    :: SolverOpt
  COMPLEX, DIMENSION(Mesh%Npanels*2**Mesh%Isym),  INTENT(IN)    :: NVel
  COMPLEX, DIMENSION(Mesh%Npanels),               INTENT(OUT)   :: ZIGB, ZIGS ! Source distribution
  COMPLEX, DIMENSION(Mesh%Npanels,Mesh%Npanels,2**Mesh%Isym),                     &
                                                  INTENT(INOUT) :: V,S,Vinv ! Influence Coef

  COMPLEX, DIMENSION(Mesh%Npanels*2**Mesh%Isym),  INTENT(OUT)   :: Potential

  INTEGER :: I, J,FLAG_CAL,ITERLID
  ! Return of GREEN_1 module
  REAL :: FSP, FSM
  REAL, DIMENSION(3) :: VSXP, VSXM

  ! Return of GREEN_2 module
  COMPLEX :: SP, SM
  COMPLEX, DIMENSION(3) :: VSP, VSM

  COMPLEX, DIMENSION(Mesh%NPanels, 2**Mesh%ISym) :: ZOL
  COMPLEX, DIMENSION(Mesh%Npanels) :: RHS ! temporary variable
  CHARACTER(LEN=*),                               INTENT(IN)  :: wd

  ! IF (ABS(Omega) <= 1.E-4) THEN
  !   WRITE(*,*)'ABS(Omega)  = ',ABS(Omega),' < 1.E-4'
  !   STOP
  ! END IF

  ! TODO: Clean up those warning messages...
  ! IF (AMH-AKH <= 1e-3) THEN
  !   PRINT*, 'Depth is almost infinite: the infinite depth solver might be more efficient'
  ! END IF
  ! IF (omega**2*depth/g <= 0.1) THEN
  !   PRINT*, 'Depth is too low for the given wavelength'
  ! END IF
  FLAG_CAL=0
  !=====================================
  ! Construction of the influence matrix
  !=====================================
  !OPEN(200, FILE=wd//'/LinearSystem_temp.txt', ACTION='WRITE')

  IF (omega /= omega_previous) THEN
      ! Do not recompute if the same frequency is studied twice
      omega_previous = omega
      FLAG_CAL=1
      ! Initialization of the Green function computations
      IF (.NOT. Env%depth == INFINITE_DEPTH) THEN
        CALL LISC(omega**2*Env%depth/Env%g, wavenumber*Env%depth,IGreen)
      END IF
      ITERLID=0
      DO I = 1, Mesh%NPanels
        DO J = 1, Mesh%NPanels

          ! First part of the Green function
          ! These output are independent of omega and computed only once in INITIALIZE_GREEN().
            FSP=IGreen%FSP1(I,J)
            FSM=IGreen%FSM1(I,J)
            VSXP=IGreen%VSP1(I,J,:)
            VSXM=IGreen%VSM1(I,J,:)

          ! Second part of the Green function
          IF ((Env%depth == INFINITE_DEPTH).OR.(omega**2*Env%depth/Env%g.GE.20)) THEN
            IF (omega**2*Env%depth/Env%g.GE.20) THEN
            ! First part of the Green function
            FSP=IGreen%FSP1_INF(I,J)
            FSM=IGreen%FSM1_INF(I,J)
            VSXP=IGreen%VSP1_INF(I,J,:)
            VSXM=IGreen%VSM1_INF(I,J,:)
            ENDIF
            CALL VNSINFD                                     &
            (I, wavenumber, Mesh%XM(:, I), J, VFace, Mesh,  &
              IGreen, SP, SM, VSP, VSM                       &
              )
          ELSE
            CALL VNSFD                                       &
            (I, wavenumber, Mesh%XM(:, I), J, VFace,  Mesh, &
                IGreen,Env%depth, SP, SM, VSP, VSM           &
              )
          END IF

          ! Store into influence matrix
          S(I, J, 1) = FSP + SP                              ! Green function
          V(I, J, 1) = DOT_PRODUCT(Mesh%N(:, I), VSXP + VSP) ! Normal Gradient of the Green function

          IF (Mesh%ISym == Y_SYMMETRY) THEN
            S(I, J, 2) = FSM + SM
            V(I, J, 2) = DOT_PRODUCT(Mesh%N(:, I), VSXM + VSM)
          ENDIF

        END DO
      END DO
    END IF
  !=========================
  ! Solve the linear problem
  !=========================
   IF (FLAG_CAL==1 .AND. (SolverOpt%ID .NE. ID_GMRES) ) THEN
      ! Invert matrix V
      IF (SolverOpt%ID== ID_GAUSS) THEN
         CALL GAUSSZ(V(:,:,1),Mesh%NPanels, Vinv(:,:,1))
         IF (Mesh%ISym == Y_SYMMETRY) THEN
              CALL GAUSSZ(V(:,:,2),Mesh%NPanels, Vinv(:,:,2))
         END IF
      ELSE
         CALL LU_INVERS_MATRIX(V(:,:,1),Mesh%NPanels, Vinv(:,:,1))
         IF (Mesh%ISym == Y_SYMMETRY) THEN
              CALL LU_INVERS_MATRIX(V(:,:,2),Mesh%NPanels, Vinv(:,:,2))
         END IF
      END IF
  END IF

  IF (Mesh%ISym == NO_Y_SYMMETRY) THEN
    IF (SolverOpt%ID .EQ. ID_GMRES) THEN
        CALL GMRES_SOLVER(V(:,:,1),NVEL(1:Mesh%NPanels), ZIGB(:), Mesh%NPanels,SolverOpt)
    ELSE
        ZIGB(:) = MATMUL(Vinv(:, :,1), NVEL(1:Mesh%NPanels))
    ENDIF
        ZIGS(:) = CZERO
  ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
    IF (SolverOpt%ID .EQ. ID_GMRES) THEN
        RHS(:)=(NVEL(1:Mesh%NPanels) + NVEL(Mesh%NPanels+1:2*Mesh%NPanels))/2
        IF (ANY(ABS(RHS)>0.)) THEN
           CALL GMRES_SOLVER(V(:,:,1),RHS(:),ZOL(:, 1),Mesh%NPanels,SolverOpt)
        ELSE
           ZOL(:, 1)=0
        ENDIF
        RHS(:)=(NVEL(1:Mesh%NPanels) - NVEL(Mesh%NPanels+1:2*Mesh%NPanels))/2
        IF (ANY(ABS(RHS)>0.)) THEN
           CALL GMRES_SOLVER(V(:,:,2),RHS(:),ZOL(:, 2),  Mesh%NPanels,SolverOpt)
         ELSE
           ZOL(:, 2)=0
        ENDIF
    ELSE
    ZOL(:, 1) = MATMUL(                                              &
      Vinv(:, :,1),                                                  &
      (NVEL(1:Mesh%NPanels) + NVEL(Mesh%NPanels+1:2*Mesh%NPanels))/2 &
      )
    ZOL(:, 2) = MATMUL(                                              &
      Vinv(:, :,2),                                                  &
      (NVEL(1:Mesh%NPanels) - NVEL(Mesh%NPanels+1:2*Mesh%NPanels))/2 &
      )
    ENDIF

    ZIGB(:) = ZOL(:, 1) + ZOL(:, 2)
    ZIGS(:) = ZOL(:, 1) - ZOL(:, 2)
  END IF

  !=============================================================
  ! Computation of potential phi = S*source on the floating body
  !=============================================================

  IF (Mesh%ISym == NO_Y_SYMMETRY) THEN
    Potential(:) = MATMUL(S(:, :, 1), ZIGB(:))

  ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
     Potential(1:Mesh%NPanels) =                    &
      MATMUL(S(:, :, 1) + S(:, :, 2), ZIGB(:))/2   &
      + MATMUL(S(:, :, 1) - S(:, :, 2), ZIGS(:))/2
     Potential(Mesh%NPanels+1:2*Mesh%NPanels) =     &
      MATMUL(S(:, :, 1) - S(:, :, 2), ZIGB(:))/2   &
      + MATMUL(S(:, :, 1) + S(:, :, 2), ZIGS(:))/2
  END IF

  RETURN
  END SUBROUTINE

END MODULE
