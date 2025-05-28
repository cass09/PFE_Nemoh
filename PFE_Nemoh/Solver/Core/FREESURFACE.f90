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
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh
!
!--------------------------------------------------------------------------------------
MODULE FREESURFACE
  ! Computation of the free surface elevation from the computed source distribution.

  USE Constants
  USE MMesh,                     ONLY: TMesh, CreateTMesh, DeleteTMesh
  USE MEnvironment,              ONLY: TEnvironment
  USE MFace,                     ONLY: TVFace
  ! Green functions
  USE M_INITIALIZE_GREEN,        ONLY: TGreen
  USE GREEN_1,                   ONLY: VAV
  USE GREEN_2,                   ONLY: VNSFD, VNSINFD

  USE OUTPUT

  IMPLICIT NONE
  ! CML - MODIFS for the Interaction Theory : 
  ! - reading of the cylindrical mesh
  ! - computation of the potential on this mesh

  PUBLIC  :: COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION, COMPUTE_AND_WRITE_CYL_SURFACE, &
  READ_CYLSURFACE_PARAMETERS
  PRIVATE :: READ_FREESURFACE_PARAMETERS, COMPUTE_POTENTIAL_AT_POINT, WRITE_FS

  ! PRIVATE
  ! TYPE(TMesh) :: MeshFS ! Mesh of the free surface

CONTAINS

  !-------------------------------------------

  SUBROUTINE COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION  &
    ! Main subroutine of the module. Called in NEMOH.f90.
    ( parameters_file,  IGreen, VFace,                 &
      Mesh, Env, omega, wavenumber, ZIGB, ZIGS,        &
      output_file                                      &
      )

    CHARACTER(LEN=*),                 INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),                      INTENT(IN) :: Mesh
    TYPE(TVFace),                     INTENT(IN) :: VFace
    TYPE(TEnvironment),               INTENT(IN) :: Env
    REAL,                             INTENT(IN) :: omega, wavenumber
    COMPLEX, DIMENSION(Mesh%NPanels), INTENT(IN) :: ZIGB, ZIGS ! Sources
    TYPE(TGREEN),                     INTENT(IN) :: IGreen

    ! Local variables
    INTEGER :: j
    COMPLEX :: PHI
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ETA
    TYPE(TMesh) :: MeshFS

    CALL READ_FREESURFACE_PARAMETERS(parameters_file, MeshFS)

    ALLOCATE(ETA(MeshFS%NPoints))

    DO j = 1, MeshFS%Npoints

      CALL COMPUTE_POTENTIAL_AT_POINT             &
      !==============================
      ( Mesh, Env, omega, wavenumber, ZIGB, ZIGS, &
        MeshFS%X(1, j), MeshFS%X(2, j), 0.0,      &
        PHI, IGreen,VFace                         &
        )

      ! Get elevation ETA from potential PHI
      ETA(j) = II*omega/Env%G*PHI
    END DO

    CALL WRITE_FS(output_file, ETA, MeshFS)

    IF (MeshFS%Npoints.GT.0) CALL DeleteTMesh(MeshFS)


    DEALLOCATE(ETA)

  END SUBROUTINE COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION

  !-------------------------------------------
  SUBROUTINE COMPUTE_AND_WRITE_CYL_SURFACE  &
    ! Main subroutine of the module. Called in NEMOH.f90.
    ( parameters_file, IGreen, VFace,                 &
      Mesh, Env, omega, wavenumber, ZIGB, ZIGS,        &
      output_file               &
      )

    CHARACTER(LEN=*),                 INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),                      INTENT(IN) :: Mesh
    TYPE(TVFace),                     INTENT(IN) :: VFace
    TYPE(TEnvironment),               INTENT(IN) :: Env
    REAL,                             INTENT(IN) :: omega, wavenumber
    COMPLEX, DIMENSION(Mesh%NPanels), INTENT(IN) :: ZIGB, ZIGS ! Sources
    TYPE(TGREEN),                     INTENT(IN) :: IGreen


    ! Local variables
    TYPE(TMesh) :: MeshCyl
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ETA
    INTEGER :: j
    COMPLEX :: PHI

    CALL READ_CYLSURFACE_PARAMETERS(parameters_file, MeshCyl)
    ALLOCATE(ETA(MeshCyl%Npoints))

    DO j = 1, MeshCyl%Npoints
      CALL COMPUTE_POTENTIAL_AT_POINT             &
      !==============================
      ( Mesh, Env, omega, wavenumber, ZIGB, ZIGS, &
        MeshCyl%X(1, j), MeshCyl%X(2, j), MeshCyl%X(3, j),      &
        PHI, IGreen,VFace                         &
        )
      ! Mesh OR MeshCyl ???????! 

      ! Get elevation ETA from potential PHI
      ETA(j) = II*omega/Env%G*PHI
    END DO
    CALL WRITE_FS(output_file, ETA, MeshCyl)

    IF (MeshCyl%Npoints.GT.0) CALL DeleteTMesh(MeshCyl)

    ! DEALLOCATE(ETA)

  END SUBROUTINE COMPUTE_AND_WRITE_CYL_SURFACE
  !-------------------------------------------

  SUBROUTINE READ_FREESURFACE_PARAMETERS(filename, MeshFS)
    ! Load mesh of the free surface

    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(TMesh), INTENT(INOUT)   :: MeshFS

    INTEGER :: j, u

    IF (.NOT. ALLOCATED(MeshFS%X)) THEN
      OPEN(NEWUNIT=u, FILE=filename, STATUS='OLD', ACTION='READ')
      READ(u, *) MeshFS%Npoints, MeshFS%Npanels
      IF (MeshFS%Npoints > 0) THEN
        CALL CreateTMesh(MeshFS, MeshFS%Npoints, MeshFS%Npanels, 1)
        DO j = 1, MeshFS%Npoints
          READ(u, *) MeshFS%X(1,j), MeshFS%X(2,j)
        END DO
        DO j = 1, MeshFS%Npanels
          READ(u, *) MeshFS%P(1,j), MeshFS%P(2,j), MeshFS%P(3,j), MeshFS%P(4,j)
        END DO
      END IF
      CLOSE(u)
    ELSE
      ! A file has already been loaded.
      ! We assume only one file should be used for each run of the program.
      ! Thus nothing happens.
    END IF

  END SUBROUTINE READ_FREESURFACE_PARAMETERS

  !-------------------------------------------
  SUBROUTINE READ_CYLSURFACE_PARAMETERS(filename, MeshCyl)
    ! Load mesh of the cylindrical surface

    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(TMesh), INTENT(INOUT)   :: MeshCyl

    INTEGER :: j, u

    IF (.NOT. ALLOCATED(MeshCyl%X)) THEN
      OPEN(NEWUNIT=u, FILE=filename, STATUS='OLD', ACTION='READ')
      READ(u, *) MeshCyl%Npoints, MeshCyl%Npanels
      IF (MeshCyl%Npoints > 0) THEN
        CALL CreateTMesh(MeshCyl, MeshCyl%Npoints, MeshCyl%Npanels, 1)
        DO j = 1, MeshCyl%Npoints
          READ(u, *) MeshCyl%X(1,j), MeshCyl%X(2,j), MeshCyl%X(3,j)
        END DO
        DO j = 1, MeshCyl%Npanels
          READ(u, *) MeshCyl%P(1,j), MeshCyl%P(2,j), MeshCyl%P(3,j), MeshCyl%P(4,j)
        END DO
      END IF
      CLOSE(u)
    ELSE
      ! A file has already been loaded.
      ! We assume only one file should be used for each run of the program.
      ! Thus nothing happens.
    END IF

  END SUBROUTINE READ_CYLSURFACE_PARAMETERS
  !-------------------------------------------

  SUBROUTINE COMPUTE_POTENTIAL_AT_POINT         &
    ( Mesh, Env, omega, wavenumber, ZIGB, ZIGS, &
      XC, YC, ZC,                               &
      PHI,IGreen,VFace                          &
      )
    ! Compute the potential PHI at the given coordinate (XC, YC, ZC) using the
    ! source distributions ZIGB abd ZIGS.

    ! Inputs
    TYPE(TMesh),                      INTENT(IN) :: Mesh ! Mesh of the floating body
    TYPE(TVFace),                     INTENT(IN) :: VFace!
    TYPE(TEnvironment),               INTENT(IN) :: Env
    REAL,                             INTENT(IN) :: omega, wavenumber
    COMPLEX, DIMENSION(Mesh%NPanels), INTENT(IN) :: ZIGB, ZIGS
    REAL,                             INTENT(IN) :: XC, YC, ZC
    TYPE(TGREEN),                     INTENT(IN) :: IGreen

    ! Output
    COMPLEX,                          INTENT(OUT) :: PHI

    ! Local variables
    INTEGER               :: J
    REAL                  :: FSP, FSM
    REAL, DIMENSION(3)    :: VSXP, VSXM
    COMPLEX               :: SP, SM
    COMPLEX, DIMENSION(3) :: VSP, VSM

    PHI = 0.0
    DO J = 1, Mesh%NPanels

      ! Compute the Greem function coefficients between the point and the face of index J.
      ! Compute also its gradient although it is not used.
      ! Values could be stored for more efficiency as in SOLVE_BEM.
      call vav(0, (/xc, yc, zc/), j,VFace, mesh, env%depth,IGREEN%EPS_ZMIN, fsp, fsm, vsxp, vsxm)
      IF ((Env%depth == INFINITE_DEPTH) .OR. (wavenumber*Env%depth >= 20)) THEN
        CALL VNSINFD(0,wavenumber, (/XC, YC, ZC/), J, VFace, Mesh, IGreen, SP, SM, VSP, VSM)
      ELSE
        CALL VNSFD  (0,wavenumber, (/XC, YC, ZC/), J, VFace, Mesh, IGreen, Env%depth, SP, SM, VSP, VSM)
      ENDIF
      
      ! Compute potential from sources and Green function.
      IF (Mesh%Isym == NO_Y_SYMMETRY) THEN
        PHI = PHI + (SP+FSP)*ZIGB(J)
      ELSE IF (Mesh%Isym == Y_SYMMETRY) THEN
        PHI = PHI + (FSP+SP+FSM+SM)/2*ZIGB(J) + (FSP+SP-FSM-SM)/2*ZIGS(J)
      ENDIF
    END DO

    RETURN
  END SUBROUTINE COMPUTE_POTENTIAL_AT_POINT

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE WRITE_FS(filename, ETA, Mesh)
    ! Write the field ETA into a text file.

    CHARACTER(LEN=*),      INTENT(IN) :: filename
    COMPLEX, DIMENSION(*), INTENT(IN) :: ETA
    TYPE(TMesh),           INTENT(IN) :: Mesh
    ! Local variables
    INTEGER          :: i, u

    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE')

    IF (output_format == TECPLOT_OUTPUT) THEN
      WRITE(u, *) 'VARIABLES="X" "Y" "Z" "abs(eta) (m)" "angle(phi) (rad)" "PRE1" "PRE2"'
      WRITE(u, '(A,I7,A,I7,A)') 'ZONE N=', Mesh%Npoints, ' , E = ', Mesh%Npanels, ' , F=FEPOINT,ET=QUADRILATERAL'
    END IF

    DO i = 1, Mesh%Npoints
      WRITE(u, '(7(X, E14.7))') Mesh%X(1, i), Mesh%X(2, i), Mesh%X(3, i), ABS(eta(i)), ATAN2(IMAG(eta(i)), REAL(eta(I))), REAL(eta(i)), IMAG(eta(i))
    END DO

    IF (output_format == TECPLOT_OUTPUT) THEN
      DO i=1,Mesh%Npanels
        WRITE(u, *) Mesh%P(1, i), Mesh%P(2, i), Mesh%P(3, i), Mesh%P(4, i)
      END DO
    END IF

    CLOSE(u)

  END SUBROUTINE WRITE_FS

  !-------------------------------------------

END MODULE FREESURFACE
