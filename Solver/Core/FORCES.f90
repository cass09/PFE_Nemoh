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
!--------------------------------------------------------------------------------------
!
! NEMOH Solver
! Deduce from the potential the forces on the floating body (for the diffraction problems)
! and the added mass and damping (for the radiation problem).
!
!--------------------------------------------------------------------------------------

MODULE FORCES


  USE Constants
  USE MMesh,     ONLY: TMesh

  IMPLICIT NONE

  PUBLIC  :: COMPUTE_AND_WRITE_FORCES
  PRIVATE :: INITIALIZE_FORCE_PARAMETERS, APPEND_FORCE_TO_FILE

  PRIVATE
  INTEGER :: Nintegration
  REAL, DIMENSION(:, :), ALLOCATABLE :: NDS

CONTAINS

  !-------------------------------------------

  SUBROUTINE COMPUTE_AND_WRITE_FORCES  &
      ! Main subroutine. Called from NEMOH.f90.
    ( parameters_file,                 &
      Mesh, rho, omega, Potential,     &
      switch_type, output_file         &
      )

    CHARACTER(LEN=*),                               INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),                                    INTENT(IN) :: Mesh
    REAL,                                           INTENT(IN) :: omega, rho
    INTEGER,                                        INTENT(IN) :: switch_type
    COMPLEX, DIMENSION(Mesh%NPanels*(Mesh%Isym+1)), INTENT(IN) :: Potential

    ! Local variables
    INTEGER :: i, j,indj
    COMPLEX, DIMENSION(:), ALLOCATABLE :: Momentum

    CALL INITIALIZE_FORCE_PARAMETERS(parameters_file, mesh, output_file)

    ALLOCATE(Momentum(Nintegration))
    Momentum(:) = 0.0

    ! Compute force coefficients
    ! momentum int_Sb phi nu dS, with nu=-NDS, NDS Normal towards fluids
    DO i = 1, Nintegration
      DO j = 1, Mesh%NPanels*(Mesh%Isym+1)
        indj=j
        IF (j> Mesh%Npanels) indj=j-Mesh%NPanels
        IF (Mesh%XM(3, indj)<0.) THEN      !dont compute at lid panels
                Momentum(i) = Momentum(i) - RHO * potential(j)*NDS(i, j)
        END IF
      END DO
    !   print*,Momentum(i)
    END DO

    CALL APPEND_FORCE_TO_FILE(output_file, switch_type, omega, Momentum)

    DEALLOCATE(Momentum)

  END SUBROUTINE COMPUTE_AND_WRITE_FORCES

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE INITIALIZE_FORCE_PARAMETERS(parameters_file, mesh, output_file)
    ! Load NDS array from the given file.

    CHARACTER(LEN=*), INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),      INTENT(IN) :: Mesh

    INTEGER :: i, j, u

    IF (.NOT. ALLOCATED(NDS)) THEN
      OPEN(NEWUNIT=u, FILE=parameters_file, STATUS='OLD', ACTION='READ')
      READ(u, *) Nintegration
      ALLOCATE(NDS(Nintegration, Mesh%NPanels*2**Mesh%Isym))
      DO i = 1, Nintegration
        READ(u, *) (NDS(i,j), j=1,Mesh%NPanels*2**Mesh%Isym)
      END DO
      CLOSE(u)

      ! Initialize output file...
      OPEN(NEWUNIT=u, FILE=output_file, ACTION='WRITE')
      WRITE(u, '(A)') ''
      CLOSE(u)
    ELSE
      ! A file has already been loaded.
      ! We assume only one file should be used for each run of the program.
      ! Thus nothing happens.
    END IF

  END SUBROUTINE INITIALIZE_FORCE_PARAMETERS

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE APPEND_FORCE_TO_FILE(filename, switch_type, omega, Momentum)
    ! Append the newly computed forces at the end of each line.
    ! By reading each line and then rewriting it.

    ! Done this way for backward compatibility with Nemoh 2.0.
    ! (Might be more efficient to add each new case on a new line.)

    CHARACTER(LEN=*),                 INTENT(IN) :: filename
    INTEGER,                          INTENT(IN) :: switch_type
    REAL,                             INTENT(IN) :: Omega
    COMPLEX, DIMENSION(Nintegration), INTENT(IN) :: Momentum

    ! Local variables
    INTEGER :: i, u
    REAL    :: line(Nintegration*2)
    COMPLEX :: Fdifrac

    ! Rewrite file with new data appended at the end of the line.
    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE',POSITION='APPEND')
    DO i = 1, Nintegration ! Loop on all lines
      IF (switch_type == DIFFRACTION_PROBLEM) THEN
        Fdifrac    =II*omega*Momentum(i) !-dmomentum/dt
        line(2*i-1)=ABS(Fdifrac)
        line(2*i)  =ATAN2(AIMAG(Fdifrac),REAL(Fdifrac))

      ELSE IF (switch_type == RADIATION_PROBLEM) THEN
        line(2*i-1)=REAL(Momentum(i))      ! added mass
        line(2*i)=omega*AIMAG(Momentum(i)) ! damping
      END IF
    END DO
    WRITE(u,*) (line(i),i=1,2*Nintegration)
    CLOSE(u)

  END SUBROUTINE APPEND_FORCE_TO_FILE

  !-------------------------------------------

END MODULE FORCES
