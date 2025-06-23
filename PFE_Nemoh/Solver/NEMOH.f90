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
!   - A. Babarit
!	  - R. Kurnia
!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - BVP solver - January 2014
!
!--------------------------------------------------------------------------------------

PROGRAM Main

  USE Constants
  USE MMesh,                ONLY: TMesh,           ReadTMesh, DeleteTMesh
  USE MEnvironment,         ONLY: TEnvironment,    ReadTEnvironment
  USE MBodyConditions,      ONLY: TBodyConditions, ReadTBodyConditions
  USE M_Solver,             ONLY: TSolver,         ReadTSolver, ID_GMRES
  USE MLogFile              !
  ! Preprocessing and initialization
  USE MFace,                ONLY: TVFace, Prepare_FaceMesh
  USE M_INITIALIZE_GREEN,   ONLY: TGREEN, INITIALIZE_GREEN
  USE Elementary_functions, ONLY: X0

  ! Resolution
  USE SOLVE_BEM_DIRECT,     ONLY: SOLVE_POTENTIAL_DIRECT
  USE SOLVE_INTERACTION_THEORY,     ONLY: SOLVE_POTENTIAL_MATRIX

  ! Post processing and output
  USE OUTPUT,               ONLY: WRITE_DATA_ON_MESH,WRITE_SOURCES
  USE FORCES,               ONLY: COMPUTE_AND_WRITE_FORCES
  USE KOCHIN,               ONLY: COMPUTE_AND_WRITE_KOCHIN
  USE FREESURFACE,          ONLY: COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION, &
                            READ_CYLSURFACE_PARAMETERS, COMPUTE_AND_WRITE_CYL_SURFACE

  IMPLICIT NONE

  CHARACTER(LEN=1000)   :: wd             ! Working directory path (max length: 1000 characters, increase if necessary)
  TYPE(TMesh)           :: Mesh         ! Mesh of the floating body
  TYPE(TBodyConditions) :: BodyConditions ! Physical conditions on the floating body
  TYPE(TEnvironment)    :: Env            ! Physical conditions of the environment
  TYPE(TSolver)         :: SolverOpt      ! Solver Option, specified by user in input_solver.txt

  INTEGER                            :: i_problem          ! Index of the current problem
  REAL                               :: omega, wavenumber  ! Wave frequency and wavenumber
  COMPLEX, DIMENSION(:), ALLOCATABLE     :: ZIGB, ZIGS     ! Computed source distribution
  COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: V, Vinv,S      ! Influece coefficients
  COMPLEX, DIMENSION(:), ALLOCATABLE :: Potential          ! Computed potential
  REAL,DIMENSION(3)                  :: NOINPUT             ! unused but has to be declared
  TYPE(TVFACE)                       :: VFace              ! Face Mesh structure variable
  TYPE(TGREEN)                       :: IGreen             ! Initial Green variables
  REAL                               :: tcpu_start
  CHARACTER(LEN=1000)                :: LogTextToBeWritten

  ! Initialization ---------------------------------------------------------------------

  WRITE(*,*) ' '
  WRITE(*,'(A,$)') '  -> Initialisation '

  ! Get working directory from command line argument
  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, wd)
  ELSE
    wd = "."
  END IF

  CALL ReadTMesh(Mesh, TRIM(wd)//'/mesh/')
  ALLOCATE(ZIGB(Mesh%NPanels), ZIGS(Mesh%NPanels))
  ALLOCATE(Potential(Mesh%NPanels*2**Mesh%Isym))

  CALL ReadTBodyConditions            &
  ( BodyConditions,                   &
    Mesh%Npanels*2**Mesh%Isym,        &
    TRIM(wd)//'/Normalvelocities.dat' &
    )
  

  CALL ReadTEnvironment(Env, file=TRIM(wd)//'/Nemoh.cal')

  CALL ReadTSolver(SolverOpt,TRIM(wd))

  CALL Prepare_FaceMesh(Mesh,SolverOpt%NP_GQ,VFace)

  CALL INITIALIZE_GREEN(VFace,Mesh,Env%depth,NOINPUT,0,SolverOpt%eps_zmin,IGreen)
  ALLOCATE(S(Mesh%NPanels,Mesh%NPanels,2**Mesh%Isym))
  ALLOCATE(V(Mesh%NPanels,Mesh%NPanels,2**Mesh%Isym))
  ALLOCATE(Vinv(Mesh%NPanels,Mesh%NPanels,2**Mesh%Isym))

  WRITE(*, *) ' '
  WRITE(LogTextToBeWritten,*) 'NP Gauss Quadrature Integ.: ', SolverOpt%NP_GQ
  CALL WRITE_LOGFILE(trim(wd)//'/logfile.txt',TRIM(LogTextToBeWritten),IdStartLog,IdprintTerm)
  WRITE(LogTextToBeWritten,*) 'EPS min z                 : ', SolverOpt%eps_zmin
  CALL WRITE_LOGFILE(trim(wd)//'/logfile.txt',TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

  WRITE(*, *) '. Done !'
  WRITE(*, *) ' '

  ! Solve BVPs and calculate forces ----------------------------------------------------
  WRITE(*, *) ' -> Solve BVPs and calculate forces '
  WRITE(LogTextToBeWritten,*) 'Linear Solver: ', SolverOpt%SNAME
  CALL WRITE_LOGFILE(trim(wd)//'/logfile.txt',TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
  CALL START_RECORD_TIME(tcpu_start,trim(wd)//'/logfile.txt',IdAppend)
  WRITE(*, *) ' '

  DO i_problem = 1, BodyConditions%Nproblems
    WRITE(*,'(A,I5,A,I5,A,A,$)') ' Problem ',i_problem,' / ',BodyConditions%Nproblems,' ',CHAR(13)
    omega = BodyConditions%omega(i_problem) ! Wave frequency
    ! Compute wave number k
    IF ((Env%depth == INFINITE_DEPTH) .OR. (omega**2*Env%depth/Env%g >= 20)) THEN
      wavenumber = omega**2/Env%g
    ELSE
      wavenumber = X0(omega**2*Env%depth/Env%g)/Env%depth
      ! X0(y) returns the solution of y = x * tanh(x)
    END IF
    !===============
    ! BEM Resolution
    !===============
      CALL SOLVE_POTENTIAL_DIRECT                                              &
      !==========================
      ( VFace, Mesh, Env, omega, wavenumber,IGreen,                            &
        BodyConditions%NormalVelocity(1:Mesh%Npanels*2**Mesh%Isym, i_problem), &
        S,V,Vinv,ZIGB, ZIGS,                                                   &
        Potential(:),SolverOpt,trim(wd))
    ! WRITE(*,*) Potential(1:10)
    !===========================
    ! Post processing and output
    !===========================
    CALL COMPUTE_AND_WRITE_FORCES            &
    !============================
    ( TRIM(wd)//'/mesh/Integration.dat',     &
      Mesh, Env%rho, omega, Potential,       &
      Bodyconditions%Switch_type(i_problem), &
      TRIM(wd)//'/results/Forces.dat'        &
      )
    IF (BodyConditions%Switch_Potential(i_problem) == 1) THEN
      ! Write pressure field on the floating body in file
      CALL WRITE_DATA_ON_MESH                                      &
      !=======================
      ( Mesh,                                                      &
        Env%rho*II*omega*Potential(:),                             &
        TRIM(wd)//'/results/pressure.'//string(i_problem)//'.dat'  &
        )
    END IF

    IF (BodyConditions%Switch_Kochin(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_KOCHIN                             &
      !============================
      ( TRIM(wd)//'/mesh/Kochin.dat',                           &
        Mesh, Env, wavenumber, ZIGB, ZIGS,                      &
        TRIM(wd)//'/results/Kochin.'//string(i_problem)//'.dat' &
        )
    END IF
    IF (BodyConditions%Switch_FreeSurface(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION                  &
      !============================================
      ( TRIM(wd)//'/mesh/Freesurface.dat', IGreen,VFace,             &
        Mesh, Env, omega, wavenumber, ZIGB, ZIGS,                    &
        TRIM(wd)//'/results/freesurface.'//string(i_problem)//'.dat' &
        )
    END IF

    IF (BodyConditions%Switch_SourceDistr(i_problem) == 1) THEN
      ! Write pressure field on the floating body in file
      CALL WRITE_SOURCES(ZIGB,ZIGS,Mesh%Npanels,                    &
       TRIM(wd)//'/results/sources/sources.'//string(i_problem)//'.dat')
    END IF
    IF (BodyConditions%ITsources(i_problem) == 1) THEN
      CALL WRITE_SOURCES(ZIGB,ZIGS,Mesh%Npanels,                    &
       TRIM(wd)//'/results/sources.'//string(i_problem)//'.dat')
    END IF
    IF (BodyConditions%Switch_Cylsurface(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_CYL_SURFACE                 &
      !============================================
      ( TRIM(wd)//'/mesh/Cylsurface.dat', IGreen,VFace,             &
        Mesh, Env, omega, wavenumber, ZIGB, ZIGS,                    &
        TRIM(wd)//'/results/cylsurface.'//string(i_problem)//'.dat')
    END IF
  END DO

  ! IF (BodyConditions%run_IT== 1) THEN
  !   CALL SOLVE_POTENTIAL_MATRIX(Env, SolverOpt, wd)
  ! END IF

  CALL END_RECORD_TIME(tcpu_start,trim(wd)//'/logfile.txt')
  WRITE(*,*) '. Done !'
  ! Finalize ---------------------------------------------------------------------------
  DEALLOCATE(ZIGB, ZIGS, Potential,S,V,Vinv)
  DEALLOCATE(IGreen%FSP1,IGreen%FSM1,IGreen%VSP1,IGREEN%VSM1)
  DEALLOCATE(IGreen%FSP1_INF,IGreen%FSM1_INF,IGreen%VSP1_INF,IGREEN%VSM1_INF)
  DEALLOCATE(IGreen%XR,IGreen%XZ)
  DEALLOCATE(IGreen%APD1X,IGREEN%APD2X,IGREEN%APD1Z,IGREEN%APD2Z)

CONTAINS

  FUNCTION string (i) result (s)
    ! For example 5 -> "00005"
    INTEGER :: i
    CHARACTER(LEN=5) :: s
    WRITE(s, '(I0.5)') i
  END FUNCTION

END PROGRAM Main
