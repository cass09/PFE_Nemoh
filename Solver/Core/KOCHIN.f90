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
!   - R. Kurnia
!--------------------------------------------------------------------------------------
MODULE KOCHIN

  USE Constants
  USE MMesh,        ONLY: TMesh
  USE MEnvironment, ONLY: TEnvironment
  USE Elementary_functions,       ONLY: DOT_PRODUCT_COMPLEX
  
  IMPLICIT NONE

  PUBLIC  :: COMPUTE_AND_WRITE_KOCHIN
  PRIVATE :: READ_KOCHIN_PARAMETERS, COMPUTE_ZS, WRITE_KOCHIN

  PRIVATE
  INTEGER                         :: Ntheta
  REAL, DIMENSION(:), ALLOCATABLE :: Thetas

CONTAINS

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE COMPUTE_AND_WRITE_KOCHIN  &
    ( parameters_file,                 &
      Mesh, Env, wavenumber, ZIGB, ZIGS,  &
      output_file)

    CHARACTER(LEN=*),                 INTENT(IN) :: parameters_file, output_file
    TYPE(TMesh),                      INTENT(IN) :: Mesh
    TYPE(TEnvironment),               INTENT(IN) :: Env
    REAL,                             INTENT(IN) :: wavenumber
    COMPLEX, DIMENSION(Mesh%NPanels), INTENT(IN) :: ZIGB, ZIGS

    ! Local variables
    INTEGER                            :: i
    COMPLEX, DIMENSION(Mesh%NPanels)   :: ZS
    COMPLEX, DIMENSION(:), ALLOCATABLE :: HKochin

    CALL READ_KOCHIN_PARAMETERS(parameters_file)

    ALLOCATE(HKochin(Ntheta))
    DO i = 1, NTheta
      CALL COMPUTE_ZS(Mesh, wavenumber, Env, Thetas(i), 1, ZS)
      Hkochin(i) = DOT_PRODUCT_COMPLEX(ZIGB, ZS,Mesh%NPanels)

      IF (Mesh%ISym == Y_SYMMETRY) THEN
        CALL COMPUTE_ZS(Mesh, wavenumber, Env, Thetas(i), -1, ZS)
        Hkochin(i) = Hkochin(i) + DOT_PRODUCT_COMPLEX(ZIGS, ZS,Mesh%NPanels)
      END IF

      HKochin(i) = HKochin(i)/(4*PI)
    END DO

    CALL WRITE_KOCHIN(output_file, HKochin)
    DEALLOCATE(HKochin)

  END SUBROUTINE COMPUTE_AND_WRITE_KOCHIN

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE READ_KOCHIN_PARAMETERS(filename)
    ! Load data from parameter files into array Thetas.

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: i, u

    IF (.NOT. ALLOCATED(Thetas)) THEN
      OPEN(NEWUNIT=u, FILE=filename, STATUS='OLD', ACTION='READ')
      READ(u, *) Ntheta
      ALLOCATE(Thetas(Ntheta))
      DO i = 1, Ntheta
        READ(u, *) Thetas(i)
      END DO
      CLOSE(u)
    ELSE
      ! A file has already been loaded.
      ! We assume only one file should be used for each run of the program.
      ! Thus nothing happens.
    END IF

  END SUBROUTINE READ_KOCHIN_PARAMETERS

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE COMPUTE_ZS(Mesh, wavenumber, Env, theta, mSym, ZS)

    ! Inputs
    TYPE(TMesh),        INTENT(IN) :: Mesh
    REAL,               INTENT(IN) :: wavenumber, theta
    TYPE(TEnvironment), INTENT(IN) :: Env
    INTEGER,            INTENT(IN) :: mSym ! for y-symetry trick

    ! Output
    COMPLEX, DIMENSION(Mesh%NPanels), INTENT(OUT) :: ZS

    ! Locals
    INTEGER                          :: I, K, L
    COMPLEX, DIMENSION(Mesh%NPoints) :: CEP, CEM, ZJ
    INTEGER, DIMENSION(5)            :: KK
    COMPLEX                          :: ZMS, ZNS, ZCS, ZC, ZHS, ZHD

    DO k = 1, Mesh%Npoints
      ZJ(k) = wavenumber*(                         &
        Mesh%X(3, k)                               &
        - II*(     Mesh%X(1, k) - Env%Xeff)*COS(Theta) &
        - II*(mSym*Mesh%X(2, k) - Env%Yeff)*SIN(Theta) &
        )

      IF (REAL(ZJ(k)) < -18.0) THEN
        CEP(k) = CZERO
      ELSE
        CEP(k) = CEXP(ZJ(k))
      END IF

      IF ((REAL(-ZJ(k)-2*wavenumber*Env%depth) < -18.0) .OR. (-wavenumber*Env%depth >= 0.0)) THEN
        CEM(k) = CZERO
      ELSE
        CEM(k) = CEXP(-ZJ(k)-2*wavenumber*Env%depth)
      END IF
    END DO

    DO i = 1, Mesh%NPanels
      KK(1:4) = Mesh%P(1:4, i)
      KK(5) = KK(1)

      ZNS =      Mesh%N(1, i) - II*Mesh%N(3, i)*COS(Theta)
      ZMS = mSym*Mesh%N(2, i) - II*Mesh%N(3, i)*SIN(Theta)

      ZS(i) = CZERO
      DO L = 1, 4

        ZCS =       ZMS*(Mesh%X(1, KK(L+1)) - Mesh%X(1, KK(L))) &
              -mSym*ZNS*(Mesh%X(2, KK(L+1)) - Mesh%X(2, KK(L)))
        ZC = ZJ(KK(L+1)) - ZJ(KK(L))

        IF ((ABS(AIMAG(ZC)) < 1.E-04) .AND. (ABS(REAL(ZC)) < 1.E-04)) THEN
          ZHS = 0.5*(CEP(KK(L+1)) + CEP(KK(L)))
          ZHD = -0.5*(CEM(KK(L+1)) + CEM(KK(L)))
        ELSE
          ZHS = (CEP(KK(L+1)) - CEP(KK(L)))/ZC
          ZHD = (CEM(KK(L+1)) - CEM(KK(L)))/ZC
        ENDIF

        ZS(i) = ZS(i) + ZCS*ZHS + CONJG(ZCS*ZHD)

      END DO

      IF ((wavenumber*Env%depth > 18.0) .OR. (wavenumber*Env%depth <= 0.0)) THEN
        ZS(i) = mSym*ZS(i)/wavenumber
      ELSE
        ZS(i) = mSym*0.5*ZS(i)*EXP(wavenumber*Env%depth)/(SINH(wavenumber*Env%depth)*wavenumber)
      END IF
    END DO

  END SUBROUTINE COMPUTE_ZS

  !-------------------------------------------
  !-------------------------------------------
  !-------------------------------------------

  SUBROUTINE WRITE_KOCHIN(filename, HKochin)

    CHARACTER(LEN=*),      INTENT(IN) :: filename
    COMPLEX, DIMENSION(*), INTENT(IN) :: HKochin

    ! Local variables
    INTEGER          :: i, u

    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE')

    DO i = 1, Ntheta
      WRITE(u, *) Thetas(i), ABS(HKochin(i)), ATAN2(AIMAG(HKochin(i)),REAL(HKochin(i)))
    END DO

    CLOSE(u)

  END SUBROUTINE WRITE_KOCHIN

  !-------------------------------------------

END MODULE KOCHIN
