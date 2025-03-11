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
! Compute the first part (independant of frequency) of the Green function.
!
!--------------------------------------------------------------------------------------
MODULE GREEN_1

  USE Constants
  USE MMesh
  USE MFace,                ONLY: TFace,TVFace,VFace_to_Face
  USE Elementary_functions, ONLY: CROSS_PRODUCT

  IMPLICIT NONE

  PUBLIC :: VAV, COMPUTE_S0, COMPUTE_ASYMPTOTIC_S0

CONTAINS

  SUBROUTINE VAV              &
    ( I, X0I, J, VFace,Mesh, depth,eps_zmin, &
      FSP, FSM, VSP, VSM)
    ! Main subroutine of the module, called in SOLVE_BEM.f90 and FREESURFACE.f90.

    ! Inputs
    INTEGER, INTENT(IN)            :: I    ! Index of the computed flow field  panel.
    REAL, DIMENSION(3), INTENT(IN) :: X0I  ! Coordinates of the computed flow field point
    INTEGER, INTENT(IN)            :: J    ! Index of the source integration panel.
    TYPE(TMesh) , INTENT(IN)       :: Mesh
    TYPE(TVFace), INTENT(IN)       :: VFace
    REAL, INTENT(IN)               :: depth
    REAL, INTENT(IN)               :: eps_zmin

    ! Outputs
    REAL, INTENT(OUT)               :: FSP, FSM ! Integral of the Green function over the panel.
    REAL, DIMENSION(3), INTENT(OUT) :: VSP, VSM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    TYPE(TFace)         :: Face
    REAL, DIMENSION(3)  :: XI
    INTEGER             :: MK
    REAL                :: S0, S1
    REAL, DIMENSION(3)  :: VS0, VS1

    IF (depth == INFINITE_DEPTH) THEN
      ! In case of infinite depth, the reflection of the problem across the free surface will be computed.
      MK = -1
      ! The contribution of the reflected problem is substracted to the solution of the original problem to ensure that φ=0 at the
      ! free surface. The full boundary condition at the interface will be obtained later in the other part of the Green function.
    ELSE
      ! In case of finite depth, the reflection of the problem around the sea bottom will be computed.
      MK = 1
      ! The contribution of the reflected problem is added to the solution of the original problem to ensure
      ! that (∇φ)_z=0 at the sea bottom.
    END IF


    CALL VFace_to_Face(VFace,Face,J) !Extract a face J from the VFace array

    XI(:) = X0I(:)
    IF (I>0) XI(3) = MIN(X0I(3), -eps_zmin*Mesh%xy_diameter) ! on body panels
                                                             !I=0 on free surface
    CALL COMPUTE_S0(XI, Face, S0, VS0)
    IF (I == J) THEN
      IF (X0I(3) >= -eps_zmin*Mesh%xy_diameter) THEN
         VS0(:) = VS0(:) + 2*PI*Face%N(:)
      ELSE
        ! print*,Face%N(1),Face%N(2),Face%N(3)
        ! print*,VS0(1),VS0(2),VS0(3)
         VS0(:) = VS0(:) -2*PI*Face%N(:)
      ENDIF
    END IF

    ! Reflected problem across the free surface/sea bottom.
    XI(:) = X0I(:)
    XI(3) = -X0I(3) - 2*depth
    CALL COMPUTE_S0(XI, Face, S1, VS1)
    VS1(3) = -VS1(3)

    IF (Mesh%Isym == NO_Y_SYMMETRY) THEN
      ! Add up the contributions of the two problems.
      FSP    = -S0     - MK*S1
      VSP(:) = -VS0(:) - MK*VS1(:)
      FSM    = ZERO
      VSM(:) = ZERO

    ELSE IF (Mesh%Isym == Y_SYMMETRY) THEN
      ! Add up the contributions of the two problems...
      FSP    = -S0     - MK*S1
      VSP(:) = -VS0(:) - MK*VS1(:)
      FSM    = -S0     - MK*S1
      VSM(:) = -VS0(:) - MK*VS1(:)

      !... and do some more.

      ! Reflected problem across the symmetry plane (xOz)

      XI(:) = X0I(:)
      XI(2) = -X0I(2)
      IF (I>0) XI(3) = MIN(X0I(3), -eps_zmin*Mesh%xy_diameter) ! on body panels

      CALL COMPUTE_S0(XI, Face, S0, VS0)
      VS0(2) = -VS0(2)

      ! Reflected problem across the symmetry plane (xOz) and across the free surface/sea bottom.
      XI(1) = X0I(1)
      XI(2) = -X0I(2)
      XI(3) = -X0I(3) - 2*depth
      CALL COMPUTE_S0(XI, Face, S1, VS1)
      VS1(2:3) = -VS1(2:3)

      ! Add up the new results.
      FSP    = FSP    - S0     - MK*S1
      VSP(:) = VSP(:) - VS0(:) - MK*VS1(:)
      FSM    = FSM    + S0     + MK*S1
      VSM(:) = VSM(:) + VS0(:) + MK*VS1(:)

    END IF

    FSP    = FSP/(4*PI)
    FSM    = FSM/(4*PI)
    VSP(:) = VSP(:)/(4*PI)
    VSM(:) = VSM(:)/(4*PI)

    RETURN
  END SUBROUTINE VAV

  !--------------------------------------------------

  SUBROUTINE COMPUTE_S0(M, Face, S0, VS0)
    ! Estimate the integral over the face S0 = ∫∫ 1/MM' dS(M')
    ! and its derivative with respect to M.

    ! Based on formulas A6.1 and A6.3 (p. 381 to 383)
    ! in G. Delhommeau thesis (referenced below as [Del]).

    ! Inputs
    REAL, DIMENSION(3), INTENT(IN) :: M
    TYPE(TFace),        INTENT(IN) :: Face

    ! Outputs
    REAL,               INTENT(OUT) :: S0
    REAL, DIMENSION(3), INTENT(OUT) :: VS0

    ! Local variables
    INTEGER               :: L
    REAL                  :: RO, GZ, DK, GY
    REAL, DIMENSION(5)    :: RR
    REAL, DIMENSION(3, 5) :: DRX
    REAL                  :: ANT, DNT, ANL, DNL, ALDEN, AT
    REAL, DIMENSION(3)    :: PJ, GYX, ANTX, ANLX, DNTX

    RO = NORM2(M(1:3) - Face%XM(1:3)) ! Distance from center of mass of the face to M.

    IF (RO > 7*Face%tDis) THEN
      ! Asymptotic value if face far away from M
      S0       = Face%A/RO
      VS0(1:3) = (Face%XM(1:3) - M)*S0/RO**2

    ELSE

      GZ = DOT_PRODUCT(M(1:3) - Face%XM(1:3), Face%N(1:3)) ! Called Z in [Del]

      DO L = 1, 5
        RR(L) = NORM2(M(1:3) - Face%X(1:3, L))       ! Distance from vertices of Face to M.
        DRX(:, L) = (M(1:3) - Face%X(1:3, L))/RR(L)  ! Normed vector from vertices of Face to M.
      END DO

      S0 = 0.0
      VS0(:) = 0.0

      DO L = 1, 4
        DK = NORM2(Face%X(:, L+1) - Face%X(:, L))    ! Distance between two consecutive points, called d_k in [Del]
        IF (DK >= 1E-3*Face%tDis) THEN
          PJ(:) = (Face%X(:, L+1) - Face%X(:, L))/DK ! Normed vector from one corner to the next
          GYX = CROSS_PRODUCT(Face%N, PJ)            ! Called (a,b,c) in [Del]
          GY = DOT_PRODUCT(M - Face%X(:, L), GYX)    ! Called Y_k in  [Del]

          ANT = 2*GY*DK                                                  ! Called N^t_k in [Del]
          DNT = (RR(L+1)+RR(L))**2 - DK*DK + 2.0*ABS(GZ)*(RR(L+1)+RR(L)) ! Called D^t_k in [Del]
          ANL = RR(L+1) + RR(L) + DK                                     ! Called N^l_k in [Del]
          DNL = RR(L+1) + RR(L) - DK                                     ! Called D^l_k in [Del]
          ALDEN = ALOG(ANL/DNL)

          IF (ABS(GZ) >= 1.E-4*Face%tDis) THEN
            AT = ATAN(ANT/DNT)
          ELSE
            AT = 0.
          ENDIF

          ANLX(:) = DRX(:, L+1) + DRX(:, L)                    ! Called N^l_k_{x,y,z} in [Del]

          ANTX(:) = 2*DK*GYX(:)                                ! Called N^t_k_{x,y,z} in [Del]
          DNTX(:) = 2*(RR(L+1) + RR(L) + ABS(GZ))*ANLX(:) &
            + 2*SIGN(1.0, GZ)*(RR(L+1) + RR(L))*Face%N(:)      ! Called D^t_k_{x,y,z} in [Del]

          S0 = S0 + GY*ALDEN - 2*AT*ABS(GZ)

          VS0(:) = VS0(:) + ALDEN*GYX(:)     &
            - 2*SIGN(1.0, GZ)*AT*Face%N(:)   &
            + GY*(DNL-ANL)/(ANL*DNL)*ANLX(:) &
            - 2*ABS(GZ)*(ANTX(:)*DNT - DNTX(:)*ANT)/(ANT*ANT+DNT*DNT)
        END IF
      END DO
    END IF

  END SUBROUTINE COMPUTE_S0

  !--------------------------------------------------

  SUBROUTINE COMPUTE_ASYMPTOTIC_S0(XI, XJ, AJ,EPS, S0, VS0)
    ! Same as above, but always use the approximate aymptotic value.

    ! Inputs
    REAL, DIMENSION(3), INTENT(IN) :: XI,XJ
    REAL                           :: AJ !area of panel J
    REAL,               INTENT(IN) :: EPS
    ! Outputs
    REAL,               INTENT(OUT) :: S0
    REAL, DIMENSION(3), INTENT(OUT) :: VS0

    ! Local variables
    REAL :: RO

    RO = NORM2(XI(1:3) - XJ(1:3)) ! Distance from XI to XJ

    IF (RO > EPS) THEN
      ! Face far away from M
      S0       = AJ/RO
      VS0(1:3) = (XJ(1:3) - XI)*S0/RO**2
    ELSE
      S0 = 0.0
      VS0(1:3) = 0.0
    END IF

    RETURN
  END SUBROUTINE COMPUTE_ASYMPTOTIC_S0

  !--------------------------------------------------

END MODULE GREEN_1
