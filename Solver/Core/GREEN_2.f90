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
!
!--------------------------------------------------------------------------------------
MODULE GREEN_2

  USE Constants
  USE MMesh
  USE MFace,              ONLY:TFace,TVFace,VFace_to_FACE
  USE Elementary_functions

  USE M_INITIALIZE_GREEN, ONLY: TGreen,FLAG_IGREEN
  USE GREEN_1,            ONLY: COMPUTE_ASYMPTOTIC_S0

  IMPLICIT NONE

  PUBLIC  :: VNSINFD, VNSFD
  PRIVATE :: COMPUTE_S2

CONTAINS

  !-------------------------------------------------------------------------------!

  SUBROUTINE VNSINFD                   &
      (I, wavenumber, X0I, J, VFace,Mesh, &
      IGreen,SP, SM, VSP, VSM )
    ! Compute the frequency-dependent part of the Green function in the infinite depth case.

    ! Inputs
    REAL,                  INTENT(IN)  :: wavenumber
    REAL, DIMENSION(3),    INTENT(IN)  :: X0I   ! Coordinates of the computed flow field point
    INTEGER,               INTENT(IN)  :: I     ! Index of the flow integration panel
    INTEGER,               INTENT(IN)  :: J     ! Index of the source integration panel
    TYPE(TMesh),           INTENT(IN)  :: Mesh
    TYPE(TGreen),          INTENT(IN)  :: IGreen ! Initial green variable
    TYPE(TVFace),          INTENT(IN)  :: VFace

    ! Outputs
    COMPLEX,               INTENT(OUT) :: SP, SM   ! Integral of the Green function over the panel.
    COMPLEX, DIMENSION(3), INTENT(OUT) :: VSP, VSM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    REAL                               :: EPS
    REAL                               :: ADPI, ADPI2, AKDPI, AKDPI2
    REAL, DIMENSION(3)                 :: XI,XJ
    COMPLEX, DIMENSION(Mesh%ISym+1)    :: FS
    COMPLEX, DIMENSION(3, Mesh%ISym+1) :: VS
    TYPE(TFace)                        :: FaceJ
    INTEGER                            :: IGQ !Index Gauss Quadrature Integration point
    COMPLEX               :: SP_IGQ, SM_IGQ   ! Integral of the Green function over the panel for Gauss point IGQ.
    COMPLEX, DIMENSION(3) :: VSP_IGQ, VSM_IGQ ! Gradient of the integral of the Green function

    ALLOCATE(FaceJ%dXdXG_WGQ_per_A(VFace%NP_GQ))
    ALLOCATE(FaceJ%XM_GQ(3,VFace%NP_GQ))
    EPS=IGREEN%EPS_ZMIN
    CALL VFace_to_FACE(VFace,FaceJ,J)    !Extract a face J from the VFace array
    !initialization
     SP=CZERO
     SM=CZERO
     VSP=CZERO
     VSM=CZERO

       DO IGQ=1, FaceJ%NP_GQ
           XI(:) = X0I(:)
           IF (I>0) XI(3) = MIN(X0I(3), -EPS*Mesh%xy_diameter) !for I on body panel, I=0 on free surface
           XJ(:) = FaceJ%XM_GQ(:,IGQ)
           XJ(3) = MIN(XJ(3), -EPS*Mesh%xy_diameter)
           CALL COMPUTE_S2(XI, XJ, INFINITE_DEPTH, wavenumber, IGreen, FS(1), VS(:, 1))

           IF (Mesh%Isym == NO_Y_SYMMETRY) THEN
             SP_IGQ   = FS(1)
             VSP_IGQ(1:3) = VS(1:3, 1)
             SM_IGQ       = CZERO
             VSM_IGQ      = CZERO

           ELSE IF (Mesh%Isym == Y_SYMMETRY) THEN
             ! Reflect the source point across the (xOz) plane and compute another coefficient
             XI(2) = -X0I(2)
             CALL COMPUTE_S2(XI, XJ, INFINITE_DEPTH, wavenumber, IGreen, FS(2), VS(:, 2))
             VS(2, 2) = -VS(2, 2) ! Reflection of the output vector

             ! Assemble the two results
             SP_IGQ       = FS(1)      + FS(2)
             VSP_IGQ(1:3) = VS(1:3, 1) + VS(1:3, 2)
             SM_IGQ       = FS(1)      - FS(2)
             VSM_IGQ(1:3) = VS(1:3, 1) - VS(1:3, 2)
           END IF

           ADPI2  = wavenumber*Mesh%A(J)/DPI2   *FaceJ%dXdXG_WGQ_per_A(IGQ)
           ADPI   = wavenumber*Mesh%A(J)/DPI    *FaceJ%dXdXG_WGQ_per_A(IGQ)
           AKDPI2 = wavenumber**2*Mesh%A(J)/DPI2*FaceJ%dXdXG_WGQ_per_A(IGQ)
           AKDPI  = wavenumber**2*Mesh%A(J)/DPI *FaceJ%dXdXG_WGQ_per_A(IGQ)

           SP  = SP+CMPLX(REAL(SP_IGQ)*ADPI2,   AIMAG(SP_IGQ)*ADPI)
           VSP = VSP+CMPLX(REAL(VSP_IGQ)*AKDPI2, AIMAG(VSP_IGQ)*AKDPI)

           IF (Mesh%ISym == Y_SYMMETRY) THEN
             SM  = SM+CMPLX(REAL(SM_IGQ)*ADPI2,   AIMAG(SM_IGQ)*ADPI)
             VSM = VSM+CMPLX(REAL(VSM_IGQ)*AKDPI2, AIMAG(VSM_IGQ)*AKDPI)
           END IF
       ENDDO

    DEALLOCATE(FaceJ%dXdXG_WGQ_per_A)
    DEALLOCATE(FaceJ%XM_GQ)

    RETURN
  END SUBROUTINE VNSINFD

  !------------------------------------------------

  SUBROUTINE VNSFD(I, wavenumber, X0I, J, VFace, Mesh,IGreen, depth, SP, SM, VSP, VSM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL,                  INTENT(IN)  :: wavenumber, depth
    REAL, DIMENSION(3),    INTENT(IN)  :: X0I   ! Coordinates of the computed flow field point
    INTEGER,               INTENT(IN)  :: I     ! Index of the flow integration panel
    INTEGER,               INTENT(IN)  :: J     ! Index of the source integration panel
    TYPE(TMesh),           INTENT(IN)  :: Mesh
    TYPE(TVFace),          INTENT(IN)  :: VFace
    TYPE(TGreen),          INTENT(IN)  :: IGreen ! Initial green variable

    ! Outputs
    COMPLEX,               INTENT(OUT) :: SP, SM   ! Integral of the Green function over the panel.
    COMPLEX, DIMENSION(3), INTENT(OUT) :: VSP, VSM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    INTEGER                                :: KE,Itemp
    TYPE(TFace)                            :: FaceJ
    REAL, DIMENSION(3)                     :: X0J   ! Coordinates of the source point
    REAL                                   :: AMH, AKH, A, COF1, COF2, COF3, COF4
    REAL                                   :: AQT, RRR
    REAL, DIMENSION(3)                     :: XI, XJ
    REAL, DIMENSION(4, 2**Mesh%Isym)       :: FTS, PSR
    REAL, DIMENSION(3, 4, 2**Mesh%Isym)    :: VTS
    COMPLEX, DIMENSION(4, 2**Mesh%Isym)    :: FS
    COMPLEX, DIMENSION(3, 4, 2**Mesh%Isym) :: VS
    INTEGER                                :: IGQ       !Index Gauss Quadrature Integration point
    COMPLEX               :: SP_IGQ, SM_IGQ  ! Integral of the Green function over the panel for Gauss point IGQ.
    COMPLEX, DIMENSION(3) :: VSP_IGQ, VSM_IGQ! Gradient of the integral of the Green function
                                             !  with respect to X0I and Gauss point IGQ.

    INTEGER                 :: NEXP
    REAL, DIMENSION(31)     :: AMBDA, AR
    REAL                    :: EPS
    ALLOCATE(FaceJ%dXdXG_WGQ_per_A(VFace%NP_GQ))
    ALLOCATE(FaceJ%XM_GQ(3,VFace%NP_GQ))

    !passing values
    NEXP =IGreen%NEXP
    AMBDA=IGreen%AMBDA(:)
    AR   =IGreen%AR(:)
    EPS  =IGREEN%EPS_ZMIN
    CALL VFace_to_FACE(VFace,FaceJ,J)    !Extract a face J from the VFace array
    !initialization
     SP=CZERO
     SM=CZERO
     VSP=CZERO
     VSM=CZERO
    !========================================
    ! Part 1: Solve 4 infinite depth problems
    !========================================

       DO IGQ=1, FaceJ%NP_GQ

           XI(:) = X0I(:)
           IF (I>0) XI(3) = MIN(X0I(3), -EPS*Mesh%xy_diameter) !for I on body panel, I=0 on free surface

           X0J(:)= FaceJ%XM_GQ(:,IGQ)
           XJ(:) = X0J
           XJ(3) = MIN(X0J(3), -EPS*Mesh%xy_diameter)

           ! Distance in xOy plane
           RRR = NORM2(XI(1:2) - XJ(1:2))

           ! 1.a First infinite depth problem
           CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(1, 1), VS(:, 1, 1))

           PSR(1, 1) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

           ! 1.b Shift and reflect XI and compute another value of the Green function
           XI(3) = -X0I(3) - 2*depth
           XJ(3) = MIN(X0J(3), -EPS*Mesh%xy_diameter)
           CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(2, 1), VS(:, 2, 1))
           VS(3, 2, 1) = -VS(3, 2, 1) ! Reflection of the output vector

           PSR(2, 1) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

           ! 1.c Shift and reflect XJ and compute another value of the Green function
           XI(3) = X0I(3)
           IF (I>0) XI(3) = MIN(X0I(3), -EPS*Mesh%xy_diameter) !for I on body panel, I=0 on free surface
           XJ(3) = -X0J(3) - 2*depth
           CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(3, 1), VS(:, 3, 1))

           PSR(3, 1) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

           ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
           XI(3) = -X0I(3)        - 2*depth
           XJ(3) = -X0J(3) - 2*depth
           CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(4, 1), VS(:, 4, 1))
           VS(3, 4, 1) = -VS(3, 4, 1) ! Reflection of the output vector

           PSR(4, 1) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

           IF (Mesh%ISym == NO_Y_SYMMETRY) THEN
             ! Add up the results of the four problems
             SP_IGQ       = -SUM(FS(1:4, 1)) - SUM(PSR(1:4, 1))
             VSP_IGQ(1:3) = -SUM(VS(1:3, 1:4, 1), 2)
             SM_IGQ       = CZERO
             VSM_IGQ(1:3) = CZERO

           ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
             ! If the y-symmetry is used, the four symmetric problems have to be solved
             XI(:) = X0I(:)
             XI(2) = -XI(2)
             IF (I>0) XI(3) = MIN(X0I(3), -EPS*Mesh%xy_diameter) !for I on body panel, I=0 on free surface
             XJ(:) = X0J
             XJ(3) = MIN(X0J(3), -EPS*Mesh%xy_diameter)

             RRR = NORM2(XI(1:2) - XJ(1:2))

             ! 1.a' First infinite depth problem
             CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(1, 2), VS(:, 1, 2))

             PSR(1, 2) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

             ! 1.b' Shift and reflect XI and compute another value of the Green function
             XI(3) = -X0I(3)        - 2*depth
             XJ(3) = MIN(X0J(3), -EPS*Mesh%xy_diameter)
             CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(2, 2), VS(:, 2, 2))
             VS(3, 2, 2) = -VS(3, 2, 2)

             PSR(2, 2) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

             ! 1.c' Shift and reflect XJ and compute another value of the Green function
             XI(3) = X0I(3)
             IF (I>0) XI(3) = MIN(X0I(3), -EPS*Mesh%xy_diameter) !for I on body panel, I=0 on free surface
             XJ(3) = -X0J(3)- 2*depth
             CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(3, 2), VS(:, 3, 2))

             PSR(3, 2) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

             ! 1.d' Shift and reflect both XI and XJ and compute another value of the Green function
             XI(3) = -X0I(3)        - 2*depth
             XJ(3) = -X0J(3)- 2*depth
             CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, IGreen, FS(4, 2), VS(:, 4, 2))
             VS(3, 4, 2) = -VS(3, 4, 2)

             PSR(4, 2) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

             ! Reflection of the four output vectors around xOz plane
             VS(2, 1:4, 2) = -VS(2, 1:4, 2)

             ! Add up the results of the 2×4 problems
             SP_IGQ       = -SUM(FS(1:4, 1))-SUM(PSR(1:4, 1))   - SUM(FS(1:4, 2))-SUM(PSR(1:4, 2))
             VSP_IGQ(1:3) = -SUM(VS(1:3, 1:4, 1), 2)            - SUM(VS(1:3, 1:4, 2), 2)
             SM_IGQ       = -SUM(FS(1:4, 1)) - SUM(PSR(1:4, 1)) + SUM(FS(1:4, 2)) + SUM(PSR(1:4, 2))
             VSM_IGQ(1:3) = -SUM(VS(1:3, 1:4, 1), 2)            + SUM(VS(1:3, 1:4, 2), 2)
           END IF
           ! Multiply by some coefficients
           AMH  = wavenumber*depth
           AKH  = AMH*TANH(AMH)
           A    = (AMH+AKH)**2/(depth*(AMH**2-AKH**2+AKH))
           COF1 = -A/(8*PI**2)*Mesh%A(J)*FaceJ%dXdXG_WGQ_per_A(IGQ)
           COF2 = -A/(8*PI)   *Mesh%A(J)*FaceJ%dXdXG_WGQ_per_A(IGQ)
           COF3 = wavenumber*COF1
           COF4 = wavenumber*COF2

           SP     = SP+CMPLX(REAL(SP_IGQ)*COF1,  AIMAG(SP_IGQ)*COF2)
           VSP(:) = VSP(:)+CMPLX(REAL(VSP_IGQ)*COF3, AIMAG(VSP_IGQ)*COF4)

           IF (Mesh%ISym == Y_SYMMETRY) THEN
             SM     = SM+CMPLX(REAL(SM_IGQ)*COF1,  AIMAG(SM_IGQ)*COF2)
             VSM(:) = VSM(:)+CMPLX(REAL(VSM_IGQ)*COF3, AIMAG(VSM_IGQ)*COF4)
           END IF

       END DO

       !=====================================================
       ! Part 2: Integrate (NEXP+1)×4 terms of the form 1/MM'
       !=====================================================
       DO IGQ=1, FaceJ%NP_GQ
           AMBDA(NEXP+1) = 0
           AR(NEXP+1)    = 2

          XJ(:) = FaceJ%XM_GQ(:,IGQ)

           DO KE = 1, NEXP+1
             XI(:) = X0I(:)

             ! 2.a Shift observation point and compute integral
             XI(3) =  X0I(3) + depth*AMBDA(KE) - 2*depth
             CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(1, 1), VTS(:, 1, 1))

             ! 2.b Shift and reflect observation point and compute integral
             XI(3) = -X0I(3) - depth*AMBDA(KE)
             CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(2, 1), VTS(:, 2, 1))
             VTS(3, 2, 1) = -VTS(3, 2, 1) ! Reflection of the output vector

             ! 2.c Shift and reflect observation point and compute integral
             XI(3) = -X0I(3) + depth*AMBDA(KE) - 4*depth
             CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(3, 1), VTS(:, 3, 1))
             VTS(3, 3, 1) = -VTS(3, 3, 1) ! Reflection of the output vector

             ! 2.d Shift observation point and compute integral
             XI(3) =  X0I(3) - depth*AMBDA(KE) + 2*depth
             CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(4, 1), VTS(:, 4, 1))

             AQT = -AR(KE)/(8*PI)

             IF (Mesh%ISym == NO_Y_SYMMETRY) THEN
               ! Add all the contributions
               SP     = SP     + AQT*SUM(FTS(1:4, 1))*FaceJ%dXdXG_WGQ_per_A(IGQ)
               VSP(:) = VSP(:) + AQT*SUM(VTS(1:3, 1:4, 1), 2)*FaceJ%dXdXG_WGQ_per_A(IGQ)
              !print*,REAL(VSP(3)),AIMAG(VSP(3))

             ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
               ! If the y-symmetry is used, the four symmetric problems have to be solved
               XI = X0I(:)
               XI(2) = -X0I(2)
               ! 2.a' Shift observation point and compute integral
               XI(3) =  X0I(3) + depth*AMBDA(KE) - 2*depth
               CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(1, 2), VTS(:, 1, 2))

               ! 2.b' Shift and reflect observation point and compute integral
               XI(3) = -X0I(3) - depth*AMBDA(KE)
               CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(2, 2), VTS(:, 2, 2))
               VTS(3, 2, 2) = -VTS(3, 2, 2) ! Reflection of the output vector

               ! 2.c' Shift and reflect observation point and compute integral
               XI(3) = -X0I(3) + depth*AMBDA(KE) - 4*depth
               CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(3, 2), VTS(:, 3, 2))
               VTS(3, 3, 2) = -VTS(3, 3, 2) ! Reflection of the output vector

               ! 2.d' Shift observation point and compute integral
               XI(3) =  X0I(3) - depth*AMBDA(KE) + 2*depth
               CALL COMPUTE_ASYMPTOTIC_S0(XI(:), XJ, FaceJ%A, EPS, FTS(4, 2), VTS(:, 4, 2))

               ! Reflection of the output vector around the xOz plane
               VTS(2, 1:4, 2) = -VTS(2, 1:4, 2)

               ! Add all the contributions
               SP     = SP     + AQT*(SUM(FTS(1:4, 1))         + SUM(FTS(1:4, 2))        )*FaceJ%dXdXG_WGQ_per_A(IGQ)
               VSP(:) = VSP(:) + AQT*(SUM(VTS(1:3, 1:4, 1), 2) + SUM(VTS(1:3, 1:4, 2), 2))*FaceJ%dXdXG_WGQ_per_A(IGQ)
               SM     = SM     + AQT*(SUM(FTS(1:4, 1))         - SUM(FTS(1:4, 2))        )*FaceJ%dXdXG_WGQ_per_A(IGQ)
               VSM(:) = VSM(:) + AQT*(SUM(VTS(1:3, 1:4, 1), 2) - SUM(VTS(1:3, 1:4, 2), 2))*FaceJ%dXdXG_WGQ_per_A(IGQ)
             END IF
           END DO
       END DO
    DEALLOCATE(FaceJ%dXdXG_WGQ_per_A)
    DEALLOCATE(FaceJ%XM_GQ)

    RETURN
  END SUBROUTINE

!------------------------------------------------------------------------

  SUBROUTINE COMPUTE_S2(XI, XJ, depth, wavenumber, IGreen, FS, VS)

    ! Inputs
    REAL, DIMENSION(3),    INTENT(IN)  :: XI, XJ
    REAL,                  INTENT(IN)  :: depth, wavenumber
   TYPE(TGREEN),           INTENT(IN)  :: IGreen

    ! Outputs
    COMPLEX,               INTENT(OUT) :: FS
    COMPLEX, DIMENSION(3), INTENT(OUT) :: VS

    ! Local variables
    INTEGER                            :: KI, KJ
    REAL                               :: RRR, AKR, ZZZ, AKZ, DD, PSURR
    REAL                               :: SIK, CSK, SQ, EPZ
    REAL                               :: PD1X, PD2X, PD1Z, PD2Z
    REAL, ALLOCATABLE,DIMENSION(:)     :: XL, ZL
    REAL                               :: EPS

    EPS=IGREEN%EPS_ZMIN

    IF (FLAG_IGREEN==1) THEN
      ALLOCATE(XL(3),ZL(3))
     ELSE
      ALLOCATE(XL(5),ZL(5))
    ENDIF

    RRR = NORM2(XI(1:2) - XJ(1:2))
    AKR = wavenumber*RRR

    ZZZ = XI(3) + XJ(3)
    AKZ = wavenumber*ZZZ

    DD  = SQRT(RRR**2 + ZZZ**2)

    IF (DD > EPS) THEN
      PSURR = PI/(wavenumber*DD)**3
    ELSE
      PSURR = 0.0
    ENDIF

    ! IF (AKZ > -1.5e-6) THEN
    !   WRITE(*,*)'AKZ < -1.5 E-6' ! Not a very explicit warning...
    ! END IF

    IF (AKZ > IGREEN%MIN_KZ) THEN             !   -16 < AKZ < -1.5e-6

      !================================================
      ! Evaluate PDnX and PDnZ depending on AKZ and AKR
      !================================================

      IF (AKR < IGREEN%MAX_KR) THEN          !     0 < AKR < 99.7
        IF (FLAG_IGREEN==1) THEN
                IF (AKZ < -1e-2) THEN       !   -16 < AKZ < -1e-2
                  KJ = INT(8*(ALOG10(-AKZ)+4.5))
                ELSE                        ! -1e-2 < AKZ < -1.5e-6
                  KJ = INT(5*(ALOG10(-AKZ)+6))
                ENDIF
                KJ = MAX(MIN(KJ, 45), 2)

                IF (AKR < 1) THEN           !     0 < AKR < 1
                  KI = INT(5*(ALOG10(AKR+1e-20)+6)+1)
                ELSE                        !     1 < AKR < 99.7
                  KI = INT(3*AKR+28)
                ENDIF
                KI = MAX(MIN(KI, 327), 2)

                XL(1) = PL2(IGreen%XR(KI),   IGreen%XR(KI+1), IGreen%XR(KI-1), AKR)
                XL(2) = PL2(IGreen%XR(KI+1), IGreen%XR(KI-1), IGreen%XR(KI),   AKR)
                XL(3) = PL2(IGreen%XR(KI-1), IGreen%XR(KI),   IGreen%XR(KI+1), AKR)
                ZL(1) = PL2(IGreen%XZ(KJ),   IGreen%XZ(KJ+1), IGreen%XZ(KJ-1), AKZ)
                ZL(2) = PL2(IGreen%XZ(KJ+1), IGreen%XZ(KJ-1), IGreen%XZ(KJ),   AKZ)
                ZL(3) = PL2(IGreen%XZ(KJ-1), IGreen%XZ(KJ),   IGreen%XZ(KJ+1), AKZ)

                PD1Z = DOT_PRODUCT(XL, MATMUL(IGreen%APD1Z(KI-1:KI+1, KJ-1:KJ+1), ZL))
                PD2Z = DOT_PRODUCT(XL, MATMUL(IGreen%APD2Z(KI-1:KI+1, KJ-1:KJ+1), ZL))

                IF (RRR > EPS) THEN
                  PD1X = DOT_PRODUCT(XL, MATMUL(IGreen%APD1X(KI-1:KI+1, KJ-1:KJ+1), ZL))
                  PD2X = DOT_PRODUCT(XL, MATMUL(IGreen%APD2X(KI-1:KI+1, KJ-1:KJ+1), ZL))
                END IF
        ELSE
                KJ=10*(ALOG10(-AKZ)+10.)
                KJ=MAX(KJ,3)
                KJ=MIN(KJ,IGREEN%JZ-2)

                IF(AKR.LT.1.)THEN
                  KI=10*(ALOG10(AKR+1.E-10)+8)+1
                ELSE
                  KI=6*AKR+75
                ENDIF
                KI=MAX(KI,3)
                KI=MIN(KI,674)

                XL(1)=PL5(IGreen%XR(KI+2),IGreen%XR(KI-1),IGreen%XR(KI  ),IGreen%XR(KI+1),IGreen%XR(KI-2),AKR)
                XL(2)=PL5(IGreen%XR(KI-2),IGreen%XR(KI  ),IGreen%XR(KI+1),IGreen%XR(KI+2),IGreen%XR(KI-1),AKR)
                XL(3)=PL5(IGreen%XR(KI-1),IGreen%XR(KI+1),IGreen%XR(KI+2),IGreen%XR(KI-2),IGreen%XR(KI  ),AKR)
                XL(4)=PL5(IGreen%XR(KI  ),IGreen%XR(KI+2),IGreen%XR(KI-2),IGreen%XR(KI-1),IGreen%XR(KI+1),AKR)
                XL(5)=PL5(IGreen%XR(KI+1),IGreen%XR(KI-2),IGreen%XR(KI-1),IGreen%XR(KI  ),IGreen%XR(KI+2),AKR)
                ZL(1)=PL5(IGreen%XZ(KJ+2),IGreen%XZ(KJ-1),IGreen%XZ(KJ  ),IGreen%XZ(KJ+1),IGreen%XZ(KJ-2),AKZ)
                ZL(2)=PL5(IGreen%XZ(KJ-2),IGreen%XZ(KJ  ),IGreen%XZ(KJ+1),IGreen%XZ(KJ+2),IGreen%XZ(KJ-1),AKZ)
                ZL(3)=PL5(IGreen%XZ(KJ-1),IGreen%XZ(KJ+1),IGreen%XZ(KJ+2),IGreen%XZ(KJ-2),IGreen%XZ(KJ  ),AKZ)
                ZL(4)=PL5(IGreen%XZ(KJ  ),IGreen%XZ(KJ+2),IGreen%XZ(KJ-2),IGreen%XZ(KJ-1),IGreen%XZ(KJ+1),AKZ)
                ZL(5)=PL5(IGreen%XZ(KJ+1),IGreen%XZ(KJ-2),IGreen%XZ(KJ-1),IGreen%XZ(KJ  ),IGreen%XZ(KJ+2),AKZ)

                PD1Z = DOT_PRODUCT(XL, MATMUL(IGreen%APD1Z(KI-2:KI+2, KJ-2:KJ+2), ZL))
                PD2Z = DOT_PRODUCT(XL, MATMUL(IGreen%APD2Z(KI-2:KI+2, KJ-2:KJ+2), ZL))

                IF (RRR > EPS) THEN
                  PD1X = DOT_PRODUCT(XL, MATMUL(IGreen%APD1X(KI-2:KI+2, KJ-2:KJ+2), ZL))
                  PD2X = DOT_PRODUCT(XL, MATMUL(IGreen%APD2X(KI-2:KI+2, KJ-2:KJ+2), ZL))
                END IF
        ENDIF
      ELSE  ! 99.7 < AKR

        EPZ  = EXP(AKZ)
        SQ   = SQRT(2*PI/AKR)
        CSK  = COS(AKR-PI/4)
        SIK  = SIN(AKR-PI/4)

        PD1Z = PSURR*AKZ - PI*EPZ*SQ*SIK
        PD2Z =                EPZ*SQ*CSK

        IF (RRR > EPS) THEN
          ! PD1X=-PSURR*AKR-PI*EPZ*SQ*(CSK-0.5/AKR*SIK) ! correction par GD le 17/09/2010
          PD1X = PI*EPZ*SQ*(CSK - 0.5*SIK/AKR) - PSURR*AKR
          PD2X =    EPZ*SQ*(SIK + 0.5*CSK/AKR)
        END IF

      ENDIF

      !====================================
      ! Deduce FS ans VS from PDnX and PDnZ
      !====================================

      FS    = -CMPLX(PD1Z, PD2Z)
      IF (depth == INFINITE_DEPTH) THEN
        VS(3) = -CMPLX(PD1Z-PSURR*AKZ, PD2Z)
      ELSE
        VS(3) = -CMPLX(PD1Z, PD2Z)
      END IF

      IF (RRR > EPS) THEN
        IF (depth == INFINITE_DEPTH) THEN
          VS(1) = (XI(1) - XJ(1))/RRR * CMPLX(PD1X+PSURR*AKR, PD2X)
          VS(2) = (XI(2) - XJ(2))/RRR * CMPLX(PD1X+PSURR*AKR, PD2X)
        ELSE
          VS(1) = (XI(1) - XJ(1))/RRR * CMPLX(PD1X, PD2X)
          VS(2) = (XI(2) - XJ(2))/RRR * CMPLX(PD1X, PD2X)
        END IF
      ELSE
        VS(1:2) = CZERO
      END IF

    ELSE ! AKZ < -16
      FS      = CMPLX(-PSURR*AKZ, 0.0)
      VS(1:3) = CZERO
    ENDIF

      DEALLOCATE(XL,ZL)
    RETURN
  END SUBROUTINE COMPUTE_S2

END MODULE GREEN_2
