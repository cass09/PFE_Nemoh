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
!   Contributors list:
!    - Ruddy Kurnia (ECN)
!--------------------------------------------------------------------------------------
Module MInfluenceMatrix

USE CONSTANTS
USE MEnvironment,               ONLY: TEnvironment
USE M_INITIALIZE_GREEN,         ONLY: TGREEN, LISC
USE GREEN_1,                    ONLY: VAV
USE GREEN_2,                    ONLY: VNSINFD, VNSFD
USE MMesh,                      ONLY: TMesh
USE MFace,                      ONLY: TVFace
IMPLICIT NONE

CONTAINS
        SUBROUTINE CONSTRUCT_INFLUENCE_MATRIX_FS(omega,wavenumber,Env,IGreen,Mesh,VFace,XM_I,S,gradS)
        !======================================
        !Construction of the influence matrix
        !=====================================
        !INPUT/OUTPUT
        REAL,                           INTENT(IN)::omega,wavenumber
        TYPE(TEnvironment),             INTENT(IN)::Env
        TYPE(TMesh),                    INTENT(IN)::Mesh
        TYPE(TVFace),                   INTENT(IN)::VFace
        TYPE(TGREEN),                   INTENT(INOUT)::IGreen
        REAL                            ,INTENT(IN):: XM_I(3)
        COMPLEX,DIMENSION(Mesh%NPanels,2**Mesh%Isym),      &
                                        INTENT(OUT)::S
        COMPLEX,DIMENSION(Mesh%NPanels,3,2**Mesh%Isym),      &
                                        INTENT(OUT)::gradS
        !LOCAL
        INTEGER   :: J
        ! Return of GREEN_1 module
        REAL :: FSP, FSM
        REAL, DIMENSION(3) :: VSXP, VSXM
        ! Return of GREEN_2 module
        COMPLEX :: SP, SM
        COMPLEX, DIMENSION(3) :: VSP, VSM

           DO J = 1, Mesh%NPanels
             IF ((Env%depth == INFINITE_DEPTH).OR.(omega**2*Env%depth/Env%g.GE.20)) THEN
               ! First part of the Green function
               ! These output are independent of omega, but always computed here to avoid
               ! large memory used for saving the matrix, particularly for the free surface panels
               CALL VAV                                                     &
               (0, XM_I, J,VFace,Mesh, INFINITE_DEPTH,IGREEN%EPS_ZMIN,&
                 FSP, FSM, VSXP, VSXM &
                 )
               ! Second part of the Green function
               CALL VNSINFD                                            &
               (0, wavenumber,XM_I , J, VFace, Mesh,                     &
                 IGreen, SP, SM, VSP, VSM                              &
                 )
             ELSE
               CALL VAV                                                &
               (0, XM_I, J,VFace,Mesh, Env%depth,IGREEN%EPS_ZMIN,&
                 FSP, FSM, VSXP, VSXM                                  &
                 )
               CALL VNSFD                                              &
               (0, wavenumber, XM_I, J, VFace,  Mesh,                    &
                   IGreen,Env%depth, SP, SM, VSP, VSM                  &
                 )
             END IF

             ! Store into influence matrix 'only for row-Ipflow
             S(J, 1) = FSP + SP                    ! Green function
             gradS(J,:, 1) = VSXP + VSP            ! Gradient of the Green function
             !print*,S(J,1)
             IF (Mesh%ISym == Y_SYMMETRY) THEN
               S(J, 2)      = FSM + SM
               gradS(J,:,2) = VSXM + VSM
             ENDIF

           END DO
         END SUBROUTINE

        SUBROUTINE CONSTRUCT_INFLUENCE_MATRIX(omega,wavenumber,Env,IGreen,Mesh,VFace,XM_add,NP_add,S,gradS)
        !======================================
        !Construction of the influence matrix
        !=====================================
        !INPUT/OUTPUT
        REAL,                           INTENT(IN)::omega,wavenumber
        TYPE(TEnvironment),             INTENT(IN)::Env
        TYPE(TMesh),                    INTENT(IN)::Mesh
        TYPE(TVFace),                   INTENT(IN)::VFace
        INTEGER,                        INTENT(IN)::NP_add
        REAL,DIMENSION(NP_add,3),       INTENT(IN)::XM_add
        TYPE(TGREEN),                   INTENT(INOUT)::IGreen
        COMPLEX,DIMENSION(Mesh%NPanels+NP_add,Mesh%NPanels,2**Mesh%Isym),      &
                                        INTENT(OUT)::S
        COMPLEX,DIMENSION(Mesh%NPanels+NP_add,Mesh%NPanels,3,2**Mesh%Isym),      &
                                        INTENT(OUT)::gradS
        !LOCAL
        INTEGER   :: I,J
        REAL      :: XM_I(3)
        ! Return of GREEN_1 module
        REAL :: FSP, FSM
        REAL, DIMENSION(3) :: VSXP, VSXM
        ! Return of GREEN_2 module
        COMPLEX :: SP, SM
        COMPLEX, DIMENSION(3) :: VSP, VSM


         IF (.NOT. Env%depth == INFINITE_DEPTH) THEN
           CALL LISC(omega**2*Env%depth/Env%g, wavenumber*Env%depth,IGreen)
         END IF
         DO I = 1, Mesh%NPanels+NP_add
           IF (I<=Mesh%NPanels) THEN
             XM_I=VFace%XM(I,:)
           ELSE
             XM_I=XM_add(I-Mesh%NPanels,:)
           ENDIF

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
               (I, wavenumber,XM_I , J, VFace, Mesh,     &
                 IGreen, SP, SM, VSP, VSM                       &
                 )
             ELSE
               CALL VNSFD                                       &
               (I, wavenumber, XM_I, J, VFace,  Mesh,    &
                   IGreen,Env%depth, SP, SM, VSP, VSM           &
                 )
             END IF

             ! Store into influence matrix
             S(I, J, 1) = FSP + SP                    ! Green function
             gradS(I, J,:, 1) = VSXP + VSP            ! Gradient of the Green function

             IF (Mesh%ISym == Y_SYMMETRY) THEN
               S(I, J, 2)      = FSM + SM
               gradS(I, J,:,2) = VSXM + VSM
             ENDIF

           END DO
         END DO

         END SUBROUTINE
END Module
