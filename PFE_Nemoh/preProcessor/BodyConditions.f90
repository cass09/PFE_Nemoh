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
!   - A. Babarit
!
!--------------------------------------------------------------------------------------
MODULE BodyConditions

IMPLICIT NONE

CONTAINS

!-- SUBROUTINE ComputeRadiationCondition

  SUBROUTINE ComputeRadiationCondition(Mesh,c,iCase,Direction,Axis,NVEL)

    USE MMesh

    IMPLICIT NONE
    TYPE(TMesh) :: Mesh
    INTEGER :: c,iCase
    REAL,DIMENSION(3) :: Direction,Axis
    COMPLEX,DIMENSION(:) :: NVEL
    REAL,DIMENSION(3) :: VEL
    INTEGER :: i

    SELECT CASE (iCase)
    CASE (1)
!       Degree of freedom is a translation
        DO i=1,Mesh%Npanels
          IF (Mesh%XM(3,i).lt.0.) THEN !if ZMN<0, dont calculate on the lid meshes (for irregular freq) by RK
           IF (Mesh%cPanel(i).EQ.c) THEN
               VEL(1)=Direction(1)
               VEL(2)=Direction(2)
               VEL(3)=Direction(3)
               NVEL(i)=CMPLX(Mesh%N(1,i)*VEL(1)+Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
           ELSE
               NVEL(i)=CMPLX(0.,0.)
           END IF
           IF (Mesh%iSym.EQ.1) THEN
               IF (Mesh%cPanel(i).EQ.c) THEN
                   VEL(1)=Direction(1)
                   VEL(2)=Direction(2)
                   VEL(3)=Direction(3)
                   NVEL(i+Mesh%Npanels)=CMPLX(Mesh%N(1,i)*VEL(1)-Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
               ELSE
                   NVEL(i+Mesh%Npanels)=CMPLX(0.,0.)
               END IF
           END IF
          ELSE
                NVEL(i)=CMPLX(0.,0.)
                IF (Mesh%iSym.EQ.1) THEN
                    NVEL(i+Mesh%Npanels)=CMPLX(0.,0.)
                END IF
          ENDIF
        END DO
    CASE (2)
!       Degree of freedom is a rotation
        DO i=1,Mesh%Npanels
          IF (Mesh%XM(3,i).lt.0.) THEN !if ZMN<0, dont calculate on the lid meshes (for irregular freq) by RK
           IF (Mesh%cPanel(i).EQ.c) THEN
               VEL(1)=Direction(2)*(Mesh%XM(3,i)-Axis(3))-Direction(3)*(Mesh%XM(2,i)-Axis(2))
               VEL(2)=Direction(3)*(Mesh%XM(1,i)-Axis(1))-Direction(1)*(Mesh%XM(3,i)-Axis(3))
               VEL(3)=Direction(1)*(Mesh%XM(2,i)-Axis(2))-Direction(2)*(Mesh%XM(1,i)-Axis(1))
               NVEL(i)=CMPLX(Mesh%N(1,i)*VEL(1)+Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
           ELSE
               NVEL(i)=CMPLX(0.,0.)
           END IF
           IF (Mesh%iSym.EQ.1) THEN
               IF (Mesh%cPanel(i).EQ.c) THEN
                   VEL(1)=Direction(2)*(Mesh%XM(3,i)-Axis(3))-Direction(3)*(-Mesh%XM(2,i)-Axis(2))
                   VEL(2)=Direction(3)*(Mesh%XM(1,i)-Axis(1))-Direction(1)*(Mesh%XM(3,i)-Axis(3))
                   VEL(3)=Direction(1)*(-Mesh%XM(2,i)-Axis(2))-Direction(2)*(Mesh%XM(1,i)-Axis(1))
                   NVEL(i+Mesh%Npanels)=CMPLX(Mesh%N(1,i)*VEL(1)-Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
               ELSE
                   NVEL(i+Mesh%Npanels)=CMPLX(0.,0.)
               END IF
           END IF
          ELSE
                NVEL(i)=CMPLX(0.,0.)
                IF (Mesh%iSym.EQ.1) THEN
                    NVEL(i+Mesh%Npanels)=CMPLX(0.,0.)
                END IF
          ENDIF
        END DO
    CASE (3)
        WRITE(*,*) 'Error: radiation case 3 not implemented yet'
        STOP
    CASE DEFAULT
        WRITE(*,*) 'Error: unknown radiation case'
        STOP
    END SELECT
    END SUBROUTINE

!-- SUBROUTINE ComputeDiffractionCondition
  SUBROUTINE ComputeDiffractionCondition(Mesh,w,beta,Environment,PRESSURE,NVEL)

    USE Constants !, only: PI
    USE MEnvironment
    USE MMEsh

    IMPLICIT NONE
!   Inputs/outputs
    TYPE(TMesh)             :: Mesh
    REAL                    :: w,beta           ! Wave period, direction and wavenumber
    TYPE(TEnvironment)      :: Environment      ! Environment
    COMPLEX,DIMENSION(*)    :: PRESSURE,NVEL    ! Pressure and normal velocities on panels
!   Locals
    REAL :: kwaveh,kwave
!    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
    REAL :: wbar
    INTEGER :: i,j
    COMPLEX :: tmp
    COMPLEX :: Phi,p,Vx,Vy,Vz

    kwave=Wavenumber(w,Environment)
!   Compute potential and normal velocities
    DO i=1,2**Mesh%Isym*Mesh%Npanels
        IF (i.LE.Mesh%Npanels) THEN
           CALL Compute_Wave(kwave,w,beta,Mesh%XM(1,i),Mesh%XM(2,i),Mesh%XM(3,i),Phi,p,Vx,Vy,Vz,Environment)
           IF (Mesh%XM(3,i).lt.0.) THEN !if ZMN<0, dont calculate on the lid meshes (for irregular freq) by RK
            PRESSURE(i)=p
            NVEL(i)=-(Vx*Mesh%N(1,i)+Vy*Mesh%N(2,i)+Vz*Mesh%N(3,i))
           ELSE
           PRESSURE(i)=0
           NVEL(i)=0                   !forcing to be zero at lid panels for the extended BIE irreg freq. removal
           END IF
        ELSE
           CALL Compute_Wave(kwave,w,beta,Mesh%XM(1,i-Mesh%Npanels),-Mesh%XM(2,i-Mesh%Npanels),Mesh%XM(3,i-Mesh%Npanels),Phi,p,Vx,Vy,Vz,Environment)
           IF (Mesh%XM(3,i-Mesh%Npanels).lt.0.) THEN !if ZMN<0, dont calculate on the lid meshes (for irregular freq) by RK
           PRESSURE(i)=p
           NVEL(i)=-(Vx*Mesh%N(1,i-Mesh%Npanels)-Vy*Mesh%N(2,i-Mesh%Npanels)+Vz*Mesh%N(3,i-Mesh%Npanels))
            ELSE
            PRESSURE(i)=0
            NVEL(i)=0
            END IF
        END IF
   END DO
   END SUBROUTINE ComputeDiffractionCondition
!--
END MODULE
