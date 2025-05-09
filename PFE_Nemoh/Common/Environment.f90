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
MODULE MEnvironment

  USE Elementary_functions, ONLY: CIH, SIH,X0,CIH_Vect,SIH_Vect

  IMPLICIT NONE

  ! Definition of TYPE Environment
  TYPE TEnvironment
    REAL :: RHO         ! Sea water density
    REAL :: G           ! Gravity constant
    REAL :: Depth       ! Water depth
    REAL :: Xeff, Yeff  ! Coordinates of point where the incident wave is measured
  END TYPE TEnvironment

CONTAINS

  SUBROUTINE ReadTEnvironment(Environment, file)
    ! Read Environment data from file

    CHARACTER(LEN=*),   INTENT(IN)  :: file
    TYPE(TEnvironment), INTENT(OUT) :: Environment

    INTEGER :: u

    OPEN(NEWUNIT=u, FILE=file, FORM='FORMATTED', STATUS='OLD')
    READ(u,*)
    READ(u,*) Environment%RHO
    READ(u,*) Environment%G
    READ(u,*) Environment%Depth
    READ(u,*) Environment%Xeff, Environment%Yeff
    CLOSE(u)

  END SUBROUTINE


  REAL FUNCTION Wavenumber(w, Environment)
    ! Calculate wave number for frequency w and given depth
    ! To be merge with X0 on Elementary_functions.

    REAL,               INTENT(IN) :: w
    TYPE(TEnvironment), INTENT(IN) :: Environment

    REAL :: k0,x0,xg,xd,xc
    INTEGER,PARAMETER :: Nitemx=10000
    INTEGER :: Nite

    Wavenumber=w*w/Environment%g
    x0=Wavenumber*Environment%Depth
    IF ((x0.LE.20.).AND.(x0.GT.0.)) THEN
      xg=0.
      xd=x0
      Nite=0
      DO WHILE ((Nite.LT.Nitemx).AND.((x0-xg*TANH(xg))*(x0-xd*TANH(xd)).GT.0.))
        xg=xd
        xd=2.*xd
        Nite=Nite+1
      END DO
      Nite=0
      IF (Nite.GE.Nitemx) THEN
        WRITE(*,*) 'Error: unable to find the wavenumber'
        STOP
      END IF
      xc=0.5*(xd+xg)
      DO WHILE ((Nite.LT.Nitemx).AND.(ABS(xd-xg)/ABS(xc).GE.1.0E-06))
        xc=0.5*(xd+xg)
        IF ((x0-xg*TANH(xg))*(x0-xc*TANH(xc)).GT.0.) THEN
          xg=xc
        ELSE
          xd=xc
        END IF
        Nite=Nite+1
      END DO
      IF (Nite.GE.Nitemx) THEN
        WRITE(*,*) 'Error: unable to find the wavenumber'
        STOP
      END IF
      Wavenumber=xc/Environment%Depth
    END IF
  END FUNCTION Wavenumber

  FUNCTION Fun_Dispersion(k,D,g)  result(ww)
        REAL, INTENT(IN) :: k,D,g
        REAL ww
        IF ((D == 0.) .OR. (k*D >= 20)) THEN
        ww=SQRT(g*k)
        ELSE
        ww=SQRT(g*k*tanh(k*D))
        ENDIF
  END FUNCTION

  FUNCTION Fun_inverseDispersion(w,D,g) result(wavenumber)
        !inverse of dispersion relation
        !should be same as the Wavenumber function above

        REAL, INTENT(IN) :: w,D,g
        REAL             :: wavenumber

        IF ((D == 0.) .OR. (w**2*D/g >= 20)) THEN
          wavenumber = w**2/g
        ELSE
          wavenumber = X0(w**2*D/g)/D
          ! X0(y) returns the solution of y = x * tanh(x)
        END IF
        RETURN
  END Function

  FUNCTION FunVect_inverseDispersion(Nw,w,D,g) result(wavenumber)
        !inverse of dispersion relation
        !should be same as the Wavenumber function above

        INTEGER,            INTENT(IN) :: Nw
        REAL,               INTENT(IN) :: D,g
        REAL,DIMENSION(Nw), INTENT(IN) :: w
        INTEGER                        :: Iw
        REAL,DIMENSION(Nw)             :: wavenumber
        DO Iw=1,Nw
             IF ((D == 0.) .OR. (w(Iw)**2*D/g >= 20)) THEN
               wavenumber(Iw) = w(Iw)**2/g
             ELSE
               wavenumber(Iw) = X0(w(Iw)**2*D/g)/D
               ! X0(y) returns the solution of y = x * tanh(x)
             END IF
        ENDDO
        RETURN
  END Function


  SUBROUTINE Compute_Wave(k,w,beta,x,y,z,Phi,p,Vx,Vy,Vz,Environment)
    ! Calculate the complex potential, pressure and fluid velocities for a regular wave eta=sin(k*wbar-wt)

    REAL :: k,w,beta,x,y,z
    COMPLEX :: Phi,p,Vx,Vy,Vz
    TYPE(TEnvironment) :: Environment
    REAL :: wbar
    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)

    wbar=(x-Environment%XEFF)*COS(Beta)+(y-Environment%YEFF)*SIN(Beta)
    Phi=-II*Environment%g/w*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    p=Environment%rho*Environment%g*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vx=Environment%g/w*k*COS(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vy=Environment%g/w*k*SIN(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vz=-II*Environment%g/w*k*SIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
  END SUBROUTINE Compute_wave

  SUBROUTINE COMPUTE_INC_POTENTIAL_VELOCITY(k,w,beta,                         &
                                  XM,Npanels,XM_ADD,NP_Add,                   &
                                  Environment,Isym,Potential,Velocity)
    ! Calculate the complex potential, and fluid velocities for a regular wave eta=sin(k*wbar-wt)
    ! INPUT/OUTPUT
    REAL,                                INTENT(IN) :: k,w,beta
    TYPE(TEnvironment),                  INTENT(IN) :: Environment
    INTEGER,                             INTENT(IN) :: Npanels,NP_Add
    INTEGER,                             INTENT(IN) :: Isym             !1=symmetry
    REAL,DIMENSION(Npanels,3),           INTENT(IN) :: XM
    REAL,DIMENSION(NP_Add,3),            INTENT(IN) :: XM_Add
        COMPLEX,DIMENSION((Npanels+NP_Add)*2**Isym),   INTENT(OUT):: Potential
    COMPLEX,DIMENSION((Npanels+NP_Add)*2**Isym,3), INTENT(OUT):: Velocity
    !Local var
    INTEGER                             :: NPTOT
    REAL,DIMENSION(Npanels+NP_Add)      :: wbar,QZ,dzQZ
    REAL,DIMENSION(Npanels+NP_Add,3)    :: XM_ALL
    COMPLEX,PARAMETER                   :: II=CMPLX(0.,1.)
    NPTOT=Npanels+NP_Add
    XM_ALL(1:Npanels,:)=XM
    IF (NP_Add>0)  THEN
    XM_ALL(Npanels+1:NPanels+NP_Add,:)=XM_Add
    ENDIF
    QZ=CIH_Vect(k,XM_ALL(:,3),Environment%Depth,Npanels+NP_Add)
    dzQZ=k*SIH_Vect(k,XM_ALL(:,3),Environment%Depth,Npanels+NP_Add)

    wbar=(XM_ALL(:,1)-Environment%XEFF)*COS(beta)+(XM_ALL(:,2)-Environment%YEFF)*SIN(beta)
    Potential(1:NPTOT)=-II*Environment%g/w*QZ*CEXP(II*k*wbar)

    Velocity(1:NPTOT,1)=Environment%g/w*k*COS(beta)*QZ*CEXP(II*k*wbar)
    Velocity(1:NPTOT,2)=Environment%g/w*k*SIN(beta)*QZ*CEXP(II*k*wbar)
    Velocity(1:NPTOT,3)=-II*Environment%g/w*dzQZ*CEXP(II*k*wbar)
    IF (ISYM==1)THEN
    wbar=(XM_ALL(:,1)-Environment%XEFF)*COS(beta)+(-XM_ALL(:,2)-Environment%YEFF)*SIN(beta)
    Potential(NPTOT+1:2*NPTOT)=-II*Environment%g/w*QZ*CEXP(II*k*wbar)
    Velocity(NPTOT+1:2*NPTOT,1)=Environment%g/w*k*COS(beta)*QZ*CEXP(II*k*wbar)
    Velocity(NPTOT+1:2*NPTOT,2)=Environment%g/w*k*SIN(beta)*QZ*CEXP(II*k*wbar)
    Velocity(NPTOT+1:2*NPTOT,3)=-II*Environment%g/w*dzQZ*CEXP(II*k*wbar)
    ENDIF

  END SUBROUTINE COMPUTE_INC_POTENTIAL_VELOCITY

END MODULE
