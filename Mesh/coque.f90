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
!
SUBROUTINE coque(X,Y,Z,NP,facettes,NF,Deplacement,Icoque,Gcoque,CG,Nsym,rho)

  IMPLICIT NONE

!	Globale
	REAL,PARAMETER :: PI=3.141592653589
!	Description du maillage
	INTEGER Np,Nf
	REAL,DIMENSION(*) :: X,Y,Z
	INTEGER,DIMENSION(4,*) :: Facettes
	REAL Deplacement
	INTEGER Nsym
	REAL :: rho
!	Caracteristiques de la coque
	REAL,DIMENSION(3) :: Gcoque,CG
	REAL mcoque
	REAL,DIMENSION(3,3) :: Icoque
	REAL decoque
!	Locales
	REAL,DIMENSION(3) :: U,V,W
	REAL,DIMENSION(3,Nf) :: CdG
	REAL :: N1,N2,ZMN
	REAL,DIMENSION(Nf) :: Aire
	REAL,DIMENSION(4,3) :: P

!	Indices
	INTEGER i,j,IMN

	mcoque=0.
	Gcoque(1)=0.
	Gcoque(2)=0.
	Gcoque(3)=0.
	DO i=1,nF
                P(1,1)=X(FACETTES(1,i))  !added by RK for not calculating lid panels
                P(1,2)=Y(FACETTES(1,i))
                P(1,3)=Z(FACETTES(1,i))
                P(2,1)=X(FACETTES(2,i))
                P(2,2)=Y(FACETTES(2,i))
                P(2,3)=Z(FACETTES(2,i))
                P(3,1)=X(FACETTES(3,i))
                P(3,2)=Y(FACETTES(3,i))
                P(3,3)=Z(FACETTES(3,i))
                P(4,1)=X(FACETTES(4,i))
                P(4,2)=Y(FACETTES(4,i))
                P(4,3)=Z(FACETTES(4,i))
                IMN=1
                ZMN=P(IMN,3)
                DO J=IMN+1,4
                        IF (P(J,3).LT.ZMN) THEN
                                ZMN=P(J,3)
                                IMN=J
                        END IF
                END DO
                IF (ZMN.LT.0.0) THEN !added by RK for not calculating lid panels
                       CdG(1,i)=0.25*(X(Facettes(1,i))+X(Facettes(2,i))+X(Facettes(3,i))+X(Facettes(4,i)))
                       CdG(2,i)=0.25*(Y(Facettes(1,i))+Y(Facettes(2,i))+Y(Facettes(3,i))+Y(Facettes(4,i)))
                       CdG(3,i)=0.25*(Z(Facettes(1,i))+Z(Facettes(2,i))+Z(Facettes(3,i))+Z(Facettes(4,i)))
                       U(1)=X(Facettes(2,i))-X(Facettes(1,i))
                       U(2)=Y(Facettes(2,i))-Y(Facettes(1,i))
                       U(3)=Z(Facettes(2,i))-Z(Facettes(1,i))
                       V(1)=X(Facettes(3,i))-X(Facettes(1,i))
                       V(2)=Y(Facettes(3,i))-Y(Facettes(1,i))
                       V(3)=Z(Facettes(3,i))-Z(Facettes(1,i))
                       CALL prdvct(U,V,W)
                       N1=0.5*SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
                       U(1)=X(Facettes(3,i))-X(Facettes(1,i))
                       U(2)=Y(Facettes(3,i))-Y(Facettes(1,i))
                       U(3)=Z(Facettes(3,i))-Z(Facettes(1,i))
                       V(1)=X(Facettes(4,i))-X(Facettes(1,i))
                       V(2)=Y(Facettes(4,i))-Y(Facettes(1,i))
                       V(3)=Z(Facettes(4,i))-Z(Facettes(1,i))
                       CALL prdvct(U,V,W)
                       N2=0.5*SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
                       Aire(i)=N1+N2
                       mcoque=mcoque+Aire(i)
                       Gcoque(1)=Gcoque(1)+Aire(i)*CdG(1,i)
                       Gcoque(2)=Gcoque(2)+Aire(i)*CdG(2,i)
                       Gcoque(3)=Gcoque(3)+Aire(i)*CdG(3,i)
              END IF
        END DO
        decoque=Deplacement*rho/mcoque
        Gcoque(1)=Gcoque(1)/mcoque
        IF (Nsym.EQ.1) THEN
            Gcoque(2)=0.
        ELSE
            Gcoque(2)=Gcoque(2)/mcoque
        END IF
        Gcoque(3)=Gcoque(3)/mcoque
        DO i=1,3
        	DO j=1,3
        		Icoque(i,j)=0.
        	END DO
        END DO
        DO i=1,nF
                P(1,1)=X(FACETTES(1,I))
                P(1,2)=Y(FACETTES(1,I))
                P(1,3)=Z(FACETTES(1,I))
                P(2,1)=X(FACETTES(2,I))
                P(2,2)=Y(FACETTES(2,I))
                P(2,3)=Z(FACETTES(2,I))
                P(3,1)=X(FACETTES(3,I))
                P(3,2)=Y(FACETTES(3,I))
                P(3,3)=Z(FACETTES(3,I))
                P(4,1)=X(FACETTES(4,I))
                P(4,2)=Y(FACETTES(4,I))
                P(4,3)=Z(FACETTES(4,I))
                IMN=1
                ZMN=P(IMN,3)
                DO J=IMN+1,4
                        IF (P(J,3).LT.ZMN) THEN
                                ZMN=P(J,3)
                                IMN=J
                        END IF
                END DO
                IF (ZMN.LT.0.0) THEN  !added by RK for not calculating lid panels
!		Icoque(1,1)=Icoque(1,1)+Aire(i)*decoque*((CdG(2,i)-Gcoque(2))**2+(CdG(3,i)-Gcoque(3))**2)
!		Icoque(1,2)=Icoque(1,2)-Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(2,i)-Gcoque(2)))
!		Icoque(1,3)=Icoque(1,3)-Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(3,i)-Gcoque(3)))
!		Icoque(2,2)=Icoque(2,2)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(3,i)-Gcoque(3))**2)
!		Icoque(2,3)=Icoque(2,3)-Aire(i)*decoque*((CdG(2,i)-Gcoque(2))*(CdG(3,i)-Gcoque(3)))
!		Icoque(3,3)=Icoque(3,3)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(2,i)-Gcoque(2))**2)
                Icoque(1,1)=Icoque(1,1)+Aire(i)*decoque*((CdG(2,i)-CG(2))**2+(CdG(3,i)-CG(3))**2)
                Icoque(1,2)=Icoque(1,2)-Aire(i)*decoque*((CdG(1,i)-CG(1))*(CdG(2,i)-CG(2)))
                Icoque(1,3)=Icoque(1,3)-Aire(i)*decoque*((CdG(1,i)-CG(1))*(CdG(3,i)-CG(3)))
                Icoque(2,2)=Icoque(2,2)+Aire(i)*decoque*((CdG(1,i)-CG(1))**2+(CdG(3,i)-CG(3))**2)
                Icoque(2,3)=Icoque(2,3)-Aire(i)*decoque*((CdG(2,i)-CG(2))*(CdG(3,i)-CG(3)))
                Icoque(3,3)=Icoque(3,3)+Aire(i)*decoque*((CdG(1,i)-CG(1))**2+(CdG(2,i)-CG(2))**2)
                END IF
	END DO
	IF (Nsym.EQ.1) THEN
	    Icoque(1,2)=0.
	    Icoque(2,3)=0.
	END IF
	!
	Icoque(2,1)=Icoque(1,2)
	Icoque(3,1)=Icoque(1,3)
	Icoque(3,2)=Icoque(2,3)
  RETURN

END SUBROUTINE
