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
SUBROUTINE calCol(NFm,Xm,Ym,Zm,Facettem,Tcol,nFacemx)
  !
  IMPLICIT NONE
  ! Maillage de la partie mouille
  INTEGER NFm,nFacemx
  REAL,DIMENSION(*) :: Xm,Ym,Zm
  INTEGER,DIMENSION(4,*) :: Facettem
  ! Hauteur de col
  REAL Tcol
  ! Locales
  INTEGER :: nFacem2,Np
  REAL,DIMENSION(4,3,nFacemx) :: Coinm2
  INTEGER :: i,j,k


  nFacem2=NFm
  DO i=1,NFm
    DO j=1,4
      Coinm2(j,1,i)=Xm(Facettem(j,i))
      Coinm2(j,2,i)=Ym(Facettem(j,i))
      Coinm2(j,3,i)=Zm(Facettem(j,i))-Tcol
    END DO
  END DO
  DO i=1,Nfm
    IF (((Coinm2(1,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(2,3,i)+Tcol).GT.-1.0E-03)) THEN
      Nfacem2=Nfacem2+1
      Coinm2(1,1,Nfacem2)=Coinm2(2,1,i)
      Coinm2(1,2,Nfacem2)=Coinm2(2,2,i)
      Coinm2(1,3,Nfacem2)=Coinm2(2,3,i)
      Coinm2(2,1,Nfacem2)=Coinm2(1,1,i)
      Coinm2(2,2,Nfacem2)=Coinm2(1,2,i)
      Coinm2(2,3,Nfacem2)=Coinm2(1,3,i)
      Coinm2(3,1,Nfacem2)=Coinm2(1,1,i)
      Coinm2(3,2,Nfacem2)=Coinm2(1,2,i)
      Coinm2(3,3,Nfacem2)=0.
      Coinm2(4,1,Nfacem2)=Coinm2(2,1,i)
      Coinm2(4,2,Nfacem2)=Coinm2(2,2,i)
      Coinm2(4,3,Nfacem2)=0.
    ELSE
      IF (((Coinm2(2,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(3,3,i)+Tcol).GT.-1.0E-03)) THEN
        Nfacem2=Nfacem2+1
        Coinm2(1,1,Nfacem2)=Coinm2(3,1,i)
        Coinm2(1,2,Nfacem2)=Coinm2(3,2,i)
        Coinm2(1,3,Nfacem2)=Coinm2(3,3,i)
        Coinm2(2,1,Nfacem2)=Coinm2(2,1,i)
        Coinm2(2,2,Nfacem2)=Coinm2(2,2,i)
        Coinm2(2,3,Nfacem2)=Coinm2(2,3,i)
        Coinm2(3,1,Nfacem2)=Coinm2(2,1,i)
        Coinm2(3,2,Nfacem2)=Coinm2(2,2,i)
        Coinm2(3,3,Nfacem2)=0.
        Coinm2(4,1,Nfacem2)=Coinm2(3,1,i)
        Coinm2(4,2,Nfacem2)=Coinm2(3,2,i)
        Coinm2(4,3,Nfacem2)=0.
      ELSE
        IF (((Coinm2(3,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(4,3,i)+Tcol).GT.-1.0E-03)) THEN
          Nfacem2=Nfacem2+1
          Coinm2(1,1,Nfacem2)=Coinm2(4,1,i)
          Coinm2(1,2,Nfacem2)=Coinm2(4,2,i)
          Coinm2(1,3,Nfacem2)=Coinm2(4,3,i)
          Coinm2(2,1,Nfacem2)=Coinm2(3,1,i)
          Coinm2(2,2,Nfacem2)=Coinm2(3,2,i)
          Coinm2(2,3,Nfacem2)=Coinm2(3,3,i)
          Coinm2(3,1,Nfacem2)=Coinm2(3,1,i)
          Coinm2(3,2,Nfacem2)=Coinm2(3,2,i)
          Coinm2(3,3,Nfacem2)=0.
          Coinm2(4,1,Nfacem2)=Coinm2(4,1,i)
          Coinm2(4,2,Nfacem2)=Coinm2(4,2,i)
          Coinm2(4,3,Nfacem2)=0.
        ELSE
          IF (((Coinm2(4,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(1,3,i)+Tcol).GT.-1.0E-03)) THEN
            Nfacem2=Nfacem2+1
            Coinm2(1,1,Nfacem2)=Coinm2(1,1,i)
            Coinm2(1,2,Nfacem2)=Coinm2(1,2,i)
            Coinm2(1,3,Nfacem2)=Coinm2(1,3,i)
            Coinm2(2,1,Nfacem2)=Coinm2(4,1,i)
            Coinm2(2,2,Nfacem2)=Coinm2(4,2,i)
            Coinm2(2,3,Nfacem2)=Coinm2(4,3,i)
            Coinm2(3,1,Nfacem2)=Coinm2(4,1,i)
            Coinm2(3,2,Nfacem2)=Coinm2(4,2,i)
            Coinm2(3,3,Nfacem2)=0.
            Coinm2(4,1,Nfacem2)=Coinm2(1,1,i)
            Coinm2(4,2,Nfacem2)=Coinm2(1,2,i)
            Coinm2(4,3,Nfacem2)=0.
          END IF
        END IF
      END IF
    END IF
  END DO
  Np=0
  NFm=Nfacem2
  DO i=1,nFacem2
    DO j=1,4
      Xm(Np+j)=Coinm2(j,1,i)
      Ym(Np+j)=Coinm2(j,2,i)
      Zm(Np+j)=Coinm2(j,3,i)
      Facettem(j,i)=Np+j
    END DO
    Np=Np+4
  END DO

END SUBROUTINE
