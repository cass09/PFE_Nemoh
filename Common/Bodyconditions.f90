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
    MODULE MBodyConditions

    PUBLIC
!
    TYPE TBodyConditions
        INTEGER :: Nproblems
        INTEGER :: Npanels
        REAL,DIMENSION(:),ALLOCATABLE :: Omega
        COMPLEX,DIMENSION(:,:),ALLOCATABLE :: NormalVelocity
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_Potential
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_FreeSurface
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_Kochin
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_SourceDistr
        INTEGER,DIMENSION(:),ALLOCATABLE :: Switch_Type
    END TYPE TBodyConditions
!
    CONTAINS
!
!       Operators for creation, copy, initialisation and destruction
!
        SUBROUTINE CreateTBodyConditions(BodyConditions,Nproblems,Npanels)
        IMPLICIT NONE
        TYPE(TBodyConditions) :: BodyConditions
        INTEGER :: Nproblems,Npanels
        BodyConditions%Nproblems=Nproblems
        BodyConditions%Npanels=Npanels
        ALLOCATE(BodyConditions%NormalVelocity(Npanels,Nproblems))
        ALLOCATE(BodyConditions%Omega(Nproblems),BodyConditions%Switch_Potential(Nproblems),&
        BodyConditions%Switch_FreeSurface(Nproblems),BodyConditions%Switch_Kochin(Nproblems),&
        BodyConditions%Switch_SourceDistr(Nproblems),BodyConditions%Switch_Type(Nproblems))
        END SUBROUTINE CreateTBodyConditions
!       ---
        SUBROUTINE CopyTBodyConditions(BodyConditionsTarget,BodyConditionsSource)
        IMPLICIT NONE
        INTEGER :: i,k
        TYPE(TBodyConditions) :: BodyConditionsTarget,BodyConditionsSource
        CALL CreateTBodyConditions(BodyConditionsTarget,BodyConditionsSource%Nproblems,BodyConditionsSource%Npanels)
        DO i=1,BodyConditionsTarget%Nproblems
            BodyConditionsTarget%Omega(i)=BodyConditionsSource%Omega(i)
            BodyConditionsTarget%Switch_Potential(i)=BodyConditionsSource%Switch_Potential(i)
            BodyConditionsTarget%Switch_FreeSurface(i)=BodyConditionsSource%Switch_FreeSurface(i)
            BodyConditionsTarget%Switch_Kochin(i)=BodyConditionsSource%Switch_Kochin(i)
            BodyConditionsTarget%Switch_SourceDistr(i)=BodyConditionsSource%Switch_SourceDistr(i)
            BodyConditionsTarget%Switch_Type(i)=BodyConditionsSource%Switch_Type(i)
            DO k=1,BodyConditionsTarget%Npanels
                BodyConditionsTarget%NormalVelocity(k,i)=BodyConditionsSource%NormalVelocity(k,i)
            END DO
        END DO
        END SUBROUTINE CopyTBodyConditions
!       ---
        SUBROUTINE ReadTBodyConditions(BodyConditions,Npanels,namefile)
        IMPLICIT NONE
        TYPE(TBodyConditions) :: BodyConditions
        CHARACTER(LEN=*) :: namefile
        INTEGER :: Nproblems,Npanels
        INTEGER :: i,k
        REAL,DIMENSION(:),ALLOCATABLE :: RBC,IBC
        OPEN(10,FILE=namefile)
        READ(10,*) Nproblems
        CALL CreateTBodyConditions(BodyConditions,Nproblems,Npanels)
        ALLOCATE(RBC(Nproblems),IBC(Nproblems))
        READ(10,*) (BodyConditions%Omega(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_Type(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_Potential(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_FreeSurface(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_Kochin(i),i=1,Nproblems)
        READ(10,*) (BodyConditions%Switch_SourceDistr(i),i=1,Nproblems)
        DO k=1,Npanels
            READ(10,*) (RBC(i),IBC(i),i=1,Nproblems)
            DO i=1,Nproblems
                BodyConditions%NormalVelocity(k,i)=CMPLX(RBC(i),IBC(i))
            END DO
        END DO
        CLOSE(10)
        DEALLOCATE(RBC,IBC)
        END SUBROUTINE ReadTBodyConditions
!       ---
        SUBROUTINE DeleteTBodyConditions(BodyConditions)
        IMPLICIT NONE
        TYPE(TBodyConditions) :: BodyConditions
        DEALLOCATE(BodyConditions%Omega,BodyConditions%Switch_potential,BodyConditions%Switch_FreeSurface,&
                BodyConditions%Switch_Kochin,BodyConditions%Switch_SourceDistr,BodyConditions%Switch_Type)
        DEALLOCATE(BodyConditions%NormalVelocity)
        END SUBROUTINE DeleteTBodyConditions
!       ---
END MODULE
