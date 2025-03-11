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
!	  - R. Kurnia
!--------------------------------------------------------------------------------------
    MODULE MIRF
!
    TYPE TIRF
        INTEGER :: Switch,Ntime,Nradiation,Nintegration,Nbeta,NtimeS
        REAL,DIMENSION(:),ALLOCATABLE :: Time,beta
        REAL,DIMENSION(:),ALLOCATABLE :: TimeS
        REAL,DIMENSION(:,:,:),ALLOCATABLE :: K
        COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: KexcForce
        REAL,DIMENSION(:,:),ALLOCATABLE :: AddedMass
    END TYPE TIRF
!
    CONTAINS
!
!       Operators for creation, copy, initialisation and destruction
!
        SUBROUTINE CreateTIRF(IRF,Ntime,Nradiation,Nintegration,Nbeta,NtimeS)
        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        INTEGER :: Ntime,Nradiation,Nintegration,Nbeta,NtimeS
        IRF%Ntime=Ntime
        IRF%NtimeS=NtimeS
        IRF%Nradiation=Nradiation
        ALLOCATE(IRF%Time(1:Ntime),IRF%AddedMass(Nradiation,Nintegration),IRF%K(1:Ntime,Nradiation,Nintegration))
        ALLOCATE(IRF%TimeS(NtimeS),IRF%KexcForce(NtimeS,Nbeta,Nintegration))
        END SUBROUTINE CreateTIRF
!       ---
        SUBROUTINE CopyTIRF(IRFTarget,IRFSource)
        IMPLICIT NONE
        INTEGER :: i,j,k
        TYPE(TIRF) :: IRFTarget,IRFSource
        CALL CreateTIRF(IRFTarget,IRFSource%Ntime,IRFSource%Nradiation,IRFSource%Nintegration,IRFSource%Nbeta,IRFSource%NtimeS)
        DO i=1,IRFTarget%Ntime
            IRFTarget%Time(i)=IRFSource%Time(i)
        END DO
        DO i=1,IRFTarget%NtimeS
            IRFTarget%TimeS(i)=IRFSource%TimeS(i)
        END DO

        DO j=1,IRFTarget%Nradiation
            DO k=1,IRFTarget%Nintegration
                IRFTarget%Addedmass(j,k)=IRFSource%Addedmass(j,k)
            END DO
        END DO
        DO i=1,IRFTarget%Ntime
            DO j=1,IRFTarget%Nradiation
                DO k=1,IRFTarget%Nintegration
                    IRFTarget%K(i,j,k)=IRFSource%K(i,j,k)
                END DO
            END DO
        END DO
        DO i=1,IRFTarget%NtimeS
            DO j=1,IRFTarget%Nbeta
                DO k=1,IRFTarget%Nintegration
                    IRFTarget%KexcForce(i,j,k)=IRFSource%KexcForce(i,j,k)
                END DO
            END DO
        END DO

        END SUBROUTINE CopyTIRF
!       ---
        SUBROUTINE Initialize_IRF(IRF,Results,namefile)
        USE MResults
        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        TYPE(TResults) :: Results
        CHARACTER(LEN=*) :: namefile
        CHARACTER(LEN=20) :: lookfor
        CHARACTER(LEN=80) :: discard
        REAL :: dt,tf
        INTEGER :: i
        OPEN(10,FILE=namefile)
        READ(10,'(A20)') lookfor
        DO WHILE (lookfor.NE.'--- Post processing ')
            READ(10,'(A20,A)') lookfor,discard
        END DO
        READ(10,*) IRF%Switch,dt,tf
        CLOSE(10)
        IF (IRF%Switch.EQ.1) THEN
            IRF%Ntime=INT(tf/dt)
            IRF%NtimeS=2*INT(tf/dt)-1
            IRF%Nradiation=Results%Nradiation
            IRF%Nintegration=Results%Nintegration
            IRF%Nbeta=Results%Nbeta
            Allocate(IRF%beta(IRF%Nbeta))
            IRF%beta=Results%beta
            CALL CreateTIRF(IRF,IRF%Ntime,IRF%Nradiation,IRF%Nintegration,IRF%Nbeta,IRF%NtimeS)
            DO i=1,IRF%Ntime
                IRF%Time(i)=i*dt
            END DO
            DO i=1,IRF%NtimeS
                IRF%TimeS(i)=(-(IRF%Ntime-1)+i)*dt
            END DO


        END IF
        END SUBROUTINE  Initialize_IRF

!       ---

        SUBROUTINE Compute_IRF(IRF,Results)

        USE Constants, only: PI,II
        USE MResults

        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        TYPE(TResults) :: Results
        REAL,DIMENSION(Results%Nw) :: CM
        INTEGER :: i,j,k,l
        COMPLEX :: ExcForce_l,ExcForce_l1

        DO i=1,IRF%Ntime
            DO j=1,Results%Nradiation
                DO k=1,Results%Nintegration
                    IRF%K(i,j,k)=0.
                    DO l=1,Results%Nw-1
!                        IRF%K(i,j,k)=IRF%K(i,j,k)-REAL(Results%Force(l,j,k))*COS(2.*PI/Results%Period(l)*IRF%Time(i))*2.*PI*ABS(1./Results%Period(l)-1./Results%Period(l+1))
                        IRF%K(i,j,k)=IRF%K(i,j,k)+0.5*(Results%RadiationDamping(l,j,k)*COS(Results%w(l)*IRF%Time(i))+Results%RadiationDamping(l+1,j,k)*COS(Results%w(l+1)*IRF%Time(i)))*(Results%w(l+1)-Results%w(l))
                    END DO
                    IRF%K(i,j,k)=IRF%K(i,j,k)*2./PI
                END DO
            END DO
        END DO

        DO j=1,Results%Nradiation
            DO k=1,Results%Nintegration
                IRF%AddedMass(j,k)=0.
                DO l=1,Results%Nw
                    CM(l)=0.
                    DO i=1,IRF%Ntime-1
 !                       CM(l)=CM(l)+IRF%K(i,j,k)*SIN(2.*PI/Results%Period(l)*IRF%Time(i))*(IRF%Time(i+1)-IRF%Time(i))
                       CM(l)=CM(l)+0.5*(IRF%K(i,j,k)*SIN(Results%w(l)*IRF%Time(i))+IRF%K(i+1,j,k)*SIN(Results%w(l)*IRF%Time(i+1)))*(IRF%Time(i+1)-IRF%Time(i))
                    END DO
                    CM(l)=Results%AddedMass(l,j,k)+CM(l)/Results%w(l)
                    IRF%AddedMass(j,k)=IRF%AddedMass(j,k)+CM(l)
                END DO
                IRF%AddedMass(j,k)=IRF%AddedMass(j,k)/Results%Nw
            END DO
        END DO

         DO i=1,IRF%NtimeS
            DO j=1,Results%Nbeta
                DO k=1,Results%Nintegration
                    IRF%KexcForce(i,j,k)=0.
                    DO l=1,Results%Nw-1
                    ExcForce_l=Results%DiffractionForce(l,j,k)+Results%FroudeKrylovForce(l,j,k)
                    ExcForce_l1=Results%DiffractionForce(l+1,j,k)+Results%FroudeKrylovForce(l+1,j,k)
                    !integration of w(1)<=w<=w(end)
                    IRF%KexcForce(i,j,k)=IRF%KexcForce(i,j,k)+                           &
                            0.5*(ExcForce_l*EXP(-II*Results%w(l)*IRF%TimeS(i))            &
                                +ExcForce_l1*EXP(-II*Results%w(l+1)*IRF%TimeS(i)))        &
                            *(Results%w(l+1)-Results%w(l))
                    !integration of -w(end)<=w<=-w(1)
                    ExcForce_l=Conjg(ExcForce_l)
                    ExcForce_l1=Conjg(ExcForce_l1)
                    IRF%KexcForce(i,j,k)=IRF%KexcForce(i,j,k)+                           &
                            0.5*(ExcForce_l*EXP(II*Results%w(l)*IRF%TimeS(i))           &
                                +ExcForce_l1*EXP(II*Results%w(l+1)*IRF%TimeS(i)))       &
                            *(Results%w(l+1)-Results%w(l))
                    END DO
                    IRF%KexcForce(i,j,k)=IRF%KexcForce(i,j,k)/2/PI
                END DO
            END DO
        END DO


        END SUBROUTINE Compute_IRF

!       ---

        SUBROUTINE Save_IRF(IRF,namefile)
        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        CHARACTER(LEN=*) :: namefile
        INTEGER :: i,j,k
        OPEN(10,FILE=namefile)
        WRITE(10,*) 'VARIABLES="Time (s)"'
        DO k=1,IRF%Nintegration
            WRITE(10,'(A,I4,A,I4,A)') '"AddedMass ',k,'" "IRF ',k,'"'
        END DO
        DO j=1,IRF%Nradiation
            WRITE(10,'(A,I4,A,I6,A)') 'Zone t="DoF ',j,'",I=',IRF%Ntime,',F=POINT'
            DO i=1,IRF%Ntime
                WRITE(10,'(80(X,E14.7))') IRF%Time(i),(IRF%AddedMass(j,k),IRF%K(i,j,k),k=1,IRF%Nintegration)
            END DO
        END DO
        CLOSE(10)
        END SUBROUTINE Save_IRF

        SUBROUTINE Save_IRF_excForce(IRF,namefile)
        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        CHARACTER(LEN=*) :: namefile
        INTEGER :: i,j,k
        OPEN(10,FILE=namefile)
        WRITE(10,*) 'VARIABLES="Time (s)"'
		    DO k=1,IRF%Nintegration
            WRITE(10,'(A,I4,A)') '"IRF ',k,'"'
        END DO
        IF ((1+IRF%Nintegration).GT.10000) THEN
          PRINT*,'Increase the value used in the following write command'
          STOP
        ENDIF
        DO j=1,IRF%Nbeta
            WRITE(10,'(A,E14.7,A,I6,A)') 'Zone t="beta ',IRF%beta(j),'",I=',IRF%NtimeS,',F=POINT'
            DO i=1,IRF%NtimeS
                WRITE(10,'(10000(X,E14.7))') IRF%TimeS(i),(REAL(IRF%KexcForce(i,j,k)),k=1,IRF%Nintegration)
            END DO
        END DO
        CLOSE(10)
        END SUBROUTINE Save_IRF_excForce

!       ---
        SUBROUTINE DeleteTIRF(IRF)
        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        DEALLOCATE(IRF%Time,IRF%K,IRF%AddedMass)
        END SUBROUTINE DeleteTIRF
!       ---
END MODULE
