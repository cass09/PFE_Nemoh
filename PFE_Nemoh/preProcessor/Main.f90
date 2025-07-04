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
!   - R. Kurnia
!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - preProcessor - January 2014
!
!--------------------------------------------------------------------------------------
!
    PROGRAM Main
!
    USE Constants
    USE Elementary_functions
    USE MEnvironment
    USE MIdentification
    USE MNemohCal,              ONLY:TNemCal,READ_TNEMOHCAL,IdFreqHz,IdPeriod,Discretized_Omega_and_Beta
    USE MMesh
    USE BodyConditions
    USE Integration
!
    IMPLICIT NONE
!
    TYPE(TID) :: ID                     ! Calculation identification data
    TYPE(TMesh) :: Mesh                 ! Mesh data
    TYPE(TEnvironment) :: Environment   ! Environment data
    TYPE(TNemCal)      :: inpNEMOHCAL
!   Wave frequencies
    INTEGER :: Nw
    REAL,DIMENSION(:),ALLOCATABLE :: w
!   Radiation cases
    INTEGER :: Nradiation
    COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: NormalVelocity
!   Diffraction cases
    INTEGER :: Nbeta
    REAL,DIMENSION(:),ALLOCATABLE :: beta
    COMPLEX,DIMENSION(:),ALLOCATABLE :: Pressure
!   Force integration cases
    INTEGER :: Switch_Potential
    INTEGER :: Nintegration
    REAL,DIMENSION(:),ALLOCATABLE :: NDS
    REAL,DIMENSION(:,:),ALLOCATABLE :: FNDS
!   Froude Krylov forces
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: FKforce
    REAL,DIMENSION(4,3) :: P
    INTEGER  :: IMN,jj
    REAL :: ZMN
!   Free surface visualisation
    INTEGER :: Switch_FreeSurface
    INTEGER :: Nx,Ny
    REAL :: Lx,Ly
!   Kochin function
    INTEGER :: Switch_Kochin
    INTEGER :: NTheta
    REAL :: Thetamin,Thetamax
!   Source distribution
    INTEGER :: Switch_SourceDistr
!   Other local variables
    INTEGER :: M
    INTEGER :: i,j,c,k,IdBody,IdMode,indsum, num_Body
!   CML - MODIF
!   Interaction Theory - Cylinder envelop
    INTEGER :: Switch_Cylsurface
    INTEGER :: cyldTheta, cyldZ
    REAL    :: cylR, cylZ
    INTEGER :: run_IT, run_BEM, ITsources
    INTEGER :: IT_Nb, IT_NDir
    REAL    :: DirMin,DirMax
    REAL,DIMENSION(:,:), ALLOCATABLE :: IT_coord
!
!   --- Initialize and read input datas ----------------------------------------------------------------------------------------
!
    CALL ReadTID(ID)
    CALL ReadTMesh(Mesh,ID)
    CALL READ_TNEMOHCAL(ID,InpNEMOHCAL)
!   ----------- passing inputs ------------------------------------------
    Environment =InpNEMOHCAL%Env
    Nradiation  =InpNEMOHCAL%Nradtot
    Nintegration=InpNEMOHCAL%Nintegtot
    Nw          =InpNEMOHCAL%waveinput%NFreq
    Nbeta          =InpNEMOHCAL%waveinput%NBeta
    CALL Discretized_Omega_and_Beta(0, InpNEMOHCAL%waveinput, w, beta)
    Switch_Potential  =InpNEMOHCAL%OptOUTPUT%Switch_POTENTIAL
    Switch_Kochin     =InpNEMOHCAL%OptOUTPUT%Kochin%Switch
    Ntheta            =InpNEMOHCAL%OptOUTPUT%Kochin%Ntheta
    thetamin          =InpNEMOHCAL%OptOUTPUT%Kochin%min_theta
    thetamax          =InpNEMOHCAL%OptOUTPUT%Kochin%max_theta
    Switch_FreeSurface=InpNEMOHCAL%OptOUTPUT%Freesurface%Switch
    Nx                =InpNEMOHCAL%OptOUTPUT%Freesurface%Nx
    Ny                =InpNEMOHCAL%OptOUTPUT%Freesurface%Ny
    Lx                =InpNEMOHCAL%OptOUTPUT%Freesurface%Lx
    Ly                =InpNEMOHCAL%OptOUTPUT%Freesurface%Ly
    Switch_SourceDistr=InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr
    Switch_Cylsurface =InpNEMOHCAL%IntTheory%Cylsurface%Switch
    cylR              =InpNEMOHCAL%IntTheory%Cylsurface%cylR
    cyldTheta         =InpNEMOHCAL%IntTheory%Cylsurface%cyldTheta
    cyldZ             =InpNEMOHCAL%IntTheory%Cylsurface%cyldZ
    ITsources         =InpNEMOHCAL%IntTheory%ITsources       
    run_IT=InpNEMOHCAL%IntTheory%run_IT
    IF (run_IT .EQ. 1) THEN 
        run_BEM=InpNEMOHCAL%IntTheory%run_BEM
        IT_Nb=InpNEMOHCAL%IntTheory%Nb
        IT_NDir=InpNEMOHCAL%IntTheory%NDir
        DirMin=InpNEMOHCAL%IntTheory%DirMin
        DirMax=InpNEMOHCAL%IntTheory%DirMax
        ALLOCATE(IT_coord(IT_Nb, 2))
        IT_coord=InpNEMOHCAL%IntTheory%Bcoord
    END IF 

! ---------------------------------------------------------------------------
!   Print summary of calculation case
    WRITE(*,*) ' '
    WRITE(*,*) ' Summary of calculation'
    WRITE(*,*) ' '
    IF (Environment%Depth.GT.0.) THEN
        WRITE(*,'(A,F7.2,A)') '  ->  Water depth = ',Environment%Depth,' m'
    ELSE
        WRITE(*,'(A)') '  ->  Infinite water depth'
    END IF
    WRITE(*,'(A,I5,A,F7.4,A,F7.4)') '  ->',Nw,' wave frequencies from ',w(1),' to ',w(Nw)
    WRITE(*,'(A,I5,A,F7.4,A,F7.4)') '  ->',Nbeta,' wave directions from  ',beta(1),' to ',beta(Nbeta)
    WRITE(*,'(A,I5,A)') '  ->',Nradiation,' radiation problems'
    WRITE(*,'(A,I5,A)') '  ->',Nintegration,' forces'
    IF (Mesh%Isym==1) THEN
    WRITE(*,'(A,I5)') '  ->  Half-body mesh (symmetric) with Npanels=', Mesh%Npanels
    ELSE
    WRITE(*,'(A,I5)') '  ->  Full-body mesh with Npanels=', Mesh%Npanels
    ENDIF

    WRITE(*,*) ' '
!
!   --- Generate force integration file ----------------------------------------------------------------------------------------
!
    ALLOCATE(FNDS(Nintegration,Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(NDS(Mesh%Npanels*2**Mesh%Isym))
    indsum=1
    DO IdBody=1,InpNEMOHCAL%Nbodies
        DO IdMode=1,InpNEMOHCAL%bodyinput(IdBody)%NIntegration
        CALL ComputeNDS(Mesh,IdBody,InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%ICase,&
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%Direction(1:3),           &
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%Axis(1:3),NDS)
                DO c=1,Mesh%Npanels*2**Mesh%Isym
                FNDS(indsum,c)=NDS(c)
                END DO
                indsum=indsum+1
        END DO
    END DO
    DEALLOCATE(NDS)
    OPEN(11,FILE=TRIM(ID%ID)//'/mesh/Integration.dat')
    WRITE(11,*) Nintegration
    DO j=1,Nintegration
        WRITE(11,*) (FNDS(j,c),c=1,Mesh%Npanels*2**Mesh%Isym)
    END DO
    CLOSE(11)
!
!   --- Generate body conditions and calculate FK forces ----------------------------------------------------------------------------------------
!
    ALLOCATE(NVEL(Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(PRESSURE(Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(FKForce(Nw,Nbeta,Nintegration))
    ALLOCATE(NormalVelocity(Mesh%Npanels*2**Mesh%Isym,(Nbeta+Nradiation)*Nw))
    DO i=1,Nw
        DO j=1,Nbeta
            CALL ComputeDiffractionCondition(ITsources, Mesh,w(i),Beta(j), j-int((Nbeta-1)/2)-1, Environment,PRESSURE,NVEL)
            DO c=1,Mesh%Npanels*2**Mesh%Isym
                NormalVelocity(c,j+(i-1)*(Nbeta+Nradiation))=NVEL(c)
            END DO
!           Calculate the corresponding FK forces
            DO k=1,Nintegration
                FKForce(i,j,k)=0.
                DO c=1,Mesh%nPanels*2**Mesh%Isym
                    FKForce(i,j,k)=FKForce(i,j,k)-PRESSURE(c)*FNDS(k,c)
                END DO
            END DO
        END DO
        indsum=1
        DO IdBody=1,InpNEMOHCAL%Nbodies
                DO IdMode=1,InpNEMOHCAL%bodyinput(IdBody)%NRadiation
                CALL ComputeRadiationCondition(Mesh,IdBody,                     &
                        InpNEMOHCAL%bodyinput(IdBody)%RadCase(IdMode)%ICase,         &
                        InpNEMOHCAL%bodyinput(IdBody)%RadCase(IdMode)%Direction(1:3),&
                        InpNEMOHCAL%bodyinput(IdBody)%RadCase(IdMode)%Axis(1:3),NVEL)
                     DO c=1,Mesh%Npanels*2**Mesh%Isym
                        NormalVelocity(c,indsum+Nbeta+(i-1)*(Nbeta+Nradiation))=NVEL(c)
                     END DO
                indsum=indsum+1
                END DO
        END DO
    END DO
    CLOSE(11)
    DEALLOCATE(PRESSURE,NVEL,FNDS)
!
    IF (run_IT .EQ. 1 .AND. InpNEMOHCAL%Nbodies .NE. 1) THEN 
        run_IT = 0
        WRITE(*,*) "BEM resolution must be on only ONE body to apply Interaction Theory"
    END IF

    IF (run_IT .EQ. 1 .AND. Switch_Cylsurface .NE. 1) THEN 
        run_IT = 0
        WRITE(*,*) "Cylindrical mesh necessary to apply Interaction Theory"
    END IF

!   --- Save body conditions ----------------------------------------------------------------------------------------
!
    OPEN(11,FILE=TRIM(ID%ID)//'/Normalvelocities.dat')
    WRITE(11,*) (Nbeta+Nradiation)*Nw
    WRITE(11,*) ((w(i),j=1,Nbeta+Nradiation),i=1,Nw)
    DO i=1,Nw
      WRITE(11,*) (/ (DIFFRACTION_PROBLEM, j=1,Nbeta), (RADIATION_PROBLEM, j=1,Nradiation) /)
    ENDDO
    WRITE(11,*) ((Switch_Potential,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) ((Switch_Freesurface,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) ((Switch_Kochin,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) ((Switch_SourceDistr,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) ((Switch_Cylsurface,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) ((ITsources,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) run_IT, run_BEM
    DO c=1,Mesh%Npanels*2**Mesh%Isym
        WRITE(11,*) (REAL(NormalVelocity(c,j)),IMAG(NormalVelocity(c,j)),j=1,(Nbeta+Nradiation)*Nw)
    END DO
    CLOSE(11)
    DEALLOCATE(NormalVelocity)
!
!   --- Save FK forces ----------------------------------------------------------------------------------------
!
    OPEN(10,FILE=TRIM(ID%ID)//'/results/FKForce.tec')
    WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
    indsum=1
    DO IdBody=1,InpNEMOHCAL%Nbodies
        DO IdMode=1,InpNEMOHCAL%bodyinput(IdBody)%NIntegration
        WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"abs(F',IdBody,indsum,')" "angle(F',IdBody,indsum,')"'
        indsum=indsum+1
        END DO
    END DO
    DO c=1,Nbeta
        WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="FKforce - beta = ',beta(c)*180./PI,'",I=',Nw,',F=POINT'
        DO i=1,Nw
            WRITE(10,'(80(X,E14.7))') w(i),(ABS(FKForce(i,c,k)),ATAN2(IMAG(FKForce(i,c,k)),REAL(FKForce(i,c,k))),k=1,Nintegration)
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=TRIM(ID%ID)//'/results/FKForce.dat')
    DO k=1,Nintegration
        WRITE(10,*) ((ABS(FKForce(i,c,k)),ATAN2(IMAG(FKForce(i,c,k)),REAL(FKForce(i,c,k))),c=1,Nbeta),(0.*c,0.*c,c=1,Nradiation),i=1,Nw)
    END DO
    CLOSE(10)
    DEALLOCATE(FKForce)
!
!   --- Generate Free Surface visualisation file ----------------------------------------------------------------------
!
    OPEN(11,FILE=TRIM(ID%ID)//'/mesh/Freesurface.dat')
    WRITE(11,*) Nx*Ny,(Nx-1)*(Ny-1)
    DO i=1,Nx
        DO j=1,Ny
            WRITE(11,'(3(X,E14.7))') -0.5*Lx+Lx*(i-1)/(Nx-1),-0.5*Ly+Ly*(j-1)/(Ny-1),0.
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(11,'(4(X,I7))') j+(i-1)*Ny,j+1+(i-1)*Ny,j+1+i*Ny,j+i*Ny
        END DO
    END DO
    CLOSE(11)
!
!   --- Generate Kochin file ----------------------------------------------------------------------------------------
!
    OPEN(11,FILE=TRIM(ID%ID)//'/mesh/Kochin.dat')
    WRITE(11,*) NTheta
    IF (Ntheta.GT.0) THEN
        IF (NTheta.GT.1) THEN
            DO j=1,NTheta
                WRITE(11,*) (Thetamin+(Thetamax-Thetamin)*(j-1)/(NTheta-1))*PI/180.
            END DO
        ELSE
            WRITE(11,*) Thetamin*PI/180.
        END IF
    END IF
    CLOSE(11)
!
! --- Generate Cylindrical control surf mesh file ----------------------------------------------------------------------
!   CML
    IF (Switch_Cylsurface .EQ. 1) THEN 
        OPEN(11,FILE=ID%ID(1:ID%lID)//'/mesh/Cylsurface.dat')
        WRITE(11,*) cyldTheta*cyldZ,cyldTheta*(cyldZ-1) 
        DO i=1,cyldZ
            DO j=1,cyldTheta
                IF (cyldZ .EQ. 1) THEN
                    cylZ = 0
                ELSE
                    cylZ = -Environment%Depth*(1.-COS(PI/2.*(i-1.)/(cyldZ-1.)))
                END IF
                WRITE(11,'(3(X,E14.6))') cylR*COS(2.*PI*(j-1)/cyldTheta),cylR*SIN(2.*PI*(j-1)/cyldTheta),cylZ
            END DO
        END DO  
        DO i=1,cyldZ-1
            DO j=1,cyldTheta-1
                WRITE(11,'(4(X,I7))') j+(i-1)*cyldTheta,j+i*cyldTheta,j+i*cyldTheta+1,j+(i-1)*cyldTheta+1
            END DO
            WRITE(11,'(4(X,I7))') i*cyldTheta,(i+1)*cyldTheta,i*cyldTheta+1,(i-1)*cyldTheta+1
        END DO
        CLOSE(11)
    END IF 
!
!   CML MODIF : transmission from preProcessor to Solver
!   --- Save Interaction Theory Inputs ----------------------------------------------------------------------------------------------
!   
    IF (run_IT .EQ. 1) THEN
        WRITE(*,*) "--------- Interaction Theory activated ------------"
        OPEN(12,FILE=TRIM(ID%ID)//'/input_IT.dat')
        WRITE(12,*) run_IT
        WRITE(12,*) run_BEM
        WRITE(12,*) IT_Nb
        DO num_Body=1,IT_Nb
            WRITE(12,*) IT_coord(num_Body, :)
        END DO
        WRITE(12,*) IT_NDir, DirMin, DirMax
        WRITE(12,*) cylR, cyldTheta, cyldZ
        CLOSE(12)
    ELSE 
        IF (Switch_Cylsurface .EQ. 1) THEN
            WRITE(*,*) "--------- Computing cylindrical mesh ------------"
        END IF
        IF (ITsources .EQ. 1) THEN
            WRITE(*,*) "--------- Computing sources for IT ------------"
        END IF
        WRITE(*,*) "--------- Interaction Theory NOT activated ------------"
    END IF
    
! 
!   --- Save index of cases ----------------------------------------------------------------------------------------------
!
    OPEN(10,FILE=TRIM(ID%ID)//'/results/index.dat')
    WRITE(10,*) Nw,Nbeta,Nradiation,Nintegration,Ntheta
    WRITE(10,*) '--- Force ---'
    indsum=1
    DO IdBody=1,InpNEMOHCAL%Nbodies
        DO IdMode=1,InpNEMOHCAL%bodyinput(IdBody)%NIntegration
        WRITE(10,*) indsum,IdBody,IdMode
        indsum=indsum+1
        END DO
    END DO
    WRITE(10,*) '--- Motion ---'
    indsum=1
    DO IdBody=1,InpNEMOHCAL%Nbodies
        DO IdMode=1,InpNEMOHCAL%bodyinput(IdBody)%NRadiation
        WRITE(10,*) indsum,IdBody,IdMode
        indsum=indsum+1
        END DO
    END DO

    WRITE(10,*) (Beta(k),k=1,Nbeta)
    WRITE(10,*) (w(k),k=1,Nw)
    WRITE(10,*) ((Thetamin+(Thetamax-Thetamin)*(k-1)/(NTheta-1))*PI/180.,k=1,Ntheta)
    CLOSE(10)
!
!   --- Finalize ----------------------------------------------------------------------------------------------------
!
    DEALLOCATE(w,Beta)
    DO IdBody=1,InpNEMOHCAL%Nbodies
      DEALLOCATE(inpNEMOHCAL%bodyinput(IdBody)%RadCase)
      DEALLOCATE(inpNEMOHCAL%bodyinput(IdBody)%IntCase)
    ENDDO
      DEALLOCATE(inpNEMOHCAL%bodyinput)
 !
    END PROGRAM Main
