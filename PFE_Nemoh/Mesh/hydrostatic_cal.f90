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
!	  - A. Babarit
!   - R Kurnia
!
!--------------------------------------------------------------------------------------
!
!  NEMOH v1.0 - Hydrostatic Calculation - July 2021
!
!--------------------------------------------------------------------------------------

 PROGRAM Hydrostatic

    USE MEnvironment
    USE MIdentification
    USE MMesh
    USE MNemohCal,              ONLY:TNemCal,READ_TNEMOHCAL

#ifndef GNUFORT
    USE iflport
#endif

    IMPLICIT NONE
    TYPE(TID) :: ID,DSCRPT              ! Calculation identification data
    TYPE(TMesh) :: Mesh                 ! Mesh data
    TYPE(TEnvironment) :: Environment   ! Environment data
    TYPE(TNemCal)      :: inpNEMOHCAL

!   Maillage proprement dit
        INTEGER,PARAMETER :: NFMX=20000 ! Nombre de facettes max
        INTEGER,PARAMETER :: NPMX=20000 ! Nombre de points max
        INTEGER :: Nmailmx              ! Nombre de facettes du maillage std max
!   Maillage du corps
        INTEGER :: NF,NP
        INTEGER,DIMENSION(4,NFMX) :: Facette
	REAL,DIMENSION(NPMX) :: X,Y,Z
!   Partie immergee du maillage
        INTEGER :: NFm,NPm
        INTEGER,DIMENSION(4,NFMX) :: Facettem
	REAL,DIMENSION(NPMX) :: Xm,Ym,Zm
!   Calcul hydrostatique
	REAL DEPLACEMENT,XF,YF,ZF,SF
	REAL,DIMENSION(6,6) :: KH
	REAL,DIMENSION(:,:), ALLOCATABLE :: KHs, Mass
	REAL :: xG,yG,zG
	REAL :: RHO,G
!   Calcul coque
	REAL,DIMENSION(3,3) :: Icoque
	REAL,DIMENSION(3) :: Gcoque,CDG

        INTEGER         :: i,j,k,Nsym
        INTEGER         :: IdMode, IdBody
        LOGICAL :: ex
        INTEGER :: cdir

        CHARACTER(len=30)  :: fmt  
        INTEGER :: Nbodies
        REAL, DIMENSION(:,:), ALLOCATABLE :: COORD
!
!   --- Initialize and read input datas ----------------------------------------------------------------------------------------
!
    CALL ReadTID(ID)
    CALL ReadTMesh(Mesh,TRIM(ID%ID)//'/mesh/')
    CALL ReadTEnvironment(Environment,TRIM(ID%ID)//'/Nemoh.cal')
    RHO=Environment%RHO
    G=Environment%G
    NP=Mesh%Npoints
    NF=Mesh%Npanels
    Nsym=Mesh%Isym
    CALL READ_TNEMOHCAL(ID,InpNEMOHCAL)
!     OPEN(10,FILE=TRIM(ID%ID)//'/Mesh.cal')
!     READ(10,*) DSCRPT%ID
!     DSCRPT%lID=LNBLNK(DSCRPT%ID)
!     READ(10,*)
!     READ(10,*)
!     READ(10,*) xG,yG,zG
!     READ(10,*)
!     CLOSE(10)
! CML modif reading CoG in Nemoh.cal for BEM computation
!                       in Hydro.txt for PIT computation
        IF (InpNEMOHCAL%IntTheory%Cylsurface%Switch==1) THEN
                OPEN(10, FILE=TRIM(ID%ID)//'/Hydro.txt')
                READ(10, *) Nbodies
                ALLOCATE(COORD(Nbodies, 3))
                DO IdBody=1, Nbodies
                        READ(10, *) COORD(IdBody, :)
                END DO
        ELSE  
                Nbodies=InpNEMOHCAL%Nbodies
        END IF
        IdMode=4
        ALLOCATE(KHs(Nbodies*6, Nbodies*6))
        ALLOCATE(Mass(Nbodies*6, Nbodies*6))
        KHs=0
        Mass=0
        ! CML modif : NP and NF in L10.dat with BEM are for Nbodies
        NP=NP/InpNEMOHCAL%Nbodies
        NF=NF/InpNEMOHCAL%Nbodies
        WRITE(fmt, '(A,I0,A)') '(', Nbodies*6, '(1X,E14.7))'
    DO IdBody=1,Nbodies
        IF (InpNEMOHCAL%IntTheory%Cylsurface%Switch==1) THEN
                xG=COORD(IdBody, 1)
                yG=COORD(IdBody, 2)
                zG=COORD(IdBody, 3)
        ELSE  
                xG=InpNEMOHCAL%bodyinput(IdBody)%RadCase(IdMode)%Axis(1)
                yG=InpNEMOHCAL%bodyinput(IdBody)%RadCase(IdMode)%Axis(2)
                zG=InpNEMOHCAL%bodyinput(IdBody)%RadCase(IdMode)%Axis(3)
        END IF
        write(*,*) "body =", IdBody, "Center of gravity: ", xG, yG, zG
        DO j=1,NP
                X(j)=Mesh%X(1,j)-xG
                Y(j)=Mesh%X(2,j)-yG
                Z(j)=Mesh%X(3,j)
        END DO
        DO j=1,NF
                DO i=1,4
                FACETTE(i,j)=Mesh%P(i,j);
                END DO
        END DO
        CALL HYDRO(X,Y,Z,NP,FACETTE,NF,DEPLACEMENT,XF,YF,ZF,SF,KH,Xm,Ym,Zm,NPm,FACETTEm,NFm,RHO,G)
        DO j=1,NP
                        X(j)=X(j)+xG
                        Y(j)=Y(j)+yG
        END DO
        IF (Nsym.EQ.1) THEN
                DEPLACEMENT=2.0*DEPLACEMENT
                YF=0.
                SF=2.0*SF
                KH(3,3)=2.*KH(3,3)
                KH(3,4)=0.
                KH(4,3)=0.
                KH(3,5)=2.*KH(3,5)
                KH(5,3)=KH(3,5)
                KH(4,4)=2.*KH(4,4)
                KH(4,5)=0.
                KH(5,4)=0.
                KH(5,5)=2.*KH(5,5)
        END IF
                KH(4,4)=KH(4,4)+deplacement*RHO*G*(ZF-ZG)
                KH(5,5)=KH(5,5)+deplacement*RHO*G*(ZF-ZG)
        DO j = 1, 6
           DO k = 1, 6
                KHs((IdBody-1) * 6 + j, (IdBody-1) * 6 + k) = KH(j, k)
           END DO
        END DO
        write(*,*) ' -> Calculate hull mass and inertia '
        WRITE(*,*) ' '
    
        CDG(1)=xG
        CDG(2)=yG
        CDG(3)=zG
        CALL coque(X,Y,Z,NP,facette,NF,Deplacement,Icoque,Gcoque,CDG,Nsym,rho)
    
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/GC_hull.dat', ACTION='WRITE',POSITION='APPEND')
        WRITE(10,'(3(1X,E14.7))') Gcoque(1),Gcoque(2),Gcoque(3)
        CLOSE(10)
    
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Inertia_hull.dat', ACTION='WRITE',POSITION='APPEND')
        DO i=1,3
                WRITE(10,'(3(1X,E14.7))') (Icoque(i,j),j=1,3)
        END DO
        CLOSE(10)

        DO j = 1, 3
                Mass((IdBody-1) * 6 + j, (IdBody-1) * 6 + j) = DEPLACEMENT*RHO
        END DO
        DO j = 1, 3
           DO k = 1, 3
                Mass((IdBody-1) * 6 + j+3, (IdBody-1) * 6 + k+3) = Icoque(i,j)
           END DO
        END DO


        WRITE(*,'(A,I3)') '   - Coordinates of buoyancy centre '
        WRITE(*,'(A,F7.3,A)') '     XB = ',XF+xG,'m'
        WRITE(*,'(A,F7.3,A)') '     YB = ',YF+yG,'m'
        WRITE(*,'(A,F7.3,A)') '     ZB = ',ZF,'m'
        WRITE(*,'(A,E14.7,A)') '    - Displacement  = ',DEPLACEMENT,' m^3'
        WRITE(*,'(A,E14.7,A)') '    - Waterplane area  = ',SF, ' m^2'
        WRITE(*,'(A,E14.7,A)') '    - Mass =',DEPLACEMENT*RHO, ' Kg'
        WRITE(*,*) ' '

        IF ((ABS(XF).GT.1E-02).OR.(ABS(YF).GT.1E-02)) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' !!! WARNING !!! '
            WRITE(*,*) ' '
            WRITE(*,'(A,I3)') ' Buoyancy center and gravity center are not vertically aligned. '
            WRITE(*,*) ' This is not an equilibrium position.'
            WRITE(*,'(A,F7.3,1X,A,F7.3)') ' XF = ',XF+xG,' XG = ',xG
            WRITE(*,'(A,F7.3,1X,A,F7.3)') ' YF = ',YF+yG,' YG = ',yG
        END IF


        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Hydrostatics.dat', ACTION='WRITE',POSITION='APPEND')
        WRITE(10,'(A,F7.3,A,F7.3)') ' XF = ',XF+xG,' - XG = ',xG
        WRITE(10,'(A,F7.3,A,F7.3)') ' YF = ',YF+yG,' - YG = ',yG
        WRITE(10,'(A,F7.3,A,F7.3)') ' ZF = ',ZF,' - ZG = ',zG
        WRITE(10,'(A,E14.7)') ' Displacement = ',DEPLACEMENT
        WRITE(10,'(A,E14.7)') ' Waterplane area = ',SF
        WRITE(10,'(A,E14.7)') ' Mass =',DEPLACEMENT*RHO
        CLOSE(10)

    END DO
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/KH.dat')
        DO i=1,6*Nbodies
                WRITE(10,fmt) (KHs(i,j),j=1,6*Nbodies)
        END DO
        CLOSE(10)
        !INQUIRE (DIRECTORY=TRIM(ID%ID)//'/Mechanics', EXIST=ex) !this is Intel-specific
        INQUIRE (FILE=TRIM(ID%ID)//'/Mechanics/.', EXIST=ex)
        IF (.NOT.ex) cdir=SYSTEM('mkdir '//TRIM(ID%ID)//'/Mechanics')

        OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Inertia.dat')
        DO i=1,6*Nbodies
                WRITE(10,fmt) (Mass(i,j),j=1,6*Nbodies)
        END DO
        ! WRITE(10,'(6(1X,E14.7))') DEPLACEMENT*RHO,0.,0.,0.,0.,0.
        ! WRITE(10,'(6(1X,E14.7))') 0.,DEPLACEMENT*RHO,0.,0.,0.,0.
        ! WRITE(10,'(6(1X,E14.7))') 0.,0.,DEPLACEMENT*RHO,0.,0.,0.
        ! DO i=1,3
        !         WRITE(10,'(6(1X,E14.7))') 0.,0.,0.,(Icoque(i,j),j=1,3)
        ! END DO
        CLOSE(10)

        OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Kh.dat')
        DO i=1,6*Nbodies
                WRITE(10,fmt) (KHs(i,j),j=1,6*Nbodies)
        END DO
        CLOSE(10)

        INQUIRE (FILE=TRIM(ID%ID)//'/Mechanics/Badd.dat', EXIST=ex)
        IF (.NOT.ex) THEN
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Badd.dat')
        DO i=1,6*Nbodies
                WRITE(10,fmt) (0.,j=1,6*Nbodies)
        END DO
        CLOSE(10)
        ENDIF

        INQUIRE (FILE=TRIM(ID%ID)//'/Mechanics/Km.dat', EXIST=ex)
        IF (.NOT.ex) THEN
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Km.dat')
        DO i=1,6*Nbodies
                WRITE(10,fmt) (0.,j=1,6*Nbodies)
        END DO
        CLOSE(10)
        ENDIF


        DEALLOCATE(KHs, Mass)

end program Hydrostatic
