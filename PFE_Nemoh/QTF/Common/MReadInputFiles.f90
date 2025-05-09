!--------------------------------------------------------------------------------------
!
!     Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
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
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr) 2014
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!
!--------------------------------------------------------------------------------------
MODULE MReadInputFiles

USE CONSTANTS,                  ONLY: II, PI,ID_BODY,ID_FREESURFACE
USE MMesh,                      ONLY: TMesh,CreateTMesh
USE MFace,                      ONLY: TWLine
USE Elementary_functions,       ONLY: CROSS_PRODUCT
USE MFileDirectoryList
USE MNemohCal,                  ONLY: Tqtfinput

IMPLICIT NONE
PUBLIC:: Read_NP_GaussQuad,Read_Mechanical_Coefs,Read_FirstOrderLoad,Read_Motion, &
         Read_SourceDistribution,Read_Eps_Zmin
!
TYPE TMech
    REAL,ALLOCATABLE,DIMENSION(:,:) :: MassMat          !Mass-Inertia Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat         !Stifness-Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat_EXT     !Additional Stifness Matrix i.e: mooring
    REAL,ALLOCATABLE,DIMENSION(:,:) :: DampCoefMat_EXT  !Additional damping coefficients
END TYPE

TYPE TLoad1
    COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: excitation
    REAL,ALLOCATABLE,DIMENSION(:,:,:)    :: addedmass
    REAL,ALLOCATABLE,DIMENSION(:,:,:)    :: dampcoef
END TYPE
TYPE TSource
    COMPLEX,ALLOCATABLE,DIMENSION(:,:)   :: ZIGB,ZIGS       !source distribution
END TYPE
! Free surface Mesh
TYPE TVFACEFS
    REAL, ALLOCATABLE    :: N(:, :)     ! Normal vectors
    REAL, ALLOCATABLE    :: XM(:, :)    ! Centre of panels
ENDTYPE
TYPE TMeshFS
  TYPE(TMesh)                       :: Mesh
  TYPE(TWLine)                      :: BdyLine       ! line properties, XM, Length
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: BdyLineNormal ! normal of the line segment
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: BdyLineP      ! the line connectivity
  TYPE(TVFACEFS)                    :: VFace         ! this transposed version of the variables
  REAL                              :: Radius_Ext    ! Exterior free surface radius (used in Asymp)
  INTEGER                           :: NpointsR      ! Npoints of discretized interior radius
  INTEGER                           :: NBessel       ! number of bessel functions
END TYPE TMeshFS
!Potentials
TYPE TPotVel
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: TotPot,RadPot ! Total&radiation pot
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:)  :: TotVel,RadVel ! Total&radiation vel
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: IncPot        ! Incident/incoming Potential
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:)  :: IncVel        ! Velocity (unperturbed)
END TYPE
!
CONTAINS
        INTEGER FUNCTION Read_NP_GaussQuad(wd)
           CHARACTER(LEN=*) :: wd
           INTEGER       :: NPGQ
           OPEN(10,file=wd//'/input_solver.txt',form='formatted',status='old')
           READ(10,*) NPGQ
           CLOSE(10)
           Read_NP_GaussQuad=NPGQ**2
           RETURN
        END FUNCTION

        FUNCTION Read_Eps_Zmin(wd) result(EPS_ZMIN)
           CHARACTER(LEN=*) :: wd
           REAL       :: EPS_ZMIN
           OPEN(10,file=wd//'/input_solver.txt',form='formatted',status='old')
           READ(10,*)
           READ(10,*) EPS_ZMIN
           CLOSE(10)
           RETURN
        END FUNCTION


        SUBROUTINE Read_Mechanical_Coefs(wd,Nradiation,MechCoef)
           !input/output
           CHARACTER(LEN=*),            INTENT(IN)::wd
           INTEGER,                     INTENT(IN)::Nradiation
           TYPE(TMech),                 INTENT(OUT)::MechCoef
           !Local
           INTEGER      ::u1,u2,u3,u4,C,I,J

           ALLOCATE(MechCoef%MassMat(Nradiation,Nradiation))
           ALLOCATE(MechCoef%StiffMat(Nradiation,Nradiation))
           ALLOCATE(MechCoef%StiffMat_EXT(Nradiation,Nradiation))
           ALLOCATE(MechCoef%DampCoefMat_EXT(Nradiation,Nradiation))
           CALL exist_file(trim(wd)//'/Mechanics/Inertia.dat')
           CALL exist_file(trim(wd)//'/Mechanics/Kh.dat')
           CALL exist_file(trim(wd)//'/Mechanics/Km.dat')
           CALL exist_file(trim(wd)//'/Mechanics/Badd.dat')

           OPEN(NEWUNIT=u1,FILE=trim(wd)//'/Mechanics/Inertia.dat',ACTION='READ')
           OPEN(NEWUNIT=u2,FILE=trim(wd)//'/Mechanics/Kh.dat',ACTION='READ')
           OPEN(NEWUNIT=u3,FILE=trim(wd)//'/Mechanics/Km.dat',ACTION='READ')
           OPEN(NEWUNIT=u4,FILE=trim(wd)//'/Mechanics/Badd.dat',ACTION='READ')

           DO I=1,Nradiation
                READ(u1,*) (MechCoef%MassMat(I,J),J=1,Nradiation)
                READ(u2,*) (MechCoef%StiffMat(I,J),J=1,Nradiation)
                READ(u3,*) (MechCoef%StiffMat_EXT(I,J),J=1,Nradiation)
                READ(u4,*) (MechCoef%DampCoefMat_EXT(I,J),J=1,Nradiation)
           ENDDO
           CLOSE(u1)
           CLOSE(u2)
           CLOSE(u3)
           CLOSE(u4)
        END SUBROUTINE

        SUBROUTINE Read_FirstOrderLoad(wd,Nw,Nbeta,Nintegration,Nradiation,Forces1)
           !Input/Output
           CHARACTER(LEN=*),            INTENT(IN) :: wd
           INTEGER,                     INTENT(IN) :: Nw,Nbeta,Nintegration,Nradiation
           TYPE(TLoad1),                INTENT(OUT):: Forces1
           !Local
           INTEGER                                 :: u1,u2,u3,I,J,K
           REAL,DIMENSION(Nintegration)            :: Amp, Phase
           REAL                                    :: w0

           ALLOCATE(Forces1%excitation(Nw,Nintegration,Nbeta))
           ALLOCATE(Forces1%addedmass(Nw,Nradiation,Nradiation))
           ALLOCATE(Forces1%dampcoef(Nw,Nradiation,Nradiation))

           CALL exist_file(trim(wd)//'/results/CA.dat')
           CALL exist_file(trim(wd)//'/results/CM.dat')
           CALL exist_file(trim(wd)//'/results/Fe.dat')

           OPEN(NEWUNIT=u1,FILE=trim(wd)//'/results/CA.dat',ACTION='READ')!Damp. Coef.
           OPEN(NEWUNIT=u2,FILE=trim(wd)//'/results/CM.dat',ACTION='READ')!Added Mass
           READ(u1,*)
           READ(u2,*)
           DO I=1,Nw
                READ(u1,*)
                READ(u2,*)
                DO J=1,Nradiation
                   READ(u1,*) ( Forces1%dampcoef(I,J,K),K=1,Nradiation )
                   READ(u2,*) ( Forces1%addedmass(I,J,K),K=1,Nradiation )
                ENDDO
           ENDDO
           CLOSE(u1)
           CLOSE(u2)

           OPEN(NEWUNIT=u3,FILE=trim(wd)//'/results/Fe.dat',ACTION='READ')!excitation
           READ(u3,*)
           READ(u3,*)
           DO K=1,Nbeta
                READ(u3,*)
                DO I=1,Nw
                   READ(u3,*) w0,(Amp(J),J=1,Nintegration),(Phase(J),J=1,Nintegration)
                   DO J=1,Nintegration
                      Phase(J)=Phase(J)*PI/180.0
                      Forces1%excitation(I,J,K)=Amp(J)*CEXP(II*Phase(J))
                   ENDDO
                ENDDO
           ENDDO
           CLOSE(u3)
        END SUBROUTINE

        SUBROUTINE Read_Motion(wd,Nw,Nbeta,Nradiation,Motion)
           !Input/Output
           CHARACTER(LEN=*),            INTENT(IN) :: wd
           INTEGER,                     INTENT(IN) :: Nw,Nbeta,Nradiation
           COMPLEX,DIMENSION(Nw,Nradiation,Nbeta),INTENT(OUT):: Motion
           !Local
           INTEGER                                 :: u1,I,J,K
           REAL,DIMENSION(Nradiation)            :: Amp, Phase
           REAL                                    :: w0

           CALL exist_file(trim(wd)//'/Motion/RAO.dat')
           OPEN(NEWUNIT=u1,FILE=trim(wd)//'/Motion/RAO.dat',ACTION='READ')!RAO
           READ(u1,*)

           DO K=1,Nbeta
              READ(u1,*)
              DO I=1,Nw
                 READ(u1,*) w0,(Amp(J),J=1,Nradiation),(Phase(J),J=1,Nradiation)
                 DO J=1,Nradiation
                     Phase(J)=Phase(J)*PI/180.0
                     IF (MODULO(J,6)>3.OR.MODULO(J,6)==0) Amp(J)=Amp(J)*PI/180.0
                     Motion(I,J,K)=Amp(J)*CEXP(II*Phase(J))
                 ENDDO
              ENDDO
           ENDDO
           CLOSE(u1)
        END

        SUBROUTINE Read_SourceDistribution(wd,Iw,Nw,Nradiation,Nbeta,Npanels,SourceDistr)
           !INPUT/OUTPUT
           CHARACTER(LEN=*),            INTENT(IN) :: wd
           INTEGER,                     INTENT(IN) :: Nw,Nbeta,Nradiation,Npanels,Iw
           TYPE(TSource),               INTENT(INOUT):: SourceDistr
           !Local
           CHARACTER*5  :: str
           REAL         :: RE,IM
           INTEGER      :: idiffrad,Pbnumber,I,K
           INTEGER      :: u1
           idiffrad=0
           !Radiation problems
           DO K=1,Nradiation
                idiffrad= idiffrad+1
                Pbnumber=(Nbeta+Nradiation)*(Iw-1)+Nbeta+K
                WRITE(str,'(I0.5)') Pbnumber
                CALL exist_file(wd//'/results/sources/sources.'//str//'.dat')
                OPEN(NEWUNIT=u1,FILE=wd//'/results/sources/sources.'//str//'.dat')
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGB(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGS(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                CLOSE(u1)
           ENDDO

           !DIffraction problems
           DO K=1,Nbeta
                idiffrad= idiffrad+1
                Pbnumber=(Nbeta+Nradiation)*(Iw-1)+K
                WRITE(str,'(I0.5)') Pbnumber
                CALL exist_file(wd//'/results/sources/sources.'//str//'.dat')
                OPEN(NEWUNIT=u1,FILE=wd//'/results/sources/sources.'//str//'.dat')
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGB(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGS(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                CLOSE(u1)
           ENDDO
        END SUBROUTINE

        SUBROUTINE READ_POTENTIALS_VELOCITIES(wd,Nw,Nbeta,NRadiation,NPFlow,      &
                        datPotVel,w,k,beta,IDDATA)
          CHARACTER(LEN=*),                         INTENT(IN)    ::wd
          INTEGER,                                  INTENT(IN)    ::Nw,Nbeta,Nradiation,NPFlow
          INTEGER,                                  INTENT(IN)    ::IDDATA
          TYPE(TPotVel),                            INTENT(INOUT) ::datPotVel
          REAL,DIMENSION(Nbeta),                    INTENT(OUT)   ::beta
          REAL,DIMENSION(Nw),                       INTENT(OUT)   ::w,k    !radfreq&wavenumber
          INTEGER                     :: Iw, Ibeta, Irad,Ipanel
          REAL, DIMENSION(NPFlow)     :: REALDATA1,IMDATA1,REALDATA2,IMDATA2,               &
                                         REALDATA3,IMDATA3
          INTEGER                     :: uF1,uF2,uF3,uF4,uF5,uF6,ILINE
          REAL,DIMENSION(3)           :: wkbeta,wkIrad
          INTEGER                     :: RECLENGTH

          ALLOCATE(datPotVel%TotPot(NPFlow,Nbeta,Nw))
          ALLOCATE(datPotVel%TotVel(NPFlow,3,Nbeta,Nw))
          ALLOCATE(datPotVel%RadPot(NPFlow,Nradiation,Nw))
          ALLOCATE(datPotVel%RadVel(NPFlow,3,Nradiation,Nw))
          INQUIRE(iolength=RECLength) (wkbeta(1))

          IF (IDDATA==ID_BODY) THEN
            !body and waterline data
            OPEN(NEWUNIT=uF1, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotPotFILE,              &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+2*NPFLOW))
            OPEN(NEWUNIT=uF2, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotVelFILE,              &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+6*NPFLOW))
            OPEN(NEWUNIT=uF3, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadPotFILE,              &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+2*NPFLOW))
            OPEN(NEWUNIT=uF4, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadVelFILE,              &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+6*NPFLOW))
          ELSE
            ALLOCATE(datPotVel%IncPot(NPFlow,Nbeta,Nw))
            ALLOCATE(datPotVel%IncVel(NPFlow,3,Nbeta,Nw))

            !free-surface data
            OPEN(NEWUNIT=uF1, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotPotFILE_FS,           &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+2*NPFLOW))
            OPEN(NEWUNIT=uF2, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotVelFILE_FS,           &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+6*NPFLOW))
            OPEN(NEWUNIT=uF3, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadPotFILE_FS,           &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+2*NPFLOW))
            OPEN(NEWUNIT=uF4, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadVelFILE_FS,           &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+6*NPFLOW))
            OPEN(NEWUNIT=uF5, FILE=TRIM(wd)//'/'//PreprocDir//'/'//IncPotFILE_FS,           &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+2*NPFLOW))
            OPEN(NEWUNIT=uF6, FILE=TRIM(wd)//'/'//PreprocDir//'/'//IncVelFILE_FS,           &
                    STATUS='UNKNOWN',ACCESS='DIRECT',RECL=RECLength*(3+6*NPFLOW))
          ENDIF

          DO Iw=1,Nw
              DO Ibeta=1,Nbeta
                 ILINE=(Iw-1)*Nbeta+Ibeta
                 READ(uF1,REC=ILINE) wkbeta,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow)
                 datPotVel%TotPot(:,Ibeta,Iw)=CMPLX(REALDATA1,IMDATA1)
                 READ(uF2,REC=ILINE) wkbeta,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow),                          &
                                (REALDATA2(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA2(Ipanel),Ipanel=1,NPFlow),                          &
                                (REALDATA3(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA3(Ipanel),Ipanel=1,NPFlow)
                 datPotVel%TotVel(:,1,Ibeta,Iw)=CMPLX(REALDATA1,IMDATA1)
                 datPotVel%TotVel(:,2,Ibeta,Iw)=CMPLX(REALDATA2,IMDATA2)
                 datPotVel%TotVel(:,3,Ibeta,Iw)=CMPLX(REALDATA3,IMDATA3)
                 IF (Iw==1) beta(Ibeta)=wkbeta(3)
                 w(Iw)=wkbeta(1)
                 k(Iw)=wkbeta(2)

                 IF (IDDATA==ID_FREESURFACE) THEN
                         READ(uF5,REC=ILINE) wkbeta,                                        &
                                        (REALDATA1(Ipanel),Ipanel=1,NPFlow),                &
                                        (IMDATA1(Ipanel),Ipanel=1,NPFlow)
                         datPotVel%IncPot(:,Ibeta,Iw)=CMPLX(REALDATA1,IMDATA1)
                         READ(uF6,REC=ILINE) wkbeta,                                        &
                                        (REALDATA1(Ipanel),Ipanel=1,NPFlow),                &
                                        (IMDATA1(Ipanel),Ipanel=1,NPFlow),                  &
                                        (REALDATA2(Ipanel),Ipanel=1,NPFlow),                &
                                        (IMDATA2(Ipanel),Ipanel=1,NPFlow),                  &
                                        (REALDATA3(Ipanel),Ipanel=1,NPFlow),                &
                                        (IMDATA3(Ipanel),Ipanel=1,NPFlow)
                         datPotVel%IncVel(:,1,Ibeta,Iw)=CMPLX(REALDATA1,IMDATA1)
                         datPotVel%IncVel(:,2,Ibeta,Iw)=CMPLX(REALDATA2,IMDATA2)
                         datPotVel%IncVel(:,3,Ibeta,Iw)=CMPLX(REALDATA3,IMDATA3)
                 ENDIF
              ENDDO

              DO Irad=1,NRadiation
                 ILINE=(Iw-1)*Nradiation+Irad
                 READ(uF3,REC=ILINE) wkIrad,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow)
                 datPotVel%RadPot(:,Irad,Iw)=CMPLX(REALDATA1,IMDATA1)
                 READ(uF4,REC=ILINE) wkIrad,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow),                          &
                                (REALDATA2(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA2(Ipanel),Ipanel=1,NPFlow),                          &
                                (REALDATA3(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA3(Ipanel),Ipanel=1,NPFlow)
                 datPotVel%RadVel(:,1,Irad,Iw)=CMPLX(REALDATA1,IMDATA1)
                 datPotVel%RadVel(:,2,Irad,Iw)=CMPLX(REALDATA2,IMDATA2)
                 datPotVel%RadVel(:,3,Irad,Iw)=CMPLX(REALDATA3,IMDATA3)
              ENDDO
          ENDDO
          CLOSE(uF1)
          CLOSE(uF2)
          CLOSE(uF3)
          CLOSE(uF4)
          IF (IDDATA==ID_FREESURFACE) THEN
          CLOSE(uF5)
          CLOSE(uF6)
          ENDIF

        END SUBROUTINE

        SUBROUTINE READ_GENERALIZED_NORMAL_BODY_dAREA(wd,Npanels,Nintegration,genNormal_dS)
          CHARACTER(LEN=*),                     INTENT(IN)   :: wd
          INTEGER,                              INTENT(IN)   :: Npanels,NIntegration
          REAL, DIMENSION(Nintegration,Npanels),INTENT(INOUT):: genNormal_dS
          !Local
          INTEGER I,J,Ninteg,u
          CALL exist_file(TRIM(wd)//'/mesh/Integration.dat')

          OPEN(NEWUNIT=u, FILE=TRIM(wd)//'/mesh/Integration.dat', STATUS='OLD', ACTION='READ')
          READ(u,*) Ninteg
          IF (Ninteg.NE.Nintegration) THEN
            CLOSE(u)
            print*,'Number of rows (Nintegration) in /mesh/Integration.dat is not correct!'
            STOP
          ENDIF
          DO I = 1, Nintegration
             READ(u, *) (genNormal_dS(I,J), J=1,Npanels)
          END DO
          CLOSE(u)

        END SUBROUTINE


        SUBROUTINE READ_PREPARE_FREESURFACE_MESH(wd,MeshFS,QTFinputNem)
          CHARACTER(LEN=*),                       INTENT(IN):: wd
          TYPE(Tqtfinput),                        INTENT(IN):: QTFinputNem
          TYPE(TMeshFS),                          INTENT(INOUT):: MeshFS
          INTEGER           :: u,Ipoint,Ip,J,Ipanel,Iline
          INTEGER           :: Npoints,Npanels,NbdyLines
          REAL,DIMENSION(3) :: X1,X2,X3,X4 !panel nodes

          OPEN(NEWUNIT=u,FILE=TRIM(wd)//'/'//QTFinputNem%FSmeshfile)
          READ(u,*) MeshFS%Mesh%ISym,Npoints,Npanels,NbdyLines
          MeshFS%Radius_Ext=QTFinputNem%FSRe
          MeshFS%NpointsR=QTFinputNem%FSNRe
          MeshFS%NBessel=QTFinputNem%FSNBessel

          !memory allocation
          CALL CreateTMesh(MeshFS%Mesh,Npoints,Npanels,1)

          DO Ipoint=1,Npoints
             READ(u,*) Ip,(MeshFS%Mesh%X(J,Ipoint),J=1,3)
          ENDDO
             READ(u,*)
          ! Free surface panels
          DO Ipanel=1,Npanels
             READ(u,*) (MeshFS%Mesh%P(J,Ipanel),J=1,4)
             !Compute Normal Vect, Area and centre of mass of the panel
             !Quadrilateral Nodes
             X1=MeshFS%Mesh%X(:,MeshFS%Mesh%P(1,Ipanel))
             X2=MeshFS%Mesh%X(:,MeshFS%Mesh%P(2,Ipanel))
             X3=MeshFS%Mesh%X(:,MeshFS%Mesh%P(3,Ipanel))
             X4=MeshFS%Mesh%X(:,MeshFS%Mesh%P(4,Ipanel))

             CALL CALC_PANEL_PROPERTIES(X1,X2,X3,X4,    &
                     MeshFS%Mesh%N(:,Ipanel),           &
                     MeshFS%Mesh%A(Ipanel),             &
                     MeshFS%Mesh%XM(:,Ipanel) )
          ENDDO

          !Free surface boundary lines
          MeshFS%BdyLine%NWlineseg=NbdyLines
          ALLOCATE(MeshFS%BdyLine%SegLength(NbdyLines))
          ALLOCATE(MeshFS%BdyLine%XM(NbdyLines,3))
          ALLOCATE(MeshFS%BdyLineNormal(NbdyLines,2))!only Nx,Ny
          ALLOCATE(MeshFS%BdyLineP(NbdyLines,2))!only Nx,Ny
          DO Iline=1,Nbdylines
             READ(u,*) (MeshFS%BdyLineP(Iline,J),J=1,2)
             !Compute Normal Vect, Area and centre of mass of the panel
             !Quadrilateral Nodes
             X1=MeshFS%Mesh%X(:,MeshFS%BdyLineP(Iline,1))
             X2=MeshFS%Mesh%X(:,MeshFS%BdyLineP(Iline,2))

             CALL CALC_BDY_LINE_PROPERTIES(X1,X2,MeshFS%Radius_Ext, &
                     MeshFS%BdyLineNormal(Iline,:),                 &
                     MeshFS%BdyLine%SegLength(Iline),               &
                     MeshFS%BdyLine%XM(Iline,:) )
          ENDDO

             READ(u,*)
             CLOSE(u)
        !transposed version
        ALLOCATE(MeshFS%VFace%XM(Npanels,3))
        ALLOCATE(MeshFS%VFace%N(Npanels,3))
        MeshFS%VFace%XM=TRANSPOSE(MeshFS%Mesh%XM)
        MeshFS%VFace%N =TRANSPOSE(MeshFS%Mesh%N)
        END SUBROUTINE

        SUBROUTINE CALC_PANEL_PROPERTIES(X1,X2,X3,X4,N,A,XM)
            REAL,DIMENSION(3),INTENT(IN)   :: X1,X2,X3,X4
            REAL,DIMENSION(3),INTENT(OUT)  :: N,XM
            REAL             ,INTENT(OUT)  :: A

            REAL,DIMENSION(3)              :: U,V,N1,N2
            REAL                           :: A1,A2
            ! Area of 1-2-4 triangle
            U(:) = X2-X1
            V(:) = X4-X2
            N1=CROSS_PRODUCT(U, V)
            A1=0.5*NORM2(CROSS_PRODUCT(U, V))

            ! Area of 2-3-4 triangle
            U(:) = X4-X3
            V(:) = X2-X3
            N2=CROSS_PRODUCT(U, V)
            A2=0.5*NORM2(CROSS_PRODUCT(U, V))

            N=N1+N2
            N=N/NORM2(N)

            A=A1+A2

            XM= (X1+X2+X4)/3*A1/A  &
               +(X2+X3+X4)/3*A2/A
        END SUBROUTINE

        SUBROUTINE CALC_BDY_LINE_PROPERTIES(X1,X2,Rext,N,L,XM)
            REAL,DIMENSION(3),INTENT(IN)   :: X1,X2
            REAL,             INTENT(IN)   :: Rext
            REAL,DIMENSION(3),INTENT(OUT)  :: XM
            REAL,DIMENSION(2),INTENT(OUT)  :: N
            REAL             ,INTENT(OUT)  :: L
            REAL                           :: Rcalc
            XM=0.5*(X1+X2) !centre point
            L=NORM2(X2-X1) !length
            Rcalc=NORM2(XM)!distance to origin (0,0)
            IF (ABS(Rcalc-Rext).LT.0.01) THEN
            N(1)=(X2(2)-X1(2))/L
            N(2)=-(X2(1)-X1(1))/L
            ELSE
            N(1)=-(X2(2)-X1(2))/L
            N(2)=(X2(1)-X1(1))/L
            ENDIF

        END SUBROUTINE


        SUBROUTINE  exist_file(filename)
          CHARACTER(LEN=*),       INTENT(IN) :: filename
          LOGICAL                            ::existfile
          INQUIRE (FILE=filename, EXIST=existfile)
          IF (.NOT.existfile) THEN
               PRINT*,filename,' data is missing!'
               STOP
          ENDIF
        END SUBROUTINE


END MODULE
