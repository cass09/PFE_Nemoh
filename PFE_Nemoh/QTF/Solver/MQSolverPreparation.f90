!-----------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------------
!   Contributors list:
!    - Ruddy Kurnia (ECN)
!--------------------------------------------------------------------------------------
MODULE MQSolverPreparation

USE CONSTANTS
USE MEnvironment,               ONLY: TEnvironment,FunVect_inverseDispersion
USE Elementary_functions,       ONLY: Fun_closest,Fun_MIN,Fun_MAX,              &
                                      CROSS_PRODUCT_COMPLEX
USE MFileDirectoryList
USE MLogFile
USE linear_interpolation_module
USE MMesh,                      ONLY: Tmesh
USE MFace,                      ONLY: TVFace,TWLine
USE MNemohCal,                  ONLY: TNemCal
USE MReadInputFiles,            ONLY: TMech,TLoad1,TpotVel,TSource,              &
                                      Read_SourceDistribution
USE MCallInterp,                ONLY: FUNVECT_INTERP1_COMPLEX

IMPLICIT NONE

TYPE TQfreq
     INTEGER                           :: NwQ
     REAL, ALLOCATABLE,DIMENSION(:,:)  :: wQ,kQ          !freq for QTF comp
     REAL, ALLOCATABLE,DIMENSION(:,:)  :: diffwQ,sumwQ   !diff and sum freq
     INTEGER, ALLOCATABLE,DIMENSION(:) :: InterpPotSwitch
END TYPE

TYPE TASYMP
     INTEGER                           :: NR,NBESSEL ! Number of points R, Nummber of bessel modes
     REAL                              :: dR         !Rf(2)-Rf(1)
     REAL,ALLOCATABLE,DIMENSION(:)     :: Rf         !finite radius of free surface domain
END TYPE

TYPE TSourceQ
COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)   :: ZIGB_Per,ZIGS_Per  !Perturbed source distribution
COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)   :: ZIGB_Rad,ZIGS_Rad  !Radiation source distribution
END TYPE

CONTAINS



SUBROUTINE PREPARE_POTENTIAL_VELOCITIES(Qfreq,Nw,w,Nbeta,beta, &
                                        NPFlow,Nradiation,datPotVel,datPotVelQ,ID_DATA)
        !input/output
        TYPE(TQfreq),                   INTENT(IN)   ::Qfreq
        INTEGER,                        INTENT(IN)   ::Nw,Nbeta
        INTEGER,                        INTENT(IN)   ::NPFlow,Nradiation
        REAL, DIMENSION(Nw),            INTENT(IN)   ::w
        REAL, DIMENSION(Nbeta),         INTENT(IN)   ::beta
        INTEGER,                        INTENT(IN)   ::ID_DATA
        TYPE(TPotVel),                  INTENT(INOUT)::datPotVel
        TYPE(TPotVel),                  INTENT(INOUT)::datPotVelQ

        !Local
        INTEGER                    :: Ibeta,Iw1,Iw2,Ipanel,iflag,IwQ,Ixyz

        ALLOCATE(datPotVelQ%TotPot(NPFlow,  Nbeta     ,Qfreq%NwQ))
        ALLOCATE(datPotVelQ%TotVel(NPFlow,3,Nbeta     ,Qfreq%NwQ))
        ALLOCATE(datPotVelQ%RadPot(NPFlow,  Nradiation,Nw))
        ALLOCATE(datPotVelQ%RadVel(NPFlow,3,Nradiation,Nw))
        IF (ID_DATA==ID_FREESURFACE) THEN
        ALLOCATE(datPotVelQ%IncPot(NPFlow,  Nbeta     ,Qfreq%NwQ))
        ALLOCATE(datPotVelQ%IncVel(NPFlow,3,Nbeta     ,Qfreq%NwQ))
        END IF

        DO Ibeta=1,Nbeta
           IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
               IF(Qfreq%NwQ.NE.Nw) THEN
                    Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                    Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                    datPotVelQ%TotPot(:,   Ibeta,:)=datPotVel%TotPot(:,Ibeta,Iw1:Iw2)
                    datPotVelQ%TotVel(:,:, Ibeta,:)=datPotVel%TotVel(:,:,Ibeta,Iw1:Iw2)
                    IF (ID_DATA==ID_FREESURFACE) THEN
                    datPotVelQ%IncPot(:,   Ibeta,:)=datPotVel%IncPot(:,Ibeta,Iw1:Iw2)
                    datPotVelQ%IncVel(:,:, Ibeta,:)=datPotVel%IncVel(:,:,Ibeta,Iw1:Iw2)
                    ENDIF
                ELSE
                    datPotVelQ%TotPot(:,  Ibeta,:)=datPotVel%TotPot(:,Ibeta,:)
                    datPotVelQ%TotVel(:,:,Ibeta,:)=datPotVel%TotVel(:,:,Ibeta,:)
                    IF (ID_DATA==ID_FREESURFACE) THEN
                    datPotVelQ%IncPot(:,  Ibeta,:)=datPotVel%IncPot(:,Ibeta,:)
                    datPotVelQ%IncVel(:,:,Ibeta,:)=datPotVel%IncVel(:,:,Ibeta,:)
                    ENDIF
                ENDIF
            ELSE
               DO Ipanel=1,NPFlow
                 !interpolating potential for the wQ rad. frequencies
                 datPotVelQ%TotPot(Ipanel,Ibeta,:)=                                   &
                         FUNVECT_INTERP1_COMPLEX(w,datPotVel%TotPot(Ipanel,Ibeta,:),      &
                                      Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
                 DO Ixyz=1,3 !vx,vy,vz
                 datPotVelQ%TotVel(Ipanel,Ixyz, Ibeta,:)=                             &
                         FUNVECT_INTERP1_COMPLEX(w,datPotVel%TotVel(Ipanel,Ixyz,Ibeta,:), &
                                      Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
                 ENDDO
                 IF (ID_DATA==ID_FREESURFACE) THEN
                    datPotVelQ%IncPot(Ipanel,Ibeta,:)=                                   &
                            FUNVECT_INTERP1_COMPLEX(w,datPotVel%IncPot(Ipanel,Ibeta,:),      &
                                         Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
                    DO Ixyz=1,3 !vx,vy,vz
                    datPotVelQ%IncVel(Ipanel,Ixyz, Ibeta,:)=                             &
                            FUNVECT_INTERP1_COMPLEX(w,datPotVel%IncVel(Ipanel,Ixyz,Ibeta,:), &
                                         Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
                    ENDDO
                 ENDIF
               ENDDO
           ENDIF
        ENDDO
        !keep radiation potential as the calculated potential in NEMOH first order
        !will be taken/interpolated for difference and sum frequency later
        datPotVelQ%RadPot=datPotVel%RadPot
        datPotVelQ%RadVel=datPotVel%RadVel

        !destroy the initial data
        DEALLOCATE(datPotVel%TotPot,datPotVel%TotVel)
        DEALLOCATE(datPotVel%RadPot,datPotVel%RadVel)

        IF (ID_DATA==ID_FREESURFACE) THEN
        DEALLOCATE(datPotVel%IncPot,datPotVel%IncVel)
        ENDIF
END SUBROUTINE

SUBROUTINE PREPARE_BODY_DISPLACEMENT(Qfreq,Nw,w,Nbeta,Nradiation,NPFlow,Nbodies,      &
                                       Mesh,WLine,InpNEMOHCAL,Motion,BdisplaceQ)
        TYPE(TQfreq),                          INTENT(IN)::Qfreq
        INTEGER,                               INTENT(IN)::Nw,Nbeta,Nradiation
        INTEGER,                               INTENT(IN)::Nbodies,NPFlow
        REAL,DIMENSION(Nw),                    INTENT(IN)::w
        TYPE(TMesh),                           INTENT(IN)::Mesh
        TYPE(TWLine),                          INTENT(IN)::WLine
        TYPE(TNemCal),                         INTENT(IN)::InpNEMOHCAL

        COMPLEX,DIMENSION(Nw,Nradiation,Nbeta),INTENT(IN)::Motion       !RAO
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT):: BdisplaceQ
        !Local
        INTEGER                     :: Ibeta,IdB,Isegline,Iw,Ipanel,Ixyz
        INTEGER                     :: Iw1,Iw2,IwQ,NPanWL
        REAL,DIMENSION(3)           :: vect_R,XCOG,vect_Rsym
        COMPLEX,DIMENSION(3)        :: vect_Theta     !complex rotation RAO
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:)  :: Bdisplace
        ALLOCATE(BdisplaceQ(Qfreq%NwQ,Nbeta,NPFlow,3))
        ALLOCATE(Bdisplace(Nw,Nbeta,NPFlow,3))
        NpanWL=Mesh%Npanels+WLine%NWLineseg

        DO Ipanel=1,Mesh%Npanels+WLine%NWLineseg
                IF (Ipanel .LE. Mesh%Npanels) THEN
                  XCOG=InpNEMOHCAL%bodyinput(Mesh%cPanel(Ipanel))%IntCase(4)%Axis(1:3)
                  IdB=6*(Mesh%cPanel(Ipanel)-1)
                  vect_R=Mesh%XM(:,Ipanel)-XCOG
                ELSE
                  Isegline=Ipanel-Mesh%Npanels
                  XCOG=InpNEMOHCAL%bodyinput(Mesh%cPanel(Wline%IndexPanel(Isegline) ))     &
                                                                 %IntCase(4)%Axis(1:3)
                  IdB=6*(Mesh%cPanel(Wline%IndexPanel(Isegline))-1)
                  vect_R=Wline%XM(Isegline,1:3)-XCOG
                ENDIF

                IF (Mesh%iSym.EQ.1) THEN
                    vect_Rsym=vect_R
                    IF (Ipanel .LE. Mesh%Npanels) THEN
                        vect_Rsym(2)=-Mesh%XM(2,Ipanel)-XCOG(2)
                    ELSE
                        Isegline=Ipanel-Mesh%Npanels
                        vect_Rsym(2)=-Wline%XM(Isegline,2)-XCOG(2)
                    ENDIF
                ENDIF

                DO Ibeta=1,Nbeta
                   DO Iw=1,Nw
                       vect_Theta=Motion(Iw,IdB+4:IdB+6,Ibeta)
                       Bdisplace(Iw,Ibeta,Ipanel,1:3)=Motion(Iw,IdB+1:IdB+3,Ibeta)              &
                            + CROSS_PRODUCT_COMPLEX(vect_Theta,CMPLX(vect_R,0))
                        IF (Mesh%iSym.EQ.1) THEN
                           Bdisplace(Iw,Ibeta,NpanWL+Ipanel,1:3)=                               &
                                   Motion(Iw,IdB+1:IdB+3,Ibeta)                                 &
                                   + CROSS_PRODUCT_COMPLEX(vect_Theta,CMPLX(vect_Rsym,0))
                        ENDIF
                   ENDDO

                    !Interpolation
                   IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
                      IF(Qfreq%NwQ.NE.Nw) THEN
                           Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                           Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                           BdisplaceQ(1:Qfreq%NwQ,Ibeta,Ipanel,1:3)=Bdisplace(Iw1:Iw2,Ibeta,Ipanel,1:3)
                           IF (Mesh%iSym.EQ.1) THEN
                           BdisplaceQ(:,Ibeta,NpanWL+Ipanel,:)=Bdisplace(Iw1:Iw2,Ibeta,NpanWL+Ipanel,:)
                           ENDIF
                       ELSE
                           BdisplaceQ(:,Ibeta,Ipanel,:)=Bdisplace(:,Ibeta,Ipanel,:)
                           IF (Mesh%iSym.EQ.1) THEN
                           BdisplaceQ(:,Ibeta,NpanWL+Ipanel,:)=Bdisplace(:,Ibeta,NpanWL+Ipanel,:)
                           ENDIF
                      ENDIF

                    ELSE
                        DO Ixyz=1,3 !for x,y,z
                            BdisplaceQ(:,Ibeta,Ipanel,Ixyz)=                                               &
                                    FUNVECT_INTERP1_COMPLEX(w,Bdisplace(:,Ibeta,Ipanel,Ixyz),              &
                                         Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
                            IF (Mesh%iSym.EQ.1) THEN
                            BdisplaceQ(:,Ibeta,NpanWL+Ipanel,Ixyz)=                                        &
                                    FUNVECT_INTERP1_COMPLEX(w,Bdisplace(:,Ibeta,NpanWL+Ipanel,Ixyz),       &
                                         Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)

                            ENDIF
                         ENDDO

                    ENDIF
                ENDDO
        ENDDO
        DEALLOCATE(Bdisplace)
END SUBROUTINE

SUBROUTINE PREPARE_INERTIA_FORCES(MechCoef,Motion,Nw,Nbeta,Nradiation,Nintegration,&
                w,Qfreq,Forces1,InertiaForceQ)
        !FInertia=[FHydrodynamic1]+[Fhydrostatic1]
        !       =[Fexc1+Ma*ddXidt2+B*dXidt]+[-K*Xi]
        INTEGER,                        INTENT(IN):: Nradiation,Nintegration
        INTEGER,                        INTENT(IN):: Nw,Nbeta
        TYPE(TLoad1),                   INTENT(IN):: Forces1
        TYPE(TMech),                    INTENT(IN):: MechCoef
        REAL,DIMENSION(Nw),             INTENT(IN):: w
        TYPE(TQfreq),                   INTENT(IN):: Qfreq
        COMPLEX,DIMENSION(Nw,Nradiation,Nbeta),INTENT(IN)::Motion       !RAO
        COMPLEX,DIMENSION(Qfreq%NwQ,Nbeta,Nintegration),                          &
                                        INTENT(OUT)::InertiaForceQ
        !local
        COMPLEX,DIMENSION(Nintegration)          ::HydroDynForce,HydroStatForce
        COMPLEX,DIMENSION(Nw,Nintegration)       ::InertiaForce
        REAL,DIMENSION(Nradiation,Nradiation)    ::StiffMat
        REAL,DIMENSION(Nradiation,Nradiation)    ::AddedMass,DampCoef
        COMPLEX,DIMENSION(Nradiation,Nradiation)    ::MATHyd
        COMPLEX,DIMENSION(Nintegration)          ::ExcitForce
        COMPLEX,DIMENSION(Nintegration)          ::Xi,ddXidt2,dXidt
        INTEGER                                  ::Iw,Ibeta,Iinteg
        INTEGER                                  ::Iw1,Iw2,IwQ,J

        StiffMat=MechCoef%StiffMat+MechCoef%StiffMat_EXT
        DampCoef=MechCoef%DampCoefMat_EXT

        DO Ibeta=1,Nbeta
           DO Iw=1,Nw
           ExcitForce   = Forces1%excitation(Iw,:,Ibeta)
           DampCoef     = MechCoef%DampCoefMat_EXT+Forces1%dampcoef(Iw,:,:)
           AddedMass    = Forces1%addedmass(Iw,:,:)
           Xi           = Motion(Iw,:,Ibeta)
           ddXidt2      = -w(Iw)*w(Iw)*Xi
           dXidt        = -II*w(Iw)*Xi
           !MATHyd       =-w(Iw)*w(Iw)*AddedMass-II*w(Iw)*DampCoef+StiffMat

           HydrodynForce= ExcitForce                            &
                          +MATMUL(-AddedMass,ddXidt2)           &
                          +MATMUL(-DampCoef,dXidt)
           HydroStatForce=-MATMUL(StiffMat,Xi)
           InertiaForce(Iw,:)=HydroDynForce+HydroStatForce
           !InertiaForce(Iw,:)=ExcitForce-MATMUL(MATHyd,Xi)
               ! DO Iinteg=1,Nintegration
                !print*,Iw,Xi(Iinteg)
                !print*,Iw,(REAL(MATHyd(Iinteg,J)),J=1,6)
                !WRITE(*,'(I3,6(X,F12.8))'),Iw,(AIMAG(MATHyd(Iinteg,J)),J=1,6)
                !print*,Iw,ExcitForce(Iinteg)
                !print*,Iw,InertiaForce(Iw,Iinteg)
               ! ENDDO
           ENDDO
            !Interpolation
            IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
               IF(Qfreq%NwQ.NE.Nw) THEN
                    Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                    Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                    InertiaForceQ(:,Ibeta,:)=InertiaForce(Iw1:Iw2,:)
                ELSE
                    InertiaForceQ(:,Ibeta,:)=InertiaForce(:,:)
                ENDIF
             ELSE
             !interpolating displacement for the wQ rad. frequencies
             DO Iinteg=1,Nintegration
                InertiaForceQ(:,Ibeta,Iinteg)=                                            &
                                        FUNVECT_INTERP1_COMPLEX(w,InertiaForce(:,Iinteg), &
                                                        Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
             ENDDO
           ENDIF
        ENDDO
END SUBROUTINE

SUBROUTINE PREPARE_ROTATION_ANGLES(Motion,Nw,Nbeta,Nradiation,&
                Nbodies,w,Qfreq,RotAnglesQ)
        INTEGER,                        INTENT(IN):: Nradiation
        INTEGER,                        INTENT(IN):: Nw,Nbeta,Nbodies
        REAL,DIMENSION(Nw),             INTENT(IN):: w
        TYPE(TQfreq),                   INTENT(IN):: Qfreq
        COMPLEX,DIMENSION(Nw,Nradiation,Nbeta),INTENT(IN)::Motion       !RAO
        COMPLEX,DIMENSION(Qfreq%NwQ,Nbeta,3*Nbodies),                    &
                                        INTENT(OUT)::RotAnglesQ
        !local
        INTEGER                                  ::Iw,Ibeta,Iinteg
        INTEGER                                  ::Itheta06,Itheta03,Itheta,Ibody
        INTEGER                                  ::Iw1,Iw2,IwQ

        DO Ibody=1,Nbodies
           DO Ibeta=1,Nbeta
              Itheta06=6*(Ibody-1)
              Itheta03=3*(Ibody-1)
              !Interpolation
              IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
                IF(Qfreq%NwQ.NE.Nw) THEN
                 Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                 Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                 RotAnglesQ(:,Ibeta,Itheta03+1:Itheta03+3)=Motion(Iw1:Iw2,Itheta06+4:Itheta06+6,Ibeta)
                ELSE
                 RotAnglesQ(:,Ibeta,Itheta03+1:Itheta03+3)=Motion(:,Itheta06+4:Itheta06+6,Ibeta)
                ENDIF
               ELSE
               !interpolating displacement for the wQ rad. frequencies
                DO Itheta=1,3 !theta_x,theta_y,theta_z
                  RotAnglesQ(:,Ibeta,Itheta03+Itheta)=                                    &
                             FUNVECT_INTERP1_COMPLEX(w,Motion(:,Itheta06+3+Itheta,Ibeta), &
                                                        Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)

                ENDDO
              ENDIF
           ENDDO
        ENDDO
END SUBROUTINE

SUBROUTINE PREPARE_TRANSLATION_MOTION(Motion,Nw,Nbeta,Nradiation,&
                Nbodies,w,Qfreq,TRANSMOTQ)
        INTEGER,                        INTENT(IN):: Nradiation
        INTEGER,                        INTENT(IN):: Nw,Nbeta,Nbodies
        REAL,DIMENSION(Nw),             INTENT(IN):: w
        TYPE(TQfreq),                   INTENT(IN):: Qfreq
        COMPLEX,DIMENSION(Nw,Nradiation,Nbeta),INTENT(IN)::Motion       !RAO
        COMPLEX,DIMENSION(Qfreq%NwQ,Nbeta,3*Nbodies),                    &
                                        INTENT(OUT)::TRANSMOTQ
        !local
        INTEGER                                  ::Iw,Ibeta,Iinteg
        INTEGER                                  ::Ind06,Ind03,Ind,Ibody
        INTEGER                                  ::Iw1,Iw2,IwQ

        DO Ibody=1,Nbodies
           DO Ibeta=1,Nbeta
               Ind06=6*(Ibody-1)
               Ind03=3*(Ibody-1)
              !Interpolation
              IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
                IF(Qfreq%NwQ.NE.Nw) THEN
                 Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                 Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                 TRANSMOTQ(:,Ibeta,Ind03+1:Ind03+3)=Motion(Iw1:Iw2,Ind06+1:Ind06+3,Ibeta)
                ELSE
                 TRANSMOTQ(:,Ibeta,Ind03+1:Ind03+3)=Motion(:,Ind06+1:Ind06+3,Ibeta)
                ENDIF
               ELSE
               !interpolating displacement for the wQ rad. frequencies
                DO Ind=1,3 !theta_x,theta_y,theta_z
                TRANSMOTQ(:,Ibeta,Ind03+Ind)=                                     &
                             FUNVECT_INTERP1_COMPLEX(w,Motion(:,Ind06+Ind,Ibeta), &
                                                  Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)

                ENDDO
              ENDIF
           ENDDO
        ENDDO
END SUBROUTINE


SUBROUTINE Discretized_omega_wavenumber_for_QTF               &
                        (Nw,w,kw,NwQ,wQminmax, Nbeta,beta,BdySpeed,    &
                          depth,g,QFreq)
           INTEGER,                         INTENT(IN) ::Nw,NwQ,Nbeta
           REAL,        DIMENSION(2),       INTENT(IN) ::wQminmax
           REAL,        DIMENSION(Nbeta),   INTENT(IN) ::beta
           REAL,                            INTENT(IN) ::BdySpeed,     &
                                                         depth,g
           REAL,        DIMENSION(Nw),      INTENT(IN) ::w,kw
           TYPE(TQfreq),                    INTENT(INOUT)::Qfreq

           INTEGER                                    ::Ibeta,Iw
           REAL                                       ::wQmax,wQmin,dwQ
           REAL,DIMENSION(NwQ)                        ::wQ_temp,kQ_temp
           INTEGER                                    ::Ind1,Ind2
           INTEGER                                    ::InterpPotSwitch

           !In case NO forward speed, wQ is not depend on beta
           !Allocation
           ALLOCATE(Qfreq%wQ(NwQ,Nbeta),Qfreq%kQ(NwQ,Nbeta))
           ALLOCATE(Qfreq%diffwQ(NwQ,Nbeta),Qfreq%sumwQ(NwQ-1,Nbeta))
           ALLOCATE(Qfreq%InterpPotSwitch(Nbeta))

           wQmin=wQminmax(1)
           wQmax=wQminmax(2)

           IF (NwQ.GT.1) THEN
                dwQ=(wQmax-wQmin)/(NwQ-1)
                IF (ABS(w(2)-w(1)-dwQ).LT.0.01) dwQ=w(2)-w(1)
                DO Iw=1,NwQ
                    wQ_temp(Iw)=wQmin+dwQ*(Iw-1)
                END DO
            ELSE
             WRITE(*,*) &
             "ERROR: QTF can not be run with 1 frequency only!"
             STOP
            END IF
           IF (w(2)-w(1)==dwQ) THEN
                Ind1=Fun_closest(Nw,w,wQmin)
                Ind2=Fun_closest(Nw,w,wQmax)
                kQ_temp=kw(Ind1:Ind2)
                InterpPotSwitch=0
           ELSE
                kQ_temp=FunVect_inverseDispersion(NwQ,wQ_temp,depth,g)!wavenumbers
                InterpPotSwitch=1
           ENDIF

           DO Ibeta=1,Nbeta
             Qfreq%wQ(:,Ibeta)=wQ_temp-kQ_temp*BdySpeed*COS(beta(Ibeta))
             IF (BdySpeed>0.) THEN
                Qfreq%kQ(:,Ibeta)=FunVect_inverseDispersion(                 &
                                NwQ,Qfreq%wQ(:,Ibeta),depth,g)!wavenumbers
                Qfreq%InterpPotSwitch(Ibeta)=1
             ELSE
                Qfreq%kQ(:,Ibeta)=kQ_temp
                Qfreq%InterpPotSwitch(Ibeta)=0
                IF (InterpPotSwitch==1) Qfreq%InterpPotSwitch(Ibeta)=1
             ENDIF
             Qfreq%diffwQ(1,Ibeta)=0
             Qfreq%diffwQ(2:NwQ,Ibeta)=Qfreq%wQ(2:NwQ,Ibeta)-Qfreq%wQ(1,Ibeta)
             Qfreq%sumwQ(1:NwQ-1,Ibeta)=Qfreq%wQ(1:NwQ-1,Ibeta)+Qfreq%wQ(1,Ibeta)
             Qfreq%NwQ=NwQ
            ! IF ((Fun_MIN(NwQ-1,Qfreq%diffwQ(2:NwQ,Ibeta))-w(1))<-0.0001) THEN
            !        print*,'INPUT ERROR: Min. diff. rad freq < w(1) in QTFpreproc data!'
            !        STOP
            ! ENDIF
            ! IF ((Fun_MAX(NwQ-1,Qfreq%sumwQ(1:NwQ-1,Ibeta))-w(Nw))>0.0001) THEN
            !        print*,'INPUT ERROR: Max. sum rad freq > w(Nw) in QTFpreproc data!'
            !        STOP
            ! ENDIF
           ENDDO
END SUBROUTINE

SUBROUTINE CALC_GENERALIZED_NORMAL_WATERLINE_dSEGMENT(Mesh,Nintegration,WLine, &
                                                InpNEMOHCAL,genNormalWLine_dGamma)
        USE MNemohCal,          ONLY:TNemCal
        IMPLICIT NONE
        TYPE(TMesh),                                  INTENT(IN):: Mesh
        TYPE(TWLINE),                                 INTENT(IN):: WLine
        TYPE(TNemCal),                                INTENT(IN):: InpNEMOHCAL
        INTEGER,                                      INTENT(IN):: Nintegration
        REAL, DIMENSION(Nintegration,Wline%NWLineseg*2**Mesh%iSym), &
                                              INTENT(INOUT) ::genNormalWLine_dGamma
        !local
        INTEGER IdBody, IDMode,c,indsum
        REAL, DIMENSION(Wline%NWLineseg*2**Mesh%iSym) :: NDGamma

        indsum=1
        DO IdBody=1,InpNEMOHCAL%Nbodies
           DO IdMode=1,InpNEMOHCAL%bodyinput(IdBody)%NIntegration
           CALL ComputeNDGamma_WLINE(Mesh,WLINE,IdBody,                                 &
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%ICase,                    &
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%Direction(1:3),           &
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%Axis(1:3),NDGamma)
                DO c=1,Wline%NWLineseg*2**Mesh%iSym
                genNormalWLine_dGamma(indsum,c)=NDGamma(c)
                END DO
                indsum=indsum+1
            END DO
        END DO
END SUBROUTINE

!-- SUBROUTINE ComputeNDGamma waterline
    SUBROUTINE ComputeNDGamma_WLINE(Mesh,WLINE,c,iCase,Direction,Axis,NDGamma)
    USE MMesh
    USE MFace
    IMPLICIT NONE
    TYPE(TMesh)  :: Mesh
    TYPE(TWline) :: Wline
    INTEGER :: c,iCase
    REAL,DIMENSION(3) :: Direction,Axis
    REAL,DIMENSION(*) :: NDGamma
    REAL,DIMENSION(3) :: VEL
    INTEGER :: i,ipanel
    REAL,DIMENSION(3) :: N,XM,DIR

    Dir=Direction
    SELECT CASE (iCase)
    CASE (1)
        DO i=1,Wline%NWLineseg
            ipanel=Wline%IndexPanel(i)
            N(1:2)=Mesh%N(1:2,ipanel)
            N(3)=0 ! Normal in vertical axis is always zero on waterline

            IF (Mesh%cPanel(ipanel).EQ.c) THEN
                VEL(1)=Dir(1)
                VEL(2)=Dir(2)
                VEL(3)=Dir(3)
                NDGamma(i)=(N(1)*VEL(1)+N(2)*VEL(2)+N(3)*VEL(3))*Wline%SegLength(i)
            ELSE
                NDGamma(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(ipanel).EQ.c) THEN
                    VEL(1)=Dir(1)
                    VEL(2)=Dir(2)
                    VEL(3)=Dir(3)
                    NDGamma(i+WLINE%NWLineseg)=(N(1)*VEL(1)-N(2)*VEL(2)+N(3)*VEL(3))*Wline%SegLength(i)
                 ELSE
                    NDGamma(i+WLINE%NWLineseg)=0.
                 END IF
            END IF
        END DO
    CASE (2)
        DO i=1,Wline%NWLineseg
            ipanel=Wline%IndexPanel(i)
            N(1:2)=Mesh%N(1:2,ipanel)
            N(3)=0 ! Normal in vertical axis is always zero on waterline

            XM(1:2)=Wline%XM(i,1:2)
            XM(3)=0
            IF (Mesh%cPanel(ipanel).EQ.c) THEN
               !Only on horizontal axes, vertical axes 0, then only sway
                VEL(1)=Dir(2)*(XM(3)-Axis(3))-Dir(3)*(XM(2)-Axis(2))
                VEL(2)=Dir(3)*(XM(1)-Axis(1))-Dir(1)*(XM(3)-Axis(3))
                VEL(3)=Dir(1)*(XM(2)-Axis(2))-Dir(2)*(XM(1)-Axis(1))
                NDGamma(i)=(N(1)*VEL(1)+N(2)*VEL(2)+N(3)*VEL(3))*Wline%SegLength(i)
            ELSE
                NDGamma(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(ipanel).EQ.c) THEN
                    !Only on horizontal axes, vertical axes 0, then only sway
                    VEL(1)=Dir(2)*(XM(3)-Axis(3))-Dir(3)*(-XM(2)-Axis(2))
                    VEL(2)=Dir(3)*(XM(1)-Axis(1))-Dir(1)*(XM(3)-Axis(3))
                    VEL(3)=Dir(1)*(-XM(2)-Axis(2))-Dir(2)*(XM(1)-Axis(1))
                    NDGamma(i+WLINE%NWLineseg)=(N(1)*VEL(1)-N(2)*VEL(2)+N(3)*VEL(3))*Wline%SegLength(i)
                ELSE
                    NDGamma(i+Wline%NWLineseg)=0.
                 END IF
            END IF
        END DO
    CASE (3)
        WRITE(*,*) 'Error: force case 3 not implemented yet'
        STOP
    CASE DEFAULT
        WRITE(*,*) 'Error: unknown radiation case'
        STOP
    END SELECT
    END SUBROUTINE

SUBROUTINE PREPARE_ASYMP_PARAM(Rext,NR,NBESSEL,ASYMP)
    INTEGER,    INTENT(IN):: NR,NBESSEL
    REAL,       INTENT(IN):: Rext
    TYPE(TASYMP)          :: ASYMP
    INTEGER               :: Ir
    REAL                  ::dR
    ALLOCATE(ASYMP%Rf(NR))
    dR=Rext/(NR-1)
    ASYMP%NBESSEL=NBESSEL
    ASYMP%NR=NR
    ASYMP%dR=dR
    DO Ir=0,NR-1
    ASYMP%Rf(Ir+1)=dR*Ir
    ENDDO
END SUBROUTINE

SUBROUTINE PREPARE_SOURCE_DISTRIBUTION(wd,Qfreq,Nw,w,Nbeta,beta,Npanels, &
                                        Nradiation,Motion,SourceDistrQ)
    !input/output
    CHARACTER(LEN=*),                      INTENT(IN)   ::wd
    TYPE(TQfreq),                          INTENT(IN)   ::Qfreq
    INTEGER,                               INTENT(IN)   ::Nw,Nbeta
    INTEGER,                               INTENT(IN)   ::NPanels,Nradiation
    REAL,   DIMENSION(Nw),                 INTENT(IN)   ::w
    REAL,   DIMENSION(Nbeta),              INTENT(IN)   ::beta
    COMPLEX,DIMENSION(Nw,Nradiation,Nbeta),INTENT(IN)   :: Motion
    TYPE(TSourceQ),                        INTENT(INOUT)::SourceDistrQ

    INTEGER                                      ::Iw,Iw1,Iw2,Ibeta,Irad,Ipanel
    TYPE(TSourceQ)                               ::SourceDistr
    TYPE(TSource)                                ::SourceDistr_Iw

    ALLOCATE(SourceDistr%ZIGB_Per(Npanels,Nbeta,Nw))
    ALLOCATE(SourceDistr%ZIGS_Per(Npanels,Nbeta,Nw))
    ALLOCATE(SourceDistr%ZIGB_Rad(Npanels,Nradiation,Nw))
    ALLOCATE(SourceDistr%ZIGS_Rad(Npanels,Nradiation,Nw))

    ALLOCATE(SourceDistr_Iw%ZIGB(Npanels,Nradiation+Nbeta))
    ALLOCATE(SourceDistr_Iw%ZIGS(Npanels,Nradiation+Nbeta))

    DO Iw=1,Nw
      CALL Read_SourceDistribution(TRIM(wd),Iw,Nw,Nradiation,Nbeta,                 &
                                   Npanels,SourceDistr_Iw)
        DO Ibeta=1,Nbeta
           !CONSTRUCT PERTURBATION(Diff+Rad) SINGULAR DISTRIBUTION FOR EACH WAVE DIRECTION
           !assign the diffraction singular source distribution for all panels
           SourceDistr%ZIGB_Per(1:Npanels,Ibeta,Iw)=                                  &
                                      SourceDistr_Iw%ZIGB(1:Npanels,Nradiation+Ibeta)
           SourceDistr%ZIGS_Per(1:Npanels,Ibeta,Iw)=                                  &
                                      SourceDistr_Iw%ZIGS(1:Npanels,Nradiation+Ibeta)
             !sum the diffraction + radiation singular source distribution for all panels
             DO Irad=1,Nradiation
              SourceDistr%ZIGB_Per(1:Npanels,Ibeta,Iw)                                &
                      =SourceDistr%ZIGB_Per(1:Npanels,Ibeta,Iw)                       &
                                 -II*w(Iw)*SourceDistr_Iw%ZIGB(1:Npanels,Irad)        &
                                  *Motion(Iw,Irad,Ibeta)
              SourceDistr%ZIGS_Per(1:Npanels,Ibeta,Iw)                                &
                      =SourceDistr%ZIGS_Per(1:Npanels,Ibeta,Iw)                       &
                                 -II*w(Iw)*SourceDistr_Iw%ZIGS(1:Npanels,Irad)        &
                                  *Motion(Iw,Irad,Ibeta)
              ENDDO
      ENDDO

      DO Irad=1,Nradiation
         SourceDistr%ZIGB_Rad(1:Npanels,Irad,Iw)= SourceDistr_Iw%ZIGB(1:Npanels,Irad)
         SourceDistr%ZIGS_Rad(1:Npanels,Irad,Iw)= SourceDistr_Iw%ZIGS(1:Npanels,Irad)
      ENDDO
    ENDDO
    DEALLOCATE(SourceDistr_Iw%ZIGB,SourceDistr_Iw%ZIGS)

    !adjust frequency as user input interval
    !interpolation applied or not
    ALLOCATE(SourceDistrQ%ZIGB_Per(Npanels,Nbeta,Qfreq%NwQ))
    ALLOCATE(SourceDistrQ%ZIGS_Per(Npanels,Nbeta,Qfreq%NwQ))
    ALLOCATE(SourceDistrQ%ZIGB_Rad(Npanels,Nradiation,Nw))
    ALLOCATE(SourceDistrQ%ZIGS_Rad(Npanels,Nradiation,Nw))
    DO Ibeta=1,Nbeta
      IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
          IF(Qfreq%NwQ.NE.Nw) THEN
               Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
               Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
               SourceDistrQ%ZIGB_Per(:,Ibeta,:)=SourceDistr%ZIGB_Per(:,Ibeta,Iw1:Iw2)
               SourceDistrQ%ZIGS_Per(:,Ibeta,:)=SourceDistr%ZIGS_Per(:,Ibeta,Iw1:Iw2)
          ELSE
               SourceDistrQ%ZIGB_Per(:,Ibeta,:)=SourceDistr%ZIGB_Per(:,Ibeta,:)
               SourceDistrQ%ZIGS_Per(:,Ibeta,:)=SourceDistr%ZIGS_Per(:,Ibeta,:)
          ENDIF
       ELSE
          DO Ipanel=1,NPanels
            !interpolating source dirstribution for the wQ rad. frequencies
            SourceDistrQ%ZIGB_Per(Ipanel,Ibeta,:)=                                   &
                    FUNVECT_INTERP1_COMPLEX(w,SourceDistr%ZIGB_Per(Ipanel,Ibeta,:),      &
                                 Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
            SourceDistrQ%ZIGS_Per(Ipanel,Ibeta,:)=                                   &
                    FUNVECT_INTERP1_COMPLEX(w,SourceDistr%ZIGS_Per(Ipanel,Ibeta,:),      &
                                 Nw,Qfreq%wQ(:,Ibeta),Qfreq%NwQ)
          ENDDO
      ENDIF
    ENDDO
    !keep radiation potential as the calculated potential in NEMOH first order
    !will be taken/interpolated for difference and sum frequency later
    SourceDistrQ%ZIGB_Rad=SourceDistr%ZIGB_Rad
    SourceDistrQ%ZIGB_Rad=SourceDistr%ZIGB_Rad

    !destroy the initial data
    DEALLOCATE(SourceDistr%ZIGB_Per)
    DEALLOCATE(SourceDistr%ZIGS_Per)
    DEALLOCATE(SourceDistr%ZIGB_Rad)
    DEALLOCATE(SourceDistr%ZIGS_Rad)
END SUBROUTINE

SUBROUTINE WRITE_QTFSOLVERLOGFILE(wd,Nbeta,beta,Qfreq)
       CHARACTER(LEN=*),             INTENT(IN)::wd
       TYPE(TQfreq),                 INTENT(IN)::Qfreq
       INTEGER,                      INTENT(IN)::Nbeta
       REAL,DIMENSION(Nbeta),        INTENT(IN)::beta
       CHARACTER(LEN=1000)                     ::LogTextToBeWritten
       INTEGER Ibeta

       WRITE(*, *) ' '
       WRITE(LogTextToBeWritten,*) '----Solver (QTF Module)---'
       CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
       DO Ibeta=1,Nbeta
       WRITE(LogTextToBeWritten,'(A,F8.3,A,I4,A,3(F8.3,A))') ' Beta= ',beta(Ibeta)*180/PI,      &
               ' deg, NFreq= ',Qfreq%NwQ, ', omega = (', Qfreq%wQ(1,Ibeta),':'                  &
               ,Qfreq%wQ(2,Ibeta)-Qfreq%wQ(1,Ibeta),':',Qfreq%wQ(Qfreq%NwQ,Ibeta),') rad/s'
       CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
       ENDDO
       WRITE(LogTextToBeWritten,*) '--------------------------'
       CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
END SUBROUTINE

FUNCTION MATMUL_COMPLEX(Mat,Vect,Nv) RESULT(Vres)
        INTEGER,                 INTENT(IN):: Nv
        COMPLEX,DIMENSION(Nv,Nv),INTENT(IN):: Mat
        COMPLEX,DIMENSION(Nv),   INTENT(IN):: Vect
        COMPLEX,DIMENSION(Nv)              :: Vres
        INTEGER                            :: I,J
        DO I=1,Nv
              Vres(I)=CMPLX(0.,0.)
              DO J=1,Nv
              Vres(I)=Vres(I)+Mat(I,J)*Vect(J)
              ENDDO
        ENDDO
END FUNCTION

END MODULE
