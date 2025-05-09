!-----------------------------------------------------------------------------------
!
!    Copyright (C) 2022 - 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
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
!   - Fabien Robaux (EDF/INNOSEA)
!   - G.DELHOMMEAU, THESE DE CHEN XIAO-BO(1988)
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
Module MQSolver
USE MFace,               ONLY:TVFace,TWLine
USE MMesh
USE MQSolverPreparation, ONLY:TPotVel,TQfreq,TASYMP,TSourceQ
USE CONSTANTS          , ONLY:II,CZERO,PI,INFINITE_DEPTH
USE Elementary_functions,ONLY:CROSS_PRODUCT_COMPLEX,Fun_closest,CIH,&
                              DOT_PRODUCT_COMPLEX,COMPLEX_CONJUGATE_VECT
USE MEnvironment,        ONLY:TEnvironment,Fun_inverseDispersion
USE MReadInputFiles,     ONLY:TMeshFS
USE MCallInterp,         ONLY:FUN_INTERP1_COMPLEX
USE MQSOLVERASYMP

IMPLICIT NONE
CONTAINS
  SUBROUTINE COMPUTATION_QTF_QUADRATIC(Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,       &
                                       Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,&
                                       rho,g,datPotVelQ,BdisplaceQ,genNormal_dS, &
                                       genNormalWLine_dGamma,wQ,beta,            &
                                       InertiaForceQ,RotAnglesQ,IntegAxis,       &
                                       StiffMat,TransMotionQ,SwitchQuadHM,       &
                                       QTFDuok)
          !Input/output
          INTEGER,                              INTENT(IN) :: Iw1,Iw2,Ibeta1,Ibeta2
          INTEGER,                              INTENT(IN) :: Nintegration,NPFlow
          INTEGER,                              INTENT(IN) :: NwQ,Nbeta,Nbodies
          REAL,                                 INTENT(IN) :: rho,g
          TYPE(TMesh),                          INTENT(IN) :: Mesh
          TYPE(TVFace),                         INTENT(IN) :: VFace
          TYPE(TWLine),                         INTENT(IN) :: WLine
          TYPE(TPotVel),                        INTENT(IN) :: datPotVelQ
          COMPLEX,DIMENSION(NwQ,Nbeta,NPFlow,3),INTENT(IN) :: BdisplaceQ
          REAL,DIMENSION(Nintegration,Mesh%Npanels*2**Mesh%Isym),                 &
                                                INTENT(IN) :: genNormal_dS
          REAL,DIMENSION(Nintegration,WLine%NWLineSeg*2**Mesh%Isym),              &
                                                INTENT(IN) :: genNormalWLine_dGamma
          REAL,DIMENSION(NwQ,Nbeta),            INTENT(IN) :: wQ
          REAL,DIMENSION(Nbeta),                INTENT(IN) :: beta
          COMPLEX,DIMENSION(NwQ,Nbeta,Nintegration),                              &
                                                INTENT(IN) :: InertiaForceQ
          COMPLEX,DIMENSION(NwQ,Nbeta,3*Nbodies),INTENT(IN):: RotAnglesQ
          COMPLEX,DIMENSION(NwQ,Nbeta,3*Nbodies),INTENT(IN):: TransMotionQ
          REAL,DIMENSION(3,Nintegration),       INTENT(IN) :: IntegAxis
          REAL,DIMENSION(Nintegration,Nintegration),                              &
                                                INTENT(IN) ::StiffMat
          INTEGER,                              INTENT(IN) :: SwitchQuadHM
          COMPLEX,DIMENSION(Nintegration,2,7),  INTENT(OUT):: QTFDuok
          !Local
          COMPLEX                             ::TotPot_Iw1,TotPot_Iw2
          COMPLEX,DIMENSION(3)                ::TotVel_Iw1,TotVel_Iw2
          COMPLEX,DIMENSION(3)                ::Bdisplace_Iw1,Bdisplace_Iw2
          COMPLEX,DIMENSION(Nintegration)     ::InertiaForce_Iw1,InertiaForce_Iw2
          COMPLEX,DIMENSION(Nintegration)     ::XiTilde_M,XiTilde_P
          COMPLEX,DIMENSION(3*Nbodies)        ::RotAngles_Iw1,RotAngles_Iw2
          COMPLEX,DIMENSION(3*Nbodies)        ::TransMotion_Iw1,TransMotion_Iw2
          COMPLEX,DIMENSION(1)                ::eta_Iw1,eta_Iw2,zeta_Iw1,zeta_Iw2
          INTEGER                             ::Npanels,Ipanel,Iinteg,NpanWLin
          COMPLEX                             ::quad_M,quad_P
          REAL                                ::w1,w2
          REAL                                ::r3
          INTEGER                             ::Iwline,Ibody,Itheta0
          INTEGER                             ::Isym,Iterm
          INTEGER,DIMENSION(2)                ::Ipanelinit,Iwlineinit

          w1=wQ(Iw1,Ibeta1)
          w2=wQ(Iw2,Ibeta2)


          InertiaForce_Iw1=InertiaForceQ(Iw1,Ibeta1,:)
          InertiaForce_Iw2=InertiaForceQ(Iw2,Ibeta2,:)

          RotAngles_Iw1=RotAnglesQ(Iw1,Ibeta1,:)
          RotAngles_Iw2=RotAnglesQ(Iw2,Ibeta2,:)

          TransMotion_Iw1=TransMotionQ(Iw1,Ibeta1,:)
          TransMotion_Iw2=TransMotionQ(Iw2,Ibeta2,:)

          Npanels=Mesh%Npanels
          NpanWlin=Npanels+WLine%NWLineSeg


          DO Iinteg=1,Nintegration
           QTFDuok(Iinteg,:,:)=0
           !integration over body panels
           DO Ipanel=1,Npanels
            IF (Mesh%XM(3, Ipanel)<0.) THEN      !dont compute at lid panels
              DO Isym=1,1+Mesh%Isym
                Ipanelinit=Function_SymIpanel(NpanWlin,Npanels,Isym)
                !prepare variables Iw1,Iw2
                TotPot_Iw1=datPotVelQ%TotPot(Ipanelinit(1)+Ipanel,Ibeta1,Iw1)
                TotVel_Iw1=datPotVelQ%TotVel(Ipanelinit(1)+Ipanel,1:3,Ibeta1,Iw1)

                TotPot_Iw2=datPotVelQ%TotPot(Ipanelinit(1)+Ipanel,Ibeta2,Iw2)
                TotVel_Iw2=datPotVelQ%TotVel(Ipanelinit(1)+Ipanel,1:3,Ibeta2,Iw2)

                Bdisplace_Iw1=BdisplaceQ(Iw1,Ibeta1,Ipanelinit(1)+Ipanel,1:3)
                Bdisplace_Iw2=BdisplaceQ(Iw2,Ibeta2,Ipanelinit(1)+Ipanel,1:3)

                ! term: product(grad_phi_1,grad_phi2)
                quad_M=DOT_PRODUCT_DIFF_BIHARM(                                 &
                                   TotVel_Iw1,TotVel_Iw2,TotVel_Iw1,TotVel_Iw2,3)
                quad_P=DOT_PRODUCT_SUM_BIHARM(                                  &
                                   TotVel_Iw1,TotVel_Iw2,TotVel_Iw1,TotVel_Iw2,3)

                QTFDuok(Iinteg,1,1)=QTFDuok(Iinteg,1,1)+0.5*rho*                &
                                quad_M*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
                QTFDuok(Iinteg,2,1)=QTFDuok(Iinteg,2,1)+0.5*rho*                &
                                quad_P*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)


                ! term: product(displacement,dtgradPhi)
                quad_M=DOT_PRODUCT_DIFF_BIHARM(                                 &
                        Bdisplace_Iw1,Bdisplace_Iw2,                            &
                        -II*w1*TotVel_Iw1,-II*w2*TotVel_Iw2,3)
                quad_P=DOT_PRODUCT_SUM_BIHARM(                                  &
                        Bdisplace_Iw1,Bdisplace_Iw2,                            &
                        -II*w1*TotVel_Iw1,-II*w2*TotVel_Iw2,3)
                QTFDuok(Iinteg,1,2)=QTFDuok(Iinteg,1,2)+rho*                    &
                                quad_M*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
                QTFDuok(Iinteg,2,2)=QTFDuok(Iinteg,2,2)+rho*                    &
                                quad_P*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)

              ENDDO
            ENDIF
           ENDDO
           !integration over waterline segments
           DO Iwline=1,WLine%NWLineSeg
               DO Isym=1,1+Mesh%Isym
                Iwlineinit=Function_SymIwline(NpanWlin,Npanels,WLine%NWLineSeg,Isym)
                !preparation var in Iw1 and Iw2
                Bdisplace_Iw1=BdisplaceQ(Iw1,Ibeta1,Iwlineinit(1)+Iwline,1:3)
                Bdisplace_Iw2=BdisplaceQ(Iw2,Ibeta2,Iwlineinit(1)+Iwline,1:3)
                TotPot_Iw1=datPotVelQ%TotPot(Iwlineinit(1)+Iwline,Ibeta1,Iw1)
                TotPot_Iw2=datPotVelQ%TotPot(Iwlineinit(1)+Iwline,Ibeta2,Iw2)
                !term: [eta-zeta]^2
                zeta_Iw1=Bdisplace_Iw1(3)    !vertical Waterline displ.
                zeta_Iw2=Bdisplace_Iw2(3)
                eta_Iw1 =II*w1*TotPot_Iw1/g !wave elevation
                eta_Iw2 =II*w2*TotPot_Iw2/g
                quad_M=DOT_PRODUCT_DIFF_BIHARM(                                    &
                     eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,                            &
                     eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,1)
                quad_P=DOT_PRODUCT_SUM_BIHARM(                                     &
                     eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,                            &
                     eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,1)
                QTFDuok(Iinteg,1,3)=QTFDuok(Iinteg,1,3)-0.5*rho*g*quad_M           &
                                *genNormalWLine_dGamma(Iinteg,Iwlineinit(2)+Iwline)
                QTFDuok(Iinteg,2,3)=QTFDuok(Iinteg,2,3)-0.5*rho*g*quad_P           &
                                *genNormalWLine_dGamma(Iinteg,Iwlineinit(2)+Iwline)

              ENDDO
           ENDDO
         ENDDO

          DO Ibody=1,Nbodies
             Itheta0=3*(Ibody-1)
             Iinteg=6*(Ibody-1)

             !term: Matrix product Rotation matrix and  Inertia Force
             !For diff freq
             !R(F_I) for translation modes (1,2,3) of FI
             QTFDuok(Iinteg+1:Iinteg+3,1,4)=QTFDuok(Iinteg+1:Iinteg+3,1,4)+   &
                      CROSS_PRODUCT_DIFF_BIHARM(                              &
                      RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
                      RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
                      InertiaForce_Iw1(Iinteg+1:Iinteg+3),                    &
                      InertiaForce_Iw2(Iinteg+1:Iinteg+3))
             !R(F_I) for rotation modes (4,5,6) of FI
             QTFDuok(Iinteg+4:Iinteg+6,1,4)=QTFDuok(Iinteg+4:Iinteg+6,1,4)+   &
                      CROSS_PRODUCT_DIFF_BIHARM(                              &
                      RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
                      RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
                      InertiaForce_Iw1(Iinteg+4:Iinteg+6),                    &
                      InertiaForce_Iw2(Iinteg+4:Iinteg+6))
             !For sum freq
             !R(F_I) for translation modes (1,2,3) of FI
             QTFDuok(Iinteg+1:Iinteg+3,2,4)=QTFDuok(Iinteg+1:Iinteg+3,2,4)+   &
                      CROSS_PRODUCT_SUM_BIHARM(                               &
                      RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
                      RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
                      InertiaForce_Iw1(Iinteg+1:Iinteg+3),                    &
                      InertiaForce_Iw2(Iinteg+1:Iinteg+3))
             !R(F_I) for rotation modes (4,5,6) of FI
             QTFDuok(Iinteg+4:Iinteg+6,2,4)=QTFDuok(Iinteg+4:Iinteg+6,2,4)+   &
                      CROSS_PRODUCT_SUM_BIHARM(                               &
                      RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
                      RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
                      InertiaForce_Iw1(Iinteg+4:Iinteg+6),                    &
                      InertiaForce_Iw2(Iinteg+4:Iinteg+6))
             IF (SwitchQuadHM==1) THEN
              !Addition term for Moment on displaced(translation) position
              !For diff freq
              !CrossProduct(X,F_I)
              QTFDuok(Iinteg+4:Iinteg+6,1,5)= QTFDuok(Iinteg+4:Iinteg+6,1,5)+  &
                       CROSS_PRODUCT_DIFF_BIHARM(                              &
                       TransMotion_Iw1(Itheta0+1:Itheta0+3),                   &
                       TransMotion_Iw2(Itheta0+1:Itheta0+3),                   &
                       InertiaForce_Iw1(Iinteg+1:Iinteg+3),                    &
                       InertiaForce_Iw2(Iinteg+1:Iinteg+3))
              !For sum freq
              QTFDuok(Iinteg+4:Iinteg+6,2,5)= QTFDuok(Iinteg+4:Iinteg+6,2,5)+  &
                       CROSS_PRODUCT_SUM_BIHARM(                               &
                       TransMotion_Iw1(Itheta0+1:Itheta0+3),                   &
                       TransMotion_Iw2(Itheta0+1:Itheta0+3),                   &
                       InertiaForce_Iw1(Iinteg+1:Iinteg+3),                    &
                       InertiaForce_Iw2(Iinteg+1:Iinteg+3))


              !term: quadratic term of second order excitation force=-[K]*Xi_Tilde
              r3=-IntegAxis(3,Iinteg+4)!r3=Zm-ZCOG, ZM=0
              XiTilde_M(Iinteg+1:Iinteg+6)=CZERO
              XiTilde_P(Iinteg+1:Iinteg+6)=CZERO
              XiTilde_M(Iinteg+3)=-0.5*r3*(DOT_PRODUCT_DIFF_BIHARM(             &
                      RotAngles_Iw1(Itheta0+1),RotAngles_Iw2(Itheta0+1),        &
                      RotAngles_Iw1(Itheta0+1),RotAngles_Iw2(Itheta0+1),1)      &
                      +DOT_PRODUCT_DIFF_BIHARM(                                 &
                      RotAngles_Iw1(Itheta0+2),RotAngles_Iw2(Itheta0+2),        &
                      RotAngles_Iw1(Itheta0+2),RotAngles_Iw2(Itheta0+2),1))
              XiTilde_M(Iinteg+4)=0.5*DOT_PRODUCT_DIFF_BIHARM(                  &
                      RotAngles_Iw1(Itheta0+2),RotAngles_Iw2(Itheta0+2),        &
                      RotAngles_Iw1(Itheta0+3),RotAngles_Iw2(Itheta0+3),1)
              XiTilde_M(Iinteg+5)=-0.5*DOT_PRODUCT_DIFF_BIHARM(                 &
                      RotAngles_Iw1(Itheta0+3),RotAngles_Iw2(Itheta0+3),        &
                      RotAngles_Iw1(Itheta0+1),RotAngles_Iw2(Itheta0+1),1)
              XiTilde_P(Iinteg+3)=-0.5*r3*(DOT_PRODUCT_SUM_BIHARM(              &
                      RotAngles_Iw1(Itheta0+1),RotAngles_Iw2(Itheta0+1),        &
                      RotAngles_Iw1(Itheta0+1),RotAngles_Iw2(Itheta0+1),1)      &
                      +DOT_PRODUCT_SUM_BIHARM(                                  &
                      RotAngles_Iw1(Itheta0+2),RotAngles_Iw2(Itheta0+2),        &
                      RotAngles_Iw1(Itheta0+2),RotAngles_Iw2(Itheta0+2),1))
              XiTilde_P(Iinteg+4)=0.5*DOT_PRODUCT_SUM_BIHARM(                   &
                      RotAngles_Iw1(Itheta0+2),RotAngles_Iw2(Itheta0+2),        &
                      RotAngles_Iw1(Itheta0+3),RotAngles_Iw2(Itheta0+3),1)
              XiTilde_P(Iinteg+5)=-0.5*DOT_PRODUCT_SUM_BIHARM(                  &
                      RotAngles_Iw1(Itheta0+3),RotAngles_Iw2(Itheta0+3),        &
                      RotAngles_Iw1(Itheta0+1),RotAngles_Iw2(Itheta0+1),1)
            ENDIF
          ENDDO

            IF (SwitchQuadHM==1) THEN
              !term: quadratic term of second order excitation force=-[K]*Xi_Tilde
              QTFDuok(:,1,6)=QTFDuok(:,1,6)-MATMUL(StiffMat,XiTilde_M)
              QTFDuok(:,2,6)=QTFDuok(:,2,6)-MATMUL(StiffMat,XiTilde_P)
            ENDIF
          DO Iterm=1,6
          !TRANSFER FUNCTION = BICHROMATIC FORCES / 2?
          QTFDuok(:,:,Iterm)=QTFDuok(:,:,Iterm)/2   !TRANSFER FUNCTION = BICHROMATIC FORCES / 2?
          QTFDuok(:,:,7)=QTFDuok(:,:,7)+QTFDuok(:,:,Iterm)!Duok total
          ENDDO

  END SUBROUTINE

  SUBROUTINE COMPUTATION_QTF_POTENTIAL_BODYFORCE                              &
                                   (Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,       &
                                    Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,&
                                    env,datPotVelQ,BdisplaceQ,genNormal_dS,   &
                                    Nw,w,Qfreq,beta,RotAnglesQ,               &
                                    QTFHasbo)
          !Input/output
          INTEGER,                              INTENT(IN) :: Iw1,Iw2,Ibeta1,Ibeta2
          INTEGER,                              INTENT(IN) :: Nintegration,NPFlow
          INTEGER,                              INTENT(IN) :: NwQ,Nbeta,Nbodies,Nw
          TYPE(TMesh),                          INTENT(IN) :: Mesh
          TYPE(TVFace),                         INTENT(IN) :: VFace
          TYPE(TWLine),                         INTENT(IN) :: WLine
          TYPE(TPotVel),                        INTENT(IN) :: datPotVelQ
          TYPE(TQfreq),                         INTENT(IN) :: Qfreq
          TYPE(TEnvironment),                   INTENT(IN) :: Env
          COMPLEX,DIMENSION(NwQ,Nbeta,NPFlow,3),INTENT(IN) :: BdisplaceQ
          COMPLEX,DIMENSION(NwQ,Nbeta,3*Nbodies),INTENT(IN):: RotAnglesQ
          REAL,DIMENSION(Nintegration,Mesh%Npanels*2**Mesh%Isym),             &
                                                INTENT(IN) :: genNormal_dS
          REAL,DIMENSION(Nw),                   INTENT(IN) :: w
          REAL,DIMENSION(Nbeta),                INTENT(IN) :: beta
          COMPLEX,DIMENSION(Nintegration,2,7),  INTENT(OUT):: QTFHasbo
          !Local
          COMPLEX                             ::TotPot_Iw1,TotPot_Iw2
          COMPLEX,DIMENSION(3)                ::TotVel_Iw1,TotVel_Iw2
          COMPLEX,DIMENSION(3)                ::Bdisplace_Iw1,Bdisplace_Iw2
          COMPLEX,DIMENSION(3)                ::dtBdisplace_Iw1,dtBdisplace_Iw2
          COMPLEX,DIMENSION(3*Nbodies)        ::RotAngles_Iw1,RotAngles_Iw2
          COMPLEX                             ::eta_Iw1,eta_Iw2,zeta_Iw1,zeta_Iw2
          INTEGER                             ::Npanels,Ipanel,Iinteg,NpanWLin
          REAL                                ::w1,w2,delw,sumw,k1,k2
          REAL                                ::rho,depth
          INTEGER                             ::Iwline,Ibody,Itheta0
          COMPLEX,DIMENSION(2)                ::Radpot
          COMPLEX,DIMENSION(3,2)              ::Radvel
          COMPLEX,Dimension(2)                ::Pressure_I
          INTEGER                             ::InterpSwitch
          REAL,DIMENSION(3)                   ::XM,Normal_Vect
          COMPLEX,DIMENSION(2)                ::PHI_I
          COMPLEX,DIMENSION(3,2)              ::GradPHI_I
          COMPLEX,DIMENSION(2)                ::dnPhi_I
          COMPLEX                             ::QB1_M,QB1_P
          COMPLEX                             ::QB2_M,QB2_P
          COMPLEX,DIMENSION(3)                ::RIw1_n0,RIw2_n0
          COMPLEX,DIMENSION(3)                ::dtDisp_GradPhi_I_Iw1, &
                                                dtDisp_GradPhi_I_Iw2
          COMPLEX,DIMENSION(3)                ::dGAMMA_Vect,QB2V_M,QB2V_P
          INTEGER                             ::Isym,Iterm
          INTEGER,DIMENSION(2)                ::Ipanelinit,Iwlineinit

          rho=Env%rho
          depth=Env%depth
          w1=Qfreq%wQ(Iw1,Ibeta1)
          w2=Qfreq%wQ(Iw2,Ibeta2)
          delw=w1-w2
          sumw=w1+w2
          k1=Qfreq%kQ(Iw1,Ibeta1)
          k2=Qfreq%kQ(Iw2,Ibeta2)


          RotAngles_Iw1=RotAnglesQ(Iw1,Ibeta1,:)
          RotAngles_Iw2=RotAnglesQ(Iw2,Ibeta2,:)


          InterpSwitch=Qfreq%InterpPotSwitch(Ibeta1)                      &
                        +Qfreq%InterpPotSwitch(Ibeta2)

          Npanels=Mesh%Npanels
          NpanWlin=Npanels+WLine%NWLineSeg
          DO Iinteg=1,Nintegration
           QTFHasbo(Iinteg,:,:)=CZERO
           !integration over body panels

          DO Ipanel=1,Npanels
            IF (Mesh%XM(3, Ipanel)<0.) THEN  !dont compute at lid panels
              XM=Mesh%XM(1:3, Ipanel)
              Normal_Vect=Mesh%N(1:3,Ipanel)

              DO Isym=1,1+Mesh%Isym
              IF (Isym==2) THEN
                 XM(2)=-XM(2)
                 Normal_Vect(2)=-Normal_Vect(2)
              ENDIF
              Ipanelinit=Function_SymIpanel(NpanWlin,Npanels,Isym)
              !prepare variables Iw1,Iw2
              TotPot_Iw1=datPotVelQ%TotPot(Ipanelinit(1)+Ipanel,Ibeta1,Iw1)
              TotVel_Iw1=datPotVelQ%TotVel(Ipanelinit(1)+Ipanel,1:3,Ibeta1,Iw1)

              TotPot_Iw2=datPotVelQ%TotPot(Ipanelinit(1)+Ipanel,Ibeta2,Iw2)
              TotVel_Iw2=datPotVelQ%TotVel(Ipanelinit(1)+Ipanel,1:3,Ibeta2,Iw2)

              Bdisplace_Iw1=BdisplaceQ(Iw1,Ibeta1,Ipanelinit(1)+Ipanel,1:3)
              Bdisplace_Iw2=BdisplaceQ(Iw2,Ibeta2,Ipanelinit(1)+Ipanel,1:3)

              dtBdisplace_Iw1=-II*w1*Bdisplace_Iw1
              dtBdisplace_Iw2=-II*w2*Bdisplace_Iw2

              !interpolating radpot and radvel at delw and sumw
              Radpot(1:2)=INTERP_RADIATION_POTENTIAL                       &
                    (Nw,w,datPotVelQ%RadPot(Ipanelinit(1)+Ipanel,Iinteg,:),&
                        InterpSwitch,delw,sumw)
              RadVel(:,1:2)=INTERP_RADIATION_VELOCITY                      &
                  (Nw,w,datPotVelQ%RadVel(Ipanelinit(1)+Ipanel,:,Iinteg,:),&
                        InterpSwitch,delw,sumw)
              !------------------------------------------------------
              !Compute second order Froude-krylov force
              Phi_I=CALC_SECONDORDER_INCOMING_POTENTIAL                     &
                   (Env,w1,w2,k1,k2,beta(Ibeta1),beta(Ibeta2),XM)
              Pressure_I(1)=II*delw*Phi_I(1) !-dtPhi_M
              Pressure_I(2)=II*sumw*Phi_I(2) !-dtPhi_P

              QTFHasbo(Iinteg,1,1)=QTFHasbo(Iinteg,1,1)                     &
                 -rho*Pressure_I(1)*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
              QTFHasbo(Iinteg,2,1)=QTFHasbo(Iinteg,2,1)                     &
                 -rho*Pressure_I(2)*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)

              !------------------------------------------------------
              !Compute second order diffraction force: body force contrib
              !-----------------------------------------------------
              !Compute iw*rho*Integ(dnPhiI*Psi)dS
              GradPhi_I=CALC_SECONDORDER_INCOMING_VELOCITY                &
                       (k1,k2,beta(Ibeta1),beta(Ibeta2),depth,XM(3),PHI_I)
              dnPhi_I(1)=DOT_PRODUCT_COMPLEX(                             &
                      GradPhi_I(:,1),CMPLX(Normal_Vect,0),3) !DnPhi^M
              dnPhi_I(2)=DOT_PRODUCT_COMPLEX(                             &
                      GradPhi_I(:,2),CMPLX(Normal_Vect,0),3) !DnPhi^P
              !
              QTFHasbo(Iinteg,1,2)=QTFHasbo(Iinteg,1,2)                   &
                    +II*delw*rho*dnPhi_I(1)*Radpot(1)*Mesh%A(Ipanel)
              QTFHasbo(Iinteg,2,2)=QTFHasbo(Iinteg,2,2)                   &
                    +II*sumw*rho*dnPhi_I(2)*Radpot(2)*Mesh%A(Ipanel)

              !-------------------------------------------------------
              ! term: product(dtBdisplace-grad_Phi,R(n0))*Psi
              ITheta0= INT((Iinteg-1)/6)*3
              RIw1_n0=CROSS_PRODUCT_COMPLEX(                              &
                         RotAngles_Iw1(Itheta0+1:Itheta0+3),              &
                         CMPLX(Normal_Vect,0.))
              RIw2_n0=CROSS_PRODUCT_COMPLEX(                              &
                        RotAngles_Iw2(Itheta0+1:Itheta0+3),               &
                        CMPLX(Normal_Vect,0.))
              dtDisp_GradPhi_I_Iw1=dtBdisplace_Iw1(:)-TotVel_Iw1(:)
              dtDisp_GradPhi_I_Iw2=dtBdisplace_Iw2(:)-TotVel_Iw2(:)

              QB1_M=DOT_PRODUCT_DIFF_BIHARM(                             &
                   dtDisp_GradPhi_I_Iw1,dtDisp_GradPhi_I_Iw2,            &
                   RIw1_n0,RIw2_n0,3)
              QB1_P=DOT_PRODUCT_SUM_BIHARM(                              &
                   dtDisp_GradPhi_I_Iw1,dtDisp_GradPhi_I_Iw2,            &
                   RIw1_n0,RIw2_n0,3)
!
              QTFHasbo(Iinteg,1,3)=QTFHasbo(Iinteg,1,3)                  &
                    -II*delw*rho*QB1_M*Radpot(1)*Mesh%A(Ipanel)
              QTFHasbo(Iinteg,2,3)=QTFHasbo(Iinteg,2,3)                  &
                    -II*sumw*rho*QB1_P*Radpot(2)*Mesh%A(Ipanel)
             !-------------------------------------------------------
              ! term: product(Grad_PHI,Grad_Psi).product(Bdisplace,n0)
              QB2_M=DOT_PRODUCT_COMPLEX(                                 &
                    COMPLEX_CONJUGATE_VECT(TotVel_Iw2,3),RadVel(:,1),3)  &
                    *DOT_PRODUCT_COMPLEX(Bdisplace_Iw1,                  &
                                        CMPLX(Normal_vect,0),3)          &
                  +DOT_PRODUCT_COMPLEX(TotVel_Iw1(:),RadVel(:,1),3)      &
                    *DOT_PRODUCT_COMPLEX(COMPLEX_CONJUGATE_VECT(         &
                    Bdisplace_Iw2(:),3),CMPLX(Normal_Vect,0.),3)
              QB2_P=DOT_PRODUCT_COMPLEX(TotVel_Iw2,RadVel(:,2),3)        &
                    *DOT_PRODUCT_COMPLEX(Bdisplace_Iw1,                  &
                                        CMPLX(Normal_vect,0),3)          &
                  +DOT_PRODUCT_COMPLEX(TotVel_Iw1,RadVel(:,2),3)         &
                    *DOT_PRODUCT_COMPLEX(Bdisplace_Iw2,                  &
                                        CMPLX(Normal_Vect,0.),3)
              QTFHasbo(Iinteg,1,4)=QTFHasbo(Iinteg,1,4)                  &
                    +II*delw*rho*(QB2_M/2)*Mesh%A(Ipanel)
              QTFHasbo(Iinteg,2,4)=QTFHasbo(Iinteg,2,4)                  &
                    +II*sumw*rho*(QB2_P/2)*Mesh%A(Ipanel)

              ! term: product(Grad_Psi,Bdisplace).product(gradPhi,n0)
              QB2_M=DOT_PRODUCT_COMPLEX(                                  &
                      COMPLEX_CONJUGATE_VECT(TotVel_Iw2,3),               &
                                        CMPLX(Normal_Vect,0.),3)          &
                     *DOT_PRODUCT_COMPLEX(RadVel(:,1),Bdisplace_Iw1,3)    &
                   +DOT_PRODUCT_COMPLEX(                                  &
                       TotVel_Iw1,CMPLX(Normal_Vect,0),3)                 &
                    *DOT_PRODUCT_COMPLEX(                                 &
                    COMPLEX_CONJUGATE_VECT(Bdisplace_Iw2,3),RadVel(:,1),3)
              QB2_P=DOT_PRODUCT_COMPLEX(TotVel_Iw2,CMPLX(Normal_Vect,0),3)&
                     *DOT_PRODUCT_COMPLEX(RadVel(:,2),Bdisplace_Iw1,3)    &
                   +DOT_PRODUCT_COMPLEX(TotVel_Iw1,CMPLX(Normal_Vect,0),3)&
                    *DOT_PRODUCT_COMPLEX(Bdisplace_Iw2,RadVel(:,2),3)

              QTFHasbo(Iinteg,1,5)=QTFHasbo(Iinteg,1,5)                   &
                    -II*delw*rho*(QB2_M/2)*Mesh%A(Ipanel)
              QTFHasbo(Iinteg,2,5)=QTFHasbo(Iinteg,2,5)                   &
                    -II*sumw*rho*(QB2_P/2)*Mesh%A(Ipanel)

            ENDDO
            ENDIF
           ENDDO
           !   print*,QTFHasbo(Iinteg,1)
           DO Iwline=1,Wline%NWLineSeg
             Normal_Vect=Mesh%N(1:3,Wline%IndexPanel(Iwline))
            DO Isym=1,1+Mesh%Isym
             Iwlineinit=Function_SymIwline(NpanWlin,Npanels,WLine%NWLineSeg,Isym)
             !preparation var in Iw1 and Iw2
                Bdisplace_Iw1=BdisplaceQ(Iw1,Ibeta1,Iwlineinit(1)+Iwline,1:3)
                Bdisplace_Iw2=BdisplaceQ(Iw2,Ibeta2,Iwlineinit(1)+Iwline,1:3)
                TotVel_Iw1=datPotVelQ%TotVel(Iwlineinit(1)+Iwline,1:3,Ibeta1,Iw1)
                TotVel_Iw2=datPotVelQ%TotVel(Iwlineinit(1)+Iwline,1:3,Ibeta2,Iw2)
             !
             IF (Isym==2) THEN
                  Normal_Vect(2)=-Normal_Vect(2)
             ENDIF

             dGAMMA_Vect(1)=CMPLX( Normal_Vect(2),0.)
             dGAMMA_Vect(2)=CMPLX(-Normal_Vect(1),0.)
             dGAMMA_Vect(3)=CMPLX( 0.,0.)
             dGAMMA_Vect=dGAMMA_Vect*Wline%SegLength(Iwline)

             !interpolating radpot and radvel at delw and sumw
             Radpot(1:2)=INTERP_RADIATION_POTENTIAL                       &
                    (Nw,w,datPotVelQ%RadPot(Iwlineinit(1)+Iwline,Iinteg,:),&
                        InterpSwitch,delw,sumw)
             !-----------------------------------------------------
               QB2V_M=CROSS_PRODUCT_COMPLEX(Radpot(1)*Bdisplace_Iw1    &
                                ,COMPLEX_CONJUGATE_VECT( TotVel_Iw2,3))&
                     +CROSS_PRODUCT_COMPLEX(Radpot(1)*                 &
                      COMPLEX_CONJUGATE_VECT(Bdisplace_Iw2,3)          &
                                ,TotVel_Iw1)
               QB2V_P=CROSS_PRODUCT_COMPLEX(Radpot(2)*Bdisplace_Iw1    &
                     ,TotVel_Iw2)              &
                     +CROSS_PRODUCT_COMPLEX(Radpot(2)*Bdisplace_Iw2    &
                     ,TotVel_Iw1)

              QTFHasbo(Iinteg,1,6)=QTFHasbo(Iinteg,1,6)                &
                     -0.5*II*delw*rho                                  &
                     *DOT_PRODUCT_COMPLEX(QB2V_M,dGAMMA_Vect,3)
              QTFHasbo(Iinteg,2,6)=QTFHasbo(Iinteg,2,6)                &
                     -0.5*II*sumw*rho                                  &
                     *DOT_PRODUCT_COMPLEX(QB2V_P,dGAMMA_Vect,3)
              !----------------------------------------------------

            ENDDO
           ENDDO
          ENDDO
          DO Iterm=1,6
          QTFHasbo(:,:,Iterm)=QTFHasbo(:,:,Iterm)/2 !TRANSFER FUNCTION = BICHROMATIC FORCES / 2?
          QTFHasbo(:,:,7)=QTFHasbo(:,:,7)+QTFHasbo(:,:,Iterm)! Hasbo total
          ENDDO
  END SUBROUTINE

  SUBROUTINE COMPUTATION_QTF_POTENTIAL_FREESURFACEFORCE                       &
                                   (Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,       &
                                    MeshFS,NwQ,Nbeta,NPFlow,Nbodies,          &
                                    env,datPotVelQ,Nw,w,Qfreq,beta,QTFHASFS)
          !Input/output
          INTEGER,                              INTENT(IN) :: Iw1,Iw2,Ibeta1,Ibeta2
          INTEGER,                              INTENT(IN) :: Nintegration,NPFlow
          INTEGER,                              INTENT(IN) :: NwQ,Nbeta,Nbodies,Nw
          TYPE(TMeshFS),                        INTENT(IN) :: MeshFS
          TYPE(TPotVel),                        INTENT(IN) :: datPotVelQ
          TYPE(TQfreq),                         INTENT(IN) :: Qfreq
          TYPE(TEnvironment),                   INTENT(IN) :: Env
          REAL,DIMENSION(Nw),                   INTENT(IN) :: w
          REAL,DIMENSION(Nbeta),                INTENT(IN) :: beta
          COMPLEX,DIMENSION(Nintegration,2,10),  INTENT(OUT):: QTFHASFS
          !Local
          COMPLEX                             ::TotPot_Iw1,TotPot_Iw2
          COMPLEX,DIMENSION(3)                ::TotVel_Iw1,TotVel_Iw2
          COMPLEX                             ::IncPot_Iw1,IncPot_Iw2
          COMPLEX,DIMENSION(3)                ::IncVel_Iw1,IncVel_Iw2
          COMPLEX                             ::PerPot_Iw1,PerPot_Iw2
          COMPLEX,DIMENSION(3)                ::PerVel_Iw1,PerVel_Iw2
          COMPLEX,DIMENSION(2)                ::Radpot
          COMPLEX,DIMENSION(3,2)              ::Radvel

          REAL                                ::w1,w2,delw,sumw,k1,k2
          REAL                                ::rho,depth,g
          INTEGER                             ::InterpSwitch
          INTEGER                             ::Npanels,Ipanel,Iinteg,NpanWLin
          REAL,DIMENSION(3)                   ::XM,Normal_Vect
          INTEGER                             ::Isym,Iterm
          INTEGER,DIMENSION(2)                ::Ipanelinit,Iwlineinit
          COMPLEX                             ::QFD1_M,QFD1_P
          COMPLEX                             ::QFD2_M,QFD2_P
          INTEGER                             ::Iwline
          REAL                                ::dGAMMA

          InterpSwitch=Qfreq%InterpPotSwitch(Ibeta1)                      &
                        +Qfreq%InterpPotSwitch(Ibeta2)

          rho=Env%rho
          depth=Env%depth
          g=Env%g
          w1=Qfreq%wQ(Iw1,Ibeta1)
          w2=Qfreq%wQ(Iw2,Ibeta2)
          delw=w1-w2
          sumw=w1+w2
          k1=Qfreq%kQ(Iw1,Ibeta1)
          k2=Qfreq%kQ(Iw2,Ibeta2)

          Npanels=MeshFS%Mesh%Npanels
          NpanWlin=Npanels+MeshFS%BdyLine%NWLineSeg

          DO Iinteg=1,Nintegration
           QTFHASFS(Iinteg,:,:)=CZERO
           !integration over body panels
           DO Ipanel=1,Npanels
              XM=MeshFS%Mesh%XM(1:3, Ipanel)
              Normal_Vect=MeshFS%Mesh%N(1:3,Ipanel)

              DO Isym=1,1+MeshFS%Mesh%Isym
                 IF (Isym==2) THEN
                    XM(2)=-XM(2)
                    Normal_Vect(2)=-Normal_Vect(2)
                 ENDIF
                 Ipanelinit=Function_SymIpanel(NpanWlin,Npanels,Isym)
                 !prepare variables Iw1,Iw2
                 TotPot_Iw1=datPotVelQ%TotPot(Ipanelinit(1)+Ipanel,Ibeta1,Iw1)
                 TotVel_Iw1=datPotVelQ%TotVel(Ipanelinit(1)+Ipanel,1:3,Ibeta1,Iw1)
                 TotPot_Iw2=datPotVelQ%TotPot(Ipanelinit(1)+Ipanel,Ibeta2,Iw2)
                 TotVel_Iw2=datPotVelQ%TotVel(Ipanelinit(1)+Ipanel,1:3,Ibeta2,Iw2)

                 IncPot_Iw1=datPotVelQ%IncPot(Ipanelinit(1)+Ipanel,Ibeta1,Iw1)
                 IncVel_Iw1=datPotVelQ%IncVel(Ipanelinit(1)+Ipanel,1:3,Ibeta1,Iw1)
                 IncPot_Iw2=datPotVelQ%IncPot(Ipanelinit(1)+Ipanel,Ibeta2,Iw2)
                 IncVel_Iw2=datPotVelQ%IncVel(Ipanelinit(1)+Ipanel,1:3,Ibeta2,Iw2)

                 PerPot_Iw1=TotPot_Iw1-IncPot_Iw1
                 PerVel_Iw1=TotVel_Iw1-IncVel_Iw1
                 PerPot_Iw2=TotPot_Iw2-IncPot_Iw2
                 PerVel_Iw2=TotVel_Iw2-IncVel_Iw2

                 !interpolating radpot and radvel at delw and sumw
                 Radpot(1:2)=INTERP_RADIATION_POTENTIAL                       &
                       (Nw,w,datPotVelQ%RadPot(Ipanelinit(1)+Ipanel,Iinteg,:),&
                           InterpSwitch,delw,sumw)
                 RadVel(:,1:2)=INTERP_RADIATION_VELOCITY                      &
                     (Nw,w,datPotVelQ%RadVel(Ipanelinit(1)+Ipanel,:,Iinteg,:),&
                           InterpSwitch,delw,sumw)

                 !------------------------------------------------------------
                 !term QFD1:(gradPhi,gradPhiPer)+(gradPhiPer,GradPhiInc)

                 QFD1_M=II*delw*(                                             &
                        DOT_PRODUCT_COMPLEX(TotVel_Iw1,CONJG(PerVel_Iw2),3)   &
                       +DOT_PRODUCT_COMPLEX(PerVel_Iw1,CONJG(IncVel_Iw2),3)   &
                        )


                 QFD1_P=II*sumw*(                                             &
                        DOT_PRODUCT_COMPLEX(TotVel_Iw1,PerVel_Iw2,3)          &
                       +DOT_PRODUCT_COMPLEX(PerVel_Iw1,IncVel_Iw2,3)          &
                        )

                 QTFHASFS(Iinteg,1,1)= QTFHASFS(Iinteg,1,1)+                  &
                                       II*delw*rho/g*QFD1_M*Radpot(1)         &
                                       *MeshFS%Mesh%A(Ipanel)
                 QTFHASFS(Iinteg,2,1)= QTFHASFS(Iinteg,2,1)+                  &
                                       II*sumw*rho/g*QFD1_P*Radpot(2)         &
                                       *MeshFS%Mesh%A(Ipanel)

                 !term QFD1:Phi*dttdzPhiPer+PhiPer*dttdzPhiInc
                 QFD1_M=-II*w1/2/g*(                                          &
                          TotPot_Iw1*(-w2**2*CONJG(PerVel_Iw2(3)))            &
                         +PerPot_Iw1*(-w2**2*CONJG(IncVel_Iw2(3)))            &
                         )
                 QFD1_P=-II*w1/2/g*(                                          &
                          TotPot_Iw1*(-w2**2*PerVel_Iw2(3))                   &
                         +PerPot_Iw1*(-w2**2*IncVel_Iw2(3))                   &
                         )

                 QTFHASFS(Iinteg,1,2)= QTFHASFS(Iinteg,1,2)+                  &
                                       II*delw*rho/g*QFD1_M*Radpot(1)         &
                                       *MeshFS%Mesh%A(Ipanel)
                 QTFHASFS(Iinteg,2,2)= QTFHASFS(Iinteg,2,2)+                  &
                                       II*sumw*rho/g*QFD1_P*Radpot(2)         &
                                       *MeshFS%Mesh%A(Ipanel)

                 QFD1_M= II*w2/2/g*(                                          &
                         CONJG(TotPot_Iw2)*(-w1**2*PerVel_Iw1(3))             &
                        +CONJG(PerPot_Iw2)*(-w1**2*IncVel_Iw1(3))             &
                         )
                 QFD1_P=-II*w2/2/g*(                                          &
                          TotPot_Iw2*(-w1**2*PerVel_Iw1(3))                   &
                         +PerPot_Iw2*(-w1**2*IncVel_Iw1(3))                   &
                         )

                 QTFHASFS(Iinteg,1,3)= QTFHASFS(Iinteg,1,3)+                  &
                                       II*delw*rho/g*QFD1_M*Radpot(1)         &
                                       *MeshFS%Mesh%A(Ipanel)
                 QTFHASFS(Iinteg,2,3)= QTFHASFS(Iinteg,2,3)+                  &
                                       II*sumw*rho/g*QFD1_P*Radpot(2)         &
                                       *MeshFS%Mesh%A(Ipanel)
                !term QFD2: direct implementation
                !dtPhi*dzzPhiPer+dtPhi_Per*dzzPhiInc
                 QFD2_M= 0.5*(                                                &
                          (-II*w1*TotPot_Iw1)*k2**2*CONJG(PerPot_Iw2)         &
                         +(-II*w1*PerPot_Iw1)*k2**2*CONJG(IncPot_Iw2)         &
                         )
                 QFD2_P= 0.5*(                                                &
                          (-II*w1*TotPot_Iw1)*k2**2*PerPot_Iw2                &
                         +(-II*w1*PerPot_Iw1)*k2**2*IncPot_Iw2                &
                         )

                 QTFHASFS(Iinteg,1,4)= QTFHASFS(Iinteg,1,4)+                  &
                                       II*delw*rho/g*QFD2_M*Radpot(1)         &
                                       *MeshFS%Mesh%A(Ipanel)
                 QTFHASFS(Iinteg,2,4)= QTFHASFS(Iinteg,2,4)+                  &
                                       II*sumw*rho/g*QFD2_P*Radpot(2)         &
                                       *MeshFS%Mesh%A(Ipanel)
                 QFD2_M= 0.5*(                                                &
                          (II*w2*CONJG(TotPot_Iw2))*k1**2*PerPot_Iw1          &
                         +(II*w2*CONJG(PerPot_Iw2))*k1**2*IncPot_Iw1          &
                         )
                 QFD2_P= 0.5*(                                                &
                          (-II*w2*TotPot_Iw2)*k1**2*PerPot_Iw1                &
                         +(-II*w2*PerPot_Iw2)*k1**2*IncPot_Iw1                &
                         )

                 QTFHASFS(Iinteg,1,5)= QTFHASFS(Iinteg,1,5)+                  &
                                       II*delw*rho/g*QFD2_M*Radpot(1)         &
                                       *MeshFS%Mesh%A(Ipanel)
                 QTFHASFS(Iinteg,2,5)= QTFHASFS(Iinteg,2,5)+                  &
                                       II*sumw*rho/g*QFD2_P*Radpot(2)         &
                                       *MeshFS%Mesh%A(Ipanel)

                !term QFD2:dtPhi*dzzPhiPer+dtPhi_Per*dzzPhiInc
                !Green theorem to express second derivative in
                !the first derivative only
                QFD2_M= 0.5*(                                               &
                        -II*w1*(                                            &
                         ( RadVel(1,1)*TotPot_Iw1+RadPot(1)*TotVel_Iw1(1))  &
                           *CONJG(PerVel_Iw2(1))                            &
                         +(RadVel(2,1)*TotPot_Iw1+RadPot(1)*TotVel_Iw1(2))  &
                           *CONJG(PerVel_Iw2(2))                            &
                              )                                             &
                        -II*w1*(                                            &
                         ( RadVel(1,1)*PerPot_Iw1+RadPot(1)*PerVel_Iw1(1))  &
                           *CONJG(IncVel_Iw2(1))                            &
                         +(RadVel(2,1)*PerPot_Iw1+RadPot(1)*PerVel_Iw1(2))  &
                           *CONJG(IncVel_Iw2(2))                            &
                              )                                             &
                        )
                QFD2_P= 0.5*(                                               &
                        -II*w1*(                                            &
                         ( RadVel(1,2)*TotPot_Iw1+RadPot(2)*TotVel_Iw1(1))  &
                           *PerVel_Iw2(1)                                   &
                         +(RadVel(2,2)*TotPot_Iw1+RadPot(2)*TotVel_Iw1(2))  &
                           *PerVel_Iw2(2)                                   &
                              )                                             &
                        -II*w1*(                                            &
                         ( RadVel(1,2)*PerPot_Iw1+RadPot(2)*PerVel_Iw1(1))  &
                           *IncVel_Iw2(1)                                   &
                         +(RadVel(2,2)*PerPot_Iw1+RadPot(2)*PerVel_Iw1(2))  &
                           *IncVel_Iw2(2)                                   &
                              )                                             &
                        )

                QTFHASFS(Iinteg,1,6)= QTFHASFS(Iinteg,1,6)+                 &
                                      II*delw*rho/g*QFD2_M                  &
                                      *MeshFS%Mesh%A(Ipanel)
                QTFHASFS(Iinteg,2,6)= QTFHASFS(Iinteg,2,6)+                 &
                                      II*sumw*rho/g*QFD2_P                  &
                                      *MeshFS%Mesh%A(Ipanel)

                 QFD2_M= 0.5*(                                              &
                        II*w2*(                                             &
                         ( RadVel(1,1)*CONJG(TotPot_Iw2)+                   &
                                 RadPot(1)*CONJG(TotVel_Iw2(1)))            &
                           *PerVel_Iw1(1)                                   &
                         +(RadVel(2,1)*CONJG(TotPot_Iw2)+                   &
                                 RadPot(1)*CONJG(TotVel_Iw2(2)))            &
                           *PerVel_Iw1(2)                                   &
                              )                                             &
                        +II*w2*(                                            &
                         ( RadVel(1,1)*CONJG(PerPot_Iw2)+                   &
                                 RadPot(1)*CONJG(PerVel_Iw2(1)))            &
                           *IncVel_Iw1(1)                                   &
                         +(RadVel(2,1)*CONJG(PerPot_Iw2)+                   &
                                 RadPot(1)*CONJG(PerVel_Iw2(2)))            &
                           *IncVel_Iw1(2)                                   &
                              )                                             &
                        )
                QFD2_P= 0.5*(                                               &
                        -II*w2*(                                            &
                         ( RadVel(1,2)*TotPot_Iw2+                          &
                                 RadPot(2)*TotVel_Iw2(1))                   &
                           *PerVel_Iw1(1)                                   &
                         +(RadVel(2,2)*TotPot_Iw2+                          &
                                 RadPot(2)*TotVel_Iw2(2))                   &
                           *PerVel_Iw1(2)                                   &
                              )                                             &
                        -II*w2*(                                            &
                         ( RadVel(1,2)*PerPot_Iw2+                          &
                                 RadPot(2)*PerVel_Iw2(1))                   &
                           *IncVel_Iw1(1)                                   &
                         +(RadVel(2,2)*PerPot_Iw2+                          &
                                 RadPot(2)*PerVel_Iw2(2))                   &
                           *IncVel_Iw1(2)                                   &
                              )                                             &
                        )

                QTFHASFS(Iinteg,1,7)= QTFHASFS(Iinteg,1,7)+                 &
                                      II*delw*rho/g*QFD2_M                  &
                                      *MeshFS%Mesh%A(Ipanel)
                QTFHASFS(Iinteg,2,7)= QTFHASFS(Iinteg,2,7)+                 &
                                      II*sumw*rho/g*QFD2_P                  &
                                      *MeshFS%Mesh%A(Ipanel)

              ENDDO

           ENDDO

           DO Iwline=1,MeshFS%BdyLine%NWLineSeg
             Normal_Vect(1:2)=MeshFS%BdyLineNormal(Iwline,:)
             Normal_Vect(3)=0.
            DO Isym=1,1+MeshFS%Mesh%Isym
              IF (Isym==2) THEN
                  Normal_Vect(2)=-Normal_Vect(2)
             ENDIF
             Iwlineinit=Function_SymIwline(NpanWlin,Npanels,                  &
                                               MeshFS%BdyLine%NWLineSeg,Isym)
             !preparation var in Iw1 and Iw2
             TotPot_Iw1=datPotVelQ%TotPot(Iwlineinit(1)+Iwline,Ibeta1,Iw1)
             TotVel_Iw1=datPotVelQ%TotVel(Iwlineinit(1)+Iwline,1:3,Ibeta1,Iw1)
             TotPot_Iw2=datPotVelQ%TotPot(Iwlineinit(1)+Iwline,Ibeta2,Iw2)
             TotVel_Iw2=datPotVelQ%TotVel(Iwlineinit(1)+Iwline,1:3,Ibeta2,Iw2)

             IncPot_Iw1=datPotVelQ%IncPot(Iwlineinit(1)+Iwline,Ibeta1,Iw1)
             IncVel_Iw1=datPotVelQ%IncVel(Iwlineinit(1)+Iwline,1:3,Ibeta1,Iw1)
             IncPot_Iw2=datPotVelQ%IncPot(Iwlineinit(1)+Iwline,Ibeta2,Iw2)
             IncVel_Iw2=datPotVelQ%IncVel(Iwlineinit(1)+Iwline,1:3,Ibeta2,Iw2)

             PerPot_Iw1=TotPot_Iw1-IncPot_Iw1
             PerVel_Iw1=TotVel_Iw1-IncVel_Iw1
             PerPot_Iw2=TotPot_Iw2-IncPot_Iw2
             PerVel_Iw2=TotVel_Iw2-IncVel_Iw2

             !interpolating radpot and radvel at delw and sumw
             Radpot(1:2)=INTERP_RADIATION_POTENTIAL                       &
                    (Nw,w,datPotVelQ%RadPot(Iwlineinit(1)+Iwline,Iinteg,:),&
                        InterpSwitch,delw,sumw)
             !-----------------------------------------------------

             dGAMMA=MeshFS%BdyLine%SegLength(Iwline)

             QFD2_M=0.5*(-II*w1*TotPot_Iw1*                                      &
                          (DOT_PRODUCT_COMPLEX(CONJG(PerVel_Iw2),                &
                                                       CMPLX(Normal_Vect,0.),3)) &
                          -II*w1*PerPot_Iw1*                                     &
                          (DOT_PRODUCT_COMPLEX(CONJG(IncVel_Iw2),                &
                                                       CMPLX(Normal_Vect,0.),3)) &
                         )
             QFD2_P=0.5*(-II*w1*TotPot_Iw1*                                      &
                          (DOT_PRODUCT_COMPLEX(PerVel_Iw2,                       &
                                                       CMPLX(Normal_Vect,0.),3)) &
                          -II*w1*PerPot_Iw1*                                     &
                          (DOT_PRODUCT_COMPLEX(IncVel_Iw2,                       &
                                                       CMPLX(Normal_Vect,0.),3)) &
                         )


             QTFHASFS(Iinteg,1,8)= QTFHASFS(Iinteg,1,8)-                         &
                                      II*delw*rho/g*QFD2_M*Radpot(1)*dGAMMA
             QTFHASFS(Iinteg,2,8)= QTFHASFS(Iinteg,2,8)-                         &
                                      II*sumw*rho/g*QFD2_P*Radpot(2)*dGAMMA

             QFD2_M=0.5*( II*w2*CONJG(TotPot_Iw2)*                               &
                          (DOT_PRODUCT_COMPLEX(PerVel_Iw1,                       &
                                                       CMPLX(Normal_Vect,0.),3)) &
                          +II*w2*CONJG(PerPot_Iw2)*                              &
                          (DOT_PRODUCT_COMPLEX(IncVel_Iw1,                       &
                                                       CMPLX(Normal_Vect,0.),3)) &
                         )
             QFD2_P=0.5*(-II*w2*TotPot_Iw2*                                      &
                          (DOT_PRODUCT_COMPLEX(PerVel_Iw1,                       &
                                                       CMPLX(Normal_Vect,0.),3)) &
                          -II*w2*PerPot_Iw2*                                     &
                          (DOT_PRODUCT_COMPLEX(IncVel_Iw1,                       &
                                                       CMPLX(Normal_Vect,0.),3)) &
                         )

             QTFHASFS(Iinteg,1,9)= QTFHASFS(Iinteg,1,9)-                         &
                                      II*delw*rho/g*QFD2_M*Radpot(1)*dGAMMA
             QTFHASFS(Iinteg,2,9)= QTFHASFS(Iinteg,2,9)-                         &
                                      II*sumw*rho/g*QFD2_P*Radpot(2)*dGAMMA


             !----------------------------------------------------
            ENDDO
           ENDDO
          ENDDO
          DO Iterm=1,9
             QTFHASFS(:,:,Iterm)=QTFHASFS(:,:,Iterm)/2 !TRANSFER FUNCTION = BICHROMATIC FORCES / 2?
             IF (Iterm.LE.5) THEN!! uses direct implementation for QFD2
                 QTFHASFS(:,:,10)=QTFHASFS(:,:,10)+QTFHASFS(:,:,Iterm)! HASFS total
             ENDIF
          ENDDO

  END SUBROUTINE
!!-------------------------- HASFS ASYMPTOTIC --------------------------------
  SUBROUTINE COMPUTATION_QTF_POTENTIAL_FREESURFACEFORCE_ASYMP                 &
                                   (Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,       &
                                    NwQ,Nbeta,env,Nw,w,Qfreq,                 &
                                    beta,Mesh,ASYMP_PARAM,SourceDistr,        &
                                    QTFHASFS_ASYMP)
  !INPUT/OUTPUT
  INTEGER,                              INTENT(IN) :: Iw1,Iw2,Ibeta1,Ibeta2
  INTEGER,                              INTENT(IN) :: Nintegration
  INTEGER,                              INTENT(IN) :: NwQ,Nbeta,Nw
  TYPE(TQfreq),                         INTENT(IN) :: Qfreq
  TYPE(TEnvironment),                   INTENT(IN) :: Env
  REAL,DIMENSION(Nw),                   INTENT(IN) :: w
  REAL,DIMENSION(Nbeta),                INTENT(IN) :: beta
  TYPE(TMesh),                          INTENT(IN) :: Mesh
  TYPE(TASYMP),                         INTENT(IN) :: ASYMP_PARAM
  TYPE(TSOURCEQ),                       INTENT(IN) :: SourceDistr
  COMPLEX,DIMENSION(Nintegration,2,3),  INTENT(OUT):: QTFHASFS_ASYMP
  !Local
  INTEGER                             ::Iinteg,Iterm,NRf,Nbessel,Npanels,InterpSwitch
  INTEGER                             ::Isym,Ipanel,Ibessel
  REAL                                ::w1,w2,delw,sumw,k1,k2,delk,sumk
  REAL                                ::rho,depth,g
  REAL,DIMENSION(ASYMP_PARAM%NR)      :: Rf !discretized finite domain free surface radius
  REAL,DIMENSION(2)                   ::DELSUMW
  COMPLEX,DIMENSION(2)                ::KAPPA1,KAPPA2
  COMPLEX,DIMENSION(2,2)              ::I_DF
  COMPLEX,DIMENSION(Mesh%Npanels*2**Mesh%Isym)   :: ZIG_Per_Iw1, ZIG_Per_Iw2
  COMPLEX,DIMENSION(Mesh%Npanels*2**Mesh%Isym,2) :: ZIG_Rad
  COMPLEX,DIMENSION(2,ASYMP_PARAM%NBESSEL+1)       :: CmSmPer_k1,CmSmPer_k2
  COMPLEX,DIMENSION(2,ASYMP_PARAM%NBESSEL+1)       :: CmSmRad_delk,CmSmRad_sumk
  COMPLEX,DIMENSION(2,ASYMP_PARAM%NBESSEL+1)       :: IR1l,IR2l
  COMPLEX,DIMENSION(2,ASYMP_PARAM%NBESSEL+3)       :: Ivartheta1l,Ivartheta2l
  INTEGER                             :: ufile

  Isym=Mesh%Isym
  Npanels=Mesh%Npanels
  Nbessel=ASYMP_PARAM%NBESSEL
  NRf=ASYMP_PARAM%NR
  Rf=ASYMP_PARAM%Rf

  rho=Env%rho
  depth=Env%depth
  g=Env%g
  w1=Qfreq%wQ(Iw1,Ibeta1)
  w2=Qfreq%wQ(Iw2,Ibeta2)
  delw=w1-w2
  sumw=w1+w2
  DELSUMW(1)=delw
  DELSUMW(2)=sumw
  k1=Qfreq%kQ(Iw1,Ibeta1)
  k2=Qfreq%kQ(Iw2,Ibeta2)
  delk=Fun_inverseDispersion(delw,depth,g)
  sumk=Fun_inverseDispersion(sumw,depth,g)
  KAPPA1(1)= II*delw*k1*k2 !for diff freq
  KAPPA1(2)=-II*sumw*k1*k2  !for sum freq
  KAPPA2   =Fun_KAPPA2_DIFFSUM(k1,k2,w1,w2,depth,g)

  InterpSwitch=Qfreq%InterpPotSwitch(Ibeta1)                      &
                        +Qfreq%InterpPotSwitch(Ibeta2)

  !perturbed source distribution
  ZIG_Per_Iw1(1:Npanels)=SourceDistr%ZIGB_Per(:,ibeta1,iw1)
  ZIG_Per_Iw2(1:Npanels)=SourceDistr%ZIGB_Per(:,ibeta2,iw2)
  IF (Isym.EQ.1) THEN !symmetric case
    ZIG_Per_Iw1(Npanels+1:2*Npanels)=SourceDistr%ZIGS_Per(:,ibeta1,iw1)
    ZIG_Per_Iw2(Npanels+1:2*Npanels)=SourceDistr%ZIGS_Per(:,ibeta2,iw2)
  ENDIF
  !Kochin coefficients for perturbed potential
  CmSmPer_k1=PREPARE_KOCHIN_COEFFICIENTS                       &
          (Isym,Npanels,Mesh%XM,Mesh%A,depth,Nbessel,k1,ZIG_Per_Iw1)
  CmSmPer_k2=PREPARE_KOCHIN_COEFFICIENTS                       &
          (Isym,Npanels,Mesh%XM,Mesh%A,depth,Nbessel,k2,ZIG_Per_Iw2)
  !Prepare Integral over radius of free surface (Rext,Infinity)
  DO Ibessel=1,Nbessel+1
   IR1l(:,Ibessel)=Fun_IR1l(Ibessel-1,k1,k2,delk,sumk,Rf,NRf)
   IR2l(:,Ibessel)=Fun_IR2l(Ibessel-1,k1,k2,delk,sumk,Rf,NRf)
   IF (delk.EQ.0)  IR1l(1,Ibessel)=0.
   IF (delk.EQ.0)  IR2l(1,Ibessel)=0.
  ENDDO

  DO Iinteg=1,Nintegration
     QTFHASFS_ASYMP(Iinteg,:,:)=CZERO

     !interpolating rad source distribution  at delw and sumw
     ZIG_Rad(:,1:2)=INTERP_RADIATION_SOURCE_DISTRIBUTION                 &
                    (Nw,w,SourceDistr%ZIGB_Rad(:,Iinteg,:),              &
                      SourceDistr%ZIGS_Rad(:,Iinteg,:),                  &
                      Npanels,Mesh%Isym,InterpSwitch,delw,sumw)
     !Kochin coefficients for radiation potential
     CmSmRad_delk=PREPARE_KOCHIN_COEFFICIENTS                           &
          (Isym,Npanels,Mesh%XM,Mesh%A,depth,Nbessel,delk,ZIG_Rad(:,1))
     CmSmRad_sumk=PREPARE_KOCHIN_COEFFICIENTS                           &
          (Isym,Npanels,Mesh%XM,Mesh%A,depth,Nbessel,sumk,ZIG_Rad(:,2))
     !Prepare Integral over vartheta (0,2pi)
     DO Ibessel=1,Nbessel+3
      Ivartheta1l(:,Ibessel)=Fun_IVartheta1l(Nbessel,beta(Ibeta1),CmSmPer_k2,   &
                                CmSmRad_delk,CmSmRad_sumk,Ibessel-2)
      Ivartheta2l(:,Ibessel)=Fun_IVartheta2l(Nbessel,beta(Ibeta2),CmSmPer_k1,   &
                                CmSmRad_delk,CmSmRad_sumk,Ibessel-2)

     ENDDO

     !all the integral terms both for diff and sum freq
     I_DF=Fun_IDF(w1,w2,k1,k2,delk,sumk,g,depth,Rf,NRf,Nbessel,           &
                                       Ivartheta1l,Ivartheta2l,IR1l,IR2l)
     ! !first integral term both for diff and sum freq
     QTFHASFS_ASYMP(Iinteg,:,1)=II*DELSUMW*rho/g*KAPPA1*I_DF(:,1)!(I_DF11+I_DF12)
     ! !second integral term both for diff and sum freq
     QTFHASFS_ASYMP(Iinteg,:,2)=II*DELSUMW*rho/g*KAPPA2*I_DF(:,2)!(I_DF21+I_DF22)
  ENDDO
  DO Iterm=1,2
     QTFHASFS_ASYMP(:,:,Iterm)=QTFHASFS_ASYMP(:,:,Iterm)/2 !TRANSFER FUNCTION = BICHROMATIC FORCES / 2?
     QTFHASFS_ASYMP(:,:,3)=QTFHASFS_ASYMP(:,:,3)+QTFHASFS_ASYMP(:,:,Iterm)! HASFS_ASYMP total
  ENDDO
  END SUBROUTINE

  !!---------------------------------------------------------------------
  FUNCTION CALC_SECONDORDER_INCOMING_POTENTIAL                         &
                  (Env,w1,w2,k1,k2,beta1,beta2,XMp) result(Phi_I)
          !input/output
          REAL,                INTENT(IN)::w1,w2,k1,k2,beta1,beta2
          REAL,DIMENSION(3),   INTENT(IN)::XMp !XM at a panel
          TYPE(TEnvironment),  INTENT(IN)::Env
          COMPLEX               ::QFI_M,QFI_P   !free surface forcing
          REAL                  ::delw,sumw,OM_M,OM_P,OM_1,OM_2
          REAL,DIMENSION(2)     ::delk_vect,sumk_vect
          REAL,DIMENSION(2)     ::k1_vect,k2_vect
          REAL,DIMENSION(3)     ::XM
          REAL                  ::abs_delk,abs_sumk
          REAL                  ::g,D
          COMPLEX,DIMENSION(2)  ::Phi_I
          g=Env%g
          D=Env%Depth
          XM(1)=XMp(1)-Env%Xeff
          XM(2)=XMp(2)-Env%Yeff
          XM(3)=XMp(3)
          delw= w1-w2
          sumw= w1+w2
          k1_vect(1)=k1*cos(beta1)
          k1_vect(2)=k1*sin(beta1)
          k2_vect(1)=k2*cos(beta2)
          k2_vect(2)=k2*sin(beta2)

          delk_vect=k1_vect-k2_vect
          sumk_vect=k1_vect+k2_vect

          OM_1=FUN_DISPERSION(k1,D,g)
          OM_2=FUN_DISPERSION(k2,D,g)

          QFI_M=II*g*g*exp(II*DOT_PRODUCT(delk_vect,XM(1:2)))*          &
             ( delw/w1/w2*                                              &
             (DOT_PRODUCT(k1_vect,k2_vect)+OM_1**2/g*OM_2**2/g)         &
             +0.5*((k1**2-(OM_1**2/g)**2)/w1-(k2**2-(OM_2**2/g)**2)/w2))
          QFI_P=II*g*g*exp(II*DOT_PRODUCT(sumk_vect,XM(1:2)))*          &
             ( sumw/w1/w2*                                              &
             (DOT_PRODUCT(k1_vect,k2_vect)-OM_1**2/g*OM_2**2/g)         &
             +0.5*((k1**2-(OM_1**2/g)**2)/w1+(k2**2-(OM_2**2/g)**2)/w2))

          abs_delk=SQRT(k1**2+k2**2-2*k1*k2*cos(beta1-beta2))
          abs_sumk=SQRT(k1**2+k2**2+2*k1*k2*cos(beta1-beta2))
          OM_M=FUN_DISPERSION(abs_delk,D,g)
          OM_P=FUN_DISPERSION(abs_sumk,D,g)
          IF (delw.GT.0.) THEN
          Phi_I(1)=QFI_M*CIH(abs_delk,XM(3),D)/(-delw**2+OM_M**2)
          ELSE
          Phi_I(1)=CMPLX(0.,0.)
          ENDIF
          Phi_I(2)=QFI_P*CIH(abs_sumk,XM(3),D)/(-sumw**2+OM_P**2)
  END FUNCTION

  FUNCTION CALC_SECONDORDER_INCOMING_VELOCITY                           &
                           (k1,k2,beta1,beta2,D,ZM,PHI_I) result(GradPhi)
          !input/output
          REAL,                 INTENT(IN)      ::k1,k2,beta1,beta2
          REAL,                 INTENT(IN)      ::D,ZM
          COMPLEX,DIMENSION(2), INTENT(IN)      ::PHI_I

          REAL                                  ::abs_delk,abs_sumk
          COMPLEX,DIMENSION(3,2)                ::GradPhi
          !dxPHI
          GradPhi(1,1)=II*(k1*cos(beta1)-k2*cos(beta2))*PHI_I(1)
          GradPhi(1,2)=II*(k1*cos(beta1)+k2*cos(beta2))*PHI_I(2)
          !dyPHI
          GradPhi(2,1)=II*(k1*sin(beta1)-k2*sin(beta2))*PHI_I(1)
          GradPhi(2,2)=II*(k1*sin(beta1)+k2*sin(beta2))*PHI_I(2)
          !dzPHI
          abs_delk=SQRT(k1**2+k2**2-2*k1*k2*cos(beta1-beta2))
          abs_sumk=SQRT(k1**2+k2**2+2*k1*k2*cos(beta1-beta2))
          IF (D .EQ. INFINITE_DEPTH) THEN
          GradPhi(3,1)=abs_delk*PHI_I(1)
          GradPhi(3,2)=abs_sumk*PHI_I(2)
          ELSE
          GradPhi(3,1)=abs_delk*tanh(abs_delk*(D+ZM))*PHI_I(1)
          GradPhi(3,2)=abs_sumk*tanh(abs_sumk*(D+ZM))*PHI_I(2)
          ENDIF
  END FUNCTION

  FUNCTION INTERP_RADIATION_SOURCE_DISTRIBUTION                                    &
              (Nw,w,ZIGB,ZIGS,Npanels,Isym,InterpSwitch,delw,sumw) RESULT(ZIG_RAD)

          INTEGER,                       INTENT(IN)::Nw,InterpSwitch
          INTEGER,                       INTENT(IN)::Npanels,Isym
          REAL, DIMENSION(Nw),           INTENT(IN)::w
          REAL,                          INTENT(IN)::delw,sumw
          COMPLEX,DIMENSION(Npanels,Nw), INTENT(IN)::ZIGB,ZIGS
          COMPLEX,DIMENSION(Npanels*2**Isym,2)     :: ZIG_RAD
          INTEGER                                  :: Ipanel

          DO Ipanel=1,Npanels
             IF (delw.GT.0.) THEN
               IF (InterpSwitch==0) THEN
                 ZIG_RAD(Ipanel,1)=ZIGB(Ipanel,Fun_closest(Nw,w,delw))
                 IF (Isym.EQ.1) THEN
                 ZIG_RAD(Npanels+Ipanel,1)=ZIGS(Ipanel,Fun_closest(Nw,w,delw))
                 ENDIF
               ELSE
                 ZIG_RAD(Ipanel,1)=FUN_INTERP1_COMPLEX(w,ZIGB(Ipanel,:),Nw,delw)
                 IF (Isym.EQ.1) THEN
                 ZIG_RAD(Npanels+Ipanel,1)=FUN_INTERP1_COMPLEX(w,ZIGS(Ipanel,:),Nw,delw)
                 ENDIF
               ENDIF
             ELSE
               ZIG_RAD(Ipanel,1)=CMPLX(0.,0.)
               IF (Isym.EQ.1) ZIG_RAD(Npanels+Ipanel,1)=CMPLX(0.,0.)
             ENDIF
             IF (sumw.LE.w(Nw)) THEN
               IF (InterpSwitch==0) THEN
                 ZIG_RAD(Ipanel,2)=ZIGB(Ipanel,Fun_closest(Nw,w,sumw))
                 IF (Isym.EQ.1) THEN
                 ZIG_RAD(Npanels+Ipanel,2)=ZIGS(Ipanel,Fun_closest(Nw,w,sumw))
                 ENDIF
               ELSE
                 ZIG_RAD(Ipanel,2)=FUN_INTERP1_COMPLEX(w,ZIGB(Ipanel,:),Nw,sumw)
                 IF (Isym.EQ.1) THEN
                 ZIG_RAD(Npanels+Ipanel,2)=FUN_INTERP1_COMPLEX(w,ZIGS(Ipanel,:),Nw,sumw)
                 ENDIF
               ENDIF
             ELSE
               ZIG_RAD(Ipanel,2)=CMPLX(0.,0.)
               IF (Isym.EQ.1) ZIG_RAD(Npanels+Ipanel,2)=CMPLX(0.,0.)
             ENDIF
          ENDDO
  END FUNCTION

  FUNCTION INTERP_RADIATION_POTENTIAL                                     &
                (Nw,w,RadPot,InterpSwitch,delw,sumw) RESULT(RadPotInterp)

          INTEGER,                   INTENT(IN)::Nw,InterpSwitch
          REAL, DIMENSION(Nw),       INTENT(IN)::w
          REAL,                      INTENT(IN)::delw,sumw
          COMPLEX,DIMENSION(Nw),     INTENT(IN)::RadPot
          COMPLEX,DIMENSION(2)                 ::RadPotInterp
          INTEGER                    :: Iprint
          IF (InterpSwitch==0) THEN
           IF (delw.GT.0.) THEN
                RadPotInterp(1)=RadPot(Fun_closest(Nw,w,delw))
           ELSE
                RadPotInterp(1)=CMPLX(0.,0.)
           ENDIF
           IF (sumw.LE.w(Nw)) THEN
                RadPotInterp(2)=RadPot(Fun_closest(Nw,w,sumw))
           ELSE
                RadPotInterp(2)=CMPLX(0.,0.)
           ENDIF
          ELSE

           IF (delw.GT.0.) THEN
              RadPotInterp(1)=FUN_INTERP1_COMPLEX(w,RadPot(:),Nw,delw) !for delta omega
           ELSE
              RadPotInterp(1)=CMPLX(0.,0.)
           ENDIF
           IF (sumw.LE.w(Nw)) THEN
              RadPotInterp(2)=FUN_INTERP1_COMPLEX(w,RadPot(:),Nw,sumw) !for sum omega
              ELSE
              RadPotInterp(2)=CMPLX(0.,0.)
           ENDIF
          ENDIF
  END FUNCTION

  FUNCTION INTERP_RADIATION_VELOCITY                                    &
                (Nw,w,RadVel,InterpSwitch,delw,sumw) RESULT(RadVelInterp)
          INTEGER,                   INTENT(IN)::Nw,InterpSwitch
          REAL, DIMENSION(Nw),       INTENT(IN)::w
          REAL,                      INTENT(IN)::delw,sumw
          COMPLEX,DIMENSION(3,Nw),   INTENT(IN)::RadVel
          COMPLEX,DIMENSION(3,2)               ::RadVelInterp
          INTEGER                    :: I

          IF (InterpSwitch==0) THEN
             IF (delw.GT.0.) THEN
                RadVelInterp(:,1)=RadVel(:,Fun_closest(Nw,w,delw))
             ELSE
                RadVelInterp(:,1)=CMPLX(0.,0.)
             ENDIF

             IF (sumw.LE.w(Nw)) THEN
                RadVelInterp(:,2)=RadVel(:,Fun_closest(Nw,w,sumw))
             ELSE
                RadVelInterp(:,2)=CMPLX(0.,0.)
             ENDIF
          ELSE
           DO I=1,3      !for Vx,Vy,Vz
             IF (delw.GT.0.) THEN
              RadVelInterp(I,1)=FUN_INTERP1_COMPLEX(w,RadVel(I,:),Nw,delw) !for delta omega
             ELSE
              RadVelInterp(I,1)=CMPLX(0.,0.) !for delta omega
             END IF
             IF (sumw.LE.w(Nw)) THEN
              RadVelInterp(I,2)=FUN_INTERP1_COMPLEX(w,RadVel(I,:),Nw,sumw) !for sum omega
             ELSE
                RadVelInterp(I,2)=CMPLX(0.,0.)
             ENDIF
           ENDDO
          ENDIF
  END FUNCTION


  FUNCTION DOT_PRODUCT_DIFF_BIHARM(var11,var12,var21,var22,Nvect) RESULT(prod)
           !Dot product for the difference frequency term of biharmonic function
           INTEGER,                  INTENT(IN):: NVect
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var21,var22
           INTEGER                             :: I
           COMPLEX                             :: prod
           prod=0.5*(DOT_PRODUCT_COMPLEX(      &
                    COMPLEX_CONJUGATE_VECT(var22,Nvect),var11,Nvect)  &
                   +DOT_PRODUCT_COMPLEX(        &
                   COMPLEX_CONJUGATE_VECT(var12,Nvect),var21,Nvect))
  END FUNCTION

  FUNCTION DOT_PRODUCT_SUM_BIHARM(var11,var12,var21,var22,Nvect) RESULT(prod)
           !Dot product for the sum frequency term of biharmonic function
           INTEGER,                  INTENT(IN):: NVect
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var21,var22
           COMPLEX                             :: prod
           prod=0.5*(DOT_PRODUCT_COMPLEX(var22,var11,Nvect)            &
                        +DOT_PRODUCT_COMPLEX(var12,var21,Nvect))
  END FUNCTION

  FUNCTION CROSS_PRODUCT_DIFF_BIHARM(var11,var12,var21,var22) RESULT(prod)
           !CROSS product for the difference frequency term of biharmonic function
           COMPLEX,DIMENSION(3), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(3), INTENT(IN):: var21,var22
           COMPLEX,DIMENSION(3)            :: prod
           prod=0.5*(CROSS_PRODUCT_COMPLEX(var11,COMPLEX_CONJUGATE_VECT(var22,3)) &
                     +CROSS_PRODUCT_COMPLEX(COMPLEX_CONJUGATE_VECT(var12,3),var21))
  END FUNCTION

  FUNCTION CROSS_PRODUCT_SUM_BIHARM(var11,var12,var21,var22) RESULT(prod)
           !CROSS product for the sum frequency term of biharmonic function
           COMPLEX,DIMENSION(3), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(3), INTENT(IN):: var21,var22
           COMPLEX,DIMENSION(3)            :: prod
           prod=0.5*(CROSS_PRODUCT_COMPLEX(var11,var22)                            &
                        +CROSS_PRODUCT_COMPLEX(var12,var21))
  END FUNCTION

  FUNCTION Function_SymIpanel(NpanWlin,Npanels,Isym) RESULT(Ipanelinit)
        INTEGER,            INTENT(IN)  ::NpanWlin,Npanels,Isym
        INTEGER,DIMENSION(2)            ::Ipanelinit
           IF (Isym==1) Ipanelinit(1:2)=0
           IF (Isym==2) THEN
                   Ipanelinit(1)=NpanWLin
                   Ipanelinit(2)=Npanels
           ENDIF
  END FUNCTION

  FUNCTION Function_SymIwline(NpanWlin,Npanels,Nwline,Isym) RESULT(Iwlineinit)
        INTEGER,            INTENT(IN)  ::NpanWlin,Npanels,Nwline,Isym
        INTEGER,DIMENSION(2)            ::Iwlineinit
           IF (Isym==1) THEN
                   Iwlineinit(1)=Npanels
                   Iwlineinit(2)=0
           ELSEIF (Isym==2) THEN
                   Iwlineinit(1)=Npanels+NpanWlin
                   Iwlineinit(2)=Nwline
           ENDIF
  END FUNCTION


END MODULE
