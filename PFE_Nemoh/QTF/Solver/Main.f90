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
!--------------------------------------------------------------------------------------
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!     THESE DE CHEN XIAO-BO(1988)
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)  Version 2014
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!   SOLVER
!
!--------------------------------------------------------------------------------------
PROGRAM MAIN
!
USE MIdentification
USE MNemohCal,          ONLY:TNemCal,READ_TNEMOHCAL,Tqtfinput,                &
                             IntegrationAXIS_FROM_MNEMOHCAL
USE MMesh
USE MFace,              ONLY:TVFace, Prepare_FaceMesh,TWLine,Prepare_Waterline
USE MReadInputFiles,    ONLY:Read_NP_GaussQuad,Read_Mechanical_Coefs,TMech,   &
                             Read_FirstOrderLoad,TLoad1,Read_Motion,TSource,  &
                             READ_POTENTIALS_VELOCITIES,TpotVel,              &
                             READ_GENERALIZED_NORMAL_BODY_dAREA,Read_Eps_Zmin,&
                             TMeshFS,Read_Prepare_FreeSurface_Mesh
USE MEnvironment,       ONLY: TEnvironment,FunVect_inverseDispersion
USE MLogFile
USE Constants,          ONLY: CZERO
USE MQSolverPreparation !CONTAINS:TQfreq,PREPARE_POTENTIAL_VELOCITIES
                        !PREPARE_BODY_DISPLACEMENT,DISCRETIZED_OMEGA_WAVENUMBER_FOR_QTF
                        !PREPARE_INERTIA_FORCE
                        !CALC_GENERALIZED_NORMAL_WATERLINE_dSEGMENT
                        !WRITE_QTFSOLVERLOGFILE
USE MQSolver            !CONTAINS:COMPUTATION_QTF_QUADRATIC,...
USE MQSolverOutputFiles !CONTAINS: OutFileDM,OutFileDP,...
                        !INITIALIZE_OUTPUT_FILE,WRITE_QTF_DATA

IMPLICIT NONE
!
INTEGER,parameter :: ID_DEBUG=0 ! for debugging, each QTFs terms will be saved
!
! ------Declaration variables
        TYPE(TID)                       :: ID
        TYPE(TMesh)                     :: Mesh
        TYPE(TMeshFS)                   :: MeshFS
        TYPE(TNemCal)                   :: inpNEMOHCAL
        TYPE(Tqtfinput)                 :: QTFinputNem
        TYPE(TEnvironment)              :: Env
        INTEGER                         :: Nw,Nbeta,Nbeta2    ! Number of Freq,direction
        REAL, ALLOCATABLE,DIMENSION(:)  :: w,kw,beta          ! vector of freq [rad/s],
                                                              ! wave numbers [rad/m],
                                                              ! direction angle [rad]
        TYPE(TWLine)  :: WLine          ! Waterline
        TYPE(TVFace)  :: VFace          ! Face of body panel
        TYPE(TVFace)  :: VFaceFS        ! Face of free surface panel
        TYPE(TMech)   :: MechCoef       ! Mechanical Coef
        TYPE(TLoad1)  :: Forces1        ! First order forces
        TYPE(TASYMP)  :: ASYMP_PARAM    ! parameters for asymptotic free surface force
        TYPE(TSourceQ):: SOURCEDISTRQ   ! perturbion and radiation source distribution
        INTEGER       :: NP_GQ          ! Number of point for Gauss Quad. Integration
        INTEGER       :: Nintegration   ! Number of excitation force integration
        INTEGER       :: Nradiation     ! Number of radiation problem
        INTEGER       :: Nbodies        ! Number of bodies
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: Motion       ! Complex RAO
        TYPE(TPotVel)                           :: datPotVel    ! data Pot & vel
        TYPE(TPotVel)                           :: datPotVelQ    ! data Pot & vel for QTF
        TYPE(TPotVel)                           :: datPotVelFS   ! data Pot & vel (FreeSurface)
        TYPE(TPotVel)                           :: datPotVelQFS  ! data Pot & vel for QTF (FreeSurface)

        INTEGER                                 :: I,Ibeta1,Ibeta2,Iinteg,Ipanel,Ibeta2temp
        INTEGER                                 :: Iw1,Iw2,IwQ
        INTEGER                                 :: NPFlow        ! Number of flow points
        INTEGER                                 :: NPFlowFS      ! Number of flow points (Free-surface)
        REAL,ALLOCATABLE,DIMENSION(:,:)         :: genNormal_dS  ! generalized Normal on panel
                                                           ! time the area of the panel
        REAL,ALLOCATABLE,DIMENSION(:,:)         :: genNormalWLine_dGamma! generalized Normal
                                                   !on wLine segment time the segm. length
        REAL,ALLOCATABLE,DIMENSION(:,:)         :: IntegAxis,StiffMat
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:)  :: BDisplaceQ!Body displacement
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: InertiaForceQ
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: RotAnglesQ,TransMotionQ
        !freq related variables For QTF computation
        TYPE(TQfreq)                            :: Qfreq
        REAL                                    :: winputQ(3), BForwardSpeed,delwiter
        INTEGER                                 :: NwQ    ! Number of wave freq
        !---------
        INTEGER                                 :: SwitchQuadHM,SwitchBiDir
        REAL                                    :: EPS_ZMIN
        COMPLEX, ALLOCATABLE,DIMENSION(:,:,:)   :: QTF_DUOK ,QTF_HASBO
        COMPLEX, ALLOCATABLE,DIMENSION(:,:,:)   :: QTF_HASFS,QTF_HASFS_ASYMP
        INTEGER                                 :: Iterm,ufile
        CHARACTER(LEN=1)                        :: strI
        !
        CHARACTER(LEN=1000)                :: LogTextToBeWritten
        REAL                               :: tcpu_start


!
!   --- Initialize and read input datas -----------------------------------------------
!
        WRITE(*,*) 'QTF Solver preparation...'
        CALL ReadTID(ID)
        CALL ReadTMesh(Mesh,TRIM(ID%ID)//'/mesh/')
        CALL READ_TNEMOHCAL(ID,InpNEMOHCAL)
        !
        QTFinputNem  =InpNEMOHCAL%qtfinput
        Env          =InpNEMOHCAL%Env
        Nw           =InpNEMOHCAL%waveinput%NFreq
        Nbeta        =InpNEMOHCAL%waveinput%NBeta
        Nbodies      =InpNEMOHCAL%Nbodies
        Nintegration =InpNEMOHCAL%Nintegtot
        Nradiation   =InpNEMOHCAL%Nradtot
        winputQ      =InpNEMOHCAL%qtfinput%omega
        SwitchQuadHM =InpNEMOHCAL%qtfinput%switch_quadHM
        SwitchBiDir  =InpNEMOHCAL%qtfinput%bidirection
        NwQ          =winputQ(1)
        BForwardSpeed=InpNEMOHCAL%qtfinput%body_forward_speed
        NP_GQ        =Read_NP_GaussQuad(TRIM(ID%ID))
        EPS_ZMIN     =Read_Eps_Zmin(TRIM(ID%ID))
        !
        CALL Prepare_FaceMesh(Mesh,NP_GQ,VFace)
        CALL Prepare_Waterline(VFace,EPS_ZMIN,Mesh%xy_diameter,Mesh%Npanels,WLine)
        !
        NPFlow   =(Mesh%Npanels+WLine%NWLineSeg)*2**Mesh%Isym !Number of flow point
        !
!       -------------------------------------
!       Prepare Free-surface mesh
        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
          CALL Read_Prepare_FreeSurface_Mesh(TRIM(ID%ID),MeshFS,QTFinputNem)
          CALL Prepare_FaceMesh(MeshFS%Mesh,NP_GQ,VFaceFS)
          NPFlowFS   =(MeshFS%Mesh%Npanels+MeshFS%BdyLine%NWLineSeg)*2**MeshFS%Mesh%Isym !Number of flow point
        ENDIF
!       --------------------------------------

        !Dynamic Memory allocation
        ALLOCATE(IntegAxis(3,Nintegration))
        ALLOCATE(StiffMat(Nintegration,Nintegration))
        ALLOCATE(Motion(Nw,Nradiation,Nbeta))

        ALLOCATE(w(Nw),kw(Nw),beta(Nbeta))
        ALLOCATE(genNormal_dS(Nintegration,Mesh%Npanels*2**Mesh%Isym))
        ALLOCATE(genNormalWLine_dGamma(Nintegration,Wline%NWLineseg*2**Mesh%Isym))
        ALLOCATE(InertiaForceQ(NwQ,Nbeta,Nintegration))
        ALLOCATE(RotAnglesQ(NwQ,Nbeta,3*Nbodies))
        ALLOCATE(TransMotionQ(NwQ,Nbeta,Nradiation))
        ALLOCATE(QTF_DUOK(Nintegration,2,7))!2 is for QTF- and QTF+, 7 for all the terms
        ALLOCATE(QTF_HASBO(Nintegration,2,7))

        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
        ALLOCATE(QTF_HASFS(Nintegration,2,10))!2 is for QTF- and QTF+, 7 for all the terms
        ALLOCATE(QTF_HASFS_ASYMP(Nintegration,2,3))
        ENDIF

        !
        IntegAxis=IntegrationAXIS_FROM_MNEMOHCAL(InpNEMOHCAL)
        !
        CALL Read_Mechanical_Coefs(TRIM(ID%ID),Nradiation,MechCoef)
        StiffMat=MechCoef%StiffMat+MechCoef%StiffMat_EXT
        !
        CALL Read_FirstOrderLoad(TRIM(ID%ID),Nw,Nbeta,Nintegration,Nradiation,Forces1)
        CALL Read_Motion(TRIM(ID%ID),Nw,Nbeta,Nradiation,Motion)!Complex RAO
        CALL READ_POTENTIALS_VELOCITIES(TRIM(ID%ID),Nw,Nbeta,NRadiation,                 &
                NPFlow,datPotVel,w,kw,beta,ID_BODY)

        CALL READ_GENERALIZED_NORMAL_BODY_dAREA(TRIM(ID%ID),Mesh%Npanels*2**Mesh%Isym,   &
                                                          Nintegration,genNormal_dS)
        CALL CALC_GENERALIZED_NORMAL_WATERLINE_dSEGMENT(Mesh,Nintegration,WLine,         &
                                                InpNEMOHCAL, genNormalWLine_dGamma)
        !READ Free-surface potentials and velocities
        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
          CALL READ_POTENTIALS_VELOCITIES(TRIM(ID%ID),Nw,Nbeta,NRadiation,               &
                NPFlowFS,datPotVelFS,w,kw,beta,ID_FREESURFACE)
        ENDIF

        CALL Discretized_omega_wavenumber_for_QTF(Nw,w,kw,NwQ,winputQ(2:3),Nbeta,beta,   &
                                                  BForwardSpeed,Env%depth,Env%g,Qfreq)


        CALL PREPARE_POTENTIAL_VELOCITIES(Qfreq,Nw,w,Nbeta,beta,NPFlow,                  &
                Nradiation,datPotVel,datPotVelQ,ID_BODY)

        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
         WRITE(*,*) 'QTF Solver preparation, Free-Surface mesh...'
         CALL PREPARE_POTENTIAL_VELOCITIES(Qfreq,Nw,w,Nbeta,beta,NPFlowFS,               &
                Nradiation,datPotVelFS,datPotVelQFS,ID_FREESURFACE)
         CALL PREPARE_ASYMP_PARAM(MeshFS%Radius_Ext,MeshFS%NpointsR,                     &
                                                        MeshFS%NBessel,ASYMP_PARAM)
         CALL PREPARE_SOURCE_DISTRIBUTION(TRIM(ID%ID),Qfreq,Nw,w,Nbeta,beta,             &
                 Mesh%Npanels,Nradiation,Motion,SOURCEDISTRQ)
        ENDIF

        CALL PREPARE_BODY_DISPLACEMENT(Qfreq,Nw,w,Nbeta,Nradiation,NPFlow,Nbodies,       &
                                        Mesh,WLine,InpNEMOHCAL,Motion,BdisplaceQ)

        CALL PREPARE_INERTIA_FORCES(MechCoef,Motion,Nw,Nbeta,Nradiation,Nintegration,    &
                                        w,Qfreq,Forces1,InertiaForceQ)
        CALL PREPARE_ROTATION_ANGLES(Motion,Nw,Nbeta,Nradiation,Nbodies,                 &
                                        w,Qfreq,RotAnglesQ)
        CALL PREPARE_TRANSLATION_MOTION(Motion,Nw,Nbeta,Nradiation,Nbodies,              &
                                        w,Qfreq,TransMotionQ)

       CALL  INITIALIZE_OUTPUT_FILES(TRIM(ID%ID),InpNEMOHCAL%qtfinput%Ncontrib,ID_DEBUG)

       WRITE(*,*) 'QTF Solver preparation, Done!'
       CALL WRITE_QTFSOLVERLOGFILE(TRIM(ID%ID),Nbeta,beta,Qfreq)
       CALL START_RECORD_TIME(tcpu_start,TRIM(ID%ID)//'/'//LogFILE,IdAppend)

        Nbeta2=Nbeta
        IF (SwitchBiDir==0) Nbeta2=1

        DO Ibeta1=1,Nbeta
           DO Ibeta2temp=1,Nbeta2
                IF (SwitchBiDir==0) Ibeta2=Ibeta1
                IF (SwitchBiDir==1) Ibeta2=Ibeta2temp

                WRITE(*,'(A,F8.3,A,F8.3,A)'),'beta1=', beta(Ibeta1)*180/PI,&
                        ', beta2=', beta(Ibeta2)*180/PI, ' [deg]'
                DO IwQ=0,NwQ-1
                    DO Iw1=IwQ+1,NwQ
                        Iw2=Iw1-IwQ
                        IF (Iw1==1 .AND. Iw2==1) THEN
                                delwiter=Qfreq%wQ(Iw1,Ibeta1)-Qfreq%wQ(Iw2,Ibeta2)
                                WRITE(*,'(A,F7.3,A,F7.3,A)')'w1-w2=',delwiter, ', w1+w2=', &
                                     Qfreq%wQ(Iw2,Ibeta2)+Qfreq%wQ(Iw1,Ibeta1), ' [rad/s]'
                        ENDIF

                        IF (Qfreq%wQ(Iw1,Ibeta1)-Qfreq%wQ(Iw2,Ibeta2).GT.1.01*delwiter) THEN
                          delwiter=Qfreq%wQ(Iw1,Ibeta1)-Qfreq%wQ(Iw2,Ibeta2)
                          IF (Qfreq%wQ(Iw2,Ibeta2)+Qfreq%wQ(Iw1,Ibeta1).LE.w(Nw)) THEN
                             WRITE(*,'(A,F7.3,A,F7.3,A)'),'w1-w2=',delwiter, ', w1+w2=', &
                                     Qfreq%wQ(Iw2,Ibeta2)+Qfreq%wQ(Iw1,Ibeta1), ' [rad/s]'
                          ELSE
                            WRITE(*,'(A,F7.3,A)'),'w1-w2=',delwiter, ', w1+w2= --NA-- [rad/s]'
                          ENDIF
                        ENDIF

                        IF (InpNEMOHCAL%qtfinput%Ncontrib.GE.1) THEN
                        CALL COMPUTATION_QTF_QUADRATIC(Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,&
                                Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,Env%rho,Env%g,  &
                                datPotVelQ,BdisplaceQ,genNormal_dS,genNormalWLine_dGamma, &
                                Qfreq%wQ,beta,InertiaForceQ,RotAnglesQ,IntegAxis,         &
                                StiffMat,TransMotionQ,SwitchQuadHM,QTF_DUOK)

                        CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileDM,OutFileDP,Nintegration, &
                                Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),                &
                                beta(Ibeta1),beta(Ibeta2), QTF_DUOK(:,:,7))
                         IF (ID_DEBUG==1) THEN
                           DO Iterm=1,6
                           WRITE(strI,'(I0.1)') Iterm
                           CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileDM_term//strI//'.dat',     &
                                   OutFileDP_term//strI//'.dat',Nintegration,                &
                                   Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2), beta(Ibeta1),  &
                                   beta(Ibeta2),QTF_DUOK(:,:,Iterm))
                           ENDDO
                         ENDIF
                        ENDIF

                        IF (InpNEMOHCAL%qtfinput%Ncontrib.GE.2) THEN
                        CALL COMPUTATION_QTF_POTENTIAL_BODYFORCE(Iw1,Iw2,Ibeta1,Ibeta2,   &
                                Nintegration,Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,   &
                                Env, datPotVelQ,BdisplaceQ,genNormal_dS,                  &
                                Nw, w,Qfreq,beta,RotAnglesQ,QTF_HASBO)

                        CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileHBM,OutFileHBP,Nintegration,&
                                Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),                &
                                beta(Ibeta1),beta(Ibeta2), QTF_HASBO(:,:,7))

                         IF (ID_DEBUG==1) THEN
                          DO Iterm=1,6
                          WRITE(strI,'(I0.1)') Iterm
                          CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileHBM_term//strI//'.dat',    &
                                  OutFileHBP_term//strI//'.dat',Nintegration,               &
                                  Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),                &
                                  beta(Ibeta1),beta(Ibeta2),QTF_HASBO(:,:,Iterm))
                          ENDDO
                         ENDIF
                        ENDIF

                        IF (InpNEMOHCAL%qtfinput%Ncontrib.EQ.3) THEN
                        ! Finite domain Free surface force
                        CALL COMPUTATION_QTF_POTENTIAL_FREESURFACEFORCE(Iw1,Iw2,Ibeta1,   &
                                Ibeta2,Nintegration,MeshFS,NwQ,Nbeta,NPFlowFS,Nbodies,    &
                                Env,datPotVelQFS,Nw,w,Qfreq,beta,QTF_HASFS)

                        CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileHFSM,OutFileHFSP,          &
                                Nintegration,Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),   &
                                beta(Ibeta1),beta(Ibeta2), QTF_HASFS(:,:,10))

                         IF (ID_DEBUG==1) THEN
                            DO Iterm=1,9
                            WRITE(strI,'(I0.1)') Iterm
                            CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileHFSM_term//strI//'.dat',   &
                                    OutFileHFSP_term//strI//'.dat',Nintegration,              &
                                    Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),                &
                                    beta(Ibeta1),beta(Ibeta2),QTF_HASFS(:,:,Iterm))
                            ENDDO
                         ENDIF
                        ! InFinite domain (Asymptotic) Free surface force
                        CALL COMPUTATION_QTF_POTENTIAL_FREESURFACEFORCE_ASYMP(Iw1,Iw2,    &
                                Ibeta1,Ibeta2,Nintegration,NwQ,Nbeta,Env,Nw,w,Qfreq,      &
                                beta,Mesh,ASYMP_PARAM,SOURCEDISTRQ,QTF_HASFS_ASYMP)

                        CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileASYM,OutFileASYP,          &
                                Nintegration,Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),   &
                                beta(Ibeta1),beta(Ibeta2),QTF_HASFS_ASYMP(:,:,3))

                         IF (ID_DEBUG==1) THEN
                            DO Iterm=1,2
                            WRITE(strI,'(I0.1)') Iterm
                            CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileASYM_term//strI//'.dat',   &
                                    OutFileASYP_term//strI//'.dat',Nintegration,              &
                                    Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),                &
                                    beta(Ibeta1),beta(Ibeta2),QTF_HASFS_ASYMP(:,:,Iterm))
                            ENDDO
                         ENDIF
                        ENDIF
                    ENDDO
                ENDDO
           ENDDO
        ENDDO

        CALL END_RECORD_TIME(tcpu_start,TRIM(ID%ID)//'/'//LogFILE)
        WRITE(LogTextToBeWritten,*) '---- DONE ---'
        CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

! ----- Finalize ---------------------------------------------------------------------------
!       DEALOCATING variables
        DEALLOCATE(w,kw,beta)
        DO I=1,InpNEMOHCAL%Nbodies
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%RadCase)
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%IntCase)
        ENDDO
        DEALLOCATE(inpNEMOHCAL%bodyinput,IntegAxis)
        DEALLOCATE(Motion,StiffMat,TransMotionQ,RotAnglesQ)
        DEALLOCATE(VFace%X,VFace%XM,VFace%N,VFace%A,VFace%tDis)
        DEALLOCATE(VFace%dXdXG_WGQ_per_A,VFace%XM_GQ)
        DEALLOCATE(datPotVelQ%TotPot,datPotVelQ%TotVel)
        DEALLOCATE(datPotVelQ%RadPot,datPotVelQ%RadVel)
        DEALLOCATE(Qfreq%wQ,Qfreq%kQ)
        DEALLOCATE(Qfreq%diffwQ,Qfreq%sumwQ)
        DEALLOCATE(Qfreq%InterpPotSwitch)
        DEALLOCATE(genNormal_dS,genNormalWLine_dGamma)
        DEALLOCATE(BDisplaceQ)
        DEALLOCATE(InertiaForceQ)
        DEALLOCATE(QTF_DUOK,QTF_HASBO)
        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
          DEALLOCATE(datPotVelQFS%IncPot)
          DEALLOCATE(datPotVelQFS%IncVel)
          DEALLOCATE(datPotVelQFS%TotPot)
          DEALLOCATE(datPotVelQFS%TotVel)
          DEALLOCATE(datPotVelQFS%RadPot)
          DEALLOCATE(datPotVelQFS%RadVel)
          DEALLOCATE(QTF_HASFS)
          DEALLOCATE(QTF_HASFS_ASYMP)
          DEALLOCATE(SOURCEDISTRQ%ZIGB_Per)
          DEALLOCATE(SOURCEDISTRQ%ZIGS_Per)
          DEALLOCATE(SOURCEDISTRQ%ZIGB_Rad)
          DEALLOCATE(SOURCEDISTRQ%ZIGS_Rad)
        ENDIF

END PROGRAM
