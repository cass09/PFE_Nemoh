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
!   NEMOH V1.0 - postProcessor - January 2014
!
!--------------------------------------------------------------------------------------
!
    PROGRAM Main
!
    USE MIdentification
    USE MEnvironment
    USE MResults
    USE MIRF
    USE MPP_ReadInputFiles,     ONLY:Read_Mechanical_Coefs,TMech
    USE MNemohCal,              ONLY:TNemCal,READ_TNEMOHCAL
    USE MPP_Compute_RAOs
#ifndef GNUFORT
    USE iflport
#endif
!
    IMPLICIT NONE
!   ID
    TYPE(TID)               :: ID
!   NEMOHCAL
    TYPE(TNemCal)      :: inpNEMOHCAL
!   Environment
    TYPE(TEnvironment) :: Environment
!   Hydrodynamic coefficients cases
    TYPE(TResults) :: Results
!   IRFs
    TYPE(TIRF) :: IRF
!   Mechanical Coef: Mass_Mat,Stiffness,...
    TYPE(TMech)   :: MechCoef
!   RAOs
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: RAOs
!   Plot Wave elevation
    REAL :: Switch_Plot_WaveElevation

!   --- Initialisation -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' '
    WRITE(*,'(A,$)') '  -> Initialisation '
!   Read case ID
    CALL ReadTID(ID)
    WRITE(*,'(A,$)') '.'

!   Read Nemoh.call
    CALL READ_TNEMOHCAL(ID,InpNEMOHCAL)
!   Read environment
    Environment =InpNEMOHCAL%Env
!   Read results
    CALL ReadTResults(Results,TRIM(ID%ID)//'/results/Forces.dat',TRIM(ID%ID)//'/results/index.dat',TRIM(ID%ID)//'/results/FKForce.tec')
    CALL SaveTResults(Results,TRIM(ID%ID)//'/results',InpNEMOHCAL)
    WRITE(*,*) '. Done !'
    WRITE(*,*) ' '
!
!   --- Compute IRFs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL Initialize_IRF(IRF,Results,TRIM(ID%ID)//'/Nemoh.cal')
    IF (IRF%Switch.EQ.1) THEN
        CALL Compute_IRF(IRF,Results)
        CALL Save_IRF(IRF,TRIM(ID%ID)//'/results/IRF.tec')
        CALL Save_IRF_ExcForce(IRF,TRIM(ID%ID)//'/results/IRF_excForce.tec')
    END IF
!
!   --- Compute RAOs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!

    ALLOCATE(RAOs(Results%Nradiation,Results%Nw,Results%Nbeta))
    IF (InpNEMOHCAL%OptOUTPUT%Switch_RAO==1 .OR. InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr==1) THEN
      CALL Read_Mechanical_Coefs(TRIM(ID%ID),Results%Nradiation,MechCoef)
      CALL Compute_RAOs(RAOs,Results,MechCoef)
      CALL SAVE_RAO(RAOs,Results%w,Results%beta,Results%Nintegration,Results%Nw,Results%Nbeta,&
              Results%IndxForce(:,3),TRIM(ID%ID)//'/Motion/','RAO.dat',InpNEMOHCAL)
    ELSE
       RAOs(:,:,:)=CMPLX(0.,0.)
    ENDIF

!
!   --- Save results -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' -> Save results '
    WRITE(*,*) ' '

!    CALL Initialize_Plot_WaveElevation(Switch_Plot_WaveElevation,TRIM(ID%ID)//'/Nemoh.cal')
!    IF (Switch_Plot_WaveElevation.GT.1 ) THEN
!!       This function is not completely develop only produce incident wave elevation
!!       Kochin coefficients for diffraction and radiation is not yet post-processed
!        CALL Plot_WaveElevation(ID,Environment,1,1,RAOs,Results)
!    END IF

!
!   --- Finalize -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CALL DeleteTResults(Results)

    IF (InpNEMOHCAL%OptOUTPUT%Switch_RAO==1 .OR. InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr==1) THEN
      DEALLOCATE(RAOs)
    ENDIF
!
    END PROGRAM Main
!
