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
!   - R. Kurnia (2022)
!
!--------------------------------------------------------------------------------------
MODULE MPP_Compute_RAOs

    USE CONSTANTS,          ONLY:II,PI
    USE MResults
    USE MPP_ReadInputFiles, ONLY:TMech
    USE M_SOLVER,           ONLY:LU_INVERS_MATRIX
    USE MNemohCal,          ONLY:TNemCal
!
    IMPLICIT NONE


CONTAINS
    SUBROUTINE Compute_RAOs(RAOS,Results,MechCoef)
!
!
!   Inputs/outputs
    TYPE(TResults),             INTENT(IN) :: Results
    TYPE(TMech),                INTENT(IN) :: MechCoef
    COMPLEX,DIMENSION(Results%Nintegration,Results%Nw,Results%Nbeta) :: RAOs
!   Locals
    INTEGER :: Iw,Ibeta
    REAL    :: w
    COMPLEX,DIMENSION(Results%Nradiation,Results%Nintegration):: MAT_A,invMAT_A
    COMPLEX,DIMENSION(Results%Nintegration)                   :: ExcitForce
!
   DO Iw=1,Results%Nw
      w= Results%w(Iw)
      MAT_A=-(MechCoef%MassMat+Results%AddedMass(Iw,:,:))*w*w                   &
            -II*w*(Results%RadiationDamping(Iw,:,:)+MechCoef%DampCoefMat_EXT)   &
            +MechCoef%StiffMat+MechCoef%StiffMat_EXT
      CALL LU_INVERS_MATRIX(MAT_A,Results%Nradiation, invMAT_A)
      DO Ibeta=1,Results%Nbeta
       ExcitForce=Results%DiffractionForce(Iw,Ibeta,:)+Results%FroudeKrylovForce(Iw,Ibeta,:)
       RAOs(:,Iw,Ibeta)=MATMUL(invMAT_A,ExcitForce)
      ENDDO
   ENDDO

    END SUBROUTINE Compute_RAOs

    SUBROUTINE SAVE_RAO(RAOs,w,beta,Nintegration,Nw,Nbeta,IndxForce,dirname,filename,InpNEMOHCAL)
    CHARACTER(LEN=*),                   INTENT(IN) :: filename,dirname
    TYPE(TNemCal),                      INTENT(IN) :: InpNEMOHCAL
    INTEGER,                            INTENT(IN) :: Nintegration,Nw, Nbeta
    REAL, DIMENSION(Nw),                INTENT(IN) :: w
    REAL, DIMENSION(Nbeta),             INTENT(IN) :: beta
    INTEGER,DIMENSION(Nintegration),    INTENT(IN) :: IndxForce
    COMPLEX,DIMENSION(Nintegration,Nw,Nbeta),&
                                        INTENT(IN) :: RAOs
    INTEGER :: Iw, Ibeta,Iinteg,u,Ninteg2
    REAL,DIMENSION(Nintegration*2) ::wRAOS
    REAL                           ::RAOSphase
    CHARACTER(LEN=23) :: FreqVar_text,FreqVar_text2
    REAL,DIMENSION(Nw)::freqVar

    IF (InpNEMOHCAL%OptOUTPUT%FreqType==1.OR.InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr==1) THEN
        FreqVar_text='VARIABLES="w (rad/s)"'
        freqVar=w
        FreqVar_text2='Number of pulsation= '
    ELSEIF (InpNEMOHCAL%OptOUTPUT%FreqType==2) THEN
        FreqVar_text='VARIABLES="f (Hz)"'
        freqVar=w/2/PI
        FreqVar_text2='Number of frequency= '
    ELSEIF (InpNEMOHCAL%OptOUTPUT%FreqType==3) THEN
        FreqVar_text='VARIABLES="T (s)"'
        freqVar=2*PI/w
        FreqVar_text2='Number of periode= '
    ENDIF


    CALL make_directory(dirname)
    Ninteg2=2*Nintegration
    IF (Ninteg2.GT.10000) THEN
      PRINT*,'Increase the value used in the following write command'
      STOP
    ENDIF
    OPEN(NEWUNIT=u, FILE=TRIM(dirname)//TRIM(filename),ACTION='WRITE')
    ! TODO: fix names of variables to match number of columns
    WRITE(u,'(13(A,X))') FreqVar_text, '"|X| (m/m)"','"|Y| (m/m)"','"|Z| (m/m)"',  &
                '"|phi| (deg)"','"|theta| (deg)"','"|psi| (deg)"', &
                 '"ang(x) (deg)"','"ang(y) (deg)"','"ang(z) (deg)"', &
                 '"ang(phi) (deg)"','"ang(theta) (deg)"','"ang(psi) (deg)"'
    DO Ibeta=1,Nbeta
    WRITE(u,'((A,X),F10.3,(A,X),I4,(A,X))') 'ZONE T="RAO - beta=',beta(Ibeta)*180/PI,' (deg)", I=', Nw, ', F=POINT'
        DO Iw=1,Nw
             DO Iinteg=1,Nintegration
                wRAOS(Iinteg)=ABS(RAOs(Iinteg,Iw,Ibeta))
                IF (IndxForce(Iinteg).GT.3) wRAOS(Iinteg)=wRAOS(Iinteg)*180/PI
                RAOSphase=ATAN2(AIMAG(RAOs(Iinteg,Iw,Ibeta)),REAL(RAOs(Iinteg,Iw,Ibeta)))
                wRAOS(Iinteg+Nintegration)=180/PI*RAOSphase
             ENDDO
             WRITE(u,'((F10.3,X),10000(E14.7,X))') freqVar(Iw),(wRAOS(Iinteg),Iinteg=1,Ninteg2)
        ENDDO
    ENDDO
    CLOSE(u)

    END SUBROUTINE

   SUBROUTINE  make_directory(dirname)
          CHARACTER(LEN=*),       INTENT(IN) ::dirname
          LOGICAL                            ::existdir
          !INQUIRE (DIRECTORY=dirname, EXIST=existdir) !this is Intel-specific
          INQUIRE (FILE=TRIM(dirname)//'/.', EXIST=existdir)
          IF (.NOT.existdir) CALL SYSTEM('mkdir '//dirname)
   END SUBROUTINE

END MODULE
