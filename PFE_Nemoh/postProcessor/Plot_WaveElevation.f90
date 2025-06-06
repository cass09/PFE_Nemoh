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
    SUBROUTINE Plot_WaveElevation(ID,Environment,iw,iBeta,RAOS,Results)
!
    USE Constants, only: PI, II
    USE MIdentification
    USE MResults
    USE MEnvironment
    USE Elementary_functions, ONLY: CIH
!
    IMPLICIT NONE
!
!   Inputs/outputs
    TYPE(TID) :: ID
    INTEGER :: iw,iBeta
    TYPE(TResults) :: Results
    COMPLEX,DIMENSION(Results%Nintegration,Results%Nw,*) :: RAOs
    TYPE(TEnvironment) :: Environment
!   Locals
    CHARACTER(LEN=20) :: lookfor
    INTEGER :: Nx,Ny
    REAL :: Lx,Ly
    REAL,DIMENSION(:),ALLOCATABLE :: X,Y
    REAL :: r,theta
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: etaI,etaP,eta
    INTEGER :: j,i,k,l
    REAL :: w,kwave
    COMPLEX :: HKleft,HKright,HKochin,Potential,p,Vx,Vy,Vz
!
!   Read data
    OPEN(10,FILE=TRIM(ID%ID)//'/Nemoh.cal')
    READ(10,'(A20)') lookfor
    DO WHILE (lookfor.NE.'--- Post processing ')
        READ(10,'(A20)') lookfor
    END DO
    DO i=1,3
        READ(10,*)
    END DO
    READ(10,*) Nx,Ny,Lx,Ly
    CLOSE(10)
!   Calculate and save wave elevations
    ALLOCATE(X(Nx),Y(Ny),etaI(Nx,Ny),etaP(Nx,Ny),eta(Nx,Ny))
    DO i=1,Nx
        X(i)=-0.5*Lx+Lx*(i-1)/(Nx-1)
    END DO
    DO i=1,Ny
        Y(i)=-0.5*Ly+Ly*(i-1)/(Ny-1)
    END DO
    w=Results%w(iw)
    kwave=Wavenumber(w,Environment)
    DO i=1,Nx
        DO j=1,Ny
            r=SQRT((X(i)-Environment%XEFF)**2+(Y(j)-Environment%YEFF)**2)
            theta=ATAN2((Y(j)-Environment%YEFF),(X(i)-Environment%XEFF))
            k=1
            DO WHILE ((k.LT.Results%Ntheta-1).AND.(Results%theta(k+1).LT.theta))
                k=k+1
            END DO
            IF (k.EQ.Results%Ntheta) THEN
                WRITE(*,*) ' Error: range of theta in Kochin coefficients is too small'
                STOP
            END IF
            CALL Compute_Wave(kwave,w,Results%beta(iBeta),X(i),Y(j),0.,Potential,p,Vx,Vy,Vz,Environment)
            EtaI(i,j)=1./Environment%G*II*w*Potential
            HKleft=0.
            HKright=0.
            ! CML w added
            DO l=1,Results%Nradiation
                HKleft=HKleft+w*RAOs(l,iw,iBeta)*Results%HKochinRadiation(iw,l,k)
                HKright=HKright+w*RAOs(l,iw,iBeta)*Results%HKochinRadiation(iw,l,k+1)
            END DO
            HKleft=-II*(HKleft+Results%HKochinDiffraction(iw,iBeta,k))
            HKright=-II*(HKright+Results%HKochinDiffraction(iw,iBeta,k+1))
            HKochin=HKleft+(HKright-HKleft)*(theta-Results%theta(k))/(Results%theta(k+1)-Results%theta(k))
            IF (r.GT.0) THEN
                ! Potential=SQRT(kwave/(2.*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
                Potential=SQRT(2./(kwave*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
            ELSE
                Potential=0.
            END IF
            EtaP(i,j)=1./Environment%G*II*w*Potential
            Eta(i,j)=EtaI(i,j)+EtaP(i,j)
        END DO
    END DO
    OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField_I.dat')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image"'
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(etaI(i,j)), ATAN2(AIMAG(etaI(i,j)), REAL(etaI(i,j))),REAL(etaI(i,j)),AIMAG(etaI(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField_S.dat')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image"'
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(etaP(i,j)), ATAN2(AIMAG(etaP(i,j)), REAL(etaP(i,j))),REAL(etaP(i,j)),AIMAG(etaP(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField.dat')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image" '
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(eta(i,j)), ATAN2(AIMAG(eta(i,j)), REAL(eta(i,j))),REAL(eta(i,j)),AIMAG(eta(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
        END DO
    END DO
    CLOSE(10)
    DEALLOCATE(X,Y,etaI,etaP,eta)
!
    END SUBROUTINE Plot_WaveElevation

    SUBROUTINE Initialize_Plot_WaveElevation(Switch_Plot_WaveElevation,namefile)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: namefile
    CHARACTER(LEN=20) :: lookfor
    CHARACTER(LEN=80) :: discard
    INTEGER :: i
    REAL :: Switch_Plot_WaveElevation
    OPEN(10,FILE=namefile)
    READ(10,'(A20)') lookfor
    DO WHILE (lookfor.NE.'--- Post processing ')
        READ(10,'(A20,A)') lookfor,discard
    END DO
    DO i=1,3
        READ(10,*)
    END DO
     READ(10,*) Switch_Plot_WaveElevation
    CLOSE(10)
    END SUBROUTINE  Initialize_Plot_WaveElevation

    SUBROUTINE Calcul_Kochin_TOT(Results,RAOs, Kochin_tot, iw, iBeta, iBody)
    USE Constants, only: PI, II
    USE MIdentification
    USE MResults
    USE MEnvironment
    USE Elementary_functions, ONLY: CIH
!
    IMPLICIT NONE
!
!   Inputs/outputs
    TYPE(TResults) :: Results
    COMPLEX,DIMENSION(Results%Nintegration,Results%Nw,Results%Nbeta) :: RAOs
    COMPLEX,DIMENSION(Results%Nw,Results%Nbeta,Results%Ntheta) :: Kochin_tot
    INTEGER :: iw, iBeta, iBody
!   Locals
    
    INTEGER :: l, itheta
    COMPLEX :: Krad

    DO itheta=1, Results%Ntheta
        Krad=0
        DO l=1+(iBody-1)*6, iBody*6
            Krad=Krad+RAOs(l,iw,iBeta)*Results%HKochinRadiation(iw,l,itheta)
        END DO
        Kochin_tot(iw,iBeta,itheta)=-II*(Results%HKochinDiffraction(iw,iBeta,itheta)+Krad*Results%w(iw))
    END DO
!
    END SUBROUTINE  Calcul_Kochin_TOT
!

    SUBROUTINE Plot_ETA(ID,Environment,RAOs, Results)
!
    USE Constants, only: PI, II
    USE MIdentification
    USE MResults
    USE MEnvironment
    USE Elementary_functions, ONLY: CIH
!
    IMPLICIT NONE
!
!   Inputs/outputs
    TYPE(TID) :: ID
    TYPE(TResults) :: Results
    COMPLEX,DIMENSION(Results%Nintegration,Results%Nw,Results%Nbeta) :: RAOs
    COMPLEX,DIMENSION(Results%Nw,Results%Nbeta, Results%Ntheta) :: Kochin, Kochin_tot
    TYPE(TEnvironment) :: Environment
!   Locals
    CHARACTER(LEN=20) :: lookfor
    INTEGER :: Nx,Ny
    REAL :: Lx,Ly
    REAL,DIMENSION(:),ALLOCATABLE :: X,Y
    REAL,DIMENSION(:,:),ALLOCATABLE :: CoG
    REAL :: r,theta
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: etaI,etaP, etaD
    COMPLEX,DIMENSION(:,:,:,:),ALLOCATABLE :: etaR
    INTEGER :: j,i,k,l
    INTEGER :: iw,iBeta, iBody, Nbodies, Np, num_pb
    REAL :: dummy_icase, dummy_dir1, dummy_dir2, dummy_dir3
    REAL :: w,kwave
    COMPLEX :: HKleft,HKright,HKochin,Potential,p,Vx,Vy,Vz
!
!   Read data
    OPEN(10,FILE=TRIM(ID%ID)//'/Nemoh.cal')
    DO i=1,6
        READ(10,*)
    END DO
    READ(10,*) Nbodies
    ALLOCATE(CoG(Nbodies, 3))
    DO I = 1, Nbodies
        DO j=1,7
        READ(10,*)
        END DO
        READ(10,*) dummy_icase, dummy_dir1, dummy_dir2, dummy_dir3, &
                            CoG(I,1), CoG(I,2), CoG(I,3)
        DO j=1,10
        READ(10,*)
        END DO        ! Skip additional line
    END DO
    READ(10,'(A20)') lookfor
    DO WHILE (lookfor.NE.'--- Post processing ')
        READ(10,'(A20)') lookfor
    END DO
    DO i=1,3
        READ(10,*)
    END DO
    READ(10,*) Nx,Ny,Lx,Ly
    CLOSE(10)
!   Calculate and save wave elevations
    ALLOCATE(X(Nx),Y(Ny),etaI(Nx,Ny),etaP(Nx,Ny), &
    & etaD(Nx, Ny), etaR(Nx, Ny, 6, Nbodies))
    DO i=1,Nx
        X(i)=-0.5*Lx+Lx*(i-1)/(Nx-1)
    END DO
    DO i=1,Ny
        Y(i)=-0.5*Ly+Ly*(i-1)/(Ny-1)
    END DO
    ! Np=Results%Nw*(Results%Nbeta+Results%Nradiation)
    DO iBeta=1, Results%Nbeta
    DO iw=1, Results%Nw
        w=Results%w(iw)
        kwave=Wavenumber(w,Environment)
        EtaP=0
        EtaD=0
        DO iBody=1, Nbodies
        CALL Calcul_Kochin_TOT(Results,RAOs, Kochin_tot, iw, iBeta, iBody)
        DO i=1,Nx
            DO j=1,Ny
                r=SQRT((X(i)-CoG(iBody, 1))**2+(Y(j)-CoG(iBody,2))**2)
                theta=ATAN2((Y(j)-CoG(iBody, 2)),(X(i)-CoG(iBody, 1)))
                k=1
                ! print*, theta
                if (theta < 0) then
                    theta = theta+2 * PI
                end if
                DO WHILE ((k.LT.Results%Ntheta-1).AND.(Results%theta(k+1).LT.theta))
                    k=k+1
                END DO
                IF (k.EQ.Results%Ntheta) THEN
                    WRITE(*,*) ' Error: range of theta in Kochin coefficients is too small'
                    STOP
                END IF
                CALL Compute_Wave(kwave,w,Results%beta(iBeta),X(i),Y(j),0.,Potential,p,Vx,Vy,Vz,Environment)
                EtaI(i,j)=1./Environment%G*II*w*Potential

                HKleft=Kochin_tot(iw, iBeta, k)
                HKright=Kochin_tot(iw, iBeta, k+1)
                HKochin=HKleft+(HKright-HKleft)*(theta-Results%theta(k))/(Results%theta(k+1)-Results%theta(k))
                IF (r.GT.0) THEN
                    ! Potential=SQRT(kwave/(2.*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
                    Potential=SQRT(2./(kwave*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
                ELSE
                    Potential=0.
                END IF
                EtaP(i,j)=EtaP(i,j)+1./Environment%G*II*w*Potential
                ! Eta(i,j)=EtaI(i,j)+EtaP(i,j)

                HKleft=-II*Results%HKochinDiffraction(iw, iBeta, k)
                HKright=-II*Results%HKochinDiffraction(iw, iBeta, k+1)
                HKochin=HKleft+(HKright-HKleft)*(theta-Results%theta(k))/(Results%theta(k+1)-Results%theta(k))
                IF (r.GT.0) THEN
                    ! Potential=SQRT(kwave/(2.*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
                    Potential=SQRT(2./(kwave*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
                ELSE
                    Potential=0.
                END IF
                EtaD(i,j)=EtaD(i,j)+1./Environment%G*II*w*Potential

                DO l=1, Results%Nradiation/Nbodies
                    HKleft=-II*w*Results%HKochinRadiation(iw, l, k)*RAOs(l, iw, iBeta)
                    HKright=-II*w*Results%HKochinRadiation(iw, l, k+1)*RAOs(l, iw, iBeta)
                    HKochin=HKleft+(HKright-HKleft)*(theta-Results%theta(k))/(Results%theta(k+1)-Results%theta(k))
                    IF (r.GT.0) THEN
                        ! Potential=SQRT(kwave/(2.*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
                        Potential=SQRT(2./(kwave*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
                    ELSE
                        Potential=0.
                    END IF
                    EtaR(i,j, l, iBody)=1./Environment%G*II*w*Potential
                END DO
            END DO
        END DO
        END DO
        CALL WRITE_ETA(ID, Nx, Ny, X, Y, etaI, etaP, etaD, etaR, iw, Results%w(iw), iBeta, Nbodies)
    END DO
    END DO

    
    DEALLOCATE(etaI, etaP, X, Y, CoG, etaD, etaR)
    END SUBROUTINE Plot_ETA

    SUBROUTINE WRITE_ETA(ID, Nx, Ny, X, Y, etaI, etaP, etaD, etaR, iw, w, iBeta, Nbodies)

    USE Constants, only: PI, II
    USE MIdentification
    USE MResults
    USE MEnvironment
    USE Elementary_functions, ONLY: CIH
    
!
    IMPLICIT NONE
!

    TYPE(TID) :: ID
    INTEGER :: Nx, Ny, iw, iBeta, Nbodies
    REAL :: w
    COMPLEX,DIMENSION(Nx,Ny) :: etaI,etaP,eta, etaD
    COMPLEX,DIMENSION(Nx,Ny,6, Nbodies) :: etaR
    REAL,DIMENSION(Nx) :: Y
    REAL,DIMENSION(Ny) :: X
    INTEGER :: i, j, k, num_pb, iBody
    CHARACTER(LEN=10) :: strPB, strW, strBeta
    WRITE(strW, '(I0)') iw  ! 'I0' = format sans espace inutile
    WRITE(strBeta, '(I0)') iBeta  ! 'I0' = format sans espace inutile
    eta=etaP+etaI

    OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField_I.dat')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image"'
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(etaI(i,j)), ATAN2(AIMAG(etaI(i,j)), REAL(etaI(i,j))),REAL(etaI(i,j)),AIMAG(etaI(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
        END DO
    END DO
    CLOSE(10)

    OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField_S.dat')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image"'
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(etaP(i,j)), ATAN2(AIMAG(etaP(i,j)), REAL(etaP(i,j))),REAL(etaP(i,j)),AIMAG(etaP(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
        END DO
    END DO
    CLOSE(10)

    OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField.dat')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image" '
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(eta(i,j)), ATAN2(AIMAG(eta(i,j)), REAL(eta(i,j))),REAL(eta(i,j)),AIMAG(eta(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
        END DO
    END DO
    CLOSE(10)
    num_pb=(iw-1)*(1+6*Nbodies)+1
    WRITE(strPb, '(I0.5)') num_pb
    OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField_'//TRIM(strPb)//'.dat')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image" '
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(etaD(i,j)), ATAN2(AIMAG(etaD(i,j)), REAL(etaD(i,j))),REAL(etaD(i,j)),AIMAG(etaD(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
        END DO
    END DO
    CLOSE(10)

    DO iBody=1, Nbodies
    DO k=1, 6
        num_pb=num_pb+1
        WRITE(strPb, '(I0.5)') num_pb
        OPEN(10,FILE=TRIM(ID%ID)//'/results/WaveField_'//TRIM(strPb)//'.dat')
        WRITE(10,'(A)') 'VARIABLES="X" "Y" "abs" "phase" "real" "image" '
        WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
        DO i=1,Nx
            DO j=1,Ny
                WRITE(10,'(10(X,E14.7))') X(i),Y(j),ABS(etaR(i,j,k, iBody)), ATAN2(AIMAG(etaR(i,j,k, iBody)), REAL(etaR(i,j,k,iBody))),REAL(etaR(i,j,k,iBody)),AIMAG(etaR(i,j,k,iBody))
            END DO
        END DO
        DO i=1,Nx-1
            DO j=1,Ny-1
                WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
            END DO
        END DO
        CLOSE(10)
    END DO
    END DO

    END SUBROUTINE