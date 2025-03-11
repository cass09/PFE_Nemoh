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
!   - R. Kurnia
!
!--------------------------------------------------------------------------------------
MODULE MPP_ReadInputFiles

IMPLICIT NONE

TYPE TMech
    REAL,ALLOCATABLE,DIMENSION(:,:) :: MassMat          !Mass-Inertia Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat         !Stifness-Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat_EXT     !Additional Stifness Matrix i.e: mooring
    REAL,ALLOCATABLE,DIMENSION(:,:) :: DampCoefMat_EXT  !Additional damping coefficients
END TYPE TMech
CONTAINS
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
