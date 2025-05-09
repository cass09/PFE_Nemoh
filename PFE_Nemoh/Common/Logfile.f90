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

MODULE MLogFile

  IMPLICIT NONE

  PUBLIC       ::  WRITE_LOGFILE,START_RECORD_TIME
  !  CPU TIME
  INTEGER,parameter, public  :: IdAppend=1      !Append writing to a file
  INTEGER,parameter, public  :: IdStartLog=0    !Starting a new file
  INTEGER,parameter, public  :: IdprintTerm=1   !print to the terminal
  INTEGER,parameter, public  :: IdNoprintTerm=0 !Not print to the terminal
CONTAINS

  SUBROUTINE WRITE_LOGFILE(logfile,textToBeWritten,FlagW,FlagWTerm)

        CHARACTER(LEN=*) logfile,textToBeWritten
        INTEGER u,FlagW,FlagWTerm
        IF (FlagW==IdStartLog) THEN
                OPEN(NEWUNIT=u, FILE=logfile, ACTION='WRITE')
        ELSE
                OPEN(NEWUNIT=u, FILE=logfile, ACTION='WRITE',POSITION='APPEND')
        ENDIF
        WRITE(u,*) textToBeWritten
        IF (FlagWTerm == IdprintTerm) THEN
        WRITE(*,*) textToBeWritten
        END IF
        CLOSE(u)
  END SUBROUTINE

  SUBROUTINE START_RECORD_TIME(tcpu_start,logfile,FlagW)
        CHARACTER(LEN=*) logfile
        INTEGER,DIMENSION(8)    :: DATETIMEVAL
        CHARACTER(LEN=1000)     :: textToBeWritten
        INTEGER                 :: FlagW
        REAL                    :: tcpu_start

        CALL CPU_TIME(tcpu_start)
        CALL DATE_AND_TIME(VALUES=DATETIMEVAL)

        CALL WRITE_LOGFILE(logfile,' STARTING TIME:',FlagW,IdprintTerm)
        WRITE(textToBeWritten,*)  ' Date: ',DATETIMEVAL(3),'-',DATETIMEVAL(2),'-',DATETIMEVAL(1)
        CALL WRITE_LOGFILE(logfile,TRIM(textToBeWritten),IDAppend,IdprintTerm)
        WRITE(textToBeWritten,*) ' Time: ',DATETIMEVAL(5),':',DATETIMEVAL(6),':',DATETIMEVAL(7)
        CALL WRITE_LOGFILE(logfile,TRIM(textToBeWritten),IDAppend,IdprintTerm)
  END SUBROUTINE

  SUBROUTINE END_RECORD_TIME(tcpu_start,logfile)
        CHARACTER(LEN=*) logfile
        INTEGER,DIMENSION(8)    :: DATETIMEVAL
        CHARACTER(LEN=1000)     :: textToBeWritten
        INTEGER                 :: FlagW
        REAL                    :: tcpu_start,tcpu_finish,comptime

        CALL CPU_TIME(tcpu_finish)
        comptime=tcpu_finish-tcpu_start
        WRITE(textToBeWritten,*) 'Computation time', comptime, ' [s]'
        CALL WRITE_LOGFILE(logfile,TRIM(textToBeWritten),IDAppend,IdprintTerm)
  END SUBROUTINE
END MODULE
