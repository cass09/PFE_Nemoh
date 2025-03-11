!--------------------------------------------------------------------------------------
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
!    - Ruddy Kurnia (ECN)
!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - Postprocessing module
!   MAIN PROGRAM
!   Description:  producing QTF output as the WAMIT QTF output format
!
!--------------------------------------------------------------------------------------
PROGRAM Main
!       Used modules
        USE MIdentification
        USE MNemohCal,          ONLY: TNemCal,READ_TNEMOHCAL
        USE MQPostproc,         ONLY: READWRITE_QTFDATA
!
        IMPLICIT NONE
!
!       Variables declaration
        TYPE(TID)       :: ID
        TYPE(TNemCal)   :: inpNEMOHCAL
!
!       Read input files
        CALL ReadTID(ID)
        CALL READ_TNEMOHCAL(ID,inpNEMOHCAL)
        CALL READWRITE_QTFDATA(inpNEMOHCAL,TRIM(ID%ID))

END PROGRAM Main
