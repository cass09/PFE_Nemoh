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
MODULE MIdentification

  ! Definition of TYPE TID
  TYPE TID
    CHARACTER(LEN=80) :: ID
    INTEGER :: lID
  END TYPE TID

#ifdef _WIN32
  CHARACTER :: FILESEP = '\'
#else
  CHARACTER :: FILESEP = '/'
#endif

CONTAINS

  SUBROUTINE ReadTID(ID)
    ! Read ID data from commad line argument

    IMPLICIT NONE
    TYPE(TID) :: ID

    IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
      CALL GET_COMMAND_ARGUMENT(1, ID%ID)
    ELSE
      ID%ID = "."
    END IF
    ID%lID = LEN(TRIM(ID%ID)) ! for legacy...

  END SUBROUTINE

END MODULE MIdentification
