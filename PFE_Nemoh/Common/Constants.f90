!--------------------------------------------------------------------------------------
!
!   NEMOH - March 2022
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
!--------------------------------------------------------------------------------------
MODULE CONSTANTS

      IMPLICIT NONE

      REAL, PARAMETER           :: PI=4.*ATAN(1.), DPI=2.*PI, DPI2=2.*PI**2
      COMPLEX, PARAMETER        :: II=CMPLX(0.,1.)
      REAL, PARAMETER           :: INFINITE_DEPTH=0.0
      INTEGER, PARAMETER        :: NO_Y_SYMMETRY=0, Y_SYMMETRY=1
      REAL, PARAMETER           :: ZERO=0.
      COMPLEX, PARAMETER        :: CZERO=CMPLX(0.0,0.0)
      INTEGER, PARAMETER        :: DIFFRACTION_PROBLEM=1, RADIATION_PROBLEM=-1
      INTEGER, PARAMETER        :: ID_BODY=0,ID_FREESURFACE=1
END MODULE CONSTANTS
