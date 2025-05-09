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
!
!   Contributors list:
!   - R. Kurnia
!--------------------------------------------------------------------------------------

MODULE MCallInterp

USE linear_interpolation_module

IMPLICIT NONE

CONTAINS

 FUNCTION FUNVECT_INTERP1_REAL(X,VAR,NX,XOUT,NXOUT) RESULT(VAROUT)
   INTEGER,            INTENT(IN) :: NX,NXOUT
   REAL,DIMENSION(NX), INTENT(IN) :: X,VAR
   REAL,DIMENSION(NXOUT)          :: XOUT,VAROUT
   Type(linear_interp_1d)         :: interp1
   INTEGER                        :: iflag,I

   CALL interp1%initialize(X,VAR,iflag)
   DO I=1,NXOUT
     CALL interp1%evaluate(XOUT(I), VAROUT(I))
   ENDDO
   CALL interp1%destroy()
 END FUNCTION

 FUNCTION FUNVECT_INTERP1_COMPLEX(X,VAR,NX,XOUT,NXOUT) RESULT(VAROUT)
   INTEGER,               INTENT(IN) :: NX,NXOUT
   REAL,DIMENSION(NX),    INTENT(IN) :: X
   COMPLEX,DIMENSION(NX), INTENT(IN) :: VAR

   REAL,DIMENSION(NXOUT)             :: XOUT
   COMPLEX,DIMENSION(NXOUT)          :: VAROUT
   Type(linear_interp_1d)            :: interp1
   INTEGER                           :: iflag
   REAL,DIMENSION(NXOUT)             :: VAROUT_R,VAROUT_I

   VAROUT_R=FUNVECT_INTERP1_REAL(X,REAL(VAR),NX,XOUT,NXOUT)
   VAROUT_I=FUNVECT_INTERP1_REAL(X,AIMAG(VAR),NX,XOUT,NXOUT)
   VAROUT  =CMPLX(VAROUT_R,VAROUT_I)
 END FUNCTION

 FUNCTION FUN_INTERP1_REAL(X,VAR,NX,XOUT) RESULT(VAROUT)
   INTEGER,            INTENT(IN) :: NX
   REAL,DIMENSION(NX), INTENT(IN) :: X,VAR
   REAL                           :: XOUT,VAROUT
   Type(linear_interp_1d)         :: interp1
   INTEGER                        :: iflag,I

   CALL interp1%initialize(X,VAR,iflag)
   CALL interp1%evaluate(XOUT, VAROUT)
   CALL interp1%destroy()
 END FUNCTION

 FUNCTION FUN_INTERP1_COMPLEX(X,VAR,NX,XOUT) RESULT(VAROUT)
   INTEGER,               INTENT(IN) :: NX
   REAL,DIMENSION(NX),    INTENT(IN) :: X
   COMPLEX,DIMENSION(NX), INTENT(IN) :: VAR

   REAL                              :: XOUT
   COMPLEX                           :: VAROUT
   Type(linear_interp_1d)            :: interp1
   INTEGER                           :: iflag
   REAL                              :: VAROUT_R,VAROUT_I

   VAROUT_R=FUN_INTERP1_REAL(X,REAL(VAR),NX,XOUT)
   VAROUT_I=FUN_INTERP1_REAL(X,AIMAG(VAR),NX,XOUT)
   VAROUT=CMPLX(VAROUT_R,VAROUT_I)
 END FUNCTION


END MODULE
