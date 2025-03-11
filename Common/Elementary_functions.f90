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
!   - G. Delhommeau
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh
!	  - R. Kurnia
!--------------------------------------------------------------------------------------
MODULE Elementary_functions

  IMPLICIT NONE

  PUBLIC :: GG, CIH, SIH, CROSS_PRODUCT,PL2,PL5

CONTAINS

  COMPLEX FUNCTION GG(Z, CEX)
    ! Estimation of ∫_z^∞ exp(-t)/t dt
    ! See p.367 of G. Delhommeau thesis (referenced as [Del]).

    COMPLEX, INTENT(IN) :: Z, CEX
    COMPLEX             :: Y

    IF (REAL(Z) < -16.0) THEN                                      ! Case 1 p. 368 in [Del]
      Y = 1./Z
      GG = Y*(1.+Y*(-1.+Y*(2.+Y*(-6.+Y*(24.+Y*(-120.))))))
    ELSE IF (ABS(AIMAG(Z)) > 10.0) THEN                            ! Case 3 p. 368 in [Del]
      GG = 0.711093/(Z+0.415775)+0.278518/(Z+2.29428)+0.010389/(Z+6.2900)
      IF (AIMAG(Z) < 0) THEN
        GG = GG-(0., 3.14159265)*CEX
      ELSE
        GG = GG+(0., 3.14159265)*CEX
      END IF
    ELSE IF (REAL(Z) > -0.5) THEN                                  ! Case 2 p. 368 in [Del]
      GG = -(CLOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+    &
        Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+                    &
        Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+                  &
        Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+                     &
        Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01 &
        +Z*(-0.8700861E-03+Z*0.2989204E-03))))))
      IF (AIMAG(Z) < 0) THEN
        GG = GG-(0., 3.14159265)*CEX
      ELSE
        GG = GG+(0., 3.14159265)*CEX
      END IF
    ELSE                                                           ! Case 4 p. 369 in [Del]
      IF (AIMAG(Z) < 0) THEN
        GG = ((((((( (1.000000, 1.3935496E-06)*Z+ (15.82958, -20.14222))  &
          *Z+ (-70.52863, -227.9511))*Z+ (-985.4221, -226.6272))*Z        &
          + (-1202.318, 1580.907))*Z+ (953.2441, 1342.447))*Z             &
          + (417.3716, -196.6665))*Z+ (-9.881266, -24.24952))/            &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, -20.14481))*Z  &
          + (-55.66969, -248.1167))*Z+ (-1068.640, -434.4707))*Z          &
          + (-2082.250, 1522.471))*Z+ (383.3455, 2730.378))*Z             &
          + (1216.791, 351.7189))*Z+ (115.3926, -161.2647))*Z             &
          + (-3.777369, -4.510900))-(0., 3.14159265)*CEX
      ELSE
        GG = ((((((( (1.000000, -1.3935496E-06)*Z+ (15.82958, 20.14222))  &
          *Z+ (-70.52863, 227.9511))*Z+ (-985.4221, 226.6272))*Z          &
          + (-1202.318, -1580.907))*Z+ (953.2441, -1342.447))*Z           &
          + (417.3716, 196.6665))*Z+ (-9.881266, 24.24952))/              &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, 20.14481))*Z   &
          + (-55.66969, 248.1167))*Z+ (-1068.640, 434.4707))*Z            &
          + (-2082.250, -1522.471))*Z+ (383.3455, -2730.378))*Z           &
          + (1216.791, -351.7189))*Z+ (115.3926, 161.2647))*Z             &
          + (-3.777369, 4.510900))+(0., 3.14159265)*CEX
      END IF
    END IF
  END FUNCTION

  !-------------------------------------------------------------------------------!

  REAL FUNCTION CIH(AK,Z,H)
    REAL, INTENT(IN) :: AK,Z,H

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      CIH=COSH(AK*(Z+H))/COSH(AK*H)
    ELSE
      CIH=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

  REAL FUNCTION SIH(AK,Z,H)
    REAL, INTENT(IN) :: AK,Z,H

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      SIH=SINH(AK*(Z+H))/COSH(AK*H)
    ELSE
      SIH=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

  FUNCTION CIH_Vect(AK,Z,H,NZ) result(CIHV)
    INTEGER           , INTENT(IN) :: NZ
    REAL              , INTENT(IN) :: AK,H
    REAL,DIMENSION(NZ), INTENT(IN) :: Z
    REAL,DIMENSION(NZ)             :: CIHV

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      CIHV=COSH(AK*(Z+H))/COSH(AK*H)
    ELSE
      CIHV=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

   FUNCTION SIH_Vect(AK,Z,H,NZ) result(SIHV)
    INTEGER           , INTENT(IN) :: NZ
    REAL              , INTENT(IN) :: AK,H
    REAL,DIMENSION(NZ), INTENT(IN) :: Z
    REAL,DIMENSION(NZ)             :: SIHV

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      SIHV=SINH(AK*(Z+H))/COSH(AK*H)
    ELSE
      SIHV=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!
  FUNCTION COMPLEX_CONJUGATE_VECT(var,Nvect) RESULT(Conj)
           INTEGER,                  INTENT(IN)::NVect
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var
           COMPLEX,DIMENSION(Nvect) :: Conj
           Conj=CMPLX(REAL(var),-AIMAG(var))
  END FUNCTION

  FUNCTION CROSS_PRODUCT(A, B)
    REAL, DIMENSION(3) :: CROSS_PRODUCT
    REAL, DIMENSION(3), INTENT(IN) :: A, B

    CROSS_PRODUCT(1) = A(2)*B(3) - A(3)*B(2)
    CROSS_PRODUCT(2) = A(3)*B(1) - A(1)*B(3)
    CROSS_PRODUCT(3) = A(1)*B(2) - A(2)*B(1)
  END FUNCTION CROSS_PRODUCT

  FUNCTION CROSS_PRODUCT_COMPLEX(A, B) RESULT(CROSS_PRODUCT)
    COMPLEX, DIMENSION(3) :: CROSS_PRODUCT
    COMPLEX, DIMENSION(3), INTENT(IN) :: A, B

    CROSS_PRODUCT(1) = A(2)*B(3) - A(3)*B(2)
    CROSS_PRODUCT(2) = A(3)*B(1) - A(1)*B(3)
    CROSS_PRODUCT(3) = A(1)*B(2) - A(2)*B(1)
  END FUNCTION CROSS_PRODUCT_COMPLEX

  FUNCTION DOT_PRODUCT_COMPLEX(var1,var2,Nvect) RESULT(prod)
       INTEGER, INTENT(IN)                  :: Nvect
       COMPLEX,DIMENSION(Nvect), INTENT(IN) :: var1,var2
       COMPLEX                              ::prod
       !NOTE: DOT_PRODUCT(A,B)=A*.B for A,B complex variables
       !we want to have DOT_PRODUCT(A,B)=A.B
       prod=DOT_PRODUCT(COMPLEX_CONJUGATE_VECT(var1,Nvect),var2)
  END FUNCTION

  !-------------------------------------------------------------------------------!

  REAL FUNCTION X0(Y)
    ! Solve x * tanh(x) = y by dichotomy

    REAL, INTENT(IN) :: Y
    REAL             :: X_LOW, X_MEAN, X_UP
    REAL             :: EPS, STEP

    X_LOW = 0.0

    ! Find upper bound for x
    X_UP = 0.0
    STEP = MAX(Y, SQRT(Y))
    DO WHILE (X_UP*TANH(X_UP) < Y)
      X_UP = X_UP + STEP
    END DO

    ! Dichotomy
    EPS = 5.E-6
    DO WHILE (X_UP - X_LOW > EPS*X_UP)
      X_MEAN = (X_LOW + X_UP)/2
      IF (X_MEAN*TANH(X_MEAN) < Y) THEN
        X_LOW = X_MEAN
      ELSE
        X_UP = X_MEAN
      END IF
    END DO

    X0 = X_MEAN

    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

  REAL FUNCTION PL2(U1,U2,U3,XU)
    REAL::U1,U2,U3,XU
    PL2=((XU-U1)*(XU-U2))/((U3-U1)*(U3-U2))
    RETURN
  END FUNCTION

  REAL FUNCTION PL5(U1,U2,U3,U4,U5,XU)
    REAL::U1,U2,U3,U4,U5,XU
    PL5=((XU-U1)*(XU-U2)*(XU-U3)*(XU-U4))/&
      ((U5-U1)*(U5-U2)*(U5-U3)*(U5-U4))
    RETURN
  END FUNCTION

  FUNCTION Fun_closest(N,w,wobs) result(Iwobs)
    INTEGER, INTENT(IN):: N
    REAL, INTENT(IN)   :: w(N),wobs
    INTEGER            :: Iwobs,Iw,Nw
    REAL               :: mindist,mindistN

    mindist=w(N)-w(1)
    DO Iw=1,N
        IF (abs(w(Iw)-wobs)<=mindist)THEN
                Iwobs=Iw
                mindist=abs(w(Iw)-wobs)
        ENDIF
    ENDDO
    RETURN
  END FUNCTION

  FUNCTION Fun_MIN(N,vect) result(minvalue)
    INTEGER,            INTENT(IN):: N
    REAL,DIMENSION(N),  INTENT(IN):: vect
    REAL                          :: minvalue
    INTEGER                       ::I
    minvalue=vect(1)
    DO I=1,N
    minvalue=MIN(vect(I),minvalue)
    ENDDO
  END FUNCTION

  FUNCTION Fun_MAX(N,vect) result(maxvalue)
    INTEGER,            INTENT(IN):: N
    REAL,DIMENSION(N),  INTENT(IN):: vect
    REAL                          :: maxvalue
    INTEGER                       :: I
    maxvalue=vect(1)
    DO I=1,N
    maxvalue=MAX(vect(I),maxvalue)
    ENDDO
  END FUNCTION

  FUNCTION Fun_KronDelta(m,n) result(delta)
        INTEGER, INTENT(IN) :: m,n
        INTEGER             :: delta
        IF (m.EQ.n) delta=1
        IF (m.NE.n) delta=0
  END FUNCTION

END MODULE Elementary_functions
