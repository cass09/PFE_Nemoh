!-------------------------------------------------------------------------------------
!
!    Copyright (C) 2022 - 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
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
!   - Fabien Robaux (EDF/INNOSEA)
!   - G.DELHOMMEAU, THESE DE CHEN XIAO-BO(1988)
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!
!--------------------------------------------------------------------------------------
MODULE MQSOLVERASYMP
!all functions used in ASYMP calculation

USE CONSTANTS          , ONLY:II,CZERO,PI
USE MEnvironment,        ONLY:Fun_Dispersion
USE Elementary_functions,ONLY:CIH,Fun_KronDelta
USE MBESSEL,             ONLY:fun_BESSJ  
USE MROMBERG,            ONLY:romberg_trap
IMPLICIT NONE

CONTAINS

  FUNCTION Fun_IDF(w1,w2,k1,k2,delk,sumk,g,D,Rf,NRf,Nbessel,             &
                  Ivartheta1l,Ivartheta2l,IR1l,IR2l)  result(IDF)

     REAL,                              INTENT(IN):: k1,k2,w1,w2,delk,sumk,g,D
     INTEGER,                           INTENT(IN):: NRf, Nbessel
     REAL, DIMENSION(NRf),              INTENT(IN):: Rf
     COMPLEX,DIMENSION(2,NBESSEL+1),    INTENT(IN):: IR1l,IR2l
     COMPLEX,DIMENSION(2,NBESSEL+3),    INTENT(IN):: Ivartheta1l,Ivartheta2l


     COMPLEX,DIMENSION(2,2)                       :: IDF
     COMPLEX,DIMENSION(2)                         :: IDF11,IDF12
     COMPLEX,DIMENSION(2)                         :: IDF21,IDF22
     INTEGER                                      :: Ibessel
     COMPLEX,DIMENSION(2)                         :: sumIRIVartheta_11,sumIRIVartheta_12
     COMPLEX,DIMENSION(2)                         :: sumIRIVartheta_21,sumIRIVartheta_22

     sumIRIVartheta_11(:)=CZERO
     sumIRIVartheta_12(:)=CZERO
     sumIRIVartheta_21(:)=CZERO
     sumIRIVartheta_22(:)=CZERO

     DO Ibessel=1,Nbessel+1
        sumIRIVartheta_11=sumIRIVartheta_11                                             &
               +0.5*IR1l(:,Ibessel)*(Ivartheta1l(:,Ibessel)+Ivartheta1l(:,Ibessel+2))
        sumIRIVartheta_12=sumIRIVartheta_12                                             &
               +0.5*IR2l(:,Ibessel)*(Ivartheta2l(:,Ibessel)+Ivartheta2l(:,Ibessel+2))
        sumIRIVartheta_21=sumIRIVartheta_21                                             &
               +IR1l(:,Ibessel)*Ivartheta1l(:,Ibessel+1)
        sumIRIVartheta_22=sumIRIVartheta_22                                             &
               +IR2l(:,Ibessel)*Ivartheta2l(:,Ibessel+1)
     ENDDO
     IDF11(1)=(-II*g*8*PI*SQRT(k2*delk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_11(1)

     IDF11(2)=(-II*g*8*PI*SQRT(k2*sumk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_11(2)

     IDF12(1)=(II*g*8*PI*SQRT(k1*delk)/w2)                                              &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_12(1)

     IDF12(2)=(-II*g*8*PI*SQRT(k1*sumk)/w2)                                             &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_12(2)

     IDF21(1)=(-II*g*8*PI*SQRT(k2*delk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_21(1)

     IDF21(2)=(-II*g*8*PI*SQRT(k2*sumk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_21(2)

     IDF22(1)=(II*g*8*PI*SQRT(k1*delk)/w2)                                              &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_22(1)

     IDF22(2)=(-II*g*8*PI*SQRT(k1*sumk)/w2)                                             &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_22(2)
     IDF(:,1)=IDF11+IDF12
     IDF(:,2)=IDF21+IDF22
     IF (delk==0.) IDF(1,:)=CZERO
  END FUNCTION

  FUNCTION Fun_KAPPA2_DIFFSUM(k1,k2,w1,w2,D,g) RESULT(KAPPA2)
     REAL,        INTENT(IN):: k1,k2,w1,w2,g,D
     COMPLEX,DIMENSION(2)   :: KAPPA2
     REAL                   :: Nu1,Nu2,SECH2_K1D,SECH2_K2D

     Nu1=Fun_Dispersion(k1,D,g)**2/g
     Nu2=Fun_Dispersion(k2,D,g)**2/g
     IF (D.EQ.0.) THEN !INFINITE DEPTH
     SECH2_K1D=0
     SECH2_K2D=0
     ELSE              !FINITE DEPTH
     SECH2_K1D=1-tanh(k1*D)**2
     SECH2_K2D=1-tanh(k2*D)**2
     ENDIF
     KAPPA2(1)=II*(w1-w2)*Nu1*Nu2+II*w1*w2/g*(k1**2*SECH2_K1D/w1-k2**2*SECH2_K2D/w2)
     KAPPA2(2)=II*(w1+w2)*Nu1*Nu2-II*w1*w2/g*(k1**2*SECH2_K1D/w1+k2**2*SECH2_K2D/w2)

  END FUNCTION

  FUNCTION Fun_Fprofile(k,D,z) result(F)
     REAL,      INTENT(IN):: k,D,z
     REAL                 :: F

     IF(D.EQ.0..OR.(k*(D+z).GT.50)) THEN !INFINITE DEPTH case
         F=exp(k*z)
     ELSE
         F=cosh(k*(D+z))/( (k*D*(1-tanh(k*D)**2)+tanh(k*D))*cosh(k*D) )
     ENDIF
  END FUNCTION

  FUNCTION Fun_IR1l(l,k1,k2,delk,sumk,Rf,NRf) result(IR1l)
     INTEGER,             INTENT(IN):: l,NRf
     REAL,                INTENT(IN):: k1,k2,delk,sumk
     REAL,DIMENSION(NRf), INTENT(IN):: Rf
     INTEGER               :: eps_l
     COMPLEX ,DIMENSION(2) :: IR1_0Rext
     COMPLEX ,DIMENSION(2) :: IR1l
     INTEGER               :: NR,Niter
     REAL                  :: rel_error,dR
     COMPLEX               :: IR_0INF,IR_RFINF


     eps_l=Fun_epsilon_l(l)
     rel_error=0.01
     NR=51;Niter=100000

     IF (delk>0) THEN
     dR=2*PI/MAX(delk-k2,k1)/NR
     IR_0INF=fun_gamma(l,delk-k2,k1)
     IR_RFINF=Fun_IntegIR_RFInf(l,rel_error,Niter,dR,Rf(NRf),NRf,k1,k2,delk,IR_0INF,11)
     IR1l(1)=eps_l*(II**l)*IR_RFINF       !for diff freq
     ELSE
        IR1l(1)=CZERO
     ENDIF


     dR=2*PI/MAX(sumk+k2,k1)/NR
     IR_0INF=II*fun_gamma(l,sumk+k2,k1)
     IR_RFINF=Fun_IntegIR_RFInf(l,rel_error,Niter,dR,Rf(NRf),NRf,k1,k2,sumk,IR_0INF,12)
     IR1l(2)=eps_l*(II**l)*IR_RFINF       !for sum freq
        !print*,l,IR1_0Rext(:)
    !print*,l,delk-k1,k1,fun_gamma(l,delk-k2,k1)
  END FUNCTION

  FUNCTION Fun_IR2l(l,k1,k2,delk,sumk,Rf,NRf) result(IR2l)
     INTEGER,             INTENT(IN):: l,NRf
     REAL,                INTENT(IN):: k1,k2,delk,sumk
     REAL,DIMENSION(NRf), INTENT(IN):: Rf
     INTEGER               :: eps_l
     COMPLEX ,DIMENSION(2) :: IR2_0Rext
     COMPLEX ,DIMENSION(2) :: IR2l
     INTEGER               :: NR,Niter
     REAL                  :: rel_error,dR
     COMPLEX               :: IR_0INF,IR_RFINF


     eps_l=Fun_epsilon_l(l)

     rel_error=0.01
     NR=51;Niter=100000

     IF (delk>0) THEN
     dR=2*PI/MAX(delk+k1,k2)/NR
     IR_0INF=II*fun_gamma(l,delk+k1,k2)
     IR_RFINF=Fun_IntegIR_RFInf(l,rel_error,Niter,dR,Rf(NRf),NRf,k1,k2,delk,IR_0INF,21)
     IR2l(1)=eps_l*CONJG(II**l)*IR_RFINF       !for diff freq
     ELSE
     IR2l(1)=CZERO
     ENDIF


     dR=2*PI/MAX(sumk+k1,k2)/NR
     IR_0INF=II*fun_gamma(l,sumk+k1,k2)
     IR_RFINF=Fun_IntegIR_RFInf(l,rel_error,Niter,dR,Rf(NRf),NRf,k1,k2,sumk,IR_0INF,22)
     IR2l(2)=eps_l*(II**l)*IR_RFINF       !for sum freq

   END FUNCTION




  FUNCTION Fun_IntegIR_RFInf(l,rel_error,Niter,dR,Rf,NRf,k1,k2,delkORsumk,IR_0INF,Iswitch) result(IR_RFINF)
     IMPLICIT NONE
     COMPLEX,  INTENT(IN) :: IR_0INF
     REAL,     INTENT(IN) :: rel_error,dR,k1,k2,delkORsumk,Rf
     INTEGER,  INTENT(IN) :: Niter,Iswitch,NRf,l
     COMPLEX              :: IR_RFINF

     COMPLEX              :: IR_0RF,IR_0RF1
     REAL                 :: comp_rel_error
     INTEGER              :: iter
     INTEGER              :: NR0RF

     REAL,DIMENSION(3)    :: paramR
     INTEGER              :: paramI
     REAL,ALLOCATABLE,DIMENSION(:)    :: R0RF
     REAL(kind=8)         :: R0,R1
     REAL(kind=8)         :: parRomberg_TOL
     INTEGER(kind=4)      :: parRomberg_M

     NR0RF  =MAX(INT(Rf/dR+1),NRf)
     ALLOCATE(R0RF(NR0RF))
     R0RF=Fun_discretized_R0R1(0.,Rf,NR0RF)

     paramR(1)=k1
     paramR(2)=k2
     paramR(3)=delkORsumk
     paramI   =l
     parRomberg_TOL=1E-8

       R0=0.
       R1=real(Rf,kind=8)

       IF (Iswitch==11) THEN
       IR_0RF=Fun_IR1M_R0_R1(l,k1,k2,delkORsumk,R0RF,NR0RF)
      ! IR_0RF=Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand1M_REAL)      &
      !       +II*Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand1M_IMAG)
      ! IR_0RF=real(romberg_trap(R0,R1,Fun_Integrand1M_REAL_R8, parRomberg_TOL,parRomberg_M))&
      !    +II*real(romberg_trap(R0,R1,Fun_Integrand1M_IMAG_R8, parRomberg_TOL,parRomberg_M))
       ELSEIF (Iswitch==12) THEN
       IR_0RF=Fun_IR1P_R0_R1(l,k1,k2,delkORsumk,R0RF,NR0RF)
     !  IR_0RF=Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand1P_REAL)      &
     !        +II*Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand1P_IMAG)
     !  IR_0RF=real(romberg_trap(R0,R1,Fun_Integrand1P_REAL_R8, parRomberg_TOL,parRomberg_M))&
     !     +II*real(romberg_trap(R0,R1,Fun_Integrand1P_IMAG_R8, parRomberg_TOL,parRomberg_M))
       ELSEIF (Iswitch==21) THEN
       IR_0RF=Fun_IR2M_R0_R1(l,k1,k2,delkORsumk,R0RF,NR0RF)
     !  IR_0RF=Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand2M_REAL)      &
     !       +II*Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand2M_IMAG)
     !  IR_0RF=real(romberg_trap(R0,R1,Fun_Integrand2M_REAL_R8, parRomberg_TOL,parRomberg_M))&
     !     +II*real(romberg_trap(R0,R1,Fun_Integrand2M_IMAG_R8, parRomberg_TOL,parRomberg_M))
       ELSEIF (Iswitch==22) THEN
       IR_0RF=Fun_IR2P_R0_R1(l,k1,k2,delkORsumk,R0RF,NR0RF)
     !  IR_0RF=Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand2P_REAL)      &
     !        +II*Fun_INTEGRATION_TRAPZ(R0RF,NR0RF,Fun_Integrand2P_IMAG)
     !  IR_0RF=real(romberg_trap(R0,R1,Fun_Integrand2P_REAL_R8, parRomberg_TOL,parRomberg_M))&
     !     +II*real(romberg_trap(R0,R1,Fun_Integrand2P_IMAG_R8, parRomberg_TOL,parRomberg_M))
       ENDIF

       IR_RFINF=IR_0INF-IR_0RF

     DEALLOCATE(R0RF)

     CONTAINS
        FUNCTION Fun_Integrand1M_REAL_R8(r) result(f1M_Re)
              REAL(kind=8),  INTENT(IN):: r
              REAL(kind=8)             :: f1M_Re
              REAL                     :: k1,k2,delk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI
              f1M_Re=real(cos((delk-k2)*real(r))*fun_BESSJ(l,k1*real(r)),kind=8)
        END FUNCTION

        FUNCTION Fun_Integrand1M_IMAG_R8(r) result(f1M_Im)
              REAL(kind=8),  INTENT(IN):: r
              REAL(kind=8)             :: f1M_Im
              REAL                     :: k1,k2,delk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI
              f1M_Im=real(sin((delk-k2)*real(r))*fun_BESSJ(l,k1*real(r)),kind=8)
        END FUNCTION

        FUNCTION Fun_Integrand1P_REAL_R8(r) result(f1P_Re)
              REAL(kind=8),  INTENT(IN):: r
              REAL(kind=8)             :: f1P_Re
              REAL                     :: k1,k2,sumk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f1P_Re=real(cos((sumk+k2)*real(r)+PI/2 )*fun_BESSJ(l,k1*real(r)),kind=8)
        END FUNCTION

        FUNCTION Fun_Integrand1P_IMAG_R8(r) result(f1P_Im)
              REAL(kind=8),   INTENT(IN):: r
              REAL(kind=8)             :: f1P_Im
              REAL                     :: k1,k2,sumk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f1P_Im=real(sin((sumk+k2)*real(r)+PI/2 )*fun_BESSJ(l,k1*real(r)),kind=8)
        END FUNCTION


        FUNCTION Fun_Integrand2M_REAL_R8(r) result(f2M_Re)
              REAL(kind=8),   INTENT(IN):: r
              REAL(kind=8)             :: f2M_Re
              REAL                     :: k1,k2,delk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI
              f2M_Re=real(cos((delk+k1)*real(r)+PI/2 )*fun_BESSJ(l,k2*real(r)),kind=8)
        END FUNCTION

        FUNCTION Fun_Integrand2M_IMAG_R8(r) result(f2M_Im)
              REAL(kind=8),   INTENT(IN):: r
              REAL(kind=8)             :: f2M_Im
              REAL                     :: k1,k2,delk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI
              f2M_Im=real(sin((delk+k1)*real(r)+PI/2 )*fun_BESSJ(l,k2*real(r)),kind=8)
        END FUNCTION

        FUNCTION Fun_Integrand2P_REAL_R8(r) result(f2P_Re)
              REAL(kind=8),   INTENT(IN):: r
              REAL(kind=8)             :: f2P_Re
              REAL                     :: k1,k2,sumk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f2P_Re=real(cos((sumk+k1)*real(r)+PI/2 )*fun_BESSJ(l,k2*real(r)),kind=8)
        END FUNCTION

        FUNCTION Fun_Integrand2P_IMAG_R8(r) result(f2P_Im)
              REAL(kind=8),   INTENT(IN):: r
              REAL(kind=8)             :: f2P_Im
              REAL                     :: k1,k2,sumk
              INTEGER                  :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f2P_Im=real(sin((sumk+k1)*real(r)+PI/2 )*fun_BESSJ(l,k2*real(r)),kind=8)
        END FUNCTION

        FUNCTION Fun_Integrand1M_REAL(r) result(f1M_Re)
              REAL,   INTENT(IN):: r
              REAL              :: f1M_Re
              REAL              :: k1,k2,delk
              INTEGER           :: l
              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI
              f1M_Re=cos((delk-k2)*r)*fun_BESSJ(l,k1*r)
        END FUNCTION

        FUNCTION Fun_Integrand1M_IMAG(r) result(f1M_Im)
              REAL,   INTENT(IN):: r
              REAL              :: f1M_Im
              REAL              :: k1,k2,delk
              INTEGER           :: l

              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI

             f1M_Im=sin((delk-k2)*r)*fun_BESSJ(l,k1*r)
        END FUNCTION

        FUNCTION Fun_Integrand1P_REAL(r) result(f1P_Re)
              REAL,   INTENT(IN):: r
              REAL              :: f1P_Re
              REAL              :: k1,k2,sumk
              INTEGER           :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f1P_Re=cos((sumk+k2)*r+PI/2 )*fun_BESSJ(l,k1*r)
        END FUNCTION

        FUNCTION Fun_Integrand1P_IMAG(r) result(f1P_Im)
              REAL,   INTENT(IN):: r
              REAL              :: f1P_Im
              REAL              :: k1,k2,sumk
              INTEGER           :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f1P_Im=sin((sumk+k2)*r+PI/2 )*fun_BESSJ(l,k1*r)
        END FUNCTION

        FUNCTION Fun_Integrand2M_REAL(r) result(f2M_Re)
              REAL,   INTENT(IN):: r
              REAL              :: f2M_Re
              REAL              :: k1,k2,delk
              INTEGER           :: l
              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI

              f2M_Re=cos((delk+k1)*r+PI/2 )*fun_BESSJ(l,k2*r)
        END FUNCTION

        FUNCTION Fun_Integrand2M_IMAG(r) result(f2M_Im)
              REAL,   INTENT(IN):: r
              REAL              :: f2M_Im
              REAL              :: k1,k2,delk
              INTEGER           :: l
              k1  =paramR(1)
              k2  =paramR(2)
              delk=paramR(3)
              l   =paramI

              f2M_Im=sin((delk+k1)*r+PI/2 )*fun_BESSJ(l,k2*r)
        END FUNCTION

        FUNCTION Fun_Integrand2P_REAL(r) result(f2P_Re)
              REAL,   INTENT(IN):: r
              REAL              :: f2P_Re
              REAL              :: k1,k2,sumk
              INTEGER           :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f2P_Re=cos((sumk+k1)*r+PI/2 )*fun_BESSJ(l,k2*r)
        END FUNCTION

        FUNCTION Fun_Integrand2P_IMAG(r) result(f2P_Im)
              REAL,   INTENT(IN) :: r
              REAL               :: f2P_Im
              REAL              :: k1,k2,sumk
              INTEGER           :: l
              k1  =paramR(1)
              k2  =paramR(2)
              sumk=paramR(3)
              l   =paramI

              f2P_Im=sin((sumk+k1)*r+PI/2 )*fun_BESSJ(l,k2*r)
        END FUNCTION


        FUNCTION Fun_IR1M_R0_R1(l,k1,k2,delk,Rf,NRf) result(IR1)
            INTEGER,             INTENT(IN):: l,NRf
            REAL,                INTENT(IN):: k1,k2,delk
            REAL,DIMENSION(NRf), INTENT(IN):: Rf
            COMPLEX                        :: IR1
            REAL                           :: dRf
            INTEGER ::Ir

            dRf=Rf(2)-Rf(1)
            IR1=CZERO
            DO Ir=1,NRf-1
               IR1=IR1+0.5*(                                                           &
                       exp(II*(delk-k2)*Rf(Ir))*fun_BESSJ(l,k1*Rf(Ir))                 &
                      +exp(II*(delk-k2)*Rf(Ir+1))*fun_BESSJ(l,k1*Rf(Ir+1)) )*dRf
            ENDDO
        END FUNCTION


        FUNCTION Fun_IR1P_R0_R1(l,k1,k2,sumk,Rf,NRf) result(IR1)
           INTEGER,             INTENT(IN):: l,NRf
           REAL,                INTENT(IN):: k1,k2,sumk
           REAL,DIMENSION(NRf), INTENT(IN):: Rf
           COMPLEX                        :: IR1
           REAL                           :: dRf
           INTEGER ::Ir

           dRf=Rf(2)-Rf(1)
           IR1=CZERO
           DO Ir=1,NRf-1
              IR1=IR1+0.5*(                                                           &
                    exp(II*( (sumk+k2)*Rf(Ir)+PI/2   ))*fun_BESSJ(l,k1*Rf(Ir))        &
                   +exp(II*( (sumk+k2)*Rf(Ir+1)+PI/2 ))*fun_BESSJ(l,k1*Rf(Ir+1)) )*dRf
           ENDDO
        END FUNCTION


        FUNCTION Fun_IR2M_R0_R1(l,k1,k2,delk,Rf,NRf) result(IR2)
           INTEGER,             INTENT(IN):: l,NRf
           REAL,                INTENT(IN):: k1,k2,delk
           REAL,DIMENSION(NRf), INTENT(IN):: Rf
           COMPLEX                        :: IR2
           REAL                           :: dRf
           INTEGER ::Ir

           dRf=Rf(2)-Rf(1)
           IR2=CZERO
           DO Ir=1,NRf-1
              IR2=IR2+0.5*(                                                           &
                   exp(II*( (delk+k1)*Rf(Ir)+PI/2   ))*fun_BESSJ(l,k2*Rf(Ir))         &
                  +exp(II*( (delk+k1)*Rf(Ir+1)+PI/2 ))*fun_BESSJ(l,k2*Rf(Ir+1)) )*dRf
           ENDDO
        END FUNCTION


        FUNCTION Fun_IR2P_R0_R1(l,k1,k2,sumk,Rf,NRf) result(IR2)
           INTEGER,             INTENT(IN):: l,NRf
           REAL,                INTENT(IN):: k1,k2,sumk
           REAL,DIMENSION(NRf), INTENT(IN):: Rf
           COMPLEX                        :: IR2
           REAL                           :: dRf
           INTEGER ::Ir

           dRf=Rf(2)-Rf(1)
           IR2=CZERO
           DO Ir=1,NRf-1
              IR2=IR2+0.5*(                                                           &
                   exp(II*( (sumk+k1)*Rf(Ir)+PI/2   ))*fun_BESSJ(l,k2*Rf(Ir))         &
                  +exp(II*( (sumk+k1)*Rf(Ir+1)+PI/2 ))*fun_BESSJ(l,k2*Rf(Ir+1)) )*dRf
           ENDDO
        END FUNCTION


  END FUNCTION

 ! FUNCTION Fun_INTEGRATION_TRAPZ(X,NX,FunX,paramR,paramI) result(INTEG)
 !       INTEGER,              INTENT(IN):: NX
 !       REAL,   DIMENSION(NX),INTENT(IN):: X
 !       REAL,                 EXTERNAL  :: FunX
 !       REAL,   DIMENSION(3), INTENT(IN):: paramR
 !       INTEGER,              INTENT(IN):: paramI
 !       INTEGER Ix
 !       REAL    INTEG
 !       INTEG=0.0
 !       DO Ix=1,NX-1
 !          INTEG=INTEG+0.5*                                             &
 !              (FunX(X(Ix),paramR(1),paramR(2),paramR(3),paramI)        &
 !              +FunX(X(Ix+1),paramR(1),paramR(2),paramR(3),paramI)      &
 !              )*(X(Ix+1)-X(Ix))
 !       ENDDO

 ! END FUNCTION

  FUNCTION Fun_INTEGRATION_TRAPZ(X,NX,FunX) result(INTEG)
        INTEGER,              INTENT(IN):: NX
        REAL,   DIMENSION(NX),INTENT(IN):: X
        REAL,                 EXTERNAL  :: FunX
        INTEGER Ix
        REAL    INTEG
        INTEG=0.0
        DO Ix=1,NX-1
           INTEG=INTEG+0.5*(FunX(X(Ix))+FunX(X(Ix+1)))*(X(Ix+1)-X(Ix))
        ENDDO

  END FUNCTION


  FUNCTION Fun_discretized_R0R1(R0,R1,NR) result(R)
    REAL,       INTENT(IN)::R0,R1
    INTEGER,    INTENT(IN)::NR
    REAL,DIMENSION(NR)    ::R
    INTEGER               ::IR
    REAL                  ::dR
    dR=(R1-R0)/(NR-1)
    DO IR=1,NR
        R(IR)=R0+(IR-1)*dR
    ENDDO

  END

  FUNCTION PREPARE_KOCHIN_COEFFICIENTS(Isym,Npanels,XM,Apanels,D,Nbessel,k,zig) &
                  result(CmSm)
   INTEGER,                             INTENT(IN) :: Isym,Npanels,Nbessel
   REAL,                                INTENT(IN) :: k,D
   REAL,DIMENSION(Npanels),             INTENT(IN) :: Apanels
   REAL,DIMENSION(3,Npanels),           INTENT(IN) :: XM
   COMPLEX,DIMENSION(Npanels*2**Isym)  ,INTENT(IN) :: zig

   COMPLEX,DIMENSION(2,Nbessel+1) :: CmSm
   INTEGER                      :: ll
   DO ll=1,Nbessel+1
   CmSm(:,ll)=Fun_KochinCoefs_l(Isym,Npanels,XM,Apanels,zig,k,D,ll-1)
   ENDDO
  END FUNCTION

  FUNCTION Fun_IVartheta1l(Nbessel,beta,CmSmPer,CnSnRad_delk,CnSnRad_sumk,ll)&
                                                         result(Ivartheta)
     INTEGER,                             INTENT(IN) :: ll,Nbessel
     REAL,                                INTENT(IN) :: beta
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CmSmPer
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_delk
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_sumk


     INTEGER                       :: mm,nn
     COMPLEX,DIMENSION(2)          :: Ivartheta
     COMPLEX,DIMENSION(2)          :: CnRad,SnRad
     COMPLEX                       :: CmPer,SmPer

     Ivartheta=CZERO
     DO mm=1,Nbessel+1
        CmPer=CmSmPer(1,mm)
        SmPer=CmSmPer(2,mm)
        DO nn=1,Nbessel+1
          CnRad(1)=CnSnRad_delk(1,nn)
          SnRad(1)=CnSnRad_delk(2,nn)
          CnRad(2)=CnSnRad_sumk(1,nn)
          SnRad(2)=CnSnRad_sumk(2,nn)
          !for diff freq
          Ivartheta(1)=Ivartheta(1)                                         &
              +CONJG(CmPer)*CnRad(1)*cos(ll*beta)*Fun_Dlmn(1,ll,mm-1,nn-1)  &
              +CONJG(CmPer)*SnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,nn-1,mm-1)  &
              +CONJG(SmPer)*CnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,mm-1,nn-1)  &
              +CONJG(SmPer)*SnRad(1)*cos(ll*beta)*Fun_Dlmn(0,nn-1,mm-1,abs(ll))
          !for sum freq
          Ivartheta(2)=Ivartheta(2)                                         &
              +CmPer*CnRad(2)*cos(ll*beta)*Fun_Dlmn(1,ll,mm-1,nn-1)         &
              +CmPer*SnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,nn-1,mm-1)         &
              +SmPer*CnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,mm-1,nn-1)         &
              +SmPer*SnRad(2)*cos(ll*beta)*Fun_Dlmn(0,nn-1,mm-1,abs(ll))
        ENDDO
     ENDDO
  END FUNCTION

  FUNCTION Fun_IVartheta2l(Nbessel,beta,CmSmPer,CnSnRad_delk,CnSnRad_sumk,ll)&
                                                         result(Ivartheta)
     INTEGER,                             INTENT(IN) :: ll,Nbessel
     REAL,                                INTENT(IN) :: beta
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CmSmPer
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_delk
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_sumk


     INTEGER                       :: mm,nn
     COMPLEX,DIMENSION(2)          :: Ivartheta
     COMPLEX,DIMENSION(2)          :: CnRad,SnRad
     COMPLEX                       :: CmPer,SmPer

     Ivartheta=CZERO
     DO mm=1,Nbessel+1
        CmPer=CmSmPer(1,mm)
        SmPer=CmSmPer(2,mm)
        DO nn=1,Nbessel+1
          CnRad(1)=CnSnRad_delk(1,nn)
          SnRad(1)=CnSnRad_delk(2,nn)
          CnRad(2)=CnSnRad_sumk(1,nn)
          SnRad(2)=CnSnRad_sumk(2,nn)

          Ivartheta(1)=Ivartheta(1)                                     &
              +CmPer*CnRad(1)*cos(ll*beta)*Fun_Dlmn(1,ll,mm-1,nn-1)     &
              +CmPer*SnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,nn-1,mm-1)     &
              +SmPer*CnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,mm-1,nn-1)     &
              +SmPer*SnRad(1)*cos(ll*beta)*Fun_Dlmn(0,nn-1,mm-1,abs(ll))
          Ivartheta(2)=Ivartheta(2)                                     &
              +CmPer*CnRad(2)*cos(ll*beta)*Fun_Dlmn(1,ll,mm-1,nn-1)     &
              +CmPer*SnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,nn-1,mm-1)     &
              +SmPer*CnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,mm-1,nn-1)     &
              +SmPer*SnRad(2)*cos(ll*beta)*Fun_Dlmn(0,nn-1,mm-1,abs(ll))
        ENDDO
     ENDDO
  END FUNCTION


  FUNCTION Fun_KochinCoefs_l(Isym,Npanels,XM,AreaPanel,sig,k,D,ll) result(ClSl)
        INTEGER,                          INTENT(IN) :: Npanels,ll,Isym
        REAL,                             INTENT(IN) :: k,D
        REAL,DIMENSION(Npanels),          INTENT(IN) :: AreaPanel
        REAL,DIMENSION(3,Npanels),        INTENT(IN) :: XM
        COMPLEX,DIMENSION(Npanels*2*Isym),INTENT(IN) :: sig
        COMPLEX,DIMENSION(2)                         :: ClSl
        COMPlEX                                      :: calcClSl
        INTEGER                                      :: Ipanel
        REAL                                         :: r,alpha
        ClSl=CZERO
        DO Ipanel=1,Npanels
         ! IF (XM(3,Ipanel)<0) THEN
          r=sqrt(XM(1,Ipanel)**2+XM(2,Ipanel)**2)
          alpha=atan2(XM(2,Ipanel),XM(1,Ipanel))
          calcClSl=sig(Ipanel)*CIH(k,XM(3,Ipanel),D)                        &
                  *Fun_epsilon_l(ll)*(-II)**ll*fun_BESSJ(ll,k*r)            &
                  *AreaPanel(Ipanel)
          ClSl(1)=ClSl(1)+calcClSl*cos(ll*alpha)
          ClSl(2)=ClSl(2)+calcClSl*sin(ll*alpha)

          IF (Isym.EQ.1) THEN
            alpha=atan2(-XM(2,Ipanel),XM(1,Ipanel))
            calcClSl=sig(Npanels+Ipanel)*CIH(k,XM(3,Ipanel),D)               &
                    *Fun_epsilon_l(ll)*(-II)**ll*fun_BESSJ(ll,k*r)           &
                    *AreaPanel(Ipanel)
            ClSl(1)=ClSl(1)+calcClSl*cos(ll*alpha)
            ClSl(2)=ClSl(2)+calcClSl*sin(ll*alpha)
          ENDIF
         ! ENDIF
        ENDDO
        ClSl=-ClSl/4/PI
  END FUNCTION

  FUNCTION Fun_Dlmn(IDEq,l,m,n) result(Dlmn)
        INTEGER,        INTENT(IN):: IDEq,l,m,n
        REAL                      :: Dlmn

        IF (IDEq.EQ.0) THEN
        !Integ sin(l*vartheta)*sin(m*vartheta)*cos(n*vartheta)dvartheta=0
        Dlmn=(Fun_KronDelta(n,abs(m-l))-Fun_KronDelta(n,m+l))
        ELSE
        !Integ sin(l*vartheta)*sin(m*vartheta)*cos(n*vartheta)dvartheta=0
        Dlmn=(Fun_KronDelta(n,abs(m-l))+Fun_KronDelta(n,m+l))
        ENDIF
        Dlmn=Dlmn*PI/Fun_epsilon_l(n)
  END FUNCTION

  FUNCTION Fun_gamma(m,beta,alpha) result(gam)
     INTEGER, INTENT(IN) :: m
     REAL,    INTENT(IN) :: beta,alpha
     COMPLEX             :: gam

     IF (0.LE.ABS(beta) .AND.ABS(beta).LE.alpha) THEN
        gam=exp(II*m*asin(ABS(beta)/alpha))/SQRT(alpha**2-beta**2)
     ELSEIF(0.LT.alpha .AND. alpha.LT.ABS(beta)) THEN
        gam=II*exp(II*m*PI/2)*(alpha/(ABS(beta)+SQRT(beta**2-alpha**2)))**m  &
              /SQRT(beta**2-alpha**2)
     ENDIF
     IF (beta<0)  gam=CONJG(gam)
  END FUNCTION


  FUNCTION Fun_epsilon_l(l) result(eps)
     INTEGER, INTENT(IN) :: l
     INTEGER             :: eps
     IF (l.EQ. 0) eps=1
     IF (l.NE. 0) eps=2
  END FUNCTION
!!---------------------------------------------------------------------
END MODULE
