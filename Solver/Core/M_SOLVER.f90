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
!   - A. Babarit
!   - R. Kurnia (2022), added LU decomp and GMRES solver
!
!--------------------------------------------------------------------------------------
MODULE M_SOLVER
  IMPLICIT NONE

  PUBLIC :: GAUSSZ,LU_INVERS_MATRIX,GMRES_SOLVER,ReadTSolver,TSolver

  INTEGER, PARAMETER, PUBLIC ::  ID_GAUSS=0
  INTEGER, PARAMETER, PUBLIC ::  ID_LU=1
  INTEGER, PARAMETER, PUBLIC ::  ID_GMRES=2

  ! Definition of TYPE TSolver
  TYPE TSolver
    INTEGER       :: NP_GQ   !0= GAUSS, 1=LU, 2=GMRES
    INTEGER       :: ID      !0= GAUSS, 1=LU, 2=GMRES
    INTEGER       :: mRestart,MaxIter
    REAL          :: Tolerance
    REAL          :: eps_zmin
    CHARACTER(20) :: SNAME
  END TYPE TSolver


  CONTAINS

  !---------------------------------------------------------------------------!

  SUBROUTINE ReadTSolver(SolverOpt,wd)

   TYPE(TSolver) :: SolverOpt
   CHARACTER(LEN=*) :: wd
   CHARACTER(20) :: SOLVER_NAME
   INTEGER       :: NPGQ
   OPEN(10,file=wd//'/input_solver.txt',form='formatted',status='old')
   READ(10,*) NPGQ
   READ(10,*) SolverOpt%eps_zmin
   READ(10,*) SolverOpt%ID
   IF (SolverOpt%ID == ID_GMRES) READ(10,*) SolverOpt%mRestart, SolverOpt%Tolerance, SolverOpt%MaxIter
   CLOSE(10)
   IF (NPGQ<1 .OR. NPGQ>4) THEN
           WRITE(*,*) ''
           WRITE(*,*) 'Specify N in input_solver.txt for Gauss Quadrature with N in (1,4)'
       STOP
   END IF
   SolverOpt%NP_GQ=NPGQ**2

   IF (SolverOpt%ID == ID_GAUSS ) THEN
            SolverOpt%SNAME='GAUSS ELIMINATION'
   ELSE IF (SolverOpt%ID == ID_LU) THEN
            SolverOpt%SNAME='LU DECOMPOSITION'
   ELSE
            SolverOpt%SNAME='GMRES'
            IF (SolverOpt%mRestart < 1) THEN
                    WRITE(*,*) 'Specify in input_solver.txt, the restart parameter for GMRES, m>0 and m<Npanels'
                    STOP
            END IF
    END IF

  END SUBROUTINE

  SUBROUTINE GAUSSZ(A,N, Ainv)

    ! Input
    INTEGER,                  INTENT(IN)   :: N
    COMPLEX, DIMENSION(N, N), INTENT(IN)   :: A
    ! Output
    COMPLEX, DIMENSION(N, N), INTENT(OUT)  :: Ainv

    ! Local variables
    INTEGER :: I, J, K, L, IL
    COMPLEX :: C, P
    COMPLEX, ALLOCATABLE :: WS(:,:) ! Working matrix
    REAL, PARAMETER :: EPS=1E-20

    ! Initialize working matrix
    ALLOCATE(WS(N,2*N))
    WS(:, 1:N)     = A
    WS(:, N+1:2*N) = CMPLX(0., 0.)
    DO I = 1, N
      WS(I, N+I) = CMPLX(1., 0.)
    END DO


    ! Gauss pivot inversion
    DO J = 1, N-1
      K = J
      DO I = J+1, N
        IF ((ABS(WS(K, J))-ABS(WS(I, J))) <= 0.0) K = I
      END DO
      IF (K <= J) THEN
        DO L = J, 2*N
          C = WS(J, L)
          WS(J, L) = WS(K, L)
          WS(K, L) = C
        END DO
      ELSE
        IF (ABS(WS(J, J)) <= EPS) THEN
          WRITE(*, '(A,E16.6)') 'PIVOT INFERIEUR A ', EPS
          STOP
        END IF
      END IF
      DO K = J+1, 2*N
        P = WS(J, K)/WS(J, J)
        DO I = J+1, N
          WS(I, K) = WS(I, K) - WS(I, J)*P
        END DO
      END DO
    END DO
    IF (ABS(WS(N, N)) < EPS) THEN
      WRITE(*, '(A,E16.6)') 'PIVOT INFERIEUR A ', EPS
      STOP
    END IF
    DO IL = N+1, 2*N
      DO J = N, 1, -1
        WS(J, IL) = WS(J, IL)/WS(J, J)
        DO I = 1, J-1
          WS(I, IL) = WS(I, IL) - WS(I, J)*WS(J, IL)
        END DO
      END DO
    END DO

    ! Extract output
     Ainv(:, :) = WS(:, N+1:2*N)
     DEALLOCATE(WS)
     RETURN

  END SUBROUTINE

   SUBROUTINE LU_INVERS_MATRIX(A,N,Ainv)

    ! Input
    INTEGER,                  INTENT(IN)   :: N
    COMPLEX, DIMENSION(N, N), INTENT(IN)   :: A
    ! Output
    COMPLEX, DIMENSION(N, N), INTENT(OUT)  :: Ainv
    !Local variables
    INTEGER:: I,J,INFO,ID_DP
    INTEGER, DIMENSION(N)    :: IPIV
    COMPLEX, DIMENSION(N)    :: WORK  ! work array for LAPACK

    CALL CHECK_DOUBLE_PRECISION(ID_DP)

    !Initialization
    Ainv(:,:)=A

    ! ZGETRF is for complex double variable, use CGETRF for complex variable
    ! Note that all of variables in this codes are defined as complex/ real (single precision)
    ! but in the compile process,it is forced to be double precision with '-r8', see the makefile
    ! if the -r8 is removed, then CGETRF and CGETRI must be used.
    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.

    IF (ID_DP.EQ.1) THEN
      CALL ZGETRF(N,N,Ainv,N,IPIV,INFO) !LAPACK function
    ELSE
      CALL CGETRF(N,N,Ainv,N,IPIV,INFO) !LAPACK function
    END IF

    IF (INFO /= 0) THEN
     stop 'Matrix is numerically singular!'
    END IF

    ! ZGETRI is for complex double variable, use CGETRI for complex variable
    ! ZGETRI computes the inverse of a matrix using the LU factorization
    ! computed by ZGETRF.

    IF (ID_DP.EQ.1) THEN
    call ZGETRI(N, Ainv, N, IPIV, WORK, N, info)
    ELSE
    call CGETRI(N, Ainv, N, IPIV, WORK, N, info)
    END IF

    IF (INFO /= 0) THEN
         STOP 'Matrix inversion failed!'
    END IF

    RETURN

    END SUBROUTINE

    SUBROUTINE LU_SOLVER(A,B,SOL,M,N,NRHS)
    !INPUT
    INTEGER,                       INTENT(IN)   :: M,N,NRHS
    COMPLEX, DIMENSION(M, N),      INTENT(IN)   :: A
    COMPLEX, DIMENSION(M,NRHS),    INTENT(IN)   :: B

    !OUTPUT
    COMPLEX, DIMENSION(M,NRHS),    INTENT(OUT)   :: SOL

    !Local
    INTEGER:: I,J,INFO,ID_DP,LDA,LDB
    COMPLEX, ALLOCATABLE,DIMENSION(:,:)  :: Ainout,Binout
    INTEGER, ALLOCATABLE,DIMENSION(:)    :: IPIV

    CALL CHECK_DOUBLE_PRECISION(ID_DP)

    LDA=MAX(1,M)
    LDB=MAX(1,M)
    ALLOCATE(Ainout(LDA,N),Binout(LDB,NRHS),IPIV(MIN(M,N)))
    !Initialization
    Ainout(1:M,1:N)=A
    Binout(1:M,1:NRHS)=B

    ! ZGETRF is for complex double variable, use CGETRF for complex variable
    ! Note that all of variables in this codes are defined as complex/ real (single precision)
    ! but in the compile process,it is forced to be double precision with '-r8', see the makefile
    ! if the -r8 is removed, then CGETRF and CGETRI must be used.
    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.

    IF (ID_DP.EQ.1) THEN
      CALL ZGETRF(M,N,Ainout,LDA,IPIV,INFO) !LAPACK function
    ELSE
      CALL CGETRF(M,N,Ainout,LDA,IPIV,INFO) !LAPACK function
    END IF

    IF (INFO /= 0) THEN
     stop 'Matrix is numerically singular!'
    END IF

    ! ZGETRS is for complex double variable, use CGETRI for complex variable
    ! ZGETRS solves the lineary system using the LU factorization
    ! computed by ZGETRF.

      IF (ID_DP.EQ.1) THEN
      call ZGETRS('N',M, NRHS, Ainout, LDA, IPIV, Binout, LDB, info)
      ELSE
      call CGETRS('N',M, NRHS, Ainout, LDA, IPIV, Binout, LDB, info)
      END IF
    IF (INFO /= 0) THEN
         STOP 'Solution procedure failed!'
    END IF
    SOL(:,:)=Binout(:,:)
    DEALLOCATE(Ainout,Binout,IPIV)
    RETURN

    END SUBROUTINE


    SUBROUTINE GMRES_SOLVER(A,B,ZOL_GMRES,nD,SolverOpt)
    ! Input
    INTEGER,                    INTENT(IN)   :: nD
    COMPLEX, DIMENSION(nD, nD), INTENT(IN)   :: A
    COMPLEX, DIMENSION(nD),     INTENT(IN)   :: B
    Type(TSolver),              INTENT(IN)   :: SolverOpt
    ! Output
    COMPLEX, DIMENSION(nD),     INTENT(OUT)  :: ZOL_GMRES

    INTEGER ID_DP

    !   For GMRES variables and parameters
    integer mD  ! dimension for othonormal basis mD, and the matrix size nD x nD
    integer i,j,lda, lwork, ldstrt
    integer revcom, colx, coly, colz, nbscal
    integer irc(5), icntl(8), info(3)
    integer matvec, precondLeft, precondRight, dotProd
    parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)
    integer nout
    complex, dimension(:), allocatable :: work
    real  cntl(5), rinfo(2)
    complex ZERO, ONE
    parameter (ZERO = (0.0e0, 0.0e0), ONE = (1.0e0, 0.0e0))

    CALL CHECK_DOUBLE_PRECISION(ID_DP)

    mD=SolverOpt%mRestart

    lda=nD
    ldstrt = mD
    !lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1
    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + 2      !if Icntl(5)=0 or 1 and Icntl(8)=1
    allocate(work(lwork))
    !write(*,*) lwork

    !if (nD.gt.lda) then
    !    STOP
    !endif

    !***************************************
    ! setting up the RHS
         work(nD+1:2*nD)=B(1:nD)
    !*****************************************
    !** Reverse communication implementation
    !*

    !*******************************************************
    !** Initialize the control parameters to default value
    !*******************************************************
    !*
    IF (ID_DP.EQ.1) THEN
          call init_zgmres(icntl,cntl)
    ELSE
          call init_cgmres(icntl,cntl)
    ENDIF
    !*
    !*************************
    !*c Tune some parameters
    !*************************
    !*The tolerance
          icntl(1)=SolverOpt%Tolerance
    !* Save the convergence history on standard output
          icntl(3) = 40
    !* Maximum number of iterations
          icntl(7) =SolverOpt%MaxIter
    !*
    !* preconditioner location
          icntl(4) = 1
    !* orthogonalization scheme
          icntl(5)=0
    !* initial guess
          icntl(6) = 0
    !* residual calculation strategy at restart
          icntl(8) = 1
    !****************************************
    !*

    IF (ID_DP.EQ.1) THEN
    10     call drive_zgmres(nD,nD,mD,lwork,work,irc,icntl,cntl,info,rinfo)
                   revcom = irc(1)
                   colx   = irc(2)
                   coly   = irc(3)
                   colz   = irc(4)
                   nbscal = irc(5)
            !*
            IF (revcom.eq.matvec) then
            !* perform the matrix vector product
            !*        work(colz) <-- A * work(colx)
                    call zgemv('N',nD,nD,ONE,A,lda,work(colx),1,ZERO,work(colz),1)
                    goto 10
            !*
            ELSE IF (revcom.eq.precondLeft) then
            !* perform the left preconditioning
            !*         work(colz) <-- M^{-1} * work(colx)

                     call zcopy(nD,work(colx),1,work(colz),1)
                     call ztrsm('L','L','N','N',nD,1,ONE,A,lda,work(colz),nD)
                     goto 10
            !*
            ELSE IF (revcom.eq.precondRight) then
            !* perform the right preconditioning

                     call zcopy(nD,work(colx),1,work(colz),1)
                     call ztrsm('L','U','N','N',nD,1,ONE,A,lda,work(colz),nD)
                     goto 10
            !*
            ELSE IF (revcom.eq.dotProd) then
            !*      perform the scalar product
            !*      work(colz) <-- work(colx) work(coly)
            !*

                     call zgemv('C',nD,nbscal,ONE,work(colx),nD,work(coly),1,ZERO,work(colz),1)
                     goto 10
            END IF
    ELSE
    12     call drive_cgmres(nD,nD,mD,lwork,work,irc,icntl,cntl,info,rinfo)
                   revcom = irc(1)
                   colx   = irc(2)
                   coly   = irc(3)
                   colz   = irc(4)
                   nbscal = irc(5)
            !*
            IF (revcom.eq.matvec) then
            !* perform the matrix vector product
            !*        work(colz) <-- A * work(colx)
                    call cgemv('N',nD,nD,ONE,A,lda,work(colx),1,ZERO,work(colz),1)
                    goto 12
            !*
            ELSE IF (revcom.eq.precondLeft) then
            !* perform the left preconditioning
            !*         work(colz) <-- M^{-1} * work(colx)

                     call ccopy(nD,work(colx),1,work(colz),1)
                     call ctrsm('L','L','N','N',nD,1,ONE,A,lda,work(colz),nD)
                     goto 12
            !*
            ELSE IF (revcom.eq.precondRight) then
            !* perform the right preconditioning

                     call ccopy(nD,work(colx),1,work(colz),1)
                     call ctrsm('L','U','N','N',nD,1,ONE,A,lda,work(colz),nD)
                     goto 12
            !*
            ELSE IF (revcom.eq.dotProd) then
            !*      perform the scalar product
            !*      work(colz) <-- work(colx) work(coly)
            !*

                     call cgemv('C',nD,nbscal,ONE,work(colx),nD,work(coly),1,ZERO,work(colz),1)
                     goto 12
            END IF
    END IF
    !*
    if (info(1).eq.0) then
        write(*,*) ' Normal exit'
        write(*,*) ' Convergence after ', info(2),' iterations'
        write(*,*) ' '
       ! write(*,*) ' Backward error - preconditioned system', rinfo(1)
       ! write(*,*) ' Backward error - unpreconditioned system', rinfo(2)
       ! write(*,*) ' Solution : '
        ZOL_GMRES(:)=work(1:nD)
       ! write(*,*) ' Optimal size for workspace ', info(3)
    else if (info(1).eq.-1) then
        write(*,*) ' Bad value of n'
    else if (info(1).eq.-2) then
        write(*,*) ' Bad value of m'
    else if (info(1).eq.-3) then
        write(*,*) ' Too small workspace. '
        write(*,*) ' Minimal value should be ', info(2)
    else if (info(1).eq.-4) then
        write(*,*) ' No convergence after ', icntl(7), ' iterations'
        write(*,*) ' switched to LU decomposition solver'
        write(*,*) ' '
        CALL LU_SOLVER(A,B,ZOL_GMRES,nD,nD,1)
     else if (info(1).eq.-5) then
        write(*,*) ' Type of preconditioner not specified'
     endif
!*******************************
       deallocate(work)
    END SUBROUTINE

    SUBROUTINE CHECK_DOUBLE_PRECISION(ID_DP)
    INTEGER ID_DP
    REAL testREAL

    IF (KIND(testREAL)==KIND(1.d0)) THEN
        ID_DP=1
    ELSE
        ID_DP=0
    END IF
    !    print*,'ID_DP=', ID_DP, '      '
    RETURN
    END SUBROUTINE

END MODULE
