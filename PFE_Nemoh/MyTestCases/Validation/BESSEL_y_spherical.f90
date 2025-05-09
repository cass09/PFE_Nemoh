PROGRAM TEST_BESSEL_Y
  IMPLICIT NONE
  INTEGER :: N, NM, K
  REAL(8) :: X
  REAL(8), ALLOCATABLE :: SY(:), DY(:)

  ! Définition des paramètres
  N = 5    ! Ordre maximal
  X = 2.0D0  ! Argument x

  ! Allocation des tableaux dynamiques
  ALLOCATE(SY(0:N), DY(0:N))

  ! Appel de la sous-routine
  CALL SPHY(N, X, NM, SY, DY)

  ! Affichage des résultats
  PRINT *, "Spherical Bessel function yn(x):"
  DO K = 0, NM
     PRINT *, "y(", K, ") =", SY(K), "   ", 
  END DO

!   PRINT *, "Derivatives y_n'(x):"
!   DO K = 0, NM
!      PRINT *, "y'(", K, ") =", DY(K)
!   END DO

  ! Libération de la mémoire
  DEALLOCATE(SY, DY)

CONTAINS

  SUBROUTINE SPHY(N, X, NM, SY, DY)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(OUT) :: NM
    REAL(8), INTENT(IN) :: X
    REAL(8), DIMENSION(0:N), INTENT(OUT) :: SY, DY
    INTEGER :: K
    REAL(8) :: F0, F1, F

    NM = N
    IF (X < 1.0D-60) THEN
       DO K = 0, N
          SY(K) = -1.0D+300
          DY(K) = 1.0D+300
       END DO
       RETURN
    END IF

    SY(0) = -DCOS(X) / X
    F0 = SY(0)
    DY(0) = (DSIN(X) + DCOS(X) / X) / X

    IF (N < 1) RETURN

    SY(1) = (SY(0) - DSIN(X)) / X
    F1 = SY(1)

    DO K = 2, N
       F = (2.0D0 * K - 1.0D0) * F1 / X - F0
       SY(K) = F
       IF (DABS(F) >= 1.0D+300) EXIT
       F0 = F1
       F1 = F
    END DO

    NM = K - 1

    DO K = 1, NM
       DY(K) = SY(K-1) - (K+1.0D0) * SY(K) / X
    END DO

  END SUBROUTINE SPHY

END PROGRAM TEST_BESSEL_Y
