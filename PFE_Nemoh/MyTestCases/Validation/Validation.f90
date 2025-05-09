PROGRAM Validation

    USE MBESSEL,            ONLY:fun_BESSJ

    IMPLICIT NONE

    REAL, PARAMETER           :: PI=4.*ATAN(1.)
    COMPLEX, PARAMETER        :: II=CMPLX(0.,1.)
    INTEGER :: N
    COMPLEX :: Hankel_2, EXP_term, correction
    REAL :: X

    N = 1
    X = 1.0

    EXP_term = CMPLX(COS(PI * N), SIN(PI * N))
    correction = EXP_term - CMPLX((-1)**N, 0.0)

    Hankel_2 = -II * fun_BESSJ(N, X) * correction / SIN(PI * N)

    print*, "Jn =", fun_BESSJ(N, X)
    print*, "EXP(II * PI * N) =", EXP_term
    print*, "Correction term =", correction
    print*, "Hn =", Hankel_2


END PROGRAM