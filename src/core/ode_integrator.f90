MODULE ODE_INTEGRATOR
!> Generic, physics-agnostic numerics: bisection table lookup
!> (LOCATE_TABLE) and an adaptive-stepsize Cash-Karp RK45 ODE integrator
!> (ODEINT, with internal RKQS/RKCK stages) -- ported from
!> ~/Desktop/EOSNS/src/locate.f and odeint.f (Numerical Recipes'
!> `locate`/`odeint`/`rkqs`/`rkck`, unmodified numerical content).
!> Deliberately has no knowledge of TOV physics or any other specific
!> ODE system -- TOV_SOLVER builds its own DERIVS callback on top of
!> this, matching this codebase's existing physics-blind-utility
!> layering (LINEAR_SOLVE is the direct precedent: "no knowledge of
!> spherical harmonics... A and B are just numbers to it").
!>
!> @warning This is a careful, faithful port, not a rewrite: every
!>   numerical constant (Cash-Karp tableau coefficients, safety factors,
!>   error exponents) is copied byte-for-byte from the F77 original.
!>   The only real changes are structural/modernizing, not numerical:
!>   `EXTERNAL derivs` becomes the `DERIVS_I` abstract interface (this
!>   codebase's existing pattern for passing procedures around, e.g.
!>   `REGIME_ADVANCE_I` in regime_interface.f90), and the F77 fixed-size
!>   `NMAX=80` work arrays become properly assumed-shape/automatic
!>   arrays sized from the caller's actual NVAR -- removing an arbitrary
!>   F77-era ceiling changes no numerics for any problem that fit under
!>   it before.
USE KINDS, ONLY: dp, i4
IMPLICIT NONE
PRIVATE
PUBLIC :: DERIVS_I, LOCATE_TABLE, ODEINT

!> Callback contract for an ODE right-hand side: dydx = f(x, y).
ABSTRACT INTERFACE
  SUBROUTINE DERIVS_I(X, Y, DYDX)
    IMPORT :: dp
    REAL(KIND=dp), INTENT(IN)  :: X
    REAL(KIND=dp), INTENT(IN)  :: Y(:)
    REAL(KIND=dp), INTENT(OUT) :: DYDX(:)
  END SUBROUTINE DERIVS_I
END INTERFACE

CONTAINS

!> Bisection search: returns J such that X lies between XX(J) and
!> XX(J+1), for a monotonic (increasing or decreasing) table XX.
!> J=0 or J=SIZE(XX) signals X is outside the table -- callers must
!> check (ported from locate.f exactly, including that convention).
SUBROUTINE LOCATE_TABLE(XX, X, J)
  REAL(KIND=dp),    INTENT(IN)  :: XX(:)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(OUT) :: J
  INTEGER(KIND=i4) :: N, JL, JU, JM

  N = SIZE(XX)
  JL = 0
  JU = N + 1
  DO WHILE (JU - JL > 1)
    JM = (JU + JL) / 2
    IF ((XX(N) > XX(1)) .EQV. (X > XX(JM))) THEN
      JL = JM
    ELSE
      JU = JM
    END IF
  END DO
  J = JL
END SUBROUTINE LOCATE_TABLE

!> Integrates dy/dx=DERIVS(x,y) from X1 to X2 with adaptive stepsize
!> control (Cash-Karp RK45, embedded error estimate via RKQS/RKCK),
!> starting from step size H1. YSTART holds the initial condition on
!> entry and the solution at X2 on return. NOK/NBAD count accepted-
!> first-try vs. accepted-after-shrinking steps (diagnostic only, ported
!> as-is). Ported from odeint.f's `odeint` verbatim in numerical content.
SUBROUTINE ODEINT(YSTART, X1, X2, H1, NOK, NBAD, DERIVS)
  REAL(KIND=dp),    INTENT(INOUT) :: YSTART(:)
  REAL(KIND=dp),    INTENT(IN)    :: X1, X2, H1
  INTEGER(KIND=i4), INTENT(OUT)   :: NOK, NBAD
  PROCEDURE(DERIVS_I) :: DERIVS

  INTEGER(KIND=i4), PARAMETER :: MAXSTP = 1000000
  REAL(KIND=dp),    PARAMETER :: TINY = 1.0E-20_dp
  REAL(KIND=dp),    PARAMETER :: EPS = 1.0E-10_dp, HMIN = 1.0E-24_dp

  INTEGER(KIND=i4) :: NVAR, NSTP
  REAL(KIND=dp) :: X, H, HDID, HNEXT
  REAL(KIND=dp), ALLOCATABLE :: Y(:), DYDX(:), YSCAL(:)

  NVAR = SIZE(YSTART)
  ALLOCATE(Y(NVAR), DYDX(NVAR), YSCAL(NVAR))

  X = X1
  H = SIGN(H1, X2 - X1)
  NOK = 0
  NBAD = 0
  Y = YSTART

  DO NSTP = 1, MAXSTP
    CALL DERIVS(X, Y, DYDX)
    YSCAL = ABS(Y) + ABS(H * DYDX) + TINY

    IF ((X + H - X2) * (X + H - X1) > 0.0_dp) H = X2 - X

    CALL RKQS(Y, DYDX, X, H, EPS, YSCAL, HDID, HNEXT, DERIVS)

    IF (HDID == H) THEN
      NOK = NOK + 1
    ELSE
      NBAD = NBAD + 1
    END IF

    IF ((X - X2) * (X2 - X1) >= 0.0_dp) THEN
      YSTART = Y
      DEALLOCATE(Y, DYDX, YSCAL)
      RETURN
    END IF

    IF (ABS(HNEXT) < HMIN) THEN
      WRITE(*,'(A)') 'ODEINT: stepsize smaller than minimum'
      STOP 1
    END IF
    H = HNEXT
  END DO

  WRITE(*,'(A)') 'ODEINT: too many steps'
  STOP 1
END SUBROUTINE ODEINT

!> One adaptive Cash-Karp step: tries HTRY, shrinks and retries until
!> the embedded error estimate is within EPS (scaled by YSCAL), then
!> reports what step size was actually used (HDID) and suggests the
!> next one (HNEXT). Ported from odeint.f's `rkqs` verbatim.
SUBROUTINE RKQS(Y, DYDX, X, HTRY, EPS, YSCAL, HDID, HNEXT, DERIVS)
  REAL(KIND=dp),    INTENT(INOUT) :: Y(:)
  REAL(KIND=dp),    INTENT(IN)    :: DYDX(:)
  REAL(KIND=dp),    INTENT(INOUT) :: X
  REAL(KIND=dp),    INTENT(IN)    :: HTRY, EPS, YSCAL(:)
  REAL(KIND=dp),    INTENT(OUT)   :: HDID, HNEXT
  PROCEDURE(DERIVS_I) :: DERIVS

  REAL(KIND=dp), PARAMETER :: SAFETY = 0.9_dp, PGROW = -0.2_dp, &
    PSHRNK = -0.25_dp, ERRCON = 1.89E-4_dp
  REAL(KIND=dp) :: ERRMAX, H, XNEW
  REAL(KIND=dp), ALLOCATABLE :: YERR(:), YTEMP(:)

  ALLOCATE(YERR(SIZE(Y)), YTEMP(SIZE(Y)))
  H = HTRY
  DO
    CALL RKCK(Y, DYDX, X, H, YTEMP, YERR, DERIVS)
    ERRMAX = MAXVAL(ABS(YERR / YSCAL)) / EPS
    IF (ERRMAX <= 1.0_dp) EXIT
    H = SAFETY * H * (ERRMAX**PSHRNK)
    IF (H < 0.1_dp * H) H = 0.1_dp * H   ! ported exactly; see header re: fidelity
    XNEW = X + H
    IF (XNEW == X) THEN
      WRITE(*,'(A)') 'RKQS: stepsize underflow'
      STOP 1
    END IF
  END DO

  IF (ERRMAX > ERRCON) THEN
    HNEXT = SAFETY * H * (ERRMAX**PGROW)
  ELSE
    HNEXT = 5.0_dp * H
  END IF
  HDID = H
  X = X + H
  Y = YTEMP
  DEALLOCATE(YERR, YTEMP)
END SUBROUTINE RKQS

!> One fifth-order Cash-Karp Runge-Kutta step with embedded fourth-order
!> error estimate (YERR = difference between the 5th- and 4th-order
!> solutions). Ported from odeint.f's `rkck` verbatim, including its
!> Butcher-tableau constants.
SUBROUTINE RKCK(Y, DYDX, X, H, YOUT, YERR, DERIVS)
  REAL(KIND=dp), INTENT(IN)  :: Y(:), DYDX(:), X, H
  REAL(KIND=dp), INTENT(OUT) :: YOUT(:), YERR(:)
  PROCEDURE(DERIVS_I) :: DERIVS

  REAL(KIND=dp), PARAMETER :: A2=0.2_dp, A3=0.3_dp, A4=0.6_dp, A5=1.0_dp, A6=0.875_dp
  REAL(KIND=dp), PARAMETER :: B21=0.2_dp
  REAL(KIND=dp), PARAMETER :: B31=3.0_dp/40.0_dp, B32=9.0_dp/40.0_dp
  REAL(KIND=dp), PARAMETER :: B41=0.3_dp, B42=-0.9_dp, B43=1.2_dp
  REAL(KIND=dp), PARAMETER :: B51=-11.0_dp/54.0_dp, B52=2.5_dp, &
    B53=-70.0_dp/27.0_dp, B54=35.0_dp/27.0_dp
  REAL(KIND=dp), PARAMETER :: B61=1631.0_dp/55296.0_dp, B62=175.0_dp/512.0_dp, &
    B63=575.0_dp/13824.0_dp, B64=44275.0_dp/110592.0_dp, B65=253.0_dp/4096.0_dp
  REAL(KIND=dp), PARAMETER :: C1=37.0_dp/378.0_dp, C3=250.0_dp/621.0_dp, &
    C4=125.0_dp/594.0_dp, C6=512.0_dp/1771.0_dp
  REAL(KIND=dp), PARAMETER :: DC1=C1-2825.0_dp/27648.0_dp, &
    DC3=C3-18575.0_dp/48384.0_dp, DC4=C4-13525.0_dp/55296.0_dp, &
    DC5=-277.0_dp/14336.0_dp, DC6=C6-0.25_dp

  REAL(KIND=dp), ALLOCATABLE :: AK2(:), AK3(:), AK4(:), AK5(:), AK6(:), YTEMP(:)
  INTEGER(KIND=i4) :: N

  N = SIZE(Y)
  ALLOCATE(AK2(N), AK3(N), AK4(N), AK5(N), AK6(N), YTEMP(N))

  YTEMP = Y + B21*H*DYDX
  CALL DERIVS(X+A2*H, YTEMP, AK2)
  YTEMP = Y + H*(B31*DYDX + B32*AK2)
  CALL DERIVS(X+A3*H, YTEMP, AK3)
  YTEMP = Y + H*(B41*DYDX + B42*AK2 + B43*AK3)
  CALL DERIVS(X+A4*H, YTEMP, AK4)
  YTEMP = Y + H*(B51*DYDX + B52*AK2 + B53*AK3 + B54*AK4)
  CALL DERIVS(X+A5*H, YTEMP, AK5)
  YTEMP = Y + H*(B61*DYDX + B62*AK2 + B63*AK3 + B64*AK4 + B65*AK5)
  CALL DERIVS(X+A6*H, YTEMP, AK6)

  YOUT = Y + H*(C1*DYDX + C3*AK3 + C4*AK4 + C6*AK6)
  YERR = H*(DC1*DYDX + DC3*AK3 + DC4*AK4 + DC5*AK5 + DC6*AK6)

  DEALLOCATE(AK2, AK3, AK4, AK5, AK6, YTEMP)
END SUBROUTINE RKCK

END MODULE ODE_INTEGRATOR
