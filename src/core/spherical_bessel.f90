!> Spherical Bessel functions of the first and second kind, j_n(x) and
!> y_n(x) (spherical Neumann), for n=0..NMAX at a single scalar x>0.
!> Physics-blind special-function utility -- no initial-condition/mode
!> semantics here, see INITIAL_CONDITIONS for the Riccati-Bessel
!> combination (a_n*x*j_n(x)+b_n*x*y_n(x)) built on top of this.
!>
!> @note The two kinds need OPPOSITE recursion directions for numerical
!>   stability (Abramowitz & Stegun; Numerical Recipes Sec 6.7): j_n
!>   DECAYS with n (for fixed x), so naive upward recursion amplifies
!>   rounding error catastrophically once n exceeds ~x -- this module
!>   uses Miller's algorithm instead (recurse DOWNWARD from an order
!>   well above NMAX using an arbitrary seed value, the stable direction
!>   since j_n GROWS as n decreases, then rescale the whole sequence
!>   against the exact closed form j_0(x)=sin(x)/x). y_n GROWS with n,
!>   so plain upward recursion from closed-form y_0/y_1 is already
!>   stable and used directly (matches Numerical Recipes' own sphbes).
!> @warning x<=0 is invalid: y_n diverges at the origin and every
!>   formula here divides by x. This project's own crust-shell radial
!>   domain never includes r=0 (R_min>0, grid_radial.f90's
!>   FULL_SPHERE=.FALSE. shell construction), so this is not expected
!>   to bite in practice; callers passing x<=0 get a hard STOP.
MODULE SPHERICAL_BESSEL
USE KINDS, ONLY: dp, i4
IMPLICIT NONE
PRIVATE
PUBLIC :: SPHERICAL_BESSEL_JN, SPHERICAL_BESSEL_YN, SPHERICAL_BESSEL_JN_YN

CONTAINS

!> j_n(x) alone, n=0..NMAX. Thin wrapper around SPHERICAL_BESSEL_JN_YN
!> for callers that don't need y_n.
SUBROUTINE SPHERICAL_BESSEL_JN(X, NMAX, J)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(IN)  :: NMAX
  REAL(KIND=dp),    INTENT(OUT) :: J(0:NMAX)
  REAL(KIND=dp) :: Y(0:NMAX)
  CALL SPHERICAL_BESSEL_JN_YN(X, NMAX, J, Y)
END SUBROUTINE SPHERICAL_BESSEL_JN

!> y_n(x) alone, n=0..NMAX. Thin wrapper, see SPHERICAL_BESSEL_JN.
SUBROUTINE SPHERICAL_BESSEL_YN(X, NMAX, Y)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(IN)  :: NMAX
  REAL(KIND=dp),    INTENT(OUT) :: Y(0:NMAX)
  REAL(KIND=dp) :: J(0:NMAX)
  CALL SPHERICAL_BESSEL_JN_YN(X, NMAX, J, Y)
END SUBROUTINE SPHERICAL_BESSEL_YN

!> Both j_n(x) and y_n(x), n=0..NMAX, in one call -- shares the
!> sin(x)/cos(x) evaluation, and is what callers building a
!> Riccati-Bessel combination (a_n*x*j_n+b_n*x*y_n) actually need.
!> @param DJ, DY Optional derivatives f_n'(x), n=0..NMAX, via the
!>   identity f_n'(x) = (n/x)*f_n(x) - f_{n+1}(x) (verified by hand at
!>   n=0 while writing this module: reduces to j_0'(x)=-j_1(x), which
!>   matches d/dx[sin(x)/x] exactly) -- computed from one extra internal
!>   order (NMAX+1), not a separate recursion.
SUBROUTINE SPHERICAL_BESSEL_JN_YN(X, NMAX, J, Y, DJ, DY)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(IN)  :: NMAX
  REAL(KIND=dp),    INTENT(OUT) :: J(0:NMAX), Y(0:NMAX)
  REAL(KIND=dp),    INTENT(OUT), OPTIONAL :: DJ(0:NMAX), DY(0:NMAX)

  INTEGER(KIND=i4) :: NEED, N
  REAL(KIND=dp), ALLOCATABLE :: JFULL(:), YFULL(:)

  IF (X <= 0.0_dp) THEN
    WRITE(*,'(A,ES16.8)') 'SPHERICAL_BESSEL_JN_YN: x must be > 0, got ', X
    STOP 1
  END IF

  NEED = NMAX
  IF (PRESENT(DJ) .OR. PRESENT(DY)) NEED = NMAX + 1_i4

  ALLOCATE(JFULL(0:NEED), YFULL(0:NEED))
  CALL YN_UPWARD(X, NEED, YFULL)
  CALL JN_MILLER(X, NEED, JFULL)

  J = JFULL(0:NMAX)
  Y = YFULL(0:NMAX)

  IF (PRESENT(DJ)) THEN
    DO N = 0, NMAX
      DJ(N) = (REAL(N,KIND=dp)/X)*JFULL(N) - JFULL(N+1)
    END DO
  END IF
  IF (PRESENT(DY)) THEN
    DO N = 0, NMAX
      DY(N) = (REAL(N,KIND=dp)/X)*YFULL(N) - YFULL(N+1)
    END DO
  END IF

  DEALLOCATE(JFULL, YFULL)
END SUBROUTINE SPHERICAL_BESSEL_JN_YN

!> y_n(x), n=0..NEED, via direct upward recursion -- stable, since y_n
!> GROWS with n (opposite of j_n). Seeded from the exact closed forms
!> y_0(x)=-cos(x)/x, y_1(x)=-cos(x)/x**2-sin(x)/x.
SUBROUTINE YN_UPWARD(X, NEED, Y)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(IN)  :: NEED
  REAL(KIND=dp),    INTENT(OUT) :: Y(0:NEED)
  INTEGER(KIND=i4) :: N

  Y(0) = -COS(X)/X
  IF (NEED >= 1_i4) THEN
    Y(1) = -COS(X)/X**2 - SIN(X)/X
    DO N = 1, NEED-1
      Y(N+1) = (REAL(2*N+1,KIND=dp)/X)*Y(N) - Y(N-1)
    END DO
  END IF
END SUBROUTINE YN_UPWARD

!> j_n(x), n=0..NEED, via Miller's algorithm: downward recursion from
!> an order well above NEED (arbitrary seed, stable direction since j_n
!> GROWS as n decreases), then rescaled against the exact closed form
!> j_0(x)=sin(x)/x. Downward recursion is stable independent of x (the
!> instability lives only in the upward direction), so no special
!> handling is needed for large x here.
SUBROUTINE JN_MILLER(X, NEED, J)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(IN)  :: NEED
  REAL(KIND=dp),    INTENT(OUT) :: J(0:NEED)
  INTEGER(KIND=i4) :: N_START, K
  REAL(KIND=dp) :: SCALE_FACTOR
  REAL(KIND=dp), ALLOCATABLE :: RAW(:)

  ! Margin above NEED: standard Miller's-algorithm heuristic (Numerical
  ! Recipes-style sqrt(40*n) growth, floored at 20), verified
  ! empirically for n<=20 by this module's own test suite
  ! (test_spherical_bessel.f90's stability sweep) rather than trusted
  ! blindly.
  N_START = NEED + MAX(20_i4, NINT(SQRT(40.0_dp*REAL(NEED,KIND=dp)), KIND=i4))

  ALLOCATE(RAW(0:N_START+1))
  RAW(N_START+1) = 0.0_dp
  RAW(N_START)   = 1.0_dp   ! arbitrary nonzero seed
  DO K = N_START-1, 0, -1
    RAW(K) = (REAL(2*K+3,KIND=dp)/X)*RAW(K+1) - RAW(K+2)
  END DO

  SCALE_FACTOR = (SIN(X)/X) / RAW(0)
  J = RAW(0:NEED) * SCALE_FACTOR

  DEALLOCATE(RAW)
END SUBROUTINE JN_MILLER

END MODULE SPHERICAL_BESSEL
