!> Regression test for FORCE_FREE_BESSEL::SOLVE_FORCE_FREE_EIGENVALUE,
!> covering both outer-BC strategies (see the module's own header).
!>
!> FF_OUTER_BC_DOUBLE_DIRICHLET:
!>  (1) l=0 closed-form check -- substituting j_0(x)=sin(x)/x,
!>      y_0(x)=-cos(x)/x into F(mu)=0 collapses it EXACTLY to
!>      -sin(mu*(R_max-R_min))/(mu**2*R_min*R_max), whose roots are
!>      mu_n=n*pi/(R_max-R_min) -- an exact, independent check of the
!>      bracket/bisect machinery itself (l=0 is not physically meaningful
!>      for a divergence-free B-field multipole, it's purely a math check,
!>      and this same bracket/bisect code path is shared by both OUTER_BC
!>      cases, so this also validates that shared machinery);
!>  (2) eigenvalue-equation residual (F(mu_solved)=0) at several l;
!>  (3) double-Dirichlet BC satisfaction -- Phi(R_min)=Phi(R_max)=0;
!>  (4) root ladder -- ROOT_INDEX=1..4 give strictly increasing mu, each
!>      independently satisfying (2)-(3).
!>
!> FF_OUTER_BC_ARBITRARY_VALUE:
!>  (5) inner Dirichlet satisfaction -- Phi(R_min)=0 -- at several l;
!>  (6) Phi(R_max) is clearly NONZERO -- the defining difference from
!>      FF_OUTER_BC_DOUBLE_DIRICHLET (confirms this case actually allows
!>      nonzero B_r at the surface, not just that it happens not to fail);
!>  (7) independent cross-check against the paper's OWN Eq. A19 form,
!>      `tan(x_min-x_max) = -beta_l(x_min)/A_l(x_min)`, with A_l,beta_l
!>      computed via the closed-form conversion from j_l,y_l derived in
!>      this module's header (A_l(x)=j_l(x)sin(x)-y_l(x)cos(x),
!>      beta_l(x)=j_l(x)cos(x)+y_l(x)sin(x)) -- genuinely independent of
!>      EIGENVALUE_RESIDUAL's own (simplified, not-A19-shaped) formula,
!>      since it's a different algebraic rearrangement;
!>  (8) root ladder, same idea as (4).
PROGRAM TEST_FORCE_FREE_BESSEL
USE KINDS,             ONLY: dp, i4
USE SPHERICAL_BESSEL,  ONLY: SPHERICAL_BESSEL_JN_YN
USE FORCE_FREE_BESSEL, ONLY: SOLVE_FORCE_FREE_EIGENVALUE, &
                              FF_OUTER_BC_DOUBLE_DIRICHLET, FF_OUTER_BC_ARBITRARY_VALUE
IMPLICIT NONE

REAL(KIND=dp),    PARAMETER :: TOL     = 1.0E-9_dp
REAL(KIND=dp),    PARAMETER :: PI      = 3.14159265358979_dp
REAL(KIND=dp),    PARAMETER :: R_MIN   = 10.3029_dp   ! matches the live TOV-derived
REAL(KIND=dp),    PARAMETER :: R_MAX   = 11.5621_dp   ! crust extent used elsewhere

INTEGER(KIND=i4) :: N_FAIL, L, ROOT_INDEX, N
REAL(KIND=dp)    :: MU, A_COEF, B_COEF, MU_PREV, MU_EXACT
REAL(KIND=dp)    :: PHI_MIN, PHI_MAX
REAL(KIND=dp)    :: AL_MIN, BETAL_MIN, A19_RESIDUAL

N_FAIL = 0

! ======================= FF_OUTER_BC_DOUBLE_DIRICHLET =======================

! --- (1) l=0 closed-form check: mu_n = n*pi/(R_max-R_min) ---
DO N = 1, 3
  CALL SOLVE_FORCE_FREE_EIGENVALUE(0_i4, R_MIN, R_MAX, FF_OUTER_BC_DOUBLE_DIRICHLET, &
    INT(N,KIND=i4), MU, A_COEF, B_COEF)
  MU_EXACT = REAL(N,KIND=dp)*PI/(R_MAX-R_MIN)
  CALL CHECK_CLOSE("dd_l0_closed_form_root", MU, MU_EXACT, TOL, N_FAIL)
END DO

! --- (2)+(3): eigenvalue residual and BC satisfaction at l=1,2,3,5,10 ---
DO L = 1, 10
  IF (L /= 1 .AND. L /= 2 .AND. L /= 3 .AND. L /= 5 .AND. L /= 10) CYCLE
  CALL SOLVE_FORCE_FREE_EIGENVALUE(INT(L,KIND=i4), R_MIN, R_MAX, FF_OUTER_BC_DOUBLE_DIRICHLET, &
    1_i4, MU, A_COEF, B_COEF)
  CALL PHI_AT(A_COEF, B_COEF, MU, INT(L,KIND=i4), R_MIN, PHI_MIN)
  CALL PHI_AT(A_COEF, B_COEF, MU, INT(L,KIND=i4), R_MAX, PHI_MAX)
  CALL CHECK_SMALL("dd_phi_vanishes_at_r_min", PHI_MIN, TOL, N_FAIL)
  CALL CHECK_SMALL("dd_phi_vanishes_at_r_max", PHI_MAX, TOL, N_FAIL)
END DO

! --- (4) root ladder: strictly increasing mu with ROOT_INDEX, each valid ---
MU_PREV = 0.0_dp
DO ROOT_INDEX = 1, 4
  CALL SOLVE_FORCE_FREE_EIGENVALUE(2_i4, R_MIN, R_MAX, FF_OUTER_BC_DOUBLE_DIRICHLET, &
    INT(ROOT_INDEX,KIND=i4), MU, A_COEF, B_COEF)
  CALL CHECK_TRUE("dd_root_ladder_increasing", MU > MU_PREV, N_FAIL)
  CALL PHI_AT(A_COEF, B_COEF, MU, 2_i4, R_MIN, PHI_MIN)
  CALL PHI_AT(A_COEF, B_COEF, MU, 2_i4, R_MAX, PHI_MAX)
  CALL CHECK_SMALL("dd_root_ladder_bc_r_min", PHI_MIN, TOL, N_FAIL)
  CALL CHECK_SMALL("dd_root_ladder_bc_r_max", PHI_MAX, TOL, N_FAIL)
  MU_PREV = MU
END DO

! ======================= FF_OUTER_BC_ARBITRARY_VALUE ========================

! --- (5)+(6)+(7): inner Dirichlet, nonzero Phi(R_max), A19 cross-check ---
DO L = 1, 10
  IF (L /= 1 .AND. L /= 2 .AND. L /= 3 .AND. L /= 5 .AND. L /= 10) CYCLE
  CALL SOLVE_FORCE_FREE_EIGENVALUE(INT(L,KIND=i4), R_MIN, R_MAX, FF_OUTER_BC_ARBITRARY_VALUE, &
    1_i4, MU, A_COEF, B_COEF)
  CALL PHI_AT(A_COEF, B_COEF, MU, INT(L,KIND=i4), R_MIN, PHI_MIN)
  CALL PHI_AT(A_COEF, B_COEF, MU, INT(L,KIND=i4), R_MAX, PHI_MAX)
  CALL CHECK_SMALL("av_phi_vanishes_at_r_min", PHI_MIN, TOL, N_FAIL)
  CALL CHECK_TRUE("av_phi_nonzero_at_r_max", ABS(PHI_MAX) > 1.0E-3_dp, N_FAIL)

  ! Independent cross-check against the paper's own Eq. A19 form, computed
  ! via a DIFFERENT algebraic path (A_l,beta_l at R_min) than
  ! EIGENVALUE_RESIDUAL's own simplified formula.
  CALL AL_BETAL_AT(MU*R_MIN, INT(L,KIND=i4), AL_MIN, BETAL_MIN)
  A19_RESIDUAL = SIN(MU*R_MIN - MU*R_MAX)*AL_MIN + COS(MU*R_MIN - MU*R_MAX)*BETAL_MIN
  CALL CHECK_SMALL("av_matches_paper_eq_A19", A19_RESIDUAL, TOL, N_FAIL)
END DO

! --- (8) root ladder ---
MU_PREV = 0.0_dp
DO ROOT_INDEX = 1, 4
  CALL SOLVE_FORCE_FREE_EIGENVALUE(2_i4, R_MIN, R_MAX, FF_OUTER_BC_ARBITRARY_VALUE, &
    INT(ROOT_INDEX,KIND=i4), MU, A_COEF, B_COEF)
  CALL CHECK_TRUE("av_root_ladder_increasing", MU > MU_PREV, N_FAIL)
  CALL PHI_AT(A_COEF, B_COEF, MU, 2_i4, R_MIN, PHI_MIN)
  CALL CHECK_SMALL("av_root_ladder_bc_r_min", PHI_MIN, TOL, N_FAIL)
  MU_PREV = MU
END DO

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_force_free_bessel"
END IF

CONTAINS

!> Phi(r) = A_COEF*x*j_l(x) + B_COEF*x*y_l(x), x=MU*r -- the same
!> construction BUILD_FORCE_FREE_BESSEL uses, re-evaluated directly here
!> (not via INITIAL_CONDITIONS) to keep this test focused purely on
!> FORCE_FREE_BESSEL's own output.
SUBROUTINE PHI_AT(A_COEF, B_COEF, MU, L, R, PHI_VAL)
  REAL(KIND=dp),    INTENT(IN)  :: A_COEF, B_COEF, MU, R
  INTEGER(KIND=i4), INTENT(IN)  :: L
  REAL(KIND=dp),    INTENT(OUT) :: PHI_VAL
  REAL(KIND=dp), ALLOCATABLE :: J(:), Y(:)
  REAL(KIND=dp) :: X
  ALLOCATE(J(0:L), Y(0:L))
  X = MU*R
  CALL SPHERICAL_BESSEL_JN_YN(X, L, J, Y)
  PHI_VAL = A_COEF*X*J(L) + B_COEF*X*Y(L)
  DEALLOCATE(J, Y)
END SUBROUTINE PHI_AT

!> A_l(x)=j_l(x)sin(x)-y_l(x)cos(x), beta_l(x)=j_l(x)cos(x)+y_l(x)sin(x)
!> (solving the 2x2 system from the Rayleigh-formula decomposition,
!> Eqs. A15-A16 -- see FORCE_FREE_BESSEL's own header), computed here
!> independently of that module's internal use of the same identity.
SUBROUTINE AL_BETAL_AT(X, L, AL, BETAL)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(IN)  :: L
  REAL(KIND=dp),    INTENT(OUT) :: AL, BETAL
  REAL(KIND=dp), ALLOCATABLE :: J(:), Y(:)
  ALLOCATE(J(0:L), Y(0:L))
  CALL SPHERICAL_BESSEL_JN_YN(X, L, J, Y)
  AL    = J(L)*SIN(X) - Y(L)*COS(X)
  BETAL = J(L)*COS(X) + Y(L)*SIN(X)
  DEALLOCATE(J, Y)
END SUBROUTINE AL_BETAL_AT

SUBROUTINE CHECK_TRUE(NAME, COND, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  LOGICAL,          INTENT(IN)    :: COND
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  IF (.NOT. COND) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A)') "FAIL  ", NAME
  ELSE
    WRITE(*,'(A,A)') "PASS  ", NAME
  END IF
END SUBROUTINE CHECK_TRUE

SUBROUTINE CHECK_CLOSE(NAME, VAL, EXPECTED, TOL, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: VAL, EXPECTED, TOL
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: ERR
  ERR = ABS(VAL-EXPECTED)/MAX(ABS(EXPECTED), 1.0E-300_dp)
  IF (ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  relerr=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  relerr=", ERR
  END IF
END SUBROUTINE CHECK_CLOSE

!> For a quantity that should be exactly zero (e.g. Phi at a Dirichlet
!> boundary) -- an absolute, not relative, tolerance check.
SUBROUTINE CHECK_SMALL(NAME, VAL, TOL, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: VAL, TOL
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  IF (ABS(VAL) > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  abs=", ABS(VAL)
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  abs=", ABS(VAL)
  END IF
END SUBROUTINE CHECK_SMALL

END PROGRAM TEST_FORCE_FREE_BESSEL
