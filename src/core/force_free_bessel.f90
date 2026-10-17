!> Force-free Bessel-Riccati eigenmode construction (ROADMAP.md programming
!> task 4), ported from Igoshev, Elfritz & Popov (2016), MNRAS 462, 3689,
!> Appendix A (arXiv:1608.08806) -- verified against the paper's actual
!> Appendix A text (fetched and independently re-derived, 2026-10-15), NOT
!> reconstructed from memory.
!>
!> The paper's force-free condition `curl(B) = mu_l * B` at a single degree
!> `l` gives (i) the toroidal potential locked directly to the poloidal one,
!> `T_l(r) = -mu_l * S_l(r)` (paper's Eq. A10; this project's
!> `Psi_l(r) = -mu_l * Phi_l(r)`, @note sign is a free chirality convention,
!> confirmed with user 2026-10-15 -- either sign is an equally valid,
!> mirror-image force-free field), and (ii) a Riccati-Bessel radial ODE for
!> `Phi_l` (already solved by SPHERICAL_BESSEL_JN_YN's `x*j_l(x)`/`x*y_l(x)`
!> combinations -- no new ODE-solving needed here, just the eigenvalue
!> `mu_l` and the `a_l`/`b_l` mix).
!>
!> **Two outer-BC strategies (selected via OUTER_BC, confirmed with user
!> 2026-10-15), since this project's real outer BC (Robin(Phi)+
!> Dirichlet(Psi)=0, `boundary_conditions.f90`) can't both be satisfied
!> exactly once `Psi=-mu*Phi` slaves them to the same function:**
!>
!>  - `FF_OUTER_BC_DOUBLE_DIRICHLET`: prioritize exact `Psi(R_out)=0`
!>    (hence `Phi(R_out)=0` too) over the Robin slope condition. Forces
!>    `B_r(R_out)=0` identically (since `B_r` ~ `Phi`) -- a special case,
!>    not the generic multipole. Cost: the exterior vacuum potential
!>    generally has a nonzero growing multipole coefficient (`A_nm`,
!>    `boundary_conditions.f90`'s header / the .tex "Boundary conditions"
!>    section), i.e. not guaranteed to decay.
!>    Eigenvalue equation (double-Dirichlet annulus, both `Phi(R_min)=0`
!>    and `Phi(R_max)=0`):
!>      `F(mu) = j_l(mu*R_max)*y_l(mu*R_min) - j_l(mu*R_min)*y_l(mu*R_max) = 0`
!>
!>  - `FF_OUTER_BC_ARBITRARY_VALUE`: the paper's OWN scheme (Eqs. A17-A19)
!>    -- normalize `Phi(R_out)` to an arbitrary nonzero value (so `B_r`
!>    at the surface is generically nonzero, the physically typical case)
!>    instead of a physical vacuum match at all; only the inner Dirichlet
!>    (`Phi(R_min)=0`) is a genuine BC. Re-derived here in closed form
!>    directly from `j_l`/`y_l` (no need for the paper's own `A_l(x)`/
!>    `beta_l(x)` polynomial decomposition, Eqs. A13-A16): writing
!>    `j_l(x)=A_l(x)sin(x)+beta_l(x)cos(x)`, `y_l(x)=-A_l(x)cos(x)+
!>    beta_l(x)sin(x)` and solving that 2x2 system for `A_l`,`beta_l`
!>    gives `beta_l(x) = j_l(x)*cos(x) + y_l(x)*sin(x)` directly from
!>    `j_l`,`y_l` -- the paper's `(a_l,b_l) ∝ (cos(x_out),sin(x_out))`
!>    choice (Eqs. A17-A18, generalized here to an arbitrary target value
!>    `V` instead of their fixed `1`) then substituted into the inner
!>    Dirichlet condition collapses to:
!>      `F(mu) = cos(mu*R_max)*j_l(mu*R_min) + sin(mu*R_max)*y_l(mu*R_min) = 0`
!>    -- numerically verified (2026-10-15, scipy cross-check) to give
!>    BIT-IDENTICAL roots to the paper's own Eq. A19
!>    (`tan[mu*(R_min-R_max)] = -beta_l(mu*R_min)/A_l(mu*R_min)`), just a
!>    simpler closed form. `V` (the target `Phi(R_out)` value) does NOT
!>    affect which `mu` are eigenvalues -- it's a pure amplitude choice,
!>    algebraically cancels out of the ratio `b/a`, so it folds into the
!>    SAME caller-supplied amplitude scale as the other outer-BC case
!>    (`SOLVE_FORCE_FREE_EIGENVALUE` always returns the `V=1`-normalized
!>    mix; the caller's own amplitude multiplies it, matching
!>    `INITIAL_CONDITIONS::BUILD_BESSEL_RICCATI`'s existing convention).
!>    Numerically confirmed (2026-10-15) to give a genuinely different
!>    eigenvalue ladder than `FF_OUTER_BC_DOUBLE_DIRICHLET` (e.g. l=2,
!>    R_min=10.3029, R_max=11.5621: fundamental mu=0.480 here vs 2.505
!>    there) -- the two strategies are not secretly equivalent.
!>
!> Both equations are oscillatory in `mu` (infinitely many positive roots,
!> a ladder of eigenvalues, loosely analogous to a particle-in-a-
!> spherical-shell problem) -- `ROOT_INDEX` selects which one (1=smallest
!> positive root, "the fundamental"). Solved by bracketing (a fine scan
!> for sign changes) then bisection, reusing SPHERICAL_BESSEL_JN_YN's
!> already-tested, stable `j_n`/`y_n` evaluator throughout -- no new
!> special-function numerics.
!>
!> @note Closed-form check (`l=0`, `FF_OUTER_BC_DOUBLE_DIRICHLET`, see
!>   test/test_force_free_bessel.f90`): substituting `j_0(x)=sin(x)/x`,
!>   `y_0(x)=-cos(x)/x` collapses that `F(mu)` to
!>   `-sin(mu*(R_max-R_min))/(mu**2*R_min*R_max)`, whose roots are EXACTLY
!>   `mu_n = n*pi/(R_max-R_min)` -- an exact, independent check of the
!>   root-finder itself (l=0 is not physically meaningful for a
!>   divergence-free B-field multipole, it's purely a math cross-check).
MODULE FORCE_FREE_BESSEL
USE KINDS,            ONLY: dp, i4
USE SPHERICAL_BESSEL, ONLY: SPHERICAL_BESSEL_JN_YN
IMPLICIT NONE
PRIVATE
PUBLIC :: SOLVE_FORCE_FREE_EIGENVALUE
PUBLIC :: FF_OUTER_BC_DOUBLE_DIRICHLET, FF_OUTER_BC_ARBITRARY_VALUE

INTEGER(KIND=i4), PARAMETER :: FF_OUTER_BC_DOUBLE_DIRICHLET = 1_i4
INTEGER(KIND=i4), PARAMETER :: FF_OUTER_BC_ARBITRARY_VALUE  = 2_i4

INTEGER(KIND=i4), PARAMETER :: N_SCAN_PER_HALFPERIOD = 8_i4
INTEGER(KIND=i4), PARAMETER :: MAX_ROOTS_SCAN         = 100000_i4
REAL(KIND=dp),    PARAMETER :: BISECT_TOL             = 1.0E-13_dp
INTEGER(KIND=i4), PARAMETER :: BISECT_MAX_ITER         = 200_i4

CONTAINS

!> Solves for the ROOT_INDEX-th positive eigenvalue `mu` of the force-free
!> annulus problem at degree `L`, under whichever outer-BC strategy
!> OUTER_BC selects (see module header), then returns the matching
!> `(A_COEF,B_COEF)` Riccati-Bessel mix, NORMALIZED (`A_COEF=1` for
!> FF_OUTER_BC_DOUBLE_DIRICHLET; the `V=1`-normalized mix for
!> FF_OUTER_BC_ARBITRARY_VALUE) -- caller applies its own amplitude scale
!> in both cases (same division of labor as
!> `INITIAL_CONDITIONS::BUILD_BESSEL_RICCATI`, which also takes raw,
!> caller-supplied `A`/`B`).
!>
!> @param L Spherical-harmonic degree, L>=0 (L=0 is unphysical for a
!>   divergence-free B-field multipole but mathematically well-posed --
!>   left unrestricted here, physical L>=1 validation is
!>   INITIAL_CONDITIONS::BUILD_FORCE_FREE_BESSEL's job, matching
!>   BUILD_BESSEL_RICCATI's own division of validation responsibility).
!> @param R_MIN,R_MAX Annulus boundaries, 0 < R_MIN < R_MAX.
!> @param OUTER_BC FF_OUTER_BC_DOUBLE_DIRICHLET or FF_OUTER_BC_ARBITRARY_VALUE.
!> @param ROOT_INDEX Which positive root of F(mu)=0 to return, 1=smallest
!>   ("fundamental"), 2=next, etc. Must be >=1.
!> @param MU Output eigenvalue (mu*R has units of radians; MU itself has
!>   units of 1/length, matching this module's R_MIN/R_MAX units).
!> @param A_COEF,B_COEF Output Riccati-Bessel mix (see above).
SUBROUTINE SOLVE_FORCE_FREE_EIGENVALUE(L, R_MIN, R_MAX, OUTER_BC, ROOT_INDEX, MU, A_COEF, B_COEF)
  INTEGER(KIND=i4), INTENT(IN)  :: L, OUTER_BC, ROOT_INDEX
  REAL(KIND=dp),    INTENT(IN)  :: R_MIN, R_MAX
  REAL(KIND=dp),    INTENT(OUT) :: MU, A_COEF, B_COEF
  REAL(KIND=dp) :: DMU, MU_LO, MU_HI, F_LO, F_HI, F_MID, MU_MID
  REAL(KIND=dp) :: JL_MIN, YL_MIN, JL_MAX, YL_MAX, X_OUT, BETA_L_OUT
  INTEGER(KIND=i4) :: N_FOUND, ISTEP, ITER

  IF (ROOT_INDEX < 1_i4) THEN
    WRITE(*,'(A,I0)') 'SOLVE_FORCE_FREE_EIGENVALUE: ROOT_INDEX must be >=1, got ', ROOT_INDEX
    STOP 1
  END IF
  IF (.NOT. (R_MIN > 0.0_dp .AND. R_MAX > R_MIN)) THEN
    WRITE(*,'(A)') 'SOLVE_FORCE_FREE_EIGENVALUE: require 0 < R_MIN < R_MAX'
    STOP 1
  END IF
  IF (OUTER_BC /= FF_OUTER_BC_DOUBLE_DIRICHLET .AND. OUTER_BC /= FF_OUTER_BC_ARBITRARY_VALUE) THEN
    WRITE(*,'(A,I0)') 'SOLVE_FORCE_FREE_EIGENVALUE: unknown OUTER_BC=', OUTER_BC
    STOP 1
  END IF

  ! F(mu) is oscillatory with characteristic period ~pi/(R_max-R_min) (exact
  ! for L=0/FF_OUTER_BC_DOUBLE_DIRICHLET, see module header; comparable
  ! otherwise -- the Riccati-Bessel phase asymptotes to a pure sinusoid for
  ! large x). Scan with several sub-steps per half-period so no root is
  ! skipped, matching this project's own "verify empirically, don't assume
  ! a heuristic" ethic (test/test_force_free_bessel.f90 checks this scan
  ! step is fine enough across several L and both OUTER_BC cases).
  DMU = 3.14159265358979_dp / (R_MAX - R_MIN) / REAL(N_SCAN_PER_HALFPERIOD, KIND=dp)

  N_FOUND = 0_i4
  MU_LO = DMU
  F_LO  = EIGENVALUE_RESIDUAL(MU_LO, L, R_MIN, R_MAX, OUTER_BC)
  DO ISTEP = 1, MAX_ROOTS_SCAN
    MU_HI = MU_LO + DMU
    F_HI  = EIGENVALUE_RESIDUAL(MU_HI, L, R_MIN, R_MAX, OUTER_BC)
    IF (F_LO*F_HI < 0.0_dp) THEN
      N_FOUND = N_FOUND + 1_i4
      IF (N_FOUND == ROOT_INDEX) THEN
        DO ITER = 1, BISECT_MAX_ITER
          MU_MID = 0.5_dp*(MU_LO+MU_HI)
          F_MID  = EIGENVALUE_RESIDUAL(MU_MID, L, R_MIN, R_MAX, OUTER_BC)
          IF (ABS(MU_HI-MU_LO) < BISECT_TOL*MU_MID) EXIT
          IF (F_LO*F_MID < 0.0_dp) THEN
            MU_HI = MU_MID; F_HI = F_MID
          ELSE
            MU_LO = MU_MID; F_LO = F_MID
          END IF
        END DO
        MU = 0.5_dp*(MU_LO+MU_HI)

        SELECT CASE (OUTER_BC)
        CASE (FF_OUTER_BC_DOUBLE_DIRICHLET)
          ! Inner Dirichlet ratio: a*j_l(mu*R_min) + b*y_l(mu*R_min) = 0.
          CALL JL_YL_AT(MU*R_MIN, L, JL_MIN, YL_MIN)
          IF (ABS(YL_MIN) < 1.0E-300_dp) THEN
            WRITE(*,'(A)') 'SOLVE_FORCE_FREE_EIGENVALUE: degenerate y_l(mu*R_min)~0 at the solved root'
            STOP 1
          END IF
          A_COEF = 1.0_dp
          B_COEF = -JL_MIN/YL_MIN
        CASE (FF_OUTER_BC_ARBITRARY_VALUE)
          ! (a,b) ~ (cos(x_out),sin(x_out)), scaled so Phi(R_out)=1 exactly
          ! (V=1 normalization -- see module header); beta_l(x_out) derived
          ! directly from j_l,y_l (no separate polynomial machinery needed).
          CALL JL_YL_AT(MU*R_MAX, L, JL_MAX, YL_MAX)
          X_OUT = MU*R_MAX
          BETA_L_OUT = JL_MAX*COS(X_OUT) + YL_MAX*SIN(X_OUT)
          IF (ABS(BETA_L_OUT) < 1.0E-300_dp) THEN
            WRITE(*,'(A)') 'SOLVE_FORCE_FREE_EIGENVALUE: degenerate beta_l(mu*R_max)~0 at the solved root'
            STOP 1
          END IF
          A_COEF = COS(X_OUT)/(X_OUT*BETA_L_OUT)
          B_COEF = SIN(X_OUT)/(X_OUT*BETA_L_OUT)
        END SELECT
        RETURN
      END IF
    END IF
    MU_LO = MU_HI; F_LO = F_HI
  END DO

  WRITE(*,'(A,I0,A,I0)') 'SOLVE_FORCE_FREE_EIGENVALUE: only found ', N_FOUND, &
    ' root(s) scanning up to MAX_ROOTS_SCAN, requested ROOT_INDEX=', ROOT_INDEX
  STOP 1
END SUBROUTINE SOLVE_FORCE_FREE_EIGENVALUE

!> Dispatches to whichever closed form OUTER_BC selects (see module
!> header for both derivations):
!>   FF_OUTER_BC_DOUBLE_DIRICHLET:
!>     F(mu) = j_l(mu*R_max)*y_l(mu*R_min) - j_l(mu*R_min)*y_l(mu*R_max)
!>   FF_OUTER_BC_ARBITRARY_VALUE:
!>     F(mu) = cos(mu*R_max)*j_l(mu*R_min) + sin(mu*R_max)*y_l(mu*R_min)
FUNCTION EIGENVALUE_RESIDUAL(MU, L, R_MIN, R_MAX, OUTER_BC) RESULT(F)
  REAL(KIND=dp),    INTENT(IN) :: MU, R_MIN, R_MAX
  INTEGER(KIND=i4), INTENT(IN) :: L, OUTER_BC
  REAL(KIND=dp) :: F, JL_MIN, YL_MIN, JL_MAX, YL_MAX, X_OUT
  CALL JL_YL_AT(MU*R_MIN, L, JL_MIN, YL_MIN)
  SELECT CASE (OUTER_BC)
  CASE (FF_OUTER_BC_DOUBLE_DIRICHLET)
    CALL JL_YL_AT(MU*R_MAX, L, JL_MAX, YL_MAX)
    F = JL_MAX*YL_MIN - JL_MIN*YL_MAX
  CASE (FF_OUTER_BC_ARBITRARY_VALUE)
    X_OUT = MU*R_MAX
    F = COS(X_OUT)*JL_MIN + SIN(X_OUT)*YL_MIN
  CASE DEFAULT
    F = 0.0_dp   ! unreachable: SOLVE_FORCE_FREE_EIGENVALUE validates OUTER_BC first
  END SELECT
END FUNCTION EIGENVALUE_RESIDUAL

!> Thin convenience wrapper: bare (not Riccati-scaled) j_l(x), y_l(x) at a
!> single degree L, via SPHERICAL_BESSEL_JN_YN's full stable 0..L array.
SUBROUTINE JL_YL_AT(X, L, JL, YL)
  REAL(KIND=dp),    INTENT(IN)  :: X
  INTEGER(KIND=i4), INTENT(IN)  :: L
  REAL(KIND=dp),    INTENT(OUT) :: JL, YL
  REAL(KIND=dp), ALLOCATABLE :: J(:), Y(:)
  ALLOCATE(J(0:L), Y(0:L))
  CALL SPHERICAL_BESSEL_JN_YN(X, L, J, Y)
  JL = J(L); YL = Y(L)
  DEALLOCATE(J, Y)
END SUBROUTINE JL_YL_AT

END MODULE FORCE_FREE_BESSEL
