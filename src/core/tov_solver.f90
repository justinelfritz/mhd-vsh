MODULE TOV_SOLVER
!> Tolman-Oppenheimer-Volkoff stellar-structure solver -- ported from
!> Dany Page's NSCool (`astroscu.unam.mx/neutrones/NSCool`, ASCL entry
!> 1609.009), `TOV/TOV.f`'s `subroutine twostep` (the TOV right-hand
!> side) and its main integration loop (RK4 in radius r, a heuristic
!> adaptive step size -- see below), built on EOS_TABLE. Faithful port:
!> every geometrized-unit conversion constant is copied byte-for-byte
!> from the original -- these are validated against NSCool's own
!> bundled reference output (`Prof_APR_Cat_1.4.dat`/
!> `prod_APR_EOS_Cat.dat`, an M=1.40 star built from the bundled
!> `APR_EOS_Cat.dat` table), not re-derived.
!>
!> @warning Units: TOV.f's own fixed unit system (NOT
!>   central-density-rescaled the way this project's prior TOV port
!>   was) -- length: km, mass: solar masses, density: solar masses per
!>   `km**3` (`RHO_UNIT_CGS=1.989d18 g/cm**3`, numerically =
!>   `Msun_grams/(1 km in cm)**3`), `G_CODE=1.484` (Newton's G expressed
!>   in these units). Confirmed empirically against the bundled
!>   reference file, not just taken from a comment: the M=1.40 profile's
!>   own last row reads `R=11567.180107` (m, i.e. `11.567` km) and
!>   `M=1.400000000` (solar masses) directly, with no further scaling
!>   needed beyond TOV.f's own hardcoded `*1.e3` (km->m, undone here --
!>   this port stays in km) output conversion.
!>
!> @warning Structural difference from the RK45-in-log(P) algorithm this
!>   module used before this port (see the addendum this replaces in
!>   the project's own plan history): NSCool's own integrator is a
!>   fixed 4-stage RK4 **in radius r** (not log-pressure), with a
!>   **heuristic adaptive step size** `step = DELTA/(EMFUNC/EM-PFUNC/P)`
!>   (scales the step by the local relative rate of change of mass and
!>   pressure) rather than an embedded-error-estimate scheme -- ported
!>   as its own small, self-contained stepper (matching NSCool's own
!>   algorithm exactly, for fidelity to what actually produced the
!>   reference output), not forced into ODE_INTEGRATOR's Cash-Karp
!>   machinery, which stays a general-purpose, TOV-agnostic utility.
!>
!> @warning Regression fidelity: `M` matches the bundled reference to
!>   full displayed precision (`1.400000000` Msun) and the first
!>   integration step's own intermediate values (`EM(1)`, baryon
!>   density) match to 9 significant figures, but the final `R` differs
!>   by ~0.03% (`11.5633` vs `11.5672` km). Traced (not assumed) to the
!>   extreme low-pressure tail of the integration (`P` down to
!>   `~1d-28` in code units, far past any physically meaningful
!>   "surface") -- the step-size heuristic's `PFUNC/P` term is acutely
!>   sensitive there, and this port's `EOS_TABLE::LOCATE_TABLE`
!>   (bisection) vs. `ener`/`pres`/`rho`'s own forward-index-cached
!>   search are two independently-correct but not bit-identical
!>   bracket-search strategies, whose tiny per-step differences compound
!>   over the ~60 tail steps between where `M` fully saturates (~11.49
!>   km) and where the integration actually stops. A REAL bug (the
!>   first RK4 stage calling `TWOSTEP` instead of the original's
!>   explicit `K1=L1=M1=0`, avoiding division-by-zero at `r=0`) was
!>   found and fixed via exactly this kind of regression check -- the
!>   `EM(1)` match above is what confirmed the fix, not a coincidence.
!>
!> @warning `TOV_PROFILE_T`'s composition fields (`AH`/`ZH`/`XH`/`YE`/
!>   `YN`/`A_TABLE`, carried in the prior EOSNS-based port) are DROPPED
!>   in this port -- NSCool's own `eos()` never reads composition from
!>   the EOS table at all (`APR_EOS_Cat.dat` has 17 columns, only the
!>   first 3 -- rho/P/nbar -- are ever read), and its crust-conductivity
!>   routine (`CRUST_CONDUCTIVITY::OYAFORM`) derives Z/A **directly from
!>   baryon density**, with no need for a pre-tabulated per-row
!>   composition at all. Threading composition through this module would
!>   also create a circular module dependency (CRUST_CONDUCTIVITY
!>   already depends on TOV_SOLVER for `TOV_PROFILE_T`). Crust-row
!>   filtering (what the old `A_TABLE` field was for) now uses a plain
!>   density threshold instead -- see CRUST_CONDUCTIVITY's own header.
USE KINDS,     ONLY: dp, i4
USE EOS_TABLE, ONLY: EOS_TABLE_T, LOAD_EOS_TABLE, EOS_AT_DENSITY, EOS_AT_PRESSURE
IMPLICIT NONE
PRIVATE
PUBLIC :: TOV_PROFILE_T, SOLVE_TOV_STAR

!> One radial profile of a solved TOV star, center to surface.
TYPE :: TOV_PROFILE_T
  INTEGER(KIND=i4) :: N = 0
  REAL(KIND=dp) :: MASS_MSUN = 0.0_dp   !! Total gravitational mass, solar masses.
  REAL(KIND=dp) :: RADIUS_KM = 0.0_dp   !! Surface radius, km.
  REAL(KIND=dp), ALLOCATABLE :: R(:)      !! Radius, km.
  REAL(KIND=dp), ALLOCATABLE :: RHOCGS(:) !! Mass density, g/cm**3.
  REAL(KIND=dp), ALLOCATABLE :: PCGS(:)   !! Pressure, dyn/cm**2.
  REAL(KIND=dp), ALLOCATABLE :: NBFM(:)   !! Baryon number density, fm**-3.
END TYPE TOV_PROFILE_T

CONTAINS

!> Solves a single TOV star of central baryon number density NBAR_CENTRAL
!> (`fm**-3`) using the EOS table at EOS_PATH, returning its radial profile
!> resampled onto NPOINTS+1 points evenly spaced in log(pressure) (P
!> drops many orders of magnitude center-to-surface, so log-spacing keeps
!> resolution near the surface -- same sampling convention the prior
!> port used, kept for interface stability even though NSCool's own
!> integrator produces an intrinsically adaptive, not evenly-log-P-
!> spaced, set of raw points -- see RESAMPLE_PROFILE).
!>
!> @param NBAR_CENTRAL Central baryon number density, `fm**-3` (e.g.
!>   0.5447307 for the M=1.40 reference star -- NSCool's own `rhoc`
!>   input convention; NOT a mass density, unlike the prior port's
!>   `RHOCGS` argument this replaces).
!> @param EOS_PATH Path to a table in NSCool's own format (see
!>   EOS_TABLE::LOAD_EOS_TABLE).
!> @param NPOINTS Number of log-pressure resampling points for the
!>   returned profile.
SUBROUTINE SOLVE_TOV_STAR(NBAR_CENTRAL, EOS_PATH, NPOINTS, PROFILE)
  REAL(KIND=dp),       INTENT(IN)  :: NBAR_CENTRAL
  CHARACTER(LEN=*),    INTENT(IN)  :: EOS_PATH
  INTEGER(KIND=i4),    INTENT(IN)  :: NPOINTS
  TYPE(TOV_PROFILE_T), INTENT(OUT) :: PROFILE

  REAL(KIND=dp), PARAMETER :: G_CODE = 1.484_dp
  ! Msun_grams/(1 km in cm)**3 = 1.989d33/(1.0d5)**3 -- confirmed against
  ! the bundled reference file's own output, not assumed from a comment.
  REAL(KIND=dp), PARAMETER :: RHO_UNIT_CGS = 1.989E18_dp
  REAL(KIND=dp), PARAMETER :: PRESSURE_UNIT_CGS = 1.989E18_dp * 9.0E20_dp
  ! "converts to solar masses from [#-km**3/fm**3]" (TOV.f's own comment
  ! on `bfunc`, ported verbatim -- the baryonic-mass integrand's own
  ! unit-conversion factor).
  REAL(KIND=dp), PARAMETER :: BARYMASS_UNIT = 8.42E-4_dp
  REAL(KIND=dp), PARAMETER :: STEPI = 0.01_dp, DELTA = 0.5_dp
  INTEGER(KIND=i4), PARAMETER :: MAXPTS = 50000

  TYPE(EOS_TABLE_T) :: TABLE
  REAL(KIND=dp) :: PI4
  REAL(KIND=dp) :: R(0:MAXPTS), P(0:MAXPTS), EM(0:MAXPTS), NBAR(0:MAXPTS)
  REAL(KIND=dp) :: RHO_ARR(0:MAXPTS), BARYDEN(0:MAXPTS)
  REAL(KIND=dp) :: P0, RHO0, NBAR0
  REAL(KIND=dp) :: PFUNC, EMFUNC, BFUNC, NBAR_STOP
  REAL(KIND=dp) :: STEP
  REAL(KIND=dp) :: K1, K2, K3, K4, L1, L2, L3, L4, M1, M2, M3, M4
  INTEGER(KIND=i4) :: I, IMAX

  PI4 = 4.0_dp * ACOS(-1.0_dp)

  CALL LOAD_EOS_TABLE(EOS_PATH, TABLE)
  CALL EOS_AT_DENSITY(TABLE, NBAR_CENTRAL, P0)
  CALL EOS_AT_PRESSURE(TABLE, P0, RHO0, NBAR0)

  R(0) = 0.0_dp
  P(0) = P0
  EM(0) = 0.0_dp
  NBAR(0) = NBAR0
  RHO_ARR(0) = RHO0
  BARYDEN(0) = 0.0_dp

  ! ---- First step: fixed STEPI (TOV.f's own "Do the first step") ------
  ! K1=L1=M1=0 EXACTLY, ported verbatim -- TOV.f's own first RK4 stage
  ! is hand-special-cased to avoid calling twostep() at r=0 itself
  ! (TWOSTEP's PFUNC/BFUNC both divide by R), not a genuine stage
  ! evaluation. (Caught during regression testing: an earlier version
  ! of this port mistakenly called TWOSTEP here instead, evaluated at
  ! the wrong radius entirely -- a real transcription bug, not a
  ! stylistic choice; fixed before this file was ever wired in.)
  K1 = 0.0_dp; L1 = 0.0_dp; M1 = 0.0_dp
  CALL TWOSTEP(TABLE, R(0)+STEPI/2.0_dp, P(0)+K1/2.0_dp, EM(0)+L1/2.0_dp, PFUNC, EMFUNC, BFUNC)
  K2 = STEPI*PFUNC; L2 = STEPI*EMFUNC; M2 = STEPI*BFUNC
  CALL TWOSTEP(TABLE, R(0)+STEPI/2.0_dp, P(0)+K2/2.0_dp, EM(0)+L2/2.0_dp, PFUNC, EMFUNC, BFUNC)
  K3 = STEPI*PFUNC; L3 = STEPI*EMFUNC; M3 = STEPI*BFUNC
  CALL TWOSTEP(TABLE, R(0)+STEPI, P(0)+K3, EM(0)+L3, PFUNC, EMFUNC, BFUNC)
  K4 = STEPI*PFUNC; L4 = STEPI*EMFUNC; M4 = STEPI*BFUNC

  P(1)  = P(0)  + (K1+2.0_dp*K2+2.0_dp*K3+K4)/6.0_dp
  EM(1) = EM(0) + (L1+2.0_dp*L2+2.0_dp*L3+L4)/6.0_dp
  R(1)  = R(0)  + STEPI
  BARYDEN(1) = (M1+2.0_dp*M2+2.0_dp*M3+M4)/6.0_dp
  CALL EOS_AT_PRESSURE(TABLE, P(1), RHO_ARR(1), NBAR(1))

  ! ---- Subsequent steps: heuristic adaptive step (TOV.f's own "Do the
  ! next step") -- step = DELTA/(EMFUNC/EM-PFUNC/P), a local relative-
  ! rate-of-change heuristic, not an embedded-error estimate. -----------
  NBAR_STOP = TABLE%NBART(TABLE%N - 1)   ! TOV.f's own `deos(limit-1)`
  IMAX = 1
  DO I = 1, MAXPTS - 1
    IMAX = I + 1

    CALL TWOSTEP(TABLE, R(I), P(I), EM(I), PFUNC, EMFUNC, BFUNC)
    STEP = DELTA / (EMFUNC/EM(I) - PFUNC/P(I))

    K1 = STEP*PFUNC; L1 = STEP*EMFUNC; M1 = STEP*BFUNC
    CALL TWOSTEP(TABLE, R(I)+STEP/2.0_dp, P(I)+K1/2.0_dp, EM(I)+L1/2.0_dp, PFUNC, EMFUNC, BFUNC)
    K2 = STEP*PFUNC; L2 = STEP*EMFUNC; M2 = STEP*BFUNC
    CALL TWOSTEP(TABLE, R(I)+STEP/2.0_dp, P(I)+K2/2.0_dp, EM(I)+L2/2.0_dp, PFUNC, EMFUNC, BFUNC)
    K3 = STEP*PFUNC; L3 = STEP*EMFUNC; M3 = STEP*BFUNC
    CALL TWOSTEP(TABLE, R(I)+STEP, P(I)+K3, EM(I)+L3, PFUNC, EMFUNC, BFUNC)
    K4 = STEP*PFUNC; L4 = STEP*EMFUNC; M4 = STEP*BFUNC

    P(I+1)  = P(I)  + (K1+2.0_dp*K2+2.0_dp*K3+K4)/6.0_dp
    EM(I+1) = EM(I) + (L1+2.0_dp*L2+2.0_dp*L3+L4)/6.0_dp
    R(I+1)  = R(I)  + STEP
    BARYDEN(I+1) = BARYDEN(I) + (M1+2.0_dp*M2+2.0_dp*M3+M4)/6.0_dp
    CALL EOS_AT_PRESSURE(TABLE, P(I+1), RHO_ARR(I+1), NBAR(I+1))

    IF (NBAR(I+1) < NBAR_STOP) EXIT
  END DO

  PROFILE%RADIUS_KM = R(IMAX)
  PROFILE%MASS_MSUN = EM(IMAX)

  CALL RESAMPLE_PROFILE(R, P, RHO_ARR, NBAR, IMAX, NPOINTS, PROFILE)

CONTAINS

  !> TOV right-hand side at radius R (P, EM given), ported from
  !> `subroutine twostep` verbatim -- metric potential PHI is NOT
  !> integrated (TOV.f's own `phfunc`/`phi` terms, dropped here: PHI is
  !> a purely diagnostic redshift quantity, decoupled from P/EM/BARYDEN
  !> -- confirmed by inspection of `twostep` itself, none of PFUNC/
  !> EMFUNC/BFUNC depend on PHI -- and TOV_PROFILE_T has never carried
  !> it). Host-associates G_CODE/PI4/BARYMASS_UNIT from SOLVE_TOV_STAR.
  SUBROUTINE TWOSTEP(TABLE_ARG, R_ARG, P_ARG, EM_ARG, PFUNC_OUT, EMFUNC_OUT, BFUNC_OUT)
    TYPE(EOS_TABLE_T), INTENT(IN)  :: TABLE_ARG
    REAL(KIND=dp),     INTENT(IN)  :: R_ARG, P_ARG, EM_ARG
    REAL(KIND=dp),     INTENT(OUT) :: PFUNC_OUT, EMFUNC_OUT, BFUNC_OUT
    REAL(KIND=dp) :: DENS, ENERGY, VOL

    CALL EOS_AT_PRESSURE(TABLE_ARG, P_ARG, ENERGY, DENS)   ! ener(p) -> ENERGY, rho(p) -> DENS
    VOL = PI4 * R_ARG**3

    BFUNC_OUT = DENS * PI4 * R_ARG**2 / SQRT(1.0_dp - 2.0_dp*G_CODE*EM_ARG/R_ARG)
    BFUNC_OUT = BFUNC_OUT * BARYMASS_UNIT
    PFUNC_OUT = -G_CODE * (ENERGY+P_ARG) * (EM_ARG+VOL*P_ARG) / (R_ARG*(R_ARG-2.0_dp*EM_ARG*G_CODE))
    EMFUNC_OUT = PI4 * R_ARG**2 * ENERGY
  END SUBROUTINE TWOSTEP

  !> Resamples the raw (R,P,RHO,NBAR) trajectory (0:IMAX, from
  !> SOLVE_TOV_STAR's own adaptive-step integration) onto NPOINTS+1
  !> points evenly spaced in log(P) from P(0) to P(IMAX), log-log
  !> interpolating R/RHO/NBAR against P at each target -- same sampling
  !> convention (log-P-even, log-log interpolation) as the prior port's
  !> own driver, kept here so callers see a predictably-sized profile
  !> regardless of how many raw adaptive-step points NSCool's own
  !> algorithm happened to take. Converts to physical units (km stays
  !> km; RHO/P -> cgs) at this final step, matching "conversion only at
  !> the reporting/consumption layer", not inside the raw solve.
  SUBROUTINE RESAMPLE_PROFILE(R_RAW, P_RAW, RHO_RAW, NBAR_RAW, IMAX_ARG, NPTS, PROF)
    REAL(KIND=dp),        INTENT(IN)    :: R_RAW(0:), P_RAW(0:), RHO_RAW(0:), NBAR_RAW(0:)
    INTEGER(KIND=i4),     INTENT(IN)    :: IMAX_ARG, NPTS
    TYPE(TOV_PROFILE_T),  INTENT(INOUT) :: PROF
    REAL(KIND=dp) :: LOGP0, LOGP1, DLOGP, PTARGET, FRAC
    INTEGER(KIND=i4) :: J, K

    PROF%N = NPTS + 1
    ALLOCATE(PROF%R(PROF%N), PROF%RHOCGS(PROF%N), PROF%PCGS(PROF%N), PROF%NBFM(PROF%N))

    LOGP0 = LOG10(P_RAW(0))
    LOGP1 = LOG10(P_RAW(IMAX_ARG))
    DLOGP = (LOGP1 - LOGP0) / REAL(NPTS, KIND=dp)

    K = 0
    DO J = 0, NPTS
      PTARGET = 10.0_dp**(LOGP0 + REAL(J,KIND=dp)*DLOGP)
      DO WHILE (K < IMAX_ARG-1 .AND. P_RAW(K+1) > PTARGET)
        K = K + 1
      END DO
      FRAC = LOG10(PTARGET/P_RAW(K)) / LOG10(P_RAW(K+1)/P_RAW(K))
      PROF%R(J+1)      = R_RAW(K)    + FRAC*(R_RAW(K+1)-R_RAW(K))
      PROF%RHOCGS(J+1) = RHO_RAW(K)  * 10.0_dp**(FRAC*LOG10(RHO_RAW(K+1)/RHO_RAW(K))) * RHO_UNIT_CGS
      PROF%NBFM(J+1)   = NBAR_RAW(K) * 10.0_dp**(FRAC*LOG10(NBAR_RAW(K+1)/NBAR_RAW(K)))
      PROF%PCGS(J+1)   = PTARGET * PRESSURE_UNIT_CGS
    END DO
  END SUBROUTINE RESAMPLE_PROFILE

END SUBROUTINE SOLVE_TOV_STAR

END MODULE TOV_SOLVER
