MODULE HALL_REGIME
!> Combined resistive+Hall regime: implicit backward-Euler diffusion
!> (DIFFUSION_REGIME, unchanged) wrapping N_SUB explicit RK4 substeps of
!> the Hall term (HALL_INDUCTION::HALL_INDUCTION_RHS), per outer
!> TIMESTEPPER step. This is the "primary interest" case -- resistive
!> and Hall effects together, not either in isolation (see
!> DIFFUSION_REGIME for pure diffusion, app/mhdvsh_hall_stability_experiment.f90
!> for pure Hall).
!>
!> Reuses DIFFUSION_STATE_T directly ({Phi,Psi} is already the right
!> shape) rather than defining a new state type -- callers USE
!> DIFFUSION_REGIME for the type and this module for HALL_INIT/
!> HALL_ADVANCE.
!>
!> Splitting: Lie (first-order), ordered [N_SUB explicit Hall RK4
!> substeps] THEN [1 implicit DIFFUSION_ADVANCE call], every outer step
!> -- diffusion always trails. This ordering is load-bearing, not
!> arbitrary: DIFFUSION_REGIME::SOLVE_FIELD unconditionally zeros the
!> boundary rows of its right-hand side before solving, and its cached
!> system matrix already encodes the true vacuum BC (Robin for Phi,
!> Dirichlet for Psi) -- so the solved boundary value depends only on
!> that matrix, never on whatever the incoming state had there. Ending
!> every outer step with a diffuse call therefore enforces the correct
!> physical BC exactly, for free, from already-tested machinery, with
!> NO new BC-on-raw-vector code needed. The explicit Hall substeps in
!> between need a cheap, numerically-stabilizing placeholder BC to keep
!> the FD stencils from being poisoned by an unconstrained boundary
!> mid-substep -- but its exact form does NOT turn out to be immaterial
!> to the final answer the way this note originally assumed. A
!> homogeneous-Dirichlet zero (matching app/mhdvsh_hall_stability_
!> experiment.f90's own placeholder, exercised there for 2000 steps on a
!> state whose boundary genuinely stayed near zero) was used here until
!> 2026-08-29, when the user's own energy-conservation investigation
!> traced a growing energy-balance residual to it: the outer boundary
!> carries a real, physically nonzero value under this regime's Robin
!> BC (this project's seed IC starts at exactly zero there, so the bug
!> was invisible early on), and by 1000yr of Hall-driven evolution that
!> boundary value had grown to ~60% of the field's own global max,
!> sitting exactly where F_HALL is ~1e7x larger than at the inner
!> boundary -- hard-zeroing it every RK4 stage, every substep, was
!> empirically confirmed (HALL_SUBSTEPS' own docstring has the numbers)
!> to pump total energy up >10x within a single outer step, silently
!> absorbed back down by the trailing diffuse call. HALL_SUBSTEPS now
!> holds the boundary FIXED at its pre-substep (last-diffusion-enforced)
!> value instead of zeroing it -- same cost, no artificial discontinuity.
!>
!> Hyperresistivity is NOT included yet (a separate, deferred piece of
!> work -- see ROADMAP.md/the Hall-regime design plan's Stage 3).
!> DT/N_SUB are fixed constants baked in at HALL_INIT time and NOT
!> touched by HALL_ADVANCE itself (it ignores its own DT/T arguments,
!> trusting only SAVED_DT_HALL, same convention as DIFFUSION_ADVANCE).
!> A caller CAN change DT between outer steps via HALL_SET_DT (added
!> 2026-08-21, matching TIMESTEPPER::REGIME_SET_DT_I) -- it just
!> re-invokes DIFFUSION_INIT (re-factorizing the implicit solve) and
!> HALL_INIT's own SAVED_DT_HALL update, using every OTHER argument
!> HALL_INIT was originally given, now cached for exactly this purpose
!> (SAVED_ETA/SAVED_ETA_PROFILE alongside the pre-existing `SAVED_*`
!> fields). HALL_COMPUTE_DT (matching REGIME_DT_I) is the paired
!> callback a driver passes to TIMESTEPPER::RUN_ADAPTIVE alongside it --
!> see app/mhdvsh_hall_adaptive.f90 for the Hall-CFL-driven dynamic-DT
!> driver built on both.
USE KINDS,             ONLY: dp, i4
USE GRID_RADIAL,       ONLY: RADIAL_GRID_T
USE RADIAL_OPERATORS,  ONLY: RADIAL_OPERATOR_T
USE FIELD_TYPES,       ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS, ONLY: HALL_COURANT_TIMESTEP
USE HALL_INDUCTION,    ONLY: HALL_INDUCTION_RHS
USE DIFFUSION_REGIME,  ONLY: DIFFUSION_STATE_T, DIFFUSION_INIT, DIFFUSION_ADVANCE
IMPLICIT NONE
PRIVATE
PUBLIC :: HALL_INIT, HALL_ADVANCE, HALL_COMPUTE_DT, HALL_SET_DT

TYPE(RADIAL_GRID_T),     SAVE :: SAVED_RGRID
TYPE(RADIAL_OPERATOR_T), SAVE :: SAVED_OPS
INTEGER(KIND=i4),        SAVE :: SAVED_LMAX, SAVED_N_SUB
REAL(KIND=dp),           SAVE :: SAVED_ETA, SAVED_F_HALL, SAVED_DT_HALL
REAL(KIND=dp), ALLOCATABLE, SAVE :: SAVED_ETA_PROFILE(:), SAVED_F_HALL_PROFILE(:)

CONTAINS

!> One-time setup: forwards to DIFFUSION_INIT unchanged (full outer DT,
!> same cached-LU-factorization machinery, untouched by the Hall term),
!> then caches everything HALL_ADVANCE's explicit substeps need that
!> REGIME_ADVANCE_I's fixed (STATE,DT,T) signature has no room for.
!>
!> @param RGRID, OPS, LMAX, ETA, DT as DIFFUSION_INIT.
!> @param F_HALL Hall prefactor (`km**2/(1e12 G)/yr`, see UNITS).
!> @param N_SUB Number of explicit Hall RK4 substeps per outer DT (so
!>   each substep has size DT/N_SUB) -- a fixed constant for the run,
!>   chosen from Stage 2's empirical stability data for a given
!>   (Lmax, F_Hall) combination, not derived automatically.
!> @param ETA_PROFILE, F_HALL_PROFILE Optional per-radial-row overrides
!>   (size RGRID%N) for ETA/F_HALL, forwarded to DIFFUSION_INIT and
!>   HALL_INDUCTION_RHS respectively (e.g. from
!>   CRUST_CONDUCTIVITY::ETA_AND_F_HALL_AT). Absent, behavior is
!>   bit-identical to the uniform-ETA/F_HALL case.
SUBROUTINE HALL_INIT(RGRID, OPS, LMAX, ETA, F_HALL, DT, N_SUB, ETA_PROFILE, F_HALL_PROFILE)
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  INTEGER(KIND=i4),        INTENT(IN) :: LMAX, N_SUB
  REAL(KIND=dp),           INTENT(IN) :: ETA, F_HALL, DT
  REAL(KIND=dp), OPTIONAL, INTENT(IN) :: ETA_PROFILE(:), F_HALL_PROFILE(:)

  IF (PRESENT(ETA_PROFILE)) THEN
    CALL DIFFUSION_INIT(RGRID, OPS, LMAX, ETA, DT, ETA_PROFILE=ETA_PROFILE)
  ELSE
    CALL DIFFUSION_INIT(RGRID, OPS, LMAX, ETA, DT)
  END IF

  SAVED_RGRID  = RGRID
  SAVED_OPS    = OPS
  SAVED_LMAX   = LMAX
  SAVED_ETA    = ETA
  SAVED_F_HALL = F_HALL
  SAVED_N_SUB  = N_SUB
  SAVED_DT_HALL = DT / REAL(N_SUB, KIND=dp)

  IF (ALLOCATED(SAVED_ETA_PROFILE)) DEALLOCATE(SAVED_ETA_PROFILE)
  IF (PRESENT(ETA_PROFILE)) THEN
    ALLOCATE(SAVED_ETA_PROFILE(SIZE(ETA_PROFILE)))
    SAVED_ETA_PROFILE = ETA_PROFILE
  END IF

  IF (ALLOCATED(SAVED_F_HALL_PROFILE)) DEALLOCATE(SAVED_F_HALL_PROFILE)
  IF (PRESENT(F_HALL_PROFILE)) THEN
    ALLOCATE(SAVED_F_HALL_PROFILE(SIZE(F_HALL_PROFILE)))
    SAVED_F_HALL_PROFILE = F_HALL_PROFILE
  END IF
END SUBROUTINE HALL_INIT

!> Matches TIMESTEPPER::REGIME_DT_I: the Hall-CFL-limited raw (safety-
!> factor-free) maximum stable step size for STATE's current field, via
!> FIELD_DIAGNOSTICS::HALL_COURANT_TIMESTEP using this module's own
!> cached RGRID/OPS/F_HALL(_PROFILE). A thin unwrap-and-forward wrapper,
!> not new physics -- see HALL_COURANT_TIMESTEP's own docstring for the
!> formula and its RMS-vs-pointwise-max caveat.
FUNCTION HALL_COMPUTE_DT(STATE) RESULT(DT_MAX)
  CLASS(*), INTENT(IN) :: STATE
  REAL(KIND=dp) :: DT_MAX
  DT_MAX = HUGE(1.0_dp)
  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    IF (ALLOCATED(SAVED_F_HALL_PROFILE)) THEN
      DT_MAX = HALL_COURANT_TIMESTEP(STATE%PHI, STATE%PSI, SAVED_OPS, SAVED_RGRID, &
        SAVED_F_HALL, F_HALL_PROFILE=SAVED_F_HALL_PROFILE)
    ELSE
      DT_MAX = HALL_COURANT_TIMESTEP(STATE%PHI, STATE%PSI, SAVED_OPS, SAVED_RGRID, SAVED_F_HALL)
    END IF
  END SELECT
END FUNCTION HALL_COMPUTE_DT

!> Matches TIMESTEPPER::REGIME_SET_DT_I: adopts a new outer step size DT
!> by re-running HALL_INIT with every OTHER argument taken from this
!> module's own cached state (SAVED_RGRID/SAVED_OPS/SAVED_LMAX/SAVED_ETA/
!> SAVED_F_HALL/SAVED_N_SUB/SAVED_ETA_PROFILE/SAVED_F_HALL_PROFILE) --
!> re-factorizes DIFFUSION_INIT's implicit solve for the new DT (the
!> expensive part) and updates SAVED_DT_HALL=DT/SAVED_N_SUB. Must only be
!> called after an initial HALL_INIT (needs its cached values); see
!> TIMESTEPPER::RUN_ADAPTIVE's own docstring re: calling this too often.
SUBROUTINE HALL_SET_DT(DT)
  REAL(KIND=dp), INTENT(IN) :: DT
  IF (ALLOCATED(SAVED_ETA_PROFILE) .AND. ALLOCATED(SAVED_F_HALL_PROFILE)) THEN
    CALL HALL_INIT(SAVED_RGRID, SAVED_OPS, SAVED_LMAX, SAVED_ETA, SAVED_F_HALL, DT, SAVED_N_SUB, &
      ETA_PROFILE=SAVED_ETA_PROFILE, F_HALL_PROFILE=SAVED_F_HALL_PROFILE)
  ELSE IF (ALLOCATED(SAVED_ETA_PROFILE)) THEN
    CALL HALL_INIT(SAVED_RGRID, SAVED_OPS, SAVED_LMAX, SAVED_ETA, SAVED_F_HALL, DT, SAVED_N_SUB, &
      ETA_PROFILE=SAVED_ETA_PROFILE)
  ELSE IF (ALLOCATED(SAVED_F_HALL_PROFILE)) THEN
    CALL HALL_INIT(SAVED_RGRID, SAVED_OPS, SAVED_LMAX, SAVED_ETA, SAVED_F_HALL, DT, SAVED_N_SUB, &
      F_HALL_PROFILE=SAVED_F_HALL_PROFILE)
  ELSE
    CALL HALL_INIT(SAVED_RGRID, SAVED_OPS, SAVED_LMAX, SAVED_ETA, SAVED_F_HALL, DT, SAVED_N_SUB)
  END IF
END SUBROUTINE HALL_SET_DT

!> Advances STATE (a DIFFUSION_STATE_T) by one outer step: N_SUB
!> explicit Hall RK4 substeps, then one implicit diffuse step. Matches
!> REGIME_INTERFACE::REGIME_ADVANCE_I; DT and T are both unused (trusts
!> HALL_INIT's baked-in values), same convention as DIFFUSION_ADVANCE.
SUBROUTINE HALL_ADVANCE(STATE, DT, T)
  CLASS(*),      INTENT(INOUT) :: STATE
  REAL(KIND=dp), INTENT(IN)    :: DT, T
  ASSOCIATE (UNUSED_DT => DT, UNUSED_T => T); END ASSOCIATE

  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    CALL HALL_SUBSTEPS(STATE)
    CALL DIFFUSION_ADVANCE(STATE, DT, T)
  END SELECT
END SUBROUTINE HALL_ADVANCE

!> N_SUB explicit classic-RK4 steps of HALL_INDUCTION_RHS at SAVED_DT_HALL,
!> holding the boundary fixed at its pre-substep value after every
!> stage/substep for numerical hygiene (see module header -- not the
!> physical BC, that's the trailing DIFFUSION_ADVANCE call's job).
!>
!> @warning Was ZERO_BOUNDARY (hard-zero the boundary), not HOLD_BOUNDARY,
!>   until 2026-08-29. Confirmed by the user's own energy-conservation
!>   investigation that this was a real bug, not harmless "numerical
!>   hygiene": the outer boundary has a Robin BC (not Dirichlet), so it
!>   legitimately carries a nonzero, physically real field value -- for
!>   this project's own seed IC that value is exactly zero at t=0 (hence
!>   zeroing looked harmless in every early check), but grows over time
!>   as Hall-driven redistribution pushes field toward the surface. By
!>   results/hall_crust_profile_1000yr/'s own t=1000yr checkpoint, the
!>   dominant mode's outer-boundary |Phi| (0.598) was already 60% of the
!>   field's global max (0.994) -- and F_HALL there is ~4e7x its value at
!>   the inner boundary. Hard-zeroing that large a value, every RK4 stage,
!>   every substep, right where F_HALL is most extreme, was empirically
!>   confirmed (a standalone dt_hall-convergence harness, isolating pure
!>   Hall substepping with no trailing diffusion) to grow total magnetic
!>   energy by 11.7x within a single outer step at t=1000yr -- a result
!>   stable across dt_hall spanning a 128x range, i.e. a real effect, not
!>   a discretization artifact. The trailing DIFFUSION_ADVANCE step (which
!>   enforces the true Robin BC) was cancelling nearly all of this every
!>   step, leaving behind the growing energy-balance residual that
!>   motivated this investigation. Holding the boundary fixed at its
!>   pre-substep (last-diffusion-enforced) value avoids the artificial
!>   discontinuity while still being just as cheap/local a placeholder as
!>   the zero it replaces.
SUBROUTINE HALL_SUBSTEPS(STATE)
  TYPE(DIFFUSION_STATE_T), INTENT(INOUT) :: STATE
  TYPE(SPECTRAL_SCALAR_T) :: K1_PHI, K1_PSI, K2_PHI, K2_PSI, K3_PHI, K3_PSI, K4_PHI, K4_PSI
  TYPE(SPECTRAL_SCALAR_T) :: TMP_PHI, TMP_PSI
  COMPLEX(KIND=dp), ALLOCATABLE :: PHI_IN(:), PHI_OUT(:), PSI_IN(:), PSI_OUT(:)
  INTEGER(KIND=i4) :: ISUB, N_R

  N_R = SAVED_RGRID%N
  ALLOCATE(PHI_IN(STATE%PHI%NLM), PHI_OUT(STATE%PHI%NLM))
  ALLOCATE(PSI_IN(STATE%PSI%NLM), PSI_OUT(STATE%PSI%NLM))
  PHI_IN = STATE%PHI%COEF(1,:); PHI_OUT = STATE%PHI%COEF(N_R,:)
  PSI_IN = STATE%PSI%COEF(1,:); PSI_OUT = STATE%PSI%COEF(N_R,:)

  DO ISUB = 1, SAVED_N_SUB
    CALL CALL_HALL_RHS(STATE%PHI, STATE%PSI, K1_PHI, K1_PSI)

    CALL ALLOC_SPECTRAL_SCALAR(TMP_PHI, N_R, SAVED_LMAX)
    CALL ALLOC_SPECTRAL_SCALAR(TMP_PSI, N_R, SAVED_LMAX)
    TMP_PHI%COEF = STATE%PHI%COEF + 0.5_dp*SAVED_DT_HALL*K1_PHI%COEF
    TMP_PSI%COEF = STATE%PSI%COEF + 0.5_dp*SAVED_DT_HALL*K1_PSI%COEF
    CALL HOLD_BOUNDARY(TMP_PHI, PHI_IN, PHI_OUT, N_R); CALL HOLD_BOUNDARY(TMP_PSI, PSI_IN, PSI_OUT, N_R)
    CALL CALL_HALL_RHS(TMP_PHI, TMP_PSI, K2_PHI, K2_PSI)

    TMP_PHI%COEF = STATE%PHI%COEF + 0.5_dp*SAVED_DT_HALL*K2_PHI%COEF
    TMP_PSI%COEF = STATE%PSI%COEF + 0.5_dp*SAVED_DT_HALL*K2_PSI%COEF
    CALL HOLD_BOUNDARY(TMP_PHI, PHI_IN, PHI_OUT, N_R); CALL HOLD_BOUNDARY(TMP_PSI, PSI_IN, PSI_OUT, N_R)
    CALL CALL_HALL_RHS(TMP_PHI, TMP_PSI, K3_PHI, K3_PSI)

    TMP_PHI%COEF = STATE%PHI%COEF + SAVED_DT_HALL*K3_PHI%COEF
    TMP_PSI%COEF = STATE%PSI%COEF + SAVED_DT_HALL*K3_PSI%COEF
    CALL HOLD_BOUNDARY(TMP_PHI, PHI_IN, PHI_OUT, N_R); CALL HOLD_BOUNDARY(TMP_PSI, PSI_IN, PSI_OUT, N_R)
    CALL CALL_HALL_RHS(TMP_PHI, TMP_PSI, K4_PHI, K4_PSI)

    STATE%PHI%COEF = STATE%PHI%COEF + (SAVED_DT_HALL/6.0_dp) * &
      (K1_PHI%COEF + 2.0_dp*K2_PHI%COEF + 2.0_dp*K3_PHI%COEF + K4_PHI%COEF)
    STATE%PSI%COEF = STATE%PSI%COEF + (SAVED_DT_HALL/6.0_dp) * &
      (K1_PSI%COEF + 2.0_dp*K2_PSI%COEF + 2.0_dp*K3_PSI%COEF + K4_PSI%COEF)
    CALL HOLD_BOUNDARY(STATE%PHI, PHI_IN, PHI_OUT, N_R); CALL HOLD_BOUNDARY(STATE%PSI, PSI_IN, PSI_OUT, N_R)
  END DO
  DEALLOCATE(PHI_IN, PHI_OUT, PSI_IN, PSI_OUT)
END SUBROUTINE HALL_SUBSTEPS

!> Forwards to HALL_INDUCTION_RHS, passing SAVED_F_HALL_PROFILE only when
!> HALL_INIT was given one (an unallocated allocatable must not be passed
!> to an OPTIONAL dummy argument -- PRESENT() would wrongly read true).
SUBROUTINE CALL_HALL_RHS(PHI, PSI, PHI_DOT, PSI_DOT)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)  :: PHI, PSI
  TYPE(SPECTRAL_SCALAR_T), INTENT(OUT) :: PHI_DOT, PSI_DOT
  IF (ALLOCATED(SAVED_F_HALL_PROFILE)) THEN
    CALL HALL_INDUCTION_RHS(PHI, PSI, SAVED_OPS, SAVED_RGRID, SAVED_F_HALL, PHI_DOT, PSI_DOT, &
      F_HALL_PROFILE=SAVED_F_HALL_PROFILE)
  ELSE
    CALL HALL_INDUCTION_RHS(PHI, PSI, SAVED_OPS, SAVED_RGRID, SAVED_F_HALL, PHI_DOT, PSI_DOT)
  END IF
END SUBROUTINE CALL_HALL_RHS

!> Pins the outer/inner radial rows of every mode's coefficient to the
!> given (pre-substep) boundary values -- HALL_INDUCTION_RHS applies no
!> BC itself (see its own docstring), this is the caller's job.
!> Numerical-hygiene placeholder only during Hall substeps (see module
!> header/HALL_SUBSTEPS' own docstring); the trailing DIFFUSION_ADVANCE
!> call establishes the true physical BC regardless of what this leaves.
!> Replaced ZERO_BOUNDARY (2026-08-29) -- see HALL_SUBSTEPS' own
!> docstring for why a hard zero was a real bug, not a harmless one.
SUBROUTINE HOLD_BOUNDARY(FIELD, BND_IN, BND_OUT, N_R)
  TYPE(SPECTRAL_SCALAR_T), INTENT(INOUT) :: FIELD
  COMPLEX(KIND=dp),        INTENT(IN)    :: BND_IN(:), BND_OUT(:)
  INTEGER(KIND=i4),        INTENT(IN)    :: N_R
  FIELD%COEF(1,:)   = BND_IN
  FIELD%COEF(N_R,:) = BND_OUT
END SUBROUTINE HOLD_BOUNDARY

END MODULE HALL_REGIME
