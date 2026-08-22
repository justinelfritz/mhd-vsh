MODULE IO_CHECKPOINT
!> Checkpoint/restart for a (Phi,Psi) run: saves everything needed to
!> resume a long-running simulation after a crash, per ROADMAP.md's own
!> design sketch (promoted to next-up 2026-08-21, after a ~2hr Hall run
!> was killed by a system crash mid-run with no way to resume it).
!>
!> Deliberately regime-agnostic (matches LINEAR_SOLVE's physics-blind
!> layering) -- takes PHI/PSI (SPECTRAL_SCALAR_T) directly, not a
!> regime's own state type (DIFFUSION_STATE_T lives in
!> src/regimes/diffusion/, outside mhdvsh_core, and every regime today
!> uses that exact {Phi,Psi} pair anyway, so there's nothing this module
!> would gain by being more generic than that).
!>
!> Only Phi/Psi's coefficient arrays plus the step/t/grid-shape metadata
!> needed to validate a resume are saved -- NOT a regime's own cached,
!> DT-dependent state (e.g. DIFFUSION_INIT's LU factorizations,
!> HALL_REGIME's saved RGRID/OPS copies). That state is cheap and
!> deterministic to rebuild from a driver's own compile-time parameters
!> via a fresh `*_INIT` call on resume, so saving it would only add
!> complexity and a second way for a checkpoint to go stale (e.g. if the
!> LAPACK/BLAS version changes between runs) for no benefit.
!>
!> @warning A checkpoint does NOT itself encode ETA/F_HALL/DT/N_SUB or
!>   any other regime-specific physics parameter -- resuming means
!>   re-running the SAME driver (same compiled-in parameters) with a
!>   `RESUME_PATH` pointing at a prior checkpoint, not resuming a
!>   different configuration. READ_CHECKPOINT returns N_R/LMAX/R_MIN/
!>   R_MAX so the caller can verify its own freshly-built RADIAL_GRID_T
!>   actually matches what produced the checkpoint before trusting the
!>   loaded Phi/Psi -- a silent mismatch there would produce silently
!>   wrong physics, not a crash.
!>
!> File format: formatted ASCII (this project's existing convention for
!> every other data file -- EOS tables, energy-budget logs -- not a new
!> binary format), full double-precision round-trip precision
!> (ES24.16). One row per (radial index, mode index) pair rather than
!> one wide row per radial index, so the format doesn't depend on NLM's
!> column count and stays simple to parse regardless of LMAX.
USE KINDS,       ONLY: dp, i4
USE GRID_RADIAL, ONLY: RADIAL_GRID_T
USE FIELD_TYPES, ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
IMPLICIT NONE
PRIVATE
PUBLIC :: WRITE_CHECKPOINT, READ_CHECKPOINT

CONTAINS

!> Writes PHI/PSI, the current step/time, and RGRID's shape (N_R, LMAX,
!> R_MIN, R_MAX -- for READ_CHECKPOINT's caller to validate on resume)
!> to PATH, overwriting any existing file there (a rolling single
!> checkpoint, not one file per call -- callers wanting a history of
!> checkpoints should vary PATH themselves, e.g. by step number).
!>
!> @param PATH Output file path.
!> @param ISTEP Current step number (whatever the caller's own step
!>   counter means -- this module has no opinion on it).
!> @param T Current simulation time.
!> @param PHI, PSI Current state.
!> @param RGRID Radial grid PHI/PSI were built on.
SUBROUTINE WRITE_CHECKPOINT(PATH, ISTEP, T, PHI, PSI, RGRID)
  CHARACTER(LEN=*),         INTENT(IN) :: PATH
  INTEGER(KIND=i4),         INTENT(IN) :: ISTEP
  REAL(KIND=dp),            INTENT(IN) :: T
  TYPE(SPECTRAL_SCALAR_T),  INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_GRID_T),      INTENT(IN) :: RGRID
  INTEGER(KIND=i4), PARAMETER :: UNIT = 51
  INTEGER(KIND=i4) :: IR, IDX

  OPEN(UNIT=UNIT, FILE=TRIM(PATH), STATUS='REPLACE', ACTION='WRITE')
  WRITE(UNIT,'(A)') '# MHD-VSH checkpoint'
  WRITE(UNIT,'(A)') '# n_r  lmax  r_min  r_max  step  t'
  WRITE(UNIT,'(2I8,2ES24.16,I8,ES24.16)') RGRID%N, PHI%LMAX, RGRID%R_MIN, RGRID%R_MAX, ISTEP, T
  WRITE(UNIT,'(A)') '# ir  idx  re(phi)  im(phi)  re(psi)  im(psi)'
  DO IR = 1, RGRID%N
    DO IDX = 1, PHI%NLM
      WRITE(UNIT,'(2I8,4ES24.16)') IR, IDX, &
        REAL(PHI%COEF(IR,IDX), KIND=dp), AIMAG(PHI%COEF(IR,IDX)), &
        REAL(PSI%COEF(IR,IDX), KIND=dp), AIMAG(PSI%COEF(IR,IDX))
    END DO
  END DO
  CLOSE(UNIT)
END SUBROUTINE WRITE_CHECKPOINT

!> Reads back a WRITE_CHECKPOINT file: allocates and fills PHI/PSI (via
!> ALLOC_SPECTRAL_SCALAR, sized from the checkpoint's own saved N_R/
!> LMAX -- NOT the caller's), and returns ISTEP/T plus the saved grid
!> shape (N_R, LMAX, R_MIN, R_MAX) for the caller to check against its
!> own freshly-built RGRID before trusting the result (see module
!> header's @warning -- this routine itself does NOT validate that,
!> since it has no independent RGRID of its own to compare against).
!>
!> @param PATH Checkpoint file path (from a prior WRITE_CHECKPOINT call).
!> Returns: ISTEP, T, PHI, PSI (freshly allocated), N_R, LMAX, R_MIN,
!>   R_MAX (the grid shape the checkpoint was written under).
SUBROUTINE READ_CHECKPOINT(PATH, ISTEP, T, PHI, PSI, N_R, LMAX, R_MIN, R_MAX)
  CHARACTER(LEN=*),         INTENT(IN)  :: PATH
  INTEGER(KIND=i4),         INTENT(OUT) :: ISTEP
  REAL(KIND=dp),            INTENT(OUT) :: T
  TYPE(SPECTRAL_SCALAR_T),  INTENT(OUT) :: PHI, PSI
  INTEGER(KIND=i4),         INTENT(OUT) :: N_R, LMAX
  REAL(KIND=dp),            INTENT(OUT) :: R_MIN, R_MAX
  INTEGER(KIND=i4), PARAMETER :: UNIT = 51
  CHARACTER(LEN=256) :: LINE
  INTEGER(KIND=i4) :: I, NROWS, IR, IDX
  REAL(KIND=dp) :: RE_PHI, IM_PHI, RE_PSI, IM_PSI

  OPEN(UNIT=UNIT, FILE=TRIM(PATH), STATUS='OLD', ACTION='READ')
  READ(UNIT,'(A)') LINE   ! '# MHD-VSH checkpoint'
  READ(UNIT,'(A)') LINE   ! '# n_r  lmax  r_min  r_max  step  t'
  READ(UNIT,*) N_R, LMAX, R_MIN, R_MAX, ISTEP, T
  READ(UNIT,'(A)') LINE   ! '# ir  idx  re(phi)  im(phi)  re(psi)  im(psi)'

  CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(PSI, N_R, LMAX)

  NROWS = N_R * PHI%NLM
  DO I = 1, NROWS
    READ(UNIT,*) IR, IDX, RE_PHI, IM_PHI, RE_PSI, IM_PSI
    PHI%COEF(IR,IDX) = CMPLX(RE_PHI, IM_PHI, KIND=dp)
    PSI%COEF(IR,IDX) = CMPLX(RE_PSI, IM_PSI, KIND=dp)
  END DO
  CLOSE(UNIT)
END SUBROUTINE READ_CHECKPOINT

END MODULE IO_CHECKPOINT
