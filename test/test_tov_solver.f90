!> Regression test for TOV_SOLVER::SOLVE_TOV_STAR against Dany Page's
!> NSCool's own bundled reference output (`~NSCool/TOV/Profile/
!> Prof_APR_Cat_1.4.dat`, `TOV/Production/prod_APR_EOS_Cat.dat`, an
!> M=1.40 star built from the bundled `APR_EOS_Cat.dat`, now
!> `data/eos/APR_EOS_Cat.dat` in this repo) -- an independent check
!> grounded in genuinely prior, non-self-referential results, not a
!> self-consistency check against this port's own code.
!>
!> @warning Radius tolerance is looser than mass: `M` matches the
!>   reference to full displayed precision, but `R` differs by ~0.03%
!>   -- traced (not assumed) to the integration's extreme low-pressure
!>   tail, where the step-size heuristic's own sensitivity amplifies
!>   tiny bracket-search differences between this port's bisection
!>   (EOS_TABLE::LOCATE_TABLE) and the original's forward-cached-index
!>   search. See TOV_SOLVER's own module header for the full writeup,
!>   including the real transcription bug (a mishandled first RK4
!>   stage) this same regression check caught and confirmed fixed.
PROGRAM TEST_TOV_SOLVER
USE KINDS,      ONLY: dp, i4
USE TOV_SOLVER, ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: NBAR_CENTRAL = 0.5447307_dp   ! fm**-3, M=1.40 reference star
INTEGER(KIND=i4), PARAMETER :: NPOINTS = 414
REAL(KIND=dp), PARAMETER :: TOL_TIGHT = 1.0E-3_dp
REAL(KIND=dp), PARAMETER :: TOL_RADIUS = 5.0E-4_dp
REAL(KIND=dp), PARAMETER :: RHOL_CGS = 2.2E14_dp

TYPE(TOV_PROFILE_T) :: PROFILE
INTEGER(KIND=i4) :: N_FAIL, I, I_FIRST_CRUST

N_FAIL = 0
CALL SOLVE_TOV_STAR(NBAR_CENTRAL, 'data/eos/APR_EOS_Cat.dat', NPOINTS, PROFILE)

CALL CHECK_CLOSE("radius_km", PROFILE%RADIUS_KM, 11.567180107_dp, TOL_RADIUS, N_FAIL)
CALL CHECK_CLOSE("mass_msun", PROFILE%MASS_MSUN, 1.400000000_dp, TOL_TIGHT,  N_FAIL)

! Center row (index 1) should reproduce the reference's own row-0
! values: central density 9.925265d14 g/cm**3, central pressure
! 1.45523d35 dyn/cm**2 (Prof_APR_Cat_1.4.dat's own first data row).
CALL CHECK_CLOSE("central_rhocgs", PROFILE%RHOCGS(1), 9.925265E14_dp, TOL_TIGHT, N_FAIL)
CALL CHECK_CLOSE("central_pcgs",   PROFILE%PCGS(1),   1.45523E35_dp,  TOL_TIGHT, N_FAIL)
CALL CHECK_CLOSE("central_nbfm",   PROFILE%NBFM(1),   NBAR_CENTRAL,   TOL_TIGHT, N_FAIL)

! Basic physical sanity: density should decrease monotonically from
! center to surface, and the last row should be at/near the surface.
CALL CHECK_TRUE("density_decreases_outward", &
  ALL(PROFILE%RHOCGS(2:PROFILE%N) <= PROFILE%RHOCGS(1:PROFILE%N-1)), N_FAIL)
CALL CHECK_CLOSE("last_row_at_surface", PROFILE%R(PROFILE%N), PROFILE%RADIUS_KM, TOL_TIGHT, N_FAIL)

! Crust extent sanity: some rows should sit below the core-crust
! boundary density (matching app/mhdvsh_tov.f90's own gate).
I_FIRST_CRUST = 0
DO I = 1, PROFILE%N
  IF (PROFILE%RHOCGS(I) <= RHOL_CGS) THEN
    I_FIRST_CRUST = I
    EXIT
  END IF
END DO
CALL CHECK_TRUE("crust_region_found", I_FIRST_CRUST > 0, N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_tov_solver"
END IF

CONTAINS

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

END PROGRAM TEST_TOV_SOLVER
