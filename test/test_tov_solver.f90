!> Regression test for TOV_SOLVER::SOLVE_TOV_STAR against the real,
!> already-solved 2015 reference output (~/Desktop/EOSNS/fort.34,
!> PL.DAT, the M=1.40 star at rhocgs=9.88d14 g/cm**3) -- an
!> independent check grounded in genuinely prior, non-self-referential
!> results, not a self-consistency check against this port's own code.
!>
!> Reference values below were read directly off fort.34's first row
!> (`10.8033325018  0.25088E-02  0.16233E+03  0.14627E+03  0.12200E+15`
!> = rad_km, n_e_fm-3, meff, tau, rhocgs) and PL.DAT's/nstot.f's own
!> printed M/R summary for this star (Radius=11.698..., M=1.40088
!> Msun) -- confirmed by rebuilding the original, unmodified F77 source
!> and re-deriving these same numbers directly (see this port's own
!> commit history/plan notes), not copied blind from file contents.
PROGRAM TEST_TOV_SOLVER
USE KINDS,      ONLY: dp, i4
USE TOV_SOLVER, ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: RHOCGS  = 9.88E14_dp
INTEGER(KIND=i4), PARAMETER :: NPOINTS = 414
REAL(KIND=dp), PARAMETER :: TOL_TIGHT = 1.0E-3_dp
REAL(KIND=dp), PARAMETER :: TOL_LOOSE = 2.0E-2_dp

TYPE(TOV_PROFILE_T) :: PROFILE
INTEGER(KIND=i4) :: N_FAIL, I, I_FIRST
REAL(KIND=dp) :: NEL_CGS_FIRST, NEL_FM3_REF

N_FAIL = 0
CALL SOLVE_TOV_STAR(RHOCGS, 'data/eos/lowd-eos.ja.tab', NPOINTS, PROFILE)

CALL CHECK_CLOSE("radius_km",       PROFILE%RADIUS_KM, 11.6982211606_dp, TOL_TIGHT, N_FAIL)
CALL CHECK_CLOSE("mass_msun",       PROFILE%MASS_MSUN, 1.40088_dp,       TOL_TIGHT, N_FAIL)

! First crust row (A_TABLE>0) should reproduce fort.34's own first row
! (10.8033325018 km, n_e=0.25088d-2 fm**-3, rhocgs=0.12200d15 g/cm**3)
! to close precision, since it sits right at the crust's own upper
! (highest-density) edge, insensitive to the two ports' differing
! NPOINTS log-pressure step count.
I_FIRST = 0
DO I = 1, PROFILE%N
  IF (PROFILE%A_TABLE(I) > 0.0_dp) THEN
    I_FIRST = I
    EXIT
  END IF
END DO
CALL CHECK_TRUE("crust_region_found", I_FIRST > 0, N_FAIL)

IF (I_FIRST > 0) THEN
  CALL CHECK_CLOSE("crust_first_row_r_km", PROFILE%R(I_FIRST), 10.8033325018_dp, TOL_TIGHT, N_FAIL)
  CALL CHECK_CLOSE("crust_first_row_rhocgs", PROFILE%RHOCGS(I_FIRST), 1.2200E14_dp, TOL_LOOSE, N_FAIL)

  NEL_FM3_REF = 0.25088E-02_dp
  NEL_CGS_FIRST = PROFILE%NEL(I_FIRST) * 1.0E39_dp   ! fm**-3 -> cm**-3, matching CRUST_CONDUCTIVITY's own conversion
  CALL CHECK_CLOSE("crust_first_row_n_e_cm3", NEL_CGS_FIRST, NEL_FM3_REF*1.0E39_dp, TOL_LOOSE, N_FAIL)
END IF

! Basic physical sanity: density should decrease monotonically from
! center to surface, and the last row should be at/near the surface.
CALL CHECK_TRUE("density_decreases_outward", &
  ALL(PROFILE%RHOCGS(2:PROFILE%N) <= PROFILE%RHOCGS(1:PROFILE%N-1)), N_FAIL)
CALL CHECK_CLOSE("last_row_at_surface", PROFILE%R(PROFILE%N), PROFILE%RADIUS_KM, TOL_TIGHT, N_FAIL)

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
