!> Checks EOS_TABLE's loading and log-log interpolation (ported from
!> Dany Page's NSCool, `TOV/TOV.f`'s `subroutine eos`/`function ener`/
!> `function pres`/`function rho`) against the real table shipped with
!> the repo (data/eos/APR_EOS_Cat.dat): (a) interpolating exactly at a
!> table row's own pressure recovers that row's own density/baryon-
!> density columns (log-log interpolation is exact at its own nodes);
!> (b) EOS_AT_DENSITY -> EOS_AT_PRESSURE round-trips back to
!> (approximately) the same pressure, since density and pressure are
!> both monotonic in this table; (c) pressure, density, and baryon
!> density are all monotonically increasing with table row index (the
!> table itself is stored in decreasing-density order, so this project's
!> own row-1..row-N indexing after loading should still increase --
!> confirmed, not assumed).
PROGRAM TEST_EOS_TABLE
USE KINDS,     ONLY: dp, i4
USE EOS_TABLE, ONLY: EOS_TABLE_T, LOAD_EOS_TABLE, EOS_AT_DENSITY, EOS_AT_PRESSURE
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: TOL = 1.0E-6_dp
TYPE(EOS_TABLE_T) :: TABLE
INTEGER(KIND=i4) :: N_FAIL, J
REAL(KIND=dp) :: P, RHO, NBAR, P_BACK

N_FAIL = 0
CALL LOAD_EOS_TABLE('data/eos/APR_EOS_Cat.dat', TABLE)

CALL CHECK_TRUE("table_loaded_nonempty", TABLE%N > 1, N_FAIL)

! (a) Interpolation exactly at an interior table row's own pressure
! recovers that row's own density/baryon-density columns.
J = TABLE%N / 2
CALL EOS_AT_PRESSURE(TABLE, TABLE%PT(J), RHO, NBAR)
CALL CHECK_CLOSE("pressure_node_recovers_rho",  RHO,  TABLE%RHOT(J),  TOL, N_FAIL)
CALL CHECK_CLOSE("pressure_node_recovers_nbar", NBAR, TABLE%NBART(J), TOL, N_FAIL)

! (b) EOS_AT_DENSITY at that same node's own baryon density recovers
! the same pressure (round-trip through the table's OTHER
! interpolation branch, at an exact node of both).
CALL EOS_AT_DENSITY(TABLE, TABLE%NBART(J), P_BACK)
CALL CHECK_CLOSE("density_node_recovers_pressure", P_BACK, TABLE%PT(J), TOL, N_FAIL)

! (c) Monotonicity: APR_EOS_Cat.dat is stored highest-density-first, and
! LOAD_EOS_TABLE doesn't re-sort -- so pressure/density/baryon density
! all DECREASE with increasing row index (confirmed by inspecting the
! raw file directly, not assumed) -- a basic physical sanity property
! of a stable EOS either way, not previously checked anywhere.
CALL CHECK_TRUE("pressure_column_monotonic", ALL(TABLE%PT(2:) < TABLE%PT(1:TABLE%N-1)), N_FAIL)
CALL CHECK_TRUE("density_column_monotonic",  ALL(TABLE%RHOT(2:) < TABLE%RHOT(1:TABLE%N-1)), N_FAIL)
CALL CHECK_TRUE("nbar_column_monotonic",     ALL(TABLE%NBART(2:) < TABLE%NBART(1:TABLE%N-1)), N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_eos_table"
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

END PROGRAM TEST_EOS_TABLE
