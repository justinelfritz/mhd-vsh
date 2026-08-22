!> Checks UNITS' derived conversion factors are internally consistent
!> with their own documented dimensional-analysis relations (energy =
!> [B]**2*[length]**3, power = energy/[time]) -- not a check on the
!> chosen numeric values themselves (10**12 G / km / Julian year are a
!> project convention, not something to "get right" against an external
!> reference), just that ENERGY_UNIT_ERG/POWER_UNIT_ERG_PER_S are what
!> the module's own header claims they are, and that none of the base
!> units were accidentally left non-positive.
PROGRAM TEST_UNITS
USE KINDS, ONLY: dp, i4
USE UNITS, ONLY: B_UNIT_GAUSS, LENGTH_UNIT_CM, TIME_UNIT_S, &
                 ENERGY_UNIT_ERG, POWER_UNIT_ERG_PER_S
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: TOL = 1.0E-9_dp
INTEGER(KIND=i4) :: N_FAIL

N_FAIL = 0
CALL CHECK("base_units_positive", &
  MIN(B_UNIT_GAUSS, LENGTH_UNIT_CM, TIME_UNIT_S), N_FAIL, POSITIVE=.TRUE.)
CALL CHECK("energy_unit_is_B2_L3", &
  ABS(ENERGY_UNIT_ERG - B_UNIT_GAUSS**2*LENGTH_UNIT_CM**3)/ENERGY_UNIT_ERG, N_FAIL)
CALL CHECK("power_unit_is_energy_over_time", &
  ABS(POWER_UNIT_ERG_PER_S - ENERGY_UNIT_ERG/TIME_UNIT_S)/POWER_UNIT_ERG_PER_S, N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_units"
END IF

CONTAINS

SUBROUTINE CHECK(NAME, ERR, N_FAIL, POSITIVE)
  CHARACTER(LEN=*),  INTENT(IN)    :: NAME
  REAL(KIND=dp),     INTENT(IN)    :: ERR
  INTEGER(KIND=i4),  INTENT(INOUT) :: N_FAIL
  LOGICAL, OPTIONAL, INTENT(IN)    :: POSITIVE
  LOGICAL :: FAILED
  IF (PRESENT(POSITIVE)) THEN
    FAILED = ERR <= 0.0_dp
  ELSE
    FAILED = ERR > TOL
  END IF
  IF (FAILED) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  err=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  err=", ERR
  END IF
END SUBROUTINE CHECK

END PROGRAM TEST_UNITS
