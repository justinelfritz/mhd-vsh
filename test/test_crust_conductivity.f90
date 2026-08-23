!> Regression test for CRUST_CONDUCTIVITY::CON_E_PHON_ION_GYP against a
!> standalone harness built directly from the ORIGINAL, unmodified
!> Dany Page's NSCool `Code/conductivity_crust.f` (not this port's own
!> code) -- an independent check of the transliteration itself, at
!> three hand-picked (T,rho) inputs spanning the crust density range.
!>
!> Reference values were produced by compiling a small F77 driver
!> against the untouched original source and running it directly (see
!> this port's own commit history/plan notes for the exact harness);
!> not hand-derived, not copied from this port's own output.
!>
!> @note All three cases use IFS=1 (finite-nuclear-size correction ON)
!>   -- matching ETA_AND_F_HALL_AT's own hardcoded choice. IFS=0 at
!>   some densities (e.g. rho=1d13, T=1d9) triggers a genuine numerical
!>   edge case in the ORIGINAL, unmodified routine itself (a `pause`
!>   inside its own `exp_int` from an invalid argument) -- a real,
!>   pre-existing limitation of that code path, not something this port
!>   introduces, and not exercised by ETA_AND_F_HALL_AT's own IFS=1
!>   usage, so not included here.
!>
!> @warning TOL is 5e-6, not full double-precision roundoff -- same
!>   root cause as the prior EOSNS-based port's identical finding (see
!>   that port's own history): some of `conductivity_crust.f`'s bare
!>   literals lack a `d0` suffix, so F77 parses them as single precision
!>   (~7 significant digits) before widening to double. This port writes
!>   every literal with `_dp` (full double-precision parsing throughout),
!>   reproducing a ~1e-6 relative difference from the original's own
!>   output -- confirmed to be exactly this effect (not a
!>   transliteration error) by rebuilding the original with
!>   `-fdefault-real-8`, which forces double-precision literal parsing
!>   throughout and reproduces this port's output bit-for-bit.
PROGRAM TEST_CRUST_CONDUCTIVITY
USE KINDS,              ONLY: dp, i4
USE CRUST_CONDUCTIVITY, ONLY: CON_E_PHON_ION_GYP
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: TOL = 5.0E-6_dp
INTEGER(KIND=i4) :: N_FAIL
REAL(KIND=dp) :: SIGMA, LAMBDA_TH, NU_E_S, NU_E_L, N_E

N_FAIL = 0

! CASE1: crust row near the top of the crust (rho~1.22d14, T=1d9 K).
CALL CON_E_PHON_ION_GYP(1.0E9_dp, 1.2200E14_dp, 1_i4, SIGMA, LAMBDA_TH, NU_E_S, NU_E_L, N_E)
CALL CHECK_CLOSE("case1_sigma",  SIGMA,     1.951752511215784E+24_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case1_lambda", LAMBDA_TH, 3.927810997515910E+20_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case1_nu_e_s", NU_E_S,    1.853567026411298E+18_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case1_nu_e_l", NU_E_L,    2.501462227264879E+18_dp, TOL, N_FAIL)

! CASE2: a lower-density crust row (rho~1d11), below neutron drip.
CALL CON_E_PHON_ION_GYP(1.0E9_dp, 1.0E11_dp, 1_i4, SIGMA, LAMBDA_TH, NU_E_S, NU_E_L, N_E)
CALL CHECK_CLOSE("case2_sigma",  SIGMA,     1.491551256092610E+22_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case2_lambda", LAMBDA_TH, 3.935645766217342E+18_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case2_nu_e_s", NU_E_S,    1.076424212120382E+19_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case2_nu_e_l", NU_E_L,    1.107942257798993E+19_dp, TOL, N_FAIL)

! CASE3: higher density (rho~1d13).
CALL CON_E_PHON_ION_GYP(1.0E9_dp, 1.0E13_dp, 1_i4, SIGMA, LAMBDA_TH, NU_E_S, NU_E_L, N_E)
CALL CHECK_CLOSE("case3_sigma",  SIGMA,     8.353293919545064E+22_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case3_lambda", LAMBDA_TH, 1.970396609311375E+19_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case3_nu_e_s", NU_E_S,    9.769079303102056E+18_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case3_nu_e_l", NU_E_L,    1.124784273794765E+19_dp, TOL, N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_crust_conductivity"
END IF

CONTAINS

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

END PROGRAM TEST_CRUST_CONDUCTIVITY
