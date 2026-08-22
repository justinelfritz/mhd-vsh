!> Regression test for CRUST_CONDUCTIVITY::POTEKHINC against a standalone
!> harness built directly from the ORIGINAL, unmodified
!> ~/Desktop/EOSNS/src/potekhinc.f (not this port's own code) -- an
!> independent check of the transliteration itself, at three hand-picked
!> (T,rho,B,composition) inputs spanning the crust density range and
!> exercising both the ZIMP>0 (impurity-scattering) and ZIMP=0 branches.
!>
!> Reference values were produced by compiling a small F77 driver
!> against the untouched original source and running it directly (see
!> this port's own commit history/plan notes for the exact harness);
!> not hand-derived, not copied from this port's own output.
!>
!> @warning TOL is 5e-6, not full double-precision roundoff: the
!>   original's `DATA AUM/1822.9/,AUD/15819.4/` and `DATA BOHR/137.036/`
!>   (potekhinc.f) -- and similar bare literals elsewhere in the file --
!>   have no `d0` suffix, so F77 parses them as single-precision
!>   constants (rounded to ~7 significant digits) BEFORE widening to
!>   double, a well-known Fortran literal-kind gotcha. This port instead
!>   writes every literal with `_dp` (full double-precision parsing
!>   throughout), which is more numerically correct but reproduces a
!>   ~1e-6 relative difference from the original's own output --
!>   confirmed to be exactly this effect (not a transliteration error)
!>   by rebuilding the original with `-fdefault-real-8`, which forces
!>   double-precision literal parsing throughout and reproduces this
!>   port's output bit-for-bit.
PROGRAM TEST_CRUST_CONDUCTIVITY
USE KINDS,              ONLY: dp, i4
USE CRUST_CONDUCTIVITY, ONLY: POTEKHINC
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: TOL = 5.0E-6_dp
INTEGER(KIND=i4) :: N_FAIL
REAL(KIND=dp) :: CKAPPA, CKAPPAT, CKAPPAH, TAU

N_FAIL = 0

! CASE1: fort.34's own first crust row (rhocgs~1.22d14, Fe-56),
! nstot.f's own hardcoded T=1d9 K, B12=10, Ximp=1d-1.
CALL POTEKHINC(1.0E9_dp, 1.2200E14_dp, 10.0_dp, 56.0_dp, 26.0_dp, 1.0_dp, 0.1_dp, &
  CKAPPA, CKAPPAT, CKAPPAH, TAU)
CALL CHECK_CLOSE("case1_ckappa",  CKAPPA,  5.630668328618629E+21_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case1_ckappat", CKAPPAT, 4.755638898667598E+21_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case1_ckappah", CKAPPAH, 2.039939344686625E+21_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case1_tau",     TAU,     7.335774440795409E+02_dp, TOL, N_FAIL)

! CASE2: a lower-density crust row (rhocgs~1d11), same composition.
CALL POTEKHINC(1.0E9_dp, 1.0E11_dp, 10.0_dp, 56.0_dp, 26.0_dp, 1.0_dp, 0.1_dp, &
  CKAPPA, CKAPPAT, CKAPPAH, TAU)
CALL CHECK_CLOSE("case2_ckappa",  CKAPPA,  5.756682995801160E+18_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case2_ckappat", CKAPPAT, 4.500912195155052E+18_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case2_ckappah", CKAPPAH, 2.377421028569920E+18_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case2_tau",     TAU,     8.457079938512022E+01_dp, TOL, N_FAIL)

! CASE3: higher density, different composition (Z=40,A=120), and no
! impurity scattering (Zimp=0, exercising CONDEGINC's ZIMP<=0 skip path).
CALL POTEKHINC(5.0E8_dp, 1.0E13_dp, 10.0_dp, 120.0_dp, 40.0_dp, 1.0_dp, 0.0_dp, &
  CKAPPA, CKAPPAT, CKAPPAH, TAU)
CALL CHECK_CLOSE("case3_ckappa",  CKAPPA,  1.175833683779469E+20_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case3_ckappat", CKAPPAT, 1.078634920016664E+20_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case3_ckappah", CKAPPAH, 3.237948705445316E+19_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case3_tau",     TAU,     1.996882377110178E+02_dp, TOL, N_FAIL)

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
