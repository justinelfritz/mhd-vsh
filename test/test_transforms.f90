!> Phase 1 check: forward(inverse(forward(f))) == forward(f) for both the
!> scalar and vector angular transforms, at every radial shell
!> independently.
!>
!> Why this property and not a direct spectral round trip: INVERSE_*
!> keeps only Re(sum_lm c_lm*Ylm), so it only reproduces its input
!> coefficients exactly when they already satisfy the real-field
!> conjugate symmetry c(l,-m) = (-1)**m*conj(c(l,m)) (scalar case; the
!> vector bases have their own, messier version of the same symmetry
!> inherited from FORTVSH's VSH_CORE). Picking coefficients that satisfy
!> it by construction risks getting that symmetry subtly wrong by hand.
!> Starting instead from an arbitrary REAL physical field sidesteps the
!> issue entirely: standard spherical-harmonic-transform theory
!> guarantees FORWARD of any real field already lands on that symmetric
!> subspace, and the non-dealiased quadrature (exact for products of two
!> Lmax-truncated fields) makes the projection exact -- so a second
!> forward/inverse pair must reproduce the first pair's coefficients
!> exactly, with no assumption on the original field's own band-limit.
PROGRAM TEST_TRANSFORMS
USE KINDS,           ONLY: dp, i4
USE GRID_ANGULAR,    ONLY: ANGULAR_GRID_T, BUILD_ANGULAR_GRID
USE TRANSFORM_PLAN,  ONLY: TRANSFORM_PLAN_T, BUILD_TRANSFORM_PLAN
USE TRANSFORMS,      ONLY: FORWARD_SCALAR, INVERSE_SCALAR, &
                            FORWARD_VECTOR, INVERSE_VECTOR
USE FIELD_TYPES,     ONLY: SPECTRAL_SCALAR_T, SPECTRAL_VECTOR3_T, &
                            PHYSICAL_SCALAR_T, PHYSICAL_VECTOR_T, &
                            ALLOC_SPECTRAL_SCALAR, ALLOC_SPECTRAL_VECTOR3, &
                            ALLOC_PHYSICAL_SCALAR, ALLOC_PHYSICAL_VECTOR
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: LMAX = 8
INTEGER(KIND=i4), PARAMETER :: N_R  = 3
REAL(KIND=dp),    PARAMETER :: TOL  = 1.0E-10_dp

TYPE(ANGULAR_GRID_T)    :: AGRID
TYPE(TRANSFORM_PLAN_T)  :: PLAN
INTEGER(KIND=i4) :: N_FAIL

CALL BUILD_ANGULAR_GRID(AGRID, LMAX, DEALIAS=.FALSE.)
CALL BUILD_TRANSFORM_PLAN(PLAN, AGRID)

N_FAIL = 0
CALL RUN_SCALAR_ROUNDTRIP(N_FAIL)
CALL RUN_VECTOR_ROUNDTRIP(N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_transforms"
END IF

CONTAINS

!> Arbitrary deterministic real field, distinct at every (r,theta,phi)
!> triple (including a radius-dependent phase, so a radial-index bug
!> such as transposed loop bounds would show up as a mismatch).
!>
!> @param PHASE Per-field phase offset, so the scalar field and the three
!>   vector components don't happen to coincide.
FUNCTION ARBITRARY_FIELD(IR, IQ, IP, PHASE) RESULT(F)
  INTEGER(KIND=i4), INTENT(IN) :: IR, IQ, IP
  REAL(KIND=dp),    INTENT(IN) :: PHASE
  REAL(KIND=dp) :: F
  F = SIN(0.7_dp*IR + 1.3_dp*IQ - 0.5_dp*IP + PHASE) &
    + 0.3_dp*COS(0.2_dp*IR*IQ - 0.4_dp*IP + PHASE)
END FUNCTION ARBITRARY_FIELD

SUBROUTINE RUN_SCALAR_ROUNDTRIP(N_FAIL)
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  TYPE(SPECTRAL_SCALAR_T) :: SPEC1, SPEC2
  TYPE(PHYSICAL_SCALAR_T) :: PHYS0, PHYS1
  INTEGER(KIND=i4) :: IR, IQ, IP
  REAL(KIND=dp) :: MAX_ERR

  CALL ALLOC_SPECTRAL_SCALAR(SPEC1, N_R, LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(SPEC2, N_R, LMAX)
  CALL ALLOC_PHYSICAL_SCALAR(PHYS0, N_R, AGRID%NTHETA, AGRID%NPHI)
  CALL ALLOC_PHYSICAL_SCALAR(PHYS1, N_R, AGRID%NTHETA, AGRID%NPHI)

  DO IP = 1, AGRID%NPHI
    DO IQ = 1, AGRID%NTHETA
      DO IR = 1, N_R
        PHYS0%F(IR,IQ,IP) = ARBITRARY_FIELD(IR, IQ, IP, 0.0_dp)
      END DO
    END DO
  END DO

  CALL FORWARD_SCALAR(SPEC1, PHYS0, PLAN)
  CALL INVERSE_SCALAR(PHYS1, SPEC1, PLAN)
  CALL FORWARD_SCALAR(SPEC2, PHYS1, PLAN)

  MAX_ERR = MAXVAL(ABS(SPEC2%COEF - SPEC1%COEF))
  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  scalar_roundtrip  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  scalar_roundtrip  max_err=", MAX_ERR
  END IF
END SUBROUTINE RUN_SCALAR_ROUNDTRIP

SUBROUTINE RUN_VECTOR_ROUNDTRIP(N_FAIL)
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  TYPE(SPECTRAL_VECTOR3_T) :: SPEC1, SPEC2
  TYPE(PHYSICAL_VECTOR_T)  :: PHYS0, PHYS1
  INTEGER(KIND=i4) :: IR, IQ, IP
  REAL(KIND=dp) :: MAX_ERR

  CALL ALLOC_SPECTRAL_VECTOR3(SPEC1, N_R, LMAX)
  CALL ALLOC_SPECTRAL_VECTOR3(SPEC2, N_R, LMAX)
  CALL ALLOC_PHYSICAL_VECTOR(PHYS0, N_R, AGRID%NTHETA, AGRID%NPHI)
  CALL ALLOC_PHYSICAL_VECTOR(PHYS1, N_R, AGRID%NTHETA, AGRID%NPHI)

  DO IP = 1, AGRID%NPHI
    DO IQ = 1, AGRID%NTHETA
      DO IR = 1, N_R
        PHYS0%R(IR,IQ,IP)  = ARBITRARY_FIELD(IR, IQ, IP, 0.0_dp)
        PHYS0%TH(IR,IQ,IP) = ARBITRARY_FIELD(IR, IQ, IP, 1.0_dp)
        PHYS0%PH(IR,IQ,IP) = ARBITRARY_FIELD(IR, IQ, IP, 2.0_dp)
      END DO
    END DO
  END DO

  CALL FORWARD_VECTOR(SPEC1, PHYS0, PLAN)
  CALL INVERSE_VECTOR(PHYS1, SPEC1, PLAN)
  CALL FORWARD_VECTOR(SPEC2, PHYS1, PLAN)

  MAX_ERR = MAX(MAXVAL(ABS(SPEC2%RAD - SPEC1%RAD)), &
                MAXVAL(ABS(SPEC2%POL - SPEC1%POL)), &
                MAXVAL(ABS(SPEC2%TOR - SPEC1%TOR)))
  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  vector_roundtrip  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  vector_roundtrip  max_err=", MAX_ERR
  END IF
END SUBROUTINE RUN_VECTOR_ROUNDTRIP

END PROGRAM TEST_TRANSFORMS
