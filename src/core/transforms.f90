MODULE TRANSFORMS
!> Physical<->spectral angular transforms: forward (physical to spectral,
!> i.e. analysis) via the 2D Gauss-Legendre(theta) x uniform(phi)
!> quadrature tabulated in a TRANSFORM_PLAN_T, inverse (spectral to
!> physical, i.e. synthesis) the matching sum over modes. Mirrors FORTVSH's own
!> vsh_decomposition example (DECOMPOSE/RECONSTRUCT) but against this
!> codebase's (N_r,Ntheta,Nphi)/(N_r,Nlm) field layout, with every radial
!> shell an independent 2D transform -- the loop OpenMP parallelizes over
!> when MHDVSH_ENABLE_OPENMP is on.
USE KINDS,          ONLY: dp, i4
USE TRANSFORM_PLAN, ONLY: TRANSFORM_PLAN_T
USE FIELD_TYPES,    ONLY: SPECTRAL_SCALAR_T, SPECTRAL_VECTOR3_T, &
                           PHYSICAL_SCALAR_T, PHYSICAL_VECTOR_T
IMPLICIT NONE
PRIVATE
PUBLIC :: FORWARD_SCALAR, INVERSE_SCALAR, FORWARD_VECTOR, INVERSE_VECTOR

CONTAINS

!> Forward transform (analysis): SPEC%COEF(r,lm) = sum_quad PHYS%F(r,theta,phi) *
!> conj(Y_lm(theta,phi)) * w_theta*w_phi, independently at every radius.
!>
!> @param SPEC Output spectral scalar (already allocated to match PHYS/PLAN).
!> @param PHYS Input physical scalar on PLAN's angular grid.
!> @param PLAN Precomputed quadrature/basis tables (BUILD_TRANSFORM_PLAN).
SUBROUTINE FORWARD_SCALAR(SPEC, PHYS, PLAN)
  TYPE(SPECTRAL_SCALAR_T), INTENT(INOUT) :: SPEC
  TYPE(PHYSICAL_SCALAR_T), INTENT(IN)    :: PHYS
  TYPE(TRANSFORM_PLAN_T),  INTENT(IN)    :: PLAN
  INTEGER(KIND=i4) :: IR, IQ, IP
  COMPLEX(KIND=dp) :: ACC(PLAN%NLM)

  !$OMP PARALLEL DO PRIVATE(IR,IQ,IP,ACC)
  DO IR = 1, PHYS%N_R
    ACC = (0.0_dp, 0.0_dp)
    DO IP = 1, PLAN%NPHI
      DO IQ = 1, PLAN%NTHETA
        ACC = ACC + (PHYS%F(IR,IQ,IP)*PLAN%QUADW(IQ)) * CONJG(PLAN%YLM(:,IQ,IP))
      END DO
    END DO
    SPEC%COEF(IR,:) = ACC
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE FORWARD_SCALAR

!> Inverse transform (synthesis): PHYS%F(r,theta,phi) = Re[ sum_lm SPEC%COEF(r,lm) *
!> Y_lm(theta,phi) ], the inverse of FORWARD_SCALAR.
!>
!> @param PHYS Output physical scalar (already allocated to match SPEC/PLAN).
!> @param SPEC Input spectral scalar.
!> @param PLAN Precomputed quadrature/basis tables (BUILD_TRANSFORM_PLAN).
SUBROUTINE INVERSE_SCALAR(PHYS, SPEC, PLAN)
  TYPE(PHYSICAL_SCALAR_T), INTENT(INOUT) :: PHYS
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)    :: SPEC
  TYPE(TRANSFORM_PLAN_T),  INTENT(IN)    :: PLAN
  INTEGER(KIND=i4) :: IR, IQ, IP

  !$OMP PARALLEL DO PRIVATE(IR,IQ,IP)
  DO IR = 1, SPEC%N_R
    DO IP = 1, PLAN%NPHI
      DO IQ = 1, PLAN%NTHETA
        PHYS%F(IR,IQ,IP) = REAL(SUM(SPEC%COEF(IR,:)*PLAN%YLM(:,IQ,IP)), KIND=dp)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE INVERSE_SCALAR

!> Forward transform (analysis): projects PHYS (physical r/theta/phi vector components) onto
!> the three polar-VSH bases at every radius, giving SPEC%RAD/POL/TOR --
!> the vector analogue of FORWARD_SCALAR (PVSH_RAD's basis vector is
!> (Y_lm,0,0), so the RAD projection reuses PLAN%YLM exactly as
!> FORWARD_SCALAR does).
!>
!> @param SPEC Output spectral vector (already allocated to match PHYS/PLAN).
!> @param PHYS Input physical vector on PLAN's angular grid.
!> @param PLAN Precomputed quadrature/basis tables (BUILD_TRANSFORM_PLAN).
SUBROUTINE FORWARD_VECTOR(SPEC, PHYS, PLAN)
  TYPE(SPECTRAL_VECTOR3_T), INTENT(INOUT) :: SPEC
  TYPE(PHYSICAL_VECTOR_T),  INTENT(IN)    :: PHYS
  TYPE(TRANSFORM_PLAN_T),   INTENT(IN)    :: PLAN
  INTEGER(KIND=i4) :: IR, IQ, IP
  REAL(KIND=dp)    :: W
  COMPLEX(KIND=dp) :: ACC_RAD(PLAN%NLM), ACC_POL(PLAN%NLM), ACC_TOR(PLAN%NLM)

  !$OMP PARALLEL DO PRIVATE(IR,IQ,IP,W,ACC_RAD,ACC_POL,ACC_TOR)
  DO IR = 1, PHYS%N_R
    ACC_RAD = (0.0_dp, 0.0_dp)
    ACC_POL = (0.0_dp, 0.0_dp)
    ACC_TOR = (0.0_dp, 0.0_dp)
    DO IP = 1, PLAN%NPHI
      DO IQ = 1, PLAN%NTHETA
        W = PLAN%QUADW(IQ)
        ACC_RAD = ACC_RAD + (W*PHYS%R(IR,IQ,IP)) * CONJG(PLAN%YLM(:,IQ,IP))
        ACC_POL = ACC_POL + W * ( PHYS%TH(IR,IQ,IP)*CONJG(PLAN%POL_TH(:,IQ,IP)) &
                                 + PHYS%PH(IR,IQ,IP)*CONJG(PLAN%POL_PH(:,IQ,IP)) )
        ACC_TOR = ACC_TOR + W * ( PHYS%TH(IR,IQ,IP)*CONJG(PLAN%TOR_TH(:,IQ,IP)) &
                                 + PHYS%PH(IR,IQ,IP)*CONJG(PLAN%TOR_PH(:,IQ,IP)) )
      END DO
    END DO
    SPEC%RAD(IR,:) = ACC_RAD
    SPEC%POL(IR,:) = ACC_POL
    SPEC%TOR(IR,:) = ACC_TOR
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE FORWARD_VECTOR

!> Inverse transform (synthesis): reconstructs PHYS (physical r/theta/phi vector components)
!> from SPEC%RAD/POL/TOR, the inverse of FORWARD_VECTOR.
!>
!> @param PHYS Output physical vector (already allocated to match SPEC/PLAN).
!> @param SPEC Input spectral vector.
!> @param PLAN Precomputed quadrature/basis tables (BUILD_TRANSFORM_PLAN).
SUBROUTINE INVERSE_VECTOR(PHYS, SPEC, PLAN)
  TYPE(PHYSICAL_VECTOR_T),  INTENT(INOUT) :: PHYS
  TYPE(SPECTRAL_VECTOR3_T), INTENT(IN)    :: SPEC
  TYPE(TRANSFORM_PLAN_T),   INTENT(IN)    :: PLAN
  INTEGER(KIND=i4) :: IR, IQ, IP

  !$OMP PARALLEL DO PRIVATE(IR,IQ,IP)
  DO IR = 1, SPEC%N_R
    DO IP = 1, PLAN%NPHI
      DO IQ = 1, PLAN%NTHETA
        PHYS%R(IR,IQ,IP)  = REAL(SUM(SPEC%RAD(IR,:)*PLAN%YLM(:,IQ,IP)), KIND=dp)
        PHYS%TH(IR,IQ,IP) = REAL(SUM(SPEC%POL(IR,:)*PLAN%POL_TH(:,IQ,IP) &
                                    + SPEC%TOR(IR,:)*PLAN%TOR_TH(:,IQ,IP)), KIND=dp)
        PHYS%PH(IR,IQ,IP) = REAL(SUM(SPEC%POL(IR,:)*PLAN%POL_PH(:,IQ,IP) &
                                    + SPEC%TOR(IR,:)*PLAN%TOR_PH(:,IQ,IP)), KIND=dp)
      END DO
    END DO
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE INVERSE_VECTOR

END MODULE TRANSFORMS
