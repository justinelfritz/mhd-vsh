MODULE TRANSFORM_PLAN
!> Precomputed angular-quadrature tables for the physical<->spectral
!> transforms in TRANSFORMS: the scalar spherical harmonic Y_l^m and the
!> horizontal (theta,phi) components of the two horizontal polar-VSH
!> basis families (FORTVSH's PVSH_POL/PVSH_TOR), evaluated once at every
!> quadrature point of an ANGULAR_GRID_T and reused every transform call
!> instead of recomputed each time. The radial polar-VSH basis member
!> (PVSH_RAD) is not tabulated separately -- FORTVSH defines it as
!> exactly Y_l^m in its own r-hat component, so the YLM table below
!> doubles as that basis. This is the same 2D Gauss-Legendre(theta) x
!> uniform(phi) quadrature FORTVSH's own vsh_decomposition example builds
!> by hand (DECOMPOSE/RECONSTRUCT), hoisted out of the per-transform hot
!> loop and reused across every radial shell and timestep.
USE KINDS,        ONLY: dp, i4
USE GLOBALS,      ONLY: pi
USE VSH,          ONLY: SSH_ALL, PVSH_POL_ALL, PVSH_TOR_ALL
USE GRID_ANGULAR, ONLY: ANGULAR_GRID_T
IMPLICIT NONE
PRIVATE
PUBLIC :: TRANSFORM_PLAN_T, BUILD_TRANSFORM_PLAN

TYPE :: TRANSFORM_PLAN_T
  INTEGER(KIND=i4) :: LMAX = 0, NLM = 0, NTHETA = 0, NPHI = 0
  REAL(KIND=dp), ALLOCATABLE    :: QUADW(:)             ! (NTHETA), w_theta*dphi
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:,:,:)            ! (NLM,NTHETA,NPHI)
  COMPLEX(KIND=dp), ALLOCATABLE :: POL_TH(:,:,:), POL_PH(:,:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: TOR_TH(:,:,:), TOR_PH(:,:,:)
END TYPE TRANSFORM_PLAN_T

CONTAINS

!> @param PLAN Output plan.
!> @param AGRID Angular grid (quadrature nodes/weights) to build the
!>   tables on; TRANSFORMS routines called with the resulting PLAN must
!>   use physical arrays sized AGRID%NTHETA x AGRID%NPHI.
SUBROUTINE BUILD_TRANSFORM_PLAN(PLAN, AGRID)
  TYPE(TRANSFORM_PLAN_T), INTENT(OUT) :: PLAN
  TYPE(ANGULAR_GRID_T),   INTENT(IN)  :: AGRID
  INTEGER(KIND=i4) :: IQ, IP
  REAL(KIND=dp)    :: DPHI
  COMPLEX(KIND=dp), ALLOCATABLE :: POLALL(:,:), TORALL(:,:)

  PLAN%LMAX   = AGRID%LMAX
  PLAN%NLM    = AGRID%NLM
  PLAN%NTHETA = AGRID%NTHETA
  PLAN%NPHI   = AGRID%NPHI

  ALLOCATE(PLAN%QUADW(PLAN%NTHETA))
  ALLOCATE(PLAN%YLM(PLAN%NLM, PLAN%NTHETA, PLAN%NPHI))
  ALLOCATE(PLAN%POL_TH(PLAN%NLM, PLAN%NTHETA, PLAN%NPHI))
  ALLOCATE(PLAN%POL_PH(PLAN%NLM, PLAN%NTHETA, PLAN%NPHI))
  ALLOCATE(PLAN%TOR_TH(PLAN%NLM, PLAN%NTHETA, PLAN%NPHI))
  ALLOCATE(PLAN%TOR_PH(PLAN%NLM, PLAN%NTHETA, PLAN%NPHI))

  DPHI = 2.0_dp*pi/REAL(PLAN%NPHI, KIND=dp)
  DO IQ = 1, PLAN%NTHETA
    PLAN%QUADW(IQ) = AGRID%WTHETA(IQ) * DPHI
  END DO

  ALLOCATE(POLALL(3,PLAN%NLM), TORALL(3,PLAN%NLM))
  DO IP = 1, PLAN%NPHI
    DO IQ = 1, PLAN%NTHETA
      CALL SSH_ALL(PLAN%YLM(:,IQ,IP), PLAN%LMAX, AGRID%THETA(IQ), AGRID%PHI(IP))
      CALL PVSH_POL_ALL(POLALL, PLAN%LMAX, AGRID%THETA(IQ), AGRID%PHI(IP))
      CALL PVSH_TOR_ALL(TORALL, PLAN%LMAX, AGRID%THETA(IQ), AGRID%PHI(IP))
      PLAN%POL_TH(:,IQ,IP) = POLALL(2,:)
      PLAN%POL_PH(:,IQ,IP) = POLALL(3,:)
      PLAN%TOR_TH(:,IQ,IP) = TORALL(2,:)
      PLAN%TOR_PH(:,IQ,IP) = TORALL(3,:)
    END DO
  END DO
  DEALLOCATE(POLALL, TORALL)
END SUBROUTINE BUILD_TRANSFORM_PLAN

END MODULE TRANSFORM_PLAN
