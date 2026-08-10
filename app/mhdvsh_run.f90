!> Phase 0 skeleton driver: builds an angular grid (spectral, via FORTVSH)
!> and a radial grid (finite-difference), and prints a summary confirming
!> the two link together. No physics yet -- that starts in Phase 3 once the
!> transform layer (Phase 1) and radial operators (Phase 2) are wired in.
PROGRAM MHDVSH_RUN
USE KINDS,        ONLY: dp, i4
USE GRID_ANGULAR,  ONLY: ANGULAR_GRID_T, BUILD_ANGULAR_GRID
USE GRID_RADIAL,   ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
IMPLICIT NONE

TYPE(ANGULAR_GRID_T) :: AGRID
TYPE(RADIAL_GRID_T)  :: RGRID
INTEGER(KIND=i4), PARAMETER :: LMAX = 8
INTEGER(KIND=i4), PARAMETER :: N_R  = 20
REAL(KIND=dp),    PARAMETER :: R_MAX = 1.0_dp

CALL BUILD_ANGULAR_GRID(AGRID, LMAX, DEALIAS=.FALSE.)
CALL BUILD_RADIAL_GRID(RGRID, N_R, 0.0_dp, R_MAX, FULL_SPHERE=.TRUE.)

WRITE(*,'(A)')            'MHD-VSH Phase 0 skeleton'
WRITE(*,'(A,I0)')         '  Lmax        = ', AGRID%LMAX
WRITE(*,'(A,I0)')         '  Nlm         = ', AGRID%NLM
WRITE(*,'(A,I0)')         '  Ntheta      = ', AGRID%NTHETA
WRITE(*,'(A,I0)')         '  Nphi        = ', AGRID%NPHI
WRITE(*,'(A,I0)')         '  N_r         = ', RGRID%N
WRITE(*,'(A,L1)')         '  full_sphere = ', RGRID%FULL_SPHERE
WRITE(*,'(A,ES12.5,A,ES12.5)') '  r(1),r(N)   = ', RGRID%R(1), ', ', RGRID%R(RGRID%N)

END PROGRAM MHDVSH_RUN
