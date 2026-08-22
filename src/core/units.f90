MODULE UNITS
!> Names this project's adopted physical unit system and provides the
!> multiplicative factors needed to report a diagnostic in physical
!> Gaussian-cgs units (erg, erg/s) instead of code units. Nothing else
!> in MHD-VSH depends on this module, and nothing else should: every
!> PDE/diagnostic formula elsewhere (DIFFUSION_REGIME, FIELD_DIAGNOSTICS,
!> the Hall induction equations in analytic_formulas/mhd-vsh-relations.tex)
!> is written in genuinely unit-system-agnostic form -- ETA/F_HALL/DT/
!> R_MIN/R_MAX are free REAL(dp) inputs with no unit baked into the code,
!> so those formulas are already dimensionally correct in ANY consistent
!> choice of units, provided ETA is given in [length]**2/[time] and
!> F_HALL in [length]**2/([B]*[time]) to match whatever [length],[time],
!> [B] the caller has chosen for R_MIN/R_MAX/DT/the state's Phi,Psi. This
!> module exists only to record what THIS project has chosen (per Justin
!> Elfritz, 2026-08-20) -- B: 10**12 Gauss, length: km, time: (Julian)
!> year -- and to convert diagnostic OUTPUTS to erg/erg-per-second at the
!> point they're printed or logged, never inside the core solve/
!> diagnostic formulas themselves.
!>
!> Energy converts by ENERGY_UNIT_ERG = B_UNIT_GAUSS**2 * LENGTH_UNIT_CM**3
!> on dimensional grounds alone (Gaussian cgs: energy = integral of
!> B**2/8pi over a volume, so [energy] = [B]**2*[length]**3 regardless of
!> the particular formula computing it) -- this holds for
!> TOTAL_(POLOIDAL/TOROIDAL/)_MAGNETIC_ENERGY without needing to re-derive
!> FIELD_DIAGNOSTICS' internal radial-integral bookkeeping. Rates
!> (JOULE_DISSIPATION_RATE, POYNTING_FLUX_RATE, HALL_POYNTING_FLUX_RATE --
!> code-unit energy per code-unit time, i.e. per year) convert the same
!> way divided by TIME_UNIT_S, to erg/s.
!>
!> @warning A code-unit numeric result times these factors is only
!>   physically meaningful if R_MIN/R_MAX/DT/ETA/F_HALL were themselves
!>   entered in km/yr/km**2 per yr/km**2 per (10**12 G) per yr consistently
!>   for that run -- this module cannot check that; it is pure
!>   documentation plus arithmetic.
USE KINDS, ONLY: dp
IMPLICIT NONE
PRIVATE
PUBLIC :: B_UNIT_GAUSS, LENGTH_UNIT_CM, TIME_UNIT_S, &
          ENERGY_UNIT_ERG, POWER_UNIT_ERG_PER_S

!> This project's field unit: 10**12 Gauss.
REAL(KIND=dp), PARAMETER :: B_UNIT_GAUSS = 1.0E12_dp
!> This project's length unit: 1 km, in cm.
REAL(KIND=dp), PARAMETER :: LENGTH_UNIT_CM = 1.0E5_dp
!> This project's time unit: 1 Julian year (365.25 d), in seconds.
REAL(KIND=dp), PARAMETER :: TIME_UNIT_S = 3.15576E7_dp

!> Multiply a code-unit energy (from TOTAL_MAGNETIC_ENERGY etc.) by this
!> to report erg.
REAL(KIND=dp), PARAMETER :: ENERGY_UNIT_ERG = &
  B_UNIT_GAUSS**2 * LENGTH_UNIT_CM**3
!> Multiply a code-unit rate (energy per code-unit time, i.e. per year --
!> from JOULE_DISSIPATION_RATE/POYNTING_FLUX_RATE/HALL_POYNTING_FLUX_RATE)
!> by this to report erg/s.
REAL(KIND=dp), PARAMETER :: POWER_UNIT_ERG_PER_S = ENERGY_UNIT_ERG / TIME_UNIT_S

END MODULE UNITS
