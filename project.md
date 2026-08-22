---
project: MHD-VSH
summary: A spherical MHD solver built on FORTVSH's vector spherical harmonics, evolving magnetic fields as poloidal/toroidal potentials (Phi,Psi) rather than raw field components.
author: Justin G. Elfritz
email: 45849387+justinelfritz@users.noreply.github.com
version: 0.1.0
year: 2026
src_dir: src
          app
output_dir: docs
display: public
         protected
         private
source: true
graph: true
coloured_edges: true
search: true
sort: permission-alpha
---

MHD-VSH advances magnetic fields as a poloidal/toroidal potential pair
(`Phi`,`Psi`) -- not the raw field components -- built on
[FORTVSH](https://github.com/justinelfritz/FORTVSH)'s vector spherical
harmonics. Every regime (pure Ohmic diffusion today; electron-MHD Hall
and, eventually, full ion-MHD dynamo physics) assembles its own PDE from
the same regime-agnostic building blocks: radial finite-difference
operators, boundary conditions, and a dense LAPACK solver.

## Where to start

- [[GRID_RADIAL]], [[GRID_ANGULAR]], [[FIELD_TYPES]] -- the grids and
  field containers every regime is built on.
- [[RADIAL_OPERATORS]] -- finite-difference `D1`/`D2`, the curvature
  term, and origin regularity.
- [[BOUNDARY_CONDITIONS]] -- the vacuum outer boundary (Robin on `Phi`,
  Dirichlet on `Psi`) and the inner core-crust boundary.
- [[LINEAR_SOLVE]] -- a physics-blind dense LU wrapper (LAPACK
  `DGETRF`/`DGETRS`); every regime assembles its own system matrix and
  hands it here only at the very last step.
- [[REGIME_INTERFACE]] and [[TIMESTEPPER]] -- the opaque per-regime
  `ADVANCE` contract and the shared fixed-step loop that drives it.
- [[FIELD_DIAGNOSTICS]] -- time-resolved diagnostics computed directly
  from `(Phi,Psi)`, including the Hall-CFL-limited timestep
  (`HALL_COURANT_TIMESTEP`); see its header for which coefficients are
  provisional, pending re-derivation.
- [[UNITS]] -- this project's adopted code-unit system (B, length,
  time) and the conversion factors to physical cgs, applied only at
  the reporting layer.
- [[TIMESTEPPER]] -- both the fixed-`DT` `RUN` loop and the
  Hall-CFL-driven `RUN_ADAPTIVE` variant.
- [[IO_CHECKPOINT]] -- checkpoint/restart for a long-running regime
  driver (`WRITE_CHECKPOINT`/`READ_CHECKPOINT`).
- `src/regimes/diffusion/diffusion_regime.f90` -- the first working
  regime (pure linear diffusion), and the reference pattern for every
  regime after it.
- `src/regimes/hall/hall_regime.f90` -- the combined resistive+Hall
  regime (implicit diffusion + explicit Hall substeps), built on
  [[HALL_INDUCTION]]'s RHS evaluator.
- [[TOV_SOLVER]], [[EOS_TABLE]], [[CRUST_CONDUCTIVITY]] -- a ported
  TOV + crust-EOS/conductivity solver supplying a real radial
  `n_e(r)`/`eta(r)`/`f_H(r)` profile (`app/mhdvsh_tov.f90`), consumed
  by the Hall regime's optional `ETA_PROFILE`/`F_HALL_PROFILE`
  arguments.

## Building

```bash
cmake -B build
cmake --build build
cd build && ctest --output-on-failure
```

See `README.md` in the repository root for the full build guide.
