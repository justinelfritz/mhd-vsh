# MHD-VSH

A spherical MHD solver built on [FORTVSH](https://github.com/justinelfritz/FORTVSH)'s
vector spherical harmonics. Magnetic fields are advanced as a
poloidal/toroidal potential pair (`Phi`, `Psi`) rather than raw field
components, matched to a vacuum exterior at the outer boundary.

See `ROADMAP.md` for current status and what's next, and `project.md`
for the FORD-generated API documentation.

## Building

Requires a Fortran compiler (gfortran), CMake >= 3.14, LAPACK, and a
built [FORTVSH](https://github.com/justinelfritz/FORTVSH) install.

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/fortvsh-install
cmake --build build
cd build && ctest --output-on-failure
```

`cmake --build build` also regenerates the diagnostic plots under
`build/artifacts/<regime>/<diagnostic>/plots/` as a normal build target
-- see `CMakeLists.txt`'s "Artifacts" section.

### Regime drivers

Each regime is its own executable (no runtime dispatch); a regime can
have more than one driver (e.g. Hall's fixed- vs. dynamic-`dt` variants):

```bash
./build/mhdvsh_diffusion [output_data_file]
./build/mhdvsh_hall [energy_data] [energy_by_l_data] [checkpoint_path] [resume_path]
./build/mhdvsh_hall_adaptive [output_data_file]
```

`mhdvsh_hall`'s checkpoint/resume args are optional and purely additive
-- see its own header comment for the exact semantics (`N_STEPS` is
always the total target step count, not "how many more to run").

### TOV / crust-EOS solver

`mhdvsh_tov` is not a regime driver -- it solves a single neutron star's
structure (ported from an existing F77 TOV+crust-EOS solver) and prints/
writes its radial `n_e(r)`/`eta(r)`/`f_H(r)` profile, which the Hall
regime's optional `ETA_PROFILE`/`F_HALL_PROFILE` arguments can consume
instead of a uniform scalar:

```bash
./build/mhdvsh_tov [eos_table_path] [output_data_file]
```

EOS tables live in `data/eos/`. `scripts/plot_eos_profile.py` plots a
run's output.

## Generating documentation

```bash
pip install ford
ford project.md
```

Writes to `docs/`; open `docs/index.html`.

## Project layout

- `src/core/` -- regime-agnostic building blocks: grids, radial
  operators, transforms, boundary conditions, the linear solver, the
  regime/timestepper contract, diagnostics, units, checkpoint/restart,
  plus a ported TOV + crust-EOS/conductivity solver (`tov_solver.f90`,
  `eos_table.f90`, `crust_conductivity.f90`, `ode_integrator.f90`).
- `src/regimes/<name>/` -- one module per physics regime, each owning
  its own PDE assembly.
- `app/` -- driver programs: one or more per regime, plus `mhdvsh_tov`
  (the TOV/EOS solver -- not itself a regime).
- `data/eos/` -- EOS table data files the TOV solver reads.
- `test/` -- CTest suite, one executable per module/regime.
- `scripts/` -- plotting scripts consuming `build/artifacts/` data.
