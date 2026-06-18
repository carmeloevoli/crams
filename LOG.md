# Changelog — `cramsmcmc` vs `CrankNicholson`

This document summarizes all differences between the previous release branch
**`CrankNicholson`** (version `1.1`) and the current branch **`cramsmcmc`**
(version `2.0`).

`cramsmcmc` is a strict descendant of `CrankNicholson`: the history is linear
(35 commits ahead, no divergence), so everything below is *added*, *removed*, or
*changed* relative to `CrankNicholson`. In aggregate the diff is ~29k insertions
and ~59k deletions across 214 files — the bulk of which is a source-tree
reorganization, a data-format overhaul, and the replacement of the C++ fitting
code with a Python MCMC suite.

---

## 1. Version

- `VERSION`: **`1.1` → `2.0`**.

## 2. Source tree reorganization

The flat `include/` and `src/` layout was replaced with a namespaced layout
under `crams/`:

- Headers moved into `include/crams/{core,physics,utils}/…`:
  - `core/`: `input.h`, `output.h`, `pid.h`, `cgs.h`
  - `physics/`: `grammage.h`, `losses.h`, `primary.h`
  - `utils/`: `csvreader.h`, `grid.h`, `logging.h`, `numeric.h`, `utilities.h`, `git_revision.h`
  - top-level: `crams.h` (umbrella header), `particle.h`, `particlelist.h`,
    `fragmentation.h`, `inelastic.h`, `secondary.h`
- Sources mirror the same structure under `src/{core,physics,utils}/…`.
- New executable entry point **`apps/main.cpp`** replaces `src/main.cpp`.
- Removed the old `include/xsecs/` headers and `src/xsecs/` implementations
  (`Evoli2019`, `Korsmeier2018`, `Tripathi99`) — cross sections are now data-driven
  (see §4).

## 3. Build system, tooling & CI

- **`CMakeLists.txt`** modernized:
  - `cmake_minimum_required` 3.2 → **3.14**.
  - Removed hard-coded compiler path and inline `-std`/`-fopenmp`/`-D…` flags;
    now uses `CMAKE_CXX_STANDARD 14`, `-Wall -Wextra -Wpedantic -march=native`,
    and exports `compile_commands.json`.
  - `find_package(GSL)` → `find_package(GSL REQUIRED)`.
  - Library renamed `CRAMS_LIB` → `crams_core`.
  - Added `enable_testing()` and a full set of `add_test(...)` targets (see §6).
- **plog** is no longer vendored in the repo. The bundled `external/plog/**` was
  removed and is fetched on demand via the new **`external/get_plog.sh`**.
- Removed the custom **`cmake/FindGSL.cmake`** (uses CMake's built-in GSL module).
- **`.github/workflows/ci.yml`** added: builds and runs the CTest suite on
  `ubuntu-latest` and `macos-latest`, installing GSL and fetching plog first.
- **`.clang-format`** modernized (deprecated keys replaced) and applied across the
  whole source tree.
- `.gitignore` updated.

## 4. Data tables & cross-section models

The old `.txt` data tables were removed wholesale and replaced with a unified,
validated CSV format, with **runtime-selectable** cross-section models:

- Removed: all `data/AMS-02_*`, `data/DAMPE_*`, `data/*_AMS-02_R.txt`,
  `data/sigProd*.txt`, `data/xsecs_*.txt`, `data/crxsecs_fragmentation_Evoli2019_*`,
  `data/nucleilist.csv`, `data/solarsystem_abundances2003.txt`, and the large
  `supplementary__XS_table_Param_II_B.txt`.
- Added CSV tables:
  - Nuclei list: `data/crams_nucleilist.csv`.
  - **Inelastic** (destruction) models: `crams_inelastic_crosec.csv`,
    `crams_inelastic_glauber.csv`, `crams_inelastic_tripathi99.csv`.
  - **Fragmentation** (production) models: `crams_fragmentation_evoli2019.csv`,
    `…_evoli2026_st99.csv`, `…_evoli2026_w93.csv`, `…_fluka4dragon.csv`,
    `…_usine_galprop17_opt12.csv`, `…_usine_galprop17_opt22.csv`,
    `…_usine_webber03+coste12.csv`.
- Cross-section model selection is now a runtime input parameter rather than a
  compile-time `-D` define; grid validation is applied on load
  (`tests/check_data_tables.cmake`).
- The numerical interpolation layer was reworked: the GSL wrapper `include/gsl.h`
  and old `grid.h` were replaced by `include/crams/utils/numeric.h` and a new
  `grid.h` (GSL is still used underneath for integration/2D splines).

## 5. Fitting: C++ chi² removed, Python MCMC added

- **Removed** the in-tree C++ optimizer and its stray copies:
  `src/chi2.cpp`, `src/chi2 copy.cpp`, `src/chi2 copy 2.cpp`, `include/chi2.h`,
  and the top-level `optimize_chi2.py`.
- **Added** a complete Python MCMC / best-fit workflow under `mcmc/`:
  - `run_mcmc.py`, `runner.py`, `fitting.py`, `find_bestfit.py`,
    `write_bestfit_ini.py`, `run_all_bestfits.sh`.
  - Plotting: `plot_bestfit.py`, `plot_model.py`, `plot_model_Be.py`,
    `plot_corner.py`.
  - `requirements.txt` for the Python dependencies.
  - Preliminary AMS-02 Be data: `mcmc/preliminary/AMS-02_preliminary_Be_isotopes.csv`,
    `…_Be_ratios.csv`.
  - Features added incrementally: free `d0`/`delta`, B/C in the fit,
    log/ratio reparametrization, selectable cross-section model and halo size,
    best-fit-start MCMC config, walker/step/autocorrelation tuning.

## 6. Test suite

A unit-test suite was introduced (previously only ad-hoc plotting sandboxes
existed). Removed `tests/apXsecsPlot.cpp`, `tests/inelasticPlot.cpp`,
`tests/sandbox.cpp`; added CTest-registered tests:

- `test_cgs`, `test_utilities`, `test_pid`, `test_numeric`, `test_csvreader`,
  `test_grid`, `test_grammage`, `test_losses`, `test_inelastic`, `test_models`,
  `test_input`, `test_output`, `test_particle`, `test_particlelist`,
  `test_primary`, `test_secondary`, `test_example_ini`.
- `tests/check_data_tables.cmake` validates the data CSV grids at configure/test time.
- `tests/test_example_ini.cpp` parses `examples/crams.ini` so the documented
  example never drifts out of sync with the parser.

## 7. Physics & runtime behavior

- Added a tabulated **fragmentation** model and wired it into the secondary source.
- Added a **Glauber** inelastic model alongside CROSEC and Tripathi99;
  Tripathi99 now throws on out-of-range kinetic energy.
- Enabled a **tertiary proton** source.
- Ongoing work on **Be isotopes** (e.g. `10Be`), including `tauhalflife` fixes and
  nuclei-list updates.
- Performance: skips per-particle dumps in quiet mode and optimizes the secondary
  source computation. A `-q` / `--quiet` flag was added to suppress per-particle logging.
- Simulation name (`simname`) is now derived from the input `.ini` basename.

## 8. Documentation & examples

- **`README.md`** rewritten for v2.0: CI/license/arXiv badges, physics overview,
  build/run instructions, and a Data tables section.
- **`examples/crams.ini`** added: a fully commented reference input file that
  documents every parameter and is exercised by the test suite.

---

*Generated 2026-06-18. Base: `CrankNicholson` (v1.1); head: `cramsmcmc` (v2.0).*
