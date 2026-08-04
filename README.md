# CRAMS

[![CI](https://img.shields.io/github/actions/workflow/status/carmeloevoli/crams/ci.yml?label=build%20%26%20tests)](https://github.com/carmeloevoli/crams/actions/workflows/ci.yml)
[![C++14](https://img.shields.io/badge/C%2B%2B-14-00599C?logo=cplusplus&logoColor=white)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.14%2B-064F8C?logo=cmake&logoColor=white)](https://cmake.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Last commit](https://img.shields.io/github/last-commit/carmeloevoli/crams)](https://github.com/carmeloevoli/crams/commits)
[![arXiv](https://img.shields.io/badge/arXiv-1904.10220-B31B1B.svg)](https://arxiv.org/abs/1904.10220)
[![arXiv](https://img.shields.io/badge/arXiv-1910.04113-B31B1B.svg)](https://arxiv.org/abs/1910.04113)

**CRAMS** — *Cosmic Rays with All Measurements from Space* — is a C++ code that computes
the spectra of galactic cosmic-ray nuclei (H through Ni) above ~1 GV rigidity in the
**weighted-slab approximation** ([Jones et al., 2000](https://arxiv.org/abs/astro-ph/0007293))
of diffusive propagation. It is aimed at reproducing the absolute fluxes
and secondary-to-primary ratios measured by space-borne detectors such as
[PAMELA](https://pamela-web.web.roma2.infn.it/), [AMS-02](https://ams02.space/),
[DAMPE](http://dpnc.unige.ch/dampe/), or [CALET](https://calet.jp/).

Most of the code is described in [Evoli et al., 2019](https://arxiv.org/abs/1904.10220),
and it was expanded to include unstable (radioactive) isotopes in
[Evoli et al., 2020](https://arxiv.org/abs/1910.04113). Please cite these works if you
use CRAMS.

## What CRAMS computes

Given a set of source (injection) and transport parameters, CRAMS solves the
weighted-slab transport equation for every nucleus in [`data/crams_nucleilist.csv`](data/crams_nucleilist.csv)
and outputs their spectra at the solar position, as functions of both rigidity and
kinetic energy per nucleon. The physics it includes:

- **Diffusive transport** with a rigidity-dependent diffusion coefficient
  $D(R) \propto R^\delta$ and a low-rigidity break ($R_b$, $\Delta\delta$), plus
  **advective transport** set by the Alfvén speed $v_A$;
- **A finite diffusive halo** of half-height $H$, with the gas concentrated in an
  infinitely thin disk of surface density $\mu$;
- **Solar modulation** in the force-field approximation (potential $\phi$);
- **Nuclear interactions**: inelastic (destruction) cross sections and
  fragmentation (production) cross sections, both selectable between several
  published models (see [Data tables](#data-tables));
- **Secondary and tertiary production** chains down to the lightest nuclei,
  including radioactive isotopes (e.g. $^{10}\mathrm{Be}$).

The full list of input parameters — transport, numerics, cross-section model
selection, injection slopes, and per-element source abundances — is documented
inline in **[`examples/crams.ini`](examples/crams.ini)**, which doubles as the
reference example and is parsed by the test suite so it never drifts out of sync
with the code.

## Requirements

- A C++ compiler with C++14 support
- [CMake](https://cmake.org/) ≥ 3.14
- The [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library)

## Installation

```bash
mkdir build
cd build
cmake ..
make
make test   # optional: run the unit-test suite to check the build
```

This builds the `crams` executable (plus the unit-test binaries) inside `build/`.
`make test` runs the CTest suite and should report all tests passing.

## Running

CRAMS takes a single input file and writes the resulting spectra to an `output/`
directory created in the current working directory:

```bash
# from the repository root, after building
./build/crams examples/crams.ini

# add -q / --quiet to suppress the per-particle logging
./build/crams examples/crams.ini -q
```

- The input file uses one `key value` pair per line; keys are case-insensitive,
  underscores/spaces in keys are ignored, unknown keys and `#` comments are
  silently skipped, and the last occurrence of a key wins.
- Running with **no input file** uses the built-in defaults.
- Output is written to `output/<simname>_spectra_R_<id>.txt`, where `<simname>`
  is derived from the input file name and `<id>` is the integer `id` key.
  Columns are rigidity `R [GV]` followed by the flux of each element Z = 1…28.

➡️ **For the meaning, units, and default of every parameter, see the fully
commented [`examples/crams.ini`](examples/crams.ini).**

## Python steering

CRAMS can be embedded into Python via SWIG-generated interface. To do this, run

```bash
swig -c++ -python -Iinclude extension/crams.i
pip install .
```

Then, from Python interpreter you can use

```python
from crams import CramsRunner, PropagationParams, InjectionParams

runner = CramsRunner(
    inelastic_model="glauber",
    fragmentation_model="fluka4dragon",
    verbose=True,
)

propagation = PropagationParams(
    H_kpc=4.0,
    ...
)
injection = InjectionParams(
    abundances=[...],
    slopes=[4.4, 4.35, 4.3],
)
rigidity_spectra = runner.compute(propagation, injection)
```

Python extension can be validated for consistency with the main executable via

```shell
python tests/test_python.py
```

**TODO**: configure and build extension Python via CMake under optional flag;
integrate testing into the same system

## Data tables

CRAMS ships the nuclear cross-section grids and the nuclei list it needs in
[`data/`](data/). Each table is a CSV with a self-describing header (model name,
generating code/version, isotope source, and a `# reference:` field for the
citation). The cross-section model used at runtime is chosen with the
`inelastic_model` and `fragmentation_model` keys in the input file.

| File | Selected by | Model | Origin / citation |
|------|-------------|-------|-------------------|
| `crams_inelastic_tripathi99.csv` | `inelastic_model = tripathi99` | Tripathi (1999) | *TODO* |
| `crams_inelastic_glauber.csv` | `inelastic_model = glauber` | Glauber | *TODO* |
| `crams_inelastic_crosec.csv` | `inelastic_model = crosec` | CROSEC | *TODO* |
| `crams_fragmentation_fluka4dragon.csv` | `fragmentation_model = fluka4dragon` | FLUKA (as used in DRAGON) | *TODO* |
| `crams_fragmentation_usine_galprop17_opt12.csv` | `fragmentation_model = usine_galprop17_opt12` | GALPROP'17 GAL12 (via USINE) | *TODO* |
| `crams_fragmentation_usine_galprop17_opt22.csv` | `fragmentation_model = usine_galprop17_opt22` | GALPROP'17 GAL22 (via USINE) | *TODO* |
| `crams_fragmentation_usine_webber03+coste12.csv` | `fragmentation_model = usine_webber03_coste12` | Webber (2003) + Coste (2012) (via USINE) | *TODO* |
| `crams_nucleilist.csv` | (always) | Nuclei list: Z, A, half-life, ISM isotopic fractions | *TODO* |

> **⚠️ Provenance is incomplete.** The cross-section tables are generated by an
> external tool (`XS4GCR`) and currently carry `reference: TODO` in their headers.
> Before publishing results, the underlying datasets/codes **must** be cited.

**How this section is meant to be maintained.** The per-file `# reference:`
header inside each CSV is the *single source of truth* for provenance; the table
above is just a human-readable index. To add or update a table:

1. Fill the `# reference:` (and, where useful, `# url:` / `# doi:`) header line
   in the CSV with the citation for the dataset or generating code.
2. Add/refresh the corresponding row in the table above, replacing *TODO* with a
   short citation (and a link to the paper or repository).
3. If a table is regenerated, keep the `# code_version:`/`# created:` header
   fields so the grid stays traceable to the tool that produced it.

## MCMC tools (`mcmc/`)

The [`mcmc/`](mcmc/) folder is an **auxiliary toolkit**, not part of the core
code: a set of Python scripts that drive the `crams` binary to fit transport and
injection parameters to cosmic-ray data (e.g. via `emcee`/`iminuit`), locate the
best-fit point, and plot the resulting fluxes, ratios, isotope predictions, and
parameter posteriors.
See the scripts and their docstrings in that directory for usage.

## Versions
* **Version 2.0 (13/06/2026):** New Crank-Nicolson integration of the transport equation, selectable inelastic/fragmentation cross-section models, and an MCMC fitting toolkit (`mcmc/`)
* **Version 1.3 (24/02/2021):** Extended to intermediate-mass and heavy nuclei up to iron, confronting the new AMS-02 measurements ([2102.12576](https://arxiv.org/abs/2102.12576))
* **Version 1.2 (20/11/2019):** Added unstable (radioactive) isotopes such as $^{10}\mathrm{Be}$ to constrain the cosmic-ray residence time and halo size from the AMS-02 beryllium data ([1910.04113](https://arxiv.org/abs/1910.04113))
* **Version 1.1 (07/02/2019):** Stable nuclei ($H–O$), with injection and diffusion parameters fit to the AMS-02 fluxes and the B/C ratio ([1904.10220](https://arxiv.org/abs/1904.10220))
* **Version 1.0 (25/01/2019):** Release version.
