# LunarResponse

[中文说明](README_zh.md)

**Lunar surface response functions from an interior model, in one command.**
Given a lunar (or planetary) interior 1D model, this package refines the model,
computes its free oscillations with a modified [MINEOS](https://github.com/geodynamics/mineos)
kernel, and sums the normal modes into frequency-dependent surface displacement
response functions — both A-type (Dyson forcing) and B-type (tidal
forcing) — ready for gravitational-wave / tidal response studies.

```bash
python lunar_response.py all --model examples/MoonModelSimp.txt --name demo --out run_demo
```

That single command produces the eigenmode catalogue and the response curves
(`response_functions.txt` + two PNG figures) in `run_demo/`.

---

## Features

- **End-to-end pipeline**: simple layered model → refined MINEOS card model →
  eigenfrequencies & eigenfunctions (spheroidal and toroidal) → response curves.
  The stages are also available separately (`model` / `eigen` / `response`).
- **Enlarged MINEOS**: the radial-knot limit is raised from the stock `mk=350`
  to `mk=50000`, so models interpolated to tens of thousands of layers can be
  solved (verified up to 55001 layers with `mk=60000`). The physics is untouched;
  results are bit-identical to stock MINEOS.
- **Prebuilt Windows binaries** in `bin/` — no compiler needed for normal use.
- **Robust mode-table parsing** (no fragile fixed-line skipping), fully
  vectorized numerics (a typical l=2 run with ~350 modes finishes in ~30 s).
- **Fluid layers** (e.g. a fluid core) are detected automatically and can also
  be set explicitly.

## Requirements

- Python 3 with `numpy` (`matplotlib` only for plotting).
- Windows: use the prebuilt executables in `bin/` as-is.
- Other platforms, or to change the layer limit: a `gfortran` compiler
  (on Windows you can install a portable one with `python install_mingw.py`,
  no admin rights required).

## Quick start

```bash
# Full pipeline (defaults reproduce a homogeneous-Moon l=2 benchmark)
python lunar_response.py all --model examples/MoonModelSimp.txt --name demo --out run_demo

# Individual stages
python lunar_response.py model    --model examples/MoonModelSimp.txt --name demo --out run_demo
python lunar_response.py eigen    --model-file run_demo/demo.txt --out run_demo
python lunar_response.py response --folder run_demo --code demo_code.txt
```

The program auto-discovers the MINEOS executables: `--exe-dir` >
`$LUNAR_MINEOS_DIR` > `bin/` next to the script.

### Input model format

Six columns, ordered from the **surface down to the centre**:

```
r(km)  vp(km/s)  vs(km/s)  rho(g/cm^3)  Q_kappa  Q_mu
```

Set `vs = 0` for fluid layers (the fluid block is detected and written into the
MINEOS header as `nic/noc` automatically).

## Key options

| Option | Meaning | Default |
|---|---|---|
| `--dr-km` / `--layers` | Interpolation accuracy: target max layer thickness (km) **or** target total layers (choose one) | legacy `--mult/--mj/--crit` scheme |
| `--nic`, `--noc` | Fluid-layer indices (ICB solid side / CMB fluid side) | auto-detected (`vs=0, rho>0`) |
| `--nmin`, `--nmax` | Radial-order range of modes (need not start from n=0) | 0, 400 |
| `--wmin`, `--wmax` | Frequency range (mHz) | 0.1, 300 |
| `--wgrav` | Frequency (mHz) above which self-gravity is neglected (~3× faster) | 0.5 |
| `--jcom` | Mode type: 1 radial, 2 toroidal, 3 spheroidal, 4 inner-core toroidal | 3 |
| `--lmin`, `--lmax` | Angular-degree range | 2, 2 |
| `--eps` | Runge–Kutta integration accuracy (toroidal: do not go below ~2e-4) | 1e-8 |
| `--save-format` | `compact` (5 columns) or `full` (7 columns, incl. raw \|T_B\|) | compact |
| `--static-rad`, `--static-hor` | B-A transform coefficients (×R; l=2 values are 0.5 / 0.25) | 0.5, 0.25 |

Run `python lunar_response.py <command> -h` for the full list.

## Output convention: A-type vs B-type

- **A-type `T_A`** — surface displacement response to a *direct body force*
  (Dyson forcing).
- **B-type `T_B`** — response to a *tidal force*. Its transform coefficients are
  `0.5·R` (radial) and `0.25·R` (horizontal) for l=2; 
  `|T_B − transform coefficient|` is comparable in magnitude to the A-type response. See our paper: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.064092

`compact` format: `freq, T_A_r, T_A_h, |T_B_r−0.5R|, |T_B_h−0.25R|`.
`full` format adds the unsubtracted `|T_B_r|`, `|T_B_h|`.
Headers in the output files state these definitions explicitly.

## Rebuilding MINEOS / changing the layer limit

```bash
python install_mingw.py            # one-time, if no gfortran is available (Windows)
python build_mineos.py             # rebuild with default mk=50000 into ./bin
python build_mineos.py --mk 100000 # raise the layer limit, one command
```

`build_mineos.py` patches every `parameter (mk=...)` in a temporary copy (the
sources in `src/` are never modified), compiles with `-O2` and semi-static
linking, bundles the runtime DLLs, and runs a smoke test.
**Warning:** if a model has more layers than `mk`, the solver *hangs* instead of
reporting an error — always keep `layers < mk`.

## Design notes

- **Normal-mode summation.** Eigenfrequencies and eigenfunctions come from
  MINEOS (`minos_bran` shooting method with Runge–Kutta integration; `eigcon`;
  `eigen2asc`), using the normalization of Woodhouse & Dahlen (1978).
  Mode coupling coefficients (tidal loading and Dyson loading, the latter via
  integration by parts of the shear-modulus gradient) are evaluated by Simpson
  quadrature on the model grid, then summed over modes with damped-oscillator
  propagators `1/(ω_n² − ω² + iωω_n/Q_n)`.
- **Source modifications vs upstream MINEOS** (kept minimal on purpose):
  `mk=350 → 50000` in three files; `fdb_eigen` helper routines inlined into
  `eigcon.f` for single-file compilation; one `loctime()` call commented out.
  Numerical results are identical to upstream.
- Verified bit-for-bit against the original notebook pipeline (356 spheroidal
  modes, 567 toroidal modes, and the final response curves all match exactly).

## Repository layout

```
lunar_response.py    end-to-end driver (model → eigenmodes → response curves)
build_mineos.py      one-command MINEOS build with adjustable layer limit mk
install_mingw.py     portable gfortran installer (Windows, no admin needed)
src/                 modified MINEOS Fortran sources (see Design notes)
bin/                 prebuilt Windows executables + runtime DLLs
examples/            example input model (homogeneous Moon)
LICENSE              GPL v2 (inherited from MINEOS)
```

## License

The Fortran sources are derived from MINEOS and remain under the
**GNU General Public License v2** (see `LICENSE`). The Python driver is
distributed under the same terms.
