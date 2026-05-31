# Anevlavi & Belibassakis (2023) — DTMB 4382 open-water reproduction

Reproduction of the open-water vortex-lattice results from

> E. Anevlavi, S. Zafeiris, G. Papadakis, K. Belibassakis,
> *A Low-Cost Vortex-Lattice Method for the Preliminary Design of Marine
> Propellers with Tip-Rake Reformation*, **J. Mar. Sci. Eng. 2023, 11, 2179.**
> doi:10.3390/jmse11112179

using the key-blade vortex-ring VLM in [`NavalArchitectToolbox`](../..)
(`openwater_vlm`, `src/vlm.jl`).

This is a **self-contained Julia 1.12 environment** (own `Project.toml`,
`NavalArchitectToolbox` brought in by path via `[sources]`) so the heavy
plotting stack (CairoMakie) stays out of the package's own dependencies.

## Run

```bash
julia +1.12 --project=papers/anevlavi_belibassakis_2023 \
    papers/anevlavi_belibassakis_2023/reproduce.jl
```

It prints the open-water table and writes four figures to `figures/`.

## Result (J = 0.889)

| coefficient | VLM (this work) | experiment | paper VLM |
|---|---|---|---|
| KT          | 0.210 | 0.208 | 0.208 |
| 10·KQ       | 0.452 | 0.445 | — |
| η           | 0.659 | 0.661 | — |

Mesh- and wake-converged. Calibration `C_LE = 0.80`, `C_Drag = 0.010`
(the paper quotes `C_LE = 0.90`, `C_Drag = 0.0050` with its exact NACA
66-mod sections; this package uses an exact NACA a=0.8 mean line and a
NACA-4-digit thickness stand-in, so it needs a slightly stronger
correction).

## Figures

- `openwater_4382.png` — KT, 10·KQ, η vs J; VLM curves + experimental
  points (paper Fig. 7).
- `convergence_4382.png` — coefficients vs mesh resolution converging to
  the experimental values.
- `loading_4382.png` — spanwise bound-circulation distribution and the
  chordwise Δc_p field on the camber surface.
- `dcp_blades_4382.png` — Δc_p on all five blades viewed down the shaft
  (parallels paper Fig. 8).

Run `discretization.jl` for two more (the method itself):

- `discretization_4382.png` — one blade: the vortex-ring lattice on the mean
  camber surface, the 1/4-chord bound (lifting) vortices, the 3/4-chord
  control points, and the surface normals (the {1/4, 3/4} rule).
- `vortex_system_4382.png` — all five blades' bound rings plus the
  geometric-pitch helical trailing wake spiralling ~2 diameters downstream.

Run `forces.jl` for the blade loading:

- `forces_4382.png` — the pressure-normal lift on each panel, and its
  decomposition into the axial component (thrust) and the tangential
  component (shaft torque).

## Interactive notebook

`reproduction.pluto.jl` is a **self-contained [Pluto](https://plutojl.org)
notebook** that walks through the whole story — table → 3-D blade → vortex
lattice → helical wake → forces → open-water curve — with an interactive
advance-ratio slider. It depends only on registered packages (StaticArrays,
CairoMakie, PlutoUI); the geometry and solver are written out inline, so it
needs nothing from this package. Open it with:

```julia
using Pluto; Pluto.run()   # then open reproduction.pluto.jl
```

## What it took to reproduce (the three subtleties)

1. **Wake = the propeller's *geometric*-pitch cylindrical helix** (winds
   with the rotation, the downstream extension of the blade surface). A
   flow-aligned helix at the advance pitch winds the opposite way, grazes
   the control points, and makes the influence matrix singular.
2. **Pressure-normal forces, not Kutta–Joukowski** — Δp·A·n̂ along the
   surface normal (paper eqs 5–7). The geometric force direction gives the
   correct η; Γ×V gives η > 4.
3. **Handedness consistency** — the relative inflow must run LE→TE with the
   camber on the suction side; getting it backwards inflates the effective
   angle of attack ~6×.
