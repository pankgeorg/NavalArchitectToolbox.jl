#!/usr/bin/env julia
#
# Visualise the vortex-lattice DISCRETIZATION and the VORTEX SYSTEM of the
# DTMB 4382 key-blade VLM: the bound vortex-ring lattice on the mean camber
# surface, the {1/4, 3/4} lifting lines and control points, and the helical
# trailing wake spiralling downstream behind all five blades.
#
#   julia +1.12 --project=papers/anevlavi_belibassakis_2023 \
#       papers/anevlavi_belibassakis_2023/discretization.jl

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate(; io=devnull)

using NavalArchitectToolbox, CairoMakie
const N = NavalArchitectToolbox
CairoMakie.activate!(type="png")

const D, Dh = 6.0, 1.2
const OUT = joinpath(@__DIR__, "figures"); mkpath(OUT)

rotz(p, ψ) = (c = cos(ψ); s = sin(ψ); (c*p[1] - s*p[2], s*p[1] + c*p[2], p[3]))
P3(p) = Point3f(p[1], p[2], p[3])

# closed ring loops (c1-c2-c3-c4-c1) for a range of elements, NaN-separated
function ringloops(elems, rng, ψ)
    pts = Point3f[]
    for k in rng
        c = elems[k].c
        for q in (1, 2, 3, 4, 1); push!(pts, P3(rotz(c[q], ψ))); end
        push!(pts, Point3f(NaN, NaN, NaN))
    end
    pts
end

# only the two streamwise (downstream) edges of each wake ring → helical filaments
function wakeedges(elems, rng, ψ)
    segs = Point3f[]
    for k in rng
        c = elems[k].c
        push!(segs, P3(rotz(c[1], ψ)), P3(rotz(c[2], ψ)))   # left  edge
        push!(segs, P3(rotz(c[4], ψ)), P3(rotz(c[3], ψ)))   # right edge
    end
    segs
end

# ─────────────────────────────────────────────────────────────────────────
# Figure A — discretization of one blade ({1/4,3/4} rule)
# ─────────────────────────────────────────────────────────────────────────
vlmA = N._build_vlm(dtmb4382, D, Dh, 0.889; nc=8, ns=12, Kw=1, dψ=deg2rad(15.0))
bound = 1:vlmA.nb

figA = Figure(size=(940, 860))
axA = Axis3(figA[1,1]; aspect=:data, azimuth=0.32π, elevation=0.30π,
            title="DTMB 4382 — vortex-ring lattice, lifting lines & control points (1 blade)")
# bound vortex-ring lattice (camber surface)
lines!(axA, ringloops(vlmA.elems, bound, 0.0); color=(:steelblue, 0.6), linewidth=1.2)
# 1/4-chord bound (lifting) segments — where Γ lives
seg = Point3f[]
for k in bound; a, b = vlmA.bseg[k]; push!(seg, P3(a), P3(b)); end
linesegments!(axA, seg; color=:navy, linewidth=2.4)
# 3/4-chord control points — where the no-penetration BC is enforced
scatter!(axA, [P3(p) for p in vlmA.cp[bound]]; color=:tomato, markersize=8)
# surface normals at the control points (a sparse subset, short arrows)
sub = 1:3:vlmA.nb
arrows!(axA, [P3(p) for p in vlmA.cp[sub]], [Vec3f(n...) for n in vlmA.nrm[sub]];
        lengthscale=0.28, color=(:seagreen, 0.85), linewidth=0.018,
        arrowsize=Vec3f(0.05, 0.05, 0.08))
# legend proxies
elems_leg = [LineElement(color=:steelblue, linewidth=2),
             LineElement(color=:navy, linewidth=3),
             MarkerElement(color=:tomato, marker=:circle, markersize=9),
             LineElement(color=:seagreen, linewidth=2)]
Legend(figA[2,1], elems_leg,
       ["vortex-ring panels", "1/4-chord bound vortices (Γ)",
        "3/4-chord control points", "surface normals (BC)"],
       orientation=:horizontal, framevisible=false)
save(joinpath(OUT, "discretization_4382.png"), figA)
println("wrote figures/discretization_4382.png")

# ─────────────────────────────────────────────────────────────────────────
# Figure B — the full vortex system: all blades + helical trailing wake
# ─────────────────────────────────────────────────────────────────────────
vlmB = N._build_vlm(dtmb4382, D, Dh, 0.889; nc=6, ns=10, Kw=44, dψ=deg2rad(12.0))
boundB = 1:vlmB.nb
wakeB  = (vlmB.nb+1):length(vlmB.elems)

figB = Figure(size=(1000, 900))
axB = Axis3(figB[1,1]; aspect=:data, azimuth=0.62π, elevation=0.10π,
            title="DTMB 4382 — bound vortex rings + geometric-pitch helical wake (Z=5)")
for b in 0:dtmb4382.Z-1
    ψ = 2π*b/dtmb4382.Z
    lines!(axB, ringloops(vlmB.elems, boundB, ψ); color=:navy, linewidth=1.3)
    linesegments!(axB, wakeedges(vlmB.elems, wakeB, ψ);
                  color=(:steelblue, 0.35), linewidth=0.7)
end
save(joinpath(OUT, "vortex_system_4382.png"), figB)
println("wrote figures/vortex_system_4382.png")

println("done — 2 figures in $OUT")
