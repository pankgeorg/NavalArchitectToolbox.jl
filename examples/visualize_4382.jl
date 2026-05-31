#!/usr/bin/env julia
#
# Visualise DTMB 4382 (D=6 m, Dh=1.2 m): radial section distributions,
# the 3-D blade geometry built from the table, and the VLM camber grid.
# CairoMakie (headless, no xvfb needed).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.add(["CairoMakie", "StaticArrays"]; io=devnull)

using NavalArchitectToolbox, CairoMakie, Printf
CairoMakie.activate!()

const D, Dh = 6.0, 1.2
t = dtmb4382
OUT = joinpath(@__DIR__, "..", "runs"); mkpath(OUT)

# ── Figure 1: radial section distributions (the table) ────────────────
fig1 = Figure(size=(1100, 640))
specs = [(:cD,"c/D",1,1), (:PD,"P/D",1,2), (:tc,"tmax/c",1,3),
         (:fc,"fmax/c",2,1), (:skew,"θs [deg]",2,2)]
for (sym, lab, i, j) in specs
    ax = Axis(fig1[i,j]; xlabel="r/R", ylabel=lab, title=lab)
    y = getfield(t, sym)
    lines!(ax, t.rR, y, color=:steelblue)
    scatter!(ax, t.rR, y, color=:tomato, markersize=7)
end
# pitch angle β(r/R) derived
axβ = Axis(fig1[2,3]; xlabel="r/R", ylabel="β [deg]", title="pitch angle β = atan(P/2πr)")
rRf = range(t.rR[1], t.rR[end], length=100)
lines!(axβ, rRf, [rad2deg(pitch_angle(t, x)) for x in rRf], color=:seagreen)
Label(fig1[0, :], "DTMB 4382 — radial section data (Z=$(t.Z))", fontsize=18)
save(joinpath(OUT, "radial_4382.png"), fig1)
println("wrote radial_4382.png")

# ── Figure 2: expanded blade outline (chord vs radius, with skew) ──────
fig2 = Figure(size=(560, 720))
ax2 = Axis(fig2[1,1]; xlabel="skewed chordwise offset [m]", ylabel="r [m]",
           title="DTMB 4382 expanded outline (LE/TE, skew)", aspect=DataAspect())
rRf = range(t.rR[1], t.rR[end], length=60)
R = D/2
# leading & trailing edge in the (mid-chord-tangential, r) plane incl. skew
le = Float64[]; te = Float64[]; rr = Float64[]
for x in rRf
    g = NavalArchitectToolbox.dimensional(t, x, D)
    # tangential offset of LE/TE relative to the shaft, including skew arc
    arc_skew = g.θs * g.r
    push!(rr, g.r)
    push!(le, arc_skew - g.c/2 * cos(g.β))   # LE
    push!(te, arc_skew + g.c/2 * cos(g.β))   # TE
end
band!(ax2, Point2f.(le, rr), Point2f.(te, rr); color=(:steelblue, 0.4))
lines!(ax2, le, rr, color=:navy); lines!(ax2, te, rr, color=:navy)
save(joinpath(OUT, "outline_4382.png"), fig2)
println("wrote outline_4382.png")

# ── Figure 3: 3-D blade geometry + VLM camber grid ────────────────────
fig3 = Figure(size=(1200, 560))
splat(g) = (getindex.(g,1), getindex.(g,2), getindex.(g,3))

# (a) all Z blades, solid surfaces
ax3 = Axis3(fig3[1,1]; aspect=:data, azimuth=0.7π, elevation=0.18π,
            title="DTMB 4382 — $(t.Z) blades (D=$D m, Dh=$Dh m)")
for b in 1:t.Z
    up, lo = blade_surface(t, D, Dh; nc=30, ns=30, blade=b)
    for skin in (up, lo)
        X, Y, Z = splat(skin)
        surface!(ax3, X, Y, Z; colormap=:viridis, color=Z, shading=NoShading)
    end
end

# (b) one blade: VLM camber grid (the lifting surface for VortexLattice)
ax4 = Axis3(fig3[1,2]; aspect=:data, azimuth=0.7π, elevation=0.18π,
            title="VLM camber grid (1 blade, nc=10×ns=16)")
grid = vlm_camber_grid(t, D, Dh; nc=10, ns=16)
Xc, Yc, Zc = splat(grid)
wireframe!(ax4, Xc, Yc, Zc; color=:tomato, linewidth=0.6)
surface!(ax4, Xc, Yc, Zc; color=fill(0.0, size(Zc)), colormap=[:lightblue,:lightblue],
         transparency=true, alpha=0.3, shading=NoShading)
save(joinpath(OUT, "blade3d_4382.png"), fig3)
println("wrote blade3d_4382.png")
