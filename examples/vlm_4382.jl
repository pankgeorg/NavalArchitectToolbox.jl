#!/usr/bin/env julia
#
# Feed the NavalArchitectToolbox camber grid into VortexLattice.jl and
# solve for the blade *loading* — the spanwise bound-circulation
# distribution — at a near-pitch-matched advance ratio. Rotation enters
# through VortexLattice's Freestream angular-velocity vector (shaft = x).
#
# Honest scope: linear, inviscid, fixed-helix-wake VLM. The *loading
# shape* (bell-shaped, zero at hub and tip) is meaningful; absolute
# CT/CQ are over-predicted for a heavily-loaded propeller (no wake
# contraction, no viscous torque ⇒ η>1) — use BEM (CCBlade) or the
# actuator-disk models for quantitative thrust. This demonstrates the
# NAT → VortexLattice pipeline.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.add(["VortexLattice", "CairoMakie", "StaticArrays"]; io=devnull)

using NavalArchitectToolbox, VortexLattice, CairoMakie, Printf, StaticArrays

const D, Dh, Z = 6.0, 1.2, 5
const R = D/2
const NC, NS = 8, 16
t = dtmb4382
OUT = joinpath(@__DIR__, "..", "runs"); mkpath(OUT)

# NAT camber grid (shaft=z) → VortexLattice grids (shaft=x), Z blades.
function vl_grids(nc, ns)
    cg = vlm_camber_grid(t, D, Dh; nc=nc, ns=ns, meanline=NACAMeanLine(0.8), tip_rR=0.95)
    ni, nj = size(cg)
    base = Array{Float64}(undef, 3, ni, nj)
    for j in 1:nj, i in 1:ni
        p = cg[i, j]
        base[1,i,j] = p[3]; base[2,i,j] = p[1]; base[3,i,j] = p[2]   # z→x, x→y, y→z
    end
    grids = [copy(base)]
    for b in 1:Z-1
        θ = 2π*b/Z; c, s = cos(θ), sin(θ); g = copy(base)
        for j in 1:nj, i in 1:ni
            y = g[2,i,j]; z = g[3,i,j]; g[2,i,j] = y*c - z*s; g[3,i,j] = y*s + z*c
        end
        push!(grids, g)
    end
    return grids
end

# Solve a SINGLE blade. The full 5-blade AIC is singular here: near the
# hub the chord (~1.04 m) exceeds the inter-blade gap (~0.75 m), so the
# blades overlap and a control point lands on another blade's vortex —
# a proper multi-blade propeller VLM needs hub treatment / finite-core
# wake. The single-blade loading *shape* is the meaningful output;
# multi-blade thrust ≈ Z× (neglecting blade-blade induction).
grids  = vl_grids(NC, NS)
system = System([grids[1]])
Sref   = π*R^2

J  = 1.0
n  = 1.0; Ω = 2π*n; Va = J*n*D
ref = Reference(Sref, R, 2R, [0.0,0.0,0.0], Va)
fs  = Freestream(Va, 0.0, 0.0, [Ω, 0.0, 0.0])
steady_analysis!(system, ref, fs; symmetric=false)
CF, CM = body_forces(system; frame=Body())
@printf("J=%.2f  per-blade CF=[%.4f, %.4f, %.4f]  CM1=%.4f  (×Z=%d, qualitative)\n",
        J, CF[1], CF[2], CF[3], CM[1], Z)

# Spanwise bound circulation of blade 1: Γ panels are column-major (nc,ns)
# per surface; sum chordwise → strip circulation Γ(r).
Γ = system.Γ
Γb1 = reshape(Γ[1:NC*NS], NC, NS)
Γspan = vec(sum(Γb1; dims=1))                 # length NS
# strip radial centres (cosine-clustered hub→tip), matching vl_grids
rR_hub = Dh/D; tip_rR = 0.95
rRedges = [rR_hub + (tip_rR-rR_hub)*(1 - cos(π*i/(2NS))) for i in 0:NS]
rRc = (rRedges[1:end-1] .+ rRedges[2:end]) ./ 2
Gnorm = Γspan ./ (Va * D)                      # non-dimensional circulation

fig = Figure(size=(760, 460))
ax = Axis(fig[1,1]; xlabel="r/R", ylabel="G = Γ/(Va·D)",
          title="DTMB 4382 spanwise loading (VLM, NACA a=0.8, J=$(J))")
lines!(ax, rRc, Gnorm; color=:steelblue, linewidth=2)
scatter!(ax, rRc, Gnorm; color=:tomato, markersize=6)
hlines!(ax, [0]; color=:grey, linestyle=:dash)
save(joinpath(OUT, "loading_4382.png"), fig)
println("wrote loading_4382.png")
@printf("peak G = %.4f at r/R = %.2f ; G(hub)=%.4f G(tip)=%.4f\n",
        maximum(Gnorm), rRc[argmax(Gnorm)], Gnorm[1], Gnorm[end])
