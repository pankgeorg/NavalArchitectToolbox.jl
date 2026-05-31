#!/usr/bin/env julia
#
# Visualise the blade FORCES from the DTMB 4382 key-blade VLM: the
# pressure-normal lift on each panel, and its decomposition into the axial
# component (thrust) and the tangential component (which produces shaft torque).
#
#   julia +1.12 --project=papers/anevlavi_belibassakis_2023 \
#       papers/anevlavi_belibassakis_2023/forces.jl

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate(; io=devnull)

using NavalArchitectToolbox, CairoMakie, Printf, LinearAlgebra
CairoMakie.activate!(type="png")

const D, Dh = 6.0, 1.2
const OUT = joinpath(@__DIR__, "figures"); mkpath(OUT)

nc, ns = 8, 14
r  = openwater_vlm(dtmb4382, D, Dh, 0.889; nc=nc, ns=ns)
cg = vlm_camber_grid(dtmb4382, D, Dh; nc=nc, ns=ns)        # camber surface

P3(p) = Point3f(p[1], p[2], p[3])
V3(p) = Vec3f(p[1], p[2], p[3])
pts   = [P3(p) for p in r.cp]
Ftot  = r.force                                            # pressure-normal lift
Fthr  = [Vec3f(0, 0, f[3]) for f in Ftot]                  # axial → thrust
function Ftan(p, f)                                        # tangential → torque
    rr = hypot(p[1], p[2]); thx, thy = -p[2]/rr, p[1]/rr   # θ̂ (in-plane)
    d = f[1]*thx + f[2]*thy
    Vec3f(d*thx, d*thy, 0)
end
Ftq = [Ftan(r.cp[k], Ftot[k]) for k in eachindex(Ftot)]

splat(g) = (getindex.(g,1), getindex.(g,2), getindex.(g,3))
Xs, Ys, Zs = splat(cg)
camber!(ax) = surface!(ax, Xs, Ys, Zs; color=fill(0.0, size(Zs)),
                       colormap=[:gray85, :gray85], shading=NoShading, transparency=true)

LS = 0.7
fig = Figure(size=(1320, 720))

# (a) the lift: total pressure-normal force per panel
axA = Axis3(fig[1,1]; aspect=:data, azimuth=0.30π, elevation=0.26π,
            title="Pressure-normal lift  (Δp·A·n̂ per panel)")
camber!(axA)
arrows!(axA, pts, [V3(f) for f in Ftot]; lengthscale=LS, color=:seagreen,
        linewidth=0.02, arrowsize=Vec3f(0.06,0.06,0.10))

# (b) decomposition: axial (thrust) + tangential (torque)
axB = Axis3(fig[1,2]; aspect=:data, azimuth=0.30π, elevation=0.26π,
            title="Decomposed:  thrust (axial, blue) + torque (tangential, red)")
camber!(axB)
arrows!(axB, pts, Fthr; lengthscale=LS, color=:royalblue,
        linewidth=0.02, arrowsize=Vec3f(0.06,0.06,0.10))
arrows!(axB, pts, Ftq;  lengthscale=LS, color=:crimson,
        linewidth=0.02, arrowsize=Vec3f(0.06,0.06,0.10))

Label(fig[0, :], @sprintf("DTMB 4382 blade loading at J=0.889   —   KT=%.3f,  10·KQ=%.3f,  η=%.3f",
                          r.KT, 10r.KQ, r.η); fontsize=18)
save(joinpath(OUT, "forces_4382.png"), fig)
println("wrote figures/forces_4382.png")
