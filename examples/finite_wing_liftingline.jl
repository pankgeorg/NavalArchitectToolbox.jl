#!/usr/bin/env julia
#
# Finite-wing lifting-line VLM (`liftingline_vlm`) — a symmetric tapered wing
# swept through angle of attack, leading-edge sweep, and geometric twist. Shows
# the pure-Julia Weissinger horseshoe lattice and cross-checks its lift against
# the multi-chordwise `Wing` (VortexLattice.jl). Headless CairoMakie.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.add(["CairoMakie", "Printf"]; io=devnull)

using NavalArchitectToolbox, CairoMakie, Printf

const CMID, CTIP, SPAN = 6.0, 3.5, 17.0      # root(mid) chord, tip chord, span
wing(; kw...) = liftingline_vlm(; chord_root=CMID, chord_tip=CTIP, span=SPAN, kw...)

# ---- (a) AoA sweep, baseline + with washout twist; cross-check vs Wing ----
aoa = 0.0:2.5:20.0
base = [wing(; alpha=a) for a in aoa]
twst = [wing(; alpha=a, twist_root=2.0, twist_tip=-2.0) for a in aoa]
wingv = [wing_forces(Wing(; chord_root=CMID, chord_tip=CTIP, span=SPAN, ns=40, nc=8),
                     deg2rad(a), 1.0) for a in aoa]

@printf("AoA   CL(ll)   CL(Wing)   CDi(ll)    e(ll)\n")
for (k,a) in enumerate(aoa)
    @printf("%4.1f  %7.4f  %7.4f   %8.5f  %5.3f\n",
            a, base[k].CL, wingv[k].CL, base[k].CDi, base[k].e)
end

fig = Figure(size=(960, 400))
ax1 = Axis(fig[1,1]; title="C_L vs AoA", xlabel="AoA (deg)", ylabel="C_L")
lines!(ax1, collect(aoa), [r.CL for r in base]; label="liftingline_vlm", linewidth=2)
lines!(ax1, collect(aoa), [r.CL for r in twst]; label="+ twist 2°/-2°", linestyle=:dash, linewidth=2)
scatter!(ax1, collect(aoa), [r.CL for r in wingv]; label="Wing (VortexLattice)", color=:black, markersize=7)
axislegend(ax1; position=:lt)
ax2 = Axis(fig[1,2]; title="Spanwise circulation Γ(y), AoA=10°", xlabel="y (m)", ylabel="Γ")
rL = wing(; alpha=10.0, N=120)
lines!(ax2, rL.y, rL.Γ; linewidth=2)
save(joinpath(@__DIR__, "finite_wing_liftingline.png"), fig; px_per_unit=2)
println("wrote examples/finite_wing_liftingline.png")
