#!/usr/bin/env julia
#
# What does "skew" mean? Two figures, built from the real DTMB modified-series
# geometry (4381 = 0° skew, 4382 = 36° skew — identical in every other respect).
# Skew is the circumferential (in-disk-plane) sweep of the blade's reference
# line as you move out along the radius; it does NOT change chord, pitch,
# thickness or camber. CairoMakie, headless.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.add(["CairoMakie"]; io=devnull)

using NavalArchitectToolbox, CairoMakie, Printf
using NavalArchitectToolbox: dimensional
CairoMakie.activate!()

const D, Dh = 6.0, 1.2
OUT = joinpath(@__DIR__, "..", "runs"); mkpath(OUT)
rRs = range(Dh/D, 1.0; length=60)          # hub → tip

# Blade outline as SEEN LOOKING DOWN THE SHAFT (project onto the disk plane).
# blade_section_point returns (X, Y, Z) with Z the shaft axis, so (X,Y) is the
# disk-plane projection. We trace leading edge (xc=0) up the span and trailing
# edge (xc=1) back down, giving a closed blade silhouette. The mid-chord line
# (xc=0.5) is the conventional skew reference line.
function outline(tab, blade)
    le = [blade_section_point(tab, rR, 0.0, D; blade) for rR in rRs]
    te = [blade_section_point(tab, rR, 1.0, D; blade) for rR in rRs]
    xs = vcat(getindex.(le,1), reverse(getindex.(te,1)))
    ys = vcat(getindex.(le,2), reverse(getindex.(te,2)))
    return xs, ys
end
# Skew reference line = the circumferential angle θs(r) of the mid-chord,
# which is exactly how skew is *defined* (θ = blade datum + θs(r), at radius r).
function midline(tab, blade)
    θ0 = 2π*(blade-1)/tab.Z
    pts = [(g = dimensional(tab, rR, D); (g.r*cos(θ0+g.θs), g.r*sin(θ0+g.θs))) for rR in rRs]
    (getindex.(pts,1), getindex.(pts,2))
end

# ── Figure 1: the concept — one blade, 4381 vs 4382 overlaid in the disk plane
fig1 = Figure(size=(760, 780))
ax = Axis(fig1[1,1]; aspect=DataAspect(),
    title="What skew is — one blade, looking down the shaft axis",
    subtitle="DTMB 4381 (0° skew) vs 4382 (36° skew): same chord/pitch/thickness, blade swept circumferentially",
    xlabel="disk-plane X [m]", ylabel="disk-plane Y [m]")

# blade 1 with skew=0 lies essentially along +X (θ = θs + tiny pitch-wrap term);
# draw that as the radial reference line the skew is measured from.
lines!(ax, [Dh/2, D/2], [0, 0]; color=(:black,0.35), linestyle=:dash, label="radial reference")

x1,y1 = outline(dtmb4381,1); x2,y2 = outline(dtmb4382,1)
poly!(ax, Point2f.(x1,y1); color=(:steelblue,0.25), strokecolor=:steelblue, strokewidth=2)
poly!(ax, Point2f.(x2,y2); color=(:tomato,0.25),   strokecolor=:tomato,    strokewidth=2)
mx1,my1 = midline(dtmb4381,1); mx2,my2 = midline(dtmb4382,1)
lines!(ax, mx1, my1; color=:steelblue, linewidth=3, label="4381 mid-chord (skew=0)")
lines!(ax, mx2, my2; color=:tomato,    linewidth=3, label="4382 mid-chord (skewed)")

# tip skew = angular sweep of the mid-chord line from hub to tip (from +X axis)
θhub = atan(my2[1],  mx2[1]); θtip_a = atan(my2[end], mx2[end])
skewdeg = rad2deg(θtip_a - θhub)
scatter!(ax, [mx2[end]], [my2[end]]; color=:tomato, markersize=11)
text!(ax, mx2[end], my2[end]; text=@sprintf("tip skew ≈ %.0f°", abs(skewdeg)),
      color=:tomato, fontsize=15, offset=(8,8), align=(:left,:bottom))
axislegend(ax; position=:lb, framevisible=true)
save(joinpath(OUT, "skew_concept.png"), fig1; px_per_unit=2)
println("wrote runs/skew_concept.png  (4382 tip mid-chord skewed $(round(abs(skewdeg),digits=1))° from radial)")

# ── Figure 2: the assembled 5-blade propellers, side by side
fig2 = Figure(size=(1150, 600))
for (k,(tab,name,col)) in enumerate(((dtmb4381,"DTMB 4381 — unskewed",:steelblue),
                                     (dtmb4382,"DTMB 4382 — 36° skew",:tomato)))
    local ax = Axis(fig2[1,k]; aspect=DataAspect(), title=name,
              xlabel="X [m]", ylabel = k==1 ? "Y [m]" : "")
    # hub circle
    φ = range(0,2π;length=120); lines!(ax, (Dh/2).*cos.(φ), (Dh/2).*sin.(φ); color=(:black,0.4))
    lines!(ax, (D/2).*cos.(φ), (D/2).*sin.(φ); color=(:black,0.15), linestyle=:dot)
    for b in 1:tab.Z
        xb,yb = outline(tab,b)
        poly!(ax, Point2f.(xb,yb); color=(col,0.35), strokecolor=col, strokewidth=1.5)
        mxb,myb = midline(tab,b)
        lines!(ax, mxb, myb; color=col, linewidth=2)
    end
    limits!(ax, -D/2*1.1, D/2*1.1, -D/2*1.1, D/2*1.1)
end
Label(fig2[0,:], "Assembled 5-blade propeller, viewed down the shaft (rotation sense identical)";
      fontsize=15, font=:bold)
save(joinpath(OUT, "skew_propellers.png"), fig2; px_per_unit=2)
println("wrote runs/skew_propellers.png")
