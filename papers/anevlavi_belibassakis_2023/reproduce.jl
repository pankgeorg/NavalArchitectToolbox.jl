#!/usr/bin/env julia
#
# Reproduction of Anevlavi, Zafeiris, Papadakis & Belibassakis,
# "A Low-Cost Vortex-Lattice Method ... Tip-Rake Reformation",
# J. Mar. Sci. Eng. 2023, 11, 2179 — open-water performance of DTMB 4382
# with the key-blade vortex-ring VLM in NavalArchitectToolbox.
#
# Activates this folder's own (Julia 1.12) environment, runs the VLM, prints
# the results table, and writes four figures to ./figures/. Headless CairoMakie.
#
#   julia +1.12 --project=papers/anevlavi_belibassakis_2023 \
#       papers/anevlavi_belibassakis_2023/reproduce.jl

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate(; io=devnull)

using NavalArchitectToolbox, CairoMakie, Printf, LinearAlgebra
CairoMakie.activate!(type="png")

const D, Dh = 6.0, 1.2          # diameter / hub (coefficients are scale-free)
const t   = dtmb4382
const OUT = joinpath(@__DIR__, "figures"); mkpath(OUT)

# Experimental open-water data for DTMB 4382, digitized (approximately) from
# the paper's Figure 7. The J=0.889 design point (★) is the firm anchor quoted
# in the text: KT=0.208, 10·KQ=0.445, η=0.661.
Jexp   = [0.6, 0.7, 0.8, 0.889, 1.0, 1.1, 1.2]
KTexp  = [0.330, 0.290, 0.245, 0.208, 0.165, 0.115, 0.070]
KQ10ex = [0.560, 0.510, 0.470, 0.445, 0.390, 0.310, 0.230]
ηexp   = [0.560, 0.605, 0.640, 0.661, 0.675, 0.650, 0.580]

# ─────────────────────────────────────────────────────────────────────────
# 1. open-water sweep
# ─────────────────────────────────────────────────────────────────────────
Js = sort(unique([collect(0.5:0.05:1.2); 0.889]))
KT = similar(Js); KQ10 = similar(Js); ηv = similar(Js)
for (i, J) in enumerate(Js)
    r = openwater_vlm(t, D, Dh, J; nc=12, ns=22)
    KT[i] = r.KT; KQ10[i] = 10r.KQ; ηv[i] = r.η
end

println("\nDTMB 4382 open-water — VLM (C_LE=0.80, C_Drag=0.010) vs experiment")
println("  J      KT_vlm  KT_exp   10KQ_vlm 10KQ_exp  η_vlm  η_exp")
for (J, kt, kq, e) in zip(Js, KT, KQ10, ηv)
    je = findfirst(≈(J; atol=1e-6), Jexp)
    if je !== nothing
        @printf("  %.3f  %.4f  %.4f   %.4f   %.4f   %.3f  %.3f\n",
                J, kt, KTexp[je], kq, KQ10ex[je], e, ηexp[je])
    end
end

# ─────────────────────────────────────────────────────────────────────────
# Figure 1 — open-water diagram (parallels paper Fig. 7)
# ─────────────────────────────────────────────────────────────────────────
fig1 = Figure(size=(820, 620))
ax = Axis(fig1[1,1]; xlabel="advance ratio  J = Vₐ/(nD)", ylabel="KT,  10·KQ,  η",
          title="DTMB 4382 open-water — key-blade vortex-lattice vs experiment")
lcKT, lcKQ, lcη = :tomato, :steelblue, :seagreen
lines!(ax, Js, KT;   color=lcKT, linewidth=2.5)
lines!(ax, Js, KQ10; color=lcKQ, linewidth=2.5)
lines!(ax, Js, ηv;   color=lcη,  linewidth=2.5)
scatter!(ax, Jexp, KTexp;  color=lcKT, marker=:utriangle, markersize=11, strokewidth=0.5)
scatter!(ax, Jexp, KQ10ex; color=lcKQ, marker=:rect,      markersize=10, strokewidth=0.5)
scatter!(ax, Jexp, ηexp;   color=lcη,  marker=:circle,    markersize=10, strokewidth=0.5)
# highlight the J=0.889 design anchor
vlines!(ax, [0.889]; color=(:gray, 0.5), linestyle=:dash)
text!(ax, 0.61, 0.10; text="KT",   color=lcKT, fontsize=18)
text!(ax, 0.61, 0.50; text="10·KQ", color=lcKQ, fontsize=18)
text!(ax, 0.61, 0.62; text="η",    color=lcη,  fontsize=18)
# legend proxies
elem_line = LineElement(color=:black, linewidth=2.5)
elem_pt   = MarkerElement(color=:black, marker=:circle, markersize=10)
axislegend(ax, [elem_line, elem_pt], ["VLM (this work)", "Exp. (Fig. 7, digitized)"];
           position=:rt, framevisible=true)
ylims!(ax, 0, 1)
save(joinpath(OUT, "openwater_4382.png"), fig1)
println("wrote figures/openwater_4382.png")

# ─────────────────────────────────────────────────────────────────────────
# 2. mesh convergence at J=0.889
# ─────────────────────────────────────────────────────────────────────────
meshes = [(6,12),(8,15),(10,18),(12,22),(14,26),(16,30)]
dof = Int[]; cKT = Float64[]; cKQ = Float64[]; cη = Float64[]
for (nc, ns) in meshes
    r = openwater_vlm(t, D, Dh, 0.889; nc=nc, ns=ns)
    push!(dof, nc*ns); push!(cKT, r.KT); push!(cKQ, 10r.KQ); push!(cη, r.η)
end

fig2 = Figure(size=(820, 560))
ax2 = Axis(fig2[1,1]; xlabel="lifting-surface panels (nc × ns)", ylabel="coefficient at J=0.889",
           title="Mesh convergence → experimental values (dashed)")
scatterlines!(ax2, dof, cKT;  color=lcKT, marker=:utriangle, label="KT")
scatterlines!(ax2, dof, cKQ;  color=lcKQ, marker=:rect,      label="10·KQ")
scatterlines!(ax2, dof, cη;   color=lcη,  marker=:circle,    label="η")
hlines!(ax2, [0.208]; color=lcKT, linestyle=:dash)
hlines!(ax2, [0.445]; color=lcKQ, linestyle=:dash)
hlines!(ax2, [0.661]; color=lcη,  linestyle=:dash)
axislegend(ax2; position=:rc)
save(joinpath(OUT, "convergence_4382.png"), fig2)
println("wrote figures/convergence_4382.png")

# ─────────────────────────────────────────────────────────────────────────
# 3. spanwise loading + chordwise Δc_p distribution (detailed solve)
# ─────────────────────────────────────────────────────────────────────────
r = openwater_vlm(t, D, Dh, 0.889; nc=12, ns=22)
nc, ns = r.nc, r.ns
lin(i,j) = (j-1)*nc + i
R = D/2
# spanwise bound circulation = cumulative chordwise ring strength at the TE
Γspan = [r.Γ[lin(nc,j)] for j in 1:ns]
rRspan = [hypot(r.cp[lin(nc,j)][1], r.cp[lin(nc,j)][2]) / R for j in 1:ns]
# non-dimensional circulation G* = Γ / (π D Vₐ)  (Vₐ at J=0.889, n=1)
Va = 0.889*1.0*D
Gstar = Γspan ./ (π*D*Va)

fig3 = Figure(size=(980, 460))
ax3a = Axis(fig3[1,1]; xlabel="r/R", ylabel="Γ / (π D Vₐ)",
            title="Spanwise bound circulation (J=0.889)")
scatterlines!(ax3a, rRspan, Gstar; color=:purple, markersize=7)
xlims!(ax3a, 0.2, 1.0)

ax3b = Axis(fig3[1,2]; xlabel="chordwise panel (LE→TE)", ylabel="r/R",
            title="Δc_p on the camber surface")
dcpM = reshape(r.dcp, nc, ns)
hm = heatmap!(ax3b, 1:nc, rRspan, dcpM; colormap=:balance,
              colorrange=(-maximum(abs, dcpM), maximum(abs, dcpM)))
Colorbar(fig3[1,3], hm; label="Δc_p")
save(joinpath(OUT, "loading_4382.png"), fig3)
println("wrote figures/loading_4382.png")

# ─────────────────────────────────────────────────────────────────────────
# 4. Δc_p on all blades, axial view (parallels paper Fig. 8)
# ─────────────────────────────────────────────────────────────────────────
rotz(p, ψ) = (c=cos(ψ); s=sin(ψ); (c*p[1]-s*p[2], s*p[1]+c*p[2], p[3]))
fig4 = Figure(size=(760, 720))
ax4 = Axis3(fig4[1,1]; xlabel="x [m]", ylabel="y [m]", zlabel="z [m]", aspect=:data,
            azimuth=-0.5π, elevation=0.5π,        # look down the shaft (+z)
            title="DTMB 4382 — Δc_p on the mean camber surface (J=0.889)")
cmax = maximum(abs, r.dcp)
hm4 = nothing
for b in 0:t.Z-1
    ψ = 2π*b/t.Z
    X = [rotz(r.cp[lin(i,j)], ψ)[1] for i in 1:nc, j in 1:ns]
    Y = [rotz(r.cp[lin(i,j)], ψ)[2] for i in 1:nc, j in 1:ns]
    Zc = [rotz(r.cp[lin(i,j)], ψ)[3] for i in 1:nc, j in 1:ns]
    global hm4 = surface!(ax4, X, Y, Zc; color=dcpM, colormap=:balance,
                          colorrange=(-cmax, cmax), shading=NoShading)
end
Colorbar(fig4[1,2], hm4; label="Δc_p")
save(joinpath(OUT, "dcp_blades_4382.png"), fig4)
println("wrote figures/dcp_blades_4382.png")

println("\ndone — 4 figures in $OUT")
