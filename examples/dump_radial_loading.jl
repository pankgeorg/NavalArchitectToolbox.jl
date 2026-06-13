#!/usr/bin/env julia
#
# Ladder 2 support: dump the DTMB 4382 open-water VLM *radial loading*
# distribution — dT/dr and dQ/dr vs r/R — at chosen advance ratios J,
# for use as the body-force source of the in-grid actuator disk in
# WaterLily (ShipFlow.jl/scripts/vlm_disk_roundtrip.jl).
#
# Method. `openwater_vlm` returns the per-panel pressure-normal force
# vector `force[k]` and collocation point `cp[k]` for the KEY blade
# (shaft = z). Per panel:
#   dThrust = Z · force[k][3]                       (axial, +z)
#   dTorque = Z · (cp[k][1]·force[k][2] − cp[k][2]·force[k][1])  (about z)
#   r       = hypot(cp[k][1], cp[k][2])
# Bin chordwise (sum over the nc panels at each spanwise strip) → a
# sectional dT/dr, dQ/dr at the strip radius. The *shape* is the VLM's;
# we then renormalise so Σ dThrust = Tt and Σ dTorque = Qt with Tt, Qt
# the VLM's validated totals (KT,KQ → which match experiment ~2%). This
# folds the friction/leading-edge-suction calibration that the integral
# coefficients carry into the radial profile consistently.
#
# Writes one CSV per J: rR, dTdr_norm, dQdr_norm (normalised so the
# trapezoidal integral over r equals 1), plus a header line carrying
# KT, KQ, eta, D, Dh, Z, nRPS so the consumer has everything.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NavalArchitectToolbox
using Printf

const D, Dh = 6.0, 1.2
const R = D/2
const Z = dtmb4382.Z
const nRPS = 1.0
const Js = [0.6, 0.889, 1.1]

OUT = joinpath(@__DIR__, "..", "runs"); mkpath(OUT)

# experimental anchors (Anevlavi-Belibassakis 2023 Fig 7, design pt firm)
Jexp  = [0.6, 0.7, 0.8, 0.889, 1.0, 1.1, 1.2]
KTexp = [0.330, 0.290, 0.245, 0.208, 0.165, 0.115, 0.070]
KQexp = [0.640, 0.600, 0.540, 0.445, 0.380, 0.290, 0.190] ./ 10  # 10·KQ digitized
nearest(J) = argmin(abs.(Jexp .- J))

@printf "DTMB 4382 radial-loading dump  (D=%.1f Dh=%.1f Z=%d n=%.1f)\n" D Dh Z nRPS
@printf "%6s | %7s %7s %7s | %7s %7s | %s\n" "J" "KT_vlm" "10KQ_v" "eta" "KT_exp" "10KQ_e" "csv"

for J in Js
    res = openwater_vlm(dtmb4382, D, Dh, J)
    nc, ns = res.nc, res.ns
    lin(i,j) = (j-1)*nc + i
    # strip radii (spanwise) and chordwise-summed sectional loads
    rR   = zeros(ns)
    dTdr = zeros(ns)   # per-strip thrust  (not yet /dr)
    dQdr = zeros(ns)   # per-strip torque
    for j in 1:ns
        rsum = 0.0
        for i in 1:nc
            k = lin(i,j)
            f = res.force[k]; pc = res.cp[k]
            dTdr[j] += Z * f[3]
            dQdr[j] += Z * (pc[1]*f[2] - pc[2]*f[1])
            rsum += hypot(pc[1], pc[2])
        end
        rR[j] = (rsum/nc)/R
    end
    # ensure positive sense (thrust +)
    if sum(dTdr) < 0; dTdr .= -dTdr; end
    if sum(dQdr) < 0; dQdr .= -dQdr; end
    # renormalise the integral to the VLM total coefficients
    Tt = res.KT * 1.0 * nRPS^2 * D^4         # ρ=1
    Qt = res.KQ * 1.0 * nRPS^2 * D^5
    sT = sum(dTdr); sQ = sum(dQdr)
    dT = dTdr .* (Tt/sT)
    dQ = dQdr .* (Qt/sQ)
    # convert per-strip → per-radius density dT/dr via strip widths
    # (trapezoidal strip edges from the cosine-clustered radii)
    redges = similar(rR, ns+1)
    redges[1] = Dh/D
    for j in 2:ns; redges[j] = 0.5*(rR[j-1]+rR[j]); end
    redges[ns+1] = 1.0
    dr = [redges[j+1]-redges[j] for j in 1:ns]  # in r/R units
    dTdr_dens = dT ./ (dr .* R)                  # per dimensional r
    dQdr_dens = dQ ./ (dr .* R)
    je = nearest(J)
    @printf "%6.3f | %7.4f %7.4f %7.4f | %7.4f %7.4f | radial_loading_J%.3f.csv\n" J res.KT 10*res.KQ res.η KTexp[je] 10*KQexp[je] J
    fn = joinpath(OUT, @sprintf("radial_loading_J%.3f.csv", J))
    open(fn, "w") do io
        @printf io "# DTMB4382 radial loading  J=%.4f KT=%.5f KQ=%.6f eta=%.5f D=%.3f Dh=%.3f Z=%d nRPS=%.3f Tt=%.5f Qt=%.5f\n" J res.KT res.KQ res.η D Dh Z nRPS Tt Qt
        println(io, "rR,dThrust_strip,dTorque_strip,dTdr,dQdr")
        for j in 1:ns
            @printf io "%.5f,%.6e,%.6e,%.6e,%.6e\n" rR[j] dT[j] dQ[j] dTdr_dens[j] dQdr_dens[j]
        end
    end
end
println("done")
