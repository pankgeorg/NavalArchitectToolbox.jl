module NavalArchitectToolbox

# Bridge from the naval-architecture propeller *section table* (chord,
# pitch, thickness, camber, skew, rake vs r/R) to 3-D blade geometry and a
# validated open-water vortex-lattice solver. We own the geometry: read the
# table, wrap each 2-D section onto its cylinder with the right
# pitch/skew/rake, and emit
#   - a vortex-lattice camber-surface grid (for VortexLattice.jl or the
#     built-in `openwater_vlm`),
#   - a full thick blade surface mesh (visualisation), and
#   - a signed-distance function (WaterLily `AutoBody` immersion;
#     ParametricBodies.jl can't represent a thick twisted blade).
#
# Section *shape* (mean line + thickness) is pluggable; the defaults are the
# propeller-standard NACA a=0.8 mean line + a closed-TE NACA-4-digit
# thickness stand-in for NACA 66-mod — see `MeanLine` / `Thickness` below.
#
# `openwater_vlm` (src/vlm.jl) reproduces the DTMB 4382 open-water curve
# (Anevlavi-Belibassakis, JMSE 2023); see the package README.

using StaticArrays

export PropellerBladeTable, read_blade_table, dtmb4382
export ParabolicMeanLine, NACAMeanLine, NACA66ish, blade_section_point
export blade_surface, vlm_camber_grid, pitch_angle, dimensional, blade_sdf
export openwater_vlm

# ----------------------------------------------------------------------------
# The section table
# ----------------------------------------------------------------------------

"""
    PropellerBladeTable{T}

Radial blade-section table for a marine propeller (one entry per `r/R`
station). Columns mirror the standard naval-architecture format:

| field  | meaning                              |
|--------|--------------------------------------|
| `rR`   | r/R — non-dimensional radius         |
| `cD`   | c/D — chord / diameter               |
| `PD`   | P/D — pitch / diameter               |
| `tc`   | tmax/c — max thickness / chord       |
| `fc`   | fmax/c — max camber / chord          |
| `skew` | θs — skew angle [degrees]            |
| `rake` | XR/R — rake / R                      |

`Z` is the blade count. Stations are assumed sorted in increasing `rR`.
"""
struct PropellerBladeTable{T}
    rR   :: Vector{T}
    cD   :: Vector{T}
    PD   :: Vector{T}
    tc   :: Vector{T}
    fc   :: Vector{T}
    skew :: Vector{T}
    rake :: Vector{T}
    Z    :: Int
end

"""
    read_blade_table(text; Z=5, T=Float64)

Parse a whitespace-delimited blade table with columns
`r/R c/D P/D tmax/c fmax/c θs XR/R` (one row per station; header lines
that don't start with a number are skipped). Returns a
[`PropellerBladeTable`](@ref).
"""
function read_blade_table(text::AbstractString; Z::Int=5, T::Type=Float64)
    cols = [T[] for _ in 1:7]
    for ln in eachline(IOBuffer(text))
        toks = split(strip(ln))
        length(toks) < 7 && continue
        v = tryparse.(T, toks[1:7])
        any(isnothing, v) && continue          # header / non-numeric row
        for j in 1:7; push!(cols[j], v[j]); end
    end
    p = sortperm(cols[1])
    PropellerBladeTable{T}((c[p] for c in cols)..., Z)
end

"DTMB 4382 (5-blade, 36°-skew DTMB modified series). Data from the standard table."
const dtmb4382 = read_blade_table("""
r/R     c/D     P/D     tmax/c  fmax/c  θs       XR/R
0.2000  0.1740  1.4550  0.2494  0.0430  0        0
0.2500  0.2020  1.4440  0.1960  0.0395  2.3280   0
0.3000  0.2290  1.4330  0.1562  0.0370  4.6550   0
0.4000  0.2750  1.4120  0.1068  0.0344  9.3630   0
0.5000  0.3120  1.3610  0.0768  0.0305  13.9480  0
0.6000  0.3370  1.2850  0.0566  0.0247  18.3780  0
0.7000  0.3470  1.2000  0.0421  0.0199  22.7470  0
0.8000  0.3340  1.1120  0.0314  0.0161  27.1450  0
0.9000  0.2800  1.0270  0.0239  0.0134  31.5750  0
0.9500  0.2100  0.9850  0.0229  0.0140  33.7880  0
1.0000  0.0100  0.9420  0.0160  0.0134  36.000   0
"""; Z=5)

# ----------------------------------------------------------------------------
# Linear interpolation of a column at an arbitrary r/R
# ----------------------------------------------------------------------------

@inline function _interp(xs::Vector{T}, ys::Vector{T}, x) where T
    x ≤ xs[1]   && return ys[1]
    x ≥ xs[end] && return ys[end]
    @inbounds for i in 1:length(xs)-1
        if xs[i] ≤ x ≤ xs[i+1]
            t = (x - xs[i]) / (xs[i+1] - xs[i])
            return ys[i] + t * (ys[i+1] - ys[i])
        end
    end
    return ys[end]
end

# ----------------------------------------------------------------------------
# Section shape: pluggable mean line + thickness (normalised, max = 1)
# ----------------------------------------------------------------------------

abstract type MeanLine end
abstract type Thickness end

"Parabolic mean line yc/fmax = 4x(1−x), x = chordwise fraction. Quick stand-in."
struct ParabolicMeanLine <: MeanLine end
@inline (::ParabolicMeanLine)(x) = 4x*(1 - x)

# NACA a-series mean-line ordinate (design cli=1), Abbott & von Doenhoff
# "Theory of Wing Sections" eq. for the uniform-load-to-a mean line.
function _camber_a(a, x)
    (x ≤ 0 || x ≥ 1) && return zero(x)
    xlnx    = x ≤ eps(x)        ? zero(x) : x*log(x)
    ax = a - x
    axterm  = abs(ax) ≤ eps(x)  ? zero(x) : (ax^2/2)*log(abs(ax))
    omx = 1 - x
    omxterm = omx ≤ eps(x)      ? zero(x) : (omx^2/2)*log(omx)
    g = -(1/(1-a))*(a^2*(log(a)/2 - 1/4) + 1/4)
    h =  (1/(1-a))*((1-a)^2*log(1-a)/2 - (1-a)^2/4) + g
    term = (1/(1-a))*(axterm - omxterm + omx^2/4 - ax^2/4) - xlnx + g - h*x
    return term / (2π*(a+1))
end

"""
    NACAMeanLine(a=0.8)

Exact NACA a-series mean line (the propeller-standard **a=0.8** uniform-
loading-to-80%-chord camber line). As a [`MeanLine`](@ref) it returns the
*shape normalised to unit max*, so it scales to the table's `fmax/c`.
"""
struct NACAMeanLine <: MeanLine
    a::Float64
    norm::Float64
end
function NACAMeanLine(a::Real=0.8)
    mx = maximum(_camber_a(a, x) for x in range(1e-4, 1-1e-4, length=2001))
    NACAMeanLine(float(a), mx)
end
@inline (m::NACAMeanLine)(x) = _camber_a(m.a, x) / m.norm

# Propeller-standard default section camber; built once (the constructor runs
# a max over the chord), shared as the default mean line across the package.
const _DEFAULT_MEANLINE = NACAMeanLine(0.8)

"Normalised NACA-4-digit half-thickness with a *closed* (sharp) trailing
edge — max = 1 at x≈0.3, zero at both LE and TE. Stand-in for NACA 66-mod."
struct NACA66ish <: Thickness end
@inline _naca4(x) = 0.2969*sqrt(x) - 0.1260x - 0.3516x^2 + 0.2843x^3 - 0.1015x^4
# Close the blunt 4-digit TE: subtract the linear term so shape(1)=0
# (propeller sections have a sharp trailing edge).
@inline _naca4_closed(x) = _naca4(x) - x*_naca4(1.0)
const _NACA4_MAX = _naca4_closed(0.3)          # ≈ 0.099 (half-thickness peak)
@inline (::NACA66ish)(x) = _naca4_closed(x) / _NACA4_MAX

# ----------------------------------------------------------------------------
# Geometry: wrap a section onto its cylinder
# ----------------------------------------------------------------------------

"""
    pitch_angle(tab, rR) -> β  (radians)

Section pitch angle `β = atan(P/(2πr)) = atan((P/D)/(π·r/R))`.
"""
pitch_angle(tab::PropellerBladeTable, rR) = atan(_interp(tab.rR, tab.PD, rR) / (π * rR))

"""
    dimensional(tab, rR, D) -> (r, c, β, tmax, fmax, θskew, rake)

Dimensional section quantities at `r/R = rR` for diameter `D`:
radius, chord, pitch angle (rad), max thickness, max camber, skew (rad),
axial rake.
"""
function dimensional(tab::PropellerBladeTable, rR, D)
    R = D/2
    r = rR * R
    c = _interp(tab.rR, tab.cD, rR) * D
    β = atan(_interp(tab.rR, tab.PD, rR) / (π * rR))
    tmax = _interp(tab.rR, tab.tc, rR) * c
    fmax = _interp(tab.rR, tab.fc, rR) * c
    θs   = deg2rad(_interp(tab.rR, tab.skew, rR))
    rake = _interp(tab.rR, tab.rake, rR) * R
    return (; r, c, β, tmax, fmax, θs, rake)
end

"""
    blade_section_point(tab, rR, xc, D; surface=:camber, blade=1,
                        meanline=ParabolicMeanLine(), thickness=NACA66ish()) -> SVector{3}

3-D point on blade `blade` (1-based, of `Z`) at radius fraction `rR` and
chordwise fraction `xc ∈ [0,1]` (0 = leading edge, 1 = trailing edge).
`surface` ∈ {`:camber`, `:upper`, `:lower`}.

The 2-D section point `(s, n)` — chordwise distance `s` from mid-chord,
normal offset `n` (camber, ± half-thickness) — is rotated by the pitch
angle into the (tangential, axial) plane, wrapped onto the cylinder of
radius `r`, then offset by skew (circumferential) and rake (axial):

```
Δtang = s·cosβ + n·sinβ      Δax = s·sinβ − n·cosβ
θ = 2π(blade−1)/Z + θskew + Δtang/r
(x, y, z) = (r·cosθ, r·sinθ, rake + Δax)
```
"""
function blade_section_point(tab::PropellerBladeTable{T}, rR, xc, D;
                             surface::Symbol=:camber, blade::Int=1,
                             meanline::MeanLine=_DEFAULT_MEANLINE,
                             thickness::Thickness=NACA66ish()) where T
    g = dimensional(tab, rR, D)
    s = (xc - one(T)/2) * g.c                       # chordwise dist from mid-chord
    yc = g.fmax * meanline(xc)                      # camber
    yt = g.tmax * thickness(xc) / 2                 # half-thickness
    n = surface === :upper ? yc + yt :
        surface === :lower ? yc - yt : yc
    sinβ, cosβ = sincos(g.β)
    # The section normal-offset (camber + thickness) sits on the suction side
    # for the rotation/inflow convention in `_inflow` (relative flow runs LE→TE
    # with the camber bulging toward suction). Δtang wraps the chord onto the
    # cylinder at the pitch angle; Δax places it along the shaft.
    Δtang = s*cosβ + n*sinβ
    Δax   = s*sinβ - n*cosβ
    θ = 2π*(blade-1)/tab.Z + g.θs + Δtang / g.r
    z = g.rake + Δax
    return SVector{3,T}(g.r*cos(θ), g.r*sin(θ), z)
end

# Cosine spacing on [0,1] clustered at both ends (chordwise: LE & TE).
_cosine(n) = [(1 - cos(π*i/n))/2 for i in 0:n]
# Cosine spacing on [a,b] clustered at the high end (spanwise: the tip).
_tip_cluster(a, b, n) = [a + (b-a)*(1 - cos(π*i/(2n))) for i in 0:n]

"""
    vlm_camber_grid(tab, D, Dh; nc=8, ns=16, blade=1, meanline=…)
        -> Matrix{SVector{3}}  of size (nc+1, ns+1)

Vortex-lattice **camber-surface** grid (no thickness): the lifting
surface for `VortexLattice.jl`. Rows are chordwise (LE→TE, cosine-
spaced), columns spanwise (hub→tip, cosine-clustered at the tip where
c/D collapses). `Dh` is the hub diameter; the root station is `Dh/D`.
"""
function vlm_camber_grid(tab::PropellerBladeTable{T}, D, Dh;
                         nc::Int=8, ns::Int=16, blade::Int=1,
                         meanline::MeanLine=_DEFAULT_MEANLINE,
                         tip_rR::Real=0.98) where T
    rR_hub = Dh / D
    # Terminate the lifting surface just inboard of the geometric tip:
    # at r/R=1 the chord collapses to ~0, giving degenerate VLM panels.
    rRs = _tip_cluster(rR_hub, T(tip_rR), ns)
    xcs = _cosine(nc)
    grid = Matrix{SVector{3,T}}(undef, nc+1, ns+1)
    for (js, rR) in pairs(rRs), (ic, xc) in pairs(xcs)
        grid[ic, js] = blade_section_point(tab, rR, xc, D;
            surface=:camber, blade=blade, meanline=meanline)
    end
    return grid
end

"""
    blade_surface(tab, D, Dh; nc=24, ns=24, blade=1, …)
        -> (upper, lower)  each a (nc+1, ns+1) Matrix{SVector{3}}

Full (thick) blade surface for visualisation or a future immersed body:
upper (suction) and lower (pressure) skins, chordwise × spanwise.
"""
function blade_surface(tab::PropellerBladeTable{T}, D, Dh;
                       nc::Int=24, ns::Int=24, blade::Int=1,
                       meanline::MeanLine=_DEFAULT_MEANLINE,
                       thickness::Thickness=NACA66ish()) where T
    rR_hub = Dh / D
    rRs = _tip_cluster(rR_hub, tab.rR[end], ns)
    xcs = _cosine(nc)
    up = Matrix{SVector{3,T}}(undef, nc+1, ns+1)
    lo = Matrix{SVector{3,T}}(undef, nc+1, ns+1)
    for (js, rR) in pairs(rRs), (ic, xc) in pairs(xcs)
        up[ic, js] = blade_section_point(tab, rR, xc, D; surface=:upper,
            blade=blade, meanline=meanline, thickness=thickness)
        lo[ic, js] = blade_section_point(tab, rR, xc, D; surface=:lower,
            blade=blade, meanline=meanline, thickness=thickness)
    end
    return up, lo
end

# ----------------------------------------------------------------------------
# Signed-distance function — for a WaterLily AutoBody / ShipShapes.TabulatedHull
# ----------------------------------------------------------------------------
#
# ParametricBodies.jl (the WaterLily ecosystem option) only models
# parametric *curves* and thin *planar membranes* (PlanarBody has a fixed
# ~grid-scale thickness), so it cannot represent a thick, twisted,
# cambered 3-D blade. The practical immersion path in this stack is an
# SDF wrapped in WaterLily's `AutoBody(sdf, map)`.
#
# Per blade we invert the section map exactly at the query's radius: at
# fixed r the (s,n)→(Δtang,Δax) map is a rotation by β, so given the
# query's (Δθ·r, Δz) we recover (s,n), hence the chordwise fraction and
# the signed normal distance to the section's thickness envelope. The
# union over Z blades is the min SDF. Approximate near the LE/TE/hub/tip
# edges (the bulk skin distance is what BDIM needs).

@inline function _sdf_one_blade(tab::PropellerBladeTable{T}, D, Dh, x, blade,
                                meanline, thickness) where T
    R = D/2; r_hub = Dh/2
    rq = hypot(x[1], x[2]); θq = atan(x[2], x[1]); zq = x[3]
    rRc = clamp(rq/R, tab.rR[1], one(T))
    g = dimensional(tab, rRc, D)
    θ0 = 2π*(blade-1)/tab.Z + g.θs
    Δθ = atan(sin(θq-θ0), cos(θq-θ0))               # wrapped to (-π,π]
    Δtang = Δθ * g.r;  Δax = zq - g.rake
    sinβ, cosβ = sincos(g.β)
    s =  Δtang*cosβ + Δax*sinβ                       # inverse of the section map
    n =  Δtang*sinβ - Δax*cosβ
    xc = s/g.c + one(T)/2
    # signed normal distance to the section thickness envelope
    if zero(T) ≤ xc ≤ one(T)
        yc = g.fmax*meanline(xc); yt = g.tmax*thickness(xc)/2
        d_sec = abs(n - yc) - yt                     # <0 inside the skin
    else                                             # past LE/TE corner
        xcl = clamp(xc, zero(T), one(T))
        yc = g.fmax*meanline(xcl); yt = g.tmax*thickness(xcl)/2
        ds = (xc < 0 ? -xc : xc - 1) * g.c           # chordwise overshoot
        d_sec = hypot(ds, max(abs(n-yc) - yt, zero(T)))
    end
    # radial overshoot beyond hub/tip
    dr = rq < r_hub ? r_hub - rq : rq > R ? rq - R : zero(T)
    return dr > 0 ? hypot(max(d_sec, zero(T)), dr) : d_sec
end

"""
    blade_sdf(tab, D, Dh; meanline=NACAMeanLine(0.8), thickness=NACA66ish())

Return a signed-distance function `sdf(x)` (negative inside, ≈0 on the
skin, positive outside) for the full `tab.Z`-blade propeller of diameter
`D`, hub `Dh`. Wrap it for WaterLily with `AutoBody((x,t)->sdf(x))`, or
sample it onto a grid for `ShipShapes.TabulatedHull`.

The SDF inverts the section map exactly at each query radius; it is
approximate near the leading/trailing edges and the hub/tip caps.
"""
function blade_sdf(tab::PropellerBladeTable{T}, D, Dh;
                   meanline::MeanLine=_DEFAULT_MEANLINE,
                   thickness::Thickness=NACA66ish()) where T
    return function (x)
        d = T(Inf)
        for b in 1:tab.Z
            d = min(d, _sdf_one_blade(tab, D, Dh, x, b, meanline, thickness))
        end
        return d
    end
end

include("vlm.jl")

end # module
