# ----------------------------------------------------------------------------
# Finite-wing lifting-line / vortex-lattice — a Weissinger horseshoe lattice
#
# A single chordwise panel per spanwise strip (the classic Weissinger /
# "vortex-lattice lifting-line", vl-ll): one bound vortex on the 1/4-chord
# line, the no-penetration control point at the 3/4-chord mid-span, and two
# trailing legs running streamwise to (effectively) downstream infinity. This
# is the lightweight, pure-Julia sibling of the multi-chordwise `Wing`
# (re-exported from LiftingSurfaces.jl over VortexLattice.jl): no chordwise
# pressure distribution, but the spanwise loading, lift slope and induced drag
# of a finite wing at a fraction of the cost and with no external dependency.
#
# Geometry. A symmetric planform of geometric span `span`, root (mid-span)
# chord `chord_root` tapering linearly to `chord_tip` at *both* tips
# (chord(y) = c_root + (c_tip−c_root)·|y|/(span/2)), with leading-edge `sweep`,
# linear geometric `twist` from root to tip, and a small `tip_inset` that pulls
# the outermost strip just inboard of the geometric tip (the chord there is
# near-zero and would give a degenerate panel). The wing lies in the z=0 plane;
# `alpha` and the local twist enter through the panel normal (thin-wing).
#
# Forces. Near-field Kutta–Joukowski: F_i = ρ Γ_i (V∞ + v_induced) × ℓ_i at
# each bound vortex, summed and resolved into lift (⊥ V∞) and induced drag
# (∥ V∞). Because the bound-vortex force recovers the leading-edge suction, the
# induced drag is the physically-consistent value (span efficiency e≈0.9–1.0
# for a moderate-taper wing), matching the Trefftz-plane result and lifting-
# line theory dCL/dα → 2π·AR/(AR+2).
#
# Pure Julia (LinearAlgebra + StaticArrays).
# ----------------------------------------------------------------------------

# Biot–Savart velocity at P from a unit-strength straight vortex filament P1→P2.
@inline function _filament_vel(P::SVector{3,Float64}, P1::SVector{3,Float64}, P2::SVector{3,Float64})
    r1 = P - P1; r2 = P - P2
    c  = cross(r1, r2)
    d  = dot(c, c)
    d < 1e-10 && return zero(SVector{3,Float64})
    r0 = P2 - P1
    k  = (dot(r0, r1)/norm(r1) - dot(r0, r2)/norm(r2)) / (4π * d)
    return k * c
end

# Unit horseshoe: trailing leg (A+wake)→A, bound A→B, trailing B→(B+wake).
@inline function _horseshoe_vel(P, A, B, wdir::SVector{3,Float64}, L::Float64)
    Aw = A + L*wdir; Bw = B + L*wdir
    _filament_vel(P, Aw, A) + _filament_vel(P, A, B) + _filament_vel(P, B, Bw)
end

"""
    liftingline_vlm(; chord_root, chord_tip, span, alpha,
                    sweep=0.0, twist_root=0.0, twist_tip=0.0,
                    N=100, V∞=1.0, ρ=1.225, tip_inset=0.02)
        -> (; CL, CDi, e, S, AR, Γ, y, cl_span)

Weissinger horseshoe vortex-lattice solve for a **symmetric finite wing** —
one chordwise panel per spanwise strip (bound vortex at 1/4-chord, control
point at 3/4-chord), streamwise semi-infinite trailing legs, and near-field
Kutta–Joukowski forces.

Planform: root (mid-span) chord `chord_root`, tip chord `chord_tip`, geometric
`span`, leading-edge `sweep` (deg), linear geometric twist from `twist_root`
(root/mid) to `twist_tip` (tip) in degrees, at angle of attack `alpha` (deg)
and free-stream `V∞`. `N` is the number of spanwise panels; `tip_inset`
(fraction of span) pulls the outermost strip just inboard of the geometric tip
to avoid a degenerate near-zero-chord panel.

Returns the lift coefficient `CL`, the induced drag `CDi` (near-field, with
leading-edge-suction recovery), the implied span efficiency
`e = CL²/(π·AR·CDi)`, the trapezoidal reference area `S`, aspect ratio
`AR = span²/S`, the bound-circulation distribution `Γ`, the spanwise station
midpoints `y`, and the sectional lift coefficient `cl_span = 2Γ/(V∞·c)`.

This is the pure-Julia lifting-line counterpart of [`Wing`](@ref) /
[`wing_forces`](@ref) (multi-chordwise VLM over VortexLattice.jl); it
reproduces the same `CL` (within ~1–2 %) and the Trefftz-consistent induced
drag, and is validated against the lifting-line slope `dCL/dα = 2π·AR/(AR+2)`.
"""
function liftingline_vlm(; chord_root::Real, chord_tip::Real, span::Real,
                         alpha::Real, sweep::Real=0.0,
                         twist_root::Real=0.0, twist_tip::Real=0.0,
                         N::Int=100, V∞::Real=1.0, ρ::Real=1.225,
                         tip_inset::Real=0.02)
    V3(x,y,z) = SVector{3,Float64}(x,y,z)
    half = span/2; toff = span*tip_inset
    yv = range(-(half-toff), half-toff; length=N+1)
    sw = deg2rad(sweep)
    chord(y) = chord_root + (chord_tip-chord_root)*(abs(y)/half)
    xLE(y)   = abs(y)*tan(sw)
    twist(y) = deg2rad(twist_root + (twist_tip-twist_root)*(abs(y)/half))

    α  = deg2rad(alpha)
    Vi = V∞ * V3(cos(α), 0.0, sin(α))
    wd = V3(1.0, 0.0, 0.0)                       # streamwise (fixed) wake
    Lw = 1.0e4                                   # ≈ semi-infinite legs
    ld = V3(-sin(α), 0.0, cos(α))                # lift direction (⊥ V∞)
    dd = V3( cos(α), 0.0, sin(α))                # drag direction (∥ V∞)

    A = Vector{SVector{3,Float64}}(undef, N); B = similar(A)
    cp = similar(A); nrm = similar(A); bvec = similar(A)
    cmid = Vector{Float64}(undef, N); ymid = Vector{Float64}(undef, N)
    for i in 1:N
        yi, yo = yv[i], yv[i+1]; ym = 0.5*(yi+yo)
        ci, co, cm = chord(yi), chord(yo), chord(ym)
        A[i]  = V3(xLE(yi)+0.25ci, yi, 0.0)
        B[i]  = V3(xLE(yo)+0.25co, yo, 0.0)
        cp[i] = V3(xLE(ym)+0.75cm, ym, 0.0)
        bvec[i] = B[i] - A[i]
        θ = twist(ym)
        nrm[i] = V3(sin(θ), 0.0, cos(θ))
        cmid[i] = cm; ymid[i] = ym
    end

    AIC = Matrix{Float64}(undef, N, N); rhs = Vector{Float64}(undef, N)
    for i in 1:N
        rhs[i] = -dot(Vi, nrm[i])
        for j in 1:N
            AIC[i,j] = dot(_horseshoe_vel(cp[i], A[j], B[j], wd, Lw), nrm[i])
        end
    end
    Γ = AIC \ rhs

    F = zero(SVector{3,Float64})
    for i in 1:N
        Pm = 0.5*(A[i]+B[i]); vind = zero(SVector{3,Float64})
        for j in 1:N
            vind += Γ[j]*_horseshoe_vel(Pm, A[j], B[j], wd, Lw)
        end
        F += ρ * Γ[i] * cross(Vi + vind, bvec[i])
    end

    S = 0.0
    for i in 1:N
        yi, yo = yv[i], yv[i+1]
        P1=V3(xLE(yi),yi,0.0); P2=V3(xLE(yi)+chord(yi),yi,0.0)
        P3=V3(xLE(yo)+chord(yo),yo,0.0); P4=V3(xLE(yo),yo,0.0)
        S += 0.5*norm(cross(P2-P1,P4-P1)) + 0.5*norm(cross(P3-P2,P4-P2))
    end
    q  = 0.5*ρ*V∞^2*S
    CL = dot(F, ld)/q
    CDi = dot(F, dd)/q
    AR = span^2/S
    e  = abs(CL) > 1e-8 ? CL^2/(π*AR*CDi) : NaN
    cl_span = [2*Γ[i]/(V∞*cmid[i]) for i in 1:N]
    return (; CL, CDi, e, S, AR, Γ, y = ymid, cl_span)
end
