# ----------------------------------------------------------------------------
# Flettner rotor — 2D potential flow around a rotating cylinder
# (Hess–Smith constant-strength source-panel method + a prescribed bound
# circulation Γ). Pure Julia (LinearAlgebra only); the lightweight sibling
# of the propeller VLM kernel in vlm.jl.
#
# Physics. A spinning cylinder of radius R in a uniform stream V∞ develops a
# bound circulation Γ (Magnus effect). For an *inviscid* rotor the rotation
# only fixes Γ — there is no Kutta point on a smooth cylinder — so we impose
# Γ directly rather than solving a trailing-edge condition. With the surface
# moving at speed ωR, the no-slip-consistent circulation is
#
#     Γ = 2π ω R²            (= surface speed ωR times perimeter 2πR / ... )
#
# i.e. solid-body rotation of the bound vortex. The closed-form surface
# pressure and lift follow from superposing the uniform stream, a doublet
# (the non-lifting cylinder) and a point vortex Γ at the centre:
#
#     V_t(θ)   = 2 V∞ sinθ + Γ/(2πR)
#     Cp(θ)    = 1 − (V_t/V∞)²                       (the analytic reference)
#     C_L      = ½ ∮ −Cp sinθ dθ  =  Γ/(R V∞)  =  4π ω R²/V∞   (with V∞ chord c=2R)
#
# The panel method reproduces these from a discrete no-penetration solve, and
# is the thing we grid-refine to an ε(C_L) < 0.02 % tolerance.
# At ω=0 it must collapse to the classic non-lifting cylinder Cp = 1 − 4sin²θ.
#
# Sign/orientation convention: panels run COUNTER-CLOCKWISE (θ increasing),
# V∞ along +x. A positive ω gives Γ = 2πωR² > 0 and **positive lift** in the
# convention C_L = ½∮ −Cp sinθ dθ (matching the analytic C_L = +4πωR²/V∞).
# Because the CCW panel ordering makes the source-resolved freestream
# surface-tangential run opposite to the analytic `Vt = 2V∞ sinθ + Γ/(2πR)`
# convention, the central point vortex is taken with the matching sign so its
# velocity reinforces the flow over the *top* of the cylinder (the lift-
# producing side); see the in-code note. Verified against `flettner_analytic`.
# ----------------------------------------------------------------------------

"""
    flettner_analytic(; R=0.5, ω=1.0, V∞=1.0, n=720) -> (; θ, Cp, CL)

Closed-form inviscid rotating-cylinder solution (the reference the panel
method is validated against). Returns the surface angle `θ` (radians,
CCW from +x, `n` samples), the analytic pressure coefficient
`Cp(θ) = 1 − ((2V∞ sinθ + Γ/(2πR))/V∞)²` with `Γ = 2π ω R²`, and the
lift coefficient `CL = ½∮ −Cp sinθ dθ` (chord `c = 2R`), which evaluates
to the exact `CL = 4π ω R² / V∞`.
"""
function flettner_analytic(; R::Real=0.5, ω::Real=1.0, V∞::Real=1.0, n::Int=720)
    Γ = 2π * ω * R^2
    θ = range(0, 2π; length=n+1)[1:n]
    vt = @. 2 * V∞ * sin(θ) + Γ / (2π * R)
    Cp = @. 1 - (vt / V∞)^2
    # CL = ½ ∮ −Cp sinθ dθ (trapezoid over the closed loop); analytic = 4πωR²/V∞
    dθ = 2π / n
    CL = 0.5 * sum(@. -Cp * sin(θ)) * dθ
    return (; θ = collect(θ), Cp = collect(Cp), CL)
end

# Influence of a constant-strength 2D source panel (Hess–Smith). Returns the
# (u, v) velocity at point `p` induced by a UNIT-strength source distributed
# over the straight panel from `a` to `b`. Local-frame closed form (Katz &
# Plotkin §11.2.1 / Kuethe & Chow), transformed back to global axes.
@inline function _source_panel_vel(px, py, ax, ay, bx, by)
    dx = bx - ax; dy = by - ay
    L = sqrt(dx*dx + dy*dy)
    sinp = dy / L; cosp = dx / L
    # point in panel-local frame (x along panel a→b, y normal)
    xt = px - ax; yt = py - ay
    xl =  xt*cosp + yt*sinp
    yl = -xt*sinp + yt*cosp
    r1 = sqrt(xl*xl + yl*yl)
    r2 = sqrt((xl - L)^2 + yl*yl)
    θ1 = atan(yl, xl)
    θ2 = atan(yl, xl - L)
    # local-frame velocities of a unit constant-strength source panel
    ul = log(r1 / r2) / (2π)
    vl = (θ2 - θ1) / (2π)
    # back to global
    u = ul*cosp - vl*sinp
    v = ul*sinp + vl*cosp
    return u, v
end

# Influence of a constant-strength 2D *vortex* panel of unit strength. By the
# source/vortex duality the vortex-panel velocity is the source-panel velocity
# rotated −90°: (u_vortex, v_vortex) = (v_source, −u_source) in the panel frame.
@inline function _vortex_panel_vel(px, py, ax, ay, bx, by)
    dx = bx - ax; dy = by - ay
    L = sqrt(dx*dx + dy*dy)
    sinp = dy / L; cosp = dx / L
    xt = px - ax; yt = py - ay
    xl =  xt*cosp + yt*sinp
    yl = -xt*sinp + yt*cosp
    r1 = sqrt(xl*xl + yl*yl)
    r2 = sqrt((xl - L)^2 + yl*yl)
    θ1 = atan(yl, xl)
    θ2 = atan(yl, xl - L)
    # vortex panel = source rotated −90° in the local frame
    ul =  (θ2 - θ1) / (2π)
    vl = -log(r1 / r2) / (2π)
    u = ul*cosp - vl*sinp
    v = ul*sinp + vl*cosp
    return u, v
end

"""
    flettner_panel(; R=0.5, ω=1.0, V∞=1.0, N=200) -> (; θ, Cp, CL, Γ)

2D Hess–Smith potential-flow solve for a rotating cylinder of radius `R`
in a uniform stream `V∞` (along +x), with a prescribed bound circulation
`Γ = 2π ω R²` (the Magnus circulation of a smooth spinning cylinder — no
Kutta condition, the rotation sets Γ directly).

The cylinder is discretized into `N` constant-strength **source** panels
arranged counter-clockwise; per-panel source strengths `σ` enforce
no-penetration at the panel control points, with the freestream **and**
the prescribed point-vortex `Γ` (placed at the centre) on the right-hand
side. Surface pressure is recovered from the tangential velocity at each
control point, `Cp = 1 − (V_t/V∞)²`, and the lift from
`CL = ½∮ −Cp sinθ dθ` over chord `c = 2R`.

Returns the control-point angles `θ` (radians, CCW), `Cp(θ)`, the lift
coefficient `CL`, and the imposed circulation `Γ`. Validate against
[`flettner_analytic`](@ref): at `ω = 0` this is the non-lifting cylinder
`Cp = 1 − 4sin²θ`; grid-refine `N` to drive `ε(CL) < 0.02 %`.
"""
function flettner_panel(; R::Real=0.5, ω::Real=1.0, V∞::Real=1.0, N::Int=200)
    Γ = 2π * ω * R^2
    # CCW panel nodes on the circle; panel j from node j to node j+1.
    ϕ = range(0, 2π; length=N+1)                       # node angles
    xn = R .* cos.(ϕ); yn = R .* sin.(ϕ)
    # panel control points (midpoints) and outward normals / CCW tangents
    xc = Vector{Float64}(undef, N); yc = similar(xc)
    nx = similar(xc); ny = similar(xc)           # outward unit normal
    tx = similar(xc); ty = similar(xc)           # CCW unit tangent (panel dir)
    θc = similar(xc)
    for j in 1:N
        xc[j] = (xn[j] + xn[j+1]) / 2
        yc[j] = (yn[j] + yn[j+1]) / 2
        dx = xn[j+1] - xn[j]; dy = yn[j+1] - yn[j]
        L = sqrt(dx*dx + dy*dy)
        tx[j] = dx / L; ty[j] = dy / L            # CCW tangent
        nx[j] = ty[j];  ny[j] = -tx[j]            # outward normal (right of tangent)
        θc[j] = atan(yc[j], xc[j])
    end

    # No-penetration system: Σ_k A[j,k] σ_k = b_j,  A = normal source-panel
    # influence, b = −(freestream + prescribed vortex) · n̂.
    A = Matrix{Float64}(undef, N, N)
    b = Vector{Float64}(undef, N)
    for j in 1:N
        # RHS: freestream (V∞, 0) plus the central point vortex Γ.
        rx = xc[j]; ry = yc[j]; r2 = rx*rx + ry*ry
        # Central point vortex of strength Γ. The CCW panel ordering makes the
        # source-resolved freestream surface-tangential run *opposite* to the
        # standard analytic convention `Vt = 2V∞ sinθ + Γ/(2πR)`; to add
        # the bound circulation on the same (lift-producing) side as the
        # analytic form, the vortex velocity is taken as Γ/(2π)·(y, −x)/r²
        # (so a positive ω/Γ speeds the flow over the top → +lift, matching
        # the closed-form CL = +4πωR²/V∞). Verified against `flettner_analytic`.
        uvx =  Γ/(2π) * ry / r2
        uvy = -Γ/(2π) * rx / r2
        b[j] = -((V∞ + uvx) * nx[j] + uvy * ny[j])
        for k in 1:N
            u, v = _source_panel_vel(xc[j], yc[j], xn[k], yn[k], xn[k+1], yn[k+1])
            if j == k
                # self-influence of a source panel is +0.5 along its own normal
                A[j, k] = 0.5
            else
                A[j, k] = u * nx[j] + v * ny[j]
            end
        end
    end
    σ = A \ b

    # Tangential surface velocity at each control point, then Cp.
    Cp = Vector{Float64}(undef, N)
    for j in 1:N
        rx = xc[j]; ry = yc[j]; r2 = rx*rx + ry*ry
        uvx =  Γ/(2π) * ry / r2                    # see RHS note on the sign
        uvy = -Γ/(2π) * rx / r2
        ut = (V∞ + uvx) * tx[j] + uvy * ty[j]      # freestream + vortex, tangential
        for k in 1:N
            if k == j
                # self-tangential influence of a constant-source panel is 0
                continue
            end
            u, v = _source_panel_vel(xc[j], yc[j], xn[k], yn[k], xn[k+1], yn[k+1])
            ut += σ[k] * (u * tx[j] + v * ty[j])
        end
        Cp[j] = 1 - (ut / V∞)^2
    end

    # CL = ½ ∮ −Cp sinθ dθ over chord c = 2R (the standard lift-coefficient definition).
    dθ = 2π / N
    CL = 0.5 * sum(@. -Cp * sin(θc)) * dθ
    return (; θ = θc, Cp, CL, Γ)
end
