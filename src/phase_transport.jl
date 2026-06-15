# ----------------------------------------------------------------------------
# 1D two-phase advection — the "phase transport equation"
#
# The volume fraction α∈[0,1] (1 = water, 0 = air) advected by a uniform
# velocity u obeys the scalar transport equation
#
#     ∂α/∂t + ∇·(u α) = 0                                  (pure advection)
#
# discretised here with **first-order upwind** in space (backward difference
# for u>0) and **Forward Euler** in time. First-order upwind is monotone and
# bounded but numerically diffusive: a top-hat profile smears as it convects.
#
# The interFoam / "phase transport equation" cure is an *artificial
# compression* term that only acts at the interface (where α(1-α)≠0) and
# pushes α back toward the sharp 0/1 step:
#
#     ∂α/∂t + ∇·(u α) + ∇·(u_r α(1-α)) = 0,   u_r = C_α |u| n̂ ≈ C_α u sign(∂α/∂x)
#
# with the compression coefficient C_α (C_α=1 is the interFoam default). The
# compression flux is discretised conservatively (flux-difference form) so the
# scheme stays mass-conserving and the bounded variant stays in [0,1].
#
# This is the 1D sibling of `VoF.step_vof_mules!` (VoF.jl): there the same
# interface-compression flux `c_α·|u_f|·n̂_j·α_f(1-α_f)` is added per face in
# 2D/3D on a WaterLily MAC grid (MULES limiter for boundedness). Here we keep
# it 1D and explicit so the mechanism — and the diffusion-vs-sharpening
# trade-off — is transparent and testable against the analytic translate.
#
# Pure Julia (no deps beyond Base).
# ----------------------------------------------------------------------------

"""
    default_patch(x) -> α₀

The reference initial condition: a water patch `α=1` on `x∈[0.2,0.4]`,
air `α=0` elsewhere. `x` is a vector of cell-centre coordinates.
"""
default_patch(x::AbstractVector) = [0.2 ≤ xi ≤ 0.4 ? 1.0 : 0.0 for xi in x]

"""
    phase_transport_1d(; N=100, u=1.0, Co=0.5, t_end=0.2, L=1.0,
                        α₀=default_patch, compression=false, Cα=1.0)
        -> (; x, α, α_analytic, t, dt, dx, mass, mass0)

Solve the 1D phase-transport (two-phase advection) equation
`∂α/∂t + ∇·(u α) = 0` on `x∈[0,L]` with `N` finite-volume cells, advecting
the volume fraction `α∈[0,1]` at uniform velocity `u>0`.

Discretisation: **first-order upwind** (backward difference, valid for
`u>0`) + **Forward Euler**, time step set by the Courant number
`Co = u·Δt/Δx` (`Δt = Co·Δx/u`). The inlet (`x=0`) is held at `α=0`
(clean air entering); the outlet uses a zero-gradient (upwind) extrapolation.

With `compression=true` an **artificial interface-compression** term
`∇·(u_r α(1-α))`, `u_r = Cα·u·sign(∂α/∂x)`, is added (the interFoam /
"phase transport equation" sharpener). It only acts where `α(1-α)≠0` (the
interface), counteracting the upwind diffusion and re-steepening the front,
while a conservative flux discretisation keeps mass conserved and α bounded
in `[0,1]`.

Returns the cell centres `x`, the final field `α`, the **pure-advection
analytic solution** `α_analytic` (the initial profile translated by `u·t`,
with `α=0` flowing in at the inlet), the reached time `t`, `dt`, `dx`, and
the total mass `∫α dx` at the end vs the start (`mass`, `mass0`) for a
conservation check.

This is the 1D explicit sibling of [`VoF.step_vof_mules!`]; the compression
flux `Cα·u·α(1-α)` is the same term as that solver's per-face
`c_α·|u_f|·n̂·α_f(1-α_f)` interface-compression flux, here in one dimension
so the diffusion-vs-sharpening trade-off is directly testable against the
analytic translate.

Reference test case: `N=100, u=1, Co=0.5, t_end=0.2`,
water patch on `[0.2,0.4]`; at `t=0.2` the analytic profile is the patch
translated to `[0.4,0.6]`.
"""
function phase_transport_1d(; N::Int=100, u::Real=1.0, Co::Real=0.5,
                            t_end::Real=0.2, L::Real=1.0,
                            α₀=default_patch, compression::Bool=false,
                            Cα::Real=1.0)
    u > 0 || throw(ArgumentError("phase_transport_1d assumes u>0 (upwind = backward diff)"))
    dx = L / N
    # cell centres
    x = [(i - 0.5) * dx for i in 1:N]
    dt = Co * dx / u
    α = α₀ isa Function ? Float64.(α₀(x)) : Float64.(collect(α₀))
    length(α) == N || throw(ArgumentError("α₀ must have length N=$N"))
    α0_save = copy(α)
    mass0 = sum(α) * dx

    # Face fluxes at the N+1 faces k=1..N+1; face k sits between cell k-1
    # (left, upwind for u>0) and cell k (right). Inlet face k=1 sees the
    # boundary value α_in=0; outlet face k=N+1 uses a zero-gradient right cell.
    #
    # Plain advection: conservative first-order upwind + Forward Euler.
    # Compression variant: a MULES / FCT scheme — the upwind flux is the
    # bounded monotone base, the interface-compression flux
    #   Φc = Cα·u·n̂·α_f(1-α_f),   n̂ = sign(∂α/∂x)
    # is the anti-diffusive correction, and a 1D Zalesak limiter scales the
    # correction per face so no cell leaves [0,1] — boundedness AND exact
    # conservation (the same construction as `VoF.step_vof_mules!` in 3D).
    Fadv = zeros(Float64, N + 1)         # monotone upwind flux
    Fc   = zeros(Float64, N + 1)         # anti-diffusive compression correction
    nsteps = max(0, round(Int, t_end / dt))
    t = 0.0
    for _ in 1:nsteps
        @inbounds for k in 1:N+1
            αup = k == 1 ? 0.0 : α[k-1]
            Fadv[k] = u * αup
            if compression
                αR = k == N + 1 ? α[N] : α[k]
                αL = k == 1 ? 0.0 : α[k-1]
                αf = 0.5 * (αL + αR)
                s  = αf * (1 - αf)
                gj = αR - αL
                nhat = gj > 0 ? 1.0 : gj < 0 ? -1.0 : 0.0
                Fc[k] = Cα * u * nhat * s     # up-gradient (interface-compacting)
            end
        end

        if !compression
            @inbounds for i in 1:N
                α[i] -= (dt / dx) * (Fadv[i+1] - Fadv[i])
            end
        else
            # Zalesak FCT in 1D. (a) low-order (upwind-only) update.
            αLD = similar(α)
            @inbounds for i in 1:N
                αLD[i] = α[i] - (dt / dx) * (Fadv[i+1] - Fadv[i])
            end
            # (b) per-cell allowable rise/fall to stay in [0,1].
            # Anti-diffusive flux entering cell i = +Fc[i] (left face) − Fc[i+1].
            Ppos = zeros(Float64, N); Pneg = zeros(Float64, N)
            @inbounds for i in 1:N
                fin  = (dt / dx) * Fc[i]        # left face adds to cell i
                fout = (dt / dx) * Fc[i+1]      # right face removes from cell i
                Ppos[i] = max(0.0,  fin) + max(0.0, -fout)   # net possible rise
                Pneg[i] = max(0.0, -fin) + max(0.0,  fout)   # net possible fall
            end
            Qpos = [1.0 - αLD[i] for i in 1:N]   # room to grow up to 1
            Qneg = [αLD[i] - 0.0 for i in 1:N]   # room to shrink down to 0
            Rpos = [Ppos[i] > 0 ? clamp(Qpos[i]/Ppos[i], 0, 1) : 1.0 for i in 1:N]
            Rneg = [Pneg[i] > 0 ? clamp(Qneg[i]/Pneg[i], 0, 1) : 1.0 for i in 1:N]
            # (c) limit each interior face's correction by the tighter of the
            #     two cells it couples (Zalesak). The face-i correction Fc[i]
            #     leaves cell i-1 (right face of it) and enters cell i (left).
            @inbounds for k in 2:N
                if Fc[k] ≥ 0      # adds to cell k, removes from cell k-1
                    λ = min(Rpos[k], Rneg[k-1])
                else              # removes from cell k, adds to cell k-1
                    λ = min(Rneg[k], Rpos[k-1])
                end
                Fc[k] *= λ
            end
            Fc[1] = 0.0; Fc[N+1] = 0.0          # no compression through boundaries
            @inbounds for i in 1:N
                α[i] = αLD[i] - (dt / dx) * (Fc[i+1] - Fc[i])
            end
        end
        t += dt
    end

    # pure-advection analytic: initial profile shifted by u·t, α=0 in from inlet
    shift = u * t
    α_analytic = map(x) do xi
        xs = xi - shift
        xs < 0 ? 0.0 : (α₀ isa Function ? α₀([xs])[1] : _sample(α0_save, x, xs))
    end

    mass = sum(α) * dx
    return (; x, α, α_analytic, t, dt, dx, mass, mass0)
end

# nearest-cell sample of a discrete field for the analytic shift fallback
@inline function _sample(field::AbstractVector, x::AbstractVector, xq)
    (xq < x[1] || xq > x[end]) && return 0.0
    i = argmin(abs.(x .- xq))
    return field[i]
end

"""
    interface_width(x, α; lo=0.05, hi=0.95) -> w

Width of the diffuse interface: the x-distance over which `α` falls from
`hi` to `lo` across the (trailing) front. A crude monotone-region estimate
used to quantify upwind diffusion vs compression sharpening — smaller is
sharper. Returns `NaN` if the band is not crossed.
"""
function interface_width(x::AbstractVector, α::AbstractVector; lo=0.05, hi=0.95)
    # locate the rising or falling edge: use the steepest contiguous front.
    # We measure across the *leading* (downstream) edge of the translated
    # top-hat by scanning from the max-α cell rightward to where α drops below lo.
    imax = argmax(α)
    xhi = nothing; xlo = nothing
    for i in imax:length(α)-1
        if xhi === nothing && α[i] ≥ hi && α[i+1] < hi
            t = (hi - α[i]) / (α[i+1] - α[i])
            xhi = x[i] + t * (x[i+1] - x[i])
        end
        if α[i] ≥ lo && α[i+1] < lo
            t = (lo - α[i]) / (α[i+1] - α[i])
            xlo = x[i] + t * (x[i+1] - x[i])
            break
        end
    end
    (xhi === nothing || xlo === nothing) && return NaN
    return abs(xlo - xhi)
end
