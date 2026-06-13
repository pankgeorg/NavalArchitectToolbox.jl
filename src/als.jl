# ----------------------------------------------------------------------------
# Air-lubrication / air-layer drag reduction — reduced 1D model (NTUA-8401 Q2α)
#
# Air lubrication reduces a ship's frictional drag by interposing a thin gas
# film between the hull and the water. The full physics (bubble breakup, film
# stability, the BDR→ALDR transition) is an active research topic; here we
# build the simplest *honest* reduced model that captures the trend and the
# order of magnitude, and is reusable as a sizing tool.
#
# MODEL (a two-layer fully-developed shear stack — the implementer's
# defensible choice, stated plainly):
#
#   • A turbulent water boundary layer of thickness δ carries the wall shear.
#     Without air, the clean wall shear is the flat-plate value
#         τ_clean = ½ ρ_w U² C_f,         C_f = ITTC-1957 friction line.
#   • With a continuous air film of thickness t_air at the wall, the film is a
#     thin laminar Couette layer: it must transmit the *same* outer-flow shear
#     stress τ (shear continuity across the air/water interface in fully-
#     developed flow), but because the air viscosity μ_a ≪ μ_w the film carries
#     that shear with a much smaller velocity gradient. The wall shear felt by
#     the hull through the air film, for a given outer boundary-layer state, is
#     reduced because part of the near-wall velocity jump is taken up across the
#     low-viscosity film instead of against the wall.
#
#   Concretely, model the near-wall region as two stacked Couette layers that
#   share the interface velocity and the shear stress τ:
#         water layer:  τ = μ_w (U_δ − U_i)/δ
#         air film:     τ = μ_a (U_i − 0)/t_air            (no-slip at the wall)
#   Eliminating the interface velocity U_i for a fixed outer driving velocity
#   U_δ gives the with-film wall shear
#         τ_air / τ_clean_couette = 1 / (1 + (μ_w/μ_a)(t_air/δ)) ,
#   i.e. the air film acts as an added low-viscosity series resistance that
#   carries the velocity jump cheaply. Because μ_w/μ_a ≈ 55 (μ_w≈1e-3,
#   μ_a≈1.8e-5 Pa·s), even a film a few percent of δ gives a large reduction
#   once it is *continuous* — the qualitative ALDR signature.
#
#   We anchor the absolute clean shear to the turbulent flat-plate ITTC C_f
#   (so the ship estimate is realistic) and apply the two-layer ratio as the
#   *fractional* reduction:
#         DR = (τ_clean − τ_air)/τ_clean = (μ_w/μ_a)(t_air/δ) /
#                                          (1 + (μ_w/μ_a)(t_air/δ)).
#
# LIMITATIONS (say so): this is a 1D fully-developed reduced model. It assumes
# a *continuous* stable film (so it really models the air-layer (ALDR) regime,
# not the dispersed-bubble (BDR) regime), ignores film breakup, wave drag at
# the interface, gravity drainage on non-horizontal surfaces, and the air
# supply power. It will NOT reproduce Elbing et al.'s exact numbers; it
# reproduces the TREND (DR rises with t_air/δ and saturates toward 100% once a
# thick continuous film forms) and the order of magnitude (tens of % at a few-%
# thickness ratio). Compare to Elbing et al., JFM 612 (2008): BDR gives modest
# (and downstream-decaying) reduction, while a continuous air layer (ALDR)
# gives ≥80% once the critical air flux is exceeded.
#
# Pure Julia (no deps beyond Base).
# ----------------------------------------------------------------------------

# Water/air properties at ~15 °C (SI). μ in Pa·s, ρ in kg/m³, ν in m²/s.
const _MU_WATER  = 1.14e-3
const _MU_AIR    = 1.81e-5
const _RHO_WATER = 999.0
const _NU_WATER  = _MU_WATER / _RHO_WATER

"""
    ittc_cf(Re) -> C_f

ITTC-1957 model–ship correlation (frictional resistance) line
`C_f = 0.075 / (log10(Re) − 2)²`. `Re = U·L/ν`.
"""
ittc_cf(Re::Real) = 0.075 / (log10(Re) - 2)^2

"""
    als_drag_reduction(; delta, t_air, mu_w=μ_water, mu_a=μ_air,
                       rho_w=ρ_water, U=10.0, Re=nothing)
        -> (; DR, tau_clean, tau_air, ratio, mu_ratio, thick_ratio, Cf)

Reduced 1D **air-layer drag-reduction** estimate for a continuous air film
of thickness `t_air` (m) under a turbulent water boundary layer of thickness
`delta` (m), at free-stream speed `U` (m/s).

Two-layer fully-developed shear model (see source notes): the air film and
water boundary layer carry the same shear; the low-viscosity film
(`μ_w/μ_a ≈ 55`) takes up the near-wall velocity jump cheaply, so the
fractional wall-shear reduction is

    DR = r / (1 + r),   r = (μ_w/μ_a)·(t_air/δ).

The absolute clean wall shear `τ_clean = ½ ρ_w U² C_f` uses the ITTC-1957
friction line (`Re = U·L/ν`; pass `Re` directly, else a 100 m length is
assumed for the friction-line anchor only — the DR fraction is
Re-independent). `τ_air = τ_clean·(1−DR)`.

Models the **air-layer (ALDR) regime** of a *continuous* film, not the
dispersed-bubble (BDR) regime; ignores film breakup, interface waves,
drainage, and air-supply power. Trend/order-of-magnitude tool, not a
quantitative match to Elbing et al. (2008) — see [`als_sweep`] and the
package notes.
"""
function als_drag_reduction(; delta::Real, t_air::Real,
                            mu_w::Real=_MU_WATER, mu_a::Real=_MU_AIR,
                            rho_w::Real=_RHO_WATER, U::Real=10.0,
                            Re::Union{Real,Nothing}=nothing)
    delta > 0 || throw(ArgumentError("delta must be > 0"))
    t_air ≥ 0 || throw(ArgumentError("t_air must be ≥ 0"))
    mu_ratio   = mu_w / mu_a
    thick_ratio = t_air / delta
    r   = mu_ratio * thick_ratio
    DR  = r / (1 + r)                      # ∈[0,1), monotone in t_air
    Re_ = Re === nothing ? U * 100.0 / (mu_w / rho_w) : Re
    Cf  = ittc_cf(Re_)
    tau_clean = 0.5 * rho_w * U^2 * Cf
    tau_air   = tau_clean * (1 - DR)
    return (; DR, tau_clean, tau_air, ratio = 1 - DR,
            mu_ratio, thick_ratio, Cf)
end

"""
    als_sweep(; delta, t_airs, kw...) -> Vector{NamedTuple}
    als_sweep(; deltas, t_air, kw...) -> Vector{NamedTuple}

Sweep helper. Provide a vector `t_airs` (vary film thickness at fixed `delta`)
or a vector `deltas` (vary boundary-layer thickness at fixed `t_air`); the
non-swept companion is a scalar keyword. Each row is the full
[`als_drag_reduction`](@ref) result with the swept value attached.
"""
function als_sweep(; delta::Union{Real,Nothing}=nothing,
                   t_air::Union{Real,Nothing}=nothing,
                   t_airs::Union{AbstractVector,Nothing}=nothing,
                   deltas::Union{AbstractVector,Nothing}=nothing, kw...)
    if t_airs !== nothing
        delta === nothing && throw(ArgumentError("sweep over t_airs needs scalar `delta`"))
        return [merge((; t_air = ta), als_drag_reduction(; delta, t_air = ta, kw...))
                for ta in t_airs]
    elseif deltas !== nothing
        t_air === nothing && throw(ArgumentError("sweep over deltas needs scalar `t_air`"))
        return [merge((; delta = d), als_drag_reduction(; delta = d, t_air, kw...))
                for d in deltas]
    else
        throw(ArgumentError("provide either `t_airs` (with `delta`) or `deltas` (with `t_air`)"))
    end
end

"""
    als_ship_saving(; L, B, T, Cb=0.6, U, frac_covered=1.0,
                    delta=nothing, t_air, kw...)
        -> (; S, Re, Cf, Rf_clean, Rf_air, DR, frac_covered, delta, dRf_kN, saving_pct)

Apply the reduced air-layer DR to a ship's **frictional** resistance.

Estimates the wetted surface `S` (Holtrop-style approximation
`S ≈ L(2T+B)·√Cb`), the length-Reynolds number `Re = U·L/ν`, the ITTC
flat-plate friction coefficient `C_f`, and the clean frictional resistance
`Rf_clean = ½ ρ_w U² S C_f`. The air layer is applied over a fraction
`frac_covered` of `S` (typically a flat-bottom patch) with film thickness
`t_air`; the boundary-layer thickness defaults to a turbulent flat-plate
estimate `δ ≈ 0.16 L Re^(-1/7)` at the mid-body if `delta` is not given.
Returns the clean/with-air friction resistances, the DR fraction, the
absolute saving (kN), and the percentage saving on the friction component.

A reduced sizing estimate (friction component only; excludes air-supply
power, wave/residuary resistance, and film-stability losses).
"""
function als_ship_saving(; L::Real, B::Real, T::Real, Cb::Real=0.6,
                         U::Real, frac_covered::Real=1.0,
                         delta::Union{Real,Nothing}=nothing, t_air::Real,
                         rho_w::Real=_RHO_WATER, nu_w::Real=_NU_WATER, kw...)
    S  = L * (2T + B) * sqrt(Cb)            # wetted-surface approximation
    Re = U * L / nu_w
    Cf = ittc_cf(Re)
    Rf_clean = 0.5 * rho_w * U^2 * S * Cf
    δ = delta === nothing ? 0.16 * L * Re^(-1/7) : delta   # turbulent BL at ~mid-body
    dr = als_drag_reduction(; delta = δ, t_air, U, Re, rho_w, kw...)
    # DR applies only over the covered fraction of the wetted area
    Rf_air = Rf_clean * (1 - frac_covered * dr.DR)
    dRf = Rf_clean - Rf_air
    return (; S, Re, Cf, Rf_clean, Rf_air, DR = dr.DR, frac_covered,
            delta = δ, dRf_kN = dRf / 1e3, saving_pct = 100 * dRf / Rf_clean)
end
