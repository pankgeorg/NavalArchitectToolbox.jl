# ----------------------------------------------------------------------------
# Wind-assisted-propulsion (WAP) required-power-change analysis — NTUA-8401 Q5
# ITTC 7.5-04-01-02 (Analysis of Speed/Power Trials), reduced for a logged
# sea-trial dataset that records SHP, speed and apparent wind.
#
# WHAT THE TOOL DOES
# ------------------
# A sea trial of a wind-assist device is run as OFF → ON → OFF: a baseline
# segment with the system stowed, an assisted segment with it deployed, then
# baseline again, all at comparable speeds. The question is the *change in
# required propulsion power* attributable to the device. The confound is that
# the OFF and ON segments do not see the same natural wind, so part of any
# SHP difference is just the superstructure aerodynamic load, not the device.
# ITTC removes that load with a wind-resistance correction.
#
# We produce ΔP in TWO passes (as the assignment asks):
#   Pass 1 (raw):       ΔP from the measured SHP directly — fit P(V) to the
#                       OFF (baseline) samples, evaluate at the ON segment's
#                       mean speed, subtract the ON segment's mean SHP.
#   Pass 2 (corrected): same, but every sample's SHP first has the apparent-
#                       wind superstructure power load removed (steps 1–3
#                       below), so OFF and ON are compared at a common (zero
#                       apparent-wind) reference. This isolates the device
#                       from the ambient wind.
#
# ITTC STEPS IMPLEMENTED (vs simplified) — stated plainly:
#   1. Per-sample wind resistance on the superstructure (IMPLEMENTED):
#         R_AA = ½ · ρ_air · CDA(AWA) · A_T · V_wr²
#      with V_wr the apparent (relative) wind speed AWS in m/s, A_T the
#      projected frontal area, and CDA linearly interpolated from the supplied
#      AWA→CDA table (a wind-tunnel/Fujiwara-type drag coefficient). CDA is
#      tabulated 0–180°; AWA is logged signed (port/starboard), so we use
#      |AWA| — the load is symmetric about the bow–stern axis for a frontal-
#      area coefficient. This is the ITTC "reference to zero wind" wind-
#      resistance increment; we do NOT add the still-air-air-drag back in
#      (a true ITTC trial subtracts the relative-wind load and adds back the
#      ship-speed-only air drag; here the baseline reference is zero apparent
#      wind, which is the simplification — documented).
#   2. Resistance → power (SIMPLIFIED): only SHP is logged, not thrust, so we
#      convert the wind-resistance increment to a shaft-power increment with a
#      stated quasi-propulsive efficiency η_D (default 0.7):
#         ΔP_wind = R_AA · V_s / η_D,   V_s = STW in m/s.
#      A full ITTC analysis would use the measured thrust/SHP relation and the
#      load-variation (thrust-identity) method; with SHP-only logging that is
#      not available, so η_D is an explicit assumption.
#   3. Corrected shaft power (IMPLEMENTED): SHP_corr = SHP_meas − ΔP_wind.
#   4. Power–speed curve + ΔP (IMPLEMENTED): fit P = a·V^b by linear least
#      squares in log space to the OFF samples, evaluate at V̄_ON, and
#      ΔP = P_baseline(V̄_ON) − P̄_ON. Per run (A,B) and combined.
#
# NOT MODELLED (say so): added wave resistance, water-temperature/density and
# shallow-water corrections, current/drift, rudder-induced drag, and the
# heading/load drift between segments beyond what speed-matching + the wind
# correction capture. This is a reduced trials-analysis tool, not the full
# ITTC 7.5-04-01-02 procedure with all environmental corrections.
#
# Pure Julia. CSV via stdlib DelimitedFiles; CDA interpolation reuses the
# package's 1-D `_interp`. No new dependencies.
# ----------------------------------------------------------------------------

using DelimitedFiles: readdlm

const _KN2MS = 0.514444          # knots → m/s

# ----------------------------------------------------------------------------
# Reading the trial CSV and the wind-coefficient table
# ----------------------------------------------------------------------------

"""
    read_trial(path) -> (; t, sog, stw, rpm, shp, rudder, aws, awa, on)

Read a Θέμα-5 sea-trial CSV (stdlib `DelimitedFiles`). Expected columns:
`DateTime, SOG [kn], STW [kn], Shaft Speed [rpm], SHP [kW], Rudder Angle [deg],
AWS [kn], AWA [deg], Status_ON`. Speeds (`sog`, `stw`, `aws`) are returned in
**knots** (as logged); `on` is a `BitVector` (`true` where `Status_ON` is the
text `"True"`). `t` is the raw `DateTime` column (strings).
"""
function read_trial(path::AbstractString)
    data, _ = readdlm(path, ','; header=true)
    col(j) = Float64.(@view data[:, j])
    statuscol = @view data[:, 9]
    on = BitVector(uppercase(strip(string(s))) == "TRUE" for s in statuscol)
    return (; t      = string.(@view data[:, 1]),
              sog    = col(2),
              stw    = col(3),
              rpm    = col(4),
              shp    = col(5),
              rudder = col(6),
              aws    = col(7),
              awa    = col(8),
              on)
end

"""
    read_wind_coef(path) -> (; awa, cda)

Read the superstructure wind-resistance coefficient table
`AWA_deg, CDA_10m` (one header row). Returns vectors sorted by `awa`,
suitable for [`wind_resistance`](@ref).
"""
function read_wind_coef(path::AbstractString)
    data, _ = readdlm(path, ','; header=true)
    awa = Float64.(@view data[:, 1]); cda = Float64.(@view data[:, 2])
    p = sortperm(awa)
    return (; awa = awa[p], cda = cda[p])
end

# ----------------------------------------------------------------------------
# Per-sample superstructure wind resistance
# ----------------------------------------------------------------------------

"""
    wind_resistance(aws_ms, awa_deg, coef; A_T=121.4, ρ=1.225) -> R_AA

Superstructure aerodynamic resistance from the apparent wind:

    R_AA = ½ · ρ · CDA(|AWA|) · A_T · aws_ms²    [N]

`aws_ms` is the apparent (relative) wind speed in **m/s**, `awa_deg` the
apparent wind angle in degrees (0° = head wind; sign/side ignored — the
frontal-area coefficient is symmetric about the bow–stern axis). `coef` is
the `(; awa, cda)` table from [`read_wind_coef`](@ref); `CDA` is linearly
interpolated (clamped outside 0–180°). Positive for head winds (`CDA>0`,
added drag), negative for following winds (`CDA<0`, thrust), and zero at
`aws_ms = 0`.
"""
function wind_resistance(aws_ms::Real, awa_deg::Real, coef;
                         A_T::Real=121.4, ρ::Real=1.225)
    cda = _interp(coef.awa, coef.cda, abs(float(awa_deg)))
    return 0.5 * ρ * cda * A_T * aws_ms^2
end

# ----------------------------------------------------------------------------
# Power–speed curve fit:  P = a · V^b  (least squares in log space)
# ----------------------------------------------------------------------------

"""
    fit_power_speed(V, P) -> (; a, b, predict, r2)

Least-squares fit of the power law `P = a·V^b` to speed/power samples by a
linear regression of `log P` on `log V` (LinearAlgebra `\\`). Returns the
coefficients `a`, `b`, a closure `predict(v) = a·v^b`, and the coefficient of
determination `r2` of the log-space fit. Requires `V, P > 0`.
"""
function fit_power_speed(V::AbstractVector, P::AbstractVector)
    length(V) == length(P) || throw(DimensionMismatch("V and P length mismatch"))
    (all(>(0), V) && all(>(0), P)) ||
        throw(ArgumentError("fit_power_speed needs strictly positive V and P"))
    x = log.(float.(V)); y = log.(float.(P))
    A = hcat(ones(length(x)), x)          # [1  log V]
    coef = A \ y                          # [log a;  b]
    loga, b = coef[1], coef[2]
    a = exp(loga)
    ŷ = A * coef
    ss_res = sum(abs2, y .- ŷ)
    ss_tot = sum(abs2, y .- sum(y)/length(y))
    r2 = ss_tot == 0 ? 1.0 : 1 - ss_res/ss_tot
    return (; a, b, predict = v -> a * v^b, r2)
end

# ----------------------------------------------------------------------------
# Segment extraction (OFF → ON → OFF) with transient trimming
# ----------------------------------------------------------------------------

# Indices of the ON block (first..last true) trimmed by `trim` each side, and
# the OFF (baseline) indices outside the ON block, also trimmed away from the
# transitions. Returns (on_idx, off_idx).
function _segments(on::BitVector, trim::Int)
    i1 = findfirst(on); i2 = findlast(on)
    (i1 === nothing || i2 === nothing) &&
        throw(ArgumentError("no ON samples (Status_ON) found in trial"))
    n = length(on)
    on_idx  = collect((i1 + trim):(i2 - trim))
    # OFF baseline: samples not in the ON block, trimmed away from the
    # transition on the inner side.
    off_lo = collect(1:(i1 - 1 - trim))
    off_hi = collect((i2 + 1 + trim):n)
    off_idx = vcat(off_lo, off_hi)
    # keep only samples that are actually OFF (guards stray flips)
    on_idx  = filter(i -> 1 ≤ i ≤ n && on[i],  on_idx)
    off_idx = filter(i -> 1 ≤ i ≤ n && !on[i], off_idx)
    return on_idx, off_idx
end

_meanstd(x) = (m = sum(x)/length(x); (; mean = m,
    std = length(x) > 1 ? sqrt(sum(abs2, x .- m)/(length(x)-1)) : 0.0))

# ----------------------------------------------------------------------------
# Top-level analysis
# ----------------------------------------------------------------------------

"""
    wap_power_analysis(run_csv, coef_csv; A_T=121.4, ρ_air=1.225,
                       η_D=0.7, trim=4)
        -> NamedTuple

Required-propulsion-power change attributable to a wind-assist system from
one OFF→ON→OFF sea-trial run (`run_csv`) and the superstructure wind-
coefficient table (`coef_csv`), per a reduced ITTC 7.5-04-01-02 analysis.

The run is split into its baseline (OFF) and assisted (ON) segments; `trim`
samples are dropped on each side of every transition to remove transients.
A power-law speed curve `P=a·V^b` is fit to the OFF samples and evaluated at
the ON segment's mean speed; the change is
`ΔP = P_baseline(V̄_ON) − P̄_ON`, computed twice:

- `ΔP_raw`  — from the measured SHP (no wind correction).
- `ΔP_corr` — after removing the per-sample apparent-wind superstructure
  power load `ΔP_wind = R_AA·V_s/η_D` (see [`wind_resistance`](@ref)),
  so OFF and ON are compared at a common zero-apparent-wind reference.

A positive ΔP means the system **reduces** required power at the matched
speed (the baseline curve sits above the assisted segment).

Returns `(; ΔP_raw, ΔP_corr, V_on, P_on_raw, P_on_corr, P_base_raw,
P_base_corr, fit_raw, fit_corr, n_on, n_off, seg)` where `seg` carries the
OFF/ON segment statistics (mean/std of speed, SHP, AWS, AWA, wind power).
Speeds are STW in **knots**; powers in **kW** (as logged). Assumes η_D for
the resistance→power conversion and a P=a·V^b fit form — see the source
header for which ITTC steps are implemented vs simplified.
"""
function wap_power_analysis(run_csv::AbstractString, coef_csv::AbstractString;
                            A_T::Real=121.4, ρ_air::Real=1.225,
                            η_D::Real=0.7, trim::Int=4)
    tr   = read_trial(run_csv)
    coef = read_wind_coef(coef_csv)

    on_idx, off_idx = _segments(tr.on, trim)
    isempty(on_idx)  && throw(ArgumentError("ON segment empty after trimming (trim=$trim too large?)"))
    isempty(off_idx) && throw(ArgumentError("OFF baseline empty after trimming"))

    # per-sample wind power load ΔP_wind [kW] = R_AA[N]·V_s[m/s]/η_D / 1000
    aws_ms = tr.aws .* _KN2MS
    stw_ms = tr.stw .* _KN2MS
    R_AA = [wind_resistance(aws_ms[i], tr.awa[i], coef; A_T, ρ=ρ_air) for i in eachindex(aws_ms)]
    dP_wind = R_AA .* stw_ms ./ η_D ./ 1000           # kW

    # corrected SHP removes the natural-wind superstructure load
    shp_corr = tr.shp .- dP_wind

    Voff = tr.stw[off_idx]
    Von  = tr.stw[on_idx]
    V̄on  = sum(Von)/length(Von)

    # --- Pass 1 (raw) ---
    fit_raw = fit_power_speed(Voff, tr.shp[off_idx])
    P_base_raw = fit_raw.predict(V̄on)
    P_on_raw   = sum(tr.shp[on_idx])/length(on_idx)
    ΔP_raw     = P_base_raw - P_on_raw

    # --- Pass 2 (corrected) ---
    fit_corr = fit_power_speed(Voff, shp_corr[off_idx])
    P_base_corr = fit_corr.predict(V̄on)
    P_on_corr   = sum(shp_corr[on_idx])/length(on_idx)
    ΔP_corr     = P_base_corr - P_on_corr

    seg = (;
        off = (; n = length(off_idx),
                 V = _meanstd(Voff), shp = _meanstd(tr.shp[off_idx]),
                 shp_corr = _meanstd(shp_corr[off_idx]),
                 aws = _meanstd(tr.aws[off_idx]), awa = _meanstd(abs.(tr.awa[off_idx])),
                 dP_wind = _meanstd(dP_wind[off_idx])),
        on  = (; n = length(on_idx),
                 V = _meanstd(Von), shp = _meanstd(tr.shp[on_idx]),
                 shp_corr = _meanstd(shp_corr[on_idx]),
                 aws = _meanstd(tr.aws[on_idx]), awa = _meanstd(abs.(tr.awa[on_idx])),
                 dP_wind = _meanstd(dP_wind[on_idx])),
    )

    return (; ΔP_raw, ΔP_corr, V_on = V̄on,
              P_on_raw, P_on_corr, P_base_raw, P_base_corr,
              fit_raw, fit_corr, n_on = length(on_idx), n_off = length(off_idx),
              seg)
end

"""
    wap_power_analysis(run_csvs::AbstractVector, coef_csv; kw...) -> NamedTuple

Combined analysis over several OFF→ON→OFF runs: pools every run's OFF samples
into one baseline `P=a·V^b` fit and every run's ON samples into one assisted
mean, then computes `ΔP_raw`/`ΔP_corr` the same way. Returns the combined
result with `per` = the vector of single-run results (one per `run_csvs`
entry). Same keyword options as the single-run method.
"""
function wap_power_analysis(run_csvs::AbstractVector{<:AbstractString},
                            coef_csv::AbstractString;
                            A_T::Real=121.4, ρ_air::Real=1.225,
                            η_D::Real=0.7, trim::Int=4)
    per = [wap_power_analysis(rc, coef_csv; A_T, ρ_air, η_D, trim) for rc in run_csvs]
    coef = read_wind_coef(coef_csv)

    Voff = Float64[]; Poff = Float64[]; Poff_c = Float64[]
    Von  = Float64[]; Pon  = Float64[]; Pon_c  = Float64[]
    for rc in run_csvs
        tr = read_trial(rc)
        on_idx, off_idx = _segments(tr.on, trim)
        aws_ms = tr.aws .* _KN2MS; stw_ms = tr.stw .* _KN2MS
        R_AA = [wind_resistance(aws_ms[i], tr.awa[i], coef; A_T, ρ=ρ_air) for i in eachindex(aws_ms)]
        dP_wind = R_AA .* stw_ms ./ η_D ./ 1000
        shp_corr = tr.shp .- dP_wind
        append!(Voff, tr.stw[off_idx]); append!(Poff, tr.shp[off_idx]); append!(Poff_c, shp_corr[off_idx])
        append!(Von,  tr.stw[on_idx]);  append!(Pon,  tr.shp[on_idx]);  append!(Pon_c,  shp_corr[on_idx])
    end

    V̄on = sum(Von)/length(Von)
    fit_raw  = fit_power_speed(Voff, Poff)
    fit_corr = fit_power_speed(Voff, Poff_c)
    P_base_raw  = fit_raw.predict(V̄on);  P_on_raw  = sum(Pon)/length(Pon)
    P_base_corr = fit_corr.predict(V̄on); P_on_corr = sum(Pon_c)/length(Pon_c)

    return (; ΔP_raw = P_base_raw - P_on_raw, ΔP_corr = P_base_corr - P_on_corr,
              V_on = V̄on, P_on_raw, P_on_corr, P_base_raw, P_base_corr,
              fit_raw, fit_corr, n_on = length(Von), n_off = length(Voff), per)
end
