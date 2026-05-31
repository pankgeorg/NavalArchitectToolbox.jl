# ----------------------------------------------------------------------------
# Key-blade vortex-ring propeller VLM (Anevlavi, Zafeiris, Papadakis,
# Belibassakis, JMSE 2023, §2.1). Shaft along +z; inflow Va along +z;
# rotation ω about +z. Solve ONE (key) blade; the other Z−1 blades enter
# through rotation symmetry.
#
# Key details (the first is the most error-prone):
#  (1) WAKE = a cylindrical helix winding the SAME sense as the rotation
#      (+ψ as z increases), so it trails off cleanly behind the blade and never
#      grazes the control points (a flow-aligned helix at the advance pitch
#      winds the *opposite* way, cuts back across the surface, and makes the
#      AIC singular). The default `:transition` wake leaves the TE at the
#      blade's geometric pitch and eases to mean(geometric, advance) downstream
#      (Kerwin transition wake); `:geometric` keeps a constant blade pitch.
#      Wake rings carry the adjacent TE ring's strength (Kutta, Kelvin).
#  (2) FORCES = pressure-normal, NOT Kutta–Joukowski. Δc_p = C_LE·2·V_m·G/V∞²
#      (steady-Bernoulli vorticity jump), and the force on each ring is
#      Δp·A·n̂ — directed along the *surface normal* (eqs 5–7). Thrust is its
#      axial projection, torque its moment about the shaft. Taking the force
#      direction from the geometry (not the local velocity) gives the correct
#      η; C_LE (0.85–0.95) restores the leading-edge suction the pressure-
#      normal model omits, C_Drag adds the friction torque.
# ----------------------------------------------------------------------------

using LinearAlgebra: ×, ⋅, norm

# Biot–Savart velocity at p from a unit-strength straight vortex segment a→b,
# with a Rankine-style vortex core `rc` (length). The core regularizes the
# 1/distance singularity (|r1×r2|² → |r1×r2|² + (rc·|r0|)²), essential here
# because the physically-correct helical wake passes near control points.
@inline function _biot_segment(p::SVector{3,T}, a, b, rc::T) where T
    r1 = p - a; r2 = p - b; r0 = b - a
    cr = r1 × r2
    n1 = norm(r1); n2 = norm(r2); l0 = norm(r0)
    (n1 < T(1e-9) || n2 < T(1e-9) || l0 < T(1e-12)) && return zero(SVector{3,T})
    denom = sum(abs2, cr) + (rc*l0)^2
    K = (one(T)/(4π)) * (r0⋅r1/n1 - r0⋅r2/n2) / denom
    return K * cr
end

# Velocity at p from a unit-strength quadrilateral vortex ring (core rc).
@inline _ring_induced(p, c, rc) =
    _biot_segment(p, c[1], c[2], rc) + _biot_segment(p, c[2], c[3], rc) +
    _biot_segment(p, c[3], c[4], rc) + _biot_segment(p, c[4], c[1], rc)

# Rotate a point about the shaft (z) axis by angle ψ.
@inline function _rotz(p::SVector{3,T}, ψ) where T
    s, c = sincos(T(ψ))
    SVector{3,T}(c*p[1] - s*p[2], s*p[1] + c*p[2], p[3])
end

# A vortex element: 4 corners + the index of the unknown circulation it
# carries (bound ring → itself; wake ring → its trailing-edge ring).
struct VElem{T}
    c::NTuple{4,SVector{3,T}}
    parent::Int
end

"VLM geometry/state for one key blade (+ its helical wake), Z blades total."
struct PropellerVLM{T}
    elems :: Vector{VElem{T}}          # bound rings + wake rings
    cp    :: Vector{SVector{3,T}}      # collocation points (one per bound ring)
    nrm   :: Vector{SVector{3,T}}      # panel normals at cp
    bseg  :: Vector{NTuple{2,SVector{3,T}}}  # leading bound-segment endpoints
    area  :: Vector{T}                 # panel areas
    nc :: Int; ns :: Int; Z :: Int
    Va :: T; ω :: T; nRPS :: T; D :: T
    rc :: T                            # bound-ring vortex-core radius (length)
    wrc :: T                           # wake-ring vortex-core radius (length)
    nb :: Int                          # number of bound rings (first nb in elems)
end

# Onset (relative) velocity at p: axial advance Va (+z) plus the rotational
# inflow. The tangential sense is set so the relative flow runs from the
# leading edge to the trailing edge of the tabulated blade (whose chord wraps
# +θ from LE→TE); i.e. the rotation/handedness convention is matched to the
# geometry. |V| = √(Va² + (ωr)²) either way, so open-water KT/KQ are unchanged.
@inline _inflow(p::SVector{3,T}, Va, ω) where T = SVector{3,T}(-ω*p[2], ω*p[1], Va)

function _build_vlm(tab::PropellerBladeTable{T}, D, Dh, J;
                    nc::Int=8, ns::Int=15, nRPS::Real=1.0,
                    meanline::MeanLine=_DEFAULT_MEANLINE, tip_rR::Real=0.97,
                    Kw::Int=32, dψ::Real=deg2rad(20.0), wake_pitch::Symbol=:transition,
                    core_frac::Real=0.1, wake_core_frac::Real=0.1,
                    pitch_corr::Real=1.9454) where T
    R = D/2; ω = T(2π*nRPS); Va = T(J*nRPS*D)
    cg = vlm_camber_grid(tab, D, Dh; nc=nc, ns=ns, meanline=meanline, tip_rR=tip_rR,
                         pitch_corr=pitch_corr)
    # quarter-chord corner lines: QL[i] at 1/4 of panel i; QL[nc+1] = TE
    QL = Matrix{SVector{3,T}}(undef, nc+1, ns+1)
    @inbounds for j in 1:ns+1
        for i in 1:nc
            QL[i,j] = cg[i,j] + T(0.25)*(cg[i+1,j] - cg[i,j])
        end
        QL[nc+1,j] = cg[nc+1,j]
    end
    elems = VElem{T}[]; cp = SVector{3,T}[]; nrm = SVector{3,T}[]
    bseg = NTuple{2,SVector{3,T}}[]; area = T[]
    lin(i,j) = (j-1)*nc + i
    @inbounds for j in 1:ns, i in 1:nc
        c = (QL[i,j], QL[i+1,j], QL[i+1,j+1], QL[i,j+1])
        push!(elems, VElem{T}(c, lin(i,j)))
        pl = cg[i,j]   + T(0.75)*(cg[i+1,j]   - cg[i,j])
        pr = cg[i,j+1] + T(0.75)*(cg[i+1,j+1] - cg[i,j+1])
        push!(cp, (pl+pr)/2)
        d1 = cg[i+1,j+1] - cg[i,j]; d2 = cg[i,j+1] - cg[i+1,j]
        nv = d1 × d2; push!(nrm, nv/norm(nv))
        push!(bseg, (QL[i,j], QL[i,j+1]))
        push!(area, norm(d1 × d2)/2)
    end
    # helical wake: per spanwise strip, march the two TE edge points
    # downstream by Δψ in azimuth, accumulating axial advance z (= ∫ r·tanβ dψ).
    # wake pitch (axial advance per revolution): :geometric = the blade pitch
    # P=(P/D)·D; :advance = the undisturbed advance Va/n; :hydro = geometric
    # mean of the two; :transition = ease the pitch *angle* from the blade
    # geometric pitch at the TE to mean(geometric, advance) downstream (a
    # Kerwin-style transition wake). Wake rings carry the TE ring's Γ.
    Padv = Va / T(nRPS)
    wpitch(rR) = (Pg = _interp(tab.rR, tab.PD, rR)*D;
                  wake_pitch === :geometric ? Pg :
                  wake_pitch === :advance   ? Padv : sqrt(Pg*Padv))
    # blade geometric pitch *angle* (DDFI-corrected, matching the blade TE)
    wβ(rR) = atan(_interp(tab.rR, tab.PD, rR)/(π*rR)) -
             T(pitch_corr)*_interp(tab.rR, tab.fc, rR)*_interp(tab.rR, tab.tc, rR)
    trans = wake_pitch === :transition
    @inbounds for j in 1:ns
        te_l = cg[nc+1,j]; te_r = cg[nc+1,j+1]
        rRl = hypot(te_l[1],te_l[2])/R; rRr = hypot(te_r[1],te_r[2])/R
        rl = rRl*R; rr = rRr*R
        βgl = wβ(rRl); βgr = wβ(rRr)                       # near-wake (geometric)
        βul = T(0.5)*(βgl + atan(Va/(ω*rl)))               # ultimate = mean(geom, advance)
        βur = T(0.5)*(βgr + atan(Va/(ω*rr)))
        cl = wpitch(rRl)/T(2π); cr = wpitch(rRr)/T(2π)     # constant-pitch advance/dψ
        prev_l = te_l; prev_r = te_r; zl = te_l[3]; zr = te_r[3]
        for m in 1:Kw
            f = Kw > 1 ? T((m-1)/(Kw-1)) : zero(T)
            zl += dψ * (trans ? rl*tan(βgl + f*(βul-βgl)) : cl)
            zr += dψ * (trans ? rr*tan(βgr + f*(βur-βgr)) : cr)
            pl = _rotz(te_l, m*dψ); pr = _rotz(te_r, m*dψ)
            nl = SVector{3,T}(pl[1], pl[2], zl); nr = SVector{3,T}(pr[1], pr[2], zr)
            push!(elems, VElem{T}((prev_l, nl, nr, prev_r), lin(nc,j)))
            prev_l = nl; prev_r = nr
        end
    end
    # vortex-core radii: small fractions of the mean panel length (√mean-area),
    # a light Rankine regularization. `wake_core_frac` is a separate knob
    # (default = the bound core): the geometric-pitch wake does not graze the
    # control points, so no extra wake smearing is needed by default.
    ℓ = sqrt(sum(area)/length(area))
    rc = T(core_frac) * ℓ
    wrc = T(wake_core_frac) * ℓ
    PropellerVLM{T}(elems, cp, nrm, bseg, area, nc, ns, tab.Z, Va, ω, T(nRPS), T(D),
                    rc, wrc, nc*ns)
end

# Total induced velocity at p from all elements, with per-element Γ given
# by Γvec[parent], summed over the Z blade images (rotation symmetry).
function _induced(vlm::PropellerVLM{T}, p, Γvec; exclude::Int=0) where T
    v = zero(SVector{3,T})
    @inbounds for (idx, e) in enumerate(vlm.elems)
        e.parent == exclude && continue          # skip the self panel (KJ)
        g = Γvec[e.parent]
        g == 0 && continue
        rc = idx <= vlm.nb ? vlm.rc : vlm.wrc
        for b in 0:vlm.Z-1
            ψ = 2π*b/vlm.Z
            cb = (_rotz(e.c[1],ψ), _rotz(e.c[2],ψ), _rotz(e.c[3],ψ), _rotz(e.c[4],ψ))
            v += g * _ring_induced(p, cb, rc)
        end
    end
    return v
end

"""
    openwater_vlm(tab, D, Dh, J; nc=10, ns=18, C_LE=0.80, C_Drag=0.010, …)

Key-blade vortex-ring VLM open-water point for the `tab.Z`-blade propeller
at advance ratio `J`. Returns a `NamedTuple`:

| field | meaning |
|-------|---------|
| `KT`, `KQ`, `η` | thrust/torque coefficients and open-water efficiency |
| `Γ`   | bound vortex-ring strengths (length `nc·ns`) |
| `cp`  | key-blade panel collocation points (3/4-chord) |
| `Vm`  | mean total velocity vector at each panel |
| `dcp` | pressure-difference coefficient Δc_p per panel |
| `force` | pressure-normal force vector per panel (key blade) |
| `nc`, `ns` | chordwise / spanwise panel counts |

The result is mesh- and wake-converged. Defaults turn on the **transition
wake** (`wake_pitch=:transition` — pitch eases from the blade's geometric
pitch at the TE to mean(geometric, advance) downstream) and a **lifting-
surface pitch correction** (`pitch_corr=1.9454`, subtracting a no-lift-angle
term ∝ (f/c)·(t/c)); see the module header. `C_LE` (leading-edge suction) and
`C_Drag` (friction) are the Anevlavi-Belibassakis calibration coefficients.
With the defaults `C_LE=0.80`, `C_Drag=0.010` the DTMB 4382 open-water point
at J=0.889 is KT≈0.207, 10·KQ≈0.444, η≈0.659 vs experimental 0.208 / 0.445 /
0.661, for this package's section approximations (exact NACA a=0.8 mean line +
a NACA-4-digit thickness stand-in for NACA 66-mod).

# Example
```julia
openwater_vlm(dtmb4382, 6.0, 1.2, 0.889)   # → (KT≈0.21, KQ≈0.045, η≈0.66, Γ)
```
"""
function openwater_vlm(tab::PropellerBladeTable{T}, D, Dh, J;
                       nc::Int=10, ns::Int=18, nRPS::Real=1.0,
                       meanline::MeanLine=_DEFAULT_MEANLINE, tip_rR::Real=0.97,
                       Kw::Int=96, dψ::Real=deg2rad(15.0), wake_pitch::Symbol=:transition,
                       core_frac::Real=0.1, wake_core_frac::Real=0.1,
                       pitch_corr::Real=1.9454,
                       C_LE::Real=0.80, C_Drag::Real=0.010, ρ::Real=1.0) where T
    vlm = _build_vlm(tab, D, Dh, J; nc, ns, nRPS, meanline, tip_rR, Kw, dψ, wake_pitch,
                     core_frac, wake_core_frac, pitch_corr)
    N = nc*ns
    # Influence matrix: A[i,k] = Σ_elements(parent=k) Σ_blades  v_ring·n_i
    A = zeros(T, N, N); rhs = zeros(T, N)
    @inbounds for i in 1:N
        ni = vlm.nrm[i]; pi = vlm.cp[i]
        for (idx, e) in enumerate(vlm.elems)
            rc = idx <= vlm.nb ? vlm.rc : vlm.wrc
            acc = zero(SVector{3,T})
            for b in 0:vlm.Z-1
                ψ = 2π*b/vlm.Z
                cb = (_rotz(e.c[1],ψ), _rotz(e.c[2],ψ), _rotz(e.c[3],ψ), _rotz(e.c[4],ψ))
                acc += _ring_induced(pi, cb, rc)
            end
            A[i, e.parent] += acc ⋅ ni
        end
        rhs[i] = -_inflow(pi, vlm.Va, vlm.ω) ⋅ ni
    end
    Γ = A \ rhs

    # Pressure-normal forces (Anevlavi-Belibassakis eqs 5–7). The chordwise
    # vorticity jump G = (Γ_i − Γ_{i−1})/δs gives Δc_p = C_LE·2·V_m·G/V∞²; the
    # force on a ring is Δp·A·n̂ along the surface normal. With Δp·A = ρ·V_m·G·A
    # and A/δs = (spanwise width) = |bound segment|, the ring force magnitude is
    # ρ·V_m·(Γ_i−Γ_{i−1})·|dl_span|, directed along n̂. Thrust = axial (z)
    # projection; torque = its moment about the shaft (z). Friction (C_Drag)
    # acts tangentially → adds torque, cuts thrust.
    lin(i,j) = (j-1)*nc + i
    Tp = zero(T); Qp = zero(T); T_fr = zero(T); Q_fr = zero(T)
    ρT = T(ρ); CLE = T(C_LE); CD = T(C_Drag)
    # per-panel post-processing fields (key blade): collocation point, Δc_p,
    # mean velocity, radius — for loading/pressure plots (paper Figs 7,8).
    cp    = Vector{SVector{3,T}}(undef, N)
    Vm_   = Vector{SVector{3,T}}(undef, N)
    dcp   = Vector{T}(undef, N)
    force = Vector{SVector{3,T}}(undef, N)         # pressure-normal force / panel
    @inbounds for j in 1:vlm.ns, i in 1:nc
        k = lin(i,j)
        γb = Γ[k] - (i==1 ? zero(T) : Γ[lin(i-1,j)])   # net chordwise vorticity jump
        a, b = vlm.bseg[k]; span = norm(b - a)         # spanwise panel width
        pc = vlm.cp[k]; ni = vlm.nrm[k]
        Vloc = _inflow(pc, vlm.Va, vlm.ω) + _induced(vlm, pc, Γ)  # mean surface vel
        Vm = norm(Vloc)
        Fp = CLE * ρT * Vm * γb * span                 # |Δp·A|; force = Fp·n̂
        Tp += Fp * ni[3]                               # axial thrust
        Qp += Fp * (pc[1]*ni[2] - pc[2]*ni[1])         # moment about z
        r  = hypot(pc[1], pc[2])
        V∞2 = vlm.Va^2 + (vlm.ω*r)^2
        cp[k] = pc; Vm_[k] = Vloc; force[k] = Fp * ni
        dcp[k] = CLE * 2 * Vm * γb / (vlm.area[k]/span) / V∞2   # Δc_p (G=γb/δs)
        # friction drag ½ρV_m²·A·C_Drag, tangential
        if Vm > 0
            Vt = hypot(Vloc[1], Vloc[2])
            qD = T(0.5)*ρT*Vm^2 * vlm.area[k] * CD
            Q_fr += qD * (Vt/Vm) * r                   # torque (resisting)
            T_fr += qD * (abs(Vloc[3])/Vm)             # thrust loss
        end
    end
    nbl = vlm.Z
    Tt = nbl * (abs(Tp) - T_fr)                        # thrust (axial)
    Qt = nbl * (abs(Qp) + Q_fr)                        # shaft torque
    n = vlm.nRPS
    KT = Tt / (ρT*n^2*D^4)
    KQ = Qt / (ρT*n^2*D^5)
    η  = KT/KQ * T(J)/(2π)
    return (; KT, KQ, η, Γ, cp, Vm=Vm_, dcp, force, nc, ns)
end
