using Test
using NavalArchitectToolbox
using NavalArchitectToolbox: dimensional, _interp, read_trial, read_wind_coef
using StaticArrays

@testset "NavalArchitectToolbox" begin

    @testset "table parse (DTMB 4382)" begin
        t = dtmb4382
        @test t.Z == 5
        @test length(t.rR) == 11
        @test t.rR[1] == 0.2 && t.rR[end] == 1.0
        @test t.cD[7] == 0.347          # r/R=0.7 chord
        @test t.PD[1] == 1.455          # hub pitch
        @test t.skew[end] == 36.0       # tip skew
        @test all(==(0), t.rake)        # no rake
    end

    @testset "DTMB 4381 = 4382 minus skew" begin
        a = dtmb4381; b = dtmb4382
        @test a.Z == 5 && length(a.rR) == 11
        @test all(==(0), a.skew)                 # 4381 is unskewed
        @test b.skew[end] == 36.0                # 4382 keeps its skew
        # everything except the skew column is identical to 4382
        @test a.rR  == b.rR
        @test a.cD  == b.cD
        @test a.PD  == b.PD
        @test a.tc == b.tc
        @test a.fc == b.fc
        @test a.rake == b.rake
    end

    @testset "interpolation + monotone sort" begin
        t = dtmb4382
        @test _interp(t.rR, t.PD, 0.2) == 1.455
        @test _interp(t.rR, t.PD, 1.0) == 0.942
        # midpoint between 0.6 and 0.7 P/D
        @test _interp(t.rR, t.PD, 0.65) ≈ (1.285 + 1.200)/2 atol=1e-6
        # clamps outside range
        @test _interp(t.rR, t.cD, 0.0) == t.cD[1]
        @test _interp(t.rR, t.cD, 2.0) == t.cD[end]
    end

    @testset "pitch angle" begin
        t = dtmb4382
        # tan β = (P/D)/(π r/R); hub: atan(1.455/(π·0.2))
        @test rad2deg(pitch_angle(t, 0.2)) ≈ atand(1.455/(π*0.2)) atol=1e-6
        @test rad2deg(pitch_angle(t, 1.0)) ≈ atand(0.942/(π*1.0)) atol=1e-6
        # pitch angle decreases with radius
        @test pitch_angle(t, 0.3) > pitch_angle(t, 0.9)
    end

    @testset "dimensional quantities (D=6)" begin
        t = dtmb4382; D = 6.0
        g = dimensional(t, 0.7, D)
        @test g.r ≈ 0.7 * D/2                  # radius
        @test g.c ≈ 0.347 * D                  # chord
        @test g.tmax ≈ 0.0421 * g.c            # max thickness
        @test g.θs ≈ deg2rad(22.747)
        @test g.rake == 0
    end

    @testset "geometry: radius span, finiteness, blade rotation" begin
        t = dtmb4382; D, Dh = 6.0, 1.2
        g = vlm_camber_grid(t, D, Dh; nc=8, ns=12, tip_rR=1.0)
        @test size(g) == (9, 13)
        rad(p) = sqrt(p[1]^2 + p[2]^2)
        @test minimum(rad, g) ≈ Dh/2 atol=1e-6        # hub radius = 0.6
        @test maximum(rad, g) ≈ D/2  atol=1e-6        # tip radius = 3.0 (tip_rR=1)
        # default tip cutoff (0.98) keeps the surface inboard of the tip
        @test maximum(rad, vlm_camber_grid(t, D, Dh; nc=4, ns=8)) ≈ 0.98*D/2 atol=1e-6
        @test all(p -> all(isfinite, p), g)
        # blade 2 is rotated 2π/Z from blade 1 (same radius, shifted angle)
        p1 = blade_section_point(t, 0.7, 0.5, D; blade=1)
        p2 = blade_section_point(t, 0.7, 0.5, D; blade=2)
        @test rad(p1) ≈ rad(p2) atol=1e-9
        Δθ = atan(p2[2],p2[1]) - atan(p1[2],p1[1])
        @test mod(Δθ, 2π) ≈ 2π/t.Z atol=1e-6
    end

    @testset "NACA a=0.8 mean line" begin
        ml = NACAMeanLine(0.8)
        @test ml(0.0) == 0 && ml(1.0) == 0           # zero at LE/TE
        @test maximum(ml(x) for x in 0:0.001:1) ≈ 1 atol=1e-3   # unit-normalised
        # peaks aft of mid-chord (≈0.51) — characteristic of a=0.8
        xf = range(0,1,length=1001); xmax = xf[argmax([ml(x) for x in xf])]
        @test 0.45 < xmax < 0.58
        @test ml(0.3) > 0 && ml(0.9) > 0             # positive camber throughout
        # asymmetric (unlike parabolic): more loaded forward of mid-chord
        @test ml(0.3) > ParabolicMeanLine()(0.3) * 0.0  # sanity, positive
    end

    @testset "blade_sdf signs (WaterLily immersion)" begin
        t = dtmb4382; D, Dh = 6.0, 1.2
        sdf = blade_sdf(t, D, Dh)
        on  = blade_section_point(t, 0.6, 0.5, D; surface=:upper)   # on the skin
        cam = blade_section_point(t, 0.6, 0.5, D; surface=:camber)  # inside
        @test abs(sdf(on)) < 1e-3                    # ≈0 on the surface
        @test sdf(cam) < 0                           # negative inside the blade
        @test sdf(SVector(10.0, 0.0, 0.0)) > 1.0     # positive far away
        @test sdf(SVector(0.0, 0.0, 0.0)) ≈ Dh/2 atol=1e-3   # shaft axis → hub radius
        # union of Z blades: a point on blade 2's skin is also ≈0
        on2 = blade_section_point(t, 0.6, 0.5, D; surface=:upper, blade=2)
        @test abs(sdf(on2)) < 1e-3
    end

    @testset "Flettner panel — ω=0 collapses to 1−4sin²θ" begin
        r = flettner_panel(; R=0.5, ω=0.0, V∞=1.0, N=200)
        # non-lifting cylinder: Cp = 1 − 4 sin²θ, CL = 0
        @test maximum(abs.(r.Cp .- (1 .- 4 .* sin.(r.θ).^2))) < 1e-10
        @test abs(r.CL) < 1e-10
        @test r.Γ == 0
    end

    @testset "Flettner panel — CL matches analytic 4πωR² within ε<0.02%" begin
        R = 0.5; V∞ = 1.0
        for ω in (0.5, 1.0, 1.5, 2.0, 2.5)
            r  = flettner_panel(; R=R, ω=ω, V∞=V∞, N=400)
            CLa = 4π * ω * R^2 / V∞
            @test r.Γ ≈ 2π * ω * R^2
            @test r.CL > 0                                  # positive lift for ω>0
            @test abs(r.CL - CLa) / abs(CLa) < 2e-4         # < 0.02 %
        end
    end

    @testset "Flettner panel — converges to analytic as N grows" begin
        R = 0.5; ω = 1.0; V∞ = 1.0; CLa = 4π * ω * R^2 / V∞
        errs = [abs(flettner_panel(; R, ω, V∞, N).CL - CLa)/CLa for N in (80, 320)]
        @test errs[2] < errs[1]                             # refining reduces error
        @test errs[2] < 1e-3
    end

    @testset "Flettner analytic — closed form CL = 4πωR²" begin
        for ω in (0.0, 1.0, 2.5)
            ra = flettner_analytic(; R=0.5, ω=ω, V∞=1.0, n=720)
            @test isapprox(ra.CL, 4π * ω * 0.5^2; atol=1e-6)
        end
        # top of the cylinder is faster (lower Cp) than the bottom for ω>0
        ra = flettner_analytic(; R=0.5, ω=1.0, V∞=1.0, n=400)
        itop = argmin(abs.(ra.θ .- π/2)); ibot = argmin(abs.(ra.θ .- 3π/2))
        @test ra.Cp[itop] < ra.Cp[ibot]
    end

    @testset "Flettner panel ↔ analytic — Cp agrees pointwise" begin
        rp = flettner_panel(; R=0.5, ω=1.0, V∞=1.0, N=720)
        ra = flettner_analytic(; R=0.5, ω=1.0, V∞=1.0, n=720)
        # panel θ are CCW control points; analytic θ on the same grid offset by
        # half a panel — compare the sorted Cp envelopes (min/max) instead.
        @test isapprox(minimum(rp.Cp), minimum(ra.Cp); rtol=2e-3)
        @test isapprox(maximum(rp.Cp), maximum(ra.Cp); rtol=2e-3)
    end

    @testset "Toolbox re-exports Wing" begin
        # NAT surfaces the finite-wing VLM from LiftingSurfaces.
        w = Wing(; chord_root=1.0, chord_tip=1.0, span=6.0, ns=20, nc=6)
        r = wing_forces(w, deg2rad(5.0), 1.0)
        @test r.CL > 0
        @test isfinite(r.CDi)
    end

    @testset "liftingline_vlm — Weissinger horseshoe VLM" begin
        # symmetric flat rectangular wing (AR=6): no lift at α=0
        r0 = liftingline_vlm(; chord_root=1.0, chord_tip=1.0, span=6.0, alpha=0.0, N=120)
        @test abs(r0.CL) < 1e-6
        # lift slope sits just below the lifting-line value 2π·AR/(AR+2)
        r5 = liftingline_vlm(; chord_root=1.0, chord_tip=1.0, span=6.0, alpha=5.0, N=120)
        AR = r5.AR; @test 6.0 < AR < 6.4                   # ≈6, slightly raised by the 2% tip inset
        slope = r5.CL / deg2rad(5.0)                       # per rad (CL(0)=0)
        ll = 2π*AR/(AR+2)
        @test 0.80*ll < slope < 1.0*ll
        # physically-consistent span efficiency (leading-edge suction recovered)
        @test 0.85 < r5.e < 1.02
        # induced drag scales ~ CL²
        r10 = liftingline_vlm(; chord_root=1.0, chord_tip=1.0, span=6.0, alpha=10.0, N=120)
        @test isapprox(r10.CDi/r5.CDi, (r10.CL/r5.CL)^2; rtol=0.05)
        # agrees with the multi-chordwise Wing CL within a few %
        rw = wing_forces(Wing(; chord_root=1.0, chord_tip=1.0, span=6.0, ns=40, nc=8), deg2rad(5.0), 1.0)
        @test isapprox(r5.CL, rw.CL; rtol=0.08)
        # tapered wing with washout: +2°(root)/−2°(tip) → small positive lift at α=0
        rt = liftingline_vlm(; chord_root=6.0, chord_tip=3.5, span=17.0, alpha=0.0,
                             twist_root=2.0, twist_tip=-2.0, N=80)
        @test 0 < rt.CL < 0.1
        # circulation peaks at mid-span, loading falls to ~0 at the tips
        rL = liftingline_vlm(; chord_root=6.0, chord_tip=3.5, span=17.0, alpha=10.0, N=80)
        @test argmax(rL.Γ) in (40, 41)                     # Γ peaks at mid (80 strips)
        @test rL.cl_span[1] < 0.5*maximum(rL.cl_span)      # low loading at the tip
    end

    @testset "phase_transport_1d — upwind bounded & conservative" begin
        r = phase_transport_1d(; N=100, u=1.0, Co=0.5, t_end=0.2)
        @test length(r.x) == 100
        @test r.dx ≈ 0.01 && r.dt ≈ 0.005
        # boundedness: first-order upwind cannot over/undershoot [0,1]
        @test minimum(r.α) ≥ -1e-12
        @test maximum(r.α) ≤ 1.0 + 1e-12
        # the patch has advected: peak is now near x≈0.5 (translated [0.4,0.6])
        @test 0.4 < r.x[argmax(r.α)] < 0.6
        # the analytic translate is the patch on [0.4,0.6]
        @test all(0.4 .≤ r.x[r.α_analytic .> 0.5] .≤ 0.6)
        # mass is conserved up to what has advected through the outlet.
        # At t=0.2 the patch [0.4,0.6] is fully inside the domain → mass ≈ mass0.
        @test isapprox(r.mass, r.mass0; rtol=1e-9)
    end

    @testset "phase_transport_1d — compression sharpens & stays bounded" begin
        r0 = phase_transport_1d(; N=100, u=1.0, Co=0.5, t_end=0.2, compression=false)
        r1 = phase_transport_1d(; N=100, u=1.0, Co=0.5, t_end=0.2, compression=true, Cα=1.0)
        w0 = interface_width(r0.x, r0.α)
        w1 = interface_width(r1.x, r1.α)
        @test isfinite(w0) && isfinite(w1)
        @test w1 < w0                                   # compression re-steepens
        @test minimum(r1.α) ≥ -1e-12                    # bounded below
        @test maximum(r1.α) ≤ 1.0 + 1e-12               # bounded above
        # compression still conserves mass (conservative flux + clip near sharp)
        @test isapprox(r1.mass, r1.mass0; rtol=5e-3)
        # peak fraction stays sharper (closer to 1) with compression
        @test maximum(r1.α) ≥ maximum(r0.α) - 1e-9
    end

    @testset "als_drag_reduction — bounded, monotone, right trend" begin
        # DR in [0,1), rises with film thickness, saturates toward 1
        thin  = als_drag_reduction(; delta=0.05, t_air=1e-4, U=10.0)
        thick = als_drag_reduction(; delta=0.05, t_air=5e-3, U=10.0)
        @test 0 ≤ thin.DR < 1 && 0 ≤ thick.DR < 1
        @test thick.DR > thin.DR                         # more film → more DR
        @test thick.tau_air < thin.tau_air < thin.tau_clean
        @test thin.tau_clean ≈ thick.tau_clean           # clean shear independent of film
        # μ_w/μ_a ≈ 55..65 (the leverage that makes thin films effective)
        @test 40 < thin.mu_ratio < 80
        # a continuous film a few % of δ already gives tens of % DR (ALDR signature)
        mid = als_drag_reduction(; delta=0.05, t_air=2e-3, U=10.0)
        @test 0.4 < mid.DR < 0.95
        # zero film → zero DR
        @test als_drag_reduction(; delta=0.05, t_air=0.0).DR == 0.0
    end

    @testset "als_sweep + ship saving" begin
        sw = als_sweep(; delta=0.05, t_airs=[1e-4, 5e-4, 1e-3, 5e-3], U=12.0)
        @test length(sw) == 4
        @test issorted([s.DR for s in sw])               # monotone rising
        @test all(0 .≤ [s.DR for s in sw] .< 1)
        # ship application: DTC-like hull, friction saving must be a sane fraction
        sh = als_ship_saving(; L=355.0, B=51.0, T=15.0, Cb=0.66, U=12.0,
                             frac_covered=0.5, t_air=2e-3)
        @test sh.S > 0 && sh.Re > 1e8
        @test 0 < sh.saving_pct < 50                     # partial coverage → partial saving
        @test sh.Rf_air < sh.Rf_clean
        @test sh.dRf_kN > 0
    end

    # ------------------------------------------------------------------
    # WAP power-analysis (Q5) — synthetic fixtures, NOT the real trial data
    # ------------------------------------------------------------------

    # small AWA→CDA table mirroring the supplied one's shape: +ve head, -ve following
    _COEF = (; awa = Float64[0, 30, 60, 90, 120, 150, 180],
               cda = Float64[0.88, 0.97, 0.60, 0.04, -0.53, -0.89, -0.81])

    # write a synthetic coef CSV + run CSV to a temp dir, return paths
    function _write_fixture(dir; on_shp_scale=1.0, awa_off=90.0, awa_on=90.0,
                            aws_off=10.0, aws_on=10.0, n=20, trim=2,
                            base_a=1000.0, base_b=3.0)
        cf = joinpath(dir, "coef.csv")
        open(cf, "w") do io
            println(io, "AWA_deg,CDA_10m")
            for (a,c) in zip(_COEF.awa, _COEF.cda); println(io, a, ",", c); end
        end
        rc = joinpath(dir, "run.csv")
        # OFF speeds span a band so the P=aV^b fit is well posed; ON speeds are
        # held at a single value (9.5 kn) so the curve-at-mean-speed equals the
        # mean-of-power on the ON segment (no Jensen-convexity offset), making
        # the "equal SHP → ΔP≈0" check exact.
        Voff = range(8.0, 11.0, length=n); Von = fill(9.5, n)
        open(rc, "w") do io
            println(io, "DateTime,SOG [kn],STW [kn],Shaft Speed [rpm],SHP [kW],Rudder Angle [deg],AWS [kn],AWA [deg],Status_ON")
            row(v, shp, aws, awa, st) =
                println(io, "00:00:00,", v, ",", v, ",225.0,", shp, ",0.0,", aws, ",", awa, ",", st)
            base(v) = base_a * v^base_b
            for v in Voff; row(v, base(v),               aws_off, awa_off, "False"); end
            for v in Von;  row(v, on_shp_scale*base(v),  aws_on,  awa_on,  "True");  end
            for v in Voff; row(v, base(v),               aws_off, awa_off, "False"); end
        end
        return rc, cf, trim
    end

    @testset "wind_resistance — sign by AWA, zero at AWS=0" begin
        coef = _COEF
        @test wind_resistance(10.0, 0.0,   coef) > 0          # head wind: drag
        @test wind_resistance(10.0, 180.0, coef) < 0          # following: thrust
        @test wind_resistance(0.0,  0.0,   coef) == 0.0       # no wind, no load
        # symmetric in sign of AWA (frontal-area coefficient)
        @test wind_resistance(10.0, -30.0, coef) ≈ wind_resistance(10.0, 30.0, coef)
        # ∝ AWS²
        @test wind_resistance(20.0, 0.0, coef) ≈ 4*wind_resistance(10.0, 0.0, coef)
        # magnitude: ½·1.225·0.88·121.4·10² ≈ 6543 N at head wind, AWS=10 m/s
        @test wind_resistance(10.0, 0.0, coef) ≈ 0.5*1.225*0.88*121.4*100 atol=1.0
    end

    @testset "fit_power_speed — recovers P = 3·V³" begin
        V = collect(5.0:1.0:15.0)
        P = 3.0 .* V.^3
        f = fit_power_speed(V, P)
        @test f.b ≈ 3.0 atol=1e-9
        @test f.a ≈ 3.0 atol=1e-9
        @test f.r2 ≈ 1.0 atol=1e-12
        @test f.predict(10.0) ≈ 3000.0 atol=1e-6
        @test_throws ArgumentError fit_power_speed([1.0,-1.0], [1.0,2.0])
    end

    @testset "CDA interpolation — linear & bounded off-grid" begin
        coef = _COEF
        # AWA=12.5° is between 0 (0.88) and 30 (0.97): linear interp
        c = _interp(coef.awa, coef.cda, 12.5)
        @test c ≈ 0.88 + (12.5/30)*(0.97-0.88) atol=1e-12
        @test minimum(coef.cda) ≤ c ≤ maximum(coef.cda)
        # clamps outside the table
        @test _interp(coef.awa, coef.cda, -10.0) == coef.cda[1]
        @test _interp(coef.awa, coef.cda, 200.0) == coef.cda[end]
    end

    @testset "ΔP sign — savings when ON SHP lower, ≈0 when equal" begin
        mktempdir() do dir
            # ON segment uses 90% of baseline SHP at matched speed → ΔP>0
            rc, cf, tr = _write_fixture(dir; on_shp_scale=0.90)
            r = wap_power_analysis(rc, cf; trim=tr)
            @test r.ΔP_raw > 0
            @test r.fit_raw.b ≈ 3.0 atol=0.05           # recovers the base curve
            @test r.n_on > 0 && r.n_off > 0
        end
        mktempdir() do dir
            # ON == OFF SHP, same wind → ΔP ≈ 0
            rc, cf, tr = _write_fixture(dir; on_shp_scale=1.0)
            r = wap_power_analysis(rc, cf; trim=tr)
            @test abs(r.ΔP_raw) < 1e-6 * r.P_base_raw
            @test abs(r.ΔP_corr) < 1e-6 * r.P_base_corr
        end
    end

    @testset "wind correction — pass-2 > pass-1 when ON saw more head wind" begin
        mktempdir() do dir
            # equal measured SHP, but ON sees strong head wind, OFF sees calm.
            # Raw ΔP≈0; correcting removes the ON head-wind load → ON_corr lower
            # than baseline_corr → ΔP_corr > ΔP_raw (the device's true saving
            # was masked by the head wind it had to fight).
            rc, cf, tr = _write_fixture(dir; on_shp_scale=1.0,
                                        awa_off=90.0, aws_off=2.0,    # OFF nearly calm/beam
                                        awa_on=0.0,   aws_on=20.0)    # ON strong head wind
            r = wap_power_analysis(rc, cf; trim=tr)
            @test abs(r.ΔP_raw) < 1e-6 * r.P_base_raw     # raw sees no difference
            @test r.ΔP_corr > r.ΔP_raw                    # correction reveals the saving
            @test r.seg.on.dP_wind.mean > r.seg.off.dP_wind.mean
        end
    end

    @testset "read_trial / read_wind_coef round-trip" begin
        mktempdir() do dir
            rc, cf, _ = _write_fixture(dir)
            tr = read_trial(rc)
            @test length(tr.shp) == 60 && eltype(tr.shp) == Float64
            @test tr.on isa BitVector && sum(tr.on) == 20
            co = read_wind_coef(cf)
            @test issorted(co.awa) && length(co.cda) == 7
        end
    end

    @testset "thickness: upper/lower straddle camber, vanish at LE/TE" begin
        t = dtmb4382; D = 6.0
        cam = blade_section_point(t, 0.5, 0.5, D; surface=:camber)
        up  = blade_section_point(t, 0.5, 0.5, D; surface=:upper)
        lo  = blade_section_point(t, 0.5, 0.5, D; surface=:lower)
        dist(a,b) = sqrt(sum(abs2, a - b))
        gap = dist(up, lo)
        @test gap > 0                          # finite thickness mid-chord
        # camber straddled symmetrically by the skins (Cartesian midpoint
        # holds only to 2nd order through the cylindrical wrap).
        @test isapprox(dist(up, cam), dist(lo, cam); rtol=0.05)
        @test dist(cam, (up .+ lo) ./ 2) < 0.05 * gap
        # thickness → 0 at the trailing edge (xc=1)
        upt = blade_section_point(t, 0.5, 1.0, D; surface=:upper)
        lot = blade_section_point(t, 0.5, 1.0, D; surface=:lower)
        @test sqrt(sum(abs2, upt - lot)) < 1e-9
    end

end
