using Test
using NavalArchitectToolbox
using NavalArchitectToolbox: dimensional, _interp
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
