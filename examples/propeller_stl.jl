#!/usr/bin/env julia
#
# Watertight STL of the DTMB 4381 propeller (5 blades + hub) for external
# meshers (OpenFOAM snappyHexMesh). Marching cubes over the analytic
# blade_sdf ∪ hub SDF — a sampled isosurface is closed by construction,
# which is far more robust than stitching parametric skins at LE/TE/tip.
#
# The SAME analytic SDF drives the WaterLily immersed body, so both codes
# see one geometry definition.
#
# Usage: julia examples/propeller_stl.jl [outfile.stl]

import Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop(path=joinpath(@__DIR__, ".."), io=devnull)
Pkg.add(["Meshing", "StaticArrays"]; io=devnull)

using NavalArchitectToolbox, Meshing, StaticArrays, Printf

# ── geometry: DTMB 4381, model scale (1 ft diameter, standard DTMB model)
const D  = 0.3048                 # propeller diameter [m]
const Dh = 0.2D                   # hub diameter = r/R 0.2 table root
const HUB_HALFLEN = 0.16D         # hub half-length [m] (compact, capped)

sdf_blades = blade_sdf(dtmb4381, D, Dh)

# exact SDF of a z-aligned capped cylinder (radius rh, half-length hl)
@inline function sdf_hub(x)
    q1 = hypot(x[1], x[2]) - Dh/2
    q2 = abs(x[3]) - HUB_HALFLEN
    return min(max(q1, q2), 0.0) + hypot(max(q1, 0.0), max(q2, 0.0))
end

sdf(x) = min(sdf_blades(x), sdf_hub(x))

# ── sample onto a grid. Blade axial extent ≈ ±(c·sinβ + rake) ≲ 0.2D, so
# z ∈ ±0.3D covers blades + hub caps with margin.
const NX = 384                    # transverse resolution (x, y)
const LO = SVector(-0.60D, -0.60D, -0.30D)
const HI = SVector( 0.60D,  0.60D,  0.30D)
NZ = round(Int, NX * (HI[3]-LO[3]) / (HI[1]-LO[1]))
dims = (NX, NX, NZ)
println("sampling SDF on $(dims[1])×$(dims[2])×$(dims[3]) grid ",
        "(h = $(round(1e3*(HI[1]-LO[1])/(NX-1); digits=4)) mm)...")

A = Array{Float32}(undef, dims)
xs = range(LO[1], HI[1]; length=dims[1])
ys = range(LO[2], HI[2]; length=dims[2])
zs = range(LO[3], HI[3]; length=dims[3])
Threads.@threads for k in 1:dims[3]
    for j in 1:dims[2], i in 1:dims[1]
        A[i,j,k] = Float32(sdf(SVector(xs[i], ys[j], zs[k])))
    end
end
lo_hi = extrema(A)
println("SDF ∈ [$(lo_hi[1]), $(lo_hi[2])]")

vts, fcs = isosurface(A, MarchingCubes(), xs, ys, zs)
verts = [SVector{3,Float64}(v...) for v in vts]
@printf("marching cubes: %d vertices, %d triangles\n", length(verts), length(fcs))

# ── weld: adjacent voxels reproduce shared edge-vertices with ~1e-9 float
# scatter, which external tools see as unconnected points. Snap to a 1e-6 m
# grid, remap faces, drop degenerate triangles.
snap(v) = SVector(round(v[1]; digits=6), round(v[2]; digits=6), round(v[3]; digits=6))
vmap = Dict{SVector{3,Float64},Int}()
newid = zeros(Int, length(verts))
wverts = SVector{3,Float64}[]
for (i, v) in pairs(verts)
    sv = snap(v)
    newid[i] = get!(vmap, sv) do
        push!(wverts, sv); length(wverts)
    end
end
wfcs = NTuple{3,Int}[]
for f in fcs
    a, b, c = newid[f[1]], newid[f[2]], newid[f[3]]
    (a == b || b == c || a == c) && continue          # degenerate after weld
    push!(wfcs, (a, b, c))
end
verts, fcs = wverts, wfcs
@printf("after weld: %d vertices, %d triangles\n", length(verts), length(fcs))

# ── keep the dominant connected component: the razor-thin tip/LE regions
# leave isolated marching-cubes specks (few-triangle shells) that break
# closedness checks in external meshers. Union-find over shared vertices.
parent = collect(1:length(verts))
findr(i) = (while parent[i] != i; parent[i] = parent[parent[i]]; i = parent[i]; end; i)
for f in fcs
    a = findr(f[1]); b = findr(f[2]); c = findr(f[3])
    parent[b] = a; parent[c] = a
end
compsize = Dict{Int,Int}()
for f in fcs
    compsize[findr(f[1])] = get(compsize, findr(f[1]), 0) + 1
end
main = argmax(compsize)
fcs = [f for f in fcs if findr(f[1]) == main]
@printf("after component filter: %d triangles (dropped %d fragments, %d specks)\n",
        length(fcs), length(compsize) - 1, sum(values(compsize)) - compsize[main])

# ── binary STL writer (little-endian; normal recomputed per facet)
function write_stl(path, verts, faces)
    open(path, "w") do io
        write(io, zeros(UInt8, 80))               # header
        write(io, UInt32(length(faces)))
        for f in faces
            a, b, c = verts[f[1]], verts[f[2]], verts[f[3]]
            n = cross(b .- a, c .- a); nn = norm(n)
            n = nn > 0 ? n ./ nn : n
            for v in (n, a, b, c), s in 1:3
                write(io, Float32(v[s]))
            end
            write(io, UInt16(0))                  # attribute byte count
        end
    end
end
using LinearAlgebra: cross, norm

out = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "runs", "dtmb4381_Z5.stl")
mkpath(dirname(out))
write_stl(out, verts, fcs)
@printf("wrote %s (%.1f MB)\n", out, filesize(out)/1e6)
