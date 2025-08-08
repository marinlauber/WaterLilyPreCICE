using WaterLilyPreCICE,StaticArrays
using GeometryBasics
using BenchmarkTools
using CUDA
using Plots
using Revise
using Adapt

import WaterLilyPreCICE: locate,normal,d²_fast,hat,load_inp,load,center

locate(tri::SMatrix{T},p::SVector{T}) where T = locate(tri[:,1],tri[:,2],tri[:,3],p)
function locate(a,b,c,p::SVector{T}) where T
    ab = b.-a
    ac = c.-a
    ap = p.-a
    d1 = sum(ab.*ap)
    d2 = sum(ac.*ap)
    # is point `a` closest?
    if ((d1 ≤ 0) && (d2 ≤ 0))
        return a
    end
    # is point `b` closest?
    bp = p.-b
    d3 = sum(ab.*bp)
    d4 = sum(ac.*bp)
    if ((d3 ≥ 0) && (d4 ≤ d3))
        return b
    end
    # is point `c` closest?
    cp = p.-c
    d5 = sum(ab.*cp)
    d6 = sum(ac.*cp)
    if ((d6 ≥ 0) && (d5 ≤ d6))
        return c
    end
    # is segment 'ab' closest?
    vc = d1*d4 - d3*d2
    if ((vc ≤ 0) && (d1 ≥ 0) && (d3 ≤ 0))
        x =  a .+ ab.*d1 ./ (d1 - d3)
        return x
    end
    #  is segment 'ac' closest?
    vb = d5*d2 - d1*d6
    if ((vb ≤ 0) && (d2 ≥ 0) && (d6 ≤ 0))
        x =  a .+ ac.*d2 ./ (d2 - d6)
        return x
    end
    # is segment 'bc' closest?
    va = d3*d6 - d5*d4
    if ((va ≤ 0) && (d4 ≥ d3) && (d5 ≥ d6))
        x =  b .+ (c .- b) .* (d4 - d3) ./ ((d4 - d3) + (d5 - d6))
        return x
    end
    # closest is interior to `abc`
    denom = one(T) / (va + vb + vc)
    v= vb*denom
    w = vc*denom
    x = a .+ ab .* v .+ ac .* w
    return x
end

# closest point and normal there
function closest(mesh,x::SVector{T};kwargs...) where T
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    n,p = normal(mesh[u]),SVector(locate(mesh[u],x))
    return p,n
end

# struct SimpleMesh{T,A<:AbstractArray,S<:AbstractVector{T}} <: AbstractBody
struct SimpleMesh{T,A<:AbstractArray,S<:AbstractVector{T},F<:Function} <: AbstractBody
    mesh :: A
    origin :: S
    width :: S
    boundary :: Bool
    map   :: F
    scale :: T
    half_thk :: T
    function SimpleMesh(mesh,origin::SVector{3,T},width::SVector{3,T},boundary=true,
                        map=(x,t)->x,scale=1,half_thk=0) where T
        new{T,typeof(mesh),typeof(origin),typeof(map)}(mesh,origin,width,boundary,map,T(scale),T(half_thk))
    end
end
Adapt.adapt_structure(to, x::SimpleMesh) = SimpleMesh(adapt(to, x.mesh), adapt(to, x.origin),
                                                      adapt(to, x.width), x.boundary, adapt(to, x.map),
                                                      x.scale, x.half_thk)

function SimpleMesh(file_name; map=(x,t)->x,boundary=true,half_thk=0,scale=1,mem=Array,T=Float32)
    # read in the mesh
    mesh = endswith(file_name,".inp") ? load_inp(file_name)[1] : load(file_name)
    points = Point3f[] # scale and map the points to the correct location
    for pnt in mesh.position
        push!(points, Point3f(SA{T}[pnt.data...]*T(scale)))
    end
    mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(mesh))
    width = SVector{3,T}(maximum(mesh.position,dims=1)...)
    origin = SVector{3,T}(minimum(mesh.position,dims=1)...)
    mesh = [hcat(vcat([mesh[i]...])...) for i in 1:length(mesh)] |> mem
    SimpleMesh(mesh,origin.-2,width-origin.+4,boundary,map,T(scale),T(half_thk))
end

# WaterLily.sdf(body::SimpleMesh,x,t;kwargs...) = distance_gpu(body.mesh,x)
WaterLily.sdf(body::SimpleMesh,x,t;kwargs...) = measure(body,x,t;kwargs...)[1]

using ForwardDiff
function WaterLily.measure(body::SimpleMesh,x::SVector{D,T},t;kwargs...) where {D,T}
    # map to correct location
    ξ = body.map(x,t)
    # we don't need to worry if the geom is a boundary or not
    outside(ξ,body.origin,body.width) && return (T(4),zeros(SVector{D,T}),zeros(SVector{D,T}))
    # locate the point on the mesh
    p,n = closest(body.mesh,ξ;kwargs...)
    # signed Euclidian distance
    s = ξ-p; d = sign(sum(s.*n))*√sum(abs2,s)
    # velocity at the mesh point
    J = ForwardDiff.jacobian(ξ->body.map(ξ,t), ξ)
    dot = ForwardDiff.derivative(t->body.map(ξ,t), t)
    # if the mesh is not a boundary, we need to adjust the distance
    !body.boundary && (d = abs(d)-body.half_thk)
    return (d,n,-J\dot)
end
outside(x::SVector,origin,width) = !(all(origin .≤ x) && all(x .≤ origin+width))

using LinearAlgebra: cross
d²_fast(tri::SMatrix,x::SVector) = sum(abs2,x-SVector(sum(tri,dims=2)/3))
normal(tri::SMatrix) = hat(SVector(cross(tri[:,2]-tri[:,1],tri[:,3]-tri[:,1])))


function make(L;U=1,mem=CuArray,T=Float32)
    α = π/10.f0 # rotation angle
    function map(x,t)
        Rx = SA[1 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]
        Ry = SA[cos(α) 0 sin(α); 0 1 0; -sin(α) 0 cos(α)]
        Rz = SA[cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
        Rx*Ry*Rz*(x.-L/2.f0).+0.5f0
    end

    body = SimpleMesh(joinpath(@__DIR__,"meshes/cube.inp");scale=L/2,map,mem)
    # Simulation((L,L,L),(U,0,0),L;body,mem,T)
end

L = 32
MEMORY = CuArray
body = make(L;mem=MEMORY)

distance = zeros(L+2,L+2,L+2)|> MEMORY;
@inside distance[I] = sdf(body,loc(0,I))
flood(distance[:,L÷2+1,:]); contour!(Array(distance)[:,L÷2+1,:]',levels=[0],lw=2)




# import WaterLilyPreCICE: @loop
# @loop new_body.mesh[I] = SMatrix{3,3,Float32}(rand(3,3)) over I in CartesianIndices(new_body.mesh)