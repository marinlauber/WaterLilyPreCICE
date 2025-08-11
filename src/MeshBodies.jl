using FileIO, MeshIO
using GeometryBasics
using WaterLily: AbstractBody,@loop,measure
using StaticArrays
using ForwardDiff

mutable struct MeshBody{T,F<:Function} <: AbstractBody
    mesh  :: GeometryBasics.Mesh
    mesh0 :: GeometryBasics.Mesh
    velocity :: GeometryBasics.Mesh
    surf_id :: NTuple
    map   :: F
    bbox  :: Rect
    scale :: T
    half_thk::T #half thickness
    boundary::Bool
end
function MeshBody(fname::String;map=(x,t)->x,scale=1.0,thk=0f0,boundary=true,T=Float32)
    if endswith(fname,".inp")
        mesh0,srf_id = load_inp(fname)
    else
        mesh0 = load(fname)
        srf_id = ntuple(i->(1,i),length(mesh0))
    end
    return MeshBody(mesh0,srf_id,map,scale,thk,boundary,T)
end
MeshBody(mesh::M;map=(x,t)->x,scale=1.0,thk=0f0,boundary=true,T=Float32) where M<:GeometryBasics.Mesh = 
        MeshBody(mesh,ntuple(i->(1,i),length(mesh)),map,scale,thk,boundary,T)
function MeshBody(mesh0::M,srf_id,map,scale,thk,boundary,T) where M<:GeometryBasics.Mesh
    points = Point3f[] # scale and map the points to the corrcet location
    for (i,pnt) in enumerate(mesh0.position)
        push!(points, Point3f(map(SA[pnt.data...]*scale,0)))
    end
    mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(mesh0))
    velocity = GeometryBasics.Mesh(zero(points),GeometryBasics.faces(mesh0))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    return MeshBody(mesh,mesh0,velocity,srf_id,map,bbox,T(scale),T(thk/2),boundary)
end
Base.copy(a::GeometryBasics.Mesh) = GeometryBasics.Mesh(a.position,GeometryBasics.faces(a));
Base.copy(b::MeshBody) = MeshBody(copy(b.mesh),copy(b.mesh0),copy(b.velocity),b.surf_id,b.map,
                                  Rect(b.bbox),b.scale,b.half_thk,b.boundary)

function load_inp(fname; facetype=GLTriangleFace, pointtype=Point3f)
    #INP file format
    @assert endswith(fname,".inp") "file type not supported"
    fs = open(fname)

    points = pointtype[]
    faces = facetype[]
    node_idx = Int[]
    srf_id = Tuple[]
    cnt = 0

    # read the first 3 lines if there is the "*heading" keyword
    line = readline(fs)
    contains(line,"*heading") && (line = readline(fs))
    BlockType = contains(line,"*NODE") ? Val{:NodeBlock}() : Val{:DataBlock}()
    
    # read the file
    while !eof(fs)
        line = readline(fs)
        contains(line,"*ELSET, ELSET=") && (cnt+=1)
        BlockType, line = parse_blocktype!(BlockType, fs, line)
        if BlockType == Val{:NodeBlock}()
            push!(node_idx, parse(Int,split(line,",")[1])) # keep track of the node index of the inp file
            push!(points, pointtype(parse.(eltype(pointtype),split(line,",")[2:4])))
        elseif BlockType == Val{:ElementBlock}()
            nodes = parse.(Int,split(line,",")[2:end])
            push!(faces, TriangleFace{Int}(facetype([findfirst(==(node),node_idx) for node in nodes])...)) # parse the face
        elseif BlockType == Val{:ElSetBlock}()
            push!(srf_id, (cnt, parse.(Int,split(line,",")[1])));
        else
            continue
        end
    end
    close(fs) # close file stream
    # case where there is no surface ID and we pass it a single tuple of all the faces
    srf_id = length(srf_id)==0 ? ntuple(i->(1,i),length(faces)) : ntuple(i->(srf_id[i].+(0,1-srf_id[1][2])),length(srf_id))
    return Mesh(points, faces), srf_id
end
function parse_blocktype!(block, io, line)
    contains(line,"*NODE") && return block=Val{:NodeBlock}(),readline(io)
    contains(line,"*ELEMENT") && return block=Val{:ElementBlock}(),readline(io)
    contains(line,"*ELSET, ELSET=") && return block=Val{:ElSetBlock}(),readline(io)
    return block, line
end

"""
    locate(p,tri)

Find the closest point `x` on the triangle `tri` to the point `p`.
"""
function locate(tri::GeometryBasics.Ngon{3},p::SVector{T}) where T #5.327 ns (0 allocations: 0 bytes)
    # is point `a` closest?
    a, b, c = tri.points
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

# inside bbox or not
outside(x::SVector,bbox::Rect) = !(all(bbox.origin .≤ x) && all(x .≤ bbox.origin+bbox.widths)) # 1.679 ns (0 allocations: 0 bytes)
# distance to box center
dist(x::SVector,bbox::Rect) = √sum(abs2,x.-bbox.origin-0.5bbox.widths) # 1.707 ns (0 allocations: 0 bytes)

WaterLily.sdf(body::MeshBody,x,t;kwargs...) = measure(body,x,t;kwargs...)[1]
"""
    measure(body::MeshBody,x,t;kwargs...)

    Measures a mesh body at point `x` and time `t`.
    If a mapping has been specified, the point `x` is first map to the new location `ξ`, and the mesh is measured there.
    We use a large bounding box around the mesh to avoid measuring  far away, where the actual distance doesn't matter as
    long as d>>O(1), this means that outside the bounding box, we return the distance of the bbox edges to the geom (d=8).
    If the mesh is not a boundary, we need to adjust the distance by the half thickness of the body.
"""
function WaterLily.measure(body::MeshBody,x::SVector{D},t;kwargs...) where D
    # if x is 2D, we need to make it 3D
    D==2 && (x=SVector{3}(x[1],x[2],0))
    # if we are outside of the bounding box, we don't even to measure
    outside(x,body.bbox) && return (max(8,2body.half_thk),zero(x),zero(x)) # we don't need to worry if the geom is a boundary or not
    d,n,v = measure(body.mesh,body.velocity,x;kwargs...)
    !body.boundary && (d = abs(d)-body.half_thk) # if the mesh is not a boundary, we need to adjust the distance
    D==2 && (n=SVector{D}(n[1],n[2]); n=n/√sum(abs2,n))
    return (d,n,v)
end
function update!(body::MeshBody,t,dt=0;kwargs...)
    # update mesh position, measure is done after
    points = Point3f[]
    for (i,pnt) in enumerate(body.mesh0.position)
        push!(points, Point3f(body.map(SA[pnt.data...]*body.scale,t)))
    end
    dt>0 && (body.velocity = GeometryBasics.Mesh((points-body.mesh.position)/dt,GeometryBasics.faces(body.mesh0)))
    body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(body.mesh0))
    # update bounding box
    bbox = Rect(points)
    body.bbox = Rect(bbox.origin.-max(4,2body.half_thk),bbox.widths.+max(8,4body.half_thk))
end

using LinearAlgebra: cross
"""
    normal(tri::GeometryBasics.Ngon{3})

Return the normal vector to the triangle `tri`.
"""
function normal(tri::GeometryBasics.Ngon{3}) #3.039 ns (0 allocations: 0 bytes)
    n = cross(SVector(tri.points[2]-tri.points[1]),SVector(tri.points[3]-tri.points[1]))
    n/(√sum(abs2,n)+eps(eltype(n))) # zero area triangles
end
@inbounds @inline d²_fast(tri::GeometryBasics.Ngon{3},x) = sum(abs2,x-center(tri)) #1.825 ns (0 allocations: 0 bytes)
@inbounds @inline d²(tri::GeometryBasics.Ngon{3},x) = sum(abs2,x-locate(tri,x)) #4.425 ns (0 allocations: 0 bytes)
@inbounds @inline center(tri::GeometryBasics.Ngon{3}) = SVector(sum(tri.points)/3.f0...) #1.696 ns (0 allocations: 0 bytes)
@inbounds @inline area(tri::GeometryBasics.Ngon{3}) = 0.5*√sum(abs2,cross(SVector(tri.points[2]-tri.points[1]),SVector(tri.points[3]-tri.points[1]))) #1.784 ns (0 allocations: 0 bytes)
@inbounds @inline dS(tri::GeometryBasics.Ngon{3}) = 0.5cross(SVector(tri.points[2]-tri.points[1]),SVector(tri.points[3]-tri.points[1]))
"""
    measure(mesh::GeometryBasics.Mesh,x,t;kwargs...)

Measure the distance `d` and normal `n` to the mesh at point `x` and time `t`.
"""
function WaterLily.measure(mesh::M,velocity::M,x::SVector{T};kwargs...) where {M<:GeometryBasics.Mesh,T}
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    n,p = normal(mesh[u]),SVector(locate(mesh[u],x))
    s = x-p # signed Euclidian distance
    d = sign(sum(s.*n))*√sum(abs2,s)
    # v = get_velocity(mesh[u],velocity[u],p)
    return d,n,zero(x) #v
end # 120.029 ns (0 allocations: 0 bytes)d # 4.266 μs (0 allocations: 0 bytes)

# linear shape function interpolation of the nodal velocity values at point `p`
function get_velocity(tri::Tr,vel::Tr,p::SVector{3,T}) where {Tr<:GeometryBasics.Ngon{3},T}
    dA = SVector{3,T}([sub_area(tri,p,Val{i}()) for i in 1:3])
    return SVector(sum(tri.points.*dA)/sum(dA))
end
@inline sub_area(t,p,::Val{1}) = area(GeometryBasics.Ngon(SA[Point(p),t[2],t[3]]))
@inline sub_area(t,p,::Val{2}) = area(GeometryBasics.Ngon(SA[t[1],Point(p),t[3]]))
@inline sub_area(t,p,::Val{3}) = area(GeometryBasics.Ngon(SA[t[1],t[2],Point(p)]))

# use divergence theorem to calculate volume of surface mesh
# F⃗⋅k⃗ = -⨕pn⃗⋅k⃗ dS = ∮(C-ρgz)n⃗⋅k⃗ dS = ∫∇⋅(C-ρgzk⃗)dV = ρg∫∂/∂z(ρgzk⃗)dV = ρg∫dV = ρgV #
volume(a::GeometryBasics.Mesh) = mapreduce(T->center(T).*dS(T),+,a)
volume(body::MeshBody) = volume(body.mesh)

import WaterLily: interp
function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,3},δ,::Val{true}) where T
    c=center(tri);n=normal(tri);ar=area(tri);
    ar.*n.*interp(c.+1.5 .+ δ.*n, p)
end
function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,3},δ,::Val{false}) where T
    c=center(tri);n=normal(tri);ar=area(tri);
    ar.*n.*(interp(c.+1.5 .+ δ.*n, p) .- interp(c.+1.5 .- δ.*n, p))
end
function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,2},δ,::Val{true}) where T
    c=center(tri)[1:2];n=normal(tri)[1:2];ar=area(tri);
    ar.*n.*interp(c.+1.5 .+ δ.*n, p)
end
forces(a::GeometryBasics.Mesh, flow::Flow, δ=2, boundary=Val{true}()) = map(T->get_p(T, flow.p, δ, boundary), a)
forces(body::MeshBody, b::Flow, δ=2) = forces(body.mesh, b, δ, Val{body.boundary}())

using Printf: @sprintf
import WaterLily
using WriteVTK

# access the WaterLily writer to save the file
function WaterLily.save!(w,a::MeshBody,t=w.count[1]) #where S<:AbstractSimulation{A,B,C,D,MeshBody}
    k = w.count[1]
    points = hcat([[p.data...] for p ∈ a.mesh.position]...)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, Base.to_index.(face)) for face in faces(a.mesh)]
    vtk = vtk_grid(w.dir_name*@sprintf("/%s_%06i", w.fname, k), points, cells)
    for (name,func) in w.output_attrib
        # point/vector data must be oriented in the same way as the mesh
        vtk[name] = ndims(func(a))==1 ? func(a) : permutedims(func(a))
    end
    vtk_save(vtk); w.count[1]=k+1
    w.collection[round(t,digits=4)]=vtk
end