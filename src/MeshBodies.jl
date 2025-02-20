using FileIO, MeshIO
using GeometryBasics
using WaterLily: AbstractBody,@loop,measure
using StaticArrays
using ForwardDiff

mutable struct MeshBody{T,F<:Function} <: AbstractBody
    mesh  :: GeometryBasics.Mesh
    srfID :: Union{Nothing,NTuple}
    map   :: F
    bbox  :: Rect
    scale :: T
    half_thk::T #half thickness
    boundary::Bool
end
function MeshBody(fname::String;scale=1.0,boundary=true,thk=0f0,T=Float32)
    if endswith(fname,".inp")
        tmp,srf_id = load_inp(fname)
    else
        tmp = load(fname)
        srf_id = nothing
    end
    points = GeometryBasics.Point.(tmp.position*scale) # can we specify types?
    mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(tmp))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    return MeshBody(mesh,srf_id,map,bbox,T(scale),T(thk/2),boundary)
end
Base.copy(b::MeshBody) = (mesh=GeometryBasics.Mesh(b.mesh.position,GeometryBasics.faces(b.mesh));
                          MeshBody(mesh,b.srfID,b.map,Rect(b.bbox),b.scale,b.half_thk,b.boundary))

function load_inp(fname; facetype=GLTriangleFace, pointtype=Point3f)
    #INP file format
    @assert endswith(fname,".inp") "file type not supported"
    fs = open(fname)

    points = pointtype[]
    faces = facetype[]
    node_idx = Int[]
    srf_id = Int[]
    tmp = Int[]
    cnt = 1

    # read the first 3 lines if there is the "*heading" keyword
    line = readline(fs)
    contains(line,"*heading") && (line = readline(fs))
    BlockType = contains(line,"*NODE") ? Val{:NodeBlock}() : Val{:DataBlock}()
    
    # read the file
    while !eof(fs)
        line = readline(fs)
        contains(line,"*ELSET, ELSET=") && (push!(tmp,cnt))
        BlockType, line = parse_blocktype!(BlockType, fs, line)
        if BlockType == Val{:NodeBlock}()
            push!(node_idx, parse(Int,split(line,",")[1])) # keep track of the node index of the inp file
            push!(points, pointtype(parse.(eltype(pointtype),split(line,",")[2:4])))
        elseif BlockType == Val{:ElementBlock}()
            nodes = parse.(Int,split(line,",")[2:end])
            push!(faces, TriangleFace{Int}(facetype([findfirst(==(node),node_idx) for node in nodes])...)) # parse the face
        elseif BlockType == Val{:ElSetBlock}()
            push!(srf_id, parse.(Int,split(line,",")[1])); cnt+=1
        else
            continue
        end
    end
    push!(tmp,cnt) # push the last element
    # reshape the surface id vector, the first ID resets the count
    srf_id = ntuple(i->srf_id[tmp[i]:tmp[i+1]-1].-srf_id[1].+1,length(tmp)-1)
    close(fs) # close file stream
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
function WaterLily.measure(body::MeshBody,x::SVector{D},t,;kwargs...) where D
    # eval d=map(x,t)-x | if ξ is 2D, we need to make it 3D
    ξ = body.map(x,t); D==2 && (ξ=SVector{3}(ξ[1],ξ[2],0))
    #  if we are outside of the bounding box, we don't even to measure
    outside(ξ,body.bbox) && return (max(8,2body.half_thk),zero(x),zero(x)) # we don't need to worry if the geom is a boundary or not
    d,n = measure(body.mesh,ξ,t;kwargs...)
    !body.boundary && (d = abs(d)-body.half_thk) # if the mesh is not a boundary, we need to adjust the distance
    D==2 && (n=SVector{D}(n[1],n[2]); n=n/√sum(abs2,n))
    # The velocity depends on the material change of ξ=m(x,t):
    #   Dm/Dt=0 → ṁ + (dm/dx)ẋ = 0 ∴  ẋ =-(dm/dx)\ṁ
    J = ForwardDiff.jacobian(x->body.map(x,t), x)
    dot = ForwardDiff.derivative(t->body.map(x,t), t)
    return (d,n,-J\dot)
end
function WaterLily.update!(body::MeshBody;update::Function=(x)->x)
    @assert((hasmethod(update,Tuple{SVector}) && typeof(update(SVector{3}(zeros(3))))==typeof(SVector{3}(zeros(3)))),
            "`update` must have a method of type ::SVector{3}->::SVector{3}")
    # update mesh position, measure is done elsewhere
    points = Point3f[]
    for (i,pnt) in enumerate(body.mesh.position)
        push!(points, Point3f(update(SA[pnt.data...])))
    end
    body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(body.mesh))
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
function WaterLily.measure(mesh::GeometryBasics.Mesh,x::SVector{T},t;kwargs...) where T
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    n = normal(mesh[u])
    v = x-locate(mesh[u],x) # signed Euclidian distance
    d = sign(sum(v.*n))*√sum(abs2,v) 
    return d,n
end # 120.029 ns (0 allocations: 0 bytes)d # 4.266 μs (0 allocations: 0 bytes)

# use divergence theorem to calculate volume of surface mesh
# F⃗⋅k⃗ = -⨕pn⃗⋅k⃗ dS = ∮(C-ρgz)n⃗⋅k⃗ dS = ∫∇⋅(C-ρgzk⃗)dV = ρg∫∂/∂z(ρgzk⃗)dV = ρg∫dV = ρgV #
volume(a::GeometryBasics.Mesh) = mapreduce(T->center(T).*dS(T),+,a)
volume(body::MeshBody) = volume(body.mesh)

import WaterLily: interp
function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,3},δ) where T
    c=center(tri);n=normal(tri);ar=area(tri);
    ar.*n.*interp(c.+1.5 .+ δ.*n, p)
end
function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,2},δ) where T
    c=center(tri)[1:2];n=normal(tri)[1:2];ar=area(tri);
    ar.*n.*interp(c.+1.5 .+ δ.*n, p)
end
forces(a::GeometryBasics.Mesh, flow::Flow, δ=2) = map(T->get_p(T,flow.p,δ), a)
forces(body::MeshBody, b::Flow, δ=2) = forces(body.mesh, b, δ)

using Printf: @sprintf
import WaterLily
using WriteVTK
# access the WaterLily writer to save the file
function WaterLily.write!(w,a::MeshBody)
    k = w.count[1]
    points = hcat([[p.data...] for p ∈ a.mesh.position]...)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, Base.to_index.(face)) for face in faces(a.mesh)]
    vtk = vtk_grid(w.dir_name*@sprintf("/%s_%06i", w.fname, k), points, cells) 
    for (name,func) in w.output_attrib
        vtk[name] = func(a)
    end
    vtk_save(vtk); w.count[1]=k+1
    w.collection[k]=vtk
end