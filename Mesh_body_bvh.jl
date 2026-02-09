using FileIO, MeshIO
using GeometryBasics
using WaterLily: AbstractBody,@loop,measure
using StaticArrays
using ForwardDiff
using WaterLily
using StaticArrays 
using KernelAbstractions
# using BiotSavartBCs
using CUDA
using LinearAlgebra

#STRUCTS ------------

struct BoundBox{T}
    lo::SVector{3, T}
    up::SVector{3, T}
end

BoundBox(lo::NTuple{3,T}, up::NTuple{3,T}) where T = BoundBox{T}(SVector{3,T}(lo), SVector{3,T}(up))
BoundBox(lo::SVector{3,T}, up::SVector{3,T}) where T = BoundBox{T}(lo, up)

struct TreeInfo
    lvl::Int
    n_nodes::Int
    n_leaves::Int
    n_meshelements::Int
end

struct BVH_simple{T}
    mesh::AbstractMesh
    nodes::AbstractVector{BoundBox{T}}
    leaves::AbstractVector{BoundBox{T}}
    all_nodes::AbstractVector{BoundBox{T}}
    entries::AbstractVector{AbstractMesh}
    info:: TreeInfo
end

mutable struct MeshBody{T,F<:Function} <: AbstractBody
    mesh  :: GeometryBasics.Mesh
    mesh0 :: GeometryBasics.Mesh
    velocity :: GeometryBasics.Mesh
    surf_id :: NTuple
    map   :: F
    bvh  :: BVH_simple{T}
    scale :: T
    half_thk::T #half thickness
    boundary::Bool
end
#-------------------

# CONSTRUCTOR FUNCTIONS ----------------------

function MeshBody(fname::String;map=(x,t)->x,scale=1.0,thk=0f0,boundary=true,T=Float32,lvl,mem,kwargs...)
    if endswith(fname,".inp")
        mesh0,srf_id = load_inp(fname)
    else
        mesh0 = load(fname)
        srf_id = ntuple(i->(1,i),length(mesh0))
    end
    return MeshBody(mesh0,srf_id,map,scale,thk,boundary,T,lvl,mem,kwargs...)
end
MeshBody(mesh::M;map=(x,t)->x,scale=1.0,thk=0f0,boundary=true,T=Float32,lvl=4,mem,kwargs...) where M<:GeometryBasics.Mesh = 
        MeshBody(mesh,ntuple(i->(1,i),length(mesh)),map,scale,thk,boundary,T,lvl,mem)
   
function MeshBody(mesh0::M,srf_id,map,scale,thk,boundary,T,lvl,mem;kwargs...) where M<:GeometryBasics.Mesh

    # SOMETHING WITH MEM
    points = Point3f[] # scale and map the points to the corrcet location
    for (i,pnt) in enumerate(mesh0.position)
        push!(points, Point3f(map(SA[pnt.data...]*scale,0)))
    end
    mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(mesh0))  #we lose the normals here?
    velocity = GeometryBasics.Mesh(zero(points),GeometryBasics.faces(mesh0))
    
    bvh=BVH_simple(mesh,lvl,T=T,mem=mem;kwargs...)
    return MeshBody(mesh,mesh0,velocity,srf_id,map,bvh,T(scale),T(thk/2),boundary)
end

Base.copy(a::GeometryBasics.Mesh) = GeometryBasics.Mesh(a.position,GeometryBasics.faces(a));
Base.copy(b::MeshBody) = MeshBody(copy(b.mesh),copy(b.mesh0),copy(b.velocity),b.surf_id,b.map,
                                  Rect(b.bvh),b.scale,b.half_thk,b.boundary)





function BVH_simple(fname::String,lvl::Int; overlap::Int =1, box0_exp::Int=2,T::DataType = Float32,mem=Array)
    mesh=load(file)
    return BVH_simple(mesh,lvl, overlap=overlap, box0_exp=box0_exp,T=T)
end

function BVH_simple(mesh::AbstractMesh, lvl::Int; overlap::Int=1, box0_exp::Int=4, T::DataType=Float32, mem=CuArray)
    
    box0 = Bounding_BBox(mesh, box0_exp, T)
    SubD = make_boxes_old(box0, lvl) |> mem
    positions= mesh.position |> mem

    @loop SubD[I] = expand_box(SubD[I], overlap) over I ∈ CartesianIndices(SubD)

    # Now using mem for arrays
    #instead of 
    entries = Vector{AbstractMesh}(undef, length(SubD))      # split mesh parts
    leaves  = mem{BoundBox{T}}(undef, length(SubD))       # bounding boxes

    @loop entries[I] = split_mesh(SubD[I], mesh) over I ∈ CartesianIndices(SubD)
    @loop leaves[I]   = Bounding_BBox(entries[I], 2, T) over I ∈ CartesianIndices(SubD)

    nodes = construct_nodes(leaves, lvl)  # make sure this is GPU-safe too
    info = TreeInfo(
        lvl,
        length(nodes),
        length(leaves),
        length(mesh.faces)
    )

    # Use vcat from GPU-compatible libraries (e.g. CUDA.jl has CUDA.vcat) if needed
    return BVH_simple{T}(mesh, nodes, leaves, vcat(nodes, leaves), entries, info)
end

# function BVH_simple( mesh::AbstractMesh, lvl ::Int ; overlap::Int =1, box0_exp::Int=4,T::DataType = Float32,mem=Array) 
#     box0=Bounding_BBox(mesh,box0_exp,T)
#     SubD=make_boxes_old(box0,lvl)
#     # box0::BoundBox{T},lvl::Int}

#     @loop SubD[I] = expand_box(SubD[I], overlap) over I ∈ CartesianIndices(SubD)  #eachindex(SubD)

#     entries = Vector{AbstractMesh}(undef, length(SubD))      # Holds the split mesh parts
#     leaves = Vector{BoundBox{T}}(undef, length(SubD))        # Bounding boxes around mesh parts

#     @loop entries[I] = split_mesh(SubD[I], mesh) over I ∈ CartesianIndices(SubD)  

#     @loop leaves[i] = Bounding_BBox(entries[i], 2, T) over i ∈ CartesianIndices(SubD) 
    
#     # SubD=[expand_box(leaf,overlap) for leaf in SubD]
#     # entries=[split_mesh(box,mesh) for box in SubD]

#     # leaves=[Bounding_BBox(mesh_part,4,T) for mesh_part in entries];
#     nodes=construct_nodes(leaves,lvl) 
#     info = TreeInfo(
#             lvl,
#         length(nodes),
#         length(leaves),
#         length(mesh.faces)
#         )

#     return BVH_simple{T}(mesh,nodes, leaves,vcat(nodes,leaves), entries,info)
# end

#BVH HELPER FUNCTIONS---------------

function Bounding_BBox(mesh_part::AbstractMesh,exp::Real,T::DataType)
    rec=Rect{3,T}(mesh_part.position) 
    bbox=boxtobox(rec) 
    box0=expand_box(bbox,exp)
    return box0
end

function construct_nodes(leaves::Vector{BoundBox{T}}, lvl::Int) where T
    n_nodes =  2^(lvl-1)-1
    all_boxes = vcat(Vector{BoundBox{T}}(undef, n_nodes),copy(leaves)) 

    for i in n_nodes:-1:1
        left = 2i
        right = 2i + 1
        parent= merge_boxes(all_boxes[left], all_boxes[right])
        all_boxes[i] = parent
    end

    return all_boxes[1:n_nodes]
end

function merge_boxes(box1::BoundBox{T}, box2::BoundBox{T}) where T
    
    is_valid(b) = all(isfinite, b.lo) && all(isfinite, b.up)

    if is_valid(box1) && !is_valid(box2)
        return box1
    elseif !is_valid(box1) && is_valid(box2)
        return box2
    elseif !is_valid(box1) && !is_valid(box2)
        # Both invalid → return a "NaN" box
        return BoundBox{T}(SVector(NaN, NaN, NaN), SVector(NaN, NaN, NaN))
    end

    lo = SVector(
        min(min(box1.lo[1], box1.up[1]), min(box2.lo[1], box2.up[1])),
        min(min(box1.lo[2], box1.up[2]), min(box2.lo[2], box2.up[2])),
        min(min(box1.lo[3], box1.up[3]), min(box2.lo[3], box2.up[3]))
    )
    up = SVector(
        max(max(box1.lo[1], box1.up[1]), max(box2.lo[1], box2.up[1])),
        max(max(box1.lo[2], box1.up[2]), max(box2.lo[2], box2.up[2])),
        max(max(box1.lo[3], box1.up[3]), max(box2.lo[3], box2.up[3]))
    )

    return BoundBox{T}(lo, up)
end


function inbox(x::AbstractVector{T}, box::BoundBox{U}) where {T,U}
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end

function boxtobox(box::HyperRectangle{3,T}) where T
    # 0 allocations
   BoundBox{T}(box.origin,box.origin .+box.widths)
end

function expand_box(box:: BoundBox{T}, overlap::Real) where T
    BoundBox{T}(box.lo .-overlap, box.up .+ overlap);
end

function split_box(box::BoundBox{T}) where T  #15 ns 0 allocations
    w= box.up .- box.lo
    max_dir=argmax(w)
    new_width= SVector(
        max_dir ==1 ? w[1]/2 : w[1],
        max_dir ==2 ? w[2]/2 : w[2],
        max_dir ==3 ? w[3]/2 : w[3]
    )
   
    if any(new_width .< 1)
        throw(ArgumentError("Width of a box can not be less then 1 cell"))
    end
    o= SVector(box.lo)
    o2 = SVector(
    max_dir ==1 ? o[1]+new_width[1] : o[1],
    max_dir ==2 ? o[2]+new_width[2] : o[2],
    max_dir ==3 ? o[3]+new_width[3] : o[3],
    )
    return BoundBox{T}(o,o .+new_width),BoundBox{T}(o2,o2 .+new_width)
end

function make_boxes_old(box0::BoundBox{T},lvl::Int) where{T}
    size=2^lvl -1
    n_leaves=2^(lvl-1)
    n_nodes = size - n_leaves

    box_array = Vector{BoundBox{T}}(undef, 2^lvl -1)
    box_array[1]=box0

    for i in 1:2^(lvl-1)-1
        child1,child2=split_box(box_array[i])
        box_array[2i]=child1
        box_array[2i+1]=child2
    end
    SubDivisions = box_array[n_nodes+1:end];
    return SubDivisions
end

function split_mesh(box::BoundBox{T}, mesh::AbstractMesh) where T
    faces = GeometryBasics.faces(mesh)
    max=length(mesh.faces)
    mask=fill(NaN32,length(mesh.position))
    count=0
    coll = Vector{eltype(faces)}(undef,max)
    @inbounds for face in faces
        for i in Tuple(face)
            v = mesh.position[i] 
            inside = all(box.lo .<= v .<= box.up)

            if inside
                count+=1
                coll[count]=face
                mask[face].=1
                break
            end
        end
    end
    # new_mesh=GeometryBasics.Mesh(mesh.position.*mask,coll[1:count],normal=mesh.normal.*mask)
    new_mesh=GeometryBasics.Mesh(mesh.position.*mask,coll[1:count])
    return new_mesh
end

# MESHBODY FUNCTIONS
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

Find the closest point x on the triangle tri to the point p.
"""
function locate(tri::GeometryBasics.Ngon{3},p::SVector{T}) where T #5.327 ns (0 allocations: 0 bytes)
    # is point a closest?
    a, b, c = tri.points
    ab = b.-a
    ac = c.-a
    ap = p.-a
    d1 = sum(ab.*ap)
    d2 = sum(ac.*ap)
    # is point a closest?
    if ((d1 ≤ 0) && (d2 ≤ 0))
        return a
    end
    # is point b closest?
    bp = p.-b
    d3 = sum(ab.*bp)
    d4 = sum(ac.*bp)
    if ((d3 ≥ 0) && (d4 ≤ d3))
        return b
    end
    # is point c closest?
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
    # closest is interior to abc
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

    Measures a mesh body at point x and time t.
    If a mapping has been specified, the point x is first map to the new location ξ, and the mesh is measured there.
    We use a large bounding box around the mesh to avoid measuring  far away, where the actual distance doesn't matter as
    long as d>>O(1), this means that outside the bounding box, we return the distance of the bbox edges to the geom (d=8).
    If the mesh is not a boundary, we need to adjust the distance by the half thickness of the body.
"""

function WaterLily.measure(body::MeshBody,x::SVector{D},t;kwargs...) where D
    # if x is 2D, we need to make it 3D
    D==2 && (x=SVector{3}(x[1],x[2],0))
    d,n,v=traverse_fsm(x,body.bvh,body.velocity)
    !body.boundary && (d = abs(d)-body.half_thk)  # if the mesh is not a boundary, we need to adjust the distance
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
    # update bvh
    body.bvh=BVH_simple(body.mesh,body.bvh.info.lvl)

end

using LinearAlgebra: cross
"""
    normal(tri::GeometryBasics.Ngon{3})

Return the normal vector to the triangle tri.
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

Measure the distance d and normal n to the mesh at point x and time t.
"""
# linear shape function interpolation of the nodal velocity values at point p
function get_velocity(tri::Tr,vel::Tr,p::SVector{3,T}) where {Tr<:GeometryBasics.Ngon{3},T}
    dA = SVector{3,T}([sub_area(tri,p,Val{i}()) for i in 1:3])
    return SVector(sum(tri.points.*dA)/sum(dA))
end
@inline sub_area(t,p,::Val{1}) = area(GeometryBasics.Ngon(SA[Point(p),t[2],t[3]]))
@inline sub_area(t,p,::Val{2}) = area(GeometryBasics.Ngon(SA[t[1],Point(p),t[3]]))
@inline sub_area(t,p,::Val{3}) = area(GeometryBasics.Ngon(SA[t[1],t[2],Point(p)]))


# TRAVERSAL AND ELEMENT DISTANCE FUNCTIONS

function get_dist(mesh::M,velocity::M,x::SVector{T};kwargs...) where {M<:GeometryBasics.Mesh,T}
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    n,p = normal(mesh[u]),SVector(locate(mesh[u],x))
    s = x-p # signed Euclidian distance
    d = sign(sum(s.*n))*√sum(abs2,s)
    # v = get_velocity(mesh[u],velocity[u],p)
    return (d,n,zero(x)) #v
end # 120.029 ns (0 allocations: 0 bytes)d # 4.266 μs (0 allocations: 0 bytes)

function traverse_fsm(x::SVector{3,T},bvh::BVH_simple,velocity::M) where {M<:GeometryBasics.Mesh,T}
    # println("trying to otraverse ", x)
    fromParent =1
    fromSibling=2
    fromChild=3

    n_nodes= bvh.info.n_nodes
    @inline function inbox_idx(current::Int)
        box=bvh.all_nodes[current]
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end

   @inline sibling(current::Int) = current%2==0 ? current+1 : current-1  #0 alloc
   @inline isLeaf(current::Int) = current > n_nodes   # 0 alloc?
   @inline parent(current::Int)= fld(current,2)  # 0 alloc
    
    sol=(8.0, zero(x),zero(x))
    state=fromParent
    current=1
    while true
        if state== fromChild
            if current ==1
                break
            elseif current==2*(parent(current))
                current=sibling(current)
                state=fromSibling
            else
                current=parent(current)
                state=fromChild
            end
        continue

        elseif state==fromSibling
            hit=inbox_idx(current)
            if !hit
                current=parent(current)
                state=fromChild
            elseif isLeaf(current)
                sol_leaf = get_dist(bvh.entries[Int(current-bvh.info.n_nodes)],velocity,x)
                abs(sol_leaf[1])<abs(sol[1]) && (sol=sol_leaf)
                current=parent(current)
                state=fromChild
            else
                current=2*current
                state=fromParent
            end
        continue

        elseif state==fromParent
            hit=inbox_idx(current)
            if !hit && current ==1
                break
            elseif !hit
                current=sibling(current)
                state=fromSibling
            elseif isLeaf(current)
                sol_leaf = get_dist(bvh.entries[Int(current-bvh.info.n_nodes)],velocity,x)
                abs(sol_leaf[1])<abs(sol[1]) && (sol=sol_leaf)
                current=sibling(current)
                state=fromSibling
            else
                current=2*current
                state=fromParent
            end
        continue

        end
    end
    return sol
end


# use divergence theorem to calculate volume of surface mesh
# F⃗⋅k⃗ = -⨕pn⃗⋅k⃗ dS = ∮(C-ρgz)n⃗⋅k⃗ dS = ∫∇⋅(C-ρgzk⃗)dV = ρg∫∂/∂z(ρgzk⃗)dV = ρg∫dV = ρgV #
volume(a::GeometryBasics.Mesh) = mapreduce(T->center(T).*dS(T),+,a)
volume(body::MeshBody) = volume(body.mesh)

import WaterLily: interp


# Changes
# negative normal (positive pressure , facing flow should be positive but n is (-1 0 0))
# if normals are outward facing ofcourse!


function get_f(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,3},δ,::Val{true}) where T
    c=center(tri);n=normal(tri);ar=area(tri);
    # n[3]<0 && n.*-1
    ar.*-n.*interp(c.+1.5 .+ δ.*n, p)
end
function get_f(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,3},δ,::Val{false}) where T
    # println("deze")
    c=center(tri);n=normal(tri);ar=area(tri);
    p1=interp(c.+1.5 .+ δ.*n, p)
    p2=interp(c.+1.5 .- δ.*n, p)
    # n[3]<0 && n.*-1
    # p1 < -1.5 && (p1 =0)
    # p2 < -1.5 && (p2=0)
    ar.*n.*(p2 .- p1)
end


# function p_test(p1,p2,n)
#     f1= p1 .* n
#     f2 = p2 .* n
#     del_p=p1-p2
#     # n .*(f)
#     # n .* (p1-p2)
#     # f1-f2
#     (p2-p1) .*n
# end
#option to get normals from mesh easier check for outwards!
#- so ensure normals face outwards!
function get_f_normals(verts::NTuple{3, <:GeometryBasics.AbstractVector}, 
    norms::NTuple{3, <:GeometryBasics.AbstractVector}, p::AbstractArray{T,3},δ,::Val{true}) where T
    tri=GeometryBasics.Ngon(verts)
    c=center(tri);ar=area(tri);
    n=normalize((norms[1] + norms[2] + norms[3]) / 3)
 
    ar.*-n.*interp(c.+1.5 .+ δ.*n, p)
end

#use U and L from simulation# 2fx/U^2 /A
function get_f_normals(verts::NTuple{3, <:GeometryBasics.AbstractVector}, 
    norms::NTuple{3, <:GeometryBasics.AbstractVector}, p::AbstractArray{T,3},δ,::Val{false}) where T
    # println("got here")
    tri=GeometryBasics.Ngon(verts)
    c=center(tri);ar=area(tri);
    n=normalize((norms[1] + norms[2] + norms[3]) / 3)
    p1= interp(c.+1.5 .+ δ.*n, p)
    p2= interp(c.+1.5 .- δ.*n, p)
    ar.*n.*(p2 .- p1 )
end


#Changed get_p to get f now get p gets a pressure


function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,3},δ,::Val{true}) where T
    c=center(tri);n=normal(tri);ar=area(tri);
    interp(c.+1.5 .+ δ.*n, p)
end
function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,3},δ,::Val{false}) where T
    c=center(tri);n=normal(tri);ar=area(tri);
    (interp(c.+1.5 .+ δ.*n, p) , interp(c.+1.5 .- δ.*n, p))
end

#option to get normals from mesh easier check for outwards!
#- so ensure normals face outwards!
function get_p_normals(verts::NTuple{3, <:GeometryBasics.AbstractVector}, 
    norms::NTuple{3, <:GeometryBasics.AbstractVector}, p::AbstractArray{T,3},δ,::Val{true}) where T
    tri=GeometryBasics.Ngon(verts)
    c=center(tri);ar=area(tri);
    n=normalize((norms[1] + norms[2] + norms[3]) / 3)
 
   interp(c.+1.5 .+ δ.*n, p)
end

function get_p_normals(verts::NTuple{3, <:GeometryBasics.AbstractVector}, 
    norms::NTuple{3, <:GeometryBasics.AbstractVector}, p::AbstractArray{T,3},δ,::Val{false}) where T
    # println("got here")
    tri=GeometryBasics.Ngon(verts)
    c=center(tri);ar=area(tri);
    n=normalize((norms[1] + norms[2] + norms[3]) / 3)
    (interp(c.+1.5 .+ δ.*n, p) , interp(c.+1.5 .- δ.*n, p))
end

@inbounds @inline normal2D(tri::GeometryBasics.Ngon{3})=  SVector{2}(normal(tri)[1:2])
@inbounds @inline center2D(tri::GeometryBasics.Ngon{3}) = SVector{2}(center(tri)[1:2])

# function get_p(tri::GeometryBasics.Ngon{3},p::AbstractArray{T,2},δ,::Val{true}) where T
#     c=center2D(tri);n=normal2D(tri);ar=area(tri);
#     p = ar.*n.*interp(c.+1.5 .+ δ.*n, p)
#     return SA[p[1],p[2],zero(T)]
# end

function get_v(tri::GeometryBasics.Ngon{3},u::AbstractArray,δ,::Val{false})
    c=center(tri);n=normal(tri);ar=area(tri);
    # println("here")
    u_p = interp(c.+1.5 .+ δ.*n,u)
    u_n = interp(c.+1.5 .- δ.*n,u)
    return u_p,u_n
end

#velocity gradient with forward difference from 2 forward points
function get_v(tri::GeometryBasics.Ngon{3},u::AbstractArray,δ, ::Val{true})
    c=center(tri);n=normal(tri);ar=area(tri);
    # println("aa")
    u_p = interp(c.+1.5 .+ δ.*n,u)
    
    return u_p
end

function get_v_norm(verts::NTuple{3, <:GeometryBasics.AbstractVector}, 
    norms::NTuple{3, <:GeometryBasics.AbstractVector}, u::AbstractArray,δ,::Val{false}) where T
    
    tri=GeometryBasics.Ngon(verts)
    c=center(tri)
    n=normalize((norms[1] + norms[2] + norms[3]) / 3)

    u_p = interp(c.+1.5 .+ δ.*n,u)
    u_n = interp(c.+1.5 .- δ.*n,u)
    return (u_p,u_n)
end

function get_v_norm(verts::NTuple{3, <:GeometryBasics.AbstractVector}, 
    norms::NTuple{3, <:GeometryBasics.AbstractVector}, u::AbstractArray,δ,::Val{true}) where T
    
    tri=GeometryBasics.Ngon(verts)
    c=center(tri);
    n=normalize((norms[1] + norms[2] + norms[3]) / 3)
  
    u_p = interp(c.+1.5 .+ δ.*n,u)
    
    return u_p
end




function get_y_plus(tri::GeometryBasics.Ngon{3},u::AbstractArray,δ, ::Val{true},ν )
    c=center(tri);n=normal(tri);ar=area(tri);
    println(n)
    println(c.+1.5 .+ δ.*n)
   
    v1=interp(c.+1.5 .+ δ.*n,u)
    v2=interp(c.+1.5 .+ (2*δ).*n,u)
    v1_sh = norm(v1 - (dot(v1,n)*n))
  
    v2_sh = norm(v2 - (dot(v2,n)*n))
    
    slope=(2*v1_sh-v2_sh)/2δ
    τ= ν*slope
    ut=√(τ)
    y_plus=ut/ν
    return y_plus
end


function get_fric_force(tri :: GeometryBasics.Ngon{3},u:: AbstractArray, δ , ::Val{true}, ν)
    c=center(tri);n=normal(tri);ar=area(tri);

    v1=interp(c.+1.5 .+ δ.*n,u)
    v2=interp(c.+1.5 .+ (2*δ).*n,u)
    v_shear=v1 .-dot(v1,n)*n
    v1_sh = norm(v1 - (dot(v1,n)*n))
    v2_sh = norm(v2 - (dot(v2,n)*n))

    slope=(2*v1_sh-v2_sh)/2δ
    τ= ν*slope

    cf_vec =  τ/0.5 .* normalize(v_shear)

    return cf_vec
end



function get_fric_force(tri :: GeometryBasics.Ngon{3},u:: AbstractArray, δ, ::Val{false}, ν)
    c=center(tri);n=normal(tri);ar=area(tri);
    
    v1=interp(c.+1.5 .+ δ.*n,u)
    v2=interp(c.+1.5 .+ (2*δ).*n,u)
    v_shear=v1 .-dot(v1,n)*n
    v1_sh = norm(v1 - (dot(v1,n)*n))
    v2_sh = norm(v2 - (dot(v2,n)*n))

    slope=(2*v1_sh-v2_sh)/2δ
    τ= ν*slope

    cf_vec =  τ/0.5 .* normalize(v_shear)

    n2=-n
    
    v1=interp(c.+1.5 .+ δ.*n2,u)
    v2=interp(c.+1.5 .+ (2*δ).*n2,u)
    v_shear=v1 .-dot(v1,n2)*n2
    v1_sh = norm(v1 - (dot(v1,n2)*n2))
    v2_sh = norm(v2 - (dot(v2,n2)*n2))

    slope=(2*v1_sh-v2_sh)/2δ
    τ= ν*slope

    cf_vec2 =  τ/0.5 .* normalize(v_shear)

    return cf_vec .+ cf_vec2

end
y_plus_coll(a::GeometryBasics.Mesh,u :: AbstractArray,ν,  δ=1,boundary=true)= map(T->get_y_plus(T,u,δ,Val{boundary}(),ν), a)
velocities(a::GeometryBasics.Mesh,u::AbstractArray,δ=1,boundary=true)  = :normal in propertynames(a) ? 
                     map(f -> get_v_norm(getindex.(Ref(a.position), Tuple(f)), getindex.(Ref(a.normal), Tuple(f)), u, δ,Val{boundary}()), a.faces) : 
                     map(T->get_v(T, u , δ,Val{boundary}()), a)


velocities(body::MeshBody,flow::Flow,δ=1) = velocities(body.mesh,flow.u,δ,body.boundary)
velocities(mesh::AbstractMesh,flow::Flow,δ=1,boundary=true) = velocities(mesh,flow.u,δ,boundary)
velocities(sim, δ=1) =velocities(sim.body.mesh,sim.flow.u,δ,sim.body.boundary)

# use copy to get parallel


fric_force(a:: GeometryBasics.Mesh, u :: AbstractArray, ν,δ=1, boundary =Val{boundary}()) =  map(T->get_fric_force(T,u,δ,boundary,ν), a)



pressures(a::GeometryBasics.Mesh,p::AbstractArray,δ=1,boundary=true)  = :normal in propertynames(a) ? 
                     map(f -> get_p_normals(getindex.(Ref(a.position), Tuple(f)), getindex.(Ref(a.normal), Tuple(f)), p, δ,Val{boundary}()), a.faces) : 
                     map(T->get_p(T, p , δ,Val{boundary}()), a)

pressures(body::MeshBody,flow::Flow,δ=2) = pressures(body.mesh,flow.u,δ,body.boundary)
pressures(mesh::AbstractMesh,flow::Flow,δ=2,boundary=true) = pressures(mesh,flow.p,δ,boundary)
pressures(sim, δ=2) =pressures(sim.body.mesh,sim.flow.p,δ,sim.body.boundary)
                     
# force on object depending on if normals are present or not
# #body.mesh does not store normals

forces(a::AbstractMesh, p::AbstractArray; δ=2, boundary=Val{true}()) = :normal in propertynames(a) ? 
                     map(f -> get_f_normals(getindex.(Ref(a.position), Tuple(f)), getindex.(Ref(a.normal), Tuple(f)), p, δ, boundary), a.faces) : 
                     map(T->get_f(T, p , δ, boundary), a)
#  map(f -> get_f_normals(getindex.(Ref(a.position), Tuple(f)), getindex.(Ref(a.normal), Tuple(f)), p|>Array, δ, boundary), a.faces)



# forces(a::AbstractMesh, p::AbstractArray; δ=2, boundary=true) = map(T->get_f(T, p, δ, Val{boundary}()), a)


forces(body::MeshBody, p::AbstractArray, δ=2) = forces(body.mesh,p,δ,body.boundary)
forces(body::MeshBody, b::Flow, δ=2) = forces(body.mesh,b.p,δ,body.boundary)
                                                                   
# forces(a::GeometryBasics.Mesh, flow::Flow, δ=2, boundary=Val{true}()) = map(T->get_f(T, flow.p, δ, boundary), a)
# forces(body::MeshBody, b::Flow, δ=2) = forces(body.mesh, b, δ, Val{body.boundary}())
# forces(body::MeshBody,p::AbstractArray, δ=2) = map(T -> get_f(T))


function pres_force(mesh::AbstractMesh,sim;δ=2,boundary=true)
    f= forces(mesh,sim.flow.p |> Array, δ=δ,boundary=Val{boundary}())
    sum(f)
end



function pres_force(body::MeshBody,sim;δ=2)
    f= forces(body.mesh,sim.flow.p |> Array, δ=δ,boundary=Val{body.boundary}())
    sum(f)
end

function CF(mesh::AbstractMesh,sim;δ=2,boundary=true)
    cfs=fric_force(mesh,sim.flow.u |> Array,sim.flow.ν,δ,Val{boundary}())
    sum(cfs)
end


function CF(body::MeshBody,sim;δ=2)
    cfs=fric_force(body.mesh,sim.flow.u |> Array,sim.flow.ν,δ,Val{body.boundary}())
    sum(cfs)
end



using Printf: @sprintf
import WaterLily
using WriteVTK

# access the WaterLily writer to save the file
function WaterLily.write!(w,a::MeshBody,t=w.count[1]) #where S<:AbstractSimulation{A,B,C,D,MeshBody}
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


#PRETTY PRINTS

function Base.show(io::IO, bvh::BVH_simple{T}) where T
    println(io, "BVH: ")
    println(io, "       ├─ type used:  $T")
    println(io, "       ├─ levels:     $(bvh.info.lvl)")
    println(io, "       ├─ Nodes:      $(bvh.info.n_nodes)")
    println(io, "       ├─ Leaves:     $(bvh.info.n_leaves)")
    println(io, "       ├─ mesh size:  $(bvh.info.n_meshelements)")
   
end



function Base.show(io::IO,body::MeshBody)
    println(io, "Body")
    # println(io,"├─surf_id: $(body.surf_id)")body
    println(io,"    ├─ boundary: $(body.boundary)")
    body.boundary==false && println(io,"    ├─ half_thk: $(body.half_thk)")
    println(io,"    ├─ scale: $(body.scale)")
    # println(io,"├─mesh size: $(body.bvh.info.n_meshelements)")
    println(io,"    ├─ $(body.bvh)")
end
