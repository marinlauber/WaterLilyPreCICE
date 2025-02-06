using StaticArrays
using WaterLily
# using BenchmarkTools
# using Plots
using GeometryBasics

abstract type AbstractBox end

mutable struct Tree
    boxes::NTuple
    points::AbstractArray
    C::SVector # these point to boxes[1]
    R::SVector # these point to boxes[1]
end
function Tree(points::AbstractArray{T,2}; δ=4, Nmax=64) where T
    b = bounding_box(points,δ=δ) 
    boxes,idx = AbstractBox[b],Int[1]
    # split the box in the longest direction
    left, right = split(b)
    # find which box each point belongs to
    points_left = filter(points, left)
    points_right = filter(points, right) # we can avoid this actually
    # downsample the boxes
    tree!(boxes, idx, 2, @view(points[:,points_left]), left, Nmax)
    tree!(boxes, idx, 3, @view(points[:,points_right]), right, Nmax)
    # sort the boxes by index and return, we need to fill some indices empty boxes
    # as we then access child of boxes[I] as boxes[2I] and boxes[2I+1]
    boxes = ntuple(i->i in idx ? boxes[findfirst(==(i),idx)] : NoBbox(T), maximum(idx))
    Tree(boxes,points,boxes[1].C,boxes[1].R)
end

Base.getindex(t::Tree, i::Int64) = t.boxes[i]
Base.length(t::Tree) = length(t.boxes)

# bbox have C (center) and R (width)
struct Bbox{S<:SVector} <: AbstractBox
    C::S
    R::S
    leaf::Bool
    indices::Union{Nothing,AbstractArray}
    function Bbox(C::S,R::S,leaf::Bool=false,indices=nothing) where S
        new{S}(C,R,leaf,indices)
    end
end
NoBbox(T) = Bbox(SA{T}[0,0],SA{T}[0,0],false)
# box = Bbox(SA[0.0,0.0,0.0],SA[0.5,0.5-eps(),0.5-eps()])

# split the box in the longest side
split_w(width::SVector{N},j) where N = SA[ntuple(i -> i==j ? width[i]/2 : width[i], N)...]
function Base.split(b::Bbox)
    # split the longest side
    w = split_w(b.R, argmax(b.R))
    return Bbox(b.C-(b.R-w),w), Bbox(b.C+(b.R-w),w)
end
# left,right = split(box)
# @btime split_w($box.R,$(argmax(box.R))) # 2.266 ns , 0 bytes
# @assert all(left.C.≈SA[-0.25,0.,0.]) && all(left.R.≈SA[0.25,0.5,0.5])
# @assert all(right.C.≈SA[0.25,0.,0.]) && all(right.R.≈SA[0.25,0.5,0.5])

# check if a point is outside the box
WaterLily.inside(x,b::Union{Tree,Bbox})::Bool = (all(b.C-b.R .≤ x) && all(x .≤ b.C+b.R))
# point = [0.,0.,0.]
# b = Bbox(SA[0.0,0.0,0.0],SA[0.5,0.5,0.5])
# @assert inside(point,b) 
# @btime inside($point,$b) # 4.985 ns , 0 bytes
# @code_warntype inside(point,b)
# @assert !inside([0.5+eps(),0.5,0.5],b) && !inside([-0.5-eps(),0.5,0.5],b)

# function Rect{N1,T1}(geometry::AbstractArray{PT}) where {N1,T1,PT<:Point}
#     N2, T2 = length(PT), eltype(PT)
#     @assert N1 >= N2
#     vmin = Point{N2,T2}(typemax(T2))
#     vmax = Point{N2,T2}(typemin(T2))
#     for p in geometry
#         vmin, vmax = _minmax(p, vmin, vmax)
#     end
#     o = vmin
#     w = vmax - vmin
#     return if N1 > N2
#         z = zero(Vec{N1 - N2,T1})
#         Rect{N1,T1}(vcat(o, z), vcat(w, z))
#     else
#         Rect{N1,T1}(o, w)
#     end
# end
# function _minmax(p::StaticVector, vmin, vmax)
#     any(isnan, p) && return (vmin, vmax)
#     return min.(p, vmin), max.(p, vmax)
# end

function bounding_box(points::AbstractArray{T,2};δ=1e-2) where T
    # vmax = SVector(ntuple(i->typemin(T),size(points,1)))
    vmax = SVector{3}(typemin(T),typemin(T),typemin(T))
    vmin = SVector{3}(typemax(T),typemax(T),typemax(T))
    for p in eachcol(points)
        p = SA[p...] # convert to static array
        vmin = min.(p, vmin)
        vmax = max.(p, vmax)
    end
    o = (vmin + vmax)/2
    r = (vmax - vmin)/2 .+ T(δ) # make it a bit bigger
    return Bbox(o,r,false)
end
function bounding_box(mesh::GeometryBasics.Mesh{3,T};δ=1e-2) where T
    r = Rect(mesh.position)
    return Bbox(SA[r.origin.+r.widths./2...],SA[r.widths...],false,nothing)
end
# function bounding_box(points::AbstractArray{T,2};δ=1e-2) where T
#     mn,mx = extrema(points,dims=2)
#     vmin = first.(mn); vmax = first.(mx)
#     o = (vmin + vmax)/2
#     r = (vmax - vmin)/2 .+ δ # make it a bit bigger
#     return Bbox(o,r,false)
# end

# some points
# points = rand(3,10)
# box = bounding_box(points)
# @btime bounding_box($points) # 2.071 μs (115 allocations: 5.19 KiB)# 25.775 μs (2594 allocations: 104.95 KiB)
# @code_warntype bounding_box(points)
# left, right = split(box)

# check if points are inside the box
Base.filter(pts,b::Bbox) = findall(x->inside(x,b),eachcol(pts))
# @btime filter($points,$left) # 273.669 ns (25 allocations: 528 bytes)
# @btime filter($points,$right)
# @assert length(p_left) + length(p_right) == size(points,2)
# @assert !(p_left in p_right) # check that there is no point in both boxes


function tree!(list, idx, Is, points::AbstractArray{T,2}, b::Bbox, Nmax) where T
    push!(idx,Is) # add the index location
    # if we are on a leaf, we push a leaf box on the list
    if (length(points)<Nmax)
        push!(list, Bbox(b.C,1.2*b.R,true,points))
        return nothing
    end
    push!(list, b);
    # we are not on a leaf node, we can downsample
    left, right = split(b)
    # find which box each point belongs to
    points_left = filter(points, left)
    points_right = filter(points, right)
    # reshape in case points are clustered on one side
    left = bounding_box(@view(points[:,points_left]))
    right = bounding_box(@view(points[:,points_right]))
    # downsample the boxes
    tree!(list, idx, Is*2, @view(points[:,points_left]), left, Nmax)
    tree!(list, idx, Is*2+1, @view(points[:,points_right]), right, Nmax)
    return nothing
end

function WaterLily.measure(tree::Tree, x::SVector)
    # are we in the bbox of the geom?
    !inside(x,tree) && return dist(tree.C,tree.R,x)
    # if we are inside the bbox, which sub-box are we in?
    # @TODO maybe we need to check if we are in the right one explicitly
    inleft = inside(x,tree[2]) # if we are not in the left box, 
    inright =inside(x,tree[3]) # we might be in the right one
    return _measure(tree, inleft ? 2 : 3, x)
end

function _measure(tree::Tree,Is,x::SVector)
    b = tree[Is] # points
    # if we are on a leaf, check that we are inside the box and measure appropriatly
    b.leaf && return !inside(x,b) ? dist(tree.C,tree.R,x) : dist(b.indices,x)
    # if not, we check which box we are in
    inleft = inside(x,tree[2Is]) # if we are not in the left box, 
    inright =inside(x,tree[2Is+1]) # we might be in the right one
    # if we are not in either, we measure to box boxes
    if !inleft && !inright # if we are in neither, return min distance to box
        return min(dist(tree[2Is].C,tree[2Is].R,x), dist(tree[2Is+1].C,tree[2Is+1].R,x))
    end
    return _measure(tree, inleft ? 2Is : 2Is+1, x)
end
# point_dist(points,indices,x::SVector) = √minimum(sum(abs2,points[:,indices].-x,dims=1))
dist(C,R,x) = √sum(abs2,max.(0,abs.(x-C)-R))+min(R...) # distance to a bbox
dist(points,x) = √minimum(sum(abs2,points.-x,dims=1)) # distance to points
brute_force(points,x) = √minimum(sum(abs2,points.-x,dims=1))


# # make a simple tree
# N = 2^8
# points = rand(2,N)
# # try a circle
# points[1,:] .= N*cos.(range(0,stop=2π,length=N)) .+ 2N
# points[2,:] .= N*sin.(range(0,stop=2π,length=N)) .+ 2N
# segments = map(i->[i,i%N+1],1:N) # connectivity of nodes

# # make a tree
# tree = Tree(points; Nmax=64);

# let
#     plot(title="KDTree",dpi=300, aspect_ratio=:equal,legend=:none)
#     for (i,b) in enumerate(tree.boxes)
#         b.C==b.R && continue # dummy box
#         p1 = b.C - b.R
#         p2 = b.C + SA[b.R[1],-b.R[2]]
#         p3 = b.C + b.R
#         p4 = b.C - SA[b.R[1],-b.R[2]]
#         plot!([p1[1],p2[1],p3[1],p4[1],p1[1]],[p1[2],p2[2],p3[2],p4[2],p1[2]],
#                color=:black, lw=2,fill=false,alpha=0.3,label=:none)
#         c = ifelse(b.leaf,:forestgreen,:brown3)
#         annotate!(b.C[1], b.C[2], ("$i", 16, c, :center))
#     end
#     plot!(points[1,:],points[2,:],seriestype=:scatter,label=:none)
# end

# # test some measure
# x = 2rand(SVector{2}).-1
# C = tree.C; R = tree.R
# @btime dist($C,$R,$x) # 2.077 ns (0 allocations: 0 bytes)
# pts_selected = tree.boxes[20].indices;
# @btime dist($pts_selected,$x) #162.201 ns (4 allocations: 768 bytes) 

# @btime measure($tree,$x) #69 ns (7 allocations: 208 bytes) # 1.212 μs (37 allocations: 2.70 KiB)
# @btime brute_force($points,$x) #1.689 μs (10 allocations: 6.25 KiB)

# # test on a grid
# M = 2^10
# σ = zeros(M,M)
# @btime @inside σ[I] = WaterLily.μ₀(brute_force(tree.points,loc(0,I)).-4,1) # 887.248 ms (13578572 allocations: 6.27 GiB)
# @btime @inside σ[I] = WaterLily.μ₀(measure(tree,loc(0,I)).-4,1) #160.230 ms (17850725 allocations: 813.04 MiB) # 254.503 ms (29340049 allocations: 1.29 GiB)
# # @btime @inside σ[I] = measure(tree,loc(0,I))
# let
#     flood(σ,levels=30)
#     plot!(title="KDTree", dpi=300, aspect_ratio=:equal,legend=:none)
#     for (i,b) in enumerate(tree.boxes)
#         b.C==b.R && continue # dummy box
#         p1 = b.C - b.R
#         p2 = b.C + SA[b.R[1],-b.R[2]]
#         p3 = b.C + b.R
#         p4 = b.C - SA[b.R[1],-b.R[2]]
#         plot!([p1[1],p2[1],p3[1],p4[1],p1[1]],[p1[2],p2[2],p3[2],p4[2],p1[2]],
#                color=:black, lw=2,fill=false,alpha=0.3,label=:none)
#         c = ifelse(b.leaf,:forestgreen,:brown3)
#     end
#     plot!()
# end
# include("MeshBodies.jl")
# sphere = load_inp("src/sphere.inp")[1]
# points_s = hcat([[a...] for a in sphere.position]...)

# @btime bounding_box($sphere) #4.948 μs (0 allocations: 0 bytes)
# @btime bounding_box($points_s) #548.199 μs (22375 allocations: 652.62 KiB)

# function Tree(mesh::GeometryBasics.Mesh{3,T}; δ=4, Nmax=64) where T
#     b = bounding_box(mesh,δ=δ) 
#     boxes,idx = AbstractBox[b],Int[1]
#     # split the box in the longest direction
#     left, right = split(b)
#     # find which box each point belongs to
#     # points_left = filter(points, left)
#     # points_right = filter(points, right) # we can avoid this actually
#     # downsample the boxes
#     # tree!(boxes, idx, 2, @view(points[:,points_left]), left, Nmax)
#     # tree!(boxes, idx, 3, @view(points[:,points_right]), right, Nmax)
#     # sort the boxes by index and return, we need to fill some indices empty boxes
#     # as we then access child of boxes[I] as boxes[2I] and boxes[2I+1]
#     boxes = ntuple(i->i in idx ? boxes[findfirst(==(i),idx)] : NoBbox(T), maximum(idx))
#     Tree(boxes,points,boxes[1].C,boxes[1].R)
# end
# Base.filter(mesh::GeometryBasics.Mesh,b::Bbox) = findall(x->inside(x,b),mesh)
# WaterLily.inside(tri::GeometryBasics.Ngon{3},b::Bbox)::Bool = (x=center(tri);(all(b.C-b.R .≤ x) && all(x .≤ b.C+b.R)))
# tree = Tree(sphere; Nmax=64)

# btime(b) = minimum(b).time
# res = 2 .^ collect(5:1:10)
# # test scaling compared to brute force
# duration = []; for N ∈ res
#     @show N
#     points = rand(2,N)
#     points[1,:] .= N*cos.(range(0,stop=2π,length=N)) .+ 2N
#     points[2,:] .= N*sin.(range(0,stop=2π,length=N)) .+ 2N
#     tree = Tree(points; Nmax=64);
#     σ = zeros(4N,4N)
#     push!(duration, btime(@benchmark @inside σ[I] = measure(tree,loc(0,I))))
#     push!(duration, btime(@benchmark @inside σ[I] = brute_force(points,loc(0,I))))
# end
# plot(res.^2,duration[1:2:end],xlabel="N×M",ylabel="Time ratio",
#      label="KDTree",lw=2, axis=:log)
# plot!(res.^2,duration[2:2:end],label="Brute force",lw=2)
# savefig("KDTree_speedup.png")