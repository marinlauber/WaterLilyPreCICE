using GLMakie
using StaticArrays
using FileIO
using MeshIO
using GeometryBasics
using BenchmarkTools
#include when plotting for debugging
 include("plot_box.jl")
# include("src/BBox.jl")

#own bvh builder, everything is 1 indexed


struct BoundBox{T}
    lo::SVector{3, T}
    up::SVector{3, T}
end
BoundBox(lo::NTuple{3,T}, up::NTuple{3,T}) where T = Boundbox{T}(SVector{3,T}(lo), SVector{3,T}(up))


struct TreeInfo
    lvl::Int
    n_nodes::Int
    n_leaves::Int
    n_meshelements::Int
end

struct BVH_simple{T}
    mesh::AbstractMesh
    nodes::Vector{BoundBox{T}}
    leaves::Vector{BoundBox{T}}
    SubD::Vector{BoundBox{T}}
    entries::Vector{AbstractMesh}
    info:: TreeInfo
end

function Base.show(io::IO, bvh::BVH_simple{T}) where T
    println(io, "BVH ")
    println(io, " ├─ type used:  $T")
    println(io, " ├─ levels:     $(bvh.info.lvl)")
    println(io, " ├─ Nodes:      $(bvh.info.n_nodes)")
    println(io, " ├─ Leaves:     $(bvh.info.n_leaves)")
    println(io, " ├─ mesh size:  $(bvh.info.n_meshelements)")
   
end


function BVH_simple( file ::String, lvl ::Int ; overlap::Int =1, T::DataType = Float64) 
    mesh=load(file)
    box0=Bounding_BBox(mesh,3,T)
    nodes,SubD=make_boxes(box0,lvl)
    
    SubD=[expand_box(leaf,overlap) for leaf in SubD]
    entries = [
    let result = split_mesh(box, mesh)
        faces, verts = result
        GeometryBasics.Mesh(verts, faces)
    end
    for box in SubD
    ]

    leaves=[Bounding_BBox(mesh_part,overlap/2,T) for mesh_part in entries];
    nodes=construct_nodes(leaves,lvl)
    info = TreeInfo(
         lvl,
        length(nodes),
        length(leaves),
        length(mesh.faces)
        )
    
    return BVH_simple{T}(mesh,nodes, leaves,SubD, entries,info)
end


function BVH_simple( mesh ::AbstractMesh, lvl ::Int ; overlap::Int =1, T::DataType = Float64) 
  
    box0=Bounding_BBox(mesh,3,T)
    nodes,SubD=make_boxes(box0,lvl)
    
    SubD=[expand_box(leaf,overlap) for leaf in SubD]
    entries = [
    let result = split_mesh(box, mesh)
        faces, verts = result
        GeometryBasics.Mesh(verts, faces)
    end
    for box in SubD
    ]

    leaves=[Bounding_BBox(mesh_part,overlap/2,T) for mesh_part in entries];
    nodes=construct_nodes(leaves,lvl)
    info = TreeInfo(
         lvl,
        length(nodes),
        length(leaves),
        length(mesh.faces)
        )
    
    return BVH_simple{T}(mesh,nodes, leaves,SubD, entries,info)
end


function Bounding_BBox(mesh_part::AbstractMesh,exp::Real,T::DataType)
    rec=Rect{3,T}(mesh_part.position)
    bbox=boxtobox(rec)
    box0=expand_box(bbox,exp)
    return box0
end

# function construct_nodes(leaves::Vector{BoundBox{T}}, lvl::Int) where T
#     total_nodes = 2^lvl - 1
#     n_leaves = 2^(lvl - 1)
#     n_internal =  length(leaves)-1

#     nodes = Vector{BoundBox{T}}(undef, n_internal)
#     all_boxes = vcat(nodes,copy(leaves))  # working array that will grow upward)

#     # Fill in parent nodes from leaves upward
#     for i in n_internal:-1:1
#         left = 2i
#         right = 2i + 1
#         parent= merge_boxes(all_boxes[left], all_boxes[right])
#         nodes[i] = parent
#         push!(all_boxes, nodes[i])  # extend the working array
#     end

#     return nodes
# end


function construct_nodes(leaves::Vector{BoundBox{T}}, lvl::Int) where T
    n_nodes =  length(leaves)-1
    all_boxes = vcat(Vector{BoundBox{T}}(undef, n_nodes),copy(leaves))  # working array that will grow upward)

    # Fill in parent nodes from leaves upward
    for i in n_nodes:-1:1
        left = 2i
        right = 2i + 1
        parent= merge_boxes(all_boxes[left], all_boxes[right])
        all_boxes[i] = parent
        # push!(all_boxes, nodes[i])  # extend the working array
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
        return BBox{T}(SVector(NaN, NaN, NaN), SVector(NaN, NaN, NaN))
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

function split_box(box::BoundBox{T}) where T
    # println(box)
    #0 allocations
    w= box.up .- box.lo
    max_dir=argmax(w)
    new_width= SVector(
        max_dir ==1 ? w[1]/2 : w[1],
        max_dir ==2 ? w[2]/2 : w[2],
        max_dir ==3 ? w[3]/2 : w[3]
    )
    # println(new_width)
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

function make_boxes(box0::BoundBox{T},lvl::Int) where{T}
    size=2^lvl -1
    n_leaves=2^(lvl-1)
    n_nodes = size - n_leaves

    # box_array=Vector{BoundBox{Float32}}
    box_array = Vector{BoundBox{T}}(undef, size)
    box_array[1]=box0

    for i in 1:2^(lvl-1)-1
        # println(i)
        child1,child2=split_box(box_array[i])
        
        box_array[2i]=child1
        box_array[2i+1]=child2
    end

    nodes = box_array[1:n_nodes]
    leaves = box_array[n_nodes+1:end];

    return nodes,leaves
end
function split_mesh(box::BoundBox{T}, mesh::AbstractMesh) where T
    faces = GeometryBasics.faces(mesh)
    positions=copy(mesh.position)
    points = SVector{3,T}.(mesh.position)

    # mesh = GeometryBasics.Mesh(SVector{3,T}.(mesh.position), GeometryBasics.faces(mesh))
    coll = Vector{eltype(faces)}()

    @inbounds for face in faces
        for i in Tuple(face)
            v = points[i]  # Now an SVector{3,T}

            # No need for two bound orders if box.lo < box.up guaranteed
            inside = all(box.lo .<= v .<= box.up)

            if inside
                push!(coll, face)
                break
            end
        end
    end

    for i in eachindex(positions)
        if !inbox(positions[i], box)
            positions[i] = Point3f(NaN32, NaN32, NaN32)
        end
    end

    return coll,positions
end






# # bbox3=BBox{Float32}((33.0f0, 3.0f0, 33.0f0), (67.0f0, 102.0f0, 73.0f0))
file="/home/raj/thesis/WaterLilyPreCICE_KD/bvh/obj/kite_large.obj" # 0 allocs
mesh = load(file) # 81745 allocations!
bvh=BVH_simple(mesh,4)


# benchmark for 4 levels
# BenchmarkTools.Trial: 2327 samples with 1 evaluation per sample.
#  Range (min … max):  1.723 ms …   7.108 ms  ┊ GC (min … max): 0.00% … 14.40%
#  Time  (median):     1.927 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   2.135 ms ± 456.916 μs  ┊ GC (mean ± σ):  5.62% ± 10.11%

#    ▁▅█▄▅▁                                                      
#   ▄██████▆▆▅▃▃▃▂▂▂▂▁▁▂▁▂▂▂▂▃▃▃▃▃▃▂▂▂▂▂▂▂▂▁▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
#   1.72 ms         Histogram: frequency by time         3.7 ms <

#  Memory estimate: 7.90 MiB, allocs estimate: 244.
 
begin
# plots for debugging
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1])
# boxes=vcat(bvh.nodes,bvh.leaves)

leaf_tocheck=2

wireframe!(ax, bvh.mesh, color = "black", ssao = true)

# lines!(ax, boxes_lines(bvh.SubD), linewidth = 2, color = "grey", linestyle=:dash)


# lines!(ax, boxes_lines([bvh.SubD[box_tocheck]]), linewidth = 2, color = "red", linestyle=:dash)
lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 10, color = "grey")
lines!(ax, boxes_lines([bvh.nodes[5]]), linewidth = 1, color = "orange",linestyle=:solid)
lines!(ax, boxes_lines(bvh.leaves[3:4]), linewidth = 5, color = "black",linestyle=:solid)
# lines!(ax, boxes_lines([bvh.nodes[6]]), linewidth = 10, color = "black",linestyle=:solid)
# lines!(ax, boxes_lines([merge_boxes(bvh.leaves[7],bvh.leaves[8])]), linewidth = 10, color = "green",linestyle=:solid)

# scatter!(ax, [Point3f(x)], color = :red, markersize = 15)
fig
end


save("entrie1.obj", bvh.entries[1])
