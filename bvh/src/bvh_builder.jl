using GLMakie
using StaticArrays
using FileIO
using MeshIO
using GeometryBasics
using BenchmarkTools
#include when plotting for debugging
#  include("plot_box.jl")
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
    all_nodes::Vector{BoundBox{T}}
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


function BVH_simple( file ::String, lvl ::Int ; overlap::Int =1, box0_exp::Int=2,T::DataType = Float64) 
    mesh=load(file)
    box0=Bounding_BBox(mesh,box0_exp,T)
    SubD=make_boxes(box0,lvl)
    
    SubD=[expand_box(leaf,overlap) for leaf in SubD]
    entries=[split_mesh(box,mesh) for box in SubD]

    leaves=[Bounding_BBox(mesh_part,1,T) for mesh_part in entries];
    nodes=construct_nodes(leaves,lvl)
    info = TreeInfo(
         lvl,
        length(nodes),
        length(leaves),
        length(mesh.faces)
        )
    
    return BVH_simple{T}(mesh,nodes, leaves,vcat(nodes,leaves), entries,info)
end



function BVH_simple( mesh::AbstractMesh, lvl ::Int ; overlap::Int =1, box0_exp::Int=2,T::DataType = Float64) 
   
    box0=Bounding_BBox(mesh,box0_exp,T)
    SubD=make_boxes_old(box0,lvl)
    # box0::BoundBox{T},lvl::Int}
    SubD=[expand_box(leaf,overlap) for leaf in SubD]
    entries=[split_mesh(box,mesh) for box in SubD]

    leaves=[Bounding_BBox(mesh_part,1,T) for mesh_part in entries];
    nodes=construct_nodes(leaves,lvl)
    info = TreeInfo(
         lvl,
        length(nodes),
        length(leaves),
        length(mesh.faces)
        )
    
    return BVH_simple{T}(mesh,nodes, leaves,vcat(nodes,leaves), entries,info)
end


function BVH_simple_inplace( mesh ::AbstractMesh, lvl ::Int ; overlap::Int =1, 
    box0_exp::Int=2,T::DataType = Float64) 
    size=2^lvl -1
    n_leaves=2^(lvl-1)
    n_nodes=size-n_leaves
    #initialize box array only once 
    box_array = Vector{BoundBox{T}}(undef, 2^lvl -1)
    box0=Bounding_BBox(mesh,box0_exp,T)  #2 alloc for some reason but only once so fine for now
    box_array[1]=box0
    make_boxes(lvl,box_array)  #changes box_array in place
    
    # box array is now array of SubDivisions
    # expand_box(SubD,overlap) for SubD in box_array[]
    box_array[n_nodes+1:end]=[expand_box(leaf,overlap) for leaf in box_array[n_nodes+1:end]] # expands only subD boxes
    entries = [split_mesh(box,mesh) for box in box_array[n_nodes+1:end]]
    
    box_array[n_nodes+1:end]=[Bounding_BBox(mesh_part,1,T) for mesh_part in entries]
    # leaves=[Bounding_BBox(mesh_part,1,T) for mesh_part in entries]
    # println(leaves)
    box_array[begin:n_nodes]=construct_nodes(box_array[n_nodes+1:end],lvl)
    # println(nodes)
    info = TreeInfo(
         lvl,
        n_nodes,
        n_leaves,
        length(mesh.faces)
        )
    
    return BVH_simple{T}(mesh,box_array[begin:n_nodes], box_array[n_nodes+1:end],box_array, entries,info)
end


function Bounding_BBox(mesh_part::AbstractMesh,exp::Real,T::DataType)
    rec=Rect{3,T}(mesh_part.position) #0
    bbox=boxtobox(rec) #0
    box0=expand_box(bbox,exp)
    return box0
end

function construct_nodes(leaves::Vector{BoundBox{T}}, lvl::Int) where T
    n_nodes =  2^(lvl-1)-1
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
    #0 alloc
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


#only 7 ns slower then function with GeometryBasics 
function split_box(box::BoundBox{T}) where T  #15 ns
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

function make_boxes_old(box0::BoundBox{T},lvl::Int) where{T}
    size=2^lvl -1
    n_leaves=2^(lvl-1)
    n_nodes = size - n_leaves

    box_array = Vector{BoundBox{T}}(undef, 2^lvl -1)
    box_array[1]=box0

    for i in 1:2^(lvl-1)-1
        # println(i)
        child1,child2=split_box(box_array[i])
        box_array[2i]=child1
        box_array[2i+1]=child2
    end
    SubDivisions = box_array[n_nodes+1:end];
    return SubDivisions
end


function make_boxes(lvl::Int, box_array::Vector{BoundBox{T}}) where{T}

    for i in 1:2^(lvl-1)-1
        # println(i)
        child1,child2=split_box(box_array[i])
        box_array[2i]=child1
        box_array[2i+1]=child2
    end
    # SubDivisions = box_array[n_nodes+1:end];
    return box_array
end


#

function split_mesh(box::BoundBox{T}, mesh::AbstractMesh) where T
    #6 allocations
    faces = GeometryBasics.faces(mesh)
    # positions=copy(mesh.position)
    # norms=copy(mesh.normal)
    max=length(mesh.faces)
    mask=fill(NaN32,length(mesh.position))
    # points = SVector{3,T}.(mesh.position)
    count=0
    # mesh = GeometryBasics.Mesh(SVector{3,T}.(mesh.position), GeometryBasics.faces(mesh))
    coll = Vector{eltype(faces)}(undef,max)

    #not fully correct yet use found faces to mask points
    @inbounds for face in faces
        for i in Tuple(face)
            v = mesh.position[i]  # Now an SVector{3,T}
            inside = all(box.lo .<= v .<= box.up)

            if inside
                count+=1
                coll[count]=face
                mask[face].=1
                break
            end
        end
    end

    new_mesh=GeometryBasics.Mesh(mesh.position.*mask,coll[1:count],normal=mesh.normal.*mask)

    return new_mesh
end





# file="obj/dragon_15k.obj" # 0 allocs
# mesh=load(file)
# # mesh = load("obj/"*file*".obj") # 81745 allocations!
# @benchmark bvh=BVH_simple($mesh,4,overlap=2)
# bvh=BVH_simple(mesh,6,overlap=2)





# # # bbox3=BBox{Float32}((33.0f0, 3.0f0, 33.0f0), (67.0f0, 102.0f0, 73.0f0))
# file="/home/raj/thesis/WaterLilyPreCICE_KD/bvh/obj/cent_dragon.obj" # 0 allocs
# mesh = load(file) # 81745 allocations!
# bvh=BVH_simple(mesh,3)

#@benchmark BVH_simple($mesh,$4)
# benchmark for 4 levels
# BenchmarkTools.Trial: 2327 samples with 1 evaluation per sample.
#  Range (min … max):  1.723 ms …   7.108 ms  ┊ GC (min … max): 0.00% … 14.40%
#  Time  (median):     1.927 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   2.135 ms ± 456.916 μs  ┊ GC (mean ± σ):  5.62% ± 10.11%

#    ▁▅█▄▅▁                                                      
#   ▄██████▆▆▅▃▃▃▂▂▂▂▁▁▂▁▂▂▂▂▃▃▃▃▃▃▂▂▂▂▂▂▂▂▁▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
#   1.72 ms         Histogram: frequency by time         3.7 ms <

#  Memory estimate: 7.90 MiB, allocs estimate: 244.
 
# begin
# # plots for debugging
# fig = Figure(size = (1200, 800))
# ax = Axis3(fig[1, 1])
# # boxes=vcat(bvh.nodes,bvh.leaves)

# leaf_tocheck=4

# # wireframe!(ax, bvh.mesh, color = "black", ssao = true)
# wireframe!(ax, bvh.entries[1], color = "green", ssao = true)
# wireframe!(ax, bvh.entries[2], color = "green", ssao = true)
# wireframe!(ax, bvh.entries[3], color = "green", ssao = true)
# wireframe!(ax, bvh.entries[4], color = "green", ssao = true)




# # lines!(ax, boxes_lines(bvh.SubD), linewidth = 2, color = "grey", linestyle=:dash)


# # lines!(ax, boxes_lines([bvh.SubD[box_tocheck]]), linewidth = 2, color = "red", linestyle=:dash)
# lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 5, color = "grey")
# lines!(ax, boxes_lines(bvh.leaves), linewidth = 1, color = "green",linestyle=:solid)
# # lines!(ax, boxes_lines([bvh.leaves[4]]), linewidth = 5, color = "black",linestyle=:solid)
# # lines!(ax, boxes_lines([bvh.nodes[6]]), linewidth = 10, color = "black",linestyle=:solid)
# # lines!(ax, boxes_lines([merge_boxes(bvh.leaves[1],bvh.leaves[2])]), linewidth = 10, color = "green",linestyle=:dash)

# # scatter!(ax, [Point3f(x)], color = :red, markersize = 15)
# # fig
# end

