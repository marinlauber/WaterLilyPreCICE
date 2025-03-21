# using ImplicitBVH
using GLMakie
using ImplicitBVH: BBox
using StaticArrays
using FileIO
using MeshIO
using GeometryBasics
using BenchmarkTools
#include when plotting for debugging
include("src/plot_box.jl")

#own bvh builder, everything is 1 indexed


struct TreeInfo
    lvl::Int
    n_nodes::Int
    n_leaves::Int
    n_meshelements::Int
end

struct BVH_simple{T}
    mesh::AbstractMesh
    nodes::Vector{BBox{T}}
    leaves::Vector{BBox{T}}
    SubD::Vector{BBox{T}}
    entries::Vector{AbstractMesh}
    info:: TreeInfo
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



function inbox(x::AbstractVector{T}, box::BBox{U}) where {T,U}
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end

function boxtobox(box::HyperRectangle{3,T}) where T
    # 0 allocations
   BBox{T}(box.origin,box.origin .+box.widths)
end


function expand_box(box:: BBox{T}, overlap::Real) where T
    BBox{T}(box.lo .-overlap, box.up .+ overlap);
end

function split_box(box::BBox{T}) where T
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
    return BBox{T}(o,o .+new_width),BBox{T}(o2,o2 .+new_width)
end

function make_boxes(box0::BBox{T},lvl::Int) where{T}
    size=2^lvl -1
    n_leaves=2^(lvl-1)
    n_nodes = size - n_leaves

    # box_array=Vector{BBox{Float32}}
    box_array = Vector{BBox{T}}(undef, size)
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
function split_mesh(box::BBox{T}, mesh::AbstractMesh) where T
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






# bbox3=BBox{Float32}((33.0f0, 3.0f0, 33.0f0), (67.0f0, 102.0f0, 73.0f0))
file="obj/kite_large.obj" # 0 allocs
# mesh = load("obj/"*file*".obj") # 81745 allocations!

bvh=BVH_simple(file,4)
begin
# plots for debugging
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1])
# boxes=vcat(bvh.nodes,bvh.leaves)

box_tocheck=3

wireframe!(ax, bvh.entries[box_tocheck], color = "black", ssao = true)

lines!(ax, boxes_lines(bvh.SubD), linewidth = 2, color = "grey", linestyle=:dash)


lines!(ax, boxes_lines([bvh.SubD[box_tocheck]]), linewidth = 2, color = "red", linestyle=:dash)
lines!(ax, boxes_lines([bvh.leaves[box_tocheck]]), linewidth = 2, color = "red")
lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 5, color = "orange",linestyle=:solid)
# lines!(ax, boxes_lines([bvh.leaves[1]]), linewidth = 5, color = "blue",linestyle=:solid)
# lines!(ax, boxes_lines([bvh.leaves[2]]), linewidth = 5, color = "blue",linestyle=:solid)

# scatter!(ax, [Point3f(x)], color = :red, markersize = 15)
fig
end


save("entrie1.obj", bvh.entries[1])
