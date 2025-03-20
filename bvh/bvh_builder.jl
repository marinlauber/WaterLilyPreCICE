using ImplicitBVH
using GLMakie
using ImplicitBVH: BBox
using StaticArrays
using FileIO
using MeshIO
using GeometryBasics
using BenchmarkTools
include("src/plot_box.jl")

#own bvh builder, everything is 1 indexed

struct BVH_simple{T}
    mesh::AbstractMesh
    nodes::Vector{BBox{T}}
    leaves::Vector{BBox{T}}
    entries::Vector{AbstractMesh}
    lvl:: Int
end

function BVH_simple( file ::String, lvl ::Int ; overlap::Int =0, T::DataType = Float64) 
    mesh=load(file)
    box0=initial_bbox(mesh,T)
    nodes,leaves=make_boxes(box0,lvl)
    leaves=[expand_box(leaf,overlap) for leaf in leaves];
    entries=[GeometryBasics.Mesh(mesh.position,split_mesh(box,mesh)) for box in leaves];
    
    return BVH_simple{T}(mesh,nodes, leaves, entries,lvl)
end

function initial_bbox(mesh::AbstractMesh,T)
    rec=Rect{3,T}(mesh.position)
    bbox=boxtobox(rec)
    box0=expand_box(bbox,3)
    return box0
end

function boxtobox(box::HyperRectangle{3,T}) where T
    # 0 allocations
   BBox{T}(box.origin,box.origin .+box.widths)
end

# bbox3=BBox{Float32}((33.0f0, 3.0f0, 33.0f0), (67.0f0, 102.0f0, 73.0f0))
file="obj/kite_large.obj" # 0 allocs
# mesh = load("obj/"*file*".obj") # 81745 allocations!


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
    #always 4 allocations!
    faces = GeometryBasics.faces(mesh)
    points = mesh.position
    coll = Vector{eltype(faces)}()#undef,6)
    # println("start  ",coll)
    for face in faces
        verts = points[face]
        for vert in verts
            # println("checking vertex: ", vert)
            # println("against box: lo = ", box.lo, ", up = ", box.up)

            if all((box.lo .<= vert) .& (vert .<= box.up)) || all((box.up .<= vert) .& (vert .<= box.lo))
                # println("→ HIT")
                push!(coll, face)
                # coll[i]=face
                # println(coll)
                
                break
            end
        end    
    end
   # coll= [coll[i] for i in eachindex(coll) if isassigned(coll, i)]
    return coll
end

function expand_box(box:: BBox{T}, overlap::Int) where T
    BBox{T}(box.lo .-overlap, box.up .+ overlap);
end


begin

fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1])
# boxes=vcat(bvh.nodes,bvh.leaves)

wireframe!(ax, bvh.mesh, color = "black", ssao = true)

lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "black", linestyle=:dash)
# lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 2, color = "gray")
# lines!(ax, boxes_lines([bvh.nodes[2]]), linewidth = 5, color = "orange",linestyle=:solid)
# lines!(ax, boxes_lines([bvh.leaves[1]]), linewidth = 5, color = "blue",linestyle=:solid)
# lines!(ax, boxes_lines([bvh.leaves[2]]), linewidth = 5, color = "blue",linestyle=:solid)

# scatter!(ax, [Point3f(x)], color = :red, markersize = 15)
fig
end