
using ImplicitBVH
using ImplicitBVH: BSphere, BBox

using MeshIO
using FileIO

using BenchmarkTools
using Profile
using PProf

using GLMakie
using GeometryBasics



"""
# Given 5 geometric elements (e.g. bounding boxes) we construct the following implicit tree
# having the 5 real leaves at implicit indices 8-12 plus 3 virtual leaves.
#         Nodes & Leaves                Tree Level
#               1                       1
#       2               3               2
#   4       5       6        7v         3
# 8   9   10 11   12 13v  14v  15v      4

ImplicitTree{Int64}
  levels: Int64 4
  real_leaves: Int64 5
  real_nodes: Int64 11
  virtual_leaves: Int64 3
  virtual_nodes: Int64 4
"""
# BBox{Float32}((-17.971245f0, -17.85344f0, -10.658194f0), (18.097673f0, 17.270058f0, 25.962446f0))
# Types used
const LeafType = BSphere{Float32}
const NodeType = BBox{Float32}
const MortonType = UInt32

mesh = load(joinpath(@__DIR__, "cone_test.obj"))

box_0 = Rect(mesh.position)
box_0= Rect(box_0.origin,box_0.widths.+3)

function constr_leaves(mesh::AbstractMesh, levels :: Int )
    box_0 = Rect(mesh.position)
    box_0= Rect(box_0.origin,bbox.widths.+3)
    leaves=[]
    splits=levels-1
    #function should output array of
    # Vector{BSphere{Float32}} (alias for Array{BSphere{Float32}, 1})
    for i in 1:splits

        #at each split split mesh in half
        println("Iteration $i")
    end


end
#box1= BBox(bbox_geom)
mutable struct populated_leaf
    node_idx ::Int
end


#rewrite geom type box to bounding box
function boxtobox(box::HyperRectangle{3,T}) where T
    up=box.origin .+box.widths
    box=BBox(box.origin,up)

end


bvh = BVH(bounding_spheres, NodeType, MortonType)


a=[0.0 ,0.0 ,0.0]
b=[1.0, 1.0 ,1.0]
test_box=BBox(a,b)

function box_lines!(lines, lo, up)
    # Write lines forming an axis-aligned box from lo to up
    @assert ndims(lines) == 2
    @assert size(lines) == (24, 3)

    lines[1:24, 1:3] .= [
        # Bottom sides
        lo[1] lo[2] lo[3]
        up[1] lo[2] lo[3]
        up[1] up[2] lo[3]
        lo[1] up[2] lo[3]
        lo[1] lo[2] lo[3]
        NaN NaN NaN

        # Vertical sides
        lo[1] lo[2] lo[3]
        lo[1] lo[2] up[3]
        NaN NaN NaN

        up[1] lo[2] lo[3]
        up[1] lo[2] up[3]
        NaN NaN NaN

        up[1] up[2] lo[3]
        up[1] up[2] up[3]
        NaN NaN NaN

        lo[1] up[2] lo[3]
        lo[1] up[2] up[3]
        NaN NaN NaN

        # Top sides
        lo[1] lo[2] up[3]
        up[1] lo[2] up[3]
        up[1] up[2] up[3]
        lo[1] up[2] up[3]
        lo[1] lo[2] up[3]
        NaN NaN NaN
    ]

    nothing
end


function boxes_lines(boxes)
    # Create contiguous matrix of lines representing boxes
    lines = Matrix{Float64}(undef, 24 * length(boxes), 3)
    for i in axes(boxes, 1)
        box_lines!(view(lines, 24 * (i - 1) + 1:24i, 1:3), boxes[i].lo, boxes[i].up)
    end
    lines
end


# lini=boxes_lines(bvh.nodes)

# Plot a wireframe of the mesh and the bounding boxes above leaf level
fig, ax = wireframe(
    mesh,
    color = [tri[1][2] for tri in mesh for i in 1:3],
    colormap=:Spectral,
    ssao=true,
)
lines!(ax, boxes_lines(test_box), linewidth=2)
fig