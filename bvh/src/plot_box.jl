# using StaticArrays
# using ImplicitBVH
# using GeometryBasics
# using ImplicitBVH: BSphere, BBox
# using MeshIO
# using FileIO
# using BenchmarkTools
# using GLMakie
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
