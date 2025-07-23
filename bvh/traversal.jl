
using StaticArrays
using GeometryBasics
using MeshIO
using FileIO
using BenchmarkTools
using GLMakie
using WriteVTK
import WriteVTK: VTKCellData

include("src/plot_box.jl")
include("src/bvh_builder.jl")


file="obj/kite_large.obj" # 0 allocs
mesh=load(file)
# mesh = load("obj/"*file*".obj") # 81745 allocations!
bvh=BVH_simple(file,4,overlap=3)

function find_box(x::SVector{3,T}, bvh::BVH_simple) where T

    function inbox(x::SVector{3,T}, box::BoundBox{U}) where U
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end
    all_boxes=vcat(bvh.nodes,bvh.leaves) #array of all boxes 1-indexed memory
    matching_leaves=Int[]
    stack=[0]
    while !isempty(stack)  
        i=pop!(stack)   #implicit index
        box= all_boxes[i] #box from memory index
        # println("we are checking box ", i," at mem ", i_mem)
        if inbox(x,box)
            if i > bvh.info.n_nodes
                push!(matching_leaves,i)
            else
                push!(matching_leaves,2i,2i+1)
            end
        end
    end
    # println(stack)
    return matching_leaves
end

function find_dist(mesh::M,x::SVector{T};kwargs...) where {M<:GeometryBasics.Mesh,T}
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

function construct_field(L::Int, bvh::BVH{U},map::AbstractVector, mesh::AbstractMesh, leaf_entries) where {U}
    dist_field=fill(25.0,L,L,L)
    for I in CartesianIndices(dist_field)
        x = SVector{3, Int}(Tuple(I)) 

        leaves_stack = find_box(x, bvh,map)
        if isempty(leaves_stack)
            # println("tested coord and skipped", x, " in ", round(elapsed / 1e6, digits=3), " ms")
            dist_field[I] = NaN
            continue
        end

        leaves = leaf_entries[leaves_stack]
        meshes = [GeometryBasics.Mesh(mesh.position, leaf) for leaf in leaves]
        try
            sols = [find_dist(sel_mesh, x)[1] for sel_mesh in meshes]
            dist_field[I] = sols[argmin(abs.(sols))]
        catch
            dist_field[I] = NaN
        end
        


        # println("tested coord ", x, " in ", round(elapsed / 1e6, digits=3), " ms")
    end
    # println("Average time per tested coord: ", round(avg_ms, digits=3), " ms")
    # println("total time for coords ", round(total_time_ns,digits=3), " ms")

    return dist_field
end
# @inside dist_field[I]=function(I,kwargs)

# use bounding box zero to start looking  @loop from waterlilly


function store_vtk(dist_field::Array{T,3}, filename::String = "dist_field_output") where T
    nx, ny, nz = size(dist_field)
    x = 1:nx
    y = 1:ny
    z = 1:nz

    vtk_grid(filename, x, y, z) do vtk
        vtk["dist"] = dist_field  # The name used in ParaView
    end

    # println("Saved VTK file: $filename.vtr")
end

begin

    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1])

    wireframe!(ax, mesh, color = "black", ssao = true)

    # lines!(ax, boxes_lines(bvh.nodes), linewidth = 2, color = "black", linestyle=:dash)
    lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "green")
    fig
end