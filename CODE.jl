
using StaticArrays
using ImplicitBVH
using GeometryBasics
using ImplicitBVH: BSphere, BBox
using MeshIO
using FileIO
using BenchmarkTools
using GLMakie
using WriteVTK
import WriteVTK: VTKCellData

include("src/plot_box.jl")
include("src/mesh_body_copy.jl")


test_box=Rect([0.0,0.0,0.0],[1.0,1.0,2.0])

# function split_box_allo(box::HyperRectangle{3,T}) where T
#     max_dir = argmax(box.widths)

#     new_width = collect(box.widths)
#     new_width[max_dir] = box.widths[max_dir] / 2  # FIXED here

#     # child one
#     o1 = collect(box.origin)
#     # child two
#     o2 = copy(o1)
#     o2[max_dir] += new_width[max_dir]

#     return Rect(SVector{3,T}(o1), SVector{3,T}(new_width)),
#            Rect(SVector{3,T}(o2), SVector{3,T}(new_width))
# end


function split_box(box::HyperRectangle{3,T}) where T
    #0 allocations
    max_dir=argmax(box.widths)
    w= box.widths
    new_width= SVector(
        max_dir ==1 ? w[1]/2 : w[1],
        max_dir ==2 ? w[2]/2 : w[2],
        max_dir ==3 ? w[3]/2 : w[3]
    )
    if any(new_width .< 1)
        throw(ArgumentError("Width of a box can not be less then 1 cell"))
    end
    o= box.origin
    o2 = SVector(
    max_dir ==1 ? o[1]+new_width[1] : o[1],
    max_dir ==2 ? o[2]+new_width[2] : o[2],
    max_dir ==1 ? o[3]+new_width[3] : o[3],
    )
    return Rect(o,new_width),Rect(o2,new_width)
end




function split_box(box::BBox{T}) where T
    #0 allocations
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
    return BBox(o,o .+new_width),BBox(o2,o2 .+new_width)
end

test_bbox=BBox{Float32}([0.0, 0.0, 0.0],[8.0,8.0,16.0])
# @btime test_splits=split_box($test_bbox)
#@btime split_box_allo($test_box)


function boxtobox(box::HyperRectangle{3,T}) where T
    # 0 allocations
   BBox(box.origin,box.origin .+box.widths)
end

#expands by inc cells in each direction and rounds to resnap to grid
function expand_box(box::HyperRectangle;inc=2)
    #0 allocations
    Rect(round.(box.origin .-inc),round.( box.widths .+ inc))
end

function expand_box(box:: BBox{T}; inc=2) where T
    BBox(round.(box.lo .-inc), round.( box.up .+ inc))
end

# @btime boxtobox($test_box)
# @btime expand_box($test_box)

# function constr_leaves(box::BBox{T}; levels::Int = 3) where T
    
#     current_level = [box]  # Start with just the root box

#     for i in 1:levels - 1
#         next_level = Vector{typeof(box)}()
#         for box in current_level
#             a, b = split_box(box)
#             push!(next_level, a)
#             push!(next_level, b)
#         end
#         current_level = next_level
#     end

#     return current_level  # Final leaf boxes
# end


function constr_leaves2(box::BBox{T}; levels::Int = 4) where T
# 1 alloc per level
    current_level = [box]  # start with root
    for _ in 1:levels-1
        n = 2 * length(current_level)
        next_level = Vector{typeof(box)}(undef, n)
        for i in 1:length(current_level)
             a, b = split_box(current_level[i])
            next_level[2i - 1] = a
            next_level[2i]     = b
        end
        current_level = next_level
    end
    return current_level
end



function collect_imp(box::BBox{T}, mesh::AbstractMesh) where T
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



#lvl leafs allos (leaves*4+1?)
#1   1       6
#2   2       9
#3   4       17
#4   8       33
#5   16      49
#6   32      81 
#7   64      113
#8   128     217
#9   256     369
#10  512     633
#11  1024    1161
#12  2048    2193
#13  4096    4258


function find_box(x::SVector{3,T}, bvh,map::AbstractVector) where T

    function inbox(x::SVector{3,T}, box::BBox{U}) where U
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end

    all_boxes=vcat(bvh.nodes,bvh.leaves) #array of all boxes 1-indexed memory
    virts_per_lvl=map[2]
    index_map=map[1]
    matching_leaves=Int[]
    stack=[0]
    while !isempty(stack)  #this will loop over implicit index
        #all logic is performed using implicit indices 
        #acces is done with mem indices
        i=pop!(stack)   #implicit index
        lv_i=Int(floor(log2(i+1)))
        i_mem=index_map[i+1]
        box= all_boxes[i_mem] #box from memory index
        
        # println("we are checking box ", i," at mem ", i_mem)
        if inbox(x,box)
            if lv_i== (bvh.tree.levels-1) # we are in a leaf
                # println("checking leaf ", i)
                leaf_idx=i_mem-(bvh.tree.real_nodes-bvh.tree.real_leaves)
                push!(matching_leaves,leaf_idx)
            elseif 2^(lv_i+1)-virts_per_lvl[lv_i+2]+2^(lv_i+1)-2 < 2i+2    #+2 to account for 1 index
                # println("hello ",2^(lv_i+1)-virts_per_lvl[lv_i+2]+2^(lv_i+1)-1)
                #the node is incomplete
                # 2^(lv+1) -virts_per_lvl  number of real nodes at next level
                # 2^(lv_i+1)-1  index setoff
                # println("we are at incomplete node ", i)
                # println("pushing ",2i+1)
                push!(stack,2i+1)
            else
                # println("we are at complete node ", i)
                # println("pushing ",2i+1," ",2i+2)
                push!(stack,2i+1,2i+2)
            end
        end
    end
    # println(stack)
    return matching_leaves

end




function inbox(x::SVector{3,T}, box::BBox{U}) where {T,U}
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end

function remove_empty_boxes(boxes, entries)
    filtered_boxes = []
    filtered_entries = []

    for i in eachindex(boxes)
        if !isempty(entries[i])
            push!(filtered_boxes, boxes[i])
            push!(filtered_entries, entries[i])
        end
    end

    return filtered_boxes, filtered_entries
end

function remove_empty_boxes2(boxes::Vector{T}, entries::Vector{Vector{F}}) where {T, F}
    filtered_boxes = Vector{T}()
    filtered_entries = Vector{Vector{F}}()

    for i in eachindex(boxes)
        if !isempty(entries[i])
            push!(filtered_boxes, boxes[i])
            push!(filtered_entries, entries[i])
        end
    end

    @assert eltype(filtered_boxes) === eltype(boxes)
    @assert eltype(filtered_entries) === eltype(entries)

    return filtered_boxes, filtered_entries
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
        sols = [find_dist(sel_mesh, x)[1] for sel_mesh in meshes]
        dist_field[I] = sols[argmin(abs.(sols))]


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

const LeafType = BBox{Float32}
const NodeType = BBox{Float32}

#setup

# μ₀ distance offset halfgrid size
#smooth jump function
# if +-1 calculate mu0


function make_map(bvh)
    l_dash=bvh.tree.levels
    Lv=bvh.tree.virtual_leaves
    function imp_to_mem(i::Int)#from 0 indexed imp to 1 indexed mem 
        l_i=floor(log2(i+1))
        Lvl=Lv >>Int(l_dash-l_i)
        Nvl=2*Lvl-count_ones(Lvl)
        i_mem=i-Nvl+1
    end
    # map=collect(0:bvh.tree.real_nodes-1)
    map=[imp_to_mem(i) for i in collect(0:2^l_dash -1)]
    virt_per_lvl=[Lv>>Int(l_dash-l_i) for l_i in collect(1:bvh.tree.levels)]
    return [map,virt_per_lvl]
end

# total_start=time_ns()

begin
file="cent_dragon" # 0 allocs
mesh = load(joinpath(@__DIR__, "obj/"*file*".obj")) # 81745 allocations!

# begin
bbox = Rect(mesh.position) # 0 allocations
bbox2=boxtobox(bbox) #0 allocations
bbox3=expand_box(bbox2, inc=1) #0 allocations
lvl=6

boundboxs_const=constr_leaves2(bbox3, levels= lvl) #1 allocation per level
println("constructed")
println(typeof(boundboxs_const))
println(eltype(boundboxs_const))
println(length(boundboxs_const))

boundboxs_exp=[expand_box(box, inc=0) for box in boundboxs_const];  #1 alloc
leaf_entries=[collect_imp(box,mesh) for box in boundboxs_exp]  #checked and works

println("expanded")
println(typeof(boundboxs_exp))
println(eltype(boundboxs_exp))
println(length(boundboxs_exp))

boundboxs,leaf_entries_pop=remove_empty_boxes2(boundboxs_exp,leaf_entries);
println("popped")
println(typeof(boundboxs))
println(eltype(boundboxs))
println(length(boundboxs))


bvh= BVH(boundboxs)
imp_map=make_map(bvh);

L=64


dist_field=construct_field(L,bvh,imp_map,mesh,leaf_entries_pop);

println("constructed bvh with $lvl levels")

store_vtk(dist_field,file*"_KD");

fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1])

wireframe!(ax, mesh, color = "black", ssao = true)

# lines!(ax, boxes_lines(bvh.nodes), linewidth = 2, color = "black", linestyle=:dash)
lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "green")
fig
end
# println("dist with KD calculated")



#comparison against previous method

@benchmark begin
    
function measure_test(;L=32,Re=5e5,U=1,mem=Array)
    
    body=MeshBody("obj/"*file*".obj";scale=1,boundary=true)

    Simulation((L,L,L), (U,0,0), L; ν=U*L/Re, body=body, mem)
end

_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
custom_attrib = Dict("d" => _body, )

sim = measure_test(L=64)

measure!(sim)

# wr=vtkWriter(file*"meshbody";attrib=custom_attrib)

# write!(wr, sim)

# close(wr)
# println("dist with waterlilly calculated")

end






begin
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1])
    
    wireframe!(ax, mesh, color = "black", ssao = true)
    
    # lines!(ax, boxes_lines(bvh.nodes), linewidth = 2, color = "black", linestyle=:dash)
    lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "green")
    fig
    end