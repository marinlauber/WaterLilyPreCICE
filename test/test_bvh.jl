using WaterLilyPreCICE
using BenchmarkTools
using WriteVTK
using GeometryBasics
using StaticArrays
using FileIO

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

function store_vtk(dist_field::Array{T,3}, filename::String = "dist_field_output") where T
    nx, ny, nz = size(dist_field)
    x = 1:nx
    y = 1:ny
    z = 1:nz

    vtk_grid(filename, x, y, z) do vtk
        vtk["dist"] = dist_field  # The name used in ParaView
    end
    println("Saved VTK file: $filename.vtr")
end

file = "/home/marin/Workspace/WaterLilyPreCICE_KD/bvh/obj/dragon_15k.obj" # 0 allocs
mesh = load(file)
lvl = 4; overlap=2
@benchmark bvh = BoundedVolumeHierarchy($mesh,$lvl,$overlap)
lvl = 6
bvh = BoundedVolumeHierarchy(mesh,lvl,overlap)

import WaterLilyPreCICE: d²_fast,normal,locate
function find_dist(mesh::M,x::SVector{T};kwargs...) where {M<:GeometryBasics.Mesh,T}
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    n,p = normal(mesh[u]),SVector(locate(mesh[u],x))
    s = x-p # signed Euclidian distance
    d = sign(sum(s.*n))*√sum(abs2,s)
    return d,n,zero(x)
end

x = SVector(40.0,30.0,20.0)
mesh_8 = bvh.entries[8] # mesh part 8
@benchmark s_t = find_dist($mesh_8,$x)
@benchmark stt = find_dist($mesh,$x)

function find_box(x::CartesianIndex{3}, bvh::BoundedVolumeHierarchy)
    x = SVector{3,Int}(Tuple(x))
    stack = [1]
    count = 0
    sol = 25
    while !isempty(stack)
        i = pop!(stack)   #implicit index
        box = bvh.all_nodes[i] #box from memory index
        if inbox(x,box)
            if i > bvh.info.n_nodes
                count += 1
                sol_leaf = find_dist(bvh.entries[Int(i-bvh.info.n_nodes)],x)[1]
                abs(sol_leaf)<abs(sol) ? sol=sol_leaf : sol=sol
            else
                push!(stack,2i,2i+1)
            end
        end
    end
   return sol
end

s_t = find_dist(bvh.entries[8],SVector(40.0,30.0,20.0))


function find_box_fsm(x::CartesianIndex{3},bvh::BoundedVolumeHierarchy)
    fromParent = 1
    fromSibling = 2
    fromChild = 3
    x = SVector{3,Int}(Tuple(x))

    n_nodes= bvh.info.n_nodes

    function inbox_idx(current::Int)
        box=bvh.all_nodes[current]
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end

    sibling(current::Int) = current%2==0 ? current+1 : current-1  #0 alloc
    isLeaf(current::Int) = current > n_nodes   # 0 alloc?
    parent(current::Int)= fld(current,2)  # 0 alloc

    sol=25
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
                sol_leaf = find_dist(bvh.entries[Int(current-bvh.info.n_nodes)],x)[1]
                abs(sol_leaf)<abs(sol) ? sol=sol_leaf : sol=sol
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
                sol_leaf = find_dist(bvh.entries[Int(current-bvh.info.n_nodes)],x)[1]
                abs(sol_leaf)<abs(sol) ? sol=sol_leaf : sol=sol
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

function inbox(x::SVector{3,T}, box::BoundingBox) where T
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end

function full_mesh_dist(x::CartesianIndex{3},mesh,box)
    x=SVector{3,Int}(Tuple(x))

    if inbox(x,box)
        sol=find_dist(mesh,x)[1]
    else
        sol=25
    end

    return sol
end


point = CartesianIndex((160,50,85))
box = BoundingBox(mesh, 4, Float64)

@benchmark find_box_fsm($point,$bvh)
@benchmark find_box($point,$bvh)
@benchmark full_mesh_dist($point,$mesh,$box)
find_box_fsm(point,bvh)
find_box(point,bvh)
full_mesh_dist(point,mesh,box)
#correct answer checking

L=210
dist_field_bvh=fill(25.0,L,L,L);
dist_field_fsm=fill(25.0,L,L,L);
dist_field_mesh=fill(25.0,L,L,L);

import WaterLilyPreCICE: expand_box
begin #box for
    box=Rect(mesh)
    box=BoundingBox(box)
    box=expand_box(box,2)
end

function make_and_fill_dist(mesh,dist_field_bvh)
    bvh = BoundedVolumeHierarchy(mesh,6,overlap=2)
    @inside dist_field_bvh[I]=find_box(I,bvh)
end

#@TODO very bad
@benchmark make_and_fill_dist($mesh,$dist_field_bvh)

@inside dist_field_bvh[I]=find_box(I,bvh); println("bvh done")
@inside dist_field_fsm[I]= find_box_fsm(I,bvh); println("fsm done")
@inside dist_field_mesh[I]=full_mesh_dist(I,mesh,box); println("mesh done")

println("benchmark bvh")
@benchmark begin
    bvh = BoundedVolumeHierarchy($mesh,6,overlap=2)
    @inside dist_field_bvh[I]=find_box(I,bvh)
end

println("benchmark fsm")
@benchmark begin
    bvh = BoundedVolumeHierarchy($mesh,6,overlap=2)
    @inside dist_field_fsm[I]=find_box_fsm(I,bvh)
end

println("benchmark mesh")
@benchmark begin
    box=Rect($mesh)
    box=boxtobox(box)
    box=expand_box(box,2)
    @inside dist_field_mesh[I]=full_mesh_dist(I,$mesh,box)
end

store_vtk(dist_field_bvh,"dragon_bvh")
store_vtk(dist_field_fsm,"dragon_fsm")
store_vtk(dist_field_mesh,"dragon_mesh")


using GLMakie
begin

    # x=SVector(50.0,20.0,62.0)
    x=CartesianIndex((160,50,85))
    # plots for debugging
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1])
    # # boxes=vcat(bvh.nodes,bvh.leaves)
    stack_boxes=find_box(x,bvh)

    # fsm_boxes,fsm_sol=find_box_fsm(x,bvh)

    # fsm2_boxes,fsm2_sol=find_box_fsm2(x,bvh)
    # box_tocheck=3
    println("stack ",stack_boxes)
    # println("fsm ", fsm_boxes,fsm_sol)

    # println("fsm2 ", fsm2_boxes,fsm2_sol)
    wireframe!(ax, bvh.mesh, color = "black", ssao = true)

    # lines!(ax, boxes_lines(bvh.SubD), linewidth = 2, color = "grey", linestyle=:dash)


    # lines!(ax, boxes_lines([bvh.SubD[box_tocheck]]), linewidth = 2, color = "red", linestyle=:dash)
    # lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "grey")
    lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "red")
    lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 1, color = "grey",linestyle=:solid)
    lines!(ax, boxes_lines([bvh.leaves[6]]), linewidth = 5, color = "blue",linestyle=:solid)
    # lines!(ax, boxes_lines([bvh.leaves[2]]), linewidth = 5, color = "blue",linestyle=:solid)
    x=SVector{3,Int}(Tuple(x))
    scatter!(ax, [Point3f(x)], color = :green, markersize = 15)
    fig
end