using WaterLily
using BenchmarkTools
using FileIO,MeshIO
using GeometryBasics
include("src/mesh_body.jl")
include("src/bvh_builder.jl")
include("src/util.jl")
include("src/plot_box.jl")


file="/home/marin/Workspace/WaterLilyPreCICE_KD/bvh/obj/dragon_15k.obj" # 0 allocs
mesh=load(file)
# mesh = load("obj/"*file*".obj") # 81745 allocations!
@benchmark bvh=BVH_simple($mesh,4,overlap=2)
bvh=BVH_simple(mesh,6,overlap=2)

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
end # 120.029 ns (0 allocations: 0 bytes)d # 4.266 μs (0 allocations: 0 bytes)


@benchmark s_t=find_dist($bvh.entries[8],$SVector(40.0,30.0,20.0))

@benchmark stt=find_dist($mesh,$SVector(40.0,30.0,20.0))


function find_box(x::CartesianIndex{3}, bvh::BVH_simple) where T
    # println(x)
    x=SVector{3,Int}(Tuple(x))
    # println(x)

    function inbox(x, box::BoundBox{U}) where U
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end
     #array of all boxes 1-indexed memory
    # matching_leaves= Vector{Int}(undef,15)
    stack=[1]
    count=0
    sol=25
    while !isempty(stack)
        i=pop!(stack)   #implicit index
        box= bvh.all_nodes[i] #box from memory index
        # println("we are checking box ", i)
        if inbox(x,box)
            if i > bvh.info.n_nodes
                count+=1
                # matching_leaves[count]=Int(i-bvh.info.n_nodes)
                sol_leaf = find_dist(bvh.entries[Int(i-bvh.info.n_nodes)],x)[1]
                abs(sol_leaf)<abs(sol) ? sol=sol_leaf : sol=sol
                # println("found leaf ", i-bvh.info.n_nodes )
                # println("with sol ", sol_leaf)
                # println("new_sol ", sol)
                # push!(matching_leaves,i-bvh.info.n_nodes)
            else
                # println(" pushing children ", 2i,"   ",2i+1)
                push!(stack,2i,2i+1)
            end
        end
    end
    # println(stack)
    # println("finished with indx ", x)
   return sol
#    return sol
end

s_t=find_dist(bvh.entries[8],SVector(40.0,30.0,20.0))


function find_box_fsm(x::CartesianIndex{3},bvh::BVH_simple) where T
    fromParent =1
    fromSibling=2
    fromChild=3
    x=SVector{3,Int}(Tuple(x))

    n_nodes= bvh.info.n_nodes
    function inbox_idx(current::Int)
        box=bvh.all_nodes[current]
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end

    sibling(current::Int) = current%2==0 ? current+1 : current-1  #0 alloc
    isLeaf(current::Int) = current > n_nodes   # 0 alloc?
    parent(current::Int)= fld(current,2)  # 0 alloc
    # matching_leaves= Vector{Int}(undef,15)
    # count=0
    sol=25
    state=fromParent
    current=1
    while true
        # println("checking box ", current)
        # println("in state ", state)
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
                # count+=1
                # println("hit on box ", current-bvh.info.n_nodes)
                # println(state)
                # matching_leaves[count]=Int(current-bvh.info.n_nodes)
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
                # count+=1
                # println(count)
                # println("hit on box ", current-bvh.info.n_nodes)
                # println(state)
                # matching_leaves[count]=Int(current-bvh.info.n_nodes)
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
    # return matching_leaves[1:count],sol
    return sol
end

function inbox(x::SVector{3,T}, box::BoundBox) where T
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end

function full_mesh_dist(x::CartesianIndex{3},mesh,box)
    x=SVector{3,Int}(Tuple(x))
    # println(x)

    # function inbox(x, box::BoundBox{U}) where U
    #     all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    # end

    if inbox(x,box)
        sol=find_dist(mesh,x)[1]
    else
        sol=25
    end

    return sol
end


point=CartesianIndex((160,50,85))
box = Bounding_BBox(mesh, 4, Float64)

@benchmark find_box_fsm($point,$bvh)
@benchmark find_box($point,$bvh)
@benchmark full_mesh_dist($point,$mesh,$box)
find_box_fsm(point,bvh)
find_box(point,bvh)
full_mesh_dist(point,mesh,box)
#correct answer checking

# file="obj/dragon_15k.obj" # 0 allocs
# mesh=load(file)
# # mesh = load("obj/"*file*".obj") # 81745 allocations!
# @benchmark bvh=BVH_simple($mesh,4,overlap=2)
# bvh=BVH_simple(mesh,6,overlap=2)

# 1 lvl 62 allocs 146 μs
# 4 lvl 249 allocs 689 μs
# 6 lvls 873 allocs 2.3 ms


L=210
dist_field_bvh=fill(25.0,L,L,L);
dist_field_fsm=fill(25.0,L,L,L);
dist_field_mesh=fill(25.0,L,L,L);

begin #box for
    box=Rect(mesh)
    box=boxtobox(box)
    box=expand_box(box,2)
end

function make_and_fill_dist(mesh,dist_field_bvh)
    bvh=BVH_simple(mesh,6,overlap=2)
    @inside dist_field_bvh[I]=find_box(I,bvh)
end

#@TODO very bad
@benchmark make_and_fill_dist($mesh,$dist_field_bvh)

@inside dist_field_bvh[I]=find_box(I,bvh); println("bvh done")
@inside dist_field_fsm[I]= find_box_fsm(I,bvh); println("fsm done")
@inside dist_field_mesh[I]=full_mesh_dist(I,mesh,box); println("mesh done")

println("benchmark bvh")
@benchmark begin
    bvh=BVH_simple($mesh,6,overlap=2)
    @inside dist_field_bvh[I]=find_box(I,bvh)
end

println("benchmark fsm")
@benchmark begin
    bvh=BVH_simple($mesh,6,overlap=2)
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
    lines!(ax, boxes_lines([bvh.leaves[4]]), linewidth = 5, color = "blue",linestyle=:solid)
    # lines!(ax, boxes_lines([bvh.leaves[2]]), linewidth = 5, color = "blue",linestyle=:solid)
    x=SVector{3,Int}(Tuple(x))
    scatter!(ax, [Point3f(x)], color = :green, markersize = 15)
    fig
end