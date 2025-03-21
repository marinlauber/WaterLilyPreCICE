

include("src/bvh_builder.jl")

function find_box(x::SVector{3,T}, bvh::BVH_simple) where T

    function inbox(x::SVector{3,T}, box::BoundBox{U}) where U
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end
    all_boxes=vcat(bvh.nodes,bvh.leaves) #array of all boxes 1-indexed memory
    matching_leaves=Int[]
    stack=[1]
    while !isempty(stack)  
        i=pop!(stack)   #implicit index
        box= all_boxes[i] #box from memory index
        println("we are checking box ", i)
        if inbox(x,box)
            if i > bvh.info.n_nodes
                println("we are in a leaf")
                push!(matching_leaves,i-bvh.info.n_nodes)
            else
                println(" pushing children ", 2i,"   ",2i+1)
                push!(stack,2i,2i+1)
            end
        end
    end
    # println(stack)
    return matching_leaves
end

function inbox(x::SVector{3,T}, box::BoundBox{U}) where {T,U}
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end



file="obj/kite_large.obj" # 0 allocs
mesh=load(file)
# mesh = load("obj/"*file*".obj") # 81745 allocations!
bvh=BVH_simple(file,4,overlap=3)


begin
    
    x=SVector(60.0,30.0,62.0)
    # plots for debugging
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1])
    # boxes=vcat(bvh.nodes,bvh.leaves)
    # foundbox=find_box(x,bvh)
    box_tocheck=3
    
    # wireframe!(ax, bvh.mesh, color = "black", ssao = true)
    
    # lines!(ax, boxes_lines(bvh.SubD), linewidth = 2, color = "grey", linestyle=:dash)
    
    
    # lines!(ax, boxes_lines([bvh.SubD[box_tocheck]]), linewidth = 2, color = "red", linestyle=:dash)
    lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "grey")
    lines!(ax, boxes_lines(bvh.leaves[foundbox]), linewidth = 2, color = "red")
    lines!(ax, boxes_lines([bvh.nodes[4]]), linewidth = 5, color = "orange",linestyle=:solid)
    lines!(ax, boxes_lines([bvh.leaves[1]]), linewidth = 5, color = "blue",linestyle=:solid)
    lines!(ax, boxes_lines([bvh.leaves[2]]), linewidth = 5, color = "blue",linestyle=:solid)
    
    scatter!(ax, [Point3f(x)], color = :red, markersize = 15)
    fig
    end
    