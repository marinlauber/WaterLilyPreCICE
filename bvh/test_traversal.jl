

include("src/bvh_builder.jl")

function find_box(x::SVector{3,T}, bvh::BVH_simple) where T

    function inbox(x::SVector{3,T}, box::BoundBox{U}) where U
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end
    all_boxes=vcat(bvh.nodes,bvh.leaves) #array of all boxes 1-indexed memory
    matching_leaves= Vector{Int}(undef,8)
    stack=[1]
    count=0
    while !isempty(stack)  
        i=pop!(stack)   #implicit index
        box= all_boxes[i] #box from memory index
        # println("we are checking box ", i)
        if inbox(x,box)
            if i > bvh.info.n_nodes
                count+=1
                matching_leaves[count]=Int(i-bvh.info.n_nodes)
                # println("we are in a leaf")
                # push!(matching_leaves,i-bvh.info.n_nodes)
            else
                # println(" pushing children ", 2i,"   ",2i+1)
                push!(stack,2i,2i+1)
            end
        end
    end
    # println(stack)
    return matching_leaves[1:count]
end



function find_box_fsm(x::SVector{3,T},bvh::BVH_simple) where T
    
    function inbox_idx(current::Int)
        box=all_boxes[current]
        all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
    end

    sibling(current::Int) = current%2==0 ? current+1 : current-1
    isLeaf(current::Int) = current > bvh.info.n_nodes ? true : false
    parent(current::Int)= Int(floor(current/2))
    all_boxes=vcat(bvh.nodes,bvh.leaves) #array of all boxes 1-indexed memory
    matching_leaves= Vector{Int}(undef,8)
    count=0
    state="fromParent"
    current=1
    while true
        println("checking box ", current)
        println("in state ", state)
        if state=="fromChild"
            if current ==1
                break
            elseif current==2*(parent(current))
                current=sibling(current)
                state="fromSibling"
            else
                current=parent(current)
                state="fromChild"
            end
        continue

        elseif state=="fromSibling"
            hit=inbox_idx(current)
            if !hit
                current=parent(current)
                state="fromChild"
            elseif isLeaf(current)
                count+=1
                println("hit on box ", current-bvh.info.n_nodes)
                println(state)
                matching_leaves[count]=Int(current-bvh.info.n_nodes)
                current=parent(current)
                state="fromChild"
            else
                current=2*current
                state="fromParent"
            end
        continue

        elseif state=="fromParent"
            hit=inbox_idx(current)
            if !hit
                current=sibling(current)
                state="fromSibling"
            elseif isLeaf(current)
                count+=1
                # println(count)
                println("hit on box ", current-bvh.info.n_nodes)
                println(state)
                matching_leaves[count]=Int(current-bvh.info.n_nodes)
                current=sibling(current)
                state="fromSibling"
            else
                current=2*current
                state="fromParent"
            end
        continue

        end


    end
return matching_leaves[1:count]

end

function inbox(x::SVector{3,T}, box::BoundBox{U}) where {T,U}
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end



file="obj/kite_large.obj" # 0 allocs
mesh=load(file)
# mesh = load("obj/"*file*".obj") # 81745 allocations!
bvh=BVH_simple(file,6,overlap=2)


begin
    
    x=SVector(50.0,28.0,62.0)
    # plots for debugging
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1])
    # boxes=vcat(bvh.nodes,bvh.leaves)
    stack_result=find_box(x,bvh)
    fsm_result=find_box_fsm(x,bvh)
    box_tocheck=3
    println("stack ",stack_result)
    println("fsm ", fsm_result)
    # wireframe!(ax, bvh.mesh, color = "black", ssao = true)
    
    # lines!(ax, boxes_lines(bvh.SubD), linewidth = 2, color = "grey", linestyle=:dash)
    
    
    # lines!(ax, boxes_lines([bvh.SubD[box_tocheck]]), linewidth = 2, color = "red", linestyle=:dash)
    lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = "grey")
    # lines!(ax, boxes_lines(bvh.leaves[foundbox]), linewidth = 2, color = "red")
    lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 5, color = "orange",linestyle=:solid)
    # lines!(ax, boxes_lines([bvh.leaves[1]]), linewidth = 5, color = "blue",linestyle=:solid)
    # lines!(ax, boxes_lines([bvh.leaves[2]]), linewidth = 5, color = "blue",linestyle=:solid)
    
    scatter!(ax, [Point3f(x)], color = :red, markersize = 15)
    fig
end
    