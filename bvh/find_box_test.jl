

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
        
        println("we are checking box ", i," at mem ", i_mem)
        if inbox(x,box)
            if lv_i== (bvh.tree.levels-1) # we are in a leaf
                println("checking leaf ", i)
                leaf_idx=i_mem-(bvh.tree.real_nodes-bvh.tree.real_leaves)
                push!(matching_leaves,leaf_idx)
            elseif 2^(lv_i+1)-virts_per_lvl[lv_i+2]+2^(lv_i+1)-2 < 2i+2    #+2 to account for 1 index
                # println("hello ",2^(lv_i+1)-virts_per_lvl[lv_i+2]+2^(lv_i+1)-1)
                #the node is incomplete
                # 2^(lv+1) -virts_per_lvl  number of real nodes at next level
                # 2^(lv_i+1)-1  index setoff
                println("we are at incomplete node ", i)
                println("pushing ",2i+1)
                push!(stack,2i+1)
            else
                println("we are at complete node ", i)
                println("pushing ",2i+1," ",2i+2)
                push!(stack,2i+1,2i+2)
            end
        end
    end
    # println(stack)
    return matching_leaves

end

function imp_to_mem(i)
    global bvh
    l_dash=bvh.tree.levels
    Lv=bvh.tree.virtual_leaves
    l_i=floor(log2(i+1))
    Lvl=Lv >>Int(l_dash-l_i)
    Nvl=2*Lvl-count_ones(Lvl)
    i_mem=i-Nvl+1

end


function inbox(x::SVector{3,T}, box::BBox{U}) where {T,U}
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end
inbox(x,bvh.nodes[3])

begin
    
x=SVector(40,70,70)
found_box=find_box(x,bvh,imp_map)
println(found_box)
end
begin
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1])
boxes=vcat(bvh.nodes,bvh.leaves)

# wireframe!(ax, mesh, color = "black", ssao = true)

# lines!(ax, boxes_lines(bvh.nodes), linewidth = 2, color = "black", linestyle=:dash)
lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 2, color = "gray")
lines!(ax, boxes_lines([bvh.nodes[4]]), linewidth = 5, color = "orange",linestyle=:solid)
# lines!(ax, boxes_lines([bvh.leaves[1]]), linewidth = 5, color = "blue",linestyle=:solid)
lines!(ax, boxes_lines([bvh.leaves[2]]), linewidth = 5, color = "blue",linestyle=:solid)

# scatter!(ax, [Point3f(x)], color = :red, markersize = 15)


fig
end
# println("dist with KD calculated")

