using GLMakie

include("/home/raj/thesis/bvh/src/plot_box.jl")
using ImplicitBVH




function find_box_scratch(x::SVector{3,T}, bvh,map::AbstractVector) where T

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

function imp_to_mem(i::Int,bvh)#from 0 indexed imp to 1 indexed mem 
l_dash=bvh.tree.levels
Lv=bvh.tree.virtual_leaves
l_i=floor(log2(i+1))
Lvl=Lv >>Int(l_dash-l_i)
Nvl=2*Lvl-count_ones(Lvl)
i_mem=Int(i-Nvl+1)
end


    begin

    # test_box=[BBox([0.0, 0., 0.],[5.,5.,5.]),BBox([0.,5.,0.],[5.,10.,10]),BBox([0.,5.,5.],[5.,10.,15])]

    # bvh=BVH(test_box)
    println("restart")
    x = SVector{3, Float32}(6.,35.,30.)

    found_box=find_box_scratch(x,bvh,imp_map)
    println(found_box)
    all_box=vcat(bvh.nodes, bvh.leaves)

    fig = Figure(size= (1200, 800),SSAO=true)
    ax = Axis3(fig[1, 1])
    # ax.scene.ssao = true  # enable SSAO like in wireframe

    # all_boxes = vcat(bvh.nodes, bvh.leaves)
    # Plot BVH leaf bounding boxes
    # lines!(ax, boxes_lines(bvh.leaves), linewidth = 2, color = :green)
    # lines!(ax, boxes_lines(all_boxes[[1]]), linewidth = 2, color = :black, linestyle=:dash)

    scatter!(ax, [Point3f(x)], markersize = 20, color = :red)
    lines!(ax, boxes_lines([bvh.nodes[1]]), linewidth = 2, color = :green ,linestyle=:solid)
    lines!(ax, boxes_lines([all_box[imp_to_mem(9,bvh)]]), linewidth = 6, color = :orange ,linestyle=:dash)
    lines!(ax, boxes_lines([all_box[imp_to_mem(20,bvh)]]), linewidth = 2, color = :black ,linestyle=:solid)
    # lines!(ax, boxes_lines([all_box[5]]), linewidth = 2, color = :red ,linestyle=:solid)
    # lines!(ax, boxes_lines([all_box[28]]), linewidth = 2, color = :orange ,linestyle=:solid)

    # lines!(ax, boxes_lines(boundboxs), linewidth = 2, color = :green ,linestyle=:solid)

    # lines!(ax, boxes_lines([boundboxs[5]]), linewidth = 4, color = :black,linestyle=:solid)
    # lines!(ax, boxes_lines([all_box[29]]), linewidth = 4, color = :blue,linestyle=:solid)
    # lines!(ax, boxes_lines([all_box[28]]), linewidth = 4, color = :black,linestyle=:solid)

    if typeof(found_box)!=Nothing
        lines!(ax, boxes_lines(bvh.leaves[found_box]), linewidth = 5, color = :red, linestyle=:solid)
    end
    int=10
    # lines!(ax, boxes_lines(bvh.nodes[[10]]), linewidth = 2, color = :red, linestyle=:dash)

    # lines!(ax,boxes_lines(bvh.leaves[[5,6]]), linewidth = 2, color = :red)


    fig

end