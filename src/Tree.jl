


struct BoundingBox{T}
    lo::SVector{3, T}
    up::SVector{3, T}
end
BoundingBox(lo::NTuple{3,T}, up::NTuple{3,T}) where T = BoundingBox{T}(SVector{3,T}(lo), SVector{3,T}(up))

function BoundingBox(mesh_part::AbstractMesh,exp::Real,T::DataType)
    rec=Rect{3,T}(mesh_part.position) #0
    bbox=BoundingBox(rec) #0
    box0=expand_box(bbox,exp)
    return box0
end
BoundingBox(box::HyperRectangle{3,T}) where T = BoundingBox{T}(box.origin,box.origin.+box.widths)

struct TreeInfo
    lvl::Int
    n_nodes::Int
    n_leaves::Int
    n_meshelements::Int
end

struct BoundedVolumeHierarchy{T}
    mesh::AbstractMesh
    nodes::Vector{BoundingBox{T}}
    leaves::Vector{BoundingBox{T}}
    all_nodes::Vector{BoundingBox{T}}
    entries::Vector{AbstractMesh}
    info:: TreeInfo
end

function Base.show(io::IO, bvh::BoundedVolumeHierarchy{T}) where T
    println(io, "BVH ")
    println(io, " ├─ type used:  $T")
    println(io, " ├─ levels:     $(bvh.info.lvl)")
    println(io, " ├─ Nodes:      $(bvh.info.n_nodes)")
    println(io, " ├─ Leaves:     $(bvh.info.n_leaves)")
    println(io, " ├─ mesh size:  $(bvh.info.n_meshelements)")
end


function BoundedVolumeHierarchy(fname::String, lvl::Int; overlap::Int=1, box0_exp::Int=2,T::DataType=Float32, mem=Array)
    mesh = load(fname)
    return BoundedVolumeHierarchy(mesh,lvl, overlap=overlap, box0_exp=box0_exp,T=T, mem=mem)
end

function BoundedVolumeHierarchy(mesh::AbstractMesh, lvl::Int; overlap::Int=1, box0_exp::Int=4, T::DataType=Float32, mem=Array)

    box0 = BoundingBox(mesh, box0_exp, T)
    SubD = make_boxes_old(box0, lvl) |> mem
    positions= mesh.position |> mem

    @loop SubD[I] = expand_box(SubD[I], overlap) over I ∈ CartesianIndices(SubD)

    # Now using mem for arrays
    entries = Vector{AbstractMesh}(undef, length(SubD))      # split mesh parts
    leaves  = mem{BoundingBox{T}}(undef, length(SubD))       # bounding boxes

    @loop entries[I] = split_mesh(SubD[I], mesh) over I ∈ CartesianIndices(SubD)
    @loop leaves[I]   = BoundingBox(entries[I], 2, T) over I ∈ CartesianIndices(SubD)

    nodes = construct_nodes(leaves, lvl)  # make sure this is GPU-safe too
    info = TreeInfo(
        lvl,
        length(nodes),
        length(leaves),
        length(mesh.faces)
        )

    # Use vcat from GPU-compatible libraries (e.g. CUDA.jl has CUDA.vcat) if needed
    return BoundedVolumeHierarchy{T}(mesh, nodes, leaves, vcat(nodes, leaves), entries, info)
end

function construct_nodes(leaves::Vector{BoundingBox{T}}, lvl::Int) where T
    n_nodes =  2^(lvl-1)-1
    all_boxes = vcat(Vector{BoundingBox{T}}(undef, n_nodes),copy(leaves))  # working array that will grow upward)

    # Fill in parent nodes from leaves upward
    for i in n_nodes:-1:1
        left = 2i
        right = 2i + 1
        parent= merge_boxes(all_boxes[left], all_boxes[right])
        all_boxes[i] = parent
        # push!(all_boxes, nodes[i])  # extend the working array
    end

    return all_boxes[1:n_nodes]
end

function merge_boxes(box1::BoundingBox{T}, box2::BoundingBox{T}) where T
    #0 alloc
    is_valid(b) = all(isfinite, b.lo) && all(isfinite, b.up)

    if is_valid(box1) && !is_valid(box2)
        return box1
    elseif !is_valid(box1) && is_valid(box2)
        return box2
    elseif !is_valid(box1) && !is_valid(box2)
        # Both invalid → return a "NaN" box
        return BoundingBox{T}(SVector(NaN, NaN, NaN), SVector(NaN, NaN, NaN))
    end

    lo = SVector(
        min(min(box1.lo[1], box1.up[1]), min(box2.lo[1], box2.up[1])),
        min(min(box1.lo[2], box1.up[2]), min(box2.lo[2], box2.up[2])),
        min(min(box1.lo[3], box1.up[3]), min(box2.lo[3], box2.up[3]))
    )
    up = SVector(
        max(max(box1.lo[1], box1.up[1]), max(box2.lo[1], box2.up[1])),
        max(max(box1.lo[2], box1.up[2]), max(box2.lo[2], box2.up[2])),
        max(max(box1.lo[3], box1.up[3]), max(box2.lo[3], box2.up[3]))
    )

    return BoundingBox{T}(lo, up)
end



function inbox(x::AbstractVector{T}, box::BoundingBox{U}) where {T,U}
    all((box.lo .<= x) .& (x .<= box.up)) || all((box.up .<= x) .& (x .<= box.lo))
end


function expand_box(box:: BoundingBox{T}, overlap::Real) where T
    BoundingBox{T}(box.lo .-overlap, box.up .+ overlap);
end


#only 7 ns slower then function with GeometryBasics 
function split_box(box::BoundingBox{T}) where T  #15 ns
    # println(box)
    #0 allocations
    w= box.up .- box.lo
    max_dir=argmax(w)
    new_width= SVector(
        max_dir ==1 ? w[1]/2 : w[1],
        max_dir ==2 ? w[2]/2 : w[2],
        max_dir ==3 ? w[3]/2 : w[3]
    )
    # println(new_width)
    if any(new_width .< 1)
        throw(ArgumentError("Width of a box can not be less then 1 cell"))
    end
    o= SVector(box.lo)
    o2 = SVector(
    max_dir ==1 ? o[1]+new_width[1] : o[1],
    max_dir ==2 ? o[2]+new_width[2] : o[2],
    max_dir ==3 ? o[3]+new_width[3] : o[3],
    )
    return BoundingBox{T}(o,o .+new_width),BoundingBox{T}(o2,o2 .+new_width)
end

function make_boxes_old(box0::BoundingBox{T},lvl::Int) where{T}
    size=2^lvl -1
    n_leaves=2^(lvl-1)
    n_nodes = size - n_leaves

    box_array = Vector{BoundingBox{T}}(undef, 2^lvl -1)
    box_array[1]=box0

    for i in 1:2^(lvl-1)-1
        # println(i)
        child1,child2=split_box(box_array[i])
        box_array[2i]=child1
        box_array[2i+1]=child2
    end
    SubDivisions = box_array[n_nodes+1:end];
    return SubDivisions
end


function make_boxes(lvl::Int, box_array::Vector{BoundingBox{T}}) where{T}

    for i in 1:2^(lvl-1)-1
        # println(i)
        child1,child2=split_box(box_array[i])
        box_array[2i]=child1
        box_array[2i+1]=child2
    end
    # SubDivisions = box_array[n_nodes+1:end];
    return box_array
end

function split_mesh(box::BoundingBox{T}, mesh::AbstractMesh) where T
    faces = GeometryBasics.faces(mesh)
    max=length(mesh.faces)
    mask=fill(NaN32,length(mesh.position))
    count=0
    coll = Vector{eltype(faces)}(undef,max)
    @inbounds for face in faces
        for i in Tuple(face)
            v = mesh.position[i] 
            inside = all(box.lo .<= v .<= box.up)

            if inside
                count+=1
                coll[count]=face
                mask[face].=1
                break
            end
        end
    end
    # new_mesh=GeometryBasics.Mesh(mesh.position.*mask,coll[1:count],normal=mesh.normal.*mask)
    new_mesh=GeometryBasics.Mesh(mesh.position.*mask,coll[1:count])
    return new_mesh
end

