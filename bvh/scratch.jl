using StaticArrays
using MeshIO
using FileIO
using GeometryBasics
using WaterLily: AbstractBody,@loop,measure
using ForwardDiff
using ImplicitBVH
using GLMakie
include("src/mesh_body_copy.jl")
include("src/plot_box.jl")
mesh = load(joinpath(@__DIR__, "obj/sphere.obj"))

bbox = Rect(mesh.position) # 0 allocations
bbox2=boxtobox(bbox) #0 allocations
bbox3=expand_box(bbox2, inc=3) #0 allocations

lvl=4

boundboxs=constr_leaves2(bbox3, levels= lvl)


@benchmark leaf_entries=[collect_imp(box,mesh) for box in boundboxs]

@benchmark leaf_entries2=[collect_imp2(box,mesh) for box in boundboxs]

# Plot a wireframe of the mesh and the bounding boxes above leaf level
function collect_imp(box::BBox{T}, mesh::AbstractMesh) where T
    faces = GeometryBasics.faces(mesh)
    points = mesh.position
    coll = Vector{eltype(faces)}()#undef,6)
    # println("start  ",coll)
    for face in faces
        verts = points[face]
        for vert in verts
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


function collect_imp2(box::BBox{T}, mesh::AbstractMesh) where T
    faces = GeometryBasics.faces(mesh)
    points = mesh.position

    coll = Vector{eltype(faces)}()

    for face in faces
        face_hit = false
        for i in Tuple(face)  # get the indices directly (no allocation)
            v = points[i]

            # Manual componentwise bounds check (avoids allocation)
            inside = (
                (box.lo[1] <= v[1] <= box.up[1] || box.up[1] <= v[1] <= box.lo[1]) &&
                (box.lo[2] <= v[2] <= box.up[2] || box.up[2] <= v[2] <= box.lo[2]) &&
                (box.lo[3] <= v[3] <= box.up[3] || box.up[3] <= v[3] <= box.lo[3])
            )

            if inside
                face_hit = true
                break
            end
        end
        if face_hit
            push!(coll, face)
        end
    end

    return coll
end


mesh_sel=GeometryBasics.Mesh(mesh.position, leaf_entries2[1])

begin
    fig, ax = wireframe(
        mesh_sel;
        color = :green,          # soft gray
        transparency = true,
        overdraw = true,          # helps with visibility of internal edges
        ssao = true
    )
    
    # scatter!(ax, mesh.position, markersize = 10, color = :red)
    
    lines!(ax, boxes_lines(boundboxs), linewidth=2)
    # lines!(ax, boxes_lines(bvh.nodes), linewidth=2)
    fig
    
    end
    