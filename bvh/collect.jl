using GeometryBasics
using ImplicitBVH
using GLMakie
using BenchmarkTools
include("src/plot_box.jl")
box=BBox([15.0,-4,22.0],[-15.0, -15.0,0])

mesh = load(joinpath(@__DIR__, "cone_10.obj"))


faces = GeometryBasics.faces(mesh)
points = mesh.position

function collect(box, mesh)
    

    faces = GeometryBasics.faces(mesh)
    points = mesh.position

    coll = []
    for face in faces
        verts = points[face]
        for vert in verts
            # println("checking vertex: ", vert)
            # println("against box: lo = ", box.lo, ", up = ", box.up)

            if all((box.lo .<= vert) .& (vert .<= box.up)) || all((box.up .<= vert) .& (vert .<= box.lo))
               # println("→ HIT")
                push!(coll, face)
                break
            end
        end    
    end
    return coll
end


function collect_imp(box, mesh)

    #always for allocations!
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


@benchmark collect($box,$mesh)
@benchmark collect_imp($box, $mesh)
# collect_imp(box,mesh)
collect(box, mesh)
collect_imp(box,mesh)
fig, ax = wireframe(
    mesh,
    # color = [tri[1][2] for tri in mesh for i in 1:3],
    color="black",
    # colormap=:Spectral,
    ssao=true,
)

lines!(ax, boxes_lines([box]), linewidth=2,color="green")
fig