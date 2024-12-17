using WaterLily
using WaterLilyPreCICE
using GeometryBasics
# using GLMakie
using Interpolations

# GLMakie.activate!(inline=false)
# # create the mesh body
MeshBody = WaterLilyPreCICE.MeshBody
body = MeshBody("/home/marin/Workspace/HHH/code/LIMO_heart/LIMO_4_contact/geom_restart.inp";boundary=false,thk=2,scale=1.f00)
# mesh(body.mesh.position, GLTriangleFace.(GeometryBasics.faces(body.mesh))) |> display
mesh = body.mesh;
# check that the mesh is loaded correctly
un = vcat(body.srfID...)
@assert length(body.mesh) == maximum(un)

T = Float32
vertices = Array{T,2}(undef, length(body.mesh.position), 3)
for i in 1:length(mes.positio.mesh)
    vertices[i,:] .= mesh.position[i]
end
# hcat(mesh.position...)

# # tri select fewer surface
# let
#     idx = body.srfID[1]
#     mesh(body.mesh.position, 
#          GLTriangleFace.(GeometryBasics.faces(body.mesh)[idx])) |> display
#     for (i,idx) in enumerate(body.srfID[2:end])
#         mesh!(body.mesh.position, 
#               GLTriangleFace.(GeometryBasics.faces(body.mesh)[idx])) |> display
#     end
# end

# forces = zeros(Float64, (length(body.mesh),3))


# let # setting local scope for dt outside of the while loop

# coupling interface
# interface = initialize!(1,1;interface=:LumpedInterface)

# while PreCICE.isCouplingOngoing()

#     # read the data from the other participant
#     readData!(interface)

#     # measure the participant
#     WaterLilyPreCICE.update!(interface)

#     # write data to the other participant
#     writeData!(interface)
    
# end
# end
# PreCICE.finalize()