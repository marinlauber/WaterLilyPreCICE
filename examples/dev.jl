using WaterLily,WriteVTK
using WaterLilyPreCICE
using GeometryBasics
using Interpolations

# # create the mesh body
body = MeshBody("/home/marin/Workspace/WaterLilyPreCICE/examples/LIMO_4/Solid/geom.inp";boundary=false,thk=2,scale=1.f00)
# check that the mesh is loaded correctly
un = vcat(body.srfID...)
@assert length(body.mesh) == maximum(un)

function static_inflation(i)
    i==1 && return  0.38
    i==6 && return -0.38
    i in [2,3] && return -0.1
    i in [4,5] && return  0.48
    i in [7,8] && return  0.1
    i in [9,10] && return -0.48
end

_srf(a::MeshBody) = getindex.(mapreduce(((i,ids),)->map(T->(i,T),ids),vcat,enumerate(a.srfID)),1)
_center(a::MeshBody) = WaterLilyPreCICE.center.(a.mesh)
_normal(a::MeshBody) = WaterLilyPreCICE.normal.(a.mesh)
_dS(a::MeshBody) = _p(a).*WaterLilyPreCICE.normal.(a.mesh)
function _p(a::MeshBody) 
    pressure = zeros(size(a.mesh))
    for (i,id) in mapreduce(((i,ids),)->map(T->(i,T),ids),vcat,enumerate(a.srfID))
        pressure[id] = static_inflation(i)
    end
    return pressure
end

custom = Dict(
    "SRF" =>_srf, "center"=>_center, "normal"=>_normal, "dS" => _dS, "p"=> _p
)
wr = vtkWriter("MeshBody"; attrib=custom)
write!(wr,body)
close(wr)

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