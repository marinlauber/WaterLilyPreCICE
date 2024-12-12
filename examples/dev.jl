using WaterLily
using WaterLilyPreCICE
using GeometryBasics
using GLMakie

GLMakie.activate!(inline=false)
# create the mesh body
MeshBody = WaterLilyPreCICE.MeshBody
body = MeshBody("/home/marin/Workspace/HHH/code/LIMO_heart/LIMO_4_contact/geom_restart.inp";boundary=false,thk=2,scale=1.f00)
mesh(body.mesh.position, GLTriangleFace.(GeometryBasics.faces(body.mesh))) |> display

# check that the mesh is loaded correctly
un = unique(vcat(body.srfID...))
@assert length(body.mesh) == maximum(un)

# tri select fewer surface
let
    idx = body.srfID[1]
    mesh(body.mesh.position, 
         GLTriangleFace.(GeometryBasics.faces(body.mesh)[idx])) |> display
    for (i,idx) in enumerate(body.srfID[2:end])
        mesh!(body.mesh.position, 
              GLTriangleFace.(GeometryBasics.faces(body.mesh)[idx])) |> display
    end
end

forces = zeros(Float64, (3,length(body.mesh)))
A1(t) =ifelse(t<10,t/10,ifelse(t<12,1,(1+(t-12))))
A2(t) = ifelse(t<10,t/10,ifelse(t<12,1,(1+(t-12))))
A3(t) = ifelse(t<10,t/10,ifelse(t<12,1,(1+(t-12))))
function static_inflation(i,t)
    i==1 && return  0.38*A1(t)
    i==6 && return -0.38*A1(t)
    i in [2,3] && return -0.1*A3(t)
    i in [4,5] && return 0.48*A2(t)
    i in [7,8] && return 0.1*A3(t)
    i in [9,10] && return -0.48*A2(t)
end
function static_pressure!(forces,body,func)
    for (i,face) in enumerate(body.srfID)
        forces[:,face] .= func(i,t)
    end
end