
using WaterLily
# using WaterLilyPreCICE
using StaticArrays
using FileIO,MeshIO
using WriteVTK
include("src/mesh_body.jl")
file="cone_10"
#To compare speedup and want to use layout
function measure_test(;L=32,Re=5e5,U=1,mem=Array)
    
    body=MeshBody("obj/"*file*".obj";scale=1,boundary=false,thk=(2+3^0.5))

    Simulation((L,L,L), (U,0,0), L; ν=U*L/Re, body=body, mem)
end

_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
custom_attrib = Dict("d" => _body, )

sim = measure_test(L=128)

measure!(sim)

wr=vtkWriter(file*"meshbody";attrib=custom_attrib)

write!(wr, sim)

close(wr)
# println("dist with waterlilly calculated")


