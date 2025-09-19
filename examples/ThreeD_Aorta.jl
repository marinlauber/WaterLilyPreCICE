using WaterLilyPreCICE,StaticArrays,WriteVTK

import WaterLily: @loop
struct ActuatorDisk{T,F<:Function}
    R::T
    map::F
    function ActuatorDisk(R=8; map=(x,t)->x)
        new{typeof(R),typeof(map)}(R, map)
    end
end
@fastmath @inline r(x) = √(x[2]^2+x[3]^2)
@fastmath @inline inside(x::SVector, disk::ActuatorDisk) = √(x[1]^2+(r(x)-min(r(x),))^2) ≤ 2.f0
function flux!(Ii::CartesianIndex,disk::ActuatorDisk{T},t) where T
    # map to the location of the disk
    x = disk.map(loc(Ii,T),t); i = Base.last(Ii)
    # if we are not inside the disk, no forcing
    !inside(x,disk) && return zero(T)
    # what is the flux there
    return one(T)
end
function force_disk!(flow, t; disk)
    @loop flow.f[Ii] += flux(Ii,disk,t) over Ii in CartesianIndices(flow.u)
end

function make_sphere(;L=32,Re=250,U=1,mem=Array)
    # make the body from the stl mesh
    body = MeshBody(joinpath(@__DIR__,"../meshes/aorta.inp");scale=L/2,
                    map=(x,t)->x+SA[L/2,L/2,L/4],boundary=false,thk=2.f0)
    # generate sim
    Simulation((L,L,L÷2), (0,0,0), L; ν=U*L/Re, body, mem)
end

# make a writer with some attributes to output to the file
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
vtk_mu0(a::Simulation) = a.flow.μ₀ |> Array;
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "μ₀"=>vtk_mu0)
# vtk attributes
vtk_srf(a::MeshBody) = Float32[el[1] for el in a.surf_id]
vtk_f(a::MeshBody) = -WaterLilyPreCICE.forces(sim.body, sim.flow)
custom = Dict("srf" =>vtk_srf, "f"=>vtk_f)

# make the sim
sim = make_sphere(L=64)

# make the paraview writer
wr = vtkWriter("Aorta";attrib=custom_attrib)
wr_mesh = vtkWriter("Aorta_mesh";attrib=custom)
save!(wr, sim); save!(wr_mesh, sim.body)
# # duration and write steps
# t₀,duration,step = 0.,10.0,0.1
# # run the sim
# @time for tᵢ in range(t₀,t₀+duration;step)
#     # update until time tᵢ in the background
#     sim_step!(sim,tᵢ;remeasure=false)
#     save!(wr, sim); save!(wr_mesh, sim.body)
#     fm = 2sum(WaterLilyPreCICE.forces(sim.body, sim.flow))/(sim.L/2)^2
#     println("Surface pressure force: ", round.(fm,digits=4))
#     pf = 2WaterLily.pressure_force(sim)/(sim.L/2)^2
#     println("Volume integrale force: ", round.(pf,digits=4))
#     # print time step
#     println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# end
close(wr); close(wr_mesh)
