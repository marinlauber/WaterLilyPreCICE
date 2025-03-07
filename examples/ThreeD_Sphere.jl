using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_sphere(;L=32,Re=250,U=1,mem=Array)
    # make the body from the stl mesh
    body = MeshBody(joinpath(@__DIR__,"../meshes/sphere.stl");scale=L/2,
                    map=(x,t)->x.+L/2)
    # generate sim
    Simulation((2L,L,L), (U,0,0), L; ν=U*L/Re, body, mem)
end

# make a writer with some attributes to output to the file
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;
mu0(a::Simulation) = a.flow.μ₀ |> Array;
custom_attrib = Dict("u" => velocity, "p" => pressure, "d" => _body, "v" => _vbody, "μ₀" => mu0,)

# make the sim
sim = make_sphere(L=64)

# make the paraview writer
wr = vtkWriter("STL_mesh_test";attrib=custom_attrib)
# duration and write steps
t₀,duration,step = 0.,10.0,0.1
# run the sim
@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    # sim_step!(sim,tᵢ;remeasure=false)
    while sim_time(sim) < tᵢ
        WaterLilyPreCICE.update!(sim.body,sum(sim.flow.Δt),sim.flow.Δt[end])
        sim_step!(sim; remeasure=true)
    end
    write!(wr, sim)
    fm = 2sum(WaterLilyPreCICE.forces(sim.body, sim.flow))/(sim.L/2)^2
    println("Surface pressure force: ", round.(fm,digits=4))
    pf = 2WaterLily.pressure_force(sim)/(sim.L/2)^2
    println("Volume integrale force: ", round.(pf,digits=4))
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)