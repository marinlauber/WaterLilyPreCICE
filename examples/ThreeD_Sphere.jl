using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_sphere(;L=32,Re=250,U=1,mem=Array)
    # make the body from the stl mesh
    meshes = ["sphere.stl","sphere.inp","sphere_srf.inp"]
    body = MeshBody(joinpath(@__DIR__,"../meshes/"*meshes[rand(1:3)]);scale=L/2,
                    map=(x,t)->x.+L/2)
    # generate sim
    Simulation((2L,L,L), (U,0,0), L; ν=U*L/Re, body, mem)
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
wr = vtkWriter("Sphere";attrib=custom_attrib)
wr_mesh = vtkWriter("Sphere_mesh";attrib=custom)

# duration and write steps
t₀,duration,step = 0.,10.0,0.1
# run the sim
@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=false)
    write!(wr, sim); write!(wr_mesh, sim.body)
    fm = 2sum(WaterLilyPreCICE.forces(sim.body, sim.flow))/(sim.L/2)^2
    println("Surface pressure force: ", round.(fm,digits=4))
    pf = 2WaterLily.pressure_force(sim)/(sim.L/2)^2
    println("Volume integrale force: ", round.(pf,digits=4))
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr); close(wr_mesh)
