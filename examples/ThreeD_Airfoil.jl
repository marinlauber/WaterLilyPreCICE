using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_sphere(;L=32,Re=1000,St=0.3,U=1,mem=Array)
    # make the body from the stl mesh
    body = MeshBody(joinpath(@__DIR__,"../meshes/naca.inp");scale=L,
                    map=(x,t)->x.+SA[L/2,L*(1+sin(2π*t*St/L)/5),0])
    # generate sim
    Simulation((3L,2L,L÷2), (U,0,0), L; ν=U*L/Re, body, mem)
end

# make a writer with some attributes to output to the file
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;) 
vtk_mu0(a::Simulation) = a.flow.μ₀ |> Array;
custom_attrib = Dict("u" => vtk_velocity, "p" => vtk_pressure, "d" =>vtk_body, "μ₀" => vtk_mu0,)

# make the sim
sim = make_sphere(L=64)

# make the paraview writer
wr = vtkWriter("Airfoil";attrib=custom_attrib)
# duration and write steps
t₀,duration,step = 0.,10.0,0.2
forces = []
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
    push!(forces,[tᵢ, fm, pf])
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)

# using Plots
# let
#     p1=plot(getindex.(forces,1),getindex.(getindex.(forces,2),1),label="Surface integral",
#             lw=2,xlabel="tU/L",ylabel="Drag",xlims=(0,100),ylims=(-5,0))
#     plot!(p1,getindex.(forces,1),getindex.(getindex.(forces,3),1),label="Volume integral",lw=2)
#     p2=plot(getindex.(forces,1),getindex.(getindex.(forces,2),2),label=:none,
#             lw=2,xlabel="tU/L",ylabel="Lift",xlims=(0,100),ylims=(-5,0))
#     plot!(p2,getindex.(forces,1),getindex.(getindex.(forces,3),2),label=:none,lw=2)
#     plot(p1, p2, layout=(1,2), dpi=300)
#     savefig("assets/Airfoil_forces.png")
# end