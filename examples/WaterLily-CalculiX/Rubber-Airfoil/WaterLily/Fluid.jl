using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_sphere(;L=32,Re=250,U=1,mem=Array)
    # make the body from the stl mesh
    # body = MeshBody(joinpath(@__DIR__,"../CalculiX/surface.inp");scale=L)
    # generate sim
    # Simulation((2L,L,L÷2), (U,0,0), L; ν=U*L/Re, body, mem, perdir=(3,))
    CoupledSimulation((2L,L,L÷2), (U,0,0), L; ν=U*L/Re, mem, perdir=(3,),
                      interface=:Interface,surface_mesh=joinpath(@__DIR__,"../CalculiX/surface.inp"),
                      scale=L,center=SA[L/2,L/2,0])
end

# make a writer with some attributes to output to the file
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
custom_attrib = Dict("u" => vtk_velocity,
                     "p" => vtk_pressure, 
                     "d" => vtk_body
)

# make the sim
sim = make_sphere(L=64)

# move to the center of the domain
# WaterLily.update!(sim.body,update=(x)->x.+SA[sim.L/2,sim.L/2,0])
# measure!(sim)

# make the paraview writer
wr = vtkWriter("Airfoil";attrib=custom_attrib)

let
    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(sim)

        # update the this participant
        sim_step!(sim)

        # write data to the other participant
        writeData!(sim)

        println("WaterLily: Time=",round(WaterLily.time(sim.flow),digits=4),
                           ", Δt=",round(sim.flow.Δt[end],digits=3))
        
        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            # save the data
            length(sim.flow.Δt)%5==0 && write!(wr, sim)
        end
    end
    close(wr)
end
PreCICE.finalize()  
println("WaterLily: Closing Julia solver...")

# # duration and write steps
# t₀,duration,step = 0.,10.0,0.1
# # run the sim
# @time for tᵢ in range(t₀,t₀+duration;step)
#     # update until time tᵢ in the background
#     sim_step!(sim,tᵢ;remeasure=false)
#     write!(wr, sim)
#     fm = 2sum(WaterLilyPreCICE.forces(sim.body, sim.flow))/(sim.L/2)^2
#     println("Surface pressure force: ", round.(fm,digits=4))
#     pf = 2WaterLily.pressure_force(sim)/(sim.L/2)^2
#     println("Volume integrale force: ", round.(pf,digits=4))
#     # print time step
#     println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# end
# close(wr)