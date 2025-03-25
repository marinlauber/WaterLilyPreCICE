using WaterLilyPreCICE,StaticArrays,WriteVTK

# make a writer with some attributes to output to the file
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body)

# make the sim
L,Re,U = 64,1000,1
sim = CoupledSimulation((3L,2L,L÷2), (U,0,0), L; ν=U*L/Re, perdir=(3,),
                         surface_mesh=joinpath(@__DIR__,"../CalculiX/surface.inp"),
                         scale=L,center=SA[L/2,L,0])

# make the paraview writer
wr = vtkWriter("Airfoil";attrib=custom_attrib)

let
    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(sim)

        # update the this participant and scale forces
        sim_step!(sim); sim.int.forces .*= 2sim.U^2/sim.L
        sim.int.forces[:,3] .= 0.0 # zero-spanwise forces 
       
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
