using WaterLilyPreCICE,StaticArrays,WriteVTK

# make a writer with some attributes to output to the file
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body)

# make the sim
L,Re,U = 128,1000,1
R = 0.1L/0.41 # radius
# velocity profile of Turek Hron
function uBC(i,x::SVector{N,T},t) where {N,T}
    i ≠ 1 && return convert(T, 0.0)
    return convert(T, 1.5*U*(x[2]/L)*(1.0-x[2]/L)/0.5^2)
end
# make a sim
sim = CoupledSimulation((6L,L), uBC, R; U, ν=U*R/Re, exitBC=true,
                         surface_mesh=joinpath(@__DIR__,"../CalculiX/surface.inp"),
                         passive_bodies=[AutoBody((x,t)->L/2-abs(x[2]-L/2)-1.5)], # wall at ±L/2
                         scale=R,center=SA[0.2L/0.41,0.2L/0.41,0])

# make the paraview writer
wr = vtkWriter("Turek-Hron";attrib=custom_attrib)
write!(wr, sim)
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
