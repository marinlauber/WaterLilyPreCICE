using WaterLilyPreCICE,StaticArrays,WriteVTK

# make a writer with some attributes
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body); a.flow.σ |> Array;)
vtk_vorticity(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u)*a.L/a.U; a.flow.σ |> Array;)
vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_mu0(a::AbstractSimulation) = a.flow.μ₀ |> Array;

custom_attrib = Dict(
    "u" => vtk_velocity, "p" => vtk_pressure,
    "d" => vtk_body,  "ω" => vtk_vorticity,
    "v" => vtk_vbody, "μ₀" => vtk_mu0
)# this maps what to write to the name in the file

# Simulation parameters and mapping
L,Re,U = 2^5,100,1
map(x,t) = x.-SA[L,L÷2,L÷2]

# construct the simulation
sim = CoupledSimulation((2L,L,L),(U,0,0),L;U,ν=U*L/Re,map,
                        fname="../CalculiX/sphere.inp",scale=L)

# writer for the sim
wr = vtkWriter("WaterLily-CalculiX"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop
    iter,every = 0,5

    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(sim)

        # update the this participant
        sim_step!(sim)

        # write data to the other participant
        writeData!(sim)

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            # save the data
            mod(iter,every)==0 && save!(wr, sim)
            iter += 1
        end
    end
    close(wr)
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")