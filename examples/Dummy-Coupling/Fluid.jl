using WaterLilyPreCICE,StaticArrays,WriteVTK

# make a writer with some attributes
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body); a.flow.σ |> Array;)
vtk_vorticity(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u)*a.L/a.U; a.flow.σ |> Array;)
vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_mu0(a::AbstractSimulation) = a.flow.μ₀ |> Array;
mesh_velocity(a::MeshBody) = [WaterLilyPreCICE.center(tri) for tri in a.velocity]

custom_attrib = Dict(
    "u" => vtk_velocity, "p" => vtk_pressure,
    "d" => vtk_body,  "ω" => vtk_vorticity,
    "v" => vtk_vbody, "μ₀" => vtk_mu0
)# this maps what to write to the name in the file
mesh_attrib = Dict("velocity"=>mesh_velocity)

# Simulation parameters and mapping
L,Re,U = 2^5,100,1

# construct the simulation
sim = CoupledSimulation((2L,L,L), (U,0,0), L; U, ν=U*L/Re,
                        surface_mesh=joinpath(@__DIR__,"../../meshes/sphere.inp"),
                        scale=L/2, center=SA[L÷2,L÷2,L÷2], T=Float32)
# writer for the sim
wr = vtkWriter("Dummy-Coupling"; attrib=custom_attrib)
mesh_wr = vtkWriter("Dummy-Mesh", attrib=mesh_attrib)
save!(wr, sim); save!(mesh_wr, sim.body, sim_time(sim))

let # setting local scope for dt outside of the while loop
    iter,every = 0,10

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
            mod(iter,every)==0 && (save!(wr, sim); save!(mesh_wr, sim.body, sim_time(sim)))
            iter += 1
        end
    end
    close(wr)
    close(mesh_wr)
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")