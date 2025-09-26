using WaterLilyPreCICE,ParametricBodies,StaticArrays,WriteVTK

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

# needed to enable set operations on boundary curves
# ParametricBodies.notC¹(l::NurbsLocator{C},uv) where C<:NurbsCurve{n,d} where {n,d} = false

# Simulation parameters
L,Re,U,Lref,Uref = 2^5,100,1,1.0,10.0
center = SA[2.95L,0]

# slow ramp up of the velocity
Ut(i,t::T) where T = i==1 ? convert(T,min(0.5*t/L,1)) : zero(T) # velocity BC

# construct the simulation
sim = CoupledSimulation((8L,4L), Ut, L; U, ν=U*L/Re,
                        Δt=0.4, interface=:GismoInterface,
                        center, dir=[1,-1,-1])

# writer for the sim
wr = vtkWriter("WaterLily-Gismo"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop

    iter,every = 0,5

    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(sim)

        # update the this participant
        sim_step!(sim)

        # access the data and change it during the run
        length(sim.flow.Δt)==2 && (sim.int.forces .= 0.0) # first time step
        sim.int.forces .*= Uref^2 # scale

        # write data to the other participant
        writeData!(sim)

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            # save the data
            mod(iter,every)==0 && save!(wr, sim)
            iter += 1
            # ...
        end
    end
end
close(wr)
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")