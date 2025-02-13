using WaterLilyPreCICE,StaticArrays,WriteVTK

# make a writer with some attributes
velocity(a::AbstractSimulation) = a.flow.u |> Array;
pressure(a::AbstractSimulation) = a.flow.p |> Array;
_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body); a.flow.σ |> Array;)
vorticity(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u)*a.L/a.U; a.flow.σ |> Array;)
_vbody(a::AbstractSimulation) = a.flow.V |> Array;
mu0(a::AbstractSimulation) = a.flow.μ₀ |> Array;

custom_attrib = Dict(
    "u" => velocity, "p" => pressure,
    "d" => _body,  "ω" => vorticity,
    "v" => _vbody, "μ₀" => mu0
)# this maps what to write to the name in the file

# Simulation parameters and mapping
L,Re,U = 2^6,250,1

# construct the simulation
sim = CoupledSimulation((8L,4L),(U,0),L;U,ν=U*L/Re,
                        surface_mesh="../CalculiX/geom.inp",scale=1.f0,center=SA[2L,2L,0],
                        boundary=false,thk=4,ϵ=1) # this one is the centerline, not the boundary

# writer for the sim
wr = vtkWriter("Inverted-flag"; attrib=custom_attrib)

# run
while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(sim)

    # update the this participant
    sim_step!(sim)
    @show size(sim.int.forces),sim_time(sim)
    sim.int.forces .= 0
    sim.int.forces[:,2] .= -5.0

    # write data to the other participant
    writeData!(sim)
    
    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        # save the data
        mod(length(sim.flow.Δt),10)==0 && write!(wr, sim)
    end
end
close(wr)
PreCICE.finalize()  
println("WaterLily: Closing Julia solver...")