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
map(x,t) = x.-SA[2L,2L]

# construct the simulation
sim = CoupledSimulation((8L,4L),(U,0),L;U,ν=U*L/Re,map,
                        fname="../CalculiX/geom.inp",scale=1.f0,
                        boundary=false,thk=4,ϵ=1) # this one is the centerline, not the boundary

# writer for the sim
wr = vtkWriter("Inverted-flag"; attrib=custom_attrib)

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
        mod(length(sim.flow.Δt),5)==0 && write!(wr, sim)
    end
end
close(wr)
PreCICE.finalize()  
println("WaterLily: Closing Julia solver...")