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
            mod(iter,every)==0 && write!(wr, sim)
            iter += 1
        end
    end
    close(wr)
end
PreCICE.finalize()  
println("WaterLily: Closing Julia solver...")