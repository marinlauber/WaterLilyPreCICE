using WaterLilyPreCICE,StaticArrays,WriteVTK

let # setting local scope for dt outside of the while loop
    
    # Simulation parameters
    L,Re,U = 2^5,100,1

    map(x,t) = x .- SA[L,L,L]

    # coupling interface
    # interface, body = initialize!(Uref,L;interface=:CalculiXInterface,map)

    # # slow ramp up of the velocity
    Ut(i,t::T) where T = i==1 ? convert(T,min(0.5*t/L,1)) : zero(T) # velocity BC

    # body = MeshBody("../CalculiX/sphere.inp";map,)

    # construct the simulation
    sim = Simulation((4L,2L,2L),Ut,L;U,ν=U*L/Re,map)
    # store = Store(sim) # allows checkpointing

    # make a writer with some attributes
    velocity(a::Simulation) = a.flow.u |> Array;
    pressure(a::Simulation) = a.flow.p |> Array;
    _body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body); 
                            a.flow.σ |> Array;)
    vorticity(a::Simulation) = (@inside a.flow.σ[I] = 
                                WaterLily.curl(3,I,a.flow.u)*a.L/a.U;
                                a.flow.σ |> Array;)
    _vbody(a::Simulation) = a.flow.V |> Array;
    mu0(a::Simulation) = a.flow.μ₀ |> Array;

    custom_attrib = Dict(
        "u" => velocity,
        "p" => pressure,
        "d" => _body,
        "ω" => vorticity,
        "v" => _vbody,
        "μ₀" => mu0
    )# this maps what to write to the name in the file

    # writer for the sim
    wr = vtkWriter("Perp-Flap"; attrib=custom_attrib)
    iter,every = 0,5

    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(interface, sim, store)

        # measure the participant
        WaterLilyPreCICE.update!(interface, sim)

        # update the this participant
        step!(sim.flow, sim.pois, sim.body, interface)

        # write data to the other participant
        writeData!(interface, sim, store)
        
        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            # save the data
            mod(iter,every)==0 && write!(wr, sim)
            iter += 1
            # ...
        end
    end
end
close(wr)
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")