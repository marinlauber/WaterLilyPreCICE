"""
The simulated flow domain is 6 units long (x) and 4 units tall (y).
The flap is located at the center of the bottom (x=0) and is 1 unit long (y) and 0.1 units thick (x).
We set the fluid density ρf = 1.0 kg/m³, the kinematic viscosity νf = 1.0 m²/s,
the solid density 3000 kg/m³, the Young’s modulus to E = 4.0×10⁶ kg/ms² and the Poisson ratio νs = 0.3.
On the left boundary a constant inflow profile in x-direction of 10 m/s is prescribed.
The right boundary is an outflow and the top and bottom of the channel as well as the surface of the flap are no-slip walls.
"""

using WaterLilyPreCICE,StaticArrays,WriteVTK

# make a writer with some attributes to output to the file
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_velbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array;)
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "ω₃"=>vtk_ω, "vᵦ"=>vtk_velbody)

# make the sim, Re=10
L,Re,U = 64,100,1.f0

# slow impulse
@inline C(t) = (tanh(4t-4.f0)+tanh(4.f0))/(1.f0+tanh(4.f0))

# velocity profile of Turek Hron
function uBC(i,x::SVector{N,T},t) where {N,T}
    i ≠ 1 && return convert(T, 0.0)
    # make sure we have no velocity outside the channel
    # return C(t)*max(convert(T, 1.5*U*((x[2]-1.5f0)/(4L-3))*(1.0-(x[2]-1.5f0)/(4L-3))/0.5^2), 0.f0)
    return 1.5f0≤x[2]≤4L-1.5f0 ? C(t) : zero(T)
end


# make a sim
sim = CoupledSimulation((6L,4L), uBC, L; U, ν=U*L/Re, exitBC=true,
                         surface_mesh=joinpath(@__DIR__,"../CalculiX/surface.inp"),
                         passive_bodies=[AutoBody((x,t)->2L - abs(x[2]-2L) - 1.5f0)], # wall at ±L/2 and cylinder
                         scale=1.0, center=SA[3L,0,0], T=Float32)

# make the paraview writer
wr = vtkWriter("Flap";attrib=custom_attrib)
save!(wr, sim)

# integrate in time
while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(sim)

    # update the this participant and scale forces
    sim_step!(sim); #sim.int.forces .*= sim.U^2/2sim.L
    sim.int.forces[:,3] .= 0.0 # zero-spanwise forces

    sim.int.forces .= 0
    sim.int.forces[:,1] .= min(sim_time(sim)/2,1.0)
    println(sum(sim.int.forces,dims=1)./size(sim.int.forces,1)) # print the average force

    # write data to the other participant
    writeData!(sim)

    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        # save the data 4 times per convective time
        length(sim.flow.Δt)%25==0 && save!(wr, sim)
    end
end
close(wr)
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")
