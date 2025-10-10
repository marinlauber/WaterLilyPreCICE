using WaterLilyPreCICE,StaticArrays,WriteVTK

# pins the outlet pressure to 0
pin!(a::AbstractArray) = (a .-= a[end-1,size(a,2)÷2])

# make a writer with some attributes to output to the file
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_V(a::AbstractSimulation) = a.flow.V |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array;)
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "ω₃"=>vtk_ω,"V"=>vtk_V)

# make the sim, Re=200
L,Re,U = 128,200,1
D = 0.1L/0.41 # radius

# slow initial ramp
@inline C(t) = (tanh(4t-4.f0)+tanh(4.f0))/(1.f0+tanh(4.f0))

# velocity profile of Turek-Hron
function uBC(i,x::SVector{N,T},t) where {N,T}
    i ≠ 1 && return convert(T, 0.0)
    # make sure we have no velocity outside the channel
    return C(t)*max(convert(T, 1.5*U*((x[2]-1.5f0)/(L-3))*(1.0-(x[2]-1.5f0)/(L-3))/0.5^2), 0.f0)
end

# make a sim
sim = CoupledSimulation((6L,L), uBC, D; U, ν=U*D/Re, exitBC=true,
                         surface_mesh=joinpath(@__DIR__,"../CalculiX/surface.inp"),
                         passive_bodies=[AutoBody((x,t)->L÷2 - abs(x[2]-L÷2) - 1.5f0),
                                         AutoBody((x,t)->√sum(abs2,x.-0.2f0L/0.41f0)-0.05f0L/0.41f0)], # wall at ±L/2 and cylinder
                         scale=1.0,center=SA[0.25L/0.41,0.2L/0.41,0],T=Float32)

# make the paraview writer
wr = vtkWriter("Turek-Hron";attrib=custom_attrib)
save!(wr, sim)

# integrate in time
while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(sim)

    # update the participant, and pin the pressure
    sim_step!(sim); pin!(sim.flow.p)

    # force computation
    sim.int.forces .=  0
    sim.int.forces[:,2] .= 0.2*sim.U^2/sim.L*sin(2π*0.19*sim_time(sim))

    # write data to the other participant
    writeData!(sim)

    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        # save the data 4 times per convective time
        length(sim.flow.Δt)%10==0 && save!(wr, sim)
    end
end
close(wr)
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")
