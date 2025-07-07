using WaterLilyPreCICE,StaticArrays,WriteVTK

# make a writer with some attributes to output to the file
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array;)
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "ω₃"=>vtk_ω)

# make the sim, Re=200
L,Re,U = 128,200,1
R = 0.1L/0.41 # radius

# velocity profile of Turek-Hron
function uBC(i,x::SVector{N,T},t) where {N,T}
    i ≠ 1 && return convert(T, 0.0)
    # make sure we have no velocity outside the channel
    return max(convert(T, 1.5*U*((x[2]-1.5f0)/(L-3))*(1.0-(x[2]-1.5f0)/(L-3))/0.5^2), 0.f0)
end

# make a sim
sim = CoupledSimulation((6L,L), uBC, R; U, ν=U*R/Re, exitBC=true,
                         surface_mesh=joinpath(@__DIR__,"../CalculiX/surface.inp"),
                         passive_bodies=[AutoBody((x,t)->L/2-abs(x[2]-L/2.f0)-1.5f0)], # wall at ±L/2
                         scale=R,center=SA[0.2L/0.41,0.2L/0.41,0])

# make the paraview writer
wr = vtkWriter("Turek-Hron";attrib=custom_attrib)
save!(wr, sim)

# integrate in time
while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(sim)

    # update the this participant and scale forces
    sim_step!(sim); sim.int.forces .*= sim.U^2/2sim.L
    sim.int.forces[:,3] .= 0.0 # zero-spanwise forces
    # println(sum(sim.int.forces,dims=1)./size(sim.int.forces,1)) # print the average force
    sim.int.forces[:,2] .= min((1-cos(2π*sim_time(sim)/2))/2,1)*sim.U^2/2sim.L*sin(2π*sim_time(sim)) # add a small perturbation in the y-direction

    # write data to the other participant
    writeData!(sim)

    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        # save the data 4 times per convective time
        # length(sim.flow.Δt)%1==0 && 
        save!(wr, sim)
    end
end
close(wr)
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")
