using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_sim(;L=128,Re=1e3,U=1)
    # move the geometry to the center of the domain
    map(x,t) = x .+ SA[0.2L/0.41,0.2L/0.41,0]

    # channel sdf
    function sdf(x,t)
        r = abs(x[2]-L/2.f0) # move to center of pipe
        return L/2 - r - 1.5f0 # remove radius and  the ghost
    end
    
    # make the body from the stl mesh
    flap = MeshBody(joinpath(@__DIR__,"../CalculiX/surface.inp");map,scale=0.1L/0.41)
    wall = AutoBody(sdf)
    body = CombinedBodies([flap,wall],[+])

    # velocity profile of Turek Hron
    function uBC(i,x::SVector{N,T},t) where {N,T}
        i ≠ 1 && return convert(T, 0.0)
        # make sure we have no velocity outside the channel
        return max(convert(T, 1.5*U*((x[2]-1.5f0)/(L-3))*(1.0-(x[2]-1.5f0)/(L-3))/0.5^2), 0.f0)
    end
    
    # generate sim
    Simulation((6L,L), uBC, L; U=U, ν=U*L/Re, body, exitBC=true)
end

# make a writer with some attributes to output to the file
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
vtk_vbody(a::Simulation) = a.flow.V |> Array;
vtk_mu0(a::Simulation) = a.flow.μ₀ |> Array;
vtk_ω(a::Simulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u);
                       a.flow.σ |> Array;)
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "v"=>vtk_vbody, "μ₀"=>vtk_mu0, "ω"=>vtk_ω)

# make the sim
sim = make_sim(L=64)
# make the paraview writer
wr = vtkWriter("test";attrib=custom_attrib)
# duration and write steps
duration,step = 10,0.1
# run the sim
@time for tᵢ in range(0,duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=true)
    save!(wr,sim)
    fm = sum(WaterLilyPreCICE.forces(sim.body[1],sim.flow))
    @show fm
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)