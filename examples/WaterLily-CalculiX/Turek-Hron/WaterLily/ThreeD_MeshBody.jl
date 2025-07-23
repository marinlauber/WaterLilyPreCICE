using WaterLilyPreCICE,StaticArrays,WriteVTK
using Plots

# cell centered velocity magnitude
mag(I,u) = √sum(ntuple(i->0.25*(u[I,i]+u[I+δ(i,I),i])^2,length(I)))

function make_sim(;L=128,Re=1e3,U=1,T=Float32)
    # move the geometry to the center of the domain
    map(x,t) = x .+ SA{T}[0.25L/0.41,0.2L/0.41,0]

    # make the body from the stl mesh and the sdf wall
    flap = MeshBody(joinpath(@__DIR__,"../CalculiX/surface.inp");map,scale=1.f0)
    wall = AutoBody((x,t)->L÷2 - abs(x[2]-L÷2) - 1.5f0)
    cylinder = AutoBody((x,t)->√sum(abs2,x.-0.2f0L/0.41f0)-0.05f0L/0.41f0)
    TurekHron = cylinder + flap + wall

    # velocity profile of Turek Hron
    function uBC(i,x::SVector{N,T},t) where {N,T}
        i ≠ 1 && return convert(T, 0.0)
        # make sure we have no velocity outside the channel
        return max(convert(T, 1.5*U*((x[2]-1.5f0)/(L-3))*(1.0-(x[2]-1.5f0)/(L-3))/0.5^2), 0.f0)
    end

    # generate sim
    Simulation((6L,L), uBC, L; U=U, ν=U*L/Re, body=TurekHron, exitBC=true)
end

# make a writer with some attributes to output to the file
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_mu0(a::AbstractSimulation) = a.flow.μ₀ |> Array;
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array;)
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "v"=>vtk_vbody, "μ₀"=>vtk_mu0, "ω"=>vtk_ω)

# make the sim
sim = make_sim(L=128)
# make the paraview writer
# wr = vtkWriter("Turek-Hron";attrib=custom_attrib)
R = inside(sim.flow.p)
# duration and write steps
duration,step = 3,0.1
# run the sim
@time @gif for tᵢ in range(0,duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=true)
    # save!(wr,sim)
    @inside sim.flow.σ[I] = mag(I,sim.flow.u)
    @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,NaN,sim.flow.σ[I])
    flood(sim.flow.σ[R],clims=(0,2), axis=([], false),
          cfill=:jet,legend=false,border=:none,size=(6*sim.L,sim.L))
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
# close(wr)