using WaterLilyPreCICE,StaticArrays,WriteVTK
using Plots

# cell centered velocity magnitude
mag(I,u) = √sum(ntuple(i->0.25*(u[I,i]+u[I+δ(i,I),i])^2,length(I)))

function make_sim(;L=64,Re=100,U=1,T=Float32)
    # move the geometry to the center of the domain
    map(x,t) = x .+ SA{T}[3L,0.,0.]

    # make the body from the stl mesh and the sdf wall
    flap = MeshBody(joinpath(@__DIR__,"../CalculiX/surface.inp");map,scale=1.f0)
    wall = AutoBody((x,t)->2L - abs(x[2]-2L) - 1.5f0)
    flap = flap + wall

    # slow impulse
    @inline C(t) = (tanh(4t-4.f0)+tanh(4.f0))/(1.f0+tanh(4.f0))

    # velocity profile of Turek Hron
    function uBC(i,x::SVector{N,T},t) where {N,T}
        i ≠ 1 && return convert(T, 0.0)
        return 1.5f0≤x[2]≤4L-1.5f0 ? C(t) : zero(T)
    end

    # generate sim
    Simulation((6L,4L), uBC, L; U=U, ν=U*L/Re, body=flap, exitBC=true)
end

# make the sim
sim = make_sim(L=64)

# duration and write steps
duration,step = 10.0,0.1
R = inside(sim.flow.p)

# run the sim
@time @gif for tᵢ in range(0,duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=false)
    @inside sim.flow.σ[I] = mag(I,sim.flow.u)
    @inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0.)<0,NaN,sim.flow.σ[I])
    flood(sim.flow.σ[R],clims=(0,2), axis=([], false),
          cfill=:jet,legend=false,border=:none,size=size(sim.flow.p))
    fm = 2sum(WaterLilyPreCICE.forces(sim.body, sim.flow))/(sim.L/2)^2
    println("Surface pressure force: ", round.(fm,digits=4))
    pf = 2WaterLily.pressure_force(sim.flow, sim.body.a)/(sim.L/2)^2
    println("Volume integrale force: ", round.(pf,digits=4))
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
