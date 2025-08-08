using WaterLilyPreCICE,StaticArrays,WriteVTK
using Plots

# cell centered velocity magnitude
mag(I,u) = √sum(ntuple(i->0.25*(u[I,i]+u[I+δ(i,I),i])^2,length(I)))

function make_sim(;L=128,Re=1e3,U=1,T=Float32)
    # move the geometry to the center of the domain
    # map(x,t) = x .+ SA{T}[0.25L/0.41,0.2L/0.41,0]
    map(x,t) = x .+ SA{T}[4L,4L,0]

    # make the body from the stl mesh and the sdf wall
    flap = MeshBody(joinpath(@__DIR__,"../CalculiX/surface.inp");map,scale=1.f0)
    wall = AutoBody((x,t)->4L - abs(x[2]-4L) - 1.5f0)
    # wall = AutoBody((x,t)->L÷2 - abs(x[2]-L÷2) - 1.5f0)
    # cylinder = AutoBody((x,t)->√sum(abs2,x.-0.2f0L/0.41f0)-0.05f0L/0.41f0)
    cylinder = AutoBody((x,t)->√sum(abs2,x.-4L)-0.05f0L/0.41f0)
    TurekHron = cylinder + flap + wall

    # impulse
    @inline C(t) = (tanh(4t-4.f0)+tanh(4.f0))/(1.f0+tanh(4.f0))

    # velocity profile of Turek Hron
    function uBC(i,x::SVector{N,T},t) where {N,T}
        i ≠ 1 && return convert(T, 0.0)
        # make sure we have no velocity outside the channel
        return C(t)*max(convert(T, 1.5*U*((x[2]-1.5f0)/(L-3))*(1.0-(x[2]-1.5f0)/(L-3))/0.5^2), 0.f0)
    end

    # generate sim (6L,L)
    Simulation((16L,8L), uBC, L; U=U, ν=U*L/Re, body=TurekHron, exitBC=false)
end

# make the sim
sim = make_sim(L=128)

# duration and write steps
duration,step = 100,0.1
R = inside(sim.flow.p)

# run the sim
@time @gif for tᵢ in range(0,duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=false)
    @inside sim.flow.σ[I] = sim.flow.p[I] #mag(I,sim.flow.u)
    @inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0.)<0,NaN,sim.flow.σ[I])
    flood(sim.flow.σ[R],clims=(-2,2), axis=([], false),
          cfill=:jet,legend=false,border=:none,size=(6*sim.L,sim.L))
    fm = 2sum(WaterLilyPreCICE.forces(sim.body, sim.flow))/(sim.L/2)^2
    println("Surface pressure force: ", round.(fm,digits=4))
    pf = 2WaterLily.pressure_force(sim.flow, sim.body.a.b)/(sim.L/2)^2
    println("Volume integrale force: ", round.(pf,digits=4))
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
