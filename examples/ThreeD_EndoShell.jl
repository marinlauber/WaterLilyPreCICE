"""
TODO:
  - [x] check velocity of center of elements
  - [x] use exact triangulation
  - [x] include in WaterLily simulation
  - [x] correct scale and motion of top
  - [x] measure on the GPU?
  - [x] add cap and valves
  - [ ] include flow rates from 0D model
  - [ ] interpolate mesh in time?
"""

using WaterLily,StaticArrays,GeometryBasics
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.25)

function make_sim(;L=64,Re=5000,U=1,mem=Array,T=Float32)
    # extract the good stuff
    data = open("/home/marin/Workspace/WaterLilyPreCICE/meshes/endo/Endo_Kasra_xyz.txt", "r") do file
        mapreduce(L->parse.(T, split(L)), hcat, readlines(file)[10:end])
    end
    # remove the reference height
    zmin,zmax = extrema(data[3,:]); scale = zmax-zmin
    # maximum radius
    R_max = max(diff([extrema(data[1,:])...])[1] , diff([extrema(data[2,:])...])[1])/2.f0/scale*L
    println("extrema of x-position: ", R_max)
    # remove 3 first columns since these are the ref config, make it unit-length and scale
    motion = [data[1+i*3:3+i*3,:]./scale.*L for i in 1:size(data)[1]÷3-1]
    # the top surface moves up and down, we must fix it
    motion = [m.+SA{T}[L,L,5L/4-maximum(m[3,:])] for m in motion]

    # read it from the connecivity file
    connectivity = open("/home/marin/Workspace/WaterLilyPreCICE/meshes/endo/connectivity.txt", "r") do f
        mapreduce(L->parse.(Int, split(L)), hcat, readlines(f))
    end

    # make points list from first snapshot, we use the last config such that
    # measuring first sets the correct velocity BC
    points = [Point3f(pnt) for pnt in eachcol(motion[end])]

    # initial mesh, create triangulation
    # connectivity must be incremented since we have 1 indexing here
    faces = [TriangleFace{Int}(GLTriangleFace((reverse(tri).+1...))) for tri in eachcol(connectivity)]

    # generate a mesh from this
    mesh = GeometryBasics.Mesh(points, faces)
    shell = MeshBody(mesh;scale=1.f0,half_thk=2f0,boundary=false,mem=mem)

    # add a top surface
    R₁,R₂ = T(1.25*R_max),T(L/9)
    println("R₁=",R₁,", R₂=",R₂)
    println("R₁=",typeof(R₁),", R₂=",typeof(R₂))
    @fastmath @inline r(x) = √(x[1]^2+x[2]^2)
    # cap = AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),R₁))^2)-2.0f0,
    #                (x,t)->x.-SA{T}[L/2.f0,L/2.f0,5.f0L÷4.f0])
    # cap -= AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),R₂))^2)-4.0f0,
    #                 (x,t)->x.-SA{T}[L/2.f0,L/2.f0,5.f0L÷4.f0])

    # make the body
    body = shell #+ cap

    # make motion usable for update!
    motion_update = []; for m in motion
        n_mesh = GeometryBasics.Mesh([Point3f(pnt) for pnt in eachcol(m)], faces)
        push!(motion_update, [hcat([n_mesh[i]...]...) for i in 1:length(n_mesh)] |> mem)
    end

    # make a simulation
    Simulation((2L,2L,2L),(0,0,0),L;U,ν=U*L/Re,body,mem,T), motion_update, shell
end

# find a body of type T in a sim
# function find(a::WaterLily.SetBody, ::Type{T}) where T
#     isa(a.a, T) && return a.a
#     isa(a.b, T) && return a.b
#     isa(a.a, WaterLily.SetBody) && return find(a.a, T)
#     isa(a.b, WaterLily.SetBody) && return find(a.b, T)
#     return nothing
# end
using CUDA
# make simulation
sim,motion,mesh_body = make_sim(L=64;mem=CuArray)

# find max displacement, and max velocity
# u = 1.f0
# dx = [√maximum(maximum.(sum.(abs2,motion[end-i].-motion[end-1-i],dims=1))) for i in 0:99]
# dt = u./dx
# println("Maximum dt: ",maximum(dt),", minimum dt=",minimum(dt))
# dt = minimum(dt) # use the maximum dt
# dt = 0.25f0

# velocity of center of elements
mesh_velocity(a) = [SVector(sum(tri,dims=2)/3) for tri in Array(a.velocity)]
mesh_wr = vtkWriter("endo_mesh", attrib=Dict("velocity"=>mesh_velocity))

# flow writer
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array);
vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_lambda(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I,a.flow.u); a.flow.σ |> Array);
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, sim_time(a)); a.flow.σ |> Array);
new_attrib = Dict("u"=>vtk_velocity, "V"=>vtk_vbody, "p"=>vtk_pressure, "d"=>vtk_body, "λ₂"=>vtk_lambda, "ω₃"=>vtk_ω)
wr = vtkWriter("endo_flow"; attrib=new_attrib)

# some time benchmarking
update_time = []
measure_time = []
step_time = []
save_time = []

# dt = 0.01f0*sim.L/sim.U

# function WaterLily.sim_step!(sim::AbstractSimulation,t_end;remeasure=true,λ=quick,max_steps=typemax(Int),verbose=false,
#                              udf=nothing,strict=false,kwargs...)
#     steps₀ = length(sim.flow.Δt)
#     while sim_time(sim) < t_end && length(sim.flow.Δt) - steps₀ < max_steps
#         strict ? enforce!(sim, t_end) : nothing # enforce t+Δt == t_end
#         sim_step!(sim; remeasure, λ, udf, kwargs...) # t→t+Δt
#         verbose && sim_info(sim)
#     end
# end
# function enforce!(a, t)
#     t_left = t - sum(a.flow.Δt)*sim.U/sim.L # what's left after next step
#     t_left < 0 && (a.flow.Δt[end] += round(t_left*sim.L/sim.U, RoundUp, digits=4)) # if we overshoot, reduce the next step
# end


# the function
function get_motion!(reference, motion, t)
    ts = collect(range(0,1,length(motion)))
    k = clamp(searchsortedfirst(ts, t)-1, firstindex(ts), lastindex(ts)-1);
    y = (clamp(t, minimum(ts), maximum(ts))-ts[k])/(ts[k+1]-ts[k])
    println("t/T=",round(t,digits=3),", k=",k,", y=",round(y,digits=3))
    reference .= (1-y).*motion[k] .+ y.*motion[k+1]
    return nothing
end

# reference = similar(motion[1])
# get_motion!(reference, motion, 0) # test the interpolation
# all(reference .≈ motion[1]) # test the first position
# get_motion!(reference, motion, 0.005) # test the interpolation
# all(reference .≈ 0.5.*motion[1] .+ 0.5.*motion[2]) # test the interpolation
# get_motion!(reference, motion, 0.01+eps()) # test the interpolation
# all(reference .≈ motion[2]) # test the interpolation

# integrate in time, one convective time in one cycle
cycle = 5.f0
reference = similar(motion[1])
for (i,t) in enumerate(range(0,cycle,step=0.02))
    while sim_time(sim) < t
        # now update the mesh, this sets the velocity BC
        get_motion!(reference, motion, (sim_time(sim)/cycle)%1) # get the reference position
        tic = time()
        sim.body = update!(sim.body, reference, sim.flow.Δt[end])
        toc = time()
        push!(update_time,toc-tic)

        tic = time()
        measure!(sim) # update the body on the flow
        toc = time()
        push!(measure_time,toc-tic)

        # update flow until t
        tic = time()
        sim_step!(sim; remasure=true)
        # push!(sim.flow.Δt, sim.flow.Δt[end])
        toc = time()
        push!(step_time,toc-tic)
    end

    tic = time()
    # measure velocity and save
    save!(mesh_wr,mesh_body,sim_time(sim))
    save!(wr,sim)
    toc = time()
    push!(save_time,toc-tic)

    println("tU/L=",round(sim_time(sim),digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(mesh_wr); close(wr)
using Plots
plot([update_time, measure_time, step_time],
     title="Time per step", xlabel="Step", ylabel="Time (s)",
     label=["Update" "Measure" "Step" "Save"], yscale=:log10)