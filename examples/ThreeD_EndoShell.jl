"""
TODO:
  - [x] check velocity of center of elements
  - [x] use exact triangulation
  - [x] include in WaterLily simulation
  - [x] correct scale and motion of top
  - [ ] measure on the GPU?
  - [x] add cap and valves
  - [ ] include flow rates from 0D model
  - [ ] interpolate mesh in time?
"""

using WaterLilyPreCICE,StaticArrays
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.5)

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
    motion = [m.+SA{T}[L/2,L/2,5L/4-maximum(m[3,:])] for m in motion]

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
    mesh = GeometryBasics.Mesh(points,faces)
    shell = MeshBody(mesh;thk=2f0,boundary=false)

    # add a top surface
    R₁,R₂ = T(1.25*R_max),T(L/9)
    println("R₁=",R₁,", R₂=",R₂)
    @fastmath r(x) = √(x[1]^2+x[2]^2)
    cap = AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),R₁))^2)-2.0f0,
                   (x,t)->x.-SA{T}[L/2.f0,L/2.f0,5.f0L÷4.f0])
    cap -= AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),R₂))^2)-4.0f0,
                   (x,t)->x.-SA{T}[L/2.f0,L/2.f0,5.f0L÷4.f0])

    # make the body
    body = shell + cap

    # make a simulation
    Simulation((L,L,2L),(0,0,0),L;U,ν=U*L/Re,body,mem,T), motion
end

# find a body of type T in a sim
function find(a::WaterLily.SetBody, ::Type{T}) where T
    isa(a.a, T) && return a.a
    isa(a.b, T) && return a.b
    isa(a.a, WaterLily.SetBody) && return find(a.a, T)
    isa(a.b, WaterLily.SetBody) && return find(a.b, T)
    return nothing
end

# make simulation
sim,motion = make_sim(L=128);
mesh_body = find(sim.body, MeshBody)

# find max displacement, and max velocity
u = 0.5
dx = [√maximum(sum(abs2,motion[end-i] .- motion[end-i-1], dims=1)) for i in 0:99]
dt = dx./u
println("Maximum dt: ",maximum(dt),", minimum dt=",minimum(dt))
dt = maximum(dt) # use the maximum dt

# velocity of center of elements
mesh_velocity(a::MeshBody) = [WaterLilyPreCICE.center(tri) for tri in a.velocity]
mesh_wr = vtkWriter("endo_mesh", attrib=Dict("velocity"=>mesh_velocity))

# flow writer
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array);
# vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_lambda(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I,a.flow.u); a.flow.σ |> Array);
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, sim_time(a)); a.flow.σ |> Array);
new_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "λ₂"=>vtk_lambda, "ω₃"=>vtk_ω)
wr = vtkWriter("endo_flow"; attrib=new_attrib)

# some time benchmarking
update_time = []
measure_time = []
step_time = []
save_time = []

# integrate in time, one convective time in one cycle
println("Run for: ", dt*length(motion))
for (i,t) in enumerate(range(0,4*dt*length(motion),step=dt))
    # now update the mesh, this sets the velocity BC
    tic = time()
    reference = motion[(i-1)%length(motion)+1] # get the reference position
    points = [Point3f(pnt) for pnt in eachcol(reference)]
    WaterLilyPreCICE.update!(mesh_body, points, dt)
    toc = time()
    push!(update_time,toc-tic)

    tic = time()
    measure!(sim) # update the body on the flow
    toc = time()
    push!(measure_time,toc-tic)

    # update flow until t
    tic = time()
    sim_step!(sim,t/sim.L)
    # push!(sim.flow.Δt,dt)
    toc = time()
    push!(step_time,toc-tic)

    tic = time()
    # measure velocity and save
    save!(mesh_wr,mesh_body,sim_time(sim))
    save!(wr,sim)
    toc = time()
    push!(save_time,toc-tic)

    println("tU/L=",round(sim_time(sim),digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(mesh_wr); close(wr)
# using Plots
# plot([update_time, measure_time, step_time, save_time],
# title="Time per step", xlabel="Step", ylabel="Time (s)",
# label=["Update" "Measure" "Step" "Save"])