using WaterLilyPreCICE,StaticArrays
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.5)

# find a body of type T in a sim
# find(a::WaterLily.AbstractBody, ::Type{T}) where T = isa(a, T) ? a : nothing
function find(a::WaterLily.SetBody, ::Type{T}) where T
    isa(a.a, T) && return a.a
    isa(a.b, T) && return a.b
    isa(a.a, WaterLily.SetBody) && return find(a.a, T)
    isa(a.b, WaterLily.SetBody) && return find(a.b, T)
end

function make_sim(;L=64,Re=5000,U=1,mem=Array,T=Float32)
    # extract the good stuff
    data = open("/home/marin/Workspace/WaterLilyPreCICE/meshes/limo/nodal_position.txt", "r") do file
        mapreduce(L->parse.(T, split(L)), hcat, readlines(file))
    end
    # remove the reference height
    ymin,ymax = extrema(data[2,2:end]); scale = ymax-ymin
    # remove first row that contains the time and then scale and move to center of domain
    motion = [data[1+i*3:3+i*3,2:end]./scale.*L.+SA{T}[3L/4,L/4,L/2] for i in 0:size(data,1)÷3-1]
    # get the times
    times = data[1:3:end,1]
    # read it from the connectivity file
    connectivity = open("/home/marin/Workspace/WaterLilyPreCICE/meshes/limo/connectivity.txt", "r") do f
        mapreduce(L->parse.(Int, split(L)), hcat, readlines(f))
    end
    # make points list from first snapshot
    points = [Point3f(pnt) for pnt in eachcol(motion[1])]
    # initial mesh, create triangulation
    # connectivity must be incremented since we have 1 indexing here
    faces = [TriangleFace{Int}(GLTriangleFace((reverse(tri).+1...))) for tri in eachcol(connectivity)]

    # generate a mesh from this
    mesh = GeometryBasics.Mesh(points,faces)
    shell = MeshBody(mesh;thk=2f0,boundary=false)

     # add a top surface
    R₁,R₂,Λ = T(L*0.35),T(L/12),T(2.0)
    println("R₁=",R₁," R₂=",R₂,", Λ=",Λ)
    @fastmath r(x,Λ=1) = √((x[1]/Λ)^2+x[3]^2)
    cap = AutoBody((x,t)->√(x[2]^2+(r(x,Λ)-min(r(x,Λ),R₁/Λ))^2)-2.0f0,
                   (x,t)->x.-SA{T}[3L/4.f0,L/4.f0,L/2.f0])
    cap -= AutoBody((x,t)->√(x[2]^2+(r(x)-min(r(x),R₂))^2)-4.0f0,
                    (x,t)->x.-SA{T}[3L/4.f0,L/4.f0,L/2.f0])

    # ful body
    body = shell + cap

    # make a simulation
    return Simulation((3L÷2,3L÷2,L),(0,0,0),L;U,ν=U*L/Re,body,mem,T), motion, times
end

# # make simulation
sim,motion,times = make_sim(L=128);
mesh_body = find(sim.body, MeshBody)
# mesh_body = sim.body

# # find max displacement, and max velocity
u = 0.5
dx = [√maximum(sum(abs2,motion[end-i] .- motion[end-i-1], dims=1)) for i in 0:length(motion)-2]
dt = dx./u
println("Maximum dt: ",maximum(dt),", minimum dt=",minimum(dt))
dt = maximum(dt) # use the maximum dt

# velocity of center of elements
mesh_velocity(a::MeshBody) = [WaterLilyPreCICE.center(tri) for tri in a.velocity]
mesh_wr = vtkWriter("limo_mesh", attrib=Dict("velocity"=>mesh_velocity))

# flow writer
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array);
# vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_lambda(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I,a.flow.u); a.flow.σ |> Array);
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, sim_time(a)); a.flow.σ |> Array);
new_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "λ₂"=>vtk_lambda, "ω₃"=>vtk_ω)
wr = vtkWriter("limo_flow"; attrib=new_attrib)

# some time benchmarking
update_time = []
measure_time = []
step_time = []
save_time = []
motion = motion[1:4:end]

# integrate in time, one convective time in one cycle
@time for (i,t) in enumerate(range(0,length(motion)-1))
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

    # # update flow until t
    tic = time()
    # sim_step!(sim,t/sim.L)
    push!(sim.flow.Δt,dt)
    toc = time()
    push!(step_time,toc-tic)

    tic = time()
    # measure velocity and save
    save!(mesh_wr,mesh_body,sim_time(sim))
    save!(wr,sim)
    toc = time()
    push!(save_time,toc-tic)

    println("tU/L=",round(sim_time(sim),digits=4),", Δt=",round(sim.flow.Δt[end],digits=3)," %=",round(i/length(motion)*100,digits=1))
end
close(mesh_wr); close(wr)

# # using Plots
# plot([update_time, measure_time, step_time, save_time],
# title="Time per step", xlabel="Step", ylabel="Time (s)",
# label=["Update" "Measure" "Step" "Save"])