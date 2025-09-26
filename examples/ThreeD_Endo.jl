"""
TODO:
  - [x] check velocity of center of elements
  - [x] use exact triangulation
  - [x] include in WaterLily simulation
  - [x] correct scale and motion of top
  - [ ] measure on the GPU?
  - [ ] add cap and valves
  - [ ] include flow rates from 0D model
  - [ ] interpolate mesh in time?
"""

using WaterLilyPreCICE,StaticArrays
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.5)

function make_sim(;L=64,Re=5000,U=1,mem=Array,T=Float32)
    # extract the good stuff
    data = open("/home/marin/Workspace/WaterLilyPreCICE/meshes/endo/Endo_Kasra_xyz.txt", "r") do file
        data = mapreduce(L->parse.(T, split(L)), hcat, readlines(file)[10:end])
    end
    # remove the reference height
    zmin,zmax = extrema(data[3,:]); scale = zmax-zmin
    # remove 3 first columns since these are the ref config, make it unit-length and scale
    motion = [data[1+i*3:3+i*3,:]./scale.*L for i in 1:size(data)[1]÷3-1]
    # the top surface moves up and down, we must fix it
    motion = [m.+SA{T}[L/2,L/2,3L/2-maximum(m[3,:])] for m in motion]

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
    # body = MeshBody(mesh)

    # add a top surface
    # # body = body ∩ AutoBody((x,t)->x[3]-3.f0L÷2.f0)
    R = T(L/3.f0)
    @fastmath r(x) = √(x[1]^2+x[2]^2)
    # cap = AutoBody((x,t)->√sum(abs2,x.-SA[clamp(x[1],-R,R),clamp(x[2],-R,R),clamp(x[3],-R,R)])-R+0.f0,
    #                (x,t)-> x.-SA{T}[L/2.f0,L/2.f0,3.f0L÷2.f0+2.f0])
    # cap -= AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),L/8.f0))^2)-4.0f0,
    #                 (x,t)->x.-SA{T}[L/2.f0,L/2.f0,3.f0L÷2.f0])
    # body += cap
    function box(x,t)
        q = abs.(x).-SA{T}[1.3R,1.3R,2R]
        return √sum(abs2,max.(q,0.f0)) + min(max(q[1],max(q[2],q[3])),0.f00);
    end
    box = AutoBody(box,(x,t)-> x.-SA{T}[L/2.f0,L/2.f0,L])

    body = box -  AutoBody((x,t)->√(max(abs(x[3])-L/4.f0,0)^2+(r(x)-min(r(x),L/8.f0))^2),
                     (x,t)->x.-SA{T}[L/2.f0,L/2.f0,3.f0L÷2.f0])
    body -= MeshBody(mesh) # remove the body from the box

    # make a simulation
    Simulation((L,L÷2,2L),(0,0,0),L;U,ν=U*L/Re,body,mem,T), motion
end

# make sim
sim,motion = make_sim(L=32);
mesh_body = sim.body.b.a;
# mesh_body = sim.body

# find max displacement, and max velocity
u = 0.5
dx = [√maximum(sum(abs2,motion[end-i] .- motion[end-i-1], dims=1)) for i in 0:99]
dt = dx./u
println("Maximum dt: ",maximum(dt),", minimum dt=",minimum(dt))
dt = maximum(dt) # use the maximum dt

# get initial bounding box
# bbox = GeometryBasics.bounding_box(body.mesh)

# velocity of center of elements
mesh_velocity(a::MeshBody) = [WaterLilyPreCICE.center(tri) for tri in a.velocity]
mesh_wr = vtkWriter("endo_mesh", attrib=Dict("velocity"=>mesh_velocity))

# flow writer
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, sim_time(a)); a.flow.σ |> Array);
new_attrib = Dict("Velocity"=>vtk_velocity, "Pressure"=>vtk_pressure, "dist"=>vtk_body, "Vbody"=>vtk_vbody)
wr = vtkWriter("endo_flow"; attrib=new_attrib)

# integrate in time, one convective time isone cycle
# Np = 4 # number of period
# dt = Float32(sim.L/length(motion)/Np) # time step in CFL
println("Run for: ", dt*length(motion))
for (i,t) in enumerate(range(0,dt*length(motion),step=dt))
    # now update the mesh, this sets the velocity BC
    reference = motion[(i-1)%length(motion)+1] # get the reference position
    points = [Point3f(pnt) for pnt in eachcol(reference)]
    WaterLilyPreCICE.update!(mesh_body,points,dt)
    measure!(sim) # update the body on the flow

    # update flow until t
    sim_step!(sim,t/sim.L)

    # measure velocity and save
    save!(mesh_wr,mesh_body,sim_time(sim))
    save!(wr,sim)

    println("tU/L=",round(sim_time(sim),digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(mesh_wr); close(wr)