"""
TODO:
  - [x] check velocity of center of elements
  - [ ] use exact triangulation
  - [ ] include in WaterLily simulation
  - [x] correct scale and motion of top
  - [ ] measure on the GPU?
  - [ ] add cap and valves
  - [ ] include flow rates from 0D model
  - [ ] interpolate mesh in time?
"""

using WaterLilyPreCICE,DelaunayTriangulation,StaticArrays,GLMakie


function make_sim(;L=64,Re=5000,U=1,mem=Array,T=Float32)
    file = open("/home/marin/Workspace/WaterLilyPreCICE/meshes/Endo_Kasra_xyz.txt", "r")
    lines = readlines(file)[10:end] # move header of the file
    # extract the good stuff
    data = mapreduce(L->parse.(T, split(L)), hcat, lines)
    # data[:,1:3:end] .+= z_ref # remove the reference height
    zmin,zmax = extrema(data[3,:]); scale = zmax-zmin
    # remove 3 first columns since these are the ref config, make it unit-length and scale
    motion = [data[1+i*3:3+i*3,:]./scale.*L for i in 1:size(data)[1]÷3-1]
    # the top surface moves up and down, we must fix it
    motion = [m.+[L/2,L/2,3L/2-maximum(m[3,:])] for m in motion]

    # triangulate, select one snapshot that is easy to triangulate
    tris = triangulate(motion[39],  predicates=ExactKernel())
    tris = delete_ghost_triangles!(tris)

    # make points list from first snapshot
    points = [Point3f(pnt) for pnt in eachcol(motion[1])]

    # initial mesh, create triangulation
    faces = [TriangleFace{Int}(GLTriangleFace(tri)) for tri in get_triangles(tris)]

    # generate a mesh from this
    mesh = GeometryBasics.Mesh(points,faces)
    body = MeshBody(mesh)

    # make a simulation
    Simulation((L,L,2L),(0,0,0),L;U,ν=U*L/Re,body,mem,T), motion
end

# make sim
sim,motion = make_sim()
body = sim.body;

# velocity of center of elements
mesh_velocity(a::MeshBody) = [WaterLilyPreCICE.center(tri) for tri in a.velocity]
mesh_wr = vtkWriter("endo_mesh", attrib=Dict("velocity"=>mesh_velocity))

# flow writer
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_vbody(a::AbstractSimulation) = a.flow.p |> Array;
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, sim_time(a)); a.flow.σ |> Array);
new_attrib = Dict("Velocity"=>vtk_velocity, "Pressure"=>vtk_pressure, "dist"=>vtk_body, "Vbody"=>vtk_vbody)
wr = vtkWriter("endo_flow"; attrib=new_attrib)

# integrate in time
for reference in motion
    push!(sim.flow.Δt,sim.L)
    points = [Point3f(pnt) for pnt in eachcol(reference)]
    WaterLilyPreCICE.update!(body,points,1)
    save!(mesh_wr,body,sim_time(sim)); save!(wr,sim)
end
close(mesh_wr); close(wr)