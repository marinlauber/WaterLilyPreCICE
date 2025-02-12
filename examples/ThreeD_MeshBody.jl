using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_sphere(;L=32,Re=250,U=1)
    # move the geometry to the center of the domain
    # map(x,t) = x .- SA[L,L,L/2]
    map(x,t) = x .- SA[L,L+L/4*sin(2π*t/L),L/2] 
    # @TODO this has a moving body if you comment this line
    # you need to change the set the flag remeasure to true in the sim_step! function
    
    # make the body from the stl mesh
    body = MeshBody(joinpath(@__DIR__,"../meshes/cube.stl");map,scale=L/2)
    # generate sim
    Simulation((4L,2L,L), (U,0,0), L; ν=U*L/Re, body)
end

# make a writer with some attributes to output to the file
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;
mu0(a::Simulation) = a.flow.μ₀ |> Array;
custom_attrib = Dict("u" => velocity, "p" => pressure, "d" => _body, "v" => _vbody, "μ₀" => mu0,)

# make the sim
sim = make_sphere(L=32)
# make the paraview writer
wr = vtkWriter("STL_mesh_test";attrib=custom_attrib)
# duration and write steps
t₀,duration,step = 0.,10,0.1
# run the sim
@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=true)
    write!(wr,sim)
    f = WaterLilyPreCICE.forces(sim.body,sim.flow)
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)