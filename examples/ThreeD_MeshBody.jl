using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_sphere(;L=32,Re=250,U=1)
    # move the geometry to the center of the domain
    map(x,t) = x .+ SA[L,L,L/2]
    # map(x,t) = x .- SA[L,L+L/4*sin(2π*t/L),L/2]
    # @TODO this has a moving body if you comment this line
    # you need to change the set the flag remeasure to true in the sim_step! function

    # make the body from the stl mesh
    body = MeshBody(joinpath(@__DIR__,"../meshes/sphere.stl");map,scale=L/2)
    # body from a S4 mesh
    # body = MeshBody(joinpath(@__DIR__,"../meshes/cube_S4.inp");map,scale=L/2)
    # generate sim
    Simulation((2L,2L,L), (U,0,0), L; ν=U*L/Re, body)
    # body
end


# make a writer with some attributes to output to the file
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
vtk_vbody(a::Simulation) = a.flow.V |> Array;
vtk_mu0(a::Simulation) = a.flow.μ₀ |> Array;
vtk_srf(a::MeshBody) = Float32[el[1] for el in a.surf_id]
vtk_normal(a::MeshBody) = [WaterLilyPreCICE.normal(el) for el in a.mesh]


# make the sim
sim = make_sphere(L=32)
# make the paraview writer
wr = vtkWriter("MeshBody";attrib=Dict("u"=>vtk_velocity,"p"=>vtk_pressure,"d"=>vtk_body,"v"=>vtk_vbody,"μ₀"=>vtk_mu0))
wr_mesh = vtkWriter("Airfoil_mesh";attrib=Dict("srf"=>vtk_srf,"normal"=>vtk_normal))
# duration and write steps
t₀,duration,step = 0.,1,0.1
# run the sim
@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=true)
    save!(wr,sim); save!(wr_mesh,sim.body,sim_time(sim))
    f = WaterLilyPreCICE.forces(sim.body,sim.flow)
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr); close(wr_mesh)