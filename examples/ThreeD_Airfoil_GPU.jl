using WaterLilyPreCICE,StaticArrays,WriteVTK

function make_airfoil(;L=32,Re=1000,U=1,mem=Array)

    # make the body from the stl mesh
    body = MeshBody(joinpath(@__DIR__,"../meshes/naca/naca.inp");scale=L,
                    map=(x,t)->x.+SA[L/2,L,0])

    # generate sim without a body
    sim = Simulation((3L,2L,L÷2), (U,0,0), L; ν=U*L/Re, mem)

    # make arrays to copy
    σ = Array(sim.flow.σ)
    V = Array(sim.flow.V)
    μ₀ = Array(sim.flow.μ₀)
    μ₁ = Array(sim.flow.μ₁)

    # measure on the CPU
    WaterLilyPreCICE.measure_CPU!(μ₀,μ₁,V,σ,body)

    # put back on the GPU
    copyto!(sim.flow.σ, σ)
    copyto!(sim.flow.V, V)
    copyto!(sim.flow.μ₀, μ₀)
    copyto!(sim.flow.μ₁, μ₁)

    # re-initialise the pressure solver
    WaterLily.update!(sim.pois)

    return sim,body
end

# make a writer with some attributes to output to the file
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_mu0(a::Simulation) = a.flow.μ₀ |> Array;
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "μ₀"=>vtk_mu0)
# vtk attributes for the MeshBody
vtk_srf(a::MeshBody) = Float32[el[1] for el in a.surf_id]
custom = Dict("srf"=>vtk_srf)

# make the sim
using CUDA
sim,body = make_airfoil(L=64;mem=CuArray)

# make the paraview writer
wr = vtkWriter("Airfoil";attrib=custom_attrib)
wr_mesh = vtkWriter("Airfoil_mesh";attrib=custom)

# duration and write steps
t₀,duration,step = 0.,10.0,0.2
forces = []
# run the sim
@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=false)
    save!(wr, sim); save!(wr_mesh, body, sim_time(sim))
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr); close(wr_mesh)