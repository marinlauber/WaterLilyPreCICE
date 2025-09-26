using WaterLilyPreCICE,StaticArrays,WriteVTK

import WaterLily: @loop
function measure_CPU!(μ₀,μ₁,V,σ::AbstractArray{T,N},body::AbstractBody;t=zero(T),ϵ=1) where {T,N}
    V .= zero(T); μ₀ .= one(T); μ₁ .= zero(T); d²=(2+ϵ)^2
    @fastmath @inline function fill!(μ₀,μ₁,V,d,I)
        d[I] = sdf(body,loc(0,I,T),t,fastd²=d²)
        if d[I]^2<d²
            for i ∈ 1:N
                dᵢ,nᵢ,Vᵢ = measure(body,loc(i,I,T),t,fastd²=d²)
                V[I,i] = Vᵢ[i]
                μ₀[I,i] = WaterLily.μ₀(dᵢ,ϵ)
                for j ∈ 1:N
                    μ₁[I,i,j] = WaterLily.μ₁(dᵢ,ϵ)*nᵢ[j]
                end
            end
        elseif d[I]<zero(T)
            for i ∈ 1:N
                μ₀[I,i] = zero(T)
            end
        end
    end
    @loop fill!(μ₀,μ₁,V,σ,I) over I ∈ inside(σ)
    BC!(μ₀,zeros(SVector{N,T}),false)
    BC!(V ,zeros(SVector{N,T}),false)
end

function make_airfoil(;L=32,Re=1000,St=0.3,U=1,mem=Array)
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
    measure_CPU!(μ₀,μ₁,V,σ,body)

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