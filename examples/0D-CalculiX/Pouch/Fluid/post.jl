using JLD2,Plots,StaticArrays
using WaterLilyPreCICE


function Elastance(t;a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)
end
function Gaussian(t;μ=0,σ=0.1)
    exp(-((t%1 - μ) / σ)^2 / 2)
end


A4(ts) = ts<=5 ? 0 : 2Gaussian((ts-5)/20;μ=0.35,σ=0.08) #Elastance((ts-10)/20)


data = jldopen("/home/marin/Workspace/WaterLilyPreCICE/examples/0D-CalculiX/Pouch/Fluid/pouch_volume.jld2")

begin
    ti = data["time"]
    Plots.plot(ti,data["intV"][1:end-1],label="V_v", ylim=(0.,0.25))
    Plots.plot!(ti,data["intP"][1:end-1],label="P_v")
    Plots.plot!(ti,data["q"],label="Q")
    Plots.plot!(0:0.01:14,0.21*A4.(0:0.01:14)/2,label="P_act/2")
end

# time = collect(0:0.1:100)

# Plots.plot(time, Elastance.(time/20))
# Plots.plot!(time, Gaussian.(time/20;μ=0.35,σ=0.08))

# Plots.plot(data["time"],data["q"], ylims=(0,0.05))
# Plots.plot(data["intdt"])
# Plots.plot(data["intP"])
# Plots.plot(data["intV"])

# Plots.plot!(time.+10, u1,label="u1")
# Plots.plot!(time.+10, u2,label="u2")
# function sub_mesh(b, srf)
#     points = Point3f[] # scale and map the points to the correct location
#     faaces = TriangleFace[]
#     for (i,pnt) in enumerate(b.mesh.position)
#         push!(points, Point3f(SA[pnt.data...]))
#     end
#     for (i,face) in enumerate(faces(b.mesh))
#         if getindex(b.surf_id[i],1) ∈ srf
#         end
#     end 
#     GeometryBasics.Mesh(points,faaces)
# end

# function volume(b::MeshBody, srf)
#     vol = zeros(SVector{3,Float64})
#     for (i,T) in enumerate(b.mesh)
#         id = getindex(b.surf_id[i],1)
#         sgn = ifelse(id ≤ getindex(maximum(body.surf_id),1)÷2, 1, -1)
#         id ∈ srf && (vol += sgn .* WaterLilyPreCICE.center(T) .* WaterLilyPreCICE.dS(T))
#     end
#     return vol
# end



# # volume(a::GeometryBasics.Mesh) = mapreduce(T->center(T).*dS(T),+,a)
# # volume(body::MeshBody) = volume(body.mesh)

# body = MeshBody("/home/marin/Workspace/WaterLilyPreCICE/examples/0D-CalculiX/Pouch/Solid/deformed.inp");
# using GLMakie
# GLMakie.mesh(body.mesh)

# v_inner = volume(body, [1,2,3,6,7,8])
# v_outer = volume(body, [1,4,5,6,9,10])

# mesh_srf_inner = sub_mesh(body, [1,2,3,6,7,8])
# mesh_srf_outer = sub_mesh(body, [1,4,5,6,9,10])

# GLMakie.mesh(mesh_srf_inner)
# GLMakie.mesh(mesh_srf_outer)