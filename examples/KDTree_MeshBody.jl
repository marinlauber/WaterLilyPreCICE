using WaterLily
using WaterLilyPreCICE
using StaticArrays

# include("../src/WaterLily/MeshBodies.jl")
# path = "/home/marin/Workspace/CardioVascularFlow.jl/example/"
function make_sphere(;L=32,Re=250,U=1)
    # make a body
    map(x,t) = x .- SA[L,L,L/2]
    # body = MeshBody("cube.stl";map,scale=L/4)
    body = MeshBody("cube.inp";map,scale=L/6)
    # body = MeshBody(path*"sphere.inp";map,scale=L/6)
    # @assert all(volume(body) .≈ (L/3)^3) # the volume is exact here!
    # body = MeshBody(path*"3D_flap.inp";map,boundary=false,thk=2,scale=1.0)
    # generate sim
    Simulation((2L,2L,L), (U,0,0), L; ν=U*L/Re, body)
end

# bb = AutoBody((x,t)->√sum(abs2,x)-32)

# make a writer with some attributes
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); 
                        a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;
mu0(a::Simulation) = a.flow.μ₀ |> Array;
norm(a::Simulation) = (a.flow.f .= 0;
                      @WaterLily.loop a.flow.f[I,:] .= measure(a.body,loc(0,I),0.0)[2] over I ∈ WaterLily.inside(a.flow.p);
                      a.flow.f |> Array;)

custom_attrib = Dict(
    "u" => velocity, "p" => pressure, "n" => norm,
    "d" => _body, "v" => _vbody, "μ₀" => mu0,
)# this

# run a sim and plot the time evolution
sim = make_sphere(L=64)

# using BenchmarkTools
body = sim.body;
mesh = body.mesh;
x = SA[20.f0,34.5f0,16.5f0]
t = 0.f0

# function test_measure(body,x,t)
#     d,n = measure(body.mesh,x,t)
# end

# function WaterLily.measure(body::MeshBody,x,t,;fastd²=Inf)
#     # eval d=map(x,t)-x, and n̂
#     # ξ = body.map(x,t)
#     d = zero(eltype(x))
#     n = zero(x)
#     # # #  if we are outside of the bounding box, we can measure approx
#     # bbox = body.bbox;
#     # outside(ξ,bbox) #&& return (dist(ξ,body.bbox),zero(x),zero(x))
#     # d,n = measure(body.mesh,ξ,t)
#     !body.boundary && (d = abs(d)-body.half_thk) # if the mesh is not a boundary, we need to adjust the distance

#     # The velocity depends on the material change of ξ=m(x,t):
#     #   Dm/Dt=0 → ṁ + (dm/dx)ẋ = 0 ∴  ẋ =-(dm/dx)\ṁ
#     J = ForwardDiff.jacobian(x->body.map(x,t), x)
#     dot = ForwardDiff.derivative(t->body.map(x,t), t)
#     return (d,n,-J\dot)
# end

import WaterLilyPreCICE: AbstractBox,Bbox,center,bounding_box
WaterLily.inside(tri::GeometryBasics.Ngon{3},b::Bbox)::Bool = (x=center(tri);(all(b.C-b.R .≤ x) && all(x .≤ b.C+b.R)))
Base.filter(mesh::GeometryBasics.Mesh,b::Bbox) = findall(x->inside(x,b),mesh)

# make a box that holds the mesh
box = bounding_box(mesh)
@assert inside(mesh[1],box) # it must be inside
@assert all(filter(mesh,box) .≈ 1:length(mesh)) # in fact, they must all be inside
# now we split that box in two
left,right = split(box)
inleft = filter(mesh,left)
inright = filter(mesh,right)

vleft = @view(mesh[inleft])

lleft,lright=split(left)
inlleft = inside(vleft,lleft)

function Tree(mesh::GeometryBasics.Mesh{3,T}; δ=4, Nmax=64) where T
    b = bounding_box(mesh,δ=δ) 
    boxes,idx = AbstractBox[b],Int[1]
    # split the box in the longest direction
    left, right = split(b)
    # find which box each point belongs to
    in_left = filter(mesh, left)
    in_right = filter(mesh, right) # we can avoid this actually
    # downsample the boxes
    tree!(boxes, idx, 2, view(mesh[in_left]), left, Nmax)
    tree!(boxes, idx, 3, view(mesh[in_right]), right, Nmax)
    # sort the boxes by index and return, we need to fill some indices empty boxes
    # as we then access child of boxes[I] as boxes[2I] and boxes[2I+1]
    boxes = ntuple(i->i in idx ? boxes[findfirst(==(i),idx)] : NoBbox(T), maximum(idx))
    Tree(boxes,mesh,boxes[1].C,boxes[1].R)
end
function tree!(list, idx, Is, sub_mesh, b::Bbox, Nmax) where T
    push!(idx,Is) # add the index location
    # if we are on a leaf, we push a leaf box on the list
    if (length(sub_mesh)<Nmax)
        push!(list, Bbox(b.C,1.2*b.R,true,sub_mesh))
        return nothing
    end
    push!(list, b);
    # we are not on a leaf node, we can downsample
    left, right = split(b)
    # find which box each point belongs to
    println()
    in_left = filter(sub_mesh, left)
    in_right = filter(sub_mesh, right)
    # reshape in case points are clustered on one side
    left = bounding_box(@view(sub_mesh[in_left]))
    right = bounding_box(@view(sub_mesh[in_right]))
    # downsample the boxes
    tree!(list, idx, Is*2, @view(sub_mesh[in_left]), left, Nmax)
    tree!(list, idx, Is*2+1, @view(sub_mesh[in_right]), right, Nmax)
    return nothing
end

# make a Tree for the mesh
tree = Tree(mesh;Nmax=10)

GeometryBasics.view(faces(mesh),2)



wr = vtkWriter("CalculiX_test"; attrib=custom_attrib)
t₀,duration,step = 0.,10,0.1
write!(wr,sim)
println("tU/L=",round(t₀,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# @time for tᵢ in range(t₀,t₀+duration;step)
#     # update until time tᵢ in the background
#     sim_step!(sim,tᵢ;remeasure=false)
#     write!(wr,sim)
    
#     # print time step
#     println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# end
close(wr)