using WaterLily
using StaticArrays
using GeometryBasics
using BenchmarkTools
using CUDA
using Plots
using Adapt

function load_inp(fname; facetype=GLTriangleFace, pointtype=Point3f)
    #INP file format
    @assert endswith(fname,".inp") "file type not supported"
    fs = open(fname)

    points = pointtype[]
    faces = facetype[]
    node_idx = Int[]
    srf_id = Tuple[]
    cnt = 0

    # read the first 3 lines if there is the "*heading" keyword
    line = readline(fs)
    contains(line,"*heading") && (line = readline(fs))
    BlockType = contains(line,"*NODE") ? Val{:NodeBlock}() : Val{:DataBlock}()

    # read the file
    while !eof(fs)
        line = readline(fs)
        contains(line,"*ELSET, ELSET=") && (cnt+=1)
        BlockType, line = parse_blocktype!(BlockType, fs, line)
        if BlockType == Val{:NodeBlock}()
            push!(node_idx, parse(Int,split(line,",")[1])) # keep track of the node index of the inp file
            push!(points, pointtype(parse.(eltype(pointtype),split(line,",")[2:4])))
        elseif BlockType == Val{:ElementBlock}()
            nodes = parse.(Int,split(line,",")[2:end])
            push!(faces, TriangleFace{Int}(facetype([findfirst(==(node),node_idx) for node in nodes])...)) # parse the face
        elseif BlockType == Val{:ElSetBlock}()
            push!(srf_id, (cnt, parse.(Int,split(line,",")[1])));
        else
            continue
        end
    end
    close(fs) # close file stream
    # case where there is no surface ID and we pass it a single tuple of all the faces
    srf_id = length(srf_id)==0 ? ntuple(i->(1,i),length(faces)) : ntuple(i->(srf_id[i].+(0,1-srf_id[1][2])),length(srf_id))
    return Mesh(points, faces), srf_id
end
function parse_blocktype!(block, io, line)
    contains(line,"*NODE") && return block=Val{:NodeBlock}(),readline(io)
    contains(line,"*ELEMENT") && return block=Val{:ElementBlock}(),readline(io)
    contains(line,"*ELSET, ELSET=") && return block=Val{:ElSetBlock}(),readline(io)
    return block, line
end

locate(tri::SMatrix{T},p::SVector{T}) where T = locate(tri[:,1],tri[:,2],tri[:,3],p)
function locate(a,b,c,p::SVector{T}) where T
    ab = b.-a
    ac = c.-a
    ap = p.-a
    d1 = sum(ab.*ap)
    d2 = sum(ac.*ap)
    # is point `a` closest?
    if ((d1 ≤ 0) && (d2 ≤ 0))
        return a
    end
    # is point `b` closest?
    bp = p.-b
    d3 = sum(ab.*bp)
    d4 = sum(ac.*bp)
    if ((d3 ≥ 0) && (d4 ≤ d3))
        return b
    end
    # is point `c` closest?
    cp = p.-c
    d5 = sum(ab.*cp)
    d6 = sum(ac.*cp)
    if ((d6 ≥ 0) && (d5 ≤ d6))
        return c
    end
    # is segment 'ab' closest?
    vc = d1*d4 - d3*d2
    if ((vc ≤ 0) && (d1 ≥ 0) && (d3 ≤ 0))
        x =  a .+ ab.*d1 ./ (d1 - d3)
        return x
    end
    #  is segment 'ac' closest?
    vb = d5*d2 - d1*d6
    if ((vb ≤ 0) && (d2 ≥ 0) && (d6 ≤ 0))
        x =  a .+ ac.*d2 ./ (d2 - d6)
        return x
    end
    # is segment 'bc' closest?
    va = d3*d6 - d5*d4
    if ((va ≤ 0) && (d4 ≥ d3) && (d5 ≥ d6))
        x =  b .+ (c .- b) .* (d4 - d3) ./ ((d4 - d3) + (d5 - d6))
        return x
    end
    # closest is interior to `abc`
    denom = one(T) / (va + vb + vc)
    v= vb*denom
    w = vc*denom
    x = a .+ ab .* v .+ ac .* w
    return x
end

# closest triangle index in the mesh
function closest(mesh,x::SVector{T};kwargs...) where T
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    return u
end

# struct MeshBody{T,A<:AbstractArray,S<:AbstractVector{T}} <: AbstractBody
struct MeshBody{T,A<:AbstractArray,S<:AbstractVector{T},F<:Function} <: AbstractBody
    mesh :: A
    origin :: S
    width :: S
    boundary :: Bool
    map   :: F
    scale :: T
    half_thk :: T
    function MeshBody(mesh,origin::SVector{3,T},width::SVector{3,T},boundary=true,
                        map=(x,t)->x,scale=1,half_thk=0) where T
        new{T,typeof(mesh),typeof(origin),typeof(map)}(mesh,origin,width,boundary,map,T(scale),T(half_thk))
    end
end
Adapt.adapt_structure(to, x::MeshBody) = MeshBody(adapt(to, x.mesh), adapt(to, x.origin),
                                                  adapt(to, x.width), x.boundary, adapt(to, x.map),
                                                  x.scale, x.half_thk)

function MeshBody(file_name; map=(x,t)->x,boundary=true,half_thk=0,scale=1,mem=Array,T=Float32)
    # read in the mesh
    mesh = endswith(file_name,".inp") ? load_inp(file_name)[1] : load(file_name)
    points = Point3f[] # scale and map the points to the correct location
    for pnt in mesh.position
        push!(points, Point3f(SA{T}[pnt.data...]*T(scale)))
    end
    mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(mesh))
    width = SVector{3,T}(maximum(mesh.position,dims=1)...)
    origin = SVector{3,T}(minimum(mesh.position,dims=1)...)
    mesh = [hcat(vcat([mesh[i]...])...) for i in 1:length(mesh)] |> mem
    MeshBody(mesh,origin.-2,width-origin.+4,boundary,map,T(scale),T(half_thk))
end

WaterLily.sdf(body::MeshBody,x,t;kwargs...) = measure(body,x,t;kwargs...)[1]

using ForwardDiff
function WaterLily.measure(body::MeshBody,x::SVector{D,T},t;kwargs...) where {D,T}
    # map to correct location
    ξ = body.map(x,t)
    # we don't need to worry if the geom is a boundary or not
    outside(ξ,body.origin,body.width) && return (T(4),zeros(SVector{D,T}),zeros(SVector{D,T}))
    # locate the point on the mesh
    #TODO this is what we replace with the BVH
    u = closest(body.mesh,ξ;kwargs...)
    # compute the normal and distance
    n,p = normal(body.mesh[u]),SVector(locate(body.mesh[u],x))
    # signed Euclidian distance
    s = ξ-p; d = sign(sum(s.*n))*√sum(abs2,s)
    # velocity at the mesh point
    J = ForwardDiff.jacobian(ξ->body.map(ξ,t), ξ)
    dot = ForwardDiff.derivative(t->body.map(ξ,t), t)
    # if the mesh is not a boundary, we need to adjust the distance
    !body.boundary && (d = abs(d)-body.half_thk)
    return (d,n,-J\dot)
end
outside(x::SVector,origin,width) = !(all(origin .≤ x) && all(x .≤ origin+width))

using LinearAlgebra: cross
@inbounds @inline d²_fast(tri::SMatrix,x::SVector) = sum(abs2,x-SVector(sum(tri,dims=2)/3))
@inbounds @inline normal(tri::SMatrix) = hat(SVector(cross(tri[:,2]-tri[:,1],tri[:,3]-tri[:,1])))
@inbounds @inline hat(v) = v/(√(v'*v)+eps(eltype(v)))

# access the WaterLily writer to save the file
function WaterLily.save!(w,a::MeshBody,t=w.count[1]) #where S<:AbstractSimulation{A,B,C,D,MeshBody}
    k = w.count[1]
    points = hcat([[p.data...] for p ∈ a.mesh.position]...)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, Base.to_index.(face)) for face in faces(a.mesh)]
    vtk = vtk_grid(w.dir_name*@sprintf("/%s_%06i", w.fname, k), points, cells)
    for (name,func) in w.output_attrib
        # point/vector data must be oriented in the same way as the mesh
        vtk[name] = ndims(func(a))==1 ? func(a) : permutedims(func(a))
    end
    vtk_save(vtk); w.count[1]=k+1
    w.collection[round(t,digits=4)]=vtk
end


function make(L;U=1,mem=CuArray,T=Float32)
    α = π/10.f0 # rotation angle
    function map(x,t)
        Rx = SA[1 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]
        Ry = SA[cos(α) 0 sin(α); 0 1 0; -sin(α) 0 cos(α)]
        Rz = SA[cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
        Rx*Ry*Rz*(x.-L/2.f0).+0.5f0
    end

    body = MeshBody(joinpath(@__DIR__,"meshes/cube.inp");scale=L/2,map,mem)
    # Simulation((L,L,L),(U,0,0),L;body,mem,T)
end

L = 32
MEMORY = CuArray
body = make(L;mem=MEMORY)

distance = zeros(L+2,L+2,L+2)|> MEMORY;
@inside distance[I] = sdf(body,loc(0,I))
flood(distance[:,L÷2+1,:]); contour!(Array(distance)[:,L÷2+1,:]',levels=[0],lw=2)




# import WaterLilyPreCICE: @loop
# @loop new_body.mesh[I] = SMatrix{3,3,Float32}(rand(3,3)) over I in CartesianIndices(new_body.mesh)