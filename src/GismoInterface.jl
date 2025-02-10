using PreCICE

# this lives here until PR:#158 of WaterLily is merged
include("Bodies.jl")

"""
    An structure to hold the coupling data for WaterLily-Gismo coupling
"""
struct GismoInterface <: AbstractInterface
    ControlPointsID::AbstractArray
    ControlPoints::AbstractArray
    quadPointID::AbstractArray
    quadPoint::AbstractArray
    forces::AbstractArray
    deformation::AbstractArray
    knots::AbstractArray
    dt::Vector{Float64}
    dir::Vector{Int16}
    N::Int16
    center::AbstractVector
end

# unpack subarray of inceasing values of an hcat
function unpack(a)
    tmp=[a[1]]; ks=Vector{Number}[]
    for i ∈ 2:length(a)
        if a[i]>=a[i-1]
            push!(tmp,a[i])
        else
            push!(ks,tmp); tmp=[a[i]]
        end
    end
    push!(ks,tmp)
    return ks
end
function knotVectorUnpack(knots)
    knots = reshape(knots,reverse(size(knots)))[1,:]
    unpack(knots)
end

function getControlPoints(points, knots)
    points = reshape(points,reverse(size(points)))
    ncp = [length(knot)-count(i->i==0,knot) for knot in knots]
    @assert sum(ncp) == size(points,2) "Number of control points does not match the number of points"
    return [points[:,1+sum(ncp[1:i])-ncp[1]:sum(ncp[1:i])] for i in 1:length(ncp)]
end
function quadPointUnpack(quadPoints)
    quadPoints = reshape(quadPoints,reverse(size(quadPoints)))
    quadPoints = [filter(!isone,filter(!iszero,quadPoints[:,i]))[1] for i in 1:size(quadPoints,2)]
    unpack(quadPoints)
end
function getDeformation(points,knots)
    ncp = [length(knot)-count(i->i==0,knot) for knot in knots]
    return [points[:,1+sum(ncp[1:i])-ncp[1]:sum(ncp[1:i])] for i in 1:length(ncp)]
end

using ParametricBodies: _pforce, _vforce
"""
    getInterfaceForces!(forces,flow::Flow,body::AbsBodies,quadPoints)

Compute the interface forces at the quadrature points `quadPoints` for each body in `body.bodies`.
"""
Index(Qs,i) = sum(length.(Qs[1:i-1]))+1:sum(length.(Qs[1:i]))
using ParametricBodies: _pforce, _vforce
function getInterfaceForces!(interface::GismoInterface,flow::Flow{T},body::AbsBodies) where T
    for (i,b) in enumerate(body.bodies[1:length(interface.quadPoint)]) # only select the active curves
        I = Index(interface.quadPoint,i)
        fi = reduce(hcat,[-1.0*_pforce(b.curve,flow.p,s,zero(T),Val{false}()) for s ∈ interface.quadPoint[i]])
        interface.dir[i] != 1 && (fi = reverse(fi;dims=2))
        interface.forces[:,I] .= fi
    end
end
import ParametricBodies: NurbsCurve, DynamicNurbsBody
function GismoInterface(; KnotMesh="KnotMesh",ControlPointMesh="ControlPointMesh",
                          ForceMesh="ForceMesh",dir=nothing,curves=nothing,center=SA[0.,0.])

    # initilise PreCICE
    PreCICE.initialize()

    # get the mesh verticies from the fluid solver
    (_, knots) = getMeshVertexIDsAndCoordinates(KnotMesh)
    knots = knotVectorUnpack(knots)
   
    # get the mesh verticies from the structural solver
    (ControlPointsID, ControlPoints) = getMeshVertexIDsAndCoordinates(ControlPointMesh)
    ControlPointsID = Array{Int32}(ControlPointsID)
    ControlPoints = getControlPoints(ControlPoints, knots)
    deformation = copy(ControlPoints)
    isnothing(dir) ? (direction = ones(length(ControlPoints))) : direction = dir
    
    # get the quad points in parameter space
    (quadPointID, quadPoint) = getMeshVertexIDsAndCoordinates(ForceMesh)
    forces = zeros(reverse(size(quadPoint))...)
    quadPointID = Array{Int32}(quadPointID)
    quadPoint = quadPointUnpack(quadPoint)

    dt = PreCICE.getMaxTimeStepSize()
    
    # construct the interface curves
    bodies = AbstractBody[]; ops = Function[]
    for (i,(cps,knot)) in enumerate(zip(ControlPoints,knots))
        direction[i] != 1 && (cps = reverse(cps;dims=2))
        cps = SMatrix{2,size(cps,2)}(cps)
        knot = SVector{length(knot)}(knot)
        weights = SA[ones(size(cps,2))...]
        push!(bodies,DynamicNurbsBody(NurbsCurve(cps.+center,knot,weights)))
        push!(ops, ∩) # always interset with the next curve
    end

    # return coupling interface
    interface = GismoInterface(ControlPointsID, ControlPoints, quadPointID, quadPoint, forces,
                                  deformation, knots, [dt], direction, length(bodies), center)
    
    # add some passive curves if we want
    !isnothing(curves) && for crv in curves
        push!(bodies,crv); push!(ops, ∪) # always union with the next curve
        println("Adding a curve to the stack...")
    end

    # return the interface and the body
    return interface, AbsBodies(bodies, ops)
end

function readData!(interface::GismoInterface)
    # Read control point displacements
    readData = transpose(PreCICE.readData("ControlPointMesh", "ControlPointData", 
                                          interface.ControlPointsID, interface.dt[end]))
    interface.deformation .= getDeformation(readData,interface.knots) # repack correctly
end
import ParametricBodies
function update!(interface::GismoInterface,sim::Simulation;kwargs...)
    # update the geom as this has not been done yet
    for i in 1:interface.N
        new = interface.ControlPoints[i].+interface.deformation[i]
        interface.dir[i] != 1 && (new = reverse(new;dims=2))
        # time step is the (numerical) time between data exchange
        new = SMatrix{size(new)...}(new.*sim.L.+interface.center)
        sim.body.bodies[i] = ParametricBodies.update!(sim.body.bodies[i],new,sim.flow.Δt[end])
    end
    # solver update
    WaterLily.measure!(sim)
end

function writeData!(interface::GismoInterface)
    # write the force at the integration points
    PreCICE.writeData("ForceMesh", "ForceData", interface.quadPointID, permutedims(interface.forces))
end