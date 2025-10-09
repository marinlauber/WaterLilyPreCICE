"""
    An structure to hold the coupling data for WaterLily-Gismo coupling
"""
struct GismoInterface <: AbstractInterface
    ControlPointsID :: AbstractArray
    ControlPoints   :: AbstractArray
    quadPointID     :: AbstractArray
    quadPoint       :: AbstractArray
    forces          :: AbstractArray
    deformation     :: AbstractArray
    knots           :: AbstractArray
    dt              :: Vector{Float64}
    dir             :: Vector{Int16}
    N               :: Int16
    center          :: AbstractVector
end

import ParametricBodies: NurbsCurve, DynamicNurbsBody
function GismoInterface(T=Float64; dir=nothing, passive_bodies=nothing, center=0, kwargs...)

    # initilise PreCICE
    PreCICE.initialize()

    # get the mesh verticies from the fluid solver
    (_, knots) = getMeshVertexIDsAndCoordinates("KnotMesh")
    knots = knotVectorUnpack(knots)

    # get the mesh verticies from the structural solver
    (ControlPointsID, ControlPoints) = getMeshVertexIDsAndCoordinates("ControlPointMesh")
    ControlPointsID = Array{Int32}(ControlPointsID)
    ControlPoints = getControlPoints(ControlPoints, knots)
    deformation = copy(ControlPoints)
    isnothing(dir) ? (direction = ones(length(ControlPoints))) : direction = dir

    # get the quad points in parameter space
    (quadPointID, quadPoint) = getMeshVertexIDsAndCoordinates("ForceMesh")
    forces = zeros(reverse(size(quadPoint))...)
    quadPointID = Array{Int32}(quadPointID)
    quadPoint = quadPointUnpack(quadPoint)

    dt = PreCICE.getMaxTimeStepSize()

    # construct the interface curves
    body = NoBody() # empty first body
    for (i,(cps,knot)) in enumerate(zip(ControlPoints,knots))
        direction[i] != 1 && (cps = reverse(cps;dims=2))
        cps = SMatrix{2,size(cps,2)}(cps)
        knot = SVector{length(knot)}(knot)
        weights = SA[ones(size(cps,2))...]
        body += DynamicNurbsBody(NurbsCurve(cps.+center, knot, weights))
    end

    # coupling interface
    interface = GismoInterface(ControlPointsID, ControlPoints, quadPointID, quadPoint, forces,
                               deformation, knots, [dt], direction, length(bodies), center)

    # add some passive_bodies if we want
    for b in passive_bodies
       body += b
    end

    # return the interface and the body
    return interface, body
end

function readData!(interface::GismoInterface)
    # Read control point displacements
    readData = transpose(PreCICE.readData("ControlPointMesh", "ControlPointData", 
                                          interface.ControlPointsID, interface.dt[end]))
    interface.deformation .= getDeformation(readData, interface.knots) # repack correctly
end
import ParametricBodies
function update!(sim::CoupledSimulation, interface::GismoInterface; kwargs...)
    # update the geom as this has not been done yet
    for i in 1:interface.N
        new = interface.ControlPoints[i].+interface.deformation[i]
        interface.dir[i] != 1 && (new = reverse(new;dims=2))
        # time step is the (numerical) time between data exchange
        new = SMatrix{size(new)...}(new.*sim.L .+ interface.center)
        sim.body.bodies[i] = ParametricBodies.update!(sim.body.bodies[i], new, sim.flow.Δt[end])
    end
end

function writeData!(interface::GismoInterface)
    # write the force at the integration points
    PreCICE.writeData("ForceMesh", "ForceData", interface.quadPointID, permutedims(interface.forces))
end

using ParametricBodies: _pforce, _vforce
"""
    getInterfaceForces!(forces,flow::Flow,body::CombinedBodies,quadPoints)

Compute the interface forces at the quadrature points `quadPoints` for each body in `body.bodies`.
"""
Index(Qs,i) = sum(length.(Qs[1:i-1]))+1:sum(length.(Qs[1:i]))
# using ParametricBodies: _pforce, _vforce
function get_forces!(interface::GismoInterface, flow::Flow{T}, body, kwargs...) where T
    # for (i,b) in enumerate(body.bodies[1:length(interface.quadPoint)]) # only select the active curves
    #     I = Index(interface.quadPoint,i)
    #     fi = reduce(hcat,[-1.0*_pforce(b.curve,flow.p,s,zero(T),Val{false}()) for s ∈ interface.quadPoint[i]])
    #     interface.dir[i] != 1 && (fi = reverse(fi;dims=2))
    #     interface.forces[:,I] .= fi
    # end
end