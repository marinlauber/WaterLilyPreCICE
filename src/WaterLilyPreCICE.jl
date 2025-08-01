module WaterLilyPreCICE

using WaterLily
using StaticArrays
using PreCICE

using Reexport
@reexport using WaterLily
@reexport using PreCICE

abstract type AbstractInterface end
export AbstractInterface

# include helpers
include("utils.jl")
export Store,store!,revert!

# include mesh bodies
include("MeshBodies.jl")
export MeshBody
@reexport using GeometryBasics

# include the KDTree
include("Tree.jl")
export BVH_simple,BoundingBox

# structure for coupled simulation
mutable struct CoupledSimulation <: AbstractSimulation
    sim :: Simulation
    int :: AbstractInterface
    store :: Store
end
# overload properties
Base.getproperty(f::CoupledSimulation, s::Symbol) = s in propertynames(f) ? getfield(f, s) : getfield(f.sim, s)

"""
    CoupledSimulation((WaterLily.Simulation inputs)...;
                      interface=:Interface,
                      surface_mesh="geom.inp", scale=1.f0,
                      boundary=true, thickness=0f0, center=0.0,
                      curve_dir=nothing, passive_bodies=nothing,
                      func=(i,t)->0, prob=nothing)

Constructor for a WaterLily.Simulation that uses PreCICE for coupling with CalculiX or G+Smo:
    - `interface`     : Interface to use, either `:Interface, :CalculiXInterface`, `:GismoInterface`, or `:LumpedInterface`.
                        default is :Interface.
    - `surface_mesh`  : Path to the surface mesh file (not used for :GismoInterface).
    - `scale`         : Scaling factor for the mesh (not used for :GismoInterface).
    - `boundary`      : Is the mesh provided the boundary of the solid (default is true)
    - `thickness`     : Thickness of the solid (default is 0, boundary must be set
                        to false for this to take effects)
    - `center`        : Center of the solid (default is 0, can be used instead of `map(x,t)`
                        to move the structure in the domain)
    - `curve_dir`     : Direction of the curve for the Gismo interface (default is nothing)
    - `passive_bodies`: Passive bodies to add to the interface, they are immersed,
                        but no FSI occurs for those
    - `func`          : Function to apply to the interface, default is no function
                        (only used for :LumpedInterface)
    - `prob`          : Problem to solve, default is nothing (only used for :LumpedInterface)
"""
function CoupledSimulation(args...; T=Float64, mem=Array, interface=:Interface, # WL specific
                           surface_mesh="geom.inp", scale=1.f0, # CalculiX specific
                           boundary=true, thickness=0f0, center=0.0, passive_bodies=nothing, # general
                           curve_dir=nothing,  # Gismo specific
                           func=(i,t)->0, prob=nothing, # Lumped specific
                           kwargs...) # args and kwargs are passed to Simulation

    # check that the interface is constructed correctly
    @assert !(:map in keys(kwargs)) "The `map` keyword argument is not allowed in the `CoupledSimulation` constructor, instead we use `center::SVector`"
    @assert interface in [:Interface,:CalculiXInterface,
                          :GismoInterface,:LumpedInterface] "The interface specified does not exist"
    if interface==:GismoInterface
        @assert length(args[1])==2 "3D simulations are not support for Gismo coupling"
    elseif interface in [:Interface,:CalculiXInterface]
        length(args[1])==2 && @warn("\nThe Interface/CalculiXInterface assumes that 2D simulation happen in the x-y plane and that the z coordinate\n"*
                                    "of the structural mesh is zero. If this is not the case, the interface will not work as expected.")
    end

     # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant("WaterLily", configFileName, 0, 1)

    # generate the coupling interface from the coupling partner
    int, body = eval(interface)(T; surface_mesh, map ,scale, curve_dir, passive_bodies,
                                center, func, prob, boundary, thickness)

    # the simulation
    sim = Simulation(args...; mem, T, body, kwargs...)

    # storage for iteration
    store =  Store(sim)

    # return coupled sim
    CoupledSimulation(sim, int, store)
end

# overload sim_step!
import WaterLily: sim_step!
function sim_step!(sim::CoupledSimulation; kwargs...)
    # update the this participant this is type-specialized
    update!(sim.int, sim; kwargs...)
    # measure and momentum step
    measure!(sim); mom_step!(sim.flow, sim.pois)
    # get forces, this is type-specialized
    get_forces!(sim.int, sim.flow, sim.body; kwargs...)
end

macro notimplemented(message="not implemented")
    quote
        error($esc(message))
    end
end
function update!(::AbstractInterface, args...; kwargs...)
    @notimplemented "`update!` not implemented for `AbstractInterface`, it must to be (interface) specialized"
end
function get_forces!(::AbstractInterface, args...; kwargs...)
    @notimplemented "`get_forces!` not implemented for `AbstractInterface`, it must to be (interface) specialized"
end
function InitializeData!(::AbstractInterface, args...; kwargs...)
    @notimplemented "`InitializeData!` not implemented for `AbstractInterface`, it must to be (interface) specialized"
end

readData!(sim::CoupledSimulation) = readData!(sim.int, sim.sim, sim.store)
function readData!(interface::AbstractInterface,sim::Simulation,store::Store)
    # set time step
    dt_precice = PreCICE.getMaxTimeStepSize()
    push!(interface.dt, min(sim.flow.Δt[end]*sim.U/sim.L, dt_precice)) # min physical time step
    sim.flow.Δt[end] = interface.dt[end]*sim.L/sim.U # numerical time step

    if PreCICE.requiresWritingCheckpoint()
        store!(store,sim)
    end

    # Read control point displacements, type specific
    readData!(interface)
end

writeData!(sim::CoupledSimulation) = writeData!(sim.int, sim.sim, sim.store)
function writeData!(interface::AbstractInterface,sim::Simulation,store::Store)
    # write the force at the integration points
    writeData!(interface)

    # do the coupling @TODO this return zeros, so don;t replace the last dt
    PreCICE.advance(interface.dt[end])

    # read checkpoint if required or move on
    if PreCICE.requiresReadingCheckpoint()
        revert!(store,sim)
    end
end
export CoupledSimulation,readData!,writeData!,InitializeData!

# the general interface type
include("Interface.jl")
export Interface

# CalculiX specific interface functions
include("CalculiXInterface.jl")
export CalculiXInterface

# G+Smo specific interface functions
include("GismoInterface.jl")
export GismoInterface

include("LumpedInterface.jl")
export LumpedInterface

end # module