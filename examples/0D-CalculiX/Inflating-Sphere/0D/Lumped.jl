using WaterLilyPreCICE,WriteVTK,StaticArrays
using Plots

# constant pressure, as in `Brown et al. 2024 CMAME 421, https://doi.org/10.1016/j.cma.2024.116764`
function constant_current(i,t,interface;p₀=0.0,Rₕ=0) # units are
    pᵢₙ = p₀ + WaterLilyPreCICE.get_Q(interface)*Rₕ # constant current
    return pᵢₙ*t
end

# vtk attributes
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces

# vtk attributes
custom = Dict("u"=>vtk_u, "f"=>vtk_f)

# zero function
zero_pressure(args...) = 0.0
# constant_pressure(i,t,args...) = 0.01*(t/1.0) # rap the pressure up

# coupling interface
# interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp", func=constant_current)


# keyword aguments might be specified
if size(ARGS, 1) < 1
    configFileName = "precice-config.xml"
else
    configFileName = ARGS[1]
end

# create coupling
PreCICE.createParticipant("LPM", configFileName, 0, 1)

# load the file
mesh0,srf_id = WaterLilyPreCICE.load_inp("../CalculiX/geom_q.inp")

# pass nodes and elements IDS to preCICE, here we use the unscaled and un-maped mesh
numberOfVertices, dimensions = length(mesh0.position), 3
vertices = Array{Float64,2}(undef, dimensions, numberOfVertices)
for i in 1:numberOfVertices
    vertices[:,i] .= mesh0.position[i].data # we need to have the same scale as in CalculiX
end
ControlPointsID = PreCICE.setMeshVertices("LPM-Nodes", reshape(vertices, (:,3)))
# mapping from center to nodes, needed for the forces
vertices .= 0 # reset since tey become the displacements
map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh0)))
forces = zeros(Float64, numberOfVertices, dimensions)
deformation = zeros(Float64, numberOfVertices, dimensions)

# make interface
interface = LumpedInterface(mesh0,deepcopy(mesh0),srf_id,deformation,ControlPointsID,
                            forces,constant_current,map_id,[0.0],[WaterLilyPreCICE.sum(volume(mesh0))/3.],[0.0],nothing,nothing,[],deepcopy(mesh0))

# check if we need to initialize the data
if PreCICE.requiresInitialData()
    println("Initial data required...")
    interface.forces .= 0 # reset the forces
    for (i,id) in interface.srf_id
        f = WaterLilyPreCICE.dS(@views(interface.mesh[id])).*interface.funct(i,sum(interface.Δt),interface)
        interface.forces[interface.map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
    end
    PreCICE.writeData("LPM-Nodes", "Forces", interface.ControlPointsID, interface.forces)
end

# initialise PreCICE
PreCICE.initialize()

# WaterLilyPreCICE.get_forces!(interface)
@show interface.V

# make the writer
wr = vtkWriter("Sphere"; attrib=custom)
v = [] # storage for the volume

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    # readData!(interface)
    dt_precice = PreCICE.getMaxTimeStepSize()
    #@TODO get max timestep from Lumped model
    push!(interface.Δt, min(1, dt_precice)) # min physical time step

    if PreCICE.requiresWritingCheckpoint()
        # save the mesh at this step
        interface.mesh_store = deepcopy(interface.mesh)
    end
    # Read control point displacements
    interface.deformation .= PreCICE.readData("LPM-Nodes", "Displacements",
                                              interface.ControlPointsID, interface.Δt[end])
    println("Displacement: ", sum(interface.deformation,dims=1))

    # compute the pressure forces
    WaterLilyPreCICE.update!(interface)
    println("Forces: ", sum(interface.forces,dims=1))

    # the volume is not correct with LumpedInterface
    pop!(interface.V); push!(interface.V, sum(WaterLilyPreCICE.volume(interface.mesh))/3)
    @show interface.V[end]
    push!(interface.P,interface.func(1,sum(interface.Δt),interface))
    @show interface.P[end]
    @show interface.Δt[end]
    push!(v, WaterLilyPreCICE.get_Q(interface))
    @show v[end]

    # write data to the other participant
    # writeData!(interface)
    PreCICE.writeData("LPM-Nodes", "Forces", interface.ControlPointsID, interface.forces)
    
    # do the coupling
    PreCICE.advance(interface.Δt[end]) # advance to t+Δt, check convergence and accelerate data
    
    # read checkpoint if required or move on
    if PreCICE.requiresReadingCheckpoint()
        # revert the mesh
        interface.mesh = deepcopy(interface.mesh_store)
        # pop the flux and pressures
        pop!(interface.Δt); pop!(interface.V); pop!(interface.P)
    end

    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        # save the data
        length(interface.Δt)%1==0 && write!(wr,interface)
    end
end
@show interface.P
@show interface.V
@show v
close(wr)
# plot(interface.V./1e-5,interface.P*1333.3,label=:none,xlabel="volume [ml]",ylabel="pressure [mmHg]"); savefig("pressure-volume.png")
PreCICE.finalize()