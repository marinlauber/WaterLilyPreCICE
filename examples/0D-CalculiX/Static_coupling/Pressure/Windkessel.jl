using GeometryBasics,WriteVTK,StaticArrays,PreCICE
using Printf,CSV,DataFrames
using LinearAlgebra: cross
using OrdinaryDiffEq,Roots

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

# oriented area of a triangle
@inbounds @inline dS(tri::GeometryBasics.Ngon{3}) = 0.5f0SVector(cross(tri.points[2]-tri.points[1],tri.points[3]-tri.points[1]))

# update the mesh with the new displacements
function update_mesh(mesh0, deformation)
    # update the mesh so that any measure on it is correct
    points = Point3f[]
    for (i,pnt) in enumerate(mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ deformation[i,:]))
    end
    return GeometryBasics.Mesh(points,GeometryBasics.faces(mesh0))
end

# compute the volume of the mesh, assumes the sphere is centered at 0 and is uniformly deformed
@inline get_volume(mesh) =4/3*π*(√sum(abs2, mesh.position[1] .- 0))^3


# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.05,Emax=2,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end # amplitude is 2

# pressure volume loop function
@inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)

# open-loop windkessel
function Windkessel!(du,u,p,t)
    # unpack
    (VLV, Pao) = u
    (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C)  = p

    #first calculate PLV from elastance and VLV
    PLV = computePLV(t,VLV)

    # calculate Qmv, Pfilling>PLV; forward transmitral flow, PLV>Pfilling - backward transmitral flow
    Qmv = Pfill ≥ PLV ? (Pfill-PLV)/Rmv_fwd : (PLV-Pfill)/Rmv_bwd

    #calculate Qao, PLV>Pao; forward aortic flow, PLV>Pao; backward aortic flow
    Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd

    # rates
    du[1] = Qmv - Qao          #dVLV/dt=Qmv-Qao
    du[2] = Qao/C - Pao/(R*C)  #dPao/dt=Qao/C-Pao/RC
end

# solver the problem
let

    # keyword argument might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # load the mesh
    mesh0,srf_id = load_inp("../Solid/geom.inp")

    # create coupling
    PreCICE.createParticipant("LPM", configFileName, 0, 1)

    # initialise PreCICE
    PreCICE.initialize()

    # we need to initialize before we can get the mesh points and coordinates
    (ControlPointsID, ControlPoints) = PreCICE.getMeshVertexIDsAndCoordinates("Solid-Mesh")
    ControlPointsID = Array{Int32}(ControlPointsID)
    vertices = Array{Float64,2}(transpose(reshape(ControlPoints,reverse(size(ControlPoints)))))
    verts = GeometryBasics.Point3f[]
    for i in 1:size(vertices,1)
        push!(verts, GeometryBasics.Point3f(vertices[i,:]))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh0))

    # mapping from center to nodes, needed for the forces
    forces = zeros(Float64, size(ControlPoints))

    # map triangle to nodes, needed for the forces
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh)))

    # solver setting
    solver_dt = 1.0
    time = [0.0]

    # initialise storage
    mesh_storage = deepcopy(mesh)
    vol0 = 100*get_volume(mesh) # convert to ml

    # iteration storage
    storage_step,pressure_iter,volume_iter = [],[],[]
    p₁ = p₀ = 0.
    iteration = step = 1
    Q_target = 60.0  # ml/s between t=0 and t=1
    Cₕ = 0.00016     # relaxation factor for the pressure
    simple = true    # use a simple fixed-point or a secant method

    # set model parameters
    # Cardiac parameters
    Emax = 2                    #mmHg/ml; slope of the ESPVR
    V0 = 20                     #ml; intercept with volume axis of the ESPVR
    Pfilling = 5                #mmHg; venous filling pressure
    EDV = 120                   #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin
    Emin = Pfilling/(EDV-V0)    #mmHg/ml
    HR = 60                     #heart rate in beats/min

    # Valve resistances
    Rmv_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rmv_bwd = 1e10              #mmHg/ml/s; leak resistance
    Rao_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rao_bwd = 1e10              #mmHg/ml/s; leak resistance

    # Arterial model parameters
    R_WK2 = 1                   #mmHg/ml/s
    C_WK2 = 2                   #ml/mmHg

    #Setup
    u₀ = [EDV, 60] # initial conditions
    tspan = (0.0, 10.0)
    params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]

    # generate a problem to solve
    prob = ODEProblem(Windkessel!, u₀, tspan, params)

    # full control over iterations
    integrator = init(prob, Tsit5(), dtmax=0.001, reltol=1e-6, abstol=1e-9,
                      save_everystep=false)
    t_sol = []
    # initial storage
    ut_store = [integrator.t, integrator.u...]

    # main time loop
    while PreCICE.isCouplingOngoing()

        if PreCICE.requiresWritingCheckpoint()
            # save state at sum(time) = t
            mesh_storage = deepcopy(mesh)
            ut_store = [integrator.t, integrator.u...]
        end

        # set time step
        precice_dt = PreCICE.getMaxTimeStepSize()
        dt = min(precice_dt, solver_dt)
        push!(time, dt) # sum(time) = t+Δt (end of this time step)

        # read data from other participant
        displacements = PreCICE.readData("Solid-Mesh", "Displacements", ControlPointsID, dt)
        # update the mesh
        mesh = update_mesh(mesh0, displacements)
        vol = 100*get_volume(mesh) # scale to ml

        # solve the ODE
        OrdinaryDiffEq.step!(integrator, dt, true) # stop exactly at t+Δt
        push!(t_sol, [integrator.t, integrator.u...])

        # volume of sphere
        push!(volume_iter, vol)

        # target volume
        # target_volume = vol0 + Q_target*sum(time) # target flow rate is Q_target
        target_volume = integrator.u[1] # volume from 0D

        # first calculate PLV and the flow rates
        PLV = computePLV(sum(time),target_volume)
        Pao = integrator.u[2]
        Qmv = Pfilling ≥ PLV ? (Pfilling-PLV)/Rmv_fwd : (PLV-Pfilling)/Rmv_bwd
        Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd

        # fixed-point for the pressure
        if simple==true
            p₁ = p₀ + Cₕ*(target_volume - volume_iter[end])
        else
            ∂p = pressure_iter[end] - pressure_iter[end-1]
            ∂v = volume_iter[end] - volume_iter[end-1]
            p₁ = p₀ + (∂p/∂v)*(target_volume - volume[end])
        end
        p₀ = p₁
        
        # actuation pressure at this time
        p_act = 0.005 * Elastance(sum(time); Emin=Emin, Emax=Emax)
        p₁ = 0.0125
        push!(pressure_iter, p₁-p_act)

        # we then need to recompute the forces with the correct volume and pressure
        forces .= 0 # reset the forces
        for (i,id) in srf_id
            f = dS(@views(mesh[id])) .* (p₁-p_act) # pressure jump
            forces[map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
        end

        # write the force at the nodes
        PreCICE.writeData("Solid-Mesh", "Forces", ControlPointsID, forces)

        # do the coupling
        PreCICE.advance(dt) # advance to t+Δt

        # save the data
        push!(storage_step, [step, iteration, sum(time), p₁, vol, target_volume,
                             integrator.u..., Qao, Qmv, PLV])

        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint()
            # revert to sum(time) = t
            mesh = deepcopy(mesh_storage)
            pop!(time) # remove last time step since we are going back
            # revert integrator
            SciMLBase.set_ut!(integrator, ut_store[2:end], ut_store[1])
            iteration += 1
        else
            # we can move on, reset counter and increment step
            iteration = 1
            step += 1
        end

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            out_data = reduce(vcat,storage_step')
             CSV.write("sphere_output.csv", DataFrame(out_data,:auto),
                       header=["timestep","iter","time","pressure","volume",
                               "target_volume","VLV", "PAO","Qao", "Qmv", "Plv"])
        end
    end
    # save the nodal force at the end of the simulation
    # open("../Solid/cload.nam", "w") do io
    #     for i in 1:size(forces,1), j in 1:3
    #         println(io, @sprintf("%d, %d, %1.8f", i, j, forces[i,j]))
    #     end
    # end
    PreCICE.finalize()
end