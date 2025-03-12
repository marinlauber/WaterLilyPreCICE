using WaterLilyPreCICE,WriteVTK,StaticArrays
using Plots

# constant pressure, as in `Brown et al. 2024 CMAME 421, https://doi.org/10.1016/j.cma.2024.116764`
function constant_current(i,t,interface;p₀=0.0,Rₕ=0.0) # units are
    pᵢₙ = p₀ + interface.Q[end]*Rₕ # constant current
    return pᵢₙ
end

# vtk attributes
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces

# vtk attributes
custom = Dict("u"=>vtk_u, "f"=>vtk_f)

# coupling interface
interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp", func=constant_current)

# make the writer
wr = vtkWriter("Sphere"; attrib=custom)
v = [] # storage for the volume

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(interface)
    @show interface.deformation[1:3,1:3]

    # compute the pressure forces
    WaterLilyPreCICE.update!(interface)
    push!(interface.P, interface.func(1,0,interface))

    # write data to the other participant
    writeData!(interface)

    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        # save the data
        length(interface.dt)%1==0 && write!(wr,interface)
        push!(v, WaterLilyPreCICE.volume(interface.mesh)[1])
    end
end
@show interface.P
@show interface.Q
@show v
close(wr)
plot(v./1e-5,interface.P[2:end]*1333.3,label=:none,xlabel="volume [ml]",ylabel="pressure [mmHg]"); savefig("pressure-volume.png")
PreCICE.finalize()