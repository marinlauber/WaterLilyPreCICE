using WaterLilyPreCICE,WriteVTK,StaticArrays,Interpolations
using Plots

# constant pressure, as in `Brown et al. 2024 CMAME 421, https://doi.org/10.1016/j.cma.2024.116764`
function constant_current(i,t,interface;p₀=100,Rₕ=0.0)
    pᵢₙ = p₀ + interface.Q[end]*Rₕ # constant current
end

# vtk attributes 
vtk_srf(a::LumpedInterface) = Float32[el[1] for el in a.srf_el]
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces

# vtk attributes
custom = Dict("u" => vtk_u, "f"=>vtk_f)

# coupling interface
interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp", func=constant_current)
WaterLilyPreCICE.get_forces!(interface)

# make the writer
wr = vtkWriter("Sphere"; attrib=custom)
v = [] # storage for the volume

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(interface)

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
        @show length(interface.P)
    end
end
close(wr)
plot(interface.P[2:end],v); savefig("pressure-volume.png")
PreCICE.finalize()