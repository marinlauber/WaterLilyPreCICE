using WaterLilyPreCICE,WriteVTK,StaticArrays
using Plots

# constant pressure, as in `Brown et al. 2024 CMAME 421, https://doi.org/10.1016/j.cma.2024.116764`
function constant_current(i,t,interface;p₀=0.01,Rₕ=0.05) # units are
    pᵢₙ = p₀ - WaterLilyPreCICE.get_Q(interface)*Rₕ # constant current
    return pᵢₙ
end

# vtk attributes
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces

# vtk attributes
custom = Dict("u"=>vtk_u, "f"=>vtk_f)

# coupling interface
interface = LumpedInterface(surface_mesh="../CalculiX/geom_q.inp", func=constant_current)
r0 = √sum(abs2,interface.mesh.position[1].-0)
pop!(interface.V); push!(interface.V,4/3*π*r0^3)

# make the writer
wr = vtkWriter("Sphere"; attrib=custom)
Qs = [] # storage for the volume
# write initial data
write!(wr,interface)

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(interface)
    
    # compute the pressure forces
    WaterLilyPreCICE.update!(interface)
   
    # the volume is not correct with this 1/8th of the sphere
    r = √sum(abs2,interface.mesh.position[1].-0) # arbitrary node to get radius
    pop!(interface.V); push!(interface.V,4/3*π*r^3)
    Qi = WaterLilyPreCICE.get_Q(interface)
    push!(interface.P, interface.func(1,0,interface))
    println("r: ", r)
    println("V: ", interface.V[max(1,length(interface.V)-2):end])
    println("Q: ", Qi)
    println("P: ", interface.P[max(1,length(interface.P)-2):end])

    # we then need to recompute the forces with the correct flow rate
    WaterLilyPreCICE.get_forces!(interface)
    
    # write data to the other participant
    writeData!(interface)
    
    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        push!(Qs, WaterLilyPreCICE.get_Q(interface))
        # save the data
        length(interface.Δt)%1==0 && write!(wr,interface)
        println("PreCICE.isTimeWindowComplete()")
        println("r: ", r)
        println("V: ", interface.V[max(1,length(interface.V)-2):end])
        println("Q: ", Qs[end])
        println("P: ", interface.P[max(1,length(interface.P)-2):end])
    end
end
@show interface.P
@show interface.V
@show Qs
plot(interface.V,interface.P,label=:none,xlabel="volume [ml]",ylabel="pressure [mmHg]");
savefig("pressure-volume.png")
close(wr)
PreCICE.finalize()