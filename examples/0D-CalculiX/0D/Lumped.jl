using WaterLilyPreCICE,WriteVTK,StaticArrays,Interpolations
using Plots,JLD2

function Elastance(t;HR,Emin,Emax,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin   
end

# parameters
Emax = 2      # mmHg/ml; slope of the ESPVR
Emin = 0.05   # mmHg/ml
HR = 60       # heart rate in beats/min

# this are the loading curves inside CalculiX
A1 = linear_interpolation([0.,10.,100], [0.,1.,1.])
A2 = linear_interpolation([0.,10.,100.], [0.,1.,10.])
A3 = linear_interpolation([0.,10.,100.], [0.,1.,-9.])
A4(t) = ifelse(t<10, t/10,  1.0 + Elastance(t/10;HR,Emin,Emax))

function static_inflation(i,t)
    i==1 && return  0.38*A1(t)
    i==6 && return -0.38*A1(t)
    i in [2,3] && return 0.38*A3(t)
    i in [4,5] && return 0.38*A2(t)
    i in [7,8] && return  -0.38*A3(t)
    i in [9,10] && return -0.38*A2(t)
end

# vtk attributes
vtk_srf(a::LumpedInterface) = getindex.(a.srfID,1)
vtk_center(a::LumpedInterface) = WaterLilyPreCICE.center.(a.mesh)
vtk_normal(a::LumpedInterface) = WaterLilyPreCICE.normal.(a.mesh)
vtk_dS(a::LumpedInterface) = WaterLilyPreCICE.dS.(a.mesh)
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces

# vtk attributes
custom = Dict("SRF" =>vtk_srf, "center"=>vtk_center, "normal"=>vtk_normal,
"dS" => vtk_dS, "u" => vtk_u, "f"=>vtk_f)

# coupling interface
interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp", func=static_inflation)

# make the writer
wr = vtkWriter("LIMO_4"; attrib=custom)
v = [] # storage for the volume

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(interface)

    # compute the pressure forces
    WaterLilyPreCICE.update!(interface)

    # write data to the other participant
    writeData!(interface)

    # do the coupling
    length(interface.dt)%5==0 && write!(wr,interface)
end
close(wr)
jldsave("LIMO_4_volume.jld2";volume=hcat(v...))
PreCICE.finalize()