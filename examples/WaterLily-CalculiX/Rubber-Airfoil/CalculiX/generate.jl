using StaticArrays#,Plots

# rotate around the z-axis
rotate_z(x,θ) = SA[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]*(x-SA[0,0,0])
# airfoil curve
NACA(s) = 0.18f0*5*(0.2969f0s-0.126f0s^2-0.3516f0s^4+0.2843f0s^6-0.1036f0s^8)
curve(s) = SA[(1-2s)^2,sign(1-2s)*NACA(abs(1-2s)),0]

# make the points
N = 32
points = [curve(s) for s ∈ 0:1/N:1]
# plot(getindex.(points,1),getindex.(points,2))

# make the template file
f = open("template.geo", "w")
println(f, "SetFactory(\"OpenCASCADE\");")
for (i,pts) in enumerate(points)
    x,y,z = pts
    i ∈ [1,N+1] && (y = sign(y)*0.005; println("Point($i) = {$x, $y, $z};"))
    x,y,z = rotate_z(SA[x,y,z],-π/10) # 18 degrees rotation
    println(f,"Point($i) = {$x, $y, $z};")
end
print(f,"Spline(1) = {1")
for i in 2:length(points)
    print(f,", $i")
end
println(f,"};\nLine(2) = {$(N+1), 1};\nCurve Loop(1) = {1, 2};")
println(f,"Plane Surface(1) = {1};\n"
         *"Extrude {0, 0, 0.5} {\n"
         *"    Surface{1};\n"
         *"}\n"
         *"Transfinite Curve {1, 5} = 32 Using Progression 1;\n"
         *"Transfinite Curve {2} = 2 Using Progression 1;\n"
         *"Transfinite Curve {3, 4} = 16 Using Progression 1;\n"
         *"// Physical Surface(\"Interface\", 7) = {2, 3}; // this line is needed when making the surface mesh"
)
close(f)
# make the mesh, second order, 2D, then 3D
run(`gmsh template.geo -order 1 -3 -o init.inp`)
println("Done!")