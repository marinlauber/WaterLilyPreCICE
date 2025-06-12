using Plots
using StaticArrays

run(`gmsh template.geo -3 -o template.msh`)
file = open("template.geo", "r")