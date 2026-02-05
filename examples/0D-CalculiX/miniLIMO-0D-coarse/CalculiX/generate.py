#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import os
import sys
import matplotlib.pyplot as plt

sys.path.append("/home/marin/Workspace/HHH/code/")
from calculix import Calculix,Mesh,replace,map_initial

# get the length scale if provided
parser = argparse.ArgumentParser(description='Generates the Finite-element mesh..')
parser.add_argument('-L','--lengthscale', help='Length scale of the body')
parser.add_argument('-N','--npoints', help='number of points along the length')
parser.add_argument('-O','--order', help='Order of the finite-element mesh')
args = parser.parse_args()

# set default variables
L = 10. if(not args.lengthscale) else float(args.lengthscale)
order = 1 if(not args.order) else int(args.order)
N = 4 if(not args.npoints) else int(args.npoints)

# for 2nd order quads, we mesh using Incomplete elements (no centroid node)
string = "-string \"Mesh.SecondOrderIncomplete = 1;\" " if (order==2) else ""
# generate mesh file
os.system("gmsh template.geo " + string + "-order %d -2 -o init.inp" % order)

# mesh generation
mesh = Mesh("init.inp")
# duplicate inside the pouch
mesh.duplicate_srf([1],[6,14,15,16,17,18])
mesh.duplicate_srf([2],[7,13,18,19])
mesh.duplicate_srf([3],[8,9,10,11,12,19])
mesh.mirror_srf([0,1,2,3,4,5,6],[0,1,2,3,4],"z",[0,0,0])
mesh.make_nodeset([0,7],lambda x: x[1]==0.0,"EDGE")
mesh.make_nodeset([0],lambda x: x[1]==0.0,"EDGE_TOP")
mesh.make_nodeset([7],lambda x: x[1]==0.0,"EDGE_BOTTOM")
mesh.save("geom.inp", eltype="S")

# boundary conditions move edge in the x-z plane
for sets in ["EDGE_TOP","EDGE_BOTTOM"]:
    xs = np.zeros(len(mesh.node_set[sets]))
    zs = np.zeros(len(mesh.node_set[sets]))
    nodes = np.zeros(len(mesh.node_set[sets]))
    for i,node in enumerate(mesh.node_set[sets]):
        nodes[i] = node
        node = mesh.get_node_from_id(node)
        xs[i] = node[1]
        zs[i] = node[2]
    # arrange the points from +x -> -x
    idx = np.argsort(xs)
    xsort = xs[idx]; zsort = zs[idx]
    nodesort = nodes[idx]
    # aspect ratio of the ellipse
    Ar = 80.2/55.2 # from Nienke
    # map to ellipse
    dxs,dzs = map_initial(xsort,zsort,Ar,1.4)
    # plt.plot(xsort,zsort,'o',label='initial')
    # plt.plot(dxs,dzs,'x',label='mapped')
    # plt.gca().set_aspect('equal')
    # plt.show()
    # put back in correct order
    dxs = dxs - xs[idx]
    dzs = dzs - zs[idx]
    if sets == "EDGE_BOTTOM":
        dzs *= -1
    # write data to BCs file
    BCs = open("BCs_%s.nam" % sets,"w")
    for node,dx,dz in zip(nodesort,dxs,dzs):
        BCs.write(" %d, 1, 1, %f\n" % (node,dx))
        BCs.write(" %d, 3, 3, %f\n" % (node,dz))
    BCs.close()
