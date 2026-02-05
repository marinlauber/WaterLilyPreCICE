#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import argparse
import sys
import os

sys.path.append("/home/marin/Workspace/HHH/code/")
from calculix import DataSet,load_srf,center,normal,area,volume,S4_to_tri


def run_analysis(folder, case, every=1, Np=4, ext="", Press=None):

    # get what pressure to use
    print(folder+"/"+case)

    # main script
    CalculiXResults = DataSet(folder+"/",case)
    # load undeformed mesh and add the displacement
    k0 = CalculiXResults.get_snapshot(10)
    data,u0,s_,e_,_,_,_,connectivity = CalculiXResults.load_snapshot(k0)

    # read the geom to get the surface IDs
    SRF = load_srf(folder+"/geom.inp")
    # fix connectivity for S8/M8 elements
    connectivity,SRF = S4_to_tri(connectivity,SRF,Np)

    # Np = (len(SRF)-2)/4 # the number of pouches is then number of (surfaces - 2)/4
    print("Number of pouches : ", Np)

    # get edge nodes to close srf
    edge = np.unique(np.where(data[:,1]==0.0))
    theta = np.arctan2(data[edge,0]+u0[0,edge],data[edge,2]+u0[2,edge])
    edge = edge[np.argsort(theta)]

    # make last node to close geom
    center_p = np.array([0.0,0.0,0.0])
    cap_connectivity = np.empty(3)

    # generate new triangles
    last_id = len(data[:,0])
    for i in range(len(edge)):
        cap_connectivity = np.vstack((cap_connectivity,np.array([last_id,edge[(i+1)%len(edge)],edge[i]])))

    # add cap triangles
    new_connectivity = np.vstack((connectivity,cap_connectivity))

    # prepare data
    data_cap = np.vstack((data, center_p))
    u_ = np.vstack((u0.T, np.zeros(3))).T

    # check the area of the cap
    defomed = data_cap + u_.T
    trian  = tri.Triangulation(defomed[:,0], defomed[:,1], triangles=new_connectivity[-len(edge):-1,:])
    print("Cap area ", area(defomed, trian))
    print("Cap area (πr²/AR) ", np.pi*np.max(abs(defomed[edge,0]))*np.max(abs(defomed[edge,2])))

    # extract surfaces
    surfaces = np.hstack((np.arange(0,Np+1),np.arange(1+2*Np,2+3*Np)))

    # for each of the surface, construct a triangulation
    tri_list = []

    for i in range(len(SRF)):
    # for i in [0]:
        pouch = np.array([el-SRF[0][0] for el in SRF[int(i)]])
        pouch = tri.Triangulation(data_cap[:,0], data_cap[:,1], triangles=new_connectivity[pouch])
        tri_list.append(pouch)
    # cap_trian = tri.Triangulation(data_cap[:,0], data_cap[:,1], triangles=cap_connectivity)
    # test for contact between the inside of pouch
    slave = np.array([el-SRF[0][0] for i in np.arange(1,Np+1) for el in SRF[int(i)]])
    slave  = tri.Triangulation(data[:,0], data[:,1], triangles=connectivity[slave])
    master = np.array([el-SRF[0][0] for i in np.arange(2*Np+2,3*Np+2) for el in SRF[int(i)]])
    master  = tri.Triangulation(data[:,0], data[:,1], triangles=connectivity[master])

    # main script
    CalculiXResults = DataSet(folder+"/",case)
    pouch_volume = []; ventricle_volume = []; area_ventricle = []; time = []; contact = []
    print("Iterating though snapshots...")

    # plot to check
    # ax = plt.figure().add_subplot(projection='3d')
    # for srf in tri_list:
    #     ax.plot_trisurf(data[:,0]+u0[0,:],data[:,1]+u0[1,:],data[:,2]+u0[2,:], linewidth=0.2,
    #                     triangles=srf.triangles,zsort='min')
    # ax.set_aspect("equal","box")
    # plt.show()

    first = True
    for k in range(k0,len(CalculiXResults.time),every):
        # load undeformed mesh and add the displacement
        data,u_,s_,e_,_,_,_,_ = CalculiXResults.load_snapshot(k)
        print("Time step ", CalculiXResults.time[k])
        # prepare data, add the cap node and the displacement
        data = np.vstack((data,center_p))
        u_ = np.vstack((u_.T,np.zeros(3))).T
        data[:,0] += u_[0,:]
        data[:,1] += u_[1,:]
        data[:,2] += u_[2,:]

        # loop over all different contributions, all the normal are no the same side but they cannot be for the volume, must be switched
        vp = np.zeros(3); vv = np.zeros(3); arrea = 0.0; a0=0

        # first set of pouches
        for srf_id,sgn in zip(np.arange(1,2*Np+1),np.where(np.arange(1,2*Np+1)<=Np,1,-1)):
            pouch = tri_list[int(srf_id)]
            vp += sgn*volume(data,pouch)
            a0 += sgn*area(data,pouch)
        # second set of pouches
        for srf_id,sgn in zip(np.arange(2*Np+2,4*Np+2),np.where(np.arange(2*Np+2,4*Np+2)<=3*Np+1,-1,1)):
            pouch = tri_list[int(srf_id)]
            vp += sgn*volume(data,pouch)
            a0 += sgn*area(data,pouch)
        # total volume of the ventricle
        for srf_id,sgn in zip(surfaces,np.where(surfaces<=Np+1,1,-1)):
            ventricle = tri_list[int(srf_id)]
            vv += sgn*volume(data, ventricle)
            arrea += area(data, ventricle)

        # add cap contribution
        # print("Cap vol contribution ", volume(data, cap_trian))
        # vv += volume(data, cap_trian)

        # test for contact between the inside of pouch
        c_slave = center(data, slave)
        c_master = center(data, master)
        cts = np.min(np.linalg.norm(c_slave-c_master,axis=1)) < 1e-4 # we can do pairwise because the surfaces are identical
        if ((cts == True) and (first==True)):
            print("Contact detected!")
            # # plot to check
            # ax = plt.figure().add_subplot(projection='3d')
            # ax.plot_trisurf(data[:,0],data[:,1],data[:,2], linewidth=0.2,
            #                     triangles=slave.triangles,zsort='min')
            # ax.plot_trisurf(data[:,0],data[:,1],data[:,2], linewidth=0.2,
            #                     triangles=master.triangles,zsort='min')
            # ax.set_aspect("equal","box")
            # plt.show()
            first = False

        print(" contact %s" % cts)
        contact.append(cts)

        # check
        print(" pouch volumes     : ", vp)
        print(" ventricle volumes : ", vv)
        pouch_volume.append(abs(np.mean(vp)))
        ventricle_volume.append(abs(np.mean(vv)))
        area_ventricle.append(arrea)
        time.append(CalculiXResults.time[k]-CalculiXResults.time[k0])
        print("")
    print("Saving...")
    # volume diff from the start state
    ventricle_volume = np.array(ventricle_volume)
    pouch_volume = np.array(pouch_volume)
    area_ventricle = np.array(area_ventricle)
    contact = np.array(contact)
    time = np.array(time)
    np.save(folder+"/"+case.split(".")[0]+".npy", np.array([ventricle_volume,pouch_volume,
                                                                              area_ventricle,contact,time],dtype=object),
            allow_pickle=True)

#
#
# main script

#
run_analysis(".","initial_condition.pvd", every=1, Np=3)
run_analysis(".","espvr.pvd", every=1, Np=3)
