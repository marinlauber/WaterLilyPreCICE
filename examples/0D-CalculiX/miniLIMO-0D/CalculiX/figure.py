#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

def add_clinician_axis(ax,scale=7.50062,label=r"$P_a$ (mmHg)"):
    # units for clinicians
    yy = ax.twiny()
    yy.set_xlabel(label)
    yy.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
    yy.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
    yy.spines['bottom'].set_position(('outward', 36))
    yy.set_xlim(ax.get_xlim()[0]*scale, ax.get_xlim()[1]*scale)
    for s in ['right', 'left', 'top']:
        yy.spines[s].set_visible(False)
    return None

# read data
fig,ax = plt.subplots(1, 2, figsize=(8,6))
scale = 10.9 # mm per unit length

# Venous afterload
data = np.load("initial_condition.npy",allow_pickle=True)
ventricle_volume,pouch_volume,area_ventricle,contact,time = data
Pa = scale*0.0738*time*40/90
print("EDV (ml) = ", ventricle_volume[0]*scale**3)
ax[0].plot(Pa, ventricle_volume*scale**3, label="Vv (Pv=6.01mmHg)")
ax[0].plot(Pa, pouch_volume*scale**3, label="Va (Pv=6.01mmHg)")
ax[0].plot(Pa, (ventricle_volume[0]-ventricle_volume)*scale**3, label="ΔVv (Pv=6.01mmHg)")
ax[1].plot(pouch_volume*scale**3, (ventricle_volume[0]-ventricle_volume)*scale**3, label="(Pv=6.01mmHg)")

# Arterial afterload
data = np.load("espvr.npy",allow_pickle=True)
ventricle_volume,pouch_volume,area_ventricle,contact,time = data
Pa = scale*0.9785*time*40/90
ax[0].plot(Pa, ventricle_volume*scale**3, label="Vv (Pv=80mmHg)")
ax[0].plot(Pa, pouch_volume*scale**3, label="Va (Pv=80mmHg)")
ax[0].plot(Pa, (ventricle_volume[0]-ventricle_volume)*scale**3, label="ΔVv (Pv=80mmHg)")
ax[1].plot(pouch_volume*scale**3, (ventricle_volume[0]-ventricle_volume)*scale**3, label="(Pv=80mmHg)")
ax[1].plot([0,300],[0,300], 'k--', label=None)
ax[0].set_xlabel("Pa (kPa)")
ax[0].set_ylabel("Volume (ml)")
ax[0].legend()
ax[0].set_xlim(0,max(Pa)+5)
add_clinician_axis(ax[0])
ax[1].set_xlim(0,150)
ax[1].set_ylim(0,150)
ax[1].set_xlabel("Vₐ (ml)")
ax[1].set_ylabel("ΔVᵥ (ml)")
plt.tight_layout()
plt.savefig("pressure_volume_test.png", dpi=300)
# plt.show()
plt.close()