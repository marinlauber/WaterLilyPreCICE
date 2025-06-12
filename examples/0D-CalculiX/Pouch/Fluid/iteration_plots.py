import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000

plt.style.use('peirlinck_lab')

cm = lambda x : x/2.56

conv = np.genfromtxt("precice-LPM-convergence.log", skip_header=1)
iter = np.genfromtxt("precice-LPM-iterations.log", skip_header=1)

x = np.array([])
resDispl = np.array([])
resForces = np.array([])

for it in range(1,iter.shape[0]+1):
    array = conv[conv[:,0]==it,2:]
    x = np.hstack((x,it+np.linspace(0,1,array.shape[0])))
    resDispl = np.hstack((resDispl,array[:,0]))
    resForces = np.hstack((resForces,array[:,1]))

fig1 = plt.figure("3D-0D iterations",figsize=(cm(18),cm(10)), constrained_layout=True)
gs = fig1.add_gridspec(nrows=1, ncols=2)
ax = fig1.add_subplot(gs[0,0])
ax.semilogy(x, resDispl, label="Displacements")
ax.semilogy(x, resForces, label="Forces")
ax.legend()
ax.set_ylabel("Relative Residuals")
ax.set_xlabel("Time-steps")
ax.set_xlim(0,iter.shape[0])
ax2 = ax.twinx()
ax2.hist(conv[:,0]-1,bins=iter.shape[0]+1,histtype='stepfilled',
         align="mid",color="gray",alpha=0.3)
ax2.set_ylabel("Iterations")
ax3 = fig1.add_subplot(gs[0,1])
ax3.hist(iter[:,2],bins=np.arange(1,max(iter[:,2]+1)),histtype='stepfilled',
        align="mid",color="gray",alpha=1.0)
ax3.set_yscale('log')
ax3.set_xlabel("Iterations")
ax3.set_ylabel("Frequency")
fig1.savefig("iterations.png", dpi=600)