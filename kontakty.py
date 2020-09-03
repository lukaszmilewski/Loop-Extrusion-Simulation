import mdtraj as md
import scipy as scp
import numpy as np
import matplotlib.pyplot as plt
#tu należy podać lokalizacje pliku dcd a w parametrze top lokalizacje
#struktury poczatkowej. Na poczatku ścierzki ZAWSZE trzeba dać "r"
#bo inaczej string jest skopany i nie działa
t = md.load(r'C:\Users\lukas\PycharmProjects\Loop-Extrusion-Simulation\1sided-trj.dcd',top=r'C:\Users\lukas\PycharmProjects\Loop-Extrusion-Simulation\initial_structure.pdb')
M=t.xyz
sum_arr=np.zeros((100,100))

for i in range(M.shape[0]):
    X=scp.spatial.distance.squareform(scp.spatial.distance.pdist(M[i],'euclidean'))

    X[X < 0.3] = 0.1
    X[X >= 0.3] = 0
    X[X >0] = 1
    sum_arr+=X

plt.imshow(sum_arr, cmap='Reds', vmin=0,vmax=2400)
plt.grid(b=None)
plt.show()
plt.savefig("1sided-heatmap.png")