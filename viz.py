import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections as mc
import fold

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
def threedimview(step, temp):
    sequence = "NLYIQWLKDGGPSSGRPPPS"
    protein = fold.Protein(sequence, temp, True)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    energies=[]
    for i in range(step):
        protein.update()
    backbones, polar_r, nonpolar_r = [], [], []
    chain, bonds = [], []
    for x,i in enumerate(protein.residues):
        print(i.backbone)
        if (x != len(protein.residues)-1):
            backbones.append(i.backbone)
            chain.append([i.backbone, protein.residues[x+1].backbone])
        cc = Line3DCollection(chain, color='black', linewidths=2, label="Chain connection")
        bc = Line3DCollection(bonds, color='red', linewidths=2, label="Bond")
        # ax.add_collection(cc)
        # ax.add_collection(bc)
        ax.scatter([x[0] for x in backbones], [x[1] for x in backbones], [x[2] for x in backbones], zorder=2, label='Backbones')
    plt.show()

#monitor(1000)
threedimview(100,1)
