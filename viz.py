import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections as mc
import fold

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

lookup = "CMFILVWYAGTSNQDEHRKP"
polar = "QNHSTYC"

def threedimview(protein: fold.Protein, step: int, temp: float):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(step):
        protein.update()
    backbones, polar_r, nonpolar_r = [], [], []
    chain, bonds = [], []
    for x,i in enumerate(protein.residues):
        print(i.backbone)
        if (x != len(protein.residues)-1):
            backbones.append(i.backbone)
            if (lookup[i.id] in polar):
                polar_r.append(i.sidechain)
            else:
                nonpolar_r.append(i.sidechain)
            chain.append([i.backbone, protein.residues[x+1].backbone])
        cc = Line3DCollection(chain, color='black', linewidths=2, label="Chain connection")
        bc = Line3DCollection(bonds, color='red', linewidths=2, label="Bond")
        # ax.add_collection(cc)
        # ax.add_collection(bc)
        ax.scatter([x[0] for x in backbones], [x[1] for x in backbones], [x[2] for x in backbones], zorder=2, label='Backbones')
        ax.scatter([x[0] for x in polar_r], [x[1] for x in polar_r], [x[2] for x in polar_r], zorder=2, label='Polar')
        ax.scatter([x[0] for x in nonpolar_r], [x[1] for x in nonpolar_r], [x[2] for x in nonpolar_r], zorder=2, label='Nonpolar')
    plt.show()

sequence = "NLYIQWLKDGGPSSGRPPPS"
temp = 1
protein = fold.Protein(sequence, temp, True)
protein.update()
threedimview(protein, 1, temp)
