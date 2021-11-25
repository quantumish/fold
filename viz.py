import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections as mc
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import fold

lookup = "CMFILVWYAGTSNQDEHRKP"
polar = "QNHSTYC"


def threedimview(protein: fold.Protein, step: int, temp: float):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    for _i in range(step):
        protein.update()
    polar_r, nonpolar_r = [], []
    chain = []
    for x, i in enumerate(protein.residues):
        if x != len(protein.residues) - 1:
            if lookup[i.id] in polar:
                polar_r.append(i.sidechain)
            else:
                nonpolar_r.append(i.sidechain)
            chain.append([i.sidechain, protein.residues[x + 1].sidechain])
        # cc = Line3DCollection(chain, color="black", linewidths=2, label="Chain connection")
        # bc = Line3DCollection(bonds, color="red", linewidths=2, label="Bond")
        # ax.add_collection(cc)
        # ax.add_collection(bc)
        # ax.scatter(
        #     [x[0] for x in backbones],
        #     [x[1] for x in backbones],
        #     [x[2] for x in backbones],
        #     zorder=2,
        #     label="Backbones",
        # )
        ax.scatter(
            [x[0] for x in polar_r],
            [x[1] for x in polar_r],
            [x[2] for x in polar_r],
            zorder=2,
            label="Polar",
        )
        ax.scatter(
            [x[0] for x in nonpolar_r],
            [x[1] for x in nonpolar_r],
            [x[2] for x in nonpolar_r],
            zorder=2,
            label="Nonpolar",
        )
    plt.show()


sequence = "NLYIQWLKDGGPSSGRPPPS"
temp = 1
protein = fold.Protein(sequence, temp, True)
for i in range(10):
    protein.update()
    threedimview(protein, 1, temp)
