import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections as mc
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import fold

seq = fold.Sequence("HPHPHP")

fig, axs = plt.subplots(8,8, subplot_kw=dict(projection='3d'))
for i in axs:
    for j in i:
        j.axis("off")
        p = fold.RawProtein.random(seq)
        residues = fold.get_residues(p)
        aminos, pos = zip(*residues)
        x, y, z = zip(*pos)
        chain = [[i, pos[x+1]] for x,i in enumerate(pos[:-1])]
        cc = Line3DCollection(chain, color="black", linewidths=2)
        j.add_collection(cc)
        j.scatter(x, y, z, s=10)        
        j.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
        
plt.show()
