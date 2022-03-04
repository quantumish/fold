import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections as mc
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import fold

seq = fold.Sequence("MSVADDDLGSLQ")

fig, axs = plt.subplots(10,8, subplot_kw=dict(projection='3d'))
gs = axs[8, 0].get_gridspec()
for i in axs[8:, 0:]:
    for ax in i:
        ax.remove()
axbig = fig.add_subplot(gs[8:, 0:])
costs = []

def animate(i):
    this_costs = []
    for a in axs[:8]:
        for j in a:
            j.clear()
            j.axis("off")
            p = fold.RawProtein.random(seq)
            this_costs.append(fold.get_cost(p))
            residues = fold.get_residues(p)
            aminos, pos = zip(*residues)
            x, y, z = zip(*pos)
            chain = [[i, pos[x+1]] for x,i in enumerate(pos[:-1])]
            cc = Line3DCollection(chain, color="black", linewidths=2)
            j.add_collection(cc)
            j.scatter(x, y, z, s=10)
            j.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
    
    costs.append(sum(this_costs)/len(this_costs))
    axbig.clear()
    axbig.plot(costs)
        
    plt.figtext(0.5, 0.9, "Initializing", wrap=True, horizontalalignment='center', fontsize=25, family="monospace", fontweight="bold")
ani = animation.FuncAnimation(fig, animate, interval=1)

plt.show()
