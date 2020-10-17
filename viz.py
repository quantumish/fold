import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections as mc
import fold

fig, axs = plt.subplots(1, 2)


protein = fold.Protein("HPPHPH")
energies = []

def animate(i):
    axs[0].clear()
    axs[1].clear()
    axs[0].set_ylim([-5, 10])
    axs[0].set_xlim([-5, 10])
    protein.update()
    energies.append(protein.score)
    h_x = []
    h_y = []
    p_x = []
    p_y = []
    chain = []
    bonds = []
    for x,i in enumerate(protein.residues):
        if i.polar == True:
            p_x.append(i.coords[0])
            p_y.append(i.coords[1])
        else:
            h_x.append(i.coords[0])
            h_y.append(i.coords[1])
            for y,j in enumerate(protein.residues):
                if (abs(y-x) > 1 and ((abs(i.coords[0]-j.coords[0]) == 1 and abs(i.coords[1]-j.coords[1]) == 0) or (abs(i.coords[0]-j.coords[0]) == 0 and abs(i.coords[1]-j.coords[1]) == 1))):
                    bonds.append([i.coords, j.coords])
        if (x != len(protein.residues)-1):
            chain.append([i.coords, protein.residues[x+1].coords])
    cc = mc.LineCollection(chain, color='black', linewidths=2, label="Chain connection")
    bc = mc.LineCollection(bonds, color='red', linewidths=2, label="Bond")
    axs[0].add_collection(cc)
    axs[0].add_collection(bc)
    axs[0].scatter(h_x, h_y, zorder=2, label='Hydrophobic')
    axs[0].scatter(p_x, p_y, zorder=2, label='Polar')
    axs[0].legend()
    axs[1].plot(energies)

ani = animation.FuncAnimation(fig, animate, interval=1000, frames=300) 
plt.show()

