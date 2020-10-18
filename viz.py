import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections as mc
import fold

fig, axs = plt.subplots(2, 2)
sequence = "HPPHPH"
protein = fold.Protein(sequence, 2, False)
energies = []
densities = []
exposures = []

def animate(i):
    axs[0,0].clear()
    axs[0,1].clear()
    axs[1,0].clear()
    axs[1,1].clear()
    axs[0,0].set_ylim([-len(sequence)/3, len(sequence)+len(sequence)/3])
    axs[0,0].set_xlim([-5, 10])
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
                if (abs(y-x) > 1 and ((abs(i.coords[0]-j.coords[0]) == 1 and abs(i.coords[1]-j.coords[1]) == 0) or (abs(i.coords[0]-j.coords[0]) == 0 and abs(i.coords[1]-j.coords[1]) == 1)) and j.polar == False):
                    bonds.append([i.coords, j.coords])
        if (x != len(protein.residues)-1):
            chain.append([i.coords, protein.residues[x+1].coords])
    cc = mc.LineCollection(chain, color='black', linewidths=2, label="Chain connection")
    bc = mc.LineCollection(bonds, color='red', linewidths=2, label="Bond")
    area = ((max(h_x) - min(h_x))*(max(h_y) - min(h_y)))
    if (area > 0): density = len(protein.residues)/area
    else: density = 1
    densities.append(density)
    exposure = 0
    for x, i in enumerate(protein.residues):
        if i.polar == False:
            exposure += fold.exposure(protein.residues, x)
    exposure /= len(h_x)
    exposures.append(exposure)
    axs[0,0].add_collection(cc)
    axs[0,0].add_collection(bc)
    axs[0,0].scatter(h_x, h_y, zorder=2, label='Hydrophobic')
    axs[0,0].scatter(p_x, p_y, zorder=2, label='Polar')
    axs[0,0].legend()
    axs[1,0].plot(energies)
    axs[1,1].plot(densities, color="#6cad50")
    axs[0,1].plot(exposures, color="#f2444f")
    axs[0,0].set_title("Lattice protein")
    axs[1,0].set_title("Energy score of protein")
    axs[1,0].set_xlabel("Iterations")
    axs[1,0].set_ylabel("Energy score")
    axs[1,1].set_title("Density of hydrophobic residues")
    axs[1,1].set_xlabel("Iterations")
    axs[1,1].set_ylabel("Density")
    axs[0,1].set_title("Average exposure of hydrophobic residues")
    axs[0,1].set_xlabel("Iterations")
    axs[0,1].set_ylabel("Avg # of exposed sides")


ani = animation.FuncAnimation(fig, animate, interval=1) 
plt.show()

