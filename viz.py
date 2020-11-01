import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import collections as mc
import fold

def draw_protein(protein, ax):
    ax.clear()
    ax.set_ylim([-len(protein.residues)/3, len(protein.residues)+len(protein.residues)/3])
    ax.set_xlim([-len(protein.residues)/3, (len(protein.residues)+len(protein.residues)/3)/2])
    h_x, h_y, p_x, p_y = [], [], [], []
    chain, bonds = [], []
    for x,i in enumerate(protein.residues):
            if i.polar == True:
                p_x.append(i.coords[0])
                p_y.append(i.coords[1])
            else:
                h_x.append(i.coords[0])
                h_y.append(i.coords[1])
            for y,j in enumerate(protein.residues):
                if (abs(y-x) > 1 and ((abs(i.coords[0]-j.coords[0]) == 1 and abs(i.coords[1]-j.coords[1]) == 0) or (abs(i.coords[0]-j.coords[0]) == 0 and abs(i.coords[1]-j.coords[1]) == 1)) and j.polar == False and i.polar == False):
                    bonds.append([i.coords, j.coords])
            if (x != len(protein.residues)-1):
                chain.append([i.coords, protein.residues[x+1].coords])
    cc = mc.LineCollection(chain, color='black', linewidths=2, label="Chain connection")
    bc = mc.LineCollection(bonds, color='red', linewidths=2, label="Bond")
    ax.add_collection(cc)
    ax.add_collection(bc)
    ax.scatter(h_x, h_y, zorder=2, label='Hydrophobic')
    ax.scatter(p_x, p_y, zorder=2, label='Polar')
    return h_x, h_y, p_x, p_y

def standard_view(updates):
    fig, axs = plt.subplots(2, 2)
    sequence = "HHPPHHPHHPPHHPHHPPHH"
    sequence = "PPHPHPPPPHPHHPHP"
    protein = fold.Protein(sequence, 2, False)
    energies, densities, exposures = [], [], []
    def animate(i):
        for i in range(updates):
            protein.update()
            energies.append(protein.score)
        h_x, h_y, p_x, p_y = draw_protein(protein, axs[0,0])
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
        axs[0,0].legend()
        axs[1,0].clear()
        axs[1,0].plot(energies)
        axs[1,1].clear()
        axs[1,1].plot(densities, color="#6cad50")
        axs[0,1].clear()
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

def landscape_view():
    sequence = "NLYIQWLKDGGPSSGRPPPS"
    protein1 = fold.Protein(sequence, 2, True)
    protein2 = fold.Protein(sequence, 2, False)
    energies_a, energies_b = [], []
    for i in range(3000):
        protein1.update()
        energies_a.append(protein1.score)
        protein2.update()
        energies_b.append(protein2.score)
    plt.plot(energies_a + list(reversed(energies_b)))
    plt.show()
    fig, ax = plt.subplots(1,1);
    draw_protein(protein1, ax)
    plt.show()

#standard_view(100)
#landscape_view()

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
def threedimview():
    sequence = "NLYIQWLKDGGPSSGRPPPS"
    protein = fold.Protein(sequence, 2, False)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    def animate(test):
        protein.update()
        ax.clear()
        #ax.set_ylim([-len(protein.residues)/3, len(protein.residues)+len(protein.residues)/3])
        #ax.set_xlim([-len(protein.residues)/3, (len(protein.residues)+len(protein.residues)/3)/2])
        h_x, h_y, h_z, p_x, p_y, p_z = [], [], [], [], [], []
        chain, bonds = [], []
        for x,i in enumerate(protein.residues):
            if i.polar == True:
                p_x.append(i.coords[0])
                p_y.append(i.coords[1])
                p_z.append(i.coords[2])
            else:
                h_x.append(i.coords[0])
                h_y.append(i.coords[1])
                h_z.append(i.coords[2])
            for y,j in enumerate(protein.residues):
                if (abs(y-x) > 1 and (np.linalg.norm(np.array(i.coords)-np.array(j.coords)) == 1) and j.polar == False and i.polar == False):
                    bonds.append([i.coords, j.coords])
            if (x != len(protein.residues)-1):
                chain.append([i.coords, protein.residues[x+1].coords])
        cc = Line3DCollection(chain, color='black', linewidths=2, label="Chain connection")
        bc = Line3DCollection(bonds, color='red', linewidths=2, label="Bond")
        ax.add_collection(cc)
        ax.add_collection(bc)
        ax.scatter(h_x, h_y, h_z, zorder=2, label='Hydrophobic')
        ax.scatter(p_x, p_y, p_z, zorder=2, label='Polar')
    ani = animation.FuncAnimation(fig, animate, interval=1) 
    plt.show()

threedimview()
