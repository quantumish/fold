import matplotlib.pyplot as plt
import fold

protein = fold.Protein("HPPHPH")
energies = []
for i in range(300):
    protein.update()
    energies.append(protein.score)
h_x = []
h_y = []
p_x = []
p_y = []
for i in protein.residues:
    if i.polar == True:
        p_x.append(i.coords[0])
        p_y.append(i.coords[1])
    else:
        h_x.append(i.coords[0])
        h_y.append(i.coords[1])

fig, axs = plt.subplots(1, 2)
axs[0].scatter(h_x, h_y)
axs[0].scatter(p_x, p_y)
axs[1].plot(energies)
axs[0].set_ylim([-5, 10])
axs[0].set_xlim([-5, 10])
plt.show()

