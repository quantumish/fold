import matplotlib.pyplot as plt

h = [(0,0), (1,0)]
p = [(0,1), (1,1)]
plt.scatter(*zip(*h))
plt.scatter(*zip(*p))
plt.ylim(-1, 2)
plt.xlim(-1, 2)
plt.show()
