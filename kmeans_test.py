import matplotlib.pyplot as plt
from kmeansclustering import Cell
from kmeansclustering import kpp_means_cluster
import random
import numpy as np

# Test datasets obtained from https://cs.joensuu.fi/sipu/datasets/

# Read the coordinates from the text file
x_coords = []
y_coords = []

cells: list[Cell] = []

with open('s1.txt', 'r') as file: # adjust k, colormap, and cell initial shortest distance references depending on dataset.
    for line in file:
        # Split each line by space and convert to float
        x, y = map(float, line.split())
        x_coords.append(x)
        y_coords.append(y)
        cells.append(Cell((x,y)))

# Plot the coordinates
plt.figure(figsize=(8, 6))
plt.scatter(x_coords, y_coords, c='blue', marker='o')
plt.title('Initial Plot')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()

clustered_cells = kpp_means_cluster(cells,15)

x = [cell.x for cell in clustered_cells]
y = [cell.w for cell in clustered_cells]
clusters = [cell.cluster.number for cell in clustered_cells]

cmap = plt.get_cmap('tab20')  # Use tab20 or any suitable colormap with at least 15 colors
colors = cmap(np.linspace(0, 1, 15))

# Plot each cluster with a different color
plt.figure(figsize=(8, 6))
for i in range(0,len(x)):
    color = colors[clusters[i]]
    plt.scatter(x[i], y[i], c=color, marker='o')

plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Cells by Cluster')
plt.legend()
plt.show()
