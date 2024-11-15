# This file implements k++ means clustering for 2D points. 
# The user is expected to use the kpp_means_cluster function and the Cell class for creating 2D points.
# The user is not expected to create Cluster objects, they may only need to get its 'number' property to find which cluster number a cell belongs to.
# The distance_between_2coords function is used internally.

import random
import math

def distance_between_2coords(a,b) -> float:
    a0 = a[0]; a1 = a[1]; b0 = b[0]; b1 = b[1]
    term1 = (b0 - a0)*(b0 - a0)
    term2 = (b1 - a1)*(b1 - a1)
    return math.sqrt(term1 + term2)

class Cluster:
    def __init__(self,number:int, x:int, w:int) -> None:
        self.number = number # for identification
        self.x = x # centroid x-coordinate
        self.w = w # centroid y-coordinate
        self.cells: list[Cell] = [] # cells mapped to this cluster
    def __str__(self) -> str:
        cells_str = ''
        for cell in self.cells:
            cells_str += '(' + str(cell.x) + ',' + str(cell.w) + '), '
        return f'Cluster {self.number}: {self.x},{self.w} | {cells_str}'
    def relocate(self, max_x=0, max_w=0) -> bool:
        """Performs the relocation step of the centroid in the k++ algorithm.
        Returns true if the location of the centroid had to be changed.

        Args:
            max_x (int/float): This was used to specify the max bounds when initiazling centroid using vanilla k-means, but since we are using k++ means now, this is not needed. Hence i set default values for it.
            max_w (int/float): This was used to specify the max bounds when initiazling centroid using vanilla k-means, but since we are using k++ means now, this is not needed. Hence i set default values for it.

        Returns:
            bool: Returns true if location of the centroid had to be changed.
        """
        x_total = 0
        w_total = 0
        count = len(self.cells)
        if count == 0: # shouldn't be reached since we are using k++ initialization of the cluster centers
            self.x = random.uniform(float(0),float(max_x))
            self.w = random.uniform(float(0),float(max_w))
            return True
        for cell in self.cells:
            x_total += cell.x
            w_total += cell.w
        new_x = x_total/count
        new_w = w_total/count
        change = distance_between_2coords((self.x,self.w),(new_x,new_w))
        self.x = new_x
        self.w = new_w
        if change != 0.0: return True
        return False

class Cell:
    def __init__(self,coord) -> None:
        self.x = coord[0]
        self.w = coord[1]
        self.cluster: Cluster = None # which cluster cell is allocated to
        self.shortest_distance = 1000000000000
    def __str__(self) -> str:
        if not self.cluster:
            return f'({self.x},{self.w}) - Cluster unallocated | distance: {self.shortest_distance}'
        return f'({self.x},{self.w}) - Cluster {self.cluster.number} | distance: {self.shortest_distance}'
    def update_cluster(self, new_cluster: Cluster, distance: float):
        """Private method used by the kpp_means_cluster function.
        This function is used to allocate a new cluster to the cell, and it also updates the old and new cluster objects.

        Args:
            new_cluster (Cluster): 
            distance (float): 
        """
        if self.cluster:
            self.cluster.cells.remove(self)
        self.cluster = new_cluster
        new_cluster.cells.append(self)
        self.shortest_distance = distance
    def get_nearest_centroid(self, clusters: list[Cluster]) -> Cluster:
        shortest_distance = 1000000000000
        nearest_centroid = None
        for cluster in clusters:
            distance = distance_between_2coords((self.x,self.w),(cluster.x,cluster.w))
            if distance <= shortest_distance:
                nearest_centroid = cluster
                shortest_distance = distance
        return nearest_centroid
    


def kpp_means_cluster(cells: list[Cell], k: int) -> list[Cell]:
    """Performs k++ means clustering on the given Cell list.
    The input cells are expected to be new.

    Args:
        cells (list[Cell]): List of NEW Cell objects.
        k (int): 

    Returns:
        list[Cell]: Same cell list, with the cluster property updated.
    """
    # get an idea of the x and w range. This is no longer needed as we use k++ means and not vanilla k-means.
    # max_x = 0; max_w = 0
    # for cell in cells:
    #     max_x = max(max_x, cell.x)
    #     max_w = max(max_w, cell.w)

    # initialize cluster centers
    clusters: list[Cluster] = []
    random_cell = random.choice(cells)
    first_cluster = Cluster(0,random_cell.x,random_cell.w)
    random_cell.update_cluster(first_cluster, 0)
    clusters.append(first_cluster)
    for n in range(1,k):
        distances: list[Cell, float] = []
        for cell in cells:
            if cell.cluster: continue
            nearest_centroid = cell.get_nearest_centroid(clusters)
            dist = distance_between_2coords((cell.x,cell.w),(nearest_centroid.x,nearest_centroid.w))
            distances.append([cell,dist])
        distances = sorted(distances, key=lambda x: x[1])
        sum_of_distances = 0.0
        for distance in distances:
            distance[1] *= distance[1] # square each distance
            sum_of_distances += distance[1]
        random_distance = random.uniform(0.0, sum_of_distances)
        sum_of_distances = 0
        for distance in distances:
            sum_of_distances += distance[1]
            if sum_of_distances >= random_distance:
                new_centroid_location_cell = distance[0]
                break
        new_cluster = Cluster(n,new_centroid_location_cell.x,new_centroid_location_cell.w)
        new_centroid_location_cell.update_cluster(new_cluster, 0)
        clusters.append(new_cluster)
        
    
    # print('Initial cluster centers:')
    # for cluster in clusters:
    #     print(cluster)
    # print('Loop begin')
    # print()

    while True:
        # assign cells to nearest cluster
        for cell in cells:
            for cluster in clusters:
                distance = distance_between_2coords((cell.x,cell.w),(cluster.x,cluster.w))
                if distance < cell.shortest_distance:
                    if cell.cluster != cluster: cluster_changed = True
                    cell.update_cluster(cluster, distance)

        # relocate cluster center
        is_change = False
        for cluster in clusters:
            change = cluster.relocate()
            if change: is_change = True

        if not is_change: break
        
        # print('Cells:')
        # for cell in cells:
        #     print(cell)
        # print('Centers:')
        # for cluster in clusters:
        #     print(cluster)
        # print()
        # print()
        # time.sleep(1)

    # I am not sure how Python deals with strong reference cycles, so I'm doing this just in case.
    for cluster in clusters:
        cluster.cells = None

    return cells


# points = [
#     (0,2),
#     (2,0),
#     (2,2),
#     (2,3),
#     (7,8),
#     (8,7),
#     (7,3),
#     (8,3),
#     (7,2)
# ]

# cells: list[Cell] = []
# for point in points:
#     cells.append(Cell(point))

# print('Initial cells:')
# for cell in cells:
#     print(cell)

# kpp_means_cluster(cells)

        


