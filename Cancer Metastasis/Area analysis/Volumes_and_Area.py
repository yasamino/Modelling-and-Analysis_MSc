import numpy as np
from scipy.spatial import ConvexHull


def volumes(f,th): 
    ''' Recieves the trajectory of nodes belonging to each cell, and calculates the Area and Volume of the cells using convex hull
    input: 
    f: The array containing the x, y, z position of nodes'''
    volume = []
    Area = []
    m = np.vstack(th.cellInd)
    c = [ f[i][:180] for i in range(len(m)) ]
    for cell in c:
        hull = ConvexHull(cell[:,0:3])
        p = cell[hull.vertices,0:3]
        Area.append(hull.area) #returns the area for 3D case
        volume.append(hull.volume) #returns the volume for 3D case
    print("Average volume of cells =  ", np.mean(volume), "Max: ", np.max(volume), "Min: ", np.min(volume))
    print("Average area of cells = ", np.mean(Area), "Max: ", np.max(Area), "Min: ", np.min(Area))