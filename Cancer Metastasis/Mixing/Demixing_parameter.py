import numpy as np
from scipy.spatial import Voronoi

def deMixing_ratio(f,th):
    m = np.vstack(th.cellInd)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]
    h = [ f[i][:180] for i in s1 ]
    len_type1 = len(h)  
    g = [ f[i][:180] for i in s2 ]
    len_type2 = len(g)

    if len_type1 == 0:
        all_coords = np.vstack(g)
    elif len_type2 == 0:
        all_coords = np.vstack(h)
    else:
        all_coords = np.concatenate((h,g), axis=0)

    all_coords = np.vstack(all_coords)
    '''create list of center of mass, form the voronoi diagram, find the neighbors of each cell'''
    all_com = []
    for mi in range(int(len(all_coords)/180)):
        all_com.append(np.mean(all_coords[mi*180:(mi+1)*180],axis=0))
    all_com = np.array(all_com)[: , 0:2]
    vor = Voronoi(all_com)
    '''create list of neighbors for each cell'''
    neighbors = [ [] for _ in range(len(all_com)) ]
    for i in range(len(vor.ridge_points)):
        neighbors[vor.ridge_points[i][0]].append(vor.ridge_points[i][1])
        neighbors[vor.ridge_points[i][1]].append(vor.ridge_points[i][0])

    Dp=0
    for cell in range(len(all_com)):
        Ns = 0
        ''' Check to see if cell index belongs in the range of type 1'''
        if cell < len_type1: 
            #it is a type 1 cell
            for neighbor in neighbors[cell]:
                if neighbor < len_type1:
                    Ns+=1
        else:
            #it is a type 2 cell
            for neighbor in neighbors[cell]:
                if neighbor >= len_type1:
                    Ns+=1

        Nt= len(neighbors[cell]) #total number of neighbors
        Dp+= (2*(Ns/Nt) - 1) #deMixing ratio

    Dp = Dp/len(all_com)
    return Dp
