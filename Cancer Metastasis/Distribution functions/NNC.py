import numpy as np


def NN_distribution(f,th,Frame,Filename):
    # Nearest neighbor distribution for soft, hard, and all cells
    m = np.vstack(th.cellInd)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]
    h = [ f[i][:180] for i in s1 ]
    g = [ f[i][:180] for i in s2 ]
    all_coords = np.vstack(f)
    all_coords = np.vstack(all_coords)
    all_com = []
    for mi in range(int(len(all_coords)/180)):
        all_com.append(np.mean(all_coords[mi*180:(mi+1)*180],axis=0))
    all_com = np.array(all_com)[: , 0:2] #2D coordinates
    print("Total number of cells : ",len(all_com))
    vor = Voronoi(all_com)
    #create list of neighbors for each cell
    neighbors = [ [] for _ in range(len(all_com)) ]
    for i in range(len(vor.ridge_points)):
        neighbors[vor.ridge_points[i][0]].append(vor.ridge_points[i][1])
        neighbors[vor.ridge_points[i][1]].append(vor.ridge_points[i][0])

    NN = []
    for Neigh in neighbors:
        NN.append(len(Neigh))
    maX = np.max(NN)
    #maX=18
    Neighbors = [0] * (maX+1)
    for i in NN:
        Neighbors[i]+=1
    Neighbor = np.array(Neighbors) / sum(Neighbors)
    x_axis = [i for i in range(maX+1)]
    
    plt.plot(x_axis, Neighbor , marker = 'o' , linestyle = '--')
    plt.xlabel('Number of Neighbors')
    plt.ylabel('Frequency')
    plt.title('CNN plot, time = {}, DC = {}%'.format(Frame*1000 , np.round(Degree_of_Confluence(f,th)*100,2)) )
    plt.savefig('{}.png'.format(Filename))
    print(Neighbor)
    plt.clf()
