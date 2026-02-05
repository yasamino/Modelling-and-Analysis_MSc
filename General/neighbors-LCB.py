import numpy as np


def z_coords_C_VORONOI_FAIL(f,th,frame_number , Filename , size , ax , z = 1.9):
    fig = plt.figure()
    ax = fig.add_subplot()
    m = np.vstack(th.cellInd)
    c = [ f[i][:180] for i in range(len(m)) ]
    all_com_in_layer = []

    hull_list = []
    # Filter cells that are in the first layer
    for index , cell in enumerate(c):
        COM_Z = np.mean(cell[:,2], axis=0) 
        if (COM_Z < 2*cell_radius + 0.5): #if it falls within the first layer
            all_com_in_layer.append(np.mean(cell[:,0:2], axis=0)) # 2D coordinates of cells in the layer
            hull_list.append(ConvexHull(cell[:,0:2]))
    vor = Voronoi(all_com_in_layer)
    number_of_neighbors = []
    patches = []
    # len([x for x in arr if x > 0])
    for points in vor.regions:
        number_of_neighbors.append(len([i for i in points if i > -1]))

    cm1 = mcol.LinearSegmentedColormap.from_list("CMAP",["b","k"])
    cnorm = mcol.Normalize(vmin=0,vmax= np.max(number_of_neighbors))
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    for index , cell in enumerate(c):
        COM_Z = np.mean(cell[:,2], axis=0) 
        if (COM_Z < 2*cell_radius + 0.5): #if it falls within the first layer
            hull = ConvexHull(cell[:,0:2])
            p = cell[hull.vertices,0:2]
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches,edgecolor="r",alpha=0.4)
            p.set_color(cpick.to_rgba(number_of_neighbors[-1]))
            ax.add_collection(p)

    print('Number of neighbors:', number_of_neighbors)
    ax.set_xlim([80,125])
    ax.set_ylim([80,125])
    ax.set_title("CNN distribution in the bottom layer, Time = {}".format(frame_number*1000))
    plt.colorbar(cpick, ax=ax , pad =0.1)
    plt.show()
    return ax.plot()

def z_coords_C_cutoff_radius(f,th,frame_number , Filename , size , ax , z = 1.9):
    Cutoff_z= 2*cell_radius +0.5

    fig = plt.figure()
    ax = fig.add_subplot()
    m = np.vstack(th.cellInd)
    c = [ f[i][:180] for i in range(len(m)) ]
    all_com_in_layer = []
    number_of_neighbors = []

    # Filter cells that are in the first layer
    COM_layer_1 = [np.mean(cell[:,0:3], axis=0) for cell in c if np.mean(cell[:,2], axis=0)< Cutoff_z]
    COM_layer_2_for_NN = [np.mean(cell[:,0:3], axis=0) for cell in c if np.mean(cell[:,2], axis=0)< (5*cell_radius)] # cutting another layer to calculate the number of nearest neighbors
    for cell in COM_layer_1:
        distances = [np.linalg.norm(cell - other_cell) for other_cell in COM_layer_2_for_NN]
        number_of_neighbors.append(len([i for i in distances if i < 2*cell_radius + 0.5]))
    cm1 = mcol.LinearSegmentedColormap.from_list("CMAP",["b","r"])
    cnorm = mcol.Normalize(vmin=np.min(number_of_neighbors),vmax= np.max(number_of_neighbors))
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    count = 0
    for index , cell in enumerate(c):
        patches =[]
        COM_Z = np.mean(cell[:,2], axis=0) 
        if (COM_Z < Cutoff_z): #if it falls within the first layer
            
            hull = ConvexHull(cell[:,0:2])
            p = cell[hull.vertices,0:2]
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches,edgecolor="r",alpha=0.3)
            p.set_color(cpick.to_rgba(number_of_neighbors[count]))
            print(number_of_neighbors[count])
            ax.add_collection(p)
            count+=1

    print('Number of neighbors:', number_of_neighbors)
    ax.set_xlim([0,200])
    ax.set_ylim([0,200])
    ax.set_title("CNN distribution in the bottom layer, Time = {}".format(frame_number*1000))
    plt.colorbar(cpick, ax=ax , pad =0.1)
    plt.show()
    return ax.plot()

def CellPlot_3D_colored(f, th, Neighcutoff , frame_number, Filename , size=200, size_z=12):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


    '''cell node coordinates soft (h) - hard(g)'''
    m = np.vstack(th.cellInd)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]

    h = [ f[i][:180] for i in s1 ]
    for cell in h:
        ax.scatter(cell[:,0], cell[:,1], cell[:,2], c='b', marker='o')

    g = [ f[i][:180] for i in s2 ]
    for cell in g:
        ax.scatter(cell[:,0], cell[:,1], cell[:,2], c='r', marker='o')


    ax.set_xlim([0,size])
    ax.set_ylim([0,size])
    ax.set_zlim([0,size_z])
    plt.title("3D projection - Cells with less than {} Neighbors - Time = {}".format( Neighcutoff, (frame_number-2)*10000))
    plt.savefig("threeDimPlot_{}_{}.png".format(Filename, frame_number))
    plt.clf()


      