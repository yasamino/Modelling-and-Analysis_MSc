def Animate_CNN(cell_radius = 0.63 , size_z = 30 , slice = "bottom"):  
    '''
    Function to animate cells and color them based on the Number of Nearest Neighbours
    '''
    start_from =0
    end = 0
    size=200
    Cutoff_black = 1
    
    
    def z_coords_C_cutoff_radius(f,th,frame_number , Filename , size , ax , z = 1.9):
        Cutoff_z= 2*cell_radius +0.5
        cut_off_neighboors = 5*cell_radius
        cut_off_NN_criteria = np.sqrt(np.power(2*cell_radius , 3))
        m = np.vstack(th.cellInd)
        c = [ f[i][:180] for i in range(len(m)) ]
        number_of_neighbors = []
        Cutoff_number_NN = 3

        # Filter cells that are in the first layer
        if slice == "bottom":
            COM_layer_1 = [np.mean(cell[:,0:3], axis=0) for cell in c if np.mean(cell[:,2], axis=0)< Cutoff_z]
            COM_layer_2_for_NN = [np.mean(cell[:,0:3], axis=0) for cell in c if np.mean(cell[:,2], axis=0)< cut_off_neighboors] # cutting another layer to calculate the number of nearest neighbors
        elif slice == "top":
            COM_layer_1 = [np.mean(cell[:,0:3], axis=0) for cell in c if np.mean(cell[:,2], axis=0)> size_z - Cutoff_z]
            COM_layer_2_for_NN = [np.mean(cell[:,0:3], axis=0) for cell in c if np.mean(cell[:,2], axis=0)> size_z - cut_off_neighboors]


        for cell in COM_layer_1:
            if COM_layer_2_for_NN:
                distances = [np.linalg.norm(cell - other_cell) for other_cell in COM_layer_2_for_NN]
            else:
                distances = [np.linalg.norm(cell - other_cell) for other_cell in COM_layer_1]
            number_of_neighbors.append(len([i for i in distances if i < Cutoff_z]))

        cm1 = mcol.LinearSegmentedColormap.from_list("CMAP",["b","r"])
        cnorm = mcol.Normalize(vmin=np.min(number_of_neighbors),vmax= np.max(number_of_neighbors))
        cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
        cpick.set_array([])
        count = 0
        for index , cell in enumerate(c):
            patches =[]
            COM_Z = np.mean(cell[:,2], axis=0) 
            if slice == "bottom":
                if (COM_Z < Cutoff_z): #if it falls within the first layer
                    
                    hull = ConvexHull(cell[:,0:2])
                    p = cell[hull.vertices,0:2]
                    patches.append(Polygon(p, closed=True))
                    p = PatchCollection(patches,edgecolor="r",alpha=0.3)
                    p.set_color(cpick.to_rgba(number_of_neighbors[count]))
                    if Cutoff_black:
                        if number_of_neighbors[count]<= Cutoff_number_NN:
                            p.set_color("k")
                    ax.add_collection(p)
                    count+=1
            elif slice == "top":
                if (COM_Z > size_z - Cutoff_z):
                    hull = ConvexHull(cell[:,0:2])
                    p = cell[hull.vertices,0:2]
                    patches.append(Polygon(p, closed=True))
                    p = PatchCollection(patches,edgecolor="r",alpha=0.3)
                    
                    if Cutoff_black and number_of_neighbors[count]<= Cutoff_number_NN:    
                        p.set_color('#000000')
                    else:
                        p.set_color(cpick.to_rgba(number_of_neighbors[count]))
                    ax.add_collection(p)
                    count+=1

        # print('Number of neighbors:', number_of_neighbors)
        ax.set_xlim([60,size-60])
        ax.set_ylim([60,size-60])
        if slice == "bottom":
            ax.set_title("CNN distribution in the bottom layer, Time = {}".format(frame_number*1000))
        elif slice == "top":
            ax.set_title("CNN distribution in the top layer, Time = {}".format(frame_number*1000))
        print(frame_number)
        return ax.plot()

    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return z_coords_C_cutoff_radius(frame, th, frame_number , filename, size, ax , size_z) #### set the size

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            fig, ax = plt.subplots(figsize = (size/10, size/10))  
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_aspect('equal')

            cm1 = mcol.LinearSegmentedColormap.from_list("CMAP",["b","r"])
            cnorm = mcol.Normalize(vmin=0,vmax= 12)
            cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
            cpick.set_array([])
            #create an animation
            for _ in range(start_from):
                th.ReadFrame(inc=nSkip)
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(start_from,th.maxFrames-end,nSkip), interval=50, repeat=False)
            cbar = plt.colorbar(cpick, ax=ax , pad =0.1)
            cbar.set_label("CNN", labelpad=10)  # Main label (on the side)
            
            # plt.show()
            # save animation as gif
            if slice == "bottom":
                ani.save('Color_Neighbors_bottom_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)
            elif slice == "top":
                ani.save('Color_Neighbors_top_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)
