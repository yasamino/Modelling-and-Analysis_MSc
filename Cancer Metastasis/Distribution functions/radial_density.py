

def radial_density_distribution(f, th, size , time , delta_r = 1  , write_to_file = False , Plot = False , Json_file = 'Radial_density_distribution.json'):
    m = np.vstack(th.cellInd)
    time = time * 10000

    s1 = np.where(m < 0)[0]
    s = [ f[i][:180] for i in s1 ] #coordinates of nodes on each soft cell
    if s:
        soft_coords = np.vstack(s)
        soft_com = []
        for mi in range(int(len(soft_coords)/180)):
            soft_com.append(np.mean(soft_coords[mi*180:(mi+1)*180],axis=0))
    soft_com = np.array(soft_com)[: , 0:2] #2D coordinates
    print("Number of soft cells: ",len(soft_com))

    s2 = np.where(m >= 0 )[0]
    r = [ f[i][:180] for i in s2 ] #coordinate of nodes on each Rigid cell
    hard_coords = np.vstack(r)
    hard_com = []
    for mi in range(int(len(hard_coords)/180)):
        hard_com.append(np.mean(hard_coords[mi*180:(mi+1)*180],axis=0))
    hard_com = np.array(hard_com)[: , 0:2] 
    print("Number of hard cells: ",len(hard_com))

    all_coords = np.vstack(f)
    all_com = []
    for mi in range(int(len(all_coords)/180)):
        all_com.append(np.mean(all_coords[mi*180:(mi+1)*180],axis=0))
    all_com = np.array(all_com)[: , 0:2] 
    print("Number of all cells: ",len(all_com))

    center = np.array([size/2, size/2])
    # find the radial distribution of the cells in intervals of delta_r
    density_soft = [0]*int(size/(delta_r))
    density_hard = [0]*int(size/(delta_r))
    density_all = [0]*int(size/(delta_r))

    for i in range(len(soft_com)):
        r = np.linalg.norm(soft_com[i] - center)
        density_soft[int(r/(delta_r))] += 1   
    density_soft = np.array(density_soft)/len(all_com)

    for i in range(len(hard_com)):
        r = np.linalg.norm(hard_com[i] - center)
        density_hard[int(r/(delta_r))] += 1
    density_hard = np.array(density_hard)/len(all_com)

    for i in range(len(all_com)):
        r = np.linalg.norm(all_com[i] - center)
        density_all[int(r/(delta_r))] += 1
    density_all = np.array(density_all)/len(all_com)

    if Plot:
        plt.plot(np.linspace(0, size, len(density_soft)), density_soft, label=r'$\rho_{Soft}$')
        plt.plot(np.linspace(0, size, len(density_hard)), density_hard, label=r'$\rho_{Hard}$')
        plt.plot(np.linspace(0, size, len(density_all)), density_all, label=r'$\rho_{All}$')
        plt.xlabel('r')
        plt.ylabel('Density')
        plt.ylim(0,0.05)
        plt.title('Radial Density Distribution - time = {} '.format(time))
        plt.legend()
        plt.savefig('Radial_density_distribution_{}.png'.format(time))
        plt.clf()

    if write_to_file:
        if os.path.isfile(Json_file):
            with open(Json_file) as file:
                data = json.load(file)
            data['Soft'].append(density_soft.tolist())
            data['Hard'].append(density_hard.tolist())
            data['All'].append(density_all.tolist())
            data['Time'].append(time)
            with open(Json_file, 'w') as file:
                json.dump(data, file)
        else:
            print(density_soft.tolist())
            data = {
                'Soft': [density_soft.tolist()],
                'Hard': [density_hard.tolist()],
                'All': [density_all.tolist()],
                'Time': [time]
            }
            with open(Json_file, 'w') as file:
                json.dump(data, file)

        with open(Json_file) as file:
            data = json.load(file)
        data['Radius'] = np.linspace(0, size, len(density_soft)).tolist()
        with open(Json_file, 'w') as file:
            json.dump(data, file)

def Animate_radius_cluster(size =200, size_z=12 , z_grid = Unit_length):
    '''
    This function reports the Radius of gyration of the cell cluster
    It reports the values using three dfferent methods:
    1) The maximum radius of cells in each z grid layer
    2) The radius of gyration of the whole cell cluster
    3) The radius of gyration of the outer layer of the cell cluster
    '''

    def z_position_report(f,th,frame_number , Filename , size , ax):
        all_coords = np.vstack(f)
        all_coords = np.vstack(all_coords)
        # Just keep the z coordinates
        COM = np.mean(all_coords, axis=0)
        Rg_grid = np.arange(0, size_z, z_grid)
        N = len(Rg_grid)
        com_coords = []
        r2_grid =[ [] for _ in range(N)]
        r2_grid_max = np.zeros(N)
        r2_coords = [[] for _ in range(N)]
        for mi in range(int(len(all_coords)/180)):
            com = np.mean(all_coords[mi*180:(mi+1)*180],axis=0)
            index = np.searchsorted(Rg_grid, com[2]) #sorts the cells into appropriate place in the grid
            # print("COM Z : {}, index:{}".format(com[2], index-1))
            r2_coords[index-1].append([com[0] - size/2 , com[1] - size/2])
            r2_grid[index-1].append((com[0] - size/2)**2+(com[1] - size/2)**2)
            r2_grid_max[index-1] = np.max( [r2_grid_max[index-1] , np.sqrt((com[0] - size/2)**2+(com[1] - size/2)**2) ] )

        Rg = np.sqrt( [np.mean(row) for row in r2_grid])
        Rg_list_only_outer_layer = [ [] for _ in range(N)]
        for i in range(N):
            if len(r2_coords[i])>2:
                hull = ConvexHull(r2_coords[i])
                r2_coords_numpy = np.array(r2_coords[i])
                Rg_list_only_outer_layer[i].append( r2_coords_numpy[hull.vertices,0]**2 + r2_coords_numpy[hull.vertices,1]**2 )

        Rg_only_outer_layer = np.sqrt( [ np.mean(row)  for row in Rg_list_only_outer_layer] )
            
        plt.title("Raduis of cell cluster, Time = {}".format(frame_number*1000))
        plt.xlabel('Z[$\\sigma$]')
        plt.ylabel('radius[$\\sigma$]')
        plt.ylim(0,60/Unit_length)
        plt.xlim(0,size_z/Unit_length)
        plt.plot(Rg_grid/Unit_length,Rg/Unit_length, label = "Rg, ROC: {}".format( np.round( Rg[N//2]/Unit_length ,2) ))
        plt.plot(Rg_grid/Unit_length, Rg_only_outer_layer/Unit_length , label = "Rg outer layer, ROC: {}".format(np.round( Rg_only_outer_layer[N//2]/Unit_length , 2)))
        plt.plot(Rg_grid/Unit_length, r2_grid_max/Unit_length , label = "Maximum radius, ROC: {}".format(np.round( r2_grid_max[N//2]/Unit_length ,2) ))
        
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        #plt.hist(z_com, bins=100)
        return plt.plot()
    
    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return z_position_report(frame, th, frame_number, filename, size, ax)

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            fig, ax = plt.subplots()
            #create an animation
            plt.xlabel('Distance from the sheets')
            plt.ylabel('Number of cells')
            plt.ylim(0,100)
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(0,th.maxFrames,nSkip), interval=30, repeat=False)
            #save animation as gif
            ani.save('Radius_{}.gif'.format(filename[:-4]), writer='Pillow', fps=10)

def Animate_distribution_of_cells_radial(size = 200, size_z = 12 , grid_size = 1 , z_report = None):
    '''
    This function reports the distribution of cells in the radial direction
    inputs:
    size: size of the box in x and y direction
    size_z: size of the box in z direction
    '''

    if (z_report == None): # if no z values are given, then report the distribution at the center and the 1/4 and 3/4 point in the box
        Z_report = [1/4 , 1/2 ,  3/4]


    def z_position_report(f,th,frame_number , Filename , size , ax):
        all_coords = np.vstack(f)
        all_coords = np.vstack(all_coords)
        COM = np.mean(all_coords, axis=0)
        COM = COM/Unit_length # convert to sigma units

        Rg_grid = np.arange(0, size_z, grid_size)
        radial = np.arange(0, size/2, grid_size)
        N = len(Rg_grid)
        R = len(radial)
        r2_coords = np.zeros((N,R))

        for mi in range(int(len(all_coords)/180)): # for each cell in the frame
            com = np.mean(all_coords[mi*180:(mi+1)*180],axis=0)
            com = np.array([com[0] - size/2 , com[1] - size/2 , com[2]])
            index_height = np.searchsorted(Rg_grid, com[2])
            radius = np.sqrt(com[0]**2 + com[1]**2)
            index_radius = np.searchsorted(radial, radius)
            r2_coords[index_height-1][index_radius-1] += 1
        
        plt.title("Arial Density distribution of the cells, Time = {}".format(frame_number*1000))
        plt.xlabel('Radius[$l^*$]')
        plt.ylabel('Cell Density(#Area) [$l^{-2}$]')
        for z_report in Z_report:
            r_density = r2_coords[ int(N*z_report) ] / (np.pi * ( (radial+grid_size)**2 - radial**2 ))
            plt.plot(radial,  r_density, label = "Z = {}".format(z_report*size_z))
        plt.legend()
        plt.ylim(0,1)
        # plt.gca().set_aspect('equal', adjustable='box')
        return plt.plot()

    
    
    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return z_position_report(frame, th, frame_number, filename, size, ax)

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            fig, ax = plt.subplots()
            #create an animation
            plt.xlabel('Distance from the sheets')
            plt.ylabel('Number of cells')
            plt.ylim(0,100)
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(0,th.maxFrames,nSkip), interval=30, repeat=False)
            #save animation as gif
            ani.save('Radial_Distribution_{}.gif'.format(filename[:-4]), writer='Pillow', fps=10)
