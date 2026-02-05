import numpy as np
import matplotlib.pyplot as plt

def Cell_Shape_Parameter(f , frame_number , Filename , size , plot_system = False , Plot_distribution = False , bin_n =50):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    order_parameter_soft = []
    order_parameter_hard = []

    '''cell node coordinates soft (h) - hard(g)'''
    m = np.vstack(th.cellInd)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]

    Order_param_shphere = 3.54
    Order_param_shphere = 3.7

    h = [ f[i][:180] for i in s1 ]
    for cell in h:
        patches=[]
        hull = ConvexHull(cell[:,0:2])
        p = cell[hull.vertices,0:2]
        patches.append(Polygon(p, closed=True))
        ord_par = hull.area/ np.sqrt(hull.volume)
        order_parameter_soft.append(ord_par)

        if plot_system:
            if ord_par > Order_param_shphere:   
                p = PatchCollection(patches, color= (1, 128/255, 0 ) ,alpha=0.8)
            else:
                p = PatchCollection(patches, color = (0.2, 0.6, 1) ,edgecolor="b",alpha=0.8)
            ax.add_collection(p)

    g = [ f[i][:180] for i in s2 ]
    for cell in g:
        patches=[]
        hull = ConvexHull(cell[:,0:2])
        p = cell[hull.vertices,0:2]
        patches.append(Polygon(p, closed=True))
        ord_par = hull.area/ np.sqrt(hull.volume)
        order_parameter_hard.append(ord_par)

        if plot_system:
            if ord_par > Order_param_shphere:
                p = PatchCollection(patches, color= (1, 128/255, 0 ) ,alpha=0.8)
            else:
                p = PatchCollection(patches, color= (222/255, 28/255, 0 ) ,alpha=0.8)
            ax.add_collection(p)
              
    if plot_system:
        ax.set_xlim([0,size])
        ax.set_ylim([0,size])
        plt.title("Cell Shape parameter p $ > ${} , Time = {}".format( Order_param_shphere , (frame_number-2)*10000))
        plt.savefig("Cell_shape_param_{}_{}.png".format(Filename, frame_number))
        plt.clf()

    if Plot_distribution:
        counts, bins, _ = plt.hist(order_parameter_soft + order_parameter_hard , bins=bin_n, edgecolor='black', alpha=0)  # alpha=0 makes the bars invisible
        fig = plt.figure()
        total_samples = np.sum(counts)
        counts = np.array(counts)
        normalized_counts = counts / total_samples *100
        plt.scatter(bins[:-1], normalized_counts)

        x= bins[:-1]
        y= normalized_counts
        print('x:',x)
        print('y:',y)

        gauss = lambda x, sigma2_times_2_rev , mu: np.sqrt(sigma2_times_2_rev/np.pi) * np.exp(-1*sigma2_times_2_rev * (x-mu)**2)
        exponential = lambda x, a, b: a*np.exp(-b*x)
        percent = 0.10
        gauss_x = x[:int(percent*bin_n)]
        gauss_y = np.array(y)[:int(percent*bin_n)]

        exp_x = x[int(percent*bin_n):]
        exp_y = np.array(y)[int(percent*bin_n):]

        parameter_gauss, covarience_gauss = curve_fit(gauss, gauss_x, gauss_y , maxfev=6000)
        parameter_exp , covarience_exp = curve_fit(exponential, exp_x, exp_y , maxfev=6000)
        #parameter_gauss_all , covarience_gauss_all = curve_fit(gauss, x, np.array(y)/len(array) *100 , maxfev=5000)

        y_data_fit_gauss = gauss(gauss_x, parameter_gauss[0], parameter_gauss[1])
        y_data_fit_exp = exponential(exp_x, parameter_exp[0], parameter_exp[1])
        #y_data_fit_gauss_all = gauss(x, parameter_gauss_all[0], parameter_gauss_all[1])

        #plt.plot(gauss_x, y_data_fit_gauss, label='$\\frac{1}{\sqrt{2 \pi  %5.3f}} e^{- \\frac{(\\Upsilon - %5.3f )^2}{%5.3f}} $' % (1/(parameter_gauss[0]*2) ,parameter_gauss[1] ,1/(parameter_gauss[0]*2)) , color = 'g')
        plt.plot(gauss_x, y_data_fit_gauss, label='$ \\sigma^2 =  %5.3f , \\mu = %5.3f $'%( 1/(parameter_gauss[0]*2) , parameter_gauss[1]) , color = 'g')
        plt.plot(exp_x, y_data_fit_exp, label='$%5.3f e^{- %5.3f \\Upsilon} $' % (parameter_exp[0], parameter_exp[1]) , color = 'orange')
        #plt.plot(x, y_data_fit_gauss_all, label='$ fit: \\frac{1}{\sqrt{2 \pi  %5.3f}} e^{- \\frac{(x - %5.3f )^2}{%5.3f}} $' % (1/(parameter_gauss_all[0]*2) ,parameter_gauss_all[1] ,1/(parameter_gauss_all[0]*2)) , color = 'k')

        plt.ylim(0.01,15)
        plt.legend()
        #plt.yscale('log') 
        plt.xlim(3.5,4.5)
        plt.xlabel('Cell Shape Parameter')
        plt.ylabel('N (%)')
        plt.title("Cell Shape parameter Distribution - time = {} ".format((frame_number-2)*10000))
        plt.savefig('Cell_shape_dist_{}_{}.png'.format(Filename, frame_number))
        plt.clf()

    return order_parameter_soft, order_parameter_hard , order_parameter_soft+order_parameter_hard


def Cell_Shape_Parameter_distribution_3D(f, frame_number, Filename, size_z=50, size = 100): # author: Chat-GPT - too tired to think
    type = "3D_N_P_Z"

    Shape_index = np.arange(4.8, 10, 0.02)
    Z_position = np.arange(0, size_z+1, 1)
    R_position = np.arange(0, 20, 1)
    grid_N_P_Z = np.zeros((len(Z_position) , len(Shape_index)))
    # make empty grid for 3D plot
    grid_P_Z_R = grid = np.empty((len(Z_position), len(R_position)), dtype=object)
    grid_P_Z_R_mean = np.zeros((len(Z_position), len(R_position)))
    for i in range(len(Z_position)):
        for j in range(len(R_position)):
            grid_P_Z_R[i][j] = []

    S, Z = np.meshgrid(Shape_index , Z_position)
    r,z = np.meshgrid(R_position , Z_position)
    p_z = np.zeros((len(Z_position)))
    p_z_err = np.zeros(len(Z_position)) 
    p_func_z = [[] for _ in range(len(Z_position))]

    order_parameter= []

    # Ensure 'th.cellInd' is properly defined before using this
    m = np.vstack(th.cellInd)  # Ensure 'th' is defined
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0)[0]

    # Process soft cells
    h = [f[i][:180] for i in s1]
    for cell in h:
        hull = ConvexHull(cell[:, 0:3])  # Switching to 3D
        ord_par = (hull.area)  / hull.volume ** (2/3)
        order_parameter.append(ord_par)

        p_func_z[int(np.mean(cell[:,2]))].append(ord_par)

        # Clip indices to prevent out-of-bounds errors
        shape_idx = np.searchsorted(Shape_index, ord_par)
        z_idx = np.searchsorted(Z_position, np.mean(cell[:,2]))
        com = np.mean(cell, axis=0)
        R_index = np.searchsorted(R_position, np.linalg.norm(com[0:2] - size/2))
        grid_P_Z_R[z_idx][int(R_index)].append(hull.volume)#append(ord_par)
        grid_N_P_Z[ z_idx , shape_idx] += 1

    # Process hard cells
    g = [f[i][:180] for i in s2]

    for cell in g:
        hull = ConvexHull(cell[:, 0:3])  # Switching to 3D
        ord_par = (hull.area)  / hull.volume ** (2/3)
        order_parameter.append(ord_par)

        p_func_z[int(np.mean(cell[:,2]))].append(ord_par)

        # Clip indices to prevent out-of-bounds errors
        shape_idx = np.searchsorted(Shape_index, ord_par)
        z_idx = np.searchsorted(Z_position, np.mean(cell[:,2]))
        com = np.mean(cell, axis=0)
        R_index = np.searchsorted(R_position, np.linalg.norm(com[0:2] - size/2))
        grid_P_Z_R[z_idx][int(R_index)].append(hull.volume)#append(ord_par)
        grid_N_P_Z[ z_idx , shape_idx] += 1


    if type == "3D_N_P_Z":
        section = int(len(Shape_index) * 1/5)  # Define how much of Shape_index you want to keep

        # Plot 3D surface
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        print("Min shape index:", np.min(order_parameter))
        print("Max shape index:", np.max(order_parameter))

        # surf = ax.plot_surface(Z, S[:int(len(Shape_index)*2/3)], grid[:][:int(len(Shape_index)*2/3)], cmap=cm.coolwarm, edgecolor='royalblue', lw=0.5, 
                            # rstride=1, cstride=1, alpha=0.3, antialiased=False)

        # Slice all Z values, but limit Shape_index range
        surf = ax.plot_surface(Z[:, :section], S[:, :section], grid[:, :section], cmap=cm.coolwarm, edgecolor='royalblue', lw=0.5, rstride=1, cstride=1, alpha=0.3, antialiased=False)

        ax.set_ylabel(r'$P^{3D} = \frac{S}{V^{2/3}}$')
        ax.set_xlabel('Z')
        ax.set_zlabel('N')
        ax.set_title(f"3D Shape Index Distribution \n frame = {frame_number}")

        # Set an initial view angle (optional)
        ax.view_init(elev=30, azim=45)  # Change elevation and azimuthal angles

        # Save figure
        plt.savefig(f'3D_shape_index_distribution_{Filename}_{frame_number}.png')

        # Enable interactive mode
        plt.show()

    if type == "2D_P_Z":
        for i in range(len(p_func_z)):
            p_z[i] = np.mean(p_func_z[i])  # Mean shape index at each Z position
            p_z_err[i] = np.std(p_func_z[i]) / np.sqrt(len(p_func_z[i]))  # Standard error (SE)

        # Plot 2D with error bars
        fig = plt.figure()
        ax = fig.add_subplot()

        ax.errorbar(Z_position, p_z, yerr=p_z_err, fmt='o', capsize=5, color="royalblue")
        ax.set_ylabel(r'$\bar{P}^{3D}$')
        ax.set_xlabel('Z')
        ax.set_title(f"Mean Shape Index vs Z \nframe = {frame_number}")

        plt.savefig(f'2D_shape_index_distribution_{Filename}_{frame_number}.png')

    if type == "3D_P_R_Z":
        for i in range(len(Z_position)):
            for j in range(len(R_position)):
                grid_P_Z_R_mean[i][j] = np.mean(grid_P_Z_R[i][j])

        # section = int(len(Shape_index) * 1/5)  # Define how much of Shape_index you want to keep

        # Plot 3D surface
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        print("Min shape index:", np.min(order_parameter))
        print("Max shape index:", np.max(order_parameter))

        # surf = ax.plot_surface(Z, S[:int(len(Shape_index)*2/3)], grid[:][:int(len(Shape_index)*2/3)], cmap=cm.coolwarm, edgecolor='royalblue', lw=0.5, 
                            # rstride=1, cstride=1, alpha=0.3, antialiased=False)

        # Slice all Z values, but limit Shape_index range
        surf = ax.plot_surface(z,r, grid_P_Z_R_mean, cmap=cm.coolwarm, edgecolor='royalblue', lw=0.5, rstride=1, cstride=1, alpha=0.3, antialiased=False)

        ax.set_ylabel("R")
        ax.set_xlabel('Z')
        ax.set_zlabel('V')
        ax.set_title(f"Volume Distribution \n frame = {frame_number}")#(f"3D Shape Index Distribution \n frame = {frame_number}")

        # Set an initial view angle (optional)
        ax.view_init(elev=30, azim=45)  # Change elevation and azimuthal angles

        # Save figure
        plt.savefig(f'Volume_3D_shape_index_distribution_P_R_Z_{Filename}_{frame_number}.png')

        # Enable interactive mode
        


    return order_parameter


def Animate_cell_packing():
    start_from = 0
    end = 0
    view = 'YZ'

    def soft_hard_animation(f,th,frame_number , Filename , size , ax , view = 'XY' , z = 1.9):
        m = np.vstack(th.cellInd)
        c = [ f[i][:180] for i in range(len(m)) ]        
        #Elong
        
        for index , cell in enumerate(c):
            hull = ConvexHull(cell[:,0:2])
            p = cell[hull.vertices,0:2]
            
            if view == 'XY':
                pass
            elif view == 'YZ':
                hull = ConvexHull(cell[:,1:3])
                p = cell[hull.vertices,1:3]

            patches = []
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches,alpha=0.8)
            ax.add_collection(p)

        if view == 'XY':
            ax.set_xlim([0,size])
            ax.set_ylim([0,size])
        elif view == 'YZ':
            ax.set_xlim([0,size])
            ax.set_ylim([0,z])
            
        ax.set_title("Cell Distribution, Time = {}".format(frame_number*1000))
        print(frame_number)
        return ax.plot()
    
    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return soft_hard_animation(frame, th, frame_number , filename, 300, ax , 'YZ' , 100) #### set the size

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            if view == 'XY':
                fig, ax = plt.subplots()
            elif view == 'YZ':
                fig, ax = plt.subplots(figsize=(21,7)) ############################################################# set the figsize for ephitelial tissue

                    
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            # cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",['#2A9B2A','r'])
            # cnorm = mcol.Normalize(vmin=0,vmax=1)
            # cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
            # cpick.set_array([])
            
            #create an animation
            for _ in range(start_from):
                th.ReadFrame(inc=nSkip)

            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(start_from,th.maxFrames-end,nSkip), interval=50, repeat=False)
            # cbar = plt.colorbar(cpick, ax=ax , pad =0.1)
            # cbar.set_label("$\\Upsilon$", labelpad=10)  # Main label (on the side)
            # cbar.ax.set_title('elongated', pad=10, fontsize=12)  # Label above the color bar
            # cbar.ax.set_xlabel('Spherical', labelpad=10 , fontsize=12)  # Label below the color bar
            # plt.show()
            #save animation as gif
            if view == 'XY':
                ani.save('cell_distribution_XY_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)
            elif view == 'YZ':
                ani.save('cell_distribution_YZ_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)

    # def soft_hard_animation(f,th,frame_number , Filename , size , ax):
    #     Area_soft = []
    #     Area_hard = []
    #     '''cell node coordinates soft (h) - hard(g)'''
    #     m = np.vstack(th.cellInd)
    #     s1 = np.where(m < 0)[0]
    #     s2 = np.where(m >= 0 )[0]

    #     h = [ f[i][:180] for i in s1 ]
    #     for cell in h:
    #         patches=[]
    #         hull = ConvexHull(cell[:,0:2])
    #         p = cell[hull.vertices,0:2]
    #         patches.append(Polygon(p, closed=True))
    #         p = PatchCollection(patches, color = (0.2, 0.6, 1) ,edgecolor="b",alpha=0.8)
    #         ax.add_collection(p)
    #         Area_soft.append(hull.volume)

    #     g = [ f[i][:180] for i in s2 ]
    #     for cell in g:
    #         patches=[]
    #         hull = ConvexHull(cell[:,0:2])
    #         p = cell[hull.vertices,0:2]
    #         patches.append(Polygon(p, closed=True))
    #         p = PatchCollection(patches, color= (222/255, 28/255, 0 ) ,alpha=0.8)
    #         ax.add_collection(p)
    #         Area_hard.append(hull.volume)

    #     ax.set_xlim([0,size])
    #     ax.set_ylim([0,size])
    #     ax.set_title("DC = {}% , Time = {}".format(int( (sum(Area_soft) + sum(Area_hard) )/(size*size)*100) , (frame_number)*10000))
    #     # collection = ax.plot()
    #     # artist.append(collection)
    #     return ax.plot()

def Animate_cell_shape_index(view='XY', size_z=12, size =200):
    start_from =0
    end = 0
    size=200

    def soft_hard_animation(f, th, frame_number, Filename, size, ax, view='XY', z=1.9):
        m = np.vstack(th.cellInd)
        s1 = np.where(m < 0)[0]
        s2 = np.where(m >= 0)[0]

        Order_param_sphere = 3.75
        Area_soft, Area_hard = [], []
        order_parameter_soft, order_parameter_hard = [], []

        def process_cells(cells, is_soft):
            patches_list = []
            order_params = []
            for cell in cells:
                if view == 'XY':
                    hull = ConvexHull(cell[:, 0:2])
                    p = cell[hull.vertices, 0:2]
                elif view == 'YZ':
                    hull = ConvexHull(cell[:, 1:3])
                    p = cell[hull.vertices, 1:3]
                elif view == 'XZ':
                    hull = ConvexHull(cell[:, ::2])
                    p = cell[hull.vertices, ::2]
                
                patches = []
                patches.append(Polygon(p, closed=True))
                ord_par = hull.area / np.sqrt(hull.volume)
                order_params.append(ord_par)

                if is_soft:
                    color = (218 / 255, 112 / 255, 214 / 255) if ord_par > Order_param_sphere else (0.2, 0.6, 1)
                    patch = PatchCollection(patches, color=color, edgecolor="b", alpha=0.8)
                else:
                    color = (1, 128 / 255, 0) if ord_par > Order_param_sphere else (222 / 255, 28 / 255, 0)
                    patch = PatchCollection(patches, color=color, edgecolor="r", alpha=0.3)

                ax.add_collection(patch)
                patches_list.append(patch)

            return order_params

        h = [f[i][:180] for i in s1]
        order_parameter_soft = process_cells(h, is_soft=True)

        g = [f[i][:180] for i in s2]
        order_parameter_hard = process_cells(g, is_soft=False)

        if view == 'XY':
            ax.set_xlim([0, size])
            ax.set_ylim([0, size])
        elif view == 'YZ':
            ax.set_xlim([0, size])
            ax.set_ylim([0, z])
        elif view == 'XZ':
            ax.set_xlim([0, size])
            ax.set_ylim([0, z])

        ax.set_title(
            "DC = {}% , Time = {}".format(
                int((sum(Area_soft) + sum(Area_hard)) / (size * size) * 100),
                frame_number * 10000,
            )
        )
        #save figure
        return ax.plot()

    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return soft_hard_animation(frame, th, frame_number, filename, size, ax, view, size_z)

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            if view == 'XY':
                fig, ax = plt.subplots(figsize=(size / 5, size / 5))
            else:  # view == 'YZ' or 'XZ'
                fig, ax = plt.subplots(figsize=(size / 5, size_z / 5))

            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            # Create animation
            for _ in range(start_from):
                th.ReadFrame(inc=nSkip)
            ani = animation.FuncAnimation(
                fig=fig,
                func=update,
                frames=range(start_from, th.maxFrames - end, nSkip),
                interval=50,
                repeat=False,
            )
            # Save animation as GIF
            if view == "XY":
                ani.save(
                    "shape_index_XY_{}.gif".format(filename[:-4]), writer="Pillow", fps=5
                )
            elif view == "YZ":
                ani.save(
                    "shape_index_YZ_{}.gif".format(filename[:-4]), writer="Pillow", fps=5
                )
            elif view == "XZ":
                ani.save(
                    "shape_index_XZ_{}.gif".format(filename[:-4]), writer="Pillow", fps=5
                )
            print("The view is: {}".format(view))

